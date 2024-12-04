#'  SCRIPT for notes around image manipulation in R

#' packages to consider
#'
#' imager
#' spatstat
#' magick
#' sf
#' sp
#' RImageJROI

# install and load required packages ####
pkg_req <- c("sf", "RImageJROI", "ggplot2", "dplyr", "pracma", "minpack.lm", "magrittr", "stringr", "circular", "svglite", "astsa", "pdftools", "scico", "ggpubr", "rstatix", "matrixStats", "gridGraphics")

# Install and load packages if not already installed
for (package in pkg_req) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# FUNCTIONS ####
#' function to import a single ROI
importROI = function(path){
  roi_path = stringr::str_replace_all(path, "\\\\", "//")
  roi_data <- RImageJROI::read.ijroi(roi_path)
  # convert to sf package compatible
  roi_coords <- roi_data$coords
  if (!all(roi_coords[1, ] == roi_coords[nrow(roi_coords), ])) {
    # Add the first coordinate at the end to close the polygon
    roi_coords <- rbind(roi_coords, roi_coords[1, ])
  }
  roi_polygon <- sf::st_polygon(list(roi_coords))
  # roi_sf <- sf::st_sfc(roi_polygon)
  return(roi_polygon)
}

#' function to import multiple ROIs
importROIs = function(path){
  rois_path = stringr::str_replace_all(path, "\\\\", "//")
  # use RImageJROI::read.ijroi to read through all the zipped rois
  tau_parts_ROI_list = read.ijzip.progr(rois_path, names = TRUE)
  print("Processing particle ROIs...")
  pb <- txtProgressBar(min = 0, max = length(tau_parts_ROI_list), style = 3)
  tau_polygons <- setNames(
    lapply(seq_along(tau_parts_ROI_list), function(i){
      coord_vec <- tau_parts_ROI_list[[i]]  # Get the current element
      # extract coordinates, contained at position 25 of the list
      coord_vals <- coord_vec[[25]]
      if (!all(coord_vals[1, ] == coord_vals[nrow(coord_vals), ])) {
        # Add the first coordinate at the end to close the polygon
        coord_vals <- rbind(coord_vals, coord_vals[1, ])
      }
      polygon <- sf::st_polygon(list(coord_vals))
      setTxtProgressBar(pb, i)
      flush.console()
      return(polygon)
    }),
    names(tau_parts_ROI_list)
  )
  close(pb)
  return(tau_polygons)
}

#' modified version of the function from the RImageJROI pkg to include progress bar
read.ijzip.progr = function(file, names = TRUE, list.files = FALSE, verbose = FALSE) {
  
  ## Read files in the zip file
  files <- unzip(file, list = TRUE)
  
  if(any(sapply(strsplit(files$Name, "\\."), '[', 2) == "roi") == FALSE) stop("The zip file contains other files than ImageJ ROIs")
  
  if(list.files == FALSE){
    ## Find a suitable location to unizip the zip file
    location <- tempfile()
    
    ## Unzip the zip file to a temporary folder
    unzip(file, exdir = location)
    
    # Read ROIs
    print("Reading ROIS...")
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    roi.dat <- sapply(seq_along(files$Name), function(i) {
      tmp <- read.ijroi(paste(location, files$Name, sep = "/")[i],
                        verbose = verbose)
      tmp2 <- list(tmp)
      names(tmp2) <- tmp[["name"]]
      setTxtProgressBar(pb, i)
      flush.console()
      return(tmp2)
    })
    close(pb)
    print("Completed")
    ## Remove the temporary folder
    unlink(location, recursive = TRUE)
    
    ## Rename elements of the returned list
    if (names == FALSE){
      rep.names <- unlist(lapply(seq_along(roi.dat), function(i) roi.dat[[i]]$strType))
      rep.names <- make.unique(rep.names)
      rep.numbers <- sapply(strsplit(rep.names, "\\."), '[', 2)
      rep.numbers[is.na(rep.numbers)] <- 0
      rep.names <- paste(sapply(strsplit(rep.names, "\\."), '[', 1), as.character(as.numeric(rep.numbers)+1), sep = "_")
      names(roi.dat) <- rep.names
    }
    class(roi.dat) <- "ijzip"
    return(roi.dat)}
  
  if (list.files == TRUE) return(files)
}

#' Import all values from kernel extraction
importGridVals = function(path){
  path_fix = stringr::str_replace_all(path, "\\\\", "//")
  #read csv and transpose for more memory efficient handling
  matrix = read.csv(path_fix) %>% t()
  #Formatting steps
  colnames(matrix) = seq(1:dim(matrix)[2])
  matrix = matrix[-1,]
  rownames(matrix) = stringr::str_remove(rownames(matrix), "Mean")
  return(matrix)
}

#' compute matrix of all distances from cells to particles
cell2particlesDistances = function(polygon_frames_list, cells_points, frame, ...){
  # select all polygons from a defined frame
  polygons_list <- polygon_frames_list[[frame]]
  
  # combine polygons from a frame into a single object
  polygons_sf <- do.call(rbind, lapply(seq_along(polygons_list), function(i) {
    polygon_sf <- st_sf(geometry = st_sfc(polygons_list[[i]]))
    polygon_sf$id <- names(polygons_list[i])  # Assign a unique ID to each polygon
    polygon_sf
  }))
  
  #' compute distances between all cells and polygons
  distance_matrix <- sf::st_distance(cells_points, polygons_sf)
  
  #' calculate minimum distance
  min_distances <- apply(distance_matrix, 1, min)
  
  #' identify the closest polygon
  closest_polygons <- apply(distance_matrix, 1, function(row) {
    polygons_sf$id[which.min(row)]
  })
  
  #' Add results to the points data frame
  cells_points$min_distance <- min_distances
  cells_points$closest_polygon_id <- closest_polygons
  cells_points$frame <- frame
  
  return(cells_points)
}

#' create plot of cells that are far/close
highlight_cells <- function(evaluation_matrices, filename, colors, all_cells = grid_coord, scn_cells = grid_points_scn, scn_col = "grey", outside_col = "black") {
  #' evaluation_matrices: list of matrices where each element represents a matrix for different conditions
  #' frame: the frame to evaluate
  #' colors: vector of colors corresponding to each matrix in `evaluation_matrices`
  #' all_cells: data frame with X and Y coordinates for all grid points
  #' scn_cells: spatial object or data frame with SCN cell coordinates
  #' scn_col: color for SCN cells
  #' outside_col: color for cells outside any specified condition
  
  # Copy `all_cells` as the base data frame
  plot_data <- all_cells
  
  # Extract coordinates from `scn_cells`
  scn_coords <- st_coordinates(scn_cells)
  
  # Add `is_SCN` column to identify SCN coordinates
  plot_data$is_SCN <- with(plot_data, paste(X, Y) %in% paste(scn_coords[, "X"], scn_coords[, "Y"]))
  
  # Initialize `color` column in `plot_data` with `outside_col` as default
  plot_data$color <- outside_col
  plot_data$alphaval <- 1
  
  # Override color for SCN cells, if needed
  plot_data$color <- ifelse(plot_data$is_SCN, scn_col, outside_col)
  
  # Loop through each matrix and corresponding color
  for (i in seq_along(evaluation_matrices)) {
    active_pixels <- evaluation_matrices[[i]]#[, frame]
    active_coords <- all_cells[as.integer(names(which(active_pixels))), ]  # Selects rows of `grid_coord` corresponding to TRUE in `evaluation_matrices`
    
    # Update color and alpha for cells that match active pixels in the current matrix
    plot_data$color <- ifelse(paste(plot_data$X, plot_data$Y) %in% paste(active_coords$X, active_coords$Y), colors[i], plot_data$color)
  }
  
  plot_data$alphaval <- ifelse(plot_data$color == "black", 0, 1)
  # Create the plot
  plot <- ggplot(plot_data, aes(x = X, y = Y, color = color, alpha = alphaval)) +
    geom_point() +
    scale_color_identity() +
    scale_alpha_identity()+
    theme_minimal() +
    labs(title = paste(filename, "- Distance groups"), x = "X Coordinate", y = "Y Coordinate") +
    coord_fixed() +
    scale_y_reverse()+
    theme(plot.title = element_text(size = 11),
          axis.title = element_text(size = 9))
  
  return(plot)
}

#' plot features of period on map
highlight_cells_period = function(all_cells = grid_coord, period_table, variable, filename = "", colorcode = "blue", sdFactor = 1){
  
  plot_data <- grid_coord %>% `colnames<-`(c("ID", "X", "Y"))
  # Merge tables to include period data
  merged_data <- left_join(plot_data, period_table, by = "ID")
  
  # Choose the variable to visualize (e.g., period, amplitude, error)
  variable <- variable
  
  mean_value <- mean(merged_data[[variable]], na.rm = TRUE)
  sd_value <- sd(merged_data[[variable]], na.rm = TRUE)
  lower_bound <- mean_value - (sdFactor * sd_value)
  upper_bound <- mean_value + (sdFactor * sd_value)
  lower_bound <- ifelse(lower_bound < 0, 0, lower_bound)
  
  colorschemes <- list("blue" = c("#DEEBF5", "#A0C9DF", "#559ECA", "#1663A5", "#003066"),
                       "green" = c("#F1F9EE", "#C0E5BB", "#6EBD71", "#07421F", "#07421F"),
                       "orange" = c("#FFEFDB", "#FCB37E", "#F5734C", "#D12B22", "#7A0006"),
                       "purple" = c("#FDEAE6", "#FAA5B0", "#F45F96", "#B21978", "#48065F"),
                       "red" = c("#FFEDE5", "#FAB29A", "#FD5C47", "#C81B21", "#640112"),
                       "violet" = c("#FAF9FC", "#B7BAD8", "#8E8CBD", "#685DA2", "#371071"),
                       "acqua" = c("#F0F9E9", "#B9E2BA", "#6EC1BD", "#217EB2", "#003B74"))
  colorscheme = colorschemes[[colorcode]]
  # Add a new variable to handle outliers and scaling
  merged_data <- merged_data %>%
    mutate(
      alphaval = ifelse(is.na(!!sym(variable)), 0, 1)
    )
  
  plot <- ggplot(merged_data, aes(x = X, y = Y, alpha = alphaval)) +
    geom_point(aes(color = !!sym(variable))) +
    scale_color_gradientn(
      colors = colorscheme,
      values = scales::rescale(c(lower_bound, mean_value, upper_bound)),
      limits = c(lower_bound, upper_bound),
      oob = scales::squish  # Ensures outliers are squished into the limits
    ) +
    scale_alpha_identity()+
    theme_minimal() +
    labs(
      title = paste0(filename, " - Distribution of ", variable),
      x = "X Coordinate",
      y = "Y Coordinate"
    ) +
    coord_fixed() +
    scale_y_reverse()+
    theme(plot.title = element_text(size = 11),
          axis.title = element_text(size = 9))
  
  return(plot)
}

#' plot features of period on map
highlight_cells_period_circ = function(all_cells = grid_coord, period_table, variable, filename = ""){
  
  plot_data <- grid_coord %>% `colnames<-`(c("ID", "X", "Y"))
  # Merge tables to include period data
  merged_data <- left_join(plot_data, period_table, by = "ID")
  
  colorscheme = c("#04305C", "#245380", "#4376A4", "#8BB09D", "#D2EA96", "#E8F4A8", "#FEFEB9", "#FEDD96", "#FDBC73", "#E08258", "#C2473D", "#8E3334", "#591F2B")
  
  merged_data <- merged_data %>%
    mutate(
      alphaval = ifelse(is.na(!!sym(variable)), 0, 1)
    )
  # Choose the variable to visualize (e.g., period, amplitude, error)
  variable <- variable
  
  # Create the plot
  plot <- ggplot(merged_data, aes(x = X, y = Y, alpha = alphaval)) +
    geom_point(aes(color = !!sym(variable))) +
    theme_minimal() +
    scale_alpha_identity()+
    scale_color_gradientn(colors = colorscheme) +
    labs(title = paste0(filename, " - Distribution of ", variable), x = "X Coordinate", y = "Y Coordinate") +
    coord_fixed() +
    scale_y_reverse()+
    theme(plot.title = element_text(size = 11),
          axis.title = element_text(size = 9))
  
  return(plot)
}

# fix path
back_to_forw <- function(x){
  stringr::str_replace_all(x, "\\\\", "//")
}

# prepare data into long table format
prep_table = function(data){
  # re-transpose table
  data = t(data)
  # add time column
  t = seq(0, length.out = length(data[, 1]), by = 0.5)
  data = cbind(t, data)
  
  # perform outliers removal, smoothening and detrending
  data_clean <- preprocess_data(data, grade = 3)
  
  #make it into long table format
  data_long <- tidyr::pivot_longer(data_clean,
                                   cols = colnames(data_clean)[-1],
                                   names_to = "ID") %>% data.table::as.data.table(key = "ID")
  
  return(data_long)
}

# function to remove outliers, smooth and detrend
preprocess_data = function(data, grade){
  # load table
  data_to_clear <- data
  # create table to store smoothened data
  clear_data <- data[,1]
  # get time values
  x_vals <- data[,1]
  print("Preprocessing data...")
  pb <- txtProgressBar(min = 2, max = length(colnames(data_to_clear)), style = 3)
  # for cycle to go through each column and smooth the data
  for (i in seq(2, length(colnames(data_to_clear)))){
    timeserie <- data_to_clear[,i]
    # create vector for cleaned data
    cleaned_ts = timeserie
    # perform first loess approximation
    smoothed_series <- loess(timeserie ~ seq_along(timeserie), span = 0.08)
    smoothed_values <- predict(smoothed_series, newdata = data.frame(x = x_vals))
    # calculate outliers based on distance from loess curve and sd
    stdev_ts <- sd(smoothed_values)
    outliers <- which(abs(smoothed_values - timeserie) > stdev_ts*0.6)
    # delete outliers from data
    cleaned_ts[outliers] <- NA
    # perform second loess approximation on clean data
    smoothed_clean_series <- loess(cleaned_ts ~ seq_along(cleaned_ts), span = 0.08)
    smoothed_clean_values <- predict(smoothed_clean_series, newdata = data.frame(x = x_vals))
    # interpolate missing values before detrending
    smoothed_clean_values <- imputeTS::na_interpolation(smoothed_clean_values, option = "spline")
    # add detrending step
    smoothed_clean_detr_values <- astsa::detrend(smoothed_clean_values, grade)
    # add smoothened data to table
    clear_data <- cbind(clear_data, smoothed_clean_detr_values)
    # update progressbar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  print("Completed")
  clear_data <- as.data.frame(clear_data)
  # re add colnames
  colnames(clear_data) <- colnames(data)
  # return new table
  return(clear_data)
}

# perform period analysis
FFT_NLLS_analyse = function(data, minPer = 15, maxPer = 32){
  # Perform FFT
  n <- length(data$value)
  dt <- mean(diff(data$t))
  fft_result <- stats::fft(data$value)
  frequencies <- seq(0, 1/dt, length.out = n)
  
  # Identify the dominant frequency (excluding the zero frequency)
  dominant_frequency <- frequencies[which.max(base::Mod(fft_result)[2:(n/2)]) + 1]
  initial_period <- 1 / dominant_frequency
  
  # Define the sinusoidal model function
  sinusoidal_model = function(params, t) {
    amplitude <- params[1]
    period <- params[2]
    phase <- params[3]
    offset <- params[4]
    return(amplitude * sin(2 * pi * t / period + phase) + offset)
  }
  
  # Define the residual function for nonlinear least squares
  residuals = function(params, t, y) {
    return(y - sinusoidal_model(params, t))
  }
  
  # Initial parameter estimates
  initial_amplitude <- (max(data$value) - min(data$value)) / 2
  initial_phase <- 0
  initial_offset <- mean(data$value)
  initial_params <- c(initial_amplitude, initial_period, initial_phase, initial_offset)
  
  # Perform nonlinear least squares fitting
  fit <- minpack.lm::nls.lm(
    par = initial_params,
    fn = residuals,
    t = data$t,
    y = data$value,
    lower = c(-Inf, minPer, -Inf, -Inf),
    upper = c(Inf, maxPer, Inf, Inf)
  )
  
  # Extract fitted parameters
  fitted_params <- fit$par
  fitted_amplitude <- fitted_params[1]
  fitted_period <- fitted_params[2]
  fitted_phase <- fitted_params[3]
  fitted_offset <- fitted_params[4]
  
  # Ensure the amplitude is positive
  if (fitted_amplitude < 0) {
    fitted_amplitude <- abs(fitted_amplitude)
    fitted_params[1] <- fitted_amplitude
    fitted_phase <- fitted_phase + pi  # Adjust the phase by 180 degrees (π radians)
    fitted_params[3] <- fitted_phase
  }
  
  # Wrap the phase into the 0-2pi radians range
  if (fitted_phase > 2*pi) {
    fitted_phase <- fitted_phase %% 2*pi
  }
  if (fitted_phase < 0) {
    fitted_phase <- fitted_phase + 2*pi
  }
  
  # Convert phase to time units (hours)
  fitted_phase <- circular::as.circular(fitted_phase, type = "angles", units = "radians", rotation = "clock")
  phase_circadian <- circular::conversion.circular(fitted_phase, units = "hours")
  phase_absolute <- phase_circadian*(fitted_period/24)
  
  # Compute residual error
  fitted_values <- sinusoidal_model(fitted_params, data$t)
  residual_error <- sqrt(mean((data$value - fitted_values)^2))
  
  # Calculate R-squared (GOF)
  ss_total <- sum((data$value - mean(data$value))^2)
  ss_res <- sum((data$value - fitted_values)^2)
  r_squared <- 1 - (ss_res / ss_total)
  
  # Return the results
  result <- list(
    keep = TRUE,
    period = fitted_period,
    amplitude = fitted_amplitude,
    phase_h = phase_absolute,
    phase_rad = fitted_phase,
    phase_circ = phase_circadian,
    offset = fitted_offset,
    error = residual_error,
    GOF = r_squared
  )
  return(result)
}

# create circular plot
circular_plot = function(data, path, saving = TRUE, ...){
  phases_rad = data
  if(!all(is.na(phases_rad))){
    real_MeanPh <- circular::mean.circular(phases_rad, Rotation = "counter", na.rm = TRUE)
    phases_rad = align_phase(phases_rad, CT = 6)
    phases = circular::circular(phases_rad, units = "radians", template = "clock24")
  }
  if(saving){
    svgdatapath = path
    svg(filename = svgdatapath, width = 5, height = 4.8, family = "sans", bg = "transparent")
  }
  # if there are no data to plot
  if(is.infinite(max(phases_rad, na.rm = TRUE))){
    phases = NA
    plot.new()
    circular::plot.circular(phases, bins = 24, axes = FALSE)
    circular::axis.circular(circular(seq(0, 2*pi, pi/2)), zero=0, rotation = 'counter', labels=c("6", "00", "18", "12", "6"), tcl.text=0.12, cex = 1.5)
    title(paste(filename), line = 0)
    dev.off()
    variance = NA
    vectorLength = NA
    RT = NA
  }else{
    #rotation = pi/2 - mean(phases)
    plot.new()
    # Create a circular plot
    circular::plot.circular(phases, bins = 24, axes = FALSE)#, col = 'black', cex = 0, pch = 0, ticks = FALSE, tcl = 0)
    # Add points to the existing plot
    circular::points.circular(phases, cex = 0.75, pch = 21, bg = '#eb0000')#, sep = 0.01)
    # Add radial axis labels to the plot
    circular::axis.circular(circular(seq(0, 2*pi, pi/2)), zero=0, rotation = 'counter', labels=c("6", "00", "18", "12", "6"), tcl.text=0.12, cex = 1.5)
    # Calculate the mean circular statistic
    meanTest <- circular::mean.circular(phases, Rotation = "counter", na.rm = TRUE)
    # Perform a Rayleigh test on RCaMP
    RT <- rayleigh.test(phases)
    # Calculate the vector length statistic
    vectorLength <- RT$statistic[1]
    # Get the p-value from the Rayleigh test
    circStats <- RT$p.value[1]
    # Round the p-value to 3 decimal places
    circStats <- round(circStats, digits = 3)
    # Add arrows to the plot based on the mean and vector length
    circular::arrows.circular(meanTest, vectorLength, length = 0.1, col = '#eb0000', lwd = 2)
    # Add an empty title (the paste("") effectively removes the default title)
    title(paste(filename), line = 0)
    #calculate angular variance
    variance = circular::angular.variance(phases, na.rm = TRUE)
    # close plot
    dev.off()
  }
  
  #capture plot and return it
  # plot_capture <- recordPlot()
  list_returns <- list(#plot = plot,
    vec_len = vectorLength,
    variance = variance,
    meanP = real_MeanPh
    #RT = RT
  )
  
  return(list_returns)
}

circStats_plot = function(data, path, colors){
  data = data
  
  svgdatapath = path
  svg(filename = svgdatapath, width = 5, height = 4.8, family = "sans", bg = "transparent")
  
  plot.new()
  # Create a circular plot
  circular::plot.circular(NA, bins = 24, axes = FALSE)
  # Add radial axis labels to the plot
  circular::axis.circular(circular(seq(0, 2*pi, pi/2)), zero=0, rotation = 'counter', labels=c("6", "00", "18", "12", "6"), tcl.text=0.12, cex = 1.5)
  # Add arrows to the plot based on the mean and vector length
  # close arrow
  vectorLength = data["close_vec_len"][[1]]
  meanTest = data["close_meanPh"][[1]]
  variance = data["close_variance"][[1]]
  color = colors[1]
  circular::arrows.circular(meanTest, vectorLength, length = 0.1, col = color, lwd = 3)
  # mid arrow
  vectorLength = data["mid_vec_len"][[1]]
  meanTest = data["mid_meanPh"][[1]]
  color = colors[2]
  circular::arrows.circular(meanTest, vectorLength, length = 0.1, col = color, lwd = 3)
  # far arrow
  vectorLength = data["far_vec_len"][[1]]
  meanTest = data["far_meanPh"][[1]]
  color = colors[3]
  circular::arrows.circular(meanTest, vectorLength, length = 0.1, col = color, lwd = 3)
  
  # Add an empty title (the paste("") effectively removes the default title)
  title(paste(filename), line = 0)
  # close plot
  dev.off()
}

# print period lengths
period_print = function(data, ...){
  df <- data
  filename
  cells = "Cells"
  
  min_y = round((min(data$period, na.rm = TRUE)-1), 0)
  max_y = round((max(data$period, na.rm = TRUE)+1), 0)
  if(is.infinite(min_y)){min_y = 20}
  if(is.infinite(max_y)){max_y = 28}
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(cells, period)) +
    ggplot2::geom_jitter(alpha=.5)+
    labs(title = filename,
         y = "Period length (h)") +                      # Add axis labels and a title
    ggpubr::theme_pubclean()+
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # Center and style the title
      axis.title.x = element_text(size = 18, face = "bold"),            # Style x-axis title
      axis.title.y = element_text(size = 18, face = "bold"),            # Style y-axis title
      axis.text.x = element_text(size = 0),                            # Style x-axis text
      axis.text.y = element_text(size = 16),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 1),
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent")
    )+
    scale_y_continuous(limits = c(min_y, max_y),
                       expand = c(0, 0),
                       breaks = seq(from = min_y, to = max_y, by = 2))
  return(plot1)
}

# function to align phases to a standard
align_phase = function(data, CT = 6, ...){
  phases = data
  mean_phase = mean(phases, na.rm = TRUE)
  CT = CT
  radCT = (CT/12)*pi
  difference = mean_phase - radCT
  if(radCT >= mean_phase){
    phases_aligned = phases + difference
  }else{
    phases_aligned = phases - difference
  }
  return(phases_aligned)
}

# foo to analyse TS and make circular plot
computePeriod = function(df, excludeNC = FALSE, top = 30, bottom = 18, ...){
  #prepare table
  data_df <- prep_table(df)
  
  # analyse period and store result in list
  unique_ids = unique(data_df$ID)
  results_list <- list()
  plot_list <- list()
  print("Computing period...")
  pb <- txtProgressBar(min = 2, max = length(unique_ids), style = 3)
  for (i in seq_len(length(unique_ids))) {
    # Filter data for the current ID
    ID = unique_ids[i]
    toget = which(data_df$ID == ID)
    data_ID <- data_df[toget,]
    
    # perform linear detrending
    fit <- lm(data_ID$value ~ data_ID$t)
    detrended_data <- residuals(fit)
    
    if(length(detrended_data == length(data_ID$values))){
      data_ID$value <- detrended_data}
    
    # Estimate period and other parameters
    result <- FFT_NLLS_analyse(data_ID, minPer = 11, maxPer = 32)
    
    # Add the ID to the result
    result$ID <- as.integer(ID)
    
    # Append the result to the results list
    results_list[[ID]] <- result
    setTxtProgressBar(pb, i)
  }
  close(pb)
  print("Completed")
  period_table = do.call(rbind, lapply(results_list, as.data.frame))
  # round all digits in the table
  period_table <- period_table %>%
    dplyr::mutate(dplyr::across(.cols = -all_of("ID"), ~ round(., 2)))
  # count total traces
  period_counts <- sum(!is.na(period_table$phase_h))
  # clean table of all values outside of limits and count valid traces
  if(excludeNC == TRUE){
    topvals <- which(period_table$period >= top)
    bottomvals <- which(period_table$period <= bottom)
    exclude = c(topvals, bottomvals)
    period_table[exclude,] <- NA
    
    print(paste("Values >", top, "h | <", bottom, "h removed", sep = ""))
  }
  non_na_count <- sum(!is.na(period_table$phase_h))
  
  return(period_table)
}

#' summarize period analysis
summarizePeriod = function(period_table){
  period_counts <- length(period_table$phase_h)
  non_na_count <- sum(!is.na(period_table$phase_h))
  period_var = stats::var(period_table$period, na.rm = TRUE)
  period_mean = mean(period_table$period, na.rm = TRUE)
  amplitude_mean = mean(abs(period_table$amplitude), na.rm = TRUE)
  error_mean = mean(period_table$error, na.rm = TRUE)
  
  # return values outside
  outreturn_list <- list(trace_no = non_na_count,
                         trace_tot = period_counts,
                         period_var = period_var,
                         period_mean = period_mean,
                         amplitude_mean = amplitude_mean,
                         error_mean = error_mean)
  
  return(outreturn_list)
}

# create pdf file with all plots
pdf_plots = function(homedir, output_filename = "combined_plots"){
  
  # Set up the output file path
  output_pdf <- file.path(homedir, paste(output_filename, "_", file_suffix, ".pdf", sep = ""))
  
  # Fetch all SVG files from the folder in alphabetical order
  svg_files <- list.files(homedir, pattern = "\\.svg$", full.names = TRUE)
  svg_files <- sort(svg_files)  # Sort files alphabetically
  
  # Set up the PDF document
  pdf(output_pdf, width = 14, height = 8.5, paper = "special", onefile = TRUE, family = "sans")
  
  # Parameters for layout
  plots_per_page <- 8
  rows <- 2
  cols <- 4
  
  # Loop through the SVG files and plot them
  for (i in seq(1, length(svg_files), by = plots_per_page)) {
    # Create a new page in the PDF
    grid::grid.newpage()
    
    # Get the current batch of SVG files
    current_batch <- svg_files[i:min(i + plots_per_page - 1, length(svg_files))]
    
    # Set up layout for multiple plots
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(rows, cols)))
    
    # Plot the SVGs in the grid
    for (j in seq_along(current_batch)) {
      
      # Calculate the row and column positions for this plot
      row_num <- ceiling(j / cols)
      col_num <- (j - 1) %% cols + 1
      
      # Read the SVG and convert it to raster
      svg_raster <- rsvg::rsvg(current_batch[j])
      
      # Push the viewport for the specific position in the layout
      grid::pushViewport(grid::viewport(layout.pos.row = row_num, layout.pos.col = col_num))
      
      # Add the SVG plot to the PDF
      grid::grid.draw(grid::rasterGrob(svg_raster, interpolate = TRUE))
      
      # Pop the viewport after drawing
      grid::popViewport()
    }
  }
  
  # Close the PDF device
  dev.off()
  
  # Return the output PDF path
  return(output_pdf)
}

#' create empty table to accumulate FFT analysis
createFFTtbl = function(){
  #define df to store variance
  variance_vec = vector()
  analysed_names_vec = vector()
  reporter_vec = vector()
  traces_no_vec = vector()
  tot_traces = vector()
  vector_len = vector()
  period_var = vector()
  period_mean = vector()
  amplitude_mean = vector()
  error_mean = vector()
  return(list(variance_vec,
              analysed_names_vec,
              reporter_vec,
              traces_no_vec,
              tot_traces,
              vector_len,
              period_var,
              period_mean,
              amplitude_mean,
              error_mean))
}

#' add values to table for FFT anaysis
updateFFTtbl = function(old_list, new_list, filename){
  # variance_vec <- c(old_list$variance_vec, new_list$variance)
  traces_no_vec <- c(old_list$traces_no_vec, new_list$trace_no)
  tot_traces <- c(old_list$tot_traces, new_list$trace_tot)
  analysed_names_vec <- c(old_list$analysed_names_vec, filename)
  period_var = c(old_list$period_var, new_list$period_var)
  period_mean = c(old_list$period_mean, new_list$period_mean)
  # vector_len = c(old_list$vector_len, new_list$vector_len)
  amplitude_mean = c(old_list$amplitude_mean, new_list$amplitude_mean)
  error_mean = c(old_list$error_mean, new_list$error_mean)
  reporter_vec <- c(old_list$reporter_vec, strsplit(filename, "_")[[1]][3])
  return(list(#variance_vec,
    analysed_names_vec,
    reporter_vec,
    traces_no_vec,
    tot_traces,
    # vector_len,
    period_var,
    period_mean,
    amplitude_mean,
    error_mean))
}

#' create boxplot to visualize differences between the groups, according to distance
metric_boxplot = function(data, metric, grouping_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, ...){
  
  # in referring to variables within the dplyr context, use the curly-curly syntax
  # when using them as part of strings, use the deparse(substitute()) syntax
  # when needed on add_xy_position, use as_label(enquo())
  
  if(norm_test){ # test for normality
    normality = data %>%  group_by({{grouping_var}}) %>% rstatix::shapiro_test({{metric}})
    print(normality)
  }
  # compute anova
  formula = as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(grouping_var))))
  within_var_sym <- rlang::sym(deparse(substitute(grouping_var)))
  res.aov = rstatix::anova_test(data, formula = formula, wid = ID, within = !!within_var_sym)
  # perform pairwise comparison
  pwc = data %>% rstatix::pairwise_t_test(formula, paired = FALSE, p.adjust.method = test)
  # add xy position of pwc
  # xvar <- as_label(enquo(grouping_var))
  pwc <- pwc %>% rstatix::add_xy_position(x = as_label(enquo(grouping_var)), fun = "max", step.increase = 0.8)
  
  # generate plot
  plot <- ggpubr::ggboxplot(data, x = deparse(substitute(grouping_var)), y = deparse(substitute(metric)), color = deparse(substitute(grouping_var))) +
    ggpubr::stat_pvalue_manual(pwc)+
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    ggplot2::labs(
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = rstatix::get_pwc_label(pwc)
    )+
    ggplot2::scale_color_manual(values = plot_colors)
  ggplot2::theme(plot.subtitle = element_text(size = 10, hjust = 0))
  
  return(plot)
}

metric_mean_plot = function(data, metric, grouping_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, ...){
  
  # Compute mean and standard error for each group
  summary_stats <- data %>%
    dplyr::group_by({{grouping_var}}) %>%
    dplyr::summarise(
      mean = mean({{metric}}, na.rm = TRUE),
      se = sd({{metric}}, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  if (norm_test) { # Test for normality
    normality = data %>% group_by({{grouping_var}}) %>% rstatix::shapiro_test({{metric}})
    print(normality)
  }
  
  # Compute ANOVA
  formula = as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(grouping_var))))
  within_var_sym <- rlang::sym(deparse(substitute(grouping_var)))
  res.aov = rstatix::anova_test(data, formula = formula, wid = ID, within = !!within_var_sym)
  
  # Perform pairwise comparison
  pwc = data %>% rstatix::pairwise_t_test(formula, paired = FALSE, p.adjust.method = test)
  
  # Add xy position of pwc
  pwc <- pwc %>% rstatix::add_xy_position(x = as_label(enquo(grouping_var)), fun = "mean_se", step.increase =0.18)
  
  # Generate plot with means and standard errors
  plot <- ggplot(summary_stats, aes(x = {{grouping_var}}, y = mean, color = {{grouping_var}})) +
    geom_point(size = 3, alpha = 0.7) +  # Add points for the means
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 1) +  # Add error bars
    ggplot2::scale_color_manual(values = plot_colors) +  # Apply custom colors
    ggplot2::labs(
      y = ylabel,
      x = xlabel,
      title = "Mean with Standard Error",
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = rstatix::get_pwc_label(pwc)
    ) +
    ylim(22.5, 24)+
    ggpubr::stat_pvalue_manual(pwc) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.subtitle = element_text(size = 10, hjust = 0)
    )
  
  return(plot)
}

concatCirc = function(table, record, filename, pos){
  # bind data
  table <- rbind.data.frame(table, record)
  # remove first empty row
  if(is.na(table[1, 1])){
    table = table[-1,]
  }
  # change rowname
  rownames(table)[dim(table)[1]] <- filename
  return(table)
}

#' function to save plots to a specific location and format
savePlots = function(obj_to_save, basepath, filename = "", extension, p.width = 560, p.height = 1100) {
  if (filename != ""){filename = paste(filename, "_", sep = "")}
  for (i in seq_along(obj_to_save)) {
    # Check if the list has names; if not, use a default numeric index
    plot_name <- if (!is.null(names(obj_to_save))) {
      names(obj_to_save)[i]
    } else {
      paste0("plot_", i)
    }
    # Create the file path
    plot_path = file.path(basepath, paste0(filename, plot_name, ".", extension))
    
    # Save the plot
    ggplot2::ggsave(plot_path, plot = obj_to_save[[i]], width = p.width, height = p.height, units = "px")
  }
}

# move plots in separate folder
pull_plots = function(base.dir, plot.name, dir.name, file.names){
  # browser()
  folder.vec <- sapply(file.names, function(filename, base.dir){
    file.path(base.dir, paste0(filename, "_results"), "spatial_analysis")
  }, base.dir)
  plot.names <- sapply(file.names, function(filename, plot.name){
    paste0(filename, "_", plot.name)
  }, plot.name)
  plot.path = file.path(folder.vec, plot.names)
  
  # create folder where to store plots if existing
  plots.dir = file.path(base.dir, "plots")
  if(!dir.exists(plots.dir)){
    dir.create(plots.dir)
  }
  plots.dir.second = file.path(plots.dir, dir.name)
  if(!dir.exists(plots.dir.second)){
    dir.create(plots.dir.second)
  }
  
  # copy file in position
  sapply(plot.path, function(file, plots.dir.second){
    destination.path = file.path(plots.dir.second, basename(file))
    file.copy(file, destination.path, overwrite = TRUE, )
  }, plots.dir.second)
}

# FILEPATHS -----

#' Base directory
wd = r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Ph_D\Science\Lab\Main_proj\4_extra\Tau_an_dev\wk5_hemislices\Left\test)"
wd = back_to_forw(wd)

files = list.files(path = wd, pattern = ".tif*$")
filenames = stringr::str_remove(files, pattern = ".tif*")
foldernames = file.path(wd, paste(filenames, "_results", sep = ""))

saving = TRUE
savingtRACE = TRUE

# initialize table for circular values
circ_summary = data.frame(matrix(nrow = 1, ncol = 9))
colnames(circ_summary) = c("close_vec_len", "close_variance", "close_meanPh", "mid_vec_len", "mid_variance", "mid_meanPh", "far_vec_len", "far_variance", "far_meanPh")

# EXECUTION ===== 

#' start for cycle here
for (i in seq(from = 1, to = length(files))){
  
  # FILE IMPORTING ####
  
  foldername = foldernames[i]
  filename = filenames[i]
  print(paste("Processing ", filename, "..."))
  #' import ROI of SCN/slice (manually annotated)
  SCN_ROI_path = file.path(foldername, paste(filename, "_SCN_nuclei.roi", sep = ""))
  SCN_roi <- importROI(SCN_ROI_path)
  
  #' import file with grid coordinates
  grid_coord_path = file.path(foldername, "grid_centroids.csv")
  grid_coord = read.csv(grid_coord_path)
  
  #' convert to sf points
  print("Converting grid coordinates")
  grid_points_sf <- sf::st_as_sf(grid_coord, coords = c("X", "Y"))
  
  #' import file with grid values ch1 (fred)
  ch1_grid_vals_path = file.path(foldername, "Ch1_1_grid_vals.csv")
  print("Importing Ch1 grid vals")
  ch1_grid_vals = importGridVals(ch1_grid_vals_path)
  
  #' import file with grid values ch1 (GCaMP)
  ch2_grid_vals_path = file.path(foldername, "ch2_2_grid_vals.csv")
  print("Importing Ch2 grid vals")
  ch2_grid_vals = importGridVals(ch2_grid_vals_path)
  
  #' import file with particle analysis of Tau
  tau_parts_ROI_path = file.path(foldername, paste("SCN_puncta_ROIs", filename, "_.zip", sep = ""))
  tau_polygons_path = file.path(foldername, paste(filename, "_particle_polygons.rds", sep = ""))
  if(file.exists(tau_polygons_path)){
    tau_polygons <- readRDS(tau_polygons_path)
  }else{
    print("Starting ROIs import")
    tau_polygons <- importROIs(tau_parts_ROI_path)
    # save tau_polygons object for faster retrieval
    saveRDS(tau_polygons, file = tau_polygons_path)
  }
  
  #' Create folder where to save all outputs
  newdir = file.path(foldername, "spatial_analysis")
  
  if(!dir.exists(newdir)){
    dir.create(newdir)
  }
  
  # MANAGE ROIs AND CELLS ####
  
  #' divide ROIs based on frame, nesting the list
  names_frames <- sapply(names(tau_polygons), strsplit, "[-]") %>% data.table::data.table(do.call(rbind, .))
  colnames(names_frames) = c("ROI", "frame", "num", "YM")
  names_frames[,1] <- seq(1:dim(names_frames)[1])
  tau_polygons_frames <- split(tau_polygons, names_frames$frame)
  
  #' compute what cells are inside the SCN
  within_scn <- sf::st_within(grid_points_sf, SCN_roi, sparse = FALSE)
  
  #' filter to only keep cells inside the SCN
  grid_points_scn =  grid_points_sf[within_scn, ] %>% purrr::set_names(., c("cell_ID", "geometry"))
  ch1_cells = ch1_grid_vals[within_scn, ]
  ch2_cells = ch2_grid_vals[within_scn, ]
  
  # PERIOD ANALYSIS ####
  
  # analyze period
  period_tbl_path = file.path(foldername, paste(filename, "_period_tbl.rds", sep = ""))
  if(file.exists(period_tbl_path)){
    period_tbl <- readRDS(period_tbl_path)
  }else{
    period_tbl <- computePeriod(ch2_cells, filename, excludeNC = TRUE, top = 27, bottom = 16)
    # save period table as RDS
    saveRDS(period_tbl, file = period_tbl_path)
  }
  # get summary of data
  period_summary_list <- summarizePeriod(period_tbl)
  
  # SPATIAL MAPS ####

  # map of red values
  # get single value for red channel
  ch1_cells_median = data.frame(ID = as.integer(rownames(ch1_cells)), 
                                   Intensity = matrixStats::rowMedians(ch1_cells))
  
  spatial_red <- highlight_cells_period(period_table = ch1_cells_median, variable = "Intensity", filename = filename, colorcode = "red")
  # map of green values
  ch2_cells_median = data.frame(ID = as.integer(rownames(ch2_cells)), 
                                SD = matrixStats::rowSds(ch2_cells))
  
  spatial_green <- highlight_cells_period(period_table = ch2_cells_median, variable = "SD", filename = filename, colorcode = "green")
  
  # plot period distribution
  spatial_period <- highlight_cells_period(period_table = period_tbl, variable = "period", filename = filename, colorcode = "purple")
  spatial_amplitude <- highlight_cells_period(period_table = period_tbl, variable = "amplitude", filename = filename, colorcode = "green", sdFactor = 3.5)
  spatial_error <- highlight_cells_period(period_table = period_tbl, variable = "error", filename = filename, colorcode = "blue", sdFactor = 1.5)
  spatial_phases <- highlight_cells_period_circ(period_table = period_tbl, variable = "phase_h", filename = filename)
  spatial_phases_circ <- highlight_cells_period_circ(period_table = period_tbl, variable = "phase_circ", filename = filename)
  
  # Save plots
  if(saving){
    savePlots(obj_to_save = list(spatial_period = spatial_period, spatial_amplitude = spatial_amplitude,
                                 spatial_error = spatial_error, spatial_phases = spatial_phases,
                                 spatial_phases_circ = spatial_phases_circ, spatial_fred = spatial_red,
                                 spatial_green = spatial_green), filename = filename, basepath = newdir, extension = "svg", p.width = 1200, p.height = 1390)}
  
  # PARTICLE-CELLS DISTANCES AND PLOT ####
  
  #'calculate min distance of all cells to particles
  grid_vals_path = file.path(foldername, paste(filename, "_grid_vals.rds", sep = ""))
  if(file.exists(grid_vals_path)){
    grid_vals <- readRDS(grid_vals_path)
  }else{
    #' initiate progress bar
    print("Calculating distances...")
    pb <- txtProgressBar(min = 0, max = length(tau_polygons_frames), style = 3)
    
    grid_vals <- do.call(rbind, lapply(seq(1, length(tau_polygons_frames)), function(frame){
      result <- cell2particlesDistances(polygon_frames_list = tau_polygons_frames,
                                        cells_points = grid_points_scn,
                                        frame = frame)
      
      # Update the progress bar
      setTxtProgressBar(pb, frame)
      return(result)
    }))
    
    close(pb)
    # save grid_vals object for faster retrieval
    saveRDS(grid_vals, file = grid_vals_path)
    print("Completed")
  }
  
  distance_mtx_path = file.path(foldername, paste(filename, "_dist_mtx.rds", sep = ""))
  if(file.exists(distance_mtx_path)){
    distance_matrix <- readRDS(distance_mtx_path)
  }else{
    #' generate matrix of distances
    distance_matrix <- cbind(cell_ID = grid_vals$cell_ID, frame = grid_vals$frame, min_distance = grid_vals$min_distance) %>%
      as.data.frame() %>%
      tidyr::pivot_wider(names_from = frame, values_from = min_distance) %>%
      select(-cell_ID) %>%
      as.matrix() %>%
      `rownames<-`(rownames(grid_points_scn))
    # save distance matrix
    saveRDS(distance_matrix, file = distance_mtx_path)
  }
  
  distance_matrix_median <- round(matrixStats::rowMedians(distance_matrix), 3)
  
  
  #' interrogate the matrix to see which cells are at less than
  limit1 = 3
  limit2 = 12
  
  close_cells_matrix = distance_matrix_median <= limit1
  far_cells_matrix = distance_matrix_median > limit2
  between_cells_matrix = distance_matrix_median >limit1 & distance_matrix_median <=limit2
  
  color_palette = c("#B23A48", "#2F9C95", "#B6DC76")
  
  #' generate plot to see which cells are close and far from particles
  groups_cells_plot = highlight_cells(evaluation_matrices = list(close_cells_matrix, between_cells_matrix, far_cells_matrix), filename = filename, colors = color_palette)
  
  
  #' save plot
  if(saving){
    savePlots(obj_to_save = list(groups_cells_plot = groups_cells_plot), filename = filename, basepath = newdir, extension = "svg", p.width = 1200, p.height = 1390)}
  
  # MERGE PERIOD AND PARTICLE ANALYSIS ####
  
  # analysis to visualize circadian parameters in relation to the distance from plaques
  distance_table <- cbind(as.numeric(rownames(distance_matrix)), as.numeric(distance_matrix_median)) %>% `colnames<-`(c("ID", "distance"))#paste("Distance_f", frame, sep = "")))
  # create a median versio of the distance table to merge with the pariod table
  merged_table <- merge(period_tbl, distance_table, by = "ID")
  
  # Define distance groups
  merged_table$distance_group <- cut(merged_table[,11],
                                     breaks = c(-Inf, limit1, limit2, Inf),
                                     labels = c(paste("<=", limit1), paste(limit1, "-",limit2 ), paste(">", limit2)))
  
  merged_table$distance_group <- as.factor(merged_table$distance_group)
  merged_table$ID <- as.factor(merged_table$ID)
  
  # save RDS image of merged table
  merged_tbl_path = file.path(foldername, paste(filename, "_merged_tbl.rds", sep = ""))
  saveRDS(merged_table, file = merged_tbl_path)
  
  # CIRCADIAN PARAMS PLOTS ####
  
  # Period vs Distance
  period_boxplot <- metric_mean_plot(data = merged_table, metric = period, grouping_var = distance_group, plot_colors = color_palette)
  
  # Amplitude vs Distance
  amplitude_boxplot <- metric_mean_plot(data = merged_table, metric = amplitude, grouping_var = distance_group, plot_colors = color_palette)
  
  # Error vs Distance
  error_boxplot <- metric_mean_plot(data = merged_table, metric = error, grouping_var = distance_group, plot_colors = color_palette)
  
  # Save plots
  if(saving){
    savePlots(obj_to_save = list(period_boxplot_w = period_boxplot, amplitude_boxplot_w = amplitude_boxplot,
                                 error_boxplot_w = error_boxplot), filename = filename, basepath = newdir, extension = "svg", p.width = 1500, p.height = 800)}
  
  # Period vs Distance
  period_dotplot_trace <- ggplot(merged_table, aes(x = distance, y = period)) +
    #geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
    geom_smooth(aes(colour = distance_group), method = "lm", se =  TRUE, level = 0.95) +
    scale_color_manual(values = color_palette)+
    ylim(22.5, 24.5) +
    theme_minimal() +
    labs(title = "Period by Distance", x = "Distance from Particle (px)", y = "Period (h)")
  
  # Amplitude vs Distance
  amplitude_dotplot_trace <- ggplot(merged_table, aes(x = distance, y = amplitude)) +
    #geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
    geom_smooth(aes(colour = distance_group), method = "lm", se = TRUE, level = 0.95) +
    scale_color_manual(values = color_palette)+
    ylim(0, 130) +
    theme_minimal() +
    labs(title = "Amplitude by Distance", x = "Distance from Particle (px)", y = "Amplitude (A.U.)")
  
  
  # Error vs Distance
  error_dotplot_trace <- ggplot(merged_table, aes(x = distance, y = error)) +
    #geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
    geom_smooth(aes(colour = distance_group), method = "lm", se =  TRUE, level = 0.95) +
    theme_minimal() +
    scale_color_manual(values = color_palette)+
    ylim(0, 30) +
    labs(title = "Error by Distance", x = "Distance from Particle (px)", y = "Error (A.U.)")
  
  # Save plots
  if(savingtRACE){
    savePlots(obj_to_save = list(period_dotplot_trace = period_dotplot_trace, amplitude_dotplot_trace = amplitude_dotplot_trace, error_dotplot_trace = error_dotplot_trace),
              filename = filename, basepath = newdir, extension = "svg", p.width = 1800, p.height = 1000)}
  
  # PHASE PLOTS ####
  
  # Phase vs distance
  groups <- unique(merged_table$distance_group)
  # add raileygh plots
  sub1 <- merged_table[merged_table$distance_group == groups[1], ]
  sub2 <- merged_table[merged_table$distance_group == groups[2], ]
  sub3 <- merged_table[merged_table$distance_group == groups[3], ]
  
  far_circ <- circular_plot(data = sub1$phase_rad, path = file.path(newdir, paste0(filename, "_far_cells.svg")), saving = saving)
  mid_circ <- circular_plot(data = sub2$phase_rad, path = file.path(newdir, paste0(filename, "_mid_cells.svg")), saving = saving)
  close_circ <- circular_plot(data = sub3$phase_rad, path = file.path(newdir, paste0(filename, "_close_cells.svg")), saving = saving)
  
  # isolate variance components and plot
  newrecord = c(close_vec_len = close_circ$vec_len,
                close_variance = close_circ$variance,
                close_meanPh = close_circ$meanP[[1]],
                mid_vec_len = mid_circ$vec_len,
                mid_variance = mid_circ$variance,
                mid_meanPh = mid_circ$meanP[[1]],
                far_vec_len = far_circ$vec_len,
                far_variance = far_circ$variance,
                far_meanPh = far_circ$meanP[[1]])
  circ_summary = concatCirc(circ_summary, newrecord, filename, i)
  
  # generate file for plot summary
  circStats_plot(newrecord, path = file.path(newdir, paste0(filename, "_phases_summary.svg")), colors = color_palette)
  
  
}

# OPS AFTER CYCLE IS FINISHED ####

#save circ_summary table
write.csv(circ_summary, file = file.path(wd, "circular_stats.csv"), row.names = TRUE)


# GROUP PLOTS IN FOLDER
pull_plots(wd, plot.name = "phases_summary.svg", dir.name = "phases_summary", file.names = filenames)
pull_plots(wd, plot.name = "groups_cells_plot.svg", dir.name = "group_maps", file.names = filenames)
pull_plots(wd, plot.name = "spatial_amplitude.svg", dir.name = "spatial_amp", file.names = filenames)
pull_plots(wd, plot.name = "spatial_period.svg", dir.name = "spatial_per", file.names = filenames)
pull_plots(wd, plot.name = "spatial_error.svg", dir.name = "spatial_err", file.names = filenames)
pull_plots(wd, plot.name = "spatial_fred.svg", dir.name = "spatial_fred", file.names = filenames)
pull_plots(wd, plot.name = "spatial_green.svg", dir.name = "spatial_green", file.names = filenames)
pull_plots(wd, plot.name = "spatial_phases_circ.svg", dir.name = "spatial_phase", file.names = filenames)


# REPORT -----

# Install and load required packages for report generation
# Install and load required packages for report generation
if (!require(rmarkdown)) install.packages("rmarkdown")
if (!require(knitr)) install.packages("knitr")

# Function to generate ordered report
generate_report <- function(filename, output_dir) {
  # Define ordered plot list
  ordered_plots <- c(
    "spatial_amplitude.svg", "spatial_error.svg", "spatial_period.svg", "spatial_phases.svg",
    "groups_cells_plot.svg", "period_boxplot.svg", "amplitude_boxplot.svg", "error_boxplot.svg",
    "period_dotplot.svg", "amplitude_dotplot.svg", "error_dotplot.svg",
    paste0(filename, "_close_cells.svg"),
    paste0(filename, "_mid_cells.svg"),
    paste0(filename, "_far_cells.svg")
  )
  
  # Define output report path
  report_path <- file.path(output_dir, paste0(filename, "_report.html"))
  
  # Generate Rmarkdown file content
  report_content <- c(
    "---",
    "title: \"Analysis Report\"",
    "output: html_document",
    "---",
    "",
    paste("# Report for file:", filename),
    "",
    "## Generated Plots",
    ""
  )
  
  # Retrieve all SVG files in the correct order and add to report
  for (plot_name in ordered_plots) {
    plot_file <- file.path(output_dir, plot_name)
    
    # Check if file exists before adding it to the report
    if (file.exists(plot_file)) {
      # Structure first 8 plots in a 2x2 grid, then dot plots in single row, then last 3 in 2 per row
      if (plot_name %in% ordered_plots[1:8]) {
        report_content <- c(report_content, paste("![](", plot_file, "){width=45%}", collapse = " "), "")
      } else if (plot_name %in% ordered_plots[9:11]) {
        report_content <- c(report_content, paste("![](", plot_file, "){width=90%}", collapse = " "), "")
      } else {
        report_content <- c(report_content, paste("![](", plot_file, "){width=45%}", collapse = " "), "")
      }
    }
  }
  
  # Write report content to an Rmd file
  rmd_file <- file.path(output_dir, paste0(filename, "_report.Rmd"))
  writeLines(report_content, rmd_file)
  
  # Render the Rmd file to HTML
  rmarkdown::render(input = rmd_file, output_file = report_path)
  
  # Clean up Rmd file after rendering
  unlink(rmd_file)
  
  message(paste("Report generated at:", report_path))
}

# Example


# Call generate_report function after each filename processing
# Example usage within your loop
for (i in seq_len(length(files))) {
  # Your existing code for processing
  
  # Call report generation function
  generate_report(filename = filenames[i], output_dir = file.path(foldernames[i], "spatial_analysis"))
}


# find SD of cells in Red channel and show the associated traces from green cells
sd_ch1 <- ch1_cells[sd(ch1_cells),]
# sort and find top 20
high.SD.cells <- names(sort(sd_ch1, decreasing = TRUE)[1500:1520])
# find names of rows in the matrix
high.SD.cells.pos <- which(rownames(ch1_cells) %in% high.SD.cells)
# plot red traces
high.SD.cells.traces.red <- ch1_cells[high.SD.cells.pos,] %>%
  t(.) %>%
  as.data.frame(.) %>%
  cbind(rownames(.), .)%>%
  tidyr::pivot_longer(., names_to = "ID", values_to = "Red.intensity", cols = colnames(.)[-1]) %>%
  `colnames<-`(., c("frame", "ID", "Red.intensity"))
plot.red <- ggplot2::ggplot(data = high.SD.cells.traces.red, ggplot2::aes(x = frame, y = Red.intensity, group = ID))+
  ggplot2::geom_line()
plot.red
# plot associated green traces
high.SD.cells.traces.green <- ch2_cells[high.SD.cells.pos,] %>%
  t(.) %>%
  as.data.frame(.) %>%
  cbind(rownames(.), .)%>%
  tidyr::pivot_longer(., names_to = "ID", values_to = "Green.intensity", cols = colnames(.)[-1]) %>%
  `colnames<-`(., c("frame", "ID", "Green.intensity"))
plot.green <- ggplot2::ggplot(data = high.SD.cells.traces.green, ggplot2::aes(x = frame, y = Green.intensity, group = ID))+
  ggplot2::geom_line()
plot.green
#'
#' compute distance between each cell centroid and Tau particles to find the
#' closest one through every frame. Order them based on distance and keep the
#' closest 10.
#'
#' Import the file with all the descriptors for each plaque and correlate
#' intensity to period disruption
#'
#' also find the density in the surrounding area to determine the "exposure"
#' classify cells based on their proximity to plaques
#' classify cells based on their location (i.e. dorsal or ventral)


# OLD CODE #####
#' Calculate cells that colocalize with particles
for (frame_index in seq_along(tau_polygons_frames)) {
  # Convert ROIs for the current frame into a single spatial object
  rois <- st_sfc(tau_polygons_frames[[frame_index]])  # Adjust CRS if needed
  
  # Check if each cell point is within any ROI in the frame
  inside <- st_intersects(grid_points_scn, rois, sparse = FALSE)
  
  # Update results matrix: any TRUE in a row indicates the cell is within an ROI
  results_matrix[, frame_index] <- apply(inside, 1, any)
}

# Import ROI
roi_path <- r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Ph_D\Science\Lab\Main_proj\4_extra\Tau_an_dev\SCN_roi.roi)"

roi_sf <- importROI(roi_path)

# import file with particle analysis
particle_path <- r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Ph_D\Science\Lab\Main_proj\4_extra\Tau_an_dev\Results.csv)"
particles_res <- read.csv(particle_path)

particle_no <- nrow(particles_res)
particle_coord <- particles_res[,c(1, 8, 9)]

# convert to sf points
particle_points_sf <- sf::st_as_sf(particle_coord, coords = c("XM", "YM"))
within_scn <- sf::st_within(particle_points_sf, roi_sf, sparse = FALSE)

# show fraction of TRUE and FALSE vals
tf_fraction <- table(within_scn)/length(within_scn)
tf_data <- data.frame(Value = names(tf_fraction),
                      Fraction = as.numeric(tf_fraction))

plot_fraction = ggplot2::ggplot(tf_data, ggplot2::aes(x = Value, y = Fraction, fill = Value))+
  ggplot2::geom_bar(stat = "identity", width = 0.5)+
  ggplot2::scale_y_continuous(limits = c(0,1))+
  ggplot2::labs(y = "Fraction")

# code removed from ComputePeriod function
# compute other period values
period_var = stats::var(period_table$period, na.rm = TRUE)
period_mean = mean(period_table$period, na.rm = TRUE)
amplitude_mean = mean(abs(period_table$amplitude), na.rm = TRUE)
error_mean = mean(period_table$error, na.rm = TRUE)

#save table
namedatapath <- stringr::str_remove(datapath, ".csv")
tabledatapath <- stringr::str_c(namedatapath, "_FFT_NLLS.csv")
write.csv(period_table, file = tabledatapath)

#create plot with period length
filename = filename
period_plot <- period_print(data = period_table, filename)

#save period plot
svgperiodpath <- stringr::str_c(namedatapath, "_period.svg")
ggplot2::ggsave(svgperiodpath, plot = period_plot, width = 2.7, height = 4.8, units = "in")

# create circular plot

svgdatapath <- stringr::str_c(namedatapath, "_circular.svg")
returns_list <- circular_plot(period_table$phase_rad, filename, path = svgdatapath)


# return values outside
outreturn_list <- list(variance = returns_list$variance,
                       trace_no = non_na_count,
                       trace_tot = period_counts,
                       period_var = period_var,
                       period_mean = period_mean,
                       amplitude_mean = amplitude_mean,
                       error_mean = error_mean,
                       vector_len = returns_list$vec_len)

return(outreturn_list)

# other multiple comparison tests
anova_model <- aov(error ~ distance_group, data = merged_table)

# Perform Tukey's HSD post-hoc test
tukey_result <- TukeyHSD(anova_model)

if (!requireNamespace("DescTools", quietly = TRUE)) {
  install.packages("DescTools")
}
library(DescTools)

# Perform Dunnett's test (specify the control group as needed)
dunnett_result <- DescTools::DunnettTest(error ~ distance_group, data = merged_table, control = "<=5")

# Print the results
print(dunnett_result)

highlight_cells = function(evaluation_matrix, frame, all_cells = grid_coord, scn_cells = grid_points_scn, active_col = "red", scn_col = "grey", outside_col = "black"){
  #' Copy `all_cells` as the base data frame
  plot_data <- all_cells
  
  #' Extract coordinates from `scn_cells`
  scn_coords <- st_coordinates(scn_cells)
  
  #' Add `is_SCN` column to identify SCN coordinates
  plot_data$is_SCN <- with(plot_data, paste(X, Y) %in% paste(scn_coords[, "X"], scn_coords[, "Y"]))
  
  #' Add `is_particle` column for coordinates where `result_matrix` is TRUE for frame `n`
  active_pixels <- evaluation_matrix[, frame]
  active_coords <- all_cells[as.integer(names(which(active_pixels))), ]  # This selects rows of `grid_coord` that correspond to TRUE in `result_matrix`
  
  # Identify rows in `plot_data` that match active pixels for the nth frame
  plot_data$is_particle <- with(plot_data, paste(X, Y) %in% paste(active_coords$X, active_coords$Y))
  
  
  # Define colors based on `is_SCN` and `is_particle`
  plot_data$color <- ifelse(plot_data$is_particle, active_col,
                            ifelse(plot_data$is_SCN, scn_col, outside_col))
  
  # Create the plot
  plot <- ggplot(plot_data, aes(x = X, y = Y, color = color)) +
    geom_point() +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Frame", frame, "- Pixel Activation"), x = "X Coordinate", y = "Y Coordinate") +
    coord_fixed()+
    scale_y_reverse()
  return(plot)
}

# close_cells_plot = highlight_cells(list(close_cells_matrix), frame = 2, colors = c("#B23A48"))
# far_cells_plot = highlight_cells(list(far_cells_matrix), frame = 40, colors = c("#2F9C95"))
# between_cells_plot = highlight_cells(list(between_cells_matrix), frame = 40, colors = c("#B6DC76"))

#' Extract IDs and TS of cells close to particles
# close_cells_IDs = which(close_cells_matrix[, 2] == TRUE) %>% names()
# close_cells_ch2_TS = ch2_cells[rownames = close_cells_IDs ,]

#' plot timeseries of cells in specific regions from particles
# plot_TS <- t(close_cells_ch2_TS) %>% as.data.frame() %>% cbind("Time" = seq(from = 0, length.out = length(plot_TS[, 1]), by = 0.5), .)
# plot_TS_long <- tidyr::pivot_longer(plot_TS, cols = c(2:length(plot_TS[1,])), values_to = "intensity", names_to = "ID") %>% order("ID")
#
# plot_TS <- ggplot2::ggplot(data = plot_TS_long, ggplot2::aes(x = Time, y = intensity, color = ID))+
# geom_line() #TODO add color or grouping based on distance from particles

