# FUNCTIONS used in Tau spatial analysis script

#' function to import a single ROI
importROI = function(path, pixelfct = 2.82, convert = TRUE){
  roi_path = stringr::str_replace_all(path, "\\\\", "//")
  roi_data <- RImageJROI::read.ijroi(roi_path)
  # convert to sf package compatible
  roi_coords <- roi_data$coords
  # convert to units of micrometers
  if(convert){
    roi_coords <- roi_coords * pixelfct
  }
  if (!all(roi_coords[1, ] == roi_coords[nrow(roi_coords), ])) {
    # Add the first coordinate at the end to close the polygon
    roi_coords <- rbind(roi_coords, roi_coords[1, ])
  }
  roi_polygon <- sf::st_polygon(list(roi_coords))
  # roi_sf <- sf::st_sfc(roi_polygon)
  return(roi_polygon)
}

#' function to import multiple ROIs
importROIs = function(path, pixelfct = 2.82, convert = TRUE){
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
      # convert to units of micrometers
      if(convert){
        coord_vals <- coord_vals * pixelfct
      }
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
highlight_cells <- function(evaluation_matrices, filename, colors, all_cells = grid_coord, scn_cells = grid_points_scn, scn_col = "grey", outside_col = "black", shape = 15, size = 2) {
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
    geom_point(shape = shape, size = size) +
    scale_color_identity() +
    scale_alpha_identity()+
    scale_size_identity()+
    theme_minimal() +
    labs(title = paste(filename, "- Distance groups"), x = "Width (µm)", y = "Height (µm)") +
    coord_fixed() +
    scale_y_reverse()+
    theme(plot.title = element_text(size = 11),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 15))
  
  return(plot)
}

#' plot features of period on map
highlight_cells_period = function(all_cells = grid_coord, period_table, variable, filename = "", colorcode = "blue", sdFactor = 1, shape=15, size = 2, ...){
  
  plot_data <- grid_coord %>% `colnames<-`(c("ID", "X", "Y"))
  # Merge tables to include period data
  merged_data <- left_join(plot_data, period_table, by = "ID")
  
  # Choose the variable to visualize (e.g., period, amplitude, error)
  variable <- variable
  
  mean_value <- mean(merged_data[[variable]], na.rm = TRUE)
  sd_value <- sd(merged_data[[variable]], na.rm = TRUE)
  
  args = list(...)
  upper_bound <- args$upper_bound
  lower_bound <- args$lower_bound
  
  if(is.null(upper_bound)){
    upper_bound <- mean_value + (sdFactor * sd_value)
  }
  if(is.null(lower_bound)){
    lower_bound <- mean_value - (sdFactor * sd_value)
    lower_bound <- ifelse(lower_bound < 0, 0, lower_bound)
  }
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
    geom_point(shape = shape, aes(colour = !!sym(variable), size = size)) +
    scale_color_gradientn(
      colors = colorscheme,
      values = scales::rescale(c(lower_bound, mean_value, upper_bound)),
      limits = c(lower_bound, upper_bound),
      oob = scales::squish  # Ensures outliers are squished into the limits
    ) +
    scale_alpha_identity()+
    scale_size_identity()+
    theme_minimal() +
    labs(
      title = paste0(filename, " - Distribution of ", variable),
      x = "Width (µm)",
      y = "Height (µm)",
    ) +
    coord_fixed() +
    scale_y_reverse()+
    theme(plot.title = element_text(size = 11),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 15))
  
  return(plot)
}

#' plot features of period on map
highlight_cells_period_circ = function(all_cells = grid_coord, period_table, variable, filename = "", shape=15, size = 2, normalized = FALSE){
  
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
    geom_point(shape = shape, aes(color = !!sym(variable)), size = size) +
    theme_minimal() +
    scale_alpha_identity()+
    scale_size_identity()+
    labs(title = paste0(filename, " - Distribution of ", variable), x = "Width (µm)", y = "Height (µm)") +
    coord_fixed() +
    scale_y_reverse()+
    theme(plot.title = element_text(size = 11),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 15))
  
  if(normalized){
    plot <- plot +
      scale_color_gradientn(colors = colorscheme, limits = c(-12.15, 12.15), breaks = c(-12, -6, 0, 6, 12))
  } else {
    plot <- plot +
      scale_color_gradientn(colors = colorscheme, limits = c(-0.15, 24.15), breaks = c(0, 6, 12, 18, 24))
  }
  
  return(plot)
}

# fix path
back_to_forw <- function(x){
  stringr::str_replace_all(x, "\\\\", "//")
}

# prepare data into long table format
prep_table = function(data, rm.start = 0, preprocess = TRUE, transp = TRUE, add_t = TRUE){

  # re-transpose table
  if(transp){
  data = t(data)
  }
  
  # remove initial data
  if(rm.start != 0){
    print("removing data")
    data = data[-c(seq(1, rm.start)), ]
  }
  
  # add time column
  if(add_t){
  t = seq(0, length.out = length(data[, 1]), by = 0.5)
  data = cbind(t, data)
  }
  
  # perform outliers removal, smoothening and detrending
  if(preprocess){
  data_clean <- preprocess_data(data, grade = 3)
  } else {
    data_clean <- as.data.frame(data)
  }
  
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
FFT_NLLS_analyse = function(data, minPer = 16, maxPer = 32){
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
  fitted_phase <- circular::as.circular(fitted_phase, type = "angles", units = "radians", 
                                        rotation = "clock", template = "none", modulo = "asis", zero = 0)
  phase_circadian <- circular::conversion.circular(fitted_phase, units = "hours")
  phase_absolute <- phase_circadian*(fitted_period/24)
  
  # Compute residual error
  fitted_values <- sinusoidal_model(fitted_params, data$t)
  residual_error <- sqrt(mean((data$value - fitted_values)^2))
  
  # Calculate R-squared (GOF)
  ss_total <- sum((data$value - mean(data$value))^2)
  ss_res <- sum((data$value - fitted_values)^2)
  r_squared <- 1 - (ss_res / ss_total)
  
  # Compute standard deviation of residuals
  residual_std <- sqrt(sum((data$value - fitted_values)^2) / (length(data$value) - length(fitted_params)))
  
  # Compute RAE
  RAE <- residual_std / fitted_amplitude
  
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
    GOF = r_squared,
    RAE = RAE
  )
  return(result)
}

# create circular plot
circular_plot = function(data, path, saving = TRUE, ...){
  phases_rad = data
  if(!all(is.na(phases_rad))){
    real_MeanPh <- circular::mean.circular(phases_rad, Rotation = "counter", na.rm = TRUE, type = "angles", units = "radians", 
                                           template = "none", modulo = "asis", zero = 0)
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
    real_MeanPh = NA
    plot.new()
    circular::plot.circular(phases, bins = 24, axes = FALSE, type = "angles", units = "radians", 
                            rotation = "clock", template = "none", modulo = "asis", zero = 0)
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
  circular::plot.circular(NA, bins = 24, axes = FALSE, type = "angles", units = "radians", 
                          rotation = "clock", template = "none", modulo = "asis", zero = 0)
  # Add radial axis labels to the plot
  circular::axis.circular(circular(seq(0, 2*pi, pi/2)), zero=0, rotation = 'counter', labels=c("6", "00", "18", "12", "6"), tcl.text=0.12, cex = 1.5, 
                          units = "radians", template = "none", modulo = "asis")
  # Add arrows to the plot based on the mean and vector length
  # close arrow
  vectorLength = data["close_vec_len"][[1]]
  meanTest = data["close_meanPh"][[1]]
  variance = data["close_variance"][[1]]
  color = colors[1]
  circular::arrows.circular(meanTest, vectorLength, length = 0.1, col = color, lwd = 3,
                            rotation = "clock")
  # mid arrow
  vectorLength = data["mid_vec_len"][[1]]
  meanTest = data["mid_meanPh"][[1]]
  color = colors[2]
  circular::arrows.circular(meanTest, vectorLength, length = 0.1, col = color, lwd = 3,
                            rotation = "clock")
  # far arrow
  vectorLength = data["far_vec_len"][[1]]
  meanTest = data["far_meanPh"][[1]]
  color = colors[3]
  circular::arrows.circular(meanTest, vectorLength, length = 0.1, col = color, lwd = 3, 
                            rotation = "clock")
  
  # Add an empty title (the paste("") effectively removes the default title)
  title(paste(filename), line = 0)
  # close plot
  dev.off()
}

# group traces plot
plot_group_traces = function(summary_data){
  plot = ggplot2::ggplot(summary_data)+
    geom_ribbon(aes(x = t, ymin = mean-se, ymax = mean+se, alpha = 0.5, colour = group, fill = group))+
    geom_line(aes(x = t, y = mean, colour = group))+
    ggplot2::scale_alpha_identity()+
    scale_x_continuous(name = "Time (h)", breaks = seq(0, as.integer(max(t)), by = 24), minor_breaks = seq(12, as.integer(max(t)), by = 24))+
    theme_minimal()+
    theme(axis.line.y = element_line(linewidth = 1),
          axis.line.x = element_line(linewidth = 1),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linewidth = 0.5, linetype = "88", colour = "black"),
          plot.subtitle = element_text(size = 10, hjust = 0),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x  = element_text(size = 20),
          axis.title.y  = element_text(size = 20))
  return(plot)
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

# function to align traces to the mean phase of the group
phase_align_trace = function(traces_table, period_table, remove_start = 0, remove_end = 0, align_to = NA, debug = FALSE){
  #' TODO check align phase value
  
  #' transform rad values in plus minus Pi
  phases = circular::minusPiPlusPi(circular(period_tbl$phase_rad, units = "rad"))
  if(is.na(align_to)){
    align_phase = circular::mean.circular(phases, na.rm = TRUE)
  } else {align_phase = circular::as.circular((align_to/12)*pi, type = "angles", units = "radians", 
                                              rotation = "clock", template = "none", modulo = "asis", zero = 0)}
  ph_diff = align_phase - phases
  adjust = circular::minusPiPlusPi(circular(ph_diff, units = "rad"))
  adjust_frames = round(circular::conversion.circular(adjust, units = "hours")*2, 0)
  # initialize new table to store vecs
  traces_aligned <- replace(traces_table, all(), NA)
  for (i in 1:dim(traces_table)[1]){
    # get vector
    trace = traces_table[i,]
    adjust_fct = adjust_frames[i]
    if(is.na(adjust_fct)){
      clean_trace = rep(NA, dim(traces_table)[2])
    } else if(adjust_fct == 0){
      clean_trace = trace
    } else if(adjust_fct < 0){
      clean_trace = c(rep_len(NA, abs(adjust_fct)), trace[1:(length(trace)-abs(adjust_fct))])
    } else if(adjust_fct > 0){
      clean_trace = c(trace[(abs(adjust_fct)+1):length(trace)], rep_len(NA, abs(adjust_fct)))
    }
    traces_aligned[i,] = clean_trace
    if(debug){
      browser()
      plot(x = 0:(dim(traces_table)[2]-1), y = trace, col = "red", type = "l")
      lines(x = 0:(dim(traces_table)[2]-1), y =clean_trace, col = "blue", type = "l")
    }
  }
  # trim start and end of table
  # evaluate if all have been shifted in one direction
  #' min_shift = min(adjust_frames, na.rm = TRUE)
  #' max_shift = max(adjust_frames, na.rm = TRUE)
  #' shift = min_shift*max_shift
  #' if(shift > 0){
  #'   #' case when all the shift happens in the same direction (need to cut only from one of the two ends)
  #'   #' TODO complete trimming in another function
  #' } else if(shift < 0){
  #'   #' case when shift happens in both direction (you can use the min-max rune simply)
  #'   traces_aligned_trim = traces_aligned[ , (max(adjust_frames, na.rm = TRUE)+1):(dim(traces_aligned)[2]+min(adjust_frames, na.rm = TRUE))]
  #' }
  
  #traces_aligned_trim = traces_aligned[ , (max(adjust_frames, na.rm = TRUE)+1):(dim(traces_aligned)[2]+min(adjust_frames, na.rm = TRUE))]
  
  return(traces_aligned)
}

# foo to analyse TS and make circular plot
computePeriod = function(df, excludeNC = FALSE, top = 30, bottom = 18, save.trace = FALSE, rm.start = 0, prep_tbl = TRUE,  preprocess = TRUE, transp = TRUE, add_t = TRUE, ...){
  
  #prepare table
  if(prep_tbl){
  data_df <- prep_table(df, rm.start, preprocess, transp, add_t)
  }else{
    data_df <- df
  }
  
  # analyse period and store result in list
  unique_ids = unique(data_df$ID)
  results_list <- list()
  traces_list <- list()
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
    
    # add the detrended trace to a list for export
    traces_list[[ID]] <- detrended_data
  }
  close(pb)
  print("Completed")
  # merge all period data and traces into tables
  period_table = do.call(rbind, lapply(results_list, as.data.frame))
  traces_table = do.call(cbind, lapply(traces_list, as.data.frame)) %>% t() %>% `rownames<-`(names(traces_list))
  # add a column for normalized phase values
  period_table$phase_norm = normalize_phase(period_table$phase_rad)
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
  
  # list to return results and traces of period analysis
  return_list = list(traces = traces_table,
                     period_table = period_table) 
  
  return(return_list)
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
    #print(normality)
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

metric_mean_plot = function(data, metric, grouping_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, user_ylims = NULL, ...){
  
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
    #print(normality)
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
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 1, ) +  # Add error bars
    ggplot2::scale_color_manual(values = plot_colors) +  # Apply custom colors
    ggplot2::labs(
      y = ylabel,
      x = xlabel,
      title = "Mean with Standard Error",
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = rstatix::get_pwc_label(pwc)
    ) +
    ggpubr::stat_pvalue_manual(pwc) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.subtitle = element_text(size = 10, hjust = 0)
    )
  
  if(!is.null(user_ylims)){
    plot = plot + ylim(user_ylims)
  }
  
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
pull_plots = function(base.dir, plot.name, dir.name, file.names, mainfold = FALSE, second_fold = ""){
  
  folder.vec <- sapply(file.names, function(filename, base.dir, mainfold, second_fold){
    if(mainfold){
      file.path(base.dir, paste0(filename, "_results"))
    }else{
    file.path(base.dir, paste0(filename, "_results"), second_fold)
      }
  }, base.dir, mainfold, second_fold)
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

# function to create 3D plots with plotly
threeDplot = function(data, x_var, y_var, z_var, color_var, col_palette = NULL, continuous = FALSE, x_lab = "x", y_lab = "y", z_lab = "z", title = "", xrange = NULL, yrange = NULL, zrange = NULL, eye = "pos1", zoom = 1){

  if(eye == "pos1"){
    eye = list(x = 0, y = 1.5, z = 0.8)
  }else if(eye == "pos2"){
    eye = list(x = 1.5, y = 0.3, z = 0.8)
  }else if(eye == "pos3"){
    eye = list(x = -1.5, y = 0.3, z = 0.8)
  }else if(length(eye != 3)){
    print("Camera values are not suitable, setting to default")
    eye = list(x = 0, y = 1.5, z = 0.8)
  }else{
    eye = eye
  }
  
  if(title == ""){
    title = paste0("3D plot of ")#, deparse(x_var), ", ", deparse(y_var), ", ", deparse(z_var))
  }
  
  if(is.null(col_palette)){
    colors = length(unique(data[,deparse(substitute(color_var))]))
  }
  
  # modify camera position according to zoom
  zoom = 1/zoom
  eye <- lapply(eye, function(x){x*zoom})

  # retrieve variable position
  x_pos <- which(colnames(data) == deparse(substitute(x_var)))[1]
  y_pos <- which(colnames(data) == deparse(substitute(y_var)))[1]
  z_pos <- which(colnames(data) == deparse(substitute(z_var)))[1]
  col_pos <- which(colnames(data) == deparse(substitute(color_var)))[1]
  data[, col_pos] <- as.factor(data[, col_pos])
  # setNames(col_palette, levels(data[, col_pos]))
  plotly <- plot_ly() %>%
    # Add scatter plot points
    add_trace(
      data = data,
      x = ~data[, x_pos],
      y = ~data[, y_pos],
      z = ~data[, z_pos],
      type = "scatter3d",
      mode = "markers",
      split = ~data[, col_pos],  # Ensures proper grouping by discrete variable
      marker = list( size = 3, opacity = 0.6#,
                     # color = ~data[, col_pos],  # Assign color based on a column
                     # colors = col_palette  # Apply custom colors
      ),
      showlegend = TRUE) %>%
    # Layout settings
    layout(
      scene = list(
        xaxis = list(title = x_lab, titlefont = list(size = 14)),
        yaxis = list(title = y_lab, titlefont = list(size = 14)),
        zaxis = list(title = z_lab, titlefont = list(size = 14)),
        camera = list(eye = eye)
      ),
      title = list(
        text = title,
        font = list(size = 16)
      )
    )
  
  if(!is.null(xrange)){
    plotly <- plotly %>% layout(scene = list(xaxis = list(range = xrange)))
  }
  if(!is.null(yrange)){
    plotly <- plotly %>% layout(scene = list(yaxis = list(range = yrange)))
  }
  if(!is.null(zrange)){
    plotly <- plotly %>% layout(scene = list(zaxis = list(range = zrange)))
  }
  return(plotly)
}

# Function to calculate AUC using the trapz function
calculate_auc <- function(y, normalize = TRUE) {
  # Trapezoidal method
  time <- seq(0, length.out = length(y), by = 0.5)
  auc_trapz <- trapz(time, y)
  if(normalize){
    auc_trapz <- auc_trapz / length(time)
  }
  
  return(AUC1 = auc_trapz)
}

#' normalize_phase
normalize_phase <- function(angles) {
  # Ensure input is a circular object
  angles_circ <- circular(angles, type = "angles", units = "radians", modulo = "asis")
  
  # Compute circular mean
  mean_angle <- mean.circular(angles_circ)
  
  # Normalize angles to be centered around 0 in the [-pi, pi] range
  normalized_angles <- (angles - as.numeric(mean_angle) + pi) %% (2 * pi) - pi
  
  # Convert radians to hours (scale π to 12)
  angles_in_hours <- normalized_angles * (12 / pi)
  
  return(angles_in_hours)
}

#' Function to create folders if not already existing
create_folder = function(base_dir, folder_name){
  folder_path = file.path(base_dir, folder_name)
  if(!dir.exists(folder_path)){
    dir.create(folder_path)
  }
  return(folder_path)
}
