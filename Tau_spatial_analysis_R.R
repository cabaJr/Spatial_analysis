# install and load required packages ####
pkg_req <- c("sf", "RImageJROI", "ggplot2", "dplyr", "pracma", "minpack.lm",
             "magrittr", "stringr", "circular", "svglite", "astsa", "pdftools",
             "scico", "ggpubr", "rstatix", "matrixStats", "gridGraphics",
             "matrixStats", "plotly", "htmlwidgets")

# Install and load packages if not already installed
for (package in pkg_req) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# FUNCTIONS import ####

source(file = r"(C:\Users\mf420\Documents\GitHub\Spatial_analysis\Tau_spatial_analysis_foo.R)")

# to clean the environment preserving functions
rm(list = setdiff(ls(), lsf.str()))

# FILEPATHS -----

#' Base directory
wd = r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Proj_Tau\Tau_an_dev\wk2_hemislices\Left)"
wd = back_to_forw(wd)

files = list.files(path = wd, pattern = ".tif*$")
filenames = stringr::str_remove(files, pattern = ".tif*")
foldernames = file.path(wd, paste(filenames, "_results", sep = ""))
pixel_fct = 2.82
# reloading options
saving = TRUE
savingtRACE = FALSE
RELOAD = TRUE
RELOADPER = TRUE
RELOAD.MTX = TRUE
# initialize table for circular values
circ_summary = data.frame(matrix(nrow = 1, ncol = 9))
colnames(circ_summary) = c("close_vec_len", "close_variance", "close_meanPh", "mid_vec_len", "mid_variance", "mid_meanPh", "far_vec_len", "far_variance", "far_meanPh")

# EXECUTION ===== 

#' start for cycle here
for (i in seq(from = 7, to = 7)){#length(files))){
  
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
  
  #' Convert centroids to micrometers
  grid_coord$X <- grid_coord$X*pixel_fct
  grid_coord$Y <- grid_coord$Y*pixel_fct
  
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
  if(!RELOAD & file.exists(tau_polygons_path)){
    tau_polygons <- readRDS(tau_polygons_path)
  }else{
    print("Starting ROIs import")
    tau_polygons <- importROIs(tau_parts_ROI_path)
    # save tau_polygons object for faster retrieval
    saveRDS(tau_polygons, file = tau_polygons_path)
  }
  
  #' Create folder where to save different outputs
  spatial_dir = create_folder(foldername, "spatial_maps")
  meanplot_dir = create_folder(foldername, "mean_plots")
  lineplot_dir = create_folder(foldername, "line_plots")
  circplot_dir = create_folder(foldername, "circ_plots")
  threeDplot_dir = create_folder(foldername, "3D_plots")
  
  #' Create folder where to save all rds outputs
  rds_dir = create_folder(foldername, "rds")
   
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
  
  #' save Ch2 cells as rds object
  ch2_cells_path = file.path(rds_dir, paste(filename, "_ch2_cells.rds", sep = ""))
  saveRDS(ch2_cells, file = ch2_cells_path)
  
  # ABSOLUTE CALCIUM LEVELS ####
  
  # Extract the area under the curve of each trace from each spot in the map
  
  # Apply function to each row and combine results into a new data frame
  auc_results <- data.frame(ID = as.integer(rownames(ch1_cells)),
                            AUC = as.data.frame(apply(ch1_cells, 1, calculate_auc))
  ) %>% `colnames<-`(c("ID", "AUC"))
  
  median_vals = data.frame(ID = as.integer(rownames(ch2_cells)), 
                           median = matrixStats::rowMedians(ch2_cells))
  
  SD_vals = data.frame(ID = as.integer(rownames(ch2_cells)), 
                           SD_green = matrixStats::rowSds(ch2_cells))
  
  # PERIOD ANALYSIS ####
  
  # analyze period
  period_tbl_path = file.path(rds_dir, paste(filename, "_period_tbl.rds", sep = ""))
  detr_traces_path = file.path(rds_dir, paste(filename, "_detr_traces.rds", sep = ""))
  if(!RELOADPER & file.exists(period_tbl_path)){
    period_tbl <- readRDS(period_tbl_path)
    detrended_traces <- readRDS(detr_traces_path)
  }else{
    results <- computePeriod(ch2_cells, filename, excludeNC = TRUE, top = 27, bottom = 16, rm.start = 12)
    period_tbl <- results$period_table
    detrended_traces <- results$traces
    
    # save period table and traces as RDS
    saveRDS(period_tbl, file = period_tbl_path)
    saveRDS(detrended_traces, file = detr_traces_path)
  }
  # get summary of data
  period_summary_list <- summarizePeriod(period_tbl)
  
  # SPATIAL MAPS ####

  # map of red values
  # get single value for red channel
  ch1_cells_median = data.frame(ID = as.integer(rownames(ch1_cells)), 
                                   Intensity = matrixStats::rowMedians(ch1_cells))
  
  spatial_red <- highlight_cells_period(period_table = ch1_cells_median, variable = "Intensity", filename = filename, colorcode = "red", shape = 15, size = 2)
  # map of green values
  
  spatial_green_median <- highlight_cells_period(period_table = median_vals, variable = "median", filename = filename, colorcode = "green", shape = 15, size = 2)
  spatial_green_AUC <- highlight_cells_period(period_table = auc_results, variable = "AUC", filename = filename, colorcode = "green", shape = 15, size = 2)
  
  ch2_cells_SD = data.frame(ID = as.integer(rownames(ch2_cells)), 
                            SD = matrixStats::rowSds(ch2_cells))
  # AUC_SD = data.frame(ID = as.integer(rownames(ch2_cells)), 
  #                     AUC_SD = ch2_cells_SD$SD * auc_results$AUC)
  # 
  # spatial_green_AUC_SD <- highlight_cells_period(period_table = AUC_SD, variable = "AUC_SD", filename = filename, colorcode = "green", shape = 15, size = 2)
  
  spatial_green_SD <- highlight_cells_period(period_table = ch2_cells_SD, variable = "SD", filename = filename, colorcode = "green", shape = 15, size = 2)
  
  # plot period distribution
  spatial_period <- highlight_cells_period(period_table = period_tbl, variable = "period", filename = filename, colorcode = "purple", shape = 15, size = 2, lower_bound = 23.2, upper_bound = 24.2)
  spatial_amplitude <- highlight_cells_period(period_table = period_tbl, variable = "amplitude", filename = filename, colorcode = "green", sdFactor = 3.5, shape = 15, size = 2, upper_bound = 130)
  spatial_error <- highlight_cells_period(period_table = period_tbl, variable = "error", filename = filename, colorcode = "blue", sdFactor = 1.5, shape = 15, size = 2)
  spatial_phases <- highlight_cells_period_circ(period_table = period_tbl, variable = "phase_h", filename = filename, shape = 15, size = 2)
  spatial_phases_circ <- highlight_cells_period_circ(period_table = period_tbl, variable = "phase_circ", filename = filename, shape = 15, size = 2)
  spatial_phases_norm <- highlight_cells_period_circ(period_table = period_tbl, variable = "phase_norm", filename = filename, shape = 15, size = 2, normalized = TRUE)
  spatial_RAE <- highlight_cells_period(period_table = period_tbl, variable = "RAE", filename = filename, colorcode = "acqua", shape = 15, size = 2)
  
  p.height = round(length(unique(grid_coord$Y))*22.64)+190
  # Save plots
  if(saving){
    savePlots(obj_to_save = list(spatial_period = spatial_period, spatial_amplitude = spatial_amplitude,
                                 spatial_error = spatial_error, spatial_phases = spatial_phases, spatial_phases_norm = spatial_phases_norm,
                                 spatial_phases_circ = spatial_phases_circ, spatial_fred = spatial_red, spatial_RAE = spatial_RAE,
                                 spatial_green_median = spatial_green_median, spatial_green_AUC = spatial_green_AUC, 
                                 spatial_green_SD = spatial_green_SD), filename = filename, basepath = spatial_dir, extension = "svg", p.width = 1500, p.height = p.height)}
  
  # PARTICLE-CELLS DISTANCES AND PLOT ####
  
  #'calculate min distance of all cells to particles
  grid_vals_path = file.path(rds_dir, paste(filename, "_grid_vals.rds", sep = ""))
  if(!RELOAD.MTX & file.exists(grid_vals_path)){
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
  
  distance_mtx_path = file.path(rds_dir, paste(filename, "_dist_mtx.rds", sep = ""))
  if(!RELOAD.MTX & file.exists(distance_mtx_path)){
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
  limit1 = 8.5
  limit2 = 34
  
  close_cells_matrix = distance_matrix_median <= limit1
  close_cells_IDs = names(close_cells_matrix[close_cells_matrix == TRUE])
  far_cells_matrix = distance_matrix_median > limit2
  far_cells_IDs = names(far_cells_matrix[far_cells_matrix == TRUE])
  between_cells_matrix = distance_matrix_median >limit1 & distance_matrix_median <=limit2
  between_cells_IDs = names(between_cells_matrix[between_cells_matrix == TRUE])
  
  color_palette = c("#B23A48", "#2F9C95", "#B6DC76")
  
  #' generate plot to see which cells are close and far from particles
  groups_cells_plot = highlight_cells(evaluation_matrices = list(close_cells_matrix, between_cells_matrix, far_cells_matrix), filename = filename, colors = color_palette, shape = 15, size = 2)
  
  
  #' save plot
  if(saving){
    savePlots(obj_to_save = list(groups_cells_plot = groups_cells_plot), filename = filename, basepath = spatial_dir, extension = "svg", p.width = 1500, p.height = p.height)}
  
  # LOOK AT TRACES IN DIFFERENT GROUPS #### 
  #' TODO comment code and create function out of it
  
  #' extract the cleaned trace from the period calculation algorithm
  detrended_traces = detrended_traces %>% `colnames<-`(seq(0, by = 0.5, length.out = dim(detrended_traces)[2]))
  #' phase align all the traces
  aligned_traces = phase_align_trace(detrended_traces, period_table = period_tbl, align_to = 6, debug = FALSE) %>% t()
    #' normalize traces in 0-1 interval
  aligned_traces_norm <- as.data.frame(aligned_traces) %>%
    mutate(across(where(is.numeric), ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))))  %>% as.matrix(.)
  t = as.numeric(rownames(aligned_traces_norm))
  aligned_traces_t = cbind(t, aligned_traces_norm) %>% as.data.frame()
  aligned_traces_long = tidyr::pivot_longer(aligned_traces_t,
                                            cols = colnames(aligned_traces_t)[-1],
                                            names_to = "ID",
                                            values_to = "intensity")
  
  #' trial with detrended values instead of normalized
  
  aligned_traces_raw_t = cbind(t, aligned_traces) %>% as.data.frame()
  aligned_traces_raw_long = tidyr::pivot_longer(aligned_traces_raw_t,
                                            cols = colnames(aligned_traces_raw_t)[-1],
                                            names_to = "ID",
                                            values_to = "intensity")
  aligned_traces_raw_long$ID = as.numeric(aligned_traces_raw_long$ID)
  aligned_traces_raw_long$group = NA
  # assign groups
  aligned_traces_raw_long$group[which(aligned_traces_raw_long$ID %in% close_cells_IDs)] = "close"
  aligned_traces_raw_long$group[which(aligned_traces_raw_long$ID %in% between_cells_IDs)] = "mid"
  aligned_traces_raw_long$group[which(aligned_traces_raw_long$ID %in% far_cells_IDs)] = "far"
  
  #' summarize groups
  traces_raw_summary = aligned_traces_raw_long %>%
    dplyr::group_by(group, t) %>% 
    dplyr::summarise(mean = mean(intensity, na.rm = TRUE),
                     se = sd(intensity, na.rm = TRUE) / sqrt(n()),
                     .groups = "drop"
    )
  
  
  aligned_traces_long$ID = as.numeric(aligned_traces_long$ID)
  # add group field and set to NA
  aligned_traces_long$group = NA
  # assign groups
  aligned_traces_long$group[which(aligned_traces_long$ID %in% close_cells_IDs)] = "close"
  aligned_traces_long$group[which(aligned_traces_long$ID %in% between_cells_IDs)] = "mid"
  aligned_traces_long$group[which(aligned_traces_long$ID %in% far_cells_IDs)] = "far"
  
  #' summarize groups
  traces_summary = aligned_traces_long %>%
    dplyr::group_by(group, t) %>% 
    dplyr::summarise(mean = mean(intensity, na.rm = TRUE),
                     se = sd(intensity, na.rm = TRUE) / sqrt(n()),
                     .groups = "drop"
    )
    
  #' save long table containing traces
  
  if(saving){
  aligned_trx_path = file.path(rds_dir, paste0(filename, "_aligned_traces.rds"))
  saveRDS(aligned_traces_long, file = aligned_trx_path)
  
  aligned_trx_wide_path = file.path(rds_dir, paste0(filename, "_aligned_traces_wide.rds"))
  saveRDS(aligned_traces_t, file = aligned_trx_wide_path)
  
  aligned_trx_raw_path = file.path(rds_dir, paste0(filename, "_raw_aligned_traces.rds"))
  saveRDS(aligned_traces_raw_long, file = aligned_trx_raw_path)
  
  aligned_trx_raw_wide_path = file.path(rds_dir, paste0(filename, "_raw_aligned_traces_wide.rds"))
  saveRDS(aligned_traces_raw_t, file = aligned_trx_raw_wide_path)
  
  summary_raw_path = file.path(rds_dir, paste0(filename, "_summary_raw_trx.rds"))
  saveRDS(traces_raw_summary, file = summary_raw_path)
  
  summary_raw_path_path = file.path(rds_dir, paste0(filename, "_summary_trx.rds"))
  saveRDS(traces_summary, file = summary_raw_path_path)
  
  print(paste0("Average trace of ", filename,  ": Saved"))
  }
  
  #' #' get traces of groups
  #' close_traces = aligned_traces[, which(colnames(aligned_traces) %in% close_cells_IDs)]
  #' between_traces = aligned_traces[, which(colnames(aligned_traces) %in% between_cells_IDs)]
  #' far_traces = aligned_traces[, which(colnames(aligned_traces) %in% far_cells_IDs)]
  #' 
  #' #' remove cells with NA
  #' close_traces_clean = close_traces[, -c(unique(which(is.na(close_traces), arr.ind=TRUE)[,2]))]
  #' write.csv(close_traces_clean, "clean.csv")
 
  # make plots
    
  summary_plot = plot_group_traces(traces_summary)
  summary_raw_plot = plot_group_traces(traces_raw_summary)
    
  # close_plot_all = ggplot2::ggplot(data = aligned_traces_long)+
  #   geom_line(aes(x = t, y = intensity, group = ID))+
  #   scale_x_discrete(t, breaks = c(seq(12, 172, by = 24)))
  
  if(saving){
    savePlots(obj_to_save = list(distance_groups_trace_plot = summary_plot, distance_groups_raw_trace_plot = summary_raw_plot), 
              filename = filename, basepath = meanplot_dir, extension = "svg", p.width = 2300, p.height = 1390)}


  # MERGE PERIOD AND PARTICLE ANALYSIS ####
  
  # analysis to visualize circadian parameters in relation to the distance from plaques
  distance_table <- cbind(as.numeric(rownames(distance_matrix)), as.numeric(distance_matrix_median)) %>% `colnames<-`(c("ID", "distance"))#paste("Distance_f", frame, sep = "")))
  
  # create a median version of the distance table to merge with the period table
  merged_table <- merge(period_tbl, distance_table, by = "ID")
  
  # merge AUC and median values
  merged_table <- merge(merged_table, auc_results, by = "ID")
  merged_table <- merge(merged_table, median_vals, by = "ID")
  merged_table <- merge(merged_table, SD_vals, by = "ID")
  
  # Define distance groups
  merged_table$distance_group <- cut(merged_table$distance,
                                     breaks = c(-Inf, limit1, limit2, Inf),
                                     labels = c(paste("<=", limit1), paste(limit1, "-",limit2 ), paste(">", limit2)))
  
  merged_table$distance_group <- as.factor(merged_table$distance_group)
  merged_table$ID <- as.factor(merged_table$ID)
  
  # save RDS image of merged table
  merged_tbl_path = file.path(rds_dir, paste(filename, "_merged_tbl.rds", sep = ""))
  saveRDS(merged_table, file = merged_tbl_path)
  
  # ABSOLUTE CALCIUM LEVELS ####
  
  # make a 3D plot that represents the median level (x) with period (y) and distance from particle (z) 
  dist_AUC_per_plot = threeDplot(data = merged_table,
                                 x_var = distance,
                                 y_var = AUC,
                                 z_var = period,
                                 color_var = distance_group, 
                                 col_palette = color_palette,
                                 x_lab = "Distance (µm)",
                                 y_lab = "AUC (a.u.)", 
                                 z_lab = "Period (h)",
                                 eye = "pos2",
                                 zoom = 0.5,
                                 xrange = c(140,0),
                                 yrange = c(0, 100),
                                 zrange = c(21.85, 25.5))
  
  phase_AUC_per_plot = threeDplot(data = merged_table,
                                  x_var = phase_norm,
                                  y_var = AUC,
                                  z_var = period,
                                  color_var = distance_group, 
                                  col_palette = color_palette,
                                  x_lab = "Phase (h)",
                                  y_lab = "AUC (a.u.)", 
                                  z_lab = "Period (h)",
                                  eye = "pos2",
                                  zoom = 0.5,
                                  xrange = c(12, -12),
                                  yrange = c(0, 100),
                                  zrange = c(21.85, 25.5))
  
  # save 3D plot object for later visualization
  if(saving){
  dist_AUC_per_plot_path = file.path(rds_dir, paste0(filename, "_3D_dist_auc_per_plot.rds"))
  saveRDS(dist_AUC_per_plot, file = dist_AUC_per_plot_path)
  dist_AUC_per_plot_path2 = file.path(threeDplot_dir, paste0(filename, "_3D_dist_auc_per_plot.html"))
  saveWidget(dist_AUC_per_plot, dist_AUC_per_plot_path2, selfcontained = TRUE)
  
  phase_AUC_per_plot_path = file.path(rds_dir, paste0(filename, "_3D_phase_auc_per_plot.rds"))
  saveRDS(phase_AUC_per_plot, file = phase_AUC_per_plot_path)
  phase_AUC_per_plot_path2 = file.path(threeDplot_dir, paste0(filename, "_3D_phase_auc_per_plot.html"))
  saveWidget(phase_AUC_per_plot, phase_AUC_per_plot_path2, selfcontained = TRUE)
  }
  # save snapshot of 3D plot
  
  # WARNING this function requires the orca command line utility
  # to install use conda: conda install -c plotly plotly-orca
  # then set the location with Sys.setenv(ORCA_PATH = "path_to_orca_exe")
  # install.packages("processx")  # Required dependency
  # install.packages("remotes")
  # remotes::install_github("plotly/orca")
  # library(orca)
  # library(processx)
  # library(remotes)
  # plotly::orca(phase_AUC_per_plot, "3d_plot_dist_AUC_per.png")
  # plotly::kaleido(phase_AUC_per_plot, "3d_plot_dist_AUC_per.png")
  
  
  # CIRCADIAN PARAMS PLOTS ####
  
  # Period vs Distance
  period_meanplot <- metric_mean_plot(data = merged_table, metric = period, grouping_var = distance_group, plot_colors = color_palette, norm_test = FALSE, user_ylims = c(22.5, 24.5))
  
  # Amplitude vs Distance
  amplitude_meanplot <- metric_mean_plot(data = merged_table, metric = amplitude, grouping_var = distance_group, plot_colors = color_palette, norm_test = FALSE, user_ylims = c(0, 90))
  
  # Error vs Distance
  error_meanplot <- metric_mean_plot(data = merged_table, metric = error, grouping_var = distance_group, plot_colors = color_palette, norm_test = FALSE, user_ylims = c(0, 20))
  
  # AUC vs Distance
  AUC_meanplot <- metric_mean_plot(data = merged_table, metric = AUC, grouping_var = distance_group, plot_colors = color_palette, norm_test = FALSE, user_ylims = c(0, 100))
  
  # RAE vs Distance
  RAE_meanplot <- metric_mean_plot(data = merged_table, metric = RAE, grouping_var = distance_group, plot_colors = color_palette, norm_test = FALSE, user_ylims = c(0, 1))
  
  # Save plots
  if(saving){
    savePlots(obj_to_save = list(period_meanplot_w = period_meanplot, 
                                 amplitude_meanplot_w = amplitude_meanplot,
                                 error_meanboxplot_w = error_meanplot,
                                 AUC_meanplot_w = AUC_meanplot,
                                 RAE_meanplot_w = RAE_meanplot), filename = filename, basepath = meanplot_dir, extension = "svg", p.width = 1500, p.height = 800)}
  
  # Period vs Distance
  period_dotplot_trace <- ggplot(merged_table, aes(x = distance, y = period)) +
    #geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
    geom_smooth(aes(colour = distance_group), method = "lm", se =  TRUE, level = 0.95) +
    scale_color_manual(values = color_palette)+
    ylim(22.5, 24.5) +
    theme_minimal() +
    labs(title = "Period by Distance", x = "Distance from Particle (µm)", y = "Period (h)")
  
  # Amplitude vs Distance
  amplitude_dotplot_trace <- ggplot(merged_table, aes(x = distance, y = amplitude)) +
    #geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
    geom_smooth(aes(colour = distance_group), method = "lm", se = TRUE, level = 0.95) +
    scale_color_manual(values = color_palette)+
    ylim(0, 130) +
    theme_minimal() +
    labs(title = "Amplitude by Distance", x = "Distance from Particle (µm)", y = "Amplitude (A.U.)")
  
  
  # Error vs Distance
  error_dotplot_trace <- ggplot(merged_table, aes(x = distance, y = error)) +
    #geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
    geom_smooth(aes(colour = distance_group), method = "lm", se =  TRUE, level = 0.95) +
    theme_minimal() +
    scale_color_manual(values = color_palette)+
    ylim(0, 30) +
    labs(title = "Error by Distance", x = "Distance from Particle (µm)", y = "Error (A.U.)")
  
  # AUC vs Distance
  AUC_dotplot_trace <- ggplot(merged_table, aes(x = distance, y = AUC)) +
    #geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
    geom_smooth(aes(colour = distance_group), method = "lm", se =  TRUE, level = 0.95) +
    theme_minimal() +
    scale_color_manual(values = color_palette)+
    ylim(0, 100) +
    labs(title = "AUC by Distance", x = "Distance from Particle (µm)", y = "AUC (A.U.)")
  
  # RAE vs Distance
  RAE_dotplot_trace <- ggplot(merged_table, aes(x = distance, y = RAE)) +
    #geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
    geom_smooth(aes(colour = distance_group), method = "lm", se =  TRUE, level = 0.95) +
    theme_minimal() +
    scale_color_manual(values = color_palette)+
    ylim(0, 1) +
    labs(title = "RAE by Distance", x = "Distance from Particle (µm)", y = "RAE (A.U.)")
  
  # Save plots
  if(savingtRACE){
    savePlots(obj_to_save = list(period_dotplot_trace = period_dotplot_trace,
                                 amplitude_dotplot_trace = amplitude_dotplot_trace,
                                 error_dotplot_trace = error_dotplot_trace,
                                 AUC_dotplot_trace = AUC_dotplot_trace,
                                 RAE_dotplot_trace = RAE_dotplot_trace),
              filename = filename, basepath = lineplot_dir, extension = "svg", p.width = 1800, p.height = 1000)}
  
  # PHASE PLOTS ####
  
  # Phase vs distance
  groups <- unique(merged_table$distance_group)
  # add raileygh plots
  sub1 <- merged_table[merged_table$distance_group == groups[1], ]
  sub2 <- merged_table[merged_table$distance_group == groups[2], ]
  sub3 <- merged_table[merged_table$distance_group == groups[3], ]
  
  far_circ <- circular_plot(data = sub1$phase_rad, path = file.path(circplot_dir, paste0(filename, "_far_cells.svg")), saving = saving)
  mid_circ <- circular_plot(data = sub2$phase_rad, path = file.path(circplot_dir, paste0(filename, "_mid_cells.svg")), saving = saving)
  close_circ <- circular_plot(data = sub3$phase_rad, path = file.path(circplot_dir, paste0(filename, "_close_cells.svg")), saving = saving)
  
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
  circStats_plot(newrecord, path = file.path(circplot_dir, paste0(filename, "_phases_summary.svg")), colors = color_palette)
  
  
}

# OPS AFTER CYCLE IS FINISHED ####

#save circ_summary table
write.csv(circ_summary, file = file.path(wd, "circular_stats.csv"), row.names = TRUE)


# GROUP PLOTS IN FOLDER -----
pull_plots(wd, plot.name = "phases_summary.svg", dir.name = "phases_summary", file.names = filenames, second_fold = "circ_plots")
pull_plots(wd, plot.name = "groups_cells_plot.svg", dir.name = "group_maps", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_amplitude.svg", dir.name = "spatial_amp", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_period.svg", dir.name = "spatial_per", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_error.svg", dir.name = "spatial_err", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_fred.svg", dir.name = "spatial_fred", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_green.svg", dir.name = "spatial_green", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_green_AUC.svg", dir.name = "spatial_green_AUC", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_green_SD.svg", dir.name = "spatial_green_SD", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_RAE.svg", dir.name = "spatial_RAE", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_phases_norm.svg", dir.name = "spatial_phase_norm", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "spatial_phases_circ.svg", dir.name = "spatial_phase", file.names = filenames, second_fold = "spatial_maps")
pull_plots(wd, plot.name = "period_dotplot_trace.svg", dir.name = "period_dotplot", file.names = filenames, second_fold = "line_plots")
pull_plots(wd, plot.name = "amplitude_dotplot_trace.svg", dir.name = "amplitude_dotplot", file.names = filenames, second_fold = "line_plots")
pull_plots(wd, plot.name = "error_dotplot_trace.svg", dir.name = "error_dotplot", file.names = filenames, second_fold = "line_plots")
pull_plots(wd, plot.name = "AUC_dotplot_trace.svg", dir.name = "AUC_dotplot", file.names = filenames, second_fold = "line_plots")
pull_plots(wd, plot.name = "RAE_dotplot_trace.svg", dir.name = "RAE_dotplot", file.names = filenames, second_fold = "line_plots")
pull_plots(wd, plot.name = "distance_groups_trace_plot.svg", dir.name = "summary_traces", file.names = filenames, second_fold = "mean_plots")
pull_plots(wd, plot.name = "period_meanplot_w.svg", dir.name = "period_meanplot", file.names = filenames, second_fold = "mean_plots")
pull_plots(wd, plot.name = "amplitude_meanplot_w.svg", dir.name = "amplitude_meanplot", file.names = filenames, second_fold = "mean_plots")
pull_plots(wd, plot.name = "error_meanboxplot_w.svg", dir.name = "error_meanplot", file.names = filenames, second_fold = "mean_plots")
pull_plots(wd, plot.name = "AUC_meanplot_w.svg", dir.name = "AUC_meanplot", file.names = filenames, second_fold = "mean_plots")
pull_plots(wd, plot.name = "RAE_meanplot_w.svg", dir.name = "RAE_meanplot", file.names = filenames, second_fold = "mean_plots")
pull_plots(wd, plot.name = "distance_groups_raw_trace_plot.svg", dir.name = "summary_traces_raw", file.names = filenames, second_fold = "mean_plots")
pull_plots(wd, plot.name = "3D_dist_auc_per_plot.html", dir.name = "3D_dist_AUC_per", file.names = filenames, second_fold = "3D_plots")
pull_plots(wd, plot.name = "3D_phase_auc_per_plot.html", dir.name = "3D_phase_AUC_per", file.names = filenames, second_fold = "3D_plots")
# extract files from each file folder
# pull_plots(wd, plot.name = "thresholded.tiff", dir.name = "thresh_imgs", file.names = filenames, mainfold = TRUE)
# pull_plots(wd, plot.name = "SCN_nuclei.roi", dir.name = "thresh_imgs", file.names = filenames, mainfold = TRUE)
#' TODO create a remove_file function




#' 
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

