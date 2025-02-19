#' script to collate the data from the different samples and run comparisons

#' objects to import and collate
#'  grid_vals or distance matrix
#'  period table

# install and load required packages
pkg_req <- c("sf", "RImageJROI", "ggplot2", "dplyr", "pracma", "minpack.lm", 
             "magrittr", "stringr", "circular", "svglite", "astsa", "pdftools", 
             "scico", "ggpubr", "rstatix", "plotly", "htmlwidgets", "reticulate")

# Install and load packages if not already installed
for (package in pkg_req) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}


 # FUNCTIONS ---------------------------------------------------------

source(file = r"(C:\Users\mf420\Documents\GitHub\Spatial_analysis\Tau_spatial_analysis_foo.R)")

source(file = r"(C:\Users\mf420\Documents\GitHub\Spatial_analysis\collate_Tau_spatial_analysis_foo.R)")

 # EXECUTION ---------------------------------------------------------

# to clean the environment preserving functions
rm(list = setdiff(ls(), lsf.str()))

 # sections to run ---------------------------------------------------------
EXCLUDE_FILES = TRUE
MULTI_FOLDER = FALSE
CELLS_COUNTING = TRUE
DISTANCE_GROUPS = TRUE
MERGE_TABLES = TRUE
ALIGNED_TRACES_PERIOD = TRUE
NESTED_PLOTS = TRUE
CORRELATIONS = TRUE
CELLS_RATIO = TRUE
LINEPLOTS_DISTANCE = TRUE
ALIGN_TRACES_COMPARE = TRUE
TOPBOTTOM = FALSE
PHASE_DIST_GROUPS = FALSE
VECTOR_SPREAD = FALSE
THREE_D_PLOTS = TRUE

 # IMPORTING ---------------------------------------------------------
#' select home dir where to pull names and data from
week = "wk2"
wd = r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Proj_Tau\Tau_an_dev\wk5_hemislices\combined_mid)"
if(!dir.exists(wd)){
  dir.create(wd)
}

wd1 = r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Proj_Tau\Tau_an_dev\wk5_hemislices\Left_size_4_trsh)"

files1 = list.files(path = wd1, pattern = ".tif*$")
filenames1 = stringr::str_remove(files1, pattern = ".tif*")
foldernames1 = file.path(wd1, paste(filenames1, "_results", sep = ""))

if(MULTI_FOLDER){
wd2 = r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Proj_Tau\Tau_an_dev\wk5_hemislices\Right_size_4_thresh)"

files2 = list.files(path = wd2, pattern = ".tif*$")
filenames2 = stringr::str_remove(files2, pattern = ".tif*")
foldernames2 = file.path(wd2, paste(filenames2, "_results", sep = ""))
}

files = c(files1)#, files2)
filenames = c(filenames1)#, filenames2)
foldernames = c(foldernames1)#, foldernames2)

dataset_list <- list(grid.vals = list(),
                     dist.mtx = list(),
                     period.tbl = list(),
                     align.trx = list())

threeD_plot_list <- list(distance = list(),
                         phase = list())

IMPORT_PLOTS <- TRUE

missing_matrices <- vector()
missing_periods <- vector()
missing_merg_tbls <- vector()
missing_traces <- vector()
missing_traces_wide <- vector()
missing_traces_raw <- vector()
missing_traces_raw_wide <- vector()
missing_puncta <- vector()

if(EXCLUDE_FILES){
#excluded data and filter filenames and foldernames vector
excluded <- c(
  "LEFT_A1_FRED", "LEFT_A4_FRED", "LEFT_A6_MUT", "LEFT_D5_MUT",
  # "RIGHT_B3_MUT", "RIGHT_C6_MUT", #low
  # "LEFT_D2_MUT", "LEFT_D3_FRED",
  # "RIGHT_B1_FRED", "RIGHT_B4_FRED", "RIGHT_B6_MUT", #mid
  "LEFT_B1_FRED", "LEFT_B3_MUT", "LEFT_B4_FRED", "LEFT_C2_MUT", "LEFT_C3_FRED",
  # "RIGHT_C1_FRED"  #high
  "LEFT_D6_FRED"
  )

filenames <- setdiff(filenames, excluded)
foldernames <- foldernames[!str_detect(foldernames, str_c(excluded, collapse = "|"))]
}
 # DATASET_LIST ---------------------------------------------------------
# period tbls, distance matrices, merged tables, and aligned traces
for (i in seq_len(length(filenames))){
  foldername = foldernames[i]
  filename = filenames[i]
  print(filename)
  #' import period tbls
  period_tbl_path = file.path(foldername, paste(filename, "_period_tbl.rds", sep = ""))
  if(file.exists(period_tbl_path)){
    period_tbl <- readRDS(period_tbl_path)
    dataset_list$period.tbl[[filename]] <- period_tbl
  }else{
    print(paste(filename, " period file missing."))
    missing_periods <- c(missing_periods, filename)
  }
  
  #' import distance matrices
  distance_mtx_path = file.path(foldername, paste(filename, "_dist_mtx.rds", sep = ""))
  if(file.exists(distance_mtx_path)){
    distance_matrix <- readRDS(distance_mtx_path)
    dataset_list$dist.mtx[[filename]] <- distance_matrix
  }else{
    print(paste(filename, " matrix file missing."))
    missing_matrices <- c(missing_matrices, filename)
  }
  
  #' import merged tables
  merg_tbl_path = file.path(foldername, paste(filename, "_merged_tbl.rds", sep = ""))
  if(file.exists(merg_tbl_path)){
    merged_table <- readRDS(merg_tbl_path)
    dataset_list$merg.tbl[[filename]] <- merged_table
  }else{
    print(paste(filename, " merged table file missing."))
    missing_merg_tbls <- c(missing_merg_tbls, filename)
  }
  
  #' import aligned traces
  align_trx_path = file.path(foldername, paste(filename, "_aligned_traces.rds", sep = ""))
  if(file.exists(align_trx_path)){
    aligned_trx <- readRDS(align_trx_path)
    dataset_list$align.trx[[filename]] <- aligned_trx
  }else{
    print(paste(filename, " aligned traces file missing."))
    missing_traces <- c(missing_traces, filename)
  }
  
  #' import aligned traces wide
  align_trx_wide_path = file.path(foldername, paste(filename, "_aligned_traces_wide.rds", sep = ""))
  if(file.exists(align_trx_wide_path)){
    aligned_trx_wide <- readRDS(align_trx_wide_path)
    dataset_list$align.trx.wide[[filename]] <- aligned_trx_wide
  }else{
    print(paste(filename, " aligned traces wide file missing."))
    missing_traces_wide <- c(missing_traces_wide, filename)
  }
  
  #' import aligned traces raw
  align_trx_raw_path = file.path(foldername, paste(filename, "_raw_aligned_traces.rds", sep = ""))
  if(file.exists(align_trx_raw_path)){
    align_trx_raw <- readRDS(align_trx_raw_path)
    dataset_list$align.trx.raw[[filename]] <- align_trx_raw
  }else{
    print(paste(filename, " raw aligned traces file missing."))
    missing_traces_raw <- c(missing_traces_raw, filename)
  }
  
  #' import aligned traces raw wide
  align_trx_raw_wide_path = file.path(foldername, paste(filename, "_raw_aligned_traces_wide.rds", sep = ""))
  if(file.exists(align_trx_raw_wide_path)){
    aligned_trx_raw_wide <- readRDS(align_trx_raw_wide_path)
    dataset_list$align.trx.raw.wide[[filename]] <- aligned_trx_raw_wide
  }else{
    print(paste(filename, " raw aligned traces wide file missing."))
    missing_traces_raw_wide <- c(missing_traces_raw_wide, filename)
  }
  
  #' Import summary of SCN_summary puncta
  SCN_puncta_summary_path = file.path(foldername, paste("SCN_Summary_puncta_", filename, ".csv", sep = ""))
  if(file.exists(SCN_puncta_summary_path)){
    SCN_puncta_summary <- read.csv (SCN_puncta_summary_path, header = TRUE)
    dataset_list$SCN_puncta_summary[[filename]] <- SCN_puncta_summary[, -1]
  }else{
    print(paste(filename, " SCN summary puncta file missing."))
    missing_puncta <- c(missing_puncta, filename)
  }
  
  #' Import 3D plot of distance period, auc
  if(IMPORT_PLOTS){
  dist_per_AUC_path = file.path(foldername, paste(filename, "_3D_dist_auc_per_plot.rds", sep = ""))
  if(file.exists(dist_per_AUC_path)){
    dist_per_AUC <- readRDS(dist_per_AUC_path)
    threeD_plot_list$distance[[filename]] <- dist_per_AUC
  }else{
    print(paste(filename, " 3D plot of Distance, period and AUC missing."))
  }
  
  #' Import 3D plot of phase period, auc
  phase_per_AUC_path = file.path(foldername, paste(filename, "_3D_phase_auc_per_plot.rds", sep = ""))
  if(file.exists(phase_per_AUC_path)){
    phase_per_AUC <- readRDS(phase_per_AUC_path)
    threeD_plot_list$phase[[filename]] <- phase_per_AUC
  }else{
    print(paste(filename, " 3D plot of Phase, period and AUC missing."))
  }
  }
  
  #clean files outside the list
  rm(list = c("period_tbl", "distance_matrix", "merged_table", "aligned_trx", "aligned_trx_wide", 
              "align_trx_raw", "aligned_trx_raw_wide", "SCN_puncta_summary", "dist_per_AUC", "phase_per_AUC"))
}

print("Datasets added!")

# create a folder to save rds files
rds_fold = file.path(wd, "rds_files")

if(!dir.exists(rds_fold)){
  dir.create(rds_fold)
}
# save dataset_list for retrieval
saveRDS(dataset_list, file = file.path(rds_fold, paste0(week, "_dataset.rds")))

 #### CELLS COUNTING ---------------------------------------------------------
if(CELLS_COUNTING){
  cells_count_fold = file.path(wd, "Cell_counting")
  if(!dir.exists(cells_count_fold)){
    dir.create(cells_count_fold)
  }
SCN_puncta_list <- dataset_list$SCN_puncta_summary
 # read tables of summary ROI, compute median and add to table
summarize_tables <- function(table_list) {
  # Get table names
  table_names <- names(table_list)
  
  # Compute median values for each table
  summary_list <- lapply(seq_along(table_list), function(i) {
    table <- table_list[[i]]
    
    # Ensure only numeric columns are considered
    numeric_cols <- table[sapply(table, is.numeric)]
    
    # Compute column-wise median
    medians <- apply(numeric_cols, 2, median, na.rm = TRUE)
    
    # Convert to data frame and add table name
    data.frame(Table = table_names[i], t(medians), check.names = FALSE)
  })
  
  # Combine summaries into a single data frame
  summary_table <- do.call(rbind, summary_list)
  
  return(summary_table)
}

SCN_puncta_res <- summarize_tables(SCN_puncta_list)

write.csv(SCN_puncta_res, file = file.path(cells_count_fold, "SCN_particle_count.csv"))

SCN_puncta_res$exclude <- FALSE
SCN_puncta_res <- SCN_puncta_res %>% mutate(exclude = ifelse(X.Area < 0.5, TRUE, exclude))
}
 #### DISTANCE GROUPS ---------------------------------------------------------
if(DISTANCE_GROUPS){
#' compare groups across conditions for only two sample being compared

limit1 = 8.5
limit2 = 34
#limit3 = 40

#' interrogate the matrix to see which cells are at less than
close_cells_matrix = lapply(dataset_list$dist.mtx,  function(x){x<=limit1})
between_cells_matrix = lapply(dataset_list$dist.mtx,  function(x){x>limit1 & x <= limit2})
far_cells_matrix = lapply(dataset_list$dist.mtx,  function(x){x>limit2})# & x <= limit3})
# veryfar_cells_matrix = lapply(dataset_list$dist.mtx,  function(x){x>limit3})

color_palette = c("#B23A48", "#2F9C95", "#B6DC76")#, "#E0D3DE")
}
 #### merge all tables from individual spatial analysis ---------------------------------------------------------
if(MERGE_TABLES){
  plot_folder = file.path(wd, "spatial_analysis")
  
  if(!dir.exists(plot_folder)){
    dir.create(plot_folder)
  }
  
  #' create list of distances and extract median value
  distance_list <- lapply(dataset_list$dist.mtx, function(matrix){
    distance_table = cbind(as.numeric(rownames(matrix)), round(matrixStats::rowMedians(matrix), 3)) %>% `colnames<-`(c("ID", "distance"))
  })
  
#' merge tables adding column for normalized phase, treatment and well
merged_list <- setNames(lapply(names(dataset_list$merg.tbl), function(name) {
  df <- dataset_list$merg.tbl[[name]]
  # calculate the mean of the phase_rad values then normalize and add new column
  # phase_norm <- normalize_and_convert(df$phase_rad)
  # df$phase_norm <- phase_norm
  # Extract part of the name after the "_"
  parts <- str_split(name, "_")[[1]]
  well = paste(parts[1], parts[2], sep = "_")
  treatment <- parts[3]  #sub("^[^_]*_", "", name)
  # Add suffix as a new column
  df$well <- as.factor(well)
  df$treatment <- as.factor(treatment)
  return(df)
}), names(distance_list))

# save merged list object
saveRDS(merged_list, file = file.path(rds_fold, paste0(week, "_merged_list.rds")))

#' merge all the tables in the merged list into a single table
merged_table = data.table::rbindlist(merged_list, idcol = "sample")

#' save merged table
saveRDS(merged_table, file = file.path(rds_fold, paste0(week, "_merged_table.rds")))
}
 #### Re-compute period based on aligned traces ---------------------------------------------------------
if(ALIGNED_TRACES_PERIOD){
# get data from dataset_list
aligned_traces_raw_list <- dataset_list$align.trx.raw.wide

# compute period
period_tbl_list <- setNames(lapply(names(aligned_traces_raw_list), function(name) {   #browser()

  print(name)
  df <- aligned_traces_raw_list[[name]]
  # trim table to cut out NA vals
  df_trim <- trim_tbl(df, rows = 24)
  
  # interpolate NA vals
  for(i in seq(2, length(colnames(df_trim)))){
    ts <- df_trim[,i]
    # interpolate missing values before detrending
    ts_clean <- imputeTS::na_interpolation(ts, option = "spline")
    # add it back to the table
    if(length(ts) == length(ts_clean)){
      df_trim[ ,i] <- ts_clean
    }
  }
  #calculate period
  period_list <- computePeriod(df = df_trim, prep_tbl = TRUE, preprocess = FALSE, transp = FALSE, add_t = FALSE)
  
  period_tbl <- period_list$period_table
  
  return(period_tbl)
  
}), names(aligned_traces_raw_list))

 #### create merged list of periods and distances on aligned detrended traces ---------------------------------------------------------

#' merge distance table with period table, then assign name to list element
merged_list_align <- setNames(
  lapply(names(distance_list), function(x, dataset_list, distance_list, period_tbl_list){
    period_tbl <- period_tbl_list[[x]]
    distance_tbl <- distance_list[[x]]
    merge(period_tbl, distance_tbl, by = "ID")
  }, dataset_list, distance_list, period_tbl_list), names(distance_list))


# Define distance groups in merged tables list
merged_list_align <- lapply(merged_list_align, function(df) {  # browser()
  
  df$distance_group <- cut(df$distance,
                           breaks = c(-Inf, limit1, limit2, Inf),# limit3, Inf),
                           labels = c(paste("<=", limit1), paste(limit1, "-",limit2 ),  paste(">", limit2)))# paste(limit2, "-",limit3), paste(">", limit3)))
  df$distance_group <- as.factor(df$distance_group)
  df$ID <- as.factor(df$ID)
  return(df)
})

#' merge tables adding column for treatment and well
merged_list_align <- setNames(lapply(names(merged_list_align), function(name) {
  df <- merged_list_align[[name]]
  # Extract part of the name after the "_"
  parts <- str_split(name, "_")[[1]]
  well = paste(parts[1], parts[2], sep = "_")
  treatment <- parts[3]  #sub("^[^_]*_", "", name)
  # Add suffix as a new column
  df$well <- as.factor(well)
  df$treatment <- as.factor(treatment)
  return(df)
}), names(distance_list))

plot_folder = file.path(wd, "spatial_analysis_all")

if(!dir.exists(plot_folder)){
  dir.create(plot_folder)
}

# save merged list object
saveRDS(merged_list_align, file = file.path(wd, paste0(week, "_merged_list_aligned.rds")))

#' merge all the tables in the merged list into a single table
merged_table_align = data.table::rbindlist(merged_list_align, idcol = "sample")
}
 #### analysis to visualize circadian parameters in relation to the distance from plaques ---------------------------------------------------------
 #### nested anova plots of period, amplitude, error ---------------------------------------------------------
if(NESTED_PLOTS){
  nest_plot_fold = file.path(wd, "Nested_plots")
  if(!dir.exists(nest_plot_fold)){
    dir.create(nest_plot_fold)
  }
#' Create nested mean plots that summarise all the samples 
   
  period_nest <- nested_anova_plot(
    data = merged_table,
    metric = period,
    distance_var = distance_group,
    treatment_var = treatment,
    plot_colors = c("black", "red"),
    ylabel = "Period (h)",
    xlabel = "Distance group (px)",
    nudge = 0.5,
    norm_test = FALSE,
    custom_levels = c("FRED", "MUT")
  )
  period_mean_plot_nest <- period_nest$plot
  
   
  amplitude_nest <- nested_anova_plot(
    data = merged_table,
    metric = amplitude,
    distance_var = distance_group,
    treatment_var = treatment,
    plot_colors = c("black", "red"),
    ylabel = "Amplitude (A.U.)",
    xlabel = "Distance group (px)",
    nudge = 8,
    step = 0.6,
    norm_test = FALSE,
    custom_levels = c("FRED", "MUT")
  )
  amplitude_mean_plot_nest <- amplitude_nest$plot
  
  # error (residual error) is calculated during period computation as follow: 
  # error <- sqrt(mean((data$value - fitted_values)^2))
  error_nest <- nested_anova_plot(
    data = merged_table,
    metric = error,
    distance_var = distance_group,
    treatment_var = treatment,
    plot_colors = c("black", "red"),
    ylabel = "Error (A.U.)",
    xlabel = "Distance group (px)",
    nudge = 4,
    step = 0.4,
    norm_test = FALSE,
    custom_levels = c("FRED", "MUT")
  )
  error_mean_plot_nest <- error_nest$plot
  
  RAE_nest <- nested_anova_plot(
    data = merged_table,
    metric = RAE,
    distance_var = distance_group,
    treatment_var = treatment,
    plot_colors = c("black", "red"),
    ylabel = "RAE (A.U.)",
    xlabel = "Distance group (px)",
    nudge = 4,
    step = 0.4,
    norm_test = FALSE,
    custom_levels = c("FRED", "MUT")
  )
  RAE_mean_plot_nest <- RAE_nest$plot
  
  AUC_nest <- nested_anova_plot(
    data = merged_table,
    metric = AUC,
    distance_var = distance_group,
    treatment_var = treatment,
    plot_colors = c("black", "red"),
    ylabel = "AUC (A.U.)",
    xlabel = "Distance group (px)",
    nudge = 4,
    step = 0.4,
    norm_test = FALSE,
    custom_levels = c("FRED", "MUT")
  )
  AUC_mean_plot_nest <- AUC_nest$plot
  
  write.csv(period_nest$well_stats, file = file.path(wd, paste0("period_wells_stats.csv")))
  savePlots(obj_to_save = list(period_mean_plot_nest = period_mean_plot_nest, amplitude_mean_plot_nest = amplitude_mean_plot_nest, error_mean_plot_nest = error_mean_plot_nest,
                               RAE_mean_plot_nest = RAE_mean_plot_nest, AUC_mean_plot_nest = AUC_mean_plot_nest), basepath = nest_plot_fold, extension = "svg", p.width = 1500, p.height = 900)
}
 #### Correlation of AUC and distance ---------------------------------------------------------
if(CORRELATIONS){

corr_folder = file.path(wd, "correlation")

if(!dir.exists(corr_folder)){
  dir.create(corr_folder)
}

 # separate data by treatment and well
 data_split <- merged_table %>%
  group_by(treatment, well, distance_group) %>%
  select(AUC, distance) %>%
  filter(!is.na(AUC), !is.na(distance))  # Remove NA values to avoid correlation issues

  # calculate the correlation for each group (treatment and well combination)
  correlations_AUC <- data_split %>%
    group_by(treatment, well, distance_group) %>%
    summarise(correlation = cor(distance, AUC), .groups = "drop")
 
  # visualize data and correlations
  linear_corr_plot <- ggplot(data_split, aes(y = AUC, x = distance, z = distance_group, color = treatment)) +
    # geom_point(size =1) +
    facet_wrap(distance_group~well) +  # Separate by well
    xlim(c(0,50))+
    geom_smooth(method = "lm", se = FALSE, aes(group = treatment)) +
    labs(title = "Correlation between AUC and Distance by Treatment and Well",
         x = "Distance",
         y = "AUC") +
    theme_minimal()
  
  correlation_plot <- ggplot(correlations_AUC)+
    geom_point(aes(x = treatment, y = correlation, color = treatment))+
    facet_wrap(~distance_group)+
    theme_pubclean()+
    theme(legend.position = "none")+
    scale_color_manual(values = c("black", "red"))
  
  # save table and plots
  write.csv(correlations_AUC, file = file.path(corr_folder, paste0("AUC_distance_corr.csv")))
  savePlots(obj_to_save = list(AUC_dist_corr = correlation_plot), 
            basepath = corr_folder, extension = "svg", p.width = 1200, p.height = 1020)
  savePlots(obj_to_save = list(AUC_dist_linear_corr = linear_corr_plot), 
            basepath = corr_folder, extension = "svg", p.width = 7200, p.height = 3200)
  
  #### Correlation of SD and distance ---------------------------------------------------------
  corr_folder = file.path(wd, "correlation")
  
  if(!dir.exists(corr_folder)){
    dir.create(corr_folder)
  }
  
  # separate data by treatment and well
  data_split <- merged_table %>%
    group_by(treatment, well, distance_group) %>%
    select(SD_green, distance) %>%
    filter(!is.na(SD_green), !is.na(distance))  # Remove NA values to avoid correlation issues
  
  # calculate the correlation for each group (treatment and well combination)
  correlations_SD <- data_split %>%
    group_by(treatment, well, distance_group) %>%
    summarise(correlation = cor(distance, SD_green), .groups = "drop")
  
  # visualize data and correlations
  linear_corr_plot <- ggplot(data_split, aes(y = SD_green, x = distance, z = distance_group, color = treatment)) +
    # geom_point(size =1) +
    facet_wrap(distance_group~well) +  # Separate by well
    xlim(c(0,50))+
    geom_smooth(method = "lm", se = FALSE, aes(group = treatment)) +
    labs(title = "Correlation between SD and Distance by Treatment and Well",
         x = "Distance",
         y = "SD") +
    theme_minimal()
  
  correlation_plot <- ggplot(correlations_SD)+
    geom_point(aes(x = treatment, y = correlation, color = treatment))+
    facet_wrap(~distance_group)+
    theme_pubclean()+
    theme(legend.position = "none")+
    scale_color_manual(values = c("black", "red"))
  
  # save table and plots
  write.csv(correlations_SD, file = file.path(corr_folder, paste0("SD_distance_corr.csv")))
  savePlots(obj_to_save = list(SD_dist_corr = correlation_plot), 
            basepath = corr_folder, extension = "svg", p.width = 1200, p.height = 1020)
  savePlots(obj_to_save = list(SD_dist_linear_corr = linear_corr_plot), 
            basepath = corr_folder, extension = "svg", p.width = 7200, p.height = 3200)
}
 #### RATIO OF CELLS IN EACH AREA ---------------------------------------------------------
#' merge together all tables containing period and spatial data 
if(CELLS_RATIO){
  ratios_fold = file.path(wd, "meta_ratios")
  
  if(!dir.exists(ratios_fold)){
    dir.create(ratios_fold)
  }
  
merged_tbl_path = file.path(plot_folder, "merged_table_all_cells.csv")
write.csv(merged_table, merged_tbl_path)

meta_merged_table <-
  merged_table %>%
  group_by(sample, distance_group, treatment) %>%
  summarise(n = n(),
            period = round(mean(period), 2),
            amplitude = round(mean(amplitude), 2))

#calculate ratio of cells in each distance group
meta_merged_ratio <- meta_merged_table
meta_merged_ratio$ratio <- with(meta_merged_table, n / ave(n, sample, treatment, FUN = sum)) %>% round(., 3)
  
meta_color <- c("black", "red")

meta_plot_ratios_period_close = meta_plot_draw(data = meta_merged_ratio,
                                              filter = "<= 3",
                                              x = ratio, y = period,
                                              fill = treatment, colour = treatment, col_scheme = c("black", "red"),
                                              shape = treatment, linetype = treatment, user_ylims = c(21.5, 25), minor_grid_y = TRUE,
                                              title = "Period by cell fraction - close", xtitle = "Cell fraction", ytitle = "Period (h)")

meta_plot_ratios_period_mid = meta_plot_draw(data = meta_merged_ratio,
                                             filter = "3 - 12",
                                             x = ratio, y = period,
                                             fill = treatment, colour = treatment, col_scheme = c("black", "red"),
                                             shape = treatment, linetype = treatment, user_ylims = c(21.5, 25), minor_grid_y = TRUE,
                                             title = "Period by cell fraction - mid", xtitle = "Cell fraction", ytitle = "Period (h)")

meta_plot_ratios_period_far = meta_plot_draw(data = meta_merged_ratio,
                                             filter = "> 12",
                                             x = ratio, y = period,
                                             fill = treatment, colour = treatment, col_scheme = c("black", "red"),
                                             shape = treatment, linetype = treatment, user_ylims = c(21.5, 25), minor_grid_y = TRUE,
                                             title = "Period by cell fraction - far", xtitle = "Cell fraction", ytitle = "Period (h)")

meta_plot_ratios_amplitude_close = meta_plot_draw(data = meta_merged_ratio,
                                                  filter = "<= 3",
                                                  x = ratio, y = amplitude,
                                                  fill = treatment, colour = treatment, col_scheme = c("black", "red"),
                                                  shape = treatment, linetype = treatment, user_ylims = c(-40, 165),
                                                  title = "Amplitude by cell fraction - close", xtitle = "Cell fraction", ytitle = "Period (h)")

meta_plot_ratios_amplitude_mid = meta_plot_draw(data = meta_merged_ratio,
                                                  filter = "3 - 12",
                                                  x = ratio, y = amplitude,
                                                  fill = treatment, colour = treatment, col_scheme = c("black", "red"),
                                                  shape = treatment, linetype = treatment, user_ylims = c(-40, 165),
                                                  title = "Amplitude by cell fraction - mid", xtitle = "Cell fraction", ytitle = "Period (h)")

meta_plot_ratios_amplitude_far = meta_plot_draw(data = meta_merged_ratio,
                                                  filter = "> 12",
                                                  x = ratio, y = amplitude,
                                                  fill = treatment, colour = treatment, col_scheme = c("black", "red"),
                                                  shape = treatment, linetype = treatment, user_ylims = c(-40, 165),
                                                  title = "Amplitude by cell fraction - far", xtitle = "Cell fraction", ytitle = "Period (h)")



savePlots(obj_to_save = list(meta_plot_ratios_period_close = meta_plot_ratios_period_close,
                             meta_plot_ratios_period_mid = meta_plot_ratios_period_mid,
                             meta_plot_ratios_period_far = meta_plot_ratios_period_far,
                             meta_plot_ratios_amplitude_close = meta_plot_ratios_amplitude_close,
                             meta_plot_ratios_amplitude_mid = meta_plot_ratios_amplitude_mid,
                             meta_plot_ratios_amplitude_far = meta_plot_ratios_amplitude_far), basepath = ratios_fold, extension = "svg", p.width = 1050, p.height = 1300)

#' plot representing how many cells are present in every distance area per sample

# Create dodged x-coordinates
meta_merged_table_dodge <- meta_merged_table %>%
  group_by(distance_group, treatment) %>%
  mutate(dodged_x = as.numeric(distance_group) + 
           (as.numeric(as.factor(treatment)) - 1.5) * 0.4 / 2)

# Plot using manually dodged x-coordinates
meta_plot <- ggplot() +
  # Points with adjusted x-coordinates
  geom_point(
    data = meta_merged_table_dodge,
    aes(x = dodged_x, y = n, color = treatment, alpha = 0.60),
    size = 3
  ) +
  # Lines with adjusted x-coordinates
  geom_line(
    data = meta_merged_table_dodge,
    aes(x = dodged_x, y = n, color = treatment, group = sample)
  ) +
  # Color scale
  scale_color_manual(values = c("black", "red")) +
  # Labels
  labs(
    title = "Cell in Every Distance Area",
    x = "Distance Group",
    y = "Cell Number"
  ) +
  # Restore original x-axis labels
  scale_x_continuous(
    breaks = unique(as.numeric(meta_merged_table$distance_group)),
    labels = unique(meta_merged_table$distance_group)
  ) +
  scale_alpha_identity()+
  theme_minimal()

savePlots(obj_to_save = list(meta_plot = meta_plot), basepath = ratios_fold, extension = "svg", p.width = 1050, p.height = 1300)
}
 #### LINEPLOT OF PERIOD/AMPLITUDE VS DISTANCE ---------------------------------------------------------
if(LINEPLOTS_DISTANCE){
#' plot showing period of all cells according to distance from particle
  lineplots_fold = file.path(wd, "lineplots_distance")
  
  if(!dir.exists(lineplots_fold)){
    dir.create(lineplots_fold)
  }
period_dotplot_trace_allcells <- ggplot()+#merged_table, aes(x = distance, y = period)) +
  # geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
  geom_smooth(data = merged_table[merged_table$treatment == "FRED", ], aes(x = distance, y = period, colour = distance_group, group = distance_group, linetype = treatment), method = "loess", se =  TRUE, level = 0.95) +
  geom_smooth(data = merged_table[merged_table$treatment == "MUT", ], aes(x = distance, y = period, colour = distance_group, group = distance_group, linetype = treatment), method = "loess", se =  TRUE, level = 0.95) +
  scale_color_manual(values = color_palette)+
  ylim(23, 24.5) +
  xlim(0, 40)+
  theme_minimal() +
  labs(title = "Period by Distance", x = "Distance from Particle (px)", y = "Period (h)")

#' plot showing amplitude of all cells according to distance from particle

amplitude_dotplot_trace_allcells <- ggplot()+#merged_table, aes(x = distance, y = period)) +
  # geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
  geom_smooth(data = merged_table[merged_table$treatment == "FRED", ], aes(x = distance, y = amplitude, colour = distance_group, group = distance_group, linetype = treatment), method = "loess", se =  TRUE, level = 0.95) +
  geom_smooth(data = merged_table[merged_table$treatment == "MUT", ], aes(x = distance, y = amplitude, colour = distance_group, group = distance_group, linetype = treatment), method = "loess", se =  TRUE, level = 0.95) +
  scale_color_manual(values = color_palette)+
  xlim(0, 40)+
  theme_minimal() +
  labs(title = "Amplitude by Distance", x = "Distance from Particle (px)", y = "Amplitude (A.U.)")

RAE_dotplot_trace_allcells <- ggplot()+#merged_table, aes(x = distance, y = period)) +
  # geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
  geom_smooth(data = merged_table[merged_table$treatment == "FRED", ], aes(x = distance, y = RAE, colour = distance_group, group = distance_group, linetype = treatment), method = "loess", se =  TRUE, level = 0.95) +
  geom_smooth(data = merged_table[merged_table$treatment == "MUT", ], aes(x = distance, y = RAE, colour = distance_group, group = distance_group, linetype = treatment), method = "loess", se =  TRUE, level = 0.95) +
  scale_color_manual(values = color_palette)+
  xlim(0, 40)+
  theme_minimal() +
  labs(title = "RAE by Distance", x = "Distance from Particle (px)", y = "RAE (A.U.)")

AUC_dotplot_trace_allcells <- ggplot()+#merged_table, aes(x = distance, y = period)) +
  # geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
  geom_smooth(data = merged_table[merged_table$treatment == "FRED", ], aes(x = distance, y = AUC, colour = distance_group, group = distance_group, linetype = treatment), method = "loess", se =  TRUE, level = 0.95) +
  geom_smooth(data = merged_table[merged_table$treatment == "MUT", ], aes(x = distance, y = AUC, colour = distance_group, group = distance_group, linetype = treatment), method = "loess", se =  TRUE, level = 0.95) +
  scale_color_manual(values = color_palette)+
  xlim(0, 40)+
  theme_minimal() +
  labs(title = "AUC by Distance", x = "Distance from Particle (px)", y = "AUC (A.U.)")

savePlots(obj_to_save = list(amplitude_vs_dist = amplitude_dotplot_trace_allcells,
                             period_vs_dist = period_dotplot_trace_allcells,
                             RAE_vs_dist = RAE_dotplot_trace_allcells,
                             AUC_vs_dist = AUC_dotplot_trace_allcells), basepath = lineplots_fold, extension = "svg", p.width = 1850, p.height = 1000)
}
 #### aligned traces all together ---------------------------------------------------------
if(ALIGN_TRACES_COMPARE){
  traces_fold = file.path(wd, "traces_comparison")
  
  if(!dir.exists(traces_fold)){
    dir.create(traces_fold)
  }
traces_list <- dataset_list$align.trx.raw
traces_annotated = list()
# merge all tables together adding info about sample name and treatment
 for(i in 1:length(traces_list)){
  #import dataset
  traces <- traces_list[[i]]
  name = names(traces_list)[i]
  # isolate treatment
  treatments = c("FRED", "MUT")
  match <- sapply(treatments, function(y) grepl(y, name))
  treatment = treatments[match]
  # attention
  # this only works if name has side_well_treatmnt_number
  # TODO swap midpoint for center of mass
  midpoint <- median(as.numeric(traces$ID))
  traces$position <- NA
  top <- which(as.numeric(traces$ID) < midpoint)
  bottom <- which(as.numeric(traces$ID) >= midpoint)
  traces$position[top] <- "top"
  traces$position[bottom] <- "bottom"
  # assign values to columns
  traces$sample = name
  traces$treatment = treatment 
  traces$ID = paste(name, traces$ID, sep = "_")
  traces_annotated[[name]] = traces
   
}

# create average of each sample and save in list
traces_summary_list = lapply(traces_annotated, function(traces){
  traces_summary <- traces %>% 
    dplyr::group_by(treatment, group, t) %>% 
    dplyr::summarise(mean = mean(intensity, na.rm = TRUE),
                     se = sd(intensity, na.rm = TRUE) / sqrt(n()),
                     size = n(),
                     .groups = "drop")
  traces_summary$sample = unique(traces$sample)
  return(traces_summary)
})

# save traces objects
saveRDS(traces_summary_list, file = file.path(rds_fold, paste0(week, "_traces_summary_list.rds")))
saveRDS(traces_annotated, file = file.path(rds_fold, paste0(week, "_traces_annotated.rds")))


# bind tables
#' averaged_traces_bound contains all the averages of traces calculated in the
#' traces_summary_list in one table. Each average is the mean trace from one
#' sample and one area 
averaged_traces_bound = do.call(rbind, traces_summary_list) %>% `rownames<-`(NULL)

avg_traces_plot4 = ggplot2::ggplot(dplyr::filter(averaged_traces_bound))+# %>% dplyr::filter(., sample == "LEFT_B3_MUT"))+
  # ggplot2::geom_ribbon(aes(x = t, ymin = average-se, ymax = average+se, alpha = 0.5, colour = group, fill = group))+
  geom_line(aes(x = t, y = mean, colour = group, linetype = sample))+
  ggplot2::scale_alpha_identity()+
  # ggplot2::scale_linetype_identity()+
  scale_x_continuous(name = "Time (h)", breaks = seq(0, 184, by = 24), minor_breaks = seq(0, 184, by = 24))+#as.integer(max(t)), by = 24), minor_breaks = seq(12, as.integer(max(t)), by = 24))+
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

#' all_cells_bound has all the traces from the traces_annotated in one table
all_cells_bound = do.call(rbind, traces_annotated) %>% `rownames<-`(NULL)

# perform averages
traces_summary_all_samples_weighted = averaged_traces_bound %>% 
  dplyr::group_by(treatment, group, t) %>% 
  dplyr::summarise(average = sum((mean*size), na.rm = TRUE)/sum(size, na.rm = TRUE),
                   se = (sd((mean*size), na.rm = TRUE)/sum(size, na.rm = TRUE))/sqrt(length(unique(sample))),
                   n = length(unique(sample)),
                   cells = sum(size, na.rm = TRUE)
                   )

traces_summary_all_samples_unweighted <- averaged_traces_bound %>% 
  dplyr::group_by(treatment, group, t) %>% 
  dplyr::summarise(average = mean(mean, na.rm = TRUE),
                   se = sd(mean, na.rm = TRUE)/sqrt(n()),
                   n = n(),
                   size = sum(size),
                   .groups = "drop")

all_cells_averaged <- all_cells_bound %>% 
  dplyr::group_by(treatment, group, t) %>% 
  dplyr::summarise(average = mean(intensity, na.rm = TRUE),
                   se = sd(intensity, na.rm = TRUE) / sqrt(n()),
                   size = n(),
                   .groups = "drop")

#' compare the same group in different treatments
close_cells_avg = dplyr::filter(all_cells_averaged, group == "close")
mid_cells_avg = dplyr::filter(all_cells_averaged, group == "mid")
far_cells_avg = dplyr::filter(all_cells_averaged, group == "far")
# veryfar_cells_avg = dplyr::filter(all_cells_averaged, group == "veryfar")

#' save tables
close_cells_path = file.path(traces_fold, "close_avgs.csv")
write.csv(close_cells_avg, file = close_cells_path)

mid_cells_path = file.path(traces_fold, "mid_avgs.csv")
write.csv(mid_cells_avg, file = mid_cells_path)

far_cells_path = file.path(traces_fold, "far_avgs.csv")
write.csv(far_cells_avg, file = far_cells_path)

#veryfar_cells_path = file.path(plot_folder, "veryfar_avgs.csv")
#write.csv(veryfar_cells_avg, file = veryfar_cells_path)

#' plot comparison



# CF_Mut_cells = dplyr::filter(all_cells_averaged, treatment == "MUT")
# MUT_cells_plot = ggplot2::ggplot(CF_Mut_cells)+
#   ggplot2::geom_ribbon(aes(x = t, ymin = average-se, ymax = average+se, alpha = 0.5, colour = group, fill = group))+
#   geom_line(aes(x = t, y = average, colour = group))+
#   ggplot2::scale_alpha_identity()+
#   scale_x_continuous(name = "Time (h)", breaks = seq(0, 184, by = 24), minor_breaks = seq(0, 184, by = 24))+#as.integer(max(t)), by = 24), minor_breaks = seq(12, as.integer(max(t)), by = 24))+
#   theme_minimal()+
#   theme(axis.line.y = element_line(linewidth = 1),
#         axis.line.x = element_line(linewidth = 1),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_line(linewidth = 0.5, linetype = "88", colour = "black"),
#         plot.subtitle = element_text(size = 10, hjust = 0),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         axis.title.x  = element_text(size = 20),
#         axis.title.y  = element_text(size = 20))

close_cells_plot = ggplot2::ggplot(close_cells_avg)+
  ggplot2::geom_ribbon(aes(x = t, ymin = average-se, ymax = average+se, alpha = 0.5, colour = treatment, fill = treatment))+
  geom_line(aes(x = t, y = average, colour = treatment))+
  ggplot2::scale_alpha_identity()+
  scale_x_continuous(name = "Time (h)", breaks = seq(0, 184, by = 24), minor_breaks = seq(12, 184, by = 24))+#as.integer(max(t)), by = 24), minor_breaks = seq(12, as.integer(max(t)), by = 24))+
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

mid_cells_plot = ggplot2::ggplot(mid_cells_avg)+
  ggplot2::geom_ribbon(aes(x = t, ymin = average-se, ymax = average+se, alpha = 0.5, colour = treatment, fill = treatment))+
  geom_line(aes(x = t, y = average, colour = treatment))+
  ggplot2::scale_alpha_identity()+
  scale_x_continuous(name = "Time (h)", breaks = seq(0, 184, by = 24), minor_breaks = seq(12, 184, by = 24))+#as.integer(max(t)), by = 24), minor_breaks = seq(12, as.integer(max(t)), by = 24))+
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

far_cells_plot = ggplot2::ggplot(far_cells_avg)+
  ggplot2::geom_ribbon(aes(x = t, ymin = average-se, ymax = average+se, alpha = 0.5, colour = treatment, fill = treatment))+
  geom_line(aes(x = t, y = average, colour = treatment))+
  ggplot2::scale_alpha_identity()+
  scale_x_continuous(name = "Time (h)", breaks = seq(0, 184, by = 24), minor_breaks = seq(12, 184, by = 24))+#as.integer(max(t)), by = 24), minor_breaks = seq(12, as.integer(max(t)), by = 24))+
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

savePlots(obj_to_save = list(close_cells_plot = close_cells_plot, mid_cells_plot = mid_cells_plot, far_cells_plot = far_cells_plot),
          filename = "", basepath = traces_fold, extension = "svg", p.width = 1800, p.height = 1000)

# savePlots(obj_to_save = list(MUT_compare = MUT_cells_plot),
# filename = "", basepath = newdir, extension = "svg", p.width = 1800, p.height = 1000)

}
 #### TOP-BOTTOM comparison ---------------------------------------------------------
if(TOPBOTTOM && ALIGN_TRACES_COMPARE){
  # compare top traces with bottom traces
  plot_trace_comparison = function(data, sample, pos1, pos2, group1, group2, col1 = position, col2 = position, p.width = 1800, p.height = 1000, basepath, extension = "svg"){
    plot = ggplot2::ggplot()+
      ggplot2::geom_ribbon(data = dplyr::filter(data[[sample]], position == pos1, group == group1), aes(x = t, ymin = mean-se, ymax = mean+se, alpha = 0.5, colour = {{col1}}, fill = {{col1}}))+
      ggplot2::geom_ribbon(data = dplyr::filter(data[[sample]], position == pos2, group == group2), aes(x = t, ymin = mean-se, ymax = mean+se, alpha = 0.5, colour = {{col2}}, fill = {{col2}}))+
      ggplot2::scale_alpha_identity()+
      scale_x_continuous(name = "Time (h)", breaks = seq(0, 184, by = 24), minor_breaks = seq(0, 184, by = 24))+#as.integer(max(t)), by = 24), minor_breaks = seq(12, as.integer(max(t)), by = 24))+
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
    
    # Create the file path
    plot_path = file.path(basepath, paste0("trace_", sample, "_", pos1, "_", pos2, "_", group1, "_", group2, ".", extension))
    
    # Save the plot
    ggplot2::ggsave(plot_path, plot = plot, width = p.width, height = p.height, units = "px")
    
    return(plot)
  }
  mutdir = file.path(traces_fold, "MUT")
  if(!dir.exists(mutdir)){
    dir.create(mutdir)
  }
  # Treated samples
  plot_trace_comparison(traces_summary_list, "LEFT_B3_MUT", "top", "bottom", "close", "close", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B3_MUT", "top", "bottom", "mid", "mid", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B3_MUT", "top", "bottom", "far", "far", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B3_MUT", "top", "top", "close", "far", group, group, basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B3_MUT", "bottom", "bottom", "close", "far", group, group, basepath = mutdir)
  
  plot_trace_comparison(traces_summary_list, "RIGHT_B6_MUT", "top", "bottom", "close", "close", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B6_MUT", "top", "bottom", "mid", "mid", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B6_MUT", "top", "bottom", "far", "far", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B6_MUT", "top", "top", "close", "far", group, group, basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B6_MUT", "bottom", "bottom", "close", "far", group, group, basepath = mutdir)
  
  plot_trace_comparison(traces_summary_list, "LEFT_C2_MUT", "top", "bottom", "close", "close", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_C2_MUT", "top", "bottom", "mid", "mid", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_C2_MUT", "top", "bottom", "far", "far", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_C2_MUT", "top", "top", "close", "far", group, group, basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_C2_MUT", "bottom", "bottom", "close", "far", group, group, basepath = mutdir)
  
  plot_trace_comparison(traces_summary_list, "LEFT_D2_MUT", "top", "bottom", "close", "close", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_D2_MUT", "top", "bottom", "mid", "mid", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_D2_MUT", "top", "bottom", "far", "far", basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_D2_MUT", "top", "top", "close", "far", group, group, basepath = mutdir)
  plot_trace_comparison(traces_summary_list, "LEFT_D2_MUT", "bottom", "bottom", "close", "far", group, group, basepath = mutdir)
  
  # Control Samples
  
  wtdir = file.path(traces_fold, "WT")
  if(!dir.exists(wtdir)){
    dir.create(wtdir)
  }
  plot_trace_comparison(traces_summary_list, "LEFT_B1_FRED", "top", "bottom", "close", "close", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B1_FRED", "top", "bottom", "mid", "mid", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B1_FRED", "top", "bottom", "far", "far", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B1_FRED", "top", "top", "close", "far", group, group, basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B1_FRED", "bottom", "bottom", "close", "far", group, group, basepath = wtdir)
  
  plot_trace_comparison(traces_summary_list, "LEFT_B4_FRED", "top", "bottom", "close", "close", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B4_FRED", "top", "bottom", "mid", "mid", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B4_FRED", "top", "bottom", "far", "far", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B4_FRED", "top", "top", "close", "far", group, group, basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_B4_FRED", "bottom", "bottom", "close", "far", group, group, basepath = wtdir)
  
  plot_trace_comparison(traces_summary_list, "LEFT_A1_FRED", "top", "bottom", "close", "close", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_A1_FRED", "top", "bottom", "mid", "mid", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_A1_FRED", "top", "bottom", "far", "far", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_A1_FRED", "top", "top", "close", "far", group, group, basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_A1_FRED", "bottom", "bottom", "close", "far", group, group, basepath = wtdir)
  
  plot_trace_comparison(traces_summary_list, "LEFT_C3_FRED", "top", "bottom", "close", "close", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_C3_FRED", "top", "bottom", "mid", "mid", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_C3_FRED", "top", "bottom", "far", "far", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_C3_FRED", "top", "top", "close", "far", group, group, basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_C3_FRED", "bottom", "bottom", "close", "far", group, group, basepath = wtdir)
  
  plot_trace_comparison(traces_summary_list, "RIGHT_B1_FRED", "top", "bottom", "close", "close", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B1_FRED", "top", "bottom", "mid", "mid", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B1_FRED", "top", "bottom", "far", "far", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B1_FRED", "top", "top", "close", "far", group, group, basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B1_FRED", "bottom", "bottom", "close", "far", group, group, basepath = wtdir)
  
  plot_trace_comparison(traces_summary_list, "LEFT_D3_FRED", "top", "bottom", "close", "close", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_D3_FRED", "top", "bottom", "mid", "mid", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_D3_FRED", "top", "bottom", "far", "far", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_D3_FRED", "top", "top", "close", "far", group, group, basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "LEFT_D3_FRED", "bottom", "bottom", "close", "far", group, group, basepath = wtdir)
  
  plot_trace_comparison(traces_summary_list, "RIGHT_C1_FRED", "top", "bottom", "close", "close", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_C1_FRED", "top", "bottom", "mid", "mid", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_C1_FRED", "top", "bottom", "far", "far", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_C1_FRED", "top", "top", "close", "mid", group, group, basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_C1_FRED", "bottom", "bottom", "close", "mid", group, group, basepath = wtdir)
  
  plot_trace_comparison(traces_summary_list, "RIGHT_B4_FRED", "top", "bottom", "close", "close", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B4_FRED", "top", "bottom", "mid", "mid", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B4_FRED", "top", "bottom", "far", "far", basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B4_FRED", "top", "top", "close", "far", group, group, basepath = wtdir)
  plot_trace_comparison(traces_summary_list, "RIGHT_B4_FRED", "bottom", "bottom", "close", "far", group, group, basepath = wtdir)
}
 #### analyse phases distances between groups ---------------------------------------------------------
if(PHASE_DIST_GROUPS){
  phase_fold = file.path(wd, "phase_dist")
  
  if(!dir.exists(phase_fold)){
    dir.create(phase_fold)
  }
  
  # read circular stats
circ_summary = read.csv(file = file.path(wd, "circular_stats.csv"), row)

  close_mid = abs(circ_summary$close_meanPh - circ_summary$mid_meanPh) %>% sapply(., function(x){if(x>pi){x = 2*pi - x}else{x = x}})
  close_far = abs(circ_summary$close_meanPh - circ_summary$far_meanPh) %>% sapply(., function(x){if(x>pi){x = 2*pi - x}else{x = x}})
  mid_far = abs(circ_summary$mid_meanPh - circ_summary$far_meanPh) %>% sapply(., function(x){if(x>pi){x = 2*pi - x}else{x = x}})
  phase_dist_table = data.frame(filename = circ_summary[,1],
                                close_far = close_far,
                                close_mid = close_mid,
                                mid_far = mid_far,
                                treatment = sapply(circ_summary[,1], function(x){strsplit(x, split = "_")[[1]][3]}))
  
  
  # create plot to show differences
  phase_dist_longtable = tidyr::pivot_longer(phase_dist_table, cols = c("close_far", "close_mid", "mid_far"), values_to = "measure", names_to = "metric")
  phase_dist_longtable$measure = ((phase_dist_longtable$measure/pi)*12)
  phase_dist_longtable$metric = as.factor(phase_dist_longtable$metric)
  phase_dist_longtable$treatment = as.factor(phase_dist_longtable$treatment)
  
  summary_phase_stats = phase_dist_longtable %>%
    group_by(treatment, metric) %>%
    summarise(phase_mean = mean(measure),
              .groups = "drop")
  
  res.aov.phase <- rstatix::anova_test(data = phase_dist_longtable, formula = measure ~ metric*treatment)
  pwc_phase <- phase_dist_longtable %>%
    group_by(metric) %>%
    rstatix::pairwise_t_test(
      formula = measure ~ treatment,
      paired = FALSE,
      p.adjust.method = "bonferroni"
    )
  
  phase_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = phase_dist_longtable, ggplot2::aes(x = metric, y = measure, color = treatment), stat = "identity", position = position_dodge(0.3), alpha = 0.5) + # scatterplot
    ggplot2::geom_point(data = summary_phase_stats, ggplot2::aes(x = metric, y = phase_mean, color = treatment), size = 3, stat = "identity", position = position_dodge(0.3))+
    labs(
      title = "Phase distance between areas",
      x = "",
      y = "Phase distance (h)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  # create table with ratio between distance close-mid and mid-far
  phase_ratio_table = data.frame(filename = phase_dist_table$filename,
                                 ratio = phase_dist_table$close_mid/phase_dist_table$mid_far,
                                 treatment = phase_dist_table$treatment)
  
  phase_ratio_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = phase_ratio_table, ggplot2::aes(x = treatment, y = ratio, color = treatment), stat = "identity", position = position_dodge(0.3), alpha = 0.5) + # scatterplot
    labs(
      title = "Phase ratio between areas",
      x = "",
      y = "Phase ratio"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  savePlots(obj_to_save = list(phase_plot, phase_ratio_plot),  savefold = phase_fold, extension = "png", p.width = 1000)

}
 #### plot vector length spread ---------------------------------------------------------
if(VECTOR_SPREAD){
  vec_folder = file.path(wd, "vector_length")
  if(!dir.exists(vec_folder)){
    dir.create(vec_folder)
  }
    mu_table = data.frame(filename = circ_summary[,1],
                          close = circ_summary$close_vec_len,
                          mid = circ_summary$mid_vec_len,
                          far = circ_summary$far_vec_len,
                          treatment = sapply(circ_summary[,1], function(x){strsplit(x, split = "_")[[1]][3]}))
    
    mu_longtable = tidyr::pivot_longer(mu_table, cols = c("close", "mid", "far"), values_to = "measure", names_to = "metric")
    
    summary_mu_stats = mu_longtable %>%
      group_by(treatment, metric) %>%
      summarise(mu_mean = mean(measure),
                .groups = "drop")
    
    mu_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = mu_longtable, ggplot2::aes(x = metric, y = measure, color = treatment), stat = "identity", position = position_dodge(0.3), alpha = 0.5) + # scatterplot
      ggplot2::geom_point(data = summary_mu_stats, ggplot2::aes(x = metric, y = mu_mean, color = treatment), size = 3, stat = "identity", position = position_dodge(0.3))+
      labs(
        title = "Vector length in areas",
        x = "",
        y = "Vector length (A.U.)"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    
    savePlots(obj_to_save = mu_plot,  savefold = vec_folder, extension = "png", p.width = 1000)
  }
 #### 3D Plots  ---------------------------------------------------------
if(THREE_D_PLOTS){
  threeD_folder = file.path(wd, "3D_plots")
  
  if(!dir.exists(threeD_folder)){
    dir.create(threeD_folder)
  }
  
  # summarize merged table to have one value per sample per distance group
  merged_table_summary <- merged_table %>% 
    group_by(treatment, well, distance_group) %>% 
    dplyr::select(period, amplitude, RAE, AUC, phase_norm, distance, median, SD_green) %>%
    filter(!is.na(period), !is.na(amplitude), !is.na(RAE), !is.na(AUC), !is.na(distance), !is.na(median), !is.na(SD_green)) %>% 
    dplyr::summarise(period = mean(period),
                     amplitude = mean(amplitude),
                     RAE = mean(RAE),
                     AUC_mean = mean(AUC),
                     AUC_SD = sd(AUC),
                     phase_mean = circular::mean.circular(circular::circular(phase_norm, type = "angles", units = "hours")),
                     phase_SD = circular::sd.circular(circular::circular(phase_norm, type = "angles", units = "hours")),
                     distance = mean(distance),
                     median = mean(median),
                     SD_green = mean(SD_green), .groups = "drop")
  
  # make a 3D plot that represents the median level (x) with period (y) and distance from particle (z) 
  dist_AUC_per_plot_summary = phase_AUC_per_plot_summary <- plot_ly() %>%
                    # Add scatter plot points
                    add_trace(
                      data = merged_table_summary %>% filter(treatment == "FRED"),
                      x = ~distance,
                      y = ~AUC_mean,
                      z = ~period,
                      type = "scatter3d",
                      mode = "markers",
                      split = ~distance_group,  # Ensures proper grouping by discrete variable
                      marker = list(symbol = "square", size = ~AUC_SD/400, opacity = 1),
                      color = ~distance_group,
                      colors = color_palette,
                      # colors = ,  # Assign custom colors to factor levels
                      showlegend = TRUE) %>%
                    # Add scatter plot points
                    add_trace(
                      data = merged_table_summary %>% filter(treatment == "MUT"),
                      x = ~distance,
                      y = ~AUC_mean,
                      z = ~period,
                      type = "scatter3d",
                      mode = "markers",
                      split = ~distance_group,  # Ensures proper grouping by discrete variable
                      marker = list(symbol = "circle", size = ~AUC_SD/400, opacity = 1),
                      color = ~distance_group,
                      colors = color_palette,
                      showlegend = TRUE) %>%
                    # Layout settings
                    layout(
                      scene = list(
                        xaxis = list(title = "Distance (px)", titlefont = list(size = 14), range = c(30, 0)),
                        yaxis = list(title = "Mean AUC (A.U.)", titlefont = list(size = 14), range = c(0,50000)),
                        zaxis = list(title = "Period (h)", titlefont = list(size = 14), range = c(22.4,24.5)),
                        camera = list(eye = list(x = 2.7, y = -1.2, z = 1))
                      ),
                      title = list(
                        text = "3D plot of phase vs AUC vs period",
                        font = list(size = 16)
                      )
                    )
  
  phase_AUC_per_plot_summary <- plot_ly() %>%
                      # Add scatter plot points
                      add_trace(
                        data = merged_table_summary %>% filter(treatment == "FRED"),
                        x = ~phase_mean,
                        y = ~AUC_mean,
                        z = ~period,
                        type = "scatter3d",
                        mode = "markers",
                        split = ~distance_group,  # Ensures proper grouping by discrete variable
                        marker = list(symbol = "square", size = ~phase_SD*10, opacity = 1),
                        color = ~distance_group,
                        colors = color_palette,
                        # colors = ,  # Assign custom colors to factor levels
                        showlegend = TRUE) %>%
                      # Add scatter plot points
                      add_trace(
                        data = merged_table_summary %>% filter(treatment == "MUT"),
                        x = ~phase_mean,
                        y = ~AUC_mean,
                        z = ~period,
                        type = "scatter3d",
                        mode = "markers",
                        split = ~distance_group,  # Ensures proper grouping by discrete variable
                        marker = list(symbol = "circle", size = ~phase_SD*10, opacity = 1),
                        color = ~distance_group,
                        colors = color_palette,
                        showlegend = TRUE) %>%
                      # Layout settings
                      layout(
                        scene = list(
                          xaxis = list(title = "Phase (h)", titlefont = list(size = 14), range = c(12, -12)),
                          yaxis = list(title = "Mean AUC (A.U.)", titlefont = list(size = 14), range = c(0,50000)),
                          zaxis = list(title = "Period (h)", titlefont = list(size = 14), range = c(22.4,24.5)),
                          camera = list(eye = list(x = -2.7, y = -1.2, z = 1))
                        ),
                        title = list(
                          text = "3D plot of phase vs AUC vs period",
                          font = list(size = 16)
                        )
                      )
  
  dist_AUC_per_plot = threeDplot_AUC(merged_table, distance, AUC, period, distance_group, 
                                     palette = color_palette, x_lab = "Distance (px)", y_lab = "AUC (a.u.)", 
                                     z_lab = "Period (h)", eye = "pos2", zoom = 0.5)
  
  phase_AUC_per_plot_FRED = threeDplot_phase(merged_table %>% filter(treatment == "FRED"), phase_circ, AUC, period, distance_group, 
                                  color_palette, x_lab = "Phase (h)", y_lab = "AUC (a.u.)", 
                                  z_lab = "Period (h)", eye = "pos2", zoom = 0.5)
  
  phase_AUC_per_plot_MUT = threeDplot_phase(merged_table %>% filter(treatment == "MUT"), phase_circ, AUC, period, distance_group, 
                                             color_palette, x_lab = "Phase (h)", y_lab = "AUC (a.u.)", 
                                             z_lab = "Period (h)", eye = "pos2", zoom = 0.5)
  
  # save 3D plot object for later visualization
  if(TRUE){
    dist_AUC_per_summary_path = file.path(threeD_folder, "3D_dist_auc_per_summary.rds")
    saveRDS(dist_AUC_per_plot_summary, file = dist_AUC_per_summary_path)
    
    phase_AUC_per_summary_path = file.path(threeD_folder, "3D_phase_auc_per_summary.rds")
    saveRDS(phase_AUC_per_plot_summary, file = phase_AUC_per_summary_path)
    
    dist_AUC_per_plot_path = file.path(threeD_folder, "All_3D_dist_auc_per_plot.rds")
    saveRDS(dist_AUC_per_plot, file = dist_AUC_per_plot_path)
    
    phase_AUC_per_plot_FRed_path = file.path(threeD_folder, "FRed_3D_phase_auc_per_plot.rds")
    saveRDS(phase_AUC_per_plot_FRED, file = phase_AUC_per_plot_FRed_path)
    
    phase_AUC_per_plot_MUT_path = file.path(threeD_folder, "MUT_3D_phase_auc_per_plot.rds")
    saveRDS(phase_AUC_per_plot_MUT, file = phase_AUC_per_plot_MUT_path)
  }
  
  if(TRUE){
    # save html files of 3D plots
    dist_AUC_per_summary_path = file.path(threeD_folder, "3D_dist_auc_per_summary.html")
    saveWidget(dist_AUC_per_plot_summary, file = dist_AUC_per_summary_path, selfcontained = TRUE)
    
    phase_AUC_per_summary_path = file.path(threeD_folder, "3D_phase_auc_per_summary.html")
    saveWidget(phase_AUC_per_plot_summary, file = phase_AUC_per_summary_path, selfcontained = TRUE)
    
    dist_AUC_per_plot_path = file.path(threeD_folder, "All_3D_dist_auc_per_plot.html")
    saveWidget(dist_AUC_per_plot, file = dist_AUC_per_plot_path, selfcontained = TRUE)
    
    phase_AUC_per_plot_FRed_path = file.path(threeD_folder, "FRed_3D_phase_auc_per_plot.html")
    saveWidget(phase_AUC_per_plot_FRED, file = phase_AUC_per_plot_FRed_path, selfcontained = TRUE)
    
    phase_AUC_per_plot_MUT_path = file.path(threeD_folder, "MUT_3D_phase_auc_per_plot.html")
    saveWidget(phase_AUC_per_plot_MUT, file = phase_AUC_per_plot_MUT_path, selfcontained = TRUE)
    
  }
}
#### not working ---------------------------------------------------------
IF(FALSE){
 # mixed model analysis of period ####
#' try mixed model effect for analysis of periods
library(lme4)
library(lmerTest)
model <- lmer(period ~ distance + distance_group + treatment +
                distance:treatment + distance_group:treatment + (1 | sample),
              data = merged_table)
model_null <- lmer(period ~ distance_group + treatment +
                     (1 | sample),
                   data = merged_table)

anovaTest <- anova(model_null, model)

library(sjPlot)
sjPlot::plot_model(model, type = "pred", terms = c("distance_group", "treatment"))

#' output parameter ready to be input in Prism for checking statistical analysis
#' TODO divide table based on parameters and order by groups
write.csv(merged_table, file = file.path(plot_folder, "merged_table.csv"))

#' compare groups across conditions with many samples grouped and add nested anova

}
# save 3d plots
#call plot
i = 7
print(filenames[i])
threeD_plot_list$distance[[i]] %>% 
  layout(scene = list(camera = list(eye = list(x = 2, y = 2, z = 1.5))))

threeD_plot_list$phase[[i]] %>% 
  layout(scene = list(camera = list(eye = list(x = 2, y = 2, z = 1.5))))

computePeriod = function(df, excludeNC = FALSE, top = 30, bottom = 18, save.trace = FALSE, rm.start = 0, ...){
  
  #prepare table
  data_df <- prep_table(df, rm.start)
  
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
  period_table = do.call(rbind, lapply(results_list, as.data.frame))
  traces_table = do.call(cbind, lapply(traces_list, as.data.frame)) %>% t() %>% `rownames<-`(names(traces_list))
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


# fix path
back_to_forw <- function(x){
  stringr::str_replace_all(x, "\\\\", "//")
}
