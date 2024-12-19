#' script to collate the data from the different samples and run comparisons

#' objects to import and collate
#'  grid_vals or distance matrix
#'  period table

# install and load required packages
pkg_req <- c("sf", "RImageJROI", "ggplot2", "dplyr", "pracma", "minpack.lm", "magrittr", "stringr", "circular", "svglite", "astsa", "pdftools", "scico", "ggpubr", "rstatix")

# Install and load packages if not already installed
for (package in pkg_req) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}


 # FUNCTIONS #####
# function to generate plot and perform statistical testing
nested_anova_plot <- function(data, metric, distance_var, well = well, treatment_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, nudge = 0.5, step = 0.2, ...) {
  
  # Compute normality test if needed
  if (norm_test) {
    normality <- data %>%
      group_by({{distance_var}}, {{treatment_var}}, {{well}}) %>%
      rstatix::shapiro_test({{metric}})
    print(normality)
  }
  
  # Nested ANOVA: Include the 'well' variable in the nested model
  formula_nested <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(distance_var)), "/", deparse(substitute(well))))
  
  # Perform the nested ANOVA
  res.aov <- rstatix::anova_test(data = data, formula = formula_nested)
  print(res.aov)
  
  # Calculate mean and SE for each well
  well_stats <- data %>%
    group_by({{distance_var}}, {{treatment_var}}, {{well}}) %>%
    summarise(
      well_mean = mean({{metric}}, na.rm = TRUE),
      well_se = sd({{metric}}, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Calculate overall mean and SE across wells
  overall_stats <- well_stats %>%
    group_by({{distance_var}}, {{treatment_var}}) %>%
    summarise(
      abs_mean = mean(well_mean, na.rm = TRUE),
      abs_se = sd(well_mean, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Perform pairwise comparisons if needed
  if (comparison) {
    formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(treatment_var))))
    pwc <- data %>%
      group_by({{distance_var}}) %>%
      rstatix::pairwise_t_test(
        formula = formula,
        paired = FALSE,
        p.adjust.method = test
      )
    
    # Add XY position for p-values
    pwc <- pwc %>%
      rstatix::add_xy_position(
        x = as_label(enquo(distance_var)),
        dodge = 0.8,
        fun = "mean_se", # Use 'max' to calculate the position relative to all values
        step.increase = 0.2
      )
    
    formula2 <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(distance_var))))
    pwc2 <- data %>%
      group_by({{treatment_var}}) %>%
      rstatix::pairwise_t_test(
        formula = formula2,
        paired = FALSE,
        p.adjust.method = test
      )
    
    # Add XY position for p-values
    pwc2 <- pwc2 %>%
      rstatix::add_xy_position(
        x = as_label(enquo(distance_var)), # Specify the distance variable for positioning
        dodge = 0.8, # Dodge ensures separation by treatment groups
        group = as_label(enquo(treatment_var)), # Include treatment grouping
        fun = "mean_se", # Position based on mean and SE
        step.increase = step # Adjust to avoid overlap with the first comparison
      )
  }
  # Generate the plot
  plot <- ggplot() +
    # Add well-level means
    geom_point(
      data = well_stats,
      aes(x = {{distance_var}}, y = well_mean, color = {{treatment_var}}),
      position = position_dodge(0.8),
      size = 2, # Smaller points for individual wells
      alpha = 0.5 # Transparency to distinguish them
    ) +
    # Add overall means with error bars
    geom_point(
      data = overall_stats,
      aes(x = {{distance_var}}, y = abs_mean, color = {{treatment_var}}),
      position = position_dodge(0.8),
      size = 4 # Larger points for overall means
    ) +
    geom_errorbar(
      data = overall_stats,
      aes(
        x = {{distance_var}},
        ymin = abs_mean - abs_se,
        ymax = abs_mean + abs_se,
        color = {{treatment_var}}
      ),
      width = 0.2,
      position = position_dodge(0.8)
    ) +
    # Add pairwise comparisons (optional)
    ggpubr::stat_pvalue_manual(
      pwc,
      dodge = 0.8
    )+
    ggpubr::stat_pvalue_manual(
      pwc2,
      dodge = 0.8,
      bracket.nudge.y = nudge
    )+
    # Customize colors and labels
    scale_color_manual(values = plot_colors) +
    labs(
      title = paste("Nested ANOVA of", deparse(substitute(metric)), "by Distance Group"),
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = { if (comparison) rstatix::get_pwc_label(pwc) else NULL },
      y = ylabel,
      x = xlabel
    ) +
    theme_minimal()
  
  return(plot)
}

# function to save plots 
savePlots = function(obj_to_save, savefold, extension, p.width = 1560, p.height = 1100) {
  for (i in seq_along(obj_to_save)) {
    # Check if the list has names; if not, use a default numeric index
    plot_name <- if (!is.null(names(obj_to_save))) {
      names(obj_to_save)[i]
    } else {
      paste0("plot_", i)
    }
    
    # Create the file path
    plot_path = file.path(savefold, paste0(plot_name, ".", extension))
    
    # Save the plot
    ggplot2::ggsave(plot_path, plot = obj_to_save[[i]], width = p.width, height = p.height, units = "px")
  }
}

# function to align traces to the mean phase of the group
phase_align_trace = function(traces_table, period_table, remove_start = 0, remove_end = 0, align_to = NA, debug = FALSE){
  browser()
  #' TODO check align phase value
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
  browser()
  return(traces_aligned)
}

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

 # EXECUTION ######

 # IMPORTING #####
#' select home dir where to pull names and data from
wd = r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Proj_Tau\Tau_an_dev\wk5_hemislices\Left)"

files = list.files(path = wd, pattern = ".tif*$")
filenames = stringr::str_remove(files, pattern = ".tif*")
foldernames = file.path(wd, paste(filenames, "_results", sep = ""))

dataset_list <- list(grid.vals = list(),
                     dist.mtx = list(),
                     period.tbl = list(),
                     align.trx = list())

missing_matrices <- vector()
missing_periods <- vector()
missing_merg_tbls <- vector()
missing_traces <- vector()

# period tbls, distance matrices, merged tables, and aligned traces
for (i in seq_len(length(files))){
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
    aligned_trx <- readRDS(align_trx_wide_path)
    dataset_list$align.trx.wide[[filename]] <- aligned_trx
  }else{
    print(paste(filename, " aligned traces wide file missing."))
    missing_traces <- c(missing_traces, filename)
  }
  
}

print("Datasets added!")

 # DISTANCE GROUS ####
#' compare groups across conditions for only two sample being compared

limit1 = 3
limit2 = 12

#' interrogate the matrix to see which cells are at less than
close_cells_matrix = lapply(dataset_list$dist.mtx,  function(x){x<=limit1})
far_cells_matrix = lapply(dataset_list$dist.mtx,  function(x){x>limit2})
between_cells_matrix = lapply(dataset_list$dist.mtx,  function(x){x>limit1 & x <= limit2})

color_palette = c("#B23A48", "#2F9C95", "#B6DC76")

 # analysis to visualize circadian parameters in relation to the distance from plaques ####

#' create list of distances and extract median value
distance_list <- lapply(dataset_list$dist.mtx, function(matrix){
  distance_table = cbind(as.numeric(rownames(matrix)), round(matrixStats::rowMedians(matrix), 3)) %>% `colnames<-`(c("ID", "distance"))
})

#' merge distance table with period table, ten assign name to list element
merged_list <- setNames(
  lapply(names(distance_list), function(x, dataset_list, distance_list){
    period_tbl <- dataset_list$period.tbl[[x]]
    distance_tbl <- distance_list[[x]]
    merge(period_tbl, distance_tbl, by = "ID")
  }, dataset_list, distance_list), names(distance_list))


# Define distance groups in merged tables list
merged_list <- lapply(merged_list, function(df) {  # browser()
  
  df$distance_group <- cut(df[,11],
                           breaks = c(-Inf, limit1, limit2, Inf),
                           labels = c(paste("<=", limit1), paste(limit1, "-",limit2 ), paste(">", limit2)))
  df$distance_group <- as.factor(df$distance_group)
  df$ID <- as.factor(df$ID)
  return(df)
})

#' merge tables adding column for treatment
merged_list <- setNames(lapply(names(merged_list), function(name) {
  df <- merged_list[[name]]
  # Extract part of the name after the "_"
  parts <- str_split(name, "_")[[1]]
  well = parts[2]
  treatment <- parts[3]  #sub("^[^_]*_", "", name)
  # Add suffix as a new column
  df$well <- as.factor(well)
  df$treatment <- as.factor(treatment)
  return(df)
}), names(distance_list))

plot_folder = file.path(wd, "spatial_analysis2")

if(!dir.exists(plot_folder)){
  dir.create(plot_folder)
}

 # nested anova plots of period, amplitude, error ####
{
  period_mean_plot_nest <- nested_anova_plot(
    data = merged_table,
    metric = period,
    distance_var = distance_group,
    treatment_var = treatment,
    plot_colors = color_palette,
    ylabel = "Period (h)",
    xlabel = "Distance group (px)",
    nudge = 0.5
  )
  
  amplitude_mean_plot_nest <- nested_anova_plot(
    data = merged_table,
    metric = amplitude,
    distance_var = distance_group,
    treatment_var = treatment,
    plot_colors = color_palette,
    ylabel = "Amplitude (A.U.)",
    xlabel = "Distance group (px)",
    nudge = 8,
    step = 0.6
  )
  
  error_mean_plot_nest <- nested_anova_plot(
    data = merged_table,
    metric = error,
    distance_var = distance_group,
    treatment_var = treatment,
    plot_colors = color_palette,
    ylabel = "Error (A.U.)",
    xlabel = "Distance group (px)",
    nudge = 4,
    step = 0.4
  )
  
  savePlots(obj_to_save = list(period_mean_plot_nest = period_mean_plot_nest, amplitude_mean_plot_nest = amplitude_mean_plot_nest,
                               error_mean_plot_nest = error_mean_plot_nest), savefold = plot_folder, extension = "png", p.width = 1000)
}

 # analyse phases distances between groups ####
circ_summary = read.csv(file = file.path(wd, "circular_stats.csv"), row)
{
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
  
  savePlots(obj_to_save = phase_plot,  savefold = plot_folder, extension = "png", p.width = 1000)
  
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
  
}

 # plot vector length spread ####
{
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
  
  savePlots(obj_to_save = mu_plot,  savefold = plot_folder, extension = "png", p.width = 1000)
}

#' merge together all tables containing period and spatial data 
merged_table = data.table::rbindlist(merged_list, idcol = "sample")

meta_merged_table =
  merged_table %>%
  group_by(sample, distance_group, treatment) %>%
  summarise(n = length(ID))

 # descriptive plots averaging all samples ####

#' plot representing how many cells are present in every distance area per sample
meta_plot = ggplot2::ggplot(data = meta_merged_table, aes(x = distance_group, y = n, color = treatment))+
  geom_point(size = 3)+
  geom_line(aes(group = sample))+
  labs(
    title = "Cell in every distance area",
    x = "Distance Group",
    y = "Cell number"
  ) +
  theme_minimal()

#' plot showing period of all cells according to distance from particle
#'  TODO add different linetype for two groups and add legend

period_dotplot_trace_allcells <- ggplot()+#merged_table, aes(x = distance, y = period)) +
  # geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
  geom_smooth(data = merged_table[merged_table$treatment == "FRED", ], aes(x = distance, y = period, colour = distance_group, group = distance_group), method = "loess", se =  TRUE, level = 0.95) +
  geom_smooth(data = merged_table[merged_table$treatment == "MUT", ], aes(x = distance, y = period, colour = distance_group, group = distance_group), method = "loess", se =  TRUE, level = 0.95) +
  scale_color_manual(values = color_palette)+
  ylim(22.5, 24.5) +
  theme_minimal() +
  labs(title = "Period by Distance", x = "Distance from Particle (px)", y = "Period (h)")

sample_weights <- merged_table %>%
  group_by(sample, distance_group, treatment) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  mutate(weight = 1 / cell_count)

#' wrong assignments of weights
merged_table_w = merged_table %>%
  left_join(sample_weights %>% select(sample, weight), by = "sample")

period_trace_w <- ggplot()+
  geom_smooth(data = merged_table_w[merged_table_w$treatment == "FRED", ], aes(x = distance, y = period, colour = distance_group, group = distance_group, weight = weight), method = "loess", se =  TRUE, level = 0.95) +
  geom_smooth(data = merged_table_w[merged_table_w$treatment == "MUT", ], aes(x = distance, y = period, colour = distance_group, group = distance_group, weight = weight), method = "loess", se =  TRUE, level = 0.95) +
  scale_color_manual(values = color_palette)+
  ylim(22.5, 24.5) +
  theme_minimal() +
  labs(title = "Period by Distance", x = "Distance from Particle (px)", y = "Period (h)")


#' plot showing period of all cells according to distance from particle
#'  TODO add different linetype for two groups and add legend

amplitude_dotplot_trace_allcells <- ggplot()+#merged_table, aes(x = distance, y = period)) +
  # geom_point(aes(color = distance_group, alpha = 0.01, stroke = NA)) +
  geom_smooth(data = merged_table[merged_table$treatment == "FRED", ], aes(x = distance, y = amplitude, colour = distance_group, group = distance_group), method = "loess", se =  TRUE, level = 0.95) +
  geom_smooth(data = merged_table[merged_table$treatment == "MUT", ], aes(x = distance, y = amplitude, colour = distance_group, group = distance_group), method = "loess", se =  TRUE, level = 0.95) +
  scale_color_manual(values = color_palette)+
  theme_minimal() +
  labs(title = "Period by Distance", x = "Distance from Particle (px)", y = "Period (h)")

 # aligned traces all together ####
 
traces_list <- dataset_list$align.trx
# merge all tables together adding info about sample
 for(i in 1:length(traces_list)){
  #import dataset
  traces <- traces_list[[i]]
  name = names(traces_list)[i]
  # isolate treatment
  treatments = c("FRED", "WT", "MUT")
  match <- sapply(treatments, function(y) grepl(y, name))
  treatment = treatments[match]
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
  
# bind tables
averaged_traces_bound = do.call(rbind, traces_summary_list) %>% `rownames<-`(NULL)

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

# plots
weighted_mean_plot = ggplot2::ggplot(traces_summary_all_samples_weighted, aes(group = group))+
  geom_ribbon(aes(x = t, ymin = average-se, ymax = average+se, alpha = 0.5, colour = group, fill = group))+
  geom_line(aes(x = t, y = average, colour = group))+
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

# perform alignment to common trace
re_aligned_traces <- 
  #setNames(
    lapply(
      names(dataset_list$align.trx.wide),
      function(x){
        browser()
        period.tbl = computePeriod(ch2_cells, filename, excludeNC = TRUE, top = 27, bottom = 16, rm.start = 0)
        phase_align_trace(dataset_list$align.trx.wide[[x]], period.tbl)
        },
      dataset_list)#,
    #names(dataset_list$align.trx.wide))

# identify length of traces
trx_len = setNames(lapply(re_aligned_traces, function(x){dim(x)[1]}), names(re_aligned_traces))

# trim all traces to common length
if(length(unique(unlist(trx_len))) > 1){
  re_aligned_traces <- 
    setNames(
      lapply(
        names(re_aligned_traces),
        function(x){
          # get difference with current trace table
          difference <- trx_len[[x]] - min(unlist(trx_len))
          if(difference > 0){
            # trim cutting half from start and half from end
            trimStart = round(difference/2)
            diff2 = difference - trimStart
            trimEnd = dim(re_aligned_traces[[x]])[2] - diff2
            traces_trim = re_aligned_traces[[x]][,trimStart:trimEnd]
            }
          return(traces_trim)
          },
        trx_len, re_aligned_traces),
      names(re_aligned_traces))
  }

# make all wide traces into long format
list_traces_long = setNames(lapply(names(dataset_list$align.trx.wide), function(x){
  traces = re_aligned_traces[[x]]
  traces_long = tidyr::pivot_longer(traces,
                                            cols = colnames(aligned_traces_t)[-1],
                                            names_to = "ID",
                                            values_to = "intensity")
  traces_long$ID = as.numeric(traces_long$ID)
  traces_long$sample = x
  
  # assign treatment
  treatments = 
  match_found <- sapply(c("FRED", "WT", "MUT"), function(y) grepl(y, x))
  matching_patterns <- treatments[match_found]
  # assign group values
  
}, re_aligned_traces), names(re_aligned_traces))

# merge all tables into a single one 
all_traces_aligned = do.call(rbind, list_traces_long)

# plot



 # import all tables 
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

 # OLD CODE ####

#' also nested anova but different

period_mean_plot_nest_all <- metric_mean_plot_nested_allpoints(
  data = merged_table,
  metric = period,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Period (h)",
  xlabel = "Distance group"
)

amplitude_mean_plot_nest_all <- metric_mean_plot_nested_allpoints(
  data = merged_table,
  metric = amplitude,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Amplitude (A.U.)",
  xlabel = "Distance group"
)

error_mean_plot_nest_all <- metric_mean_plot_nested_allpoints(
  data = merged_table,
  metric = error,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Error (A.U.)",
  xlabel = "Distance group"
)
#
# period_barplot_across_nest <- metric_barplot_across_nested(
#   data = merged_table,
#   metric = period,
#   distance_var = distance_group,
#   treatment_var = treatment,
#   plot_colors = color_palette,
#   ylabel = "Period (h)",
#   xlabel = "Treatment"
# )
#
# amplitude_barplot_across_nest <- metric_barplot_across_nested(
#   data = merged_table,
#   metric = amplitude,
#   distance_var = distance_group,
#   treatment_var = treatment,
#   plot_colors = color_palette,
#   ylabel = "Amplitude (A.U.)",
#   xlabel = "Treatment"
# )
#
# error_barplot_across_nest <- metric_barplot_across_nested(
#   data = merged_table,
#   metric = error,
#   distance_var = distance_group,
#   treatment_var = treatment,
#   plot_colors = color_palette,
#   ylabel = "Error (A.U.)",
#   xlabel = "Treatment"
# )


#' boxplot all together

period_boxplot_across <- metric_boxplot_across(
  data = merged_table,
  metric = period,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Period (h)",
  xlabel = "Distance group"
)

amplitude_boxplot_across <- metric_boxplot_across(
  data = merged_table,
  metric = amplitude,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Amplitude (A.U.)",
  xlabel = "Distance group"
)

error_boxplot_across <- metric_boxplot_across(
  data = merged_table,
  metric = error,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Error (A.U.)",
  xlabel = "Distance group"
)

#' plot mean and error across samples nested

period_mean_plot_across_nest <- metric_mean_plot_across_nested(
  data = merged_table,
  metric = period,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Period (h)",
  xlabel = "Distance group"
)

amplitude_mean_plot_across_nest <- metric_mean_plot_across_nested(
  data = merged_table,
  metric = amplitude,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Amplitude (A.U.)",
  xlabel = "Distance group"
)

error_mean_plot_across_nest <- metric_mean_plot_across_nested(
  data = merged_table,
  metric = error,
  distance_var = distance_group,
  treatment_var = treatment,
  plot_colors = color_palette,
  ylabel = "Error (A.U.)",
  xlabel = "Distance group"
)


metric_boxplot_across <- function(data, metric, distance_var, treatment_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, ...) {
  
  # Compute normality test if needed
  if (norm_test) {
    normality <- data %>%
      group_by({{distance_var}}, {{treatment_var}}) %>%
      rstatix::shapiro_test({{metric}})
    print(normality)
  }
  
  # Define formula for repeated measures ANOVA with distance and treatment variables
  formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(distance_var)), "*", deparse(substitute(treatment_var))))
  res.aov <- rstatix::anova_test(data = data, formula = formula, wid = ID, within = c({{distance_var}}, {{treatment_var}}))
  
  # Perform pairwise comparisons for treatment within each distance group
  formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(treatment_var))))
  pwc <- data %>%
    group_by({{distance_var}}) %>%
    rstatix::pairwise_t_test(formula = formula,
                             paired = FALSE,
                             p.adjust.method = test)
  
  # Add XY position to the p-values for proper placement on the plot
  pwc <- pwc %>% rstatix::add_xy_position(x = as_label(enquo(distance_var)), dodge = 0.8, fun = "max", step.increase = 0.1)
  
  # Generate the plot
  plot <- ggpubr::ggboxplot(data,
                            x = deparse(substitute(distance_var)),
                            y = deparse(substitute(metric)),
                            color = deparse(substitute(treatment_var)),
                            palette = plot_colors) +
    ggpubr::stat_pvalue_manual(pwc, dodge = 0.8, bracket.nudge.y = 1) + # Use dodge for stat_pvalue_manual only
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title = paste("Comparison of", deparse(substitute(metric)), "by Distance and Treatment"),
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = rstatix::get_pwc_label(pwc),
      y = ylabel,
      x = xlabel
    )
  
  return(plot)
}

metric_boxplot_across_nested <- function(data, metric, distance_var, well = well, treatment_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, ...) {
  
  # Compute normality test if needed
  if (norm_test) {
    normality <- data %>%
      group_by(.data[[deparse(substitute(distance_group))]], .data[[deparse(substitute(treatment))]], .data[[deparse(substitute(well))]]) %>%
      rstatix::shapiro_test({{metric}})
    print(normality)
  }
  
  # Define formula for nested ANOVA
  formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(distance_var)), "/", deparse(substitute(well))))
  
  # Run the nested ANOVA
  res.aov <- rstatix::anova_test(data = data, formula = formula)
  
  # Perform pairwise comparisons for distance groups
  formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(treatment_var))))
  pwc <- data %>%
    group_by({{distance_var}}) %>%
    rstatix::pairwise_t_test(formula = formula,
                             paired = FALSE,
                             p.adjust.method = test)
  
  # Add XY position for p-values
  pwc <- pwc %>%
    rstatix::add_xy_position(x = as_label(enquo(distance_var)), dodge = 0.8, fun = "max", step.increase = 0.4)
  
  # Generate the plot
  plot <- ggpubr::ggviolin(
    data,
    x = deparse(substitute(distance_group)),
    y = deparse(substitute(metric)),
    color = deparse(substitute(treatment)),
    palette = plot_colors
  ) +
    ggpubr::stat_pvalue_manual(pwc, dodge = 0.8, bracket.nudge.y = 1) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title = paste("Comparison of", deparse(substitute(metric)), "by Distance Group"),
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = rstatix::get_pwc_label(pwc),
      y = ylabel,
      x = xlabel
    )
  
  return(plot)
}

metric_mean_plot_across_nested = function(data, metric, distance_var, well = well, treatment_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, ...) {
  
  # Compute normality test if needed
  if (norm_test) {
    normality <- data %>%
      group_by(.data[[deparse(substitute(distance_group))]], .data[[deparse(substitute(treatment))]], .data[[deparse(substitute(well))]]) %>%
      rstatix::shapiro_test({{metric}})
    print(normality)
  }
  
  # Define formula for nested ANOVA
  formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(distance_var)), "/", deparse(substitute(well))))
  
  # Run the nested ANOVA
  res.aov <- rstatix::anova_test(data = data, formula = formula)
  
  # Perform pairwise comparisons for treatment groups within each distance group
  formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(treatment_var))))
  pwc <- data %>%
    group_by({{distance_var}}) %>%
    rstatix::pairwise_t_test(formula = formula,
                             paired = FALSE,
                             p.adjust.method = test)
  
  # Add XY position for p-values
  pwc <- pwc %>%
    rstatix::add_xy_position(x = as_label(enquo(distance_var)), dodge = 0.8, fun = "mean", step.increase = 0.1)
  
  # Calculate mean and standard error for each group
  summary_stats <- data %>%
    dplyr::group_by({{distance_var}}, {{treatment_var}}) %>%
    dplyr::summarise(
      mean = mean({{metric}}, na.rm = TRUE),
      se = sd({{metric}}, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Generate the plot with mean and standard error
  plot <- ggplot(summary_stats, aes(x = {{distance_var}}, y = mean, color = {{treatment_var}}, group = {{treatment_var}})) +
    geom_point(size = 4, position = position_dodge(0.8)) +  # Add points for means
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, size = 0.8, position = position_dodge(0.8)) +  # Add error bars
    scale_color_manual(values = plot_colors) +  # Apply custom colors
    labs(
      title = paste("Comparison of", deparse(substitute(metric)), "by Distance Group"),
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = rstatix::get_pwc_label(pwc),
      y = ylabel,
      x = xlabel
    ) +
    ggpubr::stat_pvalue_manual(pwc, dodge = 0.8) +  # Add p-values
    theme_minimal() +
    theme(
      plot.subtitle = element_text(size = 10, hjust = 0),
      legend.title = element_blank()
    )
  
  return(plot)
}

metric_mean_plot_nested_allpoints = function(data, metric, distance_var, well = well, treatment_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, ...) {
  
  # Compute normality test if needed
  if (norm_test) {
    normality <- data %>%
      group_by({{distance_var}}, {{treatment_var}}, {{well}}) %>%
      rstatix::shapiro_test({{metric}})
    print(normality)
  }
  # Calculate mean and SE for each well
  well_stats <- data %>%
    group_by({{distance_var}}, {{treatment_var}}, {{well}}) %>%
    summarise(
      well_mean = mean({{metric}}, na.rm = TRUE),
      well_se = sd({{metric}}, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Calculate absolute mean and SE across all wells
  overall_stats <- well_stats %>%
    group_by({{distance_var}}, {{treatment_var}}) %>%
    summarise(
      abs_mean = mean(well_mean, na.rm = TRUE),
      abs_se = sd(well_mean, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Perform pairwise comparisons if needed
  if (comparison) {
    formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(treatment_var))))
    pwc <- data %>%
      group_by({{distance_var}}) %>%
      rstatix::pairwise_t_test(
        formula = formula,
        paired = FALSE,
        p.adjust.method = test
      )
    # Add XY position for p-values
    pwc <- pwc %>%
      rstatix::add_xy_position(
        x = as_label(enquo(distance_var)),
        dodge = 0.8,
        fun = "mean_se",
        step.increase = 0.2
      )
  }
  
  # Generate the plot
  plot <- ggplot() +
    # Add well-level means
    geom_point(
      data = well_stats,
      aes(x = {{distance_var}}, y = well_mean, color = {{treatment_var}}),
      position = position_dodge(0.8),
      size = 2, # Smaller points for wells
      alpha = 0.5 # Transparency to make them distinct
    ) +
    # Add overall means with error bars
    geom_point(
      data = overall_stats,
      aes(x = {{distance_var}}, y = abs_mean, color = {{treatment_var}}),
      position = position_dodge(0.8),
      size = 4 # Larger points for absolute means
    ) +
    geom_errorbar(
      data = overall_stats,
      aes(
        x = {{distance_var}},
        ymin = abs_mean - abs_se,
        ymax = abs_mean + abs_se,
        color = {{treatment_var}}
      ),
      width = 0.2,
      position = position_dodge(0.8)
    ) +
    # Add pairwise comparisons (optional)
    ggpubr::stat_pvalue_manual(
      pwc,
      dodge = 0.8,
    ) +
    # Customize colors and labels
    scale_color_manual(values = plot_colors) +
    labs(
      title = paste("Comparison of", deparse(substitute(metric)), "by Distance Group"),
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = rstatix::get_pwc_label(pwc),
      y = ylabel,
      x = xlabel
    ) +
    theme_minimal()
  
  return(plot)
}

metric_barplot_across_nested <- function(data, metric, distance_var, well = well, treatment_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, ...) {
  
  # Compute normality test if needed
  if (norm_test) {
    normality <- data %>%
      group_by(.data[[deparse(substitute(distance_group))]], .data[[deparse(substitute(treatment))]], .data[[deparse(substitute(well))]]) %>%
      rstatix::shapiro_test({{metric}})
    print(normality)
  }
  
  # Define formula for nested ANOVA
  formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(distance_var)), "/", deparse(substitute(well))))
  
  # Run the nested ANOVA
  res.aov <- rstatix::anova_test(data = data, formula = formula)
  
  # Perform pairwise comparisons for distance groups
  formula <- as.formula(paste(deparse(substitute(metric)), "~", deparse(substitute(treatment_var))))
  pwc <- data %>%
    group_by({{distance_var}}) %>%
    rstatix::pairwise_t_test(formula = formula,
                             paired = FALSE,
                             p.adjust.method = test)
  
  # Add XY position for p-values
  pwc <- pwc %>%
    rstatix::add_xy_position(x = as_label(enquo(distance_var)), dodge = 0.8, fun = "max", step.increase = 0.4)
  
  # Generate the plot
  plot <- ggpubr::gghistogram(
    data,
    x = deparse(substitute(distance_group)),
    y = deparse(substitute(metric)),
    color = deparse(substitute(treatment)),
    palette = plot_colors
  ) +
    ggpubr::stat_pvalue_manual(pwc, dodge = 0.8, bracket.nudge.y = 1) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title = paste("Comparison of", deparse(substitute(metric)), "by Distance Group"),
      subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
      caption = rstatix::get_pwc_label(pwc),
      y = ylabel,
      x = xlabel
    )
  
  return(plot)
}
