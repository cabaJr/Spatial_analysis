
#' nested_anova_plot
#' 
#' @description
#' function to generate plot and perform statistical testing
#' 
#' @param data 
#' @param metric 
#' @param distance_var 
#' @param well 
#' @param treatment_var 
#' @param ylabel 
#' @param xlabel 
#' @param comparison 
#' @param test 
#' @param norm_test 
#' @param plot_colors 
#' @param nudge 
#' @param step 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
nested_anova_plot <- function(data, metric, distance_var, well = well, treatment_var, ylabel = "", xlabel = "", comparison = TRUE, test = "bonferroni", norm_test = TRUE, plot_colors, nudge = 0.5, step = 0.2, custom_levels, ...) {
  # update levels if necessary
  
  # if(!is.null(custom_levels)){
  #   levels(data$treatment) <- custom_levels
  # }
  
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
      well_n = n(),
      .groups = "drop"
    )
  
  # Calculate overall mean and SE across wells
  overall_stats <- well_stats %>%
    left_join(
      data %>% group_by({{distance_var}}, {{treatment_var}}, {{well}}) %>% summarise(n_points = n(), .groups = "drop"),
      by = c(as_label(enquo(distance_var)), as_label(enquo(treatment_var)), as_label(enquo(well)))
    ) %>%
    group_by({{distance_var}}, {{treatment_var}}) %>%
    summarise(
      abs_mean = sum(well_mean * n_points, na.rm = TRUE) / sum(n_points, na.rm = TRUE),
      abs_se = sqrt(sum((n_points * (well_mean - abs_mean)^2), na.rm = TRUE) / sum(n_points, na.rm = TRUE)),
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
      aes(x = {{distance_var}}, y = well_mean, color = {{treatment_var}}, size = well_n),
      position = position_dodge(0.8),
      alpha = 0.3 # Transparency to distinguish them
    ) +
    # Add overall means with error bars
    geom_point(
      data = overall_stats,
      aes(x = {{distance_var}}, y = abs_mean, color = {{treatment_var}}),
      position = position_dodge(0.8),
      size = 4, # Larger points for overall means
      shape = 18
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
    # ggpubr::stat_pvalue_manual(
    #   pwc2,
    #   dodge = 0.8,
    #   bracket.nudge.y = nudge
    # )+
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
  returns <- list(plot = plot,
                  well_stats = well_stats)
  return(returns)
}

# function to align traces to the mean phase of the group
#' Title
#'
#' @param traces_table 
#' @param period_table 
#' @param remove_start 
#' @param remove_end 
#' @param align_to 
#' @param debug 
#'
#' @return
#' @export
#'
#' @examples
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

#' function to save plots to a specific location and format
#'
#' @param obj_to_save 
#' @param basepath 
#' @param filename 
#' @param extension 
#' @param p.width 
#' @param p.height 
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

#' function to trim tables and remove rows with NAs
#'
#' @param df 
#' @param rows 
trim_tbl = function(df, rows){
  
  # Ensure df is a dataframe
  if (!is.data.frame(df)) {
    stop("Input must be a dataframe")
  }
  
  n <- nrow(df)
  
  if (n <= 24) {
    warning("Dataframe has 24 or fewer rows; checking all rows for NAs")
    df <- df[complete.cases(df), , drop = FALSE]
  } else {
    # Identify rows to check (first 12 and last 12)
    first_n <- 1:rows
    last_n <- (n-(rows-1)):n
    
    # Find rows containing NAs in the first and last 12 rows
    rows_to_remove <- c(
      first_n[rowSums(is.na(df[first_n, ])) > 0],
      last_n[rowSums(is.na(df[last_n, ])) > 0]
    )
    
    # Remove identified rows
    df <- df[-rows_to_remove, , drop = FALSE]
  }
  
  # Remove columns that contain only NA values
  df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]
  
  return(df)
}

#' function to generate a plot comparing the percentage of cells in each distance group and period length
#'
#' @param data 
#' @param filter 
#' @param x 
#' @param y 
#' @param fill 
#' @param colour 
#' @param col_scheme 
#' @param shape 
#' @param alpha1 
#' @param alpha2 
#' @param level1 
#' @param level2 
#' @param size1 
#' @param linetype 
#' @param user_ylims 
#' @param title 
#' @param xtitle 
#' @param ytitle 
#' @param minor_grid_x 
#' @param y_breaks 
#' @param minor_grid_y 
meta_plot_draw = function(data, filter, x, y, fill, colour, col_scheme, shape, alpha1 = 0.15, alpha2 = 0.15, 
                          level1 = 0.95, level2 = 0, size1 = 2.5, linetype, user_ylims = NULL, title, xtitle, ytitle, minor_grid_x = NULL, y_breaks = NULL, minor_grid_y = NULL){
  
  plot = ggplot2::ggplot()+
    geom_point(data = data %>% filter(distance_group == filter), aes(x = {{x}}, y = {{y}}, colour = {{colour}}, shape = {{shape}}, size = size1))+
    ggplot2::scale_size_identity()+
    stat_ellipse(data = data %>% filter(distance_group == filter), aes(x = {{x}}, y = {{y}}, fill = {{fill}}, colour = NULL, linetype = {{linetype}}), type = "t", level = level1, alpha = alpha1, geom = "polygon") +
    geom_smooth(data = data %>% filter(distance_group == filter), aes(x = {{x}}, y = {{y}}, fill = {{fill}}, colour = {{colour}}, linetype = {{linetype}}), method = lm, level = level2, alpha = alpha1)+
    scale_color_manual(values = col_scheme)+
    scale_fill_manual(values = col_scheme)+
    labs(title = title, x = xtitle, y = ytitle)+
    theme_minimal()+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))
  
  if(!is.null(user_ylims)){
    plot = plot + ylim(user_ylims)
  }
  
  if(is.null(minor_grid_x)){
    plot = plot + theme(panel.grid.minor.x = element_blank())
  }
  if(is.null(minor_grid_y)){
    plot = plot + theme(panel.grid.minor.y = element_blank())
  }
  
  if(!is.null(y_breaks)){
    plot = plot + scale_y_continuous(breaks = y_breaks)
  }
  
  
  return(plot)
}
#' function to generate a 3D plot of AUC values
#'
#' @param data 
#' @param x_var 
#' @param y_var 
#' @param z_var 
#' @param color_var 
#' @param palette 
#' @param continuous 
#' @param x_lab 
#' @param y_lab 
#' @param z_lab 
#' @param title 
#' @param eye 
#' @param zoom 
threeDplot_AUC = function(data, x_var, y_var, z_var, color_var, palette = NULL, continuous = FALSE, x_lab = "x", y_lab = "y", z_lab = "z", title = "", eye = "pos1", zoom = 1){
  
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
  
  if(is.null(palette)){
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
  
  plotly <- plot_ly() %>%
    # Add scatter plot points
    add_trace(
      data = data %>% filter(treatment == "FRED"),
      x = ~distance,
      y = ~AUC,
      z = ~period,
      type = "scatter3d",
      mode = "markers",
      split = ~distance_group,  # Ensures proper grouping by discrete variable
      marker = list(symbol = "circle", size = 4),
      colors = palette,  # Assign custom colors to factor levels
      showlegend = TRUE) %>%
    # Add scatter plot points for second group
    add_trace(
      data = data %>% filter(treatment == "MUT"),
      x = ~distance,
      y = ~AUC,
      z = ~period,
      type = "scatter3d",
      mode = "markers",
      split = ~distance_group,  # Ensures proper grouping by discrete variable
      marker = list(symbol = "diamond", size = 4),
      colors = palette,  # Assign custom colors to factor levels
      showlegend = TRUE) %>%
    # Layout settings
    layout(
      scene = list(
        xaxis = list(title = x_lab, titlefont = list(size = 14), range = c(50,0)),
        yaxis = list(title = y_lab, titlefont = list(size = 14), range = c(0,50000)),
        zaxis = list(title = z_lab, titlefont = list(size = 14), range = c(21,26)),
        camera = list(eye = eye)
      ),
      title = list(
        text = title,
        font = list(size = 16)
      )
    )
  return(plotly)
}

#' function to generate a 3D plot of phase values
#'
#' @param data 
#' @param x_var 
#' @param y_var 
#' @param z_var 
#' @param color_var 
#' @param palette 
#' @param continuous 
#' @param x_lab 
#' @param y_lab 
#' @param z_lab 
#' @param title 
#' @param eye 
#' @param zoom 
threeDplot_phase = function(data, x_var, y_var, z_var, color_var, palette = NULL, continuous = FALSE, x_lab = "x", y_lab = "y", z_lab = "z", title = "", eye = "pos1", zoom = 1){
  
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
  
  if(is.null(palette)){
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
  
  plotly <- plot_ly() %>%
    # Add scatter plot points
    add_trace(
      data = data, #%>% filter(treatment == "FRED"),
      x = ~phase_norm,
      y = ~AUC,
      z = ~period,
      type = "scatter3d",
      mode = "markers",
      split = ~distance_group,  # Ensures proper grouping by discrete variable
      marker = list(size = 3, opacity = 0.5),
      colors = palette,  # Assign custom colors to factor levels
      showlegend = TRUE) %>%
    # Add scatter plot points
    # add_trace(
    #   data = data %>% filter(treatment == "MUT"),
    #   x = ~phase_norm,
    #   y = ~AUC,
    #   z = ~period,
    #   type = "scatter3d",
    #   mode = "markers",
    #   split = ~distance_group,  # Ensures proper grouping by discrete variable
    #   marker = list(size = 3, opacity = 0.5),
    #   colors = palette,  # Assign custom colors to factor levels
    #   showlegend = TRUE) %>%
    # Layout settings
    layout(
      scene = list(
        xaxis = list(title = x_lab, titlefont = list(size = 14), range = c(12, -12)),
        yaxis = list(title = y_lab, titlefont = list(size = 14), range = c(0,50000)),
        zaxis = list(title = z_lab, titlefont = list(size = 14), range = c(21,26)),
        camera = list(eye = eye)
      ),
      title = list(
        text = title,
        font = list(size = 16)
      )
    )
  return(plotly)
}
