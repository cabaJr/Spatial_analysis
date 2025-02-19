pkg_req <- c("ggplot2", "reshape2", "gridExtra")

# Install and load packages if not already installed
for (package in pkg_req) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Example lists (replace with actual tables)
list1 <- merged_list
list2 <- merged_list_align

# Define output PDF
pdf("comparison_heatmaps.pdf", width = 18, height = 12)

# Iterate over tables that exist in both lists
common_tables <- intersect(names(list1), names(list2))
for (table_name in common_tables) {
  df1 <- list1[[table_name]]
  df2 <- list2[[table_name]]
  
  # Merge tables on ID to ensure alignment
  merged_df <- merge(df1, df2, by = "ID", suffixes = c("_1", "_2"))
  
  # Compute differences for selected columns
  diff_df <- data.frame(ID = merged_df$ID,
                        distance_group = merged_df$distance_group_1,
                        period = merged_df$period_1 - merged_df$period_2,
                        amplitude = merged_df$amplitude_1 - merged_df$amplitude_2,
                        RAE = merged_df$RAE_1 - merged_df$RAE_2)
  
  # Reshape for ggplot
  diff_melt <- melt(diff_df, id.vars = c("ID", "distance_group"))
  
  # Compute summary statistics
  summary_stats <- aggregate(value ~ variable, data = diff_melt, FUN = function(x) c(median = median(x), sd = sd(x)))
  summary_stats <- do.call(data.frame, summary_stats)
  names(summary_stats) <- c("Variable", "Median", "SD")
  
  # Create separate plots for each variable with individual legends
  plot_list <- list()
  variables <- unique(diff_melt$variable)
  
  for (var in variables) {
    plot_list[[var]] <- ggplot(subset(diff_melt, variable == var), aes(x = factor(distance_group), y = factor(ID), fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_minimal() +
      labs(title = paste("Heatmap of Differences -", table_name, "-", var), x = "Distance Group", y = "ID", fill = "Difference")
  }
  
  # Convert summary stats to grob for grid.arrange
  table_grob <- tableGrob(summary_stats)
  
  # Arrange plots in a grid with the summary table
  grid.arrange(grobs = c(plot_list, list(table_grob)), layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)), heights = c(5, 1))
}

# Close PDF
dev.off()