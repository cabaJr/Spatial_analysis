#' script to compare different datasets obtained with the collate_tau_spatial script

#### compare across weeks ####

# to clean the environment preserving functions
rm(list = setdiff(ls(), lsf.str()))

#' get location of datasets
wd_per1 = r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Proj_Tau\Tau_an_dev\wk2_hemislices\combined_all)"
period1 = "wk2"

wd_per2 = r"(C:\Users\mf420\UK Dementia Research Institute Dropbox\Brancaccio Lab\Marco F\Proj_Tau\Tau_an_dev\wk5_hemislices\combined_all)"
period2 = "wk5"

# import datasets
merged_list1 <- readRDS(file.path(wd_per1, paste0(period1, "_merged_list.rds")))
merged_list2 <- readRDS(file.path(wd_per2, paste0(period2, "_merged_list.rds")))

merged_table1 <- readRDS(file.path(wd_per1, paste0(period1, "_merged_table.rds")))
merged_table2 <- readRDS(file.path(wd_per2, paste0(period2, "_merged_table.rds")))

dataset_combined1 <- readRDS(file.path(wd_per1, paste0(period1, "_dataset.rds")))
dataset_combined2 <- readRDS(file.path(wd_per2, paste0(period2, "_dataset.rds")))

traces_annotated1 <- readRDS(file.path(wd_per1, paste0(period1, "_traces_annotated.rds")))
traces_annotated2 <- readRDS(file.path(wd_per2, paste0(period2, "_traces_annotated.rds")))

traces_summary_list1 <- readRDS(file.path(wd_per1, paste0(period1, "_traces_summary_list.rds")))
traces_summary_list2 <- readRDS(file.path(wd_per2, paste0(period2, "_traces_summary_list.rds")))

# add week column to each table in the list
merged_list1 <- lapply(merged_list1, function(x) {
  x$week <- period1
  return(x)
})

merged_list2 <- lapply(merged_list2, function(x) {
  x$week <- period2
  return(x)
})

merged_table1$week <- period1
merged_table2$week <- period2

# generate merged table
merged_table <- rbind(merged_table1, merged_table2)

# generate summary merged tables
merged_table1_summary <- merged_table1 %>% 
  group_by(treatment, well, distance_group, week) %>% 
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
                   SD_green = mean(SD_green),
                   week = unique(week), .groups = "drop")

merged_table2_summary <- merged_table2 %>% 
  group_by(treatment, well, distance_group, week) %>% 
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
                   SD_green = mean(SD_green),
                   week = unique(week), .groups = "drop")
  
merged_table_summary <- rbind(merged_table1_summary, merged_table2_summary)

# generate plots ####

# compare AUC change between weeks:
comparison_AUC <- ggplot(merged_table_summary, aes(y = factor(distance_group), x = AUC_mean, 
                                 color = factor(week), group = interaction(well, distance_group))) +
  geom_point(position = position_dodge(width = 0.5)) +  # Spread points within distance_group
  geom_line(aes(y = distance_group, x = AUC_mean, group = interaction(well, distance_group)), 
            position = position_dodge(width = 0.5), alpha = 0.5) + # Connecting lines
  facet_wrap(~ treatment) +  # Separate plots for each treatment
  labs(y = "Distance Group", x = "AUC Mean", color = "Week") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))



# compare traces between weeks

# compare period difference

# Create the plot
df <- merged_table_summary %>%
  mutate(distance_week = interaction(distance_group, week, sep = "_"))  # Create combined category

# Create the plot

# Create the plot
ggplot(merged_table_summary, aes(x = factor(week), y = AUC_mean, 
               color = treatment, group = interaction(well, distance_group), 
               shape = factor(distance_group))) +  # Different shapes for distance groups
  geom_point() +  # Spread points slightly
  geom_line(alpha = 0.5) +  # Connect values across weeks
  facet_wrap(treatment ~ distance_group) +  # Separate plots for each treatment
  labs(x = "Week", y = "AUC Mean", color = "Treatment", shape = "Distance Group") +
  scale_color_manual(values = c("black", "red")) +  # Set custom colors
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
