# script for changing the merged table object
merged_tbl_FRED = dplyr::filter(merged_table, treatment == "FRED")

merged_tbl_FRED_close_val = dplyr::filter(merged_tbl_FRED, distance_group == "<= 3") %>% 
  dplyr::select(c(sample, period, distance)) %>% 
  tidyr::pivot_wider(names_from = sample, values_from = period, values_fn = mean)

merged_tbl_FRED_close_n = dplyr::filter(merged_tbl_FRED, distance_group == "<= 3") %>% 
  dplyr::select(c(sample, period, distance)) %>% 
  tidyr::pivot_wider(names_from = sample, values_from = period, values_fn = length)

merged_tbl_FRED_close_sd = dplyr::filter(merged_tbl_FRED, distance_group == "<= 3") %>% 
  dplyr::select(c(sample, period, distance)) %>% 
  tidyr::pivot_wider(names_from = sample, values_from = period, values_fn = sd)

table_list = list(merged_tbl_FRED_close_val, merged_tbl_FRED_close_n, merged_tbl_FRED_close_sd)

merged_tbl_FRED_close<- do.call(cbind, lapply(seq_len(10), function(i) {
  do.call(rbind, lapply(table_list, function(tbl) tbl[i, , drop = FALSE]))
}))

colnames1 = colnames(merged_tbl_FRED_close_n)[-1]
fields = c("_val", "_n", "_sd")
combinations = expand.grid(fields, colnames1)
colnames2 = c("distance", paste0(combinations$Var2, combinations$Var1))
merged_tbl_FRED_close<- do.call(cbind, do.call(Map, c(f = data.frame, table_list))) %>% 
  .[-c(2, 3)] %>% `colnames<-`(colnames2) %>% arrange(distance)

arrange_tbl = function(master_tbl, sub1, sub2, var){
  filtered1 = dplyr::filter(master_tbl, treatment == sub1)
  
  filtered_val = dplyr::filter(filtered1, distance_group == sub2) %>% 
    dplyr::select(c(sample, var, distance)) %>% 
    tidyr::pivot_wider(names_from = sample, values_from = var, values_fn = mean)
  
  filtered_n = dplyr::filter(filtered1, distance_group == sub2) %>% 
    dplyr::select(c(sample, var, distance)) %>% 
    tidyr::pivot_wider(names_from = sample, values_from = var, values_fn = length)
  
  filtered_sd = dplyr::filter(filtered1, distance_group == sub2) %>% 
    dplyr::select(c(sample, var, distance)) %>% 
    tidyr::pivot_wider(names_from = sample, values_from = var, values_fn = sd)
  
  table_list = list(filtered_val, filtered_sd, filtered_n)
  
  merged_filtered<- do.call(cbind, lapply(seq_len(10), function(i) {
    do.call(rbind, lapply(table_list, function(tbl) tbl[i, , drop = FALSE]))
  }))
  
  colnames1 = colnames(filtered_n)[-1]
  fields = c("_val", "_sd", "_n")
  combinations = expand.grid(fields, colnames1)
  colnames2 = c("distance", paste0(combinations$Var2, combinations$Var1))
  merged_filtered<- do.call(cbind, do.call(Map, c(f = data.frame, table_list))) %>% 
    .[-c(2, 3)] %>% `colnames<-`(colnames2) %>% arrange(distance)
  return(merged_filtered)
}

FRED_close_period = arrange_tbl(merged_table, sub1 = "FRED", sub2 = "<= 3", var = "period")
FRED_mid_period = arrange_tbl(merged_table, sub1 = "FRED", sub2 = "3 - 12", var = "period")
FRED_far_period = arrange_tbl(merged_table, sub1 = "FRED", sub2 = "> 12", var = "period")

MUT_close_period = arrange_tbl(merged_table, sub1 = "MUT", sub2 = "<= 3", var = "period")
MUT_mid_period = arrange_tbl(merged_table, sub1 = "MUT", sub2 = "3 - 12", var = "period")
MUT_far_period = arrange_tbl(merged_table, sub1 = "MUT", sub2 = "> 12", var = "period")

FRED_close_amplitude = arrange_tbl(merged_table, sub1 = "FRED", sub2 = "<= 3", var = "amplitude")
FRED_mid_amplitude = arrange_tbl(merged_table, sub1 = "FRED", sub2 = "3 - 12", var = "amplitude")
FRED_far_amplitude = arrange_tbl(merged_table, sub1 = "FRED", sub2 = "> 12", var = "amplitude")

MUT_close_amplitude = arrange_tbl(merged_table, sub1 = "MUT", sub2 = "<= 3", var = "amplitude")
MUT_mid_amplitude = arrange_tbl(merged_table, sub1 = "MUT", sub2 = "3 - 12", var = "amplitude")
MUT_far_amplitude = arrange_tbl(merged_table, sub1 = "MUT", sub2 = "> 12", var = "amplitude")

wd = "C:\\Users\\mf420\\UK Dementia Research Institute Dropbox\\Brancaccio Lab\\Marco F\\Proj_Tau\\Tau_an_dev\\wk5_hemislices\\Left2/spatial_analysis"

tables_vec = c("FRED_close_period", "FRED_mid_period", "FRED_far_period", 
               "MUT_close_period", "MUT_mid_period", "MUT_far_period", 
               "FRED_close_amplitude", "FRED_mid_amplitude", "FRED_far_amplitude", 
               "MUT_close_amplitude", "MUT_mid_amplitude", "MUT_far_amplitude")
tables_vec_obj = c(FRED_close_period, FRED_mid_period, FRED_far_period, 
                   MUT_close_period, MUT_mid_period, MUT_far_period, 
                   FRED_close_amplitude, FRED_mid_amplitude, FRED_far_amplitude, 
                   MUT_close_amplitude, MUT_mid_amplitude, MUT_far_amplitude)
csv_path = file.path(wd, paste0(tables_vec, ".csv"))
for (i in seq(1:12)){
  print(paste0("saving: ", tables_vec[i]))
  write.csv(tables_vec_obj[i], csv_path[i])
}
write.csv(FRED_close_period, csv_path[1])
write.csv(FRED_mid_period, csv_path[2])
write.csv(FRED_far_period, csv_path[3])
write.csv(MUT_close_period, csv_path[4])
write.csv(MUT_mid_period, csv_path[5])
write.csv(MUT_far_period, csv_path[6])