# FUNCTION for creating study summary table
sum_tab = function(data) {
  # independent data
  ind_data = data %>% group_by(study, participant) %>% slice_sample(n = 1)
  # Note: there may be some subjects from different study but having same subject IDs.
  tab1 = ind_data %>% group_by(study) %>% summarize(N = length(unique(participant)))
  tab2 = data %>% group_by(study) %>% summarise(tot_obs = length(participant))
  tab = merge(tab1, tab2)
  tab = tab[order(tab$N, decreasing = TRUE),]
  rownames(tab) = as.character(1:nrow(tab))
  tab %<>% mutate(Design = ifelse(N == tot_obs, "Cross-sectional", "Longitudinal"))
  return(tab)
}