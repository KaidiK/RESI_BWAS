forest_plot = function(study_info, outcome, x_ul = 3){
  # the outcome label
  y_label = paste0(outcome, "_10000_combat")
  
  # prepare plotting data
  ES = study_info
  ES$CS_RESI_age = NA
  ES$ll_CS_RESI_age = NA
  ES$ul_CS_RESI_age = NA
  ES$L_RESI_age = NA
  ES$ll_L_RESI_age = NA
  ES$ul_L_RESI_age = NA
  
  ES$CS_RESI_sex = NA
  ES$ll_CS_RESI_sex = NA
  ES$ul_CS_RESI_sex = NA
  ES$L_RESI_sex = NA
  ES$ll_L_RESI_sex = NA
  ES$ul_L_RESI_sex = NA
  
  ES$CS_RESI_se_age = NA
  ES$CS_RESI_se_sex = NA
  
  ES$L_RESI_se_age = NA
  ES$L_RESI_se_sex = NA
  
  ES$d_sex = NA
  ES$ll_d_sex = NA
  ES$ul_d_sex = NA
  
  for (i in 1:nrow(ES)){
    study_i = ES$study[i] %>% as.character()
    study_name = gsub("-", "_", study_i)
    fit_i = readRDS(paste0("RESI_objects/age_sex/", outcome, "/", outcome, "_RESI_", study_name, ".rds"))
    if (ES$Design[i] == "Longitudinal") { 
      # age
      # CS RESI
      ES$CS_RESI_age[i] = fit_i$resi_obj_age["ns(age_years, df = 2)", "CS-RESI"]
      ES$ll_CS_RESI_age[i] = fit_i$resi_obj_age["ns(age_years, df = 2)", "CS 2.5%"]
      ES$ul_CS_RESI_age[i] = fit_i$resi_obj_age["ns(age_years, df = 2)", "CS 97.5%"]
      ES$CS_RESI_se_age[i] = fit_i$resi_obj_age["ns(age_years, df = 2)", "CS-RESI_se"]
      
      # Long RESI
      ES$L_RESI_age[i] = fit_i$resi_obj_age["ns(age_years, df = 2)", "L-RESI"]
      ES$ll_L_RESI_age[i] = fit_i$resi_obj_age["ns(age_years, df = 2)", "L 2.5%"]
      ES$ul_L_RESI_age[i] = fit_i$resi_obj_age["ns(age_years, df = 2)", "L 97.5%"]  
      ES$L_RESI_se_age[i] = fit_i$resi_obj_age["ns(age_years, df = 2)", "L-RESI_se"]
      
      # sex
      # CS RESI
      ES$CS_RESI_sex[i] = fit_i$resi_obj_sex["sex_01", "CS-RESI"]
      ES$ll_CS_RESI_sex[i] = fit_i$resi_obj_sex["sex_01", "CS 2.5%"]
      ES$ul_CS_RESI_sex[i] = fit_i$resi_obj_sex["sex_01", "CS 97.5%"]
      ES$CS_RESI_se_sex[i] = fit_i$resi_obj_sex["sex_01", "CS-RESI_se"]
      # Long RESI
      ES$L_RESI_sex[i] = fit_i$resi_obj_sex["sex_01", "L-RESI"]
      ES$ll_L_RESI_sex[i] = fit_i$resi_obj_sex["sex_01", "L 2.5%"]
      ES$ul_L_RESI_sex[i] = fit_i$resi_obj_sex["sex_01", "L 97.5%"]
      ES$L_RESI_se_sex[i] = fit_i$resi_obj_sex["sex_01", "L-RESI_se"]
      
      # Cohen's d for sex effect
      d = z_to_d(sqrt(fit_i$resi_obj_sex["sex_01", "X2"]), n = ES$N[i])
      ES$d_sex[i] = d[1, "d"]
      ES$ll_d_sex[i] = d[1, "CI_low"]
      ES$ul_d_sex[i] = d[1, "CI_high"]
      
    } else {
      # Cross-sectional studies
      # Here, longitudinal RESI = CS RESI
      # age
      ES$CS_RESI_age[i] = ES$L_RESI_age[i] = fit_i$anova["ns(age_years, df = 2)", "RESI"]
      ES$ll_CS_RESI_age[i] = ES$ll_L_RESI_age[i] = fit_i$anova["ns(age_years, df = 2)", "2.5%"]
      ES$ul_CS_RESI_age[i] = ES$ul_L_RESI_age[i] = fit_i$anova["ns(age_years, df = 2)", "97.5%"]
      # sex
      ES$CS_RESI_sex[i] = ES$L_RESI_sex[i] = fit_i$anova["sex_01", "RESI"]
      ES$ll_CS_RESI_sex[i] = ES$ll_L_RESI_sex[i] = fit_i$anova["sex_01", "2.5%"]
      ES$ul_CS_RESI_sex[i] = ES$ul_L_RESI_sex[i] = fit_i$anova["sex_01", "97.5%"]
      # se of RESI
      ES$CS_RESI_se_age[i] = ES$L_RESI_se_age[i] = fit_i$anova["ns(age_years, df = 2)", "RESI_se"]
      ES$CS_RESI_se_sex[i] = ES$L_RESI_se_sex[i] = fit_i$anova["sex_01", "RESI_se"]
      
      # Cohen's d for sex effect
      d = F_to_d(fit_i$anova["sex_01", "F"], df = fit_i$anova["sex_01", "Df"], df_error = fit_i$overall$Res.Df[2])
      ES$d_sex[i] = d[1, "d"]
      ES$ll_d_sex[i] = d[1, "CI_low"]
      ES$ul_d_sex[i] = d[1, "CI_high"]
      
    }
  }
  
  # ES = ES[order(ES$index),]
  # ES$new_index = 1:nrow(ES)
  
  plot_data = subset(analysis_data, study %in% ES$study)
  plot_data = merge(plot_data, ES[, c("study", "Design")])
  plot_data$study %<>% factor(levels = ES$study)
  
  age_dist = ggplot(plot_data, aes(x = age_years, y = study, col = Design)) +
    geom_boxplot(aes(col = Design)) + 
    scale_colour_manual(values = c("#00AFBB", "#E7B800")) +
    guides(x =  guide_axis(angle = 90)) +
    labs(title = "Age distributions", x = "Age (in years)", y = "Study")
  
  # forest plot
  # 1. for age
  title_text = outcome
  ## 1.1 plotting for Long RESI (in yellow)
  # Check whether the actual limits go beyond the specified range
  # if so, change them to the specified limit
  ES <- ES %>% mutate(corrected_ul = ifelse(ul_L_RESI_age > x_ul, x_ul, NA))
  temp_data = ES
  temp_data$L_RESI_age[temp_data$Design == "Cross-sectional"] = NA 
  temp_data$ll_L_RESI_age[temp_data$Design == "Cross-sectional"] = NA
  temp_data$ul_L_RESI_age[temp_data$Design == "Cross-sectional"] = NA
  ES_plot_age = ggplot(data = temp_data, aes(y = index_sd_age - 0.6, x = L_RESI_age, xmin = ll_L_RESI_age, xmax = ul_L_RESI_age)) +
    geom_point(colour =  "#E7B800") + geom_errorbarh(aes(y =index_sd_age - 0.6), height=.1, colour =  "#E7B800") + xlim(0, x_ul)
  # adding arrows if needed
  if (sum(!is.na(ES$corrected_ul)) > 0){
    ES_plot_age = ES_plot_age + 
      geom_errorbarh(aes(y = index_sd_age - 0.6, xmin = ll_L_RESI_age, xmax = corrected_ul), height=.1, colour = "#E7B800") + 
      # arrows
      geom_segment(data = temp_data, 
                   aes(x = corrected_ul, xend = L_RESI_age, 
                       y = index_sd_age - 0.6, yend = index_sd_age - 0.6), 
                   color =  "#E7B800",
                   arrow = arrow(length=unit(0.30,"cm"), ends="first"))
  }
  ## 1.2 plotting for CS RESI (in blue)
  ES <- ES %>% mutate(corrected_ul = ifelse(ul_CS_RESI_age > x_ul, x_ul, NA))
  
  ES_plot_age = ES_plot_age +
    geom_point(data = ES, aes(y = index_sd_age - 0.4, x = CS_RESI_age), colour = "#00AFBB") +
    geom_errorbarh(height=.1, aes(y = index_sd_age - 0.4, xmin = ll_CS_RESI_age, xmax = ul_CS_RESI_age), colour = "#00AFBB")
  
  # adding arrows if needed
  if (sum(!is.na(ES$corrected_ul)) > 0){
    ES_plot_age = ES_plot_age + 
      # geom_errorbarh(aes(xmin = ll_age_adj, xmax = corrected_ul), height=.1) +
      # arrows
      geom_segment(data = ES, 
                   aes(x = corrected_ul, xend = ll_CS_RESI_age, 
                       y = index_sd_age - 0.4, yend = index_sd_age - 0.4), 
                   colour = "#00AFBB",
                   arrow = arrow(length=unit(0.30,"cm"), ends="first"))
  }
  
  # adding other elements
  ES_plot_age = ES_plot_age +
    scale_y_continuous(breaks = (1:nrow(ES))-0.5, expand = c(0, 0), 
                       labels = ES$study,
                       limits = c(0, nrow(ES))) +
    labs(title = title_text, x = 'RESI', y = 'Study') +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    theme_classic()
  
  # 2. for Gender
  ## 2.1 plotting for Long RESI (in yellow)
  # Check whether the actual limits go beyond the specified range
  # if so, change them to the specified limit
  ES <- ES %>% mutate(corrected_ul = ifelse(ul_L_RESI_sex > x_ul, x_ul, NA))
  temp_data = ES
  temp_data$L_RESI_sex[temp_data$Design == "Cross-sectional"] = NA 
  temp_data$ll_L_RESI_sex[temp_data$Design == "Cross-sectional"] = NA
  temp_data$ul_L_RESI_sex[temp_data$Design == "Cross-sectional"] = NA
  ES_plot_sex = ggplot(data = temp_data, aes(y = index_sd_age - 0.6, x = L_RESI_sex, xmin = ll_L_RESI_sex, xmax = ul_L_RESI_sex)) +
    geom_point(colour = "#E7B800") + geom_errorbarh(height=.1, colour = "#E7B800") + xlim(0, x_ul)
  # adding arrows if needed
  if (sum(!is.na(ES$corrected_ul)) > 0){
    ES_plot_sex = ES_plot_sex + 
      geom_errorbarh(aes(xmin = ll_L_RESI_sex, xmax = corrected_ul), height=.1, colour = "#E7B800") + 
      # arrows
      geom_segment(data = temp_data, 
                   aes(x = corrected_ul, xend = L_RESI_sex, 
                       y = index_sd_age - 0.6, yend = index_sd_age - 0.6), 
                   color = "#E7B800",
                   arrow = arrow(length=unit(0.30,"cm"), ends="first"))
  }
  ## 2.2 plotting for CS RESI (in black)
  ES <- ES %>% mutate(corrected_ul = ifelse(ul_CS_RESI_sex > x_ul, x_ul, NA))
  
  ES_plot_sex = ES_plot_sex +
    geom_point(data = ES, aes(y = index_sd_age - 0.4, x = CS_RESI_sex), colour = "#00AFBB") +
    geom_errorbarh(height=.1, aes(y = index_sd_age - 0.4, xmin = ll_CS_RESI_sex, xmax = ul_CS_RESI_sex), colour = "#00AFBB")
  
  # adding arrows if needed
  if (sum(!is.na(ES$corrected_ul)) > 0){
    ES_plot_sex = ES_plot_sex + 
      # arrows
      geom_segment(data = ES, 
                   aes(x = corrected_ul, xend = ll_CS_RESI_sex, 
                       y = index_sd_age - 0.4, yend = index_sd_age - 0.4), 
                   arrow = arrow(length=unit(0.30,"cm"), ends="first"), 
                   colour = "#00AFBB")
  }
  
  # adding other elements 
  ES_plot_sex = ES_plot_sex +
    scale_y_continuous(breaks = (1:nrow(ES))-0.5, expand = c(0, 0), 
                       labels = ES$study,
                       limits = c(0, nrow(ES))) +
    labs(title = title_text, x = 'RESI', y = 'Study') +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    theme_classic()
  
  ## 2.3 Plotting Cohen's d
  # Cohen's d for sex
  # Check whether the actual limits go beyond the specified range
  # if so, change them to the specified limit
  ES <- ES %>% mutate(corrected_ul = ifelse(ul_d_sex > x_ul, x_ul, NA))
  plot_d_sex = ggplot(data = ES, aes(y = index_sd_age - 0.5, x = d_sex, xmin = ll_d_sex, xmax = ul_d_sex)) +
    geom_point() + geom_errorbarh(height=.1) + xlim(-0.5, x_ul)
  # adding arrows if needed
  if (sum(!is.na(ES$corrected_ul)) > 0){
    plot_d_sex = plot_d_sex + 
      geom_errorbarh(aes(xmin = ll_d_sex, xmax = corrected_ul), height=.1) + 
      # arrows
      geom_segment(data = ES, 
                   aes(x = corrected_ul, xend = ll_d_sex, 
                       y = index_sd_age - 0.5, yend = index_sd_age - 0.5), 
                   arrow = arrow(length=unit(0.30,"cm"), ends="first"))
  }
  
  plot_d_sex  = plot_d_sex + 
    scale_y_continuous(breaks = (1:nrow(ES))-0.5, expand = c(0, 0), 
                       labels = ES$study,
                       limits = c(0, nrow(ES))) +
    labs(title = title_text, x = "Cohen's d", y = 'Study') +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    theme_classic()
  
  
  return(list(ES = ES,
              age_dist = age_dist,
              plot_age = ES_plot_age,
              plot_sex = ES_plot_sex,
              plot_d_sex = plot_d_sex)
  )
}