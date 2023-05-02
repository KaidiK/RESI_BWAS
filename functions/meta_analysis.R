meta_resi = function(ES, variable = "age", wls = TRUE, figure = TRUE){
  
  # ES$RESI_age = ifelse(ES$Design == "Cross-sectional", ES$CS_RESI_age, ES$L_RESI_age)
  # ES$RESI_sex = ifelse(ES$Design == "Cross-sectional", ES$CS_RESI_sex, ES$L_RESI_sex)
  # 
  if (figure) {
    # plots
    if (variable == "age"){
      # RESI for age
      # 1. CS-RESI
      ## S for age vs mean(age)
      p1 = ggplot(data = ES, aes(y = CS_RESI_age, x = mean_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "age", y = 'CS RESI', x = 'Mean age')
      ## S for age vs 2nd moment of age
      p2 = ggplot(data = ES, aes(y = CS_RESI_age, x = m2_age )) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = " age", y = 'CS RESI', x = 'SD of age')
      ## S for age vs 3rd moment of age
      p3 = ggplot(data = ES, aes(y = CS_RESI_age, x = m3_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "age", y = 'CS RESI', x = 'Skewness of age')
      
      # 2. Long RESI
      ## S for age vs mean age
      p4 = ggplot(data = ES, aes(y = L_RESI_age, x =  mean_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "age", y = 'Long RESI', x = 'Mean age')
      ## S for age vs 2nd moment of age
      p5 = ggplot(data = ES, aes(y = L_RESI_age, x =  m2_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "age", y = 'Long RESI', x = 'SD of age')
      ## S for age vs 3rd moment of age
      p6 = ggplot(data = ES, aes(y = L_RESI_age, x = m3_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "age", y = 'Long RESI', x = 'Skewness of age')
      plot_age = ggarrange(p1, p2, p3, p4, p5, p6)
      print(plot_age)      
    } else {
      # RESI for sex
      # 1. CS RESI
      ## S for sex vs mean age
      p1 = ggplot(data = ES, aes(y =  CS_RESI_sex, x = mean_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "sex", y = 'CS RESI', x = 'mean age')
      ## S for sex vs SD of age
      p2 = ggplot(data = ES, aes(y = CS_RESI_sex, x = m2_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "sex", y = 'CS RESI', x = 'SD of age')
      
      ## S for sex vs mean of sex
      p3 = ggplot(data = ES, aes(y = CS_RESI_sex, x = mean_sex)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "sex", y = 'CS RESI', x = 'proportion of male')
      
      
      # 2. Long RESI
      ## S for sex vs mean age
      p4 = ggplot(data = ES, aes(y = L_RESI_sex, x =  mean_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "sex", y = 'Long RESI', x = 'mean age')
      ## S for sex vs SD of age
      p5 = ggplot(data = ES, aes(y = L_RESI_sex, x = m2_age)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "sex", y = 'long RESI', x = 'SD of age')
      ## S for sex vs mean sex
      p6 = ggplot(data = ES, aes(y = L_RESI_sex, x = mean_sex)) + geom_point() + geom_smooth(formula = y ~ x, method = loess) + labs(title = "sex", y = 'Long RESI', x = 'proportion of male')
      
      plot_sex = ggarrange(p1, p2, p3, p4, p5, p6)
      print(plot_sex)    
    }
  }
  
  # CS RESI for longitudinal studies or RESI for CS studies
  # weights
  if (wls) {
    if (variable == "age") {
      ES$w = 1/ES$CS_RESI_se_age
    } else {
      ES$w = 1/ES$CS_RESI_se_sex
    }
  } else { ES$w = 1}
  
  if (variable == "age") {
    ES_fit_csRESI = lm(CS_RESI_age ~ Design + ns(mean_age, 3) + ns(m2_age, 3) + ns(m3_age, 3), data = ES, weights = w)
  } 
  if (variable == "sex") {
    ES_fit_csRESI = lm(CS_RESI_sex ~ Design + ns(mean_age, 3) + ns(m2_age, 3) + ns(mean_sex, 3), data = ES, weights = w)
  }
  
  # using robust (sandiwich) vcov estimator
  # ES_tab_pmRESI = lmtest::coeftest(ES_fit_pmRESI, vcov. = sandwich::vcovHC)[,] %>% as.data.frame
  ES_tab_csRESI = RESI::resi(ES_fit_csRESI, data = ES)
  
  # (Longitudinal/Regular) RESI for all studies
  # weights
  if (wls) {
    if (variable == "age") {
      ES$w = 1/ES$L_RESI_se_age
    } else {
      ES$w = 1/ES$L_RESI_se_sex
    }
  } else { ES$w = 1}
  
  if (variable == "age") {
    ES_fit_lRESI = lm(L_RESI_age ~ Design + ns(mean_age, 3) + ns(m2_age, 3) + ns(m3_age, 3), data = ES, weights = w)
  } 
  if (variable == "sex") {
    ES_fit_lRESI = lm(L_RESI_sex ~ Design + ns(mean_age, 3) + ns(m2_age, 3) + ns(mean_sex, 3) , data = ES, weights = w)
  }
  
  # using robust (sandiwich) vcov estimator
  # ES_tab_RESI = lmtest::coeftest(ES_fit_RESI, vcov. = sandwich::vcovHC)[,] %>% as.data.frame
  ES_tab_lRESI = RESI::resi(ES_fit_lRESI, data = ES)
  
  
  # cat("Model: RESI", formula, ifelse(wls, "with weight = 1/sd_RESI \n\n", "\n\n"))
  
  cat("Variable:", ifelse(variable == "age", "age", "sex"), "\n")
  cat("For CS-RESI: \n")
  cat("Model: ", deparse(formula(ES_fit_csRESI$terms)), "\n")
  ES_tab_csRESI$anova %>% round(2) %>% kbl(caption = "Meta-analysis results for CS-RESI") %>% kable_classic(full_width = FALSE, html_font = "Cambria") %>% print
  ES_tab_csRESI$coefficients %>% round(2) %>% kbl(caption = "Meta-analysis results for CS-RESI") %>% kable_classic(full_width = FALSE, html_font = "Cambria") %>% print
  # summary(ES_fit_pmRESI)
  
  cat("\nFor regular RESI: \n")
  cat("Model: ", deparse(formula(ES_fit_lRESI$terms)), "\n")
  ES_tab_lRESI$anova  %>% round(2) %>% kbl(caption = "Meta-analysis results for regular RESI") %>% kable_classic(full_width = FALSE, html_font = "Cambria") %>% print
  ES_tab_lRESI$coefficients  %>% round(2) %>% kbl(caption = "Meta-analysis results for regular RESI") %>% kable_classic(full_width = FALSE, html_font = "Cambria") %>% print
  # summary(ES_fit_RESI)
}