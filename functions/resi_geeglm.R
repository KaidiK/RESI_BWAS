

resi_pe_geeglm <- function(model.full, data, anova = TRUE,
                           coefficients = TRUE, unbiased = TRUE, ...){
  # data = model.full$model
  # sample size
  N = length(summary(model.full)$clusz)
  # total num of observations
  tot_obs = nrow(data)
  # num of observations from each subject
  n_i = model.full$id %>% table 
  n_i = rep(n_i, times = n_i)
  # weight in independent model
  data$w = 1 / n_i
  
  output <- list(model.full = list(call = model.full$call, formula = formula(model.full)),
                 estimates = c())
  names.est = c()
  
  # model form
  form = formula(model.full)
  # independence model
  # w_resi = model.full$prior.weights
  # data$w_resi = w_resi
  mod_indg = glm(formula = form, family = model.full$family, data = data, w = w,
                 contrasts = model.full$contrasts)
  mod_indg$residuals = mod_indg$residuals / sqrt(mod_indg$weights) # scale the residuals to get ready for `vcovHC`
  # mod_indg = update(model.full, data = data, id = 1:nrow(data))
  # the var-cov matrix estimate from the independence model
  # Note: this is the estimate for Cov[(\hat{\beta}_ind - \beta_0)]
  cov_ind = sandwich::vcovHC(mod_indg, type = "HC0")
  
  # make copy of model.full, replace vbeta with independence
  mod_ind = model.full
  mod_ind$geese$vbeta = cov_ind
  
  # The var-cov matrix estimate from the analysis model (the one considering correlation )
  # Note: this is the estimate for Cov(\hat{\beta}_long)
  # cov_long = vcov(model.full)
  # # convert it to the estimate for \Sigma_long = Cov(\sqrt{N}(\hat{\beta}_long - \beta_0))
  # cov_long = cov_long * N
  
  # longitudinal RESI
  # coefficients (z statistics)
  if (coefficients) {
    # longitudinal resi
    coefficients.tab <- lmtest::coeftest(model.full)
    coefficients.df = data.frame(coefficients.tab[,'Estimate'],
                                 coefficients.tab[,'Std. Error'],
                                 coefficients.tab[,'z value'],
                                 coefficients.tab[,'Pr(>|z|)'],
                                 row.names = rownames(coefficients.tab))
    colnames(coefficients.df) = colnames(coefficients.tab)
    if (unbiased){
      coefficients.df[,'L-RESI'] = z2S(coefficients.df[,'z value'], N)
    } else{
      coefficients.df[,'L-RESI'] = suppressWarnings(z2S_alt(coefficients.df[,'z value'],
                                                            N))
    }
    
    # CS RESI
    # M: use independence mod
    coefficients.tabcs <- lmtest::coeftest(mod_ind, vcov. = cov_ind)
    z_cs = coefficients.tabcs[,'z value']
    if (unbiased){
      coefficients.df[,'CS-RESI'] = z2S(z_cs, N)
    } else{
      coefficients.df[,'CS-RESI'] = suppressWarnings(z2S_alt(z_cs, N))
    }
    
    output$coefficients = coefficients.df
    output$estimates = c(coefficients.df$`L-RESI`, coefficients.df[,"CS-RESI"])
    names.est = rep(rownames(coefficients.df), 2)
    names(output$estimates) = names.est
  }
  
  # anova
  if (anova) {
    # longitudinal RESI
    anova.tab <- anova(model.full)
    anova.tab[,'L-RESI'] = chisq2S(anova.tab[,'X2'], anova.tab[,'Df'], N)
    
    # CS RESI
    # M: use anova() on new mod with substituted variance
    anova.tabcs <- suppressMessages(car::Anova(mod_indg, vcov. = cov_ind,
                                               test.statistic = "Wald"))
    anova.tab[,'CS-RESI'] = chisq2S(anova.tabcs[,'Chisq'], anova.tabcs[,'Df'], N)
    
    output$anova = anova.tab
    output$estimates = c(output$estimates, anova.tab$`L-RESI`,
                         anova.tab[,"CS-RESI"])
    names.est = c(names.est, rep(rownames(anova.tab), 2))
    names(output$estimates) = names.est
    class(output$anova) = c("anova_resi", class(output$anova))
  }
  
  output$naive.var = FALSE
  return(output)
}


boot.samp <- function(data, id.var = NULL) {
  params = as.list(match.call()[-1])
  if (is.matrix(data)) data = as.data.frame(data)
  if (is.null(id.var)) {
    boot.ind = sample(1:nrow(data), replace = TRUE)
    boot.data = data[boot.ind, ]
  } else {
    boot.ind = sample(unique(data[, id.var]), replace = TRUE)
    boot.data = data[unlist(lapply(boot.ind, function(x) which(x == data[, id.var]))), ]
    boot.data$bootid = rep(1:length(unique(data[, id.var])),
                           unlist(lapply(boot.ind, function(x) length(which(x==data[,id.var])))))
  }
  return(boot.data)
}

# model.full = fit_age
# data = data_i
# anova = TRUE
# coefficients = TRUE
# nboot = 10
# alpha = 0.05
# store.boot = FALSE
# unbiased = TRUE

resi_geeglm <- function(model.full, data, anova = TRUE,
                        coefficients = TRUE, nboot = 1000,
                        alpha = 0.05, store.boot = FALSE,
                        unbiased = TRUE, ...){
  # dots = list(...)
  
  
  if (missing(data)){
    data = model.full$data
    tryCatch(update(model.full, data = data), error = function(e){
      message("Updating model fit failed. Try rerunning with providing data argument")})
  } else{
    data = as.data.frame(data)
  }
  
  # point estimation
  output <- list(alpha = alpha, nboot = nboot, boot.method = "nonparam")
  output = c(output, resi_pe_geeglm(model.full = model.full, data = data, anova = anova,
                                    coefficients = coefficients, unbiased = unbiased))
  
  # id variable name
  id_var = as.character(model.full$call$id)
  corstr_spec = model.full$corstr
  # bootstrapping
  boot.results = data.frame(matrix(nrow = nboot, ncol = length(output$estimates)))
  colnames(boot.results) = names(output$estimates)
  fail = 0
  for (i in 1:nboot){
    skip_to_next <- FALSE
    boot.data = boot.samp(data, id.var = id_var)
    # re-fit the model
    boot.model.full = update(model.full, data = boot.data, corstr = corstr_spec,
                             id = bootid)
    rv.boot = tryCatch(resi_pe_geeglm(boot.model.full, data = boot.data, anova = anova,
                                      coefficients = coefficients, unbiased = unbiased),
                       error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {
      fail = fail + 1
      next}
    #output.boot$RESI= cbind(output.boot$RESI, rv.boot$resi[, 'RESI'])
    boot.results[i,] = suppressWarnings(resi_pe_geeglm(model.full = boot.model.full,
                                                       data = boot.data, anova = anova,
                                                       coefficients = coefficients,
                                                       unbiased = unbiased)$estimates)
  }
  
  alpha.order = sort(c(alpha/2, 1-alpha/2))
  
  if (coefficients){
    lCIs = apply(boot.results[,1:nrow(output$coefficients)], 2,  quantile,
                 probs = alpha.order, na.rm = TRUE)
    lCIs = t(lCIs)
    output$coefficients[1:nrow(lCIs), paste("L ",alpha.order*100, '%', sep='')] = lCIs
    cCIs = apply(boot.results[,(nrow(output$coefficients)+1):(2*nrow(output$coefficients))],
                 2,  quantile, probs = alpha.order, na.rm = TRUE)
    cCIs = t(cCIs)
    output$coefficients[1:nrow(cCIs), paste("CS ",alpha.order*100, '%', sep='')] = cCIs
  }
  
  if (anova){
    lCIs = apply(boot.results[,(ncol(boot.results)-
                                  2*length(rownames(output$anova))+1):
                                (ncol(boot.results)-length(rownames(output$anova)))],
                 2,  quantile, probs = alpha.order, na.rm = TRUE)
    lCIs = t(lCIs)
    output$anova[1:nrow(lCIs), paste("L ", alpha.order*100, '%', sep='')] = lCIs
    cCIs = apply(boot.results[,(ncol(boot.results)-length(rownames(output$anova))+1):
                                ncol(boot.results)],
                 2,  quantile, probs = alpha.order, na.rm = TRUE)
    cCIs = t(cCIs)
    output$anova[1:nrow(cCIs), paste("CS ", alpha.order*100, '%', sep='')] = cCIs
  }
  
  if(store.boot){
    output$boot.results = boot.results
  }
  
  output$nfail = fail
  
  class(output) = "resi"
  return(output)
  
  # RESI_se = apply(output.boot$RESI, 1, sd, na.rm = TRUE)
  
  return(output)
}