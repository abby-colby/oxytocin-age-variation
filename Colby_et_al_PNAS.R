#########
## set up
#########

## packages
library(dplyr)
library(brms)
library(tidyverse)
library(tidyquant)
library(bayesplot)

## data
setwd("path/to/your/directory") # set working directory
data <- read.csv("age_ot_data.csv")

## correct for specific gravity
data$OT <- as.numeric(data$OT) # format OT
data$SG <- as.numeric(data$SG) # format SG
data$ot_sg_corrected <- (data$OT) * ((mean(data$SG) - 1.0) /
                                       ((data$SG) - 1.0)) ## control OT concentration in urine with specific gravity
sum(is.na(data$ot_sg_corrected)) # check for NAs

## format data
data$ot_sg_corrected <- as.numeric(data$ot_sg_corrected)
data$pid <- as.factor(data$pid)
data$fieldwork <- as.factor(data$fieldwork)
data$sex <- as.factor(data$sex)
data$context <- as.factor(data$context)

## center age at sample collection
data$age_centered <- scale(data$age, center = TRUE, scale = FALSE)

################################################################
## how does OT vary in females and males across the life course?
################################################################

## change reference group for the model
data$context <- relevel(factor(data$context), ref = "morning")

## age model 
age_model <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                  s(age_centered, by = sex) + sex + context + fieldwork + (1 | pid))
## run
age_model_fit <- brm(age_model,
                     data = data,
                     prior = c(
                       prior(normal(0, 2), class = "Intercept"),
                       prior(normal(0, 1), class = "b"),
                       prior(exponential(2), class = "sd")),
                     family = lognormal(),
                     chains = 4, cores = 4, warmup = 1000, iter = 4000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                     save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(age_model_fit, "age_model_pnas_2025.RDS")
## summary
summary(age_model_fit)
## conditional effects
conditional_effects(age_model_fit)
## plots
plot(age_model_fit)
## pp check
pp_check(age_model_fit) +
  coord_cartesian(xlim = c(0, 15))

## do females have higher OT than males? 
post <- as_draws_df(age_model_fit)
effect_male <- post$b_sexM
mean(effect_male)
quantile(effect_male, probs = c(0.025, 0.975))   
quantile(effect_male, probs = c(0.05, 0.95))              
mean(effect_male < 0) 

####################################################
## what drives age- and sex-related variation in OT?
####################################################

## is breastfeeding a driver of  age- and sex-related OT variation? 

## filter out NAs
data_breastfeeding <- data %>%
  filter(!is.na(breastfeeding))

## age_centered
data_breastfeeding$age_centered <- scale(data_breastfeeding$age, center = TRUE, scale = FALSE)

## change reference group for the model
data_breastfeeding$sex <- relevel(factor(data_breastfeeding$sex), ref = "F")
data_breastfeeding$context <- relevel(factor(data_breastfeeding$context), ref = "morning")
data_breastfeeding$fieldwork <- as.factor(data_breastfeeding$fieldwork)
data_breastfeeding$breastfeeding <- as.factor(data_breastfeeding$breastfeeding)
data_breastfeeding$pid <- as.factor(data_breastfeeding$pid)

## breastfeeding model
breastfeeding_model <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                            s(age_centered, by = sex) + sex + breastfeeding + context + fieldwork + (1 | pid))
## run
breastfeeding_model_fit <- brm(breastfeeding_model,
                               data = data_breastfeeding,
                               prior = c(
                                 prior(normal(0, 2), class = "Intercept"),
                                 prior(normal(0, 1), class = "b"),
                                 prior(exponential(2), class = "sd")),
                               family = lognormal(),
                               chains = 4, cores = 4, warmup = 1000, iter = 4000,
                               control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                               save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(breastfeeding_model_fit, "breastfeeding_model_pnas_2025.RDS")
## summary
summary(breastfeeding_model_fit)
## conditional effects
conditional_effects(breastfeeding_model_fit)
## plots
plot(breastfeeding_model_fit)
## pp check
pp_check(breastfeeding_model_fit) +
  coord_cartesian(xlim = c(0, 15))

## do breastfeeding women have higher OT than non-breastfeeding women?
post <- as_draws_df(breastfeeding_model_fit)
effect_breastfeeding <- post$b_breastfeeding1
mean(effect_breastfeeding)
quantile(effect_breastfeeding, probs = c(0.025, 0.975))   
quantile(effect_breastfeeding, probs = c(0.05, 0.95))              
mean(effect_breastfeeding > 0) 

## do females have higher OT than males after controlling for breastfeeding? 
post <- as_draws_df(breastfeeding_model_fit)
effect_male <- post$b_sexM
mean(effect_male)
quantile(effect_male, probs = c(0.025, 0.975))   
quantile(effect_male, probs = c(0.05, 0.95))              
mean(effect_male < 0) 

## is engagement in childcare childcare a driver of  age- and sex-related OT variation? 

## filter out NAs
data_childcare <- data %>%
  filter(!is.na(childcare) & !is.na(breastfeeding))

## age_centered
data_childcare$age_centered <- scale(data_childcare$age, center = TRUE, scale = FALSE)

## change reference group for the model
data_childcare$sex <- relevel(factor(data_childcare$sex), ref = "F")
data_childcare$childcare <- as.factor(data_childcare$childcare) 
data_childcare$fieldwork <- as.factor(data_childcare$fieldwork)
data_childcare$breastfeeding <- as.factor(data_childcare$breastfeeding)
data_childcare$pid <- as.factor(data_childcare$pid)

## childcare model
childcare_model <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                        s(age_centered, by = sex) + sex*childcare + breastfeeding + fieldwork + (1 | pid))
## run
childcare_model_fit <- brm(childcare_model,
                           data = data_childcare,
                           prior = c(
                             prior(normal(0, 2), class = "Intercept"),
                             prior(normal(0, 1), class = "b"),
                             prior(exponential(2), class = "sd")),
                           family = lognormal(),
                           chains = 4, cores = 4, warmup = 1000, iter = 4000,
                           control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                           save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(childcare_model_fit, "childcare_model_pnas_2025.RDS")
## summary
summary(childcare_model_fit)
## conditional effects
conditional_effects(childcare_model_fit)
## plots
plot(childcare_model_fit)
## pp check
pp_check(childcare_model_fit) +
  coord_cartesian(xlim = c(0, 15))

## do men that engaged in childcare have higher OT than men that did not engage in childcare?
post <- as_draws_df(childcare_model_fit)
effect_male <- post$b_childcare1 + post$`b_sexM:childcare1`
mean(effect_male)
quantile(effect_male, probs = c(0.025, 0.975))   
quantile(effect_male, probs = c(0.05, 0.95))              
mean(effect_male > 0)  

## do women that engaged in childcare have higher OT than women that did not engage in childcare?
post <- as_draws_df(childcare_model_fit)
effect_female <- post$b_childcare1 
mean(effect_female)
quantile(effect_female, probs = c(0.025, 0.975))   
quantile(effect_female, probs = c(0.05, 0.95))            
mean(effect_female > 0)  

## is being a parent a driver of age- and sex-related OT variation? 

## filter out NAs
data_parent <- data %>%
  filter(!is.na(breastfeeding) & !is.na(parent))

## age_centered
data_parent$age_centered <- scale(data_parent$age, center = TRUE, scale = FALSE)

## change reference group for the models
data_parent$sex <- relevel(factor(data_parent$sex), ref = "F")
data_parent$context <- relevel(factor(data_parent$context), ref = "morning")
data_parent$fieldwork <- as.factor(data_parent$fieldwork)
data_parent$breastfeeding <- as.factor(data_parent$breastfeeding)
data_parent$parent <- as.factor(data_parent$parent)
data_parent$pid <- as.factor(data_parent$pid)

## parent model
parent_model <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                     s(age_centered, by = sex) + sex*parent + breastfeeding + context + fieldwork + (1 | pid))
## run
parent_model_fit <- brm(parent_model,
                        data = data_parent,
                        prior = c(
                          prior(normal(0, 2), class = "Intercept"),
                          prior(normal(0, 1), class = "b"),
                          prior(exponential(2), class = "sd")),
                        family = lognormal(),
                        chains = 4, cores = 4, warmup = 1000, iter = 4000,
                        control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                        save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(parent_model_fit, "parent_model_pnas_2025.RDS")
## summary
summary(parent_model_fit)
## conditional effects
conditional_effects(parent_model_fit)
## plots
plot(parent_model_fit)
## pp check
pp_check(parent_model_fit) +
  coord_cartesian(xlim = c(0, 15))

## do mothers have higher OT than non-mothers? 
post <- as_draws_df(parent_model_fit)
effect_female <- post$b_parent1 
mean(effect_female)
quantile(effect_female, probs = c(0.025, 0.975))
quantile(effect_female, probs = c(0.05, 0.95))              
mean(effect_female > 0)  

## do fathers have higher OT than non-fathers? 
post <- as_draws_df(parent_model_fit)
effect_male <- post$b_parent1 + post$`b_sexM:parent1`
mean(effect_male)
quantile(effect_male, probs = c(0.025, 0.975))  
quantile(effect_male, probs = c(0.05, 0.95))            
mean(effect_male > 0)  

## is being a grandparent a driver of  age- and sex-related OT variation? 

## filter out NAs
data_grandparent <- data %>%
  filter(!is.na(grandparent) & !is.na(breastfeeding))

## age_centered
data_grandparent$age_centered <- scale(data_grandparent$age, center = TRUE, scale = FALSE)

## change reference group for the models
data_grandparent$sex <- relevel(factor(data_grandparent$sex), ref = "F")
data_grandparent$context <- relevel(factor(data_grandparent$context), ref = "morning")
data_grandparent$fieldwork <- as.factor(data_grandparent$fieldwork)
data_grandparent$grandparent <- as.factor(data_grandparent$grandparent) 
data_grandparent$breastfeeding <- as.factor(data_grandparent$breastfeeding) 
data_grandparent$pid <- as.factor(data_grandparent$pid)

## grandparent model
grandparent_model <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                          s(age_centered, by = sex) + sex*grandparent + breastfeeding + context + fieldwork + (1 | pid))
## run
grandparent_model_fit <- brm(grandparent_model,
                             data = data_grandparent,
                             prior = c(
                               prior(normal(0, 2), class = "Intercept"),
                               prior(normal(0, 1), class = "b"),
                               prior(exponential(2), class = "sd")),
                             family = lognormal(),
                             chains = 4, cores = 4, warmup = 1000, iter = 4000,
                             control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                             save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(grandparent_model_fit, "grandparent_model_pnas_2025.RDS")
## summary
summary(grandparent_model_fit)
## conditional effects
conditional_effects(grandparent_model_fit)
## plots
plot(grandparent_model_fit)
## pp check
pp_check(grandparent_model_fit) +
  coord_cartesian(xlim = c(0, 15))

## do grandmothers have higher OT than non-grandmothers? 
post <- as_draws_df(grandparent_model_fit)
effect_female <- post$b_grandparent1 
mean(effect_female)
quantile(effect_female, probs = c(0.025, 0.975))   
quantile(effect_female, probs = c(0.05, 0.95))              
mean(effect_female > 0) 

## do grandfathers have higher OT than non-grandfathers? 
post <- as_draws_df(grandparent_model_fit)
effect_male <- post$b_grandparent1 + post$`b_sexM:grandparent1`
mean(effect_male)
quantile(effect_male, probs = c(0.025, 0.975))   
quantile(effect_male, probs = c(0.05, 0.95))             
mean(effect_male > 0)  

## health

## filter out NAs
data_health <- data %>%
  filter(!is.na(health))

## age_centered
data_health$age_centered <- scale(data_health$age, center = TRUE, scale = FALSE)

## change reference group for the models
data_health$sex <- relevel(factor(data_health$sex), ref = "F")
data_health$context <- relevel(factor(data_health$context), ref = "morning")
data_health$health <- relevel(factor(data_health$health), ref = "Normal")
data_health$pid <- as.factor(data_health$pid)

## health model
health_model <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                     s(age_centered, by = sex) + sex*health + context + (1 | pid))
## run
health_model_fit <- brm(health_model,
                        data = data_health,
                        prior = c(
                          prior(normal(0, 2), class = "Intercept"),
                          prior(normal(0, 1), class = "b"),
                          prior(exponential(2), class = "sd")),
                        family = lognormal(),
                        chains = 4, cores = 4, warmup = 1000, iter = 4000,
                        control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                        save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(health_model_fit, "health_model_pnas_2025.RDS")
## summary
summary(health_model_fit)
## conditional effects
conditional_effects(health_model_fit_bf)
## plots
plot(health_model_fit)
## pp check
pp_check(health_model_fit) +
  coord_cartesian(xlim = c(0, 15))

## do women with good health have higher OT than women with normal health? 
post <- as_draws_df(health_model_fit)
effect_female <- post$b_healthGood
mean(effect_female)                                   
quantile(effect_female, probs = c(0.025, 0.975))   
quantile(effect_female, probs = c(0.05, 0.95))                
mean(effect_female > 0)                              

## do women with good health have higher OT than women with bad health?
post <- as_draws_df(health_model_fit)
effect_female <- post$b_healthGood - post$b_healthBad
mean(effect_female)                                  
quantile(effect_female, probs = c(0.025, 0.975))     
quantile(effect_female, probs = c(0.05, 0.95))                
mean(effect_female > 0)                              

## do men with good health have higher OT than men with normal health?
post <- as_draws_df(health_model_fit)
effect_male <- post$b_healthGood + post$`b_sexM:healthGood`
mean(effect_male)                                
quantile(effect_male, probs = c(0.025, 0.975))     
quantile(effect_male, probs = c(0.05, 0.95))              
mean(effect_male > 0)

## do men with good health have higher OT than men with bad health?
post <- as_draws_df(health_model_fit)
effect_male <- (post$b_healthGood + post$`b_sexM:healthGood`) -
  (post$b_healthBad + post$`b_sexM:healthBad`)
mean(effect_male)                                   
quantile(effect_male, probs = c(0.025, 0.975))   
quantile(effect_male, probs = c(0.05, 0.95))               
mean(effect_male > 0)    

########################
## supplemental analyses
########################

## is being a maternal grandparent a driver of  age- and sex-related OT variation? 

## filter out NAs
data_grandparent_maternal <- data %>%
  filter(!is.na(maternal_grandparent) & !is.na(breastfeeding))

## age_centered
data_grandparent_maternal$age_centered <- scale(data_grandparent_maternal$age, center = TRUE, scale = FALSE)

## change reference group for the models / format
data_grandparent_maternal$sex <- relevel(factor(data_grandparent_maternal$sex), ref = "F")
data_grandparent_maternal$context <- relevel(factor(data_grandparent_maternal$context), ref = "morning")
data_grandparent_maternal$maternal_grandparent <- as.factor(data_grandparent_maternal$maternal_grandparent) 
data_grandparent_maternal$breastfeeding <- as.factor(data_grandparent_maternal$breastfeeding) 
data_grandparent_maternal$fieldwork <- as.factor(data_grandparent_maternal$fieldwork)
data_grandparent_maternal$pid <- as.factor(data_grandparent_maternal$pid)

## maternal grandparent model
grandparent_maternal_model <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                                   s(age_centered, by = sex) + sex*maternal_grandparent + breastfeeding + context + fieldwork + (1 | pid))
## run
grandparent_maternal_model_fit <- brm(grandparent_maternal_model,
                                      data = data_grandparent_maternal,
                                      prior = c(
                                        prior(normal(0, 2), class = "Intercept"),
                                        prior(normal(0, 1), class = "b"),
                                        prior(exponential(2), class = "sd")),
                                      family = lognormal(),
                                      chains = 4, cores = 4, warmup = 1000, iter = 4000,
                                      control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                                      save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(grandparent_maternal_model_fit, "grandparent_maternal_model_pnas_2025.RDS")
## summary
summary(grandparent_maternal_model_fit)
## conditional effects
conditional_effects(grandparent_maternal_model_fit)
## plots
plot(grandparent_maternal_model_fit)
## pp check
pp_check(grandparent_maternal_model_fit) +
  coord_cartesian(xlim = c(0, 15))

## do maternal grandmothers have higher OT than non-maternal grandmothers? 
post <- as_draws_df(grandparent_maternal_model_fit)
effect_female <- post$b_maternal_grandparent1 
mean(effect_female)
quantile(effect_female, probs = c(0.025, 0.975))   
quantile(effect_female, probs = c(0.05, 0.95))              
mean(effect_female > 0) 

## do maternal grandfathers have higher OT than non-maternal grandfathers? 
post <- as_draws_df(grandparent_maternal_model_fit)
effect_male <- post$b_maternal_grandparent1 + post$`b_sexM:maternal_grandparent1`
mean(effect_male)
quantile(effect_male, probs = c(0.025, 0.975))   
quantile(effect_male, probs = c(0.05, 0.95))             
mean(effect_male > 0)  

## is being a paternal grandparent a driver of age- and sex-related OT variation? 

## filter out NAs
data_grandparent_paternal <- data %>%
  filter(!is.na(paternal_grandparent) & !is.na(breastfeeding))

## age_centered
data_grandparent$age_centered <- scale(data_grandparent$age, center = TRUE, scale = FALSE)

## change reference group for the models / format
data_grandparent_paternal$sex <- relevel(factor(data_grandparent_paternal$sex), ref = "F")
data_grandparent_paternal$context <- relevel(factor(data_grandparent_paternal$context), ref = "morning")
data_grandparent_paternal$paternal_grandparent <- as.factor(data_grandparent_paternal$paternal_grandparent) 
data_grandparent_paternal$breastfeeding <- as.factor(data_grandparent_paternal$breastfeeding) 
data_grandparent_paternal$fieldwork <- as.factor(data_grandparent_paternal$fieldwork)
data_grandparent_paternal$pid <- as.factor(data_grandparent_paternal$pid)

## paternal grandparent model
grandparent_paternal_model <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                                   s(age_centered, by = sex) + sex*paternal_grandparent + breastfeeding + context + fieldwork + (1 | pid))
## run
grandparent_paternal_model_fit <- brm(grandparent_paternal_model,
                                      data = data_grandparent_paternal,
                                      prior = c(
                                        prior(normal(0, 2), class = "Intercept"),
                                        prior(normal(0, 1), class = "b"),
                                        prior(exponential(2), class = "sd")),
                                      family = lognormal(),
                                      chains = 4, cores = 4, warmup = 1000, iter = 4000,
                                      control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                                      save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(grandparent_paternal_model_fit, "grandparent_paternal_model_pnas_2025.RDS")
## summary
summary(grandparent_paternal_model_fit)
## conditional effects
conditional_effects(grandparent_paternal_model_fit)
## plots
plot(grandparent_paternal_model_fit)
## pp check
pp_check(grandparent_paternal_model_fit) +
  coord_cartesian(xlim = c(0, 15))

## do paternal grandmothers have higher OT than non-paternal grandmothers? 
post <- as_draws_df(grandparent_paternal_model_fit)
effect_female <- post$b_paternal_grandparent1 
mean(effect_female)
quantile(effect_female, probs = c(0.025, 0.975))   
quantile(effect_female, probs = c(0.05, 0.95))              
mean(effect_female > 0) 

## do maternal grandfathers have higher OT than non-maternal grandfathers? 
post <- as_draws_df(grandparent_paternal_model_fit)
effect_male <- post$b_paternal_grandparent1 + post$`b_sexM:paternal_grandparent1`
mean(effect_male)
quantile(effect_male, probs = c(0.025, 0.975))   
quantile(effect_male, probs = c(0.05, 0.95))             
mean(effect_male > 0)  

## models for model comparison with age model
## is the age model with spline the best fit? 

## null model (no age term)
age_model_null <- bf(ot_sg_corrected | cens(censored) ~ 1 + sex + context + fieldwork + (1 | pid))
## run
age_model_null_fit <- brm(age_model_null,
                          data = data,
                          prior = c(
                            prior(normal(0, 2), class = "Intercept"),
                            prior(normal(0, 1), class = "b"),
                            prior(exponential(2), class = "sd")),
                          family = lognormal(),
                          chains = 4, cores = 4, warmup = 1000, iter = 4000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                          save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(age_model_null_fit, "age_model_null_pnas_2025.RDS")
## summary
summary(age_model_null_fit)
## conditional effects
conditional_effects(age_model_null_fit)
## plots
plot(age_model_null_fit)
## pp check
pp_check(age_model_null_fit)

## linear model (no spline function for age)
age_model_linear <- bf(ot_sg_corrected | cens(censored) ~ 1 +
                         age_centered*sex + context + fieldwork + (1 | pid))
## run
age_model_linear_fit <- brm(age_model_linear,
                            data = data,
                            prior = c(
                              prior(normal(0, 2), class = "Intercept"),
                              prior(normal(0, 1), class = "b"),
                              prior(exponential(2), class = "sd")),
                            family = lognormal(),
                            chains = 4, cores = 4, warmup = 1000, iter = 4000,
                            control = list(adapt_delta = 0.9999, max_treedepth = 15), 
                            save_pars = save_pars(all = TRUE)) 

## save RDS
saveRDS(age_model_linear_fit, "age_model_linear_pnas_2025.RDS")
## summary
summary(age_model_linear_fit)
## conditional effects
conditional_effects(age_model_linear_fit)
## plots
plot(age_model_linear_fit)
## pp check
pp_check(age_model_linear_fit)

## model comparisons
## is the age model with spline the best fit?

options(future.globals.maxSize = 2 * 1024^3) 

## model comparisons - loo
am.0 <- add_criterion(age_model_null_fit, criterion = "loo", moment_match = TRUE, reloo = TRUE) 
am.1 <- add_criterion(age_model_linear_fit, criterion = "loo", moment_match = TRUE, reloo = TRUE)
am.2 <- add_criterion(age_model_fit, criterion = "loo", moment_match = TRUE, reloo = TRUE) 
loo_compare(am.0, am.1, am.2) 
## save loo
saveRDS(am.0, "am.0.pnas.2025.RDS") 
saveRDS(am.1, "am.1.pnas.2025.RDS") 
saveRDS(am.2, "am.2.pnas.2025.RDS") 

## model comparison - kfold
kam.0 <- add_criterion(age_model_null_fit, criterion = "kfold")
kam.1 <- add_criterion(age_model_linear_fit, criterion = "kfold")
kam.2 <- add_criterion(age_model_fit, criterion = "kfold")
loo_compare(kam.0$criteria$kfold,kam.1$criteria$kfold,kam.2$criteria$kfold)  
## save kfold
saveRDS(kam.0, "kam.0.pnas.2025.RDS") 
saveRDS(kam.1, "kam.1.pnas.2025.RDS") 
saveRDS(kam.2, "kam.2.pnas.2025.RDS") 
