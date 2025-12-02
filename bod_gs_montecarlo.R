
rm(list = ls()) ### Clear environment

#-------------------------
# Load packages
#-------------------------

library(tidyverse)
library(dplyr)
library(sf)
library(ggplot2)

#-------------------------
# Load data
#-------------------------

### Population data
sofia_pop_poly <- st_read("Baseline_BOD/Clean_data/sofia_pop_poly.geojson")
sofia_pop_poly$geometry <- NULL
#sofia_pop_ct$pop_density <- NULL
names(sofia_pop_poly)[24] <- "male_70+"
names(sofia_pop_poly)[39] <- "female_70+"

# add columns with adult pop counts for disease analyses 
sofia_pop_poly$total_adults <- rowSums(sofia_pop_poly[15:24]) + rowSums(sofia_pop_poly[30:39])

### Health data 
sofia_mort <- read_csv2("Baseline_BOD/Health/Sofia_mortality.csv")
sofia_disease <- read_csv2("Baseline_BOD/Health/Sofia_diseases.csv")
sofia_mort$rate <- as.numeric(sofia_mort$count/ sofia_mort$population * 100000)
sofia_disease$rate <- as.numeric(sofia_disease$count/ sofia_disease$population * 100000)
sofia_health <- rbind(sofia_mort, sofia_disease)

### Green space data
sofia_gs_poly <- st_read("Baseline_BOD/Clean_data/sofia_green_poly.geojson")
sofia_gs_poly$geometry <- NULL
#IQR(sofia_gs_ct$ndvi) # 0.109001


#--------------------------------------------------
# First, we need to create some vectors and data 
# frames with parameters that will be used later 
# in the Monte Carlo analysis
#--------------------------------------------------

### We extract the age and sex labels
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # these need to be adjusted for each analysis based on the outcome
allages

### We create a data frame with the RRs for each pollutant and health outcome

rranalysis <- data.frame(pollutant = rep("ndvi", each = 4),
                         outcome = c("nat_mort", "cvd_mort", "hypertension", "stroke"), 
                         rr = c(0.96, 0.97, 0.98, 0.98),
                         rrlo95 = c(0.94, 0.96, 0.97, 0.96),
                         rrup95 = c(0.97, 0.99, 0.99, 0.99), 
                         scale = c(0.1, 0.1, 0.1, 0.1))


### Next, we define three functions which are needed to run the Monte Carlo simulations


#######################################
#######################################
###
### function "hia_age"
### to get one simulated value of health outcomes
### for a given age group and sex
###
#######################################
#######################################
hia_age <- function(age = allages[5],
                    pollutant = "ndvi", 
                    outcome = "nat_mort") {                 
  
  ### Uncertainty in RRs:
  # We retrieve the RR and its CIs based on the pollutant label:
  rrci <- as.numeric(rranalysis[rranalysis$pollutant == pollutant & rranalysis$outcome == outcome, c("rr", "rrlo95", "rrup95")])
  names(rrci) <- c("rr", "rrlo95", "rrup95")
  scale <- as.numeric(rranalysis[rranalysis$pollutant == pollutant & rranalysis$outcome == outcome, "scale"])
  
  # Then, we simulate a RR value assuming log(RR) is normal 
  selogrr <- log(rrci["rrup95"] / rrci["rrlo95"]) / (2 * qnorm(0.975)) # we calculate the standard deviation of log(RR)
  logrr <- log(rrci["rr"]) # we calculate log(RR)
  logrrsim <- rnorm(n = 1, mean = logrr, sd = selogrr) # we simulate one log(RR) value from a normal distribution
  rr <- exp(logrrsim) # We exponentiate back the simulated log(RR) --> the simulated RR will be then used in the HIA calculation

  # We establish the counterfactual exposure level:
  tmrel <- 0.4097499
  
  ### We simulate one error value for the counterfactual exposure level:
  # We calculate the standard error of model value (based on 95% CIs):
  sd <- (0.4136442 - 0.4058556) / (2 * qnorm(0.975))
  # We simulate the counterfactual level (assuming a normal distribution):
  tmrel_sim <- rnorm(n = 1, mean = tmrel, sd = sd)
    
  ### We retrieve the outcome rate based on the selected outcome and age group
  agesub <- sub(".*?_(.*)", "\\1", age)
  sexsub <- sub("\\_.*", "", age)
  mratecityage <- as.numeric(sofia_health[sofia_health$age == agesub & sofia_health$sex == sexsub & sofia_health$outcome == outcome, "rate"])

  ### We retrieve the city population for each neighborhood:
  popcity <- sofia_pop_poly[, c("Poly_ID", age)]
  names(popcity)[2] <- "pop"
  
  # Then, we compute the age-specific health outcome
  popcity$outcome <- popcity$pop * mratecityage / 100000
  
  ### We retrieve the pollutant values for all neighborhoods
  pollgrid <-  sofia_gs_poly[, c("Poly_ID", pollutant)]
  
  ### Next, we merge the air pollution and mortality data
  popcity <- merge(popcity, pollgrid, by = "Poly_ID")
  rm(pollgrid)
  
  ### Now that all values are simulated, we do the HIA: calculate exposure difference, scaled RR and PAF
  popcity$expdiff <- popcity$ndvi - tmrel_sim
  popcity$expdiff[popcity$expdiff > 0] <- 0  # zero if ndvi >= tmrel
  
  popcity$rr <- exp((log(rr) / scale) * popcity$expdiff)
  
  popcity$PAF <- with(popcity, (rr - 1) / rr)

  ### We calculate health outcome due to air pollution exposure in each neighborhood
  popcity$ndvi_outcome <- with(popcity, outcome * PAF)
  popcity <- popcity[c("Poly_ID", "ndvi_outcome")]
  #names(popcity)[1] <- "ct_id"
  
  return(popcity) ### We ask the function to return the health outcome estimate by neighborhood
}

### Example --> you will see that if you run this several times you get a different value each time (because of the uncertainties)
hia1 <- hia_age(age = allages[10], pollutant = "ndvi", outcome = "nat_mort")
hia1
sum(hia1$ndvi_outcome)



#######################################
#######################################
###
### function "hia_age_sim"
### to get "nsim" simulated values from the function
### "hia_age" and then
### provide both point estimate (median and mean) and
### empirical 95% CI (by sampling percentile)
### the vector of simulated values is also stored
### because it will be needed when adding the impact
### for all age groups in a given city --> this function provides point estimates + 95% CIs by age group
###
#######################################
#######################################

hia_age_sim <- function(age = allages[5],
                        pollutant = "ndvi", 
                        outcome = "nat_mort",
                        nsim = 100) {                  # number of simulations
  
  # We use the replication function to replicate the previously defined function n times (to get a sample of values):
  samp <- replicate(n = nsim,
                    expr = hia_age(age = age,
                                   pollutant = pollutant, 
                                   outcome = outcome), simplify = FALSE)
                    
  # We use the sample of values to get point (mean) and CI estimate (2.5 and 97.5 percentiles)
  summ <- bind_rows(samp, .id = "id") %>% 
    group_by(Poly_ID) %>%
    summarise(mean = mean(ndvi_outcome, na.rm = TRUE), 
              median = median(ndvi_outcome, na.rm = TRUE), 
              pct2.5 = quantile(ndvi_outcome, probs = 2.5/100, na.rm = TRUE), 
              pct97.5 = quantile(ndvi_outcome, probs = 97.5/100, na.rm = TRUE))
  
  # We ask the function to return a list with results
  res <- list(sample = samp, estimate = summ)
  return(res)
}

# example with 500 simulations --> you can try different numbers of simulations and check how long it takes to calculate
t0 <- Sys.time()
simhia1 <- hia_age_sim(age = allages[10],
                       pollutant = "ndvi", 
                       outcome = "nat_mort",
                       nsim = 500)
t1 <- Sys.time()
t1 - t0  # --> this is used to check the calculation time
names(simhia1)
simhia1$estimate # --> mean and median should be similar indicating that you have a normal distribution in your sample



#######################################
#######################################
###
### function "hia_sim"
### to get both point estimate (median and mean) and
### empirical 95% CI (by sampling percentile) of
### the HIA, for a given city, by age group as well as
### overall HIA by aggregating results for all ages.
### Results are based on using the function
### "hia_age_sim" --> this function provides point estimates + 95% CIs by age group and for all ages
###
#######################################
#######################################

hia_sim <- function(pollutant = "ndvi", 
                    outcome = "nat_mort",
                    nsim = 100, 
                    seed = NULL) {
  
  # set seed:
  if (!is.null(seed))
    set.seed(seed) ### it is important to set the seed to get replicable results every time you run the simulations
  
  # create an empty matrix to store point and CI estimates:
  nages <- length(allages) 
  est <- NULL
  
  # create an empty list to store the samples (needed to aggregate estimates for all ages):
  samp <- list()
  
  # get simulations:
  for (i in 1:nages) {
    aux <- hia_age_sim(age = allages[i],
                       pollutant = pollutant, 
                       outcome = outcome,
                       nsim = nsim)
    
    est <- rbind(est, aux$estimate)
    samp[[i]] <- aux$sample
  }
  
  # compute point and CI for all ages:
  samp <- bind_rows(samp, .id = "id")
  samp$id <- rep(1:nsim, times = nages, each = 4969)
  
  # for the total health outcome
  overallsamp <- samp %>% group_by(id, Poly_ID) %>% summarise(sum = sum(ndvi_outcome, na.rm = TRUE))
  estoverall <- overallsamp %>% group_by(Poly_ID) %>% summarise(mean = mean(sum, na.rm = TRUE), median = median(sum, na.rm = TRUE), pct2.5 = quantile(sum, probs = 2.5 / 100, na.rm = TRUE), pct97.5 = quantile(sum, probs = 97.5 / 100, na.rm = TRUE))
  
  # bind together overall results with results by age group 
  res <- as.data.frame(rbind(est, estoverall))
  # add age categories
  res$age <- rep(c(allages, "overall"), each = 4969)
  res <- res %>% arrange(Poly_ID, age)
  # add outcome
  res$outcome <- outcome
  
  # return results
  res <- res[, c("Poly_ID", "outcome", "age", "mean", "median", "pct2.5", "pct97.5")]
  return(res)
  
  
}

# example with 500 simulations
# check with different simulation numbers, the number of simulations should be enough to get robust results
# normally I start with a few simulations and increase the number until I don't see any substantial changes in the estimates
# also important to consider the balance between number of simulations and calculation time
# you can try with 100, 500, 1000 etc. 
t0 <- Sys.time()
hia2 <- hia_sim(pollutant = "ndvi",
                outcome = "nat_mort",
                nsim = 100, 
                seed = 666) ### important always to set the same seed for replicability!
t1 <- Sys.time()
t1 - t0  # --> see the calculation time
hia2




#########################################

# Apply functions to all health outcomes 

##########################################

##########################################
### NDVI
##########################################

# Natural-cause mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_ndvi_nat_mort <- hia_sim(pollutant = "ndvi",
                             outcome = "nat_mort",
                             nsim = 500, 
                             seed = 666)
sum(res_ndvi_nat_mort[res_ndvi_nat_mort$age == "overall", ]$mean) # 217 natural-cause deaths 

# CVD mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_ndvi_cvd_mort <- hia_sim(pollutant = "ndvi",
                             outcome = "cvd_mort",
                             nsim = 500, 
                             seed = 666)
sum(res_ndvi_cvd_mort[res_ndvi_cvd_mort$age == "overall", ]$mean) # 108 CVD deaths

# Hypertension
allages <- "total_adults" # define age groups to analyze
res_ndvi_hypertension <- hia_sim(pollutant = "ndvi",
                                 outcome = "hypertension",
                                 nsim = 500, 
                                 seed = 666)
sum(res_ndvi_hypertension[res_ndvi_hypertension$age == "overall", ]$mean) # 261 hypertension cases 

# Stroke
allages <- "total_adults" # define age groups to analyze
res_ndvi_stroke <- hia_sim(pollutant = "ndvi",
                           outcome = "stroke",
                           nsim = 500, 
                           seed = 666)
sum(res_ndvi_stroke[res_ndvi_stroke$age == "overall", ]$mean) # 47 stroke cases 


# Save the results
# Overall results by neighborhood
res_ndvi_rawlist <- list(res_ndvi_nat_mort, res_ndvi_cvd_mort, res_ndvi_hypertension, res_ndvi_stroke)
# Read census tracts
sofia_poly <- st_read("Baseline_BOD/Clean_data/sofia_pop_poly.geojson")[-c(7:39)]  # Sofia polygons 

for (i in 1:length(res_ndvi_rawlist)) {
  
  # Select results data frame
  res <- res_ndvi_rawlist[[i]]
  
  # Select results by age/sex
  res2 <- res[!res$age == "overall", ]
  # Merge population by age/sex
  pop <- sofia_pop_poly[c(2, 10:40)] %>% pivot_longer(cols = male_0_4:total_adults, names_to = "age", values_to = "pop")
  res2 <- merge(res2, pop, by = c("Poly_ID", "age"))
  poptot <- pop[-1] %>% group_by(age) %>% summarise(poptot = sum(pop, na.rm = T))
  
  # Select overall results
  res <- res[res$age == "overall", ]
  
  # Merge to population (adults) and health outcome rates
  aux <- sofia_pop_poly[c("Poly_ID", "total_adults")]
  res <- merge(res, aux, by.x = "Poly_ID", by.y = "Poly_ID")
  aux <- sofia_health[sofia_health$age %in% c("adults") & sofia_health$sex == "total", ]
  aux <- aux[c("outcome", "rate")]
  res <- merge(res, aux, by = "outcome")
  
  # Calculate expected cases
  res$exp_cases <- res$total_adults * (res$rate / 100000)
  
  # Calculate % of all cases
  res$per_cases <- res$mean / res$exp_cases * 100
  res$per_cases <- ifelse(is.na(res$per_cases), 0, res$per_cases)
  res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
  res$per_cases_lwr <- ifelse(is.na(res$per_cases_lwr), 0, res$per_cases_lwr)
  res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100
  res$per_cases_upr <- ifelse(is.na(res$per_cases_upr), 0, res$per_cases_upr)
  
  # Calculate attributable rate
  res2$rate <- res2$mean / res2$pop * 100000
  res2$rate <- ifelse(is.na(res2$rate), 0, res2$rate)
  res2$rate_lwr <- res2$pct2.5 / res2$pop * 100000
  res2$rate_lwr <- ifelse(is.na(res2$rate_lwr), 0, res2$rate_lwr)
  res2$rate_upr <- res2$pct97.5 / res2$pop * 100000
  res2$rate_upr <- ifelse(is.na(res2$rate_upr), 0, res2$rate_upr)
  res2 <- merge(res2, poptot, by = "age")
  res2 <- res2 %>% group_by(Poly_ID) %>% summarise(att_rate = sum(rate * poptot) / sum(poptot), att_rate_lwr = sum(rate_lwr * poptot) / sum(poptot), att_rate_upr = sum(rate_upr * poptot) / sum(poptot))
  res <- merge(res, res2, by = "Poly_ID")
  
  # res$att_rate <- res$mean / res$total_adults * 100000
  # res$att_rate <- ifelse(is.na(res$att_rate), 0, res$att_rate)
  # res$att_rate_lwr <- res$pct2.5 / res$total_adults * 100000
  # res$att_rate_lwr <- ifelse(is.na(res$att_rate_lwr), 0, res$att_rate_lwr)
  # res$att_rate_upr <- res$pct97.5 / res$total_adults * 100000
  # res$att_rate_upr <- ifelse(is.na(res$att_rate_upr), 0, res$att_rate_upr)
  
  # Export results
  res <- res[c("Poly_ID", "outcome", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
  # as csv
  write.csv(res, paste0("Baseline_BOD/Results/res_ndvi_", unique(res$outcome), ".csv"))
  # as Geojson
  res <- merge(sofia_poly, res, by.x = "Poly_ID", by.y = "Poly_ID")
  st_write(res, paste0("Baseline_BOD/Results/res_ndvi_", unique(res$outcome), ".geojson"))

}


# Overall results for the city
# Mortality outcomes
res_ndvi_rawlist <- list(res_ndvi_nat_mort, res_ndvi_cvd_mort)
res_ndvi_city_mort <- NULL # data frame to store results

for (i in 1:length(res_ndvi_rawlist)) {
  
  # Select results data frame
  res <- res_ndvi_rawlist[[i]]
  
  # Group results for the whole city
  # Select overall results
  res <- res[res$age == "overall", ]
  res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
  res$age <- ifelse(res$age == "overall", "adults", res$age)
  res$sex <- "total"
  
  # Add adult population and mortality rates by sex
  aux <- sofia_health[sofia_health$age %in% c("adults"), ]
  aux <- aux[c("outcome", "age", "sex", "rate")]
  res <- merge(aux, res, by.y = c("outcome", "age", "sex"), by.x = c("outcome", "age", "sex"))
  res$pop <- sum(sofia_pop_poly$total_adults)
  
  # Calculate expected cases
  res$exp_cases <- res$pop * (res$rate / 100000)
  
  # Calculate % of all cases
  res$per_cases <- res$mean / res$exp_cases * 100
  res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
  res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100
  
  # Calculate attributable rate
  res$att_rate <- res$mean / res$pop * 100000
  res$att_rate_lwr <- res$pct2.5 / res$pop * 100000
  res$att_rate_upr <- res$pct97.5 / res$pop * 100000
  
  # Order results
  res <- res[c("outcome", "age", "sex", "pop", "rate", "exp_cases", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
  res_ndvi_city_mort <- rbind(res_ndvi_city_mort, res)
  
  # Mortality outcomes -- by sex
  # Select results data frame
  res <- res_ndvi_rawlist[[i]]
  
  # Group results for the whole city by sex
  # Select overall results
  res <- res[!res$age == "overall", ]
  res$age <- sub("\\_.*", "", res$age)
  res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
  names(res)[2] <- "sex"
  res$age <- "adults"
  
  # Add adult population and mortality rates by sex
  aux <- sofia_mort[sofia_mort$age %in% c("adults"), ]
  aux <- aux[c("outcome", "age", "sex", "rate")]
  res <- merge(aux, res, by = c("outcome", "age", "sex"))
  sofia_pop_poly$adults_female <- rowSums(sofia_pop_poly[30:39])
  sofia_pop_poly$adults_male <- rowSums(sofia_pop_poly[15:24])
  res$pop <- c(sum(sofia_pop_poly$adults_female), sum(sofia_pop_poly$adults_male))
  
  # Calculate expected cases
  res$exp_cases <- res$pop * (res$rate / 100000)
  
  # Calculate % of all cases
  res$per_cases <- res$mean / res$exp_cases * 100
  res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
  res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100
  
  # Calculate attributable rate
  res$att_rate <- res$mean / res$pop * 100000
  res$att_rate_lwr <- res$pct2.5 / res$pop * 100000
  res$att_rate_upr <- res$pct97.5 / res$pop * 100000
  
  # Order results
  res <- res[c("outcome", "age", "sex", "pop", "rate", "exp_cases", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
  res_ndvi_city_mort <- rbind(res_ndvi_city_mort, res)
  
  # Mortality outcomes -- by age
  # Select results data frame
  res <- res_ndvi_rawlist[[i]]
  
  # Group results for the whole city by sex
  # Select overall results
  res <- res[!res$age == "overall", ]
  res$age <- sub(".*?_(.*)", "\\1", res$age)
  res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
  res$sex <- "total"
  
  # Add adult population and mortality rates by age
  aux <- sofia_mort[sofia_mort$sex %in% c("total"), ]
  aux <- aux[c("outcome", "age", "sex", "rate")]
  res <- merge(aux, res, by = c("outcome", "age", "sex"))
  res$pop <- c((sum(sofia_pop_poly$male_25_29)+sum(sofia_pop_poly$female_25_29)), (sum(sofia_pop_poly$male_30_34)+sum(sofia_pop_poly$female_30_34)), (sum(sofia_pop_poly$male_35_39)+sum(sofia_pop_poly$female_35_39)), (sum(sofia_pop_poly$male_40_44)+sum(sofia_pop_poly$female_40_44)), (sum(sofia_pop_poly$male_45_49)+sum(sofia_pop_poly$female_45_49)), (sum(sofia_pop_poly$male_50_54)+sum(sofia_pop_poly$female_50_54)), (sum(sofia_pop_poly$male_55_59)+sum(sofia_pop_poly$female_55_59)), (sum(sofia_pop_poly$male_60_64)+sum(sofia_pop_poly$female_60_64)), (sum(sofia_pop_poly$male_65_69)+sum(sofia_pop_poly$female_65_69)), (sum(sofia_pop_poly$`male_70+`)+sum(sofia_pop_poly$`female_70+`)))
  
  # Calculate expected cases
  res$exp_cases <- res$pop * (res$rate / 100000)
  
  # Calculate % of all cases
  res$per_cases <- res$mean / res$exp_cases * 100
  res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
  res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100
  
  # Calculate attributable rate
  res$att_rate <- res$mean / res$pop * 100000
  res$att_rate_lwr <- res$pct2.5 / res$pop * 100000
  res$att_rate_upr <- res$pct97.5 / res$pop * 100000
  
  # Order results
  res <- res[c("outcome", "age", "sex", "pop", "rate", "exp_cases", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
  res_ndvi_city_mort <- rbind(res_ndvi_city_mort, res)
  
}

# Disease outcomes
res_ndvi_rawlist <- list(res_ndvi_hypertension, res_ndvi_stroke)
res_ndvi_city_disease <- NULL # data frame to store results

for (i in 1:length(res_ndvi_rawlist)) {
  
  # Select results data frame
  res <- res_ndvi_rawlist[[i]]
  
  # Group results for the whole city 
  res <- res[res$age == "overall", ]
  res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
  res$age <- ifelse(res$age == "overall", "total", res$age)
 
  # Add adult population and incidence rate
  aux <- sofia_health[sofia_health$age %in% c("adults"), ]
  aux <- aux[c("outcome", "sex", "rate")]
  res <- merge(aux, res, by.y = c("outcome", "age"), by.x = c("outcome", "sex"))
  res$pop <- sum(sofia_pop_poly$total_adults)
  res$age <- "adults"
  
  # Calculate expected cases
  res$exp_cases <- res$pop * (res$rate / 100000)
  
  # Calculate % of all cases
  res$per_cases <- res$mean / res$exp_cases * 100
  res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
  res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100
  
  # Calculate attributable rate
  res$att_rate <- res$mean / res$pop * 100000
  res$att_rate_lwr <- res$pct2.5 / res$pop * 100000
  res$att_rate_upr <- res$pct97.5 / res$pop * 100000
  
  # Order results
  res <- res[c("outcome", "age", "sex", "pop", "rate", "exp_cases", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
  res_ndvi_city_disease <- rbind(res_ndvi_city_disease, res)
  
}

### Export the results
write_csv(res_ndvi_city_mort, "Baseline_BOD/Results/RESULTS_NDVI_mortality.csv")
write_csv(res_ndvi_city_disease, "Baseline_BOD/Results/RESULTS_NDVI_disease.csv")

# ### Plot the results
# # Mortality
res_ndvi_city_mort <- read_csv2("Baseline_BOD/Results/RESULTS_NDVI_mortality.csv")
res_ndvi_city_mort <- res_ndvi_city_mort[res_ndvi_city_mort$sex == "total", ]
res_ndvi_city_mort <- res_ndvi_city_mort[!res_ndvi_city_mort$age == "adults", ]

res_ndvi_city_mort$outcome <- factor(
  res_ndvi_city_mort$outcome,
  levels = c("nat_mort", "cvd_mort"), # Adjust these levels to your dataset
  labels = c("Natural-cause mortality", "Cardiovascular mortality") # Rename outcomes
)

theme_set(theme_bw())
ggplot(res_ndvi_city_mort, aes(x=age, y=mean)) + 
  geom_bar(stat = "identity", fill = "seagreen4") +
  facet_wrap(. ~ outcome, scales = "free") +
  coord_flip() + labs(title=expression("Lack of green attributable mortality by age group"),y ="Attributable deaths (n)", x = "") +
  geom_errorbar(aes(ymin = pct2.5, ymax = pct97.5), width = 0.2, color = "black")
ggsave("Baseline_BOD/Results/RESULTS_green_age.png")




