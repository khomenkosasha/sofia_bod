
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

# add columns with adult and children pop counts for disease analyses 
sofia_pop_poly$total_adults <- rowSums(sofia_pop_poly[15:24]) + rowSums(sofia_pop_poly[30:39])
sofia_pop_poly$total_children <- rowSums(sofia_pop_poly[10:12]) + rowSums(sofia_pop_poly[25:27])

### Health data 
sofia_mort <- read_csv2("Baseline_BOD/Health/Sofia_mortality.csv")
sofia_disease <- read_csv2("Baseline_BOD/Health/Sofia_diseases.csv")
sofia_mort$rate <- as.numeric(sofia_mort$count/ sofia_mort$population * 100000)
sofia_disease$rate <- as.numeric(sofia_disease$count/ sofia_disease$population * 100000)
sofia_health <- rbind(sofia_mort, sofia_disease)

### Air pollution data
sofia_airp_poly <- st_read("Baseline_BOD/Clean_data/sofia_airp_poly.geojson")
sofia_airp_poly$geometry <- NULL


#--------------------------------------------------
# First, we need to create some vectors and data 
# frames with parameters that will be used later 
# in the Monte Carlo analysis
#--------------------------------------------------

### We extract the age and sex labels
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # these need to be adjusted for each analysis based on the outcome
allages

### We create a data frame with the RRs for each pollutant and health outcome

rranalysis <- data.frame(pollutant = rep(c("pm25", "no2"), each = 12),
                         outcome = rep(c("nat_mort", "cvd_mort", "resp_mort", "lungc_mort", "ihd", "stroke", "asthma_children", "asthma_adults", "copd", "diabetes", "hypertension", "alri_children"), times = 2),
                         rr = c(1.095, 1.127, 1.136, 1.093, 1.13, 1.16, 1.34, NA, 1.18, 1.10, 1.17, NA, 1.05, 1.05, 1.05, 1.07, NA, NA, 1.10, 1.10, NA, NA, NA, 1.09),
                         rrlo95 = c(1.064, 1.102, 1.079, 1.053, 1.05, 1.12, 1.10, NA, 1.13, 1.03, 1.05, NA, 1.03, 1.03, 1.03,1.04, NA, NA, 1.05, 1.01, NA, NA, NA, 1.03),
                         rrup95 = c(1.127, 1.152, 1.197, 1.135, 1.22, 1.20, 1.63, NA, 1.23, 1.18, 1.30, NA, 1.07, 1.08, 1.07, 1.10, NA, NA, 1.18, 1.21, NA, NA, NA, 1.16), 
                         scale = c(10, 10, 10, 10, 10, 10, 10, NA, 10, 10, 10, NA, 10, 10, 10, 10, NA, NA, 10, 10, NA, NA, NA, 10))

rranalysis_ozone <- data.frame(pollutant = "o3", 
                               outcome = "resp_mort", 
                               rr = 1.05, 
                               rrlo95 = 1.02,
                               rrup95 = 1.08, 
                               scale = 10)

rranalysis <- rbind(rranalysis, rranalysis_ozone)

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
hia_age <- function(age = allages[19],
                    pollutant = "o3", 
                    outcome = "resp_mort", 
                    adjustment = "no") {  # whether to provide coefficient adjusted individual PM2.5 and NO2 results                
  
  ### Uncertainty in RRs:
  
  if(pollutant == "multiexp"){
    
    # We retrieve the RR and its CIs for both PM2.5 and NO2
    # For PM2.5:
    rrci <- as.numeric(rranalysis[rranalysis$pollutant == "pm25" & rranalysis$outcome == outcome, c("rr", "rrlo95", "rrup95")])
    names(rrci) <- c("rr", "rrlo95", "rrup95")
    scale_pm25 <- as.numeric(rranalysis[rranalysis$pollutant == "pm25" & rranalysis$outcome == outcome, "scale"])
    
    # Then, we simulate a RR value assuming log(RR) is normal 
    selogrr <- log(rrci["rrup95"] / rrci["rrlo95"]) / (2 * qnorm(0.975)) # we calculate the standard deviation of log(RR)
    logrr <- log(rrci["rr"]) # we calculate log(RR)
    logrrsim <- rnorm(n = 1, mean = logrr, sd = selogrr) # we simulate one log(RR) value from a normal distribution
    rr_pm25 <- exp(logrrsim) # We exponentiate back the simulated log(RR) --> the simulated RR will be then used in the HIA calculation
    
    # For NO2:
    rrci <- as.numeric(rranalysis[rranalysis$pollutant == "no2" & rranalysis$outcome == outcome, c("rr", "rrlo95", "rrup95")])
    names(rrci) <- c("rr", "rrlo95", "rrup95")
    scale_no2 <- as.numeric(rranalysis[rranalysis$pollutant == "no2" & rranalysis$outcome == outcome, "scale"])
    
    # Then, we simulate a RR value assuming log(RR) is normal 
    selogrr <- log(rrci["rrup95"] / rrci["rrlo95"]) / (2 * qnorm(0.975)) # we calculate the standard deviation of log(RR)
    logrr <- log(rrci["rr"]) # we calculate log(RR)
    logrrsim <- rnorm(n = 1, mean = logrr, sd = selogrr) # we simulate one log(RR) value from a normal distribution
    rr_no2 <- exp(logrrsim) # We exponentiate back the simulated log(RR) --> the simulated RR will be then used in the HIA calculation
    
  }else{
    
    # We retrieve the RR and its CIs based on the pollutant label:
    rrci <- as.numeric(rranalysis[rranalysis$pollutant == pollutant & rranalysis$outcome == outcome, c("rr", "rrlo95", "rrup95")])
    names(rrci) <- c("rr", "rrlo95", "rrup95")
    scale <- as.numeric(rranalysis[rranalysis$pollutant == pollutant & rranalysis$outcome == outcome, "scale"])
    
    # Then, we simulate a RR value assuming log(RR) is normal 
    selogrr <- log(rrci["rrup95"] / rrci["rrlo95"]) / (2 * qnorm(0.975)) # we calculate the standard deviation of log(RR)
    logrr <- log(rrci["rr"]) # we calculate log(RR)
    logrrsim <- rnorm(n = 1, mean = logrr, sd = selogrr) # we simulate one log(RR) value from a normal distribution
    rr <- exp(logrrsim) # We exponentiate back the simulated log(RR) --> the simulated RR will be then used in the HIA calculation
    
  }

  # We establish the counterfactual exposure level:
  if(pollutant == "pm25"){
    tmrel <- 5
    
  }
  if(pollutant == "no2"){
    tmrel <- 10
    
  }
  if(pollutant == "o3"){
    tmrel <- 48.7
    
    # a database in Europe documented a difference between peak season and annual ozone of 1.24
    # using this ratio the 60 ug/m2 WHO recommendation for peak exposure can be translated to 48.7 for annual exposure  
    # https://www.sciencedirect.com/science/article/abs/pii/S0160412018309759?via%3Dihub 
    # WHO guidelines  
    
  }
  if(pollutant == "multiexp"){
    tmrel_pm25 <- 5
    tmrel_no2 <- 10
    
  }

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
  if(pollutant == "multiexp"){
    
    # Retrieve both pollutants
    pollgrid <-  sofia_airp_poly[, c("Poly_ID", "pm25", "no2")]
    
  }else{
    
    # Retrieve one pollutant based on pollutant label
    pollgrid <-  sofia_airp_poly[, c("Poly_ID", pollutant)]
    names(pollgrid)[2] <- "ap"
    
  }
  
  ### Next, we merge the air pollution and mortality data
  popcity <- merge(popcity, pollgrid, by = "Poly_ID")
  rm(pollgrid)
  
  ### We simulate one error value for the pollutant level for each grid (using the reported RMSE from the models):
  # we retrieve the RMSE values:
  if(pollutant == "pm25"){ # from EXPANSE models
    rmse <- 2.31
  }
  if(pollutant == "no2"){
    rmse <- 6.33 # EXPANSE
    # rmse <- 9.98  # alternative NO2 model
  }
  if(pollutant == "o3"){
    rmse <- 6.79
  }
  if(pollutant == "multiexp"){
    rmse_pm25 <- 2.31
    rmse_no2 <- 6.33
  }
  
  # We calculate the standard error of model value:
  if(pollutant == "multiexp"){
    
    se_pm25 <- rmse_pm25 * qnorm(0.975)
    se_no2 <- rmse_no2 * qnorm(0.975)
    # We simulate the error at each grid (assuming a normal distribution):
    aperror_pm25 <- rnorm(n = dim(popcity)[1], mean = 0, sd = se_pm25)
    aperror_no2 <- rnorm(n = dim(popcity)[1], mean = 0, sd = se_no2)
    # Finally, we update the value of each pollutant by adding the error to the point estimate:
    popcity$pm25 <- popcity$pm25 + aperror_pm25
    popcity$no2 <- popcity$no2 + aperror_no2
    
  }else{
    
    se <- rmse * qnorm(0.975)
    # We simulate the error at each grid (assuming a normal distribution):
    aperror <- rnorm(n = dim(popcity)[1], mean = 0, sd = se)
    # Finally, we update the value of "ap" by adding the error to the point estimate:
    popcity$ap <- popcity$ap + aperror
    
  }
  
  ### Now that all values are simulated, we do the HIA: calculate exposure difference, scaled RR and PAF
  if(pollutant == "multiexp"){
    
    # Calculate the exposure difference
    popcity$expdiff_pm25 <- popcity$pm25 - tmrel_pm25
    popcity$expdiff_no2 <- popcity$no2 - tmrel_no2
    popcity$expdiff_pm25[popcity$expdiff_pm25 < 0] <- 0  # zero if ap <= tmrel
    popcity$expdiff_no2[popcity$expdiff_no2 < 0] <- 0  # zero if ap <= tmrel
    
    # Adjust RR based on coefficient difference
    # For PM2.5 
    rr_pm25 <- exp((log(rr_pm25) / scale_pm25) * 5) # scale PM2.5 RR to 5 ug/m3 to match the scale of the coefficient difference
    rr_pm25 <- exp(log(rr_pm25) - 0.017)
    rr_pm25 <- exp((log(rr_pm25) / 5) * scale_pm25) # scale back to 10 ug/m3
    
    # For NO2
    rr_no2 <- exp(log(rr_no2) - 0.007) # coefficient difference can be applied directly as it is calculated for 10 ug/m3 increment
    
    # Scale adjusted RR for PM2.5 and NO2 to exposure difference
    popcity$rr_pm25 <- exp((log(rr_pm25) / scale_pm25) * popcity$expdiff_pm25)
    popcity$rr_no2 <- exp((log(rr_no2) / scale_no2) * popcity$expdiff_no2)
    
    # Calculate Joint RR 
    popcity$joint_rr <- popcity$rr_pm25 * popcity$rr_no2 
    
    # Calculate joint PAF
    popcity$PAF <- with(popcity, (joint_rr - 1) / joint_rr)
    
  }else{
    
    # Calculate exposure difference
    popcity$expdiff <- popcity$ap - tmrel
    popcity$expdiff[popcity$expdiff < 0] <- 0  # zero if ap <= tmrel
    
    # If needed, adjust the RR based on coefficient difference
    if(adjustment == "yes"){
      
      if(pollutant == "pm25"){
        
        rr <- exp((log(rr) / scale) * 5) # scale PM2.5 RR to 5 ug/m3 to match the scale of the coefficient difference
        rr <- exp(log(rr) - 0.017)
        rr <- exp((log(rr) / 5) * scale) # scale back to 10 ug/m3
        
      }
      
      if(pollutant == "no2"){
        
        rr <- exp(log(rr) - 0.007) # coefficient difference can be applied directly as it is calculated for 10 ug/m3 increment
        
      }
      
    }
    
    # Scale RR to exposure difference
    popcity$rr <- exp((log(rr) / scale) * popcity$expdiff)
    
    # Calculate PAF
    popcity$PAF <- with(popcity, (rr - 1) / rr)
    
  }
  
  ### We calculate health outcome due to air pollution exposure in each neighborhood
  popcity$apoutcome <- with(popcity, outcome * PAF)
  popcity <- popcity[c("Poly_ID", "apoutcome")]
  #names(popcity)[1] <- "ct_id"
  
  return(popcity) ### We ask the function to return the health outcome estimate by neighborhood
}

### Example --> you will see that if you run this several times you get a different value each time (because of the uncertainties)
hia1 <- hia_age(age = allages[10], pollutant = "pm25", outcome = "nat_mort", adjustment = "yes")
hia1
sum(hia1$apoutcome)



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
                        pollutant = "pm25", 
                        outcome = "nat_mort",
                        adjustment = "no",
                        nsim = 100) {                  # number of simulations
  
  # We use the replication function to replicate the previously defined function n times (to get a sample of values):
  samp <- replicate(n = nsim,
                    expr = hia_age(age = age,
                                   pollutant = pollutant, 
                                   outcome = outcome, 
                                   adjustment = adjustment), simplify = FALSE)
                    
  # We use the sample of values to get point (mean) and CI estimate (2.5 and 97.5 percentiles)
  summ <- bind_rows(samp, .id = "id") %>% 
    group_by(Poly_ID) %>%
    summarise(mean = mean(apoutcome, na.rm = TRUE), 
              median = median(apoutcome, na.rm = TRUE), 
              pct2.5 = quantile(apoutcome, probs = 2.5/100, na.rm = TRUE), 
              pct97.5 = quantile(apoutcome, probs = 97.5/100, na.rm = TRUE))
  
  # We ask the function to return a list with results
  res <- list(sample = samp, estimate = summ)
  return(res)
}

# example with 500 simulations --> you can try different numbers of simulations and check how long it takes to calculate
t0 <- Sys.time()
simhia1 <- hia_age_sim(age = allages[10],
                       pollutant = "pm25", 
                       outcome = "nat_mort",
                       adjustment = "yes",
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

hia_sim <- function(pollutant = "pm25", 
                    outcome = "nat_mort",
                    adjustment = "yes",
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
                       adjustment = adjustment,
                       nsim = nsim)
    
    est <- rbind(est, aux$estimate)
    samp[[i]] <- aux$sample
  }
  
  # compute point and CI for all ages:
  samp <- bind_rows(samp, .id = "id")
  samp$id <- rep(1:nsim, times = nages, each = 4969)
  
  # for the total health outcome
  overallsamp <- samp %>% group_by(id, Poly_ID) %>% summarise(sum = sum(apoutcome, na.rm = TRUE))
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
hia2 <- hia_sim(pollutant = "pm25",
                outcome = "nat_mort",
                adjustment = "yes",
                nsim = 100, 
                seed = 666) ### important always to set the same seed for replicability!
t1 <- Sys.time()
t1 - t0  # --> see the calculation time
hia2




#########################################

# Apply functions to all health outcomes and pollutants

##########################################

##########################################
### PM2.5
##########################################

# Natural-cause mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_pm25_nat_mort <- hia_sim(pollutant = "pm25",
                             outcome = "nat_mort",
                             adjustment = "no",
                             nsim = 500, 
                             seed = 666)
sum(res_pm25_nat_mort[res_pm25_nat_mort$age == "overall", ]$mean) # 1939 natural-cause deaths

# Natural-cause mortality --> using PM2.5 RR adjusted by coefficient difference
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_pm25_nat_mort_adj <- hia_sim(pollutant = "pm25",
                             outcome = "nat_mort",
                             adjustment = "yes",
                             nsim = 500, 
                             seed = 666)
sum(res_pm25_nat_mort_adj[res_pm25_nat_mort_adj$age == "overall", ]$mean) # 1248 natural-cause deaths --> 35% less deaths

# Natural-cause mortality --> multiexposure with coefficient adjustment
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_multiexp_nat_mort <- hia_sim(pollutant = "multiexp",
                                 outcome = "nat_mort",
                                 adjustment = "no",
                                 nsim = 500, 
                                 seed = 666)
sum(res_multiexp_nat_mort[res_multiexp_nat_mort$age == "overall", ]$mean) # 2177 natural-cause deaths

# CVD mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_pm25_cvd_mort <- hia_sim(pollutant = "pm25",
                             outcome = "cvd_mort",
                             adjustment = "no",
                             nsim = 500, 
                             seed = 666)
sum(res_pm25_cvd_mort[res_pm25_cvd_mort$age == "overall", ]$mean) # 1649 CVD deaths 

# Respiratory mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_pm25_resp_mort <- hia_sim(pollutant = "pm25",
                              outcome = "resp_mort",
                              adjustment = "no",
                              nsim = 500, 
                              seed = 666)
sum(res_pm25_resp_mort[res_pm25_resp_mort$age == "overall", ]$mean) # 123 respiratory deaths 

# Lung cancer mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_pm25_lungc_mort <- hia_sim(pollutant = "pm25",
                               outcome = "lungc_mort",
                               adjustment = "no",
                               nsim = 500, 
                               seed = 666)
sum(res_pm25_lungc_mort[res_pm25_lungc_mort$age == "overall", ]$mean) # 57 lung cancer deaths 

# Hypertension
allages <- "total_adults" # define age groups to analyze
res_pm25_hypertension <- hia_sim(pollutant = "pm25",
                                 outcome = "hypertension",
                                 adjustment = "no",
                                 nsim = 500, 
                                 seed = 666)
sum(res_pm25_hypertension[res_pm25_hypertension$age == "overall", ]$mean) # 7678 hypertension cases 

# IHD
allages <- "total_adults" # define age groups to analyze
res_pm25_ihd <- hia_sim(pollutant = "pm25",
                        outcome = "ihd",
                        adjustment = "no",
                        nsim = 500, 
                        seed = 666)
sum(res_pm25_ihd[res_pm25_ihd$age == "overall", ]$mean) # 2073 IHD cases 

# Stroke
allages <- "total_adults" # define age groups to analyze
res_pm25_stroke <- hia_sim(pollutant = "pm25",
                           outcome = "stroke",
                           adjustment = "no",
                           nsim = 500, 
                           seed = 666)
sum(res_pm25_stroke[res_pm25_stroke$age == "overall", ]$mean) # 1318 stroke cases 

# Asthma children
allages <- "total_children" # define age groups to analyze
res_pm25_asthma_children <- hia_sim(pollutant = "pm25",
                                    outcome = "asthma_children",
                                    adjustment = "no",
                                    nsim = 500, 
                                    seed = 666)
sum(res_pm25_asthma_children[res_pm25_asthma_children$age == "overall", ]$mean) # 362 children asthma cases 

# COPD
allages <- "total_adults" # define age groups to analyze
res_pm25_copd <- hia_sim(pollutant = "pm25",
                         outcome = "copd",
                         adjustment = "no",
                         nsim = 500, 
                         seed = 666)
sum(res_pm25_copd[res_pm25_copd$age == "overall", ]$mean) # 1324 COPD cases 

# Diabetes
allages <- "total_adults" # define age groups to analyze
res_pm25_diabetes <- hia_sim(pollutant = "pm25",
                             outcome = "diabetes",
                             adjustment = "no",
                             nsim = 500, 
                             seed = 666)
sum(res_pm25_diabetes[res_pm25_diabetes$age == "overall", ]$mean) # 1249 diabetes cases 


# Save the results
# Overall results by neighborhood
res_pm25_rawlist <- list(res_pm25_nat_mort, res_pm25_nat_mort_adj, res_multiexp_nat_mort, res_pm25_cvd_mort, res_pm25_resp_mort, res_pm25_lungc_mort, res_pm25_hypertension, res_pm25_ihd, res_pm25_stroke, res_pm25_asthma_children, res_pm25_copd, res_pm25_diabetes)
# Read census tracts
sofia_poly <- st_read("Baseline_BOD/Clean_data/sofia_pop_poly.geojson")[-c(7:39)]  # Sofia polygons 

for (i in 1:length(res_pm25_rawlist)) {
  
  # Select results data frame
  res <- res_pm25_rawlist[[i]]
  
  # Select results by age/sex
  res2 <- res[!res$age == "overall", ]
  # Merge population by age/sex
  pop <- sofia_pop_poly[c(2, 10:41)] %>% pivot_longer(cols = male_0_4:total_children, names_to = "age", values_to = "pop")
  res2 <- merge(res2, pop, by = c("Poly_ID", "age"))
  poptot <- pop[-1] %>% group_by(age) %>% summarise(poptot = sum(pop, na.rm = T))
  
  # Select overall results
  res <- res[res$age == "overall", ]
  
  # Merge to population (adults, children) and health outcome rates
  if(unique(res$outcome) == "asthma_children"){
    
    aux <- sofia_pop_poly[c("Poly_ID", "total_children")]
    res <- merge(res, aux, by.x = "Poly_ID", by.y = "Poly_ID")
    aux <- sofia_health[sofia_health$age %in% c("children") & sofia_health$sex == "total", ]
    aux <- aux[c("outcome", "rate")]
    res <- merge(res, aux, by = "outcome")
    
    # Calculate expected cases
    res$exp_cases <- res$total_children * (res$rate / 100000)
    
    # Calculate % of all cases
    res$per_cases <- res$mean / res$exp_cases * 100
    res$per_cases <- ifelse(is.na(res$per_cases), 0, res$per_cases)
    res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
    res$per_cases_lwr <- ifelse(is.na(res$per_cases_lwr), 0, res$per_cases_lwr)
    res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100
    res$per_cases_upr <- ifelse(is.na(res$per_cases_upr), 0, res$per_cases_upr)
    
    # Calculate attributable rate
    res$att_rate <- res$mean / res$total_children * 100000
    res$att_rate <- ifelse(is.na(res$att_rate), 0, res$att_rate)
    res$att_rate_lwr <- res$pct2.5 / res$total_children * 100000
    res$att_rate_lwr <- ifelse(is.na(res$att_rate_lwr), 0, res$att_rate_lwr)
    res$att_rate_upr <- res$pct97.5 / res$total_children * 100000
    res$att_rate_upr <- ifelse(is.na(res$att_rate_upr), 0, res$att_rate_upr)
    
  }else{
    
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
    
  }
  
  # Export results
  res <- res[c("Poly_ID", "outcome", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
  # as csv
  if(i == 2){
    write.csv(res, "Baseline_BOD/Results/res_pm25_nat_mort_adj.csv")
  }else if(i == 3){
    write.csv(res, "Baseline_BOD/Results/res_pm25_nat_mort_multiexp.csv")
  }else{
    write.csv(res, paste0("Baseline_BOD/Results/res_pm25_", unique(res$outcome), ".csv"))
  }
  
  # as Geojson
  res <- merge(sofia_poly, res, by.x = "Poly_ID", by.y = "Poly_ID")
  if(i == 2){
    st_write(res, "Baseline_BOD/Results/res_pm25_nat_mort_adj.geojson")
  }else if(i == 3){
    st_write(res, "Baseline_BOD/Results/res_pm25_nat_mort_multiexp.geojson")
  }else{
    st_write(res, paste0("Baseline_BOD/Results/res_pm25_", unique(res$outcome), ".geojson"))
  }

}

# Overall results for the city
# Mortality outcomes
res_pm25_rawlist <- list(res_pm25_nat_mort, res_pm25_nat_mort_adj, res_multiexp_nat_mort, res_pm25_cvd_mort, res_pm25_resp_mort, res_pm25_lungc_mort)
res_pm25_city_mort <- NULL # data frame to store results

for (i in 1:length(res_pm25_rawlist)) {
  
  # Select results data frame
  res <- res_pm25_rawlist[[i]]
  
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
  res_pm25_city_mort <- rbind(res_pm25_city_mort, res)
  
  # Mortality outcomes -- by sex
  # Select results data frame
  res <- res_pm25_rawlist[[i]]
  
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
  res_pm25_city_mort <- rbind(res_pm25_city_mort, res)
  
  # Mortality outcomes -- by age
  # Select results data frame
  res <- res_pm25_rawlist[[i]]
  
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
  res_pm25_city_mort <- rbind(res_pm25_city_mort, res)

}

# Disease outcomes
res_pm25_rawlist <- list(res_pm25_hypertension, res_pm25_ihd, res_pm25_stroke, res_pm25_asthma_children, res_pm25_copd, res_pm25_diabetes)
res_pm25_city_disease <- NULL # data frame to store results

for (i in 1:length(res_pm25_rawlist)) {
  
  # Select results data frame
  res <- res_pm25_rawlist[[i]]
  
  # Group results for the whole city 
  res <- res[res$age == "overall", ]
  res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
  res$age <- ifelse(res$age == "overall", "total", res$age)
  
  if(res$outcome == "asthma_children"){
    
    # Add children population and incidence rate
    aux <- sofia_health[sofia_health$age %in% c("children"), ]
    aux <- aux[c("outcome", "sex", "rate")]
    res <- merge(aux, res, by.y = c("outcome", "age"), by.x = c("outcome", "sex"))
    res$pop <- sum(sofia_pop_poly$total_children)
    res$age <- "children"
    
  }else{
    
    # Add adult population and incidence rate
    aux <- sofia_health[sofia_health$age %in% c("adults"), ]
    aux <- aux[c("outcome", "sex", "rate")]
    res <- merge(aux, res, by.y = c("outcome", "age"), by.x = c("outcome", "sex"))
    res$pop <- sum(sofia_pop_poly$total_adults)
    res$age <- "adults"
  
  }

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
  res <- res[c("outcome","age", "sex", "pop", "rate", "exp_cases", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
  res_pm25_city_disease <- rbind(res_pm25_city_disease, res)
  
}

### Export the results
write_csv(res_pm25_city_mort, "Baseline_BOD/Results/RESULTS_PM25_mortality.csv")
write_csv(res_pm25_city_disease, "Baseline_BOD/Results/RESULTS_PM25_disease.csv")


# ### Plot the results
# Mortality
res_pm25_city_mort <- read_csv2("Baseline_BOD/Results/RESULTS_PM25_mortality.csv")
res_pm25_city_mort <- res_pm25_city_mort[res_pm25_city_mort$sex == "total", ]
res_pm25_city_mort <- res_pm25_city_mort[!res_pm25_city_mort$age == "adults", ]
res_pm25_city_mort <- res_pm25_city_mort[-c(11:30), ]

res_pm25_city_mort$outcome <- factor(
  res_pm25_city_mort$outcome,
  levels = c("nat_mort", "cvd_mort", "resp_mort", "lungc_mort"), # Adjust these levels to your dataset
  labels = c("Natural-cause mortality", "Cardiovascular mortality", "Respiratory mortality", "Lung cancer mortality") # Rename outcomes
)

theme_set(theme_bw())
ggplot(res_pm25_city_mort, aes(x=age, y=mean)) + 
  geom_bar(stat = "identity", fill = "skyblue2") +
  facet_wrap(. ~ outcome, scales = "free") +
  coord_flip() + labs(title=expression("PM"[2.5]*" attributable mortality by age group"),y ="Attributable deaths (n)", x = "") +
  geom_errorbar(aes(ymin = pct2.5, ymax = pct97.5), width = 0.2, color = "black")
ggsave("Baseline_BOD/Results/RESULTS_PM25_age.png")

# # Disease outcomes
# res_pm25_city_disease <- read_csv2("Baseline_BOD/Results/RESULTS_PM25_disease.csv")
# p2 <- ggplot(res_pm25_city_disease, aes(x=reorder(outcome, mean), y=mean)) + 
#              geom_bar(stat = "identity", fill = "mediumpurple4") +
#              coord_flip() + labs(title="PM2.5 attributable morbidity",y ="Attributable cases (n)", x = "") +
#              scale_x_discrete(labels=c("hypertension" = "Hypertension", "copd" = "COPD", "stroke" = "Stroke", "ihd" = "IHD", "asthma_adults" = "Asthma (adults)", "asthma_children" = "Asthma (children)", "diabetes" ="Diabetes (type II)"))
# 



##########################################
### NO2
##########################################

# Natural-cause mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_no2_nat_mort <- hia_sim(pollutant = "no2",
                            outcome = "nat_mort",
                            adjustment = "no",
                            nsim = 500, 
                            seed = 666)
sum(res_no2_nat_mort[res_no2_nat_mort$age == "overall", ]$mean) # 1172 natural-cause deaths 

# Natural-cause mortality --> using NO2 RR with coefficient adjustment
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_no2_nat_mort_adj <- hia_sim(pollutant = "no2",
                            outcome = "nat_mort",
                            adjustment = "yes",
                            nsim = 500, 
                            seed = 666)
sum(res_no2_nat_mort_adj[res_no2_nat_mort_adj$age == "overall", ]$mean) # 1013 natural-cause deaths --> 14% reduction in deaths

# CVD mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_no2_cvd_mort <- hia_sim(pollutant = "no2",
                            outcome = "cvd_mort",
                            adjustment = "no",
                            nsim = 500, 
                            seed = 666)
sum(res_no2_cvd_mort[res_no2_cvd_mort$age == "overall", ]$mean) # 776 CVD deaths 

# Respiratory mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_no2_resp_mort <- hia_sim(pollutant = "no2",
                             outcome = "resp_mort",
                             adjustment = "no",
                             nsim = 500, 
                             seed = 666)
sum(res_no2_resp_mort[res_no2_resp_mort$age == "overall", ]$mean) # 55 respiratory deaths 

# Lung cancer mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_no2_lungc_mort <- hia_sim(pollutant = "no2",
                              outcome = "lungc_mort",
                              adjustment = "no",
                              nsim = 500, 
                              seed = 666)
sum(res_no2_lungc_mort[res_no2_lungc_mort$age == "overall", ]$mean) # 47 lung cancer deaths 

# Asthma adults
allages <- "total_adults" # define age groups to analyze
res_no2_asthma_adults <- hia_sim(pollutant = "no2",
                                 outcome = "asthma_adults",
                                 adjustment = "no",
                                 nsim = 500, 
                                 seed = 666)
sum(res_no2_asthma_adults[res_no2_asthma_adults$age == "overall", ]$mean) # 207 adult asthma cases 

# Asthma children
allages <- "total_children" # define age groups to analyze
res_no2_asthma_children <- hia_sim(pollutant = "no2",
                                   outcome = "asthma_children",
                                   adjustment = "no",
                                   nsim = 500, 
                                   seed = 666)
sum(res_no2_asthma_children[res_no2_asthma_children$age == "overall", ]$mean) # 146 children asthma cases 

# ALRI children
allages <- "total_children" # define age groups to analyze
res_no2_alri <- hia_sim(pollutant = "no2",
                        outcome = "alri_children",
                        adjustment = "no",
                        nsim = 500, 
                        seed = 666)
sum(res_no2_alri[res_no2_alri$age == "overall", ]$mean) # 901 ALRI cases 


# Save the results
# Overall results by neighborhood
res_no2_rawlist <- list(res_no2_nat_mort, res_no2_nat_mort_adj, res_no2_cvd_mort, res_no2_resp_mort, res_no2_lungc_mort, res_no2_asthma_adults, res_no2_asthma_children, res_no2_alri)
# Read census tracts
sofia_poly <- st_read("Baseline_BOD/Clean_data/sofia_pop_poly.geojson")[-c(7:39)]  # Sofia polygons 

for (i in 1:length(res_no2_rawlist)) {
  
  # Select results data frame
  res <- res_no2_rawlist[[i]]
  
  # Select results by age/sex
  res2 <- res[!res$age == "overall", ]
  # Merge population by age/sex
  pop <- sofia_pop_poly[c(2, 10:41)] %>% pivot_longer(cols = male_0_4:total_children, names_to = "age", values_to = "pop")
  res2 <- merge(res2, pop, by = c("Poly_ID", "age"))
  poptot <- pop[-1] %>% group_by(age) %>% summarise(poptot = sum(pop, na.rm = T))
  
  # Select overall results
  res <- res[res$age == "overall", ]
  
  # Merge to population (adults, children) and health outcome rates
  if(unique(res$outcome) %in% c("asthma_children", "alri_children")){
    
    aux <- sofia_pop_poly[c("Poly_ID", "total_children")]
    res <- merge(res, aux, by.x = "Poly_ID", by.y = "Poly_ID")
    aux <- sofia_health[sofia_health$age %in% c("children") & sofia_health$sex == "total", ]
    aux <- aux[c("outcome", "rate")]
    res <- merge(res, aux, by = "outcome")
    
    # Calculate expected cases
    res$exp_cases <- res$total_children * (res$rate / 100000)
    
    # Calculate % of all cases
    res$per_cases <- res$mean / res$exp_cases * 100
    res$per_cases <- ifelse(is.na(res$per_cases), 0, res$per_cases)
    res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
    res$per_cases_lwr <- ifelse(is.na(res$per_cases_lwr), 0, res$per_cases_lwr)
    res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100
    res$per_cases_upr <- ifelse(is.na(res$per_cases_upr), 0, res$per_cases_upr)
    
    # Calculate attributable rate
    res$att_rate <- res$mean / res$total_children * 100000
    res$att_rate <- ifelse(is.na(res$att_rate), 0, res$att_rate)
    res$att_rate_lwr <- res$pct2.5 / res$total_children * 100000
    res$att_rate_lwr <- ifelse(is.na(res$att_rate_lwr), 0, res$att_rate_lwr)
    res$att_rate_upr <- res$pct97.5 / res$total_children * 100000
    res$att_rate_upr <- ifelse(is.na(res$att_rate_upr), 0, res$att_rate_upr)
    
  }else{
    
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
    
  }
  
  # Export results
  res <- res[c("Poly_ID", "outcome", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
  # as csv
  if(i == 2){
    write.csv(res, "Baseline_BOD/Results/res_no2_nat_mort_adj.csv")
  }else{
    write.csv(res, paste0("Baseline_BOD/Results/res_no2_", unique(res$outcome), ".csv"))
  }
  
  # as Geojson
  res <- merge(sofia_poly, res, by.x = "Poly_ID", by.y = "Poly_ID")
  if(i == 2){
    st_write(res, "Baseline_BOD/Results/res_no2_nat_mort_adj.geojson")
  }else{
    st_write(res, paste0("Baseline_BOD/Results/res_no2_", unique(res$outcome), ".geojson"))
  }
  
}

# Overall results for the city
# Mortality outcomes
res_no2_rawlist <- list(res_no2_nat_mort, res_no2_nat_mort_adj, res_no2_cvd_mort, res_no2_resp_mort, res_no2_lungc_mort)
res_no2_city_mort <- NULL # data frame to store results

for (i in 1:length(res_no2_rawlist)) {
  
  # Select results data frame
  res <- res_no2_rawlist[[i]]
  
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
  res_no2_city_mort <- rbind(res_no2_city_mort, res)
  
  # Mortality outcomes -- by sex
  # Select results data frame
  res <- res_no2_rawlist[[i]]
  
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
  res_no2_city_mort <- rbind(res_no2_city_mort, res)
  
  # Mortality outcomes -- by age
  # Select results data frame
  res <- res_no2_rawlist[[i]]
  
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
  res_no2_city_mort <- rbind(res_no2_city_mort, res)
  
}

# Disease outcomes
res_no2_rawlist <- list(res_no2_asthma_adults, res_no2_asthma_children, res_no2_alri)
res_no2_city_disease <- NULL # data frame to store results

for (i in 1:length(res_no2_rawlist)) {
  
  # Select results data frame
  res <- res_no2_rawlist[[i]]
  
  # Group results for the whole city 
  res <- res[res$age == "overall", ]
  res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
  res$age <- ifelse(res$age == "overall", "total", res$age)
  
  if(res$outcome %in% c("asthma_children", "alri_children")){
    
    # Add children population and incidence rate
    aux <- sofia_health[sofia_health$age %in% c("children"), ]
    aux <- aux[c("outcome", "sex", "rate")]
    res <- merge(aux, res, by.y = c("outcome", "age"), by.x = c("outcome", "sex"))
    res$pop <- sum(sofia_pop_poly$total_children)
    res$age <- "children"
    
  }else{
    
    # Add adult population and incidence rate
    aux <- sofia_health[sofia_health$age %in% c("adults"), ]
    aux <- aux[c("outcome", "sex", "rate")]
    res <- merge(aux, res, by.y = c("outcome", "age"), by.x = c("outcome", "sex"))
    res$pop <- sum(sofia_pop_poly$total_adults)
    res$age <- "adults"
    
  }
  
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
  res_no2_city_disease <- rbind(res_no2_city_disease, res)
  
}

### Export the results
write_csv(res_no2_city_mort, "Baseline_BOD/Results/RESULTS_NO2_mortality.csv")
write_csv(res_no2_city_disease, "Baseline_BOD/Results/RESULTS_NO2_disease.csv")


# ### Plot the results
# # Mortality
res_no2_city_mort <- read_csv2("Baseline_BOD/Results/RESULTS_NO2_mortality.csv")
res_no2_city_mort <- res_no2_city_mort[res_no2_city_mort$sex == "total", ]
res_no2_city_mort <- res_no2_city_mort[!res_no2_city_mort$age == "adults", ]
res_no2_city_mort <- res_no2_city_mort[-c(11:20), ]

res_no2_city_mort$outcome <- factor(
  res_no2_city_mort$outcome,
  levels = c("nat_mort", "cvd_mort", "resp_mort", "lungc_mort"), # Adjust these levels to your dataset
  labels = c("Natural-cause mortality", "Cardiovascular mortality", "Respiratory mortality", "Lung cancer mortality") # Rename outcomes
)

theme_set(theme_bw())
ggplot(res_no2_city_mort, aes(x=age, y=mean)) + 
  geom_bar(stat = "identity", fill = "mediumpurple4") +
  facet_wrap(. ~ outcome, scales = "free") +
  coord_flip() + labs(title=expression("NO"[2]*" attributable mortality by age group"),y ="Attributable deaths (n)", x = "") +
  geom_errorbar(aes(ymin = pct2.5, ymax = pct97.5), width = 0.2, color = "black")
ggsave("Baseline_BOD/Results/RESULTS_NO2_age.png")


# # Disease outcomes
# res_no2_city_disease <- read_csv2("Baseline_BOD/Results/RESULTS_NO2_disease.csv")
# p4 <- ggplot(res_no2_city_disease, aes(x=reorder(outcome, mean), y=mean)) + 
#             geom_bar(stat = "identity", fill = "mediumpurple4") +
#             coord_flip() + labs(title="NO2 attributable morbidity",y ="Attributable cases (n)", x = "") +
#             scale_x_discrete(labels=c("copd" = "COPD", "stroke" = "Stroke", "ihd" = "IHD", "diabetes" = "Diabetes (type II)", "asthma_adults" = "Asthma (adults)", "asthma_children" = "Asthma (children)"))
# 
# ### Export plots with PM2.5 and NO2 results
# library(ggpubr)
# ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
# ggsave("Baseline_BOD/Results/RESULTS_airpollution.png")


##########################################
### O3
##########################################

# Respiratory mortality
allages <- names(sofia_pop_poly)[c(15:24, 30:39)] # define age groups to analyze
res_o3_resp_mort <- hia_sim(pollutant = "o3",
                             outcome = "resp_mort",
                             adjustment = "no",
                             nsim = 500, 
                             seed = 666)
sum(res_o3_resp_mort[res_o3_resp_mort$age == "overall", ]$mean) # 44 respiratory deaths 

# Save the results
# Read census tracts
sofia_poly <- st_read("Baseline_BOD/Clean_data/sofia_pop_poly.geojson")[-c(7:39)]  # Sofia polygons 

# Select results data frame
res <- res_o3_resp_mort
  
# Select results by age/sex
res2 <- res[!res$age == "overall", ]
# Merge population by age/sex
pop <- sofia_pop_poly[c(2, 10:41)] %>% pivot_longer(cols = male_0_4:total_children, names_to = "age", values_to = "pop")
res2 <- merge(res2, pop, by = c("Poly_ID", "age"))
poptot <- pop[-1] %>% group_by(age) %>% summarise(poptot = sum(pop, na.rm = T))
  
# Select overall results
res <- res[res$age == "overall", ]
  
# Merge to population (adults, children) and health outcome rates
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
    
# Export results
res <- res[c("Poly_ID", "outcome", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "att_rate", "att_rate_lwr", "att_rate_upr")]
# as csv
write.csv(res, "Baseline_BOD/Results/res_o3_resp_mort.csv")
# as Geojson
res <- merge(sofia_poly, res, by.x = "Poly_ID", by.y = "Poly_ID")
st_write(res, "Baseline_BOD/Results/res_o3_resp_mort.geojson")


# Overall results for the city
# Mortality outcomes
res_o3_city_mort <- NULL # data frame to store results

# Select results data frame
res <- res_o3_resp_mort
  
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
res_o3_city_mort <- rbind(res_o3_city_mort, res)
  
# Mortality outcomes -- by sex
# Select results data frame
res <- res_o3_resp_mort
  
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
res_o3_city_mort <- rbind(res_o3_city_mort, res)
  
# Mortality outcomes -- by age
# Select results data frame
res <- res_o3_resp_mort
  
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
res_o3_city_mort <- rbind(res_o3_city_mort, res)
  
### Export the results
write_csv(res_o3_city_mort, "Baseline_BOD/Results/RESULTS_O3_mortality.csv")


# ### Plot the results
# # Mortality
res_o3_city_mort <- read_csv2("Baseline_BOD/Results/RESULTS_O3_mortality.csv")
res_o3_city_mort <- res_o3_city_mort[res_o3_city_mort$sex == "total", ]
res_o3_city_mort <- res_o3_city_mort[!res_o3_city_mort$age == "adults", ]
res_o3_city_mort <- res_o3_city_mort[-c(11:20), ]

res_o3_city_mort$outcome <- factor(
  res_o3_city_mort$outcome,
  levels = c("nat_mort", "cvd_mort", "resp_mort", "lungc_mort"), # Adjust these levels to your dataset
  labels = c("Natural-cause mortality", "Cardiovascular mortality", "Respiratory mortality", "Lung cancer mortality") # Rename outcomes
)

theme_set(theme_bw())
ggplot(res_o3_city_mort, aes(x=age, y=mean)) + 
  geom_bar(stat = "identity", fill = "darkorange") +
  facet_wrap(. ~ outcome, scales = "free") +
  coord_flip() + labs(title=expression("O"[3]*" attributable mortality by age group"),y ="Attributable deaths (n)", x = "") +
  geom_errorbar(aes(ymin = pct2.5, ymax = pct97.5), width = 0.2, color = "black")
ggsave("Baseline_BOD/Results/RESULTS_O3_age.png")



