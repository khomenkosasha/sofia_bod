
rm(list = ls()) ### Clear environment

#-------------------------
# Load packages
#-------------------------

library(tidyverse)
library(dplyr)
library(sf)
library(ggplot2)
library(dlnm)
library(mixmeta)

#-------------------------
# Load data
#-------------------------

### Population data
sofia_pop_poly <- st_read("Baseline_BOD/Clean_data/sofia_pop_poly.geojson")
sofia_pop_poly$geometry <- NULL
#sofia_pop_ct$pop_density <- NULL
names(sofia_pop_poly)[24] <- "male_70+"
names(sofia_pop_poly)[39] <- "female_70+"

# Calculate the age groups: 20–44, 45–64, 65–74, 75–84, and ≥85 years (to match ERFs)
# Proportions for 70+: 70-74 (0.427829022 M, 0.386477554 F), 75-84 (0.371877916 M, 0.377592949 F), 85+ (0.200293062 M, 0.235929497 F)
# 20–44
sofia_pop_poly$female_20_44 <- rowSums(sofia_pop_poly[29:33])
sofia_pop_poly$male_20_44 <- rowSums(sofia_pop_poly[14:18])
# 45–64
sofia_pop_poly$female_45_64 <- rowSums(sofia_pop_poly[34:37])
sofia_pop_poly$male_45_64 <- rowSums(sofia_pop_poly[19:22])
# 65-74
sofia_pop_poly$female_65_74 <- sofia_pop_poly$female_65_69 + sofia_pop_poly$`female_70+` * 0.386477554
sofia_pop_poly$male_65_74 <- sofia_pop_poly$male_65_69 + sofia_pop_poly$`male_70+` * 0.427829022
# 75-84
sofia_pop_poly$female_75_84 <- sofia_pop_poly$`female_70+` * 0.377592949
sofia_pop_poly$male_75_84 <- sofia_pop_poly$`male_70+` * 0.371877916
# 85+ 
sofia_pop_poly$`female_85+` <- sofia_pop_poly$`female_70+` * 0.235929497
sofia_pop_poly$`male_85+` <- sofia_pop_poly$`male_70+` * 0.200293062

# Select relevant age groups
sofia_pop_poly <- sofia_pop_poly[, -c(10:39)] # 945,580 adults aged > 20

### All-cause mortality (annual)
sofia_mort <- read_csv2("Baseline_BOD/Health/Sofia_allcause_mortality.csv")
sofia_mort$rate <- as.numeric(sofia_mort$count/ sofia_mort$population * 100000)

### All-cause mortality (weekly) --> to estimate daily mortality for summer months (June-August)
sofia_weekly_mort <- read_csv("Baseline_BOD/Health/Weekly_mortality.csv")[5:9]

### Summer daily temperature
sofia_summertemp_poly <- read_csv("Baseline_BOD/Clean_data/sofia_tempsummer_poly.csv", col_types = cols(Poly_ID = col_character()))

### Summer UHI
sofia_uhi_poly <- st_read("Baseline_BOD/Clean_data/sofia_uhi_poly.geojson")
sofia_uhi_poly$geometry <- NULL

# ERA5 series data for temperature adjustment
ERA5_temp <- read_csv("Baseline_BOD/Temperature/additional_data/era5series.csv")
ERA5_temp <- ERA5_temp %>% filter(URAU_CODE == "BG001C")


#--------------------------------------------------
# Estimate daily mortality by age group and sex
# For summer period (June-August)
#--------------------------------------------------

### Filter weekly mortality data
# Filter sex
sofia_weekly_mort <- sofia_weekly_mort[sofia_weekly_mort$sex %in% c("F", "M"), ]

# Recode age groups: 20–44, 45–64, 65–74, 75–84, and ≥85 years (to match ERFs)
sofia_weekly_mort$age <- ifelse(sofia_weekly_mort$age %in% c("Y20-24", "Y25-29", "Y30-34", "Y35-39", "Y40-44"), "20_44", sofia_weekly_mort$age)
sofia_weekly_mort$age <- ifelse(sofia_weekly_mort$age %in% c("Y45-49", "Y50-54", "Y55-59", "Y60-64"), "45_64", sofia_weekly_mort$age)
sofia_weekly_mort$age <- ifelse(sofia_weekly_mort$age %in% c("Y65-69", "Y70-74"), "65_74", sofia_weekly_mort$age)
sofia_weekly_mort$age <- ifelse(sofia_weekly_mort$age %in% c("Y75-79", "Y80-84"), "75_84", sofia_weekly_mort$age)
sofia_weekly_mort$age <- ifelse(sofia_weekly_mort$age %in% c("Y85-89", "Y_GE90"), "85+", sofia_weekly_mort$age)
sofia_weekly_mort <- sofia_weekly_mort %>% group_by(sex, age, TIME_PERIOD) %>% summarise(count = sum(OBS_VALUE, na.rm = T))

### Calculate total annual deaths by age group
aux <- sofia_weekly_mort %>% group_by(age, sex) %>% summarise(tot = sum(count))

### Calculate proportion of weekly summer deaths
# Filter summer weeks (22-35)
sofia_weekly_mort$TIME_PERIOD <- as.numeric(substr(sofia_weekly_mort$TIME_PERIOD, 7, 8))
sofia_weekly_mort <- sofia_weekly_mort[sofia_weekly_mort$TIME_PERIOD >= 22 & sofia_weekly_mort$TIME_PERIOD <= 35, ]

# Merge to total counts by age group and sex
sofia_weekly_mort <- sofia_weekly_mort[!sofia_weekly_mort$age == "TOTAL", ]
sofia_weekly_mort <- merge(sofia_weekly_mort, aux, by = c("age", "sex"))
sum(sofia_weekly_mort$count)/sum(unique(sofia_weekly_mort$tot)) # proportion of summer deaths among adults
rm(aux)

# Calculate proportions
sofia_weekly_mort$propmort <- sofia_weekly_mort$count / sofia_weekly_mort$tot * 100
sofia_weekly_mort$sex <- ifelse(sofia_weekly_mort$sex == "F", "female", "male")
sofia_weekly_mort[4:5] <- NULL
sofia_weekly_mort <- sofia_weekly_mort %>% pivot_wider(names_from = c(sex, age), values_from = propmort)

### Assign weekly values to each day June-August
sofia_daily_mort <- data.frame(datetime = unique(sofia_summertemp_poly$datetime))
sofia_daily_mort$week <- as.numeric(strftime(sofia_daily_mort$datetime, format = "%V"))
sofia_daily_mort <- merge(sofia_daily_mort, sofia_weekly_mort, by.x = "week", by.y = "TIME_PERIOD")
sofia_daily_mort[3:12] <- sofia_daily_mort[3:12]/7 # Assuming equal death rate in each week


#--------------------------------------------------
# Temperature adjustment
# Based on ERA 5 data
#--------------------------------------------------

################################################################################
# FUNCTION TO APPLY THE CALIBRATION PROCEDURE DESCRIBED IN:
#   Hempel et al. A trend-preserving bias correction - the ISI-MIP approach.
#      Earth Syst Dynam, 2013;4:219-236.
#
# ARGUMENTS:
#   - obs: DATE + SINGLE OBSERVED HISTORICAL SERIES
#   - mod: DATE + OPTIONALLY MULTIPLE MODELLED SERIES
#   - add: LOGICAL. IF TRUE, THE ADDITIVE CORRECTION IS APPLIED
#   - mult: LOGICAL. IF TRUE, THE MULTIPLICATIVE CORRECTION IS APPLIED
#   - output: RETURN THE SERIES ("series") OR THE PARAMETERS ("correction")
#
# Author: Antonio Gasparrini - GNU General Public License (version 3)
# Update: 24 April 2023
################################################################################

#-------------------------------------------------------------------------------

fhempel <- function(obs ,mod, add=TRUE, mult=TRUE, output="series") {
  #
  # CHECK THE OUTPUT
  output <- match.arg(output,c("series","correction"))
  #
  # CHECK CONSISTENCY
  if(ncol(obs)!=2) stop("obs must have two variables (date and series")
  if(ncol(mod)<2) stop("mod must have at least two variables (date and series")
  if(!any(class(obs[[1]])%in%"Date") | !any(class(mod[[1]])%in%"Date"))
    stop("Class of first variable of 'obs' and 'mod' must be 'Date'")
  #
  # APPLY TO EACH MODELLED SERIES
  out <- lapply(seq(ncol(mod)-1),function(j) {
    #    
    # IDENTIFY PERIOD WITH NO MISSING
    ind <- obs[[1]][obs[[1]] %in% mod[[1]]]
    notna <- complete.cases(obs[obs[[1]]%in%ind,2],mod[mod[[1]]%in%ind,j+1])
    indobs <- seq(nrow(obs))[obs[[1]]%in%ind][notna]
    indmod <- seq(nrow(mod))[mod[[1]]%in%ind][notna]
    if(length(indobs)==0) stop("no common non-missing days in 'obs' and 'mod'")
    #
    # EXTRACT MONTH AND YEAR
    month <- as.numeric(format(obs[[1]][indobs],format="%m"))
    # if(length(unique(month))!=12) 
    #   stop("some months are not reprensented in the overlapping period")
    year <- as.numeric(format(obs[[1]][indobs],format="%Y"))
    monthyear <- factor(paste(year,month,sep="-"))   
    #
    # COMPUTE ADDITIVE CORRECTION
    mavgobs <- tapply(obs[[2]][indobs],month,mean,na.rm=T)
    mavgmod <- tapply(mod[[j+1]][indmod],month,mean,na.rm=T)
    C <- mavgobs - mavgmod
    names(C) <- unique(monthyear)
    if(!add) C[] <- 0
    #
    # RESIDUALS FROM MONTHLY/YEAR AVERAGES, THEN MULTIPLICATIVE CORRECTION
    myavgobs <- tapply(obs[[2]][indobs],monthyear,mean,na.rm=T)
    myavgmod <- tapply(mod[[j+1]][indmod],monthyear,mean,na.rm=T)
    resobs <- obs[[2]][indobs] - myavgobs[monthyear]
    resmod <- mod[[j+1]][indmod] - myavgmod[monthyear]
    B <- sapply(6:8,
                function(m) coef(lm(sort(resobs[month==m])~0+sort(resmod[month==m]))))
    if(!mult) B[] <- 1 
    names(B) <- unique(monthyear)
    #
    # RETURN CORRECTION IF STATED
    if(output=="correction") return(matrix(c(as.numeric(C), as.numeric(B)),ncol=2, dimnames = list(unique(month), c("add", "mult"))))
    #
    # # OBTAIN DAY, MONTH AND YEAR FROM DATE (ORIGINAL SERIES)
    # day <- as.numeric(format(mod[[1]],format="%d"))
    # month <- as.numeric(format(mod[[1]],format="%m"))
    # year <- as.numeric(format(mod[[1]],format="%Y"))
    # monthyear <- factor(paste(year,month,sep="-"))
    #
    # # DERIVE THE QUANTITIES TO REMOVE DISCONTINUITIES FROM CORRECTION
    # # (SEE Hempel et al, PAGE 228)
    # nday <- tapply(mod[[j+1]],monthyear,length)
    # d <- (day-1)/(nday[monthyear]-1)-0.5
    # dm <- 0.5*(abs(d)-d)
    # d0 <- 1-abs(d)
    # dp <- 0.5*(abs(d)+d)    
    # #    
    # # DERIVE THE CORRECTIONS ACCOUNTING FOR DISCONTINUITIES
    # Cm <- C[c(3,1:2)] ; Bm <- B[c(3,1:2)]
    # Cp <- C[c(2:3,1)] ; Bp <- B[c(2:3,1)]
    # C <- dm*Cm[month] + d0*C[month] + dp*Cp[month]
    # B <- dm*Bm[month] + d0*B[month] + dp*Bp[month]
    #    
    # DERIVE THE CORRECTED SERIES
    myavgout <- tapply(mod[[j+1]],monthyear,mean,na.rm=T)
    resout <- mod[[j+1]] - myavgout[monthyear]
    #    
    # RETURN CORRECTED SERIES IF STATED
    return(myavgout[monthyear]+C[monthyear]+resout*B[monthyear])
  })
  #
  # RETURN CORRECTION IF STATED
  if(output=="correction") {
    
    names(out) <- names(mod)[-1]
    return(out)
  }
  #
  # ADD DATE AND RENAME
  out <- cbind(mod[,1],as.data.frame(do.call(cbind,out)))
  rownames(out) <- NULL
  #dimnames(out) <- dimnames(mod)
  #  
  # RETURN
  return(out)
}


#-------------------------------------------------------------------------------

### Apply function --> for each polygon
poly <- unique(sofia_summertemp_poly$Poly_ID)
ERA5_temp[1] <- NULL

sofia_summertemp_corrected <- NULL # data frame to store results

for (i in 1:length(poly)) {
  
  # filter for each polygon
  aux <- sofia_summertemp_poly[sofia_summertemp_poly$Poly_ID == poly[i], ]
  aux <- aux[2:3]
  
  # Apply function for correction
  corrected_temp <- fhempel(obs = ERA5_temp, mod = aux, add=TRUE, mult=TRUE, output="series")
  corrected_temp$Poly_ID <- poly[i]
  corrected_temp <- corrected_temp[c("Poly_ID", "datetime", "V1")]
  names(corrected_temp)[3] <- "temp"
  
  # save results
  sofia_summertemp_corrected <- rbind(sofia_summertemp_corrected, corrected_temp)
  
  # print progress
  print(paste0("Done ", i, "/4969"))
  
 
}

# Summary of corrected temperature values
summary(sofia_summertemp_poly$temp) # from 11.24 until 34.05
summary(sofia_summertemp_corrected$temp) # from 11.07 until 28.14
rm(sofia_summertemp_poly)


#--------------------------------------------------
# Extract ERF for temperature
# City specific for Sofia
#--------------------------------------------------

# Load data 
coefs <- read_csv("Baseline_BOD/Temperature/coefs.csv")
coefs_simu <- read_csv("Baseline_BOD/Temperature/coef_simu.csv")
vcov <- read_csv("Baseline_BOD/Temperature/vcov.csv")
tmean_distribution <- read_csv("Baseline_BOD/Temperature/tmean_distribution.csv")

# Select the city of interest: BG001C
coefs <- coefs %>% filter(URAU_CODE == "BG001C")
coefs_simu <- coefs_simu %>% filter(URAU_CODE == "BG001C")
vcov <- vcov %>% filter(URAU_CODE == "BG001C")
tmean_distribution <- tmean_distribution %>% filter(URAU_CODE == "BG001C")


#-------------------------------------------------------------------------------
### Function to estimate MMT

findmin <- function(basis,model=NULL,coef=NULL,vcov=NULL,at=NULL,from=NULL,
                    to=NULL,by=NULL,sim=FALSE,nsim=5000) {
  #
  ################################################################################
  #   ARGUMENTS:
  #   - basis: A SPLINE OR OTHER BASIS FOR AN EXPOSURE x CREATED BY DLNM FUNCTION 
  #            CROSSBASIS OR ONEBASIS
  #   - model: THE FITTED MODEL
  #   - coef AND vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
  #
  #   - at: A NUMERIC VECTOR OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT
  #   OR 
  #   - from, to: RANGE OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT.
  #   - by: INCREMENT OF THE SEQUENCES x VALUES OVER WHICH THE MINIMUM IS SOUGHT
  # 
  #   - sim: IF BOOTSTRAP SIMULATION SAMPLES SHOULD BE RETURNED
  #   - nsim: NUMBER OF SIMULATION SAMPLES
  ################################################################################
  
  
  ################################################################################
  # CREATE THE BASIS AND EXTRACT COEF-VCOV
  #
  # CHECK AND DEFINE BASIS  
  if(!any(class(basis)%in%c("crossbasis","onebasis")))
    stop("the first argument must be an object of class 'crossbasis' or 'onebasis'")
  #
  # INFO
  one <- any(class(basis)%in%c("onebasis"))
  attr <- attributes(basis)
  range <- attr(basis,"range")
  if(is.null(by)) by <- 0.1
  lag <- if(one) c(0,0) else cb=attr(basis,"lag")
  if(is.null(model)&&(is.null(coef)||is.null(vcov)))
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  name <- deparse(substitute(basis))
  cond <- if(one) paste(name,"[[:print:]]*b[0-9]{1,2}",sep="") else 
    paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
  #
  # SET COEF, VCOV CLASS AND LINK
  if(!is.null(model)) {
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- dlnm:::getlink(model,model.class)
  } else model.class <- NA
  #
  # CHECK
  if(length(coef)!=ncol(basis) || length(coef)!=dim(vcov)[1] ||
     any(is.na(coef)) || any(is.na(vcov)))
    stop("model or coef/vcov not consistent with basis")
  #
  # DEFINE at
  at <- dlnm:::mkat(at,from,to,by,range,lag,bylag=1)
  predvar <- if(is.matrix(at)) rownames(at) else at
  predlag <- dlnm:::seqlag(lag,by=1)
  #
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON TYPE)
  type <- if(one) "one" else "cb"
  Xpred <- dlnm:::mkXpred(type,basis,at,predvar,predlag,cen=NULL)
  Xpredall <- 0
  for(i in seq(length(predlag))) {
    ind <- seq(length(predvar))+length(predvar)*(i-1)
    Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]
  }
  #  
  ################################################################################
  # FIND THE MINIMUM
  #
  pred <- drop(Xpredall%*%coef)
  ind <- which.min(pred)
  min <- predvar[ind]
  #
  ################################################################################
  # SIMULATIONS
  #
  if(sim) {
    # SIMULATE COEFFICIENTS
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # COMPUTE MINIMUM
    minsim <- apply(coefsim,2,function(coefi) {
      pred <- drop(Xpredall%*%coefi)
      ind <- which.min(pred)
      return(predvar[ind])
    })
  }
  #
  ################################################################################
  #
  res <- if(sim) minsim else min
  #
  return(res)
}

#-------------------------------------------------------------------------------

### Extract RR for each temperature and age group
agegroups <- unique(coefs$agegroup)

temp_erfs_sofia <- NULL # data frame to store the results

for (i in 1:length(agegroups)) {
  
  # Create basis
  ov_basis <- onebasis(as.numeric(tmean_distribution[2:120]), fun = "bs", degree = 2, knots = quantile(as.numeric(tmean_distribution[2:120]), c(.1, .75, .9)))
  # Extract the MMT
  mmt <- findmin(ov_basis, coef = as.numeric(coefs[i, 3:7]), vcov = xpndMat(vcov[i, 3:17]), from = as.numeric(tmean_distribution[12]), to = as.numeric(tmean_distribution[110]))
  mmt_CI <- quantile(findmin(ov_basis, coef = as.numeric(coefs[i, 3:7]), vcov = xpndMat(vcov[i, 3:17]), from = as.numeric(tmean_distribution[12]), to = as.numeric(tmean_distribution[110]), sim=TRUE), c(0.025, 0.975))
  
  # Get the RR from the function
  cp <- crosspred(basis = ov_basis, coef = as.numeric(coefs[i, 3:7]), vcov = xpndMat(vcov[i, 3:17]), model.link = "log", cen = mmt, at =  seq(11, 29, by = 0.1))
  
  # Create data frame with RR values
  erfs_gasp <- data.frame(age = agegroups[i],
                          temp = round(as.numeric(cp$predvar), digits = 1),  
                          RR_FIT = cp$allRRfit,                            
                          RR_LOW = cp$allRRlow,                             
                          RR_HIGH = cp$allRRhigh,                           
                          ERF_se = cp$allse,
                          MMT = mmt, 
                          MMT_LOW = mmt_CI[1],
                          MMT_HIGH = mmt_CI[2])
  
  # Save the results
  temp_erfs_sofia <- rbind(temp_erfs_sofia, erfs_gasp)
  
  # Print progress
  print(paste0("Done ", i, " age group"))
  
  
}

# Rename variables
temp_erfs_sofia <- temp_erfs_sofia %>% mutate(age = recode(age, "20-44" = "20_44", "45-64" = "45_64", "65-74" = "65_74", "75-84" = "75_84"))

# Get maximum temperature
maxtemp <- max(as.numeric(tmean_distribution[2:120])) # 28.25 ºC

# Remove not needed datasets
rm(aux)
rm(coefs)
rm(coefs_simu)
rm(corrected_temp)
rm(cp)
rm(ERA5_temp)
rm(erfs_gasp)
rm(ov_basis)
rm(tmean_distribution)
rm(vcov)



#--------------------------------------------------
# First, we need to create some vectors and data 
# frames with parameters that will be used later 
# in the Monte Carlo analysis
#--------------------------------------------------

### We extract the age and sex labels
allages <- names(sofia_pop_poly)[10:19] # these need to be adjusted for each analysis based on the outcome
allages

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
                    pollutant = "uhi_summer", 
                    outcome = "allc_mort") {                 
  
  ### Define age and sex for the analysis: 
  agesub <- sub(".*?_(.*)", "\\1", age)
  sexsub <- sub("\\_.*", "", age)
  
  ### Extract the MMT: 
  MMT <- unique((temp_erfs_sofia[temp_erfs_sofia$age == agesub, "MMT"]))
  #MMT_LOW <-unique((temp_erfs_sofia[temp_erfs_sofia$age == agesub, "MMT_LOW"]))
  #MMT_HIGH <- unique((temp_erfs_sofia[temp_erfs_sofia$age == agesub, "MMT_HIGH"]))
  
  # simulate a MMT value assuming MMT is normally distributed:
  # semmt <- (MMT_HIGH - MMT_LOW) / (2 * qnorm(0.975))
  # mmtsim <- rnorm(n = 1, mean = MMT, sd = semmt)
  # MMT <- mmtsim
  
  ### Get daily temperature
  temp_city <- sofia_summertemp_corrected
  
  # Simulate temperature error (based on standard error from Urbclim)
  rmse_urbclim <- 1.4
  se <- rmse_urbclim * qnorm(0.975)
  gids <- length(unique(temp_city$Poly_ID))
  # simulate error at each grid:
  temperror_urbclim <- rnorm(n = (dim(temp_city)[1])/gids, mean = 0, sd = se)
  temperror_urbclim <- as.data.frame(temperror_urbclim)
  temperror_urbclim$datetime <- seq(as.Date("2017-06-01"),as.Date("2017-08-31"),by = 1)
  temp_city <- left_join(temp_city, temperror_urbclim, by = "datetime")
  temp_city <- temp_city %>% arrange(Poly_ID, datetime)
  # update the value of temperature by adding the error to the point estimate:
  temp_city$temp <- temp_city$temp + temp_city$temperror_urbclim
  
  # ### Model to adjust temp (ERA5 - Urbclim)
  # se <- regression$se
  # alfa <- regression$alfa
  # beta <- regression$beta
  # 
  # # adjust temperature (era5)  
  # temp_city$temp <- alfa + beta * temp_city$temp
  
  # # Simulate temperature value based on ERA5 error: 
  # temperror_era5 <- rnorm(n = dim(temp_city)[1], mean = 0, sd = se)
  # # update the value of temperature by adding the error to the point estimate:
  # temp_city$temp <- temp_city$temp + temperror_era5
  # 
  # ### Set maximum temperature (predvar, defined by the ERF)
  # max_temp <- as.numeric(max_predvar$max)
  temp_city$temp <- ifelse(temp_city$temp > maxtemp , maxtemp, temp_city$temp)
  
  # Round temperature for merging with ERFs
  temp_city$temp <- round(temp_city$temp, digits = 1)
  
  ### Extract the daily UHI
  df_city <- sofia_uhi_poly[c("Poly_ID", pollutant)] 
  names(df_city)[2] <- "poll"

  # Simulate UHI error based on standard error from UrbClim: 
  rmse_uhi <- 1
  se <- rmse_uhi * qnorm(0.975)
  uhierror_urbclim <- rnorm(n = (dim(df_city)[1]), mean = 0, sd = se)
  # update the value of temperature by adding the error to the point estimate:
  df_city$poll <- df_city$poll + uhierror_urbclim
  df_city$poll <- ifelse(df_city$poll < 0, 0, df_city$poll)

  ### Create data frame with counterfactual exposure values: 
  temp_count <- inner_join(temp_city, df_city, by = "Poly_ID")
  temp_count$temp <- temp_count$temp - temp_count$poll

  # Round temperature for merging with ERFs
  temp_count$temp <- round(temp_count$temp, digits = 1) 
  
  # Filter temperatures above the MMT
  temp_city <- temp_city %>% filter(temp > MMT)
  temp_count <- temp_count %>% filter(temp > MMT)
  
  ### Extract RR values: 
  # subset age group
  erfs_gasp_age <- (temp_erfs_sofia[temp_erfs_sofia$age == agesub, c("temp", "RR_FIT", "RR_LOW", "RR_HIGH")])
  
  # ASSIGN A RR TO EACH OF THE DAILY TEMPERATURES
  erfs_temp <- inner_join(temp_city, erfs_gasp_age, by=c("temp")) 
  
  # ASSIGN A RR TO EACH OF THE DAILY TEMPERATURES - COUNTERFACTUAL
  erfs_temp_count <- inner_join(temp_count, erfs_gasp_age, by=c("temp")) 
  
  # simulate a RR value assuming log(RR) is normal for baseline:
  selogrr <- as.vector(log(erfs_temp$RR_HIGH / erfs_temp$RR_LOW) / (2 * qnorm(0.975)))
  logrr <- as.vector(log(erfs_temp$RR_FIT))
  erfs_temp$logrrsim <- rnorm(length(selogrr), mean = logrr, sd = selogrr)
  erfs_temp$rr <- exp(erfs_temp$logrrsim)
  
  # simulate a RR value assuming log(RR) is normal for counterfactual:
  selogrr_count <- as.vector(log(erfs_temp_count$RR_HIGH / erfs_temp_count$RR_LOW) / (2 * qnorm(0.975)))
  logrr_count <- as.vector(log(erfs_temp_count$RR_FIT))
  erfs_temp_count$logrrsim <- rnorm(length(selogrr_count), mean = logrr_count, sd = selogrr_count)
  erfs_temp_count$rr <- exp(erfs_temp_count$logrrsim)

  ### We retrieve the outcome rate based on the selected outcome and age group
  mratecityage <- as.numeric(sofia_mort[sofia_mort$age == agesub & sofia_mort$sex == sexsub, "rate"])

  ### We retrieve the city population for each neighborhood:
  popcity <- sofia_pop_poly[, c("Poly_ID", age)]
  names(popcity)[2] <- "pop"
  
  ### We retrieve the daily mortality proportions for the age group and sex
  #agesex <- paste0(agesub, "_", sexsub)
  dailymort <- sofia_daily_mort[, c("datetime", age)]
  names(dailymort)[2] <- "prop"
  
  # Then, we compute the daily mortality 
  popcity$outcome <- popcity$pop * mratecityage / 100000
  popcity <- cross_join(popcity, dailymort)
  popcity$dailymort <- popcity$outcome * (popcity$prop / 100)
  
  ### Next, we merge the temperature and mortality data
  # Baseline
  popcity_base <- inner_join(popcity, erfs_temp[c("datetime", "Poly_ID", "temp", "rr")], by = c("Poly_ID", "datetime"))
  # Counterfactual
  popcity_counter <- inner_join(popcity, erfs_temp_count[c("datetime", "Poly_ID", "temp", "rr")], by = c("Poly_ID", "datetime"))
  
  ### Now that all values are simulated, we do the HIA: calculate the PAF and attributable deaths
  # Baseline
  popcity_base$PAF <- with(popcity_base, (rr - 1) / rr)
  popcity_base$attmort <- with(popcity_base, dailymort * PAF)
  popcity_base <- popcity_base %>% group_by(Poly_ID) %>% summarise(attmort = sum(attmort, na.rm = T))
  # Counterfactual
  popcity_counter$PAF <- with(popcity_counter, (rr - 1) / rr)
  popcity_counter$attmort_counter <- with(popcity_counter, dailymort * PAF)
  popcity_counter <- popcity_counter %>% group_by(Poly_ID) %>% summarise(attmort_counter = sum(attmort_counter, na.rm = T))
  
  ### We calculate mortality due to temperature exposure in each neighborhood
  popcity_base <- left_join(popcity_base, popcity_counter, by = "Poly_ID")
  popcity_base$attmort_counter <- ifelse(is.na(popcity_base$attmort_counter), 0, popcity_base$attmort_counter)
  popcity_base$heatmort <- popcity_base$attmort - popcity_base$attmort_counter
  
  # Save results
  #res <- data.frame(ct_id = unique(popcity$id))
  res <- popcity_base[c("Poly_ID", "heatmort")]
  res$heatmort <- ifelse(is.na(res$heatmort), 0, res$heatmort)
  #res$age <- age
  
  return(res) ### We ask the function to return the health outcome estimate by neighborhood
}

### Example --> you will see that if you run this several times you get a different value each time (because of the uncertainties)
hia1 <- hia_age(age = allages[10], pollutant = "uhi_summer", outcome = "allc_mort")
hia1
sum(hia1$heatmort)



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
                        pollutant = "uhi_summer", 
                        outcome = "allc_mort",
                        nsim = 100) {                  # number of simulations
  
  # We use the replication function to replicate the previously defined function n times (to get a sample of values):
  samp <- replicate(n = nsim,
                    expr = hia_age(age = age, 
                                   pollutant = pollutant, 
                                   outcome = outcome), simplify = FALSE)
                    
  # We use the sample of values to get point (mean) and CI estimate (2.5 and 97.5 percentiles)
  summ <- bind_rows(samp, .id = "id") %>% 
     group_by(Poly_ID) %>%
     summarise(mean = mean(heatmort, na.rm = TRUE), 
               median = median(heatmort, na.rm = TRUE), 
               pct2.5 = quantile(heatmort, probs = 2.5/100, na.rm = TRUE), 
               pct97.5 = quantile(heatmort, probs = 97.5/100, na.rm = TRUE))
  
  # We ask the function to return a list with results
  res <- list(sample = samp, estimate = summ)
  #res <- samp
  return(res)
}

# example with 500 simulations --> you can try different numbers of simulations and check how long it takes to calculate
t0 <- Sys.time()
simhia1 <- hia_age_sim(age = allages[10],
                       pollutant = "uhi_summer", 
                       outcome = "allc_mort",
                       nsim = 100)
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

hia_sim <- function(pollutant = "uhi_summer", 
                    outcome = "allc_mort",
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
    #samp[[i]] <- aux
  }
  
  # compute point and CI for all ages:
  samp <- bind_rows(samp, .id = "id")
  samp$id <- rep(1:nsim, times = nages, each = 4969)
  
  # for the total health outcome
  overallsamp <- samp %>% group_by(id, Poly_ID) %>% summarise(sum = sum(heatmort, na.rm = TRUE))
  estoverall <- overallsamp %>% group_by(Poly_ID) %>% summarise(mean = mean(sum, na.rm = TRUE), median = median(sum, na.rm = TRUE), pct2.5 = quantile(sum, probs = 2.5 / 100, na.rm = TRUE), pct97.5 = quantile(sum, probs = 97.5 / 100, na.rm = TRUE))
  #estoverall <- estoverall %>% add_column(age = "overall", .after = "Poly_ID")
  
  # by age group
  #agesamp <- samp %>% group_by(id, Poly_ID, age) %>% summarise(sum = sum(heatmort, na.rm = TRUE))
  #est <- agesamp %>% group_by(Poly_ID, age) %>% summarise(mean = mean(sum, na.rm = TRUE), median = median(sum, na.rm = TRUE), pct2.5 = quantile(sum, probs = 2.5 / 100, na.rm = TRUE), pct97.5 = quantile(sum, probs = 97.5 / 100, na.rm = TRUE))
  
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
hia2 <- hia_sim(pollutant = "uhi_summer", 
                outcome = "allc_mort",
                nsim = 100, 
                seed = 666) ### important always to set the same seed for replicability!
t1 <- Sys.time()
t1 - t0  # --> see the calculation time
hia2




#########################################

# Apply functions 

##########################################

##########################################
### UHI
##########################################

# All-cause mortality --> 45 minutes
allages <- names(sofia_pop_poly)[10:19] # define age groups to analyze
res_uhi_allc_mort <- hia_sim(pollutant = "uhi_summer", 
                             outcome = "allc_mort",
                             nsim = 500, 
                             seed = 666)
sum(res_uhi_allc_mort[res_uhi_allc_mort$age == "overall", ]$mean) # 95 all-cause deaths during summer
sum(res_uhi_allc_mort[res_uhi_allc_mort$age == "overall", ]$pct2.5)
sum(res_uhi_allc_mort[res_uhi_allc_mort$age == "overall", ]$pct97.5)

### Save the results
# Overall results by neighborhood
# Read census tracts
sofia_poly <- st_read("Baseline_BOD/Clean_data/sofia_pop_poly.geojson")[-c(7:39)]  # Sofia polygons 

# Select results data frame
res <- res_uhi_allc_mort

# Select results by age/sex
res2 <- res[!res$age == "overall", ]
# Merge population by age/sex
pop <- sofia_pop_poly[c(2, 10:19)] %>% pivot_longer(cols = female_20_44:`male_85+`, names_to = "age", values_to = "pop")
res2 <- merge(res2, pop, by = c("Poly_ID", "age"))
poptot <- pop[-1] %>% group_by(age) %>% summarise(poptot = sum(pop, na.rm = T))

# Select overall results
res <- res[res$age == "overall", ]

# Merge to population (adults) and health outcome rates
sofia_pop_poly$total_adults <- rowSums(sofia_pop_poly[10:19])
aux <- sofia_pop_poly[c("Poly_ID", "total_adults")]
res <- merge(res, aux, by.x = "Poly_ID", by.y = "Poly_ID")
aux <- sofia_mort[sofia_mort$age %in% c("adults") & sofia_mort$sex == "total", ]$rate
res$rate <- aux 

# Calculate expected cases (annual)
res$exp_cases <- res$total_adults * (res$rate / 100000)
# Calculate expected cases (summer)
res$exp_cases_summer <- res$exp_cases * 0.2564389

# Calculate % of all cases
res$per_cases <- res$mean / res$exp_cases * 100
res$per_cases <- ifelse(is.na(res$per_cases), 0, res$per_cases)
res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
res$per_cases_lwr <- ifelse(is.na(res$per_cases_lwr), 0, res$per_cases_lwr)
res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100
res$per_cases_upr <- ifelse(is.na(res$per_cases_upr), 0, res$per_cases_upr)

# Calculate % of summer cases
res$per_cases_summer <- res$mean / res$exp_cases_summer * 100
res$per_cases_summer <- ifelse(is.na(res$per_cases_summer), 0, res$per_cases_summer)
res$per_cases_lwr_summer <- res$pct2.5 / res$exp_cases_summer * 100
res$per_cases_lwr_summer <- ifelse(is.na(res$per_cases_lwr_summer), 0, res$per_cases_lwr_summer)
res$per_cases_upr_summer <- res$pct97.5 / res$exp_cases_summer * 100
res$per_cases_upr_summer <- ifelse(is.na(res$per_cases_upr_summer), 0, res$per_cases_upr_summer)

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
res <- res[c("Poly_ID", "outcome", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "per_cases_summer", "per_cases_lwr_summer", "per_cases_upr_summer", "att_rate", "att_rate_lwr", "att_rate_upr")]
# as csv
write.csv(res, "Baseline_BOD/Results/res_uhi_allcause_mort.csv")
# as Geojson
res <- merge(sofia_poly, res, by.x = "Poly_ID", by.y = "Poly_ID")
st_write(res, "Baseline_BOD/Results/res_uhi_allcause_mort.geojson")


### Overall results for the city
# Mortality outcomes -- all age groups
# Select results data frame
res <- res_uhi_allc_mort

# Group results for the whole city
# Select overall results
res <- res[res$age == "overall", ]
res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
res$age <- ifelse(res$age == "overall", "adults", res$age)
res$sex <- "total"

# Add adult population and mortality rates by sex
aux <- sofia_mort[sofia_mort$age %in% c("adults"), ]
aux <- aux[c("age", "sex", "rate")]
res <- merge(aux, res, by = c("age", "sex"))
res$pop <- sum(sofia_pop_poly$total_adults)

# Calculate expected cases
res$exp_cases <- res$pop * (res$rate / 100000)
# Calculate expected summer cases
res$exp_cases_summer <- res$exp_cases * 0.2564389

# Calculate % of all cases
res$per_cases <- res$mean / res$exp_cases * 100
res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100

# Calculate % of all summer cases
res$per_cases_summer <- res$mean / res$exp_cases_summer * 100
res$per_cases_lwr_summer <- res$pct2.5 / res$exp_cases_summer * 100
res$per_cases_upr_summer <- res$pct97.5 / res$exp_cases_summer * 100

# Calculate attributable rate
res$att_rate <- res$mean / res$pop * 100000
res$att_rate_lwr <- res$pct2.5 / res$pop * 100000
res$att_rate_upr <- res$pct97.5 / res$pop * 100000

# Order results
restot <- res[c("outcome", "age", "sex", "pop", "rate", "exp_cases", "exp_cases_summer", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "per_cases_summer", "per_cases_lwr_summer", "per_cases_upr_summer", "att_rate", "att_rate_lwr", "att_rate_upr")]


# Mortality outcomes -- by sex
# Select results data frame
res <- res_uhi_allc_mort

# Group results for the whole city by sex
# Select overall results
res <- res[!res$age == "overall", ]
res$age <- sub("\\_.*", "", res$age)
res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
names(res)[2] <- "sex"
res$age <- "adults"

# Add adult population and mortality rates by sex
aux <- sofia_mort[sofia_mort$age %in% c("adults"), ]
aux <- aux[c("age", "sex", "rate")]
res <- merge(aux, res, by = c("age", "sex"))
sofia_pop_poly$adults_female <- sofia_pop_poly$female_20_44 + sofia_pop_poly$female_45_64 + sofia_pop_poly$female_65_74 + sofia_pop_poly$female_75_84 + sofia_pop_poly$`female_85+`
sofia_pop_poly$adults_male <- sofia_pop_poly$male_20_44 + sofia_pop_poly$male_45_64 + sofia_pop_poly$male_65_74 + sofia_pop_poly$male_75_84 + sofia_pop_poly$`male_85+`
res$pop <- c(sum(sofia_pop_poly$adults_female), sum(sofia_pop_poly$adults_male))

# Calculate expected cases
res$exp_cases <- res$pop * (res$rate / 100000)
# Calculate expected summer cases
res$exp_cases_summer <- ifelse(res$sex == "female", res$exp_cases * 0.2573122, res$exp_cases * 0.2555302)   

# Calculate % of all cases
res$per_cases <- res$mean / res$exp_cases * 100
res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100

# Calculate % of all summer cases
res$per_cases_summer <- res$mean / res$exp_cases_summer * 100
res$per_cases_lwr_summer <- res$pct2.5 / res$exp_cases_summer * 100
res$per_cases_upr_summer <- res$pct97.5 / res$exp_cases_summer * 100

# Calculate attributable rate
res$att_rate <- res$mean / res$pop * 100000
res$att_rate_lwr <- res$pct2.5 / res$pop * 100000
res$att_rate_upr <- res$pct97.5 / res$pop * 100000

# Order results
ressex <- res[c("outcome", "age", "sex", "pop", "rate", "exp_cases", "exp_cases_summer", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "per_cases_summer", "per_cases_lwr_summer", "per_cases_upr_summer", "att_rate", "att_rate_lwr", "att_rate_upr")]


# Mortality outcomes -- by age
# Select results data frame
res <- res_uhi_allc_mort

# Group results for the whole city by sex
# Select overall results
res <- res[!res$age == "overall", ]
res$age <- sub(".*?_(.*)", "\\1", res$age)
res <- res %>% group_by(outcome, age) %>% summarise(mean = sum(mean, na.rm = T), pct2.5 = sum(pct2.5, na.rm = T), pct97.5 = sum(pct97.5, na.rm = T))
res$sex <- "total"

# Add adult population and mortality rates by age
aux <- sofia_mort[sofia_mort$sex %in% c("total"), ]
aux <- aux[c("age", "sex", "rate")]
res <- merge(aux, res, by = c("age", "sex"))
res$pop <- c((sum(sofia_pop_poly$male_20_44)+sum(sofia_pop_poly$female_20_44)), (sum(sofia_pop_poly$male_45_64)+sum(sofia_pop_poly$female_45_64)), (sum(sofia_pop_poly$male_65_74)+sum(sofia_pop_poly$female_65_74)), (sum(sofia_pop_poly$male_75_84)+sum(sofia_pop_poly$female_75_84)), (sum(sofia_pop_poly$`male_85+`)+sum(sofia_pop_poly$`female_85+`)))

# Calculate expected cases
res$exp_cases <- res$pop * (res$rate / 100000)
# Calculate expected summer cases
res$exp_cases_summer <- res$exp_cases * c(0.2722603, 0.2631579, 0.2593149, 0.2586207, 0.2456587)  

# Calculate % of all cases
res$per_cases <- res$mean / res$exp_cases * 100
res$per_cases_lwr <- res$pct2.5 / res$exp_cases * 100
res$per_cases_upr <- res$pct97.5 / res$exp_cases * 100

# Calculate % of all summer cases
res$per_cases_summer <- res$mean / res$exp_cases_summer * 100
res$per_cases_lwr_summer <- res$pct2.5 / res$exp_cases_summer * 100
res$per_cases_upr_summer <- res$pct97.5 / res$exp_cases_summer * 100

# Calculate attributable rate
res$att_rate <- res$mean / res$pop * 100000
res$att_rate_lwr <- res$pct2.5 / res$pop * 100000
res$att_rate_upr <- res$pct97.5 / res$pop * 100000

# Order results
resage <- res[c("outcome", "age", "sex", "pop", "rate", "exp_cases", "exp_cases_summer", "mean", "pct2.5", "pct97.5", "per_cases", "per_cases_lwr", "per_cases_upr", "per_cases_summer", "per_cases_lwr_summer", "per_cases_upr_summer", "att_rate", "att_rate_lwr", "att_rate_upr")]

### All results
resall <- rbind(restot, ressex, resage)
write_csv(resall, "Baseline_BOD/Results/RESULTS_UHI_mortality.csv")

### Plot the results (by age group)
resage <- read_csv2("Baseline_BOD/Results/RESULTS_UHI_mortality.csv")
resage <- resage[!resage$age == "adults", ]

resage$outcome <- factor(
  resage$outcome,
  levels = c("allc_mort"), # Adjust these levels to your dataset
  labels = c("All-cause mortality") # Rename outcomes
)

theme_set(theme_bw())
ggplot(resage, aes(x=age, y=mean)) + 
  geom_bar(stat = "identity", fill = "lightcoral") +
  facet_wrap(. ~ outcome, scales = "free") +
  coord_flip() + labs(title="UHI attributable mortality by age group",y ="Attributable deaths (n)", x = "") +
  geom_errorbar(aes(ymin = pct2.5, ymax = pct97.5), width = 0.2, color = "black") +
  scale_x_discrete(labels=c("nage20_44" = "20-44 years", "nage45_64" = "45-64 years", "nage65_74" = "65-74 years", "nage75_84" = "75-84 years", "nage85" = "85+ years"))
ggsave("Baseline_BOD/Results/RESULTS_UHI_age.png")




