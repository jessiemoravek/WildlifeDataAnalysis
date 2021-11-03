library(secr) ## secr package for ML spatial capture recapture models
library(raster) ## always good to have raster loaded for spatial stuff IMHO
library(devtools)

setwd("C:/Users/Owner/Desktop/WildlifeDataAnalysis/WildlifeDataAnalysis/4_SCR_SMR_uSCR") ## Set working directory to where files are stored


## Easiest way to work with secr package is to format data as tab delimited files with required columns
## and read in with 'read.capthist'

fire.capthist <- secr::read.capthist("Fire_capthist_wNAs.txt", "detectors_fire_focal.txt", covnames = "Sex",
                                     detector="count", fmt=c("trapID"))

summary(fire.capthist) ## Summarizes data nicely, we can see there are two 'occasions' - one pre-fire and one post-fire

plot(fire.capthist) ## Visualize captures and detector array

plot(fire.capthist, tracks=TRUE) ## Get a visual idea of movement with 'tracks' argument

fire.capthist.MS <- split(fire.capthist,f=factor(c(1,2)), byoccasion=TRUE)

summary(fire.capthist.MS) ## Now we have two sessions, one pre-fire and one post-fire

plot(fire.capthist.MS, tracks=TRUE) ## Two plots, one for each session

## For secr models, we need a sampling 'mask' which represents the area within which animal activity centers
## MAY be found. To accurately define, we need a good idea of space use. With no auxillary information, we could explore different masks

Fire.mask <- make.mask(traps(fire.capthist), buffer=1000, spacing=150) ## Using buffer of 1 km and spacing of 150 m, based on what we know of deer movement patterns in this area

## By default I am at least interested in sex-specific density estimates, so I will use a 'pmix' sex-proportion as my null model

secr_sex_null <- secr.fit(fire.capthist.MS, model=list( sigma ~ 1 , g0 ~ 1, D~1, pmix~1), hcov="Sex",mask = Fire.mask)

secr_sex_null ## Looks like an estimated 0.26 deer/hectare, and female/male ratio of 0.56/0.44

## I am interested in seeing if the fire had any clear effects on density
secr_sex_Dsess <- secr.fit(fire.capthist.MS, model=list( sigma ~ 1 , g0 ~ 1, D~Session, pmix~1), hcov="Sex",mask = Fire.mask)

secr_sex_Dsess ## Looks like it is estimating about a 20% decrease in density

## What about potential effects of fire on density AND sex ratio 
secr_sex_ratio_Dsess <- secr.fit(fire.capthist.MS, model=list( sigma ~ 1 , g0 ~ 1, D~Session, pmix~Session), hcov="Sex",mask = Fire.mask)

secr_sex_ratio_Dsess ## Looks like in addition to decreasing density, there was an increasing ratio of males to females

## I think that males have larger home ranges/space use than females, so we should probably fit sex-specific sigma parameters

secr_sex_sigma_ratio_Dsess <- secr.fit(fire.capthist.MS, model=list( sigma ~ h2 , g0 ~ 1, D~Session, pmix~Session), hcov="Sex",mask = Fire.mask)

secr_sex_sigma_ratio_Dsess ## Males move farther than females, as expected

## I also think that the probability of detection may vary by sex (perhaps males defecate less often?)
secr_sex_sigma_g0_ratio_Dsess <- secr.fit(fire.capthist.MS, model=list( sigma ~ h2 , g0 ~ h2, D~Session, pmix~Session), hcov="Sex",mask = Fire.mask)

secr_sex_sigma_g0_ratio_Dsess ## Females much more detectable than males

## Finally, I want to test whether movement rates changed before/after the fire
secr_sex_sigmaSess_g0_ratio_Dsess <- secr.fit(fire.capthist.MS, model=list( sigma ~ h2 + Session, g0 ~ h2, D~Session, pmix~Session), hcov="Sex",mask = Fire.mask)

secr_sex_sigmaSess_g0_ratio_Dsess ##

## We have fit models with increasing complexity. Worth performing model selection to see what the relative support of different models are

AIC(secr_sex_null,secr_sex_Dsess,secr_sex_ratio_Dsess,secr_sex_sigma_ratio_Dsess,
    secr_sex_sigma_g0_ratio_Dsess, secr_sex_sigmaSess_g0_ratio_Dsess)

## Strong evidence that the space use and detection probability varies by sex, sex ratio varies before/after
#fire, and density varies 


## Code for a goodness of fit test, takes a long time and uncertain results
# secr.test(secr_sex_sigma_g0_ratio_Dsess, nsim = 20, fit = TRUE) ## Hashed out for now

## Let's ignore fire for now, and work on looking at habitat variable effects on deer density, and mapping those results

## Read in habitat values at each mask 'cell'
mask_habitat_variables <- read.csv("Fire_mask_habitat_vars.csv")

covariates(Fire.mask) <- mask_habitat_variables

## Predict density as a function of distance to water, without looking at sex
secr_D_river <- secr.fit(fire.capthist.MS, model=list( sigma ~ 1 , g0 ~ 1, D~water.dist.clean), mask = Fire.mask)

secr_D_river ## Positive coefficient (but overlapping 0 - not particularly informative or 'significant')

secr_D_river_pred <- predictDsurface(secr_D_river, mask = Fire.mask) ## predictDsurface is a specific predict function in raster package

plot(secr_D_river_pred) ## Missing data on West side of Hopland - should

mean(covariates(secr_D_river_pred)[,14]) ## Average density of 0.266 deer/ha * 945 ha total = 283.5 deer

## To finish, let's visualize the estimated relationship between distance to water and deer

d_water_preds <- predict(secr_D_river, newdata = data.frame(water.dist.clean = seq(0,1020,10))) ## Predict to wide range of distance to water values

plot(seq(0,1020,10), sapply(d_water_preds, "[", "D","estimate"), ylim = c(0.2,0.32), ## Plot it out
     xlab = "Distance from river (m)", ylab = "Density / ha", type = "l")

################################
covariates(Fire.mask) <- mask_habitat_variables

## Predict density as a function of road distance, without looking at sex
secr_D_river <- secr.fit(fire.capthist.MS, model=list( sigma ~ 1 , g0 ~ 1, D~road.dist.clean), mask = Fire.mask)
secr_D_road <- secr_D_river

secr_D_road ## Positive coefficient (but overlapping 0 - not particularly informative or 'significant')

secr_D_road_pred <- predictDsurface(secr_D_road, mask = Fire.mask) ## predictDsurface is a specific predict function in raster package

plot(secr_D_road_pred) ## Missing data on West side of Hopland - should

mean(covariates(secr_D_road_pred)[,14]) ## Average density of 0.266 deer/ha * 945 ha total = 283.5 deer

## To finish, let's visualize the estimated relationship between distance to water and deer

d_road_preds <- predict(secr_D_road, newdata = data.frame(road.dist.clean = seq(0,1020,10))) ## Predict to wide range of distance to water values

plot(seq(0,1020,10), sapply(d_road_preds, "[", "D","estimate"), ylim = c(0.2,0.32), ## Plot it out
     xlab = "Distance from roads (m)", ylab = "Density / ha", type = "l")

