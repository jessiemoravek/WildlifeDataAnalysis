## OCCUPANCY PRACTICE

## We can now load the packaged needed for this exercise
library(dplyr)
library(unmarked)
library(ggplot2)
library(AICcmodavg)
library(spatialEco)

## Load the first dataset. These data are camera trap detection/non-detection data from the AHMbook package (Rota et al. 2016). 
mamm<-read.csv(file="Det-NonDet_mammals.csv", header=T)

##Let's look at the data
str(mamm)

## Load the second dataset
site_cov<-read.csv("Site_covs_scaled.csv", header=T)

## Let's look at this dataframe
str(site_cov)

## This second dataframe is information on each camera trap site. The ID of these sites match up with the ID of the sites in our detection dataframe. 

## Filter out a species of interest
arm<-filter(mamm, species=="Dspp")
str(arm)
## In order to run occupancy models, we need to format the data in something called an unmarked dataframe. 
#We need to specify our detection/non-detection data (i.e. the 0's and 1's for each site), which is the "y" portion of this function. 
#To do this, we tell R what columns in our data frame have these numbers (in this case it's columns 2-4). 
#We also need to specify the site covariates that we are interested in, which come from our site covariates dataframe. 
#This is done using the "siteCovs" argument in the function - everything within this portion are out site covariates, which we draw from our site_cov dataframe.
#We are using 1 covariate that we named People, one that we named Dist, one that we named HDesn, and one that we named Trail.

umf_arm <- unmarkedFrameOccu(y=arm[,c(3:14)],siteCovs=data.frame(site=as.factor(site_cov$Sampling_site), dist_road=site_cov$Dist_road, dist_water=site_cov$Dist_water, ndvi=site_cov$NDVImean_500m, pa_type=as.factor(site_cov$pa_type), CTcode2=as.factor(site_cov$CTcode2), trail = as.factor(site_cov$trail)))

## We can now check to make sure we created out unmarked frame correctly  
summary(umf_arm)

## We can also plot our unmarked dataframe to visualize detections at each site across our 3 sampling occasions 
plot(umf_arm)

## We run occupancy models using the "occu" function from the unmarked package. 
#Within this function, we first specify covariates influencing detection, then covariates influencing occupancy. 
#Let's start with a null model, meaning that nothing is influencing detection or occupancy. 
#When we put "~1" it corresponds with the null model (no covariates on that part of the model)

arm_null=occu(~1 ~1, data=umf_arm, se=TRUE)

## We can look at the model output by using the summary function
summary(arm_null)

## We can also get estimates of occupancy and detection rates, since using the summary function only gives up parameter estimates
backTransform(arm_null, "state") ## This gives us estimated occupancy
backTransform(arm_null, "det")   ## This gives us estimated detection probability

## Now that we know how to run general occupancy model, we are going to use a 2 stage approach to figuring out the top model. 
#We are first going to run models with different detection covariates while keeping occupancy constant.


## One covatiate in this dataset that we might hypothesize might influence detection is whether it's on a trail
arm_p_trail=occu(~trail ~1, data=umf_arm, se=TRUE)
summary(arm_p_trail)


## This time we can only use the backTransform function to get an estimate of occupancy (since detection probability depends on a covariate)
backTransform(arm_p_trail, "state")

## For the detection part of the model, we can look at the coefficient estimate and associated confidence intervals to see if it's a positive or negative relationship
coef(arm_p_trail, type="det")
confint(arm_p_trail, type="det")

## We can see that the detection coefficient is positive (0.13) and that the confidence interval overlaps zero (CI = -1.12 - 1.38), which means that this is not an informative covariate

## We might also think that a camera being on a road might affect detection. Let's run this model
arm_p_road=occu(~dist_road ~1, data=umf_arm, se=TRUE)
summary(arm_p_road)

coef(arm_p_road, type="det")
confint(arm_p_road, type="det")

## try pa type
arm_p_pa=occu(~pa_type ~1, data=umf_arm, se=TRUE)
summary(arm_p_pa)

coef(arm_p_pa, type="det")
confint(arm_p_pa, type="det")

## try distance to water

arm_p_water=occu(~dist_water ~1, data=umf_arm, se=TRUE)
summary(arm_p_water)

coef(arm_p_water, type="det")
confint(arm_p_water, type="det")

## Now we need to rank our detection covariate models and see how they compare the to null model. We do this by creating an AIC table. Our top model is the model with the lowest AIC and any model within 2 AIC of the top model is competitive. We first group all of our model results together using the "fitList" and then run the "modSel" function to create our model selection table

arm_det_aic<-fitList(arm_null, arm_p_trail, arm_p_road, arm_p_pa, arm_p_water)
modSel(arm_det_aic)

## Null model for detection is best. Now we move on to occupancy

## Let's first run a model where pa type affects occupancy
str(site_cov)
arm_psi_pa=occu(~1 ~pa_type, data=umf_arm, se=TRUE)
summary(arm_psi_pa)

## For the occupancy part of the model, we can look at the coefficient estimate and associated confidence intervals to see if it's a positive or negative relationship
coef(arm_psi_pa, type="state")
confint(arm_psi_pa, type="state")


## Our next model will include the effects of housing density on occupancy
arm_psi_water=occu(~1 ~dist_water, data=umf_arm, se=TRUE)
summary(arm_psi_water)


## Let's look at coefficient estimates and confidence intervals
coef(arm_psi_water, type="state")
confint(arm_psi_water, type="state")

## Looks like occupancy is not related to distnace to water

## We now want to see which model fits the data the best. We'll compare the two models we just ran with the null model for occupancy when we included our trail detection covariate

arm_occ_aic<-fitList(arm_null, arm_psi_pa, arm_psi_water)
modSel(arm_occ_aic)


## Looks like null model is best, but need to test the fit of our model. To do this, we will use the MacKenzie Bailey Goodness of Fit (GoF) test. With this procedure, ran the top model

test<-occu(~1 ~dist_water, data=umf_arm, se=TRUE)

mb.gof.test(test, nsim = 15)

## To interpret the GoF output, we look at 1) the p-value and 2) the c-hat value. A c-hat around 1 and a p-value >0.05 indicate the the data fit the model well. A c-hat >1 indicates overdispersion and <1 indicated underdispersion. A c-hat greater than ~4 suggests a lack of fit and we might want to re-think out model structure. Maybe there are other covariates we need to include or maybe our encounter length should be changed. 


## If there is evidence for dispersion, we will want to incorporate this into our model ranking. We can do this using the AICcmodav package. Here's an example based on the GoF test that we just ran

aictab(list(bob_p_trail, bob_psi_hdens, bob_psi_land), c.hat=2.28)

## If we adjust the model ranking based on c-hat, we also need to adjust the standard error associated with any predictions that we make. We aren't going to do that in this exercise since we didn't do a full GoF test, but if you needed to do that, you can also do so with the AICmodav package. For now we are going to pretend that we ran a full GoF test and didn't find evidence for dispersion.


##Let's make a predictive plots showing the occupancy relationship, as well as our detection relationship from the top model

## Starting with our predicted occupancy, we first want to make a new dataframe with values representing our range of observed values. We will do this by creating a sequence of number from our minimum observed value to our maximum observed value. It's important to name the column the same as what we called it in our model (in this case, "HDens"), Depending on the range of values, you can change the "by" argument
newdat=data.frame(dist_water=seq(from=min(site_cov$Dist_water), to=max(site_cov$Dist_water), by=0.1))

##We can now use the predict function on our new data. We use type="state" to specify that we are predicting occupancy rather than detection. If we were predicting occupancy, we would use type="det". We specify the model that we want to use to predict (bib_psi_hdens) and that our new dataframe is called "newdat" in this case. We add "append=T" to include the new data in the dataframe with the predicted values

Predicted_arm_occ=predict(arm_psi_water,type="state",newdata=newdat, append=T)


## To make a graph, we are going to make a line of our predicted values vs. our covariate values, and add a shaded region to represent our confidence intervals

ggplot(Predicted_arm_occ, aes(y=Predicted, x=dist_water))+
  geom_line()+ ## This adds the line
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+ ## This adds the shaded region. The alpha component changes how opaque the shading is
  xlab("Distance to Water (km)")+ ## This adds x label
  ylab("Occupancy probability")+ ## This adds y axis
  theme_classic()  ## This takes away the background grid

## Now we can make a similar plot for detection. First we'll make the data to predict on. Our top detection covariate is trail, which is a binary covariate (0 or 1), so this time we just need to make a dataframe with those two options
newdat2=data.frame(water=seq(from=0, to=1))
newdat2$water<-as.factor(newdat2$water)

##We can now use the predict function on our new data. This time we'll specify type="det"
Predicted_arm_det=predict(arm_psi_water,type="det",newdata=newdat2, append=T)

## Since we don't have a continuous covariate for detection, we'll make a different type of model this time - a point estimate with an error bar

ggplot(Predicted_arm_det, aes(y=Predicted, x=water))+
  geom_point()+ ## This adds the line
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.2)+
  xlab("Distance to water (km)")+ ## This adds x label
  ylab("Detection probability")+ ## This adds y axis
  theme_classic()

## The final thing we might want to do is make a predictive map based on occupancy. I don't have associated raster data for this dataset, so we are going to randomly generate a raster based on the minimum/maximum housing density values so that we have a raster to work with

rand_ras<-random.raster(n.row = 50, n.col = 50, min=min(site_cov$Dist_water), max=max(site_cov$Dist_water), distribution = "random")

## We need to make sure the name of the name of the covariate is the same as the name of the data in the model
names(rand_ras)<-"Water"

## We now have our fake raster of housing density
plot(rand_ras)

## We can use the predict function and specify this raster for the newdata argument
ras_predict<-predict(arm_psi_water, type="state",newdata=rand_ras)

## The prediction saves as a raster stack with 4 separate rasters - the predicted values (in this case predicted occupancy), the SE associated with the predicitons and the upper and lower confidence intervals. If we were making a preditive map for a publication/report we'd want to format these a little better (e.g. make sure the scales match up, add boundaries/landmarks etc) but we aren't going to spend time on that now.
plot(ras_predict)

