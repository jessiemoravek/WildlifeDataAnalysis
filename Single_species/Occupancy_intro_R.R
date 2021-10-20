## INTRO TO OCCUPANCY IN R 

## We will be using a few new packages today, so you first need to install them if you haven't used them before.
#install.packages("unmarked")
#install.packages("AICcmodavg")
#install.packages("spatialEco")

## We can now load the packaged needed for this exercise
library(dplyr)
library(unmarked)
library(ggplot2)
library(AICcmodavg)
library(spatialEco)

## Load the first dataset. These data are camera trap detection/non-detection data from the AHMbook package (Rota et al. 2016). 
meso<-read.csv(file="occupancy_intro_data.csv", header=T)

##Let's look at the data
str(meso)

## The first column is a site ID. Columns 2-4 correspond to 1 week of camera trap sampling (0=species not detected, 1=species detected). We can also see that we have data for 3 different species

## Load the second dataset
site_cov<-read.csv("occupancy_intro_sites.csv", header=T)

## Let's look at this dataframe
str(site_cov)

## This second dataframe is information on each camera trap site. The ID of these sites match up with the ID of the sites in our detection dataframe. 
#This dataframe has data on environmental conditions at each site including the proportion of disturbed land within 5km (Dist_5km), the housing density within 5 km 
#(HDens_5km), the number of people photographed at a site (People_site), and whether the site was on a trail (Trail)

## Filter out a species of interest
bobcat<-filter(meso, Species=="bobcat")

## In order to run occupancy models, we need to format the data in something called an unmarked dataframe. 
#We need to specify our detection/non-detection data (i.e. the 0's and 1's for each site), which is the "y" portion of this function. 
#To do this, we tell R what columns in our data frame have these numbers (in this case it's columns 2-4). 
#We also need to specify the site covariates that we are interested in, which come from our site covariates dataframe. 
#This is done using the "siteCovs" argument in the function - everything within this portion are out site covariates, which we draw from our site_cov dataframe.
#We are using 1 covariate that we named People, one that we named Dist, one that we named HDesn, and one that we named Trail.

umf_bobcat <- unmarkedFrameOccu(y=bobcat[,c(2:4)],siteCovs=data.frame(People=site_cov$People_site, Dist=site_cov$Dist_5km, HDens=site_cov$HDens_5km, Trail=as.factor(site_cov$Trail)))

## We can now check to make sure we created out unmarked frame correctly  
summary(umf_bobcat)

## We can also plot our unmarked dataframe to visualize detections at each site across our 3 sampling occasions 
plot(umf_bobcat)

## We run occupancy models using the "occu" function from the unmarked package. 
#Within this function, we first specify covariates influencing detection, then covariates influencing occupancy. 
#Let's start with a null model, meaning that nothing is influencing detection or occupancy. 
#When we put "~1" it corresponds with the null model (no covariates on that part of the model)

bob_null=occu(~1 ~1, data=umf_bobcat, se=TRUE)

## We can look at the model output by using the summary function
summary(bob_null)

## We can also get estimates of occupancy and detection rates, since using the summary function only gives up parameter estimates
backTransform(bob_null, "state") ## This gives us estimated occupancy
backTransform(bob_null, "det")   ## This gives us estimated detection probability


## Let's see how the calculated occupancy estimate matches up with the naive occupancy
sum(ifelse((bobcat$X1+bobcat$X2+bobcat$X3)>=1,1,0))/nrow(site_cov)


## Now that we know how to run general occupancy model, we are going to use a 2 stage approach to figuring out the top model. 
#We are first going to run models with different detection covariates while keeping occupancy constant.


## One covatiate in this dataset that we might hypothesize might influence detection is the number of people at a site 
#(e.g. maybe more people make animals more illusive). Let's run this model.

bob_p_people=occu(~People ~1, data=umf_bobcat, se=TRUE)
summary(bob_p_people)


## This time we can only use the backTransform function to get an estimate of occupancy (since detection probability depends on a covariate)
backTransform(bob_p_people, "state")

## For the detection part of the model, we can look at the coefficient estimate and associated confidence intervals to see if it's a positive or negative relationship
coef(bob_p_people, type="det")
confint(bob_p_people, type="det")

## We can see that the detection coefficient is positive (0.13) and that the confidence interval overlaps zero (CI = -1.12 - 1.38), which means that this is not an informative covariate

## We might also think that a camera being on a trail vs. not on a trail might affect detection. Let's run this model
bob_p_trail=occu(~Trail ~1, data=umf_bobcat, se=TRUE)
summary(bob_p_trail)

coef(bob_p_trail, type="det")
confint(bob_p_trail, type="det")

## In this case, the coefficient for is telling us the difference in detection between trails (1) and non trails (0). Because the coefficient is positive, it's saying that detection is higher on trails. The CI doesn't overlap 0, indicating that it's an informative covariate

## Now we need to rank our detection covariate models and see how they compare the to null model. We do this by creating an AIC table. Our top model is the model with the lowest AIC and any model within 2 AIC of the top model is competitive. We first group all of our model results together using the "fitList" and then run the "modSel" function to create our model selection table

bob_det_aic<-fitList(bob_null, bob_p_people, bob_p_trail)
modSel(bob_det_aic)

## Looks like our model with the trail covariate on detection explains the data much better than the other detection model since they have delta AIC values >2. We'll therefore retain this covariate in all models as we now model occupancy.

## Let's first run a model where disturbed land affects occupancy
bob_psi_land=occu(~Trail ~Dist, data=umf_bobcat, se=TRUE)
summary(bob_psi_land)

## For the occupancy part of the model, we can look at the coefficient estimate and associated confidence intervals to see if it's a positive or negative relationship
coef(bob_psi_land, type="state")
confint(bob_psi_land, type="state")


## Our next model will include the effects of housing density on occupancy
bob_psi_hdens=occu(~Trail ~HDens, data=umf_bobcat, se=TRUE)
summary(bob_psi_hdens)


## Let's look at coefficient estimates and confidence intervals
coef(bob_psi_hdens, type="state")
confint(bob_psi_hdens, type="state")

## Looks like occupancy is negatively related to housing density

## We now want to see which model fits the data the best. We'll compare the two models we just ran with the null model for occupancy when we included our trail detection covariate

bob_occ_aic<-fitList(bob_p_trail, bob_psi_land, bob_psi_hdens)
modSel(bob_occ_aic)


## It looks like our housing density model fits best, but before we get too far into ranking models we also need to test the fit of our model. To do this, we will use the MacKenzie Bailey Goodness of Fit (GoF) test. With this procedure, the standard practice is to look at the fit of a global model (i.e. all covariates included). Normally we'd run this with at least 1000 simulations, but it would take a while to run, so this is just to get an idea of how the code works 

test<-occu(~Trail+People ~HDens+Dist, data=umf_bobcat, se=TRUE)

mb.gof.test(test, nsim = 15)

## To interpret the GoF output, we look at 1) the p-value and 2) the c-hat value. A c-hat around 1 and a p-value >0.05 indicate the the data fit the model well. A c-hat >1 indicates overdispersion and <1 indicated underdispersion. A c-hat greater than ~4 suggests a lack of fit and we might want to re-think out model structure. Maybe there are other covariates we need to include or maybe our encounter length should be changed. 


## If there is evidence for dispersion, we will want to incorporate this into our model ranking. We can do this using the AICcmodav package. Here's an example based on the GoF test that we just ran

aictab(list(bob_p_trail, bob_psi_hdens, bob_psi_land), c.hat=2.28)

## If we adjust the model ranking based on c-hat, we also need to adjust the standard error associated with any predictions that we make. We aren't going to do that in this exercise since we didn't do a full GoF test, but if you needed to do that, you can also do so with the AICmodav package. For now we are going to pretend that we ran a full GoF test and didn't find evidence for dispersion.


##Let's make a predictive plots showing the occupancy relationship, as well as our detection relationship from the top model

## Starting with our predicted occupancy, we first want to make a new dataframe with values representing our range of observed values. We will do this by creating a sequence of number from our minimum observed value to our maximum observed value. It's important to name the column the same as what we called it in our model (in this case, "HDens"), Depending on the range of values, you can change the "by" argument
newdat=data.frame(HDens=seq(from=min(site_cov$HDens_5km), to=max(site_cov$HDens_5km), by=0.1))

##We can now use the predict function on our new data. We use type="state" to specify that we are predicting occupancy rather than detection. If we were predicting occupancy, we would use type="det". We specify the model that we want to use to predict (bib_psi_hdens) and that our new dataframe is called "newdat" in this case. We add "append=T" to include the new data in the dataframe with the predicted values

Predicted_bob_occ=predict(bob_psi_hdens,type="state",newdata=newdat, append=T)


## To make a graph, we are going to make a line of our predicted values vs. our covariate values, and add a shaded region to represent our confidence intervals

ggplot(Predicted_bob_occ, aes(y=Predicted, x=HDens))+
  geom_line()+ ## This adds the line
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+ ## This adds the shaded region. The alpha component changes how opaque the shading is
  xlab("Housing density within 5 km")+ ## This adds x label
  ylab("Occupancy probability")+ ## This adds y axis
  theme_classic()  ## This takes away the background grid

## Now we can make a similar plot for detection. First we'll make the data to predict on. Our top detection covariate is trail, which is a binary covariate (0 or 1), so this time we just need to make a dataframe with those two options
newdat2=data.frame(Trail=seq(from=0, to=1))
newdat2$Trail<-as.factor(newdat2$Trail)

##We can now use the predict function on our new data. This time we'll specify type="det"
Predicted_bob_det=predict(bob_psi_hdens,type="det",newdata=newdat2, append=T)

## Since we don't have a continuous covariate for detection, we'll make a different type of model this time - a point estimate with an error bar

ggplot(Predicted_bob_det, aes(y=Predicted, x=Trail))+
  geom_point()+ ## This adds the line
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.2)+
  xlab("Trail")+ ## This adds x label
  ylab("Detection probability")+ ## This adds y axis
  theme_classic()

## The final thing we might want to do is make a predictive map based on occupancy. I don't have associated raster data for this dataset, so we are going to randomly generate a raster based on the minimum/maximum housing density values so that we have a raster to work with

rand_ras<-random.raster(n.row = 50, n.col = 50, min=min(site_cov$HDens_5km), max=max(site_cov$HDens_5km), distribution = "random")

## We need to make sure the name of the name of the covariate is the same as the name of the data in the model
names(rand_ras)<-"HDens"

## We now have our fake raster of housing density
plot(rand_ras)

## We can use the predict function and specify this raster for the newdata argument
ras_predict<-predict(bob_psi_hdens, type="state",newdata=rand_ras)

## The prediction saves as a raster stack with 4 separate rasters - the predicted values (in this case predicted occupancy), the SE associated with the predicitons and the upper and lower confidence intervals. If we were making a preditive map for a publication/report we'd want to format these a little better (e.g. make sure the scales match up, add boundaries/landmarks etc) but we aren't going to spend time on that now.
plot(ras_predict)

