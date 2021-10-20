## Dynamic Occupancy Models in R ###

## Load packages
library(unmarked)
library(dplyr)
library(ggplot2)
library(AICcmodavg)

## Import the "crossbill" dataset from the unmarked package
data("crossbill")


## This dataset includes detection/non-detection data at 267 sites in Switzerland across 9 years (1999 - 2007). Each site was typically surveyed 3 times a year. Site covariates include: ele (elevation of the site), forest (percentage of forest at a site), and surveys (the number of times a site was monitored in a year; either 2 or 3). It also has multiple "date" columns that correspond to the Julian day of each survey
str(crossbill)

## If we look closely at the Date data, we can see some NA values. We can't have NAs in the data for the models to run, so we need to remove any observations that don't contain all the data we are interested in
crossbill2<-na.omit(crossbill)

## We know that the data were collected in different years, but we need to create a year covariate since it isn't included in the raw data.
years <- as.character(1999:2007)
years <- matrix(years, nrow(crossbill2), 9, byrow=TRUE)

## We are also interested in including date as a covariate so we need to create a date covariate from the data the crossbill dataframe
date <- as.matrix(crossbill2[,32:58])


## We can now format our data into an unmarked dataframe. Similarly to before, we specify "y" as the columns in our data frame that have the detection/non-detection data. We will also include site covariates, which are the elevation and forest columns, a yearly site covariate which is the year dataframe we made, and an observation covariate which is the date of the survey. We also have to specify how many primary occassions are in the data which is this case is 9 (i.e. the number of years)
umf <- unmarkedMultFrame(y=crossbill2[,c(5:31)], siteCovs=crossbill2[,2:3], yearlySiteCovs=list(year=years),obsCovs=list(date=date), numPrimary=9)

summary(umf)
plot(umf)


## To run a dynamic occupancy model, we will be using the "colext" function in the unmarked package. The order of arguments in the function are 1) covariates for probability of occupancy, 2) covariates for probability of colonization, 3) covariates for probability of extinction, 4) covariates for detection probability

# First let's run a null model
null <- colext(~1, ~1, ~1, ~1, dat=umf)
summary(null)

## We can now look at estimates from the model
backTransform(null, type="psi")   ## Estimate of occupancy probability
backTransform(null, type="col")   ## Estimate of colonization probability
backTransform(null, type="ext")   ## Estimate of extinction probability
backTransform(null, type="det")   ## Estimate of detection probability


## We are going to employ a similar strategy as before to avoid model dredging - first we are going to investigate detection models while holding everything else constant

## First we'll try a model where detection differs by year (maybe there are different people doing the surveys), followed by a model where detection differs by date of survey
p_year<-colext(~1,~1, ~1, ~year, dat=umf)

p_date<-colext(~1, ~1, ~1, ~date, dat=umf)


## Let's rank these models and compare them to the null model
det_aic<-fitList(null, p_year, p_date)
modSel(det_aic)

## Looks like detection varies by year. We'll retain this covariate and start running models that match up with hypotheses related to occupancy, colonization, and extinction. There are many different hypotheses that we might be interested in, but for the sake of time we are only going to run models where forest or elevation affect occupancy, colonization, and extinction 

forest_mod<-colext(~forest, ~forest, ~forest, ~year, dat=umf, se=TRUE)
elev_mod<-colext(~ele, ~ele, ~ele, ~year, dat=umf, se=TRUE)


## Let's rank these models and compare them to the null model
all_aic<-fitList(p_year, forest_mod, elev_mod)
modSel(all_aic)

## We can also run a similar GoF test like we did for the first exercise. We are going to skip that for this exercise (bad practice, but this is just an example exercise!)

## Now we want to plot our predictions. We'll start with occupancy, colonization, and extinction
newdat=data.frame(ele=seq(from=min(crossbill2$ele), to=max(crossbill2$ele), by=0.5))

Pred_occ<-predict(elev_mod, type="psi", newdata=newdat, appendData=T)

ggplot(Pred_occ, aes(y=Predicted, x=ele))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  theme_classic()+
  ylab("Occupancy probability")+
  xlab("Elevation (m)")


Pred_col<-predict(elev_mod, type="col", newdata=newdat, appendData=T)

ggplot(Pred_col, aes(y=Predicted, x=ele))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  theme_classic()+
  ylab("Colonization probability")+
  xlab("Elevation (m)")

Pred_ext<-predict(elev_mod, type="ext", newdata=newdat, appendData=T)

ggplot(Pred_ext, aes(y=Predicted, x=ele))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  theme_classic()+
  ylab("Extinction probability")+
  xlab("Elevation (m)")

## Finally out detection model will have a different dataset to predict on since we found that detection varied by year
newdat2=data.frame(year=seq(from=1999, to=2007))
newdat2$year<-as.factor(newdat2$year)

Pred_det<-predict(elev_mod, type="det", newdata=newdat2, appendData=T)

ggplot(Pred_det, aes(y=Predicted, x=year))+
  geom_point()+
  geom_errorbar(aes(ymax=upper, ymin=lower))+
  theme_classic()+
  ylab("Extinction probability")+
  xlab("Year")
