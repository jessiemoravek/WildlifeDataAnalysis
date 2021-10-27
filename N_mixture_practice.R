## N-mixture models in R ####

## Load in the librarys for this exercise
library(unmarked)
library(AHMbook)
library(dplyr)
library(AICcmodavg)
library(ggplot2)
library(raster)
library(tidyverse)

## Load in the dataset
birds <- read.csv("hubbard_brook_birds_2018.csv")
str(birds)

## Pull a species
BHVI <- birds %>% filter(species == "BHVI")

## Next well pull out the count data
counts<-BHVI[,c(2:4)]

## We might also be interested in the time of day that the points were surveyed 
time<-as.data.frame(BHVI[,c(9:11)])

## Make a data frame for site covariates
elev <- as.data.frame(BHVI[,14])
str(elev)
colnames(elev)[1] <- "elev"

## Pull forest data
forest <- as.data.frame(BHVI[,17])

## We can now make our unmarked dataframe. Similar to occupany, we specify our count data (y), our site covariates, and our observation covariates
umf<-unmarkedFramePCount(y=counts, siteCovs = elev, obsCovs = list(time=time))


## First we are going to figure out the mixture most suited for the data. We will try Poisson vs. Negative Binomial
null_p<-pcount(~1 ~1, data=umf, mixture="P")                         
null_nb<-pcount(~1 ~1, data=umf, mixture="NB")

mix_aic<-fitList(null_p, null_nb)
modSel(mix_aic)

## The negative binomial model fit better, so we will use that moving forward and next test models that affect detection while holding abundance constant. For this type of model, we specify detection covariates first in the formula, followed by abundance covariates

## A model where time of day affects detection
p_time<-pcount(~time ~1, data=umf, mixture="NB")
summary(p_time)

## A model where the elevation affects detection
elev <- sapply(elev[1],as.numeric)
p_elev<-pcount(~elev ~1, data=umf, mixture="NB")
summary(p_elev)

## We can now rank the detection models
p_aic<-fitList(null_nb, p_time, p_elev)
modSel(p_aic)

### Keep null detection model and now run models with covariates on abundance

## elevation affects abundance
abun_elev<-pcount(~elev ~1, data=umf, mixture="NB")
summary(abun_elev)

## forest affects abundance
forest <- sapply(forest[1],as.character)
abun_forest<-pcount(~forest ~1, data=umf, mixture="NB")
summary(abun_forest)

## Rank the abundance models
abun_aic<-fitList(p_time, abun_area, abun_edge, abun_height, abun_dens)
modSel(abun_aic)

## The model with willow height has the lowest AIC. We can run a goodness of fit test on this model, this time using the NMix.gof.test function. We'll only do 20 simulations but ideally you'd want 1000+
Nmix.gof.test(abun_height, nsim=20, report=T)

#C-hat should be close ish to 1
#red line of simulated statistic should be in the middle (slightly under dispersed but bc only 20 sim)
## Even with a small number of simulations, there's nothing obvious to suggest a poor fit

## We can now make some predictive plots starting with abundance
newdat=data.frame(height=seq(from=min(sites$height), to=max(sites$height), by=0.1))

Predicted_abun=predict(abun_height,type="state",newdata=newdat, append=T)

ggplot(Predicted_abun, aes(y=Predicted, x=height))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  ylab("Predicted abundance")+
  xlab("Willow height (cm)")+
  theme_classic()

## Now a plot for detection
newdat2<-data.frame(time=seq(from=min(time_2006, na.rm = T), to=max(time_2006, na.rm=T), by=1))

Predicted_det<-predict(abun_height, type="det", newdata=newdat2, append=T)

## The times are in military time expresses as thousands, so we are going to convert the number to make it in a more understandable format
Predicted_det$time2<-Predicted_det$time/100

ggplot(Predicted_det, aes(y=Predicted, x=time2))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.5)+
  ylab("Detection probability")+
  xlab("Hour of day")+
  theme_classic()


## The last thing we might want to do is get a total abundance number across our entire study area, or map out how abundance varies spatially. 

## The first thing we need to do is figure out our effective sampling area. These point counts were conducted in 100 m radius circles. This means the estimated abundance is associated with an area of 31415.9 square m (area of a circle = pi*r^2) or 0.03 km2

## Import raster. Note - this is fake data which is really elevation from my PhD study area. The scale of numbers matched up closely enough to make it work for this example
height_ras<-raster("fake_height_data.tif")

## Make sure the name is the same as what was in the model
names(height_ras)<-"height"

plot(height_ras)

## We can check to see how large our raster cells are
res(height_ras)

## We'll round and say that they are 30 by 30. This is an issue because each cell is only 900 sq m compared to the 31415 sq m of the sample point area. One options is to resample the raster to the correct cell size

## We are first going to aggregate our raster. The factor that we are using is the length of a cell equal to the plot area divide by the current length of cell
height2<-aggregate(height_ras, fact=sqrt(31415)/30)
plot(height2)

##Let's check
res(height2)
sqrt(31415)

## Close enough for this exercise - the rounding is a bit off but we'll ignore that

## In order for the prediction to run, the raster needs to be saved as a raster stack
height_stack<-stack(height2)

## We can now run the predict function on the raster
height_predict<-predict(abun_height, newdata=height_stack, type="state")
plot(height_predict)

## We can now calculate the predicted total population size in our fake study area, as well as the lower and upper confidence intervals
cellStats(height_predict$Predicted, stat="sum") ## This is the total estimate
cellStats(height_predict$lower, stat="sum")
cellStats(height_predict$upper, stat="sum")
