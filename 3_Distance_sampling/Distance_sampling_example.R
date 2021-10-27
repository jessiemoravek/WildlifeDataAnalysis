## Distancing Sampling to Estimate Abundance ## 

## Load packages
library(unmarked)
library(AHMbook)
library(AICcmodavg)

## One useful thing to know how to do is format data you collected in the field to get it to the correct format for distance sampling - your distances will need to be in distance bins. In some cases you might have collected the data in this format (e.g. count the number of bird within 50 m and between 50-100m). In other cases, you might have measured distances that you need to put into bins

## We are going to make some fake data to demonstrate how to format data. You don't have to worry too much about this line of code but it basically makes a fake dataframe with a column with transect IDs and a column of distances
dist_dat<-data.frame(transect=gl(6,5, labels=letters[1:6]),distance=rpois(30, 15))
summary(dist_dat)

## We are going to pretend that there is another transect where we didn't detect any animals
levels(dist_dat$transect) <- c(levels(dist_dat$transect), "g")
levels(dist_dat$transect)

## Let's take a look at a histogram of our distances.
hist(dist_dat$distance)

## This isn't completely realistic because we randomly generated the data, but we'd want to look at the distribution of the data to consider different bin sizes. We want at least 3 bins of equal sizes, but can play around with the breakpoints. 

## This function is how we actually format the data. We specify the raw dataframe name, the name of the column containing the distances, the name of the column containing the transect/point IDs, and we specify the breakpoints for our bins
y_format<-formatDistData(dist_dat, distCol="distance",transectNameCol="transect", dist.breaks=c(0, 10, 20, 30))

y_format


## Now that you know how to format your own data, we are going to use data included with the unmarked package from the island scrub jay that was collected in a distance sampling framework. This data is from the Sillett et al. 2012 paper.

data(issj)

## Let's look at the data structure
str(issj)

## These data were collected based on point counts and each row of data was a survey point. The first 3 columns are the count of birds in each of 3 distance categories (0-100m, 100-200m, and 200-300m). The remaining columns are environmental covariates associated with the survey location.

## Just like with occupancy, we are going to be using the unmarked package and need to format our data into an unmarked dataframe. A few difference here: 1) the y data (observed counts) needs to be in a matrix format, 2) we have to specify the distances we binned the data into, 3) we have to specify the units associated with the distance, 4) we have to specify if the data was collected at point or a transect. If we did have data collected along a transect, we'd also add the "tlength" argument here to specify the total length of the transects

umf<-unmarkedFrameDS(y=as.matrix(issj[,1:3]), siteCovs = issj[,6:8], dist.breaks = c(0,100,200,300), unitsIn = "m", survey = "point")

summary(umf)

plot(umf)

## Similar to before, we are first going to determine the best detection model. We will be using a function called "distsamp." In this function, we first specify the detection covariates followed by the abundance covariates. THe first thing we are going to do is test different detection functions using a null model (no covariates on detection or abundance)

null_halfnorm<-distsamp(~1 ~1, data=umf, keyfun="halfnorm", output="abund")

null_hazard<-distsamp(~1 ~1, data=umf, keyfun = "hazard", output="abund")

null_exp<-distsamp(~1 ~1, data=umf, keyfun = "exp", output="abund")
## The warning message tells us there are convergence issues with the expoential model.

null_uniform<-distsamp(~1 ~1, data=umf, keyfun = "uniform", output="abund")

## Rank the detection function models
det_fun<-fitList(null_halfnorm, null_hazard, null_exp, null_uniform)
modSel(det_fun)

## Let's take a look at some estimates from the top model to see if they are reasonable. Distance functions can sometimes be finicky and generate unreasonable estimates depending on the data they are trying to fit
backTransform(null_hazard, type="state")

## This estimates the abundance at a site. We can see that the SE value is nearly as large as the estimate, which suggests there is an issue with the model.

## We can compare it to the estimates from the second best model
backTransform(null_halfnorm, type="state")


##The SE estimates for the halfnormal function are more reasonable, so we are going to move forward using it for all our models. Next, we'll look at covariates affecting detection while holding abundance constant. We'll try model with forest and chaparral affecting detection (maybe there is something about the habitats that affects how well we hear or see the birds)

det_forest<-distsamp(~forest ~1, data=umf, keyfun = "halfnorm", output="abund")
det_chap<-distsamp(~chaparral ~1, data=umf, keyfun = "halfnorm", output="abund")

## We'll now rank these with the null model included
det_aic<-fitList(null_halfnorm, det_forest, det_chap)
modSel(det_aic)

## Looks like chaparral is the best detection covariate, so we'll retain it and continue on with abunance covariates

abun_elev<-distsamp(~chaparral ~elevation, data=umf, keyfun = "halfnorm", output="abund")
abun_for<-distsamp(~chaparral ~forest, data=umf, keyfun = "halfnorm", output="abund")
abun_chap<-distsamp(~chaparral ~chaparral, data=umf, keyfun = "halfnorm", output="abund")


abun_aic<-fitList(null_halfnorm, abun_elev, abun_for, abun_chap)
modSel(abun_aic)

## We can now run a goodness of fit test on our top model. The GoF function we used last week doesn't work for this time of model, so we'll use this "fitstats" function in conjunction with the parboot function

fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2)
  chisq <- sum((observed - expected)^2 / expected)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

gof<-parboot(abun_chap,statistic=fitstats,  nsim=25)

## We can look at the results here. This time there are 3 GoF statistics (SEE, Chi Square, and Freeman Tukey). THe p-values for all of them suggests a lack of fit (P<0.05)
gof

## This gives us the chi-square c-hat value
gof@t0[2]/mean(gof@t.star[,2])


## Ideally we'd want to use 1000+ simulations for the GoF test. If you read the Sillett et al paper or the AHM book, you'll know that there are some issue with overdispersion with the dataset that we are using. You can read the paper to learn how they adapted their models to account for this, but it's beyond what we are learning today. In short - if you are running distance sampling models with your own data and the GoF test indicates overdispersion, there are ways to account for it.


## We are going to move on and make predictive graphs based on the top model
newdat=data.frame(chaparral=seq(from=min(issj$chaparral), to=max(issj$chaparral), by=0.01))

Pred_abun<-predict(abun_chap, type="state", newdata=newdat, appendData=T)

ggplot(Pred_abun, aes(y=Predicted, x=chaparral))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  theme_classic()+
  ylab("Predicted abundance")+
  xlab("Proportion chaparral")

## We might also want to express the estimate as a density estimate. You'll remember when we ran the models that we specified output=abundance. We can change that if we want estimates as densities and the units that we want it in (ha or kmsq)

abun_chap2<-distsamp(~chaparral ~chaparral, data=umf, keyfun = "halfnorm", output="density", unitsOut = "ha")

Pred_dens<-predict(abun_chap2, type="state", newdata=newdat, appendData=T)

ggplot(Pred_dens, aes(y=Predicted, x=chaparral))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  theme_classic()+
  ylab("Predicted density (birds per ha)")+
  xlab("Proportion chaparral")

## Our predictive figure for detection works a little bit differently this time. Our top model tells us that chaparral affect detection, but what this really means is that the relationship between distance and detection probability depends on the proportion of chaparral in a plot. This makes it a bit hard to visualize - we could either make a 3-d plot or see what the detection function looks like at different values of the chaparral covariate. We are going to use this second option.

## We'll first use the predict function on the new chaparral dataframe but specify type=det
Pred_dens_chap<-predict(abun_chap, type="det", newdata=newdat, appendData=T)

## The predicted values are the sigma values associated with the detection function, which correspond with the shape of the detection function. We'll pull out specific chaparral values and the associated sigmas. We then use the gxhn function which corresponds with the half normal detection function. We could use gxhaz for the hazard function. We also specify the lowest and highest distance, as well as the predicted sigma from our prediction dataframe

#First is when chaparral is 0%
det_graph0<-as.data.frame(gxhn(c(0:300), 150.82366))
colnames(det_graph0)<-"Prob"
det_graph0$Distance<-rep(0:300)
det_graph0$Chaparral<-"0%"

## We can do the same thing when chaparral is 50%
det_graph50<-as.data.frame(gxhn(c(0:300), 88.62775))
colnames(det_graph50)<-"Prob"
det_graph50$Distance<-rep(0:300)
det_graph50$Chaparral<-"50%"

## Last we'll do it at the highest chaparral percent
det_graph94<-as.data.frame(gxhn(c(0:300), 55.51085))
colnames(det_graph94)<-"Prob"
det_graph94$Distance<-rep(0:300)
det_graph94$Chaparral<-"94%"

## We can now combine them and plot
det_graph_all<-rbind(det_graph0, det_graph50, det_graph94)

ggplot(det_graph_all,aes(y=Prob, x=Distance, color=Chaparral))+
  geom_line()+
  ylab("Detection probability")+
  xlab("Distance (m)")+
  theme_classic()
