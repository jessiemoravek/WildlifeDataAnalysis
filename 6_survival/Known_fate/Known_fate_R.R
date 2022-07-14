## Known-fate Survival with RMark ###

##Install some new packages 
#install.packages("RMark")
#install.packages("msm")

## Load packages
library(RMark)
library(dplyr)
library(ggplot2)
library(msm)
setwd('C:/Users/Owner/Desktop/WildlifeDataAnalysis/WildlifeDataAnalysis/6_survival/Known_fate')

## We are going to use a known-fate dataset called Blackduck that is included with the RMark package. We can import it using the import.chdata function. We can use the field.types argument to specify if columns follwing the first column are a factor (f) or numeric (n)
Blackduck<-import.chdata("blackduck.txt", header=T, field.types = c("f", "n","n","n","f",rep("n", 8)))


## Let's take a look at the data
str(Blackduck)

## The column called "ch" is the known-fate capture history. The other columns are individual covariates (age, weight, wing length, condition index), or time-varying covriates (Temperature). For the time varying covariate, there is one value for each of the 8 occassions (we will pretend that each occassion is a week)


## We need to format the data in the correct format for the models to run. We first use the "process.data" function and specify type=known to tell the package that the data are known-fate data. We also need to specify that we have a grouping variable in the data (sex) 
duck_processed<-process.data(Blackduck, model="Known", groups=c("Sex"))

## Next we need to make a design matrix with the data. This creates indexed parameters based on our processed data. For example, if we look at the design matrix, we can see that we can estimate a survival covariate for each of the occasions and for males vs females
duck_ddl<-make.design.data(duck_processed)
duck_ddl

## Now we can run some models. For each model, we need to specify the formula with covariates for survival in a list. We'll start with the null model, which is often called the "dot" model
S_dot<-list(formula=~1)
S_dot_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_dot))

S_dot_mod

## If we look at the model, we can see a covariate for survival as well as an estimated survival rate. Because this is the null model, this is telling us the survival in any given month for any bird (male or female)

## Next we can run a model where survival depends on the bird's weight
S_weight<-list(formula=~Weight)
S_weight_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_weight))
S_weight_mod

## Next, a model where wing length affects survival
S_wing<-list(formula=~Wing_Len)
S_wing_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_wing))
S_wing_mod

## A model where both weight and wing length affect survival
S_weight_wing<-list(formula=~Weight+Wing_Len)
S_weight_wing_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_weight_wing))

## A model where survival differs between males and females
S_sex<-list(formula=~Sex)
S_sex_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_sex))

## This time if we look at the model output, we see we have a covariate for Sex:M. This means that female is our reference cateogory and the covariate tell us if male survival is different than females. The covariate is positive but the 95% CI overlaps zero, meaning that males have higher survival but not significantly different from females
S_sex_mod

## Next we might be interested in seeing if survival differs by time. Right now we are treating survival constant across months. The design matrix has a built in "time" covariate (note: There is also a built in "Time" covariate - the capital one specifies time as a continuous variable whereas the lowercase is discrete. We are assuming that there isn't a trend by time so we are using the lowercase one)
S_time<-list(formula=~time)
S_time_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_time))

S_time_mod

## If we look at the output, we can see there are issues. The last survival period has a super high SE and confidence intervals. This indicates there are convergence issues. There are a few things we could do- we could decide to re-group the occasions (maybe estimate survival during 2 month periods), or we can look at the raw data and see if there's anything going on. 

#When we look at the data, we see that there were no deaths in the last month. This means that survival during this period is 1, which is why the model is having trouble (i.e. there's no error and the CI should be from 1 to 1). In this case since we know survival in 1 during this period, we can fix the value to 1 so that the model doesn't have to estimate it. (NOTE: ONLY DO THIS IF YOU KNOW FOR SURE THAT YOUR SURVIVAL IS 1!)

## To fix the data, we need to know the index of the covariates that we want to fix in our design matrix
fixed<-as.numeric(row.names(duck_ddl$S[which(duck_ddl$S$time==9),]))

## Let's re-run the model but this time fix the final month's survival to 1
S_time<-list(formula=~time, fixed=list(index=fixed, value=1))
S_time_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_time))

S_time_mod

## Now let's run a model where temperature affects survival.
S_temp<-list(formula=~Temp)
S_temp_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_temp))                  

S_temp_mod
## We can see that we have a covariate for temperature. The estimates are based on average temperature values in each month

## Finally, maybe we think that temperature affects males vs. female differently
S_temp_sex<-list(formula=~Temp*Sex)
S_temp_sex_mod<-mark(data=duck_processed, ddl=duck_ddl, model.parameters = list(S=S_temp_sex))

S_temp_sex_mod

## There are lots more models that might match up with hypotheses we are interested in, but we are going to move forward and rank models. There are a few ways to do this. This is the easiest:
collect.models()

## This is fine if you want to include all the models in your AIC table. If you only want to select certain ones you can remove the ones you don't want to include.  We'll pretend that the weight+wing length model was misspcified and don't want to include it in our table. The number needs to match the model number in the lefthand column
aic<-collect.models()
aic
aic2<-remove.mark(aic,c(7))
aic2

## Looks like we have several competitive models. However the null model is the top model, which suggests that nothing is influencing survival more than considering constant survival. You could decide to model average the top models using the model.average() function. If you are in a situation like this, you should read more of the literature on the debates about model averaging

## As an example, we are still going to make some predictive plots. We will first look at how survival varies by monthly temperature. First we'll make our dataframe based on the min and max temp in any month. Although we don't have an interaction between time and temp, we need the column name to match the name of the original data. We will call it Temp2 even though the relationship with temperature is irrespective of month

newdat<-data.frame(Temp2=seq(from=min(Blackduck[,7:14]), to=max(Blackduck[,7:14]), by=1))

## This package uses a slightly different predict function. We still specify the model and new data, but we also need to specify the parameter index that we are using to predict - this is more important if we have a time varying or group varying model. In this case we only have the 1 covariate (temperature) so the index is 1 
temp_pred<-covariate.predictions(S_temp_mod, data=newdat, indices=c(1))

temp_graph<-temp_pred$estimates

ggplot(temp_graph, aes(y=estimate, x=covdata))+
  geom_line()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.3)+
  theme_classic()+
  ylab("Monthly survival")+
  xlab("Temperature")

## Now we will make a figure of how survival varies by week. This one is slightly easier because the model out put for our time model already has the estimates by week. 
time_graph<-S_time_mod$results$real

## We'll add a "week" label to make graphing easier
time_graph$month<-rep(1:8)

ggplot(time_graph, aes(y=estimate, x=month))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  theme_classic()+
  ylab("Monthly survival")+
  xlab("Month")

## The last thing we might need to know is cumulative survival, or survival during defined intervals. In this case, we estimated weekly survival but we might want to know annual survival

## If we look at our null model for example, we see a monthly survival of 0.937. To convert this to an annual survival, we'd raise it to the 12 power (12 months in a year) 
0.937^12

## Getting the error/confidence intervals is trickier because we just can't multiply standard errors. Instead we need to use the delta method. First we need to save the real estimates and the variance-covariance matrix

rr=get.real(S_dot_mod,"S",se=TRUE,vcv=TRUE)

## This gives us the SE associated with the annual survival
SE<-deltamethod(~x1^12, mean=rr$estimates$estimate[1], cov=rr$vcv.real, ses=T)

## We can then multiple the SE by 1.96 and subtract/add to get the CI
((rr$estimates$estimate[1])^12)-1.96*SE
((rr$estimates$estimate[1])^12)+1.96*SE

## So our annual survival is 0.46 (95% CI = 0.29 to 0.62)

## The delta methods gets trickier when we have different variances from different covariates. Let's calculate annual survival of males and females separately
rr2<-get.real(S_sex_mod,"S", se=TRUE, vcv=TRUE)

## If we look at our vcv matrix, it's now a 2x2 matrix
rr2$vcv.real

## Because females were the first level of our sex factor, the number in the [1,1] cell of the matrix is the variance associated with females. The [2,2] number is the variance associated with males and the other 2 numbers are related to covariance

SE_female<-deltamethod(~x1^12, mean=rr2$estimates$estimate[1], cov=rr2$vcv.real[1,1], ses=T)
((rr2$estimates$estimate[1])^12)-1.96*SE_female
((rr2$estimates$estimate[1])^12)+1.96*SE_female
((rr2$estimates$estimate[1])^12)


SE_male<-deltamethod(~x1^12, mean=rr2$estimates$estimate[9], cov=rr2$vcv.real[2,2], ses=T)
((rr2$estimates$estimate[9])^12)-1.96*SE_male
((rr2$estimates$estimate[9])^12)+1.96*SE_male
((rr2$estimates$estimate[9])^12)

## Using this method we can see that females have an annual survival rate of 0.42 (95% CI = 0.20 to 0.65) and males have an annual survival rate of 0.49 (95% CI = 0.25 to 0.73)


######################################################
#Practice#############################################
######################################################


condors<-import.chdata("known_fate_practice.txt", header=T, field.types = c("f", "n","n","n","n","f",rep("n", 11)))
## Let's take a look at the data
str(condors)
condors_processed<-process.data(condors, model="Known", groups=c("Sex"))
condors_ddl<-make.design.data(condors_processed)
condors_ddl

## Now we can run some models. For each model, we need to specify the formula with covariates for survival in a list. We'll start with the null model, which is often called the "dot" model
S_dot<-list(formula=~1)
S_dot_mod<-mark(data=condors_processed, ddl=condors_ddl, model.parameters = list(S=S_dot))
S_dot_mod

## Next we can run a model where survival depends on the bird's weight
S_weight<-list(formula=~Mass)
S_weight_mod<-mark(data=condors_processed, ddl=condors_ddl, model.parameters = list(S=S_weight))
S_weight_mod

## A model where survival differs between males and females
S_sex<-list(formula=~Sex)
S_sex_mod<-mark(data=condors_processed, ddl=condors_ddl, model.parameters = list(S=S_sex))

## This time if we look at the model output, we see we have a covariate for Sex:M. This means that female is our reference cateogory and the covariate tell us if male survival is different than females. The covariate is positive but the 95% CI overlaps zero, meaning that males have higher survival but not significantly different from females
S_sex_mod

