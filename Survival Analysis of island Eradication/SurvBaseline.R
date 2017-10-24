# Basic test of the Product-Limit Estimation Method with the given eradication data

# Analysis only includes uncensored or right censored data due to the nature of the dataset (uncensored if eradication is known, 
#right-censored if: eradication has yet to occur or if natural extinction occured)

#The product of this dataset is only practice:

setwd("H:/Time to Event Analysis/Time to Event in R/Survival Analysis of island Eradication")
library(survival)

erad<-read.csv("SurvBaseline.csv")
yrs<-erad$timeYrs
event<-erad$Delta
area<-erad$Area_ha_Cov1
# With an original skew of 4.47 and a kurtosis of 19.99, the area frequency was far from normal (skewed very var right: the median is 13 with a max of 168537) --> refitted the area on a log_10 scale to get a new skew of 1.62 and a kurtosis of 2.42. Its not great but its much better than previous. No need to calculate with reflection either due to positive skew
area.log<-erad$Area_Cov1_LogTrans
latitude<-erad$Latitude_Cov2
Eradnumb<-erad$Erad_Numb_Cov3
Eradnumb.logx1<-erad$EradNumb_Cov2_logx.1

#############################################################################################
#Survival curve estimation and plot:

Surv(yrs,event)
#resultkm.erad<-survfit(Surv(yrs,event)~1, conf.type="log-log") ## testing of log-log confidence intervals
resultkm.erad<-survfit(Surv(yrs,event)~1) #the "~" represents a single curve for all erads within the datset
summary(resultkm.erad)
resultkm.erad #the result shows that there is no median value or 0.95LCL/UCL -> this is because the survival confidence never goes below 0.5 for the probability 

plot(resultkm.erad, xlab="years", ylab="Invasive Pest Survival Probability (%)", ymin=0.7, conf.int=F)
#axis(2, at=seq(0,1,by=.2),labels=paste(100*seq(0,1,by=.2), "%") ) -> change axis to percentage

#Plotting the inverse of the survival curve to show the probability of eradication occurring-> very similar to the amount of eradicated area graph produced for BHNSC workshop
plot(resultkm.erad, fun="event", xlab="years", ylab="Invasive Pest Historical Eradication Probability (%)", conf.int=F)

#############################################################################################
# Obtaining a smoothed hazard and survival function estimate:

#Hazard:
library(muhaz)
result.eradpe11<-pehaz(yrs, event, width=11, max.time=100) #ENSURE THAT THE BIN WIDTH IS CORRECT
plot(result.eradpe11, xlab="Time (years)", ylab="Hazard of Invasive Pest Eradication (offshore NZ islands)", main= "Hazard Plot")

#result.eradpe1<-pehaz(yrs, delta, width=5, max.time=100)
#lines(result.eradpe1)

eradHazSmooth<-muhaz(yrs, event, bw.method="local", b.cor="none", max.time=100)
lines(eradHazSmooth, xlab="years", col="red")

#Obtaining Survival plot from Hazard function:
haz<-eradHazSmooth$haz.est
times<-eradHazSmooth$est.grid
surv<-exp(-cumsum(haz[1:length(haz)-1]*diff(times)))

#A smooth survival plot will allow for better identification of the distribution type that we may have to fit to the data
result.km<-survfit(Surv(yrs,event)~1, conf.type="plain")
plot(result.km, xlab= "time (years)", ylab="Invasive Pest Survival Probability (offshore NZ islands)", main="Survival Plot",ymin=0.6)
lines(surv~times[1:(length(times)-1)], col="red")

############################################################################################
#Comparing the effect of different covariates on the survival distribution using partial likelihood functions (the Cox proportional hazard model):
#EVENTUALLY RUN A SENSITIVITY ANALYSIS TO SEE WHAT COVARIATES ARE MOST EFFECTED BY CHANGE

#Inquiry of the effect of area (log_10 transformed) on the survival dataset:
erad.area.cox<-coxph(Surv(yrs,event)~area.log, data=erad)
summary(erad.area.cox)


#The positive coef value means that the hazard (the instantaneous risk of failure, where failure is eradication in this case) is higher with the inclusion of area; that is invasive pest survival slightly decreases with the inclusion of eradicated area -> as demonstrated by the slightly higher survival output below --> the probability of eradication slightly increases
#In other words, increasing the area coefficient by one unit (one log unit) increases the the the log hazard ratio by 0.3296.
#Interpretation of the coef value-> even 


##THE COVARIATE ANALYSIS IS SHOWING AN INCREASE IN IN THE HAZARD RATIO -> THAT IS, AN INCREASE IN THE 



#increasing the hazard increases the eradication probability

#Summary results estimate a log hazard ratio (coef value) of 0.32962. The positive coef means that the hazard (risk of invasive pest survival) is higher; That is, invasive pest survival increases with the inclusion of eradicated area.
#In other words, increasing the area coefficient by one unit (one log unit b/c "area" was log transformed) will increase the log hazard ratio by 0.3296.
#-->Interpretation of coef value: increasing island size makes it more difficult to successfully eradicate invasive pests, as a result inclusion of the coefficient "area" will increase pest survival probability.

#exp(coef)--> the hazard ratio --> the effect size of the covariate: with a hazard ratio of 1.39, there is a 39% increase in hazard through an increase in the area coefficient by one (log) unit.
#--> interpreation of hazard ratio: there is a 39% increase in pest survival with an increase in one log area unit --> is there a way to transform this back to get sensical results? Understanding how increasing island area is affecting eradication success would be really interesting.

#The p-value is statistically significant (p=0.0004, alpha=0.05): we reject the null hypothesis that there is no difference between the calculated coefficient and 0--> there is a difference: So this model is "moving" the survival curve.


#logrank test tests a null hypothesis where there is no difference between the populations in the probability of an event at any time point. Tests the null hypothesis of no difference in survival between two (or more) independent groups.
# The log rank test statistic is approx distributed as chi-square with 1 degree of freedom -> thus critical value for the test can be found in a table of critical values of the x^2 distribution (test value with alpha=0.05 is 3.841)--> we reject the null and have significant evidence to show that the survival curves are different from one another.
#if the survival curves cross, the log rank test may become null--> these curves do not cross, though they decrease in the same stepwise fashion.

cox.zph(erad.area.cox)
plot.cox.logarea<-survfit(erad.area.cox)
plot(plot.cox.logarea, ymin=0.7, col="red",conf.int=F)

#Inquiry of the effect of latitude on the survival dataset:
#latitude skew of 0.9 and kurtosis of 0.05 --> no need to transform the data
erad.latitude.cox<-coxph(Surv(yrs,event)~latitude, data=erad)
summary(erad.latitude.cox)
cox.zph(erad.latitude.cox)
plot.cox.latitude<-plot(survfit(erad.latitude.cox), col="green", ymin=0.7, conf.int=F)


#Inquiry of the effect of the number of eradicated Species for each island (includes repeat eradications):
erad.Eradnumb.cox<-coxph(Surv(yrs,event)~Eradnumb.logx1, data=erad)
summary(erad.Eradnumb.cox)
plot.cox.Eradnumb<-plot(survfit(erad.Eradnumb.cox), col="purple", ymin=0.7, conf.int=F)

#graph all survival lines together:
lines(resultkm.erad, conf.int=F)
lines(plot.cox.logarea, col="red", conf.int=F)
lines(plot.cox.Eradnumb, col="purple",conf.int=F)

############################################################################################
###Using AIC to Determine the best model:

#Adding all covariates in to the same model for AIC analysis:
eradAll.coxph<-coxph(Surv(yrs,event)~latitude + area.log+ Eradnumb.logx1, data=erad)
eradAll.surv<-plot(survfit(eradAll.coxph), ymin=0.75, col="red",conf.int=F)

#Step function to determine the best model via AIC for non-nested models:

#Step function currently determines that removing both latitude and island area are the best choice for the model (assuming that area or latitude covariates are fixed within the model)
#Remember that the AIC formula considers the number of parameters to be a penalty term- the larger the number of parameters the more penalised the model is.
#Output determines that only the number of eradicated species affects the model output
erad.step<-step(eradAll.coxph, scope=list(upper= ~area.log +latitude+Eradnumb.logx1, lower= ~Eradnumb.logx1),direction="both")
erad.step<-step(eradAll.coxph, scope=eradAll.coxph, direction="both") #there is an error message but its no big deal- the output is the same as the previous step function approach

#Visualising the Model selection with a Forest Plot: