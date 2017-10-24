setwd("H:/Time to Event Analysis/Time to Event in R/Survival Analysis of island Eradication")
library(survival)

#Load in dataset:
erad<-read.csv("SurvBaselineV2.csv")

#Survival data and covariates:
yrs<-erad$TimeYrs #Data time (in years) from the year 1900 to present
event<-erad$Delta #Delta for event: 1= event(eradication), 0= non-event (eradication has not yet occurred)
area<-erad$Cov1_Area #covariate 1
latitude<-erad$Cov2_Latitude #covariate 2
erad.number<-(erad$Cov3_Erad_Number) #covariate 3
rodent.presence<- erad$Cov4_Rodent_pres #covariate 4
possum.presence<- erad$Cov5_Possum_pres #covariate 5
mustelid.presence<- erad$Cov6_Mustelid_pres #covariate 6
erad.method<- erad$Cov7_Erad_method #covariate 7

#Obtain a kaplan-meier estimate of the survival distribution:
resultkm.erad<-survfit(Surv(yrs,event)~1)

#Exrtract the survival estimates and the time variables:
survEst <-resultkm.erad$surv
survTime <-resultkm.erad$time

logLogsurvEst <- log(-log(survEst)) #log-log may not be necessary
logSurvTime <- log(survTime) #lopg transform may not be necessary
plot(logLogsurvEst ~ logSurvTime)
result.weibull.lm <- lm(logLogsurvEst ~logSurvTime)
abline(result.weibull.lm) #weibull distribution may not be appropriate for these data- points dont follow a linear relationship

## parametric survival distribution selection:
model.lognormal <- survreg(Surv(yrs, event)~erad.number + area + possum.presence + erad.method, dist="loglogistic")
summary(model.lognormal)
plot(model.lognormal)

erad.number + area + possum.presence + erad.method,  data=erad)