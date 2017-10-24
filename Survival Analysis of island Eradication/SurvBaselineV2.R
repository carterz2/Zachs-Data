setwd("H:/Time to Event Analysis/Time to Event in R/Survival Analysis of island Eradication")
library(survival)
library(moments)
library(muhaz)
library(forestplot)
library(ggplot2)
library(survminer)
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

#Relevel Covariates for proper comparisons
erad.method<- relevel(erad.method, "none")


###################################################################################################################
## Survival Curve estimation and plot:
Surv(yrs,event)
resultkm.erad<-survfit(Surv(yrs,event)~1)
summary(resultkm.erad)

#Survival curve:
plot(resultkm.erad, xlab="years", ylab="Invasive Pest Survival Probability (%)", ymin=0.7, conf.int=FALSE)

#Survival curve inverse to show eradication probability:
plot(resultkm.erad, fun="event", xlab="years", ylab="Invasive Pest Historical Eradication Probability (%)", conf.int=F)

###################################################################################################################
## Obtaining a smoothed hazard and survival function estimate:
#Hazard:
result.eradpe11<-pehaz(yrs, event, width=11, max.time=100) #ENSURE THAT THE BIN WIDTH IS CORRECT
plot(result.eradpe11, xlab="Time (years)", ylab="Hazard of Invasive Pest Eradication (offshore NZ islands)", main= "Hazard Plot")
result.eradpe11<-pehaz(yrs, event, width=11, max.time=100)
plot(result.eradpe11, xlab="Time (years)", ylab="Hazard of Invasive Pest Eradication (offshore NZ islands)", main= "Hazard Plot")
eradHazSmooth<-muhaz(yrs, event, bw.method="local", b.cor="none", max.time=100)
lines(eradHazSmooth, xlab="years", col="red")

#Obtaining Survival plot from Hazard function:
haz<-eradHazSmooth$haz.est
times<-eradHazSmooth$est.grid
surv<-exp(-cumsum(haz[1:length(haz)-1]*diff(times)))

#A smooth survival plot will allow for better identification of the distribution type that we may have to fit to the data
result.km<-survfit(Surv(yrs,event)~1, conf.type="plain")
plot(result.km, xlab= "time (years)", ylab="Invasive Pest Survival Probability (offshore NZ islands)", main="Survival Plot",ymin=0.7, conf.int=F)
lines(surv~times[1:(length(times)-1)], col="red")

###################################################################################################################
## Comparing the effect of individual covariates on the survival distribution using partial likelihood functions (the Cox proportional hazard model):
#Area
erad.area.cox<-coxph(Surv(yrs,event)~area, data=erad)
summary(erad.area.cox)
plot.cox.area<-survfit(erad.area.cox)
plot(plot.cox.area, ymin=0.7, col="red",conf.int=F)

#Latitude
erad.latitude.cox<-coxph(Surv(yrs,event)~latitude, data=erad)
summary(erad.latitude.cox)
plot.cox.latitude<-plot(survfit(erad.latitude.cox), col="green", ymin=0.7, conf.int=F)

#Number of eradication regimes
erad.number.cox<-coxph(Surv(yrs,event)~erad.number, data=erad)
summary(erad.number.cox)
plot.cox.number<-plot(survfit(erad.number.cox), col="purple", ymin=0.7, conf.int=F)

#Rodent presence
erad.rodent.cox<-coxph(Surv(yrs,event)~rodent.presence, data=erad)
summary(erad.rodent.cox)
plot.cox.rodent<-plot(survfit(erad.rodent.cox), col="brown", ymin=0.7, conf.int=F)

#Possum presence
erad.possum.cox<-coxph(Surv(yrs,event)~possum.presence, data=erad)
summary(erad.possum.cox)
plot.cox.possum<-plot(survfit(erad.possum.cox), col="grey", ymin=0.7, conf.int=F)

#Mustelid presence
erad.mustelid.cox<-coxph(Surv(yrs,event)~mustelid.presence, data=erad)
summary(erad.mustelid.cox)
plot.cox.mustelid<-plot(survfit(erad.mustelid.cox), col="orange", ymin=0.7, conf.int=F)

#Eradication Method
#Erad method really fucks with the model- it increases survival probability substantially
erad.method.cox<-coxph(Surv(yrs,event)~erad.method, data=erad)
summary(erad.method.cox)
plot.cox.eradmethod<- plot(survfit(erad.method.cox), col="cyan", ymin=0.7, conf.int=F)

#Eradication Method and area interaction: hypothesized that the size of the island effects the eradication method:
erad.method_area.cox<-coxph(Surv(yrs,event)~area:erad.method, data=erad)
summary(erad.method_area.cox)
#there is an effect- larger islands generally use mixed eradication methods as opposed to a single method

###################################################################################################################
## Assessing the porportional hazards assumption:
#Schoenfeld residuals for continuous variables (if the proportional hazards assumption is met, the residuals should look like a random scatter around zero):

# Area
erad.resid.schoenfeld <- cox.zph(coxph(Surv(yrs, event) ~area)) #the cox.zph function produces the schoenfeld residuals
erad.resid.schoenfeld #produces a p-value of 0.306- which is statistically insignificant. the covariate area does not fulfill the proportional hazards assumption.
plot(erad.resid.schoenfeld)

# Log10 Area
erad.resid.schoenfeld <- cox.zph(coxph(Surv(yrs, event) ~log10(area)))
erad.resid.schoenfeld 
plot(erad.resid.schoenfeld)

#Latitude and log(latitude)
erad.resid.schoenfeld <- cox.zph(coxph(Surv(yrs, event) ~log10(latitude)))
erad.resid.schoenfeld #produces a p-value of 0.219- which is statistically insignificant. Log(latitude) is also statistically insignificant (p=0.319) which means that latitude likely cannot be used as a predictor
plot(erad.resid.schoenfeld)


# parallel log cumulative hazard plots for categorical variables (curves should be approximately parallel and should not intersect after time apart):
    #because we are comparing the difference between categorical groups, we need to either a): compare the groups through the contruction of a stratified log-rank test
    # or b) include the other covariate (or multiple covariates) as regression terms for the hazard function

#I am confused here: I do not believe that there will be a confounding variable in this instance 
#(such as having two treatments- 0 and 1 - and we want to control (adjust) for a categorical covariate (such as gender)
#how to assess 

#Rodent presence:
rodent.surv.present <- survfit(Surv(yrs,event) ~rodent.presence, subset = {rodent.presence=="1"})
time.present <- rodent.surv.present$time
surv.present <- rodent.surv.present$surv
cloglog.present <- log10(-log10(surv.present)) #complimentary log-log survival: the y-output
logtime.present <-log10(time.present)
rodent.surv.non <- survfit(Surv(yrs,event) ~rodent.presence, subset = {rodent.presence=="0"})
time.non <- rodent.surv.non$time
surv.non <- rodent.surv.non$surv
cloglog.non <- log10(-log10(surv.non))
logtime.non <- log10(time.non)
plot (cloglog.present ~ logtime.present, type="s", col="blue", lwd=2)
lines (cloglog.non ~ logtime.non, type= "s", col="red", lwd=2)
legend("topleft", legend=c("present rodents", "no rodents"), col=c("blue", "red"), lwd=2) #doesnt seem to meet the the proportional hazards assumtion- the lines cross

#possum presence:
possum.surv.present <- survfit(Surv(yrs,event) ~possum.presence, subset = {possum.presence=="1"})
time.present <- possum.surv.present$time
surv.present <- possum.surv.present$surv
cloglog.present <- log10(-log10(surv.present)) #complimentary log-log survival: the y-output
logtime.present <-log10(time.present)
possum.surv.non <- survfit(Surv(yrs,event) ~possum.presence, subset = {possum.presence=="0"})
time.non <- possum.surv.non$time
surv.non <- possum.surv.non$surv
cloglog.non <- log10(-log10(surv.non))
logtime.non <- log10(time.non)
plot (cloglog.present ~ logtime.present, type="s", col="blue", lwd=2)
lines (cloglog.non ~ logtime.non, type= "s", col="red", lwd=2)
legend("bottomright", legend=c("present possums", "no possums"), col=c("blue", "red"), lwd=2)

#mustelid presence:
mustelid.surv.present <- survfit(Surv(yrs,event) ~mustelid.presence, subset = {mustelid.presence=="1"})
time.present <- mustelid.surv.present$time
surv.present <- mustelid.surv.present$surv
cloglog.present <- log10(-log10(surv.present)) #complimentary log-log survival: the y-output
logtime.present <-log10(time.present)
mustelid.surv.non <- survfit(Surv(yrs,event) ~mustelid.presence, subset = {mustelid.presence=="0"})
time.non <- mustelid.surv.non$time
surv.non <- mustelid.surv.non$surv
cloglog.non <- log10(-log10(surv.non))
logtime.non <- log10(time.non)
plot (cloglog.present ~ logtime.present, type="s", col="blue", lwd=2)
lines (cloglog.non ~ logtime.non, type= "s", col="red", lwd=2)
legend("bottomright", legend=c("present mustelids", "no mustelids"), col=c("blue", "red"), lwd=2) #doesnt seem to meet the the proportional hazards assumtion- the lines cross

#Erad methodology (none, hunt/trap, other, toxicant, toxicant and hunt/trap):
method.surv.none <- survfit(Surv(yrs,event) ~erad.method, subset = {erad.method=="none"})
time.none <- method.surv.none$time
surv.none <- method.surv.none$surv
cloglog.none <- log10(-log10(surv.none)) #cloglog= complimentary log-log survival: the y-output
logtime.none <-log10(time.none)
method.surv.huntTrap <- survfit(Surv(yrs,event) ~erad.method, subset = {erad.method=="hunt/trap"})
time.huntTrap <- method.surv.huntTrap$time
surv.huntTrap <- method.surv.huntTrap$surv
cloglog.huntTrap <- log10(-log10(surv.huntTrap))
logtime.huntTrap <- log10(time.huntTrap)
method.surv.other <- survfit(Surv(yrs,event) ~erad.method, subset = {erad.method=="other"})
time.other <- method.surv.other$time
surv.other <- method.surv.other$surv
cloglog.other <- log10(-log10(surv.other))
logtime.other <- log10(time.other)
method.surv.toxicant <- survfit(Surv(yrs,event) ~erad.method, subset = {erad.method=="toxicant"})
time.toxicant <- method.surv.toxicant$time
surv.toxicant <- method.surv.toxicant$surv
cloglog.toxicant <- log10(-log10(surv.toxicant))
logtime.toxicant <- log10(time.toxicant)
method.surv.both <- survfit(Surv(yrs,event) ~erad.method, subset = {erad.method=="toxicant and hunt/trap"})
time.both <- method.surv.both$time
surv.both <- method.surv.both$surv
cloglog.both <- log10(-log10(surv.both))
logtime.both <- log10(time.both)
#plot (cloglog.none ~ logtime.none, type="s", col="blue", lwd=2) #screws the graph up
plot (cloglog.huntTrap ~ logtime.huntTrap, type= "s", col="red", lwd=2)
# plot (cloglog.other ~ logtime.other, type= "s", col="green", lwd=2) #Not enough data to make a good plot
lines(cloglog.toxicant ~ logtime.toxicant, type= "s", col="purple", lwd=2)
lines (cloglog.both ~ logtime.both, type= "s", col="orange", lwd=2)
legend("bottomright", legend=c("none", "hunt/trap", "other", "toxicant"), col=c("blue", "red", "green", "purple", "orange"), lwd=2) #produces some odd results that definitely do not meet the proportional hazards assumption

#Number of eradication regimes:
eradNumber.surv.0 <- survfit(Surv(yrs,event) ~erad.number, subset = {erad.number=="0"})
time.0 <- eradNumber.surv.0$time
surv.0 <- eradNumber.surv.0$surv
cloglog.0 <- log10(-log10(surv.0)) 
logtime.0 <-log10(time.0)
eradNumber.surv.1 <- survfit(Surv(yrs,event) ~erad.number, subset = {erad.number=="1"})
time.1 <- eradNumber.surv.1$time
surv.1 <- eradNumber.surv.1$surv
cloglog.1 <- log10(-log10(surv.1)) 
logtime.1 <-log10(time.1)
eradNumber.surv.2 <- survfit(Surv(yrs,event) ~erad.number, subset = {erad.number=="2"})
time.2 <- eradNumber.surv.2$time
surv.2 <- eradNumber.surv.2$surv
cloglog.2 <- log10(-log10(surv.2))
logtime.2 <-log10(time.2)
plot (cloglog.0 ~ logtime.0, type= "s", col="red", lwd=2)
plot (cloglog.1 ~ logtime.1, type= "s", col="blue", lwd=2) #Can tell with just three of the different eradication regime numbers that it is not proportional
lines (cloglog.2 ~ logtime.2, type= "s", col="purple", lwd=2)

###################################################################################################################
# Use of the stratified log-rank test to determine if the proportional hazards assumpion is met with the categorical predictors:
# Needs to be done on our categorical predictors because we want to make statements compared to a baseline value (like male compared to female, for example)

survdiff(Surv(yrs,event)~rodent.presence)
rodent.surv.1 <- survfit(Surv(yrs,event) ~rodent.presence, subset = {rodent.presence=="1"})
rodent.surv.0 <- survfit(Surv(yrs,event) ~rodent.presence, subset = {rodent.presence=="0"})
plot(rodent.surv.1, col="red", conf.int=F, ylab= "Island pest survival probability (%)")
lines(rodent.surv.0, col="blue", conf.int= F)
legend ("bottomright", legend =c("present", "not present"), col= c("red", "blue"), lwd=2)

rodent.surv.strata <- survdiff(Surv(yrs,event)~rodent.presence +strata(erad.method))
lines(rodent.surv.strata, col="purple", conf.int=F)

###################################################################################################################

#A demonstration of the effect of island area on eradication method:

###################################################################################################################
##Using AIC to Determine the best model:
#Adding all covariates in to the same model for AIC analysis:
eradAll.coxph<-coxph(Surv(yrs,event)~ log10(area) + latitude + erad.number + rodent.presence + possum.presence + mustelid.presence +erad.method, data=erad)
summary(eradAll.coxph)
eradAll.surv<-plot(survfit(eradAll.coxph), ymin=0.75, col="red",conf.int=F)
erad.step<-step(eradAll.coxph, scope=eradAll.coxph, direction="both") #The step function reveals that latitude is a non-factor for predicting island pest survival in NZ

summary(erad.step)

#Visualising coefficient selection with a forest plot:
forestdata.rows <- matrix(c("Area", "Latitude", "#erads", "rodent presence", "possum presence", "mustelid presence"), ncol=1)
forestdata.coef.est <- data.frame(coef=c(-0.000056, -0.0283, 0.3560, 1.779, -0.8198, -0.6376),
                                low=c(0.9999, 0.9301, 1.2698, 2.7356, 0.1776, 0.2843),
                                high=c(1.0001, 1.0271, 1.6051, 12.8246, 1.0923, 0.9825))
forestplot(forestdata.rows, forestdata.coef.est,zero=0, lineheight="auto", clip=c(-3,3), boxsize=0.1, xticks= c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5),
          xlab="log hazard ratio",cex.lab=1.25, txt_gp=fpTxtGp(label=gpar(cex=0.85, lwd=1.5)))

#The forest plot demonstrates that, although latitude has a larger log hazard ratio (coef value) than some of the other covariates (namely mustelid/possum presence),
#its narrow confidence bands that trend toward zero indicate a small and ultimately marginal effect on predicting pest survival probability in NZ.


###################################################################################################################
## Model Diagnostics- Assessing goodness of fit using residuals:
# Plotting the martingale residuals against continuous covariates is a common approach used to determine nonlinearity
# or in other words, to asses the functional form of a covariate. Patterns in the plot may suggest that the variable is not properly fit

# Using Martingale residuals and a fit of the null Cox proportional hazards model:
result.null.coxph <- coxph(Surv(yrs, event)~ 1, data=erad)
residuals.null<- residuals(result.null.coxph, type = "martingale")


# Next, an assessment of the potential relationship between individual continuous covariates and survival.
# plot the null model residuals versus each variable used in the model:
# the below function (provided by Dirk Moore) fits a "loess" smooth curve through each set of residuals the better assess the shapes of the distribution (using the loes function)
#To get a 95% confidence interval for thew smooth curve the funtion "SmoothSEcurve" is generated:

smoothSEcurve<- function(yy,xx) {
  
  xx.list <- min(xx) + ((0:100)/100) * (max(xx)-min(xx))
  yy.xx <-predict(loess(yy ~ xx), se=T, newdata=data.frame(xx=xx.list))
  lines(yy.xx$fit ~ xx.list, lwd=2)
  lines(yy.xx$fit - qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
  lines(yy.xx$fit + qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
}

par (mfrow=c(1,2))
plot(residuals.null~log10(area))
smoothSEcurve(residuals.null, log10(area))
title("Martingale residuals\nversus log area")
#Another way using the "ggcoxfunctional" function in the survminer package
ggcoxfunc.area <- (ggcoxfunctional(Surv(yrs,event) ~ log10(area), data=erad))
ggcoxfunc.area
residuals.area <- plot(residuals.null~log10(area))
ggcoxfunctional(ggcoxfunc.area,  data = erad, point.col = "blue", point.alpha = 0.5)
plot(residuals.area,ggcoxfunc.area)
lines(residuals.area)

plot(residuals.null ~ latitude)
smoothSEcurve(residuals.null, latitude)
title("Martingale residuals\nversus latitude")

## An assessment of the final model by creating residuals of the complete model and then plotting them versus the final selected predictors:
residual.final <- residuals(erad.step, type="martingale")
par(mfrow=c(2,3))

plot (residual.final ~log10(area))
smoothSEcurve (residual.final, log10(area))
title("Martingale residuals \nvs area (log trans)")
#When compared to the null model output, the log transformed area variable demonstrates a much improved fit for the model.
#It still shows deviation from a linear-type model (hence the curvature) but the fit is improved over the plot for the null model (above).

#with eradication number considered as a continuous variable:
#plot(residual.final ~ erad.number)
#smoothSEcurve (residual.final, erad.number)
#title("Martingale residuals \nvs # of erads")

#with eradication number considered as a categorical variable:
boxplot( residual.final ~ erad.number, xlab="number of island eradications", ylab="residual.final", cex=0.6)
title("Martingale residuals \nvs # of erads") #try plotting with ggplot2 at somepoint

boxplot( residual.final ~ rodent.presence, xlab="Rodent presence \n(binary, 0=none, 1=present)", ylab="residual.final", cex=0.6)
title("Martingale residuals \nvs rodent presence")
#the residuals of rodent presence are not very evenly distributed- does this indicate a non-linear relationship?

boxplot( residual.final ~ possum.presence, xlab="possum presence \n(binary, 0=none, 1=present)", ylab="residual.final", cex=0.6)
title("Martingale residuals \nvs possum presence")

boxplot( residual.final ~ mustelid.presence, xlab="mustelid presence \n(binary, 0=none, 1=present)", ylab="residual.final", cex=0.6)
title("Martingale residuals \nvs mustelid presence")

boxplot(residual.final ~ erad.method, xlab="eradication method", ylab="residual.final", cex=0.6)
title("Martingale residuals \nvs eradication method")

# Analysis of case deletion residuals (jackknife residuals):
# The jackknife residual plot will determine if any subject within the covariate is having an especially large influence on the parameter estimate:

coef.all<- erad.complete.coxph$coef[2]

n.obs<- length(yrs)
jkbeta.vec <- rep(NA, n.obs) 
for ( i in 1:n.obs) {
  yrs.i <- yrs [-i]
  delta.i <-event [-i]
  area.i <- area [-i]
  erad.number.i <- erad.number [-i]
  rodent.i <- rodent.presence [-i]
  possum.i <- possum.presence [-i]
  mustelid.i <- mustelid.presence [-i]
  result.coxph.i <-coxph (Surv(yrs.i, delta.i) ~area.i + erad.number.i +rodent.i +possum.i +mustelid.i)
  coef.i <- result.coxph.i$coef[2]
  jkbeta.vec[i] <- (coef.all -coef.i)
}

index.obs <- 1:n.obs
plot(jkbeta.vec ~index.obs, type="h", xlab = "observation number", ylab = "change in coefficient for the number of erads", cex.axis=1.3, cex.lab=1.3)
abline(h=0)
#identify(jkbeta.vec ~index.obs) #identified observation numbers 160, 194, and 207 as possibly influencing the results because of their large magnitude

###################################################################################################################
#Plotted demonstration on the effect of eradication method on pest survival probability:
#Toxicant survival curve:
#Hunt/trap survival curve:
#combination eradication method survival curve:

###################################################################################################################
#Assessment of the proportional hazards assumption