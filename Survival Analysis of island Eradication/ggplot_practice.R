install.packages("tidyverse")
library(tidyverse)

#ggplot with aesthetics mapping
ggplot(data=mpg) + geom_point (mapping = aes(x = displ, y= hwy, color = class)) #geom_point is used to create scatter plots in R- this works best with continuous variables

#ggplot with facets
ggplot(data=mpg) + geom_point(mapping = aes(x = displ, y= hwy)) + facet_wrap(~class,nrow =2) 

#overlaying Geom objects:
ggplot(data=mpg) + geom_point(mapping = aes(x = displ, y= hwy)) + geom_smooth(mapping = aes(x = displ, y= hwy))
#shorthand version with added scatterplot color
ggplot(data=mpg, mapping= aes(x = displ, y = hwy)) + geom_point(mapping = aes(color= class)) + geom_smooth()

#To select only subcompact cars:
ggplot(data=mpg, mapping= aes(x = displ, y = hwy)) + geom_point(mapping = aes(color= class)) + geom_smooth(data=filter(mpg, class=="subcompact"), se=F)

#Geom practice for line chart, boxplot, histogram, area chart:
#Boxplot:
ggplot(data=mpg) + geom_boxplot(mapping = aes(x = displ, y= hwy, color = class)) #

###################################################################################################################
# Example with my data:
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

#ggplot with aesthetics mapping
ggplot(data=erad) + geom_smooth (mapping = aes(x = log(area), y = erad.method))
