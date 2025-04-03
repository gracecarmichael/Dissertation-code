
################################################################################
# R code used for the interrupted time series analysis of the DHIS headcount data 
# author: Grace Carmichael
################################################################################

# load the packages
library(foreign)
library(tsModel)
library(Epi)
library(splines)
library(vcd)
library(nlme)
library(car)
library(AER)
library(tseries)
library(DHARMa)
library(Metrics)

# read in dataset (this is a simulated dataset, the actual dataset is stored
# on Box)

clinic_headcount <- read.csv("clinic_headcount.csv")

# This dataset includes the following variables:
# year
# month
# time = time elapsed since start of dataset
# headcount = country-wide (excluding Western Cape province) clinic headcounts
# ccmdd = indicator for when the ccmdd programme was in place (0 = not in place,
# 1 = in place)
# covid = indicator for pre (0) versus post (1) the start of the COVID-19 
# pandemic

clinic_headcount_nocovid <- clinic_headcount[-(88:132),]

#-------------------------------------------------------------------------------
# Exploratory Data Analysis
#-------------------------------------------------------------------------------

### Scatter plot

## With COVID-19 period
# plot axes
plot(clinic_headcount$headcount,type="n",
     ylim=c(min(clinic_headcount$headcoun)-50000,
            max(clinic_headcount$headcoun)+50000),
     xlab="Year", ylab="Headcount", bty="l", xaxt="n", xlim = c(0,132))
# shade the post ccmdd initiation period
rect(36,5000000,132,10306000,col=grey(0.9),border=F)
# add axis labels
axis(1,at=0:11*12,labels=F)
axis(1,at=0:10*12+6,tick=F,labels=2013:2023)
# add data
lines(0:131,clinic_headcount$headcount)

## Without COVID-19 period
# plot axes
plot(clinic_headcount_nocovid$headcount,type="n",
     ylim=c(min(clinic_headcount_nocovid$headcount)-50000,
            max(clinic_headcount_nocovid$headcount)+50000),
     xlab="Year", ylab="Headcount",bty="l",xaxt="n",xlim = c(0,96))
# shade the post ccmdd initiation period
rect(36,5000000,96,10306000,col=grey(0.9),border=F)
# add axis labels
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020)
# add data
lines(0:86,clinic_headcount_nocovid$headcount)



# Produce summary statistics
summary(clinic_headcount)

# look at headcount summaries before and after ccmdd introduction
summary(clinic_headcount$headcount[clinic_headcount$ccmdd==0])
summary(clinic_headcount$headcount[clinic_headcount$ccmdd==1])

# look at headcount summaries before and after ccmdd introduction (before 
# beginning of covid-19 pandemic)
summary(clinic_headcount$headcount[clinic_headcount$ccmdd==0])
summary(clinic_headcount$headcount[clinic_headcount$ccmdd==1 & 
                                     clinic_headcount$covid == 0])

# check for trend stationarity (before ccmdd introduction)
kpss.test(clinic_headcount[1:36,1], null = "Trend")

#-------------------------------------------------------------------------------
# Poisson regression model (with COVID-19 period)
#-------------------------------------------------------------------------------

# data is count data (headcounts) so start with Poisson regression model

# Poisson model
model_pois <- glm(headcount ~  ccmdd + time + covid, family=poisson, 
                  data = clinic_headcount)
summary(model_pois)

# create a new dataframe with 0.1 time units to improve the graph
datanew <- data.frame(ccmdd=rep(c(0,1),c(360,960)), 
                      covid = rep(c(0,1),c(870,450)),time= 1:1320/10,
                      month=rep(1:120/10,11))

# generate predicted values for plotting
pred1 <- predict(model_pois,type="response",datanew)

# recreate scatter plot
plot(clinic_headcount$headcount,type="n",
     ylim=c(min(clinic_headcount$headcoun)-50000,
            max(clinic_headcount$headcoun)+50000),
     xlab="Year", ylab="Headcount",bty="l",xaxt="n")
rect(36,5000000,132,9306000,col=grey(0.9),border=F)
points(clinic_headcount$headcount,cex=0.7)
axis(1,at=0:11*12,labels=F)
axis(1,at=0:10*12+6,tick=F,labels=2013:2023)
# add fitted line
lines((1:1320/10),pred1,col=2)


# plot the counterfactual scenario by creating a data frame as if the ccmdd
# (the intervention) was never introduced
datanew_count <- data.frame(ccmdd=rep(0,1320), covid = rep(c(0,1),c(870,450)),
                      time= 1:1320/10,month=rep(1:120/10,11))

# generate predictions under the counterfactual scenario and add line to the 
# scatterplot
pred_count <- predict(model_pois,datanew_count,type="response")
lines(datanew$time,pred_count,col=2,lty=2)


#-------------------------------------------------------------------------------
# Checking Model Fit and Possible Data Issues
#-------------------------------------------------------------------------------

# 1. Overdispersion: Quasi-Poisson model 

# conduct overdispersion test
dispersiontest(model_pois)

# generate plots to check residuals
sim_model_pois <- simulateResiduals(model_pois, refit=T)
plotSimulatedResiduals(sim_model_pois)

# fit quasipoisson model to deal with overdispersion
model_quasi <- glm(headcount ~  ccmdd + time + covid, family=quasipoisson, 
              data = clinic_headcount)
summary(model_quasi)



# 2. Checking residuals for autocorrelation

# Check the residuals by plotting against time
res2 <- residuals(model_quasi,type="deviance")
plot(clinic_headcount$time,res2,pch=19,cex=0.7,col=grey(0.6),
     main="Residuals over time",ylab="Deviance residuals",xlab="Time")
abline(h=0,lty=2,lwd=2)

hist(res2, main="", breaks = 10)

# Further check for autocorrelation by examining the autocorrelation and
# partial autocorrelation functions
acf(res2)
pacf(res2)

# conduct Durbin Watson test for autocorrelation
durbinWatsonTest(model_quasi, max.lag=5)
dwtest(model_quasi, alternative = "two.sided")


# 3. Adjusting for seasonality
# There are a number of ways to adjusting for seasonality - here harmonic
# terms are used 
# the length of the period is 12 months

model3a <- glm(headcount ~  ccmdd + time + covid + harmonic(month,2,12), 
               family=quasipoisson, data = clinic_headcount)
model3b <- glm(headcount ~  ccmdd + time + covid + harmonic(month,3,12), 
               family=quasipoisson, data = clinic_headcount)
model3c <- glm(headcount ~  ccmdd + time + covid + harmonic(month,4,12), 
               family=quasipoisson, data = clinic_headcount)
model3d <- glm(headcount ~  ccmdd + time + covid + harmonic(month,5,12), 
               family=quasipoisson, data = clinic_headcount)

anova(model3a,model3b,model3c,model3d)

summary(model3c)


# Check residuals for autocorrelation
res3 <- residuals(model3c,type="deviance")
plot(res3,pch=19,cex=0.7,col=grey(0.6),main="Residuals over time",
     ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2)
acf(res3)
pacf(res3)

# conduct Durbin Watson test for autocorrelation
durbinWatsonTest(model3c, max.lag=5)
dwtest(model3c, alternative = "two.sided")


#-------------------------------------------------------------------------------
# Adding a change in slope
#-------------------------------------------------------------------------------

# add a change-in-slope
model4 <- glm(headcount ~  ccmdd * time + covid*time + harmonic(month,4,12), 
              family=quasipoisson, data = clinic_headcount)
summary(model4)


# Check residuals for autocorrelation
res4 <- residuals(model4,type="deviance")
plot(res4,pch=19,cex=0.7,col=grey(0.6),main="Residuals over time",
     ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2)
acf(res4)
pacf(res4)

# conduct Durbin Watson test for autocorrelation
durbinWatsonTest(model4, max.lag=5)
dwtest(model4, alternative = "two.sided")


# generate fitted values from
pred4 <- predict(model4,type="response",transform(datanew,month=6))
pred3c <- predict(model3c,type="response",transform(datanew,month=6))
# recreate scatter plot
plot(clinic_headcount$headcount,type="n",
     ylim=c(min(clinic_headcount$headcoun)-50000,
            max(clinic_headcount$headcoun)+50000),
     xlab="Year", ylab="Headcount",bty="l",xaxt="n")
rect(36,5000000,132,9306000,col=grey(0.9),border=F)
points(clinic_headcount$headcount,cex=0.7)
axis(1,at=0:11*12,labels=F)
axis(1,at=0:10*12+6,tick=F,labels=2013:2023)
lines(1:1320/10,pred4,col=4)
lines(1:1320/10,pred3c,col=2)
legend("topright",c("Step-change only","Step-change + change-in-slope"),lty=1,
       col=c(2,4),inset=0.05,bty="n",cex=0.7)

# test if the change-in-slope improved the fit
anova(model3c,model4,test="F")

# recreate scatter plot
plot(clinic_headcount$headcount,type="n",ylim=c(5000000,9306000),xlab="Year", 
     ylab="Headcount",bty="l",xaxt="n")
rect(36,5000000,132,9306000,col=grey(0.9),border=F)
axis(1,at=0:11*12,labels=F)
axis(1,at=0:10*12+6,tick=F,labels=2013:2023)
lines(1:132,clinic_headcount$headcount)

# create counterfactual data (under scenario that CCMDD not introduced)
datanew_count <- data.frame(ccmdd=rep(0,1320), covid = rep(c(0,1),c(870,450)),
                            time= 1:1320/10,month=rep(1:120/10,11))

# generate predictions under the counterfactual scenario and add line to the 
# scatterplot
pred_count <- predict(model4,datanew_count,type="response")
lines(datanew$time,pred_count,col=2,lty=2)


#-------------------------------------------------------------------------------
# Model without COVID-19 period
#-------------------------------------------------------------------------------

# fit same model as model 4 but without the COVID-19 period in the data
model5 <- glm(headcount ~  ccmdd * time + harmonic(month,4,12), 
              family=quasipoisson, data = clinic_headcount_nocovid)
summary(model5)


# Check residuals for autocorrelation
res5 <- residuals(model5,type="response")
names(res5) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res5,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.2)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(res5, main="", cex.lab=1.7, cex.axis=1.5)
pacf(res5, main="", cex.lab=1.7, cex.axis=1.5)
hist(res5, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)

# conduct durbin watson test for autocorrelation
durbinWatsonTest(model5, max.lag=5)
dwtest(model5, alternative = "two.sided")

# create data for prediction
datanew <- data.frame(ccmdd=rep(c(0,1),c(360,510)),time= 1:870/10,
                      month=c(rep(1:120/10,7),1:30/10))

# recreate scatter plot
plot(clinic_headcount_nocovid$headcount,type="n",
     ylim=c(min(clinic_headcount_nocovid$headcount)-50000,
            max(clinic_headcount_nocovid$headcount)+50000),
     xlab="Year", ylab="Headcount",bty="l",xaxt="n",xlim = c(0,96))
rect(36,5000000,96,10306000,col=grey(0.9),border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020)
lines(0:86,clinic_headcount_nocovid$headcount)
# add fitted line
lines(0:86,fitted(model5),col="#C11C06")

# create counterfactual data (under scenario that CCMDD not introduced)
datanew_count <- data.frame(ccmdd=rep(0,87),
                            time= 1:87,month=c(rep(1:12,7),1:3))

# generate predictions under the counterfactual scenario and add line to the 
# scatterplot
pred_count <- predict(model5,datanew_count,type="response")
lines(36:86,pred_count[-c(1:36)],col="#C11C06")


#-------------------------------------------------------------------------------
# GLS model
#-------------------------------------------------------------------------------

# fit same model as model 5 but including correlation structure
mod.gls.har <- gls(headcount ~  ccmdd * time + harmonic(month,4,12),
         data = clinic_headcount_nocovid,
         correlation=corARMA(p=1, form=~time), method = "ML")
summary(mod.gls.har)

# include month as a factor to avoid needing to specify harmonic terms
clinic_headcount_nocovid$month <- as.factor(clinic_headcount_nocovid$month)
mod.gls <- gls(headcount ~  ccmdd * time + month,
               data = clinic_headcount_nocovid,
               correlation=corARMA(p=1, form=~time), method = "ML")
summary(mod.gls)

# check residuals for autocorrelation
res6 <- residuals(mod.gls,type="response")
names(res6) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res6,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n")
axis(1,at=0:7*12,labels = 2013:2020)
abline(h=0,lty=2,lwd=2)
acf(res6, main="")
pacf(res6, main="")

# recreate scatter plot
plot(clinic_headcount_nocovid$headcount,type="n",
     ylim=c(min(clinic_headcount_nocovid$headcount)-50000,
            max(clinic_headcount_nocovid$headcount)+50000),
     xlab="Year", ylab="Headcount",bty="l",xaxt="n",xlim = c(0,96))
rect(36,5000000,96,10306000,col=grey(0.9),border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020)
lines(0:86,clinic_headcount_nocovid$headcount)
# add fitted line
lines(0:86,fitted(mod.gls),col="#C11C06")


# Check different correlation structures by varying p and q
AIC_vec <- matrix(rep(0,49), ncol = 7, nrow = 7)

for(p in 1:7){
  for(q in 1:7){
    tryCatch({mod <-update(mod.gls, correlation=corARMA(p=p-1,q=q-1,form=~time))
    AIC_vec[p,q] <- AIC(mod)}, error=function(e){AIC_vec[p,q] <- NA})
    
  }
}

# look at AIC vector to see which correlation structure gave smallest AIC
AIC_vec

mod.gls.4 <- update(mod.gls, correlation=corARMA(p=4,q=4,form=~time))

# Check residuals 
res6 <- residuals(mod.gls.4,type="response")
names(res6) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res6,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.2)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(res6, main="", cex.lab=1.7, cex.axis=1.5)
pacf(res6, main="", cex.lab=1.7, cex.axis=1.5)
hist(res6, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)


# recreate scatter plot
plot(clinic_headcount_nocovid$headcount,type="n",
     ylim=c(min(clinic_headcount_nocovid$headcount)-50000,
            max(clinic_headcount_nocovid$headcount)+50000),
     xlab="Year", ylab="Headcount",bty="l",xaxt="n",xlim = c(0,96))
rect(36,5000000,96,10306000,col=grey(0.9),border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020)
lines(0:86,clinic_headcount_nocovid$headcount)
# add fitted line
lines(0:86,fitted(mod.gls.4),col="#C11C06")

# create counterfactual data
datanew_count$month <- as.factor(datanew_count$month)
pred_count.gls <- predict(mod.gls.4 ,datanew_count,type="response")
attr(pred_count.gls,"label") <- NULL 
# add counterfactual line to plot
lines(36:86,pred_count.gls[-c(1:36)],col="#C11C06")


# Get MAE and RMSE for final model
mae(clinic_headcount_nocovid$headcount,fitted(mod.gls.4))
rmse(clinic_headcount_nocovid$headcount,fitted(mod.gls.4))

