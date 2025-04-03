################################################################################
# R code used for the interrupted time series analysis of the DHIS headcount 
# data using ARIMA model
# author: Grace Carmichael
################################################################################

# read in packages needed for analysis
library(astsa)
library(forecast)
library(dplyr)
library(zoo)

# read in data (this is simulated data not the real dataset)
clinic_headcount <- read.csv("clinic_headcount.csv")

# remove COVID years from dataset
clinic_headcount_nocovid <- clinic_headcount[-(88:132),]

# convert data to time series object
clinic_ts <- ts(clinic_headcount_nocovid[,1], frequency=12, start=c(2013,1))
clinic_ts

# Plot data to visualise time series
plot(clinic_headcount_nocovid$headcount,type="n",
     ylim=c(min(clinic_headcount_nocovid$headcount)-50000,
            max(clinic_headcount_nocovid$headcount)+50000),
     xlab="Year", ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96),
     cex.lab=1.7, cex.axis=1.2)
# add shaded area to indicate when intervention started
rect(36,5000000,96,10306000,col=grey(0.9),border=F)
# add axis labels
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
# add data
lines(0:86,na.omit(clinic_ts))

# look at lag plot
lag.plot(clinic_ts, lags = 12, do.lines = FALSE)

# View ACF/PACF plots
acf2(clinic_ts, max.lag=24)
acf(clinic_ts, main="")
pacf(clinic_ts, main="")

# View ACF/PACF plots of differenced/seasonally differenced data
acf2(diff(diff(clinic_ts,12)), max.lag=24)
acf2(diff(clinic_ts,12), max.lag=24)

# Create variable representing step change and view
step <- as.numeric(as.yearmon(time(clinic_ts))>='Jan 2016')
step

# Create variable representing ramp (change in slope) and view
ramp <- append(rep(0,36), seq(1,51,1))
ramp  


# Use automated algorithm to identify p,q and d parameters
# Specify first difference = 1 and seasonal difference = 1
model_full <- auto.arima(clinic_ts, seasonal=TRUE, xreg=cbind(step,ramp), 
                         max.d=1, max.D=1, stepwise=FALSE, trace=TRUE)

# Check residuals
checkresiduals(model_full)
Box.test(model_full$residuals, lag = 24, type = "Ljung-Box")

res_ar <- as.matrix(model_full$residuals)
names(res_ar) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res_ar,pch=19,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date", type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.5)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(model_full$residuals, main="", cex.lab=1.7, cex.axis=1.5)
pacf(model_full$residuals, main="", cex.lab=1.7, cex.axis=1.5)
hist(model_full$residuals, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)


# Estimate parameters and confidence intervals
summary(model_full)
confint(model_full)

# get p-values for estimates
library(lmtest)
coeftest(model_full)


# get counterfactual by modelling data excluding the post-intervention time 
# period
model_preint <- Arima(window(clinic_ts, end=c(2015,12)), order=c(0,0,4), 
                      seasonal=list(order=c(1,1,0), period=12))

# forecast 51 months post-intervention and create ts object with forecast
fc <- forecast(model_preint, h=51)
fc.ts <- ts(as.numeric(fc$mean), start=c(2016,1), frequency=12)

# combine forcasted data with observed data
clinic_ts_2 <- ts.union(clinic_ts, fc.ts)
clinic_ts_2

# Check model fit
# recreate scatter plot
plot(clinic_headcount_nocovid$headcount,type="n",
     ylim=c(min(clinic_ts_2[,1])-1000,max(clinic_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96),cex.lab=1.7, cex.axis=1.2)
rect(36,min(clinic_ts_2[,1])-200000,96,max(clinic_ts_2[,1])+200000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(clinic_ts_2[,1]))
# add fitted line
lines(0:86,fitted(model_full), col = "#C11C06")


# plot counterfactual prediction
# recreate scatter plot
par(mar=c(5,5,3,2))
plot(clinic_headcount_nocovid$headcount,type="n",
     ylim=c(min(clinic_ts_2[,1])-1000,max(clinic_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96),cex.lab=1.7, cex.axis=1.2)
rect(36,min(clinic_ts_2[,1])-200000,96,max(clinic_ts_2[,1])+200000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020, cex.axis=1.4)
lines(0:86,na.omit(clinic_ts_2[,1]))
# add counterfactual line
lines(36:86,na.omit(clinic_ts_2[,2]), col="#C11C06")

