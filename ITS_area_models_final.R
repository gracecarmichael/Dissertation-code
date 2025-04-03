################################################################################
# R code used for a zoom in on the differences in headcounts for different areas
# author: Grace Carmichael
################################################################################


library(lmtest)
library(nlme)

#-------------------------------------------------------------------------------
# Headcount rural vs peri-urban vs urban
#-------------------------------------------------------------------------------

# seperate headcounts into the 3 area types
rural_clinics <- DHIS_clinics[which(DHIS_clinics$OrgUnitRuralUrban=="Rural"),]
periurban_clinics <- DHIS_clinics[which(DHIS_clinics$OrgUnitRuralUrban=="Peri-Urban"),]
urban_clinics <- DHIS_clinics[which(DHIS_clinics$OrgUnitRuralUrban=="Urban"),]

rural_headcount <- clinics_all[which(clinics_all$organisationunitname %in% rural_clinics$`Facility Name`),]
periurban_headcount <- clinics_all[which(clinics_all$organisationunitname %in% periurban_clinics$`Facility Name`),]
urban_headcount <- clinics_all[which(clinics_all$organisationunitname %in% urban_clinics$`Facility Name`),]

# get total headcounts for each month for each of the 3 area types
rural_headcount_total <- apply(rural_headcount[-1],2,sum, na.rm=TRUE)
periurban_headcount_total <- apply(periurban_headcount[-1],2,sum, na.rm=TRUE)
urban_headcount_total <- apply(urban_headcount[-1],2,sum, na.rm=TRUE)

# plot headcounts over time for each area type
par(mar=c(8,4,4,2))
plot(y=rural_headcount_total,x=1:132, type = "l", ylab = "headcount", xaxt="n", xlab = "")
axis(1, at = 1:132, labels = names(clinics_total), cex.axis = 0.75, las=2)

par(mar=c(8,4,4,2))
plot(y=periurban_headcount_total,x=1:132, type = "l", ylab = "headcount", xaxt="n", xlab = "")
axis(1, at = 1:132, labels = names(clinics_total), cex.axis = 0.75, las=2)

par(mar=c(8,4,4,2))
plot(y=urban_headcount_total,x=1:132, type = "l", ylab = "headcount", xaxt="n", xlab = "")
axis(1, at = 1:132, labels = names(clinics_total), cex.axis = 0.75, las=2)


#-------------------------------------------------------------------------------
# Rural model
#-------------------------------------------------------------------------------

### ARIMA model ###

# create dataset using monthly headcounts
rural_headcount <- data.frame(headcount = rural_headcount_total, 
                               year = rep(c(2013,2014,2015,2016,2017,2018,2019,2020,2021,
                                            2022,2023), each = 12), 
                               month = rep(c(1,2,3,4,5,6,7,8,9,10,11,12), times = 11),
                               time = 1:132,
                               ccmdd = c(rep(0, times = 36), rep(1, times=96)),
                               covid = c(rep(0, times = 87), rep(1, times = 45)),
                               row.names = 1:132)

# remove covid years from dataset
rural_nocovid <- rural_headcount[-(88:132),]

# convert data to time series object
rural_ts <- ts(rural_nocovid[,1], frequency=12, start=c(2013,1))
rural_ts

# check lag plot
lag.plot(rural_ts, lags = 12, do.lines = FALSE)

# look at ACF/PACF plots
acf2(rural_ts, max.lag=24)

# look at ACF/PACF plots of differenced/seasonally differenced data
acf2(diff(diff(rural_ts,12)), max.lag=24)

# create variable representing step change
step <- as.numeric(as.yearmon(time(rural_ts))>='Jan 2016')
step

# create variable representing change in slope
ramp <- append(rep(0,36), seq(1,51,1))
ramp  


# Use automated algorithm to identify p,q and d parameters
rural_full <- auto.arima(rural_ts, seasonal=TRUE, xreg=cbind(step,ramp), 
                         max.d=1, max.D=6, max.P = 6, stepwise=FALSE, trace=TRUE)


# check residuals for autocorrelation
checkresiduals(rural_full)
Box.test(rural_full$residuals, lag = 24, type = "Ljung-Box")

res_rural <- as.matrix(rural_full$residuals)
names(res_rural) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                   rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))

par(mar=c(5,5,3,2))
plot(res_rural,pch=19,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date", type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.5)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(rural_full$residuals, main="", cex.lab=1.7, cex.axis=1.5)
pacf(rural_full$residuals, main="", cex.lab=1.7, cex.axis=1.5)
hist(rural_full$residuals, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)


# check model fit
par(mar=c(5,5,3,2))
plot(rural_nocovid$headcount,type="n",ylim=c(min(rural_ts)-1000,max(rural_ts)+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96),cex.lab=1.7, cex.axis=1.2)
rect(36,min(rural_ts)-200000,96,max(rural_ts)+200000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(rural_ts))
lines(0:86,fitted(rural_full), col = "#C11C06")


# get parameter estimates and confidence intervals
summary(rural_full)
round(confint(rural_full),2)
coeftest(rural_full)

# get counterfactual by modelling data excluding the post-intervention time 
# period
rural_preint <- Arima(window(rural_ts, end=c(2015,12)), order=c(0,0,4), 
                      seasonal=list(order=c(1,1,0), period=12))

# forecast 51 months post-intervention and create ts object with forecast
rural_fc <- forecast(rural_preint, h=51)
rural_fc_ts <- ts(as.numeric(rural_fc$mean), start=c(2016,1), frequency=12)

# combine forcasted data with observed data
rural_ts_2 <- ts.union(rural_ts, rural_fc_ts)
rural_ts_2

# cretae scatterplot of headcounts
par(mar=c(5,5,3,2))
plot(rural_nocovid$headcount,type="n",
     ylim=c(min(rural_ts_2[,1])-1000,max(rural_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(rural_ts_2[,1])-100000,96,max(rural_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(rural_ts_2[,1]))
# add counterfactual line
lines(36:86,na.omit(rural_ts_2[,2]), col="#C11C06")



### GLS regression model ###

# convert month to a factor variable
rural_nocovid$month <- as.factor(rural_nocovid$month)

# fit model with gls
mod.gls <- gls(headcount ~  ccmdd * time + month,
               data = rural_nocovid,
               correlation=corARMA(p=1, form=~time), method = "ML")
summary(mod.gls)

# check residuals
res <- residuals(mod.gls,type="response")
names(res) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n")
axis(1,at=0:7*12,labels = 2013:2020)
abline(h=0,lty=2,lwd=2)
acf(res, main="")
pacf(res, main="")

# check model fit
par(mar=c(5,5,3,2))
plot(rural_nocovid$headcount,type="n",
     ylim=c(min(rural_ts_2[,1])-1000,max(rural_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7)
rect(36,min(rural_ts_2[,1])-100000,96,max(rural_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(rural_ts_2[,1]))
lines(0:86,fitted(mod.gls), col="#C11C06")


# try different correlation structures
AIC_vec <- matrix(rep(0,49), ncol = 7, nrow = 7)

for(p in 1:7){
  for(q in 1:7){
    tryCatch({mod <-update(mod.gls, correlation=corARMA(p=p-1,q=q-1,form=~time))
    AIC_vec[p,q] <- AIC(mod)}, error=function(e){AIC_vec[p,q] <- NA})
    
  }
}

# fit best model
mod.gls.22 <- update(mod.gls, correlation=corARMA(p=2,q=2,form=~time))

# check residuals of new correlation structures
res <- residuals(mod.gls.22,type="response")
names(res) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.2)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(res, main="", cex.lab=1.7, cex.axis=1.5)
pacf(res, main="", cex.lab=1.7, cex.axis=1.5)
hist(res, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)

# check model fit
par(mar=c(5,5,3,2))
plot(rural_nocovid$headcount,type="n",
     ylim=c(min(rural_ts_2[,1])-1000,max(rural_ts_2[,1])+100000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(rural_ts_2[,1])-100000,96,max(rural_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(rural_ts_2[,1]))
lines(0:86,fitted(mod.gls.22), col="#C11C06")

# create counterfactual predictions
datanew_count$month <- as.factor(datanew_count$month)
pred_count.gls <- predict(mod.gls.22 ,datanew_count,type="response")
attr(pred_count.gls,"label") <- NULL 

# create plot of headcounts
par(mar=c(5,5,3,2))
plot(rural_nocovid$headcount,type="n",
     ylim=c(min(rural_ts_2[,1])-1000,max(rural_ts_2[,1])+100000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(rural_ts_2[,1])-100000,96,max(rural_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(rural_ts_2[,1]))
# add counterfactual line
lines(36:86,pred_count.gls[-c(1:36)],col="#C11C06")

# get model summary
summary(mod.gls.22)
coef(mod.gls.22)
confint(mod.gls.22)


#-------------------------------------------------------------------------------
# Per-urban model
#-------------------------------------------------------------------------------

### ARIMA model ###

# create dataset using monthly headcounts
periurban_headcount <- data.frame(headcount = periurban_headcount_total, 
                               year = rep(c(2013,2014,2015,2016,2017,2018,2019,2020,2021,
                                            2022,2023), each = 12), 
                               month = rep(c(1,2,3,4,5,6,7,8,9,10,11,12), times = 11),
                               time = 1:132,
                               ccmdd = c(rep(0, times = 36), rep(1, times=96)),
                               covid = c(rep(0, times = 87), rep(1, times = 45)),
                               row.names = 1:132)

# remove covid years from dataset
periurban_nocovid <- periurban_headcount[-(88:132),]

# convert data to time series object
periurban_ts <- ts(periurban_nocovid[,1], frequency=12, start=c(2013,1))
periurban_ts

# look at lage plot
lag.plot(periurban_ts, lags = 12, do.lines = FALSE)

# look at ACF/PACF plots
acf2(periurban_ts, max.lag=24)

# look at ACF/PACF plots of differenced/seasonally differenced data
acf2(diff(diff(periurban_ts,12)), max.lag=24)

# create variable representing step change and view
step <- as.numeric(as.yearmon(time(periurban_ts))>='Jan 2016')
step

# create variable representing ramp (change in slope) and view
ramp <- append(rep(0,36), seq(1,51,1))
ramp  


# Use automated algorithm to identify p,q and d parameters
periurban_full <- auto.arima(periurban_ts, seasonal=TRUE, xreg=cbind(step,ramp), 
                         max.d=1, max.D=1, max.P = 3, stepwise=FALSE, trace=TRUE)

# check residuals
checkresiduals(periurban_full)
Box.test(periurban_full$residuals, lag = 24, type = "Ljung-Box")

res_periurban <- as.matrix(periurban_full$residuals)
names(res_periurban) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                      rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))

par(mar=c(5,5,3,2))
plot(res_periurban,pch=19,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date", type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.5)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(periurban_full$residuals, main="", cex.lab=1.7, cex.axis=1.5)
pacf(periurban_full$residuals, main="", cex.lab=1.7, cex.axis=1.5)
hist(periurban_full$residuals, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)


# Check model fit
plot(periurban_nocovid$headcount,type="n",
     ylim=c(min(periurban_ts)-1000,max(periurban_ts)+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96),cex.lab=1.7, cex.axis=1.2)
rect(36,min(periurban_ts)-200000,96,max(periurban_ts)+200000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(periurban_ts))
lines(0:86,fitted(periurban_full), col = "#C11C06")

# get parameters estimates and confidence intervals
summary(periurban_full)
round(confint(periurban_full),2)
coeftest(periurban_full)

# get counterfactual by modelling data excluding the post-intervention time 
# period
periurban_preint <- Arima(window(periurban_ts, end=c(2015,12)), order=c(2,1,0), 
                      seasonal=list(order=c(1,1,0), period=12))

# forecast 51 months post-intervention and create ts object with forecast
periurban_fc <- forecast(periurban_preint, h=51)
periurban_fc_ts <- ts(as.numeric(periurban_fc$mean), start=c(2016,1), frequency=12)

# combine forcasted data with observed data
periurban_ts_2 <- ts.union(periurban_ts, periurban_fc_ts)
periurban_ts_2

# plot headcount data
par(mar=c(5,5,3,2))
plot(periurban_nocovid$headcount,type="n",
     ylim=c(min(periurban_ts_2[,1])-1000,max(periurban_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(periurban_ts_2[,1])-100000,96,max(periurban_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020, cex.axis=1.4)
lines(0:86,na.omit(periurban_ts_2[,1]))
# add counterfactual line
lines(36:86,na.omit(periurban_ts_2[,2]), col="#C11C06")


### GLS regression model ###

# convert month to a factor variable
periurban_nocovid$month <- as.factor(periurban_nocovid$month)
# fit model with GLS
mod.gls <- gls(headcount ~  ccmdd * time + month,
               data = periurban_nocovid,
               correlation=corARMA(p=1, form=~time), method = "ML")
summary(mod.gls)

# check residuals
res <- residuals(mod.gls,type="response")
names(res) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n")
axis(1,at=0:7*12,labels = 2013:2020)
abline(h=0,lty=2,lwd=2)
acf(res, main="")
pacf(res, main="")

# check model fit
par(mar=c(5,5,3,2))
plot(periurban_nocovid$headcount,type="n",
     ylim=c(min(periurban_ts_2[,1])-1000,max(periurban_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7)
rect(36,min(periurban_ts_2[,1])-100000,96,max(periurban_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(periurban_ts_2[,1]))
lines(0:86,fitted(mod.gls), col="#C11C06")

# check AIC of different correlation structures
AIC_vec <- matrix(rep(0,49), ncol = 7, nrow = 7)
for(p in 1:7){
  for(q in 1:7){
    tryCatch({mod <-update(mod.gls, correlation=corARMA(p=p-1,q=q-1,form=~time))
    AIC_vec[p,q] <- AIC(mod)}, error=function(e){AIC_vec[p,q] <- NA})
    
  }
}

# fit best model
mod.gls.42 <- update(mod.gls, correlation=corARMA(p=4,q=2,form=~time))

# check residuals
res <- residuals(mod.gls.42,type="response")
names(res) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.2)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(res, main="", cex.lab=1.7, cex.axis=1.5)
pacf(res, main="", cex.lab=1.7, cex.axis=1.5)
hist(res, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)

# check model fit
par(mar=c(5,5,3,2))
plot(periurban_nocovid$headcount,type="n",
     ylim=c(min(periurban_ts_2[,1])-1000,max(periurban_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(periurban_ts_2[,1])-100000,96,max(periurban_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(periurban_ts_2[,1]))
lines(0:86,fitted(mod.gls.42), col="#C11C06")

# get counterfactual prediction
datanew_count$month <- as.factor(datanew_count$month)
pred_count.gls <- predict(mod.gls.42 ,datanew_count,type="response")
attr(pred_count.gls,"label") <- NULL

# plot headcount data
par(mar=c(5,5,3,2))
plot(periurban_nocovid$headcount,type="n",
     ylim=c(min(periurban_ts_2[,1])-1000,max(periurban_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(periurban_ts_2[,1])-100000,96,max(periurban_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(periurban_ts_2[,1]))
# add counterfactual line
lines(36:86,pred_count.gls[-c(1:36)],col="#C11C06")

# get model summary
summary(mod.gls.42)
coef(mod.gls.42)
confint(mod.gls.42)


#-------------------------------------------------------------------------------
# Urban model
#-------------------------------------------------------------------------------

### ARIMA model ###

# create dataset using monthly headcounts
urban_headcount <- data.frame(headcount = urban_headcount_total, 
                               year = rep(c(2013,2014,2015,2016,2017,2018,2019,2020,2021,
                                            2022,2023), each = 12), 
                               month = rep(c(1,2,3,4,5,6,7,8,9,10,11,12), times = 11),
                               time = 1:132,
                               ccmdd = c(rep(0, times = 36), rep(1, times=96)),
                               covid = c(rep(0, times = 87), rep(1, times = 45)),
                               row.names = 1:132)

# remove covid years from dataset
urban_nocovid <- urban_headcount[-(88:132),]

# convert data to time series object
urban_ts <- ts(urban_nocovid[,1], frequency=12, start=c(2013,1))
urban_ts

lag.plot(urban_ts, lags = 12, do.lines = FALSE)

# look at ACF/PACF plots
acf2(urban_ts, max.lag=24)

# look at ACF/PACF plots of differenced/seasonally differenced data
acf2(diff(diff(urban_ts,12)), max.lag=24)

# Create variable representing step change
step <- as.numeric(as.yearmon(time(urban_ts))>='Jan 2016')
step

# Create variable representing change in slope
ramp <- append(rep(0,36), seq(1,51,1))
ramp  


# Use automated algorithm to identify p,q and d parameters
urban_full <- auto.arima(urban_ts, seasonal=TRUE, xreg=cbind(step,ramp), 
                         max.d=1, max.D=1, stepwise=FALSE, trace=TRUE)

# check residuals
checkresiduals(urban_full)
Box.test(urban_full$residuals, lag = 24, type = "Ljung-Box")

res_urban <- as.matrix(urban_full$residuals)
names(res_urban) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                          rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))

par(mar=c(5,5,3,2))
plot(res_urban,pch=19,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date", type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.5)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(urban_full$residuals, main="", cex.lab=1.7, cex.axis=1.5)
pacf(urban_full$residuals, main="", cex.lab=1.7, cex.axis=1.5)
hist(urban_full$residuals, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)


# check model fit
plot(urban_nocovid$headcount,type="n",
     ylim=c(min(urban_ts)-1000,max(urban_ts)+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96),cex.lab=1.7, cex.axis=1.2)
rect(36,min(urban_ts)-200000,96,max(urban_ts)+200000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(urban_ts))
lines(0:86,fitted(urban_full), col = "#C11C06")

# get parameter estimates and confidence intervals
summary(urban_full)
round(confint(urban_full),2)
coeftest(urban_full)

# get counterfactual by modelling data excluding the post-intervention time 
# period
urban_preint <- Arima(window(urban_ts, end=c(2015,12)), order=c(1,0,3), 
                      seasonal=list(order=c(0,1,1), period=12))

# forecast 51 months post-intervention and create ts object with forecast
urban_fc <- forecast(urban_preint, h=51)
urban_fc_ts <- ts(as.numeric(urban_fc$mean), start=c(2016,1), frequency=12)

# combine forcasted data with observed data
urban_ts_2 <- ts.union(urban_ts, urban_fc_ts)
urban_ts_2

# plot headcount data
par(mar=c(5,5,3,2))
plot(urban_nocovid$headcount,type="n",
     ylim=c(min(urban_ts_2[,1])-1000,max(urban_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(urban_ts_2[,1])-100000,96,max(urban_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020, cex.axis=1.4)
lines(0:86,na.omit(urban_ts_2[,1]))
# add counterfactual line
lines(36:86,na.omit(urban_ts_2[,2]), col="#C11C06")


### GLS regression model ###

# convert month to a factor variable
urban_nocovid$month <- as.factor(urban_nocovid$month)
# fit model with GLS
mod.gls <- gls(headcount ~  ccmdd * time + month,
               data = urban_nocovid,
               correlation=corARMA(p=1, form=~time), method = "ML")
summary(mod.gls)

# check residuals
res <- residuals(mod.gls,type="response")
names(res) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n")
axis(1,at=0:7*12,labels = 2013:2020)
abline(h=0,lty=2,lwd=2)
acf(res, main="")
pacf(res, main="")

# check model fit
par(mar=c(5,5,3,2))
plot(urban_nocovid$headcount,type="n",
     ylim=c(min(urban_ts_2[,1])-1000,max(urban_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(urban_ts_2[,1])-100000,96,max(urban_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(urban_ts_2[,1]))
lines(0:86,fitted(mod.gls), col="#C11C06")

# get AIC for different correlation structures
AIC_vec <- matrix(rep(0,49), ncol = 7, nrow = 7)
for(p in 1:7){
  for(q in 1:7){
    tryCatch({mod <-update(mod.gls, correlation=corARMA(p=p-1,q=q-1,form=~time))
    AIC_vec[p,q] <- AIC(mod)}, error=function(e){AIC_vec[p,q] <- NA})
    
  }
}

# fit best model
mod.gls.64 <- update(mod.gls, correlation=corARMA(p=6,q=4,form=~time))

# check residuals
res <- residuals(mod.gls.64,type="response")
names(res) <- c(rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
                 rep(2017,12),rep(2018,12),rep(2019,12),rep(2020,3))
plot(res,pch=19,cex=0.7,col=grey(0.6),main="",
     ylab="Deviance residuals",xlab="Date",type = "l", xaxt="n", cex.lab=1.7, cex.axis=1.2)
axis(1,at=0:7*12,labels = 2013:2020, cex.axis=1.5)
abline(h=0,lty=2,lwd=2)
acf(res, main="", cex.lab=1.7, cex.axis=1.5)
pacf(res, main="", cex.lab=1.7, cex.axis=1.5)
hist(res, breaks = 15, xlab = "Residuals", main = "", cex.lab=1.7, cex.axis=1.5)

# check model fit
par(mar=c(5,5,3,2))
plot(urban_nocovid$headcount,type="n",
     ylim=c(min(urban_ts_2[,1])-1000,max(urban_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(urban_ts_2[,1])-100000,96,max(urban_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(urban_ts_2[,1]))
lines(0:86,fitted(mod.gls.64), col="#C11C06")

# get counterfactual prediction
datanew_count$month <- as.factor(datanew_count$month)
pred_count.gls <- predict(mod.gls.64 ,datanew_count,type="response")
attr(pred_count.gls,"label") <- NULL 

# plot headcount data
par(mar=c(5,5,3,2))
plot(urban_nocovid$headcount,type="n",
     ylim=c(min(urban_ts_2[,1])-1000,max(urban_ts_2[,1])+1000),xlab="Year", 
     ylab="Headcount", bty="l",xaxt="n", xlim = c(0,96), cex.lab=1.7, cex.axis=1.2)
rect(36,min(urban_ts_2[,1])-100000,96,max(urban_ts_2[,1])+100000,col=grey(0.9),
     border=F)
axis(1,at=0:8*12,labels=F)
axis(1,at=0:7*12+6,tick=F,labels=2013:2020,cex.axis=1.4)
lines(0:86,na.omit(urban_ts_2[,1]))
# add counterfactual line
lines(36:86,pred_count.gls[-c(1:36)],col="#C11C06")

# get model summary
summary(mod.gls.64)
coef(mod.gls.64)
confint(mod.gls.64)
