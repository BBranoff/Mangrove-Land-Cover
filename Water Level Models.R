## This subroutine calculates an estuarine water level model for a given set of water level observations.
## The model includes tidal harmonics and precipitation inputs.
## Required Inputs are:
##    A set of water level observations (WTRLVLS) taken with at least an hour frequency and over at least a few months. This data set has a "time" column and 
##        an "h" column, which is the water level in meters
##    A set of precipitation observations (Weather) matching the time-period of the water levels. In this case the observations are hourly, and this will be the
##        resulting frequency of the water level model. The precipitation dataset must contain the full desired length of the study, in this case five years.
##        This is because rainfall must be known to calculate its contribution to waterlevels. Tides, however can be calculated with only a few months of observations.
##    


library(data.table)
library(zoo)
library(TideHarmonics)

Mods <- read.csv("Mods.csv")  ## This is the complete list of observations and models. We will subset it to contain only one set of waterlevel observations from PONMID
WTRLVLS <- Mods[,c("time","Ponhmid")]   ## only the time column and the observstions from PONMID
names(WTRLVLS) <- c("time","h")

Weather$time <- as.POSIXct(Weather$time,format="%m/%d/%Y %H:%M",tz="america/caracas") ## Make sure the date column of the precipitation data is in POSIXct format and in the proper timezone
WTRLVLS$time <- as.POSIXct(WTRLVLS$time,format="%m/%d/%Y %H:%M",tz="america/caracas")
setDT(WTRLVLS)[setDT(Weather), Rain := i.Rain, on = "time", roll="nearest"]  ## bind the precipitation data to the nearest time of the water levels data
WTRLVLS$Rain[is.na(WTRLVLS$Rain)] <- 0   ## any gaps or na values are set to 0 precipitation

#####  Set aside 10% for validation   #################
WTRLVLS_val <- WTRLVLS  ## save the complete observations for later
###  set the last 10% of observations to NA. These wont be included in the calibration, so we can test the model later on these for validation
WTRLVLS$h[(max(which(!is.na(WTRLVLS$h)))-0.1*(max(which(!is.na(WTRLVLS$h)))-min(which(!is.na(WTRLVLS$h))))):max(which(!is.na(WTRLVLS$h)))] <- NA
write.csv(WTRLVLS,"WTRLVLS.csv",row.names=F)  ## save the observations for later use

#############             Tidal Models    ###############
###   Before calculating tidal harmonics, we need to exclude observations where rain might be an influence
WTRLVLS$RainMS <- rollapply(WTRLVLS$Rain,100,sum,fill="NA",align="right")  ## calculate the 100 hour (~4 day) moving sum of rainfall
WTRLVLS2 <- subset(WTRLVLS,!WTRLVLS$RainMS>0)  ## create a separate data frame that does not include any observations where the 4 day rainfall was greater than 0
WTRLVLS2$h <- WTRLVLS2$h - mean(WTRLVLS2$h,na.rm=T)  ## normalize the data to be centered on 0
WTRLVLSFit <- ftide(WTRLVLS2$h, WTRLVLS2$time, hc60) ## from the tidal harmonics package, calculate the tidal constituents from the water level data
saveRDS(WTRLVLSFit,"WTRLVLSFit.rds")  ## save the fitted tidal constituents for later use

## create a new data frame where we will create tidal predictions for each timestep
## The weather data set must contain the full length of teh study period
WTRLVLPreds <- as.data.frame(seq.POSIXt(as.POSIXct("2012-06-01 01:00:00"),max(WTRLVLS$time),by="hour" ))
attr(WTRLVLPreds[,1],"tzone") = "America/Caracas"  ## set the timeznone of the times
colnames(WTRLVLPreds) <- "time"
WTRLVLPreds$hPred <- predict(WTRLVLSFit,from=min(WTRLVLPreds$time), to=max(WTRLVLPreds$time),by=1) ## calculate the predicted water levels based on tides only
setDT(WTRLVLPreds)[setDT(WTRLVLS), h := i.h, on = "time", roll="nearest"]  ## Join the actual water level observations to the data frame
setDT(WTRLVLPreds)[setDT(Weather), Rain := i.Rain, on = "time", roll="nearest"]  ## Join the rainfall observations to the data frame
WTRLVLPreds$Rain[is.na(WTRLVLPreds$Rain)] <- 0



######################################################
##########   Rain model  #############################
######################################################
###   first the longterm precipitation input
###  Not all sites will require this longterm input. Usually only sites with poor tidal connectivity will respond to long-term rainfall
###  If there is no long-term input the coefficients calculated below will max out and you can ignore this long-term input.
###  The optimization uses a bobyqa algorithm to find the best rain coeffifients that minimize the error between the observations and the model
###  It will likely require some trial and erro to find the best model. You will need to play around with the moving average window as well
###  as with the coefficients. 


library(minqa)
library(TTR)   ## for zlema
min.RSS <- function(par) {
  x =  ZLEMA(as.data.frame(WTRLVLPreds)[!is.na(WTRLVLPreds$h),"h"],24*96)  ##  The long term rainfall should match a moving average of long-term water levels, in this case 20 days. This timeframe might change for each site.
  y =  as.data.frame(WTRLVLPreds$Rain[!is.na(WTRLVLPreds$h)])
  y =  rollapply(y,par[1],sum,fill=NA,align="right")/par[2]   ## the par[1] and par[2] are the variables we are trying to optimize.
  sum((x - y)^2,na.rm=TRUE)  ## optimization depends on minimizing the residual sum of squares between the moving average of water levels (x) and the predicted long-term rainfall model
}
lower=c(1,.00001)  ## set the lower limits for the two parameters par[1] and par[2]
upper=c(8000,8000)  ## set the upper limits for the two parameters par[1] and par[2]
opt <- bobyqa(c(500,300), min.RSS, lower, upper)  ## run the bobyqa optimization to find the ideal parameters that match the ong-term water level average.
pars <- opt$par ## save the parameters as a variable

WTRLVLPreds$Rain_long <- rollapply(WTRLVLPreds$Rain,pars[1],sum,fill=NA,align="right")/pars[2]  ## create the long-term precipitation input in the water level data frame

#### look at the results and mosify the bobyqa parameters as needed
plot(WTRLVLS$time,WTRLVLS$h, pch=19,cex=0.5)     ###  first plot the actual observations
points(WTRLVLS$time[!is.na(WTRLVLS$h)],ZLEMA(WTRLVLS$h[!is.na(WTRLVLS$h)],24*96),col="red",pch=19,cex=0.5)   ## then the moving average. This is what the long-term model is trying to replicate
points(WTRLVLPreds$time,WTRLVLPreds$Rain_long,col="blue",cex=0.5,pch=19)   ###  now the model

#########################################
######   Now the short term model
#########################################
## this should be detectable in most systems, except open ocean or other large water bodies that wont respond to most rain events
##  Again, you will need to do some trial and error to find the best parameters. 
## Much of the success of this step depends on how accurate the long-term rain model is.

min.RSS <- function(par) {
  x = as.data.frame(WTRLVLPreds)[!is.na(WTRLVLPreds$Rain),"h"]   ### the actual water level observations
  y =  as.data.frame(WTRLVLPreds)[!is.na(WTRLVLPreds$Rain),"hPred"]+
    (EMA(WTRLVLPreds$Rain[!is.na(WTRLVLPreds$Rain)],n=par[1])/par[2])  ## The exponential moving average of rainfall added to the tidal mode. Again, 
##                                                                        the parameters par[1] and par[2] will be optimized to match the observations
  sum((x - y)^2,na.rm=TRUE)
}
lower=c(1,0)  ## set the lower limits for the two parameters par[1] and par[2]
upper=c(10000,10000)  ## set the upper limits for the two parameters par[1] and par[2]
opt <- bobyqa(c(20,12), min.RSS, lower, upper)  ## run the bobyqa optimization to find the ideal parameters that match the short-term water level fluctuations.
pars <- opt$par ## save the parameters as a variable
WTRLVLPreds$Rain_short <-(EMA(WTRLVLPreds$Rain[!is.na(WTRLVLPreds$Rain)],n=pars[1])/pars[2])


#### look at the results and mosify the bobyqa parameters as needed
plot(WTRLVLS$time,WTRLVLS$h, pch=19,cex=0.5)     ###  first plot the actual observations
points(WTRLVLPreds$time,WTRLVLPreds$Rain_long+WTRLVLPreds$Rain_short,col="blue",cex=0.5,pch=19)   ###  now the two combined rain models
points(WTRLVLPreds$time,WTRLVLPreds$Rain_long+WTRLVLPreds$Rain_short+WTRLVLPreds$hPred,col="red",cex=0.5,pch=19)   ###  now the two combined rain models and the tidla model


###############################################
####    Test with 10% for validation    ######
############################################
WTRLVLS_val <- WTRLVLS_val[WTRLVLS_val$time %in% WTRLVLPreds$time,]   ### get the observations included in the model
WTRLVLS_val$h[1:(max(which(!is.na(WTRLVLS_val$h)))-0.1*(max(which(!is.na(WTRLVLS_val$h)))-min(which(!is.na(WTRLVLS_val$h)))))] <- NA  ### remove the calibration observations
setDT(WTRLVLPreds)[WTRLVLS_val, h_val := i.h, on = "time", roll="nearest"]  ##  add in the extra observations 

plot(WTRLVLPreds$h_val[!is.na(WTRLVLPreds$h_val)],rowSums(WTRLVLPreds[!is.na(WTRLVLPreds$h_val),c("hPred","Rain_long","Rain_short")]))  ## this is a plot of the observations vs the model, including tides and both rainfalls. A line of slope 1 would be perfect
abline(lm(rowSums(WTRLVLPreds[!is.na(WTRLVLPreds$h_val),c("hPred","Rain_long","Rain_short")])~WTRLVLPreds$h_val[!is.na(WTRLVLPreds$h_val)]),col="red")  ## the line describing the reationship between the observations and the predictions, again, this should be a line of slope 1

######   if satisfied, save the results for later processing   ########
write.csv(WTRLVLPreds,"WTRLVLPreds.csv",row.names=F)
