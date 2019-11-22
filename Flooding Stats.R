##############     Calculate inundation statistics      #####################
##############                                            ##########################
## this code uses a preconstructed list of dataframes (as output from the "Mangrove Points.R script"), each containing a column for every location and a row for each water level observation
## the code cycles through these dataframes and calculates the flooding statistics for each location
## these statistics are saved to the computer after the completion of each dataframe

library(lubridate)  # this package is necessary for the timing statistics
library(plyr) # this package is necessary for aggregations

Mang <- readOGR("D:/Dropbox/Dissertation/GIS files/MangroveDistribution/PRMangrovesCombined/PRMangroves_HOBOS.shp")

mflood_data_t <- as.data.frame(rep(seq(1:70),16))  ## create the monthly statistics data set, with the number of months in the time period (70) by the number of sites (16)
colnames(mflood_data_t)<-"month" ## this first column is the consecutive month number
mflood_data_t$loc <- rep(unique(Mang$HOBO),each=70) ## the second column is the site name
mflood_data <- mflood_data_t  # a temporary data set that will be used within the loop
listofdfs <- readRDS("D:/Documents/PhD/Water Level dfs/dfs_100.rds")  ## read in the list of dataframes created in the previous step
listofpoints <- readRDS("D:/Documents/PhD/Water Level dfs/points_100.rds") ## read in the list of coordinates corresponding to the columns in the list of data frames
d=i=c=1  ## these counters will be used in the loop
total <- sum(sapply(listofdfs, length))  # the total number of repitions, which is used by the progress bar
pb <- txtProgressBar(min = 1, max = total, style = 3) # create progress bar
repeat{  ### go through each set of waterlevel values at each point within each site
  df <- listofdfs[[1]]  ### start with the first site. This df will be removed from the list at the end of the loop and we will come back here to continue with the next
  points <- listofpoints[[1]]  ###  get the list of coordinates corresponding to each site in the df
  c = 1  ## set the point counter to one
  for (C in colnames(df)[1:(ncol(df)-7)]){  ## Each column in the df is a random location and each row is the water level at the location for every hour from 2012 to 2017. For each column in the df, we will calculate the flooding statistics. Except for the last 7 columns, which are identifiers
    flood_dur_temp <- rle(df[,C]>=0) ### all flood events
    dry_dur_temp <- rle(df[,C]<0)  ### all dry events
    flood_freq_temp <- rep(flood_dur_temp$values,flood_dur_temp$lengths)  ## all dry and flood conditions
    if (length(flood_dur_temp$lengths)==1){  ## if the length of the flood_dur_temp data frame is one, this means the water never changed from flood to dry throughout the time period. So we can easily calculate the statistics
      if (flood_dur_temp$values==F){ ## If the value for the 1 observation is FALSE (F), this means that location was always dry
        flood_dur_temp <- dflood_depth_temp <- dflood_freq_temp <- mflood_depth_temp <- mflood_freq_temp <-as.data.frame(0) ## so the flooding depth and frequency are 0
        dry_dur_temp <- as.data.frame(nrow(df)) ## the dry duration is the total length of the time period
        dflood_freq_temp$freq <- dflood_depth_temp$depth <- mflood_depth_temp$depth <- mflood_freq_temp$floods <- 0  ## set the frequency and depths to 0
        colnames(mflood_depth_temp)[1] <- colnames(mflood_freq_temp)[1] <- "V2"  ## to be consistent with later steps and to ensure the cbind works correctly
      }else{  ## otherwise the site is permanently flooded
        flood_dur_temp <- as.data.frame(nrow(df))  ## the flood duration is the length of the analysis, which is the number of rows
        dry_dur_temp <- dflood_depth_temp <- as.data.frame(0) ##  set the dry and flood durations dataframes
        mflood_depth_temp <- mflood_freq_temp <- as.data.frame(c(NA))  ## the monthly flood depth and freq are set to NA
        dflood_freq_temp <- as.data.frame(flood_dur_temp[,1]) ## set the daily flood duration to the length of the data frame
        dflood_freq_temp$freq <- mflood_depth_temp$depth <- mflood_freq_temp$floods <- Inf ## set the monthly flood duration, depth, and number of floods to infinity
        dflood_depth_temp$depth <- mean(df[,C])  ## the daily flood depth is the average depth at the site
        colnames(mflood_depth_temp)[1] <- colnames(mflood_freq_temp)[1] <- "V2" ## to be consistent with later steps and to ensure the cbind works correctly
      }
    }else{  ##otherwise, the site is periodically flooded, so calculate the monthly and daily statistics
      flood_dur_temp <- as.data.frame(flood_dur_temp$lengths[flood_dur_temp$values==T]) ## only flood events
      dry_dur_temp <- as.data.frame(dry_dur_temp$lengths[dry_dur_temp$values==T])  ## only dry events
      mflood_freq_temp <- as.data.frame(cbind(flood_freq_temp,df[,"month"])) ## set up the monthly stats table, starting with the column representing the consecutive month number of each observation
      dflood_freq_temp <- as.data.frame(cbind(flood_freq_temp,df[,"day"])) ## set up the daily stats table, starting with the column representing the consecutive day number of each observation
      mflood_depth_temp <- cbind(mflood_freq_temp,df[,C]) # join the water level to the month column
      dflood_depth_temp <- cbind(dflood_freq_temp,df[,C]) # join the water level to the day column
      mflood_depth_temp <- mflood_depth_temp[mflood_depth_temp$flood_freq_temp==1,] ## extract the depths only when flooded in the monthly dataset
      dflood_depth_temp <- dflood_depth_temp[dflood_depth_temp$flood_freq_temp==1,]  ## extract the depths only when flooded in the daily dataset
      colnames(mflood_depth_temp)[3] <- "depth"  ## the new column is the depth column
      colnames(dflood_depth_temp)[3] <- "depth"  ## the new column is the depth column
      mflood_freq_temp <- ddply(mflood_freq_temp, .(V2), summarise, 
                                floods = sum(rle(flood_freq_temp)$values==1))  ### aggregate the values into monthly flood frequency
      dflood_freq_temp <- ddply(dflood_freq_temp, .(V2), summarise, 
                                floods = sum(rle(flood_freq_temp)$values==1))  ### aggregate the values into daily flood frequency
      mflood_depth_temp <- ddply(mflood_depth_temp, .(V2), summarise, 
                                 depth = mean(depth,na.rm=T))  ### aggregate the values into monthly flood depth
      dflood_depth_temp <- ddply(dflood_depth_temp, .(V2), summarise, 
                                 depth = mean(depth,na.rm=T))  ### aggregate the values into daily flood depth
    }
    # with the basic flooding stats calculated, create a data table with them all 
    dflood_depth_temp <- as.data.frame(cbind(nrow(dflood_depth_temp),mean(dflood_depth_temp$depth,na.rm=T),max(dflood_depth_temp$depth,na.rm=T))) ## the daily data frame is the number of days in which flooding ocurred nrow(dflood_depth_temp), the average flooding depth mean(dflood_depth_temp$depth,na.rm=T), and the maximum flood depth max(dflood_depth_temp$depth,na.rm=T)
    colnames(dflood_depth_temp) <- c("nFloodDays","meanFloodDepth","maxFloodDepth") ## set the column names for the daily stats
    mflood_depth_temp$loc <- unique(df$loc) ## assign the site name to the depth data
    mflood_freq_temp$loc <- unique(df$loc) ## assign the site name to the frequency data
    mflood_dur_temp <- aggregate(df[,C]>=0~month,data=df,FUN=sum) ### aggregate the data into the total flood duration per month
    mdry_dur_temp <- aggregate(df[,C]<0~month,data=df,FUN=sum) ### aggregate the data into the total dry duration per month
    mFlood_prop_temp <- 100*mflood_dur_temp[,2]/(mflood_dur_temp[,2]+mdry_dur_temp[,2]) ## the proportion of time flooded as a percentage
    colnames(mflood_dur_temp) <- c("month","flood_length_hr")  # set the names of the columns
    mFlood_prop_temp <- as.data.frame(cbind(mflood_dur_temp$month,mFlood_prop_temp)) ## attach the months to the flood proportion data
    mflood_dur_temp$loc <- unique(df$loc)  # set the site name to the data set
    colnames(mdry_dur_temp) <- c("month","dry_length_hr") # set the names of the columns
    mdry_dur_temp$loc <- unique(df$loc) # set the site name to the data set
    colnames(mFlood_prop_temp) <- c("month","flood_length_hr") # set the names of the columns
    mFlood_prop_temp$loc <- unique(df$loc)# set the site name to the data set
    mflood_range_temp <- aggregate(get(C)~month,data=df,FUN=min,na.rm=T) ### get the minimum depth for each month
    mflood_range_temp$max <- aggregate(get(C)~month,data=df,FUN=max,na.rm=T)[,2] ### get the maximum depth for each month
    colnames(mflood_range_temp)[1:2] <- c("month","min") # set the names of the columns
    mflood_range_temp$range <- mflood_range_temp$max-mflood_range_temp$min ## calcualate the full depth range for each month
    mflood_range_temp$loc <- unique(df$loc) # set the site name to the data set
    dflood_freq_temp <- as.data.frame(cbind(nrow(dflood_freq_temp[dflood_freq_temp$floods>0,]),mean(dflood_freq_temp$floods,na.rm=T),max(dflood_freq_temp$floods,na.rm=T)))  ## the daily data frame is the number of times in which flooding ocurred ina a day nrow(dflood_freq_temp[dflood_freq_temp$floods>0,]), the average daily flooding frequency mean(dflood_freq_temp$floods,na.rm=T), and the maximum daily flooding frequency max(dflood_freq_temp$floods,na.rm=T)
    colnames(dflood_freq_temp) <- c("nFloodDays","meandFloodfreq","maxdFloodfreq") # set the names of the columns
    dflood_freq_temp$loc <- unique(df$loc) # set the site name to the data set
    flood_dur_temp <- as.data.frame(mean(flood_dur_temp[,1],na.rm=T)) ## the daily flood duration
    dry_dur_temp <- as.data.frame(mean(dry_dur_temp[,1],na.rm=T)) ## the daily dry duration
    Flood_prop_temp <- 100*flood_dur_temp/(flood_dur_temp+dry_dur_temp) ## the daily flood proportion as a percentage
    colnames(flood_dur_temp) <- "FloodLength_hours" # set the names of the columns
    colnames(dry_dur_temp) <- "DryLength_hours" # set the names of the columns
    colnames(Flood_prop_temp) <- "PropFlooded" # set the names of the columns
    flood_dur_temp$loc <- unique(df$loc) # set the site name to the data set
    dry_dur_temp$loc <- unique(df$loc) # set the site name to the data set
    Flood_prop_temp$loc <- unique(df$loc) # set the site name to the data set
    ### now join all statistics to one data set, this is only for this particular location
    flooddata_temp <- as.data.frame(c(flood_dur_temp[,2],flood_dur_temp[,1],dry_dur_temp[,1],Flood_prop_temp[,1],dflood_freq_temp[,c(1:3)],dflood_depth_temp[,c(2:3)],
                                      mean(df[,C],na.rm=T),points$x[c],points$y[c]))
    colnames(flooddata_temp) <- c("loc","FloodLength_hours","DryLength_hours","FloodProp","nFloodDays","meandFloodfreq",
                                  "maxdFloodfreq","meanFloodDepth","maxFloodDepth","Avg_depth","x","y")
    flooddata_temp[flooddata_temp=="NaN"] <- NA
    ###########   store values for all  locations   ##########
    ## take the datasets for each location and add them to a data set that will contain all the statistics for each location
    if (d==1&C==colnames(df)[1]){ ## if this is the first location, we have to first create the data set, which will be equal to the above created "temp" data sets
      flooddata <- flooddata_temp
      mflood_depth <- merge(mflood_data_t,mflood_depth_temp,by.x=c("month","loc"),by.y=c("V2","loc"),all.x=T)
      mflood_freq <- merge(mflood_data_t,mflood_freq_temp,by.x=c("month","loc"),by.y=c("V2","loc"),all.x=T)
      mflood_dur <- merge(mflood_data_t,mflood_dur_temp,by=c("month","loc"),all.x=T)
      mdry_dur <- merge(mflood_data_t,mdry_dur_temp,by=c("month","loc"),all.x=T)
      mFlood_prop <- merge(mflood_data_t,mFlood_prop_temp,by=c("month","loc"),all.x=T)
      mflood_range <- merge(mflood_data_t,mflood_range_temp[,c("month","loc","range")],by=c("month","loc"),all.x=T)
    }else{  ## othwerwise, add the individual location statistics to the data set containing all locations
      flooddata <- rbind(flooddata,flooddata_temp)
      mflood_depth <-  merge(mflood_depth,mflood_depth_temp,by.x=c("month","loc"),by.y=c("V2","loc"),all.x=T)
      mflood_freq <-  merge(mflood_freq,mflood_freq_temp,by.x=c("month","loc"),by.y=c("V2","loc"),all.x=T)
      mflood_dur <-  merge(mflood_dur,mflood_dur_temp,by=c("month","loc"),all.x=T)
      mdry_dur <-  merge(mdry_dur,mdry_dur_temp,by=c("month","loc"),all.x=T)
      mFlood_prop <- merge(mFlood_prop,mFlood_prop_temp,by=c("month","loc"),all.x=T)
      mflood_range <- merge(mflood_range,mflood_range_temp[,c("month","loc","range")],by=c("month","loc"),all.x=T)
    }
    i=i+1 ## a running counter for the progress bar
    c = c+1 ## move to the next point
    Sys.sleep(0.05) ## to reset the progress bar, sleep the progra for 0.05 seconds
    setTxtProgressBar(pb,i) # update progress bar
  }
  ## when all columns and all locations have been analyzed, join the data to a master dataset of all sites
  flooddata_means <- aggregate(cbind(FloodLength_hours,DryLength_hours,FloodProp,nFloodDays,meandFloodfreq,meanFloodDepth,Avg_depth) ~ loc, 
                               data = flooddata, mean, na.rm = TRUE)  ## find the mean values for that site 
  mflood_depth <- cbind(mflood_depth[,(1:2)],rowMeans(mflood_depth[,3:ncol(mflood_depth)],na.rm=T))   # bind the flood depth data to the previously created dataset with month and site columns
  colnames(mflood_depth)[3] <- "flood_depth"  # name the new column 
  setDT(mflood_data)[setDT(mflood_depth[mflood_depth$loc==df$loc,]), flood_depth := i.flood_depth, on = .(month,loc)]  ## bind this new column to the master dataset 
  mflood_freq <- cbind(mflood_freq[,(1:2)],rowMeans(mflood_freq[,3:ncol(mflood_freq)],na.rm=T))  # bind the flood frequency data to the previously created dataset with month and site columns
  colnames(mflood_freq)[3] <- "flood_freq" # name the new column 
  setDT(mflood_data)[setDT(mflood_freq[mflood_freq$loc==df$loc,]), flood_freq := i.flood_freq, on = .(month,loc)] ## bind this new column to the master dataset 
  mflood_dur <- cbind(mflood_dur[,(1:2)],rowMeans(mflood_dur[,3:ncol(mflood_dur)],na.rm=T)) # bind the flood duration data to the previously created dataset with month and site columns
  colnames(mflood_dur)[3] <- "flood_dur"  # name the new column  
  setDT(mflood_data)[setDT(mflood_dur[mflood_dur$loc==df$loc,]), flood_dur := i.flood_dur, on = .(month,loc)] ## bind this new column to the master dataset 
  mdry_dur <- cbind(mdry_dur[,(1:2)],rowMeans(mdry_dur[,3:ncol(mdry_dur)],na.rm=T)) # bind the dry duration data to the previously created dataset with month and site columns
  colnames(mdry_dur)[3] <- "dry_dur" # name the new column  
  setDT(mflood_data)[setDT(mdry_dur[mdry_dur$loc==df$loc,]), dry_dur := i.dry_dur, on = .(month,loc)] ## bind this new column to the master dataset 
  mFlood_prop <- cbind(mFlood_prop[,(1:2)],rowMeans(mFlood_prop[,3:ncol(mFlood_prop)],na.rm=T)) # bind the flood proportion data to the previously created dataset with month and site columns
  colnames(mFlood_prop)[3] <- "flood_prop"  # name the new column  
  setDT(mflood_data)[setDT(mFlood_prop[mFlood_prop$loc==df$loc,]), flood_prop := i.flood_prop, on = .(month,loc)] ## bind this new column to the master dataset 
  mflood_range <- cbind(mflood_range[,(1:2)],rowMeans(mflood_range[,3:ncol(mflood_range)],na.rm=T)) # bind the depth range data to the previously created dataset with month and site columns
  colnames(mflood_range)[3] <- "range"  # name the new column  
  setDT(mflood_data)[setDT(mflood_range[mflood_range$loc==df$loc,]), range := i.range, on = .(month,loc)] ## bind this new column to the master dataset 
  write.csv(mflood_data,"D:/Documents/PhD/Water Level dfs/mflood_data_partial_100.csv",row.names = F) ## save the monthly data
  write.csv(flooddata_means,"D:/Documents/PhD/Water Level dfs/flooddata_means_partial_100.csv",row.names = F) ## save the means data
  write.csv(flooddata,"D:/Documents/PhD/Water Level dfs/flooddata_partial_100.csv",row.names = F) ## save the daily data
  d=d+1  ## to make sure we dont create another master data set when we finish analyzing the next site.
  listofdfs[[1]] <- NULL ## erase the current data frame from the list of data frames, this will make the next data fram available for analysis
  listofpoints[[1]] <- NULL ## erase the current points from the list of points, this will make the points available for analysis
  remove(mflood_freq_temp,mflood_depth_temp,dflood_depth_temp,flood_freq_temp,dflood_freq_temp,dry_dur_temp,flood_dur_temp,
         Avg_depth,flooddata_temp,flooddata_temp2,mflood_range_temp,mflood_dur_temp,mdry_dur_temp,mFlood_prop_temp,Flood_prop_temp) ## remove all of the temporary datasets
  gc() ## clean up the environment
  if (length(listofdfs) ==0){ ## if we've finished with all data frames, remove the datasets and break from the repeat loop
    remove(mflood_depth,mflood_dur,mflood_freq,mflood_range,mdry_dur,mFlood_prop)
    break
  }
}

###   Save the resulting datasets. This can also be done in the loop to make sure progress is saved in case the computer shuts down or other errors. 
flooddata_means <- aggregate(cbind(FloodLength_hours,DryLength_hours,nFloodDays,meandFloodfreq,maxdFloodfreq,meanFloodDepth,maxFloodDepth,Avg_depth) ~ loc, 
                             data = flooddata, mean, na.rm = TRUE)
write.csv(flooddata_means,"flooddata_means_partial_100.csv",row.names = F)
write.csv(flooddata,"flooddata_partial_100.csv",row.names = F)
write.csv(mflood_data,"mflood_data_partial_100.csv",row.names = F)
