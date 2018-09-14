##################################################################
########    Create Water Level Data Frames at Random Mangrove Locations 
####################################################################
### This subroutine creates random points within the mangrove habitat. It then calculates the water levels at each point based on the water level models
### and the elevation of the random point and the piezometer. It saves the points as a list of dataframes, one dataframe for each location. These are used 
### later to calculate urbanization. It also saves the water levels for each points, again one dataframe for each location. Depending upon the number of 
##  locations and the number of points for each location, this list can get very big. For 16 locations and 1000 points at each location, the list of water level
###  data frames is around 3 GB. 
### Required Inputs for the process are:
###         A data frame of modeled water levels for each location for the study time period (Mods). 
##          A shapefile of the mangrove habitat (Mang), which has already been divided into separate locations. The locations are saved within the shapefile as
###               an attribute with the name "loc". The dataframes will be saved with a column that matches each location within the mangrove shapefile.
###         A raster of elevations (ElevPon and Elev), or a digital elevation model. The  vertical datum for the DEM is not important, because the water levels are calculated
##                relative to the elevation of the piezometer.In this case, there is one DEM for the north of PR and one for the south (Ponce)
###         A shapefile of piezometer locations. This shapefile must have an attribute called "HOBOS", which must match the name of the mangrove locations, "loc".
library(raster)   ## spatial processing package
library(rgdal)   ## spatial processing package
library(lubridate)  ## this package is necessary for calculating the day, week, month, and year numbers of the timeseries

listofdfs <- list()
listofdfs_100 <-list()
listofpoints <- list()
listofpoints_100 <-list()

Mods <- read.csv("Mods.csv")      ## this is the set of waterlevel observations and models at each piezometer. It has the tidal and rainfall models as well as both of them combined into a full prediction. It also has rainfall observations for each stusy area
                                  ## The names of the columns of Mods is important, they must contain the three letter location identifiers that will be matched to the appropriate elevations and piezometer
                                  ## For example: "PonRain" is the rainfall observations at Ponce. "Ponhmin" is the water levels observations at PONMIN site, "PonhminPred" is the tidal model at PONMIN, "PONRainMSmin" is the shorterm rainfall model, "Ponhminmod" is the combination of tides and rainfall and is the complete model.
                                  ## Each waterlevel model from Mods will be matched with the appropriate piezometer and mangrove elevations based on the names.
Mang <- readOGR("D:/Dropbox/Dissertation/GIS files/MangroveDistribution/PRMangrovesCombined/PRMangroves_HOBOS.shp")  ## the pre-delimited mangrove areas. Each is names as "HOBO" which corresponds to the piezometer and the water level models for that piezometer
Elev <- raster("D:/Documents/GIS DataBase/San_Juan_DEM_1602/San_Juan_DEM_1602/DEM.tif")  #DEMs are referenced to the vertical tidal datum of Mean High Water (MHW)
ElevPon <- raster("D:/Documents/GIS DataBase/PonceDEM/PonceMosaicNOAASLR_DEM.tif") 
ElevPon <- crop(ElevPon, extent(Mang))
ElevPon <- aggregate(ElevPon,fact = res(Elev)/res(ElevPon),fun="mean",na.rm=T)
Mods$time <- as.POSIXct(Mods$time)  ##  Set the time column of the Mods dataframe to a POSIXct data type
for (i in 1:length(Mang$HOBO)){  ### For each location within the mangrove (each piezometer), we will create a dataframe of waterlevels at random locations
  if (Mang$HOBO[i] %in% c("PONMAX","PONMIN","PONMID")){  ## If the mangroves are in Ponce and include "PON" in their name, we use the Ponce DEM
    rec <- crop(ElevPon, extent(Mang[i,]))  ## rec will be the dem used for each location
  }else{  ## If the location is not in Ponce, we will use the San Juan DEM, which includes Levittown
    rec <- crop(Elev, extent(Mang[i,]))
  }
  rec <- mask(rec, Mang[i,])  ## remove elevation values outside of mangrove habitat  
  points <- sampleRandom(rec,1000, na.rm=TRUE, xy=TRUE)   #sample elevations from 1000 random points inside the mangroves and save the x and y coordinates
  points <- points[!duplicated(points[,c(1,2)]),]  ## make sure any duplicated points are not include
  points <- points[points[,3] < quantile(points[,3], 0.90)&points[,3]>quantile(points[,3], 0.10),]  ## exclude points if they are not within the 10-90th percentile of elevation. This eleminates points that are likely placed outside of the mangrove zone due to misalignment between the two datasets
  df<- data.frame(matrix(NA, ncol = nrow(points), nrow = nrow(Mods)))  ## create an empty data frame for this location. There will be a column for each of the random points and a row for each timestep in the water level model
  elev <- extract(rec,HOBOloc[HOBOloc$HOBOS==Mang$HOBO[i],])  ## extract the elevation of the piezometer
  for (j in 1:nrow(points)){  ## for each of the random points, calculate the water level at each time step
    df[,j] <- Mods[,gsub("Pred","mod",Mang$Preds[i])]+(elev-points[j,3])  ## the water level is the water level of the piezometer plus the elevation of the piezometer minus the elevation fo the point
  }
  df$time <- Mods$time   ## save a time column to the datafram, which matches the time column from Mods
  df$hour <- as.numeric(interval(df$time[1], df$time) %/% hours(1)) + 1 ## save an hour column to the dataframe
  df$day <- as.numeric(floor(julian(df$time, origin = df$time[1]))+ 1) ## save an day column to the dataframe
  df$week <- as.numeric(interval(df$time[1], df$time) %/% weeks(1)) + 1 ## save a week column to the dataframe
  df$month <- as.numeric(interval(df$time[1], df$time) %/% months(1)) + 1 ## save a month column to the dataframe
  df$year <- as.numeric(interval(df$time[1], df$time) %/% years(1)) + 1  ## save an year column to the dataframe
  df$loc <- as.character(Mang$HOBO[i]) ## save the location name to the dataframe
  rndm <- c(sample(seq(1,ncol(df)-7),100),c((ncol(df)-6):ncol(df)))  ## randomly sample 100 of the points. This can be used for test runs before running all 1000 points
  listofdfs_100[[i]] <- df[,rndm]  ## save the 100 sampled waterlevel points to a list of dataframes, one for each location
  listofdfs[[i]] <- df ## save the water levels to a list of dataframes, one for each location
  listofpoints[[i]] <- as.data.frame(points)    ###  save the locations of all points
  listofpoints_100[[i]] <- as.data.frame(points)[rndm[1:(length(rndm)-7)],]   ###  save the locations of the 100 sampled points
  listofpoints_100[[i]]$loc <- as.character(Mang$HOBO[i])
  remove(df)  ## remove the current dataframe and start again for the next location
}
### Save all of the dataframes for later use
saveRDS(listofdfs,"dfs.rds")
saveRDS(listofdfs_100,"dfs_100.rds")
saveRDS(listofpoints,"points.rds")
saveRDS(listofpoints_100,"points_100.rds")
gc()  ## clean the memory


####   You now have lists of dataframes, one for each location in the study. Each dataframe contains the waterlevels at 1000 points for each timestep in the
###     study period. You also have a list of x and y coordinates for each point. This is a less heavy version that can be used later to calculate the landcover
###     around each points. Additionally, you have smaller versions of each, with only 100 random locations instead of 10. This can be used for test runs.


