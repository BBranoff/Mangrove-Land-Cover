BufferStats <- function(p,L,P,M,St,High,Cen){
  ## This routine calculates surrounding land cover at a given spatial point. 
  ## Inputs are:
  ##        A spatial point (p), this is created before running this subroutine as a random point within the mangroves
  ##        A land cover raster (L), this is a NOAA 2010 Land Cover raster for the coastal US and Puerto Rico: Office for Coastal Management, 2017. C-CAP Land Cover, Puerto Rico, 2010. Charleston, SC. https://coast.noaa.gov/dataviewer/#/imagery/search/where:ID=6201%0A
  ##        A mangrove cover raster (M), this is a data set from: Hamilton, S.; and Casey, D., 2016. Creation of a high spatiotemporal resolution global database of continuous mangrove forest cover for the 21st Century (CGMFC-21). Glob. Ecol. Biogeogr. 25.6, 729-738.
  ##        A network of roads within the study area: U.S. Census Bureau, 2015. 2015 TIGER/Line Shapefiles roads.
  ##        A list of total population numbers for census blocks (P): U.S. Census Bureau, 2010. Total population at census block level.
  ##        A shapefile of census blocks (Cen) : U.S. Census Bureau, 2010. Total population at census block level.
  
  stopifnot(length(p) == 1)  ## There must only be one point in the object
  cust <- sprintf("+proj=gnom +lat_0=%s +lon_0=%s +x_0=0 +y_0=0", 
                  p@coords[[2]], p@coords[[1]])  ## create a custom gnomic projection centered on the point. All analyses will be conducted in this CRS to ensure equal lengths/areas
  projected <- spTransform(p, CRS(cust))  #re project the point to the custom CRS
  ### Save the datasets to global environment, so that they can be used in all sub routines
  LCType <<- L  ## land cover
  Pop <<- P   #population density raster
  Mangr <<- M  #Mangrove coverage
  Streets <<- St   ## Street network
  Hwy <<- High   ### highway network
  Census <<- Cen   ## Census table
  id <- as.character(p$HOBO)  ## Site of the current point
  ### Run the LandUse subroutine to on the projected point to gather the surrounding land cover statistics
  ## the second argument to LandUse is the circle radius to create the sampling circle. The other arguments are the previously created datasets
  LandUse(projected,l=500,LCType,Mangr,p,Pop,cust,Streets,Hwy,Census)  
  
  ###  Save the collected data from LandUse to the statistical data table "pdata"
  ## the order of columns is: x coordinate, y coordinate, elevation, Site name, the number of land cover cells within the sampling circle, 
  ## the number of mangrove cells within the sampling circle, the area of the sampling circle, the area of mangroves in the sampling circle,
  ## the urban area within the sampling circle, the green and blue are within the sampling circle, the crop land within the sampling circle,
  ## the population density within the sampling circle, the total road length within the sampling circle, and the total highway length within the sampling circle.
  pdata <<- rbind(pdata,c(p@coords[[2]],p@coords[[1]],p$elev,as.character(id),
                              Count500,CellsMang500,Area500,
                              Mang500/1000000,Urban500*Area500/Count500/1000000,GB500*Area500/Count500/1000000,
                              Crop500*Area500/Count500/1000000,PopDens500/Area500/1000000,StLength,HwyLength))   
}


###  THe LandUse subroutine is the primary subroutine to gather surrounding land cover statistics at each point
###  A number of other subroutines are called from within LandUse and are given below
LandUse <- function(projected,l,LCType,Mang,p,Pop,cust,Streets,Hwy,Census){
  buffered <- gBuffer(projected, width=l, quadsegs=50)  ## create the sampling circle from the projected point, given the circle radius "l"
  assign(paste("Area",l,sep=""), as.numeric(gArea(buffered)),.GlobalEnv)  ## save the circle area to the global environment for later. Units depend on the CRS and are meters for the custom gnomic projection from the BufferStats function
  buffered <- spTransform(buffered, LCType@crs) ## With the accurate area saved, reproject the sampling circle to the native CRS of the datasets
  roads2 <- gIntersection(Streets, buffered, byid=TRUE)  ## Cut out the roads within the sampling circle
  highways2 <- gIntersection(Highways, buffered, byid=TRUE)  ## cut out the highways form within the sampling circle
  PopDens <- PopDensity(buffered,Census,Pop,cust)  ## Run the PopDensity subroutine to calculate population density
  assign(paste("PopDens",l,sep=""), as.numeric(PopDens),.GlobalEnv) ## Save the previously calculated population density to the global environment
   MangAg <- crop(Mang,spTransform(gBuffer(projected, width=l, quadsegs=50),Mang@crs),snap="out")  ## crop the mangrove coverage to the sampling circle, this speeds up the following step of masking
  MangAg <- mask(MangAg,spTransform(gBuffer(projected, width=l, quadsegs=50),Mang@crs))  ## mask the mangroves values from the sampling circle, keeping only those within the circle
  LCTypeAg <- crop(LCType,spTransform(gBuffer(projected, width=l, quadsegs=50),LCType@crs),snap="out")  ## crop the land coverage to the sampling circle, this speeds up the following step of masking
  LCTypeAg <- mask(LCTypeAg,spTransform(gBuffer(projected, width=l, quadsegs=50),LCType@crs))  ## mask the land cover values from the sampling circle, keeping only those within the circle
  MangAg <-as.data.frame(getValues(MangAg))  ###  get all of the mangrove values within the sampling circle, this includes values partially outside of the circle and is only recommended for high resolution datasets. Low resolution datasets (relative to the sampling circle size) would be much less accurate
  LCTypeAg <-as.data.frame(getValues(LCTypeAg))    ###  get all of the land cover values within the sampling circle, this includes values partially outside of the circle and is only recommended for high resolution datasets. Low resolution datasets (relative to the sampling circle size) would be much less accurate
  assign(paste("Mang",l,sep=""), as.numeric(sum(MangAg,na.rm=T)),.GlobalEnv)  ## sum the mangrove values to get the total mangrove area within the circle and save this value to the global environment
  assign(paste("CellsMang",l,sep=""), as.numeric(length(MangAg[!is.na(MangAg$`getValues(MangAg)`),])),.GlobalEnv)  ## sum the total number of mangrove cells within the circle, and save this value to the global environment. This is useful for error checking later and for calculating the spatial area of each cell
  assign(paste("Urban",l,sep=""), as.numeric(length(LCTypeAg[LCTypeAg[,1] %in% c(2,5),])),.GlobalEnv) ## sum the total number of urban classified land cover cells within the circle, and save this value to the global environment. Urban from the original dataset is defined by "Urban Developed" and "Urban Open Space" from the original dataset
  assign(paste("GB",l,sep=""), as.numeric(length(LCTypeAg[LCTypeAg[,1] %in% c(6:18,21),])),.GlobalEnv)## sum the total number of vegetation classified land cover cells within the circle, and save this value to the global environment. These are all vegetated habitats, including croplands
  assign(paste("Crop",l,sep=""), as.numeric(length(LCTypeAg[LCTypeAg[,1] %in% c(6),])),.GlobalEnv) ## sum the total number of cropland classified land cover cells within the circle, and save this value to the global environment. 
  assign(paste("Count",l,sep=""),as.numeric(length(LCTypeAg[!is.na(LCTypeAg$`getValues(LCTypeAg)`),])),.GlobalEnv) ## sum the total number of land cover cells within the circle,and save this value to the global environment. This is useful for error checking later and for calculating the spatial area of each cell
  assign("StLength",as.numeric(roadlengths(roads2,cust)),.GlobalEnv)  ## run the roadlengths subroutine to count the total road length within the sampling circle, and save this value to the global environment
  assign("HwyLength",as.numeric(roadlengths(highways2,cust)),.GlobalEnv)  ## run the roadlengths subroutine to count the total highway length within the sampling circle, and save this value to the global environment
}


PopDensity <- function(buff,Blocks,Poptemp,cust){
  ## This subroutine calculates population density within the sampling circle
  ## The algorithm assumes people only live in urban classified land that is not a road or highway
  ## this is done using the US Census block shapefile with population counts for each block saved in a separate dataset
  ## as well as the landcover dataset and the road dataset to find the total amount of non-road urban land
  ## The algorithm first cycles through each census block that falls within a sampling circle
  ## for each block, it calcualtes the population density per area of non-road developed land
  ## it saves this density to a dataset, which I recommend saving to the computer and using later, as it takes a realtively long time to calculate this density for each census block
  ## with this information saved to the data frame, it then uses the density to calculate the total population within each portion of the census blocks contained in the sampling area
  
  popcount<-0  ## start out with a population of zero. the total population of each census block will be added to this number
  buff <- spTransform(buff,CRS(as.character(Blocks@proj4string))) ### re-project the sampling circle to the CRS of the census data
  Blocks <- intersect(Blocks, buff)  ## extract the census block portions within the sampling circle
  Blocks@data$ALAND <- as.numeric(levels(Blocks@data$ALAND))[Blocks@data$ALAND]  ## get the native stored land area from the census block dataset
  Blocks@data$AWATER <- as.numeric(levels(Blocks@data$AWATER))[Blocks@data$AWATER]  ## get the native stored water area from the census block dataset
  Blocks@data$Area <- gArea(spTransform(Blocks,cust),byid = T)  ## calculate the total area of the sampling circle
  Blockdata <- Poptemp[which(Poptemp$GEOID %in% Blocks$GEOID),]  ## Get the total population within all of the census blocks inside of the sampling circle, with blocks identified by their GEOID
  Blockdata <- merge(Blockdata,Blocks@data[,c("GEOID","Area","ALAND","AWATER")])  ## merge the total population data with the previously calculated land and water areas
  for (b in 1:length(Blocks)){  ## for each of the census blocks in the native dataset...
    if (is.na(Poptemp$DensityPerDeveloped[which(Poptemp$GEOID %in% Blocks$GEOID[b])])){  ## if the block is not already in the stats table, then calculate the necessary statistics as follows. If all of the points to be analyzed are close to eachother, its likely they will include some of the same census blocks, so there's no need to calculate the stats twice 
      LCType2 <- crop(LCType,Census[Census$GEOID==Blocks$GEOID[b],], snap="out") ## Crop the landcover data to the census block 
      clip <- rasterize(Census[Census$GEOID==Blocks$GEOID[b],],LCType2, mask=TRUE)  ## extract all the land cover points within the census block 
      class<-as.data.frame(table(getValues(clip)))  ### count all of the land cover classes within the census block
      dev <- sum(class$Freq[class$Var1==2])*(Blockdata$ALAND[Blockdata$GEOID==Blocks$GEOID[b]]+Blockdata$AWATER[Blockdata$GEOID==Blocks$GEOID[b]])/sum(class$Freq,na.rm=T)  ## find the developed land area (when the land cover class equals 2) within the sampling block by counting the number of urban developed cells and multiplying by the cell area
      roads.sub <- gIntersects(spTransform(Census[Census$GEOID==Blocks$GEOID[b],],CRS(as.character(Streets@proj4string))), Streets, byid=TRUE)  ## extract all the roads in the census block
      roads.sub <- Streets[which(roads.sub > 0),]  ## get rid of zero length roads
      roads.sub <- gIntersection(roads.sub, spTransform(Census[Census$GEOID==Blocks$GEOID[b],],CRS(as.character(St@proj4string))))  ## take the network of roads and create one object, a collection of roads
      if (is(roads.sub,"SpatialCollections")){roads.sub <- roads.sub@lineobj}  ## make sure the new roads collection is a line object
      rlength <- roadlengths(roads.sub,cust)  ## calculate the total road length 
      roads.sub <- gIntersects(spTransform(Census[Census$GEOID==Blocks$GEOID[b],],CRS(as.character(Highways@proj4string))), Highways, byid=TRUE)  ## extract all the highways in the census block
      roads.sub <- Highways[which(roads.sub > 0),] ## get rid of zero length highways
      roads.sub <- gIntersection(roads.sub, spTransform(Census[Census$GEOID==Blocks$GEOID[b],],CRS(as.character(Highways@proj4string))))  ## take the network of highways and create one object, a collection of highways
      if (is(roads.sub,"SpatialCollections")){roads.sub <- roads.sub@lineobj}  ## make sure the new highways collection is a line object
      hlength <- roadlengths(roads.sub,cust)   ## calculate the total highways length 
      dev <- dev - (rlength-hlength)*6 - hlength*16  ## calculate the are of roads (assumed 6 meters wide) and highways (assumed 16 meters wide) and subtract the total roads and highway areas from the developed area in the census block, to give the total area of non-roads developed land
      if (dev == 0|dev <0){ ## if there was no developed area int he census block, or if there was only roads in the cenus block (because our assumption of road and highway widths is slightly off), the non-road developed land is 0
        Poptemp[Poptemp$GEOID ==Blocks$GEOID[b] ,"DensityPerDeveloped"] <-0  ## save the value to the popualtion density dataset for each census block
      }else{   ## Otherwise, the popualtion density is the total number of people in the census block divided by the non road developed land area
        Poptemp[Poptemp$GEOID ==Blocks$GEOID[b] ,"DensityPerDeveloped"] <- Poptemp[Poptemp$GEOID == Blocks$GEOID[b],"PopTotal"]/(dev)  ## save the value to the popualtion density dataset for each census block
      }
      stopifnot(Poptemp[Poptemp$GEOID ==Blocks$GEOID[b] ,"DensityPerDeveloped"] < 0|Poptemp[Poptemp$GEOID ==Blocks$GEOID[b] ,"DensityPerDeveloped"]=="NaN")  ## if the population density was less than zero or if it results in "NaN", something is wrong and an error is thrown
      Pop <<- Poptemp  ## save the calculate population density to the global environment
    }
    ### with the population density aready calculated for that block, calculate the population density pertaining to the portion of the block within the sampling circle 
    LCType2 <- crop(LCType,Blocks[b,], snap="out")  ## crop the land cover dataset to the sampling circle
    dens <- Poptemp$DensityPerDeveloped[Pop$GEOID==Blocks$GEOID[b]]  ## get the previously calcualted density for the current censusblock
    clip <- rasterize(Blocks[b,],LCType2, mask=TRUE)  ## Extract the landcover cells within the portion of the census block within the sampling circle
    class<-as.data.frame(table(getValues(clip)))  ## calculate the frequency of each land cover class within that portion of census block
    dev <- sum(class$Freq[class$Var1==2])  ##  sum the number of urban developed cells within the portion of census block 
    roads.sub <- gIntersects(spTransform(Blocks[b,],CRS(as.character(Streets@proj4string))), Streets, byid=TRUE)  ## extract the roads within the portion of census block
    roads.sub <- Streets[which(roads.sub > 0),] ## get rid of roads with zero length
    roads.sub <- gIntersection(roads.sub, spTransform(Blocks[b,],CRS(as.character(Streets@proj4string))))   ## take the network of roads and create one object, a collection of roads
    if (is(roads.sub,"SpatialCollections")){roads.sub <- roads.sub@lineobj}  ## make sure the new roads collection is a line object
    rlength <- roadlengths(roads.sub,cust) ## calculate the total road length
    roads.sub <- gIntersects(spTransform(Blocks[b,],CRS(as.character(Highways@proj4string))), Highways, byid=TRUE) ## extract the highways within the portion of census block
    roads.sub <- Highways[which(roads.sub > 0),]  ## get rid of highways with zero length
    roads.sub <- gIntersection(roads.sub, spTransform(Blocks[b,],CRS(as.character(Highways@proj4string))))  ## take the network of highways and create one object, a collection of highways
    if (is(roads.sub,"SpatialCollections")){roads.sub <- roads.sub@lineobj}  ## make sure the new roads collection is a line object
    hlength <- roadlengths(roads.sub,cust)  ## calculate the total highway length
    dev <- dev - (rlength-hlength)*6 - hlength*16   ## calculate the area of roads (assumed 6 meters wide) and highways (assumed 16 meters wide) and subtract the total roads and highway areas from the developed area in the portion of the census block, to give the total area of non-roads developed land
    if (dev <0){dev=0}  ## if the total non-road developed area is less than 0 (because our assumption of road and highway widths is slightly off), then set it to 0
    popcount <- popcount+dens*dev  ## The total population count is then the previous population count plus the population density per non-road developed area times the amount of non-road developed area
  }
  popcount  ## return the total population count claculated for each portion of census block within the sampling circle.
}



roadlengths <- function(r,cust){
  ## this subroutine calculates road lengths
  ## it is basically a check to make sure their are not 0 roads
  ## because this has to be done often, it is delegated to a subroutine to be easily repeated
  if (is.null(r)){
    return(0)
  }else{
    r <-  spTransform(r, cust)#p@proj4string
    r <- as(r, "SpatialLines")
    return(as.numeric(sum(SpatialLinesLengths(r))))
  }
}

