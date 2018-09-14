# Mangrove-Land-Cover
A collection of R language routines from calculating land cover surrounding mangrove locations

This collection is meant to accompany the manuscript entitled "Urbanization plays a minor role in the flooding and surface water chemistry of Puerto Rico’s mangroves", currently in review. The manuscript describes the modeling of surface water levels in mangroves across Puerto Rico, and the use of these models to calculate flooding metrics within the mangroves. It also describes how surrounding land cover at each site is assessed to quantify the level of urbanization around the mangroves. The flooding metrics are then used along with surface water chemistry measurements in correlations with the sampled land cover. The manuscript describes how urbanization influences the flooding and surface water chemistry of the mangroves of Puerto Rico.

As Figure 2 of the manuscript portrays, there are four basic steps to the process. The first is to calculate the water level models, then distribute points randomly throughout the mangroves of Puerto Rico and calculate their water levels based on elevation and the model, then calculate the flooding metrics from the water levels and finally sample the land cover surrounding each site. Subroutines are included for each of these processes as described below:

Water Level Models: This routine takes as an input water level observations from a piezometer in the mangroves, as well as rainfall observations from a nearby weather station. The algorithm first calculates the tidal harmonics from the water leve observations, then the contributions from rainfall. It creates a final water level model from the addition of the tidal and rainfall contributions.

Mangrove Points: This routine takes as inputs shapefiles of mangrove zones in Puerto Rico, each zone corresponding to a piezometer, as well as the point locations of the piezometers. It also takes digitial elevation models (DEMS) for the area and randomly samples elevations within the mangrove zones. It then uses the observed water levels at the piezometers and the difference in elevation between the random point and the piezometer to calculate the water level at the random point.

Flooding Stats:  With the water levels calculated at all points, this routine returns flooding metrics for each point. Metrics include hydroperiod (flood duration) flood depth, flood frequency etc.

Land Cover Stats:  This routine uses a number of spatial inputs to sample the surrounding land cover at each points. It extracts urban and green & blue (vegetation + open water) coverage, population density, road lengths, and mangrove coverage surrounding each point. 

I have included files small enough for GitHub, and others are available for free online following the commented sources in the R scripts. 
