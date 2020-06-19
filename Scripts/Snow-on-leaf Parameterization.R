##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script estimates autumn phenology in 2011 to determine vegetation status during the 2011 Snowtober 
#snow-on-leaf event. 

#The phenophase calculated is for 90% leaf coloration, or peak leaf coloration, using a
#cooling degree day model. Model parameters were based on Jeong and Medvigy, who used leaf coloration data from 
#the Harvard Forest and from the National Phenology Network and daily average temperature as model inputs

#Input data are from DAYMET (https://daymet.ornl.gov/)

#scripts to download data and calculate daily average temperature from tmin and tmax use the
#daymetr package (https://github.com/bluegreen-labs/daymetr)

library(daymetr)
library(lubridate)
library(raster)
library(rgdal)
library(rgeos)
library(stringr)

#############################
#Define extent of study area#
#############################

#read in shapefile
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects")
study.region = readOGR(dsn = "Spatial_Domain/", layer = "NA_Eco_East")

ex = extent(study.region)

######################################################
#determine daymet tiles needed that fit within extent#
######################################################

download_daymet_tiles(location = c(ex[3], ex[1], ex[4], ex[2]),
                     start = 2011,
                     end = 2011,
                     param = c("tmin"),
                     path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Tile_Names",
                     silent = TRUE)

#note that the code will only download data for one tile. That's OK. We are using it to identify
#tiles within the study area for which we will will download daymet data
files.1 = list.files("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Tile_Names")
files.2 = data.frame(do.call(rbind, str_split(files.1, "_")))
files.3 = data.frame(do.call(rbind, str_split(files.2$X3, ".nc")))
files.4 = files.3$X1

#now download daymet tiles one at a time

#tmin
for(i in 1:length(files.4)){
download_daymet_tiles(tiles =  files.4[i],
                      start = 2011,
                      end = 2011,
                      param = c("tmin"),
                      path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMIN",
                      silent = TRUE)

}


#tmax
for(i in 1:length(files.4)){
  download_daymet_tiles(tiles =  files.4[i],
                        start = 2011,
                        end = 2011,
                        param = c("tmax"),
                        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMAX",
                        silent = TRUE)
  
}

########################################################
#calculate average daily temperature for each grid cell# 
#in each tile using daymet_grid_tmean                  #
########################################################

#list files in tmin and tmax
tmin.files = list.files("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMIN")
tmax.files = list.files("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMAX")

#run loop
for(i in 1:length(files.4)){
  
  tile = as.character(unique(files.4[i]))
    
  #copy files from tmin and tmax into a temporary directory
  setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMIN")
  file.copy(tmin.files[i], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Temp_Dir")
  setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMAX")
  file.copy(tmax.files[i], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Temp_Dir")
  
  #calculate average temperature
  tmean = daymet_grid_tmean(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Temp_Dir",
                               product = tile,
                               year = 2011,
                               internal = TRUE)

  out_name = paste("tmean_2011_", tile, ".tif", sep = "")
  
  #set working directory and export average temperature data
  setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMEAN")
  writeRaster(tmean, out_name, format = "GTiff")
 
  #set working directory and remove files to make space for the next pair
  setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Temp_Dir")
  
  file.remove(tmin.files[i])
  file.remove(tmax.files[i])
  
    }

#############################
#run autumn phenology models# 
#############################

#list files in TMEAN folder
tmean.files = list.files("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMEAN")
#split file names to separate .aux.xml files from .tif files
tmean.files.2 = data.frame(do.call(rbind, str_split(tmean.files, ".aux")))
#subset to only include .tif tiles
tmean.files.3 = tmean.files.2[tmean.files.2$X2 != ".xml", ]

#extract the names of individual tiles
tiles.1 = data.frame(do.call(rbind, str_split(tmean.files.3$X1, "_")))
tiles.2 = data.frame(do.call(rbind, str_split(tiles.1$X3, ".tif")))

#begin loop

for(i in 1:nrow(tmean.files.3)){

  #define model parametera for three species groups as delieated by Jeong and Medvidgy (2014)
  #group 1 = A rubrum, A saccharum, F grandifolia
  #group 2 = Q rubra, Q velutina
  #group 3 = B lenta, B papyrifera
  
  #define intial day from which to calculate cdd
  g1doy = 244:365
  g2doy = 203:365
  g3doy = 191:365
  
  #define temperature threshold
  g1temp = 25
  g2temp = 30
  g3temp = 30
  
  #define CDD threshold for peak foliage coloration
  g1thresh = -640
  g2thresh = -672
  g3thresh = -732
  
  #define functions for determining whether tmean is below the degree day threshold
  ftemp1 = function(x, ...) {ifelse(x < g1temp, x - g1temp, 0)}
  ftemp2 = function(x, ...) {ifelse(x < g2temp, x - g2temp, 0)}
  ftemp3 = function(x, ...) {ifelse(x < g3temp, x - g3temp, 0)}
  
  #define functions for determining whether cdd is below the cumulative degree day threshold
  #fcdd1 = function(x, ...) {ifelse(x <= g1thresh, 0, 1)}
  #fcdd2 = function(x, ...) {ifelse(x <= g2thresh, 0, 1)}
  #fcdd3 = function(x, ...) {ifelse(x <= g3thresh, 0, 1)}
  
  #define g1 layer names. Scrip sometimes get stuck on g1 layer names that are not always consistent
  nums = c(58:179)
  g1names = paste("layer.", nums, sep = "")

#set directory for importing files
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/TMEAN/")

#name the tile to import
tilename = (as.character(unique(tmean.files.3$X1[i])))
    
#extract the tile name
tile = (as.character(unique(tiles.2$X1[i])))

#import the tile as a raster stack  
rs = stack(tilename)

#subset raster stack for g1, g2, and g3. The [[]] notation pulls invidual raster bands that 
#correspond to days of year
g1sub = rs[[g1doy]]
g2sub = rs[[g2doy]]
g3sub = rs[[g3doy]]

#replace names in g1sub
names(g1sub) = g1names
  
#calculate the difference between tmean and temperature thresholds
g1temp = calc(g1sub, fun = ftemp1)
g2temp = calc(g2sub, fun = ftemp2)
g3temp = calc(g3sub, fun = ftemp3)

#calculate the cumulative cooling degree days (cumulative difference between tmean and 
#temperature threhsold)
g1cdd = calc(g1temp, fun = cumsum)
g2cdd = calc(g2temp, fun = cumsum)
g3cdd = calc(g3temp, fun = cumsum)

#extract bands for October 28, 2011, just prior to the arrivial of the October nor'easter.
#this corresponds to layer.58 for g1 (i.e., 58 days from 9/1 or DOY 244), layer.99 for g2
#and layer.111

#
g1eve = g1cdd$layer.58
g2eve = g2cdd$layer.99
g3eve = g3cdd$layer.111

#reclassify cdd vales to be either below (1) or above (0) cdd threshold
#g1stat = calc(g1eve, fun = fcdd1)
#g2stat = calc(g2eve, fun = fcdd2)
#g3stat = calc(g3eve, fun = fcdd3)

#define file output names for each group
g1out_name = paste("g1eve_", tile, ".tif", sep = "")
g2out_name = paste("g2eve_", tile, ".tif", sep = "")
g3out_name = paste("g3eve_", tile, ".tif", sep = "")

#write reclassified rasters to folders for g1, g2, and g3

#set working directory for g1
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Maple_Beech")
writeRaster(g1eve, g1out_name, format = "GTiff")

#set working directory for g2
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Oak")
writeRaster(g2eve, g2out_name, format = "GTiff")

#set working directory for g3
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/DAYMET/Birch")
writeRaster(g3eve, g3out_name, format = "GTiff")

}

