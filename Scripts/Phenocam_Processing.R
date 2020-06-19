##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script estimates model parameters used to simulate phenology from daily temperature
#using the phenor and phenocamr packages (https://khufkens.github.io/phenor/) developed
#by K. Hufkens et al. (2018) and modified for identifying winter weather whiplash events

####################################################################################
#Initial set up                                                                   
####################################################################################

#install packages and load libraries

install.packages("devtools")
library(devtools)

install.packages("phenocamr")
library(phenocamr)

#note that the phenor package here is the most recent (under current development)
devtools::install_github("khufkens/phenor@refactoring_v2")
library(phenor)                

library(maps)
library(plyr)
library(raster)
library(rgdal)
library(RJSONIO)
library(stringr)


####################################################################################
#Select phenocams to use in the model parameterization                                                                 
####################################################################################

#Import land cover classifications for each weather station
station.lc = read.table("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Land_Cover/station_lc.csv", head = T, sep = ",")

##########################
#Import phenocam metadata 
##########################

#File is in .json format and has to be converted to data.frame
#Be sure to have file path to main project directory in place

#set working directory
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data")

#import file
con = file('site_information.json', "r")

#use ldply to split list, apply function, and return results in dataframe
pcmeta  <- ldply(fromJSON(con), data.frame)

#close connection to file
close(con)

##################
#Select phenocams
##################

#define spation extent

#use level II ecoregion classification to select northern and eastern temperate forest types in 
#North America. Defines spatial extent.

#read in shapefile 
study.region = readOGR(dsn = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Spatial_Domain/", layer = "NA_Eco_Level2_Clip_SSC")

#set the extent of the whole study region
ex = extent(study.region)

#select stations that fit within extent (ex)

pcsub = pcmeta[pcmeta$lat >= ex[3] &
                 pcmeta$lat <= ex[4] &
                 pcmeta$lon >= ex[1] &
                 pcmeta$lon <= ex[2] , ]

#examine IGBP land cover for pcsub and compare of IGBP land cover for station.lc
table(pcsub$landcover_igbp)
table(station.lc$stalcc)

#both have lc classes 4 (deciduous broadleaf forest), 5 (mixed forest), 10 (grassland), 
#13 (urban), and 14 (cropland / vegetation mosaics). The weather stations also include 
#8 (woody savannas), 9 (savannas), and 11 (wetlands)

#some of these (8,9, 10, and 13) wil be reclassified as mixed forest for the purpose of 
#phenology modeling. Grasslands and cropland / vegetation mosaics will be reclassified as cropland
station.lc$lc = ifelse(station.lc$stalcc == 8 | station.lc$stalcc == 9 | station.lc$stalcc == 11 | station.lc$stalcc == 13, 5,
                       station.lc$stalcc)

#subset station list into IGBP land cover classification
#"DB" = deciduous broadleaf (4), "MX" = mixed forest (5), "GL" = grassland (10), "UR" = urban (13),
#"AG" = mixed cropland and natural vegetation

pcdbf = pcsub[pcsub$landcover_igbp == 4, ]
pcmxf = pcsub[pcsub$landcover_igbp == 5, ]
pcgrl = pcsub[pcsub$landcover_igbp == 10, ]
pcurb = pcsub[pcsub$landcover_igbp == 13, ]
pcagr = pcsub[pcsub$landcover_igbp == 14, ]

#extract station names from subsets
dbf = as.character(pcdbf$site)
mxf = as.character(pcmxf$site)
grl = as.character(pcgrl$site)
urb = as.character(pcurb$site)
agr = as.character(pcagr$site)

####################################################################################
#Download phenocam data for each veg type  
#Placing each veg type in a separate folder
#Using the code (out_dir = ".") causes the downloaded data, both the 3-day time series 
#and the calculated transition dates, to be stored in current working directory.
####################################################################################

##########################
#Deciduous Broadleaf (dbf)
##########################

setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/DBF")

for(i in 1:length(dbf)){
  try(
    phenocamr::download_phenocam(
      frequency = 3,
      veg_type = "DB",
      roi_id = "1000",
      site = dbf[i],
      phenophase = TRUE,
      out_dir = "."
    )  
  )
  if(inherits(phenocamr::download_phenocam, "try-error"))next
}

##########################
#Mixed Forest (mxf)
##########################

setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/MXF")

for(i in 1:length(mxf)){
  try(
    phenocamr::download_phenocam(
      frequency = 3,
      veg_type = "DB",
      roi_id = "1000",
      site = mxf[i],
      phenophase = TRUE,
      out_dir = "."
    )  
  )
  if(inherits(phenocamr::download_phenocam, "try-error"))next
}

##########################
#Urban (urb)
##########################

setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/URB")

for(i in 1:length(urb)){
  try(
    phenocamr::download_phenocam(
      frequency = 3,
      veg_type = "DB",
      roi_id = "1000",
      site = urb[i],
      phenophase = TRUE,
      out_dir = "."
    )  
  )
  if(inherits(phenocamr::download_phenocam, "try-error"))next
}

##########################
#Agriculture (agr)
##########################

setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/AGR")

for(i in 1:length(agr)){
  try(
    phenocamr::download_phenocam(
      frequency = 3,
      veg_type = "AG",
      roi_id = "1000",
      site = agr[i],
      phenophase = TRUE,
      out_dir = "."
    )  
  )
  if(inherits(phenocamr::download_phenocam, "try-error"))next
}

####################################################################################
#Format phenocam transition date data, specifying the direction of the curve (rising 
#or falling), the gcc percentile (10, 50, or 90), an the temporal offset
#must be done separately for each vegetation type, direction, and gcc threshold (n = 18) 
####################################################################################

##########################
#Deciduous Broadleaf (dbf)
##########################

#set working directory
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/DBF")

#read files names in working directory
files.dbf = list.files(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/DBF", full.names = F)

#split file names to ID the ones that have raw greenness data
#this will allow for the sequential selecting of stations in the loop below
#R will throw a warning because not all items in data frame have same length (no worries about this)
dbf_spl.1 = data.frame(do.call(rbind, str_split(files.dbf, "_")))

#create empty containers to hold the lists that will be generated for each threshold (10, 50, 90) and
#each time period (rising, falling)
dbf_rising_10 = list()
dbf_rising_50 = list()
dbf_rising_90 = list()

dbf_falling_10 = list()
dbf_falling_50 = list()
dbf_falling_80 = list()

#loop through pairs of files from each station in DBF. Each pair contains the raw gcc data and 
#predetermined transition dates. The code will pull these pairs and place them in a temporary directory.
#Once in the directory, the pr_fm_phenocam will try to format the data to prepare it for downstream modeling.
#the "try" argument will keep the loop going even if the function throws and error, allowing it
#to move to the next pair of files in DBF. Some of the code will alter the file with transition dates
#so instead of the predetermined ones of 10, 25, and 50 percentile gcc, it does 10, 50, and 90 gcc for the "rising"
#or spring phenophases, and 10, 50, and 80 for the "falling" or fall phenophases

#begin outer loop 

#########################################
#select files and place them in directory
#########################################

for(i in 1:nrow(dbf_spl.1)){
  #choose files from deciduous broadleaf veg_type and 
  if(dbf_spl.1$X6[i] == "DB"){
  #set working directory to get these files
  setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/DBF")
 
  #copy the first two files from the DBF directory to the temporary directory. This pair
  #contains the raw gcc data and predetermined transition dates. 
  file.copy(files.dbf[i], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
  file.copy(files.dbf[i+1], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")

  #record the file names in the temporary directory
  files.temp = list.files(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir/", full.names = F)
  
  ###########################################
  #format files to alter the transition dates
  ###########################################
  
  #set working directory to temporary directory. IMPORTANT to do this because many of the functions in 
  #phenor use the "." syntax to refer to the home directory or have this embedded in the background code
  #to enable the processing of files. Things will end up in the wrong place and the loop will
  #likely break of the working directory is not specified properly!
  setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
  
  # read in the transition date file (the second file in the temporary directory)
  td = read.table(files.dbf[i+1],
                    header = TRUE,
                     sep = ",")
  
  # smooth data using an optimization approach; forcing (force = TRUE) allows for the skipping of the 
  # procedure if the smoothed data are already available
  smooth_ts(as.character(files.dbf[i]),
            out_dir = ".",
            force = TRUE)
  
  # generate phenological transition dates (10, 50, and 90 percentile thresholds). Results in a list
  # with two items: a list of gcc thresholds with the transition dates for the rising and falling limbs
  # of the phenophases. For this, we are interested in the rising limb
  td.1 = phenophases(as.character(files.dbf[i]),
                     internal = TRUE,
                     lower_thresh = 0.1,
                     middle_thresh = 0.5,
                     upper_thresh = 0.9)
  
  #take the rising part of the list, convert to a data.frame, and remove row names
  rising.1 = do.call(rbind.data.frame, td.1[1])
  rownames(rising.1) = NULL
  
  #NOTE: not all of the phenocams have data for both rising and falling phenophases. The next few lines of
  #code enable the delineation of rising as "NA" if no data exist
  
  #determine if data.frame is empty
  if(is.data.frame(rising.1) && nrow(rising.1)!=0){
    #if false
    #convert columns 2 to 10 (estimated phenophases with 95% CIs) to dates
    rising.1[2:10] = lapply(rising.1[2:10], as.Date, origin = "1970-01-01")
    #create columns containing site, veg_type, roi, and phenophase direction info
    #these are the columns that appear in the original transition date (td) file
    #that phenor will be expecting when it runs pr_fit to estimate phenophases below.
    #the list (td.1) that we created does not do this automatically
    rinfo = data.frame("site" = unique(td$site), "veg_type" = unique(td$veg_type), "roi_id" = 
                         unique(td$roi_id), "direction" = "rising")
    #make rising.2 by combining rinfo with rising.1
    rising.2 = cbind(rinfo, rising.1)}else
    
    #if data.frame is empty, make rising.2 NA
    rising.2 = NA   
  
  # generate phenological transition dates (10, 50, and 80 percentile thresholds). Results in a list
  # with two items: a list of gcc thresholds with the transition dates for the rising and falling limbs
  # of the phenophases. For this, we are interested in the falling limb
  td.2 = phenophases(as.character(files.dbf[i]),
                     internal = TRUE,
                     lower_thresh = 0.1,
                     middle_thresh = 0.5,
                     upper_thresh = 0.8)
  
  #take the falling part of the list, convert to a data.frame, and remove row names
  falling.1 = do.call(rbind.data.frame, td.2[2])
  rownames(falling.1) = NULL
  
  #NOTE: not all of the phenocams have data for both falling and falling phenophases. The next few lines of
  #code enable the delineation of falling as "NA" if no data exist
  
  #determine if data.frame is empty
  if(is.data.frame(falling.1) && nrow(falling.1)!=0){
    #if false
    #convert columns 2 to 10 (estimated phenophases with 95% CIs) to dates
    falling.1[2:10] = lapply(falling.1[2:10], as.Date, origin = "1970-01-01")
    #create columns containing site, veg_type, roi, and phenophase direction info
    #these are the columns that appear in the original transition date (td) file
    #that phenor will be expecting when it runs pr_fit to estimate phenophases below.
    #the list (td.1) that we created does not do this automatically
    rinfo = data.frame("site" = unique(td$site), "veg_type" = unique(td$veg_type), "roi_id" = 
                         unique(td$roi_id), "direction" = "falling")
    #make falling.2 by combining rinfo with falling.1
    falling.2 = cbind(rinfo, falling.1)
    #force names of falling.2 to be the same as rising.2, otherwise they cannot be combined
    names(falling.2) = names(rising.2)
    }else
      
      #if data.frame is empty, make falling.2 NA
      falling.2 = NA   
  
  #combine rising.2 and falling.2 into a new data.frame only if both are present. 
  #otherwise, make the new data.frame out of either rising.2 (falling.2 is NA)
  #or falling.2 (rising.2 is NA)
  #note that R throws a warning saying there are multiple conditions and it will only use the first one
  #there are not multiple conditions. R is just being cranky. Don't worry about the warnings
  
  #determine if rising.2 = NA
    if(is.na(rising.2) == T) {
      #if yes, td.3 = falling.2
      td.3 = falling.2}
      #determine if falling.2 = NA
     else if(is.na(falling.2) == T) {
        #if yes, td.3 = rising.2
        td.3 = rising.2
    } else  
    #otherwise td.2 is the combination of rising.2 and falling.2
     td.3 = rbind(rising.2, falling.2)
  #write td.3 to the temporary directory, using the name name as the transition date file
  #this will overwrite the original file, which is what we want to have the user-defined
  #transition dates
  write.csv(td.3, files.dbf[i+1])
  
  ###################################################################################
  #format phenocam data to include DAYMET variables and user-defined transition dates
  ###################################################################################
 
  for(j in 1:length(files.temp)){
    
  # formate phenocam
  pc_form_rising_10 = try(pr_fm_phenocam(
    path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
    direction = "rising",
    gcc_value = "gcc_90",
    threshold = 10,
    offset = 264,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_rising_10, "try-error"))
  
    next

  pc_form_rising_50 = try(pr_fm_phenocam(
    path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
    direction = "rising",
    gcc_value = "gcc_90",
    threshold = 50,
    offset = 264,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_rising_50, "try-error"))
    
    next
  
  pc_form_rising_90 = try(pr_fm_phenocam(
    path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
    direction = "rising",
    gcc_value = "gcc_90",
    threshold = 90,
    offset = 264,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_rising_90, "try-error"))
    
    next
  
  pc_form_falling_10 = try(pr_fm_phenocam(
    path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
    direction = "falling",
    gcc_value = "gcc_90",
    threshold = 10,
    offset = 365,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_falling_10, "try-error"))
    
    next
  
  pc_form_falling_50 = try(pr_fm_phenocam(
    path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
    direction = "falling",
    gcc_value = "gcc_90",
    threshold = 50,
    offset = 365,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_falling_50, "try-error"))
    
    next
  
  pc_form_falling_80 = try(pr_fm_phenocam(
    path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
    direction = "falling",
    gcc_value = "gcc_90",
    threshold = 90,
    offset = 365,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_falling_80, "try-error"))
    
    next
  
  #append formatted phenocam files to lists
  
  dbf_rising_10 = append(dbf_rising_10, pc_form_rising_10)
  dbf_rising_50 = append(dbf_rising_50, pc_form_rising_50)
  dbf_rising_90 = append(dbf_rising_90, pc_form_rising_90)
  dbf_falling_10 = append(dbf_falling_10, pc_form_falling_10)
  dbf_falling_50 = append(dbf_falling_50, pc_form_falling_50)
  dbf_falling_80 = append(dbf_falling_80, pc_form_falling_80)
  
  
  }
  
  #set working directory and remove files to make space for the next pair
  setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
  
  file.remove(files.dbf[i])
  file.remove(files.dbf[i+1])
  
  #end loop
  }else
  next
}

#remove null elements from list (in case stations had null values for both rising and falling
#during the process of assigning phenophases)
#dbf_ris_10 = dbf_rising_10[-which(lapply(dbf_rising_10,is.null) == T)]
#dbf_ris_50 = dbf_rising_50[-which(lapply(dbf_rising_50,is.null) == T)]
#dbf_ris_90 = dbf_rising_90[-which(lapply(dbf_rising_90,is.null) == T)]
 
#dbf_fal_10 = dbf_falling_10[-which(lapply(dbf_falling_10,is.null) == T)]
#dbf_fal_50 = dbf_falling_50[-which(lapply(dbf_falling_50,is.null) == T)]
#dbf_fal_90 = dbf_falling_90[-which(lapply(dbf_falling_90,is.null) == T)]
  
#determine the length of the list in order to remove every other item
#loop writes each station twice and needs to be halved
#not super elegant programming solution, but it works, and pr_fit will not
#operate with duplicates!

len_lis_dbf = seq(2, length(dbf_rising_10), by = 2)

#remove every other item (indexed as len_list_dbf) by making it NULL
#rename each list prior to doing this so as not to lose list completely (paranoid move, but 
#safer choice than overwriting each list and having to rerun again)

dbf_ris_10 = dbf_rising_10
dbf_ris_10[len_lis_dbf] = NULL

dbf_ris_50 = dbf_rising_50
dbf_ris_50[len_lis_dbf] = NULL

dbf_ris_90 = dbf_rising_90
dbf_ris_90[len_lis_dbf] = NULL

dbf_fal_10 = dbf_falling_10
dbf_fal_10[len_lis_dbf] = NULL

dbf_fal_50 = dbf_falling_50
dbf_fal_50[len_lis_dbf] = NULL

dbf_fal_80 = dbf_falling_80
dbf_fal_80[len_lis_dbf] = NULL

#######################################################################
#use the pr_fit function to select the parameter set for each threshold
#(10, 50, 90) and time period (rising and falling)
#######################################################################

dbf_par_ris_10 = pr_fit(data = dbf_ris_10,
                   model = "TT",
                   method = "GenSA",
                   par_ranges = system.file("extdata",
                  "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                   plot = TRUE)

dbf_par_ris_50 = pr_fit(data = dbf_ris_50,
                        model = "TT",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                        "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

dbf_par_ris_90 = pr_fit(data = dbf_ris_90,
                        model = "TT",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                        "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

dbf_par_fal_10 = pr_fit(data = dbf_fal_10,
                        model = "CDD",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                        "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

dbf_par_fal_50 = pr_fit(data = dbf_fal_50,
                        model = "CDD",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                        "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

dbf_par_fal_80 = pr_fit(data = dbf_fal_80,
                        model = "CDD",
                        method = "GenSA",
                        par_ranges = system.file("extdata",
                        "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

#combine phenology parameters into a single list
dbf_pars = append(list(dbf_par_ris_10), list(dbf_par_ris_50))
dbf_pars = append(dbf_pars, list(dbf_par_ris_90))
dbf_pars = append(dbf_pars, list(dbf_par_fal_10))
dbf_pars = append(dbf_pars, list(dbf_par_fal_50))
dbf_pars = append(dbf_pars, list(dbf_par_fal_80))

#setwd and export list
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling")
list.save(dbf_pars, "dbf_pars.rdata")

##########################
##########################
#Mixed Forest (mxf)
##########################
##########################

#set working directory
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/MXF")

#read files names in working directory
files.mxf = list.files(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/MXF", full.names = F)

#split file names to ID the ones that have raw greenness data
#this will allow for the sequential selecting of stations in the loop below
#R will throw a warning because not all items in data frame have same length (no worries about this)
mxf_spl.1 = data.frame(do.call(rbind, str_split(files.mxf, "_")))

#create empty containers to hold the lists that will be generated for each threshold (10, 50, 90) and
#each time period (rising, falling)
mxf_rising_10 = list()
mxf_rising_50 = list()
mxf_rising_90 = list()

mxf_falling_10 = list()
mxf_falling_50 = list()
mxf_falling_80 = list()

#loop through pairs of files from each station in MXF. Each pair contains the raw gcc data and 
#predetermined transition dates. The code will pull these pairs and place them in a temporary directory.
#Once in the directory, the pr_fm_phenocam will try to format the data to prepare it for downstream modeling.
#the "try" argument will keep the loop going even if the function throws and error, allowing it
#to move to the next pair of files in MXF. Some of the code will alter the file with transition dates
#so instead of the predetermined ones of 10, 25, and 50 percentile gcc, it does 10, 50, and 90 gcc for the "rising"
#or spring phenophases, and 10, 50, and 80 for the "falling" or fall phenophases

#begin outer loop 

#########################################
#select files and place them in directory
#########################################

for(i in 1:nrow(mxf_spl.1)){
  #choose files from deciduous broadleaf veg_type and 
  if(mxf_spl.1$X6[i] == "DB"){
    #set working directory to get these files
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/MXF")
    
    #copy the first two files from the MXF directory to the temporary directory. This pair
    #contains the raw gcc data and predetermined transition dates. 
    file.copy(files.mxf[i], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    file.copy(files.mxf[i+1], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    #record the file names in the temporary directory
    files.temp = list.files(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir/", full.names = F)
    
    ###########################################
    #format files to alter the transition dates
    ###########################################
    
    #set working directory to temporary directory. IMPORTANT to do this because many of the functions in 
    #phenor use the "." syntax to refer to the home directory or have this embedded in the background code
    #to enable the processing of files. Things will end up in the wrong place and the loop will
    #likely break of the working directory is not specified properly!
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    # read in the transition date file (the second file in the temporary directory)
    td = read.table(files.mxf[i+1],
                    header = TRUE,
                    sep = ",")
    
    # smooth data using an optimization approach; forcing (force = TRUE) allows for the skipping of the 
    # procedure if the smoothed data are already available
    smooth_ts(as.character(files.mxf[i]),
              out_dir = ".",
              force = TRUE)
    
    # generate phenological transition dates (10, 50, and 90 percentile thresholds). Results in a list
    # with two items: a list of gcc thresholds with the transition dates for the rising and falling limbs
    # of the phenophases. For this, we are interested in the rising limb
    td.1 = phenophases(as.character(files.mxf[i]),
                       internal = TRUE,
                       lower_thresh = 0.1,
                       middle_thresh = 0.5,
                       upper_thresh = 0.9)
    
    #take the rising part of the list, convert to a data.frame, and remove row names
    rising.1 = do.call(rbind.data.frame, td.1[1])
    rownames(rising.1) = NULL
    
    #NOTE: not all of the phenocams have data for both rising and falling phenophases. The next few lines of
    #code enable the delineation of rising as "NA" if no data exist
    
    #determine if data.frame is empty
    if(is.data.frame(rising.1) && nrow(rising.1)!=0){
      #if false
      #convert columns 2 to 10 (estimated phenophases with 95% CIs) to dates
      rising.1[2:10] = lapply(rising.1[2:10], as.Date, origin = "1970-01-01")
      #create columns containing site, veg_type, roi, and phenophase direction info
      #these are the columns that appear in the original transition date (td) file
      #that phenor will be expecting when it runs pr_fit to estimate phenophases below.
      #the list (td.1) that we created does not do this automatically
      rinfo = data.frame("site" = unique(td$site), "veg_type" = unique(td$veg_type), "roi_id" = 
                           unique(td$roi_id), "direction" = "rising")
      #make rising.2 by combining rinfo with rising.1
      rising.2 = cbind(rinfo, rising.1)}else
        
        #if data.frame is empty, make rising.2 NA
        rising.2 = NA   
    
    # generate phenological transition dates (10, 50, and 80 percentile thresholds). Results in a list
    # with two items: a list of gcc thresholds with the transition dates for the rising and falling limbs
    # of the phenophases. For this, we are interested in the falling limb
    td.2 = phenophases(as.character(files.mxf[i]),
                       internal = TRUE,
                       lower_thresh = 0.1,
                       middle_thresh = 0.5,
                       upper_thresh = 0.8)
    
    #take the falling part of the list, convert to a data.frame, and remove row names
    falling.1 = do.call(rbind.data.frame, td.2[2])
    rownames(falling.1) = NULL
    
    #NOTE: not all of the phenocams have data for both falling and falling phenophases. The next few lines of
    #code enable the delineation of falling as "NA" if no data exist
    
    #determine if data.frame is empty
    if(is.data.frame(falling.1) && nrow(falling.1)!=0){
      #if false
      #convert columns 2 to 10 (estimated phenophases with 95% CIs) to dates
      falling.1[2:10] = lapply(falling.1[2:10], as.Date, origin = "1970-01-01")
      #create columns containing site, veg_type, roi, and phenophase direction info
      #these are the columns that appear in the original transition date (td) file
      #that phenor will be expecting when it runs pr_fit to estimate phenophases below.
      #the list (td.1) that we created does not do this automatically
      rinfo = data.frame("site" = unique(td$site), "veg_type" = unique(td$veg_type), "roi_id" = 
                           unique(td$roi_id), "direction" = "falling")
      #make falling.2 by combining rinfo with falling.1
      falling.2 = cbind(rinfo, falling.1)
      #force names of falling.2 to be the same as rising.2, otherwise they cannot be combined
      names(falling.2) = names(rising.2)
    }else
      
      #if data.frame is empty, make falling.2 NA
      falling.2 = NA   
    
    #combine rising.2 and falling.2 into a new data.frame only if both are present. 
    #otherwise, make the new data.frame out of either rising.2 (falling.2 is NA)
    #or falling.2 (rising.2 is NA)
    #note that R throws a warning saying there are multiple conditions and it will only use the first one
    #there are not multiple conditions. R is just being cranky. Don't worry about the warnings
    
    #determine if rising.2 = NA
    if(is.na(rising.2) == T) {
      #if yes, td.3 = falling.2
      td.3 = falling.2}
    #determine if falling.2 = NA
    else if(is.na(falling.2) == T) {
      #if yes, td.3 = rising.2
      td.3 = rising.2
    } else  
      #otherwise td.2 is the combination of rising.2 and falling.2
      td.3 = rbind(rising.2, falling.2)
    #write td.3 to the temporary directory, using the name name as the transition date file
    #this will overwrite the original file, which is what we want to have the user-defined
    #transition dates
    write.csv(td.3, files.mxf[i+1])
    
    ###################################################################################
    #format phenocam data to include DAYMET variables and user-defined transition dates
    ###################################################################################
    
    for(j in 1:length(files.temp)){
      
      # formate phenocam
      pc_form_rising_10 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 10,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_10, "try-error"))
        
        next
      
      pc_form_rising_50 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 50,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_50, "try-error"))
        
        next
      
      pc_form_rising_90 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 90,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_90, "try-error"))
        
        next
      
      pc_form_falling_10 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 10,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_10, "try-error"))
        
        next
      
      pc_form_falling_50 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 50,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_50, "try-error"))
        
        next
      
      pc_form_falling_80 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 90,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_80, "try-error"))
        
        next
      
      #append formatted phenocam files to lists
      
      mxf_rising_10 = append(mxf_rising_10, pc_form_rising_10)
      mxf_rising_50 = append(mxf_rising_50, pc_form_rising_50)
      mxf_rising_90 = append(mxf_rising_90, pc_form_rising_90)
      mxf_falling_10 = append(mxf_falling_10, pc_form_falling_10)
      mxf_falling_50 = append(mxf_falling_50, pc_form_falling_50)
      mxf_falling_80 = append(mxf_falling_80, pc_form_falling_80)
      
      
    }
    
    #set working directory and remove files to make space for the next pair
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    file.remove(files.mxf[i])
    file.remove(files.mxf[i+1])
    
    #end loop
  }else
    next
}

#remove null elements from list (in case stations had null values for both rising and falling
#during the process of assigning phenophases)
#mxf_ris_10 = mxf_rising_10[-which(lapply(mxf_rising_10,is.null) == T)]
#mxf_ris_50 = mxf_rising_50[-which(lapply(mxf_rising_50,is.null) == T)]
#mxf_ris_90 = mxf_rising_90[-which(lapply(mxf_rising_90,is.null) == T)]

#mxf_fal_10 = mxf_falling_10[-which(lapply(mxf_falling_10,is.null) == T)]
#mxf_fal_50 = mxf_falling_50[-which(lapply(mxf_falling_50,is.null) == T)]
#mxf_fal_90 = mxf_falling_90[-which(lapply(mxf_falling_90,is.null) == T)]

#determine the length of the list in order to remove every other item
#loop writes each station twice and needs to be halved
#not super elegant programming solution, but it works, and pr_fit will not
#operate with duplicates!

len_lis_mxf = seq(2, length(mxf_rising_10), by = 2)

#remove every other item (indexed as len_list_mxf) by making it NULL
#rename each list prior to doing this so as not to lose list completely (paranoid move, but 
#safer choice than overwriting each list and having to rerun again)

mxf_ris_10 = mxf_rising_10
mxf_ris_10[len_lis_mxf] = NULL

mxf_ris_50 = mxf_rising_50
mxf_ris_50[len_lis_mxf] = NULL

mxf_ris_90 = mxf_rising_90
mxf_ris_90[len_lis_mxf] = NULL

mxf_fal_10 = mxf_falling_10
mxf_fal_10[len_lis_mxf] = NULL

mxf_fal_50 = mxf_falling_50
mxf_fal_50[len_lis_mxf] = NULL

mxf_fal_80 = mxf_falling_80
mxf_fal_80[len_lis_mxf] = NULL

#######################################################################
#use the pr_fit function to select the parameter set for each threshold
#(10, 50, 90) and time period (rising and falling)
#######################################################################

mxf_par_ris_10 = pr_fit(data = mxf_ris_10,
                        model = "TT",
                        method = "GenSA",
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

mxf_par_ris_50 = pr_fit(data = mxf_ris_50,
                        model = "TT",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

mxf_par_ris_90 = pr_fit(data = mxf_ris_90,
                        model = "TT",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

mxf_par_fal_10 = pr_fit(data = mxf_fal_10,
                        model = "CDD",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

mxf_par_fal_50 = pr_fit(data = mxf_fal_50,
                        model = "CDD",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

mxf_par_fal_80 = pr_fit(data = mxf_fal_80,
                        model = "CDD",
                        method = "GenSA",
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

#combine phenology parameters into a single list
mxf_pars = append(list(mxf_par_ris_10), list(mxf_par_ris_50))
mxf_pars = append(mxf_pars, list(mxf_par_ris_90))
mxf_pars = append(mxf_pars, list(mxf_par_fal_10))
mxf_pars = append(mxf_pars, list(mxf_par_fal_50))
mxf_pars = append(mxf_pars, list(mxf_par_fal_80))

#setwd and export list
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling")
list.save(mxf_pars, "mxf_pars.rdata")

##########################
##########################
#Urban (urb)
##########################
##########################

#set working directory
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/URB")

#read files names in working directory
files.urb = list.files(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/URB", full.names = F)

#split file names to ID the ones that have raw greenness data
#this will allow for the sequential selecting of stations in the loop below
#R will throw a warning because not all items in data frame have same length (no worries about this)
urb_spl.1 = data.frame(do.call(rbind, str_split(files.urb, "_")))

#create empty containers to hold the lists that will be generated for each threshold (10, 50, 90) and
#each time period (rising, falling)
urb_rising_10 = list()
urb_rising_50 = list()
urb_rising_90 = list()

urb_falling_10 = list()
urb_falling_50 = list()
urb_falling_80 = list()

#loop through pairs of files from each station in URB. Each pair contains the raw gcc data and 
#predetermined transition dates. The code will pull these pairs and place them in a temporary directory.
#Once in the directory, the pr_fm_phenocam will try to format the data to prepare it for downstream modeling.
#the "try" argument will keep the loop going even if the function throws and error, allowing it
#to move to the next pair of files in URB. Some of the code will alter the file with transition dates
#so instead of the predetermined ones of 10, 25, and 50 percentile gcc, it does 10, 50, and 90 gcc for the "rising"
#or spring phenophases, and 10, 50, and 80 for the "falling" or fall phenophases

#begin outer loop 

#########################################
#select files and place them in directory
#########################################

for(i in 1:nrow(urb_spl.1)){
  #choose files from deciduous broadleaf veg_type and 
  if(urb_spl.1$X6[i] == "DB"){
    #set working directory to get these files
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/URB")
    
    #copy the first two files from the URB directory to the temporary directory. This pair
    #contains the raw gcc data and predetermined transition dates. 
    file.copy(files.urb[i], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    file.copy(files.urb[i+1], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    #record the file names in the temporary directory
    files.temp = list.files(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir/", full.names = F)
    
    ###########################################
    #format files to alter the transition dates
    ###########################################
    
    #set working directory to temporary directory. IMPORTANT to do this because many of the functions in 
    #phenor use the "." syntax to refer to the home directory or have this embedded in the background code
    #to enable the processing of files. Things will end up in the wrong place and the loop will
    #likely break of the working directory is not specified properly!
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    # read in the transition date file (the second file in the temporary directory)
    td = read.table(files.urb[i+1],
                    header = TRUE,
                    sep = ",")
    
    # smooth data using an optimization approach; forcing (force = TRUE) allows for the skipping of the 
    # procedure if the smoothed data are already available
    smooth_ts(as.character(files.urb[i]),
              out_dir = ".",
              force = TRUE)
    
    # generate phenological transition dates (10, 50, and 90 percentile thresholds). Results in a list
    # with two items: a list of gcc thresholds with the transition dates for the rising and falling limbs
    # of the phenophases. For this, we are interested in the rising limb
    td.1 = phenophases(as.character(files.urb[i]),
                       internal = TRUE,
                       lower_thresh = 0.1,
                       middle_thresh = 0.5,
                       upper_thresh = 0.9)
    
    #take the rising part of the list, convert to a data.frame, and remove row names
    rising.1 = do.call(rbind.data.frame, td.1[1])
    rownames(rising.1) = NULL
    
    #NOTE: not all of the phenocams have data for both rising and falling phenophases. The next few lines of
    #code enable the delineation of rising as "NA" if no data exist
    
    #determine if data.frame is empty
    if(is.data.frame(rising.1) && nrow(rising.1)!=0){
      #if false
      #convert columns 2 to 10 (estimated phenophases with 95% CIs) to dates
      rising.1[2:10] = lapply(rising.1[2:10], as.Date, origin = "1970-01-01")
      #create columns containing site, veg_type, roi, and phenophase direction info
      #these are the columns that appear in the original transition date (td) file
      #that phenor will be expecting when it runs pr_fit to estimate phenophases below.
      #the list (td.1) that we created does not do this automatically
      rinfo = data.frame("site" = unique(td$site), "veg_type" = unique(td$veg_type), "roi_id" = 
                           unique(td$roi_id), "direction" = "rising")
      #make rising.2 by combining rinfo with rising.1
      rising.2 = cbind(rinfo, rising.1)}else
        
        #if data.frame is empty, make rising.2 NA
        rising.2 = NA   
    
    # generate phenological transition dates (10, 50, and 80 percentile thresholds). Results in a list
    # with two items: a list of gcc thresholds with the transition dates for the rising and falling limbs
    # of the phenophases. For this, we are interested in the falling limb
    td.2 = phenophases(as.character(files.urb[i]),
                       internal = TRUE,
                       lower_thresh = 0.1,
                       middle_thresh = 0.5,
                       upper_thresh = 0.8)
    
    #take the falling part of the list, convert to a data.frame, and remove row names
    falling.1 = do.call(rbind.data.frame, td.2[2])
    rownames(falling.1) = NULL
    
    #NOTE: not all of the phenocams have data for both falling and falling phenophases. The next few lines of
    #code enable the delineation of falling as "NA" if no data exist
    
    #determine if data.frame is empty
    if(is.data.frame(falling.1) && nrow(falling.1)!=0){
      #if false
      #convert columns 2 to 10 (estimated phenophases with 95% CIs) to dates
      falling.1[2:10] = lapply(falling.1[2:10], as.Date, origin = "1970-01-01")
      #create columns containing site, veg_type, roi, and phenophase direction info
      #these are the columns that appear in the original transition date (td) file
      #that phenor will be expecting when it runs pr_fit to estimate phenophases below.
      #the list (td.1) that we created does not do this automatically
      rinfo = data.frame("site" = unique(td$site), "veg_type" = unique(td$veg_type), "roi_id" = 
                           unique(td$roi_id), "direction" = "falling")
      #make falling.2 by combining rinfo with falling.1
      falling.2 = cbind(rinfo, falling.1)
      #force names of falling.2 to be the same as rising.2, otherwise they cannot be combined
      names(falling.2) = names(rising.2)
    }else
      
      #if data.frame is empty, make falling.2 NA
      falling.2 = NA   
    
    #combine rising.2 and falling.2 into a new data.frame only if both are present. 
    #otherwise, make the new data.frame out of either rising.2 (falling.2 is NA)
    #or falling.2 (rising.2 is NA)
    #note that R throws a warning saying there are multiple conditions and it will only use the first one
    #there are not multiple conditions. R is just being cranky. Don't worry about the warnings
    
    #determine if rising.2 = NA
    if(is.na(rising.2) == T) {
      #if yes, td.3 = falling.2
      td.3 = falling.2}
    #determine if falling.2 = NA
    else if(is.na(falling.2) == T) {
      #if yes, td.3 = rising.2
      td.3 = rising.2
    } else  
      #otherwise td.2 is the combination of rising.2 and falling.2
      td.3 = rbind(rising.2, falling.2)
    #write td.3 to the temporary directory, using the name name as the transition date file
    #this will overwrite the original file, which is what we want to have the user-defined
    #transition dates
    write.csv(td.3, files.urb[i+1])
    
    ###################################################################################
    #format phenocam data to include DAYMET variables and user-defined transition dates
    ###################################################################################
    
    for(j in 1:length(files.temp)){
      
      # formate phenocam
      pc_form_rising_10 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 10,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_10, "try-error"))
        
        next
      
      pc_form_rising_50 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 50,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_50, "try-error"))
        
        next
      
      pc_form_rising_90 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 90,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_90, "try-error"))
        
        next
      
      pc_form_falling_10 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 10,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_10, "try-error"))
        
        next
      
      pc_form_falling_50 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 50,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_50, "try-error"))
        
        next
      
      pc_form_falling_80 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 90,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_80, "try-error"))
        
        next
      
      #append formatted phenocam files to lists
      
      urb_rising_10 = append(urb_rising_10, pc_form_rising_10)
      urb_rising_50 = append(urb_rising_50, pc_form_rising_50)
      urb_rising_90 = append(urb_rising_90, pc_form_rising_90)
      urb_falling_10 = append(urb_falling_10, pc_form_falling_10)
      urb_falling_50 = append(urb_falling_50, pc_form_falling_50)
      urb_falling_80 = append(urb_falling_80, pc_form_falling_80)
      
      
    }
    
    #set working directory and remove files to make space for the next pair
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    file.remove(files.urb[i])
    file.remove(files.urb[i+1])
    
    #end loop
  }else
    next
}

#remove null elements from list (in case stations had null values for both rising and falling
#during the process of assigning phenophases)
#urb_ris_10 = urb_rising_10[-which(lapply(urb_rising_10,is.null) == T)]
#urb_ris_50 = urb_rising_50[-which(lapply(urb_rising_50,is.null) == T)]
#urb_ris_90 = urb_rising_90[-which(lapply(urb_rising_90,is.null) == T)]

#urb_fal_10 = urb_falling_10[-which(lapply(urb_falling_10,is.null) == T)]
#urb_fal_50 = urb_falling_50[-which(lapply(urb_falling_50,is.null) == T)]
#urb_fal_90 = urb_falling_90[-which(lapply(urb_falling_90,is.null) == T)]

#determine the length of the list in order to remove every other item
#loop writes each station twice and needs to be halved
#not super elegant programming solution, but it works, and pr_fit will not
#operate with duplicates!

len_lis_urb = seq(2, length(urb_rising_10), by = 2)

#remove every other item (indexed as len_list_urb) by making it NULL
#rename each list prior to doing this so as not to lose list completely (paranoid move, but 
#safer choice than overwriting each list and having to rerun again)

urb_ris_10 = urb_rising_10
urb_ris_10[len_lis_urb] = NULL

urb_ris_50 = urb_rising_50
urb_ris_50[len_lis_urb] = NULL

urb_ris_90 = urb_rising_90
urb_ris_90[len_lis_urb] = NULL

urb_fal_10 = urb_falling_10
urb_fal_10[len_lis_urb] = NULL

urb_fal_50 = urb_falling_50
urb_fal_50[len_lis_urb] = NULL

urb_fal_80 = urb_falling_80
urb_fal_80[len_lis_urb] = NULL

#######################################################################
#use the pr_fit function to select the parameter set for each threshold
#(10, 50, 90) and time period (rising and falling)
#######################################################################

urb_par_ris_10 = pr_fit(data = urb_ris_10,
                        model = "TT",
                        method = "GenSA",
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

urb_par_ris_50 = pr_fit(data = urb_ris_50,
                        model = "TT",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

urb_par_ris_90 = pr_fit(data = urb_ris_90,
                        model = "TT",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

urb_par_fal_10 = pr_fit(data = urb_fal_10,
                        model = "CDD",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

urb_par_fal_50 = pr_fit(data = urb_fal_50,
                        model = "CDD",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

urb_par_fal_80 = pr_fit(data = urb_fal_80,
                        model = "CDD",
                        method = "GenSA",
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

#combine phenology parameters into a single list
urb_pars = append(list(urb_par_ris_10), list(urb_par_ris_50))
urb_pars = append(urb_pars, list(urb_par_ris_90))
urb_pars = append(urb_pars, list(urb_par_fal_10))
urb_pars = append(urb_pars, list(urb_par_fal_50))
urb_pars = append(urb_pars, list(urb_par_fal_80))

#setwd and export list
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling")
list.save(urb_pars, "urb_pars.rdata")

##########################
##########################
#Agricultural (agr)
##########################
##########################

#set working directory
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/AGR")

#read files names in working directory
files.agr = list.files(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/AGR", full.names = F)

#split file names to ID the ones that have raw greenness data
#this will allow for the sequential selecting of stations in the loop below
#R will throw a warning because not all items in data frame have same length (no worries about this)
agr_spl.1 = data.frame(do.call(rbind, str_split(files.agr, "_")))

#create empty containers to hold the lists that will be generated for each threshold (10, 50, 90) and
#each time period (rising, falling)
agr_rising_10 = list()
agr_rising_50 = list()
agr_rising_90 = list()

agr_falling_10 = list()
agr_falling_50 = list()
agr_falling_80 = list()

#loop through pairs of files from each station in AGR. Each pair contains the raw gcc data and 
#predetermined transition dates. The code will pull these pairs and place them in a temporary directory.
#Once in the directory, the pr_fm_phenocam will try to format the data to prepare it for downstream modeling.
#the "try" argument will keep the loop going even if the function throws and error, allowing it
#to move to the next pair of files in AGR. Some of the code will alter the file with transition dates
#so instead of the predetermined ones of 10, 25, and 50 percentile gcc, it does 10, 50, and 90 gcc for the "rising"
#or spring phenophases, and 10, 50, and 80 for the "falling" or fall phenophases

#begin outer loop 

#########################################
#select files and place them in directory
#########################################

for(i in 1:nrow(agr_spl.1)){
  #choose files from deciduous broadleaf veg_type and 
  if(agr_spl.1$X6[i] == "AG"){
    #set working directory to get these files
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/AGR")
    
    #copy the first two files from the AGR directory to the temporary directory. This pair
    #contains the raw gcc data and predetermined transition dates. 
    file.copy(files.agr[i], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    file.copy(files.agr[i+1], "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    #record the file names in the temporary directory
    files.temp = list.files(path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir/", full.names = F)
    
    ###########################################
    #format files to alter the transition dates
    ###########################################
    
    #set working directory to temporary directory. IMPORTANT to do this because many of the functions in 
    #phenor use the "." syntax to refer to the home directory or have this embedded in the background code
    #to enable the processing of files. Things will end up in the wrong place and the loop will
    #likely break of the working directory is not specified properly!
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    # read in the transition date file (the second file in the temporary directory)
    td = read.table(files.agr[i+1],
                    header = TRUE,
                    sep = ",")
    
    # smooth data using an optimization approach; forcing (force = TRUE) allows for the skipping of the 
    # procedure if the smoothed data are already available
    smooth_ts(as.character(files.agr[i]),
              out_dir = ".",
              force = TRUE)
    
    # generate phenological transition dates (10, 50, and 90 percentile thresholds). Results in a list
    # with two items: a list of gcc thresholds with the transition dates for the rising and falling limbs
    # of the phenophases. For this, we are interested in the rising limb
    td.1 = phenophases(as.character(files.agr[i]),
                       internal = TRUE,
                       lower_thresh = 0.1,
                       middle_thresh = 0.5,
                       upper_thresh = 0.9)
    
    #take the rising part of the list, convert to a data.frame, and remove row names
    rising.1 = do.call(rbind.data.frame, td.1[1])
    rownames(rising.1) = NULL
    
    #NOTE: not all of the phenocams have data for both rising and falling phenophases. The next few lines of
    #code enable the delineation of rising as "NA" if no data exist
    
    #determine if data.frame is empty
    if(is.data.frame(rising.1) && nrow(rising.1)!=0){
      #if false
      #convert columns 2 to 10 (estimated phenophases with 95% CIs) to dates
      rising.1[2:10] = lapply(rising.1[2:10], as.Date, origin = "1970-01-01")
      #create columns containing site, veg_type, roi, and phenophase direction info
      #these are the columns that appear in the original transition date (td) file
      #that phenor will be expecting when it runs pr_fit to estimate phenophases below.
      #the list (td.1) that we created does not do this automatically
      rinfo = data.frame("site" = unique(td$site), "veg_type" = unique(td$veg_type), "roi_id" = 
                           unique(td$roi_id), "direction" = "rising")
      #make rising.2 by combining rinfo with rising.1
      rising.2 = cbind(rinfo, rising.1)}else
        
        #if data.frame is empty, make rising.2 NA
        rising.2 = NA   
    
    # generate phenological transition dates (10, 50, and 80 percentile thresholds). Results in a list
    # with two items: a list of gcc thresholds with the transition dates for the rising and falling limbs
    # of the phenophases. For this, we are interested in the falling limb
    td.2 = phenophases(as.character(files.agr[i]),
                       internal = TRUE,
                       lower_thresh = 0.1,
                       middle_thresh = 0.5,
                       upper_thresh = 0.8)
    
    #take the falling part of the list, convert to a data.frame, and remove row names
    falling.1 = do.call(rbind.data.frame, td.2[2])
    rownames(falling.1) = NULL
    
    #NOTE: not all of the phenocams have data for both falling and falling phenophases. The next few lines of
    #code enable the delineation of falling as "NA" if no data exist
    
    #determine if data.frame is empty
    if(is.data.frame(falling.1) && nrow(falling.1)!=0){
      #if false
      #convert columns 2 to 10 (estimated phenophases with 95% CIs) to dates
      falling.1[2:10] = lapply(falling.1[2:10], as.Date, origin = "1970-01-01")
      #create columns containing site, veg_type, roi, and phenophase direction info
      #these are the columns that appear in the original transition date (td) file
      #that phenor will be expecting when it runs pr_fit to estimate phenophases below.
      #the list (td.1) that we created does not do this automatically
      rinfo = data.frame("site" = unique(td$site), "veg_type" = unique(td$veg_type), "roi_id" = 
                           unique(td$roi_id), "direction" = "falling")
      #make falling.2 by combining rinfo with falling.1
      falling.2 = cbind(rinfo, falling.1)
      #force names of falling.2 to be the same as rising.2, otherwise they cannot be combined
      names(falling.2) = names(rising.2)
    }else
      
      #if data.frame is empty, make falling.2 NA
      falling.2 = NA   
    
    #combine rising.2 and falling.2 into a new data.frame only if both are present. 
    #otherwise, make the new data.frame out of either rising.2 (falling.2 is NA)
    #or falling.2 (rising.2 is NA)
    #note that R throws a warning saying there are multiple conditions and it will only use the first one
    #there are not multiple conditions. R is just being cranky. Don't worry about the warnings
    
    #determine if rising.2 = NA
    if(is.na(rising.2) == T) {
      #if yes, td.3 = falling.2
      td.3 = falling.2}
    #determine if falling.2 = NA
    else if(is.na(falling.2) == T) {
      #if yes, td.3 = rising.2
      td.3 = rising.2
    } else  
      #otherwise td.2 is the combination of rising.2 and falling.2
      td.3 = rbind(rising.2, falling.2)
    #write td.3 to the temporary directory, using the name name as the transition date file
    #this will overwrite the original file, which is what we want to have the user-defined
    #transition dates
    write.csv(td.3, files.agr[i+1])
    
    ###################################################################################
    #format phenocam data to include DAYMET variables and user-defined transition dates
    ###################################################################################
    
    for(j in 1:length(files.temp)){
      
      # formate phenocam
      pc_form_rising_10 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 10,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_10, "try-error"))
        
        next
      
      pc_form_rising_50 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 50,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_50, "try-error"))
        
        next
      
      pc_form_rising_90 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "rising",
        gcc_value = "gcc_90",
        threshold = 90,
        offset = 264,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_rising_90, "try-error"))
        
        next
      
      pc_form_falling_10 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 10,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_10, "try-error"))
        
        next
      
      pc_form_falling_50 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 50,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_50, "try-error"))
        
        next
      
      pc_form_falling_80 = try(pr_fm_phenocam(
        path = "/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir",
        direction = "falling",
        gcc_value = "gcc_90",
        threshold = 90,
        offset = 365,
        internal = TRUE), silent = T)
      
      if(inherits(pc_form_falling_80, "try-error"))
        
        next
      
      #append formatted phenocam files to lists
      
      agr_rising_10 = append(agr_rising_10, pc_form_rising_10)
      agr_rising_50 = append(agr_rising_50, pc_form_rising_50)
      agr_rising_90 = append(agr_rising_90, pc_form_rising_90)
      agr_falling_10 = append(agr_falling_10, pc_form_falling_10)
      agr_falling_50 = append(agr_falling_50, pc_form_falling_50)
      agr_falling_80 = append(agr_falling_80, pc_form_falling_80)
      
      
    }
    
    #set working directory and remove files to make space for the next pair
    setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/Phenocam_Data/Temp_Dir")
    
    file.remove(files.agr[i])
    file.remove(files.agr[i+1])
    
    #end loop
  }else
    next
}

#remove null elements from list (in case stations had null values for both rising and falling
#during the process of assigning phenophases)
#agr_ris_10 = agr_rising_10[-which(lapply(agr_rising_10,is.null) == T)]
#agr_ris_50 = agr_rising_50[-which(lapply(agr_rising_50,is.null) == T)]
#agr_ris_90 = agr_rising_90[-which(lapply(agr_rising_90,is.null) == T)]

#agr_fal_10 = agr_falling_10[-which(lapply(agr_falling_10,is.null) == T)]
#agr_fal_50 = agr_falling_50[-which(lapply(agr_falling_50,is.null) == T)]
#agr_fal_90 = agr_falling_90[-which(lapply(agr_falling_90,is.null) == T)]

#determine the length of the list in order to remove every other item
#loop writes each station twice and needs to be halved
#not super elegant programming solution, but it works, and pr_fit will not
#operate with duplicates!

len_lis_agr = seq(2, length(agr_rising_10), by = 2)

#remove every other item (indexed as len_list_agr) by making it NULL
#rename each list prior to doing this so as not to lose list completely (paranoid move, but 
#safer choice than overwriting each list and having to rerun again)

agr_ris_10 = agr_rising_10
agr_ris_10[len_lis_agr] = NULL

agr_ris_50 = agr_rising_50
agr_ris_50[len_lis_agr] = NULL

agr_ris_90 = agr_rising_90
agr_ris_90[len_lis_agr] = NULL

agr_fal_10 = agr_falling_10
agr_fal_10[len_lis_agr] = NULL

agr_fal_50 = agr_falling_50
agr_fal_50[len_lis_agr] = NULL

agr_fal_80 = agr_falling_80
agr_fal_80[len_lis_agr] = NULL

#######################################################################
#use the pr_fit function to select the parameter set for each threshold
#(10, 50, 90) and time period (rising and falling)
#######################################################################

agr_par_ris_10 = pr_fit(data = agr_ris_10,
                        model = "TT",
                        method = "GenSA",
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

agr_par_ris_50 = pr_fit(data = agr_ris_50,
                        model = "TT",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

agr_par_ris_90 = pr_fit(data = agr_ris_90,
                        model = "TT",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

agr_par_fal_10 = pr_fit(data = agr_fal_10,
                        model = "CDD",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

agr_par_fal_50 = pr_fit(data = agr_fal_50,
                        model = "CDD",
                        method = "GenSA",
                        control = list(max.call = 2000),
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

agr_par_fal_80 = pr_fit(data = agr_fal_80,
                        model = "CDD",
                        method = "GenSA",
                        par_ranges = system.file("extdata",
                                                 "parameter_ranges.csv", package = "phenor", mustWork = TRUE),
                        plot = TRUE)

#combine phenology parameters into a single list
agr_pars = append(list(agr_par_ris_10), list(agr_par_ris_50))
agr_pars = append(agr_pars, list(agr_par_ris_90))
agr_pars = append(agr_pars, list(agr_par_fal_10))
agr_pars = append(agr_pars, list(agr_par_fal_50))
agr_pars = append(agr_pars, list(agr_par_fal_80))

#setwd and export list
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling")
list.save(agr_pars, "agr_pars.rdata")

##################################
#Export phenology model parameters
##################################

#combine parameters from all lists into single dataframe

#create empy dataframe
all_pars = data.frame(matrix(vector(), nrow = 4, ncol = 19))
names(all_pars) = c("lc", "T0_R10", "Tb_R10", "Fc_R10", 
                    "T0_R50", "Tb_R50", "Fc_R50", 
                    "T0_R90", "Tb_R90", "Fc_R90", 
                    "T0_F10", "Tb_F10", "Fc_F10",
                    "T0_F50", "Tb_F50", "Fc_F50",
                     "T0_F80", "Tb_F80", "Fc_F80")
              
#assing land cover. 4 = deciduous broadleaf, 5 = mixed forest, 10 = urban, 14 = mixed cropland / 
#natural vegetation mosaic
all_pars$lc = c(4, 5, 10, 14)

#add values for T0, Tbase, and Fcrit for each of the land cover types and phenological phases

#for dbf
all_pars[1, 2:4] = dbf_pars[[1]]$par
all_pars[1, 5:7] = dbf_pars[[2]]$par
all_pars[1, 8:10] = dbf_pars[[3]]$par
all_pars[1, 11:13] = dbf_pars[[4]]$par
all_pars[1, 14:16] = dbf_pars[[5]]$par
all_pars[1, 17:19] = dbf_pars[[6]]$par

#for mxf
all_pars[2, 2:4] = mxf_pars[[1]]$par
all_pars[2, 5:7] = mxf_pars[[2]]$par
all_pars[2, 8:10] = mxf_pars[[3]]$par
all_pars[2, 11:13] = mxf_pars[[4]]$par
all_pars[2, 14:16] = mxf_pars[[5]]$par
all_pars[2, 17:19] = mxf_pars[[6]]$par

#for urb
all_pars[3, 2:4] = urb_pars[[1]]$par
all_pars[3, 5:7] = urb_pars[[2]]$par
all_pars[3, 8:10] = urb_pars[[3]]$par
all_pars[3, 11:13] = urb_pars[[4]]$par
all_pars[3, 14:16] = urb_pars[[5]]$par
all_pars[3, 17:19] = urb_pars[[6]]$par

#for agr
all_pars[4, 2:4] = agr_pars[[1]]$par
all_pars[4, 5:7] = agr_pars[[2]]$par
all_pars[4, 8:10] = agr_pars[[3]]$par
all_pars[4, 11:13] = agr_pars[[4]]$par
all_pars[4, 14:16] = agr_pars[[5]]$par
all_pars[4, 17:19] = agr_pars[[6]]$par

#merge all_pars with station.lc to have a list of phenology parameters for each station
station_pars = merge(station.lc, all_pars, by.x = "lc", by.y = "lc")

#write table
write.table(station_pars, "station_pars.csv", col.names = T, row.names = F, sep = ",")
