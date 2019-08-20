##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script estimates model parameters used to simulate phenology from daily temperature
#using the phenor and phenocamr packages (https://khufkens.github.io/phenor/) developed
#by K. Hufkens et al. (2018)

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

##########################
#Import phenocam metadata 
##########################

#File is in .json format and has to be converted to data.frame
#Be sure to have file path to main project directory in place

#set working directory
setwd("/research-home/acontosta/noaadata")

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
study.region = readOGR(dsn = "Spatial_Domain/", layer = "NA_Eco_Level2_Clip_SSC")

#set the extent of the whole study region
ex = extent(study.region)

#select stations that fit within extent (ex)

pcsub = pcmeta[pcmeta$lat >= ex[3] &
                 pcmeta$lat <= ex[4] &
                 pcmeta$lon >= ex[1] &
                 pcmeta$lon <= ex[2] , ]

#subset station list into vegetation type for ROI (region of interest) at which phenocam points
#"EN" = evergreen; "DB" = deciduous broadleaf", "AG" = agricultural
pcenf = pcsub[pcsub$veg_type == "EN", ]
pcdbf = pcsub[pcsub$veg_type == "DB", ]
pcagr = pcsub[pcsub$veg_type == "AG", ]

#extract station names from subsets
dbf = as.character(pcdbf$site)
enf = as.character(pcenf$site)
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

setwd("/research-home/acontosta/noaadata/Phenocam_Data/DBF")

for(i in 1:length(dbf)){
  try(
  phenocamr::download_phenocam(
    frequency = 3,
    veg_type = "DB",
    roi_id = 1000,
    site = dbf[i],
    phenophase = TRUE,
    out_dir = "."
  )  
  )
  if(inherits(phenocamr::download_phenocam, "try-error"))next
}

##########################
#Evergreen (enf)
##########################

setwd("/research-home/acontosta/noaadata/Phenocam_Data/ENF")

for(i in 1:length(enf)){
  try(
    phenocamr::download_phenocam(
      frequency = 3,
      veg_type = "EN",
      roi_id = 1000,
      site = enf[i],
      phenophase = TRUE,
      out_dir = "."
    )  
  )
  if(inherits(phenocamr::download_phenocam, "try-error"))next
}


##########################
#Agriculture (agr)
##########################

setwd("/research-home/acontosta/noaadata/Phenocam_Data/AGR")

for(i in 1:length(agr)){
  try(
    phenocamr::download_phenocam(
      frequency = 3,
      veg_type = "AG",
      roi_id = 1000,
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
setwd("/research-home/acontosta/noaadata/Phenocam_Data/DBF")

#read files names in working directory
files.dbf = list.files(path = "/research-home/acontosta/noaadata/Phenocam_Data/DBF/", full.names = F)

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
  #choose files from deciduous broadleaf veg_type
  if(dbf_spl.1$X6[i] == "DB"){
  #set working directory to get these files
  setwd("/research-home/acontosta/noaadata/Phenocam_Data/DBF")
 
  #copy the first two files from the DBF directory to the temporary directory. This pair
  #contains the raw gcc data and predetermined transition dates. 
  file.copy(files.dbf[i], "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir")
  file.copy(files.dbf[i+1], "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir")

  #record the file names in the temporary directory
  files.temp = list.files(path = "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir/", full.names = F)
  
  ###########################################
  #format files to alter the transition dates
  ###########################################
  
  #set working directory to temporary directory. IMPORTANT to do this because many of the functions in 
  #phenor use the "." syntax to refer to the home directory or have this embedded in the background code
  #to enable the processing of files. Things will end up in the wrong place and the loop will
  #likely break of the working directory is not specified properly!
  setwd("/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir")
  
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
    path = "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir",
    direction = "rising",
    gcc_value = "gcc_90",
    threshold = 10,
    offset = 264,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_rising_10, "try-error"))
  
    next

  pc_form_rising_50 = try(pr_fm_phenocam(
    path = "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir",
    direction = "rising",
    gcc_value = "gcc_90",
    threshold = 50,
    offset = 264,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_rising_50, "try-error"))
    
    next
  
  pc_form_rising_90 = try(pr_fm_phenocam(
    path = "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir",
    direction = "rising",
    gcc_value = "gcc_90",
    threshold = 90,
    offset = 264,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_rising_90, "try-error"))
    
    next
  
  pc_form_falling_10 = try(pr_fm_phenocam(
    path = "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir",
    direction = "falling",
    gcc_value = "gcc_90",
    threshold = 10,
    offset = 365,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_falling_10, "try-error"))
    
    next
  
  pc_form_falling_50 = try(pr_fm_phenocam(
    path = "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir",
    direction = "falling",
    gcc_value = "gcc_90",
    threshold = 50,
    offset = 365,
    internal = TRUE), silent = T)
  
  if(inherits(pc_form_falling_50, "try-error"))
    
    next
  
  pc_form_falling_80 = try(pr_fm_phenocam(
    path = "/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir",
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
  setwd("/research-home/acontosta/noaadata/Phenocam_Data/Temp_Dir")
  
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

#extract parameters and write to data.frame


