##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script estimates phenology using thermal time models for the rising and falling limbs of deciduous 
#vegetation. Model parameters were generated from phenocam imagery

####################################################################################
#Initial set up                                                                   
####################################################################################

#call libraries

library(dplyr)
library(zoo)


#set wd and read in files
setwd("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Weather_Station_Data/Daily_Data/")

#call file names in Daily_modSWE subdirectory
files.modSWE = list.files(path = "Daily_modSWE/",
                          full.names = T)

#read table with station, land cover, and phenology modeling parameters
station_pars = read.table("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Phenology-Modeling/station_pars.csv",
                          head = T, sep = ",")

#split the ID names in station_pars to get just the station ID
station.1 = data.frame(do.call(rbind, str_split(station_pars$ID, ":")))
names(station.1) = c("dat", "station")

#add to data.frame
station_pars$station = station.1$station

####################################################################################
#Model phenology                                                                  
####################################################################################

#loop through stations in files.modSWE that match stations in station_pars

st.PHEN = c()

for(i in 1:nrow(station_pars)){
  
  st.PHEN[i] = paste("Daily_modSWE//GHCND_", station_pars$station[i], "_Daily_modSWE.csv", sep = "")
  station = as.character(station_pars$station[i])
  
  #read in weather data
  phen = read.table(st.PHEN[i], header = T, sep = ",")
  
  #remove duplicate dates
  phen = phen[!duplicated(phen$DATE), ]

  #create new column for modeling the rising limb of the phenological year, which begins September 1 and ends August 30
  phen$pyear = ifelse(phen$doy < 244, phen$year - 1, phen$year)
  
  #remove the first phenological year to avoid calculating spurious greenness values
  phendat = phen[phen$pyear > min(phen$pyear), ]
  
  #create new column in which September 1 is doy 1 for modeling the rising limb
  phendat$pdoy = ifelse(phendat$doy < 244, phendat$doy + 122, phendat$doy-243)
  
  #gap fill missing temperature and precipitation data because the model will fail without complete data
  #temperature data can be interpolated
  phendat$TMAXgap = ifelse(is.na(phendat$TMAX) == T, na.approx(phendat$TMAX, rule = 2, na.rm = F), phendat$TMAX)
  phendat$TMINgap = ifelse(is.na(phendat$TMIN) == T, na.approx(phendat$TMIN, rule = 2, na.rm = F), phendat$TMIN)
  
  #calcuate daily average temperature
  phendat$TMEAN = (phendat$TMAXgap + phendat$TMINgap) / 2
  
  #convert into units that the themal time phenology model expects (degrees C)
  phendat$TMEANc = phendat$TMEAN / 10
 
  #calculate degree days
  
  #rising (heating degree days)
  phendat$R10_hdd = ifelse(phendat$pdoy >= station_pars$T0_R10[i] & phendat$TMEANc > station_pars$Tb_R10[i], 
                           phendat$TMEANc - station_pars$Tb_R10[i], 0)
  
  phendat$R50_hdd = ifelse(phendat$pdoy >= station_pars$T0_R50[i] & phendat$TMEANc > station_pars$Tb_R50[i], 
                           phendat$TMEANc - station_pars$Tb_R50[i], 0)
  
  phendat$R90_hdd = ifelse(phendat$pdoy >= station_pars$T0_R90[i] & phendat$TMEANc > station_pars$Tb_R90[i], 
                           phendat$TMEANc - station_pars$Tb_R90[i], 0)
  
  #falling (cooling degree days)
  phendat$F10_cdd = ifelse(phendat$doy >= station_pars$T0_F10[i] & phendat$TMEANc < station_pars$Tb_F10[i], 
                           phendat$TMEANc - station_pars$Tb_F10[i], 0)
  
  phendat$F50_cdd = ifelse(phendat$doy >= station_pars$T0_F50[i] & phendat$TMEANc < station_pars$Tb_F50[i], 
                           phendat$TMEANc - station_pars$Tb_F50[i], 0)
  
  phendat$F80_cdd = ifelse(phendat$doy >= station_pars$T0_F80[i] & phendat$TMEANc < station_pars$Tb_F80[i], 
                           phendat$TMEANc - station_pars$Tb_F80[i], 0)
  
  
  #calculate cumulative sums of hdds and cdds within years
  phendat$R10_csum = ave(phendat$R10_hdd, phendat$pyear,FUN=cumsum)
  phendat$R50_csum = ave(phendat$R50_hdd, phendat$pyear,FUN=cumsum)
  phendat$R90_csum = ave(phendat$R90_hdd, phendat$pyear,FUN=cumsum)
  phendat$F10_csum = ave(phendat$F10_cdd, phendat$year,FUN=cumsum)
  phendat$F50_csum = ave(phendat$F50_cdd, phendat$year,FUN=cumsum)
  phendat$F80_csum = ave(phendat$F80_cdd, phendat$year,FUN=cumsum)
  
  #flag days when the cumulative sum of hdd > F_crit
  phendat$R10_crit = ifelse(phendat$R10_csum > station_pars$Fc_R10[i], phendat$doy, NA)
  phendat$R50_crit = ifelse(phendat$R50_csum > station_pars$Fc_R50[i], phendat$doy, NA)
  phendat$R90_crit = ifelse(phendat$R90_csum > station_pars$Fc_R90[i], phendat$doy, NA)
  phendat$F10_crit = ifelse(phendat$F10_csum < station_pars$Fc_F10[i], phendat$doy, NA)
  phendat$F50_crit = ifelse(phendat$F50_csum < station_pars$Fc_F50[i], phendat$doy, NA)
  phendat$F80_crit = ifelse(phendat$F80_csum < station_pars$Fc_F80[i], phendat$doy, NA)
  
  #select the first days of the year when the cumsum reaches the F_crit threshold
  #then subset each data frame to include year and doy that threshold is reached
  #add column identifying the threshold (needed when merging with phendat below and 
  #estimating %gcc on days between thresholds)
  
  R10 = phendat %>%
    group_by(year) %>%
    filter(R10_crit == min(R10_crit, na.rm = T))
  R10_sub = R10[ , c("year", "doy")]
  R10_sub$thresh = "R10"
  
  R50 = phendat %>%
    group_by(year) %>%
    filter(R50_crit == min(R50_crit, na.rm = T))
  R50_sub = R50[ , c("year", "doy")]
  R50_sub$thresh = "R50"
  
  R90 = phendat %>%
    group_by(year) %>%
    filter(R90_crit == min(R90_crit, na.rm = T))
  R90_sub = R90[ , c("year", "doy")]
  R90_sub$thresh = "R90"
  
  F10 = phendat %>%
    group_by(year) %>%
    filter(F10_crit == min(F10_crit, na.rm = T))
  F10_sub = F10[ , c("year", "doy")]
  F10_sub$thresh = "F10"
  
  F50 = phendat %>%
    group_by(year) %>%
    filter(F50_crit == min(F50_crit, na.rm = T))
  F50_sub = F50[ , c("year", "doy")]
  F50_sub$thresh = "F50"
  
  F80 = phendat %>%
    group_by(year) %>%
    filter(F80_crit == min(F80_crit, na.rm = T))
  F80_sub = F80[ , c("year", "doy")]
  F80_sub$thresh = "F80"
  
  #combine into single dataframe with year and doy that each phenological threshold was reached
  allphen =  data.frame(rbind(R10_sub, R50_sub, R90_sub, F80_sub, F50_sub, F10_sub))
  
  #omit September 1 for the first full phenological year in the dataset
  allphen = allphen[allphen$year >= min(allphen$year) & allphen$doy != 244, ]
  
  
  #merge with phendat
  #make ID columns in each data.frame 
  allphen$ID = paste(allphen$year, allphen$doy, sep = " ")
  phen$ID = paste(phen$year, phen$doy, sep = " ")
  
  phenest = merge(phen, allphen, by.x = "ID", by.y = "ID", all.x = T, all.y = T)
  
  #order by date
  phenest = phenest[order(phenest$year.x, phenest$doy.x), ]
  
  #make new column that contains the phenophase (rising: 10, 50, 90, or falling: 80, 50, 10) as represented
  #by percent greenness. 
  phenest$pergreen = ifelse(phenest$thresh == "R10", 10, 
                            ifelse(phenest$thresh == "R50", 50,
                                   ifelse(phenest$thresh == "R90", 90,
                                          ifelse(phenest$thresh == "F80", 80,
                                                 ifelse(phenest$thresh == "F50", 50,
                                                        ifelse(phenest$thresh == "F10", 10, NA))))))
  
  #fill in pergreen between phenophases using na.approx
  phenest$pergreen = ifelse(is.na(phenest$pergreen) == T, na.approx(phenest$pergreen, na.rm = F), phenest$pergreen)
  
  #omit columns no longer needed
  phenest_fin=phenest[c("ID", "year.x", "month", "ym", "doy.x", "DATE", "station",
                      "PRCPfin", "SNOWfin", "SNWDfin", "TMAXfin", "TMINfin", "mod_ID", 
                      "mSWE", "mSNOW", "mLIQ", "thresh", "pergreen")]
  names(phenest_fin) = c("ID", "year", "month", "ym", "doy", "DATE", "station",
                         "PRCPfin", "SNOWfin", "SNWDfin", "TMAXfin", "TMINfin", "mod_ID", 
                         "mSWE", "mSNOW", "mLIQ", "thresh", "pergreen")

  #write table to folder
  #create output table name
  out_name = paste("/nfs/WinterWeatherWhiplash-data/Contosta_Projects/Weather_Station_Data/Daily_Data/Daily_Phenology/", 
                   sub(":", "_", station), 
                   "_", 
                   "Daily_Phenology", 
                   ".csv", 
                   sep="")
  
  write.csv(phenest_fin, out_name)
  
}

#additional code for plotting phenology across years

#omit duplicate years from the beginning and the end of the dataframe to determaxe when to start / end model subsets

start.rows <- !duplicated(phenest_fin$year)
end.rows <- !duplicated(phenest_fin$year, fromLast = TRUE)

#extract original row indexes associated with the start and end of an SiteYr

sr.ind <- seq_along(phenest_fin$year)[start.rows]
er.ind <- seq_along(phenest_fin$year)[end.rows]

plot(phenest_fin$doy, phenest_fin$pergreen, ylim = c(0, 100),type = "n", xlab = "Day of Year", ylab = "Estimated Greenness (%)")

cl <- rainbow(61)

#cl = brewer.pal(61, "Blues")

for (i in seq(1,length(sr.ind))) {
  
  lines(phenest_fin$doy[sr.ind[i]:er.ind[i]], phenest_fin$pergreen[sr.ind[i]:er.ind[i]], col = cl[i], type = "l")
}



