## Combining the PREDICTS database with the fire atlas database
rm(list = ls())

## Load the required packages:
library(ggplot2)
library(dplyr)
library(tidyr) 
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lwgeom)
library(predictsFunctions) 
library(StatisticalModels)
library(yarg)
library(roquefort)
library(rmarkdown)


## Load in PREDICTS data
diversity <- readRDS("~/Imperial/PROJECT/PREDICTS data/database.rds")
sites <- readRDS("~/Imperial/PROJECT/PREDICTS data/sites.rds")

## Subset PREDICTS so that is 2003 onwards 
sites_post_2003 <- subset(sites, sites$Sample_start_earliest > "2003-01-01")
## Remove NAs
sites_post_2003_NA <- subset(sites_post_2003, (! is.na(Longitude) | ! is.na(Latitude)))

## Make the PREDICTS data a spatial feature
sites.sf.point <- st_as_sf(x = sites_post_2003_NA, 
                           coords = c("Longitude", "Latitude"),
                           crs = "+proj=longlat +datum=WGS84")

## Creating the function which takes each perimeter spatial feature, selects the attributes needed
#converts the start date of the fire to a date, intersects the two spatial features, changes it to a 
#dataframe and then subsets the data where the fire has occured after sampling. 
fires <- function(data.frame){
  data.frame %>% 
    select(fire_ID, size, start_date, start_DOY, duration) -> perimetershort
  perimetershort$start_date <- as.Date(perimetershort$start_date)
  out<- st_intersection(sites.sf.point, perimetershort)
  out <- as.data.frame(out)
  out[!out$Sample_start_earliest < out$start_date, ]
}

setwd("~/Imperial/PROJECT/FIRE ATLAS DATA/CMS_Global_Fire_Atlas_1642/CMS_Global_Fire_Atlas_1642/data")
## 2013
perimeter2013 <- st_read('Global_fire_atlas_V1_perimeter_2013.shp')
PREDICTS13 <- fires(perimeter2013)

## 2012
perimeter2012 <- st_read('Global_fire_atlas_V1_perimeter_2012.shp')
PREDICTS12 <- fires(perimeter2012)

## 2011
perimeter2011 <- st_read('Global_fire_atlas_V1_perimeter_2011.shp')
PREDICTS11 <- fires(perimeter2011)

## Save the dateframe we're getting as csv's
setwd("~/Imperial/PROJECT/Testing the fire atlas")
write.csv(PREDICTS11, 'PREDICTS_fires_2011.csv')
write.csv(PREDICTS12, 'PREDICTS_fires_2012.csv')
write.csv(PREDICTS13, 'PREDICTS_fires_2013.csv')

## Free up space in r memory
rm(perimeter2013)
rm(perimeter2012)
rm(perimeter2011)

## carrying on...
setwd("~/Imperial/PROJECT/FIRE ATLAS DATA/CMS_Global_Fire_Atlas_1642/CMS_Global_Fire_Atlas_1642/data")

## 2010
perimeter2010 <- st_read('Global_fire_atlas_V1_perimeter_2010.shp')
PREDICTS10 <- fires(perimeter2010)

## 2009
perimeter2009 <- st_read('Global_fire_atlas_V1_perimeter_2009.shp')
PREDICTS09 <- fires(perimeter2009)

## 2008
perimeter2008 <- st_read('Global_fire_atlas_V1_perimeter_2008.shp')
PREDICTS08 <- fires(perimeter2008)

rm(perimeter2010)
rm(perimeter2009)
rm(perimeter2008)

setwd("~/Imperial/PROJECT/Testing the fire atlas")
write.csv(PREDICTS10, 'PREDICTS_fires_2010.csv')
write.csv(PREDICTS09, 'PREDICTS_fires_2009.csv')
write.csv(PREDICTS08, 'PREDICTS_fires_2008.csv')

setwd("~/Imperial/PROJECT/FIRE ATLAS DATA/CMS_Global_Fire_Atlas_1642/CMS_Global_Fire_Atlas_1642/data")

## 2007
## error in 2007 - see end of script

## 2006
perimeter2006 <- st_read('Global_fire_atlas_V1_perimeter_2006.shp')
PREDICTS06 <- fires(perimeter2006)

## 2005
perimeter2005 <- st_read('Global_fire_atlas_V1_perimeter_2005.shp')
PREDICTS05 <- fires(perimeter2005)

setwd("~/Imperial/PROJECT/Testing the fire atlas")
write.csv(PREDICTS06, 'PREDICTS_fires_2006.csv')
write.csv(PREDICTS05, 'PREDICTS_fires_2005.csv')

rm(perimeter2006)
rm(perimeter2005)

setwd("~/Imperial/PROJECT/FIRE ATLAS DATA/CMS_Global_Fire_Atlas_1642/CMS_Global_Fire_Atlas_1642/data")

## 2004
perimeter2004 <- st_read('Global_fire_atlas_V1_perimeter_2004.shp')
PREDICTS04 <- fires(perimeter2004)

setwd("~/Imperial/PROJECT/Testing the fire atlas")
write.csv(PREDICTS04, 'PREDICTS_fires_2004.csv')

rm(perimeter2004)

## 2003
## also error in 2003 

setwd("~/Imperial/PROJECT/FIRE ATLAS DATA/CMS_Global_Fire_Atlas_1642/CMS_Global_Fire_Atlas_1642/data")

## Fixing issues with 2007 and 2003
fires <- function(data.frame){
  data.frame %>% 
    select(fire_ID, size, start_date, start_DOY, duration) -> perimetershort
  perimetershort$start_date <- as.Date(perimetershort$start_date)
  out<- st_intersection(st_make_valid(perimetershort), sites.sf.point)
  out <- as.data.frame(out)
  out[!out$Sample_start_earliest < out$start_date, ]
}
## Need to change the intersection line to make sure that all geometries within the perimeter data 
# are valid. Increases the time it takes to run. 

perimeter2007 <- st_read('Global_fire_atlas_V1_perimeter_2007.shp')
PREDICTS07 <- fires(perimeter2007)
perimeter2003 <- st_read('Global_fire_atlas_V1_perimeter_2003.shp')
PREDICTS03 <- fires(perimeter2003)

## Merge all the dataframes together...
Fire_in_PREDICTS <- rbind(PREDICTS03, PREDICTS04, PREDICTS05, PREDICTS06, PREDICTS07, PREDICTS08, 
                          PREDICTS09, PREDICTS10, PREDICTS11, PREDICTS12, PREDICTS13)
# Can now add in column of time since fire (in days)
Fire_in_PREDICTS$days_since_fire <- Fire_in_PREDICTS$Sample_start_earliest - Fire_in_PREDICTS$start_date
## convert geometry to lat and lon so i can plot the sites on a map
Fire_in_PREDICTS_latlon <- Fire_in_PREDICTS %>% extract(geometry, c('lon', 'lat'), '\\((.*), (.*)\\)', convert = TRUE)
## Adding a years column to the dataframe..
Fire_in_PREDICTS_latlon[, "year"] <- format(Fire_in_PREDICTS_latlon[,"start_date"], "%Y")

## save the new dataframe
setwd("~/Imperial/PROJECT/Testing the fire atlas")
write.csv(Fire_in_PREDICTS_latlon, 'PREDICTS_Fires.csv') ## Dataframe containing all details of 
#PREDICTS sites which have had a fire. 

###################################################################################################
## Editting of the data and dataframe

## Cutting down the number of columns and looking at how many unique sites we have and adding on a 
#frequency column.

## Working off the saved df "PREDICTS_Fires.csv"
PREDICTS_Fires <- read.csv("PREDICTS_Fires.csv")
## Checking for fire frequencies in PREDICTS 
tab <- table(PREDICTS_Fires$SSBS)
tab_df <- as.data.frame(tab) ## We have 905 sites which have atleast 1 fire
## Renamed the variable column to the same name as in the firePREDICTS data table. 
names(tab_df)[names(tab_df) == "Var1"] <- "SSBS"
## Merge the two dataframes together
new_df <- merge(x=PREDICTS_Fires, y=tab_df, by = "SSBS")
## Making the df shorter
new_df %>% 
  dplyr::select(SSBS, Reference, Sample_start_earliest, Sample_end_latest, Predominant_land_use, Use_intensity,
         Country, UN_subregion, Biome, lon, lat, fire_ID, size, start_date, start_DOY, duration, 
         days_since_fire, year, Freq) -> Fire_Predicts_Short_withFreq

## Now have a cut down table of all the sites and the key information!
setwd("~/Imperial/PROJECT/Testing the fire atlas")
write.csv(Fire_Predicts_Short_withFreq, 'PREDICTS_and_Fire_cut.csv')
## Now working off the saved df "PREDICTS_Fire_cut"
PREDICTS_Fires_cut <- read.csv("PREDICTS_and_Fire_cut.csv")

##Fire histories 

## Looking at which sites have had a fire within the last 3 years and have a three year period 
#available before sampling (so from 2006 onwards)
## we can just exclude based on the days e.g. 3 years is 1095 days - therefore we can exlude all 
#rows where the time since fire column has a number greater than 1095. 
PREDICTS_fire_3yrs <- PREDICTS_Fires[!PREDICTS_Fires$days_since_fire > 1095, ]
PREDICTS_fire_3yrs$Sample_start_earliest <- as.Date(PREDICTS_fire_3yrs$Sample_start_earliest)
PREDICTS_fire_3yrs<- subset(PREDICTS_fire_3yrs, PREDICTS_fire_3yrs$Sample_start_earliest > "2006-01-01")

PREDICTS_fire_3yrs%>% 
  dplyr::select(SSBS, Reference, Sample_start_earliest, Sample_end_latest, Predominant_land_use, Use_intensity,
         Country, UN_subregion, Biome, lon, lat, fire_ID, size, start_date, start_DOY, duration, 
         days_since_fire, year) -> PREDICTS_fire_3yrs

## Extracting the number of sites which have had a fire. 
tab <- table(PREDICTS_fire_3yrs$SSBS)
tab_df <- as.data.frame(tab)
View(tab)
tab_df_no0 <- tab_df[!(tab_df$Freq== 0),]
## 594 sites with at least 1 fire
## Renamed the variable column to the same name as in the firePREDICTS data table. 
names(tab_df_no0)[names(tab_df_no0) == "Var1"] <- "SSBS" 

PREDICTS_fire_3yrs_freq <- merge(x=PREDICTS_fire_3yrs, y=tab_df_no0, by = "SSBS")
write.csv(PREDICTS_fire_3yrs_freq, "PREDICTS_fire_3yrs_freq.csv")
## this dataframe show all PREDICTS sites which have had a fire sampled after 2006 with details on 
#individual fires + freq

sites_post_2003_NA %>% 
  dplyr::select(SSBS, Reference, Sample_start_earliest, Sample_end_latest, Predominant_land_use, Use_intensity,
         Country, UN_subregion, Biome, Longitude, Latitude) -> sites_short

PREDICTS_sites_3yrhist <- merge(x=sites_short, y=tab_df_no0, by = "SSBS")
## THis data frame shows all PREDICTS sites that ahve a three year fire history with atleast 1 fire. 
write.csv(PREDICTS_sites_3yrhist, "PREDICTS_sites_3yrhist.csv")
## Data frame contains the 594 sites (and details of these sites) which have had at least 1 fire in 
#a three year period before sampling. 

## We now want to extract the details of the most recent fire that has effected each of the 594 sites
#which have fire histories. 

## Summaries of PREDICTS sites and their most recent fires:
## Can generate dataframes of sites which have had one fire and sites which have had more than one
## First need to create blank dataframes that have the same column headings
## Working off the PREDICTS_Fires_cut dataframe 
PREDICTS_fire_freq <- PREDICTS_Fires_cut
PREDICTS_fire_freq <- PREDICTS_fire_freq[0,] ## Dataframe for more than one fire
PREDICTS_fire_once <- PREDICTS_Fires_cut
PREDICTS_fire_once <- PREDICTS_fire_once[0,] ## Dataframe for one fire

for (i in 1:nrow(PREDICTS_Fires_cut)) {
  if (PREDICTS_Fires_cut$Freq[i] != 1) {
    PREDICTS_fire_freq = rbind(PREDICTS_fire_freq, PREDICTS_Fires_cut[i,])}
  else {
    PREDICTS_fire_once = rbind(PREDICTS_fire_once, PREDICTS_Fires_cut[i,])}
}

##Now trying to filter out all the rows with duplicated SSBS and keep only the most recent. 
PREDICTS_fire_freq_copy <- PREDICTS_fire_freq ## creating a copy to do the editing on
## order the df by date
PREDICTS_fire_freq_copy <- PREDICTS_fire_freq_copy[order(PREDICTS_fire_freq_copy$start_date),]
## Use the rev() and duplicated() functions to remove all replicates earlier in the dataframe.
PREDICTS_fire_freq_copy <- PREDICTS_fire_freq_copy[!rev(duplicated(rev(PREDICTS_fire_freq_copy$SSBS))),]
##Merging with the single fire sites
PREDICTS_most_recent_fire <- rbind(PREDICTS_fire_freq_copy, PREDICTS_fire_once)
write.csv(PREDICTS_most_recent_fire, "PREDICTS_most_recent_fire.csv") ## Saved as a csv
PREDICTS_most_recent_fire <- read.csv("PREDICTS_most_recent_fire.csv")
## Now have a dataframe containing all the sites which have had fires and the details of the most 
#recent fire at that site. 

## We now want to generate the fire histories with the information on the most recent fires at that 
#site

Fire_history <- merge(PREDICTS_sites_3yrhist, PREDICTS_most_recent_fire
                      [,c("SSBS", "fire_ID", "size", "start_date", "duration", "days_since_fire", "year")],
                      by="SSBS", all.x=TRUE)
## Creating a column of years since fire (categories not exact numbers)
Fire_history$years_since_fire <- ifelse(Fire_history$days_since_fire<365, 1,
                                        ifelse(Fire_history$days_since_fire>365& Fire_history$days_since_fire<730, 2,
                                               ifelse(Fire_history$days_since_fire>730 & Fire_history$days_since_fire<1095, 3,
                                                      NA)))
Fire_history <- subset(Fire_history, select = -c(X))
## We now have a table which gives the detail of sites with fires in the last 3years, how many 
#in the last 3 years and details of the most recent fire in the last 3 years
write.csv(Fire_history, "Fire_history.csv")

## Merging with PREDICTS for analysis:
## First need to make sure the PREDICTS database is ready and has the correct columns including 
#measures of biodiversity

## Making Measurement and Effort_corrected_measurement the same thing and Sampling_effort and 
#Rescaled_sampling_effort the same. 
diversity <- mutate(diversity, 
                    Measurement = Effort_corrected_measurement,
                    Sampling_effort = Rescaled_sampling_effort)

## Calculate diversity level metrics include rescaled abundance. 
sites <- diversity %>%
  # add Diversity_metric_is_valid column
  mutate(Diversity_metric_is_valid = TRUE) %>%
  # calculate SiteMetrics  
  predictsFunctions::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Biome", "Country")) %>%
  # calculate the total abundance within each study
  group_by(SS) %>%
  mutate(MaxAbundance = ifelse(Diversity_metric_type == "Abundance",
                               max(Total_abundance),
                               NA)) %>%
  ungroup() %>%
  # now calculate the rescaled abundance (abundance divided by the maximum within each study)
  mutate(RescaledAbundance = ifelse(Diversity_metric_type == "Abundance",
                                    Total_abundance/MaxAbundance,
                                    NA))

## Tidying up the predominant land use category 
sites <- sites %>%
  mutate(Predominant_land_use = recode_factor(Predominant_land_use, 
                                         "Primary forest" = "Primary vegetation", 
                                         "Primary non-forest" = "Primary vegetation"),
    Predominant_land_use = na_if(Predominant_land_use, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use = na_if(Predominant_land_use, "Cannot decide"),
    Use_intensity = na_if(Use_intensity, "Cannot decide"),
    Predominant_land_use = factor(Predominant_land_use),
    Predominant_land_use = relevel(Predominant_land_use, ref = "Primary vegetation"),
    Use_intensity = factor(Use_intensity),
    Use_intensity = relevel(Use_intensity, ref = "Minimal use")
  )

## Subset so that all sites are after 2006 (so that we can have 3yrs of fire data)
sites_post_2006 <- subset(sites, sites$Sample_start_earliest > "2006-01-01")
## Removing rows that have NAs in latitude and longitude as these won't have fire data
sites_post_2006_NA <- subset(sites_post_2006, (! is.na(Longitude) | ! is.na(Latitude)))

## Merge with fire data
PREDICTS_fire_hist <- merge(sites_post_2006_NA, Fire_history
                            [,c("SSBS", "Freq", "fire_ID", "size", "start_date", "duration", "days_since_fire",
                                "year", "years_since_fire")], by="SSBS", all.x=TRUE)

## Change freq from NA to 0
PREDICTS_fire_hist$Freq[is.na(PREDICTS_fire_hist$Freq)] <- 0
## Adding a column of fire presence
PREDICTS_fire_hist$fire_presence <- ifelse(PREDICTS_fire_hist$Freq==0, "No",
                                           ifelse(PREDICTS_fire_hist$Freq!=0, "Yes",
                                                  NA))
PREDICTS_fire_hist$fire_presence <- as.factor(PREDICTS_fire_hist$fire_presence)

## Creating a map of all sites (abundance and richness)
map_data <- drop_na(PREDICTS_fire_hist, Predominant_land_use, Use_intensity, fire_presence)
table(map_data$fire_presence)
## 9579 no fire, 542 fire (total = 10121)

world <- ne_countries(scale = "medium", returnclass = "sf")
model_data_fire <- subset(map_data, map_data$fire_presence == "Yes")
model_data_nofire <- subset(map_data, map_data$fire_presence == "No")
ggplot(data = world) +
  geom_sf()+
  geom_point(data = model_data_nofire, aes(x = Longitude, y = Latitude, fill = "#56B4E9"), 
             size = 4, shape = 21, alpha = 0.5) +
  geom_point(data = model_data_fire, aes(x = Longitude, y = Latitude, fill = "#D55E00"), 
             size = 4, shape = 21) +
  scale_fill_identity(name = NULL,
                      breaks = c("#D55E00", "#56B4E9"),
                      labels = c("Sites with fires", "Sites without fires"), 
                      guide = "legend")+
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE)+
  xlab('Longitude')+
  ylab('Latitude')+
  theme_bw()+
  theme(legend.position = c(0.1, 0.2))

## Looking at the savanna:
sites_savanna <- subset(PREDICTS_fire_hist, PREDICTS_fire_hist$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands")
table(sites_savanna$fire_presence)
