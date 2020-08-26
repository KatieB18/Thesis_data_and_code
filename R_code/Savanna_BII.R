## Savanna BII models:
rm(list = ls())

library(dplyr) 
library(tidyr) 
library(magrittr) 
library(lme4) 
library(car)
library(raster) 
library(geosphere) 
library(foreach) 
library(doParallel)
library(gower)
library(sjPlot)
library(ggplot2)


## PREDICTS data
diversity <- readRDS("~/Imperial/PROJECT/PREDICTS data/database.rds")
## Fire atlas data
fire_hist <- read.csv("~/Imperial/PROJECT/Combined data/Fire_history.csv")

#################################################################################################
## Data preparation:
diversity <- mutate(diversity, 
                    Measurement = Effort_corrected_measurement,
                    Sampling_effort = Rescaled_sampling_effort)
#Subset the data to the fire atlas range
diversity_post_2006 <- subset(diversity, diversity$Sample_start_earliest > "2006-01-01")
## Removing rows that have NAs in latitude and longitude as these won't have fire data
diversity_post_2006_NA <- subset(diversity_post_2006, (! is.na(Longitude) | ! is.na(Latitude)))

## Merge with fire data
diversity_fire_hist <- merge(diversity_post_2006_NA, fire_hist
                             [,c("SSBS", "Freq", "fire_ID", "size", "start_date", "duration", "days_since_fire",
                                 "year", "years_since_fire")], by="SSBS", all.x=TRUE)
## Change freq from NA to 0
diversity_fire_hist$Freq[is.na(diversity_fire_hist$Freq)] <- 0
## Adding a column of fire presence
diversity_fire_hist$fire_presence <- ifelse(diversity_fire_hist$Freq==0, "No",
                                            ifelse(diversity_fire_hist$Freq!=0, "Yes",
                                                   NA))
diversity_fire_hist$fire_presence <- as.factor(diversity_fire_hist$fire_presence)

## Subset to just the savanna
savanna_diversity <- subset(diversity_fire_hist, diversity_fire_hist$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands")

### Exploring data
table(savanna_diversity$Predominant_land_use, savanna_diversity$Use_intensity)

## Initially we're going to use just primary minimal sites as the baseline
savanna_diversity <- savanna_diversity %>%
  # make a level of Primary minimal. Everything else gets the coarse land use
  mutate(
    LandUse = ifelse(Predominant_land_use == "Primary vegetation" & Use_intensity == "Minimal use",
                     "Primary minimal",
                     paste(Predominant_land_use)),
    # collapse the secondary vegetation classes together
    LandUse = ifelse(grepl("secondary", tolower(LandUse)),
                     "Secondary vegetation",
                     paste(LandUse)),
    # change cannot decide into NA
    LandUse = ifelse(Predominant_land_use == "Cannot decide",
                     NA, 
                     paste(LandUse)),
    # relevel the factor so that Primary minimal is the first level (so that it is the intercept term in models)
    LandUse = factor(LandUse),
    LandUse = relevel(LandUse, ref = "Primary minimal")
  )
table(savanna_diversity$LandUse)

table(savanna_diversity$LandUse, savanna_diversity$fire_presence) 

## As we're interested in fires I am going to remove urban and plantation forests due to lack of 
#data (as in models of abundance). Also going to drop the rows containing NAs

savanna_diversity<- drop_na(savanna_diversity, LandUse)
savanna_diversity<-savanna_diversity[!(savanna_diversity$LandUse == "Plantation forest"),]
savanna_diversity<-savanna_diversity[!(savanna_diversity$LandUse == "Urban"),]
savanna_diversity$LandUse <- as.factor(as.character(savanna_diversity$LandUse))
savanna_diversity <- savanna_diversity %>%
  mutate( LandUse = factor(LandUse),
          LandUse = relevel(LandUse, ref = "Primary minimal"))
levels(savanna_diversity$LandUse)

## Total abundance
abundance_data <- savanna_diversity %>%
  # pull out just the abundance measures
  filter(Diversity_metric_type == "Abundance") %>%
  # group by SSBS (each unique value corresponds to a unique site)
  group_by(SSBS) %>%
  # now add up all the abundance measurements within each site
  mutate(TotalAbundance = sum(Effort_corrected_measurement)) %>%
  # ungroup
  ungroup() %>%
  # pull out unique sites
  distinct(SSBS, .keep_all = TRUE) %>%
  # now group by Study ID
  group_by(SS) %>%
  # pull out the maximum abundance for each study
  mutate(MaxAbundance = max(TotalAbundance)) %>%
  # ungroup
  ungroup() %>%
  # now rescale total abundance, so that within each study, abundance varies from 0 to 1.
  mutate(RescaledAbundance = TotalAbundance/MaxAbundance)

## Compositional similarity (using the asymmetric Jaccard Index)

## data set up
cd_data_input <- savanna_diversity %>%
  # drop any rows with unknown LandUse
  filter(!is.na(LandUse)) %>%
  # pull out only the abundance data
  filter(Diversity_metric_type == "Abundance") %>%
  # group by Study
  group_by(SS) %>%
  # calculate the number of unique sampling efforts within that study
  mutate(n_sample_effort = n_distinct(Sampling_effort)) %>%
  # calculate the number of unique species sampled in that study
  mutate(n_species = n_distinct(Taxon_name_entered)) %>%
  # check if there are any Primary minimal sites in the dataset
  mutate(n_primin_records = sum(LandUse == "Primary minimal")) %>%
  # ungroup
  ungroup() %>%
  # now keep only the studies with one unique sampling effort
  filter(n_sample_effort == 1) %>%
  # and keep only studies with more than one species 
  # as these studies clearly aren't looking at assemblage-level diversity
  filter(n_species > 1) %>%
  # and keep only studies with at least some Primary minimal data
  filter(n_primin_records > 0) %>%
  # drop empty factor levels
  droplevels()

table(cd_data_input$LandUse, cd_data_input$fire_presence)

## Now set up a function to calculate compositional similarity between a single pair of sites 
# in a study.
getJacAbSym <- function(s1, s2, data){
  s1species <- data %>%
    filter(SSBS == s1) %>%
    filter(Measurement > 0) %>%
    distinct(Taxon_name_entered) %>%
    pull
  s2abundance_s1species <- data %>%
    filter(SSBS == s2) %>%
    filter(Taxon_name_entered %in% s1species) %>%
    pull(Measurement) %>%
    sum()
  s2_sum <- data %>%
    filter(SSBS == s2) %>%
    pull(Measurement) %>%
   sor <- s2abundance_s1species / s2_sum
   return(sor)
}

# vector of studies to loop over
studies <- distinct(cd_data_input, SS) %>%
  pull()

## Now need to loop over each element in studies and calculate the compositional similarity 
# between pairs of sites, the geographical distance between pairs and land uses for each pair.

registerDoParallel(cores = 2)
cd_data <- foreach(s = studies, 
                   .combine = rbind,
                   .packages = c("dplyr", "magrittr", "geosphere")) %dopar% {
                     data_ss <- filter(cd_data_input, SS == s)
                     site_data <- data_ss %>%
                       dplyr::select(SSBS, LandUse, Sample_start_earliest, Longitude, Latitude, Freq,
                                     fire_presence, size, start_date, duration, days_since_fire) %>%
                       distinct(SSBS, .keep_all = TRUE)
                     baseline_sites <- site_data %>%
                       filter(LandUse == "Primary minimal") %>%
                       pull(SSBS)
                     site_list <- site_data %>%
                       pull(SSBS)
                     site_comparisons <- expand.grid(baseline_sites, site_list) %>%
                       rename(s1 = Var1, s2 = Var2) %>%
                       filter(s1 != s2)
                     sor <- apply(site_comparisons, 1, function(y) getJacAbSym(data = data_ss, s1 = y['s1'], s2 = y['s2']))
                     s1LatLong <- as.matrix(data_ss[match(site_comparisons$s1, data_ss$SSBS), c('Longitude','Latitude')])
                     s2LatLong <- as.matrix(data_ss[match(site_comparisons$s2, data_ss$SSBS), c('Longitude','Latitude')])
                     dist <- distHaversine(s1LatLong, s2LatLong)
                     Contrast <- paste(site_data$LandUse[match(site_comparisons$s1, site_data$SSBS)],
                                       site_data$LandUse[match(site_comparisons$s2, site_data$SSBS)], 
                                       sep = "-")
                     study_results <- data.frame(site_comparisons,
                                                 sor,
                                                 dist,
                                                 Contrast,
                                                 SS = s,
                                                 stringsAsFactors = TRUE)
                     
                   }

# stop running things in parallel
registerDoSEQ()

## Calculating environmental distance
r <- raster::getData('worldclim', var='bio', res=10) # dowload climate variables
alt <- raster::getData("worldclim", var="alt", res=10) #download elevation
#want BIO5(max temp of warmest) BIO6(min temp of coldest) BIO13(precip wettest) and 
#BIO14(precip driest)
r <- r[[c(5,6,13,14)]]
names(r) <- c("max_temp", "min_temp", "max_precip", "min_precip")

## Extracting data for PREDICTS sites
cd_data_climate <- raster::extract(r, cd_data_input[,37:38])
cd_data_alt <- raster::extract(alt, cd_data_input[,37:38])
cd_data_envi <- cbind(cd_data_climate,cd_data_alt)
cd_data_input <- cbind(cd_data_input, cd_data_envi)

s1Enviro <- as.data.frame(cd_data_input[match(cd_data$s1, cd_data_input$SSBS), 
                                        c('max_temp', 'min_temp', 'max_precip', 'min_precip','cd_data_alt')])
s2Enviro <- as.data.frame(cd_data_input[match(cd_data$s2, cd_data_input$SSBS), 
                                        c('max_temp', 'min_temp', 'max_precip', 'min_precip', 'cd_data_alt')])
enviro_dist <- gower_dist(s1Enviro, s2Enviro)

cd_data <- cbind(cd_data, enviro_dist)

## Want to add on information on fire presence:
s1_fire <- as.data.frame(cd_data_input[match(cd_data$s1, cd_data_input$SSBS), 
                                       c('fire_presence')])
s2_fire <- as.data.frame(cd_data_input[match(cd_data$s2, cd_data_input$SSBS), 
                                       c('fire_presence')])
cd_fire <- cbind(s1_fire, s2_fire)
names(cd_fire) <- c("s1_fire", "s2_fire")

cd_data <- cbind(cd_data, cd_fire)

cd_data$fire_contrast <- paste(cd_data$s1_fire, cd_data$s2_fire, sep = "-")

## Editting the comp sim data for modelling
cd_data <- cd_data %>%
  separate(Contrast, c("s1_LandUse", "s2_LandUse"), sep = "-", remove = FALSE) %>%
  filter(s1_LandUse == "Primary minimal") %>%
  mutate(logitCS = logit(sor, adjust = 0.001, percents = FALSE)) %>%
  mutate(log10geo = log10(dist + 1)) %>%
  mutate(cuberoot_envi = enviro_dist^(1/3)) %>%
  mutate(Contrast = factor(Contrast), 
         Contrast = relevel(Contrast, ref = "Primary minimal-Primary minimal"))


## Setting up the abundance data and model:
## Check complete cases:
ab_model_data <- drop_na(abundance_data, RescaledAbundance, LandUse, fire_presence)
## Transform abundance
ab_model_data <- mutate(ab_model_data, 
                        logAbundance = log(RescaledAbundance + 1),
                        sqrtAbundance = sqrt(RescaledAbundance))
## Set up time since fire data:
str(ab_model_data$days_since_fire)
ab_model_data$days_since_fire <- as.numeric(ab_model_data$days_since_fire)
ab_model_data$days_since_fire[is.na(ab_model_data$days_since_fire)] <- 1277
ab_model_data$time_since_fire <- ab_model_data$days_since_fire/365

source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(ab_model_data[ , c("LandUse", "time_since_fire")])

## Model
ab_m <- lmer(sqrtAbundance ~ LandUse + time_since_fire + LandUse:time_since_fire +
               (1|SS) + (1|SSB), data = ab_model_data)
Anova(ab_m)
summary(ab_m) 

## Model simplification - remove the interaction
ab_m.1 <- update(ab_m,~.- LandUse:time_since_fire)
anova(ab_m, ab_m.1)
Anova(ab_m.1)
## Need to leave time since fire in the model. 
summary(ab_m.1)

theme_set(theme_classic())
plot_model(ab_m.1, type = "pred", terms = c("time_since_fire", "LandUse"),
           legend.title = c("Land use"),
           axis.title = c("Time since fire (years)",  "Square root rescaled abundance"))

## Compositional similarity:
str(cd_data_input$days_since_fire)
cd_data_input$days_since_fire <- as.numeric(cd_data_input$days_since_fire)
cd_data_input$days_since_fire[is.na(cd_data_input$days_since_fire)] <- 1277
cd_data_input$time_since_fire <- cd_data_input$days_since_fire/365

s1_time <- as.data.frame(cd_data_input[match(cd_data$s1, cd_data_input$SSBS),
                                       c('time_since_fire')])
s2_time <- as.data.frame(cd_data_input[match(cd_data$s2, cd_data_input$SSBS),
                                       c('time_since_fire')])
cd_time <- cbind(s1_time, s2_time)
names(cd_time) <- c("s1_time", "s2_time")

cd_data <- cbind(cd_data, cd_time)
cd_data$time_since_fire <- cd_data$s2_time

cd_data <- drop_na(cd_data, Contrast, log10geo, cuberoot_envi, time_since_fire)
corvif(cd_data[ , c("Contrast", "log10geo", "cuberoot_envi", "time_since_fire")])

## Check normality of logitCS
hist(cd_data$logitCS)
qqnorm(cd_data$logitCS)

table(cd_data$Contrast, cd_data$time_since_fire)

cd_m <- lmer(logitCS ~ Contrast + log10geo + cuberoot_envi + time_since_fire +
               (1|SS), data = cd_data)
Anova(cd_m)
summary(cd_m)

cd_m.1 <- update(cd_m,~.- cuberoot_envi)
anova(cd_m, cd_m.1)
Anova(cd_m.1)
summary(cd_m.1)

plot_model(cd_m.1, type = "pred", terms = c("time_since_fire"),
           axis.title = c("Time since fire (years)",  "Square root rescaled abundance"))

## Projecting models and BII

## Abundance:
new_data_ab <- expand.grid(time_since_fire = seq(0,3.5, by=0.1), 
                           LandUse = levels(ab_model_data$LandUse))

new_data_ab$ab_predicted <- predict(ab_m.1, newdata = new_data_ab, re.form=NA) ^ 2

## Compositional similarity:
inv_logit <- function(f, a){
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}

new_data_cd <- expand.grid(time_since_fire = seq(0,3.5, by=0.1), 
                           Contrast = levels(cd_data$Contrast),
                           log10geo = 0,
                           cuberoot_envi = 0)

new_data_cd$cd_predicted <-predict(cd_m.1, newdata = new_data_cd, re.form=NA) %>%
  inv_logit(a = 0.001)

## BII
bii <- new_data_ab$ab_predicted*new_data_cd$cd_predicted

bii_df <- new_data_ab
bii_df <- cbind(bii_df, new_data_cd$cd_predicted)
colnames(bii_df)
names(bii_df)[names(bii_df) == "new_data_cd$cd_predicted"] <- "cd_predicted"
## Want to divide the predicted values by a reference value - reference value wants to be
#primary minimal with little fire disturbance so 3.5 years? 
ref_value <- subset(bii_df, bii_df$LandUse == "Primary minimal" & bii_df$time_since_fire == "3.5")

bii_df$ab_predicted <- bii_df$ab_predicted/ref_value$ab_predicted
bii_df$cd_predicted <- bii_df$cd_predicted/ref_value$cd_predicted

bii_df$bii <- (bii_df$ab_predicted*bii_df$cd_predicted)*100

plot(bii_df$time_since_fire, bii_df$bii)

ggplot(data=bii_df, aes(x=time_since_fire, y=bii, group=LandUse)) +
  geom_line(aes(col=LandUse))+
  geom_point(aes(col=LandUse))+
  xlab(label = "Time since fire (years)") +
  ylab(label = "BII (%)")+
  theme_classic()

## The lack of studies is making the estimates quite unreliable. 

################################################################################################
## Trying a coarser baseline...
savanna_diversity <- subset(diversity_fire_hist, diversity_fire_hist$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands")

## Create a baseline of all primary plus mature secondary. 
levels(savanna_diversity$Predominant_land_use)

savanna_diversity_v2 <- savanna_diversity %>%
  mutate(Landuse2 = ifelse(Predominant_land_use == "Primary vegetation", "Baseline vegetation", 
                           paste(Predominant_land_use)),
         Landuse2 = ifelse(Predominant_land_use == "Mature secondary vegetation", "Baseline vegetation",
                           paste(Landuse2)),
         Landuse2 = ifelse(Predominant_land_use == "Cropland", "Agriculture",
                           paste(Landuse2)),
         Landuse2 = ifelse(Predominant_land_use == "Pasture", "Agriculture",
                           paste(Landuse2)),
         Landuse2 = ifelse(grepl("secondary", tolower(Landuse2)), "Secondary vegetation", 
                           paste(Landuse2)),
         Landuse2 = ifelse(Predominant_land_use == "Cannot decide",
                           NA,
                           paste(Landuse2)),
         Landuse2 = factor(Landuse2),
         Landuse2 = relevel(Landuse2, ref = "Baseline vegetation"))

levels(savanna_diversity_v2$Landuse2)

table(savanna_diversity_v2$Landuse2)
## remove urban and plantation

savanna_diversity_v2<- drop_na(savanna_diversity_v2, Landuse2)
savanna_diversity_v2<-savanna_diversity_v2[!(savanna_diversity_v2$Landuse2 == "Plantation forest"),]
savanna_diversity_v2<-savanna_diversity_v2[!(savanna_diversity_v2$Landuse2 == "Urban"),]
savanna_diversity_v2$Landuse2 <- as.factor(as.character(savanna_diversity_v2$Landuse2))
savanna_diversity_v2 <- savanna_diversity_v2 %>%
  mutate( Landuse2 = factor(Landuse2),
          Landuse2 = relevel(Landuse2, ref = "Baseline vegetation"))
levels(savanna_diversity_v2$Landuse2)

## Total abundance
abundance_data <- savanna_diversity_v2 %>%
  filter(Diversity_metric_type == "Abundance") %>%
  group_by(SSBS) %>%
  mutate(TotalAbundance = sum(Effort_corrected_measurement)) %>%
  ungroup() %>%
  distinct(SSBS, .keep_all = TRUE) %>%
  group_by(SS) %>%
  mutate(MaxAbundance = max(TotalAbundance)) %>%
  ungroup() %>%
  mutate(RescaledAbundance = TotalAbundance/MaxAbundance)


## data set up
cd_data_input <- savanna_diversity_v2 %>%
  filter(!is.na(Landuse2)) %>%
  filter(Diversity_metric_type == "Abundance") %>%
  group_by(SS) %>%
  mutate(n_sample_effort = n_distinct(Sampling_effort)) %>%
  mutate(n_species = n_distinct(Taxon_name_entered)) %>%
  mutate(n_primin_records = sum(Landuse2 == "Baseline vegetation")) %>%
  ungroup() %>%
  filter(n_sample_effort == 1) %>%
  filter(n_species > 1) %>%
  filter(n_primin_records > 0) %>%
  droplevels()

table(cd_data_input$Landuse2, cd_data_input$fire_presence) ## Good levels of data
table(cd_data_input$Landuse2)

## Compositional similarity function already set up. 

# Now need the dataset. First need a vector of studies to loop over 
studies <- distinct(cd_data_input, SS) %>%
  pull() ##Now have 12 studies

## Calculate compositional similarity 
registerDoParallel(cores = 2)
cd_data <- foreach(s = studies,
                   .combine = rbind,
                   .packages = c("dplyr", "magrittr", "geosphere")) %dopar% {
                     data_ss <- filter(cd_data_input, SS == s)
                     site_data <- data_ss %>%
                       dplyr::select(SSBS, Landuse2, Sample_start_earliest, Longitude, Latitude, Freq,
                                     fire_presence, size, start_date, duration, days_since_fire) %>%
                       distinct(SSBS, .keep_all = TRUE)
                     baseline_sites <- site_data %>%
                       filter(Landuse2 == "Baseline vegetation") %>% 
                       pull(SSBS)
                     site_list <- site_data %>%
                       pull(SSBS)
                     Site_comparisons <- expand.grid(baseline_sites, site_list) %>%
                       rename(s1 = Var1, s2 = Var2) %>%
                       filter(s1 != s2)
                     sor <- apply(Site_comparisons, 1, function(y) getJacAbSym(data = data_ss, s1 = y['s1'], s2 = y['s2']))
                     s1LatLong <- as.matrix(data_ss[match(Site_comparisons$s1, data_ss$SSBS), c('Longitude','Latitude')])
                     s2LatLong <- as.matrix(data_ss[match(Site_comparisons$s2, data_ss$SSBS), c('Longitude','Latitude')])
                     dist <- distHaversine(s1LatLong, s2LatLong)
                     Contrast <- paste(site_data$Landuse2[match(Site_comparisons$s1, site_data$SSBS)],
                                       site_data$Landuse2[match(Site_comparisons$s2, site_data$SSBS)],
                                       sep = "-")
                     study_results <- data.frame(Site_comparisons,
                                                 sor,
                                                 dist,
                                                 Contrast,
                                                 SS = s,
                                                 stringsAsFactors = TRUE)
                   }


# stop running things in parallel
registerDoSEQ()

## Environmental distance
r <- raster::getData('worldclim', var='bio', res=10) # dowload climate variables
alt <- raster::getData("worldclim", var="alt", res=10) #download elevation
#want BIO5(max temp of warmest) BIO6(min temp of coldest) BIO13(precip wettest) and 
#BIO14(precip driest)
r <- r[[c(5,6,13,14)]]
names(r) <- c("max_temp", "min_temp", "max_precip", "min_precip")

## Extracting data for PREDICTS sites
cd_data_climate <- raster::extract(r, cd_data_input[,37:38])
cd_data_alt <- raster::extract(alt, cd_data_input[,37:38])
cd_data_envi <- cbind(cd_data_climate,cd_data_alt)
cd_data_input <- cbind(cd_data_input, cd_data_envi)

s1Enviro <- as.data.frame(cd_data_input[match(cd_data$s1, cd_data_input$SSBS), 
                                        c('max_temp', 'min_temp', 'max_precip', 'min_precip','cd_data_alt')])
s2Enviro <- as.data.frame(cd_data_input[match(cd_data$s2, cd_data_input$SSBS), 
                                        c('max_temp', 'min_temp', 'max_precip', 'min_precip', 'cd_data_alt')])
enviro_dist <- gower_dist(s1Enviro, s2Enviro)

cd_data <- cbind(cd_data, enviro_dist)

## Want to add on information on fire presence:
s1_fire <- as.data.frame(cd_data_input[match(cd_data$s1, cd_data_input$SSBS), 
                                       c('fire_presence')])
s2_fire <- as.data.frame(cd_data_input[match(cd_data$s2, cd_data_input$SSBS), 
                                       c('fire_presence')])
cd_fire <- cbind(s1_fire, s2_fire)
names(cd_fire) <- c("s1_fire", "s2_fire")

cd_data <- cbind(cd_data, cd_fire)

cd_data$fire_contrast <- paste(cd_data$s1_fire, cd_data$s2_fire, sep = "-")

## Editting the comp sim data for modelling
cd_data <- cd_data %>%
  separate(Contrast, c("s1_LandUse", "s2_LandUse"), sep = "-", remove = FALSE) %>%
  filter(s1_LandUse == "Baseline vegetation") %>%
  mutate(logitCS = logit(sor, adjust = 0.001, percents = FALSE)) %>%
  mutate(log10geo = log10(dist + 1)) %>%
  mutate(cuberoot_envi = enviro_dist^(1/3)) %>%
  mutate(Contrast = factor(Contrast),
         Contrast = relevel(Contrast, ref = "Baseline vegetation-Baseline vegetation"))


## Setting up the abundance data and model:
## Check complete cases:
ab_model_data <- drop_na(abundance_data, RescaledAbundance, Landuse2, fire_presence)
## Transform abundance
ab_model_data <- mutate(ab_model_data, 
                        logAbundance = log(RescaledAbundance + 1),
                        sqrtAbundance = sqrt(RescaledAbundance))
## Set up time since fire data:
str(ab_model_data$days_since_fire)
ab_model_data$days_since_fire <- as.numeric(ab_model_data$days_since_fire)
ab_model_data$days_since_fire[is.na(ab_model_data$days_since_fire)] <- 1277
ab_model_data$time_since_fire <- ab_model_data$days_since_fire/365

source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(ab_model_data[ , c("Landuse2", "time_since_fire")])

## Model
ab_m <- lmer(sqrtAbundance ~ Landuse2 + time_since_fire + Landuse2:time_since_fire +
               (1|SS) + (1|SSB), data = ab_model_data)
Anova(ab_m)
summary(ab_m) 

## Model simplification - remove the interaction
ab_m.1 <- update(ab_m,~.- Landuse2:time_since_fire)
anova(ab_m, ab_m.1)
Anova(ab_m.1)
## Need to leave time since fire in the model. 
summary(ab_m.1)
## Because we cant simplify the compositional similarity model we're going to use the maximal
#models for each. 

## Diagnostics:
ab_output <- simulateResiduals(ab_m)
plot(ab_output)
ab_resids <- residuals(ab_m)
summary(ab_resids)
hist(ab_resids)

par(mfrow=c(1,2))
qqnorm(resid(ab_m), main="normal qq-plot, residuals")
qqline(resid(ab_m)) ## not bad
plot(fitted(ab_m), resid(ab_m), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)
ran <- ranef(ab_m, condVar=TRUE)
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
grid.arrange(a,b, nrow=2, ncol=1)

## Plot:
theme_set(theme_classic())
plot_model(ab_m, type = "pred", terms = c("time_since_fire", "Landuse2"),
           legend.title = c("Land use"),
           axis.title = c("Time since fire (years)",  "Square root rescaled abundance"))

r.squaredGLMM(ab_m)

## Compositional similarity:
str(cd_data_input$days_since_fire)
cd_data_input$days_since_fire <- as.numeric(cd_data_input$days_since_fire)
cd_data_input$days_since_fire[is.na(cd_data_input$days_since_fire)] <- 1277
cd_data_input$time_since_fire <- cd_data_input$days_since_fire/365

s1_time <- as.data.frame(cd_data_input[match(cd_data$s1, cd_data_input$SSBS),
                                       c('time_since_fire')])
s2_time <- as.data.frame(cd_data_input[match(cd_data$s2, cd_data_input$SSBS),
                                       c('time_since_fire')])
cd_time <- cbind(s1_time, s2_time)
names(cd_time) <- c("s1_time", "s2_time")

cd_data <- cbind(cd_data, cd_time)
cd_data$time_since_fire <- cd_data$s2_time

cd_data <- drop_na(cd_data, Contrast, log10geo, cuberoot_envi, time_since_fire)
corvif(cd_data[ , c("Contrast", "log10geo", "cuberoot_envi", "time_since_fire")])

## Check normality of logitCS
hist(cd_data$logitCS)
qqnorm(cd_data$logitCS)

cd_m <- lmer(logitCS ~ Contrast + log10geo + cuberoot_envi + time_since_fire +
               Contrast:time_since_fire +
               (1|SS), data = cd_data)
Anova(cd_m) ## All significant but due to pseudoreplication they can't be assumed to be correct. 
summary(cd_m)

## Diagnostics
cd_m_output <- simulateResiduals(cd_m)
plot(cd_m_output)
cd_resids <- residuals(cd_m)
summary(cd_resids)
hist(cd_resids)

par(mfrow=c(1,2))
qqnorm(resid(cd_m), main="normal qq-plot, residuals")
qqline(resid(cd_m)) ## not bad
plot(fitted(cd_m), resid(cd_m), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

ran <- ranef(cd_m, condVar=TRUE)
dotplot(ran)[["SS"]]

## R2
MuMIn::r.squaredGLMM(cd_m)

## Plot:
plot_model(cd_m, type = "pred", terms = c("time_since_fire", "Contrast"),
           axis.title = c("Time since fire (years)",  "Square root rescaled abundance"))

###############################################################################################
## Projecting models and BII

## Abundance:
new_data_ab <- expand.grid(time_since_fire = seq(0,3.5, by=0.1), 
                           Landuse2 = levels(ab_model_data$Landuse2))

new_data_ab$ab_predicted <- predict(ab_m, newdata = new_data_ab, re.form=NA) ^ 2

## Compositional similarity:
inv_logit <- function(f, a){
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}

new_data_cd <- expand.grid(time_since_fire = seq(0,3.5, by=0.1), 
                           Contrast = levels(cd_data$Contrast),
                           log10geo = 0,
                           cuberoot_envi = 0)

new_data_cd$cd_predicted <-predict(cd_m, newdata = new_data_cd, re.form=NA) %>%
  inv_logit(a = 0.001)

## BII
bii <- new_data_ab$ab_predicted*new_data_cd$cd_predicted

bii_df <- new_data_ab
bii_df <- cbind(bii_df, new_data_cd$cd_predicted)
colnames(bii_df)
names(bii_df)[names(bii_df) == "new_data_cd$cd_predicted"] <- "cd_predicted"
## Want to divide the predicted values by a reference value - reference value wants to be
#primary minimal with little fire disturbance so 3.5 years? 
ref_value <- subset(bii_df, bii_df$Landuse2 == "Baseline vegetation" & bii_df$time_since_fire == "3.5")

bii_df$ab_predicted <- bii_df$ab_predicted/ref_value$ab_predicted
bii_df$cd_predicted <- bii_df$cd_predicted/ref_value$cd_predicted

bii_df$bii <- (bii_df$ab_predicted*bii_df$cd_predicted)*100

plot(bii_df$time_since_fire, bii_df$bii)

ggplot(data=bii_df, aes(x=time_since_fire, y=bii, group=Landuse2)) +
  geom_line(aes(col=Landuse2))+
  xlab(label = "Time since fire (years)") +
  ylab(label = "BII (%)")+
  labs(col = "Land use")+
  theme_classic()
