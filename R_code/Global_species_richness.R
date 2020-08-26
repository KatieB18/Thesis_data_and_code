rm(list = ls())

## Load the required packages
library(predictsFunctions) 
library(StatisticalModels)
library(raster)
library(dplyr) 
library(tidyr) 
library(lme4)
library(car) 
library(DHARMa) 
library(MuMIn) 
library(rmarkdown)
library(sjPlot)
library(yarg)
library(roquefort)
library(performance)
library(optimx)
library(dfoptim)

## Load data
diversity <- readRDS("~/Imperial/PROJECT/PREDICTS data/database.rds")
fire_hist <- read.csv("~/Imperial/PROJECT/Combined data/Fire_history.csv")

##############################################################################################
## PREDICTS data preparation - calculating site level biodiversity and tidying up land use. 
diversity <- mutate(diversity, 
                    Measurement = Effort_corrected_measurement,
                    Sampling_effort = Rescaled_sampling_effort)

sites <- diversity %>%
  mutate(Diversity_metric_is_valid = TRUE) %>%
  predictsFunctions::SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use", "Biome", "Country")) %>%
  group_by(SS) %>%
  mutate(MaxAbundance = ifelse(Diversity_metric_type == "Abundance",
                               max(Total_abundance),
                               NA)) %>%
  ungroup() %>%
  mutate(RescaledAbundance = ifelse(Diversity_metric_type == "Abundance",
                                    Total_abundance/MaxAbundance,
                                    NA))

sites <- sites %>%
  mutate(
    Predominant_land_use = recode_factor(Predominant_land_use, 
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

## For my analysis I need to subset so that all sites are after 2006 (so that we can have 3yrs of
#fire data)
sites_post_2006 <- subset(sites, sites$Sample_start_earliest > "2006-01-01")
## Removing rows that have NAs in latitude and longitude as these won't have fire data
sites_post_2006_NA <- subset(sites_post_2006, (! is.na(Longitude) | ! is.na(Latitude)))

###############################################################################################
## Merge with fire data
PREDICTS_fire_hist <- merge(sites_post_2006_NA, fire_hist
                            [,c("SSBS", "Freq", "fire_ID", "size", "start_date", "duration", "days_since_fire",
                                "year", "years_since_fire")], by="SSBS", all.x=TRUE)

## editing of data
## Change freq from NA to 0
PREDICTS_fire_hist$Freq[is.na(PREDICTS_fire_hist$Freq)] <- 0
## Adding a column of fire presence
PREDICTS_fire_hist$fire_presence <- ifelse(PREDICTS_fire_hist$Freq==0, "No",
                                           ifelse(PREDICTS_fire_hist$Freq!=0, "Yes",
                                                  NA))
PREDICTS_fire_hist$fire_presence <- as.factor(PREDICTS_fire_hist$fire_presence)

###############################################################################################
## Data exploration
model_data <- drop_na(PREDICTS_fire_hist, Species_richness, Predominant_land_use, Use_intensity, fire_presence)

table(model_data$fire_presence)
## 9579 sites without fire and 542 sites with fire

table(model_data$Biome, model_data$fire_presence)
## There are huge differences in data availability across biomes. 

## remove the levels with no data/Can't have fires 
model_data<-model_data[!(model_data$Biome == "Rock & Ice"),]
model_data<-model_data[!(model_data$Biome == "Flooded Grasslands & Savannas"),]
model_data<-model_data[!(model_data$Biome == "Mangroves"),]
model_data<-model_data[!(model_data$Biome == "Inland Water"),]
model_data<-model_data[!(model_data$Biome == "Tundra"),]
model_data$Biome <- as.factor(as.character(model_data$Biome))
table(model_data$Biome)

## Could try merging biomes together..
## Forest vs non forest:
model_data$adj.biome <- ifelse(model_data$Biome=="Boreal Forests/Taiga", "Forest",
                               ifelse(model_data$Biome=="Temperate Conifer Forests", "Forest",
                                      ifelse(model_data$Biome=="Temperate Broadleaf & Mixed Forests", "Forest",
                                             ifelse(model_data$Biome=="Montane Grasslands & Shrublands", "Non Forest",
                                                    ifelse(model_data$Biome=="Temperate Grasslands, Savannas & Shrublands", "Non Forest",
                                                           ifelse(model_data$Biome=="Mediterranean Forests, Woodlands & Scrub", "Forest",
                                                                  ifelse(model_data$Biome=="Deserts & Xeric Shrublands", "Non Forest",
                                                                         ifelse(model_data$Biome=="Tropical & Subtropical Grasslands, Savannas & Shrublands", "Non Forest",
                                                                                ifelse(model_data$Biome=="Tropical & Subtropical Coniferous Forests", "Forest",
                                                                                       ifelse(model_data$Biome=="Tropical & Subtropical Dry Broadleaf Forests", "Forest",
                                                                                              ifelse(model_data$Biome=="Tropical & Subtropical Moist Broadleaf Forests", "Forest",
                                                                                                     NA)))))))))))
model_data$adj.biome <- as.factor(model_data$adj.biome)
table(model_data$Biome, model_data$adj.biome)
table(model_data$adj.biome, model_data$fire_presence)
## we have a lot more non forest fires than forest fires

## More representative biome spilts
model_data$new.biome <- model_data$Biome
model_data <- model_data %>%
  mutate(
    new.biome = recode_factor(new.biome, 
                              "Deserts & Xeric Shrublands" = "Deserts, Shrublands and Grasslands", 
                              "Montane Grasslands & Shrublands" = "Deserts, Shrublands and Grasslands",
                              "Temperate Grasslands, Savannas & Shrublands" = "Deserts, Shrublands and Grasslands",
                              "Temperate Broadleaf & Mixed Forests" = "Temperate Forest",
                              "Temperate Conifer Forests" = "Temperate Forest",
                              "Tropical & Subtropical Coniferous Forests" = "Tropical & Subtropical Forest",
                              "Tropical & Subtropical Dry Broadleaf Forests" = "Tropical & Subtropical Forest",
                              "Tropical & Subtropical Moist Broadleaf Forests" = "Tropical & Subtropical Forest")
  )
model_data$new.biome <- as.factor(as.character(model_data$new.biome))
table(model_data$new.biome, model_data$fire_presence)
table(model_data$new.biome, model_data$Predominant_land_use) #some land use classes need to be
#removed/merged


# Use intensity ok.
table(model_data$Predominant_land_use, model_data$Use_intensity)
table(model_data$Use_intensity, model_data$fire_presence)

table(model_data$Predominant_land_use, model_data$fire_presence)
## Need to remove NAs and plantation and urban, and then merge mature and intermediate secondary
#vegetation. 
model_data<- drop_na(model_data, Predominant_land_use)
model_data<-model_data[!(model_data$Predominant_land_use == "Plantation forest"),]
model_data<-model_data[!(model_data$Predominant_land_use == "Urban"),]
model_data$Predominant_land_use <- as.factor(as.character(model_data$Predominant_land_use))
model_data <- model_data %>%
  mutate( Predominant_land_use = recode_factor(Predominant_land_use, 
                                               "Intermediate secondary vegetation" = "Older secondary vegetation", 
                                               "Mature secondary vegetation" = "Older secondary vegetation"),
          Predominant_land_use = factor(Predominant_land_use),
          Predominant_land_use = relevel(Predominant_land_use, ref = "Primary vegetation")
  )
model_data$Predominant_land_use <- factor(model_data$Predominant_land_use, levels = c("Primary vegetation", "Young secondary vegetation", "Older secondary vegetation",
                                                                                      "Cropland", "Pasture"))
table(model_data$Predominant_land_use, model_data$fire_presence)

###############################################################################################
## Assumption checking:
## Checking collinearity 
source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence")])
## GVIF values low

## complete cases
model_data <- drop_na(model_data, Species_richness, Predominant_land_use, Use_intensity, fire_presence)

## Looking at distribution of species richness (models assume normality)
par(mfrow=c(1,2))
hist(model_data$Species_richness, xlab="Species richness", main="Histogram of species richness") 
## At the moment it shows a skewed distribution 
qqnorm(model_data$Species_richness) ## very curved.

# Need to model with poisson errors

## Initial model:
m1 <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + fire_presence +
              Predominant_land_use:Use_intensity + Predominant_land_use:fire_presence +
              (1|SS) + (1|SSB) + (1|Biome), data = model_data, family = poisson)
Anova(m1)
summary(m1)

## Diagnostics:
m1_output <- simulateResiduals(fittedModel = m1, n = 250)
plot(m1_output)

## Data is overdispersed:
testDispersion(m1_output)
check_overdispersion(m1)

## going to try fitting a site level random effect to account for this.
m2 <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + fire_presence +
              Predominant_land_use:Use_intensity + Predominant_land_use:fire_presence +
              (1|SS) + (1|SSB) + (1|Biome)+ (1|SSBS), data = model_data, family = poisson)
## Model fails to converge so checking for a false positive

m2_check <- with(m2@optinfo$derivs,solve(Hessian,gradient))
max(abs(m2_check)) # = 0.0003968403 - I think this implies there is a false positive here
## Checking different optimisers
m2_allfit <- allFit(m2, data = model_data)
m2_sum <- summary(m2_allfit)
m2_sum
m2_sum$which.OK #shows which optimisers worked 
## Model converges when using different optimiser

m2_refit <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + fire_presence +
                     Predominant_land_use:Use_intensity + Predominant_land_use:fire_presence +
                     (1|SS) + (1|SSB) + (1|Biome)+ (1|SSBS), data = model_data, family = poisson, 
                   control = glmerControl(optimizer = "bobyqa"))
## No convergence warning...
Anova(m2_refit)
summary(m2_refit)

m2refit_output <- simulateResiduals(fittedModel = m2_refit, n = 250)
plot(m2_refit)
plot(m2refit_output)
check_overdispersion(m2_refit) #No overdispersion detected 

## Random effects structure already optimised (want to keep the same as the abundance model)

## Fixed effects structure:
max <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + fire_presence +
              Predominant_land_use:Use_intensity + Predominant_land_use:fire_presence +
              Use_intensity:fire_presence +
              (1|SS) + (1|SSB) + (1|Biome)+ (1|SSBS), data = model_data, family = poisson, 
            control = glmerControl(optimizer = "bobyqa"))

Anova(max) # Use intensity fire presence interaction isn't significant so we can remove this from 
#the model and just use the previous version (m2_refit2)

anova(max, m2_refit) ## No significant explanatory power lost. 

## Looking at selected model:
richness_model <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + fire_presence +
                          Predominant_land_use:Use_intensity + Predominant_land_use:fire_presence +
                          (1|SS) + (1|SSB) + (1|Biome)+ (1|SSBS), data = model_data, family = poisson, 
                        control = glmerControl(optimizer = "bobyqa"))
Anova(richness_model)
summary(richness_model)

## R2:
R2GLMER(richness_model)

## Diagnostics
par(mfrow=c(1,2))
qqnorm(resid(richness_model), main="normal qq-plot, residuals")
qqline(resid(richness_model)) ## not bad
plot(fitted(richness_model), resid(richness_model), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

ran <- ranef(richness_model, condVar=TRUE)
#double square bracket access the lists in d1
a <- dotplot(ran)[["SSBS"]]
b <- dotplot(ran)[["SSB"]]
c <- dotplot(ran)[["SS"]]
d <- dotplot(ran)[["Biome"]]
grid.arrange(a,b,c,d, nrow=2, ncol=2) 

# Plot: 
set_theme(base = theme_classic(), legend.item.backcol = "white")
plot_model(richness_model, type = "pred", terms = c("Predominant_land_use", "fire_presence"), 
           legend.title = "Fire presence", show.intercept = TRUE, 
           axis.title = c("Predominant land use", "Species richness")) +
  theme(axis.text.x = element_text(angle = 10, hjust = 0.5, vjust = 0.75))

### testing hypothesis:
richness_fire <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + fire_presence +
                         Predominant_land_use:Use_intensity + Predominant_land_use:fire_presence +
                         (1|SS) + (1|SSB) + (1|Biome)+ (1|SSBS), data = model_data, family = poisson, 
                       control = glmerControl(optimizer = "bobyqa"))

richness_nofire <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + 
                           Predominant_land_use:Use_intensity + 
                           (1|SS) + (1|SSB) + (1|Biome)+ (1|SSBS), 
                         data = model_data, family = poisson, 
                         control = glmerControl(optimizer = "bobyqa"))
anova(richness_fire, richness_nofire)
## Highly significant effect of fire on richness. 


###########################################################################################
## Comparing savanna vs rest of the world.
## savanna/not savanna level
model_data$savanna <- ifelse(model_data$Biome== "Tropical & Subtropical Grasslands, Savannas & Shrublands","Savanna",
                             ifelse(model_data$Biome != "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                    "Not_savanna", NA))
model_data$savanna <- as.factor(model_data$savanna)
levels(model_data$savanna)
table(model_data$Biome, model_data$savanna)

## Making sure we dont have rank deficiency
model_data <- model_data %>%
  mutate(
    Use_intensity_2 = recode_factor(Use_intensity, 
                                    "Light use" = "Greater use", 
                                    "Intense use" = "Greater use"),
    Use_intensity_2 = factor(Use_intensity_2),
    Use_intensity_2 = relevel(Use_intensity_2, ref = "Minimal use")
  )
table(model_data$Predominant_land_use, model_data$Use_intensity_2)

## Creating a second land use classs. 
model_data$Landuse <- model_data$Predominant_land_use
model_data <- model_data %>%
  mutate(Landuse = recode_factor(Landuse, "Young secondary vegetation" = "Secondary vegetation",
                                 "Older secondary vegetation" = "Secondary vegetation"),
         Landuse = factor(Landuse),
         Landuse = relevel(Landuse, ref = "Primary vegetation"))
model_data$Landuse <- factor(model_data$Landuse, levels = c("Primary vegetation",
                                                            "Secondary vegetation",
                                                            "Cropland",
                                                            "Pasture"))
table(model_data$Landuse, model_data$Use_intensity_2, model_data$savanna)

savanna_comp <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + fire_presence +
                        Predominant_land_use:Use_intensity + Predominant_land_use:fire_presence +
                        fire_presence:savanna +
                        (1|SS) + (1|SSB) + (1|Biome)+ (1|SSBS), data = model_data, family = poisson, 
                      control = glmerControl(optimizer = "bobyqa"))
Anova(savanna_comp)
## Not significant therefore the effects of fire do not vary significiant between savanna and not

no_savanna <- glmer(Species_richness ~ Predominant_land_use + Use_intensity + fire_presence +
                      Predominant_land_use:Use_intensity + Predominant_land_use:fire_presence +
                      (1|SS) + (1|SSB) + (1|Biome)+ (1|SSBS), data = model_data, family = poisson, 
                    control = glmerControl(optimizer = "bobyqa"))
anova(savanna_comp, no_savanna)
