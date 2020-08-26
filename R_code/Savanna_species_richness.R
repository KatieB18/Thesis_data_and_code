## Savanna species richness models

rm(list = ls())
## Packages:
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
library(roquefort)
library(yarg)
library(sjstats)
library(gt)
library(performance)


## Data: 
diversity <- readRDS("~/Imperial/PROJECT/PREDICTS data/database.rds")
fire_hist <- read.csv("~/Imperial/PROJECT/Combined data/Fire_history.csv")

#####################################################################################################
## Data preparation: 

diversity <- mutate(diversity, 
                    Measurement = Effort_corrected_measurement,
                    Sampling_effort = Rescaled_sampling_effort)
## Calculating diversity metrics including rescaled abundance. 
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
## Tidying up the predominant land use category 
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

## Merge with fire data
PREDICTS_fire_hist <- merge(sites_post_2006_NA, fire_hist
                            [,c("SSBS", "Freq", "fire_ID", "size", "start_date", "duration", "days_since_fire",
                                "year", "years_since_fire")], by="SSBS", all.x=TRUE)
## Creating fire presence column (need to change NAs in freq to 0 as these sites have no fires in the
#last 3 years)
PREDICTS_fire_hist$Freq[is.na(PREDICTS_fire_hist$Freq)] <- 0
PREDICTS_fire_hist$fire_presence <- ifelse(PREDICTS_fire_hist$Freq==0, "No",
                                           ifelse(PREDICTS_fire_hist$Freq!=0, "Yes",
                                                  NA))
PREDICTS_fire_hist$fire_presence <- as.factor(PREDICTS_fire_hist$fire_presence)

## Need to remove data which has NAs in categories required for our analysis
PREDICTS_fire_hist2 <- drop_na(PREDICTS_fire_hist, Species_richness, Predominant_land_use, Use_intensity, fire_presence)

## Subset to just the Tropical & Subtropical Grasslands, Savannas & Shrublands
sites_savanna <- subset(PREDICTS_fire_hist2, PREDICTS_fire_hist2$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands")

## Looking at the data
## We need to make sure we have sufficient data for analysis
table(sites_savanna$Predominant_land_use, sites_savanna$fire_presence)
# Removing land use level which lack data
sites_savanna<-sites_savanna[!(sites_savanna$Predominant_land_use == "Plantation forest"),]
sites_savanna<-sites_savanna[!(sites_savanna$Predominant_land_use == "Urban"),] 
sites_savanna$Predominant_land_use <- as.factor(as.character(sites_savanna$Predominant_land_use))
## Merge together intermediate and mature secondary vegetation into one class
sites_savanna <- sites_savanna %>%
  mutate( Predominant_land_use = recode_factor(Predominant_land_use, 
                                               "Young secondary vegetation" = "Secondary vegetation",
                                               "Intermediate secondary vegetation" = "Secondary vegetation", 
                                               "Mature secondary vegetation" = "Secondary vegetation"),
          Predominant_land_use = factor(Predominant_land_use),
          Predominant_land_use = relevel(Predominant_land_use, ref = "Primary vegetation")
  )
sites_savanna$Predominant_land_use <- factor(sites_savanna$Predominant_land_use, levels = c("Primary vegetation", "Secondary vegetation",
                                                                                            "Cropland", "Pasture"))

table(sites_savanna$Predominant_land_use, sites_savanna$fire_presence)

## Also want to be able to use interaction between land use and use intensity
table(sites_savanna$Predominant_land_use, sites_savanna$Use_intensity)
## If we use the data as is we will have rank deficiency and therefore it may be better to combine 
#use intensities into coarser land uses. 

sites_savanna <- sites_savanna %>%
  mutate(
    Use_intensity_2 = recode_factor(Use_intensity, 
                                    "Light use" = "Greater use", 
                                    "Intense use" = "Greater use"),
    Use_intensity_2 = factor(Use_intensity_2),
    Use_intensity_2 = relevel(Use_intensity_2, ref = "Minimal use")
  )
table(sites_savanna$Predominant_land_use, sites_savanna$Use_intensity_2)

###############################################################################################
## Checking asusmptions:

# Colinearity
source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(sites_savanna[ , c("Predominant_land_use", "Use_intensity_2", "fire_presence")])

# Complete cases
model_data <- drop_na(sites_savanna, Species_richness, Predominant_land_use, Use_intensity_2, fire_presence)

# Normality of the response variable
par(mfrow=c(1,2))
hist(model_data$Species_richness) ## At the moment it shows a skewed distribution 
qqnorm(model_data$Species_richness)
## Data is not normally distributed so we need to model with Poisson errors.

## Starting model:
m1 <- glmer(Species_richness ~ Predominant_land_use + Use_intensity_2 + fire_presence +
              Predominant_land_use:Use_intensity_2 + Predominant_land_use:fire_presence +
              (1|SS) + (1|SSB), data = model_data, family = poisson)
## Convergence warning: Need to check if this is a false positive
m1_check <- with(m1@optinfo$derivs,solve(Hessian,gradient))
max(abs(m1_check)) #0.000536
## Trying different optimisers
m1_all <- allFit(m1)
m1_sum <- summary(m1_all)
m1_sum$which.OK

##refit with diff optimiser
m1.1 <- glmer(Species_richness ~ Predominant_land_use + Use_intensity_2 + fire_presence +
                Predominant_land_use:Use_intensity_2 + Predominant_land_use:fire_presence +
                (1|SS) + (1|SSB), data = model_data, family = poisson,
              control = glmerControl(optimizer ="bobyqa"))
## Model now converges
Anova(m1.1)
summary(m1.1)

## Model diagnostics
m1_output <- simulateResiduals(fittedModel = m1.1, n = 250)
plot(m1_output)

check_overdispersion(m1)
## Some overdispersion

## Account for overdispersion using a site level random effect:
m2 <- glmer(Species_richness ~ Predominant_land_use + Use_intensity_2 + fire_presence + 
              Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2 +
              (1|SS) + (1|SSB) + (1|SSBS), data = model_data, family = poisson)
## Also failed to converge
m2_check <- with(m2@optinfo$derivs,solve(Hessian,gradient))
max(abs(m2_check)) #0.0049
## Trying a new optimiser
m2.1 <- glmer(Species_richness ~ Predominant_land_use + Use_intensity_2 + fire_presence + 
                Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2 +
                (1|SS) + (1|SSB) + (1|SSBS), data = model_data, family = poisson,
              control = glmerControl(optimizer ="bobyqa"))
# Now converges
Anova(m2.1)
summary(m2.1)

m2_output <- simulateResiduals(fittedModel = m2.1, n = 250)
plot(m2_output)

check_overdispersion(m2.1) 
## No longer overdispersed 

## Model selection:
## Want to keep the random effects structure used in the maximal model
## selecting fixed effects:

## Maximal model:
max <- glmer(Species_richness ~ Predominant_land_use + Use_intensity_2 + fire_presence + 
               Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2 +
               (1|SS) + (1|SSB) + (1|SSBS), data = model_data, family = poisson,
             control = glmerControl(optimizer ="bobyqa"))
Anova(max)
# Predominant land use and its interaction with fire presence is signficant! In the savanna fire
#presence does have a significant effect on species richness
# This implies that we can remove use intensity from this model, but because we want to compare 
#this model to the savanna abundance model we want to keep the elements the same?
max.1 <- update(max,~.- Predominant_land_use:Use_intensity_2)
anova(max, max.1)
## Haven't lost significant explanatory power
Anova(max.1)
## can also try removing use intensity entirely
max.2 <- update(max.1,~.- Use_intensity_2)
anova(max.1, max.2)
## Haven't lost significant explanatory power therefore best fixed effects structure for this model
# = just land use and fire presence

## Selected model =
richness_model <- glmer(Species_richness ~ Predominant_land_use + fire_presence + 
                          Predominant_land_use:fire_presence +
                          (1|SS) + (1|SSB) + (1|SSBS), data = model_data, family = poisson,
                        control = glmerControl(optimizer ="bobyqa"))

Anova(richness_model)
summary(richness_model)

## Diagnostics:
## Residuals
residuals <- resid(richness_model)
summary(residuals)  
hist(residuals) #ok 

par(mfrow=c(1,2))
qqnorm(resid(richness_model), main="normal qq-plot, residuals")
qqline(resid(richness_model)) ## not bad
plot(fitted(richness_model), resid(richness_model), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

ran <- ranef(richness_model, condVar=TRUE)
a <- dotplot(ran)[["SSBS"]]
b <- dotplot(ran)[["SSB"]]
c <- dotplot(ran)[["SS"]]
grid.arrange(a,b,c, nrow=2, ncol=2) 

## Plotting:
set_theme(base = theme_classic(), legend.item.backcol = "white")
plot_model(richness_model, type = "pred", terms = c("Predominant_land_use", "fire_presence"), 
           legend.title = "Fire presence", show.intercept = TRUE, 
           axis.title = c("Predominant land use", "Species richness"))

## Hypothesis testing.. Does fire matter for species richness
## compare max against model with no fire
fire <- glmer(Species_richness ~ Predominant_land_use + fire_presence + 
                Predominant_land_use:fire_presence +
                (1|SS) + (1|SSB) + (1|SSBS), data = model_data, family = poisson,
              control = glmerControl(optimizer ="bobyqa"))
no_fire <- glmer(Species_richness ~ Predominant_land_use + 
                   (1|SS) + (1|SSB) + (1|SSBS), data = model_data, family = poisson, 
                 control = glmerControl(optimizer ="bobyqa"))
anova(max, no_fire)
AIC(max, no_fire)
## Significant effect of fire on species richness in the savanna 