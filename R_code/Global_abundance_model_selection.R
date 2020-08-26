## Selection of global abundance model
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
library(ggplot2)
library(ggpubr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lwgeom)
library(tidyverse)

## Load data
diversity <- readRDS("~/Imperial/PROJECT/PREDICTS data/database.rds")
fire_hist <- read.csv("~/Imperial/PROJECT/Combined data/Fire_history.csv")

##############################################################################################
## PREDICTS data preparation 
## merging measurements and sampling efforts
diversity <- mutate(diversity, 
                    Measurement = Effort_corrected_measurement,
                    Sampling_effort = Rescaled_sampling_effort)

## Calculate diversity level metrics including calculated rescaled abundance 
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

###############################################################################################
## Merge with fire data
PREDICTS_fire_hist <- merge(sites_post_2006_NA, fire_hist
                            [,c("SSBS", "Freq", "fire_ID", "size", "start_date", "duration", "days_since_fire",
                                "year", "years_since_fire")], by="SSBS", all.x=TRUE)

## Change frequency from NA to 0
PREDICTS_fire_hist$Freq[is.na(PREDICTS_fire_hist$Freq)] <- 0
## Adding a column of fire presence
PREDICTS_fire_hist$fire_presence <- ifelse(PREDICTS_fire_hist$Freq==0, "No",
                                           ifelse(PREDICTS_fire_hist$Freq!=0, "Yes",
                                                  NA))
PREDICTS_fire_hist$fire_presence <- as.factor(PREDICTS_fire_hist$fire_presence)

###############################################################################################
## Data exploration
model_data <- drop_na(PREDICTS_fire_hist, Total_abundance, Predominant_land_use, Use_intensity, fire_presence)

table(model_data$fire_presence)
## 8235 sites without fire and 531 sites with fire
table(model_data$Predominant_land_use, model_data$fire_presence)
## Plantation forest and urban are lacking in sites so these will probably need to be excluded 
#as in the savanna only models. 
table(model_data$Biome, model_data$fire_presence)
## There are huge differences in data availability across biomes. 

## remove the levels with no data/Can't have fires 
model_data<-model_data[!(model_data$Biome == "Rock & Ice"),]
model_data<-model_data[!(model_data$Biome == "Flooded Grasslands & Savannas"),]
model_data<-model_data[!(model_data$Biome == "Mangroves"),]
model_data<-model_data[!(model_data$Biome == "Inland Water"),]
model_data<-model_data[!(model_data$Biome == "Tundra"),]
model_data$Biome <- as.factor(as.character(model_data$Biome))
table(model_data$Biome, model_data$fire_presence)

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
table(model_data$Biome, model_data$adj.biome)
table(model_data$adj.biome, model_data$fire_presence)
## we have a lot more non forest fires than forest fires. There is still a large imbalance globally
model_data$adj.biome <- as.factor(model_data$adj.biome)

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
table(model_data$new.biome, model_data$Predominant_land_use)

## Looking at land use categories, and aggregating categories to avoid data imbalance
## Need to remove NAs and plantation and urban
model_data<- drop_na(model_data, Predominant_land_use)
model_data<-model_data[!(model_data$Predominant_land_use == "Plantation forest"),]
model_data<-model_data[!(model_data$Predominant_land_use == "Urban"),]
model_data$Predominant_land_use <- as.factor(as.character(model_data$Predominant_land_use))
## merge intermediate and mature secondary together to avoid serious data imbalances. 
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

################################################################################################
## Assumption checking before modelling:
## Checking collinearity 
source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence")])
## GVIF values are all quite low. 

## complete cases
model_data <- drop_na(model_data, Total_abundance, Predominant_land_use, Use_intensity, fire_presence)

## Need to transform abundance
par(mfrow=c(1,2))
hist(model_data$RescaledAbundance, xlab="Rescaled abundance", main="Histogram of rescaled abundance") 
## At the moment it shows a skewed distribution 
qqnorm(model_data$RescaledAbundance) ## curve (skewed data) and long tails
model_data <- mutate(model_data, 
                     logAbundance = log(RescaledAbundance + 1),
                     sqrtAbundance = sqrt(RescaledAbundance)
)
dev.off()
par(mfrow=c(2,2))
hist(model_data$logAbundance, xlab="Log+1 abundance", main="Histogram of log+1 abundance") ## still shows a skew
qqnorm(model_data$logAbundance, main = "Normal Q-Q plot log+1 abundance") ## shows quite strong tails on either end also a slight curve indicating
#skewed data. 
hist(model_data$sqrtAbundance, xlab="Square root abundance", main="Histogram of square root abundance") ## Shows a normal distribution
qqnorm(model_data$sqrtAbundance, main = "Normal Q-Q plot sqaure root abundance") ## also shows quite distinct tails

#################################################################################################
## Model selection:

## Random effects:
## Models with and without biome as a random intercept and biome as a random slope.

## First going to create a maximal model:
max <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
              Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
              Use_intensity:fire_presence +
              (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(max) ## all elements significant except fire presence on its own
summary(max)

# Biome intercept
max_ref1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                   Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                   Use_intensity:fire_presence +
                   (1|SS) + (1|SSB) + (1|Biome), data = model_data)

# Aggregated biomes intercept
max_ref2 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                   Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                   Use_intensity:fire_presence +
                   (1|SS) + (1|SSB) + (1|new.biome), data = model_data)

#Biome slope
max_ref3 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                   Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                   Use_intensity:fire_presence +
                   (1|SS) + (1|SSB) + (1+fire_presence|Biome), data = model_data)
## Singular fit warning

# Aggregated biomes slope
max_ref4 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                   Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                   Use_intensity:fire_presence +
                   (1|SS) + (1|SSB) + (1+fire_presence|new.biome), data = model_data)
## Singular fit warning

# Without biome
max_ref5 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                   Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                   Use_intensity:fire_presence +
                   (1|SS) + (1|SSB), data = model_data)

AIC(max_ref1, max_ref2, max_ref3, max_ref4, max_ref5)
# lowest AIC = random effects structure 1 (biome intercept)
summ.table <- do.call(rbind, lapply(list(max_ref1, max_ref2, max_ref3, max_ref4, max_ref5), broom::glance))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid.Df", "Resid.Dev", "AIC")
reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
model.names <- c("SS, SSB and biome intercepts", "SS, SSB and aggregated biomes intercepts", 
                 "Biome slope", "Aggregated biome slope", "SS and SSB only")
row.names(reported.table) <- model.names
View(reported.table)

## Fixed effects
Anova(max)
## All elements expect fire presence have significant effect on abundance. Fire presence 
#stays in the model as it is involved in two significant interactions. 

### Hypothesis testing whether fire has an effect on biodiversity
m_with_fire <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                      Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                      Use_intensity:fire_presence +
                      (1|SS) + (1|SSB) + (1|Biome), data = model_data)

m_no_fire <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + 
                    Predominant_land_use:Use_intensity +
                    (1|SS) + (1|SSB) + (1|Biome), data = model_data)

anova(m_with_fire, m_no_fire)
## Lost significant explanatory power therefore fire presence and its interaction with other factors
# is having an important influence on abundance. 

###### Looking at the selected model in more detail:
max <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
              Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
              Use_intensity:fire_presence +
              (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(max) ## all elements significant except fire presence on its own
summary(max)
# Model suggest that the presence of fire has different effects in different land uses compared to
#primary

## Diagnostics
## USing DHARMa
max_output <- simulateResiduals(fittedModel = max, n = 250)
plot(max_output)
# Signidicant deviation on QQ plot, and some outliers in residuals vs predictors

## Plotting residuals against each predictor
par(mfrow=c(1,3))
plotResiduals(max_output, model_data$Predominant_land_use, main="Predominant land use")
plotResiduals(max_output, model_data$Use_intensity, main="Use intensity")
plotResiduals(max_output, model_data$fire_presence, main="fire presence")
dev.off()
## all apears ok

## Residuals
residuals <- resid(max)
summary(residuals) # seem ok (small disparity in max and min but nothing too serious and mean
#is around 0 indicated normality)
hist(residuals) ## Residuals look ok

## residuals qq plot and residuals vs fitted plot
par(mfrow=c(1,2))
qqnorm(resid(max), main="normal qq-plot, residuals")
qqline(resid(max)) ## not bad
plot(fitted(max), resid(max), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

## Normality of random effects
ran <- ranef(max, condVar=TRUE)
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
c <- dotplot(ran)[["Biome"]]
grid.arrange(a,b,c, nrow=2, ncol=2) 

## Plots:
set_theme(base = theme_classic(), legend.item.backcol = "white")

## Land use fire presence interaction plot
LU_fp <- plot_model(max, type = "pred", terms = c("Predominant_land_use", "fire_presence"), legend.title = "Fire presence",
                    show.intercept = TRUE, axis.title = c("Predominant land use", "Square root of rescaled abundance")) +
  theme(axis.text.x = element_text(angle = 10, hjust = 0.5, vjust = 0.75))
LU_fp

## Land use use intensity interaction plot
LU_UI <- plot_model(max, type = "pred", terms = c("Predominant_land_use", "Use_intensity"), legend.title = "Use intensity",
                    show.intercept = TRUE, axis.title = c("Predominant land use", "Square root of rescaled abundance")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
LU_UI

## Use intensity fire presence plot
UI_fp <- plot_model(max, type = "pred", terms = c("Use_intensity", "fire_presence"), legend.title = "Fire presence",
                    show.intercept = TRUE, axis.title = c("Land use intensity", "Square root of rescaled abundance")) 
UI_fp

## Land use plot
LU <- plot_model(max, type = "pred", terms = "Predominant_land_use", show.intercept = TRUE,
                 axis.title = c("Predominant land use", "Square root of rescaled abundance"))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
LU

ggarrange(LU_fp, UI_fp, nrow = 2, ncol = 1, labels = c("A", "B"), 
          common.legend = TRUE,legend = "right")

## Generating a summary table:
Global_abundance_table <- tab_model(max, show.stat = TRUE, show.se = TRUE, 
                                    show.df = TRUE, string.stat = c("t statistic"))
Global_abundance_table

r.squaredGLMM(max)
# fixed effects = 0.036, random effects = 0.607

##############################################################################################
## Do the savannas behave differently?
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

savanna_comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                       Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                       Use_intensity:fire_presence + fire_presence:savanna +
                       (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(savanna_comp)
## Not significant therefore the effects of fire do not vary significiant between savanna and not

no_savanna <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                     Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                     Use_intensity:fire_presence + 
                     (1|SS) + (1|SSB) + (1|Biome), data = model_data)
anova(savanna_comp, no_savanna)

###############################################################################################
##### The following sections were not fully completed/did not make it into the final thesis

#################################################################################################
## Investigating which sites are pulling pasture down
Pasture_sites <- subset(model_data, model_data$Predominant_land_use == "Pasture")
## 1400 pasture sites 
table(Pasture_sites$fire_presence) ## In these sites only 38 have had fires
table(Pasture_sites$Biome, Pasture_sites$fire_presence) ## 26 are in savannas

## 4 sites in tropical and subtropical moist broadleaf forests and 7 in montane grasslands and 
#shrublands (1 in desert)
table(Pasture_sites$new.biome, Pasture_sites$fire_presence)

## Subset to just tropical forests 
pasture_tropforest <- subset(Pasture_sites, Pasture_sites$Biome == "Tropical & Subtropical Moist Broadleaf Forests")

plot(pasture_tropforest$fire_presence, pasture_tropforest$sqrtAbundance)

## Subset to just desrts, shrublands and grasslands
pasture_montane <- subset(Pasture_sites, Pasture_sites$Biome == "Montane Grasslands & Shrublands")

plot(pasture_montane$fire_presence, pasture_montane$sqrtAbundance)

## I believe the large decline in pasture is being partly caused by savanna sites and partly caused
#by tropical forest sites! 

#################################################################################################
## fire prone/fire sensitive biomes 
## Archibald et al., 2013 identifies 5 pyromes FIL, FCS, RIL, RCS, ICS. 
## Defines them based off biomes
table(model_data$Biome)

model_data$pyrome <- ifelse(model_data$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands", 
                            "Frequent", "Rare")
table(model_data$pyrome, model_data$fire_presence)

corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence", "pyrome")])
#Fab

pyrome_model <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                       Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                       Use_intensity:fire_presence + fire_presence:pyrome +
                       (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(pyrome_model)
summary(pyrome_model)
## Fire presence:pyrome interaction is not significant... 

## Trying frequent, intermediate and rare...
model_data$pyrome2 <- ifelse(model_data$Biome == "Boreal Forests/Taiga", "Rare",
                             ifelse(model_data$Biome == "Deserts & Xeric Shrublands", "Rare",
                                    ifelse(model_data$Biome == "Mediterranean Forests, Woodlands & Scrub", "Rare",
                                           ifelse(model_data$Biome == "Montane Grasslands & Shrublands", "Intermediate",
                                                  ifelse(model_data$Biome == "Temperate Broadleaf & Mixed Forests", "Intermediate",
                                                         ifelse(model_data$Biome == "Temperate Conifer Forests", "Rare",
                                                                ifelse(model_data$Biome == "Temperate Grasslands, Savannas & Shrublands", "Intermediate",
                                                                       ifelse(model_data$Biome == "Tropical & Subtropical Coniferous Forests", "Intermediate",
                                                                              ifelse(model_data$Biome == "Tropical & Subtropical Dry Broadleaf Forests", "Intermediate",
                                                                                     ifelse(model_data$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Frequent",
                                                                                            ifelse(model_data$Biome == "Tropical & Subtropical Moist Broadleaf Forests", "Intermediate",
                                                                                                   NA)))))))))))


table(model_data$pyrome2, model_data$fire_presence)

pyrome_model2 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                        Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                        Use_intensity:fire_presence + pyrome2:fire_presence +
                        (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(pyrome_model2)
summary(pyrome_model2)

## Trying a different biome separation...
model_data$fire_dependence <- ifelse(model_data$Biome == "Boreal Forests/Taiga", "Dependent",
                                     ifelse(model_data$Biome == "Deserts & Xeric Shrublands", "Dependent",
                                            ifelse(model_data$Biome == "Mediterranean Forests, Woodlands & Scrub", "Dependent",
                                                   ifelse(model_data$Biome == "Montane Grasslands & Shrublands", "Dependent",
                                                          ifelse(model_data$Biome == "Temperate Broadleaf & Mixed Forests", "Sensitive",
                                                                 ifelse(model_data$Biome == "Temperate Conifer Forests", "Dependent",
                                                                        ifelse(model_data$Biome == "Temperate Grasslands, Savannas & Shrublands", "Dependent",
                                                                               ifelse(model_data$Biome == "Tropical & Subtropical Coniferous Forests", "Sensitive",
                                                                                      ifelse(model_data$Biome == "Tropical & Subtropical Dry Broadleaf Forests", "Sensitive",
                                                                                             ifelse(model_data$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Dependent",
                                                                                                    ifelse(model_data$Biome == "Tropical & Subtropical Moist Broadleaf Forests", "Sensitive",
                                                                                                           NA)))))))))))
table(model_data$Biome, model_data$fire_dependence)
table(model_data$fire_dependence, model_data$fire_presence)
