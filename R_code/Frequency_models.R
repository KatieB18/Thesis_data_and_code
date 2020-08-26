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
library(performance)
library(ape)
library(gridExtra)
library(gamm4)
library(mgcv)

## Load data
diversity <- readRDS("~/Imperial/PROJECT/PREDICTS data/database.rds")
fire_hist <- read.csv("~/Imperial/PROJECT/Combined data/Fire_history.csv")

##############################################################################################
## PREDICTS data preparation 
## Move Measurement and Effort_corrected_measurement into the same thing, and Sampling_effort and 
#Rescaled_sampling_effort 
diversity <- mutate(diversity, 
                    Measurement = Effort_corrected_measurement,
                    Sampling_effort = Rescaled_sampling_effort)

## Calculate diversity level metrics - calculated rescaled abundance 
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
  mutate(
    # collapse primary forest and non-forest together into primary vegetation as these aren't well 
    #distinguished
    Predominant_land_use = recode_factor(Predominant_land_use, 
                                         "Primary forest" = "Primary vegetation", 
                                         "Primary non-forest" = "Primary vegetation"),
    # indeterminate secondary veg and cannot decide get NA
    Predominant_land_use = na_if(Predominant_land_use, "Secondary vegetation (indeterminate age)"),
    Predominant_land_use = na_if(Predominant_land_use, "Cannot decide"),
    Use_intensity = na_if(Use_intensity, "Cannot decide"),
    # set reference levels
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

#################################################################################################
## Data exploration

model_data <- drop_na(PREDICTS_fire_hist, Predominant_land_use, Use_intensity, Freq)
## 10121 sites in total

## Remove biomes which lack data
table(model_data$Biome, model_data$Freq)
## There are huge differences in data availability across biomes. 

## remove the levels with no data/Can't have fires 
model_data<-model_data[!(model_data$Biome == "Rock & Ice"),]
model_data<-model_data[!(model_data$Biome == "Flooded Grasslands & Savannas"),]
model_data<-model_data[!(model_data$Biome == "Mangroves"),]
model_data<-model_data[!(model_data$Biome == "Inland Water"),]
model_data<-model_data[!(model_data$Biome == "Tundra"),]
model_data$Biome <- as.factor(as.character(model_data$Biome))

## Land use classes:
table(model_data$Predominant_land_use, model_data$Freq)
## Remove plantation and urban and merge intermediate and mature secondary to avoid low data levels
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
table(model_data$Predominant_land_use, model_data$Freq)

## Use intensity:
table(model_data$Use_intensity, model_data$Freq)

## Checking frequency is a continuous numeric vector.
str(model_data$Freq)

model_data <- model_data %>%
  mutate(Use_intensity_2 = recode_factor(Use_intensity, "Light use" = "Greater use",
                                         "Intense use" = "Greater use"),
         Use_intensity_2 = factor(Use_intensity_2),
         Use_intensity_2 = relevel(Use_intensity_2, ref = "Minimal use"))

model_data <- model_data %>%
  mutate(Landuse = recode_factor(Predominant_land_use, "Young secondary vegetation" = "Secondary vegetation",
                                 "Older secondary vegetation" = "Secondary vegetation"),
         Landuse = factor(Landuse),
         Landuse = relevel(Landuse, ref = "Primary vegetation"))
model_data$Landuse <- factor(model_data$Landuse, levels = c("Primary vegetation",
                                                            "Secondary vegetation",
                                                            "Cropland",
                                                            "Pasture"))


################################################################################################
## Modelling:

## check assumptions:
# collinearity:
source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(model_data[ , c("Landuse", "Use_intensity", "Freq")])
## all ok. 

## Normality of the response variable:
par(mfrow=c(1,2))
hist(model_data$Freq, xlab="Frequency", main="Histogram of fire frequency") 
## At the moment it shows a skewed distribution 
qqnorm(model_data$Freq)
## not great, way more 0 than anything else

m1 <- glmer(Freq ~ Landuse + Use_intensity +
              (1|SS) + (1|SSB) ,
            data = model_data, family = poisson)
## Worked!!!
Anova(m1) ## Both significant
summary(m1) ## V similar to savanna model

m1_output <- simulateResiduals(m1)
plot(m1_output)

## Model is underdispersed
check_overdispersion(m1) ## No over dispersion so I think this model is underdispered. 
## Underdispersion is not caused by zero inflation:
check_zeroinflation(m1)

#######
## Fitting spatial models

## glmmTMB
library(glmmTMB)
## Creat a numeric factor of coordinates of locations
model_data$pos <- numFactor(scale(model_data$Longitude), scale(model_data$Latitude))
# create a group factor to be used as a random term
model_data$ID <- factor(rep(1, nrow(model_data)))

# Model:
m_tmb <- glmmTMB(Freq ~ Landuse + Use_intensity + 
                   (1|SS) + (1|SSB)+
                   mat(pos + 0 | ID), model_data)
## Hasn't worked due to an error in memory allocation

## spaMM
library(spaMM)

m_spamm <- fitme(Freq ~ Landuse + Use_intensity + 
                   (1|SS) + (1|SSB) +
                   Matern(1 | Longitude + Latitude), 
                 data = model_data, family = "gaussian")
## Also has an error with memory allocation. 

## gamm
library(gamm4)
library(mgcv)

## Need to fixed mixed effects model as poisson
gamm_model <- gamm4(Freq ~ s(Longitude, Latitude) + Landuse + Use_intensity, 
                    random = ~(1|SS) + (1|SSB), data = model_data, family = poisson)
summary(gamm_model$mer)
summary(gamm_model$gam)
anova(gamm_model$mer)
Anova(gamm_model$gam)

check_overdispersion(gamm_model$mer) ## Model is still underdispersed. 

## Model selection
## Create a maximum model including an interaction between land use, use intensity and savanna/not 
#to test if the savanna behaves differently to the rest of the world or not. 

## Creating a savanna/not factor
model_data$savanna <- ifelse(model_data$Biome== "Tropical & Subtropical Grasslands, Savannas & Shrublands","Savanna",
                             ifelse(model_data$Biome != "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                    "Not_savanna", NA))
model_data$savanna <- as.factor(model_data$savanna)
levels(model_data$savanna)

## Looking at interaction effect between use intensity and land use
m1 <- gamm4(Freq ~ s(Longitude, Latitude) + Landuse + Use_intensity_2 + Landuse:Use_intensity_2, 
            random = ~(1|SS) + (1|SSB), data = model_data, family = poisson)
summary(m1$mer)
summary(m1$gam)
Anova(m1$mer)
## it seems that there is not significant interaction effect between land use and use intensity. 

## Interaction with savanna or not
m2 <- gamm4(Freq ~ s(Longitude, Latitude) + Landuse + Use_intensity_2 + Landuse:Use_intensity_2:savanna, 
            random = ~(1|SS) + (1|SSB), data = model_data, family = poisson)
summary(m2$mer)
summary(m2$gam)
Anova(m2$mer)

## Model comparison/selection:
## Settiing up functions required for MuMIn selection
gamm4 <- function(...) structure(c(gamm4::gamm4(...), list(call = match.call())),
                                 class = c("gamm", "list"))
logLik.gamm <- function(object, ...) logLik(object[[if (is.null(object$lme)) "mer" else "lme"]], ...)
formula.gamm <- function(x, ...) formula(x$gam, ...)
nobs.gamm <- function(object, ...) nobs(object$gam, ...)

nobs.gam <- function(object, ...) stats:::nobs.glm(object, ...)
coeffs.gamm <- function(model) coef(model$gam)
getAllTerms.gamm <- function(x, ...) getAllTerms(x$gam)
tTable.gamm <- function(model, ...) tTable(model$gam)

## Create maximal model
m3 <- uGamm(Freq ~ s(Longitude, Latitude) + Landuse + Use_intensity_2 + Landuse:Use_intensity_2:savanna, 
            random = ~(1|SS) + (1|SSB), data = model_data, family = poisson)
## all variations of the model compared together: 
(dd <- dredge(m3))
summary(model.avg(dd, subset = delta <= 8))
## Best model according to this uses only land use and use intensity as fixed effects - does not
#consider lat/lon or the interaction with savanna. 
## However, I want to retain the lat/lon in the model because it is important for accounting for 
#underdispersion. 

## Looking into models 
m4 <- gamm4(Freq ~ s(Longitude, Latitude) + Landuse + Use_intensity_2, 
            random = ~(1|SS) + (1|SSB), data = model_data, family = poisson)
m5 <- gamm4(Freq ~ s(Longitude, Latitude) + Use_intensity, 
            random = ~(1|SS) + (1|SSB), data = model_data, family = poisson)

summary(m4$mer)
summary(m4$gam)
anova(m4$gam)

summary(m5$mer)
summary(m5$gam)
anova(m5$gam)

AIC(m4, m5)
## Going to select m5 as the model (lowest AIC that contains a spatial term)

## Looking into the selected model:
check_overdispersion(m5$mer) ## still underdispersed

## Diagnostics:
par(mfrow=c(1,2))
qqnorm(resid(m5$mer), main="normal qq-plot, residuals")
qqline(resid(m5$mer)) ## not bad
plot(fitted(m5$mer), resid(m5$mer), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

ran <- ranef(m5$mer, condVar=TRUE)
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
grid.arrange(a,b, nrow=2, ncol=1) 

m5_mer <- m5$mer
summary(m5_mer)

## PLotting
library(mgcViz)
m5.2 <- gamm4V(Freq ~ s(Longitude, Latitude) + Use_intensity, 
               random = ~(1|SS) + (1|SSB), data = model_data, family = poisson)

plot.gamViz(m5.2, select = 2, xlab="Land use intensity") + 
  labs(title = NULL) +
  xlab("Land use intensity") +
  ylab("Change in fire frequency") +
  theme_classic()

Anova(m5_mer)

################################################################################################
## Savanna frequency models
Savanna_data <- subset(model_data, model_data$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands")

## Savanna data set up
## Land use
Savanna_data$Landuse <- Savanna_data$Predominant_land_use
Savanna_data <- Savanna_data %>%
  mutate(Landuse = recode_factor(Landuse, "Young secondary vegetation" = "Secondary vegetation",
                                 "Older secondary vegetation" = "Secondary vegetation"),
         Landuse = factor(Landuse),
         Landuse = relevel(Landuse, ref = "Primary vegetation"))
Savanna_data$Landuse <- factor(Savanna_data$Landuse, levels = c("Primary vegetation",
                                                                "Secondary vegetation",
                                                                "Cropland",
                                                                "Pasture"))
table(Savanna_data$Landuse, Savanna_data$Freq)

Savanna_data <- Savanna_data %>%
  mutate(
    Use_intensity_2 = recode_factor(Use_intensity, 
                                    "Light use" = "Greater use", 
                                    "Intense use" = "Greater use"),
    Use_intensity_2 = factor(Use_intensity_2),
    Use_intensity_2 = relevel(Use_intensity_2, ref = "Minimal use")
  )

table(Savanna_data$Landuse, Savanna_data$Use_intensity_2)
table(Savanna_data$Landuse, Savanna_data$Use_intensity) ## Will have rank deficiency if we use this 
#interaction. 

## Maximal model:
max <- glmer(Freq ~ Landuse + Use_intensity_2 + Landuse:Use_intensity_2 +
                (1|SS) + (1|SSB), data = Savanna_data, family = poisson)
## Doesn't converge
## try combining cropland and pasture together into one category to make the model more simple
#to plot.
Savanna_data <- Savanna_data %>%
  mutate(Landuse2 = recode_factor(Landuse, "Cropland" = "Agriculture",
                                  "Pasture" = "Agriculture"),
         Landuse2 = factor(Landuse2),
         Landuse2 = relevel(Landuse2, ref = "Primary vegetation"))
Savanna_data$Landuse2 <- factor(Savanna_data$Landuse2, levels = c("Primary vegetation",
                                                                  "Secondary vegetation",
                                                                  "Agriculture"))

max.1 <- glmer(Freq ~ Landuse2 + Use_intensity_2 + Landuse2:Use_intensity_2 +
                (1|SS) + (1|SSB), data = Savanna_data, family = poisson)
Anova(max.1) ## interaction is not significant nor is land use. 
summary(max.1)

## Neither land use or the interaction is significant.

## removing the interaction and not using the aggregated land use and use intensity classes
max.2 <- glmer(Freq ~ Landuse + Use_intensity + 
                 (1|SS) + (1|SSB), data = Savanna_data, family = poisson)
anova(max.2, max.1)
Anova(max.2) ## Land use and use intensity significant
summary(max.2)

## Model selection of the fixed effects therefore results in the selected model:
savanna_mod <- glmer(Freq ~ Landuse + Use_intensity +
                (1|SS) + (1|SSB), data = Savanna_data, family = poisson)
Anova(savanna_mod)
summary(savanna_mod)
## Increased land use intensity reduced fire frequency. 

## Diagnostics
savanna_output <- simulateResiduals(savanna_mod)
plot(savanna_output)
## No over dispersion but do have issues with normality and outliers
check_overdispersion(savanna_mod) #No overdispersion.

## Residuals 
residuals <- resid(savanna_mod)
summary(residuals)
hist(residuals)

par(mfrow=c(1,2))
qqnorm(resid(savanna_mod), main="normal qq-plot, residuals")
qqline(resid(savanna_mod)) ## not bad
plot(fitted(savanna_mod), resid(savanna_mod), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

ran <- ranef(savanna_mod, condVar=TRUE)
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
grid.arrange(a,b, nrow=2, ncol=1)

## Plots
set_theme(base = theme_classic(), axis.angle.x = 5, axis.textsize.x = 0.9)
LU <- plot_model(savanna_mod, type = "pred", terms = c("Landuse"), show.intercept = TRUE, 
                 axis.title = c("Predominant land use", "Frequency of fires"),
                 title = "", axis.lim = c(0,0.9)) 
UI <- plot_model(savanna_mod, type = "pred", terms = c("Use_intensity"), show.intercept = TRUE, 
                 axis.title = c("Land use intensity", "Frequency of fires"),
                 title = "", axis.lim = c(0,0.9))
ggarrange(LU, UI + 
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank()), 
          labels = c("A", "B"))

## R2
r.squaredGLMM(savanna_mod)

tab <- tab_model(savanna_mod, show.stat = TRUE, show.se = TRUE, 
                 show.df = TRUE, string.stat = c("z statistic"))
tab