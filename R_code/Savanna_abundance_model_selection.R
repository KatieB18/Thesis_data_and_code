## Choosing a savanna model to use.
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
library(sjmisc)
library(sjlabelled)
library(roquefort)


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
PREDICTS_fire_hist2 <- drop_na(PREDICTS_fire_hist, Total_abundance, Predominant_land_use, Use_intensity, fire_presence)

## Subset to just the Tropical & Subtropical Grasslands, Savannas & Shrublands
sites_savanna <- subset(PREDICTS_fire_hist2, PREDICTS_fire_hist2$Biome == "Tropical & Subtropical Grasslands, Savannas & Shrublands")

## Looking at the data to make sure we have sufficient data for analysis
table(sites_savanna$Predominant_land_use, sites_savanna$fire_presence)
# Removing land use level which lack data
sites_savanna<-sites_savanna[!(sites_savanna$Predominant_land_use == "Plantation forest"),]
sites_savanna<-sites_savanna[!(sites_savanna$Predominant_land_use == "Urban"),] 
sites_savanna$Predominant_land_use <- as.factor(as.character(sites_savanna$Predominant_land_use))
## Merge together intermediate and mature secondary vegetation into one class
sites_savanna <- sites_savanna %>%
  mutate( Predominant_land_use = recode_factor(Predominant_land_use, 
                                               "Intermediate secondary vegetation" = "Older secondary vegetation", 
                                               "Mature secondary vegetation" = "Older secondary vegetation"),
          Predominant_land_use = factor(Predominant_land_use),
          Predominant_land_use = relevel(Predominant_land_use, ref = "Primary vegetation")
  )
sites_savanna$Predominant_land_use <- factor(sites_savanna$Predominant_land_use, levels = c("Primary vegetation", "Young secondary vegetation",
                                                                                            "Older secondary vegetation", "Cropland", "Pasture"))
table(sites_savanna$Predominant_land_use, sites_savanna$fire_presence)

## May also be best to merge the two age classes of secondary vegetation together due to low sample 
#sizes, especially when combined with use intensity. 

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

## Creating a second land use classs. 
sites_savanna$Landuse <- sites_savanna$Predominant_land_use
sites_savanna <- sites_savanna %>%
  mutate(Landuse = recode_factor(Landuse, "Young secondary vegetation" = "Secondary vegetation",
                                 "Older secondary vegetation" = "Secondary vegetation"),
         Landuse = factor(Landuse),
         Landuse = relevel(Landuse, ref = "Primary vegetation"))
sites_savanna$Landuse <- factor(sites_savanna$Landuse, levels = c("Primary vegetation",
                                                                  "Secondary vegetation",
                                                                  "Cropland",
                                                                  "Pasture"))
table(sites_savanna$Landuse, sites_savanna$fire_presence)
table(sites_savanna$Landuse, sites_savanna$Use_intensity_2)

#####################################################################################################
## Assumption checking before modelling

## Checking collinearity 
source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(sites_savanna[ , c("Predominant_land_use", "Use_intensity", "fire_presence")])
corvif(sites_savanna[ , c("Predominant_land_use", "Use_intensity_2", "fire_presence")])
corvif(sites_savanna[ , c("Landuse", "Use_intensity_2", "fire_presence")])
## GVIF values are all quite low. Using the recoded land use also lowers the GVIF values 

## complete cases have all ready been checked earlier on too
model_data <- drop_na(sites_savanna, Total_abundance, Predominant_land_use, Landuse, Use_intensity_2, fire_presence)

table(model_data$fire_presence)
## Need to transform abundance due to non normality
par(mfrow=c(1,2))
hist(model_data$RescaledAbundance) ## At the moment it shows a skewed distribution 
qqnorm(model_data$RescaledAbundance) ## curve (skewed data) and long tails
model_data <- mutate(model_data, 
                     logAbundance = log(RescaledAbundance + 1),
                     sqrtAbundance = sqrt(RescaledAbundance))
dev.off()
par(mfrow=c(2,2))
hist(model_data$logAbundance) ## still shows a skew
qqnorm(model_data$logAbundance) ## shows quite strong tails on either end also a slight curve indicating
#skewed data. 
hist(model_data$sqrtAbundance) ## Shows a normal distribution
qqnorm(model_data$sqrtAbundance) ## also shows quite distinct tails
## Use the sqrt abundance

###################################################################################################
## Model selection

## Choosing the maximal/starting model 
m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence +
             Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2 +
             (1|SS) + (1|SSB), data = model_data)

## Using recoded or original use intensity
m1.F <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence +
               Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2 +
               (1|SS) + (1|SSB), data = model_data, REML = FALSE)
m1.F2 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                (1|SS) + (1|SSB), data = model_data, REML = FALSE)
## Then test using AIC
AIC(m1.F, m1.F2)
## the model which uses the recoded use intensity has the lower AIC value therefore we shall use that

## Choosing whether to use Predominant_land_use or Landuse 
## Predominant land use model = m1.F
m1.F3 <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                Landuse:fire_presence + Landuse:Use_intensity_2 +
                (1|SS) + (1|SSB), data = model_data, REML = FALSE)
AIC(m1.F, m1.F3)
## Lowest AIC is m1.F3 therefore we should use the model with the recoded land use categories. 

## Choosing the random effects structure:
m1.re1 <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                 Landuse:fire_presence + Landuse:Use_intensity_2 +
                 (1|SS) + (1|SSB), data = model_data)
m1.re2 <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                 Landuse:fire_presence + Landuse:Use_intensity_2 +
                 (1+Landuse|SS) + (1|SSB), data = model_data)
## Error message warning of singular fit
m1.re3 <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                 Landuse:fire_presence + Landuse:Use_intensity_2 +
                 (1+Use_intensity_2|SS) + (1|SSB), data = model_data)
m1.re4 <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                 Landuse:fire_presence + Landuse:Use_intensity_2 +
                 (1+fire_presence|SS) + (1|SSB), data = model_data)

AIC(m1.re1, m1.re2, m1.re3, m1.re4)
## Model with the lowest AIC is m1.re1 therefore we shall use this random effects structure
## Better AIC table 
summ.table <- do.call(rbind, lapply(list(m1.re1, m1.re2, m1.re3, m1.re4), broom::glance))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid.Df", "Resid.Dev", "AIC")
reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
model.names <- c("Intercepts only", "Land use slope", "Use intensity slope", "fire presence slope")
row.names(reported.table) <- model.names
View(reported.table)

## Choosing the fixed effects

## Backwards stepwise selection:
max.ms <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                 Landuse:fire_presence + fire_presence:Use_intensity_2 + 
                 Landuse:Use_intensity_2 +
                 (1|SS) + (1|SSB), data = model_data)
Anova(max.ms)
## Can remove the use intensity fire presence interaction as it is the least significant
## But because we are particulary interested in using fire in these models going to see which model
## that included fire performs best. 

## Creating models with different combinations of fixed effects. To compare different fixed effect
#we need REML = false
m1.ms <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                Landuse:fire_presence + Landuse:Use_intensity_2 +
                (1|SS) + (1|SSB), data = model_data, REML = FALSE)
m2.ms <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                Landuse:fire_presence +
                (1|SS) + (1|SSB), data = model_data, REML = FALSE)
m3.ms <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                Landuse:fire_presence + Use_intensity_2:fire_presence +
                (1|SS) + (1|SSB), data = model_data, REML = FALSE)
m4.ms <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
                Landuse:fire_presence + fire_presence:Use_intensity_2 + 
                Landuse:Use_intensity_2 +
                (1|SS) + (1|SSB), data = model_data)

m.null <- lmer(sqrtAbundance ~ 
                 (1|SS) + (1|SSB), data = model_data, REML = FALSE)

AIC(m1.ms, m2.ms, m3.ms, m4.ms, m.null)
#m1.ms comes out with the lowest AIC. 

summ.table2 <- do.call(rbind, lapply(list(m1.ms, m2.ms, m3.ms, m4.ms, m.null), broom::glance))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table2 <- summ.table2[table.cols]
names(reported.table2) <- c("Resid.Df", "Resid.Dev", "AIC")
reported.table2[['dAIC']] <-  with(reported.table2, AIC - min(AIC))
reported.table2[['weight']] <- with(reported.table2, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table2$AIC <- NULL
reported.table2$weight <- round(reported.table2$weight, 2)
reported.table2$dAIC <- round(reported.table2$dAIC, 1)
model.names2 <- c("m1", "m2", "m3", "m4", "null model")
row.names(reported.table2) <- model.names2
reported.table2 ## Based off this table we should use model 1!

###################################################################################################
## Model to use = 
m1 <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
             Landuse:fire_presence + Landuse:Use_intensity_2+
             (1|SS) + (1|SSB), data = model_data)

Anova(m1)
summary(m1)

## Model diagnostics:
m1_output <- simulateResiduals(fittedModel = m1, n = 250)
plot(m1_output) 
## Looking at each predictor seperately
par(mfrow=c(1,3))
plotResiduals(m1_output, model_data$Landuse, main="Predominant land use")
plotResiduals(m1_output, model_data$Use_intensity_2, main="Use intensity")
plotResiduals(m1_output, model_data$fire_presence, main="fire presence")

## Residuals
residuals <- resid(m1)
summary(residuals) 
hist(residuals) ## Residuals look ok

par(mfrow=c(1,2))
qqnorm(resid(m1), main="normal qq-plot, residuals")
qqline(resid(m1)) ## not bad
plot(fitted(m1), resid(m1), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)
ran <- ranef(m1, condVar=TRUE)
#double square bracket access the lists in d1
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
grid.arrange(a,b, nrow=2, ncol=1)

## Using the MuMIn package to get r2 values
r.squaredGLMM(m1)

## Table of results:
Savanna_abundance_table <- tab_model(m1, show.stat = TRUE, show.se = TRUE, show.df = TRUE, 
                                     string.stat = c("t statistic"))
Savanna_abundance_table

## Plots:
plot_model(m1, type = "pred", terms = c("Landuse", "fire_presence"), legend.title = "Fire presence",
           show.intercept = TRUE, axis.title = c("Predominant land use", "Square root of rescaled abundance")) 

## Hypothesis testing with and fire out fire presence
fire <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + fire_presence +
               Landuse:fire_presence + Landuse:Use_intensity_2+
               (1|SS) + (1|SSB), data = model_data, REML = FALSE)
no_fire <- lmer(sqrtAbundance ~ Landuse + Use_intensity_2 + Landuse:Use_intensity_2+
                  (1|SS) + (1|SSB), data = model_data, REML = FALSE)

anova(fire, no_fire) #anova with maximum liklihood = chisq=2.4947, df=4, p=0.6456