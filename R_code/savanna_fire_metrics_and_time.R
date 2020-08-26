## Savanna abundance models with fire metrics and time since fire. 

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
library(gridExtra)
library(ggplot2)
library(roquefort)
library(ggpubr)

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

## Looking at the data
## We need to make sure we have sufficient data for analysis
table(sites_savanna$Predominant_land_use, sites_savanna$fire_presence)
# Removing land use level which lack data
sites_savanna<-sites_savanna[!(sites_savanna$Predominant_land_use == "Plantation forest"),]
sites_savanna<-sites_savanna[!(sites_savanna$Predominant_land_use == "Urban"),] 
sites_savanna$Predominant_land_use <- as.factor(as.character(sites_savanna$Predominant_land_use))
## Merge secondary vegetation categories together
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

## Ensuring complete cases and normality in the response variable. 
model_data <- drop_na(sites_savanna, Total_abundance, Predominant_land_use, Use_intensity_2, fire_presence)
model_data <- mutate(model_data, 
                     logAbundance = log(RescaledAbundance + 1),
                     sqrtAbundance = sqrt(RescaledAbundance))

###################################################################################################
## fire presence model
m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence +
             Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2+
             (1|SS) + (1|SSB), data = model_data)

Anova(m1)
summary(m1)

## Model selection for the basic model is in the savanna model selection r file. 
## The same fixed and random effects structures will be used here. 

## Looking at our years and days since fire:
days_hist <- ggplot(data=model_data, aes(days_since_fire)) + 
  geom_histogram(binwidth = 100, fill = 'light blue', colour = 'Black')+
  xlab("Days since fire")+
  ylab("Number of sites")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

years_hist <- ggplot(data=model_data, aes(years_since_fire)) + 
  geom_histogram(binwidth = 1, fill = 'light blue', colour = 'Black')+
  xlab("Years since fire")+
  ylab("Number of sites")+
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))
grid.arrange(days_hist, years_hist, nrow=2)

## Fire metrics:

# Frequency of fires
ggplot(data=model_data, aes(Freq)) + 
  geom_histogram(binwidth = 1, fill = 'light blue', colour = 'Black')+
  xlab("Number of fires")+
  ylab("Number of sites")+
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

# Fire size
ggplot(data=model_data, aes(size)) + 
  geom_histogram(binwidth = 100, fill = 'light blue', colour = 'Black')+
  xlab("Fire size (km2)")+
  ylab("Number of sites")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

# Fire duration
gplot(data=model_data, aes(duration)) + 
  geom_histogram(binwidth = 5, fill = 'light blue', colour = 'Black')+
  xlab("Duration of fire (days)")+
  ylab("Number of sites")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

## Looking at using different years of data - 3yr window vs 2yr vs 1yr
#ie any fires within 1 yr of sampling, any fires within 2 years and any within 1 (so I think we can 
#just subset this from the years since fire column)

model_data_2yr <- subset(model_data, years_since_fire < 3 | is.na(years_since_fire))
model_data_1yr <- subset(model_data, years_since_fire < 2 | is.na(years_since_fire))

## Now we can try the model with this data
m_2yr <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence +
                Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2+
                (1|SS) + (1|SSB), data = model_data_2yr)
m_1yr <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence +
                Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2+
                (1|SS) + (1|SSB), data = model_data_1yr)

Anova(m1)
Anova(m_2yr)
Anova(m_1yr)
## The 1yr data has a higher significant of fire in the model (interesting maybe?- possibly because 
#the sampling would have been done closer to the disturbance and so the vegetation has had a chance 
#to recover yet and so will be in that initially damaged state still) The 2yr data doesn't suggest 
#that this is a linear trend. 
## The effect of land use also appears to become less important as the sampling window is shortened. 

summary(m_1yr)
summary(m_2yr)
## It appears that the main cause of the increased signficant of fire is possibly the older secondary 
#vegetation. In the 2yr model estimate is 0.088 (0.35 T value), whereas in the 1yr model it is -0.41
#(-2.886 t value). THis suggests a rapid recovery of secondary vegetation in the space of a year. 

## I think this supports the idea of trying to look at how the effect of fire in different land use 
#categories changes with time since the fire. 

####################################################################################################
## Modelling with different fire metrics...'
## M1 is just the original model with fire presence in it
m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence +
             Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2+
             (1|SS) + (1|SSB), data = model_data)
m1.comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence +
                  Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity_2+
                  (1|SS) + (1|SSB), data = model_data, REML = FALSE)

## using fire frequency
str(sav_ab_model_data$Freq) ## Make sure frequency is numeric variable
m1.freq <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + Freq +
                  Predominant_land_use:Freq + Predominant_land_use:Use_intensity_2+
                  (1|SS) + (1|SSB), data = model_data)
Anova(m1.freq)
summary(m1.freq)

m1.freq.comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + Freq +
                       Predominant_land_use:Freq + Predominant_land_use:Use_intensity_2+
                       (1|SS) + (1|SSB), data = model_data, REML=FALSE)
## Get a convergence warning in this model? confused as it doesn't happen in REML=TRUE version.


## To look at size and duration we need to recode the factors so that we are dealing with the same 
#number of rows of data.
model_data$size[is.na(model_data$size)] <- 0
model_data$duration[is.na(model_data$duration)] <- 0


## using fire size
m1.size <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + size +
                  Predominant_land_use:size + Predominant_land_use:Use_intensity_2+
                  (1|SS) + (1|SSB), data = model_data)

Anova(m1.size)
summary(m1.size)

## using fire duration
m1.duration <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + duration +
                      Predominant_land_use:duration + Predominant_land_use:Use_intensity_2+
                      (1|SS) + (1|SSB), data = model_data)

Anova(m1.duration)
summary(m1.duration)

## Compare all models
AIC(m1.comp, m1.freq.comp, m1.duration.comp, m1.size.comp)
## all very similar AIC values not much in it really as to which is a better predictor if any. 
## implies that duration is the better model? 

## Trying modelling with multiple metrics

## Checking for collinearity
source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "fire_presence", "Freq", "size", "duration")])
## really high colinearity 
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "fire_presence", "size", "duration")])
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "fire_presence", "size")])
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "fire_presence", "duration")])
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "Freq", "duration")])
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "Freq", "size")])
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "duration", "size")])
## fire presence and size don't have much collinearity 


## model with fire presence and size
m1.sizpres <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence + size +
                     Predominant_land_use:fire_presence + Predominant_land_use:size +
                     Predominant_land_use:Use_intensity_2+
                     (1|SS) + (1|SSB), data = model_data)
Anova(m1.sizpres)
summary(m1.sizpres)

m1.sizpres.comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence + size +
                          Predominant_land_use:fire_presence + Predominant_land_use:size +
                          Predominant_land_use:Use_intensity_2+
                          (1|SS) + (1|SSB), data = model_data, REML=FALSE)

## model with fire presence and duration
m1.durpres <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence + duration +
                     Predominant_land_use:fire_presence + Predominant_land_use:duration +
                     Predominant_land_use:Use_intensity_2+
                     (1|SS) + (1|SSB), data = model_data)
Anova(m1.durpres) 
summary(m1.durpres)

m1.durpres.comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + fire_presence + duration +
                          Predominant_land_use:fire_presence + Predominant_land_use:duration +
                          Predominant_land_use:Use_intensity_2+
                          (1|SS) + (1|SSB), data = model_data, REML=FALSE)

## model with fire frequency and size
m1.sizfreq <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + Freq + size +
                     Predominant_land_use:Freq + Predominant_land_use:size +
                     Predominant_land_use:Use_intensity_2+
                     (1|SS) + (1|SSB), data = model_data)
Anova(m1.sizfreq)
summary(m1.sizfreq)

m1.sizfreq.comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + Freq + size +
                          Predominant_land_use:Freq + Predominant_land_use:size +
                          Predominant_land_use:Use_intensity_2+
                          (1|SS) + (1|SSB), data = model_data, REML=FALSE)

## model with fire frequency and duration
m1.durfreq <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + Freq + duration +
                     Predominant_land_use:Freq + Predominant_land_use:duration +
                     Predominant_land_use:Use_intensity_2+
                     (1|SS) + (1|SSB), data = model_data)
Anova(m1.durfreq)
summary(m1.durfreq)

m1.durfreq.comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + Freq + duration +
                          Predominant_land_use:Freq + Predominant_land_use:duration +
                          Predominant_land_use:Use_intensity_2+
                          (1|SS) + (1|SSB), data = model_data, REML=FALSE)


## size and duration

m1.sizdur <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + duration + size +
                    Predominant_land_use:duration + Predominant_land_use:size + duration:size +
                    Predominant_land_use:Use_intensity_2+
                    (1|SS) + (1|SSB), data = model_data)

Anova(m1.sizdur)
summary(m1.sizdur)

m1.sizdur.comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + duration + size +
                         Predominant_land_use:duration + Predominant_land_use:size + duration:size +
                         Predominant_land_use:Use_intensity_2+
                         (1|SS) + (1|SSB), data = model_data, REML = FALSE)

######### COMPARING ALL MODELS 

AIC(m1.comp, m1.freq.comp, m1.duration.comp, m1.size.comp, m1.sizpres.comp, m1.sizfreq.comp, 
    m1.durpres.comp, m1.durfreq.comp, m1.sizdur.comp)

###################################################################################################
## DAYS SINCE FIRE

## Linear models using days since fire... Issue as it contains NAs for values where we don't know if
#there has been a fire. can we code this as 0 or at least x number of days?
str(model_data$days_since_fire) # not numeric
model_data$days_since_fire <- as.numeric(model_data$days_since_fire)
str(model_data$days_since_fire)

## for now will just change the na's to be 1277 (3.5 years)
model_data$time_since_fire <- model_data$days_since_fire
model_data$time_since_fire[is.na(model_data$time_since_fire)] <- 1277
str(model_data$time_since_fire)

corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "fire_presence", "time_since_fire")])
## Very high correlation with fire presence and time need to just include one
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "time_since_fire")])

## Days model
m1_days <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + time_since_fire +
                  Predominant_land_use:time_since_fire +
                  Predominant_land_use:Use_intensity_2+
                  (1|SS) + (1|SSB), data = model_data)
Anova(m1_days)
summary(m1_days)

## May aid in understanding the output to use years since fire. 

model_data$years_since_int <- model_data$time_since_fire
model_data$years_since_int <- model_data$time_since_fire/365

corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "fire_presence", "years_since_int")])
## Very high correlation with fire presence and time need to just include one
corvif(model_data[ , c("Predominant_land_use", "Use_intensity_2", "years_since_int")])

m1_years <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + years_since_int +
                   Predominant_land_use:years_since_int+
                   Predominant_land_use:Use_intensity_2+
                   (1|SS) + (1|SSB), data = model_data, REML=TRUE)
Anova(m1_years)
summary(m1_years)

## Model diagnostics:
par(mfrow=c(1,2))
qqnorm(resid(m1_years), main="normal qq-plot, residuals")
qqline(resid(m1_years)) ## not bad
plot(fitted(m1_years), resid(m1_years), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)
ran <- ranef(m1_years, condVar=TRUE)
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
grid.arrange(a,b, nrow=2, ncol=1) 

## R2:
r.squaredGLMM(m1_years)

## Plotting and looking at time since fire model - m1_years
set_theme(base = theme_classic())
Time_plot <- plot_model(m1_years, type = "pred", terms = c("years_since_int", "Predominant_land_use"), show.intercept = TRUE, 
                        axis.title = c("Time since fire (years))", "Square root of rescaled abundance"),
                        legend.title = c("Predominant land use")) 
Time_plot

### Testing if this relationship is linear:

## Checking for misspecifications in time model:
years_output <- simulateResiduals(fittedModel = m1_years)
plot(years_output, quantreg = T)
## Now want to plot residuals against the predictors
dev.off()
par(mfrow = c(1,2))
plotResiduals(years_output1, model_data$Predominant_land_use)
plotResiduals(years_output1, model_data$years_since_int)

## Transforming time:
## -1/time 
model_data$time_transformed <- model_data$years_since_int
model_data$time_transformed <- -1/model_data$time_transformed

m1_years3 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + time_transformed +
                    Predominant_land_use:time_transformed + 
                    Predominant_land_use:Use_intensity_2 +
                    (1|SS) + (1|SSB), data = model_data)
Anova(m1_years3)
summary(m1_years3)
years_output2 <- simulateResiduals(fittedModel = m1_years3)
plot(years_output2, quantreg = T)
par(mfrow = c(1,2))
plotResiduals(years_output2, model_data$Predominant_land_use)
plotResiduals(years_output2, model_data$time_transformed)
## This transformation hasn't improved the model 

##log(time)
model_data$logtime <- model_data$years_since_int
model_data$logtime <- log(model_data$logtime)

m1_years4 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + logtime +
                    Predominant_land_use:logtime+
                    Predominant_land_use:Use_intensity_2 +
                    (1|SS) + (1|SSB), data = model_data)
Anova(m1_years4)
summary(m1_years4)
years_output3 <- simulateResiduals(fittedModel = m1_years4)
plot(years_output3, quantreg = T)
par(mfrow = c(1,2))
plotResiduals(years_output, model_data$Predominant_land_use)
plotResiduals(years_output3, model_data$logtime)

## No indication that the relationship is not linear

################################################################################################
## Comparing models:
m1_years.comp <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + years_since_int +
                        Predominant_land_use:years_since_int+
                        Predominant_land_use:Use_intensity_2+
                        (1|SS) + (1|SSB), data = model_data, REML=FALSE)

AIC(m1.comp, m1.freq.comp, m1.duration.comp, m1.size.comp, m1.sizpres.comp, m1.sizfreq.comp, 
    m1.durpres.comp, m1.durfreq.comp, m1.sizdur.comp, m1_years.comp)

## Best model is m1 so the model that only contains fire presence. 
summ.table <- do.call(rbind, lapply(list(m1.comp, m1.freq.comp, m1.duration.comp, m1.size.comp, m1.sizpres.comp, m1.sizfreq.comp, 
                                         m1.durpres.comp, m1.durfreq.comp, m1.sizdur.comp, m1_years.comp), broom::glance))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid.Df", "Resid.Dev", "AIC")
reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
model.names <- c("Fire presence", "Frequency", "Duration", "Size", "Presence and size", 
                 "Frequency and size", "Presence and duration", "Frequency and duration",
                 "Size and duration", "Time since fire")
row.names(reported.table) <- model.names
reported.table

##############################################################################################
## Looking at the duration model in more detail:

m1.duration <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + duration +
                      Predominant_land_use:duration + Predominant_land_use:Use_intensity_2+
                      (1|SS) + (1|SSB), data = model_data)

Anova(m1.duration)
# we have a significant interaction effect between land use and use intensity
summary(m1.duration)
r.squaredGLMM(m1.duration)

## Diagnostics
par(mfrow=c(1,2))
qqnorm(resid(m1.duration), main="normal qq-plot, residuals")
qqline(resid(m1.duration)) ## not bad
plot(fitted(m1.duration), resid(m1.duration), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)
ran <- ranef(m1.duration, condVar=TRUE)
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
grid.arrange(a,b, nrow=2, ncol=1)

## Model criticism/diagnostics
duration_output <- simulateResiduals(fittedModel = m1.duration, n = 250)
plot(duration_output) ## Plot not great. 
## Looking at each predictor seperately

par(mfrow=c(1,3))
plotResiduals(duration_output, model_data$Predominant_land_use, main="Predominant land use")
plotResiduals(duration_output, model_data$Use_intensity_2, main="Use intensity")
plotResiduals(duration_output, model_data$duration, main="Duration")
## Some issues with duration in this model. 

## Residuals
residuals <- resid(m1.duration)
summary(residuals)
hist(residuals) ## Residuals look ok

## Hypothesis testing:
m1.duration.2 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity_2 + 
                        Predominant_land_use:Use_intensity_2+
                        (1|SS) + (1|SSB), data = model_data)
anova(m1.duration, m1.duration.2)