## Global fire metrics and time since fire

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
library(ggpubr)
library(roquefort)
library(gt)

## Data: 
diversity <- readRDS("~/Imperial/PROJECT/PREDICTS data/database.rds")
fire_hist <- read.csv("~/Imperial/PROJECT/Combined data/Fire_history.csv")

## Data preparation: Same process and category aggregation as the model selection r script. 
diversity <- mutate(diversity, 
                    Measurement = Effort_corrected_measurement,
                    Sampling_effort = Rescaled_sampling_effort)
## Calculating diversity metrics
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
## Subset to 2006
sites_post_2006 <- subset(sites, sites$Sample_start_earliest > "2006-01-01")
## Removing rows that have NAs in latitude and longitude as these won't have fire data
sites_post_2006_NA <- subset(sites_post_2006, (! is.na(Longitude) | ! is.na(Latitude)))

## Merge with fire data
PREDICTS_fire_hist <- merge(sites_post_2006_NA, fire_hist
                            [,c("SSBS", "Freq", "fire_ID", "size", "start_date", "duration", "days_since_fire",
                                "year", "years_since_fire")], by="SSBS", all.x=TRUE)
PREDICTS_fire_hist$Freq[is.na(PREDICTS_fire_hist$Freq)] <- 0
PREDICTS_fire_hist$fire_presence <- ifelse(PREDICTS_fire_hist$Freq==0, "No",
                                           ifelse(PREDICTS_fire_hist$Freq!=0, "Yes",
                                                  NA))
PREDICTS_fire_hist$fire_presence <- as.factor(PREDICTS_fire_hist$fire_presence)

## Need to remove data which has NAs in categories required for our analysis
model_data <- drop_na(PREDICTS_fire_hist, Total_abundance, Predominant_land_use, Use_intensity, fire_presence)

table(model_data$Biome)

## Biome alteration/aggregation
table(model_data$Biome, model_data$fire_presence)
## remove the levels with no data/Can't have fires 
model_data<-model_data[!(model_data$Biome == "Rock & Ice"),]
model_data<-model_data[!(model_data$Biome == "Flooded Grasslands & Savannas"),]
model_data<-model_data[!(model_data$Biome == "Mangroves"),]
model_data<-model_data[!(model_data$Biome == "Inland Water"),]
model_data<-model_data[!(model_data$Biome == "Tundra"),]
model_data$Biome <- as.factor(as.character(model_data$Biome))
table(model_data$Biome)

## Forest/non forest
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

## Aggregated biomes
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

## Predominant land use and use intensity
table(model_data$Use_intensity, model_data$fire_presence)
table(model_data$Predominant_land_use, model_data$fire_presence)
table(model_data$Use_intensity, model_data$Predominant_land_use)
## Need to remove plantation and urban due to lack of data and merge intermediate and mature together
## Use intensity can stay as it is.

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

## Assumption checking:
## Collinearity
source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence")])

## complete cases
model_data <- drop_na(model_data, Total_abundance, Predominant_land_use, Use_intensity, fire_presence)

## Need to transform abundance
model_data <- mutate(model_data, 
                     logAbundance = log(RescaledAbundance + 1),
                     sqrtAbundance = sqrt(RescaledAbundance))

## basic fire presence model
m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
             Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
             Use_intensity:fire_presence +
             (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(m1)
summary(m1)

################################################################################################
## Visual inspection of metrics:
## Size
hist <- ggplot(data=model_data, aes(size)) + 
  geom_histogram(binwidth = 100, fill = 'light blue', colour = 'Black')+
  xlab("Fire size (km2)")+
  ylab("Number of sites")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

box <- ggplot(model_data, aes(x = Predominant_land_use, y = size)) + 
  geom_boxplot()+
  xlab("Predominant land use")+
  ylab("Size of fires")+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))
ggarrange(hist, box, ncol = 2, nrow = 1)

## Duration
hist <- ggplot(data=model_data, aes(duration)) + 
  geom_histogram(binwidth = 5, fill = 'light blue', colour = 'Black')+
  xlab("Duration of fire (days)")+
  ylab("Number of sites")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

box <- ggplot(model_data, aes(x = Predominant_land_use, y = duration)) + 
  geom_boxplot()+
  xlab("Predominant land use")+
  ylab("Duration of fire (days)")+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

ggarrange(hist, box, ncol = 2, nrow = 1)

## Frequency
graph_data <- subset(model_data, model_data$fire_presence == "Yes")
ggplot(data=graph_data, aes(Freq)) + 
  geom_histogram(binwidth = 1, fill = 'light blue', colour = 'Black')+
  xlab("Number of fires")+
  ylab("Number of sites")+
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

ggplot(model_data, aes(x = Predominant_land_use, y = Freq)) + 
  geom_boxplot()+
  xlab("Predominant land use")+
  ylab("Frequency of fires in a three year period")+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

## Days since fire
ggplot(data=model_data, aes(days_since_fire)) + 
  geom_histogram(binwidth = 100, fill = 'light blue', colour = 'Black')+
  xlab("Days since fire")+
  ylab("Number of sites")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

ggplot(model_data, aes(x = Predominant_land_use, y = days_since_fire)) + 
  geom_boxplot()+
  xlab("Predominant land use")+
  ylab("Days since fire")+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1), axis.title.x = element_text(size = 11), axis.text.y = element_text(size = 10, color = "black"))

###################################################################################################
## Models with the different fire metrics:

## Fire presence:
m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
             Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
             Use_intensity:fire_presence +
             (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(m1)

######
## Fire frequency:
str(model_data$Freq)
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "Freq")])
m_freq <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + Freq +
                 Predominant_land_use:Freq + Predominant_land_use:Use_intensity +Use_intensity:Freq +
                 (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(m_freq) #freq is significant
summary(m_freq)
theme_set(theme_classic())
frequency_plot <- plot_model(m_freq, type = "pred", terms = c("Freq", "Predominant_land_use"), 
                             axis.title = c("Fire frequency in a three year period", "Squareroot rescaled abundance"),
                             legend.title = c("Predominant land use"), title = "") +
  guides(col=guide_legend(nrow=5,byrow=TRUE))

frequency_UI_plot <- plot_model(m_freq, type = "pred", terms = c("Freq", "Use_intensity"), 
                                axis.title = c("Fire frequency in a three year period", "Squareroot rescaled abundance"),
                                legend.title = c("Land use intensity"), title = "")+
  guides(col=guide_legend(nrow=3,byrow=TRUE))

ggarrange(frequency_plot, frequency_UI_plot, nrow = 1, ncol = 2, legend = "bottom", 
          labels = c("A", "B"))

r.squaredGLMM(m_freq)

## Diagnostics for report
par(mfrow=c(1,2))
qqnorm(resid(m_freq), main="normal qq-plot, residuals")
qqline(resid(m_freq)) ## not bad
plot(fitted(m_freq), resid(m_freq), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)
ran <- ranef(m_freq, condVar=TRUE)
#double square bracket access the lists in d1
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
c <- dotplot(ran)[["Biome"]]
grid.arrange(a,b,c, nrow=2, ncol=2)

######
## Fire size:
model_data$size[is.na(model_data$size)] <- 0
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "size")])
m_size <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + size +
                 Predominant_land_use:size + Predominant_land_use:Use_intensity + Use_intensity:size +
                 (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(m_size) ## size significant
summary(m_size)
r.squaredGLMM(m_size)
theme_set(theme_classic())
size_plot <- plot_model(m_size, type = "pred", terms = c("size", "Predominant_land_use"), 
                        axis.title = c("Fire size (km2)", "Squareroot rescaled abundance"),
                        legend.title = c("Predominant land use"), title = "")+
  guides(col=guide_legend(nrow=5,byrow=TRUE))
size_UI_plot <- plot_model(m_size, type = "pred", terms = c("size", "Use_intensity"), 
                           axis.title = c("Fire size (km2)", "Squareroot rescaled abundance"),
                           legend.title = c("Land use intensity"), title = "")+
  guides(col=guide_legend(nrow=3,byrow=TRUE))
ggarrange(size_plot, size_UI_plot, nrow = 1, ncol = 2, legend = "bottom", 
          labels = c("A", "B"))
plot_model(m_size, type = "pred", terms = "size")

## Diagnostics
par(mfrow=c(1,2))
qqnorm(resid(m_size), main="normal qq-plot, residuals")
qqline(resid(m_size)) ## not bad
plot(fitted(m_size), resid(m_size), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

######
## Fire duration:
model_data$duration[is.na(model_data$duration)] <- 0
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "duration")])
m_duration <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + duration +
                     Predominant_land_use:duration + Predominant_land_use:Use_intensity +
                     Use_intensity:duration +
                     (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(m_duration)
summary(m_duration)
r.squaredGLMM(m_duration)
theme_set(theme_classic())
duration_plot <- plot_model(m_duration, type = "pred", terms = c("duration", "Predominant_land_use"), 
                            axis.title = c("Fire duration (days)", "Squareroot rescaled abundance"),
                            legend.title = c("Predominant land use"), title = "")+
  guides(col=guide_legend(nrow=5,byrow=TRUE))
duration_UI_plot <- plot_model(m_duration, type = "pred", terms = c("duration", "Use_intensity"), 
                               axis.title = c("Fire duration (days)", "Squareroot rescaled abundance"),
                               legend.title = c("Land use intensity"), title = "")+
  guides(col=guide_legend(nrow=3,byrow=TRUE))

ggarrange(duration_plot, duration_UI_plot, nrow = 1, ncol = 2, legend = "bottom", 
          labels = c("A", "B"))

par(mfrow=c(1,2))
qqnorm(resid(m_duration), main="normal qq-plot, residuals")
qqline(resid(m_duration)) ## not bad
plot(fitted(m_duration), resid(m_duration), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

## Comparing all these models 
m1.ml <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence +
                Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                Use_intensity:fire_presence +
                (1|SS) + (1|SSB) + (1|Biome), data = model_data, REML=FALSE)
m_freq.ml <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + Freq +
                    Predominant_land_use:Freq + Predominant_land_use:Use_intensity +Use_intensity:Freq +
                    (1|SS) + (1|SSB) + (1|Biome), data = model_data, REML=FALSE)
m_size.ml <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + size +
                    Predominant_land_use:size + Predominant_land_use:Use_intensity + Use_intensity:size +
                    (1|SS) + (1|SSB) + (1|Biome), data = model_data, REML=FALSE)
m_duration.ml <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + duration +
                        Predominant_land_use:duration + Predominant_land_use:Use_intensity +
                        Use_intensity:duration +
                        (1|SS) + (1|SSB) + (1|Biome), data = model_data, REML=FALSE)
AIC(m1.ml, m_freq.ml, m_size.ml, m_duration.ml)
## Best model with a single metric is basic fire presence/absence model

## New AIC table
summ.table <- do.call(rbind, lapply(list(m1.ml, m_freq.ml, m_size.ml, m_duration.ml), broom::glance))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid.Df", "Resid.Dev", "AIC")
reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
reported.table

## Combinations of metrics
## Trying models with combinations, It is important not to run into issues with collinearity
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence", "Freq", "size", "duration")])
## really high colinearity, so can't use all 
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence", "size", "duration")])
# v high GVIF for duration and fire presence
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence", "size")])
#ok
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence", "duration")])
# quite high collinearity here between fire presence and duration
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "Freq", "duration")])
# high but not as bad as before
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "Freq", "size")])
#ok
corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "duration", "size")])
# relatively high

# Combinations to model = fire presence and size, frequency and size
## Fire presence and size
m.pressize <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence + size +
                     Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                     Use_intensity:fire_presence + Predominant_land_use:size + Use_intensity:size +
                     (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(m.pressize) # Non of the size aspects come out as significant
summary(m.pressize)

## Fire frequency and size
m.freqsize <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + Freq + size +
                     Predominant_land_use:Freq + Predominant_land_use:Use_intensity +
                     Use_intensity:Freq + Predominant_land_use:size + Use_intensity:size +
                     (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(m.freqsize)
summary(m.freqsize)

## Comparing all models:
m.pressize.ml <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + fire_presence + size +
                        Predominant_land_use:fire_presence + Predominant_land_use:Use_intensity +
                        Use_intensity:fire_presence + Predominant_land_use:size + Use_intensity:size +
                        (1|SS) + (1|SSB) + (1|Biome), data = model_data, REML=FALSE)
m.freqsize.ml <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + Freq + size +
                        Predominant_land_use:Freq + Predominant_land_use:Use_intensity +
                        Use_intensity:Freq + Predominant_land_use:size + Use_intensity:size +
                        (1|SS) + (1|SSB) + (1|Biome), data = model_data, REML=FALSE)

AIC(m1.ml, m_freq.ml, m_size.ml, m_duration.ml, m.pressize.ml, m.freqsize.ml)
## Best model is m1 so the model that only contains fire presence. 
summ.table <- do.call(rbind, lapply(list(m1.ml, m_freq.ml, m_size.ml, m_duration.ml, m.pressize.ml, m.freqsize.ml), broom::glance))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid.Df", "Resid.Dev", "AIC")
reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
reported.table

#################################################################################################
## Time since fire:
str(model_data$days_since_fire) # Need to make it numeric
model_data$days_since_fire <- as.numeric(model_data$days_since_fire)

## Need to change the NAs to a value in order to model them. Going to use 3.5 years or 1277 days as
#an arbitary value. 
model_data$days_since_fire[is.na(model_data$days_since_fire)] <- 1277

## for ease of understanding it may be best to convert days since fire into years.
model_data$time_since_fire <- model_data$days_since_fire
model_data$time_since_fire <- model_data$time_since_fire/365
str(model_data$time_since_fire)

corvif(model_data[ , c("Predominant_land_use", "Use_intensity", "fire_presence", "time_since_fire")])
# High collinearity in fire presence and days since fire

m1_time <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + time_since_fire + 
                  Predominant_land_use:Use_intensity + Predominant_land_use:time_since_fire + 
                  Use_intensity:time_since_fire + (1|SS) + (1|SSB) + (1|Biome), data = model_data)
Anova(m1_time)
summary(m1_time)

## Diagnostics:
par(mfrow=c(1,2))
qqnorm(resid(m1_time), main="normal qq-plot, residuals")
qqline(resid(m1_time)) ## not bad
plot(fitted(m1_time), resid(m1_time), main = "Residuals vs fitted") #residuals vs fitted
abline(h=0)

ran <- ranef(m1_time, condVar=TRUE)
#double square bracket access the lists in d1
a <- dotplot(ran)[["SSB"]]
b <- dotplot(ran)[["SS"]]
c <- dotplot(ran)[["Biome"]]
grid.arrange(a,b,c, nrow=2, ncol=2) 

table <- tab_model(m1_time, show.stat = TRUE, show.se = TRUE, 
                   show.df = TRUE, string.stat = c("t statistic"))
r.squaredGLMM(m1_time)

## Hypothesis testing:
m1_notime <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity +  
                    Predominant_land_use:Use_intensity + 
                    (1|SS) + (1|SSB) + (1|Biome), data = model_data)
anova(m1_time, m1_notime)
## chisq = 36.256, df=7. p<0.001

## Plots
theme_set(theme_classic())
a <- plot_model(m1_time, type = "pred", terms = c("time_since_fire", "Predominant_land_use"), 
                legend.title = "Predominant land use", axis.title = c("Time since fire (years)", 
                                                                      "Square root rescaled abundance"),
                title = "")+
  font_size(axis_title.x = 10, axis_title.y = 10)
b <- plot_model(m1_time, type = "pred", terms = c("time_since_fire", "Use_intensity"), 
                legend.title = "Use intensity", axis.title = c("Time since fire (years)", 
                                                               "Square root rescaled abundance"),
                title = "")+
  font_size(axis_title.x = 10, axis_title.y = 10)

ggarrange(a, b, nrow = 2, ncol=1, legend = "right", labels = c("A", "B"))

### Checking linearity:
## Investigating if the relationship is truely linear
plot(m1_time) 
## Clear bounding occuring
## checking diagnostics
time_output <- simulateResiduals(fittedModel = m1_time)
plot(time_output, quantreg = T)
plotResiduals(time_output)
## This plot is not great so we can try plotting the specific predictor
plotResiduals(time_output, model_data$time_since_fire)
# Plot implies there are no real issues happening here so the relationship is likely to  be linear

## -1/time
model_data$minustime <- model_data$time_since_fire
model_data$minustime <- -1/model_data$time_since_fire
m1_minustime <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + minustime + 
                       Predominant_land_use:Use_intensity + Predominant_land_use:minustime + 
                       Use_intensity:minustime + (1|SS) + (1|SSB) + (1|Biome), data = model_data)
plot(m1_minustime)
minustime_output <- simulateResiduals(fittedModel = m1_minustime)
plot(minustime_output, quantreg = T) ## More issues in this plot
plotResiduals(minustime_output, model_data$minustime)

## Log time
model_data$logtime <- model_data$time_since_fire
model_data$logtime <- log(model_data$time_since_fire)
m1_logtime <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + logtime + 
                     Predominant_land_use:Use_intensity + Predominant_land_use:logtime + 
                     Use_intensity:logtime + (1|SS) + (1|SSB) + (1|Biome), data = model_data)
plot(m1_logtime)
logtime_output <- simulateResiduals(fittedModel = m1_logtime)
plot(logtime_output, quantreg = T) ## More issues in this plot
plotResiduals(logtime_output, model_data$logtime)

## Transformation don't improve time therefore can assume that the relationship is linear. 

## Comparing to other models
m1_time.ml <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + time_since_fire + 
                     Predominant_land_use:Use_intensity + Predominant_land_use:time_since_fire + 
                     Use_intensity:time_since_fire + (1|SS) + (1|SSB) + (1|Biome), 
                   data = model_data, REML = FALSE)
AIC(m1.ml, m_freq.ml, m_size.ml, m_duration.ml, m.pressize.ml, m.freqsize.ml, m1_time.ml)
## m1 remains the best model. 