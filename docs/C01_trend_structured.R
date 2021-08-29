##########################################################################
#### Workshop CSD-Toolbox, 29.08.2021                                 ####
#### Lionel Hertzog, Thünen Institute for Biodiversity                ####
#### Fitting spatial and temporal GLMs to classical data              ####
#### using lme4 and DHARMa for checking                               ####
##########################################################################

## set wd where the script and the data are located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load libraries
library(tidyverse)
library(lme4)
library(DHARMa)

## load temporal data
temp_dat <- read.csv("C01_temporal_data.csv")

## plot the data
temp_dat$int_year <- as.numeric(gsub("Year", "", temp_dat$Year))

ggplot(temp_dat, aes(x=int_year, y=Occupied)) +
  geom_jitter(height = 0.1, width = 0.1) +
  scale_x_continuous(breaks = 1:10) +
  theme_classic()
  
# count number of occupied plots per year
temp_dat %>%
  group_by(Year, int_year) %>%
  summarise(Occupied = sum(Occupied)) -> temp_year

ggplot(temp_year, aes(x=int_year, y=Occupied)) +
  geom_point() +
  scale_x_continuous(breaks = 1:10) +
  theme_classic()

## fit a binomial model
m_temp <- glmer(Occupied ~ Year + (1 | Site), temp_dat,
              family = "binomial")  

## model summary
summary(m_temp)

## look at residuals
resid_temp <- simulateResiduals(m_temp, integerResponse = TRUE)
plot(resid_temp)

## derive model predictions,
## see http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions
temp_year$p_occ_link <- predict(m_temp, newdata = temp_year, re.form = ~ 0, type = "link")

## derive CIs around the predictions
modmat <- model.matrix(~Year, temp_year)
pvar1 <- diag(modmat %*% tcrossprod(vcov(m_temp),modmat))
temp_year$lci_link <- temp_year$p_occ_link - 1.96 * sqrt(pvar1)
temp_year$uci_link <- temp_year$p_occ_link + 1.96 * sqrt(pvar1)
## now transform back to original scale
## and put into expected number of occupied sites
temp_year$pred_occ <- 100 * arm::invlogit(temp_year$p_occ_link)
temp_year$lci_pred_occ <- 100 * arm::invlogit(temp_year$lci_link)
temp_year$uci_pred_occ <- 100 * arm::invlogit(temp_year$uci_link)

## plot this
ggplot(temp_year, aes(x = int_year)) +
  geom_point(aes(y = Occupied), color = "blue",
             shape = 4, size = 2) +
  geom_linerange(aes(ymin = lci_pred_occ, 
                     ymax = uci_pred_occ),
                 color = "red") +
  geom_point(aes(y = pred_occ), color = "red") +
  scale_x_continuous(breaks = 1:10) +
  theme_classic()

## load spatial data
spat_dat <- read.csv("C01_spatial_data.csv")
spat_pred <- read.csv("C01_spatial_data_prediction.csv") # contains data for prediction

## look at the data
summary(spat_dat)

## make a plot
ggplot(spat_dat, aes(x = elev, y = obs)) +
  geom_point() 


## fit the model
m_spat <- glm(obs ~ elev + elev2, data = spat_dat,
         family = "binomial")

## model summary
summary(m_spat)

## look at residuals
dd <- DHARMa::simulateResiduals(m_spat, integerResponse = TRUE)
plot(dd)

## derive model prediction over the whole domain
spat_pred$p_occ <- predict(m_spat, newdata = spat_pred, type = "response")

# plot the predicted spatial distribution
ggplot(spat_pred, aes(x = x, y = y, fill = p_occ)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) +
  labs(title = "Estimated probability\nof occurence (glm)")


