## Coding session 4: Spatial analysis

## set wd
# DEV !!!
setwd("~/PostDoc_Thunen/Stat-stuff/Gfoe_workshop/")

## load libraries
library(tidyverse)
library(INLA)
library(inlabru)
library(DHARMa)

## load georeferenced data from C01
dat <- readRDS("C01_spatial_data.rds")
dat_spat <- dat$data
dat_pred <- dat$prediction

## create the mesh
mesh <- inla.mesh.2d(loc = dat_spat[,c("x", "y")],
                     max.edge = c(0.1, 1))

## plot the mesh
plot(mesh)

## create the SPDE
spat <- inla.spde2.pcmatern(mesh, 
                               prior.range = c(1, 0.05),
                               prior.sigma = c(5, 0.05))

## turn the data into a SPDF
dat_spat_spdf <- dat_spat
coordinates(dat_spat_spdf) <- c("x", "y")

## set up the model formula
form <- obs ~ Intercept(1) + elev + elev2 + spatial(coordinates, model = spat)

## fit the model
m_inla <- bru(form, data = dat_spat_spdf, family = "binomial")

## check model summary
summary(m_inla)

## check spatial correlation
plot(spde.posterior(m_inla, "spatial", what = "matern.correlation"))

## turn the prediction df into a SPDF
dat_pred_spdf <- dat_pred
coordinates(dat_pred_spdf) <- c("x", "y")

## project into space
pred <- predict(
  m_inla, dat_pred_spdf,
  ~ arm::invlogit(Intercept + elev + elev2 + spatial))

dat_pred$p_pred <- pred$mean

## the color scale
colsc <- color_scale(dat_pred$p_pred, dat_pred$real_p_occ)

# plot the predicted spatial distribution
gg_p_inla <- ggplot(dat_pred, aes(x = x, y = y)) +
  geom_raster(interpolate = TRUE, aes(fill = p_pred)) +
  geom_point(data = dat_spat, aes(shape = factor(obs))) +
  colsc +
  labs(title = "Estimated probability\nof occurence (inla)")

# plot the real predicted spatial distribution
gg_real_glm <- ggplot(dat_pred, aes(x = x, y = y, fill = real_p_occ)) +
  geom_raster(interpolate = TRUE) +
  colsc +
  labs(title = "Real probability\nof occurence")



gg_spat_comp <- gridExtra::grid.arrange(gg_p_inla, gg_real_glm, ncol = 2)

## point pattern analysis

# load the data
data("gorillas", package = "inlabru")

# look at the data
ggplot() +
  gg(gorillas$gcov$waterdist, alpha = 0.5) +
  gg(gorillas$boundary, size = 2) +
  gg(gorillas$nests) +
  scale_fill_continuous(type = "viridis")
  
# set up the SPDE
spde <- inla.spde2.pcmatern(gorillas$mesh,
                            prior.range = c(1, 0.05),
                            prior.sigma = c(5, 0.05))
# get te covariate
water <- gorillas$gcov$waterdist
elev <- gorillas$gcov$elevation

# set up the model formula
form <- coordinates ~ Intercept(1) + water(waterdist, model = "linear") +
  spatial(coordinates, model = spde)

# fit the model
m1 <- lgcp(form, data = gorillas$nests,
           samplers = gorillas$boundary,
           domain = list(coordinates = gorillas$mesh))

# explore model summary
summary(m1)

# derive prediction over space
pred <- predict(m1, pixels(gorillas$mesh, nx = 50, ny = 50),
                ~ list(spatial = spatial,
                       water = water,
                       all = spatial + water))

# put the predictions in one df
pred_df <- data.frame(x = rep(coordinates(pred[[1]])[,1], 3),
                      y = rep(coordinates(pred[[1]])[,2], 3),
                      pred = as.numeric(sapply(pred, function(x) x$mean)),
                      type = rep(c("spatial", "water", "all"), each = 2148))

ggplot(pred_df, aes(x = x, y = y)) +
  geom_raster(interpolate = TRUE, aes(fill = pred)) +
  gg(gorillas$nests, alpha = 0.3) +
  facet_wrap(~type) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd")) +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

# derive prediction across water gradient
water_pred <- predict(
  m1,
  data = data.frame(waterdist = seq(min(water$waterdist), 
                                    max(water$waterdist),
                                    length.out = 10)),
  formula = ~ water_eval(waterdist)
)

ggplot(water_pred, aes(x = waterdist, y = mean, ymin = q0.025, ymax = q0.975)) +
  geom_line() +
  geom_ribbon(alpha = 0.2)

# let's do another model with elevation as a quadratic term
elev_sq <- elev
elev_sq$elevation_sq <- elev_sq$elevation ** 2
# set up the model formula
form <- coordinates ~ Intercept(1) + elev(elevation, model = "linear") +
  elev_sq(elevation_sq, model = "linear") +  spatial(coordinates, model = spde)

# fit the model
m2 <- lgcp(form, data = gorillas$nests,
           samplers = gorillas$boundary,
           domain = list(coordinates = gorillas$mesh))

# compare the IC
deltaIC(m1, m2)
