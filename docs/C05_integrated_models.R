##########################################################################
#### Workshop CSD-Toolbox, 29.08.2021                                 ####
#### Lionel Hertzog, Thünen Institute for Biodiversity                ####
#### Fitting spatial models (georeferenced data plus point patterns)  ####
#### using INLA via inlabru                                           ####
##########################################################################


## set wd where the script and the data are located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load the libraries
library(tidyverse)
library(INLA)
library(inlabru)

# load the data
intdat <- readRDS("~/PostDoc_Thunen/Stat-stuff/Gfoe_workshop/site2/minimal-master/docs/C05_integrated_data.rds")

# look at the data
env_df <- as.data.frame(cbind(coordinates(intdat$env), intdat$env$env))

ggplot() +
  geom_raster(data = env_df, aes(x = x, y = y, fill = V3), interpolate = TRUE) +
  gg(data = intdat$str, aes(shape = factor(presence)), color = "grey70", size = 1.5)

ggplot() +
  geom_raster(data = env_df, aes(x = x, y = y, fill = V3), interpolate = TRUE) +
  gg(data = intdat$unstr, color = "grey70")





# set up the mesh
mesh <- inla.mesh.2d(loc.domain = raster::geom(intdat$bnd)[,c("x", "y")], 
                     max.edge =c(7, 14),   
                     cutoff = 2, 
                     offset = c(2, 7)) 

# set up the spde
spatial <- inla.spde2.pcmatern(mesh,
                               prior.range = c(0.1, 0.05),
                               prior.sigma = c(5, 0.05))

# set up the formulas
form_common <- ~ -1 + intercept_str + intercept_unstr + spat(coordinates, model = spatial)
form_str <- presence ~ intercept_str + spat
form_unstr <- coordinates ~ intercept_unstr + spat

# set up the integration points
ips <- ipoints(intdat$bnd, mesh)

# set up the likelihood
lik_str <- like(formula = form_str,
                family = "binomial",
                data = intdat$str)

lik_unstr <- like(formula = form_unstr,
                  family = "cp",
                  data = intdat$unstr,
                  ips =ips)

# fit the model
m_spat <- bru(form_common, lik_str, lik_unstr)

# plot the resulting spatial field
pix <- SpatialPixels(points = SpatialPoints(expand.grid(x = 1:100,
                                               y = 1:100)))
pred <- predict(m_spat, pix,
                ~ list(spatial = spat))
# turning the prediction / list object into a data.frame
pred_df <- data.frame(x = rep(coordinates(pred[[1]])[,1], 3),
                      y = rep(coordinates(pred[[1]])[,2], 3),
                      pred = as.numeric(sapply(pred, function(x) x$mean)))

gg_1 <- ggplot(pred_df, aes(x=x, y=y, fill=pred)) +
  geom_raster(interpolate = TRUE)

## new models with covariate

# create a function to interpolate the covariate values throughout the field
f.env <- function(x, y){
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y))
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, intdat$env)
  if (any(is.na(v$env))) {
    v$env <- inlabru:::bru_fill_missing(intdat$env, spp, v$env)
  }
  return(v$env)
}

# set up the formulas
form_common <- ~ -1 + intercept_str + intercept_unstr + spat(coordinates, model = spatial) +
  env(f.env(x, y), model = "linear")
form_str <- presence ~ intercept_str + spat + env
form_unstr <- coordinates ~ intercept_unstr + spat + env

# set up the likelihood
lik_str <- like(formula = form_str,
                family = "binomial",
                data = intdat$str)

lik_unstr <- like(formula = form_unstr,
                  family = "cp",
                  data = intdat$unstr,
                  ips =ips)

# fit the model
m_spatenv <- bru(form_common, lik_str, lik_unstr)

#new spatial field
pred_2 <- predict(m_spatenv, pix,
                ~ list(all = spat + env))
# turning the prediction / list object into a data.frame
pred2_df <- data.frame(x = rep(coordinates(pred_2[[1]])[,1], 3),
                      y = rep(coordinates(pred_2[[1]])[,2], 3),
                      pred = as.numeric(sapply(pred_2, function(x) x$mean)))

gg_2 <- ggplot(pred2_df, aes(x=x, y=y, fill=pred)) +
  geom_raster(interpolate = TRUE)

# make the two plots together
ggpubr::ggarrange(gg_1, gg_2, common.legend = TRUE)

## plot the unique effect of the covariate
pred_cov <- predict(m_spatenv, 
                    data = data.frame(env = seq(0, 1, length.out = 25)),
                  ~ env_eval(env))

ggplot(pred_cov, aes(x=env, y = mean, ymin = q0.025, ymax = q0.975)) +
  geom_ribbon(alpha = 0.1) +
  geom_line()


# compare the ICs
deltaIC(m_spat, m_spatenv)

# compare the matern correlation
p1 <- rbind(spde.posterior(m_spat, "spat", what = "matern.correlation"),
            spde.posterior(m_spatenv, "spat", what = "matern.correlation"))
p1$type <- rep(c("without env", "with env"), each = 200)

ggplot(p1, aes(x=x, y=median, ymin=q0.025, ymax=q0.975, color = type, fill = type)) +
  geom_ribbon(alpha = 0.3) +
  geom_line()
