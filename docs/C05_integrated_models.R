# integrated model of the form
# str ~ beta0 + beta1 * env + spat
# unstr ~ lambda0 + beta1 * env + spat


str_dat <- structured_data$whole
unstr_dat <- unstructured_data

coordinates(str_dat) <- ~ x + y
coordinates(unstr_dat) <- ~ x + y

# set up the mesh
mesh <- inla.mesh.2d(loc.domain = biasfield[,c(1,2)], 
                     max.edge = mesh.edge,   
                     cutoff = 2, 
                     offset = mesh.offset) 

# set up the spde
spatial <- inla.spde2.pcmatern(mesh,
                               prior.range = c(0.1, 0.05),
                               prior.sigma = c(5, 0.05))

# set up the formulas
form_common <- ~ -1 + intercept_str + intercept_unstr + spat(coordinates, model = spatial)
form_str <- presence ~ intercept_str + spat
form_unstr <- coordinates ~ intercept_unstr + spat

# set up the integration points
bnd_coords <- matrix(c(1, 1, 1, 100, 100, 100, 100, 1, 1, 1), 
                     ncol = 2, byrow = TRUE)
bnd <- SpatialPolygons(list(Polygons(list(Polygon(bnd_coords)),
                                     ID = "bnd")))
ips <- ipoints(bnd, mesh)

# set up the likelihood
lik_str <- like(formula = form_str,
                family = "binomial",
                data = str_dat)

lik_unstr <- like(formula = form_unstr,
                  family = "cp",
                  data = unstr_dat,
                  ips =ips)

# fit the model
m <- bru(form_common, lik_str, lik_unstr)

# plot the resulting spatial field
pix <- SpatialPixels(points = SpatialPoints(expand.grid(x = 1:100,
                                               y = 1:100)))
pred <- predict(m, pix,
                ~ list(spatial = spat))

pred_df <- data.frame(x = rep(coordinates(pred[[1]])[,1], 3),
                      y = rep(coordinates(pred[[1]])[,2], 3),
                      pred = as.numeric(sapply(pred, function(x) x$mean)))
ggplot(pred_df, aes(x=x, y=y, fill=pred)) +
  geom_raster(interpolate = TRUE)

# real spatial field
field <- expand.grid(x = 1:100,
                     y = 1:100)
field$spat <- as.numeric(dat1$Lam)
ggplot(field, aes(x=x, y=y, fill=spat)) +
  geom_raster(interpolate = TRUE)

## new models with covariate
# turn the env covariate into a SPDF
envdata <- SpatialPointsDataFrame(coords = expand.grid(x = 1:100, y = 1:100),
                       data = data.frame(env = as.numeric(dat1$gridcov)))

# fit a model with a common environmental covariate
f.env <- function(x, y){
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y))
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, envdata)
  if (any(is.na(v$env))) {
    v$env <- inlabru:::bru_fill_missing(envdata, spp, v$env)
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
                data = str_dat)

lik_unstr <- like(formula = form_unstr,
                  family = "cp",
                  data = unstr_dat,
                  ips =ips)

# fit the model
m <- bru(form_common, lik_str, lik_unstr)

# separate env variables
# set up the formulas
form_common <- ~ -1 + intercept_str + intercept_unstr + spat(coordinates, model = spatial) +
  env1(f.env(x, y), model = "linear") + env2(f.env(x, y), model = "linear")
form_str <- presence ~ intercept_str + spat + env1
form_unstr <- coordinates ~ intercept_unstr + spat + env2

# set up the likelihood
lik_str <- like(formula = form_str,
                family = "binomial",
                data = str_dat)

lik_unstr <- like(formula = form_unstr,
                  family = "cp",
                  data = unstr_dat,
                  ips =ips)

# fit the model
m <- bru(form_common, lik_str, lik_unstr)
