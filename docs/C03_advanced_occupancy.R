#### useful R packages ######################

### some useful R packages for citizen science data

#sampBias - coordhttps://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05102
#CoordinateCleaner - https://github.com/ropensci/CoordinateCleaner
#https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_GBIF_data_with_CoordinateCleaner.html
#bRactus - https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13629
#sparta - https://github.com/BiologicalRecordsCentre/sparta
#rinat - https://cran.r-project.org/web/packages/rinat/index.html
#occAssess - https://www.biorxiv.org/content/10.1101/2021.04.19.440441v3.full
#Bio-geo -  https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.02118
#see SI for https://www.sciencedirect.com/science/article/pii/S235198941930633X
#naturaList - https://www.biorxiv.org/content/10.1101/2020.05.26.115220v1.full
#scrubr - https://rdrr.io/github/ropensci/scrubr/

#lots here:
#https://ropensci.org/packages/

#### practical ##############################

#`  Code largely taken from:`
#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   by Marc Kery & J. Andy Royle

# Read in libraries that we will need

library(AHMbook)
library(unmarked)
library(AICcmodavg)
library(boot)
library(tidyverse)

#use install.packages() if you dont have any of these packages

### Dynamic models: inferring colonization and extinction ####

set.seed(13973)

# Lets start with a simulation
M <- 50                                 # Number of sites 
J <- 3                                  # num secondary sample periods (i.e. visits)
T <- 10                                 # num primary sample periods (i.e., years)

#parameters
psi <- rep(NA, T)                       # Occupancy probability
muZ <- z <- array(dim = c(M, T))        # Expected and realized occurrence
y <- array(NA, dim = c(M, J, T))        # Detection histories
psi[1] <- 0.4                           # Initial occupancy probability
p <- runif(n=T, min=0.2, max=0.3)     # Detection probability
phi <- runif(n=T-1, min=0.6, max=0.8)   # Persistence probability (1-epsilon, extinction)   
gamma <- runif(n=T-1, min=0.1, max=0.2) # Colonization probability

# Generate latent states of occurrence

# First year
z[,1] <- rbinom(M, 1, psi[1])           # Initial occupancy state

# Later years
for(i in 1:M){                          # Loop over sites
  for(k in 2:T){                        # Loop over years
    muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1]
    z[i,k] <- rbinom(1, 1, muZ[k])}}


#Observation process:

# Generate detection/non-detection data
for(i in 1:M){
  for(k in 1:T){
    prob <- z[i,k] * p[k]
    for(j in 1:J){
      y[i,j,k] <- rbinom(1, 1, prob)
    }
  }
}
dim(y)

#plot simulated data
propOcc <- colMeans(z)
obsOcc <- colMeans(apply(y, c(1,3), max))

ggplot(data.frame(Year=1:T, propOcc, obsOcc)) +
  geom_line(aes(x = Year, y = propOcc)) +
  geom_line(aes(x = Year, y = obsOcc), color = "red", linetype = "dashed")+
  theme_classic()

#Lets fit the model with unmarked
library(unmarked)


#first, we need to format the data for unmarked

#organise the reponse long format
yy <- matrix(y, M, J*T)

#compare
y[1,,]
yy[1,]

#make spatio-temporal time variable
year <- matrix(c('01','02','03','04','05','06','07','08','09','10'),
               nrow(yy), T, byrow=TRUE)

simUMF <- unmarkedMultFrame(y = yy,
                            yearlySiteCovs = list(year = year),
                            numPrimary=T)
summary(simUMF)

# Model with all constant parameters
m0 <- colext(psiformula= ~ 1,       # First-year occupancy
             gammaformula = ~ 1,    # Colonization
             epsilonformula = ~ 1,  # Extinction
             pformula = ~ 1,        # Detection
             data = simUMF)
summary(m0)

names(m0)#we can extract all these elements
backTransform(m0, type="psi")
confint(backTransform(m0, type="psi"))

#model with year effects
m1 <- colext(psiformula = ~1,   
             gammaformula = ~ year-1,    
             epsilonformula = ~ year-1,  
             pformula = ~ year-1,        
             data = simUMF)
summary(m1)

#get predictions for each year
predict(m1, type='ext')
predict(m1, type='ext')
predict(m1, type='col')
predict(m1, type='det')

#Lets look at predicted occupancy each year
predictOccupancy <- function(fm) {
  
  psi.hat <- plogis(coef(fm, type="psi"))
  T <- fm@data@numPrimary
  gamma.hat <- plogis(coef(fm, type="col"))
  phi.hat <- 1 - plogis(coef(fm, type="ext"))
  
  for(t in 2:T) {
    psi.hat[t] <- psi.hat[t-1]*phi.hat[t-1] +(1-psi.hat[t-1])*gamma.hat[t-1]
  }
  return(psi.hat)
}

#also see: projected(m1)

pb <- parboot(m1, statistic = predictOccupancy, nsim=50)#increase to 1000
predsCI <- data.frame(cbind(pb@t0,t(apply(pb@t.star, 2, quantile, probs=c(0.025, 0.975)))))
colnames(predsCI) <- c("estimate", "lower", "upper")

#plot data
ggplot(data.frame(Year=1:T, propOcc, obsOcc, 
                  preds = predsCI$estimate,
                  lower = predsCI$lower, 
                  upper = predsCI$upper)) +
  geom_line(aes(x = Year, y = propOcc)) +
  geom_line(aes(x = Year, y = obsOcc), color = "red", linetype = "dashed")+
  geom_ribbon(aes(x = Year, y = preds, ymin = lower, ymax = upper), fill = "blue", alpha=0.4)+
  theme_classic()


#goodness of fit - experimental approach!!! check out the latest research
#https://rdrr.io/cran/AICcmodavg/man/mb.gof.test.html
mb.gof.test(m1)

#comparing different models
fl <- fitList(Null=m0, Annual=m1)
modSel(fl, nullmod="Null")


### Dynamics model in JAGS ####

data(frogs)
frogData <- masspcru
head(frogData)
str(frogData)

sort(unique(frogData$SurveyYear))
sort(unique(frogData$RouteNumStopNum))

#make response binary
frogData$Pcru <- ifelse(frogData$Pcru>0, 1, 0)

#make short to keep the computation times quick
frogData <- subset(frogData, RouteNumStopNum %in% 
                     unique(frogData$RouteNumStopNum)[1:20])
frogData$RouteNumStopNum <- factor(frogData$RouteNumStopNum)
table(frogData$RouteNumStopNum,frogData$SurveyYear)
#between 2 and 4 visits per year

#fit using unmarked
umf <- formatMult(frogData)
obsCovs(umf) <- scale(obsCovs(umf))

## constant transition rates
(fm <- colext(psiformula = ~ 1,
              gammaformula = ~ 1,
              epsilonformula = ~ 1,
              pformula = ~ JulianDate + MinAfterSunset, umf, control = list(trace=1, maxit=1e4)))

#we can plot the effect of each covariate, like, usual, with a predict function
predict(fm, type="p", newdata, appendData=TRUE)

##estimates of population parameters
projected(fm)

# same example in JAGS

sink("dynOccu.txt")
cat("
model {
  
  #priors
  
  #year 1 prior on psi
  psi1 ~ dunif(0,1)
  
  #prior on mean phi
  muphi.prob ~ dunif(0,1)
  lphi <- logit(muphi.prob)
  
  #prior on mean gamma
  mugam.prob ~ dunif(0,1)
  lgam <- logit(mugam.prob)
  
  #priors on detection probability
  p.prob ~ dunif(0,1)
  int.p <- logit(p.prob)
  beta.date ~ dnorm(0,0.01)
  beta.time ~ dnorm(0,0.01)
    
  #state model
  for(i in 1:nsite){
  
    #for year one
    z[i,1 ] ~ dbern(psi1)
    
    for(t in 2:nyear){
    
      #for subsequent years
      z[i,t] ~ dbern(psi[i,t])
      
      psi[i,t] <- z[i,t-1]*phi[i,t] + (1-z[i,t-1])*gamma[i,t]
      
      #survival model
      logit(phi[i,t]) <- lphi
      
      #colonization model
      logit(gamma[i,t]) <- lgam
      
    }
  }
  
  #detection model - written in long form (best for CS data with variable number of visits per site)
    for(j in 1:nvisit) {
    
    y[j] ~ dbern(Py[j]) #data is Y 
    Py[j]<- z[siteIndex[j],yearIndex[j]]*p[j] 
    
    #detection model:
    logit(p[j]) <-int.p + beta.date * date[j] + beta.time * time[j]
    
  }
  
  #metapopulation summaries
  for(t in 1:nyear){
    psi.fs[t] <- sum(z[1:nsite,t])/nsite  
  }
  
  }
        ",fill = TRUE)
sink()

# Bundle and summarize data set

#JAGS like indicies
frogData$RouteNumStopNum <- as.numeric(as.factor(frogData$RouteNumStopNum))
frogData$SurveyYear <- as.numeric(as.factor(frogData$SurveyYear))

jags.data <- list(nsite = n_distinct(frogData$RouteNumStopNum), 
                  nyear = n_distinct(frogData$SurveyYear), 
                  nvisit = nrow(frogData),
                  y = frogData$Pcru,
                  siteIndex = as.numeric(as.factor(frogData$RouteNumStopNum)),
                  yearIndex = as.numeric(as.factor(frogData$SurveyYear)),
                  date = as.numeric(scale(frogData$JulianDate)),
                  time = as.numeric(scale(frogData$MinAfterSunset)))


# Initial values
zst <- reshape2::acast(frogData, RouteNumStopNum ~ SurveyYear, value.var = "Pcru", max) # Observed occurrence as starting values for z
zst
inits <- function() list(z = zst)

# Parameters monitored
params <- c("muphi.prob","mugam.prob","beta.time","beta.date","psi.fs")

# MCMC settings
ni <- 2000   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3

# Call JAGS 
library(rjags)
library(jagsUI)
fm_dynamic <- jags(jags.data, inits, params, "dynOccu.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

fm_dynamic

plot(fm_dynamic)

### False presence ##############################################

#model of:
#Miller, D.A., J.D. Nichols, B.T. McClintock, E.H.C. Grant, L.L. Bailey, and L.A. #Weir. 2011. Improving occupancy estimation when two types of observational error #occur: non-detection and species misidentification. Ecology 92:1422-1428.

#3 types of data
#type 1 - no false positives
#type 2 - false negative and false positive detection probabilities
#type 3 - certain detections and uncertain detections that may include false positive detections. Uncertain detections are given a value of 1 and certain detections a value of 2.

n = 100 # number of sites
o = 10 #total number of occasions 
o1 = 5 #number of occasions with data type 1

#simulate each data type
y = matrix(NA,n,o)

#data type 1
r = 0.5 
y[1:(n*0.5),1:(o-o1)] <- rbinom((o-o1)*n*0.5,1,r) #first 5 occassions - data type 1

#data type 2
p = 0.7
y[1:(n*0.5),(o-o1+1):o] <- rbinom((n*o1*0.5),1,p) # second 5 occasions - data type 2

#add false positives to data type 2
fp = 0.05
y[(n*0.5+1):n,(o-o1+1):o] <- rbinom((n*o1*0.5),1,fp)

#create a covariate telling R of the amount of each data type
type <- c((o-o1),o1,0)  

#create observation-level data frame
occ <- matrix(c(rep(0,n*(o-o1)),rep(1,n*o1)),n,o)
occocc <- list(METHOD = occ)

#create site-level data
site <- data.frame(habitat = c(rep(1,n*0.5*.8),rep(0,n*0.5*.2),rep(1,n*0.5*.2),rep(0,n*0.8*0.5)))

#package up for unmarked
umf1 <- unmarkedFrameOccuFP(y, 
                            siteCovs = site, 
                            obsCovs = occ, 
                            type = type)

m1 <- occuFP(detformula = ~ METHOD, 
             FPformula = ~1,#covariates of false positive detection probability
             Bformula = ~1, #covariates of probability detections are certain
             stateformula = ~ 1, 
             data = umf1)

#summary
m1

#we can also:
coef(m1)
confint(m1, type = 'det')

#components of model
names(m1)
predict(m1, type = 'fp') 
predict(m1, type = 'state') 
predict(m1, type = 'det')


### Metacommunity data from the Swiss breeding bird survey MHB ####

#taken from Kery and Royle AHM book

data(MHB2014)
?MHB2014
str(MHB2014)

# Check the detection data in 3D array MHB2014$count: site x rep x species
nsite <- nrow(MHB2014$sites)    # number of sites in Swiss MHB
nrep <- 3                           # maximum number of replicate surveys per season
nspec <- nrow(MHB2014$species)   # 158 species occur in the 2014 data
dim(MHB2014$count) == c(nsite, nrep, nspec) # check

# Create the detection/nondetection (1/0) array
y <- MHB2014$count ; y[y > 1] <- 1  ## 'Y' replaced with 'y'
str(y)

# Check data for one species, here chaffinch, and pull them out from 3D array
(tmp <- y[, , "Common Chaffinch"])

# Frequency distribution of number of surveys actually carried out per site in 2014
# NB MHB2014$sites$nsurvey gives the number of surveys *planned*.
table(nsurveys <- apply(!is.na(y[,,1]), 1, sum))

# Which site has all NA data in 2014 ?
(NAsites <- which(nsurveys == 0) )

# Observed number of occupied sites
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
# For the 'all NA' site, max returns -Inf with a warning
tmp[tmp == -Inf] <- NA         # Change -Inf to NA
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))

# Plot species 'occurrence frequency' distribution (not shown)
plot(sort(obs.occ), xlab = "Species number", ylab = "Number of quads with detections")

# Drop data from species that were not observed in 2014
toss.out <- which(obs.occ == 0)
y <- y[,,-toss.out]
obs.occ <- obs.occ[-toss.out]

# Redefine nspec as the number of species observed in 2014: 145
(nspec <- dim(y)[3])

str(y)

# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum))     # Compute and print sorted species counts

plot(table(C), xlim = c(0, 60), xlab = "Observed number of species", ylab = "Number of quadrats", frame = FALSE)
abline(v = mean(C, na.rm = TRUE), col = "blue", lwd = 3)


### Independent species-specific random effects ####

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2], R = matrix(c(5,0,0,1), ncol = 2), df = 3) )

# Specify model in BUGS language
sink("model6.txt")
cat("
         model {
        
        # Priors
        for(k in 1:nspec){          # Loop over species
          psi[k] ~ dnorm(mu.psi, tau.psi)
          p[k] ~ dnorm(mu.p, tau.p)
        }
        prob.psi ~ dunif(0,1)
        prob.p ~ dunif(0,1)
        mu.psi <- logit(prob.psi)
        mu.p <- logit(prob.p)
        tau.psi ~ dgamma(0.001,0.001)
        tau.p ~ dgamma(0.001,0.001)
        
        # Ecological model for latent occurrence z (process model)
        for(k in 1:nspec){          # Loop over species
            for (i in 1:M) {         # Loop over sites
              z[i,k] ~ dbern(psi[k])
          }
        }
        
        # Observation model for observed data y
        for(k in 1:nspec){          # Loop over species
          for (i in 1:M) {
            mup[i,k] <- z[i,k] * p[k]
            ysum[i,k] ~ dbin(mup[i,k], J[i])
          }
        }
        
        # Derived quantities
        for(k in 1:nspec){          # Loop over species
          Nocc.fs[k] <- sum(z[,k]) # Add up number of occupied sites among the 267
        }
        for (i in 1:M) {            # Loop over sites
          Nsite[i] <- sum(z[i,])   # Add up number of occurring species at each site
          }
        }
        ",fill = TRUE)
sink()


# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as starting values for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, Omega = matrix(c(1,0,0,1), ncol = 2), eta = matrix(0, nrow = nspec, ncol = 2))

# Parameters monitored
params <- c("mu.eta", "probs", "psi", "p", "Nsite", "Nocc.fs", "Sigma", "rho")

# MCMC settings
ni <- 20000   ;   nt <- 15   ;   nb <- 5000   ;   nc <- 3

# Call JAGS from R (ART 12 min), check traceplots and summarize posteriors
out6 <- jags(win.data, inits, params, "model6.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(3,3))   ;   traceplot(out6, c('mu.eta', 'probs', 'Sigma', 'rho'))
print(out6, 3)


### The Dorazio-Royle (DR) community occupancy model with data augmentation ####

# Augment data set (DA part)
nz <- 150                # Number of potential species in superpopulation
M <- nspec + nz          # Size of augmented data set ('superpopulation')
yaug <- cbind(ysum, array(0, dim=c(nsite, nz))) # Add all zero histories

# Bundle and summarize data set
str( win.data <- list(yaug = yaug, nsite = nrow(ysum), nrep = MHB2014$sites$nsurvey, M = M, nspec = nspec, nz = nz) )

# Specify model in BUGS language
sink("model9.txt")
cat("
        model {
        
        # Priors to describe heterogeneity among species in community
        for(k in 1:M){                  # Loop over all species in augmented list
        lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
        lp[k] ~ dnorm(mu.lp, tau.lp)
        }
        
        # Hyperpriors to describe full community
        omega ~ dunif(0,1)              # Data augmentation or 'occupancy' parameter
        mu.lpsi ~ dnorm(0,0.001)        # Community mean of occupancy (logit)
        mu.lp ~ dnorm(0,0.001)          # Community mean of detection (logit)
        tau.lpsi <- pow(sd.lpsi, -2)
        sd.lpsi ~ dunif(0,5)            # Species heterogeneity in logit(psi)
        tau.lp <- pow(sd.lp, -2)
        sd.lp ~ dunif(0,5)              # Species heterogeneity in logit(p)
        
        # Superpopulation process:this is the 'paramater expansion' part of PX-DA
        for(k in 1:M){
        w[k] ~ dbern(omega)           # Metacommunity membership indicator
        }                               # (or data augmentation variable)
        
        # Ecological model for latent occurrence z (process model)
        for(k in 1:M){
        mu.psi[k] <- w[k] * psi[k]    # species not part of community zeroed out for z
        logit(psi[k]) <- lpsi[k]
        for (i in 1:nsite) {
        z[i,k] ~ dbern(mu.psi[k])
        }
        }
        
        # Observation model for observed detection frequencies
        for(k in 1:M){
        logit(p[k]) <- lp[k]
        for (i in 1:nsite) {
        mu.p[i,k] <- z[i,k] * p[k]  # non-occurring species are zeroed out for p
        yaug[i,k] ~ dbin(mu.p[i,k], nrep[i])
        }
        }
        
        # Derived quantities
        for(k in 1:M){
        Nocc.fs[k] <- sum(z[,k])     # Number of occupied sites among the 267
        }
        for (i in 1:nsite) {
        Nsite[i] <- sum(z[i,])       # Number of occurring species at each site
        }
        n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species in metacommunity
        Ntotal <- sum(w[])              # Total metacommunity size (= nspec + n0)
        }
        ",fill = TRUE)
sink()

# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at 'occurring'
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto for z
inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz))

# Parameters monitored
params <- c("mu.lpsi", "sd.lpsi", "mu.lp", "sd.lp", "psi", "p", "Nsite", "Ntotal", "omega", "n0")

# MCMC settings
ni <- 22000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call JAGS from R (ART 62 min), check convergence and summarize posteriors
out9 <- jags(win.data, inits, params, "model9.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2,2)) ; traceplot(out9, c('mu.lpsi', 'sd.lpsi', 'mu.lp', 'sd.lp'))
print(out9, dig = 3)


# Plot posterior distribution of site-specific species richness (Nsite)
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in 1:267){
  plot(table(out9$sims.list$Nsite[,i]), main = paste("Quadrat", i),
       xlab = "Local species richness", ylab = "", frame = F,
       xlim = c((min(C[i], out9$sims.list$Nsite[,i], na.rm = T)-2),
                max(out9$sims.list$Nsite[,i]) ))
  abline(v = C[i], col = "grey", lwd = 4)
  browser()
}

# Plot it only for a selection of sites
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
  plot(table(out9$sims.list$Nsite[,i]), main = paste("Quadrat", i),
       xlab = "Local species richness", ylab = "", frame = F,
       xlim = c((min(C[i], out9$sims.list$Nsite[,i], na.rm = T)-2),
                max(out9$sims.list$Nsite[,i]) ))
  abline(v = C[i], col = "grey", lwd = 4)
}


# Plot posterior distribution of total species richness (Ntotal)
plot(table(out9$sims.list$Ntotal), main = "", ylab = "", xlab = "Avian metacommunity size in Swiss MHB survey (267 1km2 quadrats)", frame = F, xlim = c(144, 245))
abline(v = nspec, col = "grey", lwd = 4)





