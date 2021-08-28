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
#   https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/

# Read in libraries that we will need

library(AHMbook)
library(unmarked)
library(AICcmodavg)
library(boot)
library(tidyverse)

#use install.packages() if you dont have any of these packages

### Dynamic models: inferring colonization and extinction ####

set.seed(20)

# Lets start with a simulation
M <- 50                                 # Number of sites 
J <- 3                                  # num secondary sample periods (i.e. visits)
T <- 10                                 # num primary sample periods (i.e., years)

#empty frames for data
psi <- rep(NA, T)                       # Occupancy probability
muZ <- z <- array(dim = c(M, T))        # Expected and realized occurrence
y <- array(NA, dim = c(M, J, T))        # Detection histories

#parameters
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
    z[i,k] <- rbinom(1, 1, muZ[k])
  }
}


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
trueOcc <- colMeans(z)
obsOcc <- colMeans(apply(y, c(1,3), max))

ggplot(data.frame(Year=1:T, trueOcc, obsOcc)) +
  geom_line(aes(x = Year, y = trueOcc)) +
  geom_line(aes(x = Year, y = obsOcc), color = "red", linetype = "dashed")+
  theme_classic()

#Lets fit the model with unmarked

#first, we need to format the data for unmarked

#organise the reponse long format
y_long <- matrix(y, M, J*T)

#compare
y[1,,]
y_long[1,]

#make spatio-temporal time variable
year <- matrix(c('01','02','03','04','05','06','07','08','09','10'),
               nrow(y_long), T, byrow=TRUE)

simUMF <- unmarkedMultFrame(y = y_long,
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

#extract components separately
names(m0)#we can extract all these elements
backTransform(m0, type="col")
confint(backTransform(m0, type="col"))

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

#Lets look at predicted occupancy each year and get 95%CI
predictOccupancy <- function(fm) {
  
  psi.hat <- plogis(coef(fm, type="psi"))
  gamma.hat <- plogis(coef(fm, type="col"))
  phi.hat <- 1 - plogis(coef(fm, type="ext"))
  
  T <- fm@data@numPrimary
  for(t in 2:T) {
    psi.hat[t] <- psi.hat[t-1]*phi.hat[t-1] +(1-psi.hat[t-1])*gamma.hat[t-1]
  }
  return(psi.hat)
}

#also see: projected(m1)

pb <- parboot(m1, statistic = predictOccupancy, nsim=20)#increase to 100
predsCI <- data.frame(cbind(pb@t0,t(apply(pb@t.star, 2, quantile, probs=c(0.025, 0.975)))))
colnames(predsCI) <- c("estimate", "lower", "upper")

#plot data
ggplot(data.frame(Year=1:T, trueOcc, obsOcc, 
                  preds = predsCI$estimate,
                  lower = predsCI$lower, 
                  upper = predsCI$upper)) +
  geom_line(aes(x = Year, y = trueOcc)) +
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
length(unique(frogData$RouteNumStopNum))

#make data frame shorter to keep the computation times quick
frogData <- subset(frogData, RouteNumStopNum %in% 
                     unique(frogData$RouteNumStopNum)[1:20])
frogData$RouteNumStopNum <- factor(frogData$RouteNumStopNum)

#scale data
frogData$JulianDate <- as.numeric(scale(frogData$JulianDate))
frogData$MinAfterSunset <- as.numeric(scale(frogData$MinAfterSunset))  

#do we have repeat visits per year?
table(frogData$RouteNumStopNum,frogData$SurveyYear)
#between 2 and 4 visits per year

#make response binary (just for teaching purposes!)
frogData$Pcru <- ifelse(frogData$Pcru>0, 1, 0)


#fit using unmarked for comparison
umf <- formatMult(frogData)

## constant transition rates
(fm <- colext(psiformula = ~ 1,
              gammaformula = ~ 1,
              epsilonformula = ~ 1,
              pformula = ~ JulianDate + MinAfterSunset, umf, control = list(trace=1, maxit=1e4)))

#we can plot the effect of each covariate, like, usual, with a predict function
predDF <- predict(fm, type="det", newdata = data.frame(JulianDate = mean(frogData$JulianDate), 
                                           MinAfterSunset = seq(min(frogData$MinAfterSunset),
                                                                max(frogData$MinAfterSunset),
                                                                length.out=100)), appendData=TRUE)

ggplot(predDF)+
  geom_line(aes(x = MinAfterSunset, y = Predicted))+
  geom_ribbon(aes(x = MinAfterSunset, 
                  ymin = Predicted - 1.96 * SE,
                  ymax = Predicted + 1.96 * SE),alpha=0.2)+
  theme_classic() +
  ylab("Predicted Detection Probability")


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
                  date = frogData$JulianDate,
                  time = frogData$MinAfterSunset)


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
fm_dynamic <- jags(jags.data, inits, params, "dynOccu.txt", n.chains = nc, n.thin = nt, n.iter = ni,
                   n.burnin = nb, parallel = TRUE)

#look at model coefficients
fm_dynamic

#check model convergence
plot(fm_dynamic)

### False presence ##############################################

#model of:
#Miller, D.A., J.D. Nichols, B.T. McClintock, E.H.C. Grant, L.L. Bailey, and L.A. #Weir. 2011. Improving occupancy estimation when two types of observational error #occur: non-detection and species misidentification. Ecology 92:1422-1428.

#3 types of data
#type 1 - no false positives
#type 2 - false negative and false positive detection probabilities
#type 3 - certain detections and uncertain detections that may include false positive detections. Uncertain detections are given a value of 1 and certain detections a value of 2.

#sampling
nsites = 100 # number of sites
nsurveys_type1 = 5 #total number of occasions (i.e. repeat visits) with type 1 data
nsurveys_type2 = 5 #total number of occasions (i.e. repeat visits) with type 2 data

#parameters
psi <- 0.5 #occupancy probability
p <- c(0.5,0.3)   #detection probability of method 1 and method 2
fp <- 0.05 #false positive probability

#simulate data
z <- rbinom(nsites, 1, psi)
y <- matrix(NA, nrow = nsites, ncol = nsurveys_type1 + nsurveys_type2)
for(i in 1:nsites){
  p1 <- p[1] * z[i]
  p2 <- p[2] * z[i]
  y[i,1:5] <- rbinom(nsurveys_type1, 1, p1)
  y[i, 6:10] <- rbinom(nsurveys_type2, 1, p2)
}

#method covariate 
Method = matrix(c(rep("1",nsurveys_type1),rep("2",nsurveys_type2)), nrow = nsites,
                ncol = nsurveys_type1 + nsurveys_type2, byrow = TRUE)
type = c(nsurveys_type1, nsurveys_type2, 0) 

#package data for unmarked
umf1 <- unmarkedFrameOccuFP(y = y, 
                            obsCovs = list(Method = Method),
                            type = type)

#fit model
m1 <- occuFP(detformula = ~ Method, 
             FPformula = ~ 1,#covariates of false positive detection probability
             stateformula = ~ 1, 
             data = umf1)

#summary
m1

#predicted probability
names(m1)
predict(m1, type="det", newdata = data.frame(Method=c("1","2")))
predict(m1, type="state", newdata = data.frame(Method=c("1","2")))

### Example code for community models ##############

#see Kery and Royle for full example

#could be analysed with e.g.:
data(MHB2014)
?MHB2014
head(MHB2014$counts)

### Independent species-specific random effects ####

# Specify model in BUGS language
sink("model.txt")
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


### The Dorazio-Royle (DR) community occupancy model with data augmentation ####

# Augment data set (DA part)
nz <- 150                # Number of potential species in superpopulation
M <- nspec + nz          # Size of augmented data set ('superpopulation')
yaug <- cbind(ysum, array(0, dim=c(nsite, nz))) # Add all zero histories

# Bundle and summarize data set
str( win.data <- list(yaug = yaug, nsite = nrow(ysum), nrep = MHB2014$sites$nsurvey, M = M, nspec = nspec, nz = nz) )

# Specify model in BUGS language
sink("model.txt")
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







