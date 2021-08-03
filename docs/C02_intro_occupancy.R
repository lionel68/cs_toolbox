##########################################################################
#### Workshop CSD-Toolbox, 29.08.2021                                 ####
#### Swantje Loebel, Inst. für Geoökologie, TU Braunschweig           ####
#### Fitting Single-Season (static) Occupancy models                  ####
#### using unmarked and Jags                                          ####
##########################################################################

rm(list=ls())

# Read in libraries that we will need
library(unmarked)
library(boot)

###############################################################################
## Let's create some "typical" CSD data                                       #
###############################################################################

#We will (1) simulate a state process and then (2) the observation process.

# Choose sample sizes and prepare observed data array y
set.seed(100)                       # So we all get same data set
M <- 180                            # Number of sites
J <- 3                              # Number of presence/absence measurements (visits)
y <- matrix(NA, nrow = M, ncol = J) # to contain the observation data

#-------------------------------------------------------------------------------
#State process - occupancy is affected by vegetation height
#-------------------------------------------------------------------------------

# Create a covariate called vegHt (scaled)
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for occupancy model and compute occupancy
beta0 <- 0                    # Logit-scale intercept
beta1 <- 3                    # Logit-scale slope for vegHt
psi <- plogis(beta0 + beta1 * vegHt) # Occupancy probability (psi)
plot(vegHt, psi, ylim = c(0,1), type = "l", lwd = 3, col="chartreuse3") # Plot psi relationship

# And draw some random realization of this:
z <- rbinom(M, 1, psi)        # True presence/absence

# Look at data so far
table(z)

# Plot the true system state
plot(vegHt, z, xlab="Vegetation height", ylab="True presence/absence (z)", frame = F, cex = 1.5)
plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "chartreuse3")

#-------------------------------------------------------------------------------
# Detection process - detection assumed to be affected by wind and observer ID
#-------------------------------------------------------------------------------

# Create a numeric covariate called wind(speed), scaled
wind <- array(runif(M * J, -1, 1), dim = c(M, J))
wind

# Create a factor covariate called observed ID (3 different observers)
observer <- c(rep(1,3),rep(2,3),rep(3,3))

obsID <- matrix(observer, nrow=M, ncol=J,byrow = TRUE)
summary(as.factor(obsID))

# assume different detection rates for the three observers:
alpha0.1 <- -3                  # Logit-scale intercept, observer 1
alpha0.2 <- -1.5                  # Logit-scale intercept, observer 2
alpha0.3 <-  0               # Logit-scale intercept, observer 3

alpha0 <- c(ifelse(obsID==1,alpha0.1,ifelse(obsID==2,alpha0.2,alpha0.3)))
alpha1 <- -3

p <- plogis(alpha0 + alpha1 * wind) # Detection probability

# Plot relationships
plot(p ~ wind, ylim = c(0,1),col=c("green","red","blue")[as.factor(obsID)])    

# Take J = 3 presence/absence measurements at each site
for(j in 1:J) {
  y[,j] <- rbinom(M, z, p[,j])
}

sum(apply(y, 1, max))               # Number of sites with observed presences

#-------------------------------------------------------------------------------
# Plot observed data and true effects of wind and observer ID on detection probability
#-------------------------------------------------------------------------------
plot(wind, y, xlab="Wind", ylab="Observed det./nondetection data (y)", frame = F, cex = 1.5)
plot(function(x) plogis(alpha0.1 + alpha1*x), -1, 1, add=T, lwd=3, col = "green")
plot(function(x) plogis(alpha0.2 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")
plot(function(x) plogis(alpha0.3 + alpha1*x), -1, 1, add=T, lwd=3, col = "blue")
legend(x=0.3,y=0.9,c("Observer 1","Observer 2", "Observer 3"),
       col=c("green","red","blue"),lwd=3,bty="n")

# Look at the data: occupancy, true presence/absence (z), and measurements (y)
cbind(psi=round(psi,2), z=z, y1=y[,1], y2=y[,2], y3=y[,3])


################################################################################
#             Fitting occupancy models in unmarked                             #
################################################################################

# Load unmarked, format data and summarize
library(unmarked)

#create unmarked data frame
umf <- unmarkedFrameOccu(
  y = y,                                              # Pres/Abs measurements
  siteCovs = data.frame(vegHt = vegHt),               # site-specific covs.
  obsCovs = list(wind = wind, obsID = obsID))         # obs-specific covs.
summary(umf)

#---- Intercepts-Only model----------------------------------------------------- 
# Detection covariates follow first tilde, then occupancy covariates
summary(fm.occ0 <- occu(~1 ~ 1, data=umf))
summary(fm.occ0)

#use backTransform to go from logit scale to probability scale
backTransform(fm.occ0, "state")      # Get estimates on probability scale
backTransform(fm.occ0, "det")

#---- Add Occupancy covariate "Vegetation Height" (Detection constant) --------
summary(fm.occ1 <- occu(~ 1 ~ vegHt, data=umf))

#---- Occupancy covariate and numeric detection covariate "Wind" --------------
summary(fm.occ2 <- occu(~ wind ~ vegHt, data=umf))

#---- Occupancy covariate and detection covariates "Wind" & "Observer ID" (factor)
summary(fm.occ3 <- occu(~ as.factor(obsID) + wind ~ vegHt, data=umf))

# Compare AICs
data.frame(Model=c("fm.occ0","fm.occ1","fm.occ2","fm.occ2"),
           AIC=c(fm.occ0@AIC, fm.occ1@AIC, fm.occ2@AIC, fm.occ3@AIC),
           deltaAIC=c(fm.occ0@AIC, fm.occ1@AIC, fm.occ2@AIC, fm.occ3@AIC)-fm.occ0@AIC)

# Likelihood ratio test
LRT(fm.occ0, fm.occ1)
LRT(fm.occ1, fm.occ2)
LRT(fm.occ2, fm.occ3)  # fm.occ3 (full model) = best model

#-------------------------------------------------------------------------------
# Predict occupancy and detection as function of covs (with 95% CIs)
# and compare model predictions
#-------------------------------------------------------------------------------

# Fit detection-naive GLM to observed occurrence and plot comparison
summary(fm.glm <- glm(apply(y, 1, max) ~ vegHt, family=binomial))

# Plot occupancy probabilities
plot(vegHt, apply(y, 1, max), xlab="Vegetation height", ylab="Occurrence probability", frame = F, cex = 1.5)
plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red") #truth
lines(vegHt, predict(fm.glm,type="response"), type = "l", lwd = 3)
lines(vegHt, predict(fm.occ1, type="state")[,1], col = "seagreen", lwd = 3, lty=3)
lines(vegHt, predict(fm.occ2, type="state")[,1], col = "cornflowerblue", lwd = 3, lty=2)
lines(vegHt, predict(fm.occ3, type="state")[,1], col = "orange", lwd = 3, lty=4)

legend(-1, 0.95, c("Truth", "Occ Det constant","Occ Det Wind", "Occ Det Wind & Obs. ID", "naive GLM"), 
       col=c("red", "seagreen", "cornflowerblue","orange","black"), lty = c(1,3,2,4,1), lwd=3, bty="n",cex=1)

# conclusion: Model fm.occ3 is quite close to the "truth"

# Plot detection probabilities
newdat <- data.frame(wind=seq(-1, 1, 0.1))
pred.det <- predict(fm.occ1, type="det", newdata=newdat)
pred.det1 <- predict(fm.occ2, type="det", newdata=newdat)

layout(matrix(c(1,2),nrow=1)) 

plot(newdat$wind, pred.det$Predicted,type="l",lwd=3, col = "seagreen",ylim=c(0,1), 
     xlab="Wind", ylab="Detection probability")
lines(newdat$wind, pred.det$lower, lty=3, col="seagreen")
lines(newdat$wind, pred.det$upper, lty=3, col="seagreen")
lines(newdat$wind, pred.det1$Predicted, lwd=3, col="cornflowerblue")
lines(newdat$wind, pred.det1$upper, lty=3, col="cornflowerblue")
lines(newdat$wind, pred.det1$lower, lty=3, col="cornflowerblue")
legend("topright", c("Occ Det constant","Occ Det Wind"),col=c("seagreen","cornflowerblue"),lwd=3,bty="n")

# Wind speed & Observer ID
newdat2 <- data.frame(wind=rep(seq(-1, 1, 0.1),3),obsID=as.factor(rep(1:3,each=21)))
pred.det2 <- predict(fm.occ3, type="det", newdata=newdat2)

plot(newdat2$wind[newdat2$obsID==1], pred.det2$Predicted[newdat2$obsID==1], lwd=3, type="l",col="green",ylim=c(0,1), 
     xlab="wind", ylab="detection probability")
lines(newdat2$wind[newdat2$obsID==1], pred.det2$upper[newdat2$obsID==1], lty=3, col="green")
lines(newdat2$wind[newdat2$obsID==1], pred.det2$lower[newdat2$obsID==1], lty=3, col="green")
lines(newdat2$wind[newdat2$obsID==2], pred.det2$Predicted[newdat2$obsID==2], lwd=3, col="red")
lines(newdat2$wind[newdat2$obsID==2], pred.det2$upper[newdat2$obsID==2], lty=3, col="red")
lines(newdat2$wind[newdat2$obsID==2], pred.det2$lower[newdat2$obsID==2], lty=3, col="red")
lines(newdat2$wind[newdat2$obsID==3], pred.det2$Predicted[newdat2$obsID==3], lwd=3, col="blue")
lines(newdat2$wind[newdat2$obsID==3], pred.det2$upper[newdat2$obsID==3], lty=3, col="blue")
lines(newdat2$wind[newdat2$obsID==3], pred.det2$lower[newdat2$obsID==3], lty=3, col="blue")

legend("topright",title="Occ Det Wind & Observer ID", c("Observer 1","Observer 2", "Observer 3"),
       col=c("green","red","blue"),lwd=3,bty="n")

layout(1) 

# Predict detection probability for Observer factor levels at average 
# covariate values (Wind speed)
newdat3 <- data.frame(wind=0, obsID = c("1", "2", "3"))
predict(fm.occ3, type="det", newdata = newdat3, appendData = TRUE)

#Extract predicted z - true occurrence state
ranef(fm.occ3)

#compare with actually observations
obs <- apply(y,1,max)
cbind(obs, predicted_z = bup(ranef(fm.occ3), stat="mean"))[1:20,]

#When species was observed, then z is 1

#When species was never observed, then z is calculated by:
#example for site in row 1
(psi1 <- predict(fm.occ3, type="state")[1,1])
(p1 <- predict(fm.occ3, type="det")[c(1:3),1])
(z1 <- (psi1 * prod(1-p1)) / ((1 - psi1) + psi1 * prod(1-p1))) #this equation

# Define function for finite-sample number and proportion of occupied sites
fs.fn <- function(fm){
  Nocc <- sum(ranef(fm)@post[,2,])    #also could use sum(bup(ranef(fm.occ3), stat="mean"))
  psi.fs <- Nocc / nrow(fm@data@y)
  out <- c(Nocc = Nocc, psi.fs = psi.fs)
  return(out)
}

# Bootstrap the function
fs.hat <- fs.fn(fm.occ3)           # Point estimate
pb.fs <- parboot(fm.occ3, fs.fn, nsim=1000, report=10) # run for 1000 for final answer

# Summarize bootstrap distributions
summary(pb.fs@t.star)

# Get 95% bootstrapped confidence intervals
(tmp1 <- quantile(pb.fs@t.star[,1], prob = c(0.025, 0.975)))
(tmp2 <- quantile(pb.fs@t.star[,2], prob = c(0.025, 0.975)))

# Plot bootstrap distribution of number of occupied sites
hist(pb.fs@t.star[,1], col = "grey",border="darkgrey", breaks = 20, xlim = c(40, 140), 
     main = "Number of occupied sites",xlab = "Number of occupied sites", freq = F)
abline(v = fs.hat[1], col = "blue", lwd = 3)    # add point estimate
abline(v = tmp1, col = "blue", lty = 2)         # add 95% CI
abline(v = sum(apply(y, 1, max)), lwd = 3)      # observed #occ sites
abline(v = sum(z), col = "red", lwd = 3)        # true #occ sites

# Conclusion: our "full" occupancy model matches the "truth" quite well 


################################################################################
#             Fitting occupancy models in Jags                                 #
################################################################################

library(rjags)
library(jagsUI)
#library(coda)

setwd('...')  #set working directory (text files will be stored there)
#setwd("C:/Users/swantjeloebel/Desktop/Vorbereitung_CSD_Workshop/Occupancy_Models/Simulated_Data")

#-------------------------------------------------------------------------------
# Simple intercepts-only model
#-------------------------------------------------------------------------------

# Some remarks regarding the choice of priors for the intercept parameter in logistic models

# Problem: The choice of prior in occupancy (and other logistic) models is crucial, 
# especially when parameters are near the probability boundaries, 0 and 1 (Northrup & Gerber 2018).
# The common practice of using 'vague' normal priors (e.g. with a large standard deviation or 
# small precision (tau)) is problematic: the logit transformation is non-linear such that as values 
# become more negative or positive, the transformed probability values approach zero and one, respectively. 
# This non-linearity in the transformation leads to priors with large standard deviations becoming informative
# on the probability scale and strongly bimodal with large standard deviations. 

# Possible priors:
# (1) beta0 ~ dlogis(0,1)  --> recommended by Northrup & Gerber (2018) for intercept-parameters of logistic models,
#                              if models only include intercepts (and standardized covariates)
# (2) beta0 ~ dunif(0,1)
#     l.beta0 <-  logit(beta0) --> suggested by Kéry and Royle (2016), i.e. uniform prior that is than logit transformed 

# (3) beta0 ~ dnorm(0,1/2)     --> suggested Hobbs & Hooten (2015) and Northrup & Gerber (2018) if researchers want to use normal prior
#     beta0 ~ dnorm(0,1/2.71)  --> suggested by Lunn et al. (2013) 

# Literature
# Northrup & Gerber (2018). A comment on priors for Bayesian occupancy models. 
#                            PLoS One, 13, e0192819, doi: 10.1371/journal.pone.0192819
# Kéry & Royle (2016). Applied Hierarchical Modelling in Ecology. Volume 1.
# Hobbs & Hooten (2015).  Bayesian Models. A Statistical Primer for Ecologists.
# Lunn et al. (2013). The BUGS book.


# Specify model in BUGS language
cat(file = "StaticOcc_InterceptsOnly.txt","
model{ 
    
    ## Specify priors
    beta0.occ ~  dlogis(0,1)       # Intercept Occurrence probability 
    alpha0.det ~ dlogis(0,1)       # Intercept Detection probability    
                                                                              
    #----------------------------------------------------------------------
    # True state model for the partially observed true state   
    #----------------------------------------------------------------------
    
    for (i in 1:nSites) {           
    
      z[i] ~ dbern(psi[i])        ## True occupancy z at site i  
      logit(psi[i]) <- beta0.occ 
      
       #-------------------------------------------------------------------
       # Correction for imperfect detection
       #-------------------------------------------------------------------
    
      for (j in 1:nVisits) {                          # Loop over Visits
    
        # Observation model for the actual observations
         y[i,j] ~ dbern(muy[i,j])	# Detection-nondetection at i and j
         muy[i,j] <- z[i] * p[i,j]
         logit(p[i,j]) <- alpha0.det 

        # Computation of fit statistic (for Bayesian p-value)
         Presi[i,j] <- abs(y[i,j]-p[i,j])	 # Absolute residual
         y.new[i,j]~dbern(muy[i,j])
         Presi.new[i,j] <- abs(y.new[i,j]-p[i,j])
        
      } #j
    } #i
    
    
  fit <- sum(Presi[,])            # Discrepancy for actual data set
  fit.new <- sum(Presi.new[,]) 		# Discrepancy for replicate data set

  # Derived quantities
  occ.fs <- sum(z[])			    # Number of occupied sites
  mean.det.prob <- mean(p[,]) # Mean detectability  (alt: ilogit(alpha0.det))
}

")

# Inits function
zst <- apply(y, 1, max)			# Starting values for latent states
inits.fn <- function(){list(z = zst, beta0.occ=runif(1, -5, 5), alpha0.det = runif(1, -5, 5))}

# Parameters to estimate
para.names <- c("beta0.occ","alpha0.det", "mean.det.prob","mean.det.prob1","occ.fs", "fit", "fit.new")

#Data preparation
Data <- list(y=y, nSites=dim(y)[1], nVisits=dim(y)[2])

jagsModel0 <- jags(model.file="StaticOcc_InterceptsOnly.txt", data=Data, inits = inits.fn, n.chains = 3,
                      n.iter=5000, n.burnin=2500,n.thin=1,parameters.to.save=para.names, DIC=T)
print(jagsModel0)

# compare estimates with unmarked intercept-only model
data.frame(Mean=c(jagsModel0$mean$beta0.occ,jagsModel0$mean$alpha0.det),
           sd=c(jagsModel0$sd$beta0.occ,jagsModel0$sd$alpha0.det),
           row.names=c("beta0.occ","alpha0.det"))
fm.occ0@estimates


#-------------------------------------------------------------------------------
# Jags Model including vegetation height (occupancy),
#                             wind speed & Observer ID (detectability)
#-------------------------------------------------------------------------------

# Specify model in BUGS language
cat(file = "StaticOcc_Veg_WindObsID.txt","
model{ 
    
    ## Specify priors
    beta0.occ ~  dlogis(0,1)       # Intercept Occurrence probability 
    beta1.occ  ~  dnorm(0,0.01)    # Slope occurence covariate, note that in BUGS language
                                   # the normal distribution is defined by mean and tau = (1/sd^2)  
    for (a in 1:3) {
       alpha0.det[a] ~ dlogis(0,1)    }   # Intercepts Observer ID (categorial)
    
    alpha1.det ~ dnorm(0,0.01)       # Slope detectability covariate (Wind)
      
    #------------------------------------------------------------------
    # True state model for the partially observed true state   
    #------------------------------------------------------------------
    
    for (i in 1:nSites) {           
    
      z[i] ~ dbern(psi[i])        ## True occupancy z at site i  
      logit(psi[i]) <- beta0.occ + beta1.occ * vegHt[i]
      
       #-------------------------------------------------------------------
       # Correction for imperfect detection
       #-------------------------------------------------------------------
    
      for (j in 1:nVisits) {                          # Loop over Visits
    
        # Observation model for the actual observations
         y[i,j] ~ dbern(muy[i,j])	# Detection-nondetection at i and j
         muy[i,j] <- z[i] * p[i,j]
         logit(p[i,j]) <- alpha0.det[obsID[i,j]] + alpha1.det * wind[i,j]
    
        # Computation of fit statistic (for Bayesian p-value)
         Presi[i,j] <- abs(y[i,j]-p[i,j])	 # Absolute residual
         y.new[i,j]~dbern(muy[i,j])
         Presi.new[i,j] <- abs(y.new[i,j]-p[i,j])
        
      } #j
    } #i
    
    
  fit <- sum(Presi[,])            # Discrepancy for actual data set
  fit.new <- sum(Presi.new[,]) 		# Discrepancy for replicate data set
  
  # Derived quantities
  occ.fs <- sum(z[])			      # Number of occupied sites
  mean.det.prob <- mean(p[,])   # Mean detectability
  
}
")

# Inits function
zst <- apply(y, 1, max)			# Good starting values for latent states essential !
inits.fn <- function(){list(z = zst, beta0.occ=runif(1, -5, 5), 
                            beta1.occ=runif(1, -5, 5),  alpha0.det=runif(3, -5, 5), alpha1.det=runif(1, -5, 5))}

# Parameters to estimate
para.names <- c("beta0.occ","beta1.occ","alpha0.det","alpha1.det","mean.det.prob","occ.fs", "fit", "fit.new")
Data <- list(y=y,nSites=dim(y)[1],nVisits=dim(y)[2],vegHt=vegHt, wind=wind, obsID=obsID)

jagsModel3<- jags(model.file="StaticOcc_Veg_WindObsID.txt", data=Data, inits = inits.fn, n.chains = 3,
                   n.iter=5000, n.burnin=2500,n.thin=1,parameters.to.save=para.names, DIC=T)

print(jagsModel3)

Posterior0 <- as.data.frame(jagsModel0$sims.list)  #Intercepts-Only
Posterior1 <- as.data.frame(jagsModel3$sims.list)  #Vegetation height (occ) + Wind, ObserverID (obs)

# Plot posterior distribution of number of occupied sites
hist(Posterior1[,c("occ.fs")],freq=F,xlim=c(40,140),border="darkgrey")
polygon(density(Posterior1[,c("occ.fs")]),border="blue", lwd=2)
abline(v=sum(z[]),lwd=3,col="red")
abline(v=jagsModel3$mean$occ.fs,lwd=3,col="blue")
abline(v=jagsModel3$q2.5$occ.fs,lwd=2,lty=2,col="blue")
abline(v=jagsModel3$q97.5$occ.fs,lwd=2,lty=2,col="blue")
# add boots-trapped point estimates and 95% CIs for unmarked model
abline(v = fs.hat[1])             # add point estimate
abline(v = tmp1, lty = 2)         # add 95% CI


#-------------------------------------------------------------------------------
# Predict occupancy as function of covs (with 95% CIs)
# and compare model predictions
#-------------------------------------------------------------------------------

# Predictions of occupancy
vegHt <- matrix(seq(-1, 1, 0.1))
Posterior <- as.matrix(Posterior1[c("beta0.occ","beta1.occ")])

predOcc <- matrix(NA, nrow = dim(Posterior)[1], ncol = length(vegHt),dimnames=list(c(1:7500),c(seq(-1, 1, 0.1)))) # to contain the observation data
for (h in 1:length(vegHt)) {
  occupancy <- inv.logit(Posterior[,1] + Posterior[,2] * vegHt[h])
  predOcc[,h] <- occupancy  }

summary(predOcc)

# Plot occupancy probabilities
plot(vegHt, apply(predOcc,2,"mean"),type="l",lwd=3, col="blue",ylim=c(0,1),xlab="Vegetation height", ylab="Occurrence probability", frame = F, cex = 1.5)
lines(vegHt, apply(predOcc,2,quantile,probs =0.025),type="l",lty=2,col="blue", lwd=2)
lines(vegHt, apply(predOcc,2,quantile,probs =0.975),type="l",lty=2,col="blue", lwd=2)
#truth
plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red") #"truth"
#for comparison: unmarked model
lines(vegHt,predict(fm.occ3, type="state", newdata=as.data.frame(vegHt))[,"Predicted"], col="black")
lines(vegHt,predict(fm.occ3, type="state", newdata=as.data.frame(vegHt))[,"lower"],col="black",lty=3)
lines(vegHt,predict(fm.occ3, type="state", newdata=as.data.frame(vegHt))[,"upper"],col="black",lty=3)
legend(-1, 0.95, c("Truth", "Jags Model","Unmarked Model"),
       lwd=c(3,3,1),col=c("red","blue","black"),bty="n",cex=1)


## -----------------------------------------------------------------------------
## Some typical checks for JagsModels
##------------------------------------------------------------------------------

# Trace-plots (--> convergence (stationarity reached), chain mixing)
xyplot(jagsModel3$samples[,1:6], cex=2,cex.lab=2,cex.axis=2) 

traceplot(jagsModel3)

# Gelman-Diagnostics plots (--> convergence)
gelman.plot(jagsModel3$samples[,1:6],cex=2,cex.lab=1.3, cex.axis=1.3,cex.main=2)

# Posterior predictive check
par(mfrow=c(1,1))
plot(jagsModel3$sims.list$fit, jagsModel3$sims.list$fit.new, main = "jagsModel3",
     xlab = "Discrepancy for actual data set", 
     ylab = "Discrepancy for perfect data sets", las = 1)
abline(0,1, lwd = 2)

# Bayes p-value  # should be around 0.5, proportion of points above line
mean(jagsModel3$sims.list$fit.new > jagsModel3$sims.list$fit) 

# Correlation among parameters
Samples<-as.matrix(jagsModel3$samples)
Samples<-as.data.frame(Samples[seq(1,2500,25),])  #thin 

plot(Samples[,c("beta0.occ","beta1.occ","alpha0.det[1]","alpha0.det[2]",
        "alpha0.det[3]","alpha1.det")],col="blue")


#save.image("Session_OccupancyModels.RData")

