###Step 1### 
##Load workspace with Stan model fit.##
###

Iter<-#Number of MCMC interations
Ns<-#Number of site (aka points)
NySim<-#Number of simulation years in each scenario. 
X<-#Array of covariates structured as years by site by covariate  
  
###Step 2###
##Extract Parameters
alpha<-extract(fit,"alpha")$alpha ###Knot Random Effect
beta<-extract(fit,"b")$b
sigma2<-extract(fit,"sigma")$sigma ###Gaussian Process Error SD
sigma2Y<-extract(fit,"sigmaYr")$sigmaYr ###Year Random Effect Error SD


###Notes on error structures in simulation###
###Eta (spatial random effects) are fixed for each point across the simulation.
###Year random effects are randomly drawn for each simulation year, but are fixed across all points within a year. 
###Remaining Gaussian error (sigma) is completely IID for everysite/year within a given parameter set
### Parameter uncertainty we have thousands of different parameter values. Each climate change scenario is run with every set of parameters


###Step 3###
###Because Eta (spatial random effects) are fixed across all years, it is fastest to simulate these once up front because matrix manipuations in R can be slow
eta<-matrix(NA,Iter,Ns) ###Create a matrix to store site random effects  

for (i in 1:Iter){eta[i,]<-K%*%alpha[i,]} ###This multiplies the knot effects from each itteration (alpha[i,]) by the distance decay matrix to calculate the random effect for each point.
###

###Step 4###
##We now want to simulate over parameter posteriors, years, and site to simualte the full future dataset for one climate scenario
##There are two ways to do this: 1) Vectorize years and sites and simulate the entire dataset in a single parameter loop, 2)Nest a time loop inside the parameter loop and vectorize only the sites
##The first approach will probably be faster, but let's start with #2 becuase it is easier to understand. 

Pred.out<-array(NA,c(Iter,NySim,Ns)) ###This is an array to store the predictions from each 

for (p in 1:Iter) { ##This is a loop over all of the different parameter sets from the posterior simualtion. This captures parameter uncertainty. Iter is the number of iterations during the model fitting
  
YrRand<-rnorm(NySim,0,sigma2Y[p,]) ##Within each parameter set we want to draw a random effect of each year. NySim is the number of years the simualtion is for.  

for (t in 1:NySim) { ###Within a given parameter set this simulates each year for each site
  mu<-X[t,,]%*%beta[p,]+eta[p,] ###This is the key step it take a Ns*Np Matrix with covariates (X) and multiplies it by a Np vector of parameters (beta[p,]). We then add on the Spatial random effect for that parameter set (eta[p,]). 
  ###mu is a vector of Ns length
  ### If we are only interested in the mean NPP then we can save mu, but this doesn't included the process error
  Pred.out[p,t,]<-rnorm(Ns,mu+YrRand[t],sigma2[p,]) ###Predicts the future NPP including the proccess uncertaintiy (i.e. sigma2 and sigma2Y)

} ###End Time loop  
} ###End Parameter loop  

###If we run multiple scenarios we can add an extra loop outside of the parameter loop for scenarios. At the end of each scenario we can run summary stats (i.e. mean and 95% CI) for each scenario and only save that to reduce memory required.




