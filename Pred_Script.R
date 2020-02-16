###Step 1### 
library(rstan)

# Set working directory to source file location 
# (this method only works on Rstudio)
#library(rstudioapi)
#current_path <- getActiveDocumentContext()$path 
#setwd(dirname(current_path )) # set working directory to location of this file

##Load workspace with Stan model fit.##
load("./../parameters/........")

# MAKE FILENAMES
scenars <- c(1,seq(25,45,2)) # 1 is current, the others are all the RCP85, late century scenarios

ecoregion <- "CA" #name of ecoregion: CA, SGS, CD, NM, HD

###These lines make sure the model fit data match what the region specified is
if (ecoregion=="CA"){HistData=CA;fit=CAfit}
if (ecoregion=="CD"){HistData=CD;fit=CDfit}
if (ecoregion=="SGS"){HistData=SGS;fit=SGSfit}
if (ecoregion=="NM"){HistData=NM;fit=NMfit}
if (ecoregion=="HD"){HistData=HD;fit=HDfit}

###Read in sitematch table 
matchtable<-read.csv("./../data/site_coordinate_match.csv")
HistData<-merge(HistData,matchtable)

for(i in scenars){
  
  #import covariates
  infile <- paste0("./../data/future/by_scenario/covariates_sc",i,".csv")
  future_covars <- read.csv(infile)
  

chains<-3
Iter<-1500*chains #Number of MCMC interations* the number of chains
Ns<-Nsite #Number of site (aka points)
FYrs<-unique(future_covars$year)
NySim<-length(FYrs)#Number of simulation years in each scenario. 

###########################
###Add command to subset covariates by region and match sites####
###########################
###Remove sites not in ecoregion
####need vector to match sites in X with random effects
unqsite<-HistData$site[which(HistData$year==1986)] ###Take the name of each site by extracting 1 year of data. 
SiteMatch<-match(as.character(future_covars$site),unqsite) ###Match site names of sites with historic data locations in the historic data

future_covars<-future_covars[which(is.na(SiteMatch)==F),] ###remove locations that do not match with region sites
future_covars<-merge(future_covars,matchtable) ###as a backup merge x and y with future data. 

SiteMatch<-SiteMatch[which(is.na(SiteMatch)==F)]####This is the index needed for the site random effect below #remove NAs

###need to find mean for each site, and create annual deviates for each pixel

sitemeanPHist<-HistData$mm.mean[which(HistData$year==1986)] ###Get historic site mean precip and temp by just extracting one year
sitemeanTHist<-HistData$WatYrTEMP.mean[which(HistData$year==1986)]

sitemeanPFut<-aggregate((future_covars[,3]*10)~future_covars[,1],FUN=mean)[,2]
sitemeanTFut<-aggregate(future_covars[,4]~future_covars[,1],FUN=mean)[,2]

future_covars$mm.mean<-rep(sitemeanPFut,each=NySim)
future_covars$mm.dev<-(future_covars$WatYrPRECIP*10)/sitemeanP[SiteMatch] ###get future deviations from the historic mean ###convert future values from cm to mm
future_covars$WatYrTEMP.mean<-rep(sitemeanTFut,each=NySim)
future_covars$WatYrTEMP.dev<-(future_covars$WatYrTEMP)/sitemeanT[SiteMatch] ###get future deviations from the historic mean


Xfutuns<-future_covars[,c('mm.mean','mm.dev','WatYrTEMP.mean','WatYrTEMP.dev')]###Extract future covs in a way that matches historic data
  
  
Xfut<-t((t(Xfutuns)-apply(Xuns,2,mean))/apply(Xuns,2,sd))
Xfut<-cbind(Xfut,  Xfut[,'mm.mean']*Xfut[,'mm.dev']  ,   Xfut[,'mm.dev']*Xfut[,'mm.dev']   ,   Xfut[,'mm.dev']*Xfut[,'mm.dev']*Xfut[,'mm.mean']  ,   Xfut[,'WatYrTEMP.mean']*Xfut[,'WatYrTEMP.dev']   ,  Xfut[,'mm.dev']*Xfut[,'WatYrTEMP.dev']  )
Xfut<-cbind(1,Xfut)




  
  
###Step 2###
##Extract Parameters
alpha<-extract(fit,"alpha")$alpha ###Knot Random Effect
beta<-extract(fit,"b")$b
#sigma2<-extract(fit,"sigma")$sigma ###Gaussian Process Error SD
#sigma2Y<-extract(fit,"sigmaYr")$sigmaYr ###Year Random Effect Error SD


###Notes on error structures in simulation###
###Eta (spatial random effects) are fixed for each point across the simulation.
###Year random effects are randomly drawn for each simulation year, but are fixed across all points within a year. 
###Remaining Gaussian error (sigma) is completely IID for everysite/year within a given parameter set
### Parameter uncertainty we have thousands of different parameter values. Each climate change scenario is run with every set of parameters





###Step 4###
##We now want to simulate over parameter posteriors, years, and site to simualte the full future dataset for one climate scenario
##There are two ways to do this: 1) Vectorize years and sites and simulate the entire dataset in a single parameter loop, 2)Nest a time loop inside the parameter loop and vectorize only the sites
##The first approach will probably be faster, but let's start with #2 becuase it is easier to understand. 

Pred.out<-matrix(NA,NySim*Ns,Iter) ###This is an array to store the predictions from each 

for (p in 1:Iter) { ##This is a loop over all of the different parameter sets from the posterior simualtion. This captures parameter uncertainty. Iter is the number of iterations during the model fitting
  
  ### I've removed year random effects from simualtion since these will be averaged over###
  #YrRand<-rnorm(NySim,0,sigma2Y[p,]) ##Within each parameter set we want to draw a random effect of each year. NySim is the number of years the simualtion is for.  
  ####
  eta<-K%*%alpha[p,]
  mu<-Xfut%*%beta[p,]+eta[SiteMatch] ###This is the key step it take a Ns*Np Matrix with covariates (X) and multiplies it by a Np vector of parameters (beta[p,]). We then add on the Spatial random effect for that parameter set (eta[p,]). 
  ###mu is a vector of Ns length
  ### If we are only interested in the mean NPP then we can save mu, but this doesn't included the process error
  Pred.out[,p]<-mu #rnorm(Ns,mu+YrRand[t],sigma2[p,]) ###Predicts the future NPP including the proccess uncertaintiy (i.e. sigma2 and sigma2Y)


} ###End Parameter loop  

###Add colnames to file for site and year ####
Pred.out<-rbind(FYrs,Pred.out)###append the year to the first row 
rownames(Pred.out)<-future_covars$site ####rep site names for 


# write output to file...though maybe it will be .rds (arrays) instead of .csv?
outfile <- paste0("./../output/",ecoregion,"/projected_npp_sc",i,".csv")
write.csv(Pred.out, outfile, row.names = F)

} ###End Scenario loop

###If we run multiple scenarios we can add an extra loop outside of the parameter loop for scenarios. At the end of each scenario we can run summary stats (i.e. mean and 95% CI) for each scenario and only save that to reduce memory required.




