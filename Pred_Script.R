###Step 1### 
##Load workspace with Stan model fit.##
###

###Step 2###
##Extract Parameters
alpha<-extract(fit,"alpha")$alpha ###Knot Random Effect
beta<-extract(fit,"b")$b
sigma2<-extract(fit,"sigma")$sigma ###Gaussian Process Error SD
sigma2Y<-extract(fit,"sigmaYr")$sigmaYr ###Year Random Effect Error SD


###Notes on error structures in simulation###
###Eta (spatial random effects) are fixed for each point across the simulation.
###Year random effects are randomly drawn for each simulation year, but are fixed across all points within a year. 
###Remaining Gaussian error (sigma) is completely IID for everysite/year. 


###Step 3###
###Because Eta (spatial random effects) are fixed across all years, it is fastest to simulate these once up front 
eta<-matrix(NA,It,Ns) ###Create a matrix to store site random effects  

for (i in 1:It){eta[i,]<-K%*%alpha[i,]} 


