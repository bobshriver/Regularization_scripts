###Set Path### 
setwd('Google Drive/range-resilience/Sensitivity/projections/data/historical/cleaned')

HistData<-readRDS(file.path("historical_npp_climate_for_projections.rds")) 

SGS<-subset.data.frame(HistData,region=="shortgrass_steppe")
SGS<-SGS[,-6]
head(SGS)


attach(SGS)

#Y<-npp
X<-cbind(mm.mean,mm.dev,WatYrTEMP.mean,WatYrTEMP.dev)
Xuns<-X
X<-t((t(Xuns)-apply(Xuns,2,mean))/apply(Xuns,2,sd))
X<-cbind(X,  X[,'mm.mean']*X[,'mm.dev']  ,   X[,'mm.dev']*X[,'mm.dev']   ,   X[,'mm.dev']*X[,'mm.dev']*X[,'mm.mean']  ,   X[,'WatYrTEMP.mean']*X[,'WatYrTEMP.dev']   ,  X[,'mm.dev']*X[,'WatYrTEMP.dev']  )
X<-cbind(1,X)

colnames(X)<-c('Intcpt','PMean','Pdev','Tmean','Tdev','Pmean*Pdev','Pdev2','Pdev2*Pmean','Tmean*Tdev','Pdev*Tdev')

b<-solve(t(X)%*%X)%*%t(X)%*%Y 
b

Ypred<-X%*%b
plot(Y,Ypred, cex=.01)
abline(0,1)

Yresid<-(Y-Ypred)
VarioData<-as.data.frame(cbind(Yresid[which(as.character(year)=="2015")],x[which(as.character(year)=="2015")],y[which(as.character(year)=="2015")]))
colnames(VarioData)<-c('Yresid','lat','long') 
library(gstat) ###Needed for the variogram function
Vout<-variogram(Yresid~1,loc=~lat+long, width=.1,data=VarioData) ###this will take a few minutes
plot(Vout)


plot(x[which(as.character(year)=="2015")],y[which(as.character(year)=="2015")], ylab="Lat", xlab="Long")
knotloc<-expand.grid(seq(-105,-100,1.5),seq(31,43,1.5))
points(knotloc[,1],knotloc[,2], col='red',pch=8)

Nsite<-length(which(as.character(year)=="2015"))
Nknot<-length(knotloc[,1])
siteloc<-cbind(x[which(as.character(year)=="2015")],y[which(as.character(year)=="2015")])
colnames(siteloc)<-c("Long","Lat")
colnames(knotloc)<-c("Long","Lat")
d<-as.matrix(dist(rbind(siteloc,knotloc),upper=T,diag=T)) ### Find distance between all points and knots
d<-d[-((Nsite+1):(Nsite+Nknot)),-(1:Nsite)]

remKnot<-which((apply(d,2,min)>.5))
knotloc<-knotloc[-remKnot,]
Nknot<-length(knotloc[,1])
d<-d[,-remKnot]

plot(x[which(as.character(year)=="2015")],y[which(as.character(year)=="2015")], ylab="Lat", xlab="Long")
points(knotloc[,1],knotloc[,2], col='red',pch=8)



phi<-1.5/3 ###kernal bandwidth
W<-exp(-d/phi)

K<-t(t(W)/apply(W,2,sum)) ##normalize kernel ###"t" is the transpose which flips matrix to make sure division is applied across columns. R does row-wise by default
apply(K,2,sum)###check that all columns sum to 1


STmodel="

data{
int Nsy;//This is the number of site-years
int Nk; //This is the number of knots
int Ns; //This is the number of sites
int Ny; //This is the number of years
int par; //This is the number of covariates

vector[Nsy] P; //This is all of the NPP measurments as a vector, length=Nsy
matrix[Nsy,par] X; //This is a matrix of covariates Nsy by 2, 1 column is intercept and the second precip

matrix[Ns,Nk] K; //This is the K matrix that translates knots to sites. Remember this contains paramters, but we fixed thes based on our variogram analysis. 

int sitevec[Nsy] ; //This is a vector that will be helpful in assigning random effects to each site. Each element will tell us which element of P corresponds to which site. Since each site has multiple measurments. 
int yrvec[Nsy] ; //This is a vector that will be helpful in assigning random effects to each site. Each element will tell us which element of P corresponds to which site. Since each site has multiple measurments. 


}

transformed data {
  real m0 = 3;           // Expected number of large slopes 
real slab_scale = 3;    // Scale for large slopes //Threshold needed to pass to be considered important//
real slab_scale2 = square(slab_scale);
real slab_df = 25;      // Effective degrees of freedom for large slopes
real half_slab_df = 0.5 * slab_df;
}

parameters{
vector[par] b_tilde; //This is a 2 length vector of all 
real<lower=0> eps; // This is the random effects SD
real<lower=0> sigma; //This is the SD on the process model
vector[Nk] alpha; //This is a vector of random effects for each knots
vector[Ny] alphaYr; //This is a vector of random effects for each year
real<lower=0> sigmaYr; //This is the SD on the year random effects

vector<lower=0>[par] lambda;
  real<lower=0> c2_tilde;
real<lower=0> tau_tilde;
}

transformed parameters {
  vector[par] b;

    real tau0 = (m0 / (par - m0)) * (sigma / sqrt(1.0 * Nsy));
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0);

    real c2 = slab_scale2 * c2_tilde;

    vector[par] lambda_tilde =
    sqrt( c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)) );


    b = tau * lambda_tilde .* b_tilde;
}

model{
vector[Nsy] mu;
vector[Ns] Eta; 

mu=X*b ;
Eta=K*alpha;


alpha~normal(0,eps);

P~normal(mu+Eta[sitevec]+(alphaYr[yrvec]*sigmaYr),sigma);

b_tilde ~ normal(0, 1);
lambda ~ cauchy(0, 1);
tau_tilde ~ cauchy(0, 1);
c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

alphaYr ~ normal(0, 1);
sigmaYr ~ cauchy(0, 1);


}



"


data=list("par"=4,"Nsy"=length(Y), 'yrvec'=as.numeric(as.character(year))-1985,"Ny"=length(unique(year)), "Ns"=Nsite, "Nk"=Nknot, "P"=Y, "X"=X, "K"=K, "sitevec"=rep(1:Nsite,each=30))
library(rstan)
fit = stan(model_code=STmodel, iter=1000, data=data,chains=1,refresh = 1,control = list(max_treedepth = 10), sample_file = 'horseshoe.csv')
save.image("Testfit.Rdata")

