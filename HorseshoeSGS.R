SGS<-readRDS(file.path("sgs_covariates_for_regularization.rds")) 

SGS<-as.data.frame(SGS)
head(SGS)


attach(SGS)

Y<-npp.x
X<-cbind(mm.y,mm.dev,day_of_50_total_transp.y,transp.dev)
Xuns<-X
X<-t((t(Xuns)-apply(Xuns,2,mean))/apply(Xuns,2,sd))
X<-cbind(X,X[,1]*X[,2],X[,3]*X[,4])
X<-cbind(1,X)

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
knotloc<-expand.grid(seq(-105,-100,1),seq(31,43,1))
points(knotloc[,1],knotloc[,2], col='red',pch=8)

Nsite<-length(which(as.character(year)=="2015"))
Nknot<-length(knotloc[,1])
siteloc<-cbind(x[which(as.character(year)=="2015")],y[which(as.character(year)=="2015")])
colnames(siteloc)<-c("Long","Lat")
colnames(knotloc)<-c("Long","Lat")
d<-as.matrix(dist(rbind(siteloc,knotloc),upper=T,diag=T)) ### Find distance between all points and knots
d<-d[-((Nsite+1):(Nsite+Nknot)),-(1:Nsite)]

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

vector[Nsy] P; //This is all of the NPP measurments as a vector, length=Nsy
matrix[Nsy,7] X; //This is a matrix of covariates Nsy by 2, 1 column is intercept and the second precip

matrix[Ns,Nk] K; //This is the K matrix that translates knots to sites. Remember this contains paramters, but we fixed thes based on our variogram analysis. 

int sitevec[Nsy] ; //This is a vector that will be helpful in assigning random effects to each site. Each element will tell us which element of P corresponds to which site. Since each site has multiple measurments. 
int yrvec[Nsy] ; //This is a vector that will be helpful in assigning random effects to each site. Each element will tell us which element of P corresponds to which site. Since each site has multiple measurments. 


}

transformed data {
  real m0 = 4;           // Expected number of large slopes
real slab_scale = 3;    // Scale for large slopes
real slab_scale2 = square(slab_scale);
real slab_df = 25;      // Effective degrees of freedom for large slopes
real half_slab_df = 0.5 * slab_df;
}

parameters{
vector[7] b_tilde; //This is a 2 length vector of all 
real<lower=0> eps; // This is the random effects SD
real<lower=0> sigma; //This is the SD on the process model
vector[Nk] alpha; //This is a vector of random effects for each knots
vector[Ny] alphaYr; //This is a vector of random effects for each year
real<lower=0> sigmaYr; //This is the SD on the year random effects

vector<lower=0>[7] lambda;
  real<lower=0> c2_tilde;
real<lower=0> tau_tilde;
}

transformed parameters {
  vector[7] b;

    real tau0 = (m0 / (7 - m0)) * (sigma / sqrt(1.0 * Nsy));
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0);

    real c2 = slab_scale2 * c2_tilde;

    vector[7] lambda_tilde =
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


data=list("Nsy"=length(Y), 'yrvec'=as.numeric(as.character(year))-1985,"Ny"=length(unique(year)), "Ns"=Nsite, "Nk"=Nknot, "P"=Y, "X"=X, "K"=K, "sitevec"=rep(1:Nsite,each=30))
library(rstan)
fit = stan(model_code=STmodel, data=data,chains=1,refresh = 1,control = list(max_treedepth = 12), sample_file = 'horseshoe.csv')
