
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


parameters{
vector[par] b;
real<lower=0> eps; // This is the random effects SD
real<lower=0> sigma; //This is the SD on the process model
vector[Nk] alpha; //This is a vector of random effects for each knots
vector[Ny] alphaYr; //This is a vector of random effects for each year
real<lower=0> sigmaYr; //This is the SD on the year random effects


}



model{
vector[Nsy] mu;
vector[Ns] Eta; 

mu=X*b ;
Eta=K*(alpha*eps);


alpha~normal(0, 1);

P~normal(mu+Eta[sitevec]+(alphaYr[yrvec]*sigmaYr),sigma);

b ~ normal(0, 100);

alphaYr ~ normal(0, 1);
sigmaYr ~ cauchy(0, 1);
sigma ~ cauchy(0, 1);
eps ~ cauchy(0, 1);

}

