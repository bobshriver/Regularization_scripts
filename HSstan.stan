
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
real m0 = 6;           // Expected number of large slopes 
real slab_scale = 3;    // Scale for large slopes //Threshold needed to pass to be considered important//
real slab_scale2 = square(slab_scale);
real slab_df = 1;      // Effective degrees of freedom for large slopes
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
sigmaYr ~ cauchy(0, 5);
sigma ~ cauchy(0, 5);
eps ~ cauchy(0, 5);

}

