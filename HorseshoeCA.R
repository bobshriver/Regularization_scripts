###Set Path### 

HistData<-readRDS(file.path("Google Drive/range-resilience/Sensitivity/projections/data/historical/cleaned/historical_npp_climate_for_projections.rds")) 

CA<-subset.data.frame(HistData,region=="california_annuals")
CA<-CA[,-6]
head(CA)


attach(CA)

Y<-npp
X<-cbind(mm.mean,mm.dev,WatYrTEMP.mean,WatYrTEMP.dev)
Xuns<-X
X<-t((t(Xuns)-apply(Xuns,2,mean))/apply(Xuns,2,sd))
X<-cbind(X,  X[,'mm.mean']*X[,'mm.dev']  ,   X[,'mm.dev']*X[,'mm.dev']   ,   X[,'mm.dev']*X[,'mm.dev']*X[,'mm.mean']  ,   X[,'WatYrTEMP.mean']*X[,'WatYrTEMP.dev']   ,  X[,'mm.dev']*X[,'WatYrTEMP.dev']  )
X<-cbind(1,X)

colnames(X)<-c('Intcpt','PMean','Pdev','Tmean','Tdev','Pmean*Pdev','Pdev2','Pdev2*Pmean','Tmean*Tdev','Pdev*Tdev')

b<-solve(t(X)%*%X)%*%t(X)%*%Y ###Solve for th ML estimate of the regression parameters
b

Ypred<-X%*%b ###Predict the data using the MLE
plot(Y,Ypred, cex=.01) ###Plot Predicted Vs. Observed
abline(0,1)

Yresid<-(Y-Ypred) ###Calculate the residuals 
VarioData<-as.data.frame(cbind(Yresid,x,y))
colnames(VarioData)<-c('Yresid','lat','long') 
library(gstat) ###Needed for the variogram function
Vout<-variogram(Yresid~1,loc=~lat+long, width=.1,data=VarioData) ###calculate a variogram using the data. this will take a few minutes if the full dataset
plot(Vout)

phi<-1/3 ###kernal bandwidth---distance at which variogram stabalizes



plot(x[which(as.character(year)=="2015")],y[which(as.character(year)=="2015")], ylab="Lat", xlab="Long")###plot locations of all datapoints
knotloc<-expand.grid(seq(-122,-118,1),seq(38,33,-1))
points(knotloc[,1],knotloc[,2], col='red',pch=8)

Nsite<-length(which(as.character(year)=="2015"))
Nknot<-length(knotloc[,1])
siteloc<-cbind(x[which(as.character(year)=="2015")],y[which(as.character(year)=="2015")])
colnames(siteloc)<-c("Long","Lat")
colnames(knotloc)<-c("Long","Lat")
d<-as.matrix(dist(rbind(siteloc,knotloc),upper=T,diag=T)) ### Find distance between all points and knots
d<-d[-((Nsite+1):(Nsite+Nknot)),-(1:Nsite)]

remKnot<-which((apply(d,2,min)>1))
knotloc<-knotloc[-remKnot,]
Nknot<-length(knotloc[,1])
d<-d[,-remKnot]

plot(x[which(as.character(year)=="2015")],y[which(as.character(year)=="2015")], ylab="Lat", xlab="Long")
points(knotloc[,1],knotloc[,2], col='red',pch=8)


W<-exp(-d/phi)

K<-t(t(W)/apply(W,2,sum)) ##normalize kernel ###"t" is the transpose which flips matrix to make sure division is applied across columns. R does row-wise by default
apply(K,2,sum)###check that all columns sum to 1



data=list("par"=ncol(X),"Nsy"=length(Y), 'yrvec'=as.numeric(as.character(year))-1985,"Ny"=length(unique(year)), "Ns"=Nsite, "Nk"=Nknot, "P"=Y, "X"=X, "K"=K, "sitevec"=rep(1:Nsite,each=30))
library(rstan)
fit = stan('Google Drive/research/NPP_Model/Regularization_scripts/HSstan.stan', iter=1000, data=data,chains=1,refresh = 1,control = list(max_treedepth = 10), sample_file = 'horseshoe.csv')
save.image("Testfit.Rdata")

