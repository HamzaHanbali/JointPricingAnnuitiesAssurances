
## R code for replicating Section 5 in "Pricing insurance contracts with offsetting liabilities"

## The user should have downloaded the number of deaths and exposures for England and Whales (EW) and the United States (US).
## EWD and USD : matrices containing the number of deaths in EW and US, respectively, with 59 rows corresponding to years 1950:2008 and 61 columns correponding to ages 30:90.

##########################
### Fitting the Li-Lee ###
##########################

## Log crude rates
DR_E = log(EWD/EWE) ; DR_U = log(USD/USE) ; DR_T = log((EWD+USD)/(EWE+USE))
DRn_E = log(EWD/EWE) ; DRn_U = log(USD/USE) ; DRn_T = log((EWD+USD)/(EWE+USE))

## Step 1 - Estimating the common components AlphaT, BetaT and KappaT using singular value decomposition
AlphaT = colMeans(DRn_T) ; for(j in 1:ncol(DRn_T)){DRn_T[,j] = DRn_T[,j] - AlphaT[j]}
SVD_T = svd(DRn_T, 1, 1)
BetaT = SVD_T$v/sum(SVD_T$v) ; KappaT = SVD_T$u * sum(SVD_T$v) * SVD_T$d[1]

## Step 2 - Estimating the individual components BetaE and KappaE for EW, and BetaU and KappaU for US
for(j in 1:ncol(DRn_E)){DRn_E[,j] = DRn_E[,j] - BetaT[j]*KappaT - AlphaT[j]}
for(j in 1:ncol(DRn_U)){DRn_U[,j] = DRn_U[,j] - BetaT[j]*KappaT - AlphaT[j]}
SVD_E = svd(DRn_E, 1, 1) ; SVD_U = svd(DRn_U, 1, 1)

BetaE = SVD_E$v/sum(SVD_E$v) ; KappaE = SVD_E$u * sum(SVD_E$v) * SVD_E$d[1]
BetaU = SVD_U$v/sum(SVD_U$v) ; KappaU = SVD_U$u * sum(SVD_U$v) * SVD_U$d[1]


####################################
###   Simulation of the Kappas   ###
####################################
library(mvtnorm)

## Number of simulations and maturities of the contracts
NSim = 10000
Horizon = 30

## Means and standard deviations of the kappa's
M1 = mean(diff(KappaE))  ;SD1 = sd(diff(KappaE))
M2 = mean(diff(KappaU))  ;SD2 = sd(diff(KappaU))
Mt = mean(diff(KappaT));SDt = sd(diff(KappaT))

## Initializing the matrices of simulations of the kappa's
SimKappaE = matrix(0,ncol=Horizon,nrow=NSim)
SimKappaU = matrix(0,ncol=Horizon,nrow=NSim)
SimKappaT = matrix(0,ncol=Horizon,nrow=NSim)

## Multivariate normal random vector
ZZ = rmvnorm(NSim,sigma=cbind(c(1,cor(KappaE,KappaU),cor(KappaE,KappaT)),c(cor(KappaE,KappaU),1,cor(KappaU,KappaT)),c(cor(KappaE,KappaT),cor(KappaU,KappaT),1)))
Z1 = ZZ[,1] ; Z2 = ZZ[,2] ; Zt = ZZ[,3]

## Forecast
SimKappaE[,1] = KappaE[length(KappaE)] + M1 + SD1*Z1
SimKappaU[,1] = KappaU[length(KappaU)] + M2 + SD2*Z2
SimKappaT[,1] = KappaT[length(KappaT)] + Mt + SDt*Zt

for (k in 2:Horizon){
ZZ = rmvnorm(NSim,sigma=cbind(c(1,cor(KappaE,KappaU),cor(KappaE,KappaT)),c(cor(KappaE,KappaU),1,cor(KappaU,KappaT)),c(cor(KappaE,KappaT),cor(KappaU,KappaT),1)))
Z1 = ZZ[,1] ; Z2 = ZZ[,2] ; Zt = ZZ[,3]
SimKappaE[,k] = SimKappaE[,k-1] + M1 + SD1*Z1
SimKappaU[,k] = SimKappaU[,k-1] + M2 + SD2*Z2
SimKappaT[,k] = SimKappaT[,k-1] + Mt + SDt*Zt
}


## Transforing the kappa's into forces of mortality 
N0 = length(KappaU)
MuUS = matrix(0,ncol=Horizon,nrow=NSim)
MuEW = matrix(0,ncol=Horizon,nrow=NSim)
for (j in 1:30){
MuUS[,j] = exp(DR_U[N0,j]) * exp( BetaU[j]* (SimKappaU[,j] - KappaU[N0]) + BetaT[j]*(SimKappaT[,j] - KappaT[N0]) )
MuEW[,j] = exp(DR_E[N0,30+j]) * exp( BetaE[30+j]* (SimKappaE[,j] - KappaE[N0]) + BetaT[30+j]*(SimKappaT[,j] - KappaT[N0]) )
}

## Transforming the forces of mortality into survival probabilities (for annuity using EW data) and into survival and death probabilities (for assurances using US data)
PUS = matrix(0,ncol=Horizon,nrow=NSim) ; PEW = matrix(0,ncol=Horizon,nrow=NSim)
for (i in 1:NSim){PUS[i,] = exp( - cumsum( MuUS[i,]) ) ; PEW[i,] = exp( - cumsum( MuEW[i,]) ) }
QUS = 1 - exp(-MuUS)


## Simulations of assurances and annuities using 1% interest rate
Discount = 1/1.01
Assurances = rowSums((Discount^t(matrix(rep((seq(1,Horizon)),NSim),ncol=NSim))) * QUS * (cbind(rep(1,NSim),PUS[,-Horizon])))
Annuities = rowSums((Discount^t(matrix(rep((seq(1,Horizon)),NSim),ncol=NSim))) * PEW)

#####################################################################
###### Premium Loading under <<standard deviation principle>> #######
#####################################################################

##
## Figure 1, Section 5
##


## Function for conditional value-at-risk (RM_TV), value-at-risk (RM_V) and mean-standard deviation (RM_sd)
RM_TV = function(X,level) {as.numeric(mean(X[which(X>=quantile(X,level))]))}
RM_V = function(X,level) {as.numeric(quantile(X,level))}
RM_sd = function(X,level) {mean(X) + level*sd(X)}

## Function to calculate the joint loading. The function returns the individual loadings Psi_60 (annuities) and Psi_30 (assurances)
##   as well as the joint loading Psi_Ptf. The function also returns the corresponding Premium under join loading (NH for natural hedging) and stand-alone (SA)
Psis = function(Level,Target,b,B,n,RM){
B_60 = B ; B_30 = b*B
CVaR_30 = RM(B_30*Assurances,Level)
Psi_30 = Target * (CVaR_30 - mean(B_30*Assurances)) / mean(B_30*Assurances)
Premium_30_SA = mean(B_30*Assurances)*(1 + Psi_30)
CVaR_60 = RM(B_60*Annuities,Level)
Psi_60 = Target * (CVaR_60 - mean(B_60*Annuities)) / mean(B_60*Annuities)
Premium_60_SA = mean(B_60*Annuities)*(1 + Psi_60)
## Premiums with Natural Hedging ##
Value_Ptf = n*B_30*Assurances + (1-n)*B_60*Annuities
CVaR_Ptf = RM(Value_Ptf,Level)
Psi_Ptf = Target*(CVaR_Ptf - mean(Value_Ptf)) / mean(Value_Ptf)
Premium_30_NH = mean(B_30*Assurances)*(1 + Psi_Ptf)
Premium_60_NH = mean(B_60*Annuities)*(1 + Psi_Ptf)
return( c(Psi_60,Psi_30,Psi_Ptf,Premium_30_NH,Premium_30_SA,Premium_60_NH,Premium_60_SA) )}

## Vectors of outpout, under different assumptions
pA_V = c() ; pB_V = c() ; pT_V = c()
pA_V2 = c() ; pB_V2 = c() ; pT_V2 = c()
pA_V3 = c() ; pB_V3 = c() ; pT_V3 = c()
pA_V4 = c() ; pB_V4 = c() ; pT_V4 = c()

pA_S = c() ; pB_S = c() ; pT_S = c()
pA_S2 = c() ; pB_S2 = c() ; pT_S2 = c()
pA_S3 = c() ; pB_S3 = c() ; pT_S3 = c()
pA_S4 = c() ; pB_S4 = c() ; pT_S4 = c()

Pols = seq(0,1,length.out=500)
for (n in 1:length(Pols)){
print(n)
PPT = Psis(0.95,0.5,1,1,Pols[n],RM_V) ; pA_V[n] = PPT[1] ; pB_V[n] = PPT[2] ; pT_V[n] = PPT[3]
PPT = Psis(0.95,0.5,10,1,Pols[n],RM_V) ; pA_V2[n] = PPT[1] ; pB_V2[n] = PPT[2] ; pT_V2[n] = PPT[3]
PPT = Psis(0.95,0.5,100,1,Pols[n],RM_V) ; pA_V3[n] = PPT[1] ; pB_V3[n] = PPT[2] ; pT_V3[n] = PPT[3]
PPT = Psis(0.95,0.75,10,1,Pols[n],RM_V) ; pA_V4[n] = PPT[1] ; pB_V4[n] = PPT[2] ; pT_V4[n] = PPT[3]

PPT = Psis(1.686,0.5,1,1,Pols[n],RM_sd) ; pA_S[n] = PPT[1] ; pB_S[n] = PPT[2] ; pT_S[n] = PPT[3]
PPT = Psis(1.686,0.5,10,1,Pols[n],RM_sd) ; pA_S2[n] = PPT[1] ; pB_S2[n] = PPT[2] ; pT_S2[n] = PPT[3]
PPT = Psis(1.686,0.5,100,1,Pols[n],RM_sd) ; pA_S3[n] = PPT[1] ; pB_S3[n] = PPT[2] ; pT_S3[n] = PPT[3]
PPT = Psis(1.686,0.75,10,1,Pols[n],RM_sd) ; pA_S4[n] = PPT[1] ; pB_S4[n] = PPT[2] ; pT_S4[n] = PPT[3]
}


## Visualization
mm = min(pT_V,pT_V2,pT_V3,pT_V4)
MM = max(pT_V,pT_V2,pT_V3,pT_V4)


xx = c(1,seq(20,length(Pols),20))
par(mfrow=c(2,2))
plot(cbind(Pols,pT_S),ylim=c(mm,MM),type='l',xlab='Proportion of term assurances',ylab=' ',main=expression('$\\zeta=0.5$ and $C_B=C_A$'))
points(cbind(Pols[xx],pT_V[xx]),col='blue',cex=0.5) ; lines(cbind(Pols,pA_S),lty=2) ; lines(cbind(Pols,pB_S),lty=2)
plot(cbind(Pols,pT_S2),ylim=c(mm,MM),type='l',xlab='Proportion of term assurances',ylab=' ',main=expression('$\\zeta=0.5$ and $C_B=10C_A$'))
points(cbind(Pols[xx],pT_V2[xx]),col='blue',cex=0.5) ; lines(cbind(Pols,pA_S2),lty=2) ; lines(cbind(Pols,pB_S2),lty=2)
plot(cbind(Pols,pT_S3),ylim=c(mm,MM),type='l',xlab='Proportion of term assurances',ylab=' ',main=expression('$\\zeta=0.5$ and $C_B=100C_A$'))
points(cbind(Pols[xx],pT_V3[xx]),col='blue',cex=0.5) ; lines(cbind(Pols,pA_S3),lty=2) ; lines(cbind(Pols,pB_S3),lty=2)
plot(cbind(Pols,pT_S4),ylim=c(mm,MM),type='l',xlab='Proportion of term assurances',ylab=' ',main=expression('$\\zeta=0.75$ and $C_B=10C_A$'))
points(cbind(Pols[xx],pT_V4[xx]),col='blue',cex=0.5) ; lines(cbind(Pols,pA_S4),lty=2) ; lines(cbind(Pols,pB_S4),lty=2)


################################
### Total collected premiums ###
################################

##
## Figure 2, Section 5
##


## Critical thresholds at portfolio level Nct, and at market level Wct
Nct = (2*Lambda2*piA)/(2*Lambda2*piA + (Lambda1 - 2*Lambda2)*piB)
Wct1 = Nct/(Nct + (1-Nct)*(1 + ((kB-1)/kB)*((PsiB-PsiA)/(1+PsiB)) *QB1 )*(kA/kB))
Wct2 = Nct/(Nct + (1-Nct)*(1 + ((kB-1)/kB)*((PsiB-PsiA)/(1+PsiB)) *QB2 )*(kA/kB))

## Function for the number of policyholders
Np = function(w,qq,kk,Pstar,PP,pp){
cc = ((1+Pstar)*pp-PP)/PP
Output = ((w)/kk)*(1 - ((kk-1)/kk)*qq*cc )
return( Output )}

## Main function to obtain the output, where 'bb' is the ratio of benefits. All other parameters are consistent with the notations in the paper
## The function returns the difference in total collected premiums between the two strategies

TotCol = function(bb,zeta,Gamma,w,qA,qB,kA,kB){

## Stand-alone loadings
PsiA = Gamma*zeta*sd(Annuities)/mean(Annuities)
PsiB = Gamma*zeta*sd(Assurances)/mean(Assurances)

## Correlation
Rho = cor(Annuities,Assurances)

## Pure premiums
piA = mean(Annuities)
piB = bb*mean(Assurances)

## Ratio of risk (denoted 'b' in the paper)
B = PsiB/PsiA

## Lambda's
Lambda1 = 1 + B^2 - 2*B*Rho
Lambda2 = 1 - B*Rho

## Minimum joint loading
Psimin = PsiA*sqrt((Lambda1 - Lambda2^2)/Lambda1)

## Equation (4.1) to determine Psi^star
psiF = function(ps){
nstar = Np(w,qB,kB,ps,(1+PsiB)*piB,piB)/(Np(w,qB,kB,ps,(1+PsiB)*piB,piB)+Np(1-w,qA,kA,ps,(1+PsiA)*piA,piA))
nts = (nstar*piB)/((nstar*piB)+((1-nstar)*piA))
Ps = PsiA*sqrt( Lambda1*nts*nts - 2*Lambda2*nts +1 )
return(Ps-ps)}

## Solving the equation
psol = uniroot(psiF,interval=c(0,PsiB))$root

## Stand-alone (sa) and joint (st) premiums for business lines A (annuities) and B (assurances)
PA_sa = (1+PsiA)*piA
PB_sa = (1+PsiB)*piB
PA_st = (1+psol)*piA
PB_st = (1+psol)*piB

## Stand-alone (sa) and joint (st) number of policyholders for business lines A (annuities) and B (assurances)
NA_st = Np(1-w,qA,kA,psol,(1+PsiA)*piA,piA)
NB_st = Np(w,qB,kB,psol,(1+PsiB)*piB,piB)
NA_sa = (1-w)/kA
NB_sa = w/kB

## Difference in total collected premiums
Diff = (NA_st*PA_st + NB_st*PB_st)/(NA_sa*PA_sa + NB_sa*PB_sa) - 1
return(Diff)}

## Calculating the output for different input values. 
KA1 = 10
KA2 = 5
KA3 = 15
KB = 10
QB1 = 0.5
QB2 = 3
QA1 = 0.5
QA2 = 3

Psis1_1 = c() ; Psis2_1 = c() ; Psis3_1 = c() ; Psis4_1 = c()
Psis1_2 = c() ; Psis2_2 = c() ; Psis3_2 = c() ; Psis4_2 = c()
Psis1_3 = c() ; Psis2_3 = c() ; Psis3_3 = c() ; Psis4_3 = c()


for (n in 1:length(Pols)){
print(n)
Psis1_1[n] = TotCol(10,0.5,1.686,Pols[n],QA1,QB1,KA1,KB)
Psis2_1[n] = TotCol(10,0.5,1.686,Pols[n],QA1,QB2,KA1,KB)
Psis3_1[n] = TotCol(10,0.5,1.686,Pols[n],QA2,QB1,KA1,KB)
Psis4_1[n] = TotCol(10,0.5,1.686,Pols[n],QA2,QB2,KA1,KB)

Psis1_2[n] = TotCol(10,0.5,1.686,Pols[n],QA1,QB1,KA2,KB)
Psis2_2[n] = TotCol(10,0.5,1.686,Pols[n],QA1,QB2,KA2,KB)
Psis3_2[n] = TotCol(10,0.5,1.686,Pols[n],QA2,QB1,KA2,KB)
Psis4_2[n] = TotCol(10,0.5,1.686,Pols[n],QA2,QB2,KA2,KB)

Psis1_3[n] = TotCol(10,0.5,1.686,Pols[n],QA1,QB1,KA3,KB)
Psis2_3[n] = TotCol(10,0.5,1.686,Pols[n],QA1,QB2,KA3,KB)
Psis3_3[n] = TotCol(10,0.5,1.686,Pols[n],QA2,QB1,KA3,KB)
Psis4_3[n] = TotCol(10,0.5,1.686,Pols[n],QA2,QB2,KA3,KB)
}


## Visualization
YY = 100*c(min(Psis1_1,Psis2_1,Psis3_1,Psis4_1,Psis1_2,Psis2_2,Psis3_2,Psis4_2,Psis1_3,Psis2_3,Psis3_3,Psis4_3),max(Psis1_1,Psis2_1,Psis3_1,Psis4_1,Psis1_2,Psis2_2,Psis3_2,Psis4_2,Psis1_3,Psis2_3,Psis3_3,Psis4_3))

plot(cbind(Pols,100*Psis1_1),type='l',ylim=YY,xlab='Proportion of market demand for term assurances',ylab=' ',cex.lab=2)
lines(cbind(Pols,100*Psis2_1),col='red')
lines(cbind(Pols,100*Psis3_1),col='magenta')
lines(cbind(Pols,100*Psis4_1),col='blue')
abline(h=0,lty=2)
points(cbind(Pols[seq(1,250,length(Pols)/6)],100*Psis1_1[seq(1,250,length(Pols)/6)]))
points(cbind(Pols[seq(15,250,length(Pols)/6)],100*Psis2_1[seq(15,250,length(Pols)/6)]),col='red',pch=2)
points(cbind(Pols[seq(30,250,length(Pols)/6)],100*Psis3_1[seq(30,250,length(Pols)/6)]),col='magenta',pch=4)
points(cbind(Pols[seq(45,250,length(Pols)/6)],100*Psis4_1[seq(45,250,length(Pols)/6)]),col='blue',pch=7)

legend('topleft',legend=c(expression('$q_A=0.5$ and $q_B=0.5$'),
				expression('$q_A=0.5$ and $q_B=3$'),
				expression('$q_A=3$ and $q_B=0.5$'),
				expression('$q_A=3$ and $q_B=3$')),
		    ,col=c('black','red','magenta','blue'),lty=c(1,1,1,1),pch=c(1,2,4,7),cex=1.5)


