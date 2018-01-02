######################################################
#Example 1
library("Rcpp")
library("RcppArmadillo")
sourceCpp("rss_bvsr.cpp")


betahat=matrix(c(0.05,1.2,0.3),1,3)
se=matrix(c(0.01,0.02,0.03),1,3)
Nsnp = c(100,200,300)
Ndraw = 1000
Nburn = 500
Nthin = 50

#RSS-BVSR
#INPUT:
   # betahat: 1 by p, effect size estimates under single-SNP model
   # se: 1 by p, standard errors of betahat
   # Nsnp: 1 by p, the number of individuals for each SNP
   # Ndraw: scalar, the total number of MCMC draws
   # Nburn: scalar, the length of burn-in
   # Nthin: scalar, the length of thinning interval
#OUTPUT:
    #logpi: Nsam by 1, the output chain of logpi
    #h: Nsam by 1, the output chain of h
    #output: p by 1, the relative PIP estimation
rss_bvsr(betahat,se,Nsnp,Ndraw,Nburn,Nthin)