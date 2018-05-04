########################
## LATENT SPACE MODEL ##
########################

library(latentnet)
load("DTI-connectome.RData") # load the DTI dataset
covariates <-  read.table("SUBJ-covariates.txt", header = T) # subject specific covariates

P = dim(D)[3] # number of subjects. No data available for 4 subjects -> they won't be considered in the analysis 
n = dim(D)[1] # number of brain regions
DD = 2 # number of dimensions in the latent space
Z_hoff <- array(NA, dim = c(n, DD, P-4)) # array where to store the estimated latent coordinates for each one of the 20 available subjects
beta_hoff <- matrix( NA, nrow = 4000, ncol = P-4) # matrix to store the chains of the beta intercept (4000 mcmc iterations)

missingSubj <- c(6, 17, 20, 22) # id of the missing subjects
avaSubj <- seq( 1, P )[ -missingSubj ] # id of the available subjects

for( p in (1:(P-4)) ) {  # estimation of the latent coordinates
  dati <- D[,, avaSubj[p], 1] # one subject at the time, considering the first scan 
  diag(dati) <- 0 # no self loop allowed
  mod = ergmm( as.network(dati) ~ euclidean(d = DD), verbose = TRUE, tofit = "mkl")
  Z_hoff[, , p] <- mod$mkl$Z # storing the estimated latent coordinates
  beta_hoff[, p] <- mod$sample$beta # storing chain of beta 
}

rm(list = ls())

###########################
## DYNAMIC NETWORK MODEL ##
###########################

load("fMRI-ROI-time-series.RData") # load the fMRI dataset
n <- 70 # number of brain regions
TT <- 400 # number of time points
subj <- 9 # choose a subject                           
scan <- 1 # chose if scan (1) or rescan (2)
weight.fun <- function(t,t.star,h){ # function to compute the weights. Provide: t.star, scalar; t, vector; h, scalar.
  K=exp((-(t-t.star)^2)/h)
  return(K/sum(K))
}

#--------------------------------------#

y <- Y[,1:TT, subj, scan] # data for subject subj and scan 1
A.rescan2 <- array(0, c(n, n, (TT)) )
conv.count <- matrix(0, nrow = n,ncol = TT)
h <- 2365
lambda <- 0.001

w.matrix=sapply(seq(2,TT),function(x) weight.fun(seq(2,TT),x,h)) # weight matrix (t.star: column, t: row)
S <- b <- rep(0, n)
y.tilde=rep(0,TT)

#--- BEGIN ---#

A.rescan2[,,1] <- diag(rep(1,n))
for( i in 1:n ){
  for( t.star in 2:TT){
    temp.A <- A.rescan2[i,, (t.star-1)]
    y.tilde <- sqrt(w.matrix[, t.star-1])*y[i, (2:TT)]
    prev.y.tilde <- sapply(1:(TT-1), function(k) y[, k]*sqrt(w.matrix[k, t.star-1]) )
    iter <- 0
    cond <- 1
    while(( cond > 0.001) & ( iter<=5 ) ){
      for( j in 1:n ){
        S[j] <- 2*mean((temp.A[-j]%*%prev.y.tilde[-j, ]-y.tilde)*prev.y.tilde[j, ] )
        b[j] <- 2*mean(prev.y.tilde[j, ]*prev.y.tilde[j, ] )
        temp.A[j] <- ifelse(abs(S[j] ) > lambda, (sign(S[j]-lambda )*lambda-S[j] )/b[j], 0 )
      }# end n (j)
      cond <- mean(abs(A.rescan2[i,, (t.star-1)] -temp.A) )
      iter <- iter+1
    }# end while
    A.rescan2[i,, t.star] <- temp.A
    conv.count[i, t.star] <- iter 
  }# end TT
}#end n (i)

