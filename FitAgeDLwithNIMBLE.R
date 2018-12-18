##
## NIMBLE Code to fit Age DL effects
## You need to have run the initialization code first
##

## Libraries
library(nimble)

## Load the data
load("./SimulatedAgeDLPois.RData")
int <- beta0
DL <- DL.df$DL

## Write out model code
ageDL <- nimbleCode({
  for(d in 1:D){
    for(a in 1:A){
      dl.effect[d,a] <- inprod(Temps[d,1:(L+1)],(DLBasis[((a-1)*(L+1)+1):((a-1)*(L+1)+1+L),1:n.DLBasis]%*%beta.DL[1:n.DLBasis]))
      y[d,a] ~ dpois(lambda=mu[d,a])
      mu[d,a] <- exp(intcpt+dl.effect[d,a])
    }
  }
  
  ## Prior Distributions
  intcpt ~ dnorm(mean=0,sd=1)
  beta.DL[1:n.DLBasis] ~ dmnorm(mean=zero[1:n.DLBasis],
                                prec=DLprec[1:n.DLBasis,1:n.DLBasis])
})

## Setup the nimble model
modelConsts <- list(D=nrow(LagTemp),A=diff(range(DL.df$Age))+1,
                    Temps=LagTemp,
                    DLBasis=B,L=ncol(LagTemp)-1,
                    n.DLBasis=ncol(B),zero=rep(0,ncol(B)),
                    DLprec=pri.R.inv)
modelData <- list(y=y)
logmuInit <- array(log(mean(y)),dim=dim(y))
muInit <- exp(logmuInit)
modelInitVals <- list(intcpt=int,
                      beta.DL=DL,
                      mu=muInit)
ageDLmodel <- nimbleModel(ageDL,constants=modelConsts,data=modelData,
                          inits=modelInitVals)
c.ageDLmodel <- compileNimble(ageDLmodel)

## Configure the MCMC
ageDL.MCMC <- configureMCMC(ageDLmodel)
ageDL.MCMC$printSamplers()

## Build the MCMC algorithm and compile it
ageDL.MCMCbuild <- buildMCMC(ageDL.MCMC)
c.ageDL.MCMC <- compileNimble(ageDL.MCMCbuild,project=c.ageDLmodel)

## Run the MCMC
system.time(result <- runMCMC(c.ageDL.MCMC,samplesAsCodaMCMC=TRUE,
                  nburnin=1000000,niter=2000000,thin=1000,inits=modelInitVals))
save(file="./FitAgeDLResults.RData",list=c("result"))




