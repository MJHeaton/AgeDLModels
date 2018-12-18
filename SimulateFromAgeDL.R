##
## Fit the AgeDL model with NB liklihood
## using NIMBLE code
##

## Libraries
library(magrittr)
library(splines)
library(DiceDesign)
library(LatticeKrig)
library(dplyr)
library(lubridate)
library(ggplot2)

###################################
## Simulated Data Specifications ##
###################################
min.age <- 65
max.age <- 98
n.age <- length(min.age:max.age)
L <- 3 #Maximum non-zero lag
M <- 10 #Maximum lag
r <- 0.5 #Reduction factor for number of basis functions K = ceiling((1-r)*(L+1)*n.age)
nu <- 3.1 #Smoothness parameter of correlation function
alpha.age <- 0.5
alpha.lag <- 1.575
sigma2 <- 0.001 # Variance of DL coefficients
beta0 <- log(5) # Intercept for counts

##############################
## Load in Temperature Data ##
##############################
load("./AvgLagTemp.RData")
LagTemp <- as.matrix(LagTemp)[,1:(L+1)]

####################################
## Pred Proc basis knot Locations ##
####################################
n.knots <- ceiling((1-r)*(L+1)*n.age)
zero.locs <- expand.grid(L:M,min.age:max.age) %>% as.matrix()
nonzero.locs <- expand.grid(0:L,min.age:max.age) %>% as.matrix()
al.knots <- lhsDesign(n.knots,2)$design %>% maximinSA_LHS()
al.knots <- cbind(al.knots$design[,1]*L,min.age+al.knots$design[,2]*(max.age-min.age))
plot(zero.locs[,1], zero.locs[,2], pch=19, col="black", ylim=c(min.age,max.age), xlim=c(0,M))
points(nonzero.locs[,1], nonzero.locs[,2], pch=19, col="red")
points(al.knots[,1], al.knots[,2], pch=19, col="blue")
legend("topright", legend=c("Zero Effect", "Nonzero Effect", "Knots"), pch=19,
       col=c("black","red","blue"))

####################################
## Pred Proc Basis Functions from ##
## anisotropic corr function      ##
####################################
R.lag <- rdist(rbind(nonzero.locs,al.knots,zero.locs)[,1]) %>% Matern(alpha=alpha.lag,nu=nu)
R.age <- rdist(rbind(nonzero.locs,al.knots,zero.locs)[,2]) %>% Matern(alpha=alpha.age,nu=nu)
R <- R.lag*R.age
zero <- (nrow(nonzero.locs)+nrow(al.knots)+1):(nrow(R))
klocs <- nrow(nonzero.locs)+(1:nrow(al.knots))
R <- R[-zero,-zero]-R[-zero,zero]%*%solve(R[zero,zero],R[zero,-zero])
B <- R[-klocs,klocs]%*%solve(R[klocs,klocs])
pri.R <- sigma2*R[klocs,klocs]
pri.R.inv <- chol2inv(chol(R[klocs,klocs]))
TAind <- matrix(1:nrow(B),nrow=n.age,byrow=TRUE)

########################################
## Simulate and Plot a DLxAge Surface ##
########################################
DL <- matrix(B%*%t(chol(pri.R))%*%rnorm(nrow(pri.R)),nrow=L+1)
DL.df <- data.frame(Age=rep(min.age:max.age,each=L+1),Lag=rep(0:L,ncol(DL)),
                                        DL=c(DL))
DL.df %>% filter(Age%in%seq(65,max(Age),by=5)) %>% ggplot(aes(x=Lag,y=DL,color=Age))+geom_path()+theme_bw()+
  scale_color_distiller(palette="Spectral")+geom_hline(yintercept=0,lwd=1,lty=2)+
  ylab(expression(theta(a,l)))

##################################
## Simulate Counts from Poisson ##
##################################
mu <- exp(beta0+LagTemp[,1:(L+1)]%*%DL)
y <- matrix(rpois(length(mu), lambda=mu), nrow=nrow(mu))
save(file="./SimulatedAgeDLPois.RData",list=c("y","DL.df","beta0","B","n.knots", "LagTemp", "pri.R.inv"))
