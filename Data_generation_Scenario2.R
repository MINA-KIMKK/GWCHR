############################################################
############################################################
### Scenario 2
### Spatially varying coefficients
### Data generation
### Mina Kim
############################################################
############################################################

library(MASS)
library(foreach)
library(doParallel)
library(cvTools)
library(dplyr)
library(survival)


data.generator <- function(N, lam1, lam2, alp, beta11, beta12, beta21, beta22, probC, tau){
  getdata.f <- function(id, tau, lam1, lam2, alp, beta11, beta12, beta21, beta22, X1, X2, coor) {
    
    tbeta11 <- beta11+0.05*(coor[1]+coor[2])
    tbeta12 <- beta12+0.05*(coor[1]+coor[2])
    tbeta21 <- beta21+0.05*sqrt(coor[1]^2+coor[2]^2)
    tbeta22 <- beta22+0.05*sqrt(coor[1]^2+coor[2]^2)
    
    lam1 <- lam1 * exp(tbeta11 * X1 + tbeta12 * X2)
    lam2 <- lam2 * exp(tbeta21 * X1 + tbeta22 * X2)
    
    cur.t1 <- rexp(1, rate=lam1)
    cur.t2 <- rexp(1, rate=lam2)
    
    if (min(tau, cur.t1, cur.t2)==tau) {
      estart <- 0
      estop <- tau
      estatus <- 0
    } else if(min(tau, cur.t1, cur.t2)==cur.t1) {
      estart <- 0
      estop <- cur.t1
      estatus <- 1
    }else{
      estart <- 0
      estop <- cur.t2
      estatus <- 2
    }
    
    tmp <- data.frame(id = id, loc1 = coor[1], loc2 = coor[2], loc = coor[3], 
                      estart = estart, estop = estop, estatus = estatus,
                      tau = tau, X1 = X1, X2 = X2,  
                      beta11 = beta11, beta12 = beta12,
                      beta21 = beta21, beta22 = beta22)
    return(tmp)
  }
  
  
  if (probC == 0) {
    CC <- rep(tau, N)
  } else {
    CC <- rexp(N, rate = ((-1) * log(1 - probC))) 
    CC <- ifelse(CC > tau, tau, CC)
  }
  
  # covariates
  X1 <- runif(N, 0, 1)
  X2 <- rbinom(N, 1, 0.5)
  
  
  # coordinates representing the locations of a dataset 
  griddf <- cbind(as.matrix(expand.grid(latcoords = seq(1,10,1),
                                        lngcoords = seq(1,10,1))),loc=1:100)
  idx <- sample(1:100,N,replace=TRUE)
  coor <- griddf[idx,]
  
  
  event <- lapply(1:N, function(i) getdata.f(id = i, coor = coor[i,], X1 = X1[i], X2 = X2[i], 
                                             tau = CC[i], lam1 = lam1, lam2 = lam2, alp = alp, 
                                             beta11 = beta11, beta12 = beta12,
                                             beta21 = beta21, beta22 = beta22))
  data <- do.call(rbind, event)
  
  return(data)
}