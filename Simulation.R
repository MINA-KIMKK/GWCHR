############################################################
### Simulation
### Mina Kim
############################################################

library(MASS)
library(foreach)
library(doParallel)
library(cvTools)
library(dplyr)
library(survival)


registerDoParallel(cores = 8)
optimal.m <- data.frame()
rr <- 1000
griddf <- cbind(as.matrix(expand.grid(latcoords = seq(1,10,1),
                                      lngcoords = seq(1,10,1))),loc=1:100)


for (r in 1:rr) {
  set.seed(r)
  data2 <- data.generator(N = 5000,
                           lam1 = 1,
                           lam2 = 0.5,
                           beta11 = log(1.2),
                           beta12 = log(0.8),
                           beta21 = log(1.3),
                           beta22 = log(0.4),
                           probC = 0.3,
                           tau = 1)
  
  dist.mat2 <- as.matrix(dist(griddf[,c(1,2)], method = "euclidean", diag = TRUE, upper = TRUE))
  m <- c(seq(1.5,max(dist.mat2),0.5))
  wt.c1 <- wt.c2 <-c()
  k=1

####################################################################
### Optimal bandwidth
  
  for (i in m) {
    print(paste("r =",r,",","i =",i))
    m.cvpl <-  foreach(j = 1:10,.combine=rbind,.packages=c("doParallel","cvTools","dplyr","survival")) %dopar% {
      set.seed(1)
      cv10 <- cvTools::cvFolds(n=nrow(data2), K=10, R=1)
      idx <- cv10$subsets[cv10$which==j,]
      data2.new <- data2[-idx, ]
      
      coef <- foreach(x = 1:nrow(dist.mat2),.combine=rbind,.packages=c("doParallel","cvTools","dplyr","survival")) %dopar% {
        wt <- ifelse(dist.mat2[x, ]<i, (1-(dist.mat2[x, ]/i)^2)^2, 0)
        datwt <- data.frame(wts = wt, loc = griddf[,3])
        data2.new2 <- data2.new %>% left_join(datwt, by="loc") %>% filter(wts!=0)
        model1 <- coxph(Surv(estop, estatus == 1) ~ X1 + X2, data=data2.new2, method="breslow", weights=wts)
        model2 <- coxph(Surv(estop, estatus == 2)~X1 + X2, data=data2.new2, method="breslow", weights=wts)
        
        return(c(x,model1$coef,model2$coef))
        
      }
      coef <- as.data.frame(coef)
      colnames(coef) <- c("loc","coef11","coef12","coef21","coef22")
      tmp <- data2 %>% left_join(coef, by = 'loc') %>%
        mutate(X1beta11 = X1*coef11,
               X2beta12 = X2*coef12,
               X1beta21 = X1*coef21,
               X2beta22 = X2*coef22)
      
      m1 <- coxph(Surv(estop, estatus == 1)~X1beta11 + X2beta12,
                  data=tmp[idx,], init=c(1,1), control=coxph.control(iter.max=0))
      m2 <- coxph(Surv(estop, estatus == 2)~X1beta21 + X2beta22,
                  data=tmp[idx,], init=c(1,1), control=coxph.control(iter.max=0))
      return(data.frame(m1=m1$loglik[2], m2=m2$loglik[2]))
    }
    wt.c1[k] = sum(m.cvpl$m1)
    wt.c2[k] = sum(m.cvpl$m2)
    k=k+1
  }
  opt1=m[which(wt.c1==max(wt.c1))]
  opt2=m[which(wt.c2==max(wt.c2))]
  
  optimal.bw <- rbind(optimal.m,c(r,opt1,opt2))
}
names(optimal.bw) <- c("iter","bw1","bw2")


####################################################################
### Spatially varying coefficients & summary

data.all <- foreach(rr = 1:1000,.combine=rbind,.packages=c("doParallel","dplyr","survival")) %dopar% {
  
  set.seed(rr)
  data2 <- data.generator2(N = 5000,
                           lam1 = 1,
                           lam2 = 0.5,
                           beta11 = log(1.2),
                           beta12 = log(0.8),
                           beta21 = log(1.3),
                           beta22 = log(0.4),
                           probC = 0.3,
                           tau = 1)
  
gw.coef <- foreach(j = 1:nrow(dist.mat2),.combine=rbind,.packages=c("doParallel","dplyr","survival")) %dopar% {
  ## GWCRR
  # Cause 1
  wt <- ifelse(dist.mat2[j, ]<optimal.bw$bw1[rr], (1-(dist.mat2[j, ]/optimal.bw$bw1[rr])^2)^2, 0)
  datwt <- data.frame(wts = wt, loc = griddf[,3])
  data3 <- data2 %>% left_join(datwt, by="loc") %>% filter(wts!=0)
  gw.cox1 <- coxph(Surv(estop, estatus == 1) ~ X1 + X2, data=data3, method="breslow", weights=wts)
  gw.coef11 <- gw.cox1$coefficients[1]
  gw.coef12 <- gw.cox1$coefficients[2]
  
  # Cause 2
  wt <- ifelse(dist.mat2[j, ]<optimal.bw$bw2[rr], (1-(dist.mat2[j, ]/optimal.bw$bw2[rr])^2)^2, 0)
  datwt <- data.frame(wts = wt, loc = griddf[,3])
  data3 <- data2 %>% left_join(datwt, by="loc") %>% filter(wts!=0)
  gw.cox2 <- coxph(Surv(estop, estatus == 2) ~ X1 + X2, data=data3, method="breslow", weights=wts)
  gw.coef21 <- gw.cox2$coefficients[1]
  gw.coef22 <- gw.cox2$coefficients[2]
  
  return(data.frame(gw.coef11,gw.coef12,gw.coef21,gw.coef22,
                    t(sqrt(diag(vcov(gw.cox1)))), t(sqrt(diag(vcov(gw.cox2))))))
}

# Global
glb.cox1 <- coxph(Surv(estop, estatus == 1)~X1 + X2, data2)
glb.cox2 <- coxph(Surv(estop, estatus == 2)~X1 + X2, data2)
glb.coef11 <- glb.cox1$coefficients[1]
glb.coef12 <- glb.cox1$coefficients[2]
glb.coef21 <- glb.cox2$coefficients[1]
glb.coef22 <- glb.cox2$coefficients[2]


return(data.frame(rr, gw.coef,
                  glb.coef11,glb.coef12,glb.coef21,glb.coef22,
                  t(sqrt(diag(vcov(glb.cox1)))), t(sqrt(diag(vcov(glb.cox2))))))
}

names(data.all) <- c("rr", "gw.coef11", "gw.coef12", "gw.coef21", "gw.coef22", 
                     "se.gw.coef11","se.gw.coef12","se.gw.coef21","se.gw.coef22", 
                     "glb.coef11","glb.coef12","glb.coef21","glb.coef22",
                     "se.glb.coef11","se.glb.coef12","se.glb.coef21","se.glb.coef22")

data.cp <- cbind(data.all, true.beta)

data.cp <- data.cp %>%
  mutate(bias.gw.beta11 = abs(gw.coef11 - beta11),
         bias.gw.beta12 = abs(gw.coef12 - beta12),
         bias.gw.beta21 = abs(gw.coef21 - beta21),
         bias.gw.beta22 = abs(gw.coef22 - beta22),
         cp.gw.beta11 = ((gw.coef11 - se.gw.coef11*1.96) < beta11) & (beta11 < (gw.coef11 + se.gw.coef11*1.96)),
         cp.gw.beta12 = ((gw.coef12 - se.gw.coef12*1.96) < beta12) & (beta12 < (gw.coef12 + se.gw.coef12*1.96)),
         cp.gw.beta21 = ((gw.coef21 - se.gw.coef21*1.96) < beta21) & (beta21 < (gw.coef21 + se.gw.coef21*1.96)),
         cp.gw.beta22 = ((gw.coef22 - se.gw.coef22*1.96) < beta22) & (beta22 < (gw.coef22 + se.gw.coef22*1.96)),
         bias.glb.beta11 = abs(glb.coef11 - beta11),
         bias.glb.beta12 = abs(glb.coef12 - beta12),
         bias.glb.beta21 = abs(glb.coef21 - beta21),
         bias.glb.beta22 = abs(glb.coef22 - beta22),
         cp.glb.beta11 = ((glb.coef11 - se.glb.coef11*1.96) < beta11) & (beta11 < (glb.coef11 + se.glb.coef11*1.96)),
         cp.glb.beta12 = ((glb.coef12 - se.glb.coef12*1.96) < beta12) & (beta12 < (glb.coef12 + se.glb.coef12*1.96)),
         cp.glb.beta21 = ((glb.coef21 - se.glb.coef21*1.96) < beta21) & (beta21 < (glb.coef21 + se.glb.coef21*1.96)),
         cp.glb.beta22 = ((glb.coef22 - se.glb.coef22*1.96) < beta22) & (beta22 < (glb.coef22 + se.glb.coef22*1.96)))


data.cp.sum <- data.cp %>% group_by(loc) %>% 
  summarise(mean(gw.coef11),mean(gw.coef12),
            mean(gw.coef21),mean(gw.coef22),
            mean(bias.gw.beta11), mean(bias.gw.beta12),
            mean(bias.gw.beta21), mean(bias.gw.beta22), 
            mean(bias.gw.beta11^2), mean(bias.gw.beta12^2),
            mean(bias.gw.beta21^2), mean(bias.gw.beta22^2), 
            sd(gw.coef11), sd(gw.coef12), 
            sd(gw.coef21), sd(gw.coef22), 
            mean(se.gw.coef11), mean(se.gw.coef12), 
            mean(se.gw.coef21), mean(se.gw.coef22), 
            mean(cp.gw.beta11), mean(cp.gw.beta12), 
            mean(cp.gw.beta21), mean(cp.gw.beta22), 
            mean(glb.coef11),mean(glb.coef12),
            mean(glb.coef21),mean(glb.coef22),
            mean(bias.glb.beta11), mean(bias.glb.beta12), 
            mean(bias.glb.beta21), mean(bias.glb.beta22), 
            mean(bias.glb.beta11^2), mean(bias.glb.beta12^2), 
            mean(bias.glb.beta21^2), mean(bias.glb.beta22^2), 
            sd(glb.coef11), sd(glb.coef12), 
            sd(glb.coef21), sd(glb.coef22), 
            mean(se.glb.coef11), mean(se.glb.coef12),  
            mean(se.glb.coef21), mean(se.glb.coef22), 
            mean(cp.glb.beta11), mean(cp.glb.beta12), 
            mean(cp.glb.beta21), mean(cp.glb.beta22))


