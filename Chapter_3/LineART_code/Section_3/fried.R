library(coda)

{
  set.seed(123)
  
  n <- 100 #sample size
  npred <- 500 #sample size
  sigma <- 1
  p <- 5 #number of covariates
  x <- matrix(runif(n*p), ncol = p)
  xpred <- matrix(runif(npred*p), ncol = p)
  
  #Friedman
  y <- 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2 + 10*x[,4] + 5*x[,5] + rnorm(n, sd = sigma)
  ypred <- 10*sin(pi*xpred[,1]*xpred[,2]) + 20*(xpred[,3]-0.5)^2 + 10*xpred[,4] + 5*xpred[,5]
  
  
  
  set.seed(328645)
  
  
  burn <- 50000
  th <- 1
  post <- 5000
  
  model_lbart_default <- list()
  model_bart_default <- list()
  
  rep <- 20
  
  for(i in 1:rep){
    set.seed(2135+i)
    nt <- 50
    alpha <- 0.5
    beta <- 1
    lbart_default <- lbart:::bcf_new(y, #response
                                     rep(1,n), #Every z equals 1 because we want regular BART
                                     x_control = as.matrix(x), #matrix of covariates
                                     pihat = rep(1,n), #in case of causal inference
                                     z_est = rep(1,npred), x_control_est = as.matrix(xpred), pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
                                     include_pi = "none", #can be "none", "control", "moderate", or "both"
                                     linear = "control", #can be "none", "control", "moderate", or "both"
                                     base_control = alpha, #base_moderate = 0.95, #hyperparameters
                                     power_control = beta, #power_moderate = 2, #hyperparameters
                                     nburn = burn, nthin = th, nsim = post, #draws
                                     ntree_control = nt, #control trees
                                     ntree_moderate = 0, #moderate trees (treatment effect)
                                     dart = F, #Linero's Prior
                                     save_trees = T) #Not saving the trees
    model_lbart_default[[i]] <- lbart_default
    
    nt <- 50
    alpha <- 0.5
    beta <- 1
    bart_default <- lbart:::bcf_new(y, #response
                                    rep(1,n), #Every z equals 1 because we want regular BART
                                    x_control = as.matrix(x), #matrix of covariates
                                    pihat = rep(1,n), #in case of causal inference
                                    z_est = rep(1,npred), x_control_est = as.matrix(xpred), pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
                                    include_pi = "none", #can be "none", "control", "moderate", or "both"
                                    linear = "control", #can be "none", "control", "moderate", or "both"
                                    base_control = alpha, #base_moderate = 0.95, #hyperparameters
                                    power_control = beta, #power_moderate = 2, #hyperparameters
                                    nburn = burn, nthin = th, nsim = post, #draws
                                    ntree_control = nt, #control trees
                                    ntree_moderate = 0, #moderate trees (treatment effect)
                                    dart = F, #Linero's Prior
                                    gcon = 0, #1 to use the gprior, 0 to regular bart prior on drmu_linear
                                    gmod = 0, #1 to use the gprior, 0 to regular bart prior on drmu_linear
                                    save_trees = T) #Not saving the trees
    model_bart_default[[i]] <- bart_default
  }
}


{bcf_mu <- NULL
  lbcf_mu <- NULL}

for(i in 1:rep){
  bcf_mu <- rbind(bcf_mu, model_bart_default[[i]]$mu_post)
  lbcf_mu <- rbind(lbcf_mu, model_lbart_default[[i]]$mu_post)
}


{
  cd <- mcmc(lbcf_mu)
  cd2 <- mcmc(bcf_mu)
  ef1 <- effectiveSize(cd)
  ef2 <- effectiveSize(cd2)
  par(mfrow=c(1,1))
  pdf("ess_fried.pdf")
  plot(ef1, ef2, pch=19, ylim = c(0,8000), xlim = c(0,8000),
       col = "blue", cex = 0.7,
       xlab = "Proposed Prior", ylab = "Simple prior",
       main = "Effective size comparison")
  abline(0,1, col = "red",lty = 2, lwd = 2)
  dev.off()
}
# 
# {
#   k <- 87
#   plot.ts(lbart_default$mu_est_post[,k])
#   plot.ts(bart_default$mu_est_post[,k])
#   acf(lbart_default$mu_est_post[,k])
#   acf(bart_default$mu_est_post[,k])
#   
#   ypred[k]
# }
# 
# cd <- mcmc(lbart_default$mu_est_post)
# cd2 <- mcmc(bart_default$mu_est_post)
# ef1 <- effectiveSize(cd)
# ef2 <- effectiveSize(cd2)
# 
# par(mfrow=c(1,1))
# plot(ef1, ef2,
#      xlab = "Proposed Prior", ylab = "Simple prior",
#      main = "Effective size comparison")
# abline(0,1)
# {
#   cd <- mcmc(lbart_default$mu_est_post)
#   cd2 <- mcmc(bart_default$mu_est_post)
#   ef1 <- effectiveSize(cd)
#   ef2 <- effectiveSize(cd2)
#   par(mfrow=c(1,2))
#   plot(xpred[,1], ef1, ylim=c(0,500), 
#        ylab = "Effective size", main = "Proposed prior")
#   plot(xpred[,1], ef2, ylim=c(0,500), 
#        ylab = "Effective size", main = "Simple prior")
#   plot(xpred[,2], ef1, ylim=c(0,500), 
#        ylab = "Effective size", main = "Proposed prior")
#   plot(xpred[,2], ef2, ylim=c(0,500), 
#        ylab = "Effective size", main = "Simple prior")
#   plot(xpred[,3], ef1, ylim=c(0,500), 
#        ylab = "Effective size", main = "Proposed prior")
#   plot(xpred[,3], ef2, ylim=c(0,500), 
#        ylab = "Effective size", main = "Simple prior")
#   plot(xpred[,4], ef1, ylim=c(0,500), 
#        ylab = "Effective size", main = "Proposed prior")
#   plot(xpred[,4], ef2, ylim=c(0,500), 
#        ylab = "Effective size", main = "Simple prior")
#   plot(xpred[,5], ef1, ylim=c(0,500), 
#        ylab = "Effective size", main = "Proposed prior")
#   plot(xpred[,5], ef2, ylim=c(0,500), 
#        ylab = "Effective size", main = "Simple prior")
# }


sqrt(mean((ypred-apply(lbart_default$mu_est_post, 2, mean))^2))
sqrt(mean((ypred-apply(bart_default$mu_est_post, 2, mean))^2))

