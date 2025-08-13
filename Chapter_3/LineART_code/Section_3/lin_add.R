library(coda)

{
  set.seed(1234)
  
  n <- 100 #sample size
  sigma <- 1
  npred <- 500 #sample size
  p <- 5 #number of covariates
  x <- matrix(runif(n*p), ncol = p)
  xpred <- matrix(runif(npred*p), ncol = p)
  
  #Linear additive
  y <- 3*x[,1] - 2*x[,2] + 0.5*x[,3] + x[,4] - x[,5]  + 1.5*x[,1]*x[,2] - 2*x[,3]*x[,4] + 0.8*x[,2]*x[,5] + rnorm(n, sd = sigma)
  ypred <- 3*xpred[,1] - 2*xpred[,2] + 0.5*xpred[,3] + xpred[,4] - xpred[,5]  + 1.5*xpred[,1]*xpred[,2] - 2*xpred[,3]*xpred[,4] + 0.8*xpred[,2]*xpred[,5]
  
  
  
  set.seed(328645678)
  
  
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


{bcf_mu_linear <- NULL
  lbcf_mu_linear <- NULL}

for(i in 1:rep){
  bcf_mu_linear <- rbind(bcf_mu_linear, model_bart_default[[i]]$mu_post)
  lbcf_mu_linear <- rbind(lbcf_mu_linear, model_lbart_default[[i]]$mu_post)
}
options(scipen=20)

{
  cd_linear <- mcmc(lbcf_mu_linear)
  cd2_linear <- mcmc(bcf_mu_linear)
  ef1_linear <- effectiveSize(cd_linear)
  ef2_linear <- effectiveSize(cd2_linear)
  par(mfrow=c(1,1))
  pdf("ess_linear.pdf")
  plot(ef1_linear, ef2_linear, pch=19, ylim = c(0,100000), xlim = c(0,100000),
       col = "blue", cex = 0.7,
       xlab = "Proposed Prior", ylab = "Simple prior",
       main = "Effective size comparison")
  abline(0,1, col = "red",lty = 2, lwd = 2)
  dev.off()
}


