#Cross validation

{
  setwd("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart/cv_anxiety")
  library(foreach)
  library(doParallel)
  library(dplyr)
  library(magrittr)
  
  #reading data
  load("C:/Users/pedro/Downloads/study6.rdta")
  
  #Cleaning up NA in treatment and anxiety
  study6_clean <- study6[!is.na(study6$treatment),]
  study6_clean <- study6_clean[!is.na(study6_clean$anxiety),]
  
  #control n=172; treatment n=179, n=351 
  study6_clean <- study6_clean %>% dplyr::select(treatment, anxiety,
                                                 testanxiety, sex, 
                                                 stressmindset_baseline, 
                                                 fixedmindset_baseline, pss_impute, 
                                                 bfi_c, bfi_o, bfi_e, bfi_n,
                                                 fgen_moth, fgen_fath, ses)
  
  #now n is smaller, only 331
  study6_clean <- study6_clean %>% filter(complete.cases(.))
  
  Xc <- study6_clean %>% dplyr::select(testanxiety, sex, stressmindset_baseline, fixedmindset_baseline, pss_impute, 
                                       bfi_c, bfi_o, bfi_e, bfi_n,
                                       fgen_moth, fgen_fath, ses) %>% dbarts::makeModelMatrixFromDataFrame()
  
  
  Xm <- study6_clean %>% dplyr::select(sex, stressmindset_baseline, fixedmindset_baseline, pss_impute )%>% dbarts::makeModelMatrixFromDataFrame()
  
  
  z_main = study6_clean$treatment # "Treatment" variable - Actually a measured variable
  y_main = study6_clean$anxiety
  
}

#Accounting for all the possible hyperparameters possibilities

iter <- 0
hyper <- list(NULL)

for(mtree in c(25, 50, 100)){
  for(ctree in c(50, 100, 200)){
    for(cbase in c(0.5, 0.8, 0.95)){
      for(mbase in c(0.25, 0.5, 0.8)){
        for(cpower in c(1, 2, 3)){
          for(mpower in c(1, 2, 3)){
            for(muscale in c(T, F)){
              for(tauscale in c(T, F)){
                for(linear in c("control", "moderate", "both", "none")){
                  iter <- iter + 1
                  df <- data.frame(mtree, ctree, cbase, mbase, cpower, mpower, muscale, tauscale, linear)
                  hyper[[iter]] <- df
                }
              }
            }
          }
        }
      }
    }
  }
}


#5 fold cross validation 
npred <- length(y_main)
id_aux <- rep(1:5, times = ceiling(npred/5))[1:npred]

set.seed(1234567)
nsamp <- sample(1:npred, size = npred, replace = F)

folds <- 5

burn <- 1000
th <- 1
post <- 1000

cl <- parallel::makeCluster(20)
doParallel::registerDoParallel(cl)
foreach(k = 1:length(hyper), .packages=c('lbart')) %dopar%  {
  #for(k in 1:length(hyper)){
  df <- hyper[[k]]
  
  mtree <- df$mtree
  ctree <- df$ctree
  cbase <- df$cbase
  mbase <- df$mbase
  cpower <- df$cpower
  mpower <- df$mpower
  muscale <- df$muscale
  tauscale <- df$tauscale
  linear <- df$linear
  
  #for each possibility, do it 5 times
  for(i in 1:5){
    idx <- (1:npred)[(id_aux[nsamp]==i)]
    y <- y_main[-idx]
    yp <- y_main[idx]
    z <- z_main[-idx]
    zp <- z_main[idx]
    xc <- Xc[-idx,]
    xcp <- Xc[idx,]
    xm <- Xm[-idx,]
    xmp <- Xm[idx,]
    np <- length(idx)
    n <- npred-np
    
    
    set.seed(1)
    #simtime <- Sys.time()
    lbcf_anx1 <- lbart:::bcf_new(y, #response
                                    z = z, z_est = zp, #Z
                                    x_control = xc, #matrix of covariates
                                    x_control_est = xcp,
                                    x_moderate = xm,
                                    x_moderate_est = xmp,
                                    pihat = rep(1,n), #in case of causal inference
                                    pihat_est = rep(1,np), #for predictions (does not matter by now, not implemented)
                                    include_pi = "none", #can be "none", "control", "moderate", or "both"
                                    linear = linear, #can be "none", "control", "moderate", or "both"
                                    base_control = cbase, base_moderate = mbase, #hyperparameters for control
                                    power_control = cpower, power_moderate = mpower, #hyperparameters for moderate
                                    nburn = burn, nthin = th, nsim = post, #draws
                                    ntree_control = ctree, #control trees
                                    ntree_moderate = mtree, #moderate trees (treatment effect)
                                    dart = F, #Linero's Prior
                                    save_trees = T, use_muscale = muscale, use_tauscale = tauscale) #Not saving the trees
    
    #time <- difftime(Sys.time(), simtime, units = "secs")[[1]]
    r1_train <- sqrt(mean((apply(lbcf_anx1$mu_post+t(t(lbcf_anx1$b_post)*z), 2, mean)-y)^2))
    r1 <- sqrt(mean((apply(lbcf_anx1$mu_est_post+t(t(lbcf_anx1$b_est_post)*zp), 2, mean)-yp)^2))
    
    #Save stuff after
    if(muscale){
      mm <- "true"
    }else{
      mm <- "false"
    }
    if(tauscale){
      tt <- "true"
    }else{
      tt <- "false"
    }
    
    #write.table(time, paste0("times_cv_anx__a_", cbase, "__am_",mbase, "__b_",cpower,"__bm_",mpower,"__nt_",ctree,"__ntm_",mtree,"__mu_",mm,"__tau_",tt,"__linear_",linear,"__fold_",i,"_partial.csv"), sep = ";")
    write.table(r1_train, paste0("rmsetrain_cv_anx__a_", cbase, "__am_",mbase, "__b_",cpower,"__bm_",mpower,"__nt_",ctree,"__ntm_",mtree,"__mu_",mm,"__tau_",tt,"__linear_",linear,"__fold_",i,"_partial.csv"), sep = ";")
    write.table(r1, paste0("rmse_cv_anx__a_", cbase, "__am_",mbase, "__b_",cpower,"__bm_",mpower,"__nt_",ctree,"__ntm_",mtree,"__mu_",mm,"__tau_",tt,"__linear_",linear,"__fold_",i,"_partial.csv"), sep = ";")
    
    #Clean memory
    rm(lbcf_anx1)
    gc()
    
  }
  
}
parallel::stopCluster(cl)
