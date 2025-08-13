#Cross validation

{
  setwd("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart/cv_tpr_small")
  
  library(dplyr)
  library(ggplot2)
  library(latex2exp)
  library(foreach)
  library(doParallel)
  
  source("C:/Users/pedro/OneDrive/Área de Trabalho/Important_so_doing_it_again/Nature/Project_Stress_Final/Study_1/Code/Auxiliary/gp_approx_fun.R")
  
  data_1 <- read.csv("C:/Users/pedro/OneDrive/Área de Trabalho/Important_so_doing_it_again/Nature/Project_Stress/New Dataset/bcfstudy1_long.csv")
  
  #Adding a variable that is the observed mean for "math" for each individual
  
  data_1$tpr[data_1$tpr<650] = NA
  data_1$tpr[data_1$tpr>5000] = NA
  
  mean_math <- data_1[data_1$epoch=="math",] %>%
    group_by(id) %>%
    summarise(mean_tpr_math = mean(tpr, na.rm = T),
              mean_bps_math = mean(bps, na.rm = T))
  mean_speech <- data_1[data_1$epoch=="speech",] %>%
    group_by(id) %>%
    summarise(mean_tpr_speech = mean(tpr, na.rm = T),
              mean_bps_speech = mean(bps, na.rm = T))
  
  #Adding new columns in the data
  data_1 <- left_join(data_1, mean_math, by ="id")
  data_1 <- left_join(data_1, mean_speech, by ="id")
  
  
  data_1$sex <- (data_1$sex-1) %>% as.factor()
  
  # Person-level dataframe for everyone
  person_df <- data_1 %>% group_by(id) %>% filter(row_number()==1) %>%group_by()
  
  
  get_controls = function(x) {
    
    ret = x %>% dplyr::select(bothmindsets,
                              fixedmindset_baseline,
                              stressmindset_baseline,
                              selfesteem_baseline, 
                              sex) %>%
      dbarts::makeModelMatrixFromDataFrame(drop = FALSE)
  }
  
  get_moderators = function(x) {
    x %>% dplyr::select(bothmindsets,
                        fixedmindset_baseline,
                        stressmindset_baseline,
                        sex) %>%
      dbarts::makeModelMatrixFromDataFrame(drop = FALSE)
  }
  
  
  #Get response (still not clean)
  y_math <- person_df$mean_tpr_math
  y_speech <- person_df$mean_tpr_speech
  
  #Get covariates
  x_con_person <- get_controls(person_df)
  x_mod_person <- get_moderators(person_df)
  #Remove where y have NA's
  x_con_math <- x_con_person[!is.na(y_math),]
  x_mod_math <- x_mod_person[!is.na(y_math),]
  x_con_speech <- x_con_person[!is.na(y_speech),]
  x_mod_speech <- x_mod_person[!is.na(y_speech),]
  
  #The treatment
  z_math <- person_df$treatment[!is.na(y_math)]
  z_speech <- person_df$treatment[!is.na(y_speech)]
  
  #The response without NA's
  y_math <- y_math[!is.na(y_math)]
  y_speech <- y_speech[!is.na(y_speech)]
  
  
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
npred <- length(y_speech)
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
    idx <- nsamp[(29*(i-1)+1):(29*i)] 
    y <- y_speech[-idx]
    yp <- y_speech[idx]
    z <- z_speech[-idx]
    zp <- z_speech[idx]
    xc <- x_con_speech[-idx,]
    xcp <- x_con_speech[idx,]
    xm <- x_mod_speech[-idx,]
    xmp <- x_mod_speech[idx,]
    np <- length(idx)
    n <- npred-np
    
    
    set.seed(1)
    #simtime <- Sys.time()
    lbcf_speech1 <- lbart:::bcf_new(y, #response
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
    r1_train <- sqrt(mean((apply(lbcf_speech1$mu_post+t(t(lbcf_speech1$b_post)*z), 2, mean)-y)^2))
    r1 <- sqrt(mean((apply(lbcf_speech1$mu_est_post+t(t(lbcf_speech1$b_est_post)*zp), 2, mean)-yp)^2))
    
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
    
    #write.table(time, paste0("times_cv_speech__a_", cbase, "__am_",mbase, "__b_",cpower,"__bm_",mpower,"__nt_",ctree,"__ntm_",mtree,"__mu_",mm,"__tau_",tt,"__linear_",linear,"__fold_",i,"_partial.csv"), sep = ";")
    write.table(r1_train, paste0("rmsetrain_cv_speech__a_", cbase, "__am_",mbase, "__b_",cpower,"__bm_",mpower,"__nt_",ctree,"__ntm_",mtree,"__mu_",mm,"__tau_",tt,"__linear_",linear,"__fold_",i,"_partial.csv"), sep = ";")
    write.table(r1, paste0("rmse_cv_speech__a_", cbase, "__am_",mbase, "__b_",cpower,"__bm_",mpower,"__nt_",ctree,"__ntm_",mtree,"__mu_",mm,"__tau_",tt,"__linear_",linear,"__fold_",i,"_partial.csv"), sep = ";")
    
    #Clean memory
    rm(lbcf_speech1)
    gc()
    
  }
  
}
parallel::stopCluster(cl)
