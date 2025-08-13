library(magrittr)

setwd("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart/cv_tpr_small")

folds <- 5
df_speech <- NULL

for(mtree in c(25, 50, 100)){
  for(ctree in c(50, 100, 200)){
    for(cbase in c(0.5, 0.8, 0.95)){
      for(mbase in c(0.25, 0.5, 0.8)){
        for(cpower in c(1, 2, 3)){
          for(mpower in c(1, 2, 3)){
            for(muscale in c(T, F)){
              for(tauscale in c(T, F)){
                for(linear in c("control", "moderate", "both", "none")){
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
                  
                  rmse_train <- numeric(0)
                  rmse_test <- numeric(0)
                  for(i in 1:folds){
                    rmse_train[i] <- read.csv(paste0("rmsetrain_cv_speech__a_", cbase, "__am_",mbase, "__b_",cpower,"__bm_",mpower,"__nt_",ctree,"__ntm_",mtree,"__mu_",mm,"__tau_",tt,"__linear_",linear,"__fold_",i,"_partial.csv"), sep = ";") %>% as.numeric()
                    rmse_test[i] <- read.csv(paste0("rmse_cv_speech__a_", cbase, "__am_",mbase, "__b_",cpower,"__bm_",mpower,"__nt_",ctree,"__ntm_",mtree,"__mu_",mm,"__tau_",tt,"__linear_",linear,"__fold_",i,"_partial.csv"), sep = ";") %>% as.numeric()
                    
                    
                  }
                  rmse_avg_train <- mean(rmse_train)
                  rmse_sd_train <- sd(rmse_train)
                  rmse_avg_test <- mean(rmse_test)
                  rmse_sd_test <- sd(rmse_test)
                  df_speech <- rbind(df_speech, data.frame(rmse_train = rmse_avg_train, rmse_test = rmse_avg_test, rmse_train_sd = rmse_sd_train, rmse_test_sd = rmse_sd_test, a_con = cbase, a_mod = mbase, b_con = cpower, b_mod = mpower, nt_con = ctree, nt_mod = mtree, muscale = mm, tauscale = tt, linear = linear))
                  
                  
                }
              }
            }
          }
        }
      }
    }
  }
}

#save(df_speech, file = "df_cvspeech.RData")

load("df_cvspeech.RData")

