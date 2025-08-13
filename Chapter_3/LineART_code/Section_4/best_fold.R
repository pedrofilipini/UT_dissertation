library(magrittr)
library(ggplot2)
library(purrr)
library(dplyr)
setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/5fold")

#The datasets
current <- 1:41
#Reps
reps <- 20
#folds
folds <- 5
#Model settings
settings <- list()
#Results
results <- list()

df_final <- NULL

for(i in current){
  #Auxiliary variable so I know which columns are which
  aux <- 0
  for(hh in 1:reps){
    df <- NULL
    for(alpha in c(0.5, 0.8, 0.95)){
      for(beta in c(1, 2, 3)){
        for(nt in c(10, 50, 100, 150, 200)){
          #Updates auxiliary variable
          aux <- aux+1  
          
          if(i==1){
            #This creates a list of settings, just so I know which columns are which
            settings[[aux]] <- c(alpha, beta, nt)
          }
          #RMSE
          
          
          
          vec_aux <- NULL
          for(k in 1:folds){
            val_aux <- read.csv(paste0("rmse_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"__k_",k,".csv"), sep = ";", header = T) %>% as.numeric()
            vec_aux <- c(vec_aux, val_aux)
          }
          rmse_avg <- mean(vec_aux)
          
          df <- rbind(df, data.frame(rmse = rmse_avg, dataset = i, rep = hh, a = alpha, b = beta, nt = nt))
          
        }
      }
    }
    df_aux <- df %>% slice(which.min(rmse))
    df_final <- rbind(df_final, df_aux)
    
  }
}

#Saving the object so I do not have to run this again, since there are 200k files
setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/complete")
save(df_final, file = "5fold_cv.RData")


load("5fold_cv.RData")

df_cv <- NULL
for(i in 1:dim(df_final)[1]){
  aux <- df_final[i,] 
  j <- aux$dataset
  nt <- aux$nt
  alpha <- aux$a
  beta <- aux$b
  k <- aux$rep
  vec_aux <- read.csv(paste0("rmse_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",j,"_total.csv"), sep = ";", header = F)[-1,2] %>% as.numeric()
  df_cv <- rbind(df_cv, data.frame(rmse = vec_aux[k], dataset = j, rep = k, a = alpha, b = beta, nt = nt))
}

save(df_cv, file = "df_cv.RData")











