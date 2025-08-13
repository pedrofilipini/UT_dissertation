#setwd("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart")

#Softbart broke at dataset 24, sim 13
#dataset3, sim 18 too


{
  library(foreach)
  library(magrittr)
  library(SoftBart)
  #Where the datasets are located
  setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/nipsdata")
  
  #Repetitions of the simulations
  simrep <- 20
  
  
  #These are the large datasets that take "forever" to run (a few days)
  #Best thing is to use 20 threads to parallelize these separately from everything else
  large <- c(22,1,3,4,10,24,26,41)
  
  
  #Reading files
  x_files <- list.files(path=getwd(), pattern = "x.txt", full.names=TRUE, recursive=TRUE)
  y_files <- list.files(path=getwd(), pattern = "y.txt", full.names=TRUE, recursive=TRUE)
  
  #Saving those in lists so its gonna be easier to make a for loop
  x <- list()
  y <- list()
  x <- lapply(x_files, function(x){read.table(x, header=F)})
  y <- lapply(y_files, function(x){read.table(x, header=F)})
  
  #Just to compare sizes and variables with BART paper (LABOR is missing, so 41 datasets)
  size_x <- matrix(unlist(lapply(x, dim)), ncol = 2, byrow = T)
  size_y <- matrix(unlist(lapply(y, dim)), ncol = 2, byrow = T)[,1]
  
  
  #Fixing categorical variables
  sizes <- readr::read_csv("sizes_paper.csv")
  
  sz <- sizes$sz
  
  string_to_vector <- function(input_string) {
    # Remove parentheses and split the string by commas
    numbers <- strsplit(gsub("[()]", "", input_string), ",")[[1]]
    
    # Convert the character vector to numeric
    result_vector <- as.numeric(numbers)
    
    return(result_vector)
  }
  
  lista_sz <- lapply(sz, function(x){string_to_vector(x)})
  
  set.seed(9465)
  
  times_table <- list()
  
  for(i in c(large,(1:length(lista_sz))[-large])){
    
    aa <- apply(x[[i]],2,function(x){sum(length(unique(x)))})
    bb <- lista_sz[[i]]
    cc <- match(aa, bb, nomatch = 0)
    ii <- which(cc!=0)
    
    if(i==25){
      xx <- x[[i]] %>% dbarts::makeModelMatrixFromDataFrame(drop = FALSE)
    }
    
    for(j in ii){
      x[[i]] <- cbind(x[[i]] ,dbarts::makeModelMatrixFromDataFrame(data.frame(j=as.factor(x[[i]][,j])), drop = F))
    }
    
    for(j in ii){
      x[[i]][,j] <- as.factor(x[[i]][,j])
    }
    
    if(i==15){
      j <- 13
      x[[i]] <- cbind(x[[i]] ,dbarts::makeModelMatrixFromDataFrame(data.frame(j=as.factor(x[[i]][,j])), drop = F))
      x[[i]][,j] <- as.factor(x[[i]][,j])
    }
    
    if(length(ii)>0){
      if(i==15){
        x[[i]] <- x[[i]][,-c(ii,13)]
      }else{
        x[[i]] <- x[[i]][,-ii]
      }
    }
    
    for(jjj in 1:dim(x[[i]])[2]){
      if((length(unique(x[[i]][,jjj]))>6)|((length(unique(x[[i]][,jjj]))>6)&(max(abs(unique(x[[i]][,jjj])))>100))){
        x[[i]][,jjj] <- scale(x[[i]][,jjj])
      }
    }
    
    #These are the datasets
    if(i!=25){
      xx <- x[[i]] %>% dbarts::makeModelMatrixFromDataFrame(drop = FALSE)
    }
    
    yy <- as.numeric(unlist(y[[i]]))
    
    #This will be the size of the training
    n_train <- round(length(yy)*5/6,0)
    n_test <- length(yy)-n_train
    
    times_table[[i]] <- matrix(0, nrow = simrep, ncol = 3)
    
    write.table(0, paste0("start",i,".csv"), sep = ";")
    
    #On BART paper is hh 1:20
    for(alpha in c(0.5, 0.8, 0.95)){ 
      for(beta in c(1, 2, 3)){
        for(nt in c(10,50,100,150,200)){
          cl <- parallel::makeCluster(20)
          doParallel::registerDoParallel(cl)
          foreach(hh = 1:simrep, .packages=c('lbart','dbarts','dplyr','magrittr')) %dopar%  {
            
            
            #indexes of training
            set.seed(1236*hh+356*i)
            idx <- sample(1:length(yy), size = n_train, replace = F)
            
            #This creates a vector with the samples for k-fold
            kfold <- 5
            kidx <- rep(1:kfold, times = ceiling(length(idx)/kfold))[1:length(idx)]
            
            for(k in 1:kfold){
              #Getting the indexes
              idx_train <- idx[kidx!=k]
              idx_val <- idx[kidx==k]
              
              #Train
              x_train <- xx[idx_train,]
              y_train <- yy[idx_train]
              
              #Test
              x_test <- xx[idx_val,]
              y_test <- yy[idx_val]
              
              
              x_train_min <- apply(x_train,2,min)
              x_train_max <- apply(x_train,2,max)
              
              for(kkk in 1:dim(x_test)[2]){
                x_test[(x_test[,kkk]<x_train_min[kkk]),kkk] <- x_train_min[kkk]
                x_test[(x_test[,kkk]>x_train_max[kkk]),kkk] <- x_train_max[kkk]
              }
              
              
              #From BART paper
              burn <- 10000
              post <- 5000
              th <- 3

              #They only use the RRMSE to measure the performance
              times_vec <- numeric(0)
              rmse_train <- numeric(0)
              rmse <- numeric(0)
              coverage <- numeric(0)
              coveragep <- numeric(0)
              range <- numeric(0)
              rangep <- numeric(0)
              cl_names <- character(0) #This is just to identify whatever
              
              print(paste0("Dataset: ", i, "; Simulation: ", hh, "; k: ", k, "; a: ", alpha, "; b: ", beta, "; nt: ", nt, "; DEFAULT LBART;"))
              
              simtime <- Sys.time()
              set.seed(1236*hh+356*i)
              #default lbart (a=0.5, b=1, m=50)
              lbart_default <- lbart:::bcf_new(y_train, #response
                                               rep(1,length(idx_train)), #Every z equals 1 because we want regular BART
                                               x_control = x_train, #matrix of covariates
                                               pihat = rep(1,length(idx_train)), #in case of causal inference
                                               z_est = rep(1,length(idx_val)), x_control_est = x_test, pihat_est = rep(1,length(idx_val)), #for predictions (does not matter by now, not implemented)
                                               include_pi = "none", #can be "none", "control", "moderate", or "both"
                                               linear = "control", #can be "none", "control", "moderate", or "both"
                                               base_control = alpha, #base_moderate = 0.95, #hyperparameters
                                               power_control = beta, #power_moderate = 2, #hyperparameters
                                               nburn = burn, nthin = th, nsim = post, #draws
                                               ntree_control = nt, #control trees
                                               ntree_moderate = 0, #moderate trees (treatment effect)
                                               dart = F, #Linero's Prior
                                               save_trees = F) #Not saving the trees
              
              times_lbart <- difftime(Sys.time(), simtime, units = "secs")[[1]]
              
              sqrt(mean((apply(lbart_default$mu_est_post, 2, mean)-y_test)^2))
              
              times1 <- times_lbart
              
              r1_train <- sqrt(mean((apply(lbart_default$mu_post, 2, mean)-y_train)^2))
              
              r1 <- sqrt(mean((apply(lbart_default$mu_est_post, 2, mean)-y_test)^2))
              
              cove1 <- mean(dplyr::between(y_train,apply(lbart_default$mu_post, 2, quantile, 0.025), apply(lbart_default$mu_post, 2, quantile, 0.975)))
              
              covep1 <- mean(dplyr::between(y_test,apply(lbart_default$mu_est_post, 2, quantile, 0.025), apply(lbart_default$mu_est_post, 2, quantile, 0.975)))
              
              rangep1 <- mean(apply(lbart_default$mu_est_post, 2, quantile, 0.975)-apply(lbart_default$mu_est_post, 2, quantile, 0.025))
              
              range1 <- mean(apply(lbart_default$mu_post, 2, quantile, 0.975)-apply(lbart_default$mu_post, 2, quantile, 0.025))
              
              times_vec <- c(times_vec, times1)
              rmse_train <- c(rmse_train, r1_train)
              rmse <- c(rmse, r1)
              coverage <- c(coverage, cove1)
              coveragep <- c(coveragep, covep1)
              range <- c(range, range1)
              rangep <- c(rangep, rangep1)
              
              
              #Clean space
              rm(lbart_default)
              gc()
              
              
              #For large datasets it is faster to parallelize the repetitions
              #So I save partial results and concatenate everything after
              
              
              write.table(times_vec, paste0("times_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"__k_",k,".csv"), sep = ";")
              write.table(rmse_train, paste0("rmse_train_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"__k_",k,".csv"), sep = ";")
              write.table(rmse, paste0("rmse_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"__k_",k,".csv"), sep = ";")
              write.table(coverage, paste0("coverage_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"__k_",k,".csv"), sep = ";")
              write.table(coveragep, paste0("coveragep_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"__k_",k,".csv"), sep = ";")
              write.table(range, paste0("range_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"__k_",k,".csv"), sep = ";")
              write.table(rangep, paste0("rangep_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"__k_",k,".csv"), sep = ";")
              
            } 
          }
          parallel::stopCluster(cl)
          
        }
      }
    }
    
  }
  
}













