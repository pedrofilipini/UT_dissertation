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
  
  
  #set.seed(23864)
  set.seed(9465)
  
  #i is the dataset i
  #for(i in 1:length(lista_sz))
  
  
  
  times_table <- list()
  
  #cl <- parallel::makeCluster(28)
  #doParallel::registerDoParallel(cl)
  
  #foreach(i = 1:length(lista_sz), .packages=c('lbart','dbarts','dplyr','magrittr')) %dopar%  {
    
    #Don't do 1, 3, 4, 10, 22 and 42 for now.
  for(i in c(large,(1:length(lista_sz))[-large])){
    #for(i in 3)
    #       i <- 3
    #{
    
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
      #for(j in 1:4){
      x[[i]][,j] <- as.factor(x[[i]][,j])
    }
    
    if(i==15){
      j <- 13
      x[[i]] <- cbind(x[[i]] ,dbarts::makeModelMatrixFromDataFrame(data.frame(j=as.factor(x[[i]][,j])), drop = F))
      x[[i]][,j] <- as.factor(x[[i]][,j])
    }
    
    # if(i==35){
    #   j <- 1
    #   x[[i]] <- cbind(x[[i]] ,dbarts::makeModelMatrixFromDataFrame(data.frame(j=as.factor(x[[i]][,j])), drop = F))
    #   x[[i]][,j] <- as.factor(x[[i]][,j])
    #   x[[i]] <- x[[i]][,-c(ii,1)]
    # }
    
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
            
            #Train
            x_train <- xx[idx,]
            y_train <- yy[idx]
            
            #Test
            x_test <- xx[-idx,]
            y_test <- yy[-idx]
            
            
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
            #nt <- 200 #50
            
            #They only use the RRMSE to measure the performance
            times_vec <- numeric(0)
            rmse_train <- numeric(0)
            rmse <- numeric(0)
            coverage <- numeric(0)
            coveragep <- numeric(0)
            range <- numeric(0)
            rangep <- numeric(0)
            cl_names <- character(0) #This is just to identify whatever
            
            print(paste0("Dataset: ", i, "; Simulation: ", hh, "; DEFAULT LBART;"))
            
            simtime <- Sys.time()
            set.seed(1236*hh+356*i)
            #default lbart (a=0.5, b=3, m=50)
            lbart_default <- lbart:::bcf_new(y_train, #response
                                             rep(1,n_train), #Every z equals 1 because we want regular BART
                                             x_control = x_train, #matrix of covariates
                                             pihat = rep(1,n_train), #in case of causal inference
                                             z_est = rep(1,n_test), x_control_est = x_test, pihat_est = rep(1,n_test), #for predictions (does not matter by now, not implemented)
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
            
            # print(paste0("Dataset: ", i, "; Simulation: ", hh, "; DEFAULT BART;"))
            # 
            # 
            # simtime <- Sys.time()
            # set.seed(1236*hh+356*i)
            # #default bart (a=0.95, b=2, m=200)
            # bart_default <- lbart:::bcf_new(y_train, #response
            #                                 rep(1,n_train), #Every z equals 1 because we want regular BART
            #                                 x_control = x_train, #matrix of covariates
            #                                 pihat = rep(1,n_train), #in case of causal inference
            #                                 z_est = rep(1,n_test), x_control_est = x_test, pihat_est = rep(1,n_test), #for predictions (does not matter by now, not implemented)
            #                                 include_pi = "none", #can be "none", "control", "moderate", or "both"
            #                                 linear = "none", #can be "none", "control", "moderate", or "both"
            #                                 base_control = 0.95, #base_moderate = 0.95, #hyperparameters
            #                                 power_control = 2, #power_moderate = 2, #hyperparameters
            #                                 nburn = burn, nthin = th, nsim = post, #draws
            #                                 ntree_control = 200, #control trees
            #                                 ntree_moderate = 0, #moderate trees (treatment effect)
            #                                 dart = F, #Linero's Prior
            #                                 save_trees = T) #Not saving the trees
            # 
            # times_table[[i]][hh,2] <- difftime(Sys.time(), simtime, units = "secs")[[1]]
            # sqrt(mean((apply(bart_default$mu_est_post, 2, mean)-y_test)^2))
            
            # print(paste0("Dataset: ", i, "; Simulation: ", hh, "; DEFAULT SOFTBART;"))
            
            # simtime <- Sys.time()
            # if((i==24) & (hh==13)){
            #   set.seed(1236*hh+356*i+4367)
            # }else{
            #   set.seed(1236*hh+356*i)
            # }
            # #default softbart
            # soft_default <- softbart(X = x_train, Y = y_train, X_test = x_test,
            #                         opts = Opts(num_burn = burn, num_thin = th, num_save = post))
            # 
            # times_table[[i]][hh,3] <- difftime(Sys.time(), simtime, units = "secs")[[1]]
            
            
            #Add softbart, randomforests, boosting, cv-lbart,
            #default bart and default lbart
            
            times1 <- times_lbart
            times2 <- 0#times_bart
            times3 <- 0#times_soft
            
            r1_train <- sqrt(mean((apply(lbart_default$mu_post, 2, mean)-y_train)^2))
            r2_train <- 0#sqrt(mean((apply(bart_default$mu_post, 2, mean)-y_train)^2))
            r3_train <- 0#sqrt(mean((apply(soft_default$y_hat_train, 2, mean)-y_train)^2))
            
            r1 <- sqrt(mean((apply(lbart_default$mu_est_post, 2, mean)-y_test)^2))
            r2 <- 0#sqrt(mean((apply(bart_default$mu_est_post, 2, mean)-y_test)^2))
            r3 <- 0#sqrt(mean((apply(soft_default$y_hat_test, 2, mean)-y_test)^2))
            
            cove1 <- mean(dplyr::between(y_train,apply(lbart_default$mu_post, 2, quantile, 0.025), apply(lbart_default$mu_post, 2, quantile, 0.975)))
            cove2 <- 0#mean(dplyr::between(y_train,apply(bart_default$mu_post, 2, quantile, 0.025), apply(bart_default$mu_post, 2, quantile, 0.975)))
            cove3 <- 0#mean(dplyr::between(y_train,apply(soft_default$y_hat_train, 2, quantile, 0.025), apply(soft_default$y_hat_train, 2, quantile, 0.975)))
            
            covep1 <- mean(dplyr::between(y_test,apply(lbart_default$mu_est_post, 2, quantile, 0.025), apply(lbart_default$mu_est_post, 2, quantile, 0.975)))
            covep2 <- 0#mean(dplyr::between(y_test,apply(bart_default$mu_est_post, 2, quantile, 0.025), apply(bart_default$mu_est_post, 2, quantile, 0.975)))
            covep3 <- 0#mean(dplyr::between(y_test,apply(soft_default$y_hat_test, 2, quantile, 0.025), apply(soft_default$y_hat_test, 2, quantile, 0.975)))
            
            rangep1 <- mean(apply(lbart_default$mu_est_post, 2, quantile, 0.975)-apply(lbart_default$mu_est_post, 2, quantile, 0.025))
            rangep2 <- 0#mean(apply(bart_default$mu_est_post, 2, quantile, 0.975)-apply(bart_default$mu_est_post, 2, quantile, 0.025))
            rangep3 <- 0#mean(apply(soft_default$y_hat_test, 2, quantile, 0.975)-apply(soft_default$y_hat_test, 2, quantile, 0.025))
            
            range1 <- mean(apply(lbart_default$mu_post, 2, quantile, 0.975)-apply(lbart_default$mu_post, 2, quantile, 0.025))
            range2 <- 0#mean(apply(bart_default$mu_post, 2, quantile, 0.975)-apply(bart_default$mu_post, 2, quantile, 0.025))
            range3 <- 0#mean(apply(soft_default$y_hat_train, 2, quantile, 0.975)-apply(soft_default$y_hat_train, 2, quantile, 0.025))
            
            times_vec <- c(times_vec, times1, times2, times3)
            rmse_train <- c(rmse_train, r1_train, r2_train, r3_train)
            rmse <- c(rmse, r1, r2, r3)
            coverage <- c(coverage, cove1, cove2, cove3)
            coveragep <- c(coveragep, covep1, covep2, covep3)
            range <- c(range, range1, range2, range3)
            rangep <- c(rangep, rangep1, rangep2, rangep3)
            
            
            #Clean space
            rm(lbart_default)
            # rm(bart_default)
            #rm(soft_default)
            gc()
            
            
            #For large datasets it is faster to parallelize the repetitions
            #So I save partial results and concatenate everything after
            
            
            write.table(times_vec, paste0("times_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")
            write.table(rmse_train, paste0("rmse_train_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")
            write.table(rmse, paste0("rmse_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")
            write.table(coverage, paste0("coverage_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")
            write.table(coveragep, paste0("coveragep_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")
            write.table(range, paste0("range_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")
            write.table(rangep, paste0("rangep_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")
            
          }
          parallel::stopCluster(cl)
          
          #Now I gather all the parts
          
          #First I create all the matrices and name the models
          cl_names <- NULL
          cl_names <- c(cl_names, paste0("LBART; a = ", alpha, "; b = ", beta, "; nt = ", nt, ";"))
          cl_names <- c(cl_names, paste0("BART; a = 0.95; b = 2; nt = 200;"))
          cl_names <- c(cl_names, paste0("SOFTBART; a = 0.95; b = 2; nt = 20;"))
          
          times_matrix <- matrix(0, ncol = length(cl_names), nrow = simrep)
          colnames(times_matrix) <- cl_names
          
          rmse_matrix_train <- matrix(0, ncol = length(cl_names), nrow = simrep)
          colnames(rmse_matrix_train) <- cl_names
          
          rmse_matrix <- matrix(0, ncol = length(cl_names), nrow = simrep)
          colnames(rmse_matrix) <- cl_names
          
          coverage_matrix <- matrix(0, ncol = length(cl_names), nrow = simrep)
          colnames(coverage_matrix) <- cl_names
          
          coveragep_matrix <- matrix(0, ncol = length(cl_names), nrow = simrep)
          colnames(coveragep_matrix) <- cl_names
          
          range_matrix <- matrix(0, ncol = length(cl_names), nrow = simrep)
          colnames(range_matrix) <- cl_names
          
          rangep_matrix <- matrix(0, ncol = length(cl_names), nrow = simrep)
          colnames(rangep_matrix) <- cl_names
          
          
          for(hh in 1:simrep){
            times_matrix[hh,] <- read.csv(paste0("times_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";") %>% t() %>% as.numeric()
            rmse_matrix_train[hh,] <- read.csv(paste0("rmse_train_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")%>% t() %>% as.numeric()
            rmse_matrix[hh,] <- read.csv(paste0("rmse_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";") %>% t() %>% as.numeric()
            coverage_matrix[hh,] <- read.csv(paste0("coverage_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";") %>% t() %>% as.numeric()
            coveragep_matrix[hh,] <- read.csv(paste0("coveragep_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";") %>% t() %>% as.numeric()
            range_matrix[hh,] <- read.csv(paste0("range_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")%>% t() %>% as.numeric()
            rangep_matrix[hh,] <- read.csv(paste0("rangep_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"__hh_",hh,"_partial.csv"), sep = ";")%>% t() %>% as.numeric()
          }
          
          #Now I can save in the same format
          rownames(times_matrix) <- NULL
          rownames(rmse_matrix_train) <- NULL
          rownames(rmse_matrix) <- NULL
          rownames(coverage_matrix) <- NULL
          rownames(coveragep_matrix) <- NULL
          rownames(range_matrix) <- NULL
          rownames(rangep_matrix) <- NULL
          write.table(times_matrix, paste0("times_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"_total.csv"), sep = ";")
          write.table(rmse_matrix_train, paste0("rmse_train_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"_total.csv"), sep = ";")
          write.table(rmse_matrix, paste0("rmse_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"_total.csv"), sep = ";")
          write.table(coverage_matrix, paste0("coverage_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"_total.csv"), sep = ";")
          write.table(coveragep_matrix, paste0("coveragep_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"_total.csv"), sep = ";")
          write.table(range_matrix, paste0("range_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"_total.csv"), sep = ";")
          write.table(rangep_matrix, paste0("rangep_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",i,"_total.csv"), sep = ";")
          
        }
      }
    }
    
  }
  
}













