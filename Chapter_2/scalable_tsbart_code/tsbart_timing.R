setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Basis/multibart/timing")

#This is one of Jennifer's examples adapted for timing
library(parallel)
library(doParallel)
rep <- 10
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
foreach(i = 1:rep, .packages=c('tsbart','multibart','MASS')) %dopar%  {
  for(j in c(8,16,32,64,128,256)){
    set.seed(123*i)
    basis_size <- j
    n <- 100
    t <- seq(1,8, length.out = basis_size)

    x12 <- MASS::mvrnorm(n, mu=c(0,0), Sigma=matrix(c(1,0.2,0.2,1), nrow=2))
    x34 <- MASS::mvrnorm(n, mu=c(0,0), Sigma=matrix(c(1,0.2,0.2,1), nrow=2))

    xa <- NULL
    xb <- NULL
    tt <- NULL
    for(k in 1:basis_size){
      xa <- rbind(xa, x12)
      xb <- rbind(xb, x34)
      tt <- c(tt, rep(t[k], times = n))
    }

    y <- (xa[,1]+xa[,2])*cos(tt+2*pi*(xb[,1]+xb[,2])) + rnorm(n*basis_size)

    #Seeing the truee value for patient 5
    #plot(t,(xa[5,1]+xa[5,2])*cos(t+2*pi*(xb[5,1]+xb[5,2])))
    #curve((xa[5,1]+xa[5,2])*cos(x+2*pi*(xb[5,1]+xb[5,2])), add = T)


    x_train <- cbind(xa,xb)
    y_train <- y
    t_train <- tt
    x_pred <- x_train
    t_pred <- t_train

    time1 <- Sys.time()
    fit_1 <- tsbart:::tsbart(y_train,
                             tgt = t_train,
                             x = x_train,
                             tpred = t_pred,
                             xpred = x_pred,
                             nburn = 1000,
                             nsim = 1000,
                             ntree = 200,
                             ecross = 1,
                             base_tree = 0.95,
                             power_tree = 200,
                             verbose = T,
                             use_fscale = TRUE)
    ftime1 <- difftime(Sys.time(), time1, units = "secs")[[1]]
    rmse1 <- sqrt(mean((apply(fit_1$mcmcdraws,2,mean)-y)^2))

    rm(fit_1)
    gc()

    time2 <- Sys.time()
    fit_2 <- multibart:::bcf_core(y_train,
                                  z_train = rep(1, times = length(y_train)),
                                  z_out = rep(1, times = length(y_train)),
                                  tvar_con = t_train,
                                  x_control = x_train,
                                  x_control_out = x_pred,
                                  tvar_con_out = t_pred,
                                  vanilla = F, dart = F,
                                  enc_con = 1,
                                  ntree_moderate = 0,
                                  ntree_control = 200,
                                  nburn = 1000, nsim = 1000, nthin = 1)
    ftime2 <- difftime(Sys.time(), time2, units = "secs")[[1]]
    rmse2 <- sqrt(mean((apply(fit_2$mu_post,2,mean)-y)^2))

    rm(fit_2)
    gc()

    times_vec <- c(ftime1, ftime2, rmse1, rmse2)
    write.table(times_vec, paste0("tsbart_timing_i_",i,"_j_",j,".csv"), sep = ";")

  }
}
parallel::stopCluster(cl)


ts1 <- NULL
ts2 <- NULL
for(i in i:10){
  for(j in c(8,16,32,64,128,256)){
    ts1 <- c(ts1,read.table(paste0("tsbart_timing_i_",i,"_j_",j,".csv"), sep = ";")[1,])
    ts2 <- c(ts2,read.table(paste0("tsbart_timing_i_",i,"_j_",j,".csv"), sep = ";")[2,])
  }
}

tsm1 <- matrix(ts1, ncol = 6, byrow = T)
tsm2 <- matrix(ts2, ncol = 6, byrow = T)

time_v <- NULL
log2_v <- NULL
Model_v <- NULL
for(i in 1:10){
  for(j in 1:6){
    time_v <- c(time_v,tsm1[i,j])
    log2_v <- c(log2_v,j+2)
    Model_v <- c(Model_v,"tsBART")
  }
}
time_v2 <- NULL
log2_v2 <- NULL
Model_v2 <- NULL
for(i in 1:10){
  for(j in 1:6){
    time_v2 <- c(time_v2,tsm2[i,j])
    log2_v2 <- c(log2_v2,j+2)
    Model_v2 <- c(Model_v2,"Scalable tsBART")
  }
}
time_df <- data.frame(time=c(time_v,time_v2), log2=c(log2_v,log2_v2), Model=c(Model_v,Model_v2))


library(ggplot2)
ranges <- function(x){
  return(data.frame(ymin = min(x), y = mean(x), ymax = max(x)))
}
pdf("timing_ts_ggplot.pdf")
ggplot(time_df, aes(log2, time, colour = Model)) +
  stat_summary(fun.data=ranges, geom="pointrange") +
  xlab(expression(log[2]*m^'*')) +
  ylab("Time (seconds)") +
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(color = "gray80", size = 0.5)
  )
dev.off()



