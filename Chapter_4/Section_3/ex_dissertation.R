set.seed(123)

n <- 500

#Covariates
X1 <- runif(n, 0, 1)
X2 <- rnorm(n, 0, 1)

#Treatment Z in [-1, 1]
Z <- runif(n, -1, 1)

#True treatment effect curve (can be any smooth function)
# true_tau <- function(z) z + sin(pi * z)
true_tau <- function(z, x1) {
  (1 + x1) * z + sin(pi * z * (1 + 2 * x1))
}

#Define true response function
f <- function(x1, x2, z) {
  2 * x2 + true_tau(z, x1)
}

#Generate outcome
Y<- f(X1, X2, Z) + rnorm(n, 0, 1)

# Create data frame
dat <- data.frame(Y, Z, X1, X2)

#The three values I will try to estimate (using monte carlo)
mean(true_tau(0.5, seq(0,1, length.out=10000))-true_tau(0, seq(0,1, length.out=10000)))
mean(true_tau(0, seq(0,1, length.out=10000))-true_tau(-0.5, seq(0,1, length.out=10000)))
mean(true_tau(0.5, seq(0,1, length.out=10000))-true_tau(-0.5, seq(0,1, length.out=10000)))

#Run a default bcf+ model
source("C:/Users/pedro/OneDrive/Área de Trabalho/Important_so_doing_it_again/Nature/Project_Stress/Dataset_1/Code/gp_approx_fun.R")
set.seed(68763)
burn <- 50000
post <- 1000
th <- 5
fit <- conbart::bcf_core(as.matrix(Y),
                         z_train = rep(1, times = n), #placeholder
                         z_out = rep(1, times = n), #placeholder
                         tvar_mod = Z, #actual treatment
                         x_control = as.matrix(cbind(X2)),
                         x_moderate = as.matrix(cbind(X1)),
                         x_control_out = as.matrix(cbind(rep(0,times=n))),
                         x_moderate_out = as.matrix(cbind(X1)),
                         tvar_mod_out = rep(0, times = n),
                         vanilla = T, #Keeping the first part of the function a regular BART
                         dart = F, enc_mod = 3,
                         ntree_moderate = 50,
                         ntree_control = 200, #Same as bart
                         nburn = burn, nsim = post, nthin = th)

#Get time grid so I can calculate the curves using the betas
time_grid = c(-0.5,0,0.5)
{
  for(j in 1:dim(fit$coefs_mod_est)[3]){
    omega_mod_grid = (Omega((time_grid), dim(fit$coefs_mod_est)[1], fit$L_mod_post[j]) %*% sqrt(SDiag(1, fit$l_mod_post[j], dim(fit$coefs_mod_est)[1], fit$L_mod_post[j])))
    if(j==1){
      ite = lapply(1:dim(fit$coefs_mod_est)[2],
                   function(i) t(t(omega_mod_grid%*%fit$coefs_mod_est[,i,j])))
    }else{
      ite_aux = lapply(1:dim(fit$coefs_mod_est)[2],
                       function(i) t(t(omega_mod_grid%*%fit$coefs_mod_est[,i,j])))
      
      for(i in 1:dim(fit$coefs_mod_est)[2]){
        ite[[i]] <- cbind(ite[[i]],ite_aux[[i]])
      }
    }
  }
}


#Calculate the predictions
auxmat <- matrix(0,ncol=post,nrow=n)
auxmat2 <- matrix(0,ncol=post,nrow=n)
auxmat3 <- matrix(0,ncol=post,nrow=n)

for(j in 1:n){
  auxmat[j,] <- ite[[j]][2,]-ite[[j]][1,]
  auxmat2[j,] <- ite[[j]][3,]-ite[[j]][2,]
  auxmat3[j,] <- ite[[j]][3,]-ite[[j]][1,]
}



library(ggplot2)
library(tibble)
library(dplyr)
library(gridExtra)

#Compare some graphs
df1 <- tibble(estimate = apply(auxmat, 2, mean))
true1 <- mean(true_tau(0.5, seq(0,1, length.out=10000)) - true_tau(0, seq(0,1, length.out=10000)))

df2 <- tibble(estimate = apply(auxmat2, 2, mean))
true2 <- mean(true_tau(0, seq(0,1, length.out=10000)) - true_tau(-0.5, seq(0,1, length.out=10000)))

df3 <- tibble(estimate = apply(auxmat3, 2, mean))
true3 <- mean(true_tau(0.5, seq(0,1, length.out=10000)) - true_tau(-0.5, seq(0,1, length.out=10000)))

plot_hist <- function(df, true_value, title) {
  ggplot(df, aes(x = estimate)) +
    geom_density(fill = "lightblue", color = "black") +
    geom_vline(xintercept = true_value, color = "red", linetype = "dashed", size = 1) +
    labs(title = title, x = "Estimated ATE", y = "Count") +
    theme_minimal(base_size = 14)
}

p1 <- plot_hist(df1, true1, "ATE: t=0.5 vs t=0")
p2 <- plot_hist(df2, true2, "ATE: t=0 vs t=-0.5")
p3 <- plot_hist(df3, true3, "ATE: t=0.5 vs t=-0.5")

# Arrange them in a single row
grid.arrange(p1, p2, p3, nrow = 1)


#Now for BART
set.seed(12453)
fit2 <- conbart::bcf_core(as.matrix(Y),
                         z_train = rep(1, times = n), #placeholder
                         z_out = rep(1, times = 3*n), #placeholder
                         tvar_mod = Z, #actual treatment
                         x_control = as.matrix(cbind(X2,X1,Z)),
                         x_moderate = as.matrix(cbind(X1)),
                         x_control_out = as.matrix(cbind(rep(X2, each = 3),rep(X1, each = 3),rep(c(-0.5,0,0.5), times = n))),
                         x_moderate_out = as.matrix(cbind(rep(X1, each = 3))),
                         tvar_mod_out = Z,
                         vanilla = T, #Keeping the first part of the function a regular BART
                         dart = F, enc_mod = 3,
                         ntree_moderate = 0,
                         ntree_control = 200, #Same as bart
                         nburn = burn, nsim = post, nthin = th)

#version with another prediction (predict function need tweaking)
#Seed is the same to guarantee the same model
set.seed(12453)
fit2b <- conbart::bcf_core(as.matrix(Y),
                          z_train = rep(1, times = n), #placeholder
                          z_out = rep(1, times = 3*n), #placeholder
                          tvar_mod = Z, #actual treatment
                          x_control = as.matrix(cbind(X2,X1,Z)),
                          x_moderate = as.matrix(cbind(X1)),
                          x_control_out = as.matrix(cbind(rep(rep(0, times = n), each = 3),rep(X1, each = 3),rep(c(-0.5,0,0.5), times = n))),
                          x_moderate_out = as.matrix(cbind(rep(X1, each = 3))),
                          tvar_mod_out = Z,
                          vanilla = T, #Keeping the first part of the function a regular BART
                          dart = F, enc_mod = 3,
                          ntree_moderate = 0,
                          ntree_control = 200, #Same as bart
                          nburn = burn, nsim = post, nthin = th)


#Calculate the predictions
bauxmat <- matrix(0,ncol=post,nrow=n)
bauxmat2 <- matrix(0,ncol=post,nrow=n)
bauxmat3 <- matrix(0,ncol=post,nrow=n)

bauxmat <- fit2b$mu_est_post[,seq(2, ncol(fit2b$mu_est_post), by = 3)]-fit2b$mu_est_post[,seq(1, ncol(fit2b$mu_est_post), by = 3)]
bauxmat2 <- fit2b$mu_est_post[,seq(3, ncol(fit2b$mu_est_post), by = 3)]-fit2b$mu_est_post[,seq(2, ncol(fit2b$mu_est_post), by = 3)]
bauxmat3 <- fit2b$mu_est_post[,seq(3, ncol(fit2b$mu_est_post), by = 3)]-fit2b$mu_est_post[,seq(1, ncol(fit2b$mu_est_post), by = 3)]

#And create the graphs
df1b <- tibble(estimate = apply(bauxmat, 1, mean))
true1 <- mean(true_tau(0.5, seq(0,1, length.out=10000)) - true_tau(0, seq(0,1, length.out=10000)))

df2b <- tibble(estimate = apply(bauxmat2, 1, mean))
true2 <- mean(true_tau(0, seq(0,1, length.out=10000)) - true_tau(-0.5, seq(0,1, length.out=10000)))

df3b <- tibble(estimate = apply(bauxmat3, 1, mean))
true3 <- mean(true_tau(0.5, seq(0,1, length.out=10000)) - true_tau(-0.5, seq(0,1, length.out=10000)))

p1b <- plot_hist(df1b, true1, "ATE: t=0.5 vs t=0")
p2b <- plot_hist(df2b, true2, "ATE: t=0 vs t=-0.5")
p3b <- plot_hist(df3b, true3, "ATE: t=0.5 vs t=-0.5")

# Arrange them in a single row
grid.arrange(p1b, p2b, p3b, nrow = 1)


#Again, bcf+, but with a different kind of prediction
#Setting X2 to be zero, so I can analyze X1
set.seed(68763)
burn <- 50000
post <- 1000
th <- 5
fit3 <- conbart::bcf_core(as.matrix(Y),
                         z_train = rep(1, times = n), #placeholder
                         z_out = rep(1, times = n), #placeholder
                         tvar_mod = Z, #actual treatment
                         x_control = as.matrix(cbind(X2)),
                         x_moderate = as.matrix(cbind(X1)),
                         x_control_out = as.matrix(cbind(rep(0, times = n))),
                         x_moderate_out = as.matrix(cbind(X1)),
                         tvar_mod_out = rep(0, times = n),
                         vanilla = T, #Keeping the first part of the function a regular BART
                         dart = F, enc_mod = 3,
                         ntree_moderate = 100,
                         ntree_control = 200, #Same as bart
                         nburn = burn, nsim = post, nthin = th)

time_grid = c(-0.5,0,0.5)
{
  for(j in 1:dim(fit3$coefs_mod_est)[3]){
    omega_mod_grid = (Omega((time_grid), dim(fit3$coefs_mod_est)[1], fit3$L_mod_post[j]) %*% sqrt(SDiag(1, fit3$l_mod_post[j], dim(fit3$coefs_mod_est)[1], fit3$L_mod_post[j])))
    if(j==1){
      ite = lapply(1:dim(fit3$coefs_mod_est)[2],
                   function(i) t(t(omega_mod_grid%*%fit3$coefs_mod_est[,i,j])))
    }else{
      ite_aux = lapply(1:dim(fit3$coefs_mod_est)[2],
                       function(i) t(t(omega_mod_grid%*%fit3$coefs_mod_est[,i,j])))
      
      for(i in 1:dim(fit3$coefs_mod_est)[2]){
        ite[[i]] <- cbind(ite[[i]],ite_aux[[i]])
      }
    }
  }
  
}



cauxmat <- matrix(0,ncol=post,nrow=n)
cauxmat2 <- matrix(0,ncol=post,nrow=n)
cauxmat3 <- matrix(0,ncol=post,nrow=n)

for(j in 1:n){
  cauxmat[j,] <- ite[[j]][2,]-ite[[j]][1,]
  cauxmat2[j,] <- ite[[j]][3,]-ite[[j]][2,]
  cauxmat3[j,] <- ite[[j]][3,]-ite[[j]][1,]
}


df1c <- tibble(estimate = c(apply(bauxmat, 1, mean),
                            apply(cauxmat, 2, mean)),
               Method = rep(c("BART", "bcf-con"), each = ncol(auxmat)))
true1 <- mean(true_tau(0.5, seq(0,1, length.out=10000)) - true_tau(0, seq(0,1, length.out=10000)))

df2c <- tibble(estimate = c(apply(bauxmat2, 1, mean),
                            apply(cauxmat2, 2, mean)),
               Method = rep(c("BART", "bcf-con"), each = ncol(auxmat)))
true2 <- mean(true_tau(0, seq(0,1, length.out=10000)) - true_tau(-0.5, seq(0,1, length.out=10000)))

df3c <- tibble(estimate = c(apply(bauxmat3, 1, mean),
                            apply(cauxmat3, 2, mean)),
               Method = rep(c("BART", "bcf-con"), each = ncol(auxmat)))
true3 <- mean(true_tau(0.5, seq(0,1, length.out=10000)) - true_tau(-0.5, seq(0,1, length.out=10000)))


plot_hist_compare <- function(df, true_value, title) {
  ggplot(df, aes(x = estimate, fill = Method)) +
    geom_density(alpha = 0.6, position = "identity", color = "black") +
    geom_vline(xintercept = true_value, color = "red", linetype = "dashed", size = 1) +
    labs(title = title, x = "Estimated ATE", y = "Count") +
    scale_fill_manual(values = c("skyblue", "orange")) +
    theme_minimal(base_size = 14)
}


p1c <- plot_hist_compare(df1c, true1, "ATE: t=0.5 vs t=0.0")
p2c <- plot_hist_compare(df2c, true2, "ATE: t=0.0 vs t=-0.5")
p3c <- plot_hist_compare(df3c, true3, "ATE: t=0.5 vs t=-0.5")

grid.arrange(p1c, p2c, p3c, nrow = 1)

library(cowplot)
#Removing the legends, because I don't need them
plot_density_compare_nolegend <- function(df, true_value, title, xleg) {
  ggplot(df, aes(x = estimate, color = Method, fill = Method)) +
    geom_density(alpha = 0.4, size = 1) +
    geom_vline(xintercept = true_value, color = "red", linetype = "dashed", size = 1) +
    labs(title = title, x = xleg, y = "Density") +
    ylim(0,2.2)+
    scale_fill_manual(values = c("skyblue", "orange")) +
    scale_color_manual(values = c("skyblue", "orange")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
}

# Create plots without legends
p1 <- plot_density_compare_nolegend(df1c, true1, "", expression(hat(ATE)[(0 * "," * 0.5)])) + scale_x_continuous(breaks = c(-1, 0, 1))
p2 <- plot_density_compare_nolegend(df2c, true2, "", expression(hat(ATE)[(-0.5 * "," * 0)])) + scale_x_continuous(breaks = c(-1, 0, 1, 2))
p3 <- plot_density_compare_nolegend(df3c, true3, "", expression(hat(ATE)[(-0.5 * "," * 0.5)])) + scale_x_continuous(breaks = c(-1, 0, 1, 2))

# Combine with one legend
pdf("con_simulation.pdf", height = 4, width = 8)
plot_grid(
  plot_grid(p1, p2, p3, nrow = 1), 
  ncol = 1,
  rel_heights = c(1, 0.1)
)
dev.off()

#Curve on X
idx1 <- order(X1)

pdf("cate.pdf", width = 12, height = 9)
par(mfrow=c(2,3))
plot(X1[idx1], true_tau(0.5, X1[idx1]) - true_tau(0, X1[idx1]), col = "purple", 
     ylim = c(-2,3), type = "l", lwd = 2, xlab = "X1", ylab = expression(CATE[(0 * "," * 0.5)](x)))
lines(X1[idx1],apply(bauxmat,2,mean)[idx1], col = "blue", lwd = 2)
lines(X1[idx1],apply(bauxmat,2,quantile, 0.025)[idx1], col = "blue", lty = 2, lwd = 2)
lines(X1[idx1],apply(bauxmat,2,quantile, 0.975)[idx1], col = "blue", lty = 2, lwd = 2)
legend("topright", legend = c("True function", "BART"), col = c("purple", "blue"), lty = c(1,1), lwd = 2)

plot(X1[idx1], true_tau(0, X1[idx1]) - true_tau(-0.5, X1[idx1]), col = "purple", 
     ylim = c(-2,3), type = "l", lwd = 2, xlab = "X1", ylab = expression(CATE[(0 * "," * 0.5)](x)))
lines(X1[idx1],apply(bauxmat2,2,mean)[idx1], col = "blue", lwd = 2)
lines(X1[idx1],apply(bauxmat2,2,quantile, 0.025)[idx1], col = "blue", lty = 2, lwd = 2)
lines(X1[idx1],apply(bauxmat2,2,quantile, 0.975)[idx1], col = "blue", lty = 2, lwd = 2)
legend("topright", legend = c("True function", "BART"), col = c("purple", "blue"), lty = c(1,1), lwd = 2)


plot(X1[idx1], true_tau(0.5, X1[idx1]) - true_tau(-0.5, X1[idx1]), col = "purple", 
     ylim = c(-2,4), type = "l", lwd = 2, xlab = "X1", ylab = expression(CATE[(0 * "," * 0.5)](x)))
lines(X1[idx1],apply(bauxmat3,2,mean)[idx1], col = "blue", lwd = 2)
lines(X1[idx1],apply(bauxmat3,2,quantile, 0.025)[idx1], col = "blue", lty = 2, lwd = 2)
lines(X1[idx1],apply(bauxmat3,2,quantile, 0.975)[idx1], col = "blue", lty = 2, lwd = 2)
legend("topright", legend = c("True function", "BART"), col = c("purple", "blue"), lty = c(1,1), lwd = 2)


plot(X1[idx1], true_tau(0.5, X1[idx1]) - true_tau(0, X1[idx1]), col = "purple", 
     ylim = c(-2,3), type = "l", lwd = 2, xlab = "X1", ylab = expression(CATE[(0 * "," * 0.5)](x)))
lines(X1[idx1],apply(cauxmat,1,mean)[idx1], col = "red", lwd = 2)
lines(X1[idx1],apply(cauxmat,1,quantile, 0.025)[idx1], col = "red", lty = 2, lwd = 2)
lines(X1[idx1],apply(cauxmat,1,quantile, 0.975)[idx1], col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("True function", "bcf+"), col = c("purple", "red"), lty = c(1,1), lwd = 2)

plot(X1[idx1], true_tau(0, X1[idx1]) - true_tau(-0.5, X1[idx1]), col = "purple", 
     ylim = c(-2,3), type = "l", lwd = 2, xlab = "X1", ylab = expression(CATE[(0 * "," * 0.5)](x)))
lines(X1[idx1],apply(cauxmat2,1,mean)[idx1], col = "red", lwd = 2)
lines(X1[idx1],apply(cauxmat2,1,quantile, 0.025)[idx1], col = "red", lty = 2, lwd = 2)
lines(X1[idx1],apply(cauxmat2,1,quantile, 0.975)[idx1], col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("True function", "bcf+"), col = c("purple", "red"), lty = c(1,1), lwd = 2)

plot(X1[idx1], true_tau(0.5, X1[idx1]) - true_tau(-0.5, X1[idx1]), col = "purple", 
     ylim = c(-2,4), type = "l", lwd = 2, xlab = "X1", ylab = expression(CATE[(0 * "," * 0.5)](x)))
lines(X1[idx1],apply(cauxmat3,1,mean)[idx1], col = "red", lwd = 2)
lines(X1[idx1],apply(cauxmat3,1,quantile, 0.025)[idx1], col = "red", lty = 2, lwd = 2)
lines(X1[idx1],apply(cauxmat3,1,quantile, 0.975)[idx1], col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("True function", "bcf+"), col = c("purple", "red"), lty = c(1,1), lwd = 2)
dev.off()




#Conditional dose-response curve
X1pred

#I need a predict function working so I don't need to keep running the models again
X1pred <- c(0.2, 0.5, 0.8)
set.seed(68763)
burn <- 50000
post <- 1000
th <- 5
fit5 <- conbart::bcf_core(as.matrix(Y),
                          z_train = rep(1, times = n), #placeholder
                          z_out = rep(1, times = n), #placeholder
                          tvar_mod = Z, #actual treatment
                          x_control = as.matrix(cbind(X2)),
                          x_moderate = as.matrix(cbind(X1)),
                          x_control_out = as.matrix(cbind(rep(0, times = 3))),
                          x_moderate_out = as.matrix(cbind(X1pred)),
                          tvar_mod_out = rep(0, times = n),
                          vanilla = T, #Keeping the first part of the function a regular BART
                          dart = F, enc_mod = 3,
                          ntree_moderate = 100,
                          ntree_control = 200, #Same as bart
                          nburn = burn, nsim = post, nthin = th)


time_grid = seq(-1,1, length.out = 100)
{
  for(j in 1:dim(fit5$coefs_mod_est)[3]){
    omega_mod_grid = (Omega((time_grid), dim(fit5$coefs_mod_est)[1], fit5$L_mod_post[j]) %*% sqrt(SDiag(1, fit5$l_mod_post[j], dim(fit5$coefs_mod_est)[1], fit5$L_mod_post[j])))
    if(j==1){
      ite = lapply(1:dim(fit5$coefs_mod_est)[2],
                   function(i) t(t(omega_mod_grid%*%fit5$coefs_mod_est[,i,j])))
    }else{
      ite_aux = lapply(1:dim(fit5$coefs_mod_est)[2],
                       function(i) t(t(omega_mod_grid%*%fit5$coefs_mod_est[,i,j])))

      for(i in 1:dim(fit5$coefs_mod_est)[2]){
        ite[[i]] <- cbind(ite[[i]],ite_aux[[i]])
      }
    }
  }

}

#Now BART
#Seed is the same to guarantee the same model
set.seed(12453)
fit5b <- conbart::bcf_core(as.matrix(Y),
                           z_train = rep(1, times = n), #placeholder
                           z_out = rep(1, times = 100*3), #placeholder
                           tvar_mod = Z, #actual treatment
                           x_control = as.matrix(cbind(X2,X1,Z)),
                           x_moderate = as.matrix(cbind(X1)),
                           x_control_out = as.matrix(cbind(rep(rep(0, times = 3), each = 100),rep(X1pred, each = 100),rep(seq(-1,1, length.out = 100), times = 3))),
                           x_moderate_out = as.matrix(cbind(rep(X1pred, each = 100))),
                           tvar_mod_out = Z,
                           vanilla = T, #Keeping the first part of the function a regular BART
                           dart = F, enc_mod = 3,
                           ntree_moderate = 0,
                           ntree_control = 200, #Same as bart
                           nburn = burn, nsim = post, nthin = th)


#Now the graphs
pdf("doseres.pdf", width = 12, height = 9)
par(mfrow=c(2,3))
plot(time_grid, apply(fit5b$mu_est_post[,1:100],2, mean), type = "l", lwd = 2, col ="blue", ylim = c(-4,4), xlab = "z", ylab = expression(xi[z](x[1]* "=" *0.2 * "," * x[2]* "=" *0)))
lines(time_grid, true_tau(time_grid, X1pred[1]), col = "purple", lwd = 2)
lines(time_grid, apply(fit5b$mu_est_post[,1:100],2, quantile, 0.025), col = "blue", lwd = 2, lty = 2)
lines(time_grid, apply(fit5b$mu_est_post[,1:100],2, quantile, 0.975), col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("True function", "BART"), col = c("purple", "blue"), lty = c(1,1), lwd = 2)

plot(time_grid, apply(fit5b$mu_est_post[,101:200],2, mean), type = "l", lwd = 2, col ="blue", ylim = c(-4,4), xlab = "z", ylab = expression(xi[z](x[1]* "=" *0.5 * "," * x[2]* "=" *0)))
lines(time_grid, true_tau(time_grid, X1pred[2]), col = "purple", lwd = 2)
lines(time_grid, apply(fit5b$mu_est_post[,101:200],2, quantile, 0.025), col = "blue", lwd = 2, lty = 2)
lines(time_grid, apply(fit5b$mu_est_post[,101:200],2, quantile, 0.975), col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("True function", "BART"), col = c("purple", "blue"), lty = c(1,1), lwd = 2)

plot(time_grid, apply(fit5b$mu_est_post[,201:300],2, mean), type = "l", lwd = 2, col ="blue", ylim = c(-4,4), xlab = "z", ylab = expression(xi[z](x[1]* "=" *0.8 * "," * x[2]* "=" *0)))
lines(time_grid, true_tau(time_grid, X1pred[3]), col = "purple", lwd = 2)
lines(time_grid, apply(fit5b$mu_est_post[,201:300],2, quantile, 0.025), col = "blue", lwd = 2, lty = 2)
lines(time_grid, apply(fit5b$mu_est_post[,201:300],2, quantile, 0.975), col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("True function", "BART"), col = c("purple", "blue"), lty = c(1,1), lwd = 2)

plot(time_grid, apply(ite[[1]], 1, mean), type = "l", lwd = 2, col ="red", ylim = c(-4,4), xlab = "z", ylab = expression(xi[z](x[1]* "=" *0.2 * "," * x[2]* "=" *0)))
lines(time_grid, true_tau(time_grid, X1pred[1]), col = "purple", lwd = 2)
lines(time_grid, apply(ite[[1]], 1, quantile, 0.025), col = "red", lwd = 2, lty = 2)
lines(time_grid, apply(ite[[1]], 1, quantile, 0.975), col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("True function", "bcf+"), col = c("purple", "red"), lty = c(1,1), lwd = 2)

plot(time_grid, apply(ite[[2]], 1, mean), type = "l", lwd = 2, col ="red", ylim = c(-4,4), xlab = "z", ylab = expression(xi[z](x[1]* "=" *0.5 * "," * x[2]* "=" *0)))
lines(time_grid, true_tau(time_grid, X1pred[2]), col = "purple", lwd = 2)
lines(time_grid, apply(ite[[2]], 1, quantile, 0.025), col = "red", lwd = 2, lty = 2)
lines(time_grid, apply(ite[[2]], 1, quantile, 0.975), col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("True function", "bcf+"), col = c("purple", "red"), lty = c(1,1), lwd = 2)

plot(time_grid, apply(ite[[3]], 1, mean), type = "l", lwd = 2, col ="red", ylim = c(-4,4), xlab = "z", ylab = expression(xi[z](x[1]* "=" *0.8 * "," * x[2]* "=" *0)))
lines(time_grid, true_tau(time_grid, X1pred[3]), col = "purple", lwd = 2)
lines(time_grid, apply(ite[[3]], 1, quantile, 0.025), col = "red", lwd = 2, lty = 2)
lines(time_grid, apply(ite[[3]], 1, quantile, 0.975), col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("True function", "bcf+"), col = c("purple", "red"), lty = c(1,1), lwd = 2)
dev.off()


