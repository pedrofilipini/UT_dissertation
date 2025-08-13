{
library(dplyr)
library(ggplot2)
library(latex2exp)
library(patchwork)
  
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
  
  ret = x %>% dplyr::select(fixedmindset_baseline,
                            stressmindset_baseline, bothmindsets,
                            selfesteem_baseline, sex) %>%
    dbarts::makeModelMatrixFromDataFrame(drop = FALSE)
}

get_moderators = function(x) {
  x %>% dplyr::select(bothmindsets,
                      stressmindset_baseline,
                      fixedmindset_baseline,
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


#The math model
npred <- length(y_math)
burn <- 15000
th <- 3
post <- 5000

}



#The best models
load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart/cv_tpr_small/df_cvmath.RData")
load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart/cv_tpr_small/df_cvspeech.RData")

#Best for bcf
df1 <- df_math %>% filter(linear=="none")
df2 <- df_speech %>% filter(linear=="none")

#Best math
df_math[which.min(df_math$rmse_test),]

rep <- 10
model_lbcf_math <- list()
model_bcf_math <- list()
model_lbcf_speech <- list()
model_bcf_speech <- list()

# for(i in 1:rep){
# set.seed(2135+i)
# 
# lbcf_math <- lbart:::bcf_new(y_math, #response
#                              z = z_math, #Z
#                              z_est = rep(1,npred), #Z pred
#                              x_control = x_con_math, #matrix of covariates
#                              x_control_est = x_con_math,
#                              x_moderate = x_mod_math,
#                              x_moderate_est = x_mod_math,
#                              pihat = rep(1,npred), #in case of causal inference
#                              pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
#                              include_pi = "none", #can be "none", "control", "moderate", or "both"
#                              linear = "both", #can be "none", "control", "moderate", or "both"
#                              base_control = 0.5, base_moderate = 0.25, #hyperparameters for control
#                              power_control = 3, power_moderate = 3, #hyperparameters for moderate
#                              nburn = burn, nthin = th, nsim = post, #draws
#                              ntree_control = 200, #control trees
#                              ntree_moderate = 100, #moderate trees (treatment effect)
#                              dart = F, #Linero's Prior
#                              save_trees = T, 
#                              use_muscale = T, 
#                              use_tauscale = T) #Not saving the trees
# model_lbcf_math[[i]] <- lbcf_math
# 
# #Best bcf
# df1[which.min(df1$rmse_test),]
# 
# set.seed(2135+i)
# bcf_math <- lbart:::bcf_new(y_math, #response
#                              z = z_math, #Z
#                              z_est = rep(1,npred), #Z pred
#                              x_control = x_con_math, #matrix of covariates
#                              x_control_est = x_con_math,
#                              x_moderate = x_mod_math,
#                              x_moderate_est = x_mod_math,
#                              pihat = rep(1,npred), #in case of causal inference
#                              pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
#                              include_pi = "none", #can be "none", "control", "moderate", or "both"
#                              linear = "none", #can be "none", "control", "moderate", or "both"
#                              base_control = 0.5, base_moderate = 0.5, #hyperparameters for control
#                              power_control = 1, power_moderate = 1, #hyperparameters for moderate
#                              nburn = burn, nthin = th, nsim = post, #draws
#                              ntree_control = 50, #control trees
#                              ntree_moderate = 100, #moderate trees (treatment effect)
#                              dart = F, #Linero's Prior
#                              save_trees = T, 
#                              use_muscale = F, 
#                              use_tauscale = T) #Not saving the trees
# model_bcf_math[[i]] <- bcf_math
# 
# }




# b_og <- bcf::bcf(y=y_math,
#     z=z_math,
#     pihat = rep(1,npred),
#     x_control = as.matrix(x_con_math),n_chains = 1,
#     x_moderate = as.matrix(x_mod_math), use_muscale = F, use_tauscale = F, n_threads = 1,
#     nburn = 1000, nthin = 1, nsim = 5000)
# 


#The speech model
npred <- length(y_speech)
burn <- 15000
th <- 3
post <- 5000

df_speech[which.min(df_speech$rmse_test),]


# for(i in 1:rep){
# set.seed(2135+i)
# lbcf_speech <- lbart:::bcf_new(y_speech, #response
#                                z = z_speech, #Z
#                                z_est = rep(1,npred), #Z pred
#                                x_control = x_con_speech, #matrix of covariates
#                                x_control_est = x_con_speech,
#                                x_moderate = x_mod_speech,
#                                x_moderate_est = x_mod_speech,
#                                pihat = rep(1,npred), #in case of causal inference
#                                pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
#                                include_pi = "none", #can be "none", "control", "moderate", or "both"
#                                linear = "both", #can be "none", "control", "moderate", or "both"
#                                base_control = 0.5, base_moderate = 0.5, #hyperparameters for control
#                                power_control = 3, power_moderate = 3, #hyperparameters for moderate
#                                nburn = burn, nthin = th, nsim = post, #draws
#                                ntree_control = 200, #control trees
#                                ntree_moderate = 100, #moderate trees (treatment effect)
#                                dart = F, #Linero's Prior
#                                save_trees = T, 
#                                use_muscale = T, 
#                                use_tauscale = T) #Not saving the trees
# 
# model_lbcf_speech[[i]] <- lbcf_speech
# #Best bcf
# df2[which.min(df2$rmse_test),]
# #I got one with more trees, the best seems to be too small of a model
# set.seed(2135+i)
# bcf_speech <- lbart:::bcf_new(y_speech, #response
#                               z = z_speech, #Z
#                               z_est = rep(1,npred), #Z pred
#                               x_control = x_con_speech, #matrix of covariates
#                               x_control_est = x_con_speech,
#                               x_moderate = x_mod_speech,
#                               x_moderate_est = x_mod_speech,
#                               pihat = rep(1,npred), #in case of causal inference
#                               pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
#                               include_pi = "none", #can be "none", "control", "moderate", or "both"
#                               linear = "none", #can be "none", "control", "moderate", or "both"
#                               base_control = 0.8, base_moderate = 0.25, #hyperparameters for control
#                               power_control = 2, power_moderate = 3, #hyperparameters for moderate
#                               nburn = burn, nthin = th, nsim = post, #draws
#                               ntree_control = 50, #control trees
#                               ntree_moderate = 100, #moderate trees (treatment effect)
#                               dart = F, #Linero's Prior
#                               save_trees = T, 
#                               use_muscale = F, 
#                               use_tauscale = T) #Not saving the trees
# model_bcf_speech[[i]] <- bcf_speech
# }


# save(model_lbcf_math, file = "lbcf_math.RData")
# save(model_bcf_math, file = "bcf_math.RData")
# save(model_lbcf_speech, file = "lbcf_speech.RData")
# save(model_bcf_speech, file = "bcf_speech.RData")


load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/lbcf_math.RData")
load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/bcf_math.RData")
load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/lbcf_speech.RData")
load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/bcf_speech.RData")



{bcf_mu_math <- NULL
  bcf_tau_math <- NULL
  bcf_eta_con_math <- NULL
  bcf_eta_mod_math <- NULL
  bcf_sigma_math <- NULL
  lbcf_mu_math <- NULL
  lbcf_tau_math <- NULL
  lbcf_eta_con_math <- NULL
  lbcf_eta_mod_math <- NULL
  lbcf_sigma_math <- NULL}

for(i in 1:rep){
  bcf_mu_math <- rbind(bcf_mu_math, model_bcf_math[[i]]$mu_est_post)
  bcf_tau_math <- rbind(bcf_tau_math, model_bcf_math[[i]]$b_est_post)
  bcf_eta_con_math <- c(bcf_eta_con_math, model_bcf_math[[i]]$eta_con)
  bcf_eta_mod_math <- c(bcf_eta_mod_math, model_bcf_math[[i]]$eta_mod)
  bcf_sigma_math <- c(bcf_sigma_math, model_bcf_math[[i]]$sigma)
  lbcf_mu_math <- rbind(lbcf_mu_math, model_lbcf_math[[i]]$mu_est_post)
  lbcf_tau_math <- rbind(lbcf_tau_math, model_lbcf_math[[i]]$b_est_post)
  lbcf_eta_con_math <- c(lbcf_eta_con_math, model_lbcf_math[[i]]$eta_con)
  lbcf_eta_mod_math <- c(lbcf_eta_mod_math, model_lbcf_math[[i]]$eta_mod)
  lbcf_sigma_math <- c(lbcf_sigma_math, model_lbcf_math[[i]]$sigma)
}

{bcf_mu_speech <- NULL
  bcf_tau_speech <- NULL
  bcf_eta_con_speech <- NULL
  bcf_eta_mod_speech <- NULL
  bcf_sigma_speech <- NULL
  lbcf_mu_speech <- NULL
  lbcf_tau_speech <- NULL
  lbcf_eta_con_speech <- NULL
  lbcf_eta_mod_speech <- NULL
  lbcf_sigma_speech <- NULL}

for(i in 1:rep){
  bcf_mu_speech <- rbind(bcf_mu_speech, model_bcf_speech[[i]]$mu_est_post)
  bcf_tau_speech <- rbind(bcf_tau_speech, model_bcf_speech[[i]]$b_est_post)
  bcf_eta_con_speech <- c(bcf_eta_con_speech, model_bcf_speech[[i]]$eta_con)
  bcf_eta_mod_speech <- c(bcf_eta_mod_speech, model_bcf_speech[[i]]$eta_mod)
  bcf_sigma_speech <- c(bcf_sigma_speech, model_bcf_speech[[i]]$sigma)
  lbcf_mu_speech <- rbind(lbcf_mu_speech, model_lbcf_speech[[i]]$mu_est_post)
  lbcf_tau_speech <- rbind(lbcf_tau_speech, model_lbcf_speech[[i]]$b_est_post)
  lbcf_eta_con_speech <- c(lbcf_eta_con_speech, model_lbcf_speech[[i]]$eta_con)
  lbcf_eta_mod_speech <- c(lbcf_eta_mod_speech, model_lbcf_speech[[i]]$eta_mod)
  lbcf_sigma_speech <- c(lbcf_sigma_speech, model_lbcf_speech[[i]]$sigma)
}

#Traceplots to include in the supplemental material
plot(lbcf_sigma_speech, type = "l")
plot(lbcf_eta_mod_speech, type = "l")
plot(lbcf_eta_con_speech, type = "l")
set.seed(12354)
samp <- sample(dim(lbcf_tau_speech)[2], size = 10, replace = F)
for(i in samp){
  plot(lbcf_tau_speech[,i], type = "l")
  plot(lbcf_mu_speech[,i], type = "l")
}



#ITEs and ATEs
ite_math <- data.frame(mean = apply(bcf_tau_math,2, mean), 
                       lb = apply(bcf_tau_math,2, quantile, 0.025), 
                       ub = apply(bcf_tau_math,2, quantile, 0.975))
cate_math <- data.frame(mean = apply(bcf_tau_math,1, mean), 
                        lb = apply(bcf_tau_math,1, quantile, 0.025), 
                        ub = apply(bcf_tau_math,1, quantile, 0.975))

lite_math <- data.frame(mean = apply(lbcf_tau_math,2, mean), 
                       lb = apply(lbcf_tau_math,2, quantile, 0.025), 
                       ub = apply(lbcf_tau_math,2, quantile, 0.975))
lcate_math <- data.frame(mean = apply(lbcf_tau_math,1, mean), 
                        lb = apply(lbcf_tau_math,1, quantile, 0.025), 
                        ub = apply(lbcf_tau_math,1, quantile, 0.975))


#The story: got those two models, taking the average of the epoch
#Did 5-fold cross validation to select the best model

#Results are align, in general, but linear seems to have more magnitude

#ITE and CATE for math
boxplot(lite_math$mean, ite_math$mean,
        main = "Posterior Mean ITE", names= c("lbcf","bcf"))

pdf("density_ite_math.pdf")
plot(density(lite_math$mean), ylim = c(0,0.025),
     main = "Density of ITE posterior mean - Math epoch")
lines(density(ite_math$mean), col = "red")
legend("topleft", legend = c("lbcf", "bcf"), fill = c("black", "red"))
dev.off()

hist(lite_math$mean, breaks = 70)
hist(ite_math$mean, breaks = 70)

boxplot(lcate_math$mean, cate_math$mean,
        main = "Posterior Mean ATE", names= c("lbcf","bcf"))

pdf("density_ate_math.pdf")
plot(density(lcate_math$mean), ylim = c(0,0.008),
     main = "Density of ATE posterior - Math epoch")
lines(density(cate_math$mean), col = "red")
legend("topleft", legend = c("lbcf", "bcf"), fill = c("black", "red"))
dev.off()

#Results, in general, aligned
plot(lite_math$mean, ite_math$mean, col = x_mod_math[,4]+1, pch = 19,
     xlab = "lbcf", ylab = "bcf", main = "Posterior Mean ITE")
abline(0,1)


#ITEs and ATEs
ite_speech <- data.frame(mean = apply(bcf_tau_speech,2, mean), 
                       lb = apply(bcf_tau_speech,2, quantile, 0.025), 
                       ub = apply(bcf_tau_speech,2, quantile, 0.975))
cate_speech <- data.frame(mean = apply(bcf_tau_speech,1, mean), 
                        lb = apply(bcf_tau_speech,1, quantile, 0.025), 
                        ub = apply(bcf_tau_speech,1, quantile, 0.975))

lite_speech <- data.frame(mean = apply(lbcf_tau_speech,2, mean), 
                        lb = apply(lbcf_tau_speech,2, quantile, 0.025), 
                        ub = apply(lbcf_tau_speech,2, quantile, 0.975))
lcate_speech <- data.frame(mean = apply(lbcf_tau_speech,1, mean), 
                         lb = apply(lbcf_tau_speech,1, quantile, 0.025), 
                         ub = apply(lbcf_tau_speech,1, quantile, 0.975))


#Same goes for speech
#Same thing, linear has more magnitude, mostly because
#it is separating sex more

#ITE and CATE for speech
boxplot(lite_speech$mean, ite_speech$mean,
        main = "Posterior Mean ITE", names= c("lbcf","bcf"))
boxplot(lcate_speech$mean, cate_speech$mean,
        main = "Posterior ATE", names= c("lbcf","bcf"))


pdf("density_ite_speech.pdf")
plot(density(lite_speech$mean), ylim = c(0,0.016),
     main = "Density of ITE posterior mean - Speech epoch")
lines(density(ite_speech$mean), col = "red")
legend("topleft", legend = c("lbcf", "bcf"), fill = c("black", "red"))
dev.off()

pdf("density_ate_speech.pdf")
plot(density(lcate_speech$mean), ylim = c(0,0.008),
     main = "Density of ATE posterior - Speech epoch")
lines(density(cate_speech$mean), col = "red")
legend("topleft", legend = c("lbcf", "bcf"), fill = c("black", "red"))
dev.off()

mean(lcate_speech$mean<0)
mean(cate_speech$mean<0)

#Again, results aligned
plot(lite_speech$mean, ite_speech$mean, col = x_mod_speech[,4]+1, pch = 19,
     xlab = "lbcf", ylab = "bcf")
abline(0,1)



#Looking for subgroups using Spencers package
library(possum)

#Getting the data
fhat <- rowMeans(t(lbcf_tau_math))
fhatbcf <- rowMeans(t(bcf_tau_math))
# fhatog <- rowMeans(t(b_og$tau))

#Additive summary for math
ad_sum <- additive_summary(
  fhat~s(stressmindset_baseline)+s(bothmindsets)+s(fixedmindset_baseline)+(sex.1),
  fhatSamples=t(lbcf_tau_math),
  fhat = rowMeans(t(lbcf_tau_math)),
  df = as.data.frame(x_mod_math),
  alpha = 0.05,
  #fast = TRUE,
  #quants = seq(0, 1, by = 0.005),
  #grid_size = 100,
  verbose = FALSE,
  return_samples = TRUE,
  meta = NA
)

resid_sum <- t(lbcf_tau_math)-ad_sum$fittedValues
pdf("sum_res_math_1_lbcf.pdf")
plot(x_mod_math[,1], apply(resid_sum, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Bothmindsets", main = "lbcf summary residuals - Math epoch")
dev.off()
pdf("sum_res_math_2_lbcf.pdf")
plot(x_mod_math[,2], apply(resid_sum, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Stress mindset", main = "lbcf summary residuals - Math epoch")
dev.off()
pdf("sum_res_math_3_lbcf.pdf")
plot(x_mod_math[,3], apply(resid_sum, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Fixed mindset", main = "lbcf summary residuals - Math epoch")
dev.off()

plot(x_mod_math[,4], apply(resid_sum, 1, mean))


# gm <- gam(fhat ~ s(stressmindset_baseline)+s(fixedmindset_baseline)+factor(sex.1), data = as.data.frame(x_mod_math))
# summary(gm)
# 
# gm2 <- gam(fhatbcf ~ s(stressmindset_baseline)+s(fixedmindset_baseline)+factor(sex.1), data = as.data.frame(x_mod_math))
# summary(gm2)
# 
# gm3 <- gam(fhatog ~ s(stressmindset_baseline)+s(fixedmindset_baseline)+factor(sex.1), data = as.data.frame(x_mod_math))
# summary(gm3)

#Compared with bcf additive summary for math
ad_sum2 <- additive_summary(
  rowMeans(t(bcf_tau_math))~s(bothmindsets)+s(stressmindset_baseline)+s(fixedmindset_baseline)+(sex.1),
  fhatSamples=t(bcf_tau_math),
  fhat = rowMeans(t(bcf_tau_math)),
  df = as.data.frame(x_mod_math),
  alpha = 0.05,
  #fast = TRUE,
  #quants = seq(0, 1, by = 0.005),
  grid_size = 100,
  verbose = FALSE,
  return_samples = TRUE,
  meta = NA
)


resid_sum2 <- t(bcf_tau_math)-ad_sum2$fittedValues

pdf("sum_res_math_1_bcf.pdf")
plot(x_mod_math[,1], apply(resid_sum2, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Bothmindsets", main = "bcf summary residuals - Math epoch")
dev.off()
pdf("sum_res_math_2_bcf.pdf")
plot(x_mod_math[,2], apply(resid_sum2, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Stress mindset", main = "bcf summary residuals - Math epoch")
dev.off()
pdf("sum_res_math_3_bcf.pdf")
plot(x_mod_math[,3], apply(resid_sum2, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Fixed mindset", main = "bcf summary residuals - Math epoch")
dev.off()


#Additive summary for speech
ad_sum3 <- additive_summary(
  rowMeans(t(lbcf_tau_speech))~s(bothmindsets)+s(stressmindset_baseline)+s(fixedmindset_baseline)+sex.1,
  fhatSamples=t(lbcf_tau_speech),
  fhat = rowMeans(t(lbcf_tau_speech)),
  df = as.data.frame(x_mod_speech),
  alpha = 0.05,
  # fast = TRUE,
  # quants = seq(0, 1, by = 0.005),
  # grid_size = 100,
  verbose = FALSE,
  return_samples = TRUE,
  meta = NA
)


resid_sum3 <- t(lbcf_tau_speech)-ad_sum3$fittedValues

plot(x_mod_speech[,1], apply(resid_sum3, 1, mean))
plot(x_mod_speech[,2], apply(resid_sum3, 1, mean))
plot(x_mod_speech[,3], apply(resid_sum3, 1, mean))
plot(x_mod_speech[,4], apply(resid_sum3, 1, mean))

pdf("sum_res_speech_1_lbcf.pdf")
plot(x_mod_speech[,1], apply(resid_sum3, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Bothmindsets", main = "lbcf summary residuals - Speech epoch")
dev.off()
pdf("sum_res_speech_2_lbcf.pdf")
plot(x_mod_speech[,2], apply(resid_sum3, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Stress mindset", main = "lbcf summary residuals - Speech epoch")
dev.off()
pdf("sum_res_speech_3_lbcf.pdf")
plot(x_mod_speech[,3], apply(resid_sum3, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Fixed mindset", main = "lbcf summary residuals - Speech epoch")
dev.off()


#Compared with bcf additive summary for speech
ad_sum4 <- additive_summary(
  rowMeans(t(bcf_tau_speech))~s(bothmindsets)+s(stressmindset_baseline)+s(fixedmindset_baseline)+sex.1,
  fhatSamples=t(bcf_tau_speech),
  fhat = rowMeans(t(bcf_tau_speech)),
  df = as.data.frame(x_mod_speech),
  alpha = 0.05,
  # fast = TRUE,
  # quants = seq(0, 1, by = 0.005),
  # grid_size = 100,
  verbose = FALSE,
  return_samples = TRUE,
  meta = NA
)

resid_sum4 <- t(bcf_tau_speech)-ad_sum4$fittedValues

plot(x_mod_speech[,1], apply(resid_sum4, 1, mean))
plot(x_mod_speech[,2], apply(resid_sum4, 1, mean))
plot(x_mod_speech[,3], apply(resid_sum4, 1, mean))
plot(x_mod_speech[,4], apply(resid_sum4, 1, mean))

pdf("sum_res_speech_1_bcf.pdf")
plot(x_mod_speech[,1], apply(resid_sum4, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Bothmindsets", main = "bcf summary residuals - Speech epoch")
dev.off()
pdf("sum_res_speech_2_bcf.pdf")
plot(x_mod_speech[,2], apply(resid_sum4, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Stress mindset", main = "bcf summary residuals - Speech epoch")
dev.off()
pdf("sum_res_speech_3_bcf.pdf")
plot(x_mod_speech[,3], apply(resid_sum4, 1, mean), 
     ylab = "Summary residuals mean", xlab = "Fixed mindset", main = "bcf summary residuals - Speech epoch")
dev.off()


#Fixed has upward trend
#Sex has downward trend
#stress is more or less stable across all values

#lbcf math
pdf("gams_math_1_lbcf.pdf")
additive_summary_plot(ad_sum) + ylim(c(-600,500)) + 
  plot_annotation(title = "GAMs summary plots - Math epoch - lbcf")
dev.off()
#bcf math
pdf("gams_math_1_bcf.pdf")
additive_summary_plot(ad_sum2) + ylim(c(-1000,1000)) + 
  plot_annotation(title = "GAMs summary plots - Math epoch - bcf")
dev.off()
#lbcf speech
pdf("gams_speech_1_lbcf.pdf")
additive_summary_plot(ad_sum3) + ylim(c(-600,500)) + 
  plot_annotation(title = "GAMs summary plots - Speech epoch - lbcf")
dev.off()
#bcf speech
pdf("gams_speech_1_bcf.pdf")
additive_summary_plot(ad_sum4) + ylim(c(-1400,1400)) + 
  plot_annotation(title = "GAMs summary plots - Speech epoch - bcf")
dev.off()


hist(ad_sum$summaryRsq, ylim=c(0,38000), xlim = c(0,1), breaks = 70)
hist(ad_sum2$summaryRsq, ylim=c(0,38000), xlim = c(0,1), breaks = 70)

mean(ad_sum$summaryRsq>0.95)
mean(ad_sum2$summaryRsq>0.95)


pdf("summary_r2_math.pdf")
plot(density(ad_sum$summaryRsq), col = "blue", lwd = 2,
     main = "")
lines(density(ad_sum2$summaryRsq), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf-CV summary", "bcf-CV summary"), fill = c("blue", "red"))
dev.off()

pdf("summary_r2_speech.pdf")
plot(density(ad_sum3$summaryRsq), col = "blue", lwd = 2,
     main = "")
lines(density(ad_sum4$summaryRsq), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf-CV summary", "bcf-CV summary"), fill = c("blue", "red"))
dev.off()

mean(ad_sum3$summaryRsq>0.95)
mean(ad_sum4$summaryRsq>0.95)




#Try the doublebad and doublegood subgroups
#and separate by sex

##################################
#Doublebad and doublegood for math
##################################
x_mod_math <- as.data.frame(x_mod_math)

groups1 <- ifelse(x_mod_math$stressmindset_baseline>3&x_mod_math$fixedmindset_baseline>3, 1, 0)
groups2 <- ifelse(x_mod_math$stressmindset_baseline<=3&x_mod_math$fixedmindset_baseline<=3, 1, 0)
#Getting the interaction
intr <- interaction(x_mod_math$stressmindset_baseline>3,x_mod_math$fixedmindset_baseline>3)

#Now getting the subgroups
lsubgroup_ates = subgroup_average_posterior(lbcf_tau_math, intr)
subgroup_ates = subgroup_average_posterior(bcf_tau_math, intr)

#Does not look like much
boxplot(lsubgroup_ates)

#Let us take the difference between the two groups that matter
difl <- lsubgroup_ates$FALSE.FALSE-lsubgroup_ates$TRUE.TRUE
difbcf <- subgroup_ates$FALSE.FALSE-subgroup_ates$TRUE.TRUE

#Seems like they are pretty close to zero
#Does not look like anything interesting is happening
plot(density(difbcf), ylim = c(0,0.015))
lines(density(difl), col = "red")

#Calculate the probability
#Not much is happening either
mean(difl>0)
mean(difbcf>0)

#It does not look like there is something here for doublebad and doublegood
quantile(lsubgroup_ates$FALSE.FALSE, c(0.025, 0.975))
quantile(subgroup_ates$FALSE.FALSE, c(0.025, 0.975))

quantile(lsubgroup_ates$TRUE.TRUE, c(0.025, 0.975))
quantile(subgroup_ates$TRUE.TRUE, c(0.025, 0.975))
#Quantiles tending to be leaning towards negative side
#Either way, too much uncertainty around zero to get something from this


#######################
#For sex
#######################
groups3 <- x_mod_math$sex.1

#Now getting the subgroups
lsubgroup_ates = subgroup_average_posterior(lbcf_tau_math, groups3)
subgroup_ates = subgroup_average_posterior(bcf_tau_math, groups3)

#Does not look like much
boxplot(lsubgroup_ates)

#Let us take the difference
difl <- lsubgroup_ates$X0-lsubgroup_ates$X1
difbcf <- subgroup_ates$X0-subgroup_ates$X1

#Seems like they are pretty close to zero
pdf("math_sex_diff.pdf")
plot(density(difbcf), lwd = 2, main = "", col = "blue")
lines(density(difl), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf", "bcf"), fill = c("red", "blue"))
dev.off()

#Interesting density curve for linear
#spreading way more
#Calculate the probabilities
mean(difl>0)
mean(difbcf>0)

quantile(difbcf, c(0.025, 0.975))
quantile(difl, c(0.025, 0.975))


##################################
#Doublebad and doublegood for speech
##################################
x_mod_speech <- as.data.frame(x_mod_speech)

groups1 <- ifelse(x_mod_speech$stressmindset_baseline>3&x_mod_speech$fixedmindset_baseline>3, 1, 0)
groups2 <- ifelse(x_mod_speech$stressmindset_baseline<=3&x_mod_speech$fixedmindset_baseline<=3, 1, 0)
#Getting the interaction
intr <- interaction(x_mod_speech$stressmindset_baseline>3,x_mod_speech$fixedmindset_baseline>3)

#Now getting the subgroups
lsubgroup_ates = subgroup_average_posterior(lbcf_tau_speech, intr)
subgroup_ates = subgroup_average_posterior(bcf_tau_speech, intr)

#Does not look like much
boxplot(lsubgroup_ates)

#Let us take the difference
difl <- lsubgroup_ates$FALSE.FALSE-lsubgroup_ates$TRUE.TRUE
difbcf <- subgroup_ates$FALSE.FALSE-subgroup_ates$TRUE.TRUE

#Seems like they are pretty close to zero
#Again, not much is happening
plot(density(difbcf), ylim = c(0,0.015))
lines(density(difl), col = "red")

#Calculate the probability
mean(difl>0)
mean(difbcf>0)

quantile(difl, c(0.025, 0.975))
quantile(difbcf, c(0.025, 0.975))

#It does not look like there is something here for doublebad and doublegood
quantile(lsubgroup_ates$FALSE.FALSE, c(0.025, 0.975))
quantile(subgroup_ates$FALSE.FALSE, c(0.025, 0.975))

quantile(lsubgroup_ates$TRUE.TRUE, c(0.025, 0.975))
quantile(subgroup_ates$TRUE.TRUE, c(0.025, 0.975))

#######################
#For sex
#######################
groups3 <- x_mod_speech$sex.1

#Now getting the subgroups
lsubgroup_ates = subgroup_average_posterior(lbcf_tau_speech, groups3)
subgroup_ates = subgroup_average_posterior(bcf_tau_speech, groups3)

#Does not look like much
boxplot(lsubgroup_ates)

pdf("sex_speech_lbcf.pfd")
plot(density(lsubgroup_ates$X0), lwd = 2, xlim = c(-800,1000),
     main = "", col = "blue")
lines(density(lsubgroup_ates$X1), col = "red", lwd = 2)
legend("topleft", legend = c("Male", "Female"), fill = c("blue", "red"))
dev.off()

pdf("sex_speech_bcf.pfd")
plot(density(subgroup_ates$X0), lwd = 2, xlim = c(-800,1000),
     main = "", col = "blue")
lines(density(subgroup_ates$X1), col = "red", lwd = 2)
legend("topleft", legend = c("Male", "Female"), fill = c("blue", "red"))
dev.off()


#Let us take the difference
difl <- lsubgroup_ates$X0-lsubgroup_ates$X1
difbcf <- subgroup_ates$X0-subgroup_ates$X1

#Seems like they are pretty close to zero
pdf("speech_sex_diff.pdf")
plot(density(difbcf), lwd = 2, main = "", col = "blue")
lines(density(difl), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf", "bcf"), fill = c("red", "blue"))
dev.off()


#Interesting density curve for linear

#Calculate the probabilities
mean(difl>0)
mean(difbcf>0)

quantile(difl, c(0.025, 0.975))
quantile(difbcf, c(0.025, 0.975))

#Comparing variables with ITEs for math
plot(x_mod_math$bothmindsets, apply(lbcf_tau_math,2,mean), ylim = c(-100,70))
plot(x_mod_math$bothmindsets, apply(bcf_tau_math,2,mean), ylim = c(-100,70))

plot(x_mod_math$stressmindset_baseline, apply(lbcf_tau_math,2,mean), ylim = c(-100,70))
plot(x_mod_math$stressmindset_baseline, apply(bcf_tau_math,2,mean), ylim = c(-100,70))

plot(x_mod_math$fixedmindset_baseline, apply(lbcf_tau_math,2,mean), ylim = c(-100,70))
plot(x_mod_math$fixedmindset_baseline, apply(bcf_tau_math,2,mean), ylim = c(-100,70))

plot(x_mod_math$sex.1, apply(lbcf_tau_math,2,mean), ylim = c(-100,70))
plot(x_mod_math$sex.1, apply(bcf_tau_math,2,mean), ylim = c(-100,70))



#Comparing variables with ITEs for speech
plot(x_mod_speech$bothmindsets, apply(lbcf_tau_speech,2,mean), ylim = c(-150,120))
plot(x_mod_speech$bothmindsets, apply(bcf_tau_speech,2,mean), ylim = c(-150,120))

plot(x_mod_speech$stressmindset_baseline, apply(lbcf_tau_speech,2,mean), ylim = c(-150,120))
plot(x_mod_speech$stressmindset_baseline, apply(bcf_tau_speech,2,mean), ylim = c(-150,120))

plot(x_mod_speech$fixedmindset_baseline, apply(lbcf_tau_speech,2,mean), ylim = c(-150,120))
plot(x_mod_speech$fixedmindset_baseline, apply(bcf_tau_speech,2,mean), ylim = c(-150,120))

plot(x_mod_speech$sex.1, apply(lbcf_tau_speech,2,mean), ylim = c(-150,120))
plot(x_mod_speech$sex.1, apply(bcf_tau_speech,2,mean), ylim = c(-150,120))



#########################
#Finding the best cuts
#########################
search_dual <- function(y, data, grid_idx, grid1, grid2, prop = 0.2){
  n <- nrow(data)
  min_n <- round(n*prop,0)
  
  fixed_obs <- data$fixedmindset_baseline
  stress_obs <- data$stressmindset_baseline
  
  dif <- numeric(0)
  for(i in 1:nrow(grid_idx)){
    grid_aux <- as.numeric(grid_idx[i,])
    
    x1_aux <- grid1[grid_aux[1],]
    x2_aux <- grid2[grid_aux[2],]
    
    #The groups
    g1 <- numeric(0)
    g2 <- numeric(0)
    
    #Separating the groups
    g1 <- (data$fixedmindset_baseline<(x1_aux[1]))&((data)$stressmindset_baseline<(x2_aux[1]))
    g2 <- (data$fixedmindset_baseline>=(x1_aux[2]))&((data)$stressmindset_baseline>=(x2_aux[2]))
    
    g_test <- (sum(g1)>=min_n)&(sum(g2)>=min_n)
    
    #calculating the difference
    if(g_test==T){
      #Calculate the difference of means the groups
      dif[i] <- (mean(y[g1])-mean(y[g2]))
    }else{
      dif[i] <- 0
    }
    
  }
  
  return(c(dif))
}

#Function to find available cutpoints
middle <- function(x){
  aux <- numeric(0)
  y <- unique(sort(x))
  for(i in 1:(length(y)-1)){
    aux[i] <- (y[i]+y[(i+1)])/2
  }
  return(aux)
}

#Unique expand grid (https://stackoverflow.com/questions/17171148/non-redundant-version-of-expand-grid)
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#The available cutpoints
cuts_fixed <- middle(x_mod_math$fixedmindset_baseline)
cuts_stress <- middle(x_mod_math$stressmindset_baseline)

#Create grid of a<b
grid_x1 <- expand.grid.unique(cuts_fixed,cuts_fixed)
grid_x2 <- expand.grid.unique(cuts_stress,cuts_stress)

#And include the case where a=b
grid_x1 <- rbind(grid_x1,cbind(cuts_fixed,cuts_fixed))
grid_x2 <- rbind(grid_x2,cbind(cuts_stress,cuts_stress))

#Creating a grid of cutpoints
#In this case a have matrices, and I want to try combinations of the lines
#So I will create a grid of indexes and use those

ind_1 <- 1:(nrow(grid_x1))
ind_2 <- 1:(nrow(grid_x2))

grid_index <- expand.grid(ind_1,ind_2)


#Separete by using only math epoch
dif_vec <- search_dual(y = apply(lbcf_tau_math,2,mean), 
                       data = x_mod_math, grid_idx = grid_index, grid1 = grid_x1, grid2 = grid_x2, 
                       prop = 0.15)
indexes_max <- dif_vec==max(dif_vec)
fixed_cut = grid_x1[as.numeric((grid_index[indexes_max,])[,1]), , drop = F]
stress_cut = grid_x2[as.numeric((grid_index[indexes_max,])[,2]), , drop = F]

cuts_list <- list()
for(i in 1:sum(indexes_max)){
  cuts_list[[i]] <- data.frame(fixed = fixed_cut[i,], stress = stress_cut[i,])
}

cuts <- cuts_list[[1]]

doublegood <- (x_mod_math$fixedmindset_baseline<(cuts[1,1]))&(x_mod_math$stressmindset_baseline<(cuts[1,2]))
doublebad <- (x_mod_math$fixedmindset_baseline>=(cuts[2,1]))&(x_mod_math$stressmindset_baseline>=(cuts[2,2]))

#########################
#########################

intr <- interaction(doublegood,doublebad)

library(possum)
lsubgroup_ates = subgroup_average_posterior(lbcf_tau_math, intr)
subgroup_ates = subgroup_average_posterior(bcf_tau_math, intr)
boxplot(lsubgroup_ates)
boxplot(subgroup_ates)
difl <- lsubgroup_ates$TRUE.FALSE-lsubgroup_ates$FALSE.TRUE
difbcf <- subgroup_ates$TRUE.FALSE-subgroup_ates$FALSE.TRUE

#Doublebad, both mindsets are high
hist(lsubgroup_ates$FALSE.TRUE, breaks = 70)
hist(subgroup_ates$FALSE.TRUE, breaks = 90, add = T, col = "red")

#Doublegood, both mindsets are low
hist(lsubgroup_ates$TRUE.FALSE, breaks = 70)
hist(subgroup_ates$TRUE.FALSE, breaks = 70, add = T, col = "red")

pdf("math_priormind_diff.pdf")
plot(density(difbcf), lwd = 2, main = "", col = "blue")
lines(density(difl), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf", "bcf"), fill = c("red", "blue"))
dev.off()

mean(difbcf>0)
mean(difl>0)

quantile(difbcf, c(0.1, 0.9))
quantile(difl, c(0.1, 0.9))


#######################################

#The available cutpoints
cuts_fixed <- middle(x_mod_speech$fixedmindset_baseline)
cuts_stress <- middle(x_mod_speech$stressmindset_baseline)

#Create grid of a<b
grid_x1 <- expand.grid.unique(cuts_fixed,cuts_fixed)
grid_x2 <- expand.grid.unique(cuts_stress,cuts_stress)

#And include the case where a=b
grid_x1 <- rbind(grid_x1,cbind(cuts_fixed,cuts_fixed))
grid_x2 <- rbind(grid_x2,cbind(cuts_stress,cuts_stress))

#Creating a grid of cutpoints
#In this case a have matrices, and I want to try combinations of the lines
#So I will create a grid of indexes and use those

ind_1 <- 1:(nrow(grid_x1))
ind_2 <- 1:(nrow(grid_x2))

grid_index <- expand.grid(ind_1,ind_2)


#Separete by using only speech epoch
dif_vec <- search_dual(y = apply(lbcf_tau_speech,2,mean), 
                       data = x_mod_speech, grid_idx = grid_index, grid1 = grid_x1, grid2 = grid_x2, 
                       prop = 0.2)
indexes_max <- dif_vec==max(dif_vec)
fixed_cut = grid_x1[as.numeric((grid_index[indexes_max,])[,1]), , drop = F]
stress_cut = grid_x2[as.numeric((grid_index[indexes_max,])[,2]), , drop = F]

cuts_list <- list()
for(i in 1:sum(indexes_max)){
  cuts_list[[i]] <- data.frame(fixed = fixed_cut[i,], stress = stress_cut[i,])
}

cuts <- cuts_list[[1]]

doublegood <- (x_mod_speech$fixedmindset_baseline<(cuts[1,1]))&(x_mod_speech$stressmindset_baseline<(cuts[1,2]))
doublebad <- (x_mod_speech$fixedmindset_baseline>=(cuts[2,1]))&(x_mod_speech$stressmindset_baseline>=(cuts[2,2]))

#########################
#########################

intr <- interaction(doublegood,doublebad)

library(possum)
lsubgroup_ates = subgroup_average_posterior(lbcf_tau_speech, intr)
subgroup_ates = subgroup_average_posterior(bcf_tau_speech, intr)
boxplot(lsubgroup_ates)
boxplot(subgroup_ates)
difl <- lsubgroup_ates$TRUE.FALSE-lsubgroup_ates$FALSE.TRUE
difbcf <- subgroup_ates$TRUE.FALSE-subgroup_ates$FALSE.TRUE

#Doublebad, both mindsets are high
hist(lsubgroup_ates$FALSE.TRUE, breaks = 70)
hist(subgroup_ates$FALSE.TRUE, breaks = 90, add = T, col = "red")

#Doublegood, both mindsets are low
hist(lsubgroup_ates$TRUE.FALSE, breaks = 70)
hist(subgroup_ates$TRUE.FALSE, breaks = 70, add = T, col = "red")


plot(density(difbcf))
lines(density(difl), col = "red")

pdf("speech_priormind_diff.pdf")
plot(density(difbcf), lwd = 2, main = "", col = "blue")
lines(density(difl), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf", "bcf"), fill = c("red", "blue"))
dev.off()

#Calculate the probabilities
mean(difl>0)
mean(difbcf>0)

quantile(difbcf, c(0.1, 0.9))
quantile(difl, c(0.1, 0.9))

hist(apply(lbcf_tau_speech,2,mean), col = x_mod_speech$sex.1, breaks = 100)
hist(apply(lbcf_tau_speech,2,mean), col = doublegood, breaks = 100)
hist(apply(lbcf_tau_speech,2,mean), col = doublebad, breaks = 100)
hist(apply(lbcf_tau_math,2,mean), col = x_mod_speech$sex.1, breaks = 100)
hist(apply(lbcf_tau_math,2,mean), col = doublegood, breaks = 100)
hist(apply(lbcf_tau_math,2,mean), col = doublebad, breaks = 100)



plot(density(apply(lbcf_tau_speech[,x_mod_speech$sex.1==0],1,mean)))
lines(density(apply(lbcf_tau_speech[,x_mod_speech$sex.1==1],1,mean)), col = "red")

plot(density(apply(lbcf_tau_speech[,doublegood],1,mean)))
lines(density(apply(lbcf_tau_speech[,doublebad],1,mean)), col = "red")

plot(density(apply(lbcf_tau_speech[,as.logical(doublebad*(x_mod_speech$sex.1==1))],1,mean)))
lines(density(apply(lbcf_tau_speech[,as.logical(doublegood*(x_mod_speech$sex.1==1))],1,mean)), col = "red")

plot(density(apply(lbcf_tau_speech[,as.logical(doublebad*(x_mod_speech$sex.1==0))],1,mean)))
lines(density(apply(lbcf_tau_speech[,as.logical(doublegood*(x_mod_speech$sex.1==0))],1,mean)), col = "red")


##################################
#RMSE and 5-fold CV

testando <- df_speech %>% arrange(rmse_test)
bcfdefault <- testando %>% filter(nt_con==200, nt_mod==50, 
                                  a_con==0.95, a_mod==0.8, 
                                  b_con==2, b_mod==3, 
                                  muscale=="true", tauscale=="true", 
                                  linear=="none")
lbcfdefault <- testando %>% filter(nt_con==100, nt_mod==25, 
                                   a_con==0.5, a_mod==0.8, 
                                   b_con==1, b_mod==3, 
                                   muscale=="true", tauscale=="true", 
                                   linear=="both")
pdf("rmse_speech_full.pdf")
par(mfrow=c(1,1))
plot(testando$rmse_test, type = "p", cex = 0.1, ylim = c(520,830),
     ylab = "RMSE", xlab = "")
points((1:length(testando$rmse_train)),
       testando$rmse_train,
       pch = 19, cex = 0.2, col = "purple")
points(which(testando$rmse_test==bcfdefault$rmse_test),bcfdefault$rmse_test, 
       col = "red", pch=19)
points(which(testando$rmse_test==bcfdefault$rmse_test),bcfdefault$rmse_train, 
       col = "red", pch=19)
points(which(testando$rmse_test==lbcfdefault$rmse_test),lbcfdefault$rmse_test, 
       col = "blue", pch=19)
points(which(testando$rmse_test==lbcfdefault$rmse_test),lbcfdefault$rmse_train, 
       col = "blue", pch=19)
legend("topleft", legend = c("lbcf-default", "bcf-default"), fill = c("blue", "red"))
dev.off()

pdf("rmse_speech_partial.pdf")
par(mfrow=c(2,2))
plot(testando$rmse_test, type = "l", ylim = c(520,830),
     ylab = "RMSE", xlab = "")
points((1:length(testando$rmse_train))[testando$linear=="none"],
       testando$rmse_train[testando$linear=="none"],
       pch = 19, cex = 0.2, col = "red")
legend("topleft", legend = c("bcf"), fill = c("red"))
plot(testando$rmse_test, type = "l", ylim = c(520,830),
     ylab = "RMSE", xlab = "")
points((1:length(testando$rmse_train))[testando$linear=="moderate"],
       testando$rmse_train[testando$linear=="moderate"],
       pch = 19, cex = 0.2, col = "purple")
legend("topleft", legend = c("partial lbcf (moderate)"), fill = c("purple"))
plot(testando$rmse_test, type = "l", ylim = c(520,830),
     ylab = "RMSE", xlab = "")
points((1:length(testando$rmse_train))[testando$linear=="control"],
       testando$rmse_train[testando$linear=="control"],
       pch = 19, cex = 0.2, col = "blue")
legend("topleft", legend = c("partial lbcf (control)"), fill = c("blue"))
plot(testando$rmse_test, type = "l", ylim = c(520,830),
     ylab = "RMSE", xlab = "")
points((1:length(testando$rmse_train))[testando$linear=="both"],
       testando$rmse_train[testando$linear=="both"],
       pch = 19, cex = 0.2, col = "orange")
legend("topleft", legend = c("lbcf"), fill = c("orange"))
dev.off()


testando2 <- df_math %>% arrange(rmse_test)

bcfdefault2 <- testando2 %>% filter(nt_con==200, nt_mod==50, 
                                  a_con==0.95, a_mod==0.8, 
                                  b_con==2, b_mod==3, 
                                  muscale=="true", tauscale=="true", 
                                  linear=="none")
lbcfdefault2 <- testando2 %>% filter(nt_con==100, nt_mod==25, 
                                   a_con==0.5, a_mod==0.8, 
                                   b_con==1, b_mod==3, 
                                   muscale=="true", tauscale=="true", 
                                   linear=="both")

pdf("rmse_math_full.pdf")
par(mfrow=c(1,1))
plot(testando2$rmse_test, type = "p", cex = 0.1, ylim = c(500,750),
     ylab = "RMSE", xlab = "")
points((1:length(testando2$rmse_train)),
       testando2$rmse_train,
       pch = 19, cex = 0.2, col = "purple")
points(which(testando2$rmse_test==bcfdefault2$rmse_test),bcfdefault2$rmse_test, 
       col = "red", pch=19)
points(which(testando2$rmse_test==bcfdefault2$rmse_test),bcfdefault2$rmse_train, 
       col = "red", pch=19)
points(which(testando2$rmse_test==lbcfdefault2$rmse_test),lbcfdefault2$rmse_test, 
       col = "blue", pch=19)
points(which(testando2$rmse_test==lbcfdefault2$rmse_test),lbcfdefault2$rmse_train, 
       col = "blue", pch=19)
legend("topleft", legend = c("lbcf-default", "bcf-default"), fill = c("blue", "red"))
dev.off()




pdf("rmse_math_partial.pdf")
par(mfrow=c(2,2))
plot(testando2$rmse_test, type = "l", ylim = c(500,750),
     ylab = "RMSE", xlab = "")
points((1:length(testando2$rmse_train))[testando2$linear=="none"],
       testando2$rmse_train[testando2$linear=="none"],
       pch = 19, cex = 0.2, col = "red")
legend("topleft", legend = c("bcf"), fill = c("red"))
plot(testando2$rmse_test, type = "l", ylim = c(500,750),
     ylab = "RMSE", xlab = "")
points((1:length(testando2$rmse_train))[testando2$linear=="moderate"],
       testando2$rmse_train[testando2$linear=="moderate"],
       pch = 19, cex = 0.2, col = "purple")
legend("topleft", legend = c("partial lbcf (moderate)"), fill = c("purple"))
plot(testando2$rmse_test, type = "l", ylim = c(500,750),
     ylab = "RMSE", xlab = "")
points((1:length(testando2$rmse_train))[testando2$linear=="control"],
       testando2$rmse_train[testando2$linear=="control"],
       pch = 19, cex = 0.2, col = "blue")
legend("topleft", legend = c("partial lbcf (control)"), fill = c("blue"))
plot(testando2$rmse_test, type = "l", ylim = c(500,750),
     ylab = "RMSE", xlab = "")
points((1:length(testando2$rmse_train))[testando2$linear=="both"],
       testando2$rmse_train[testando2$linear=="both"],
       pch = 19, cex = 0.2, col = "orange")
legend("topleft", legend = c("lbcf"), fill = c("orange"))
dev.off()



