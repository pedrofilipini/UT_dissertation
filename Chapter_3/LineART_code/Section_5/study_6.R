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

#bfi is personality
#ses is social class
#test anxiety is Highly test-anxious adolescents could be lower-performing and could be more likely to show negative stress reactivity

Xm <- study6_clean %>% dplyr::select(sex, stressmindset_baseline, fixedmindset_baseline, pss_impute )%>% dbarts::makeModelMatrixFromDataFrame()


z = study6_clean$treatment # "Treatment" variable - Actually a measured variable
y = study6_clean$anxiety

npred <- length(y)
burn <- 15000
th <- 3
post <- 5000

cate <- list()
ite <- list()
lcate <- list()
lite <- list()
rep <- 10

model_lbcf <- list()
model_bcf <- list()

load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart/cv_anxiety/df_cvanx.RData")

# for(i in 1:rep){
#   set.seed(1+i)
#   lbcf_default <- lbart:::bcf_new(y, #response
#                                   z = z, z_est = rep(1,npred), #Z
#                                   x_control = Xc, #matrix of covariates
#                                   x_control_est = Xc,
#                                   x_moderate = Xm,
#                                   x_moderate_est = Xm,
#                                   pihat = rep(1,npred), #in case of causal inference
#                                   pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
#                                   include_pi = "none", #can be "none", "control", "moderate", or "both"
#                                   linear = "moderate", #can be "none", "control", "moderate", or "both"
#                                   base_control = 0.5, base_moderate = 0.25, #hyperparameters
#                                   power_control = 3, power_moderate = 2, #hyperparameters
#                                   nburn = burn, nthin = th, nsim = post, #draws
#                                   ntree_control = 50, #control trees
#                                   ntree_moderate = 25, #moderate trees (treatment effect)
#                                   dart = F, #Linero's Prior
#                                   save_trees = T, 
#                                   use_muscale = T, 
#                                   use_tauscale = F)
#   model_lbcf[[i]] <- lbcf_default
# 
# 
#   set.seed(1+i)
#   bcf_default <- lbart:::bcf_new(y, #response
#                                  z = z, z_est = rep(1,npred), #Z
#                                  x_control = Xc, #matrix of covariates
#                                  x_control_est = Xc,
#                                  x_moderate = Xm,
#                                  x_moderate_est = Xm,
#                                  pihat = rep(1,npred), #in case of causal inference
#                                  pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
#                                  include_pi = "none", #can be "none", "control", "moderate", or "both"
#                                  linear = "none", #can be "none", "control", "moderate", or "both"
#                                  base_control = 0.5, base_moderate = 0.5, #hyperparameters
#                                  power_control = 2, power_moderate = 3, #hyperparameters
#                                  nburn = burn, nthin = th, nsim = post, #draws
#                                  ntree_control = 50, #control trees
#                                  ntree_moderate = 50, #moderate trees (treatment effect)
#                                  dart = F, #Linero's Prior
#                                  save_trees = T, 
#                                  use_muscale = F, 
#                                  use_tauscale = F)
#   model_bcf[[i]] <- bcf_default
# 
# }


# save(model_lbcf, file="lbcf_study6.RData")
# save(model_bcf, file="bcf_study6.RData")

load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/lbcf_study6.RData")
load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/bcf_study6.RData")


{bcf_mu <- NULL
bcf_tau <- NULL
bcf_eta_con <- NULL
bcf_eta_mod <- NULL
bcf_sigma <- NULL
lbcf_mu <- NULL
lbcf_tau <- NULL
lbcf_eta_con <- NULL
lbcf_eta_mod <- NULL
lbcf_sigma <- NULL}

for(i in 1:rep){
  bcf_mu <- rbind(bcf_mu, model_bcf[[i]]$mu_est_post/sd(y))
  bcf_tau <- rbind(bcf_tau, model_bcf[[i]]$b_est_post/sd(y))
  bcf_eta_con <- c(bcf_eta_con, model_bcf[[i]]$eta_con)
  bcf_eta_mod <- c(bcf_eta_mod, model_bcf[[i]]$eta_mod)
  bcf_sigma <- c(bcf_sigma, model_bcf[[i]]$sigma)
  lbcf_mu <- rbind(lbcf_mu, model_lbcf[[i]]$mu_est_post/sd(y))
  lbcf_tau <- rbind(lbcf_tau, model_lbcf[[i]]$b_est_post/sd(y))
  lbcf_eta_con <- c(lbcf_eta_con, model_lbcf[[i]]$eta_con)
  lbcf_eta_mod <- c(lbcf_eta_mod, model_lbcf[[i]]$eta_mod)
  lbcf_sigma <- c(lbcf_sigma, model_lbcf[[i]]$sigma)
}

#Comparing some stuff
par(mfrow=c(1,2))
plot(lbcf_sigma/sd(y), ylim = c(0.5,1), type = "l")
plot(bcf_sigma/sd(y), ylim = c(0.5,1), type = "l")

par(mfrow=c(1,2))
plot(lbcf_eta_con, type = "l")
plot(bcf_eta_con, type = "l")

par(mfrow=c(1,2))
plot(lbcf_eta_mod, type = "l")
plot(bcf_eta_mod, type = "l")

#Across chains seems to be the same, also similar values for bcf and lbcf

#Getting all the chains together
lite <- data.frame(mean = apply(lbcf_tau,2, mean),
                       lb = apply(lbcf_tau,2, quantile, 0.025),
                       ub = apply(lbcf_tau,2, quantile, 0.975))
lcate <- data.frame(mean = apply(lbcf_tau,1, mean),
                        lb = apply(lbcf_tau,1, quantile, 0.025),
                        ub = apply(lbcf_tau,1, quantile, 0.975))
ite <- data.frame(mean = apply(bcf_tau,2, mean),
           lb = apply(bcf_tau,2, quantile, 0.025),
           ub = apply(bcf_tau,2, quantile, 0.975))
cate <- data.frame(mean = apply(bcf_tau,1, mean),
                       lb = apply(bcf_tau,1, quantile, 0.025),
                       ub = apply(bcf_tau,1, quantile, 0.975))



par(mfrow=c(1,2))
boxplot(lcate$mean, cate$mean,
        names = c("lbcf", "bcf"), main = "Posterior ATE")
boxplot(lite$mean, ite$mean,
        names = c("lbcf", "bcf"), main = "Posterior Mean ITE")

par(mfrow=c(1,1))
#pdf("density_ite_anx.pdf")
plot(density(lite$mean), ylim = c(0,3),
     main = "Density of ITE posterior mean - Anxiety")
lines(density(ite$mean), col = "red")
legend("topleft", legend = c("lbcf", "bcf"), fill = c("black", "red"))
#dev.off()

#pdf("density_ate_anx.pdf")
plot(density(lcate$mean), ylim = c(0,6),
     main = "Density of ATE posterior - Anxiety")
lines(density(cate$mean), col = "red")
legend("topleft", legend = c("lbcf", "bcf"), fill = c("black", "red"))
#dev.off()



Xm <- as.data.frame(Xm)
plot(Xm$fixedmindset_baseline, lite$mean,
     main = "Posterior Mean ITE")
points(Xm$fixedmindset_baseline, ite$mean, col ="red",
     main = "Posterior Mean ITE")

plot(Xm$stressmindset_baseline, lite$mean,
     main = "Posterior Mean ITE")
points(Xm$stressmindset_baseline, ite$mean, col ="red",
       main = "Posterior Mean ITE")

plot(Xm$sex, lite$mean,
     main = "Posterior Mean ITE")
points(Xm$sex, ite$mean, col ="red",
       main = "Posterior Mean ITE")

plot(Xm$pss_impute, lite$mean,
     main = "Posterior Mean ITE")
points(Xm$pss_impute, ite$mean, col ="red",
       main = "Posterior Mean ITE")


diff_t <- lite$mean- ite$mean
plot(Xm$pss_impute, diff_t)
plot(Xm$stressmindset_baseline, diff_t)
plot(Xm$fixedmindset_baseline, diff_t)
plot(Xm$sex, diff_t)

plot(lite$mean, ite$mean,
     xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
abline(0,1)

plot(apply(lbcf_tau,2, mean), apply(bcf_tau,2, mean),
     xlim=c(-1.5,1.5),ylim=c(-1.5,1.5), col = study6_clean$sex+1, pch = 19)
abline(0,1)


boxplot(apply(bcf_tau,1, mean)[study6_clean$sex==1],apply(lbcf_tau,1, mean)[study6_clean$sex==1])
boxplot(apply(bcf_tau,1, mean)[study6_clean$sex==0],apply(lbcf_tau,1, mean)[study6_clean$sex==0])

boxplot(apply(bcf_tau,1, mean),apply(lbcf_tau,1, mean))
abline(h=0, lwd = 3, col = "red")



groups1 <- ifelse(study6_clean$stressmindset_baseline>3&study6_clean$fixedmindset_baseline>3, 1, 0)
groups2 <- ifelse(study6_clean$stressmindset_baseline<=3&study6_clean$fixedmindset_baseline<=3, 1, 0)

intr <- interaction(study6_clean$stressmindset_baseline>3,study6_clean$fixedmindset_baseline>3)

library(possum)
lsubgroup_ates = subgroup_average_posterior(lbcf_tau, intr)
subgroup_ates = subgroup_average_posterior(bcf_tau, intr)
boxplot(lsubgroup_ates)
difl <- lsubgroup_ates$FALSE.FALSE-lsubgroup_ates$TRUE.TRUE
difbcf <- subgroup_ates$FALSE.FALSE-subgroup_ates$TRUE.TRUE


hist(lsubgroup_ates$FALSE.FALSE, breaks = 70)
hist(subgroup_ates$FALSE.FALSE, breaks = 70)

#Negative prior mindsets
#negative fixed mindset beliefs, 
#or the belief that intellectual ability is fixed and cannot change, 
#which can lead to the appraisal that negative events are uncontrollable 
#and harmful
#stress-is-debilitating mindset, which is the belief that stress is 
#inherently negative and compromises performance, health and well-being
#Negative prior mindsets
#At baseline, participants in all experiments except study 2 
#completed standard measures of negative event-focused mindsets 
#(fixed mindset of intelligence; that is, “Your intelligence is 
#something about you that you can't change very much”) and 
#response-focused mindsets (the stress-is-debilitating mindset; 
#that is, “The overall effect of stress on my life is negative”) 
#(for both, 1=strongly disagree, 6=strongly agree). 
hist(lsubgroup_ates$TRUE.TRUE, breaks = 70)
hist(subgroup_ates$TRUE.TRUE, breaks = 70)


cor(Xm$pss_impute, Xm$fixedmindset_baseline)
plot(Xm$pss_impute, Xm$fixedmindset_baseline)


plot(density(difbcf))
lines(density(difl), col = "red")

mean(difbcf>0)
mean(difl>0)

quantile(difbcf, c(0.1, 0.9))
quantile(difl, c(0.1, 0.9))



groups <- study6_clean$sex

lsubgroup_ates = subgroup_average_posterior(lbcf_tau, groups)
subgroup_ates = subgroup_average_posterior(bcf_tau, groups)


dif1 <- subgroup_ates$X0-subgroup_ates$X1
dif2 <- lsubgroup_ates$X0-lsubgroup_ates$X1

quantile(dif1, c(0.025, 0.5, 0.975))

#pdf("anx_sex_diff.pdf")
par(mfrow=c(1,1))
plot(density(dif1), lwd = 2, main = "", col = "blue")
lines(density(dif2), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf", "bcf"), fill = c("red", "blue"))
#dev.off()

mean(dif1>0)
mean(dif2>0)


boxplot(subgroup_ates)
abline(h=0, col="red", lwd = 2)

l_subgroup_ates = subgroup_average_posterior(lbcf_tau, groups)

boxplot(l_subgroup_ates)
abline(h=0, col="red", lwd = 2)





ad_sum <- additive_summary(
  rowMeans(t(lbcf_tau))~s(stressmindset_baseline)+
    s(fixedmindset_baseline)+sex+s(pss_impute),
  fhatSamples=t(lbcf_tau),
  fhat = rowMeans(t(lbcf_tau)),
  df = study6_clean,
  alpha = 0.05,
  fast = TRUE,
  quants = seq(0, 1, by = 0.005),
  grid_size = 100,
  verbose = FALSE,
  return_samples = TRUE,
  meta = NA
)

ad_sum2 <- additive_summary(
  rowMeans(t(bcf_tau))~s(stressmindset_baseline)+
    s(fixedmindset_baseline)+sex+s(pss_impute),
  fhatSamples=t(bcf_tau),
  fhat = rowMeans(t(bcf_tau)),
  df = study6_clean,
  alpha = 0.05,
  fast = TRUE,
  quants = seq(0, 1, by = 0.005),
  grid_size = 100,
  verbose = FALSE,
  return_samples = TRUE,
  meta = NA
)




hist(ad_sum$summaryRsq, ylim=c(0,38000))
hist(ad_sum2$summaryRsq, ylim=c(0,38000))

#pdf("gams_anx_1_lbcf.pdf")
additive_summary_plot(ad_sum) + ylim(c(-1,1.3))
#dev.off()


temp <- ad_sum$gamDf[-c(1:201),]
ribbonFill = "lightsalmon"
pdf("gams_anx_1_lbcf.pdf", width = 8, height = 5)
temp %>% distinct() %>% ggplot() + geom_hline(yintercept = 0) + 
  geom_ribbon(aes(x_j, ymin = fx_j_lo, ymax = fx_j_hi), 
              fill = ribbonFill, alpha = 0.5) + 
  geom_line(aes(x_j,fx_j_mean), col = "firebrick3") + 
  geom_rug(aes(x_j, fx_j_mean),sides = "b", alpha = 0.25) + 
  facet_wrap(~term, scale = "free_x") + 
  labs(x = (""), y = ("Partial effect")) +
  ylim(c(-1,1.3))
dev.off()

temp <- ad_sum2$gamDf[-c(1:201),]
ribbonFill = "lightsalmon"
pdf("gams_anx_1_bcf.pdf", width = 8, height = 5)
temp %>% distinct() %>% ggplot() + geom_hline(yintercept = 0) + 
  geom_ribbon(aes(x_j, ymin = fx_j_lo, ymax = fx_j_hi), 
              fill = ribbonFill, alpha = 0.5) + 
  geom_line(aes(x_j,fx_j_mean), col = "firebrick3") + 
  geom_rug(aes(x_j, fx_j_mean),sides = "b", alpha = 0.25) + 
  facet_wrap(~term, scale = "free_x") + 
  labs(x = (""), y = ("Partial effect")) +
  ylim(c(-1,1.3))
dev.off()


#pdf("gams_anx_1_bcf.pdf")
additive_summary_plot(ad_sum2) + ylim(c(-1,1.3))
#dev.off()


#pdf("summary_r2_anx.pdf")
plot(density(ad_sum$summaryRsq), col = "blue", lwd = 2,
     main = "", ylim = c(0,12))
lines(density(ad_sum2$summaryRsq), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf-CV summary", "bcf-CV summary"), fill = c("blue", "red"))
#dev.off()

res <- (t(bcf_tau)-ad_sum$fittedValues)


plot(Xm$stressmindset_baseline, apply(res,1,mean))
plot(Xm$fixedmindset_baseline, apply(res,1,mean))
plot(Xm$sex, apply(res,1,mean))
plot(Xm$pss_impute, apply(res,1,mean))

res2 <- (t(lbcf_tau)-ad_sum$fittedValues)
plot(Xm$stressmindset_baseline, apply(res2,1,mean))
plot(Xm$fixedmindset_baseline, apply(res2,1,mean))
plot(Xm$sex, apply(res2,1,mean))
plot(Xm$pss_impute, apply(res2,1,mean))



# #What if pss is "removed"?
# model_lbcfp <- list()
# model_bcfp <- list()
# 
# #Predict everyone with the same level of pss
# Xmp <- Xm %>% as.data.frame() %>% mutate(pss_impute=mean(pss_impute)) %>% as.matrix()
# Xcp <- Xc %>% as.data.frame() %>% mutate(pss_impute=mean(pss_impute)) %>% as.matrix()
# 
# for(i in 1:rep){
#   set.seed(1+i)
#   lbcf_default <- lbart:::bcf_new(y, #response
#                                   z = z, z_est = rep(1,npred), #Z
#                                   x_control = Xc, #matrix of covariates
#                                   x_control_est = Xcp,
#                                   x_moderate = Xm,
#                                   x_moderate_est = Xmp,
#                                   pihat = rep(1,npred), #in case of causal inference
#                                   pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
#                                   include_pi = "none", #can be "none", "control", "moderate", or "both"
#                                   linear = "moderate", #can be "none", "control", "moderate", or "both"
#                                   base_control = 0.5, base_moderate = 0.25, #hyperparameters
#                                   power_control = 3, power_moderate = 2, #hyperparameters
#                                   nburn = burn, nthin = th, nsim = post, #draws
#                                   ntree_control = 50, #control trees
#                                   ntree_moderate = 25, #moderate trees (treatment effect)
#                                   dart = F, #Linero's Prior
#                                   save_trees = T,
#                                   use_muscale = T,
#                                   use_tauscale = F)
#   model_lbcfp[[i]] <- lbcf_default
# 
# 
#   set.seed(1+i)
#   bcf_default <- lbart:::bcf_new(y, #response
#                                  z = z, z_est = rep(1,npred), #Z
#                                  x_control = Xc, #matrix of covariates
#                                  x_control_est = Xcp,
#                                  x_moderate = Xm,
#                                  x_moderate_est = Xmp,
#                                  pihat = rep(1,npred), #in case of causal inference
#                                  pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
#                                  include_pi = "none", #can be "none", "control", "moderate", or "both"
#                                  linear = "none", #can be "none", "control", "moderate", or "both"
#                                  base_control = 0.5, base_moderate = 0.5, #hyperparameters
#                                  power_control = 2, power_moderate = 3, #hyperparameters
#                                  nburn = burn, nthin = th, nsim = post, #draws
#                                  ntree_control = 50, #control trees
#                                  ntree_moderate = 50, #moderate trees (treatment effect)
#                                  dart = F, #Linero's Prior
#                                  save_trees = T,
#                                  use_muscale = F,
#                                  use_tauscale = F)
#   model_bcfp[[i]] <- bcf_default
# 
# }
# 
# 
# save(model_lbcfp, file="lbcf_study6p.RData")
# save(model_bcfp, file="bcf_study6p.RData")

# load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/lbcf_study6p.RData")
# load("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/bcf_study6p.RData")
# 
# {bcf_taup <- NULL
#   lbcf_taup <- NULL}
# 
# for(i in 1:rep){
#   bcf_taup <- rbind(bcf_taup, model_bcfp[[i]]$b_est_post)
#   lbcf_taup <- rbind(lbcf_taup, model_lbcfp[[i]]$b_est_post)
# }
# 
# 
# ad_sum3 <- additive_summary(
#   rowMeans(t(lbcf_taup))~s(stressmindset_baseline)+
#     s(fixedmindset_baseline)+sex+s(pss_impute),
#   fhatSamples=t(lbcf_taup),
#   fhat = rowMeans(t(lbcf_taup)),
#   df = study6_clean,
#   alpha = 0.05,
#   fast = TRUE,
#   quants = seq(0, 1, by = 0.005),
#   grid_size = 100,
#   verbose = FALSE,
#   return_samples = TRUE,
#   meta = NA
# )
# 
# ad_sum4 <- additive_summary(
#   rowMeans(t(bcf_taup))~s(stressmindset_baseline)+
#     s(fixedmindset_baseline)+sex+s(pss_impute),
#   fhatSamples=t(bcf_taup),
#   fhat = rowMeans(t(bcf_taup)),
#   df = study6_clean,
#   alpha = 0.05,
#   fast = TRUE,
#   quants = seq(0, 1, by = 0.005),
#   grid_size = 100,
#   verbose = FALSE,
#   return_samples = TRUE,
#   meta = NA
# )
# 
# 
# 
# 
# hist(ad_sum$summaryRsq, ylim=c(0,10000), breaks = 70)
# hist(ad_sum2$summaryRsq, ylim=c(0,10000), breaks = 70)
# 
# additive_summary_plot(ad_sum) + ylim(c(-8,8))
# additive_summary_plot(ad_sum2) + ylim(c(-8,8))
# 
# 
# res <- (t(bcf_taup)-ad_sum$fittedValues)
# 
# Xm <- Xm %>% as.data.frame()
# plot(Xm$stressmindset_baseline, apply(res,1,mean))
# plot(Xm$fixedmindset_baseline, apply(res,1,mean))
# plot(Xm$sex, apply(res,1,mean))
# 
# res2 <- (t(lbcf_taup)-ad_sum$fittedValues)
# plot(Xm$stressmindset_baseline, apply(res2,1,mean))
# plot(Xm$fixedmindset_baseline, apply(res2,1,mean))
# plot(Xm$sex, apply(res2,1,mean))




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
cuts_fixed <- middle(study6_clean$fixedmindset_baseline)
cuts_stress <- middle(study6_clean$stressmindset_baseline)

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
dif_vec <- search_dual(y = apply(lbcf_tau,2,mean), 
                       data = study6_clean, grid_idx = grid_index, grid1 = grid_x1, grid2 = grid_x2, 
                       prop = 0.20)
indexes_max <- dif_vec==max(dif_vec)
fixed_cut = grid_x1[as.numeric((grid_index[indexes_max,])[,1]), , drop = F]
stress_cut = grid_x2[as.numeric((grid_index[indexes_max,])[,2]), , drop = F]

cuts_list <- list()
for(i in 1:sum(indexes_max)){
  cuts_list[[i]] <- data.frame(fixed = fixed_cut[i,], stress = stress_cut[i,])
}

cuts <- cuts_list[[1]]

doublegood <- (study6_clean$fixedmindset_baseline<(cuts[1,1]))&(study6_clean$stressmindset_baseline<(cuts[1,2]))
doublebad <- (study6_clean$fixedmindset_baseline>=(cuts[2,1]))&(study6_clean$stressmindset_baseline>=(cuts[2,2]))

#########################
#########################

intr <- interaction(doublegood,doublebad)

library(possum)
lsubgroup_ates = subgroup_average_posterior(lbcf_tau, intr)
subgroup_ates = subgroup_average_posterior(bcf_tau, intr)
boxplot(lsubgroup_ates)
boxplot(subgroup_ates)
difl <- lsubgroup_ates$TRUE.FALSE-lsubgroup_ates$FALSE.TRUE
difbcf <- subgroup_ates$TRUE.FALSE-subgroup_ates$FALSE.TRUE

#Doublebad, both mindsets are high
hist(lsubgroup_ates$FALSE.TRUE, breaks = 70)
hist(subgroup_ates$FALSE.TRUE, breaks = 70, add = T, col = "red")

#Doublegood, both mindsets are low
hist(lsubgroup_ates$TRUE.FALSE, breaks = 70)
hist(subgroup_ates$TRUE.FALSE, breaks = 70, add = T, col = "red")


cor(Xm$pss_impute, Xm$fixedmindset_baseline)
plot(Xm$pss_impute, Xm$fixedmindset_baseline)



#pdf("anx_priormind_diff.pdf")
plot(density(difbcf), lwd = 2, main = "", col = "blue")
lines(density(difl), col = "red", lwd = 2)
legend("topleft", legend = c("lbcf", "bcf"), fill = c("red", "blue"))
#dev.off()

mean(difbcf>0)
mean(difl>0)

quantile(difbcf, c(0.1, 0.9))
quantile(difl, c(0.1, 0.9))



testando <- df_anx %>% arrange(rmse_test)

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
#pdf("rmse_anx_full.pdf")
par(mfrow=c(1,1))
plot(testando$rmse_test, type = "l", 
     ylim = c(1,5), xlab = "", ylab = "RMSE")
points(testando$rmse_train, pch = 19, cex = 0.1, col = "purple")
points(which(testando$rmse_test==bcfdefault$rmse_test),bcfdefault$rmse_test, 
       col = "red", pch=19)
points(which(testando$rmse_test==bcfdefault$rmse_test),bcfdefault$rmse_train, 
       col = "red", pch=19)
points(which(testando$rmse_test==lbcfdefault$rmse_test),lbcfdefault$rmse_test, 
       col = "blue", pch=19)
points(which(testando$rmse_test==lbcfdefault$rmse_test),lbcfdefault$rmse_train, 
       col = "blue", pch=19)
legend("topleft", legend = c("lbcf-default", "bcf-default"), fill = c("blue", "red"))
#dev.off()

# points((1:length(testando$rmse_train))[testando$linear=="none"],
#        testando$rmse_train[testando$linear=="none"],
#        pch = 19, cex = 0.2, col = "red")
# points((1:length(testando$rmse_train))[testando$linear=="moderate"],
#        testando$rmse_train[testando$linear=="moderate"],
#        pch = 19, cex = 0.2, col = "purple")
# points((1:length(testando$rmse_train))[testando$linear=="control"],
#        testando$rmse_train[testando$linear=="control"],
#        pch = 19, cex = 0.2, col = "blue")
# points((1:length(testando$rmse_train))[testando$linear=="both"],
#        testando$rmse_train[testando$linear=="both"],
#        pch = 19, cex = 0.2, col = "orange")


#pdf("rmse_speech_full.pdf")
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
#dev.off()



##

library(rpart)

# Example dataset: use built-in 'iris' dataset
# We will predict 'Species' using only 'Petal.Length'
df_tr <- data.frame(x=study6_clean$fixedmindset_baseline, y=apply(lbcf_tau,2,mean))
tree_model <- rpart(y ~ x, data = df_tr)
plot(tree_model)



groups1 <- ifelse(study6_clean$fixedmindset_baseline>3.166667, 1, 0)

library(possum)
lsubgroup_ates = subgroup_average_posterior(lbcf_tau, groups1)
subgroup_ates = subgroup_average_posterior(bcf_tau, groups1)
boxplot(lsubgroup_ates)
difl <- lsubgroup_ates$X0-lsubgroup_ates$X1
difbcf <- subgroup_ates$X0-subgroup_ates$X1

hist(lsubgroup_ates$X1, breaks = 70)
hist(subgroup_ates$X1, breaks = 70)
hist(lsubgroup_ates$X0, breaks = 70)
hist(subgroup_ates$X0, breaks = 70)

plot(density(difbcf))
lines(density(difl), col = "red")

mean(difbcf>0)
mean(difl>0)

quantile(difbcf, c(0.1, 0.9))
quantile(difl, c(0.1, 0.9))

