# Imai, Kosuke, and David A. van Dyk. (2004).
#``Causal Inference With General Treatment Regimes: Generalizing the Propensity Score,
#'' Journal of the American Statistical Association, Vol. 99, No. 467 (September), pp. 854-866.
# 
#https://imai.fas.harvard.edu/research/pscore.html

#Libraries
library(dplyr)
#Smoking data
###
data<-read.csv("C:/Users/pedro/OneDrive/Documentos/Code/Conbart/smoke/pscore/nmesdata.txt", header=T)
#Must have some medical expenditure
data<-subset(data, TOTALEXP>0)
#Must be a smoker 
data<-subset(data, packyears>0)
plot(data$packyears,data$LASTAGE)
#But not too old to avoid overlap
data<-subset(data, packyears<=70)

#Complete case analysis
data<-subset(data, select=c(packyears, #treatment 
                            AGESMOKE, #Numeric, when started smoking
                            LASTAGE, #Numeric, when stopped smoking
                            MALE, RACE3, beltuse, educate, marital, #factors
                            SREGION, POVSTALB, #factors
                            HSQACCWT, #Weight to account for population
                            TOTALEXP)) #Response
#Selected this based on the observed packyears vs lastage
plot(data$packyears,data$LASTAGE)
#Removing NAs
data<-na.omit(data)
#Factor variables
data$RACE3<-factor(data$RACE3, ordered=F)
data$marital<-factor(data$marital, ordered=F)
data$SREGION<-factor(data$SREGION, ordered=F)
data$educate<-factor(data$educate, ordered=F)
data$beltuse<-factor(data$beltuse, ordered=F)
data$POVSTALB<-factor(data$POVSTALB, ordered=F)

#Keeping without scale for plots later
data_og <- data

data$AGESMOKE <- scale(data$AGESMOKE)
data$LASTAGE <- scale(data$LASTAGE)
data$HSQACCWT <- scale(data$HSQACCWT)

# View(data)

#Data for the models
x_con <- data %>% select(AGESMOKE, LASTAGE, MALE, RACE3, beltuse,
                         educate, marital, SREGION, POVSTALB, HSQACCWT) %>%
  dbarts::makeModelMatrixFromDataFrame(drop = FALSE)
x_mod <- data %>% select(AGESMOKE, LASTAGE, MALE, RACE3, beltuse,
                         educate, marital, SREGION, POVSTALB, HSQACCWT) %>%
  dbarts::makeModelMatrixFromDataFrame(drop = FALSE)

#The treatment
z <- log(data$packyears)
npred <- length(z)

#The response
y <- log(data$TOTALEXP)


burn <- 50000
post <- 5000
th <- 5
set.seed(23454)


#I can use lineart to estimate the pscore, then apply
lbcf_gps <- lbart:::bcf_new(z, #response
                            z = rep(1,npred), #placeholder
                            z_est = rep(1,npred), #placeholder
                            x_control = x_con, #matrix of covariates
                            x_control_est = x_con,
                            x_moderate = x_mod,
                            x_moderate_est = x_mod,
                            pihat = rep(1,npred), #in case of causal inference
                            pihat_est = rep(1,npred), #for predictions (does not matter by now, not implemented)
                            include_pi = "none", #can be "none", "control", "moderate", or "both"
                            linear = "both", #can be "none", "control", "moderate", or "both"
                            base_control = 0.5, base_moderate = 0.25, #hyperparameters
                            power_control = 1, power_moderate = 3, #hyperparameters
                            nburn = burn, nthin = th, nsim = post, #draws
                            ntree_control = 50, #control trees
                            ntree_moderate = 0, #moderate trees (treatment effect)
                            dart = F, #Linero's Prior
                            save_trees = T,
                            use_muscale = T,
                            use_tauscale = T)

#Posterior mean
zhat <- apply(lbcf_gps$mu_post, 2, mean)

#Sigma
sigma_hat <- mean(lbcf_gps$sigma)

#GPS (density assuming normality)
gps_hat <- numeric(0)
for(i in 1:length(z)){
  #gps_hat[i] <- dnorm(z[i], mean = zhat[i], sd = sigma_hat)
  gps_hat[i] <- zhat[i] #Not using the GPS now
}

rm(lbcf_gps)
gc()

#Using GPS now as a variable in conbart
burn <- 50000
post <- 5000
th <- 5

#Chaging z to be between -1 and 1
z2 <- (z-(max(z)+min(z))/2)/((max(z)-min(z))/2)


#Change expected number of crossings to 1
set.seed(32859123)
source("C:/Users/pedro/OneDrive/Área de Trabalho/Important_so_doing_it_again/Nature/Project_Stress/Dataset_1/Code/gp_approx_fun.R")
fit <- conbart::bcf_core(as.matrix(y), #response
                         pihat = as.matrix(gps_hat), #estimated gps
                         pihat_out = as.matrix(gps_hat), #estimated gps
                         include_pi = "control", #include gps as covariate in con
                         z_train = rep(1, times = npred), #placeholder
                         z_out = rep(1, times = npred), #placeholder
                         tvar_mod = as.matrix(z2), #treatment variable
                         x_control = x_con, #control
                         x_moderate = x_mod, #mod 
                         x_control_out = x_con, #For prediction
                         x_moderate_out = x_mod, #For prediction
                         tvar_mod_out = as.matrix(z2), #treatment variable
                         vanilla = T, #first part is a bart, second is a stsbart
                         dart = F, #No variable selection
                         enc_mod = 3,
                         ntree_moderate = 50, #Same default as bcf
                         ntree_control = 200, #Same default as bcf
                         nburn = burn, nsim = post, nthin = th) #MCMC sample



time_grid = (log(c(2,4,8,16,32,64))-(max(z)+min(z))/2)/((max(z)-min(z))/2)
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
  
  #Then calculate the Average Treatment Effect curve posterior draws
  #  ate = Reduce("+", ite)/length(ite)
}



auxmat <- matrix(0,ncol=post,nrow=npred)
auxmat2 <- matrix(0,ncol=post,nrow=npred)
auxmat3 <- matrix(0,ncol=post,nrow=npred)
auxmat4 <- matrix(0,ncol=post,nrow=npred)
auxmat5 <- matrix(0,ncol=post,nrow=npred)
#auxmat6 <- matrix(0,ncol=post,nrow=npred)
#auxmat7 <- matrix(0,ncol=post,nrow=npred)

for(j in 1:npred){
  auxmat[j,] <- ite[[j]][2,]-ite[[j]][1,]
  auxmat2[j,] <- ite[[j]][3,]-ite[[j]][2,]
  auxmat3[j,] <- ite[[j]][4,]-ite[[j]][3,]
  auxmat4[j,] <- ite[[j]][5,]-ite[[j]][4,]
  auxmat5[j,] <- ite[[j]][6,]-ite[[j]][5,]
}

dim(auxmat)
hist(exp(apply(auxmat, 1, mean)))

library(ggplot2)
library(dplyr)

ate_t1_t2 <- apply(auxmat, 2, mean) %>% exp()
ate_t2_t3 <- apply(auxmat2, 2, mean) %>% exp()
ate_t3_t4 <- apply(auxmat3, 2, mean) %>% exp()
ate_t4_t5 <- apply(auxmat4, 2, mean) %>% exp()
ate_t5_t6 <- apply(auxmat5, 2, mean) %>% exp()


posterior_df <- bind_rows(
  data.frame(ATE = ate_t1_t2, Packyears = "02 vs 04"),
  data.frame(ATE = ate_t2_t3, Packyears = "04 vs 08"),
  data.frame(ATE = ate_t3_t4, Packyears = "08 vs 16"),
  data.frame(ATE = ate_t4_t5, Packyears = "16 vs 32"),
  data.frame(ATE = ate_t5_t6, Packyears = "32 vs 64")
)

#Plot
pdf("smoke_cate1.pdf", width = 7, height = 5)
ggplot(posterior_df, aes(x = ATE, fill = Packyears)) +
  geom_density(alpha = 0.5) +
  labs(title = "",
       x = "Multiplicative factor for total medical expenses",
       y = "Density") +
  theme_minimal()
dev.off()

auxmat8 <- matrix(0,ncol=post,nrow=npred)
for(j in 1:npred){
  auxmat8[j,] <- ite[[j]][6,]-ite[[j]][3,]
}


ate_t3_t6 <- apply(auxmat8, 2, mean) %>% exp()
posterior_df2 <- data.frame(ATE = ate_t3_t6, Packyears = "08 vs 64")

summary_stats <- quantile(ate_t3_t6, probs = c(0.025, 0.5, 0.975))
pdf("smoke_cate2.pdf")
ggplot(posterior_df2, aes(x = ATE)) +
  geom_density(fill = "#69b3a2", alpha = 0.6, color = "darkgreen") +
  geom_vline(xintercept = summary_stats[2], color = "black", linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = summary_stats[c(1, 3)], color = "darkred", linetype = "dashed", size = 1.5) +
  labs(
    title = "",
    x = "Multiplicative factor for total medical expenses",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)
dev.off()

#
pdf("smoke_itevsage.pdf")
ite_t3_t6 <- apply(auxmat8, 1, mean) %>% exp()
plot(data_og$LASTAGE, ite_t3_t6, cex = 0.5, pch = 19,
     xlab = "Age", ylab = "Estimated treatment effect")
dev.off()

library(rpart)
#Creating the dummy variables
data_tree <- data_og %>% select(AGESMOKE, LASTAGE, MALE, RACE3, beltuse,
                         educate, marital, SREGION, POVSTALB, HSQACCWT) %>%
  dbarts::makeModelMatrixFromDataFrame(drop = FALSE) %>% as.data.frame()


ite_tree <- rpart(ite_t3_t6 ~ ., data = data_tree,
                  control = rpart.control(maxdepth = 2))

library(partykit)
ite_tree_party <- as.party(ite_tree)
pdf("smoke_agetree.pdf", width = 12, height = 5)
plot(ite_tree_party,ip_args = list(id=F), tp_args = list(id=F))
dev.off()



g1 <- auxmat8[data_og$LASTAGE<48.5 & data_og$MALE==1,]
g2 <- auxmat8[data_og$LASTAGE>=48.5 & data_og$MALE==0,]

dim(g1)
dim(g2)

dif_g <- (apply(exp(g1), 2, mean)-apply(exp(g2), 2, mean))

#Subgroup difference
posterior_df3 <- data.frame(ATE = dif_g, Comparison = "Young men vs older women")

summary_stats_2 <- quantile(dif_g, probs = c(0.025, 0.5, 0.975))
pdf("smoke_catediff.pdf")
ggplot(posterior_df3, aes(x = ATE)) +
  geom_density(fill = "#69b3a2", alpha = 0.6, color = "darkgreen") +
  #geom_vline(xintercept = summary_stats_2[2], color = "black", linetype = "dashed", size = 1.5) +
  #geom_vline(xintercept = summary_stats_2[c(1, 3)], color = "darkred", linetype = "dashed", size = 1.5) +
  labs(
    title = "",
    x = "Subgroup difference",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)
dev.off()

mean(posterior_df3>0)



lmod_df <- data.frame(lmod=fit$l_mod_post)
summary_stats_3 <- quantile(lmod_df$lmod, probs = c(0.025, 0.5, 0.975))
pdf("smoke_ls.pdf")
ggplot(lmod_df, aes(x = lmod)) +
  geom_density(fill = "blue", alpha = 0.6, color = "blue") +
  geom_vline(xintercept = summary_stats_3[2], color = "black", linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = summary_stats_3[c(1, 3)], color = "darkred", linetype = "dashed", size = 1.5) +
  labs(
    title = "",
    x = "Length-scale posterior sample",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)
dev.off()




