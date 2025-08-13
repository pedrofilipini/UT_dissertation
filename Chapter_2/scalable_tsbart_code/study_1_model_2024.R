{library(dplyr)
library(ggplot2)
library(latex2exp)

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

#Auxiliary functions
source("C:/Users/pedro/OneDrive/Área de Trabalho/Important_so_doing_it_again/Nature/Project_Stress/Dataset_1/Code/gp_approx_fun.R")

#Add a baseline variable that is the average of first 5 min
use_mean_baseline = TRUE

#The original dataset
data_1 <- read.csv("C:/Users/pedro/OneDrive/Área de Trabalho/Important_so_doing_it_again/Nature/Project_Stress/New Dataset/bcfstudy1_long.csv")

#Adding a variable that is the observed mean for "base" for each individual
mean_base <- data_1[data_1$epoch=="base",] %>%
  group_by(id) %>%
  summarise(mean_tpr = mean(tpr, na.rm = T),
            mean_bps = mean(bps, na.rm = T))

#Adding new columns in the data
data_1 <- left_join(data_1, mean_base, by ="id")
if(use_mean_baseline) data_1  = filter(data_1, epoch!="base")
if(use_mean_baseline) data_1  = data_1[!is.na(data_1$mean_tpr),]

#Remove all the NAs and technical issues
data_1$tpr[data_1$tpr<650] = NA
data_1$tpr[data_1$tpr>5000] = NA
data_1_clean <- data_1[!is.na(data_1$tpr),]
data_1_clean <- data_1_clean %>% mutate(reactivity_tpr = tpr - mean_tpr)


#Separating epochs depending on each case
epoch_train <- as.numeric(factor(data_1_clean$epoch,
                                 ordered = T,
                                 levels = c("prep","speech","math","recov")))


# Person-level dataframe for everyone, useful for prediction later
person_df <- data_1_clean %>% group_by(id) %>% filter(row_number()==1) %>%group_by()
epochs = sort(unique(data_1_clean$epoch))


#Functions to get the covariates
get_controls = function(x) {

  ret = x %>% dplyr::select(fixedmindset_baseline,
                            stressmindset_baseline, bothmindsets,
                            selfesteem_baseline, sex,
                            mean_tpr
                            ) %>%
    {if (use_mean_baseline) . else dplyr::select(., -mean_tpr) } %>%
    dbarts::makeModelMatrixFromDataFrame(drop = FALSE)

}

get_moderators = function(x) {
  x %>% dplyr::select(bothmindsets,
                      stressmindset_baseline,
                      fixedmindset_baseline,
                      sex)%>%
    dbarts::makeModelMatrixFromDataFrame(drop = FALSE)
}


#The treatment and the t variable, scale it between -1 and 1
t_train <- (data_1_clean$time-(max(data_1_clean$time)+min(data_1_clean$time))/2)/((max(data_1_clean$time)-min(data_1_clean$time))/2)
z_train <- data_1_clean$treatment


#Now for the covariates
x_train <- get_controls(data_1_clean)
x_train_person <- get_controls(person_df)

x_mod <-   get_moderators(data_1_clean)
x_mod_person <- get_moderators(person_df)


#Random effect matrices, same structure as Yeager's syntax (random intercepts)
library(lme4)
fit = lmer(reactivity_tpr~(1|id), data=data_1_clean)
randeff_design = getME(fit, "Z")
Gp = getME(fit, "Gp")
Q = matrix(0, nrow=ncol(randeff_design), ncol = length(Gp)-1)
for(j in 2:length(Gp)) Q[(1+Gp[j-1]):Gp[j],j-1]=1

#Model - TPR
set.seed(34764)
}
fit_2_long <- multibart:::bcf_core(data_1_clean$reactivity_tpr,
                                   z_train = z_train,
                                   z_out = rep(1, times = dim(x_train_person)[1]),
                                   tvar_con = t_train,
                                   tvar_mod = t_train,
                                   x_control = x_train,
                                   x_moderate = x_mod,
                                   x_control_out = x_train_person,
                                   x_moderate_out = x_mod_person,
                                   tvar_con_out = rep(0, times = dim(x_train_person)[1]),
                                   tvar_mod_out = rep(0, times = dim(x_train_person)[1]),
                                   randeff_design = as.matrix(randeff_design),
                                   randeff_variance_component_design = Q,
                                   randeff_scales = 300,
                                   vanilla = FALSE, dart = F,
                                   ntree_moderate = 100,
                                   ntree_control = 250,
                                   nburn = 15000, nsim = 3000, nthin = 5)


mean((data_1_clean$reactivity_tpr[(data_1_clean$time<9 & data_1_clean$treatment==1)]))
mean((data_1_clean$reactivity_tpr[(data_1_clean$time<9 & data_1_clean$treatment==0)]))

#Quick check just to see if everything is more or less in line before deep analysis
plot.ts(fit_2_long$eta_con)
plot.ts(fit_2_long$l_con_post)
plot.ts(fit_2_long$sigma)

plot.ts(fit_2_long$eta_mod)
plot.ts(fit_2_long$l_mod_post)
plot.ts(fit_2_long$sigma)

acf(fit_2_long$l_con_post)
acf(fit_2_long$l_mod_post)

plot(density(fit_2_long$l_con_post))
plot(density(fit_2_long$l_mod_post))


setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Basis")
save(fit_2_long, data_1_clean, person_df,
     file=paste0("study_1_long_randeff_reactivity_tau.Rdata"))

#}}
