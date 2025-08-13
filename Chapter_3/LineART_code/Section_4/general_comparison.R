library(magrittr)
library(ggplot2)
library(purrr)
library(dplyr)
setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/complete")
current <- 1:41

#MSE
#My default
alpha <- 0.5
beta <- 1
nt <- 50
alldata <- NULL

for(j in current){
  tbl_lbart <- read.csv(paste0("rmse_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",j,"_total.csv"), sep = ";", header = F)[-1,2] %>% as.numeric()
  tbl_bart <- read.csv(paste0("rmse_cv_bart__a_0.95__b_2__nt_200__i_",j,"_total.csv"), sep = ";", header = F)[-1,3] %>% as.numeric()
  tbl_soft <- read.csv(paste0("rmse_cv_soft__s_off__i_",j,"_full.csv"), sep = ";", header = F)[-1,4] %>% as.numeric()
  
  result_table <- NULL
  for(i in 1:length(tbl_lbart)){
    row_i <- t(tbl_lbart)[i]
    row_i2 <- t(tbl_bart)[i]
    row_i3 <- t(tbl_soft)[i]
    
    alldata <- rbind(alldata, data.frame(bart = row_i2^2, lbart = row_i^2, soft = row_i3^2,
                                           a=alpha, b=beta, nt=nt,
                                           dataset=j, rep = i))
    
  }
  
}

load("C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/complete/df_cv.RData")





df_full <- df_cv %>% mutate(lbartcv=rmse^2) %>% 
  select(lbartcv, dataset, rep) %>% 
  left_join(alldata, by = c("dataset", "rep"))

#Let's do RRMSE now
d1 <- cbind(sqrt(df_full$lbart), sqrt(df_full$bart), sqrt(df_full$soft), sqrt(df_full$lbartcv)) %>% apply(1, min)
d1 <- cbind(sqrt(df_full$lbart), sqrt(df_full$bart), sqrt(df_full$soft)) %>% apply(1, min)
d2 <- data.frame(RRMSE=sqrt(df_full$lbart)/d1, name="LineART")
d3 <- data.frame(RRMSE=sqrt(df_full$bart)/d1, name="BART")
d4 <- data.frame(RRMSE=sqrt(df_full$soft)/d1, name="SoftBART")
d5 <- data.frame(RRMSE=sqrt(df_full$lbartcv)/d1, name="LineART-CV")

dfn5 <- rbind(d2,d3,d4,d5)
dfn5 <- rbind(d2,d3,d4)

#Full graphs
ggplot(dfn5, aes(name, RRMSE)) +
  geom_boxplot() +
  xlab(NULL)

#Closer
#setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/nipsdata")
#pdf("RRMSE.pdf")
ggplot(dfn5, aes(name, RRMSE)) +
  geom_boxplot() +
  xlab(NULL)  +
  coord_cartesian(ylim=c(1,1.3))
#dev.off()

# idx <- which(dfn5$RRMSE<1)-820*3
# 
# idx2 <- c(idx, idx+820)
# dfn6 <- dfn5[idx2,]
# 
# ggplot(dfn6, aes(name, RRMSE)) +
#   geom_boxplot() +
#   xlab(NULL)  +
#   coord_cartesian(ylim=c(1,1.3))


#How many times was the best one
dfn5 %>% 
  filter(RRMSE==1) %>% 
  group_by(name) %>% 
  summarise(sum=sum(RRMSE))

#Quantiles
dfn5 %>% 
  group_by(name) %>% 
  summarise(q1=quantile(RRMSE, 0.25), 
            q2=quantile(RRMSE, 0.5), 
            q3=quantile(RRMSE, 0.75))

#% of times larger than 1.3
dfn5 %>% 
  group_by(name) %>% 
  filter(RRMSE>1.3) %>% 
  summarise(n=n()/820)


#clean it up a little bit
df_full <- df_full %>% select(-a, -b, -nt)

graphs <- list()
for(i in 1:41){
  df_aux <- df_full %>% filter(dataset == i) %>% select(-dataset, -rep)
  df_aux1 <- data.frame(MSE = df_aux$lbart, name="LineART")
  df_aux2 <- data.frame(MSE = df_aux$lbartcv, name="LineART-CV")
  df_aux3 <- data.frame(MSE = df_aux$bart, name="BART")
  df_aux4 <- data.frame(MSE = df_aux$soft, name="SoftBART")
  df_aux <- rbind(df_aux1, df_aux2, df_aux3, df_aux4)
  graphs[[i]] <- ggplot(df_aux, aes(MSE, name)) +
    ggtitle(paste0("MSE - Dataset ",i)) +
    xlab(NULL) + 
    ylab(NULL) +
    geom_boxplot()
}

for(i in 1:41){
  print(graphs[[i]])
}


df_aux1 <- data.frame(MSE = df_full$lbart, name="LineART", dataset = df_full$dataset, rep = df_full$rep)
df_aux2 <- data.frame(MSE = df_full$lbartcv, name="LineART-CV", dataset = df_full$dataset, rep = df_full$rep)
df_aux3 <- data.frame(MSE = df_full$bart, name="BART", dataset = df_full$dataset, rep = df_full$rep)
df_aux4 <- data.frame(MSE = df_full$soft, name="SoftBART", dataset = df_full$dataset, rep = df_full$rep)
df_aux <- rbind(df_aux1, df_aux2, df_aux3, df_aux4)


datasets <- list.files(path="C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/nipsdata", 
                       pattern = "x.txt", full.names=F, recursive=TRUE)
data_names <- unlist(strsplit(datasets, split='/', fixed=TRUE))[as.logical(1:82%%2)]

for(i in 1:41){
  df_aux$dataset[df_aux$dataset==i] <- data_names[i]
}


pdf("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/LineART_code/all_datasets.pdf")
ggplot(df_aux, aes(MSE, name)) +
  xlab("MSE") + 
  ylab("") +
  geom_boxplot() +
  facet_wrap(~dataset, scales = "free_y", nrow = 6) +
  coord_flip() +
  scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()



tbl_aux <- df_aux %>% group_by(dataset, name) %>% 
  summarise(avgmse=mean(MSE))
tbl_aux2 <- df_aux %>% group_by(dataset, name) %>% 
  summarise(sdmse=sd(MSE))
tbl_aux3 <- tbl_aux
tbl_aux3$avgmse <- paste(round(tbl_aux$avgmse,5), "(", round(tbl_aux2$sdmse,5), ")", sep = "")

library(tidyverse)
tbl_avg <- tbl_aux3 %>%
  pivot_wider(names_from = name, values_from = avgmse)


write.csv(tbl_avg, "C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/LineART_code/datasets_avgmse.csv", row.names = F)


