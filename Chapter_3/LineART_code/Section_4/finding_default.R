library(magrittr)
library(ggplot2)
library(purrr)
library(dplyr)
setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/complete")
current <- 1:41
settings <- list()
counts <- list()
aux <- 0
for(alpha in c(0.5, 0.8, 0.95)){
  for(beta in c(1, 2, 3)){
    for(nt in c(10, 50, 100, 150, 200)){
      aux <- aux+1
      settings[[aux]] <- c(alpha, beta, nt)
      
      results <- NULL
      #RMSE
      lpairt <- numeric(0)
      lpv <- numeric(0)
      lspairt <- numeric(0)
      lspv <- numeric(0)
      pairt <- numeric(0)
      pvaluesp <- numeric(0)
      softp <- numeric(0)
      softp2 <- numeric(0)
      bmeans <- numeric(0)
      
      #alpha <- 0.95
      #beta <- 2
      #nt <- 200
      
      for(j in current){
        tbl_aux <- read.csv(paste0("rmse_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",j,"_total.csv"), sep = ";", header = F)[-1,2] %>% as.numeric()
        tb_i <- read.csv(paste0("rmse_cv_bart__a_0.95__b_2__nt_200__i_",j,"_total.csv"), sep = ";", header = F)[-1,]
        tb_i3 <- read.csv(paste0("rmse_cv_soft__s_off__i_",j,"_full.csv"), sep = ";", header = F)[-1,]
        
        result_table <- NULL
        tab_list <- tb_i[,3]
        tab_list3 <- tb_i3[,4]
        for(i in 1:length(tab_list)){
          row_i <- t(as.numeric(tab_list))[i]
          row_i2 <- t(tbl_aux)[i]
          row_i3 <- t(as.numeric(tab_list3))[i]
          
          if(aux==1 & j == current[1] & i ==1){
            alldata <- data.frame(lbart = row_i2^2, bart = row_i^2, soft = row_i3^2,
                                  a=alpha, b=beta, nt=nt,
                                  dataset=j, rep = i)
          }else{
            alldata <- rbind(alldata, data.frame(lbart = row_i2^2, bart = row_i^2, soft = row_i3^2,
                                                 a=alpha, b=beta, nt=nt,
                                                 dataset=j, rep = i))
          }
          
          if((i==1)){
            result_table <- c(row_i^2, row_i2^2, row_i3^2)
          }else{
            result_table <- rbind(result_table, c(row_i^2, row_i2^2, row_i3^2))
          }
        }
        
        if(j == current[1]){
          lpairt[j] <- t.test(log(result_table[,2]/result_table[,1]))$estimate
          lpv[j] <- t.test(log(result_table[,2]/result_table[,1]))$p.value
          lspairt[j] <- t.test(log(result_table[,3]/result_table[,1]))$estimate
          lspv[j] <- t.test(log(result_table[,3]/result_table[,1]))$p.value
          pairt[j] <- t.test(result_table[,2], result_table[,1], paired = T)$estimate
          pvaluesp[j] <- t.test(result_table[,2], result_table[,1], paired = T)$p.value
          bmeans[j] <- mean(result_table[,1])
          softp[j] <- t.test(result_table[,3], result_table[,1], paired = T)$estimate/bmeans[j]
          softp2[j] <- t.test(result_table[,3], result_table[,1], paired = T)$p.value
          
          
          pairs2 <- numeric(0)
          pairs2 <- pairt/bmeans
          pairs3 <- ifelse(pvaluesp<0.05, pairs2, NA)
          softp3 <- ifelse(softp2<0.05, softp, NA)
          counts[[aux]] <- c(sum(pairs3<0, na.rm = T), sum(pairs3>0, na.rm = T))
          
          if(aux==1 & j == current[1]){
            maindf <- data.frame(ldiv=log(mean(result_table[,2])/mean(result_table[,1])),lsoft=log(mean(result_table[,3])/mean(result_table[,1])),
                                 ptlbart=pairt[j], ptbart= bmeans[j], div=pairs2[j], pv=pvaluesp[j],
                                 ptsoft=softp[j], pvsoft= softp2[j], a=alpha, b=beta, nt=nt,
                                 dataset=j, a1=paste0("a=",alpha),b1=paste0("b=",beta),nt1=paste0("nt=",nt),
                                 lpairt = lpairt[j], lpv = lpv[j], lspairt=lspairt[j], lspv=lspv[j])
          }else{
            maindf <- rbind(maindf, data.frame(ldiv=log(mean(result_table[,2])/mean(result_table[,1])),lsoft=log(mean(result_table[,3])/mean(result_table[,1])),
                                               ptlbart=pairt[j], ptbart= bmeans[j], div=pairs2[j], pv=pvaluesp[j],
                                               ptsoft=softp[j], pvsoft= softp2[j], a=alpha, b=beta, nt=nt,
                                               dataset=j, a1=paste0("a=",alpha),b1=paste0("b=",beta),nt1=paste0("nt=",nt),
                                               lpairt = lpairt[j], lpv = lpv[j], lspairt=lspairt[j], lspv=lspv[j]))
            
          }
          
          
        }else{
          lpairt[j] <- t.test(log(result_table[,2]/result_table[,1]))$estimate
          lpv[j] <- t.test(log(result_table[,2]/result_table[,1]))$p.value
          lspairt[j] <- t.test(log(result_table[,3]/result_table[,1]))$estimate
          lspv[j] <- t.test(log(result_table[,3]/result_table[,1]))$p.value
          pairt[j] <- t.test(result_table[,2], result_table[,1], paired = T)$estimate
          pvaluesp[j] <- t.test(result_table[,2], result_table[,1], paired = T)$p.value
          bmeans[j] <- mean(result_table[,1])
          softp[j] <- t.test(result_table[,3], result_table[,1], paired = T)$estimate/bmeans[j]
          softp2[j] <- t.test(result_table[,3], result_table[,1], paired = T)$p.value
          
          
          pairs2 <- numeric(0)
          pairs2 <- pairt/bmeans
          pairs3 <- ifelse(pvaluesp<0.05, pairs2, NA)
          softp3 <- ifelse(softp2<0.05, softp, NA)
          counts[[aux]] <- c(sum(pairs3<0, na.rm = T), sum(pairs3>0, na.rm = T))
          
          if(aux==1 & j == current[1]){
            maindf <- data.frame(ldiv=log(mean(result_table[,2])/mean(result_table[,1])),lsoft=log(mean(result_table[,3])/mean(result_table[,1])),
                                 ptlbart=pairt[j], ptbart= bmeans[j], div=pairs2[j], pv=pvaluesp[j],
                                 ptsoft=softp[j], pvsoft= softp2[j], a=alpha, b=beta, nt=nt,
                                 dataset=j, a1=paste0("a=",alpha),b1=paste0("b=",beta),nt1=paste0("nt=",nt),
                                 lpairt = lpairt[j], lpv = lpv[j], lspairt=lspairt[j], lspv=lspv[j])
          }else{
            maindf <- rbind(maindf, data.frame(ldiv=log(mean(result_table[,2])/mean(result_table[,1])),lsoft=log(mean(result_table[,3])/mean(result_table[,1])),
                                               ptlbart=pairt[j], ptbart= bmeans[j], div=pairs2[j], pv=pvaluesp[j],
                                               ptsoft=softp[j], pvsoft= softp2[j], a=alpha, b=beta, nt=nt,
                                               dataset=j, a1=paste0("a=",alpha),b1=paste0("b=",beta),nt1=paste0("nt=",nt),
                                               lpairt = lpairt[j], lpv = lpv[j], lspairt=lspairt[j], lspv=lspv[j]))
            
          }
          
        }
        
      }
    }
  }
}

#Organizing the factors
maindf$a1 = factor(maindf$a1, levels=c('a=0.5','a=0.8','a=0.95'))
maindf$b1 = factor(maindf$b1, levels=c('b=1','b=2','b=3'))
maindf$nt1 = factor(maindf$nt1, levels=c('nt=10','nt=50','nt=100','nt=150', 'nt=200'))

alldata <- alldata %>% mutate(nt1 = paste0("nt=",nt), a1 = paste0("a=",a), b1 = paste0("b=",b))
alldata$a1 = factor(alldata$a1, levels=c('a=0.5','a=0.8','a=0.95'))
alldata$b1 = factor(alldata$b1, levels=c('b=1','b=2','b=3'))
alldata$nt1 = factor(alldata$nt1, levels=c('nt=10','nt=50','nt=100','nt=150','nt=200'))

#Saving graphs from all datasets with all settings
graphs <- list()
graphs2 <- list()
for(i in current){
  dfn3 <- alldata %>% 
    filter(dataset==i) %>% 
    group_by(a1, b1, nt1) %>% 
    summarise(mean = mean(lbart),
              median = median(lbart),
              lb2 = quantile(lbart, probs = 0.1),
              ub2 = quantile(lbart, probs = 0.90),
              lb = quantile(lbart, probs = 0.025),
              ub = quantile(lbart, probs = 0.975),
              meanb = mean(bart),
              medianb = median(bart),
              lb2b = quantile(bart, probs = 0.1),
              ub2b = quantile(bart, probs = 0.90),
              lbb = quantile(bart, probs = 0.025),
              ubb = quantile(bart, probs = 0.975)
    )
  
  default <- dfn3 %>% filter(a1 == "a=0.5", b1 == "b=1", nt1=="nt=50")
  dfn3_aux <- dfn3
  dfn3_aux <- dfn3_aux %>% group_by(a1, b1) %>% summarise()
  dfn3_aux$nt1 <- "Default"
  dfn3_aux$mean <- default$mean
  dfn3_aux$median <- default$median
  dfn3_aux$ub <- default$ub
  dfn3_aux$lb <- default$lb
  dfn3_aux$ub2 <- default$ub2
  dfn3_aux$lb2 <- default$lb2
  
  dfn3 <- rbind(dfn3, dfn3_aux)
  
  dfn3$nt1 = factor(dfn3$nt1, levels=c('nt=10','nt=50','nt=100','nt=150','nt=200','Default'))
  
  
  graphs[[i]] <- ggplot(dfn3, aes(nt1,mean)) +
    ggtitle(paste0("Dataset ",i, "; 2.5%, Mean, 97.5%")) +
    geom_point() + 
    geom_point(size = 0.9, color = "red", aes(nt1, ub)) +
    geom_point(size = 0.9, color = "red", aes(nt1, lb)) +
    geom_abline(slope=0,intercept=dfn3$ubb, color = "green") +
    geom_abline(slope=0,intercept=dfn3$lbb, color = "green") +
    geom_abline(slope=0,intercept=dfn3$meanb, color = "green4") +
    geom_point(size = 0.9, color = "purple", aes(nt1, ub), data = ~subset(., nt1 == "Default")) +
    geom_point(size = 0.9, color = "purple", aes(nt1, lb), data = ~subset(., nt1 == "Default")) +
    geom_point(size = 1.5, color = "purple", aes(nt1, mean), data = ~subset(., nt1 == "Default")) +
    xlab(NULL) + 
    ylab(NULL) +
    facet_grid(b1~a1) +
    theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))
  
  graphs2[[i]] <- ggplot(dfn3, aes(nt1,median)) +
    ggtitle(paste0("Dataset ",i, "; 10%, Median, 90%")) +
    geom_point() + 
    geom_point(size = 0.9, color = "firebrick", aes(nt1, ub2)) +
    geom_point(size = 0.9, color = "firebrick", aes(nt1, lb2)) +
    geom_abline(slope=0,intercept=dfn3$ub2b, color = "green") +
    geom_abline(slope=0,intercept=dfn3$lb2b, color = "green") +
    geom_abline(slope=0,intercept=dfn3$medianb, color = "green4") +
    geom_point(size = 0.9, color = "darkmagenta", aes(nt1, ub2), data = ~subset(., nt1 == "Default")) +
    geom_point(size = 0.9, color = "darkmagenta", aes(nt1, lb2), data = ~subset(., nt1 == "Default")) +
    geom_point(size = 1.5, color = "purple", aes(nt1, median), data = ~subset(., nt1 == "Default")) +
    xlab(NULL) + 
    ylab(NULL) +
    facet_grid(b1~a1) +
    theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))
}

for(i in current){
  par(mfrow=c(1,2))
  print(graphs[[i]])
  print(graphs2[[i]])
}



#The lowest RMSE for every simulation
winner <- list()
for(i in unique(alldata$dataset)){
  for(j in unique(alldata$rep)){
    if(j==1){
      winner[[i]] <- list()
    }
    winner[[i]][[j]] <- alldata %>% 
      filter(dataset==i, rep==j) %>% 
      arrange(lbart) %>% 
      slice(1) %>% 
      select(a, b, nt)
  }
}

#Find the one that wins most datasets
winner2 <- list()
for(i in unique(alldata$dataset)){
  winner2[[i]] <- map_df(winner[[i]], ~as.data.frame(.)) %>% 
    add_count(a, b, nt) %>% 
    filter(n==max(n))
}


#Reorganize in a dataframe
winner3 <- list()
for(i in unique(alldata$dataset)){
  winner3[[i]] <- map_df(winner[[i]], ~as.data.frame(.))
}
winner4 <- map_df(winner3, ~as.data.frame(.))

#This is just to add the dataset and rep columns
winner4_aux <- alldata %>% 
  filter(a==0.5,b==1,nt==50) %>% 
  select(dataset, rep)

winner4 <- cbind(winner4,winner4_aux)
winner5 <- left_join(winner4, alldata, by = c("a", "b", "nt", "dataset", "rep"))


#0.5, 1, 10
#This does not make much sense because 
#I am not considering how bad it was in the other datasets
map_df(winner2, ~as.data.frame(.)) %>% 
  add_count(a, b, nt) %>% 
  filter(nn==max(nn))


#0.5, 1, 50
#So I will count how many times it was better than BART in general
#better than bart 448 times
#0.5463415
bet <- alldata %>% 
  mutate(better=ifelse(lbart<bart, T, F)) %>% 
  group_by(a, b, nt) %>% 
  summarise(sbetter = sum(better)) #%>% 
  filter(sbetter==max(sbetter))

#And how many times it was better than SoftBART
#better than softbart 359 times
#0.4378049
alldata %>% 
  mutate(better=ifelse(lbart<soft, T, F)) %>% 
  group_by(a, b, nt) %>% 
  summarise(sbetter = sum(better)) %>% 
  filter(sbetter==max(sbetter))

#As a measure of comparison, that's how many times SoftBART was better than BART
#better than bart 493 times
#0.6012195
alldata %>% 
  filter(a==0.5, b==1, nt==50) %>% #I am just selecting some feature so I only count soft once
  mutate(better=ifelse(soft<bart, T, F)) %>% 
  group_by(a, b, nt) %>% 
  summarise(sbetter = sum(better)) %>% 
  filter(sbetter==max(sbetter))


#Ok, seems decent, let's take a closer look
df1 <- maindf %>% filter(a==0.5, b==1, nt==50)

#par(mfrow=(c(1,1)))
plot(df1$dataset, df1$lpairt, main = "Linear", ylim = c(-2.5, 1.5))
inda <- which(df1$lpv<0.05 & df1$lpairt<0)
indb <- which(df1$lpv<0.05 & df1$lpairt>0)
points(df1$dataset[inda], df1$lpairt[inda], col = "red", pch = 19)
points(df1$dataset[indb], df1$lpairt[indb], col = "blue", pch = 19)
abline(h=0)
abline(h=.1, col ="purple")
abline(h=-.1, col ="purple")

plot(df1$dataset, df1$lspairt, main = "softBART", ylim = c(-2.5, 1.5))
indsa <- which(df1$lspv<0.05 & df1$lspairt<0)
indsb <- which(df1$lspv<0.05 & df1$lspairt>0)
points(df1$dataset[indsa], df1$lspairt[indsa], col = "red", pch = 19)
points(df1$dataset[indsb], df1$lspairt[indsb], col = "blue", pch = 19)
abline(h=0)
abline(h=.1, col ="purple")
abline(h=-.1, col ="purple")



#Let's do RRMSE now
dfn4 <- alldata %>% 
  filter(a==0.5, b==1, nt==50)

d4 <- cbind(dfn4$lbart, dfn4$bart, dfn4$soft) %>% apply(1, min)
#d4 <- cbind(dfn4$lbart, dfn4$bart, dfn4$soft, winner5$lbart) %>% apply(1, min)

d1 <- data.frame(RRMSE=dfn4$lbart/d4, name="LineART")
d2 <- data.frame(RRMSE=dfn4$bart/d4, name="BART")
d3 <- data.frame(RRMSE=dfn4$soft/d4, name="SoftBART")
#d5 <- data.frame(RRMSE=winner5$lbart/d4, name="LineART-CV")

dfn5 <- rbind(d1,d2,d3)
#dfn5 <- rbind(d1,d2,d3,d5)

#Full graphs
g1 <- ggplot(dfn5, aes(name, RRMSE)) +
  geom_boxplot() +
  xlab(NULL)

#Closer
g2 <- ggplot(dfn5, aes(name, RRMSE)) +
  geom_boxplot() +
  xlab(NULL)  +
  coord_cartesian(ylim=c(1,1.3))

#Percentage over 1.3
dfn5 %>% group_by(name) %>% summarise(sum(RRMSE>1.3)/820)

library(gridExtra)
pdf("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/LineART_code/plots_rrmse_default.pdf", width = 14, height = 7)
grid.arrange(g1, g2, ncol = 2)
dev.off()

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

#% of times larger than 1.5
dfn5 %>% 
  group_by(name) %>% 
  filter(RRMSE>1.5) %>% 
  summarise(n=n()/820)

