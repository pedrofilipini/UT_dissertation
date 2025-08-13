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
#d1 <- data.frame(RRMSE=sqrt(df_full$bart))
d2 <- data.frame(RRMSE=sqrt(df_full$lbart)/d1, name="LineART")
d3 <- data.frame(RRMSE=sqrt(df_full$bart)/d1, name="BART")
d4 <- data.frame(RRMSE=sqrt(df_full$soft)/d1, name="SoftBART")
d5 <- data.frame(RRMSE=sqrt(df_full$lbartcv)/d1, name="LineART-CV")

dfn5 <- rbind(d2,d3,d4,d5)
dfn5 <- rbind(d2,d3,d4)

#Full graphs
ggplot(dfn5, aes(name, RRMSE)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(1,1.5))+
  xlab(NULL)

#Comparing CV with default
t.test(dfn5$RRMSE[dfn5$name=="LineART"])
t.test(dfn5$RRMSE[dfn5$name=="LineART-CV"])
t.test(dfn5$RRMSE[dfn5$name=="SoftBART"])

mean((dfn5$RRMSE[dfn5$name=="LineART"])<1)
mean((dfn5$RRMSE[dfn5$name=="LineART-CV"])<1)
mean((dfn5$RRMSE[dfn5$name=="SoftBART"])<1)


d7 <- data.frame(RRMSE=sqrt(df_full$bart), name="BART")
d8 <- data.frame(RRMSE=sqrt(df_full$soft)/d7$RRMSE, name="SoftBART")
d9 <- data.frame(RRMSE=sqrt(df_full$lbart)/d7$RRMSE, name="LineART")
d10 <- data.frame(RRMSE=sqrt(df_full$lbartcv)/d7$RRMSE, name="LineART-CV")
dfn6 <- rbind(d8,d9,d10)

plot((dfn6$RRMSE[dfn6$name=="LineART"]),(dfn6$RRMSE[dfn6$name=="LineART-CV"]),
     xlab = "LineART", ylab = "LineART-CV",
     main = "RMSE divided by BART RMSE")
abline(0,1)

for(i in 1:41){
#i <- 41
  points((dfn5$RRMSE[dfn5$name=="LineART"])[(20*(i-1)+1):(20*(i-1)+20)],(dfn5$RRMSE[dfn5$name=="LineART-CV"])[(20*(i-1)+1):(20*(i-1)+20)], col = i, pch = 19)
}


pv <- NULL
est <- NULL
min1 <- NULL
max1 <- NULL
for(i in 1:41){
  pv <- c(pv,t.test(log(d9$RRMSE)[(20*(i-1)+1):(20*i)])$p.value)
  est <- c(est,t.test(log(d9$RRMSE)[(20*(i-1)+1):(20*i)])$estimate)
  min1 <- c(min1, min(log(d9$RRMSE)[(20*(i-1)+1):(20*i)]))
  max1 <- c(max1, max(log(d9$RRMSE)[(20*(i-1)+1):(20*i)]))
}

datasets <- list.files(path="C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/nipsdata", 
                       pattern = "x.txt", full.names=F, recursive=TRUE)
data_names <- unlist(strsplit(datasets, split='/', fixed=TRUE))[as.logical(1:82%%2)]


# plot(est, col =(pv<0.05)+2, pch=19,
#      xlab = "Dataset", ylab = "log(RMSE/RMSE_BART)",
#      main = "t-test for log(RMSE/RMSE_BART)")
# legend("bottomright", legend = c("p<0.05", "p>=0.05"), 
#        col = c("green","red"), pch = c(19,19))
# abline(h=0, col = "blue")
# abline(h=0.1, col = "blue")
# abline(h=-0.1, col = "blue")

df_tgraph <- data.frame(est=est, names=data_names, pv=pv, Min=min1, Max=max1)
df_tgraph$Significance <- ifelse(df_tgraph$pv < 0.05 & df_tgraph$est > 0, "p-value<0.05 and positive",
                          ifelse(df_tgraph$pv < 0.05 & df_tgraph$est < 0, "p-value<0.05 and negative",
                                 "p-value>0.05"))
pdf("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/LineART_code/pvalue_lbart.pdf")
ggplot(df_tgraph, aes(x = names, y = est, color = Significance)) +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "black") +
  geom_point(size = 3) +  # Scatter plot points
  scale_color_manual(values = c("p-value<0.05 and positive" = "red", 
                                "p-value<0.05 and negative" = "green", 
                                "p-value>0.05" = "blue")) +  # Custom colors
  xlab("") +
  ylab("log of RMSE ratio") +
  scale_y_continuous(limits = c(-2, 1))+
  geom_hline(yintercept = log(1/1.1), linetype = "dashed", color = "black") +  # First horizontal line
  geom_hline(yintercept = log(1/0.9), linetype = "dashed", color = "black") +  # Second horizontal line
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25),  # Rotate x-axis labels
    legend.position = c(0.8, 0.22),  # Move legend inside the plot (adjust as needed)
    legend.background = element_rect(fill = "white", color = "black")  # Add background to legend
  )
dev.off()

pv2 <- NULL
est2 <- NULL
min2 <- NULL
max2 <- NULL
for(i in 1:41){
  pv2 <- c(pv2,t.test(log(d8$RRMSE)[(20*(i-1)+1):(20*i)])$p.value)
  est2 <- c(est2,t.test(log(d8$RRMSE)[(20*(i-1)+1):(20*i)])$estimate)
  min2 <- c(min2, min(log(d8$RRMSE)[(20*(i-1)+1):(20*i)]))
  max2 <- c(max2, max(log(d8$RRMSE)[(20*(i-1)+1):(20*i)]))
}

df_tgraph2 <- data.frame(est=est2, names=data_names, pv=pv2, Min=min2, Max=max2)
df_tgraph2$Significance <- ifelse(df_tgraph2$pv < 0.05 & df_tgraph2$est > 0, "p-value<0.05 and positive",
                                  ifelse(df_tgraph2$pv < 0.05 & df_tgraph2$est < 0, "p-value<0.05 and negative",
                                         "p-value>0.05"))
pdf("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/LineART_code/pvalue_softbart.pdf")
ggplot(df_tgraph2, aes(x = names, y = est, color = Significance)) +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "black") +
  geom_point(size = 3) +  # Scatter plot points
  scale_color_manual(values = c("p-value<0.05 and positive" = "red", 
                                "p-value<0.05 and negative" = "green", 
                                "p-value>0.05" = "blue")) +  # Custom colors
  xlab("") +
  ylab("log of RMSE ratio") +
  scale_y_continuous(limits = c(-2, 1))+
  geom_hline(yintercept = log(1/1.1), linetype = "dashed", color = "black") +  # First horizontal line
  geom_hline(yintercept = log(1/0.9), linetype = "dashed", color = "black") +  # Second horizontal line
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25),  # Rotate x-axis labels
    legend.position = c(0.8, 0.22),  # Move legend inside the plot (adjust as needed)
    legend.background = element_rect(fill = "white", color = "black")  # Add background to legend
  )
dev.off()




plot((log(dfn6$RRMSE)[dfn6$name=="LineART"]),(log(dfn6$RRMSE)[dfn6$name=="LineART-CV"]),
     xlab = "LineART", ylab = "LineART-CV",
     main = "log of RMSE divided by BART RMSE")
abline(0,1)
abline(h=0, col = "red", lwd=2)
abline(v=0, col = "red", lwd=2)
legend("topleft", legend = c("2","14","20"), pch=c(19,19,19),
       col=c(2,14,20))

for(i in c(2,14,20)){
  #i <- 41
  points((log(dfn6$RRMSE)[dfn6$name=="LineART"])[(20*(i-1)+1):(20*(i-1)+20)],(log(dfn6$RRMSE)[dfn6$name=="LineART-CV"])[(20*(i-1)+1):(20*(i-1)+20)], col = i, pch = 19)
}


plot((log(dfn6$RRMSE)[dfn6$name=="LineART"]),(log(dfn6$RRMSE)[dfn6$name=="LineART-CV"]),
     xlab = "LineART", ylab = "LineART-CV",
     main = "log of RMSE divided by BART RMSE")
abline(0,1)
abline(h=0, col = "red", lwd=2)
abline(v=0, col = "red", lwd=2)
legend("topleft", legend = c("15","17","23"), pch=c(19,19,19),
       col=c("red","blue","purple"))

aux <- 0
for(i in c(15,17,23)){
  #i <- 41
  aux <- aux+1
  col_aux <- c("red","blue","purple")
  points((log(dfn6$RRMSE)[dfn6$name=="LineART"])[(20*(i-1)+1):(20*(i-1)+20)],(log(dfn6$RRMSE)[dfn6$name=="LineART-CV"])[(20*(i-1)+1):(20*(i-1)+20)], col = col_aux[aux], pch = 19)
}

df_data <- data.frame(dataset=df_full$dataset)
dfn7 <- data.frame(lineart=log(d9$RRMSE),
                   softbart=log(d8$RRMSE),
                   dataset=df_data) %>% group_by(dataset) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

pdf("rmse_ratio.pdf")
plot(dfn7$lineart, dfn7$softbart, xlim=c(-1.2,0.5),ylim=c(-1.2,0.5),
     ylab="Average Log of RMSE ratio (SoftBART/BART)",
     xlab="Average Log of RMSE ratio (LineART/BART)",
     cex = 1.5, col = "black", pch = 1)
points(dfn7$lineart[2], dfn7$softbart[2], cex = 1.5, col = "darkgreen", pch = 19)
points(dfn7$lineart[14], dfn7$softbart[14], cex = 1.5, col = "green", pch = 19)
points(dfn7$lineart[17], dfn7$softbart[17], cex = 1.5, col = "blue", pch = 19)
points(dfn7$lineart[23], dfn7$softbart[23], cex = 1.5, col = "purple", pch = 19)
abline(h=0, col = "red", lwd=2)
abline(v=0, col = "red", lwd=2)
abline(v=log(1/1.1), col = "red", lwd=2, lty = 2)
abline(v=log(1/0.9), col = "red", lwd=2, lty = 2)
legend("topleft", legend = c("Ais","Cpu","Diamond","Hatco"), 
       pch=c(19,19,19,19),
       col=c("darkgreen", "green", "blue","purple"))
dev.off()


