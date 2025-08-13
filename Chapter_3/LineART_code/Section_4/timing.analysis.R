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
  tbl_lbart <- read.csv(paste0("times_cv_lbart__a_", alpha, "__b_",beta,"__nt_",nt,"__i_",j,"_total.csv"), sep = ";", header = F)[-1,2] %>% as.numeric()
  tbl_bart <- read.csv(paste0("times_cv_bart__a_0.95__b_2__nt_200__i_",j,"_total.csv"), sep = ";", header = F)[-1,3] %>% as.numeric()
  tbl_soft <- read.csv(paste0("times_cv_soft__s_off__i_",j,"_full.csv"), sep = ";", header = F)[-1,4] %>% as.numeric()
  
  result_table <- NULL
  for(i in 1:length(tbl_lbart)){
    row_i <- t(tbl_lbart)[i]
    row_i2 <- t(tbl_bart)[i]
    row_i3 <- t(tbl_soft)[i]
    
    alldata <- rbind(alldata, data.frame(bart = row_i2, lbart = row_i, soft = row_i3,
                                         a=alpha, b=beta, nt=nt,
                                         dataset=j, rep = i))
    
  }
  
}


#Where the datasets are located
setwd("C:/Users/pedro/OneDrive/Área de Trabalho/Datasets/data/nipsdata")

#Reading files
x_files <- list.files(path=getwd(), pattern = "x.txt", full.names=TRUE, recursive=TRUE)
y_files <- list.files(path=getwd(), pattern = "y.txt", full.names=TRUE, recursive=TRUE)

#Saving those in lists so its gonna be easier to make a for loop
y <- list()
y <- lapply(y_files, function(x){read.table(x, header=F)})

#Just to compare sizes and variables with BART paper (LABOR is missing, so 41 datasets)
size_y <- matrix(unlist(lapply(y, dim)), ncol = 2, byrow = T)[,1]
n_train <- round(size_y*5/6,0)


for(i in 1:41){
  alldata$n[alldata$dataset==i] <- n_train[i]
}

df_full <- alldata

#df_aux1 <- data.frame(Time = df_full$lbart, name="LineART", dataset = df_full$dataset, rep = df_full$rep, n = df_full$n)
df_aux3 <- data.frame(Time = df_full$bart/df_full$lbart, name="BART", dataset = df_full$dataset, rep = df_full$rep, n = df_full$n)
df_aux4 <- data.frame(Time = df_full$soft/df_full$lbart, name="SoftBART", dataset = df_full$dataset, rep = df_full$rep, n = df_full$n)
df_aux <- rbind(df_aux3, df_aux4)
df_aux <- df_aux %>% select(-rep, -dataset)


# ggplot(df_aux, aes(Time, group = n)) +
#   geom_boxplot() +
#   coord_flip()

b_df <- data.frame(logn=log(df_aux$n[df_aux$name=="BART"]), time=df_aux$Time[df_aux$name=="BART"])
s_df <- data.frame(logn=log(df_aux$n[df_aux$name=="SoftBART"]),  time=df_aux$Time[df_aux$name=="SoftBART"])


# pdf("timing.pdf")
# plot(b_df$logn, b_df$time, ylim = c(0,35), col = "purple",
#      xlab = "log(n)", ylab = "Time Ratio")
# points(s_df$logn, s_df$time, col = "red")
# abline(h=1, lwd = 2, col = "blue")
# legend("topleft", legend = c("BART", "SoftBART"), pch = c(19,19), col = c("purple", "red"))
# dev.off()

b_df$name <- "BART"
s_df$name <- "SoftBART"
time_df <- rbind(b_df, s_df)



ranges <- function(x){
  return(data.frame(ymin = min(x), y = mean(x), ymax = max(x)))
}
colnames(time_df)[3] <- "Competitor"
pdf("C:/Users/pedro/OneDrive/Documentos/Code/Linear/lineart_20250117/LineART_code/timing_ggplot.pdf")
ggplot(time_df, aes(logn, time, colour = Competitor)) +
  stat_summary(fun.data=ranges, geom="pointrange") +
  xlab("log(n)") +
  ylab("Time Ratio") +
  geom_abline(slope = 0, intercept = 1, colour = "purple") +
  theme(
    legend.position = c(0.13, 0.85),  # Adjust (x, y) for placement inside the plot
    legend.background = element_rect(fill = "white", color = "black")  # Optional: add background
  )
dev.off()












