#Libraries
library(ggplot2)  #For graphs
library(ggpubr) #For theme
library(dplyr) #For data cleaning

#To calculate the Phi and Omega matrices
source("C:/Users/pedro/OneDrive/Área de Trabalho/Important_so_doing_it_again/Nature/Project_Stress/Dataset_1/Code/gp_approx_fun.R")

# #Load dataset
load(file=paste0("C:/Users/pedro/OneDrive/Área de Trabalho/Basis/study_1_long_randeff_reactivity_tau.Rdata"))

#To separate the variable t into epochs
time_graph = seq(min(data_1_clean$time), max(data_1_clean$time), length.out=100) #To be used in the graphs
end_times = data_1_clean %>% group_by(epoch) %>% summarize(beg = min(time), end=max(time))
epoch_times = end_times[order(end_times$end),]
epoch_grid = epoch_times$epoch[colSums(sapply(time_graph, function(x) x>=epoch_times$beg))]

#The grid that will be used in predictions (-1,1)
time_grid = seq(6, 21, length.out=100)

{
#Making the smooth curve using the time grid and the estimated coefs for control
for(j in 1:dim(fit_2_long$coefs_con_est)[3]){
  omega_con_grid = (Omega((time_grid-13.5)/7.5, dim(fit_2_long$coefs_con_est)[1], fit_2_long$L_con_post[j]) %*% sqrt(SDiag(1, fit_2_long$l_con_post[j], dim(fit_2_long$coefs_con_est)[1], fit_2_long$L_con_post[j])))
  if(j==1){
    control = lapply(1:dim(fit_2_long$coefs_con_est)[2],
                     function(i) t(t(omega_con_grid%*%fit_2_long$coefs_con_est[,i,j]) + fit_2_long$random_effects[j,i]))
  }else{
    control_aux = lapply(1:dim(fit_2_long$coefs_con_est)[2],
                      function(i) t(t(omega_con_grid%*%fit_2_long$coefs_con_est[,i,j]) + fit_2_long$random_effects[j,i]))

    for(i in 1:dim(fit_2_long$coefs_con_est)[2]){
      control[[i]] <- cbind(control[[i]],control_aux[[i]])
    }
  }
}

#Then calculate the Average Control curve posterior draws
acontrol = Reduce("+", control)/length(control)


#Same thing for ITE
#Making the smooth curve using the time grid and the estimated coefs for ite
for(j in 1:dim(fit_2_long$coefs_mod_est)[3]){
  omega_mod_grid = (Omega((time_grid-13.5)/7.5, dim(fit_2_long$coefs_mod_est)[1], fit_2_long$L_mod_post[j]) %*% sqrt(SDiag(1, fit_2_long$l_mod_post[j], dim(fit_2_long$coefs_mod_est)[1], fit_2_long$L_mod_post[j])))
  if(j==1){
    ite = lapply(1:dim(fit_2_long$coefs_mod_est)[2],
                     function(i) t(t(omega_mod_grid%*%fit_2_long$coefs_mod_est[,i,j])))
  }else{
    ite_aux = lapply(1:dim(fit_2_long$coefs_mod_est)[2],
                         function(i) t(t(omega_mod_grid%*%fit_2_long$coefs_mod_est[,i,j])))

    for(i in 1:dim(fit_2_long$coefs_mod_est)[2]){
      ite[[i]] <- cbind(ite[[i]],ite_aux[[i]])
    }
  }
}

#Then calculate the Average Treatment Effect curve posterior draws
ate = Reduce("+", ite)/length(ite)
}


#Aggregating stuff for the boxplot graphs
epoch_aggregate = function(draws, epoch_grid) {
  epochs = unique(epoch_grid)
  sapply(epochs, function(x) colMeans(draws[epoch_grid==x,]))
}
control_epoch = epoch_aggregate(acontrol, epoch_grid)
ate_epoch = epoch_aggregate(ate, epoch_grid)
trt_epoch = epoch_aggregate(acontrol+ate, epoch_grid)

ite_epoch = lapply(ite, function(x) epoch_aggregate(x, epoch_grid))
icontrol_epoch = lapply(control, function(x) epoch_aggregate(x, epoch_grid))

epoch_cate = sapply(ite_epoch, function(x) mean(x[,2]))
epoch_control = sapply(icontrol_epoch, function(x) mean(x[,2]))

#########
#The Graphs
#########

theme_set(theme_pubr())

get_bands = function(post_samples, dim=1, q = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975) ){
  out = t(apply(post_samples, dim, function(x) quantile(x, q)))
  colnames(out) = paste0("q",q)
  out = cbind(out, mean = apply(post_samples, dim, mean))
  out
}


#Graph A
##TPR vs Time (Control and Treated)

all_control = Reduce("+", control)/length(control)
all_ate = Reduce("+", ite)/length(ite)


control_band = data.frame(get_bands(all_control+mean(data_1_clean$reactivity_tpr)))
treat_band = data.frame(get_bands(all_ate+all_control+mean(data_1_clean$reactivity_tpr)))


control_band$Group = "Control"
control_band$epoch = epoch_grid
treat_band$Group = "Mindset Treatment"
treat_band$epoch = epoch_grid


#Creating the dataframe for the data
plotdf = data.frame(rbind(control_band,treat_band),
                    time=c(time_grid, time_grid))
plotdf$Group <- factor(plotdf$Group, levels = c("Control", "Mindset Treatment"))
plotdf$epoch <- factor(plotdf$epoch, levels=c("prep", "speech", "math", "recov"))

y_lim_min = min(plotdf$q0.025)
y_lim_max = max(plotdf$q0.975)+150

#Creating the labels to separate the epochs on the graph
label_df = data.frame(
  label_y = y_lim_max,
  label_x = ((epoch_times$beg-c(0,0.5,0.5,0.5)) + (epoch_times$end+c(0.5,0.5,0.5,0)))/2,
  label  = c("Prep.", "Speech", "Math", "Recov.")
)

color1 <- "blue"
color2 <- "red"


graphA <- ggplot(aes(x=time, y=mean), data=plotdf) +
  geom_vline(xintercept = epoch_times$beg[-1]-0.5, linetype="dotted", col="gray30") +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), fill="gray70", data = plotdf[plotdf$Group=="Control",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), fill="gray70", data = plotdf[plotdf$Group!="Control",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.1, ymax=q0.9, fill=Group), alpha = 0.2) +
  geom_line(aes(color=Group)) +
  geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  scale_colour_manual(values = c(color1, color2)) +
  scale_fill_manual(values = c(color1, color2)) +
  scale_y_continuous(limits = c(y_lim_min-50,y_lim_max+50) #, expand = c(0, 0)
  ) +
  scale_x_continuous(limits = c(min(data_1_clean$time),max(data_1_clean$time))
                     #, expand = c(0, 0)
  ) +
  xlab("Time") +
  ylab("Reactivity - Total Peripheral Resistance") + #labs(title=paste0("Study 1")) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "Reactivity - Total peripheral resistance (R-TPR) levels")
graphA


#Treatment Effect vs Time

ate_band = data.frame(get_bands(all_ate))

ate_band$epoch = epoch_grid

#Creating the dataframe for the data
plotdf = data.frame(rbind(ate_band),
                    time=c(time_grid))
plotdf$epoch <- factor(plotdf$epoch, levels=c("Prep.", "Speech", "Math", "Recov."))

y_lim_min = min(plotdf$q0.025)
y_lim_max = max(plotdf$q0.975)+25

#Creating the labels to separate the epochs on the graph
label_df = data.frame(
  label_y = y_lim_max,
  label_x = epoch_times$beg + (epoch_times$end + 1 - epoch_times$beg)/2,
  label  = c("Prep.", "Speech", "Math", "Recov.")
)


graphB <- ggplot(aes(x=time, y=mean), data=plotdf) +
  geom_vline(xintercept = epoch_times$beg[-1], linetype="dotted", col="gray30") +
  geom_hline(yintercept = 0, linetype="dotted", col="gray30") +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha = 0.2, fill = "gray70") +
  geom_ribbon(aes(ymin=q0.1, ymax=q0.9), alpha = 0.2) +
  geom_line() +
  geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  ylim(y_lim_min,y_lim_max) +
  xlab("Time") +
  ylab("CATE") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "CATE over time")
graphB


####################################################################

#Getting the indexes of males and females
male_idx <- which(person_df$sex==1)
female_idx <- which(person_df$sex==2)

#Now I get the control and ate for each sex
all_control_male = Reduce("+", control[male_idx])/length(control[male_idx])
all_ate_male = Reduce("+", ite[male_idx])/length(ite[male_idx])

all_control_female = Reduce("+", control[female_idx])/length(control[female_idx])
all_ate_female = Reduce("+", ite[female_idx])/length(ite[female_idx])

#Getting the quantiles for each one of the 100 points of the curve
ate_band_male = data.frame(get_bands(all_ate_male))
ate_band_female = data.frame(get_bands(all_ate_female))

#Now I add the grid of epoch to help in the graph
#Creating the groups
ate_band_male$Group = "Male"
ate_band_male$epoch = epoch_grid
ate_band_female$Group = "Female"
ate_band_female$epoch = epoch_grid

#Creating the dataframes for the data
plotdf = data.frame(rbind(ate_band_male, ate_band_female),
                      time=c(time_grid, time_grid))
plotdf$epoch <- factor(plotdf$epoch, levels=c("prep", "speech", "math", "recov"))

#Just adjusting the graph limits here (in this case is male, but can change if needed)
y_lim_min = min(plotdf$q0.025)
y_lim_max = max(plotdf$q0.975)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

#Creating the labels to separate the epochs on the graph
label_df = data.frame(
  label_y = y_lim_max,
  label_x = epoch_times$beg + (epoch_times$end + 1 - epoch_times$beg)/2,
  label  = c("Prep.", "Speech", "Math", "Recov.")
)


graphC <- ggplot(aes(x=time, y=mean), data=plotdf) +
  geom_vline(xintercept = epoch_times$beg[-1], linetype="dotted", col="gray") +
  geom_hline(yintercept = 0, linetype="dotted", col="gray") +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), fill="gray70", data = plotdf[plotdf$Group=="Male",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), fill="gray70", data = plotdf[plotdf$Group!="Male",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.1, ymax=q0.9, fill=Group), alpha = 0.2) +
  geom_line(aes(color=Group)) +
  geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  ylim(y_lim_min,y_lim_max) +
  xlab("Time") +
  ylab("CATE") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "CATE by sex over time")
graphC



#Do the difference in Cates
all_atediff_sex <- all_ate_male-all_ate_female
#Same thing for previous graph, but for this dataframe
atesexdif_band = data.frame(get_bands(all_atediff_sex))

atesexdif_band$epoch = epoch_grid

#Creating the dataframe for the data
plotdf = data.frame(rbind(atesexdif_band),
                    time=c(time_grid))
plotdf$epoch <- factor(plotdf$epoch, levels=c("Prep.", "Speech", "Math", "Recov."))

y_lim_min = min(plotdf$q0.025)
y_lim_max = max(plotdf$q0.975)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

#Creating the labels to separate the epochs on the graph
label_df = data.frame(
  label_y = y_lim_max,
  label_x = epoch_times$beg + (epoch_times$end + 1 - epoch_times$beg)/2,
  label  = c("Prep.", "Speech", "Math", "Recov.")
)


graphD <- ggplot(aes(x=time, y=mean), data=plotdf) +
  geom_vline(xintercept = epoch_times$beg[-1], linetype="dotted", col="gray30") +
  geom_hline(yintercept = 0, linetype="dotted", col="gray30") +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha = 0.2, fill = "gray70") +
  geom_ribbon(aes(ymin=q0.1, ymax=q0.9), alpha = 0.2) +
  geom_line() +
  geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  ylim(y_lim_min,y_lim_max) +
  xlab("Time") +
  ylab("CATE") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "Difference of CATE by sex over time")
graphD


#Worth it to do the violin plots here to get the probability
ate_epoch = epoch_aggregate(all_atediff_sex, epoch_grid)
ate_epoch_bp <- data.frame((get_bands(t(ate_epoch))))
ate_epoch_bp$epoch <- factor(c("Prep.", "Speech", "Math", "Recov.") , levels=c("Prep.", "Speech", "Math", "Recov.") )

ate_df <- data.frame(mean = as.vector(ate_epoch), epoch = factor(rep(c("Prep.", "Speech", "Math", "Recov."), each = nrow(ate_epoch)),levels = c("Prep.", "Speech", "Math", "Recov.")))

#The probabilities are pretty high tho
probs_diff <- apply((ate_epoch)>0,2,function(x){sum(x)/length(x)})

y_lim_min = min(ate_df$mean)
y_lim_max = max(ate_df$mean)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

graphE <- ggplot(aes(x=0, y=mean), data=ate_epoch_bp) +
  geom_violin(data = ate_df)+
  geom_point(size=1.2) + facet_grid(~epoch) +
  theme_pubr() +
  geom_segment(aes(x=0, xend=0, y=q0.025, yend=q0.975), size=0.5, color = "gray70") +
  geom_segment(aes(x=0, xend=0, y=q0.1, yend=q0.9), size=1.2) +
  ylim(y_lim_min,y_lim_max) +
  geom_hline(yintercept = 0, linetype="dotted", col="gray30")  +
  ylab("CATE") + xlab("") +
  geom_text(data =  data.frame(epoch = factor("Prep."),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[1],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Speech"),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[2],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Math"),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[3],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Recov."),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[4],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  theme(axis.ticks.x=element_blank(), legend.title = element_blank(),
        axis.text.x=element_blank(),
        legend.position="bottom") +
  labs(title = "Diff. of CATE by sex for each epoch")
graphE



#For doublebad and doublegood I need to separate the subgroups by the person_df dataframe
#Separating the mindsets
fixedmind <- person_df$fixedmindset_baseline
stressmind <- person_df$stressmindset_baseline

#Function to find available cutpoints
#In this case I'm doing a CART-like search, so I need the cutpoints
#Cutpoints will be just the average between two points
middle <- function(x){
  aux <- numeric(0)
  y <- unique(sort(x))
  for(i in 1:(length(y)-1)){
    aux[i] <- (y[i]+y[(i+1)])/2
  }
  return(aux)
}

#I need a grid with all available choices of cuts that can be performed
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
cuts_fixed <- middle(fixedmind)
cuts_stress <- middle(stressmind)

#Create grid of a<b, since I have two groups, I will need two cutpoints for each variable
grid_x1 <- expand.grid.unique(cuts_fixed,cuts_fixed)
grid_x2 <- expand.grid.unique(cuts_stress,cuts_stress)

#And include the case where a=b
grid_x1 <- rbind(grid_x1,cbind(cuts_fixed,cuts_fixed))
grid_x2 <- rbind(grid_x2,cbind(cuts_stress,cuts_stress))

#Creating a grid of cutpoints
#In this case I have matrices, and I want to try combinations of the lines
#So I will create a grid of indexes and use those

ind_1 <- 1:(nrow(grid_x1))
ind_2 <- 1:(nrow(grid_x2))

#Since I'm combining the two grids, I can use the expand.grid
grid_index <- expand.grid(ind_1,ind_2)

#Finding the best cuts
#The rule here is to exhaust every possible combination of points and test all of them
#I want to get the two most distinct subgroups I can have
#I need to set a minimum size to that group, otherwise I will get a single value
#The way the difference is calculated is by simply taking the average of y between subgroups
search_dual <- function(y, data, grid_idx, grid1, grid2, prop = 0.1){
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

    #Separating the groups (the equal should never happen, but I include it anyway)
    g1 <- (data$fixedmindset_baseline<(x1_aux[1]))&((data)$stressmindset_baseline<(x2_aux[1]))
    g2 <- (data$fixedmindset_baseline>=(x1_aux[2]))&((data)$stressmindset_baseline>=(x2_aux[2]))

    #Here I test if the groups are big enough (depending in the prop argumnt)
    g_test <- (sum(g1)>=min_n)&(sum(g2)>=min_n)

    #calculating the difference
    if(g_test==T){
      #Calculate the absolute value of the difference of means the groups
      dif[i] <- abs(diff(c(mean(y[g1]),mean(y[g2]))))
    }else{
      #Otherwise, set it to be zero, since I'm woring with abs values, this is enough
      dif[i] <- 0
    }

  }

  return(c(dif))
}

#Separete by using only speech epoch
#Notice that here we use the estimated control, so we are using whatever our model gave us as responses
#But we do NOT include the treatment effect, so we are trying to find subgroups of treatment effects using only the characteristics found on control
dif_vec <- search_dual(y = sapply(lapply(control, function(x) epoch_aggregate(x, epoch_grid)[,2]), mean),
                       data = person_df, grid_idx = grid_index, grid1 = grid_x1, grid2 = grid_x2,
                       prop = 0.1)
if(sum(dif_vec>0)){
  print("At least one good cutpoint was found")
}else{
  print("Something is wrong, no valid cutpoints were chosen")
}

#Getting the actual cutpoints
indexes_max <- dif_vec==max(dif_vec)
fixed_cut = grid_x1[as.numeric((grid_index[indexes_max,])[,1]), , drop = F]
stress_cut = grid_x2[as.numeric((grid_index[indexes_max,])[,2]), , drop = F]

#Save inside a list to keep track of the cutpoints selected
cuts_list <- list()
for(i in 1:sum(indexes_max)){
  cuts_list[[i]] <- data.frame(fixed = fixed_cut[i,], stress = stress_cut[i,])
}

#Now I create the subgroups
cuts <- cuts_list[[1]]
#write.table(cuts, "cuts.txt")
doublegood <- (person_df$fixedmindset_baseline<(cuts[1,1]))&(person_df$stressmindset_baseline<(cuts[1,2]))
doublebad <- (person_df$fixedmindset_baseline>=(cuts[2,1]))&(person_df$stressmindset_baseline>=(cuts[2,2]))

#Just counting the number of observation in each subgroup compared to the total num. of obs.
c(sum(doublegood),sum(doublebad),nrow(person_df))

#Now I get the control and ate for each sex
all_control_db = Reduce("+", control[doublebad])/length(control[doublebad])
all_ate_db = Reduce("+", ite[doublebad])/length(ite[doublebad])

all_control_dg = Reduce("+", control[doublegood])/length(control[doublegood])
all_ate_dg = Reduce("+", ite[doublegood])/length(ite[doublegood])


#Getting the quantiles for each one of the 100 points of the curve
ate_band_db = data.frame(get_bands(all_ate_db))
ate_band_dg = data.frame(get_bands(all_ate_dg))

#Now I add the grid of epoch to help in the graph
ate_band_db$Group = "Negative prior mindsets"
ate_band_db$epoch = epoch_grid
ate_band_dg$Group = "Positive prior mindsets"
ate_band_dg$epoch = epoch_grid

#Creating the dataframes for the data
plotdf = data.frame(rbind(ate_band_db,ate_band_dg),
                      time=c(time_grid,time_grid))
plotdf$epoch <- factor(plotdf$epoch, levels=c("prep", "speech", "math", "recov"))

#Just adjusting the graph limits here
y_lim_min = min(plotdf$q0.025)
y_lim_max = max(plotdf$q0.975)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

#Creating the labels to separate the epochs on the graph
label_df = data.frame(
  label_y = y_lim_max,
  label_x = epoch_times$beg + (epoch_times$end + 1 - epoch_times$beg)/2,
  label  = c("Prep.", "Speech", "Math", "Recov.")
)


graphF <- ggplot(aes(x=time, y=mean), data=plotdf) +
  geom_vline(xintercept = epoch_times$beg[-1], linetype="dotted", col="gray") +
  geom_hline(yintercept = 0, linetype="dotted", col="gray") +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), fill="gray70", data = plotdf[plotdf$Group=="Negative prior mindsets",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), fill="gray70", data = plotdf[plotdf$Group!="Negative prior mindsets",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.1, ymax=q0.9, fill=Group), alpha = 0.2) +
  geom_line(aes(color=Group)) +
  geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  ylim(y_lim_min,y_lim_max) +
  xlab("Time") +
  ylab("CATE") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "CATE by subgroup over time")
graphF



#Do the difference in Cates
all_atediff_sg <- all_ate_dg-all_ate_db
#Same thing for previous graph, but for this dataframe
atesgdif_band = data.frame(get_bands(all_atediff_sg))

atesgdif_band$epoch = epoch_grid

#Creating the dataframe for the data
plotdf = data.frame(rbind(atesgdif_band),
                    time=c(time_grid))
plotdf$epoch <- factor(plotdf$epoch, levels=c("Prep.", "Speech", "Math", "Recov."))

y_lim_min = min(plotdf$q0.025)
y_lim_max = max(plotdf$q0.975)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

#Creating the labels to separate the epochs on the graph
label_df = data.frame(
  label_y = y_lim_max,
  label_x = epoch_times$beg + (epoch_times$end + 1 - epoch_times$beg)/2,
  label  = c("Prep.", "Speech", "Math", "Recov.")
)


graphG <- ggplot(aes(x=time, y=mean), data=plotdf) +
  geom_vline(xintercept = epoch_times$beg[-1], linetype="dotted", col="gray30") +
  geom_hline(yintercept = 0, linetype="dotted", col="gray30") +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha = 0.2, fill = "gray70") +
  geom_ribbon(aes(ymin=q0.1, ymax=q0.9), alpha = 0.2) +
  geom_line() +
  geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  ylim(y_lim_min,y_lim_max) +
  xlab("Time") +
  ylab("CATE") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "Difference of CATE by subgroup over time")
graphG


#Worth it to do the violin plots here to get the probability

ate_epoch = epoch_aggregate(all_atediff_sg, epoch_grid)
ate_epoch_bp <- data.frame((get_bands(t(ate_epoch))))
ate_epoch_bp$epoch <- factor(c("Prep.", "Speech", "Math", "Recov.") , levels=c("Prep.", "Speech", "Math", "Recov.") )

ate_df <- data.frame(mean = as.vector(ate_epoch), epoch = factor(rep(c("Prep.", "Speech", "Math", "Recov."), each = nrow(ate_epoch)),levels = c("Prep.", "Speech", "Math", "Recov.")))

#The probabilities are pretty high tho
probs_diff <- apply((ate_epoch)>0,2,function(x){sum(x)/length(x)})

y_lim_min = min(ate_df$mean)
y_lim_max = max(ate_df$mean)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

graphH <- ggplot(aes(x=0, y=mean), data=ate_epoch_bp) +
  geom_violin(data = ate_df)+
  geom_point(size=1.2) + facet_grid(~epoch) +
  theme_pubr() +
  geom_segment(aes(x=0, xend=0, y=q0.025, yend=q0.975), size=0.5, color = "gray70") +
  geom_segment(aes(x=0, xend=0, y=q0.1, yend=q0.9), size=1.2) +
  ylim(y_lim_min,y_lim_max) +
  geom_hline(yintercept = 0, linetype="dotted", col="gray30")  +
  ylab("CATE") + xlab("") +
  geom_text(data =  data.frame(epoch = factor("Prep."),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[1],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Speech"),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[2],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Math"),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[3],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Recov."),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[4],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  theme(axis.ticks.x=element_blank(), legend.title = element_blank(),
        axis.text.x=element_blank(),
        legend.position="bottom") +
  labs(title = "Diff. of CATE by subgroup for each epoch")
graphH


##########################################################################
#Now for the self esteem thing
ptr = rpart::rpart(epoch_control
                   ~stressmindset_baseline+fixedmindset_baseline+sex+selfesteem_baseline,
                   data=person_df,
                   model = T, control=list(maxdepth=4, cp = 0.0001))
pdf(paste0("tree_study_1.pdf"), width=6, height=4)
rpart.plot::prp(ptr, extra = 100, main = paste0("Study 1 - Decision Tree - Speech - R-TPR"))
dev.off()

#Creating the self-esteem subgroups
lowse <- (person_df$selfesteem_baseline<3.9)
highse <- (person_df$selfesteem_baseline>=3.9)

#Just counting the number of observation in each subgroup compared to the total num. of obs.
c(sum(lowse),sum(highse),nrow(person_df))

#Now I get the control and ate for each sex
all_control_lse = Reduce("+", control[lowse])/length(control[lowse])
all_ate_lse = Reduce("+", ite[lowse])/length(ite[lowse])

all_control_hse = Reduce("+", control[highse])/length(control[highse])
all_ate_hse = Reduce("+", ite[highse])/length(ite[highse])


#Getting the quantiles for each one of the 100 points of the curve
ate_band_lse = data.frame(get_bands(all_ate_lse))
ate_band_hse = data.frame(get_bands(all_ate_hse))

#Now I add the grid of epoch to help in the graph
ate_band_lse$Group = "Low self-esteem"
ate_band_lse$epoch = epoch_grid
ate_band_hse$Group = "High self-esteem"
ate_band_hse$epoch = epoch_grid

#Creating the dataframes for the data
plotdf = data.frame(rbind(ate_band_lse,ate_band_hse),
                    time=c(time_grid,time_grid))
plotdf$epoch <- factor(plotdf$epoch, levels=c("prep", "speech", "math", "recov"))

#Just adjusting the graph limits here
y_lim_min = min(plotdf$q0.025)
y_lim_max = max(plotdf$q0.975)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

#Creating the labels to separate the epochs on the graph
label_df = data.frame(
  label_y = y_lim_max,
  label_x = epoch_times$beg + (epoch_times$end + 1 - epoch_times$beg)/2,
  label  = c("Prep.", "Speech", "Math", "Recov.")
)


graphI <- ggplot(aes(x=time, y=mean), data=plotdf) +
  geom_vline(xintercept = epoch_times$beg[-1], linetype="dotted", col="gray") +
  geom_hline(yintercept = 0, linetype="dotted", col="gray") +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), fill="gray70", data = plotdf[plotdf$Group=="Low self-esteem",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), fill="gray70", data = plotdf[plotdf$Group!="Low self-esteem",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.1, ymax=q0.9, fill=Group), alpha = 0.2) +
  geom_line(aes(color=Group)) +
  geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  ylim(y_lim_min,y_lim_max) +
  xlab("Time") +
  ylab("CATE") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "CATE by subgroup over time")
graphI



#Do the difference in Cates
all_atediff_se <- all_ate_hse-all_ate_lse
#Same thing for previous graph, but for this dataframe
atesedif_band = data.frame(get_bands(all_atediff_se))

atesedif_band$epoch = epoch_grid

#Creating the dataframe for the data
plotdf = data.frame(rbind(atesedif_band),
                    time=c(time_grid))
plotdf$epoch <- factor(plotdf$epoch, levels=c("Prep.", "Speech", "Math", "Recov."))

y_lim_min = min(plotdf$q0.025)
y_lim_max = max(plotdf$q0.975)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

#Creating the labels to separate the epochs on the graph
label_df = data.frame(
  label_y = y_lim_max,
  label_x = epoch_times$beg + (epoch_times$end + 1 - epoch_times$beg)/2,
  label  = c("Prep.", "Speech", "Math", "Recov.")
)


graphJ <- ggplot(aes(x=time, y=mean), data=plotdf) +
  geom_vline(xintercept = epoch_times$beg[-1], linetype="dotted", col="gray30") +
  geom_hline(yintercept = 0, linetype="dotted", col="gray30") +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha = 0.2, fill = "gray70") +
  geom_ribbon(aes(ymin=q0.1, ymax=q0.9), alpha = 0.2) +
  geom_line() +
  geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  ylim(y_lim_min,y_lim_max) +
  xlab("Time") +
  ylab("CATE") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "Difference of CATE by self-esteem levels over time")
graphJ


#Worth it to do the violin plots here to get the probability

ate_epoch = epoch_aggregate(all_atediff_se, epoch_grid)
ate_epoch_bp <- data.frame((get_bands(t(ate_epoch))))
ate_epoch_bp$epoch <- factor(c("Prep.", "Speech", "Math", "Recov.") , levels=c("Prep.", "Speech", "Math", "Recov.") )

ate_df <- data.frame(mean = as.vector(ate_epoch), epoch = factor(rep(c("Prep.", "Speech", "Math", "Recov."), each = nrow(ate_epoch)),levels = c("Prep.", "Speech", "Math", "Recov.")))

#The probabilities are pretty high tho
probs_diff <- apply((ate_epoch)>0,2,function(x){sum(x)/length(x)})

y_lim_min = min(ate_df$mean)
y_lim_max = max(ate_df$mean)
y_lim_range <- y_lim_max - y_lim_min
y_lim_min <- y_lim_min - y_lim_range*0.1
y_lim_max <- y_lim_max + y_lim_range*0.1

graphK <- ggplot(aes(x=0, y=mean), data=ate_epoch_bp) +
  geom_violin(data = ate_df)+
  geom_point(size=1.2) + facet_grid(~epoch) +
  theme_pubr() +
  geom_segment(aes(x=0, xend=0, y=q0.025, yend=q0.975), size=0.5, color = "gray70") +
  geom_segment(aes(x=0, xend=0, y=q0.1, yend=q0.9), size=1.2) +
  ylim(y_lim_min,y_lim_max) +
  geom_hline(yintercept = 0, linetype="dotted", col="gray30")  +
  ylab("CATE") + xlab("") +
  geom_text(data =  data.frame(epoch = factor("Prep."),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[1],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Speech"),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[2],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Math"),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[3],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  geom_text(data =  data.frame(epoch = factor("Recov."),
                               mean = y_lim_max),
            label = paste0("P(Y>0| . )=",round(probs_diff[4],3)*100,"%"),
            color = "gray30",
            show.legend = FALSE)+
  theme(axis.ticks.x=element_blank(), legend.title = element_blank(),
        axis.text.x=element_blank(),
        legend.position="bottom") +
  labs(title = "Diff. of CATE by self-esteem levels for each epoch")
graphK




#Saving graphs
{
  ggsave(paste0("graphA.pdf"),
         graphA, width=6, height=4)
  ggsave(paste0("graphB.pdf"),
         graphB, width=6, height=4)
  ggsave(paste0("graphC.pdf"),
         graphC, width=6, height=4)
  ggsave(paste0("graphD.pdf"),
         graphD, width=6, height=4)
  ggsave(paste0("graphE.pdf"),
         graphE, width=6, height=4)
  ggsave(paste0("graphF.pdf"),
         graphF, width=6, height=4)
  ggsave(paste0("graphG.pdf"),
         graphG, width=6, height=4)
  ggsave(paste0("graphH.pdf"),
         graphH, width=6, height=4)
  ggsave(paste0("graphI.pdf"),
         graphI, width=6, height=4)
  ggsave(paste0("graphJ.pdf"),
         graphJ, width=6, height=4)
  ggsave(paste0("graphK.pdf"),
         graphK, width=6, height=4)
}



