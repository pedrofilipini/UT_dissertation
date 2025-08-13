set.seed(123456)

#sq_exponential
sqe_kernel <- function(x1, x2, ls) {
  exp(-(x1-x2)^2/(2*ls^2))
}
par(mfrow=c(1,1))
plot(seq(-1,1, length.out = 100), sqe_kernel(-1, seq(-1,1, length.out = 100),0.13), type = "l")


library(magrittr)
#GP example
n <- 100

set.seed(1234)
t <- seq(-1,1, length.out = n) #my t variable
t_train_person <- c(t,t)

x_train <- runif(n, 0,1) %>% as.matrix(ncol = 1)
x_train <- ifelse(x_train<0.5,0.25,0.75)
x_train_person <- matrix(c(rep(0.25, times = n),rep(0.75, times = n)),ncol=1)


#Using beta instead
t1 <- rbeta(n, shape1 = 2, shape2 = 5)
#curve(dbeta(x, shape1 = 2, shape2 = 5))
t_train <- 2*((t1-(max(t1)+min(t1))/2)/(max(t1)-min(t1)) ) %>% sort()

t_train_person <- c(seq(-1,1, length.out=100),seq(-1,1, length.out=100))

curve(sin(x*pi), from = 0, to = 1)
curve(sin(-0.1*x*pi), from = 0, to = 1, add = T, col = "blue")
y1 <- sin(1*t_train*pi)
y2 <- sin(-0.1*t_train*pi)+0.3


y <- y1*(x_train<0.5) + y2*(x_train>=0.5) + rnorm(n, sd = 0.1)



fit_1 <- multibart:::bcf_core(y,
                              z_train = c(rep(0, times = (dim(x_train)[1]-1)),1),
                              z_out = rep(1, times = dim(x_train_person)[1]),
                              tvar_con = t_train,
                              tvar_mod = t_train,
                              x_control = x_train,
                              x_moderate = x_train, #placeholder
                              x_control_out = x_train_person,
                              x_moderate_out = x_train_person, #placeholder
                              tvar_con_out = t_train_person,
                              tvar_mod_out = t_train_person,
                              vanilla = FALSE, dart = F,
                              ntree_moderate = 0,
                              ntree_control = 100,
                              update_ls = T, enc_con = 5,
                              nburn = 1000, nsim = 3000, nthin = 1)

fit_2 <- multibart:::bcf_core(y,
                              z_train = c(rep(0, times = (dim(x_train)[1]-1)),1),
                              z_out = rep(1, times = dim(x_train_person)[1]),
                              tvar_con = t_train,
                              tvar_mod = t_train,
                              x_control = x_train,
                              x_moderate = x_train, #placeholder
                              x_control_out = x_train_person,
                              x_moderate_out = x_train_person, #placeholder
                              tvar_con_out = t_train_person,
                              tvar_mod_out = t_train_person,
                              vanilla = FALSE, dart = F,
                              ntree_moderate = 0,
                              ntree_control = 100,
                              update_ls = F, enc_con = 5,
                              nburn = 1000, nsim = 3000, nthin = 1)


plot(density(fit_1$l_con_post))
plot.ts((fit_1$l_con_post))
# fit_2$l_con_post
#
# dim(fit_1$mu_est_post)

#Variable ls

p1 <- apply(fit_1$mu_est_post[,1:100],2,mean)
p2 <- apply(fit_1$mu_est_post[,1:100],2,quantile, 0.025)
p3 <- apply(fit_1$mu_est_post[,1:100],2,quantile, 0.975)

p4 <- apply(fit_1$mu_est_post[,101:200],2,mean)
p5 <- apply(fit_1$mu_est_post[,101:200],2,quantile, 0.025)
p6 <- apply(fit_1$mu_est_post[,101:200],2,quantile, 0.975)



#Fixed ls

p1b <- apply(fit_2$mu_est_post[,1:100],2,mean)
p2b <- apply(fit_2$mu_est_post[,1:100],2,quantile, 0.025)
p3b <- apply(fit_2$mu_est_post[,1:100],2,quantile, 0.975)

p4b <- apply(fit_2$mu_est_post[,101:200],2,mean)
p5b <- apply(fit_2$mu_est_post[,101:200],2,quantile, 0.025)
p6b <- apply(fit_2$mu_est_post[,101:200],2,quantile, 0.975)

par(mfrow=c(2,2))

curve(sin(x*pi), from = -1, to = 1, ylim =c(-2,2), main = "Variable")
lines(t, p1, col = "red")
lines(t, p2, col = "purple")
lines(t, p3, col = "purple")


curve(sin(x*pi), from = -1, to = 1, ylim =c(-2,2), main = "Fixed")
lines(t, p1b, col = "red")
lines(t, p2b, col = "purple")
lines(t, p3b, col = "purple")


curve(sin(-0.1*x*pi)+0.3, from = -1, to = 1, ylim =c(-2,2), main = "Variable")
lines(t, p4, col = "red")
lines(t, p5, col = "purple")
lines(t, p6, col = "purple")


curve(sin(-0.1*x*pi)+0.3, from = -1, to = 1, ylim =c(-2,2), main = "Fixed")
lines(t, p4b, col = "red")
lines(t, p5b, col = "purple")
lines(t, p6b, col = "purple")

# plot(t_train,y)
# plot(t_train,y1)
# plot(t_train,y2)

library(ggplot2)

plotdf <- data.frame(mean=c(p1,p1b), time = c(t,t), Group = c(rep("Variable", times = n),rep("Fixed", times = n)),
           q0.025=c(p2,p2b),q0.975=c(p3,p3b))
graphA_ls <- ggplot(aes(x=time, y=mean), data=plotdf) +
  stat_function(fun = function(x) {sin(1 * pi * x)},
                color = "black", linetype = "longdash", size = 1.5) +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975, fill=Group), data = plotdf[plotdf$Group=="Variable",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975, fill=Group), data = plotdf[plotdf$Group!="Variable",], alpha = 0.2) +
  #geom_ribbon(aes(ymin=q0.1, ymax=q0.9, fill=Group), alpha = 0.2) +
  geom_line(aes(color=Group)) +
  #geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  scale_colour_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  #scale_y_continuous(limits = c(y_lim_min-50,y_lim_max+50) #, expand = c(0, 0)
  #) +
  #scale_x_continuous(limits = c(min(data_1_clean$time),max(data_1_clean$time))
  #, expand = c(0, 0)
  #) +
  xlab("t") +
  ylab("Y") + #labs(title=paste0("Study 1")) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "Response over targeted smoothing variable: X<0.5")
graphA_ls


plotdf2 <- data.frame(mean=c(p4,p4b), time = c(t,t), Group = c(rep("Variable", times = n),rep("Fixed", times = n)),
                     q0.025=c(p5,p5b),q0.975=c(p6,p6b))
graphB_ls <- ggplot(aes(x=time, y=mean), data=plotdf2) +
  stat_function(fun = function(x) {sin(-0.1 * pi * x)+0.3},
                color = "black", linetype = "longdash", size = 1.5) +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975, fill=Group), data = plotdf2[plotdf2$Group=="Variable",], alpha = 0.2) +
  geom_ribbon(aes(ymin=q0.025, ymax=q0.975, fill=Group), data = plotdf2[plotdf2$Group!="Variable",], alpha = 0.2) +
  #geom_ribbon(aes(ymin=q0.1, ymax=q0.9, fill=Group), alpha = 0.2) +
  geom_line(aes(color=Group)) +
  #geom_label(aes(x=label_x, y=label_y, label = label), data=label_df) +
  scale_colour_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  #scale_y_continuous(limits = c(y_lim_min-50,y_lim_max+50) #, expand = c(0, 0)
  #) +
  #scale_x_continuous(limits = c(min(data_1_clean$time),max(data_1_clean$time))
  #, expand = c(0, 0)
  #) +
  xlab("t") +
  ylab("Y") + #labs(title=paste0("Study 1")) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(title = "Response over targeted smoothing variable: X>0.5")
graphB_ls

{ ggsave(paste0("graphA_ls.pdf"),
         graphA_ls, width=6, height=4)
  ggsave(paste0("graphB_ls.pdf"),
         graphB_ls, width=6, height=4)}




