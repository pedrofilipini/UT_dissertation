Omega <- function(x,j,L){
  phi <- function(x,j,L){
    sqrt(1/L)*sin(pi*j*(x+L)/(2*L))
  }
  omega <- matrix(0, ncol = j, nrow = length(x))
  for(i in 1:j){
    for(k in 1:length(x)){
      omega[k,i] <- phi(x[k], i, L)
    }
  }
  return(omega)
}


#Select m. First m to cross e, then c that minimizes that error
m_select <- function(a=1, mmin=5, mmax=100,cmin=1, cmax=10, step = 0.01, softmin,xmin,xmax,e=0.01){
  omega <- function(x,a,l,L,m){
    aux <- 0
    for(i in 1:m){
      aux <- aux + sqrt(1/L)*sin((i*pi/(2*L))*(x+L))*sqrt(1/L)*sin((i*pi/(2*L))*(0+L))*(a*sqrt(2*pi)*l*exp(-0.5*l^2*((i*pi/(2*L))^2)))
    }
    return(aux)
  }

  sqexp <- function(x,a,l){
    a*exp(-0.5*(x/l)^2)
  }

  integrand <- function(x,a,l,L,m){
    return(abs(sqexp(x,a,l)-omega(x,a,l,L,m)))
  }

  #m to be tested
  m <- mmin:mmax
  #c to be tested
  c <- seq(cmin,cmax, by = step)
  #The softmin
  l <- softmin

  #Lasterror to pick the minimum
  lasterror <- Inf
  for(mm in 1:length(m)){

    for(cc in 1:length(c)){

      #Update L
      L <- ((xmax-xmin)/2)*c[cc]
      #Integrate
      aux <- integrate(integrand, lower = 0, upper = (xmax-xmin),
                       a=a,l=l,L=L,m=m[mm])$value/integrate(sqexp, lower = 0, upper = (xmax-xmin), a=a,l=l)$value

      #If it is less than 1% and minimum, return
      if((aux<e)&&(aux>lasterror)){
        return(c(m[mm],c[(cc-1)]))
      }
      lasterror <- aux
    }
  }
}

#Return c for the minimum error given m and l
c_select <- function(a=1, m, cmin=1, cmax=10, step = 0.01, l, xmin, xmax){

  #The approximation
  omega <- function(x,a,l,L,m){
    aux <- 0
    for(i in 1:m){
      aux <- aux + sqrt(1/L)*sin((i*pi/(2*L))*(x+L))*sqrt(1/L)*sin((i*pi/(2*L))*(0+L))*(a*sqrt(2*pi)*l*exp(-0.5*l^2*((i*pi/(2*L))^2)))
    }
    return(aux)
  }

  #The true covariance function
  sqexp <- function(x,a,l){
    a*exp(-0.5*(x/l)^2)
  }

  #The absolute difference
  integrand <- function(x,a,l,L,m){
    return(abs(sqexp(x,a,l)-omega(x,a,l,L,m)))
  }

  #Taking the vector of possible c for testing
  c <- seq(cmin,cmax, by = step)

  #Keeps track of the lowest error
  lasterror <- Inf

  #The loop to find c
  for(cc in 1:length(c)){
    #Half range times c
    L <- ((xmax-xmin)/2)*c[cc]

    #Calculate the integral
    aux <- integrate(integrand, lower = 0, upper = (xmax-xmin),
                     a=a,l=l,L=L,m=m)$value/integrate(sqexp, lower = 0, upper = (xmax-xmin), a=a,l=l)$value

    #If the error starts increasing, stop
    if(aux>lasterror){
      break;
    }
    #Else, update lasterror with current error
    lasterror <- aux

  }
  #Return the previous value of c and its error
  return(c(c[(cc-1)],lasterror))
}

#Select the l that would be required for minimizing the error for a given m and c
lmax_select <- function(a=1, m=10, c=10, step = 0.01, lmin, xmin=0, xmax=1){

  #The approximation
  omega <- function(x,a,l,L,m){
    aux <- 0
    for(i in 1:m){
      aux <- aux + sqrt(1/L)*sin((i*pi/(2*L))*(x+L))*sqrt(1/L)*sin((i*pi/(2*L))*(0+L))*(a*sqrt(2*pi)*l*exp(-0.5*l^2*((i*pi/(2*L))^2)))
    }
    return(aux)
  }

  #The true covariance function
  sqexp <- function(x,a,l){
    a*exp(-0.5*(x/l)^2)
  }

  #The absolute difference
  integrand <- function(x,a,l,L,m){
    return(abs(sqexp(x,a,l)-omega(x,a,l,L,m)))
  }

  #Taking the vector of possible c for testing
  l <- lmin

  #Keeps track of the lowest error
  lasterror <- Inf
  L <- ((xmax-xmin)/2)*c
  #The loop to find c
  ll=1
  while(T){
    #Calculate the integral
    aux <- integrate(integrand, lower = 0, upper = (xmax-xmin),
                     a=a,l=l[ll],L=L,m=m)$value/integrate(sqexp, lower = 0, upper = (xmax-xmin), a=a,l=l[ll])$value

    ll <- ll+1
    l[ll] <- l[(ll-1)]+step

    #If the error starts increasing, stop
    if(aux>lasterror){
      #print(paste(c(l[(ll-1)],lasterror)))
      if(c_select(a=1, m=m, cmin=1, cmax=c*2, step = step, l=l[(ll-1)], xmin=xmin,xmax=xmax)[1]>=c){
        break;
      }
    }
    #Else, update lasterror with current error
    lasterror <- aux

  }
  #Return the previous value of c and its error
  return(c(l[(ll-1)],lasterror))
}

graph_fun <- function(a, mmax,cmax,lmax, xmax, e){
  omega_2 <- function(x,a,l,L,m){
    aux <- 0
    for(i in 1:m){
      aux <- aux + sqrt(1/L)*sin((i*pi/(2*L))*(x+L))*sqrt(1/L)*sin((i*pi/(2*L))*(0+L))*(a*sqrt(2*pi)*l*exp(-0.5*l^2*((i*pi/(2*L))^2)))
    }
    return(aux)
  }

  sqexp <- function(x,a,l){
    a*exp(-0.5*(x/l)^2)
  }

  integrand2 <- function(x,a,l,L,m){
    return(abs(sqexp(x,a,l)-omega_2(x,a,l,L,m)))
  }


  m <- 5:mmax
  c <- seq(1,cmax, by = 1)
  ls <-seq(0.05,lmax, by = 0.05)
  xmin <- 0
  l <- ls*xmax/2

  for(cc in 1:length(c)){

    L <- (xmax/2)*c[cc]

    result <- NULL
    m_result <- NULL
    l_result <- NULL
    c_result <- NULL

    for(ll in 1:length(l)){
      for(mm in 1:length(m)){
        aux <- integrate(integrand2, lower = xmin, upper = xmax,
                         a=a,l=l[ll],L=L,m=m[mm])$value/integrate(sqexp, lower = xmin, upper = xmax, a=a,l=l[ll])$value
        if(aux<e){
          result <- c(result, aux)
          m_result <- c(m_result, m[mm])
          c_result <- c(c_result, c[cc])
          l_result <- c(l_result, l[ll])
          break
        }
      }
    }
    if(cc==1){
      plot(l_result*2/xmax, m_result, type = 'l', col = cc, lwd = 2, ylim = c(0,mmax), xlim = c(0,lmax),
           ylab = "m", xlab = "2*lengthscale/Xmax", main = paste0("Error = ",e))
    }else{
      lines(l_result*2/xmax, m_result, col = cc, lwd = 2)
    }

  }
}

par(mfrow = c(1,1))
pdf("cand_1.pdf")
graph_fun(a = 1, mmax = 100,cmax = 6,lmax = 3,xmax = 1, e = 0.01)
dev.off()


#Fix m and l, plot error as f(L)
c_select_2 <- function(a=1, m, cmin=1, cmax=10, step = 0.01, l, xmin, xmax){

  #The approximation
  omega <- function(x,a,l,L,m){
    aux <- 0
    for(i in 1:m){
      aux <- aux + sqrt(1/L)*sin((i*pi/(2*L))*(x+L))*sqrt(1/L)*sin((i*pi/(2*L))*(0+L))*(a*sqrt(2*pi)*l*exp(-0.5*l^2*((i*pi/(2*L))^2)))
    }
    return(aux)
  }

  #The true covariance function
  sqexp <- function(x,a,l){
    a*exp(-0.5*(x/l)^2)
  }

  #The absolute difference
  integrand <- function(x,a,l,L,m){
    return(abs(sqexp(x,a,l)-omega(x,a,l,L,m)))
  }

  #Taking the vector of possible c for testing
  c <- seq(cmin,cmax, by = step)

  #Keeps track of the lowest error
  lasterror <- Inf

  error <- NULL
  #The loop to find c
  for(cc in 1:length(c)){
    #Half range times c
    L <- ((xmax-xmin)/2)*c[cc]

    #Calculate the integral
    aux <- integrate(integrand, lower = 0, upper = (xmax-xmin),
                     a=a,l=l,L=L,m=m)$value/integrate(sqexp, lower = 0, upper = (xmax-xmin), a=a,l=l)$value
    error <- c(error, aux)
    #If the error starts increasing, stop
    #if(aux>lasterror){
    #  break;
    #}
    #Else, update lasterror with current error
    #lasterror <- aux

  }
  #Return the previous value of c and its error
  return(data.frame(c=c,e=error))
}

#test <- c_select(a=1,m=10,l=0.1, xmin = 0, xmax = 1)
#test2 <- c_select(a=1,m=20,l=0.1, xmin = 0, xmax = 1)
#test3 <- c_select(a=1,m=50,l=0.1, xmin = 0, xmax = 1)

test <- list()
count <- 0
for(m in 5:25){
  count <- count + 1
  test[[count]] <- c_select_2(a=1,m=m,l=0.1, xmin = 0, xmax = 1)
}

pdf("cand_21.pdf")
plot(test[[1]]$c, test[[1]]$e, ylim = c(0,1.5),
     type = "l", lwd = 2, xlab = "c", ylab = "e",
     main = expression(paste(theta,"=0.1")))
abline(v=test[[1]]$c[which.min(test[[1]]$e)], lty = 3, lwd = 2, col = "orange")
lines(test[[3]]$c, test[[3]]$e, col = 2, lwd = 2)
lines(test[[5]]$c, test[[5]]$e, col = 3, lwd = 2)
lines(test[[7]]$c, test[[7]]$e, col = 4, lwd = 2)
lines(test[[9]]$c, test[[9]]$e, col = 5, lwd = 2)
lines(test[[11]]$c, test[[11]]$e, col = 6, lwd = 2)
lines(test[[13]]$c, test[[13]]$e, col = 7, lwd = 2)
lines(test[[15]]$c, test[[15]]$e, col = 8, lwd = 2)
#abline(h = 0.05, lty = 2, lwd = 3)
legend("bottomright", legend = c("m=5","m=7","m=9","m=11","m=13","m=15","m=17","m=19"),
       fill=c(1:8), inset = 0.08)
dev.off()


pdf("cand_22.pdf")
plot(test[[1]]$c, test[[1]]$e, ylim = c(0,1.5),
     type = "l", lwd = 2, xlab = "c", ylab = "e",
     main = expression(paste(theta,"=0.1")))
lines(test[[3]]$c, test[[3]]$e, col = 2, lwd = 2)
lines(test[[5]]$c, test[[5]]$e, col = 3, lwd = 2)
lines(test[[7]]$c, test[[7]]$e, col = 4, lwd = 2)
lines(test[[9]]$c, test[[9]]$e, col = 5, lwd = 2)
lines(test[[11]]$c, test[[11]]$e, col = 6, lwd = 2)
lines(test[[13]]$c, test[[13]]$e, col = 7, lwd = 2)
lines(test[[15]]$c, test[[15]]$e, col = 8, lwd = 2)
abline(h = 0.05, lty = 2, lwd = 3)
legend("bottomright", legend = c("m=5","m=7","m=9","m=11","m=13","m=15","m=17","m=19"),
       fill=c(1:8), inset = 0.08)
dev.off()

test2 <- list()
count <- 0
for(m in 5:25){
  count <- count + 1
  test2[[count]] <- c_select_2(a=1,m=m,l=0.25, xmin = 0, xmax = 1)
}

pdf("cand_3.pdf")
plot(test2[[1]]$c, test2[[1]]$e, ylim = c(0,1.6),
     type = "l", lwd = 2, xlab = "c", ylab = "e",
     main = expression(paste(theta,"=0.25")))
abline(v=test[[1]]$c[which.min(test[[1]]$e)], lty = 3, lwd = 2, col = "orange")
abline(v=test2[[1]]$c[which.min(test2[[1]]$e)], lty = 3, lwd = 2, col = "red")
lines(test2[[3]]$c, test2[[3]]$e, col = 2, lwd = 2)
lines(test2[[5]]$c, test2[[5]]$e, col = 3, lwd = 2)
lines(test2[[7]]$c, test2[[7]]$e, col = 4, lwd = 2)
lines(test2[[9]]$c, test2[[9]]$e, col = 5, lwd = 2)
lines(test2[[11]]$c, test2[[11]]$e, col = 6, lwd = 2)
lines(test2[[13]]$c, test2[[13]]$e, col = 7, lwd = 2)
lines(test2[[15]]$c, test2[[15]]$e, col = 8, lwd = 2)
legend("topleft", legend = c("m=5","m=7","m=9","m=11","m=13","m=15","m=17","m=19"),
       fill=c(1:8))
#abline(h = 0.05, lty = 2)
dev.off()

pdf("cand_31.pdf")
plot(test[[3]]$c, test[[3]]$e, ylim = c(0,1.6),
     type = "l", lwd = 2, xlab = "c", ylab = "e",
     main = "m=9")
lines(test2[[3]]$c, test2[[3]]$e, col = 2, lwd = 2)
legend("topleft", legend = c(expression(paste(theta,"=0.1")),
                             expression(paste(theta,"=0.25"))),
       fill=c(1:2))
#abline(h = 0.05, lty = 2)
dev.off()

test3 <- list()
count <- 0
for(m in 5:25){
  count <- count + 1
  test3[[count]] <- c_select_2(a=1,m=m,l=0.01, xmin = 0, xmax = 1)
}

pdf("cand_4.pdf")
plot(test3[[1]]$c, test3[[1]]$e, ylim = c(0,4),
     type = "l", lwd = 2, xlab = "c", ylab = "Error",
     main = expression(paste("Changing m, fixed",theta,"=0.01")))
lines(test3[[3]]$c, test3[[3]]$e, col = 2, lwd = 2)
lines(test3[[5]]$c, test3[[5]]$e, col = 3, lwd = 2)
lines(test3[[7]]$c, test3[[7]]$e, col = 4, lwd = 2)
lines(test3[[9]]$c, test3[[9]]$e, col = 5, lwd = 2)
lines(test3[[11]]$c, test3[[11]]$e, col = 6, lwd = 2)
lines(test3[[13]]$c, test3[[13]]$e, col = 7, lwd = 2)
lines(test3[[15]]$c, test3[[15]]$e, col = 8, lwd = 2)
abline(h = 0.05, lty = 2)
dev.off()


pdf("cand_5.pdf")
l <- seq(0.01, 1, length.out = 100)
c <- numeric(0)
for(ll_con in 1:length(l)){
  c[ll_con] <- c_select(a=1, m=9, cmin=1, cmax=10, step = 0.01, l=l[ll_con], xmin=0, xmax=1)[1]
  #print(paste0("C_con: ", c_con_aux[ll_con], "; l: ", l_con_aux[ll_con]))
}
pdf("cand_51.pdf")
plot(l, c, ylab = "c", xlab = expression(theta),
     type = "l", lwd = 2, col = 1)
dev.off()
count <- 1
for(m in c(7,9,11,13,15,17)){
  count <- count+1
  l <- seq(0.01, 1, length.out = 100)
  c <- numeric(0)
  for(ll_con in 1:length(l)){
    c[ll_con] <- c_select(a=1, m=m, cmin=1, cmax=10, step = 0.01, l=l[ll_con], xmin=0, xmax=1)[1]
    #print(paste0("C_con: ", c_con_aux[ll_con], "; l: ", l_con_aux[ll_con]))
  }
  lines(l, c, lwd = 2, col = count)
}
dev.off()

#Error curve, fixed L, fixed m, changing l in the curve

#trying to create a story
#Why picking l is important? Show flat case, show wiggly case


library(extraDistr)

pdf("hnorm.pdf")
curve(dhnorm(x), 0, 5, ylab = "Density",
      main = "HNorm(0,1)", lwd = 2, col = "blue")
dev.off()








