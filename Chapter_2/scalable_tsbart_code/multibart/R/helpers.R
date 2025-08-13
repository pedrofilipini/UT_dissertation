.ident <- function(...){
  # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
  args <- c(...)
  if( length( args ) > 2L ){
    #  recursively call ident()
    out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
  }else{
    out <- identical( args[1] , args[2] )
  }
  return( all( out ) )
}

.cp_quantile = function(x, num=10000, cat_levels=8){
  nobs = length(x)
  nuniq = length(unique(x))

  if(nuniq==1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/nobs))
    ind = seq(min(x),max(x),length.out=num)
    ret = q(ind)
  }

  return(ret)
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

SDiag <- function(sigma, l, j, L){
  eigen <- function(j,L){
    (pi*j/(2*L))^2
  }
  s <- function(sigma, l, w){
    (sigma^2)*sqrt(2*pi)*l*exp(-0.5*(l^2)*(w^2))
  }
  return(diag(s(sigma,l,sqrt(eigen(1:j,L)))))
}

