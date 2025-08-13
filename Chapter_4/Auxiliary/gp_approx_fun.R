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

# antiderivatve of the phi function above
Phi = function(x, j, L) {
  -2*sqrt(L)*cos(pi*j*(x+L)/(2*L))/(pi*j) #+C hahaha
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
