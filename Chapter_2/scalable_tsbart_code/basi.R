phi <- function(x,j,L){
  sin(pi*j*(x+L)/(2*L))
}

sqexp <- function(x,a,l){
  a*exp(-0.5*(x/l)^2)
}

s <- function(i,l,c){(l/sqrt(c))*exp(-0.5*l^2*(i*pi/(2*c))^2)}

x <- seq(0,10, length.out = 1000)

l <- 0.1
aa <- s(x,l,1)
l <- 0.25
bb <- s(x,l,1)
l <- 0.5
cc <- s(x,l,1)

L <- 1
sqrt((pi*1)/(2*L))

pdf("s_graph.pdf")
plot(x,aa, type = "l",lwd=2, xlab = "q", ylab = expression(delta[theta](q)), xlim = c(0,10), ylim = c(0,0.5))
lines(x,bb, col = "red",lwd=2)
lines(x,cc, col = "blue",lwd=2)
legend("topright", legend = c(expression(paste(theta,"=0.10")),
                              expression(paste(theta,"=0.25")),
                              expression(paste(theta,"=0.50"))),
       fill=c("black","red", "blue"))
dev.off()

x <- seq(0,1, length.out = 1000)
l <- 0.1
aa <- sqexp(x,1,l)
l <- 0.25
bb <- sqexp(x,1,l)
l <- 0.5
cc <- sqexp(x,1,l)

L <- 1

pdf("s2_graph.pdf")
plot(x,aa, type = "l",lwd=2, xlab = "t", ylab = expression(C[theta](t)), xlim = c(0,1), ylim = c(0,1))
lines(x,bb, col = "red",lwd=2)
lines(x,cc, col = "blue",lwd=2)
legend("topright", legend = c(expression(paste(theta,"=0.10")),
                              expression(paste(theta,"=0.25")),
                              expression(paste(theta,"=0.50"))),
       fill=c("black","red", "blue"))
dev.off()

(L^(1/2))*sin(pi*j*(x+L)/(2*L))


sqrt(1/L)*sin((i*pi/(2*L))*(x+L))*sqrt(1/L)*sin((i*pi/(2*L))*(0+L))*(a*sqrt(2*pi)*l*exp(-0.5*l^2*((i*pi/(2*L))^2)))

i <- 1
L <- 0.5

x <- seq(-1,1, length.out = 1000)
plot(x,phi(x,1,1), type = "l", lwd = 2)
lines(x,phi(x,2,1), col = 2, lwd = 2)
lines(x,phi(x,3,1), col = 3, lwd = 2)
lines(x,phi(x,4,1), col = 4, lwd = 2)
lines(x,phi(x,5,1), col = 5, lwd = 2)

pdf("basis_curves.pdf")
par(mfrow = c(2,2))
plot(x,phi(x,1,1), type = "l", lwd = 2, col = 1,
     ylab = expression(phi[1](t)), xlab = "t")
plot(x,phi(x,2,1), type = "l", lwd = 2, col = 2,
     ylab = expression(phi[2](t)), xlab = "t")
plot(x,phi(x,3,1), type = "l", lwd = 2, col = 3,
     ylab = expression(phi[3](t)), xlab = "t")
plot(x,phi(x,4,1), type = "l", lwd = 2, col = 4,
     ylab = expression(phi[4](t)), xlab = "t")
dev.off()
