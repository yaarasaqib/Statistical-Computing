####### Standard Cauchy ROU #####

a <- 1/sqrt(pi)
b <- (-1/sqrt(2*pi))
c <- (1/sqrt(2*pi))

draw_rect <- function(a,b,c){
  U1 <- a*runif(1)
  U2 <- b+(c-b)*runif(1)
  return(c(U1,U2))
}

sqrt.C <- function(x){
  y <- sqrt(1/(pi*(1+x^2)))
  return(y)
}

N <- 1e3  # number of sample required
counter <- 0
n<- 0
samp <- numeric(N)
while(n<N){
  prop <- draw_rect(a,b,c)
  u <- prop[1]
  v <- prop[2]
  if(u <= sqrt.C(v/u)){
    n<- n+1
    samp[n]<- v/u
    }
}

samp[1:5]

plot(density(samp),main="Estimated density for Exp(1)")
lines(density(rcauchy(1e4,location = 0,scale=1)), col = "red")
legend("topright", col = c("black", "red"), lty = 1, legend = c("RoU", "Truth"))

################################
### ROU of f(x)=1/x^2  x>0 #####
################################
a <- 1
b <- 0
c <- 1

draw_rect <- function(a,b,c){
  U1 <- a*runif(1)
  U2 <- b+(c-b)*runif(1)
  return(c(U1,U2))
}

sqrt.f <- function(x){1/x
}

N <- 1e3  # number of sample required
counter <- 0
n<- 0
samp <- numeric(N)
while(n<N){
  prop <- draw_rect(a,b,c)
  u <- prop[1]
  v <- prop[2]
  if(u <= sqrt.f(v/u)){
    n<- n+1
    samp[n]<- v/u
  }
}

samp[1:5]

x<- seq(0.00001,5,length=1e3)
f <- 1/x^2
plot(density(samp),main="Estimated density for 1/x^2",xlim=c(0,5))
lines(x,f, col = "red",xlim=c(0,5))
legend("topright", col = c("black", "red"), lty = 1, legend = c("RoU", "Truth"))



