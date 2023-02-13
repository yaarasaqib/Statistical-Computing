########################
#### LAB 4  SC  #################
#######################

#### Qsn2 |Moments  of Gamma Distribution using Imp Sampling ####
##### part A ##############
k <- 2  # for second moment
N <- 1e3  # number of sample
alpha <- 2
beta <- 5

(truth <- alpha/beta^2+(alpha/beta)^2)   # True Second moment

#Simple monte carlo

samples <- rgamma(N,shape=alpha,rate=beta)
mean(samples^k)

## Importance Sampling
set.seed(1)
lambda<- 3 # proposal
N <- 1e4
samp <- rexp(N,rate=lambda)

# evaluate Inside the sum
funcs <- samp^k*dgamma(samp,shape=alpha,rate=beta)/dexp(samp,rate=lambda)

mean(funcs)
# sigmaˆ2_g
var(funcs)

########### Repeating for different lambda ##############
num <- 15
mean <- numeric(num)
var.g <- numeric(num)
for( c in 1:num){
set.seed(1)
lambda<- c # proposal
N <- 1e4
samp <- rexp(N,rate=lambda)

# evaluate Inside the sum
funcs <- samp^k*dgamma(samp,shape=alpha,rate=beta)/dexp(samp,rate=lambda)

mean[c]<-mean(funcs)

# sigmaˆ2_g
var.g[c]<-var(funcs)
}
mean
var.g   ## lambda=1 giving minimum var.g

### Qsn 2 | part B ####

##repeating for proposal gamma(3,3)
k <- 2  # for second moment
N <- 1e3  # number of sample
alpha <- 2
beta <- 5

samples <- rgamma(N,shape=alpha,rate=beta)
mean(samples^k)

## Importance Sampling
set.seed(1)
theta<- 3 # proposal
N <- 1e4
samp <- rgamma(N,shape=theta,rate=theta)

# evaluate Inside the sum
funcs <- samp^k*dgamma(samp,shape=alpha,rate=beta)/rgamma(N,shape=theta,rate=theta)
mean(funcs)
# sigmaˆ2_g
var(funcs)

##############################

##### Question 3 - Law of Large Number #####

# Checking the Convergence
N <- 1e5  # N is large
lambda <- 3
samp <- rexp(N,rate=lambda)    # importance sample
func <- samp^2*dgamma(samp,shape=alpha,rate=beta) /dexp(samp,rate=lambda)
x.axis<- 1:N   # sample size on the x-axis
y.axis <- cumsum(func)/(1:N)  # IS estimator for each N

# Plotting the running average
plot(x.axis, y.axis , type = 'l', xlab = "N", ylab = "Running average")
abline(h = truth, col = "red")

###################################
###### Question 4  ################

# Estimate MGF using importance sampling

# choosing 50 values of t between (-5, 5)
t <- seq(-5, 5, length = 50)

mu <- 1
lambda <- 3

#, using importance sampling with importance distribution
#Gamma(10, 3)    proposal

h <- function(t,x){exp(t*x)}
f <- function(x){sqrt(lambda/(2*pi*x^3)*exp((-lambda*(x-mu)^2))/(2*mu^2*x))}
##f target 

N <-1e5
est <- numeric(length(t))
for(i in t){
  samp <- rgamma(N,shape=10,rate=3)  # IS
  s <- h(i,samp)*(f(samp)/dgamma(samp,shape=10,rate=3))
  est[which(t==i)]= mean(s)
  }
plot(t,est,type="l",xlab="t",ylab=expression(M_x(t)),main="Plot of mgf")
#############################################
## ma,am solution

par(mfrow=c(1,1))
# Inv Gaussian density function
dinvGaus <- function(x, lam, mu)
{
  # assuming x > 0 always
  den <- sqrt(lam/(2 * pi * x^3)) * exp(-lam * (x-mu)^2/(2 * mu^2 * x))
  return(den)
}

momInvG <- function(t, lam = 3, mu = 1, N = 1e4)
{
  samp <- rgamma(N, shape = 10, rate = 3)
  #using the same sample for each t
  ests <- numeric(length = length(t))
  for(k in 1:length(t))
  {
    foo <- exp(t[k]*samp) * dinvGaus(samp, lam, mu)/dgamma(samp, shape = 10, rate = 3)
    ests[k] <- mean(foo)
  }
  return(ests)
}
# the original grid won't work
# using this grid
t <- seq(-4, .6, length = 50)
lam <- 3
mu <- 1
ests <- momInvG(t = t)
truth <- exp(lam/mu* (1 - sqrt(1 - 2*mu^2 *t/lam)) )
plot(t, ests, type = "l")
lines(t, truth, col = "red")
legend("topleft", legend = c("estimated", "true"), col = c("black", "red"), lty = 1)

### different samples foe every t

momInvG <- function(t, lam = 3, mu = 1, N = 1e4)
{
  #using the same sample for each t
  ests <- numeric(length = length(t))
  for(k in 1:length(t))
  {
    samp <- rgamma(N, shape = 10, rate = 3)
    foo <- exp(t[k]*samp) * dinvGaus(samp, lam, mu)/dgamma(samp, shape = 10, rate = 3)
    ests[k] <- mean(foo)
  }
  return(ests)
}
# the original grid won't work
# using this grid
t <- seq(-4, .6, length = 50)
lam <- 3
mu <- 1
ests <- momInvG(t = t)
truth <- exp(lam/mu* (1 - sqrt(1 - 2*mu^2 *t/lam)) )
plot(t, ests, type = "l")
lines(t, truth, col = "red")
legend("topleft", legend = c("estimated", "true"), col = c("black", "red"), lty = 1)


