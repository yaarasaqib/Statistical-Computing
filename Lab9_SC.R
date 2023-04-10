    ######### Lab 9  #########
##### Class Demonstration #######
    

################################################
## Old Faithful Geyser data
################################################
data(faithful)
head(faithful)

x <- faithful$eruptions
hist(x, breaks = 30, main = "Eruptions")


################################################
## EM Algorithm for the Old Faithful Geyser data
################################################
## Estep gamma calculation
## this is a general function that we are not actuall
gamma_ick <- function(x, mu, sig2, pis, C = 2)
{
  gamma_prob <- matrix(0, nrow = length(x), ncol = C)
  for(c in 1:C)
  {
    gamma_prob[ ,c]  <- dnorm(x, mean = mu[c], sd = sqrt(sig2[c]))* pis[c]
  }
  
  gamma_prob <- gamma_prob/(rowSums(gamma_prob))
  return(gamma_prob)
}


# Starting values
pis <- c(.6, .4) 
mu <- c(1, 5)
sig2 <- c(1, 2)
diff <- 100
tol <- 1e-5
iter <- 0

# just for visualizing
current <- c(pis, mu, sig2)
store <- current
C <- 2
while(diff > tol)
{
  previous <- current
  iter <- iter + 1
  
  # E step: find gamma_{i,c,k} for just c = 1, since for c = 2 is just 1-Ep
  # Ep <- current[1]*dnorm(x, current[2], sqrt(current[4]))/
  #   (current[1]*dnorm(x, current[2], sqrt(current[4])) + (1 - current[1])*dnorm(x, current[3], sqrt(current[5])))
  # 
  Ep <- gamma_ick(x, mu, sig2, pis, C = 2)
  
  # M-step
  pis <- colMeans(Ep)
  mu <- colSums(Ep*x) / colSums(Ep)
  for(c in 1:C)
  {
    sig2[c] <- sum(Ep[,c]*(x - mu[c])^2) / sum(Ep[,c])
  }
  current <- c(pis, mu, sig2)
  
  diff <- norm(previous - current, "2")
  store <- rbind(store, current)
}

current # final estimates


# Final estimates of the probability
# that each observation is in Class C.
Prob.Z <- current[1]*dnorm(x, current[2], sqrt(current[4]))/
  (current[1]*dnorm(x, current[2], sqrt(current[4])) + (1 - current[1])*dnorm(x, current[3], sqrt(current[5])))

head(round(Prob.Z, 10))


# Make plot of iterative model fits
hist(x, breaks = 30, main = "Eruptions", freq = FALSE)
for(i in 1:dim(store)[1])
{
  test.x <- seq(min(x), max(x), length = 1000)
  test.y <- store[i,1]* dnorm(test.x, mean = store[i,3], sd = sqrt(store[i,5])) + (store[i,2]) *dnorm(test.x, mean = store[i,4], sd = sqrt(store[i,6]))
  lines(test.x, test.y, col = rgb(1,0,0, alpha = .5))
  i <- i + 1
}
lines(test.x, test.y, col = rgb(0,0,1, alpha = 1))

# add color
color <- 1*(Ep < .5) + 3*(Ep >= .5)
points(x, rep(.0, length(x)), pch = 16, col = color)

##########################
##  Problem 2 ###

## 1 dim data to  C Cluster  ###

GMMoneDim <- function(x, C = 2)
{
  # Starting values
  pis <- rep(1/C, C)
  mu <- mean(x) + rnorm(C)
  sig2 <- rep(var(x), C)
  
  diff <- 100
  tol <- 1e-5
  iter <- 0
  # just for visualizing
  current <- c(pis, mu, sig2)
  while(diff > tol)
  {
    previous <- current
    iter <- iter + 1
    # E step
    Ep <- gamma_ick(x, mu, sig2, pis, C = C)
    # M-step
    pis <- colMeans(Ep)
    mu <- colSums(Ep*x) / colSums(Ep)
    for(c in 1:C)
    {
      sig2[c] <- sum(Ep[,c]*(x - mu[c])^2) / sum(Ep[,c])
    }
    current <- c(pis, mu, sig2)
    diff <- norm(previous - current, "2")
  }
  Ep <- gamma_ick(x, mu = mu, sig2 = sig2, pis = pis, C = C)
  Zs <- apply(Ep, 1, which.max)
  rtn <- list(Zs, mu, sig2, pis)
  return(rtn)
}

###  Runninf for c= 2 and 3
set.seed(1)
erup2 <- GMMoneDim(faithful$eruptions, C = 2)

erup3 <- GMMoneDim(faithful$eruptions, C = 3)
n <- length(faithful$eruptions)
par(mfrow = c(1,2))
hist(faithful$eruptions, breaks = 30, main = "Two clusters", xlab = "x")
points(faithful$eruptions, rep(0,n), col = erup2[[1]], pch = 16)
hist(faithful$eruptions, breaks = 30, main = "Three clusters", xlab = "x")
points(faithful$eruptions, rep(0,n), col = erup3[[1]], pch = 16)


######################################
####  Problem 3 ###

## BIC For a given C=c

# First, we complete the function to calculate the log-likelihood.

log.like <- function(x, mu, sig2, pis)
{
  n <- length(x)
  track <- 0
  # calculate for each data point
  for(i in 1:n)
  {
    track <- track + log(sum(dnorm(x[i], mean = mu, sd = sqrt(sig2))*pis))
  }
  return(track)
}


# first finding the log-lielihood values for both models
x <- faithful$eruptions
loglike2 <- log.like(x, mu = erup2[[2]], sig2 = erup2[[3]], pis = erup2[[4]])
loglike3 <- log.like(x, mu = erup3[[2]], sig2 = erup3[[3]], pis = erup3[[4]])

# number of parameters
C <- 2
K2 <- 3*C - 1
C <- 3
K3 <- 3*C - 1
BIC.2 <- 2*loglike2 - K2*log(n)
BIC.3 <- 2*loglike3 - K3*log(n)

# comparing BICs
c(BIC.2, BIC.3)

#############################
##   Problem 4 ###

# Repeat for waiting Componenet of Data set

# Fitting the EM:
set.seed(1)
x <- faithful$waiting
erup2 <- GMMoneDim(x, C = 2)
erup3 <- GMMoneDim(x, C = 3)
n <- length(x)
par(mfrow = c(1,2))
hist(x, breaks = 30, main = "Two clusters", xlab = "x")
points(x, rep(0,n), col = erup2[[1]], pch = 16)
hist(x, breaks = 30, main = "Three clusters", xlab = "x")
points(x, rep(0,n), col = erup3[[1]], pch = 16)

loglike2 <- log.like(x, mu = erup2[[2]], sig2 = erup2[[3]], pis = erup2[[4]])
loglike3 <- log.like(x, mu = erup3[[2]], sig2 = erup3[[3]], pis = erup3[[4]])
# number of parameters
C <- 2
K2 <- 3*C - 1
C <- 3
K3 <- 3*C - 1
BIC.2 <- 2*loglike2 - K2*log(n)
BIC.3 <- 2*loglike3 - K3*log(n)

# comparing BICs
c(BIC.2, BIC.3)

## BIC2 is clearly better



###############################################

########## Problem 6   #############

# function that calculates the
# the log_likelihood
# for this multivariate setup
library(mvtnorm)
log_like <- function(X, pi.list, mu.list, Sigma.list, C)
{
  foo <- 0
  for(c in 1:C)
  {
    foo <- foo + pi.list[[c]]*dmvnorm(X, mean = mu.list[[c]], sigma = Sigma.list[[c]])
  }
  return(sum(log(foo)))
}
# Now I recommend the following
# mu is a list
# Sigma is a list
GMMforD <- function(X, C = 2, tol = 1e-5, maxit = 1e3)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  ######## Starting values ###################
  ## pi are equally split over C
  pi.list <- rep(1/C, C)
  mu <- list() # list of all the means
  Sigma <- list() # list of all variances
  # The means for each C cannot be the same,
  # since then the three distributions overlap
  # Hence adding random noise to colMeans(X)
  # the variance of the noise depends on the components
  for(i in 1:C)
  {
    mu[[i]] <- rnorm(p) + colMeans(X)
    Sigma[[i]] <- var(X) # same covariance matrix
  }
  # Choosing good starting values is important since
  # The GMM likelihood is not concave, so the algorithm
  # may converge to a local optima.
  ######## EM algorithm steps ###################
  iter <- 0
  diff <- 100
  epsilon <- 1e-05 # postive-def of Sigma
  Ep <- matrix(0, nrow = n, ncol = C) # gamma_{i,c}
  current <- c(unlist(mu), unlist(Sigma), pi.list)
  while((diff > tol) && (iter < maxit) )
  {
    iter <- iter + 1
    previous <- current
    ## E step: find gammas
    for(c in 1:C)
    {
      Ep[ ,c] <- pi.list[c]*apply(X, 1, dmvnorm , mu[[c]], Sigma[[c]])
    }
    Ep <- Ep/rowSums(Ep)
    ### M-step
    pi.list <- colMeans(Ep)
    for(c in 1:C)
    {
      mu[[c]] <- colSums(Ep[ ,c] * X )/sum(Ep[,c])
    }
    for(c in 1:C)
    {
      foo <- 0
      for(i in 1:n)
      {
        foo <- foo + (X[i, ] - mu[[c]]) %*% t(X[i, ] - mu[[c]]) * Ep[i,c]
      }
      Sigma[[c]] <- foo/sum(Ep[,c])
      if(min(eigen(Sigma[[c]])$values) <=0)
      {
        # To ensure the estimator is positive definite
        # otherwise next iteration gamma_i,c,k cannot be calculated
        Sigma[[c]] <- Sigma[[c]] + diag(epsilon, p)
        print("Matrix not positive-definite")
      }
    }
    # Difference in the log-likelihoods as the difference criterion
    current <- c(unlist(mu), unlist(Sigma), pi.list)
    diff <- norm(previous - current, "2")
  }
  # Final allocation updates
  for(c in 1:C)
  {
    Ep[ ,c] <- pi.list[c]*apply(X, 1, dmvnorm , mu[[c]], Sigma[[c]])
  }
  Ep <- Ep/rowSums(Ep)
  # calculate the loglikelihood of the final est
  save.loglike <- log_like(X = X, pi.list = pi.list, mu.list = mu, Sigma.list = Sigma, C = C)
  return(list("pi" = pi.list, "mu" = mu, "Sigma" = Sigma, "Ep" = Ep,
              "log.like" = save.loglike))
}
    
########## Problem 7  #######

# Run for Bivariate faithful data
X <- as.matrix(faithful)
C <- 2
class2.1 <- GMMforD(X = X, C = C)
class2.2 <- GMMforD(X = X, C = C)
class2.3 <- GMMforD(X = X, C = C)
class2.4 <- GMMforD(X = X, C = C)

library(ellipse)

# making 4 plots, one for each fit
par(mfrow = c(2,2))
for(i in 1:4){
  # get() allows to access the 4 different models
  model <- get(paste("class2.",i, sep = ""))
  # find which cluster is assigned to which data point
  allot <- apply(model$Ep, 1, which.max)
  # plot the data
  plot(faithful[,1], faithful[,2], col = allot, pch = 16,
       main = paste("Log-Like = ", round(model$log.like,3))) # plot allotment
  # make an ellipse for each cluster, to see the shape of clusters.
  ell <- list()
  for(c in 1:C)
  {
    ell[[c]] <- ellipse(model$Sigma[[c]], centre = as.numeric(model$mu[[c]]))
    lines(ell[[c]], col = c)
  }
}

### Problem 8

#Use BIC to obtain the best value of 𝐶 among 2, 3, 4 for this dataset.

Kpar <- function(C, p)
{
  (C-1) + C*p + p*(p+1)/2*C
}

class3.1 <- GMMforD(X = X, C = 3)
class3.2 <- GMMforD(X = X, C = 3)
class3.3 <- GMMforD(X = X, C = 3)
class3.4 <- GMMforD(X = X, C = 3)
class3.5 <- GMMforD(X = X, C = 3)
class3.6 <- GMMforD(X = X, C = 3)

#Plotting them:

  C <- 3
par(mfrow = c(2,3))
for(i in 1:6)
{
  model <- get(paste("class3.",i, sep = ""))
  allot <- apply(model$Ep, 1, which.max) ## Final allotment of classification
  plot(faithful[,1], faithful[,2], col = allot, pch = 16,
       main = paste("Log-Like = ", round(model$log.like,3))) # plot allotment
  ell <- list()
  for(c in 1:C)
  {
    ell[[c]] <- ellipse(model$Sigma[[c]], centre = as.numeric(model$mu[[c]]))
    lines(ell[[c]], col = c)
  }
}

### for c=4

class4.1 <- GMMforD(X = X, C = 4)
class4.2 <- GMMforD(X = X, C = 4)
class4.3 <- GMMforD(X = X, C = 4)
class4.4 <- GMMforD(X = X, C = 4)
class4.5 <- GMMforD(X = X, C = 4)
class4.6 <- GMMforD(X = X, C = 4)

## For Plot
C <- 4
par(mfrow = c(2,3))
for(i in 1:6)
{
  
  model <- get(paste("class4.",i, sep = ""))
  allot <- apply(model$Ep, 1, which.max) ## Final allotment of classification
  plot(faithful[,1], faithful[,2], col = allot, pch = 16,
       main = paste("Log-Like = ", round(model$log.like,3))) # plot allotment
  ell <- list()
  for(c in 1:C)
  {
    ell[[c]] <- ellipse(model$Sigma[[c]], centre = as.numeric(model$mu[[c]]))
    lines(ell[[c]], col = c)
  }
}

# To compare which one is better
BIC.2 <- 2*class2.1$log.like - log(n)*Kpar(C = 2, p = 2)
BIC.3 <- 2*class3.4$log.like - log(n)*Kpar(C = 3, p = 2)
BIC.4 <- 2*class4.5$log.like - log(n)*Kpar(C = 4, p = 2)
c(BIC.2, BIC.3, BIC.4)


