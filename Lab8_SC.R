########### Lab 8 Solution #######

# Class Demontration ###

################################################
## MLE for location Cauchy distribution using
## MM-algorithm
################################################
set.seed(1)
mu.star <- 5  # Setting true mu
n <- 4  # sample size
X <- rt(n, df = 1) + mu.star

## Function calculates  the log-likelihood
log.like <- function(mu, X)
{
  n <- length(X)
  rtn <- -n*log(pi) - sum( log(1 + (X - mu)^2) )
  return(rtn)
}

mu.x <- seq(-10, 40, length = 1e3)  # A sequence of mu's 
ll.est <- sapply(mu.x, log.like, X)  # evaluating log-likelihood at the mus
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu))  # plotting log-likelihood. Not concave, so we need to choose good starting values.


## Starting Newton-Raphson method
tol <- 1e-5  # tolerance level for when to stop algorithm


## Returns derivate of log-likelihood
f.prime <- function(X, mu)
{
  rtn <- sum(2* (X - mu)/(1 + (X-mu)^2))  #f.prime
  return(rtn)
}

## Loop below stops when |mu_(k+1) - mu_(k)| < tol

current <- -10 # Good starting value
diff <- 100  # inital large value for difference
iter <- 0    # counting the number of iterations

mu.k <- current
while( (diff > tol) && iter < 1000)
{
  iter <- iter + 1
  update <- current - f.prime(X, current)/(-2*n)  # NR update
  mu.k <- c(mu.k, update)
  diff <- abs(current - update)
  current <- update
}
current  # final approximation to MLE
evals <- sapply(mu.k, log.like, X)
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu))
points(mu.k, evals, pch = 16, col = rgb(0,0,1, alpha = .5))
abline(v = current, lty = 2)

############################################

## problem 2 

# Bridge regression estimate using MM Algorithm

# y = response
# X = covariate matrix
# lambda = penlaty term
# alpha = bridge term
# max.ter = maximum iterations for the MM algorithm. if max iteration has been reached, the function should print: "Maximum iterations reached"
# tol = tolerance level for when to stop MM

bridgeReg <- function(y, X, alpha, lambda, max.iter, tol)
{
  # a value larger than tol
  distance <- tol + 1
  iter <- 0
  p <- dim(X)[2]
  # starting from the ridge solution
  current <- matrix(1, ncol = 1, nrow = p)
  #(t(X) %*% X + lambda*diag(p)) %*% t(X) %*%y
  while(distance > tol){
    iter <- iter + 1
    if(iter > max.iter)
    {
      print("Maximum iterations reached")
      stop
    }
    # MM steps
    previous <- current
    mjs <- alpha* abs(current)^(alpha - 2)
    ## using qr.solve since that is more stable than solve
    current <- qr.solve(t(X)%*%X + lambda/alpha * diag(mjs) ) %*% t(X) %*% y
    distance <- norm(previous - current, "2")
  }
  #returning the last iterate of the
  return(current)
}

## Problem 3
set.seed(1)
n <- 100
p <- 5
beta.star <- c(3, 1, .01, -2, -.003)
# X matrix with intercept
X <- cbind(1, matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1)))
# generating response
y <- X %*% beta.star + rnorm(n)

alpha.vec <- seq(1, 2, length = 5)
lambda <- c(.01, .1, 1, 10, 100)
for(i in 1:length(alpha.vec))
{
  for(j in 1:length(lambda))
  {
    bridge.est <- bridgeReg(y, X, alpha = alpha.vec[i], lambda = lambda[j], max.iter = 100, tol = 1e-5)
    print(bridge.est)
  }
}
 #####################

## Problem 4

# Expected value of ||beta_{k,lamda}-beta||

alpha.vec <- seq(1, 2, length = 5)
reps <- 100
dist <- matrix(0, nrow = reps, ncol = length(alpha.vec))
for(r in 1:reps)
{
  # generate X and y again here
  # X matrix with intercept
  X <- cbind(1, matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1)))
  # generating response
  y <- X %*% beta.star + rnorm(n)
  for(i in 1:length(alpha.vec))
  {
    bridge.est <- bridgeReg(y, X, alpha = alpha.vec[i], lambda = 10, max.iter = 100, tol = 1e-5)
    dist[r, i] <- norm(bridge.est - beta.star, "2")
  }
}

alpha.vec <- seq(1, 2, length = 5)
reps <- 100
dist <- matrix(0, nrow = reps, ncol = length(alpha.vec))
for(r in 1:reps)
{
# generate X and y again here
# X matrix with intercept
X <- cbind(1, matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1)))
# generating response
y <- X %*% beta.star + rnorm(n)
for(i in 1:length(alpha.vec))
{
bridge.est <- bridgeReg(y, X, alpha = alpha.vec[i], lambda = 10, max.iter = 100, tol = 1e-5)
dist[r, i] <- norm(bridge.est - beta.star, "2")
}
}

# naming the columns of dist
colnames(dist) <- alpha.vec
# expectation
colMeans(dist)

boxplot(dist, xlab = expression(alpha))
