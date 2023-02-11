########## LAB_3, Statistical Computing#########

### Question 1 ######
##################################################
### Ratio of Uniforms for Exp(1)
##################################################
set.seed(1)
# function to sample from the rectangle
drawFromRect <- function(a, b, c)
{
  u <- runif(1, min = 0, max = a)
  v <- runif(1, min = b, max = c)
  return(c(u,v))
}

# sqrt f function
sqrt.f <- function(x) {
  return(exp(-x/2))
  }

# Starting the process for Exp(1)
  a <-   1#fill from notes
  b <-  0 #fill from notes
  c <- 2/exp(1) #fill from notes
  prob.of.acceptance <- 1/(2*a*(c-b))  # true prob. of acceptance for AR

N <- 1e4 # number of samples
samp <- numeric(length = N)
i <- 1
counter <- 0  # to check acceptance
while(i <= N)
{
  counter <- counter + 1
  prop <- drawFromRect(a = a, b = b, c = c)
  vbyu <- prop[2]/prop[1]
  if( prop[1] <  sqrt.f( vbyu) )
  {
    samp[i] <- vbyu
      i <- i + 1
  }
}


plot(density(samp), main = "Estimated density for Exp(1)")
lines(density(rexp(1e4, 1)), col = "red")
legend("topright", col = c("black", "red"), lty = 1, legend = c("RoU", "Truth"))


(prob.of.acceptance)
# [1] 0.6795705

N/counter  # very close
# [1] 0.6796248

##################################################
 
### Question 2 ##########


##############################################
#### RoU region for exponetial
##############################################

# accept-reject from C
N <- 1e3
samples <- matrix(0, nrow = N, ncol = 2)
n <- 0
while(n < N)
{
  prop.u <- runif(1, min = 0, max = 1)
  prop.v <- runif(1, min = 0, max = 2/exp(1))
  if(prop.v <= -2*prop.u * log(prop.u))
  {
    n <- n+1
    samples[n, ] <- c(prop.u,prop.v)
  }
}

# define color based on regions
color <- 0
for(i in 1:10)
{
  color <- color + i*(samples[,2] > i*samples[,1] & samples[,2] < (i+1)*samples[,1])
}



# to plot the regions
u <- seq(0, 1, length = 1e4)  # b = 1
v <- -2*u * log(u)   #v = -2*u log(u)

# plot C with colors
par(mfrow = c(1,2))
plot(u, v, type = 'l', main = "C region for Exp")
abline(h = c(0, .8), v = c(0,1))
points(samples, col = color + 1, pch = 16)

x <- samples[,2]/samples[,1] # samples from Exp(1)
y <- dexp(x)
plot(x, y, col = color + 1, pch = 16, main = "Exponential samples with density")


##############################################
#### RoU region for N(theta, s2)
##############################################

theta <- 5
s2 <- 1

a <- 1/(2*pi*s2)^(.25)
foo <- (theta - sqrt(theta^2 + 8*s2))/(2*s2)
b <- foo * sqrt(dnorm(foo, theta, sqrt(s2)))

foo <- (theta + sqrt(theta^2 + 8*s2))/(2*s2)
c <- foo * sqrt(dnorm(foo, theta, sqrt(s2)))

N <- 5e3
samples <- matrix(0, nrow = N, ncol = 2)
n <- 0
while(n < N)
{
  prop.u <- runif(1, min = 0, max = a)
  prop.v <- runif(1, min = b, max = c)
  if(abs(prop.v - theta*prop.u) <= sqrt(-4* s2* prop.u^2 *(log(prop.u) + log(2*pi*s2)/4) ) )
  {
    n <- n+1
    samples[n, ] <- c(prop.u,prop.v)
  }
}

x <- samples[,2]/samples[,1]
z <- (x - theta)/sqrt(s2)
# define color based on regions
color <- 0
for(i in 0:3)
{
  color <- color + (i+1)*(z > i & z < (i+1) )
}
for(i in (-1:-3) )
{
  color <- color + (i+8)*(z > (i) & z < (i+1))
}


# to plot the regions
u <- seq(0.00001, (2*pi*s2)^(-.25), length = 1e3)
foo <-  sqrt(-4*s2*u^2 *(log(u) + log(2*pi*s2)/4))  

par(mfrow = c(1,2))
plot(u, theta * u + foo, type = 'l', main = "C region for N(theta,s2)", ylim = range(c(theta * u + foo, theta * u - foo)))
lines(u,theta * u - foo )

points(samples, col = color + 1, pch = 16)

y <- dnorm(x, mean = theta, sd = sqrt(s2))
plot(x, y, col = color + 1, pch = 16, main = "Normal samples with density")


####################################################
############### Question 3##########

####################################################################
## Generate from multivariate normal distribution
####################################################################
set.seed(1)
par(mfrow = c(2,2))

# Function produces samples from a multivariate normal
multinorm <- function(mu, Sigma, N = 5e2)
{
  # Eigenvalue (spectral) decomposition
  decomp <- eigen(Sigma)
  
  # Finding matrix square-root
  Sig.sq <- decomp$vectors %*% diag(decomp$values^(1/2)) %*% solve(decomp$vectors)
  
  samp <- matrix(0, nrow = N, ncol = 2)
  for(i in 1:N)
  {
    U1=runif(1)
    U2=runif(1)
    Z=c(sqrt(-2*log(U1))*cos(2*pi*U2),sqrt(-2*log(U1))*sin(2*pi*U2))
    samp[i,]=mu+Sig.sq%*%Z 
  }
  return(samp)
}

###
# First: Mean (-5, 10) and .5 correlation
mu <- c(-5,10)
Sigma <- matrix(c(1, .5, .5, 1), nrow = 2, ncol = 2)
samp <- multinorm(mu = mu, Sigma = Sigma)

par(mfrow = c(1,3))
plot(samp, asp = 1, main = "Correlation = .5", xlab = "x_1", ylab = "x_2")
plot(density(samp[,1]), main = "Marginal density for X1")
plot(density(samp[,2]), main = "Marginal density for X2")


par(mfrow = c(2,2))
plot(samp, asp = 1, main = "Correlation = .5", xlab = "x_1", ylab = "x_2")



### 
# Second: Mean (-5, 10) and .99 correlation
mu <- c(-5,10)
Sigma <- matrix(c(1, .99, .99, 1), nrow = 2, ncol = 2)
plot(samp, asp = 1, main = "Correlation = .99", xlab = "x_1", ylab = "x_2")

### 
# Third: Mean (-5, 10) and -.8 correlation
mu <- c(-5,10)
Sigma <- matrix(c(1, -.8, -.8, 1), nrow = 2, ncol = 2)
samp <- multinorm(mu = mu, Sigma = Sigma)
plot(samp, asp = 1, main = "Correlation = -.8", xlab = "x_1", ylab = "x_2")

### 
# Fourth: Mean (-5, 10) and no correlation
mu <- c(-5,10)
Sigma <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
samp <- multinorm(mu = mu, Sigma = Sigma)
plot(samp, asp = 1, main = "Correlation = 0", xlab = "x_1", ylab = "x_2")
