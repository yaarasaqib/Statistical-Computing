#########################
## Accept-reject for
## Beta(4,3) distribution
## Using U(0,1) proposal
#########################
set.seed(1)
beta_ar <- function() 
{
  c <- 60 *(3/5)^3 * (2/5)^2 
  accept <- 0
  counter <- 0   # count the number of loop
  while(accept == 0)
  {
    counter <- counter + 1
    U <- runif(1)
    prop=runif(1)
    ratio=dbeta(prop,shape1=4,shape2=3)/c
    
    if(U <= ratio)
    {
      accept <- 1
      return(c(prop, counter))
    }
  }
}


### Obtaining 10^4 samples from Beta() distribution
N <- 1e4
samp <- numeric(length = N)
counts <- numeric(length = N)
for(i in 1:N)
{
  rep <-  beta_ar()
    samp[i] <-  rep[1]
    counts[i] <- rep[2]
}

# Make a plot of the estimated density from the samples
# versus the true density
x <- seq(0, 1, length = 500)
plot(density(samp), main = "Estimated density from 1e4 samples")
lines(x,dbeta(x,shape1=4,shape2=3) , col = "red", lty = 2) ## Complete this
legend("topleft", lty = 1:2, col = c("black", "red"), legend = c("AR", "truth"))

# This is c
(c <- 60 *(3/5)^3 * (2/5)^2)

# This is the mean number of loops required
mean(counts)

#They should be almost the same!


#### ACcept Reject for Beta(2,0.1)

beta_2 <- function(){
  c <- 1.1
  accept <- 0
  count <- 0
  
  while(accept==0){
    count <- count+1
    U <- runif(1)
    x <- runif(1)
    prop <- 1-(1-x)^10
    ratio <- log(dbeta(prop,shape1=2,shape2=0.1))-log(c*prop)
    ratio <- exp(ratio)
    if(U <= ratio)
    {
      accept <- 1
      return(c(prop, count))
    }
  }
}


### Obtaining 10^4 samples from Beta() distribution
N <- 1e4
samp <- numeric(length = N)
counts <- numeric(length = N)
for(i in 1:N)
{
  rep <-  beta_2()
  samp[i] <-  rep[1]
  counts[i] <- rep[2]
}

# Make a plot of the estimated density from the samples
# versus the true density
x <- seq(0, 1, length = 500)
plot(density(samp), main = "Estimated density from 1e4 samples")
lines(x,dbeta(x,shape1=2,shape2=0.1) , col = "red", lty = 2) ## Complete this
legend("topleft", lty = 1:2, col = c("black", "red"), legend = c("AR", "truth"))

# This is c
(c <- 1.1)

# This is the mean number of loops required
mean(counts)
