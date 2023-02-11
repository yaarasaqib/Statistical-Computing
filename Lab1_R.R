##### Lab1  ####


#### Ex1 PRNG ###

a= 7^5
m=2^30-1
x <- numeric(length=1e3)
x[1] <- 7

for(i in 2:1e3){
  x[i]= a*x[i-1] %% m
}

par(mfrow=c(1,2))
hist(x/m)
plot.ts(x/m)

###################################
### Ex2 Mixed Congruential method  ###
set.seed(2)
a= 7^5
m=2^30-1
c=7

x <- numeric(length=1e3)
for( i in 2:1e3){
  x[i]= (x[i-1]*a+c) %% m
}

hist(x/m)
plot.ts(x/m)
####################################

## Accept Reject Binom(n,p)

AR_Binom <-function(n,p){
  x <- seq(0,n,by=1)
  c <- max(choose(n,x)*p^(x-1)*(1-p)^(n-2*x))
  try <- 0
  accept <- 0
  while(accept==0){
    try <- try+1
    U <- runif(1)
    prop <- rgeom(1,p)
    ratio <- dbinom(prop,n,p)/(c*dgeom(prop,prob=p))
    if(U <= ratio){
      accept <- accept+1
      return(c(prop,try,c))
    }
  }
}

AR_Binom(10,0.25)

N <- 1e3
samp <-rep(0,N)
counts <- rep(0,N)
for(i in 1:N){
  foo <- AR_Binom(10,0.25)
  samp[i]<- foo[1]
  counts[i]<- foo[2]
  c <- foo[3]
  }

mean(samp)
mean(counts)
c

############################


(c=1/ppois(20,lambda=20))

trunc_pois <- function(m = 20, lam = 20)
{
  accept <- 0
  try <- 0
  while(accept==0)
  {
    try <- try + 1
    prop <- rpois(1, lambda = 20)
    
    if(prop <= m)
    {
      accept <- 1
    }
  }
  return(c(prop, try))
}

trunc_pois()
N <- 1e3
out <- replicate(N, trunc_pois())
# out is 2 x 1000 matrix. First row are the samples
# second row are the number of loops
mean(out[1, ]) # mean of trunc pois
mean(out[2, ]) # mean of number loops, similar to c
hist(out[1, ], main = "Hist of Truncated Poisson")
