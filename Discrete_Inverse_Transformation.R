### Inverse Transformation ####

## Negative Binomial

rm(list=ls())
set.seed(221348)

Inv_Nbin <- function(k,p){
  U=runif(1)
  p0 = p^k
  i=k
  A= p0
  repeat{
    if(U <=A){
      break
    }else{
      p0 = (i/(i-k+1))*(1-p)*p0
      A = A+ p0
      i=i+1
    }
  }
  return(i)
}

Inv_Nbin(10,.6)