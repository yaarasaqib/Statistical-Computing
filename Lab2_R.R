##### Statistical Computing ####
### Lab 2 ####

### Accept reject Beta(4,3)

AR_Beta <- function(){
  try <- 0
  accept <- 0
  c  <- 60 *(3/5)^3 * (2/5)^2
  
  while(accept==0){
    try <- try+1
    U <- runif(1)
    prop <- runif(1)
    ratio <- dbeta(prop, 4, 3)/c
    if (U <= ratio){
      accept <- 1
      return(c(prop,try))
    }
    
  }
}
AR_Beta()


