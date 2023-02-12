 #############################
########## Lab 5 | SC ###############
 ######  based on weighted and likelihood #########
 
 
 #### Question 1
 
 # importance sampling 
 
 N <- 1e3   # size of importance sample
 h <- function(x){exp(-x)}

 Zs <- runif(N,min=0,max=1)   # importance sample from prop U(0,1)
 
 mean(h(Zs)*dbeta(Zs,shape1 = 4,shape2 = 5))

 # Weighted Importance Sampling Estimator
 N <- 1e3
 Zs <- runif(N, min = 0, max = 1)
 numerator <- sum( h(Zs)*dbeta(Zs,shape1 = 4,shape2 = 5))
 denominator <- sum(dbeta(Zs,shape1 = 4,shape2 = 5))
 WIS <- numerator/denominator
 print(WIS) 

############################################
 
 ##### Question 2 | compare the estimator variance #####
 
 ## variance of simple importance sampling estimator
 
 N <- 1e3
 Zs <- runif(N, min = 0, max = 1)
 func <- exp(-Zs) * dbeta(Zs, 4, 5)/1
 var(func) # sigma^2_g estimate
 
 # variance of weighted sampling estimator
 
 N <- 1e3
 reps <- 1e2
 ests <- numeric(length = reps)
 for(t in 1:reps)
 {
   Zs <- runif(N, min = 0, max = 1)
   num <- sum(h(x)*dbeta(Zs,4,5))
   den <- sum(dbeta(Zs,4,5))
   ests[t] <- num/den
 }
 var(ests)*N ## var of IS is less than weighted sampling
 #######################################
 
 ############ Question3  ##################
 ##### weighted importance sampling bivariate ##
 
 
 N <- 1e3
 samp <- matrix(0,nrow=N,ncol=2)
 f <- function(x,y)(exp(sin(x*y)))
 h <- function(x,y){x*y}
 
 for(i in 1:N){
   samp[1,i]<- pi*runif(1)
   samp[2,i]<- pi*runif(1)
 }
 num <- sum(h(samp[1,],samp[2,])*f(samp[1,],samp[,2]))
 den <- sum(f(samp[1,],samp[,2]))

 (WIS<- num/den) 
############################################
 
 ######## Question 4 , plotting likelihood graph ########
 
 x1 <- 2
 x2 <- 3  # sample values
 f <- function(theta){
   out <- 1/(2*pi)*exp(-(2-theta)^2/2)*exp(-(3-theta)^2/2)
   return(out)
  }
 x <- seq(-5,5,length=1e3)
 
 plot(x,f(x),xlab=expression(theta),ylab=expression(L(theta)),
      type="l", main="likelihood vs theta")
 abline(h=f(2.5),v=2.5,col=c("red"))
 
        
