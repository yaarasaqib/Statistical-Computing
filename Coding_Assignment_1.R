# Coding Assignment 1 
# sampling from ellipse

### Paste your function here


ellipse<- function(a,b){
  count <- 0
  accept <- 0
  
  a= sqrt(a)
  b= sqrt(b)
  
  while(accept == 0){
    count <- count+1
    
    U1 <- runif(1)
    U2 <- runif(1)
    
    U1 <- 2*a*U1- a
    U2 <- 2*b*U2- b
    
    if((U1/a)^2 +(U2/b)^2 <=1 ){
      accept <- accept+1
      return(c(U1,U2,count))
    }
  }
} 

N = 1e4 
samp <- matrix(0,ncol=2,nrow=N)
counts<- numeric(N)

for(i in 1:N){
  foo <- ellipse(3,4)
  samp[i,]<- foo[1:2]
  counts[i]<- foo[3]
}

plot(samp[,1],samp[,2],pch=16,main= "estimate of sample from ellipse ")

mean(counts)
