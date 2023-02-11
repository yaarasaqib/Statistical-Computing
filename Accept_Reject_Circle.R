##################################
## Accept-reject for obtaining
## sample uniformlyfrom a standard circle
## using a box as a proposal
##############################
set.seed(1)
circle_ar <- function() 
{
  accept <- 0
  counter <- 0   # count the number of loop
  while(accept == 0)
  {
    counter <- counter + 1
    U1 <- runif(1)
    U2 <- runif(1)
    
    U1 <- 2*U1-1
    U2 <- 2*U2-1 
    
    if( U1^2+U2^2 <=1 ) # fill condition
    {
      accept <- 1
      prop <- c(U1,U2)
      return(c(prop, counter))
    }
  }
}

# Simulation 10^4 samples from circle
N <- 1e4
samp <- matrix(0, ncol = 2, nrow = N)
counts <- numeric(length = N)
for(i in 1:N)
{
  foo <- circle_ar()  # I use foo as a dummy name
  samp[i,] <- foo[1:2]
  counts[i] <- foo[3]
}


4/pi
# [1] 1.27324
mean(counts)  # should be very close

# Plotting the obtained samples
# no paritcular part of the circle is favored more
# than any other part.
plot(samp[,1], samp[,2], xlab = "x", ylab = "y", 
     main = "Uniform samples from a circle", asp = 1)


