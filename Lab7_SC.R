##########################
#### Lab 7 Worksheet #####
##########################


##### Lab Question 7.2 ####
## newton raphson on cosx

f.grad <- function(x) -sin(x)
f.hessian <- function(x) -cos(x)

tol <- 1e-10
compare <- 100
iter <- 1

xk <- c() # will store sequence here
xk[1] <- .5 # starting value

while(compare > tol){
  iter <- iter + 1 # tracking iterations
  gradient <- f.grad(xk[iter - 1])
  hessian <- f.hessian(xk[iter - 1])
  xk[iter] <- xk[iter - 1] - gradient/hessian
  compare <- abs(gradient)
}
iter

xk[iter] 

#### Gradient Ascent

tol <- 1e-10
compare <- 100
iter <- 1
xk <- c() # will store sequence here
gradient <- c() # will store gradient sequence here
xk[1] <- .5 # starting value
gradient[1] <- f.grad(xk[1])
t <- 1

while(compare > tol && iter < 100){
  iter <- iter + 1 # tracking iterations
  gradient[iter] <- f.grad(xk[iter - 1])
  xk[iter] <- xk[iter - 1] + t*gradient[iter]
  compare <- abs(gradient[iter])
}

xk[iter] # GA last iterate  

plot.ts(gradient) # wrt time


#################################
## 3 Titanic data logistic regression

titanic <- read.csv("https://dvats.github.io/assets/titanic.csv")
head(titanic)

f.gradient <- function(y, X, beta){
  # converting beta to compatible matrix form
  beta <- matrix(beta, ncol = 1)
  pi.vec <- 1 / (1 + exp(-X%*% beta))
  rtn <- colSums(X* as.numeric(y - pi.vec))
  return(rtn)
}


f.hessian <- function(y, X, beta)
{
  beta <- matrix(beta, ncol = 1)
  W_i <- exp(X%*%beta) / (1 + exp(X%*%beta))^2
  W <- diag(as.numeric(W_i))
  rtn <- -t(X) %*% W %*% X
}

y <- titanic$Survived
X <- as.matrix(titanic[, -1]) # everything but the first column is the X
# will need these later
p <- dim(X)[2]
n <- length(y)

tol <- 1e-10
compare <- 100
iter <- 1
# starting from the zero-vector
grad.vec <- c() # will store gradients here
beta.current <- rep(0, p)
beta.new <- beta.current

while(compare > tol)
{
  iter <- iter + 1 # tracking iterations
  gradient <- f.gradient(y, X, beta.current)
  hessian <- f.hessian(y, X, beta.current)
  beta.new <- beta.current - qr.solve(hessian) %*% gradient
  grad.vec[iter] <- norm(gradient, "2")
  beta.current <- beta.new
  compare <- grad.vec[iter]
}

iter

beta.new # N-R last iterate.

plot.ts(grad.vec)


### lab 7.4 ####
## Probability of jack and  Rose dying #


# 1 for intercept, 1 for male,
jack.x <- c(1, 1, 20, 0, 0, 7.5)
rose.x <- c(1, 0, 19, 1, 1, 512)
# estimate from logistic reg is in beta.new
pi.jack <- 1/ (1 + exp( - sum(jack.x * beta.new)))
pi.rose <- 1/ (1 + exp( - sum(rose.x * beta.new)))

pi.jack
pi.rose



## Implementing Gradient Ascent

tol <- 1e-10
compare <- 100
iter <- 1

# starting from the zero-vector
grad.vec <- c() # will store gradients here
beta.current <- rep(0, p)

beta.new <- beta.current

t <- 1
while(compare > tol && iter < 1000){
  iter <- iter + 1 # tracking iterations
  gradient <- f.gradient(y, X, beta.current)
  beta.new <- beta.current + t * gradient
  grad.vec[iter] <- norm(gradient, "2")
  beta.current <- beta.new
  compare <- grad.vec[iter]
}
iter

beta.new # GA last iterate.

plot.ts(grad.vec)

#In the above the learning rate is clearly too large, since ‖∇𝑙(𝛽)‖ oscillated the whole time. When this
#happens, it typically means that the learning rate 𝑡 is too large.

tol <- 1e-5
compare <- 100
iter <- 1

# starting from the zero-vector
grad.vec <- c() # will store gradients here
beta.current <- rep(0, p)
beta.new <- beta.current
t <- .000007
while(compare > tol && iter < 3e5)
{
  iter <- iter + 1 # tracking iterations
  gradient <- f.gradient(y, X, beta.current)
  beta.new <- beta.current + t * gradient
  grad.vec[iter] <- norm(gradient, "2")
  beta.current <- beta.new
  compare <- grad.vec[iter]
}
iter
beta.new # GA last iterate

plot.ts(grad.vec) 


## 7.6

#install.packages("faraway")
library(faraway)
?motorins    # to know about the data


set.seed(1)
n <- 50
p <- 5
# generate covariates
X <- cbind(1, matrix( rnorm(n*(p-1)), ncol = p-1, nrow = n))
beta.star <- c(.1, .2, .1, .6, -.3)
pi <- exp(X %*% beta.star)
y <- rpois(n, pi) # generate response

tol <- 1e-10
compare <- 100
iter <- 1
# starting from the zero-vector
grad.vec <- c() # will store gradients here
foo <- glm(y ~ X - 1, family = poisson)
beta.current <- rep(0, p)
beta.new <- beta.current


while(compare > tol)
{
  iter <- iter + 1 # tracking iterations
  gradient <- f.gradient(y, X, beta.current)
  hessian <- f.hessian(y, X, beta.current)
  beta.new <- beta.current - qr.solve(hessian) %*% gradient
  grad.vec[iter] <- norm(gradient, "2")
  beta.current <- beta.new
  compare <- grad.vec[iter]
}
iter

beta.new # N-R last iterate.

plot.ts(grad.vec)

#Gradient Ascent algorithm.
tol <- 1e-8
compare <- 100
iter <- 1



