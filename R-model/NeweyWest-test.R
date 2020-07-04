# Test out the NeweyWest function from the "sandwich" package.

library(sandwich)


cat("AUTOCORRELATION TEST\n")

set.seed(1)

n <- 1000
f <- c(0.5,0.7,1.0,0.6,0.3,0.2,0.1)

true_beta <- c(2,3)

reps <- 200
beta <- matrix(NA,reps,2)

for (r in 1:reps)
{
  e <- filter (rnorm(n+length(f)-1), f)
  e <- e[!is.na(e)]
  
  x <- rep(rnorm(n/10),each=10)
  y <- true_beta[2]*x + true_beta[1] + e
  
  m <- lm(y~x)

  beta[r,] <- coef(m)

  if (r>3) next

  cat("\n\nREPETITION",r,"\n\n")
  
  print (as.vector (round (acf(e,lag.max=length(f)+4,plot=FALSE) $ acf, 2)))
  
  cat("\nLM STANDARD ERRORS:\n")  
  print (summary(m))
  
  cat("NW STANDARD ERRORS:\n")
  NW <- NeweyWest(m,adjust=TRUE,prewhite=TRUE,verbose=TRUE)
  print(NW)
  cat("\n")
  print(cbind(Estimate=coef(m),Std.Error=sqrt(diag(NW))))

}

cat ("\nActual standard deviations of estimates:", 
     round (sqrt (mean ((beta[,1]-true_beta[1])^2)), 5),
     round (sqrt (mean ((beta[,2]-true_beta[2])^2)), 5),
     "\n")


cat("\n\nHETEROSKEDASTICITY TEST\n")

set.seed(1)

n <- 1000

true_beta <- c(2,3)

reps <- 200
beta <- matrix(NA,reps,2)

for (r in 1:reps)
{
  e <- rnorm(n)
  e[1:(n/2)] <- e[1:(n/2)] * 10
  
  x <- rnorm(n)
  x[1:(n/2)] <- x[1:(n/2)]-2
  y <- true_beta[2]*x + true_beta[1] + e
  
  m <- lm(y~x)

  beta[r,] <- coef(m)

  if (r>3) next

  cat("\n\nREPETITION",r,"\n")

  cat("\nLM STANDARD ERRORS:\n")  
  print (summary(m))
  
  cat("NW STANDARD ERRORS:\n")
  NW <- NeweyWest(m,lag=0,adjust=TRUE,prewhite=TRUE,verbose=TRUE)
  print(NW)
  cat("\n")
  print(cbind(Estimate=coef(m),Std.Error=sqrt(diag(NW))))

}

cat ("\nActual standard deviations of estimates:", 
     round (sqrt (mean ((beta[,1]-true_beta[1])^2)), 5),
     round (sqrt (mean ((beta[,2]-true_beta[2])^2)), 5),
     "\n")


cat("\n\nHETEROSKEDASTICITY AND AUTOCORRELATION TEST\n")

set.seed(1)

n <- 1000
f <- c(0.5,0.7,1.0,0.6,0.3,0.2,0.1)

true_beta <- c(2,3)

reps <- 200
beta <- matrix(NA,reps,2)

for (r in 1:reps)
{
  e <- filter (rnorm(n+length(f)-1), f)
  e <- e[!is.na(e)]
  e[1:(n/2)] <- e[1:(n/2)] * 10
  
  x <- rep(rnorm(n/10),each=10)
  x[1:(n/2)] <- x[1:(n/2)]-2
  y <- true_beta[2]*x + true_beta[1] + e
  
  m <- lm(y~x)

  beta[r,] <- coef(m)

  if (r>3) next

  cat("\n\nREPETITION",r,"\n\n")
  
  print (as.vector (round (acf(e,lag.max=length(f)+4,plot=FALSE) $ acf, 2)))
  
  cat("\nLM STANDARD ERRORS:\n")  
  print (summary(m))
  
  cat("NW STANDARD ERRORS:\n")
  NW <- NeweyWest(m,adjust=TRUE,prewhite=TRUE,verbose=TRUE)
  print(NW)
  cat("\n")
  print(cbind(Estimate=coef(m),Std.Error=sqrt(diag(NW))))

}

cat ("\nActual standard deviations of estimates:", 
     round (sqrt (mean ((beta[,1]-true_beta[1])^2)), 5),
     round (sqrt (mean ((beta[,2]-true_beta[2])^2)), 5),
     "\n")
