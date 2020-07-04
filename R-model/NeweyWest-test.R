# Test out the NeweyWest function from the "sandwich" package.

library(sandwich)

set.seed(1)

n <- 1000
f <- c(0.5,0.7,1.0,0.6,0.3,0.2,0.1)

for (r in 1:3)
{
  cat("\n\nREPETITION",r,"\n")
  
  e <- filter (rnorm(n+length(f)-1), f)
  e <- e[!is.na(e)]
  
  acf(e,lag.max=length(f)+5)
  
  x <- rep(rnorm(n/10),each=10)
  y <- 3*x + 2 + e
  
  m <- lm(y~x)

  cat("\nLM STANDARD ERRORS:\n")  
  print (summary(m))
  
  cat("NW STANDARD ERRORS:\n")
  NW <- NeweyWest(m,lag=length(f)-1,adjust=TRUE,prewhite=TRUE,verbose=TRUE)
  print(NW)
  cat("\n")
  print(cbind(Estimate=coef(m),Std.Error=sqrt(diag(NW))))

}
