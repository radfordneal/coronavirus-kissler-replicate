# Tests of immunity decline schemes.

par(mfcol=c(4,2))

n <- 400
x <- 0:n

a <- 0.01
a2 <- 0.03
a3 <- 0.005

b <- 0.001

plot (x,exp(-a*x),
      ylim=c(0,1),yaxs="i",xaxs="i",type="l")
abline(h=seq(0.2,1.0,by=0.2),v=seq(100,n,by=100))

# From maple:
#
# > limit ((1/b)*(exp(-ax) - (1-b)*exp(-ax/(1-b))), b=0);
#                                     1 + ax
#                                     -------
#                                     exp(ax)
#
# > limit((a2/(a2-a+a*b)) * (exp(-a*x) - (a/a2)*(1-b)*exp(-a2*x/(1-b))),b=0);
#                           -a2 exp(a2 x) + a exp(a x)
#                  ---------------------------------------------
#                  -exp(a x) a2 exp(a2 x) + exp(a x) exp(a2 x) a
#
# > simplify(limit((a2/(a2-a+a*b)) * (exp(-a*x) 
#                    - (a/a2)*(1-b)*exp(-a2*x/(1-b))),b=0));
#                               (-a2 exp(a2 x) + a exp(a x)) exp(-x (a + a2))
#                               ---------------------------------------------
#                                                  -a2 + a


plot (x, (1/b) * (exp(-a*x) - (1-b)*exp(-a*x/(1-b))),
      ylim=c(0,1),yaxs="i",xaxs="i",type="l")
abline(h=seq(0.2,1.0,by=0.2),v=seq(100,n,by=100))

plot (x,(1+a*x)*exp(-a*x),
      ylim=c(0,1),yaxs="i",xaxs="i",type="l")
abline(h=seq(0.2,1.0,by=0.2),v=seq(100,n,by=100))

y <- c(1,numeric(n))
z <- c(0,numeric(n))
z2 <- c(0,numeric(n))
z3 <- c(0,numeric(n))

for (i in 1:n)
{ y[i+1] <- y[i] - a*y[i]
  z[i+1] <- z[i] + a*y[i] - a*z[i]
  z2[i+1] <- z2[i] + a*y[i] - a2*z2[i]
  z3[i+1] <- z3[i] + a*y[i] - a3*z3[i]
}

plot (x, y+z,
      ylim=c(0,1),yaxs="i",xaxs="i",type="l")
abline(h=seq(0.2,1.0,by=0.2),v=seq(100,n,by=100))

plot (x, (a2/(a2-a+a*b)) * (exp(-a*x) - (a/a2)*(1-b)*exp(-a2*x/(1-b))),
      ylim=c(0,1),yaxs="i",xaxs="i",type="l")
abline(h=seq(0.2,1.0,by=0.2),v=seq(100,n,by=100))

plot (x, y+z2,
      ylim=c(0,1),yaxs="i",xaxs="i",type="l")
abline(h=seq(0.2,1.0,by=0.2),v=seq(100,n,by=100))

plot (x, (a3/(a3-a+a*b)) * (exp(-a*x) - (a/a3)*(1-b)*exp(-a3*x/(1-b))),
      ylim=c(0,1),yaxs="i",xaxs="i",type="l")
abline(h=seq(0.2,1.0,by=0.2),v=seq(100,n,by=100))

plot (x, y+z3,
      ylim=c(0,1),yaxs="i",xaxs="i",type="l")
abline(h=seq(0.2,1.0,by=0.2),v=seq(100,n,by=100))

