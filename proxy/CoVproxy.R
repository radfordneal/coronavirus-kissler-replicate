# Code to reproduce production of proxies for common-cold coronavirus
# incidence in the Kissler, et al paper, and to investigate alternative
# proxies.  Reads the ILI proxies from ILIproxy.csv.  Produces various
# plots that are written to CoVproxy.pdf.  Writes data to CoVproxy.csv, 
# with all proxies produced included, with column names having the form 
# <virus>_<proxy>[<var>], with <proxy> being an ILI proxy, and the optional
# <var> as described below (eg, "OC43_proxyW" or "OC43_proxyAs").
#
# The "o" variant has (manually-identified) outliers adjusted, and
# zeros for the fraction testing positive for each coronavirus amongst 
# all positive tests replaced by non-zero values.
#
# The "s" variant has outliers and zeros adjusted and is then smoothed; 
# the "ss" variant is smoothed more. 

# The "m" variant has outliers and zeros adjusted, and then replaces
# percentages amongst coronavirus tests with values from a smoothing
# spline. The "n" variant is like "m" except that values for the total
# for all coronaviruses are also from a spline fit, rather than the
# actual data.
#
# Copyright 2020 by Radford M. Neal
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, a copy is available at
# https://www.gnu.org/licenses/old-licenses/gpl-2.0.html


library(splines)

pdf("CoVproxy.pdf",height=8,width=6)
par(mar=c(1.5,2.3,3,0.5),mgp=c(1.4,0.3,0),tcl=-0.22)

yrcols <- c("red","green","blue","orange","darkcyan","darkmagenta")

min_pct <- 0.007      # Minimum positive test percentage (for logarithmic plots)

zero_below <- 1.3     # Factor to put zeros under minimum non-zero (or other
                      # reference) value


# READ THE FILE OF ILI PROXIES.  As produced by ILIproxy.R.  Dates are
# for the start of the week (Sunday).

ILIproxy <- read.csv ("ILIproxy.csv",
                      header=TRUE, stringsAsFactors=FALSE)

ILIproxy <- ILIproxy[-1,]  # First week not used by Kissler, et al

ILIproxy$start <- as.Date(ILIproxy$start)


# INCLUDE UTILITY FUNCTIONS.  They need "start" and "week" to be defined.

start <- ILIproxy$start    # Start dates of weeks being analysed
year  <- ILIproxy$year     # Year for each week
week  <- ILIproxy$week     # Number of each week in its year

source("../util/util.R")


# READ THE RAW DATA ON CORONAVIRUS TESTS.  Data was obtained by
# Kissler, et al from NREVSS, and is weekly (date is last day of week)
# from 2014-07-05 (the week starting on 2014-06-29) to 2019-06-29 (the
# week starting on 2019-06-23).
#
# (The file has some fields computed by Kissler, et al as well as raw data.)

CoV_raw <- read.csv ("../data/nrevssCDC_ILI.csv", 
                     header=TRUE, stringsAsFactors=FALSE)

CoV_raw <- CoV_raw[-1,]    # First week not used by Kissler, et al


# CREATE A MORE CONVENIENT VERSION FROM THE RAW COV DATA.

CoV <- data.frame (start = as.Date(CoV_raw$WEEKEND) - 6)

stopifnot (nrow(ILIproxy)==nrow(CoV) && all(ILIproxy$start==CoV$start))

for (virus in viruses)
{ CoV[,virus] <- 100 * CoV_raw[,virus]  # convert to %
}

CoV$total <- 0
for (virus in viruses)
{ CoV$total <- CoV$total + CoV[,virus]
}


# PRINT SUMMARIES OF DATA.

cat("\nILI PROXY SUMMARY:\n\n")
print(summary(ILIproxy))

cat("\nCORONAVIRUS DATA SUMMARY:\n\n")
print(summary(CoV))


# CREATE TWEAKED PERCENTAGES AMONGST POSITIVE CORONAVIRUS TESTS.

pct <- data.frame(start=start)    # original data
pcto <- data.frame(start=start)   # with outliers and zeros adjusted
pcts <- data.frame(start=start)   #  - then smoothed once
pctss <- data.frame(start=start)  #  - and again

# Adjust outliers.  Manually identified.

adjCoV <- CoV

adjCoV$NL63[8]   <- mean(adjCoV$NL63[c(6,7,9,10)]) / zero_below
adjCoV$HKU1[61]  <- mean(adjCoV$HKU1[c(59,60,62,63)])
adjCoV$HKU1[127] <- mean(adjCoV$HKU1[c(125,126,128,129)])

# Recompute totals after outliers adjusted.

adjCoV$total[] <- 0
for (virus in viruses) 
{ adjCoV$total <- adjCoV$total + adjCoV[,virus]
}

# Adjust zeros.  Replace with minimum non-zero value for virus divided by 1.5.

for (virus in viruses)
{ nonzero <- adjCoV[,virus] != 0
  min_nonzero <- min(adjCoV[,virus][nonzero])
  adjCoV[,virus] <- pmax (min_nonzero/zero_below, adjCoV[,virus])
}

# Recompute totals after zero adjustment.

adjCoV$total[] <- 0
for (virus in viruses) 
{ adjCoV$total <- adjCoV$total + adjCoV[,virus]
}


for (virus in viruses)
{ 
  # Raw version.

  pct[,virus] <- 100 * CoV[,virus] / CoV$total 

  # With adjustment for outliers and zeros.

  pcto[,virus] <- 100 * adjCoV[,virus] / adjCoV$total

  # Smooth once.

  x <- pcto[,virus]

  x <- c(x[1],x,x[length(x)])     # extend so edge values won't be NA
  x <- filter (x, c(0.2,0.6,0.2)) # filter with three-point convolution
  x <- x[-c(1,length(x))]         # shrink to original length
  pcts[,virus] <- x

  # Smooth a second time.

  x <- c(x[1],x,x[length(x)])     # extend so edge values won't be NA
  x <- filter (x, c(0.2,0.6,0.2)) # filter (again) with three-point convolution
  x <- x[-c(1,length(x))]         # shrink to original length
  pctss[,virus] <- x
}


# PRODUCE VARIOUS PLOTS.

par(mfrow=c(4,1))

for (virus in c(viruses,"total"))
{ plot (start, CoV[,virus], pch=20, ylim=c(0,7+7*(virus=="total")), yaxs="i",
        ylab="percent")
  week_lines()
  title(paste("Percent positive tests,",virus))
  plot (start, pct_logit(pmax(min_pct,CoV[,virus])),
               pch=20, ylim=c(-9.5,-2.5)+(virus=="total"), ylab="logit")
  week_lines()
}

for (virus in c("total"))
{ plot (start, adjCoV[,virus], pch=20, ylim=c(0,7+7*(virus=="total")), yaxs="i",
        ylab="percent")
  week_lines()
  title(paste("Percent positive tests after outlier/zero adjustment,",virus))
  plot (start, pct_logit(adjCoV[,virus]),
               pch=20, ylim=c(-9.5,-2.5)+(virus=="total"), ylab="logit")
  week_lines()
}

par(mfrow=c(2,1))

plot (start, rep(0,length(start)), ylim=c(0,100), 
      type="n", yaxs="i", ylab="", xlab="")
for (virus in viruses)
{ lines (start, pct[,virus], lwd=2,
         col=c("orange","green","blue","red")[virus==viruses])
}
title("Pct HKU1(r), OC43(b), NL63(g), 229E(o) in positive tests")
week_lines()

plot (start, rep(0,length(start)), ylim=c(0,100), 
      type="n", yaxs="i", ylab="", xlab="")
for (virus in viruses)
{ lines (start, pcto[,virus], lwd=2,
         col=c("orange","green","blue","red")[virus==viruses])
}
title ("Percent in positive tests, adjusting for outliers")
week_lines()

plot (start, rep(0,length(start)), ylim=c(0,100), 
      type="n", yaxs="i", ylab="", xlab="")
for (virus in viruses)
{ lines (start, pcts[,virus], lwd=2,
         col=c("orange","green","blue","red")[virus==viruses])
}
title ("Pct in positive tests, adjusting for outliers, smoothing once")
week_lines()

plot (start, rep(0,length(start)), ylim=c(0,100), 
      type="n", yaxs="i", ylab="", xlab="")
for (virus in viruses)
{ lines (start, pctss[,virus], lwd=2,
         col=c("orange","green","blue","red")[virus==viruses])
}
title ("Pct in positive tests, adjusting for outliers, smoothing twice")
week_lines()

par(mfrow=c(2,1))

for (virus in c(viruses,"total"))
{ plot_vs_doy (start, adjCoV[,virus], pch=20, 
               ylim=c(0,7+7*(virus=="total")), yaxs="i")
  title(paste("Pct positive vs day of year, outliers & zeros adjusted,",virus))
  plot_vs_doy (start, pct_logit(adjCoV[,virus]), 
               pch=20, ylim=c(-9.5,-2.5)+(virus=="total"), ylab="logit")
}


# FIT SPLINE TO VALUES FOR LOGIT PERCENT AMONGST POSITIVE CORONAVIRUS TESTS.
# Fits done after outlier / zero adjustment. Two fake data points are added at 
# the beginning, equal to the average of the first four points, to reduce the 
# tendency of the spline to fit noise in this high-variance initial part.

par(mfrow=c(4,1))
pct_model <- list()
pctm <- data.frame(start=start)

for (virus in viruses)
{
  d <- pct_logit(pcto[,virus])
  t <- start
  d <- c(rep(mean(d[1:4]),2),d)
  t <- c(t[1]-c(14,7),t)
  pct_model[[virus]] <- lm (d ~ bs(t,df=35))
  pctm[,virus] <- predict(pct_model[[virus]]) [-c(1,2)]

  plot (start, pct_logit(pcto[,virus]), pch=20, ylab="logit")
  lines (start, pctm[,virus], col="red")
  week_lines()
  title (paste (
    "Spline fit for logit percent of",virus,"among coronavirus positive tests"))
}


# CREATE SPLINE FIT FOR LOGIT TOTAL PERCENT POSITIVE.  Fit is done after 
# outlier/zero adjustment.

par(mfrow=c(2,1))

tot_model <- lm (pct_logit(adjCoV$total) ~ bs(start,df=35))

plot (start, pct_logit(adjCoV$total), pch=20, ylab="log percent")
lines (start, predict(tot_model), col="red")
week_lines()
title ("Spline fit for logit total percent positive")

plot (start, residuals(tot_model), pch=20, ylab="")
abline(h=0)
week_lines()
title ("Residuals from spline fit for logit total percent positive")


# COMPUTE AND PLOT ALL THE PROXIES.

CoVproxy <- data.frame (start = as.character(start), year=year, week=week)

proxies <- names(ILIproxy) [substring(names(ILIproxy),1,5)=="proxy"]

for (virus in viruses)
{ par(mfrow=c(3,2))
  for (proxy in proxies)
  { 
    CoVproxy[,paste0(virus,"_",proxy)] <- 
      ILIproxy[,proxy] * CoV[,virus]
    CoVproxy[,paste0(virus,"_",proxy,"o")] <- 
      ILIproxy[,proxy] * adjCoV$total * pcto[,virus] / 100
    CoVproxy[,paste0(virus,"_",proxy,"s")] <- 
      ILIproxy[,proxy] * adjCoV$total * pcts[,virus] / 100
    CoVproxy[,paste0(virus,"_",proxy,"ss")] <- 
      ILIproxy[,proxy] * adjCoV$total * pctss[,virus] / 100
    CoVproxy[,paste0(virus,"_",proxy,"m")] <- 
      ILIproxy[,proxy] * adjCoV$total * pct_logit_inv(pctm[,virus])  /  100
    CoVproxy[,paste0(virus,"_",proxy,"n")] <- 
      ILIproxy[,proxy] * pct_logit_inv(predict(tot_model)) * 
                         pct_logit_inv(pctm[,virus]) / 100

    plot (start, log(pmax(min_pct,CoVproxy[,paste0(virus,"_",proxy)])),
          ylim=c(-5.3,3.4), pch=20,
          ylab=paste("Log",proxy,"proxy,",virus))
    # week_lines()
    if (proxy==proxies[1]) title(paste("Proxies for",virus))
    for (w in c("o","s","ss","m","n"))
    { plot (start, log(CoVproxy[,paste0(virus,"_",proxy,w)]),
            ylim=c(-4.8,3.4), pch=20,
            ylab=paste("Log",paste0(proxy,w),"proxy,",virus))
      # week_lines()
    }
  }
}

cat("\nCORONAVIRUS PROXY SUMMARY:\n\n")
print(summary(CoVproxy))


# PLOT SOME COMPARISONS.

par(mfrow=c(4,1))

for (pair in list (c("proxyW","proxyWXss"), 
                   c("proxyAXo","proxyDn"),
                   c("proxyDss","proxyDn")))
{ for (virus in viruses)
  { plot_two_with_lines (start, 
      log(CoVproxy[,paste0(virus,"_",pair[1])]), 
      log(CoVproxy[,paste0(virus,"_",pair[2])]),
      pch=19, ylab=paste("Line goes to log",pair[2]))
    week_lines()
    title(paste(pair[1],"versus",pair[2],"for",virus))
  }
}


# WRITE A FILE WITH THE VARIOUS PROXIES FOR CORONAVIRUS INCIDENCE.

write.table (CoVproxy, "CoVproxy.csv", sep=",",
             quote=FALSE, row.names=FALSE, col.names=TRUE)


# ALL DONE.

dev.off()
