# Code to reproduce production of proxies for common-cold coronavirus
# incidence in the Kissler, et al paper, and to investigate alternative
# proxies.  Reads the ILI proxies from ILIproxy.csv, and for proxies that
# have daily interpolations, from ILIproxy-daily.csv.  Produces various
# plots that are written to CoVproxy.pdf.  Writes data to CoVproxy.csv, 
# with all weekly proxies produced included, with column names having the 
# form <virus>_<proxy>[<var>], with <proxy> being an ILI proxy, and the 
# optional <var> as described below (eg, "OC43_proxyW" or "OC43_proxyAs").
# Daily proxies from daily ILI proxies are written to CoVproxy-daily.csv,
# with only the "n" variant currently produced.
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


# READ THE FILE OF WEEKLY ILI PROXIES.  As produced by ILIproxy.R.  Dates are
# for the start of the week (Sunday).

ILIproxy <- read.csv ("ILIproxy.csv",
                      header=TRUE, stringsAsFactors=FALSE)

ILIproxy <- ILIproxy[-1,]  # First week not used by Kissler, et al

ILIproxy$start <- as.Date(ILIproxy$start)


# READ THE FILE OF DAILY-INTERPOLATED ILI PROXIES.  As produced by ILIproxy.R.

ILIproxy_daily <- read.csv ("ILIproxy-daily.csv",
                             header=TRUE, stringsAsFactors=FALSE)

ILIproxy_daily <- ILIproxy_daily[-(1:7),]# First week not used by Kissler, et al
ILIproxy_daily$date <- as.Date(ILIproxy_daily$date)


# INCLUDE UTILITY FUNCTIONS.  They need "start" and "week" to be defined.

start <- ILIproxy$start    # Start dates of weeks being analysed
year  <- ILIproxy$year     # Year for each week
week  <- ILIproxy$week     # Number of each week in its year

all_days <- rep(start,each=7) + (0:6)        # Dates for every day of the period
#all_days <- all_days[1:(length(all_days)-6)]

yrcont <- (0:(length(start)-1)) / (365.24/7) # Continuous year over whole period
sncont <- yrcont %% (365.24/7)               # Continuous (0,1) value for season

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

# Compute total percent positive for all coronaviruses. Note that this could
# conceivably be greater than 100, though in fact it never is.

CoV$total <- 0
for (virus in viruses)
{ CoV$total <- CoV$total + CoV[,virus]
}


# PRINT SUMMARIES OF DATA.

cat("\nILI PROXY SUMMARY:\n\n")
print(summary(ILIproxy))

cat("\nCORONAVIRUS DATA SUMMARY:\n\n")
print(summary(CoV))


# PLOT ORIGINAL DATA.

par(mfrow=c(4,1))

for (virus in viruses)
{ plot (start, CoV[,virus], pch=20, ylim=c(0,7.1), yaxs="i", ylab="percent")
  week_lines()
  title(paste("Percent positive tests,",virus))
  plot (start, log(pmax(min_pct,CoV[,virus])),
               pch=20, ylim=c(-5,2), ylab="log, zeros at bottom")
  week_lines()
}


# CREATE TWEAKED PERCENTAGES AMONGST POSITIVE CORONAVIRUS TESTS.

pct <- data.frame(start=start)    # original data
pcto <- data.frame(start=start)   # with outliers and zeros adjusted
pcts <- data.frame(start=start)   #  - then smoothed once
pctss <- data.frame(start=start)  #  - and again

# Adjust outliers.  Manually identified.  Make plots too.

par(mfrow=c(4,1))

adjCoV <- CoV

adjCoV$NL63[8]   <- mean(adjCoV$NL63[c(6,7,9,10)]) / zero_below

plot (start, log(pmax(min_pct,CoV$NL63)), ylab="NL63", 
      ylim=c(-5,2), pch=20, col="gray")
points (start[8], log(pmax(min_pct,CoV$NL63[8])), pch=19, col="red")
lines (rep(start[8],2), log(pmax(min_pct,c(CoV$NL63[8],adjCoV$NL63[8]))))
title (
 "Adjustment of outliers in percent positive tests for each virus (log scale)")

adjCoV$E229[10]  <- mean(adjCoV$E229[c(8,9,11,12)])

plot (start, log(pmax(min_pct,CoV$E229)), ylab="E229", 
      ylim=c(-5,2), pch=20, col="gray")
points (start[10], log(pmax(min_pct,CoV$E229[10])), pch=19, col="red")
lines (rep(start[10],2), log(pmax(min_pct,c(CoV$E229[10],adjCoV$E229[10]))))

plot (start, log(pmax(min_pct,CoV$OC43)), ylab="OC43", pch=20, 
      ylim=c(-5,2), col="gray")

adjCoV$HKU1[61]  <- mean(adjCoV$HKU1[c(59,60,62,63)])
adjCoV$HKU1[127] <- mean(adjCoV$HKU1[c(125,126,128,129)])

plot (start, log(pmax(min_pct,CoV$HKU1)), ylab="HKU1", pch=20,
      ylim=c(-5,2), col="gray")
points (start[c(61,127)], log(pmax(min_pct,CoV$HKU1[c(61,127)])), 
        pch=19, col="red")
lines (rep(start[61],2), log(pmax(min_pct,c(CoV$HKU1[61],adjCoV$HKU1[61]))))
lines (rep(start[127],2), log(pmax(min_pct,c(CoV$HKU1[127],adjCoV$HKU1[127]))))

# Recompute totals after outliers adjusted.

adjCoV$total[] <- 0
for (virus in viruses) 
{ adjCoV$total <- adjCoV$total + adjCoV[,virus]
}

# Adjust zeros.  Replace with minimum non-zero value for virus divided 
# by zero_below.

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

plot (start, CoV$total, pch=20, ylim=c(0,14), yaxs="i", ylab="percent")
week_lines()
title("Total percent positive tests")
plot(start, log(pmax(min_pct,CoV$total)), pch=20, ylim=c(-1.5,3), ylab="log")
week_lines()

plot (start, adjCoV$total, pch=20, ylim=c(0,14), yaxs="i", ylab="percent")
week_lines()
title("Total Percent positive tests after outlier/zero adjustment")
plot (start, log(adjCoV$total), pch=20, ylim=c(-1.5,3), ylab="log")
week_lines()

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

for (virus in viruses)
{ plot_vs_doy (start, adjCoV[,virus], pch=20, ylim=c(0,7), yaxs="i")
  title(paste("Pct positive vs day of year, outliers & zeros adjusted,",virus))
  plot_vs_doy (start, pct_logit(adjCoV[,virus]), 
               pch=20, ylim=c(-9.5,-2.5), ylab="logit")
}


# FIT SPLINE TO VALUES FOR LOGIT PERCENT AMONGST POSITIVE CORONAVIRUS TESTS.
# Fits done after outlier / zero adjustment. Two fake data points are added at 
# the beginning, equal to the average of the first four points, to reduce the 
# tendency of the spline to fit noise in this high-variance initial part.

par(mfrow=c(4,1))
pct_model <- list()
pctm <- data.frame(start=start)

t <- start
t <- c(t[1]-c(14,7),t+7/2)
pctm_spline <- bs(t,df=35)

for (virus in viruses)
{
  cat ("\nSpline model for logit percent positive for",virus,"\n\n")

  d <- pct_logit(pcto[,virus])
  d <- c(rep(mean(d[1:4]),2),d)
  pct_model[[virus]] <- lm (d ~ pctm_spline)

  print(summary(pct_model[[virus]]))

  pctm[,virus] <- predict(pct_model[[virus]]) [-c(1,2)]

  plot (start, pct_logit(pcto[,virus]), pch=20, ylab="logit")
  lines (start, pctm[,virus], col="red")
  week_lines()
  title (paste (
    "Spline fit for logit percent of",virus,"among coronavirus positive tests"))
}

par(mfrow=c(4,1))

for (g in seq_along(virus_groups))
{ 
  virus_group <- virus_groups[[g]]
  plot (start, pct_logit(pcto[,virus_group[1]]+pcto[,virus_group[2]]),
        pch=20, ylab="logit")
  week_lines()
  title (paste ("Logit total percent of",names(virus_groups)[g],
                "viruses among coronavirus positive tests"))
}


# CREATE SPLINE FIT FOR LOG TOTAL PERCENT POSITIVE.  Fit is done after 
# outlier/zero adjustment.

par(mfrow=c(2,1))

cat("\nModels for log total percent\n\n")

tot_model_spline <- bs(start,df=40)
tot_model <- lm (log(adjCoV$total) ~ tot_model_spline)

print(summary(tot_model))

plot (start, log(adjCoV$total), pch=20, ylab="log percent")
lines (start, predict(tot_model), col="red")
week_lines()
title ("Spline fit for log total percent positive")

plot (start, residuals(tot_model), pch=20, ylab="")
abline(h=0)
week_lines()
title ("Residuals from spline fit for log total percent positive")

# tot_model2 <- lm (log(adjCoV$total) ~ bs(start,df=20)
#                                     + sin(1*2*pi*sncont) + cos(1*2*pi*sncont))
# 
# print(summary(tot_model2))
# 
# plot (start, log(adjCoV$total), pch=20, ylab="log percent")
# lines (start, predict(tot_model2), col="red")
# week_lines()
# title ("Spline + sine fit for log total percent positive")
# 
# plot (start, residuals(tot_model2), pch=20, ylab="")
# abline(h=0)
# week_lines()
# title ("Residuals from spline + sine fit for log total percent positive")


# COMPUTE AND PLOT ALL THE WEEKLY PROXIES.

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
      ILIproxy[,proxy] * exp(predict(tot_model)) * 
                         pct_logit_inv(pctm[,virus]) / 100

    plot (start, log(pmax(min_pct,CoVproxy[,paste0(virus,"_",proxy)])),
          ylim=c(-5.3,3.4), pch=20,
          ylab=paste0("Log ",proxy,", ",virus))
    # week_lines()
    if (proxy==proxies[1]) title(paste("Weekly proxies for",virus))
    for (w in c("o","s","ss","m","n"))
    { plot (start, log(CoVproxy[,paste0(virus,"_",proxy,w)]),
            ylim=c(-4.8,3.4), pch=20,
            ylab=paste0("Log ",paste0(proxy,w),", ",virus))
      # week_lines()
    }
  }
}

cat("\nCORONAVIRUS WEEKLY PROXY SUMMARY:\n\n")
print(summary(CoVproxy))


# COMPUTE AND PLOT ALL THE DAILY PROXIES.

CoVproxy_daily <- data.frame (date = as.character(all_days))

proxies <- names(ILIproxy_daily) [substring(names(ILIproxy_daily),1,5)=="proxy"]

par(mfrow=c(3,2))

tot_model_predict_all_days <- 
  as.vector (cbind (1, predict(tot_model_spline,all_days-7/2)) 
               %*% coef(tot_model))
                                
for (virus in viruses)
{ for (proxy in proxies)
  { 
    CoVproxy_daily[,paste0(virus,"_",proxy,"n")] <-
       ILIproxy_daily[,proxy] * exp(tot_model_predict_all_days) *
         pct_logit_inv (as.vector (cbind (1, predict(pctm_spline,all_days))
                                     %*% coef (pct_model[[virus]]))) / 100
    for (w in c("n"))
    { plot (all_days, log(CoVproxy_daily[,paste0(virus,"_",proxy,w)]),
            ylim=c(-4.8,3.4), pch=20,
            ylab=paste0("Log ",paste0(proxy,w),", ",virus))
      if (proxy==proxies[1]) title(paste("Daily proxies for",virus))
    }
  }
}

cat("\nCORONAVIRUS DAILY PROXY SUMMARY:\n\n")
print(summary(CoVproxy_daily))


# MORE WEEKLY PROXY PLOTS, GROUPED BY PROXY RATHER THAN VIRUS.

proxies_to_show <- c("proxyW","proxyAXo","proxyAXss","proxyDn","proxyEn")
for (proxy in proxies_to_show)
{ 
  par(mfrow=c(3,2))
  for (virus in viruses)
  { plot (start, CoVproxy[,paste0(virus,"_",proxy)],
          ylim=c(0,35), yaxs="i", pch=20,
          ylab=paste("Weekly",proxy,"for",virus))
    if (virus==viruses[1]) title(proxy)
    # week_lines()
  }

  par(mfrow=c(3,2))
  for (virus in viruses)
  { plot (start, log(pmax(min_pct,CoVproxy[,paste0(virus,"_",proxy)])),
          ylim=c(-5.3,3.4), pch=20,
          ylab=paste("Log weekly",proxy,"for",virus))
    if (virus==viruses[1]) title(paste("Log",proxy))
    # week_lines()
  }
}


# DAILY PROXY PLOTS, GROUPED BY PROXY RATHER THAN VIRUS.

proxies_to_show <- c("proxyEn")
for (proxy in proxies_to_show)
{ 
  par(mfrow=c(3,2))
  for (virus in viruses)
  { plot (all_days, CoVproxy_daily[,paste0(virus,"_",proxy)],
          ylim=c(0,35), yaxs="i", pch=20,
          ylab=paste("Daily",proxy,"for",virus))
    if (virus==viruses[1]) title(proxy)
  }

  par(mfrow=c(3,2))
  for (virus in viruses)
  { plot (all_days, log(pmax(min_pct,CoVproxy_daily[,paste0(virus,"_",proxy)])),
          ylim=c(-5.3,3.4), pch=20,
          ylab=paste("Log daily",proxy,"for",virus))
    if (virus==viruses[1]) title(paste("Log",proxy))
  }
}


# PLOT SOME COMPARISONS.

par(mfrow=c(4,1))

for (pair in list (c("proxyW","proxyWXss"), 
                   c("proxyDn","proxyW"),
                   c("proxyAXo","proxyDn"),
                   c("proxyDss","proxyDn"),
                   c("proxyDm","proxyDn"),
                   c("proxyDn","proxyEn")))
{ for (virus in viruses)
  { plot_two_with_lines (start, 
      log(CoVproxy[,paste0(virus,"_",pair[1])]), 
      log(CoVproxy[,paste0(virus,"_",pair[2])]),
      pch=19, ylab=paste("Line goes to log",pair[2]))
    week_lines()
    title(paste(pair[1],"versus",pair[2],"for",virus))
  }
}


# WRITE FILES WITH THE VARIOUS PROXIES FOR CORONAVIRUS INCIDENCE.

write.table (CoVproxy, "CoVproxy.csv", sep=",",
             quote=FALSE, row.names=FALSE, col.names=TRUE)


write.table (CoVproxy_daily, "CoVproxy-daily.csv", sep=",",
             quote=FALSE, row.names=FALSE, col.names=TRUE)


# ALL DONE.

warnings()

dev.off()
