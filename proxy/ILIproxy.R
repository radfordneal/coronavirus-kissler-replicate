# Code to reproduce production of a proxy for Influenza-Like-Illness (ILI)
# incidence in the Kissler, et al paper, and to investigate alternative
# proxies.  Produces various plots that are written to ILIproxy.pdf.
# Writes data for all proxies to ILIproxy.csv, with names the same as
# the proxy name (eg, "proxyA").  ProxyW is the one used by Kissler, et al.
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

pdf("ILIproxy.pdf",height=8,width=6)
par(mar=c(1.5,2.3,3,0.5),mgp=c(1.4,0.3,0),tcl=-0.22)
yrcols <- c("red","green","blue","orange","darkcyan","darkmagenta")


# READ THE RAW DATA ON ILI VISITS.  Data was obtained by Kissler, et al 
# from FluView, and is weeky from 1997 week 40 through 2020 week 7.
#
# Weeks start on Sunday, with first week of a year being the first
# with at least four days in that year.

ILI_raw <- read.csv ("../data/ILINet.csv", 
                     header=TRUE, stringsAsFactors=FALSE)


# CREATE A MORE CONVENIENT VERSION FROM THE RAW ILI DATA.  Note that the 
# column names were tweaked from the original version, to reduce special
# characters.

ILI <- data.frame (start = as.Date("1997-09-28") + 7 * (0:(nrow(ILI_raw)-1)))
ILI$year  <- ILI_raw$YEAR
ILI$week  <- ILI_raw$WEEK
ILI$syear <- ifelse (ILI$week<27, ILI$year-1, ILI$year)

ILI$providers <- ILI_raw$NUM.OF.PROVIDERS

ILI$all_visits <- ILI_raw$TOTAL
ILI$ili_visits <- ILI_raw$ILITOTAL
ILI$age_0_4    <- ILI_raw$AGE.0.4
ILI$age_5_24   <- ILI_raw$AGE.5.24
ILI$age_65_up  <- ILI_raw$AGE.65
ILI$age_25_64  <- ILI$ili_visits - ILI$age_0_4 - ILI$age_5_24 -ILI$age_65_up

ILI$weighted_pct_ili   <- ILI_raw$PCT.WEIGHTED.ILI
ILI$unweighted_pct_ili <- ILI_raw$PCT.UNWEIGHTED.ILI

ILI$nonili_visits <- ILI$all_visits - ILI$ili_visits

# Verify that the unweighted_pct_ili is just ili_vists/all_visits.
# (Weighted_pc_ili adjusts for census population of regions.)

stopifnot (
 max (abs (100 * (ifelse (ILI$all_visits==0 ,0, ILI$ili_visits/ILI$all_visits))
           - ILI$unweighted_pct_ili)) < 1e-5)


# REMOVE THE DATA FOR WEEKS BEFORE AND AFTER THE CORONAVIRUS DATA. That is,
# keep from the week starting 2014-06-29 through the week starting 2019-06-23.

ILI <- ILI [ILI$start >= "2014-06-29" & ILI$start <= "2019-06-23", ]


# INCLUDE UTILITY FUNCTIONS.  They need "start" and "week" to be defined.

start <- ILI$start  # Start dates of weeks being analysed
year  <- ILI$year   # Year for each week
week  <- ILI$week   # Number of each week in its year

source("../util/util.R")


# PLOT NUMBER OF PROVIDERS.

par(mfrow=c(2,1))

plot (start, ILI$providers, pch=20); week_lines()
title("Number of providers")

plot (start, log(ILI$providers), pch=20); week_lines()
title("Log number of providers")

plot_vs_doy (start, ILI$providers, pch=20)
title("Number of providers vs day of year from July 1")

plot_vs_doy (start, log(ILI$providers), pch=20)
title("Log number of providers vs day of year from July 1")

plot_vs_doy (start, ILI$all_visits/ILI$providers, pch=20)
title("All visits per provider vs day of year from July 1")

plot_vs_doy (start, log(ILI$all_visits/ILI$providers), pch=20)
title("Log all visits per provider vs day of year from July 1")


# PLOT VISITS PER PROVIDER, ILI AND NON-ILI.

plot_vs_doy (start, log(ILI$nonili_visits/ILI$providers), ylim=c(5,7.5), pch=20)
title("Log non-ILI visits per provider vs day of year from July 1")

plot_vs_doy (start, log(ILI$ili_visits/ILI$providers), ylim=c(1.2,3.7), pch=20)
title("Log ILI visits per provider vs day of year from July 1")

plot (start, log(ILI$nonili_visits/ILI$providers), pch=20); week_lines()
title("Log non-ILI visits per provider")

plot (start, log(ILI$ili_visits/ILI$providers), pch=20); week_lines()
title("Log ILI visits per provider")


# PLOT ILI AND NON-ILI VISITS. Done with linear and log forms, and
# with significant weeks marked. Log scales have same size divisions
# for both log plots.

par(mfrow=c(2,1))

plot (start, ILI$nonili_visits, pch=20); week_lines()
title("Non-ILI physician visits")
plot (start, log(ILI$nonili_visits), ylim=c(13,16), pch=20); week_lines()
title("Log non-ILI physician visits")

plot (start, ILI$ili_visits, pch=20); week_lines()
title("ILI physician visits")
plot (start, log(ILI$ili_visits), ylim=c(8.5,11.5), pch=20); week_lines()
title("Log ILI physician visits")

par(mfrow=c(2,1))

plot_vs_doy (start, log(ILI$nonili_visits), pch=20)
abline(h=c(13.0,13.5,14.0,14.5),col="gray")
title("Log non-ILI visits vs day of year from July 1")

plot_vs_doy (start, log(ILI$ili_visits), pch=20)
abline(h=c(8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0),col="gray")
title("Log ILI visits vs day of year from July 1")


# PLOT ILI VISITS OF DIFFERENT AGE GROUPS.

par(mfrow=c(4,1))

plot (start, ILI$age_0_4, pch=20, ylab="age 0-4"); week_lines()
title("ILI visits for different age groups")
plot (start, ILI$age_5_24, pch=20, ylab="age 5-24"); week_lines()
plot (start, ILI$age_25_64, pch=20, ylab="age 25-64"); week_lines()
plot (start, ILI$age_65_up, pch=20, ylab="age 65 up"); week_lines()

plot (start, log(ILI$age_0_4), pch=20, ylab="age 0-4, log"); week_lines()
title("Log ILI visits for different age groups")
plot (start, log(ILI$age_5_24), pch=20, ylab="age 5-24, log"); week_lines()
plot (start, log(ILI$age_25_64), pch=20,ylab="age 25-64, log");week_lines()
plot (start, log(ILI$age_65_up), pch=20,ylab="age 65 up, log");week_lines()


# PLOT WEIGHTED AND UWEIGHTED ILI PERCENTAGE, AND RATIO.  Also repeat the
# plots of numbers of providers, for reference.

par(mfrow=c(2,1))

plot (start, ILI$weighted_pct_ili, pch=20); week_lines()
title("Weighted percent ILI")
plot (start, log(ILI$weighted_pct_ili), pch=20); week_lines()

plot (start, ILI$unweighted_pct_ili, pch=20); week_lines()
title("Unweighted percent ILI")
plot (start, log(ILI$unweighted_pct_ili), pch=20); week_lines()

plot (start, ILI$weighted_pct_ili/ILI$unweighted_pct_ili, pch=20)
week_lines(); abline(h=1)
title("Weighted / Unweighted percent ILI ratio")
plot (start, log(ILI$weighted_pct_ili/ILI$unweighted_pct_ili), pch=20)
week_lines(); abline(h=0)
title("Log Weighted / Unweighted percent ILI ratio")

plot(start,ILI$providers,pch=20); week_lines()
title("Number of providers")
plot(start,log(ILI$providers),pch=20); week_lines()
title("Log number of providers")


# MODEL THE NUMER OF NON-ILI VISITS.

nonili_visits_spline_df <- 7
nonili_visits_spline <- bs(start,df=nonili_visits_spline_df)

par(mfrow=c(2,1))

# spline + providers

visits_mod1 <- lm (log(nonili_visits) 
                     ~ nonili_visits_spline + log(providers)
                   , data=ILI)
print(summary(visits_mod1))

plot (start, log(ILI$nonili_visits), pch=20); week_lines()
lines (start, predict(visits_mod1), col="blue", lwd=2)
title("Non-ILI model, spline+providers")

plot (start, residuals(visits_mod1), pch=20, ylab="Residuals"); week_lines()
abline(h=0)

# spline + providers + holidays

visits_mod2 <- lm (log(nonili_visits) 
                     ~ nonili_visits_spline + log(providers) 
                     + July4th_indicator
                     + Labor_Day_indicator 
                     + Thanksgiving_indicator
                     + Christmas_indicator
                     + New_Year_indicator
                   , data=ILI)
print(summary(visits_mod2))

par(mfrow=c(2,1))

plot (start, log(ILI$nonili_visits), pch=20); week_lines()
lines (start, predict(visits_mod2), col="blue", lwd=2)
title("Non-ILI model, spline+providers+holidays")

plot (start, residuals(visits_mod2), pch=20, ylab="Residuals"); week_lines()
abline(h=0)

# spline + providers + holidays + seasonality

xtra_df_prov <- 1  # df from providers
xtra_df_ind <- 5   # df from indicators
xtra_df <- 1 + xtra_df_prov + xtra_df_ind

visits_mod3 <- lm (log(nonili_visits) 
                     ~ nonili_visits_spline + log(providers)
#                    + I(log(providers)^2) + I(log(providers)^3)
                     + July4th_indicator
                     + Labor_Day_indicator 
                     + Thanksgiving_indicator
                     + Christmas_indicator
                     + New_Year_indicator
                     + sin(2*pi*week/52) + cos(2*pi*week/52)
                     + sin(2*2*pi*week/52) + cos(2*2*pi*week/52)
                     + sin(3*2*pi*week/52) + cos(3*2*pi*week/52)
                     + sin(4*2*pi*week/52) + cos(4*2*pi*week/52)
                     + sin(5*2*pi*week/52) + cos(5*2*pi*week/52)
#                    + sin(6*2*pi*week/52) + cos(6*2*pi*week/52)
#                    + sin(7*2*pi*week/52) + cos(7*2*pi*week/52)
                   , data=ILI, x=TRUE)
print(summary(visits_mod3))

residuals_nonili_visits3 <- residuals(visits_mod3)

smooth_nonili_visits3 <- as.vector (exp (cbind (1,
   nonili_visits_spline, log(ILI$providers), 
#  log(ILI$providers)^2, log(ILI$providers)^3, 
   July4th_indicator,
   Labor_Day_indicator,
   Thanksgiving_indicator,
   Christmas_indicator,
   New_Year_indicator
  ) %*% coef(visits_mod3) [1:(nonili_visits_spline_df+xtra_df)]))

par(mfrow=c(2,1))

plot (start, log(ILI$nonili_visits), pch=20); week_lines()
lines (start, predict(visits_mod3), col="blue", lwd=2)
lines (start, log(smooth_nonili_visits3), col="orange", lwd=2)
title("Non-ILI model, spline+providers+holidays+season")

plot (start, residuals(visits_mod3), pch=20, ylab="Residuals"); week_lines()
abline(h=0)
title("Residuals")

plot_vs_doy (start, residuals(visits_mod3), pch=20)
abline(h=0)
title("Residuals vs DOY")

plot (start, log(ILI$nonili_visits) - log(smooth_nonili_visits3), pch=20, ylab="")
week_lines(); abline(h=0)
title("Residuals without seasonal component")

xx <- visits_mod3$x
xx[,1:(nonili_visits_spline_df+xtra_df)] <- 0

seasonal_visits3 <- as.vector (xx %*% coef(visits_mod3))

plot (start, seasonal_visits3, type="l"); week_lines()
title("Seasonal component of model (log domain)")

xx <- visits_mod3$x
xx[,-(2:(nonili_visits_spline_df+1))] <- 0

trend_visits3 <- as.vector (xx %*% coef(visits_mod3))

plot (start, trend_visits3, type="l"); week_lines()
title("Trend component of model (log domain)")

xx <- visits_mod3$x
xx[,1:(1+nonili_visits_spline_df+xtra_df_prov)] <- 0
xx[,-(1:(nonili_visits_spline_df+xtra_df))] <- 0

holiday_visits3 <- as.vector (xx %*% coef(visits_mod3))

plot (start, holiday_visits3, pch=20); week_lines()
title("Holiday component of model (log domain)")


# MODEL THE NUMER OF ILI VISITS.

ili_bound <- range(start)

ili_knots_within_year <- 
  as.Date(c("2000-07-01","2000-08-01","2000-09-01","2000-10-15","2000-12-01",
            "2001-01-01","2001-02-01","2001-03-15","2001-05-01","2001-06-01")) -
  as.Date("2000-07-01")

ili_knots <- rep (seq (start[1],start[ILI$syear==2018][1],length=5),
                  each = length(ili_knots_within_year)) +
             ili_knots_within_year

ili_knots <- ili_knots[-1]

ili_visits_spline_df <- length(ili_knots) + 3
ili_visits_spline <- bs (start, knots=ili_knots, Bound=ili_bound)

# spline + providers + holidays + seasonality

ili_visits_mod3 <- lm (log(ili_visits) 
                     ~ log(providers)
#                    + I(log(providers)^2) + I(log(providers)^3)
                     + July4th_indicator
                     + Labor_Day_indicator 
                     + Thanksgiving_indicator
                     + Christmas_indicator
                     + New_Year_indicator
                     + ili_visits_spline
                     + sin(2*pi*week/52) + cos(2*pi*week/52)
                     + sin(2*2*pi*week/52) + cos(2*2*pi*week/52)
                     + sin(3*2*pi*week/52) + cos(3*2*pi*week/52)
#                    + sin(4*2*pi*week/52) + cos(4*2*pi*week/52)
#                    + sin(5*2*pi*week/52) + cos(5*2*pi*week/52)
#                    + sin(6*2*pi*week/52) + cos(6*2*pi*week/52)
#                    + sin(7*2*pi*week/52) + cos(7*2*pi*week/52)
#                    + sin(8*2*pi*week/52) + cos(8*2*pi*week/52)
#                    + sin(9*2*pi*week/52) + cos(9*2*pi*week/52)
                  , data=ILI, x=TRUE)
print(summary(ili_visits_mod3))

par(mfrow=c(2,1))

plot (start, log(ILI$ili_visits), pch=20); week_lines()
lines (start, predict(ili_visits_mod3), col="blue", lwd=2)
title("ILI model, spline+providers+holidays+season")

plot (start, residuals(ili_visits_mod3), pch=20, ylab="Residuals")
week_lines(); abline(h=0)
title("Residuals")

plot_vs_doy (start, residuals(ili_visits_mod3), pch=20)
abline(h=0)
title("Residuals vs DOY")

# spline + providers + holidays + seasonality + non-ili residuals

ili_visits_mod4 <- lm (log(ili_visits) 
                     ~ log(providers)
#                    + I(log(providers)^2) + I(log(providers)^3)
                     + July4th_indicator
                     + Labor_Day_indicator 
                     + Thanksgiving_indicator
                     + Christmas_indicator
                     + New_Year_indicator
                     + residuals_nonili_visits3
                     + ili_visits_spline
                     + sin(2*pi*week/52) + cos(2*pi*week/52)
                     + sin(2*2*pi*week/52) + cos(2*2*pi*week/52)
                     + sin(3*2*pi*week/52) + cos(3*2*pi*week/52)
#                    + sin(4*2*pi*week/52) + cos(4*2*pi*week/52)
#                    + sin(5*2*pi*week/52) + cos(5*2*pi*week/52)
#                    + sin(6*2*pi*week/52) + cos(6*2*pi*week/52)
#                    + sin(7*2*pi*week/52) + cos(7*2*pi*week/52)
#                    + sin(8*2*pi*week/52) + cos(8*2*pi*week/52)
#                    + sin(9*2*pi*week/52) + cos(9*2*pi*week/52)
                  , data=ILI, x=TRUE)
print(summary(ili_visits_mod4))

par(mfrow=c(2,1))

plot (start, log(ILI$ili_visits), pch=20); week_lines()
lines (start, predict(ili_visits_mod4), col="blue", lwd=2)
title("ILI model, spline+providers+holidays+season+noniliresid")

plot (start, residuals(ili_visits_mod4), pch=20, ylab="Residuals")
week_lines(); abline(h=0)
title("Residuals")

plot_vs_doy (start, residuals(ili_visits_mod4), pch=20)
abline(h=0)
title("Residuals vs DOY")

xx <- ili_visits_mod4$x
xx[,2:8] <- 0
ili_spline_season <- as.vector (xx %*% coef(ili_visits_mod4))

plot(start,ili_spline_season,type="l",col="darkgray",lwd=3)
points(start,ili_spline_season+residuals(ili_visits_mod4),pch=20,ylab="")
week_lines()
title("Spline + seasonality, and that plus residuals")


# CREATE PROXYA FROM RATIO OF ILI VISITS TO NON-ILI VISITS.

proxyA <- 100 * ILI$ili_visits / ILI$nonili_visits

par(mfrow=c(4,1))

plot (start, proxyA, pch=20); week_lines()
title("ProxyA: 100 x ratio of ILI visits to non-ILI visits")

plot (start, log(proxyA), pch=20, ylim=c(-0.5,2.5))
week_lines()


# CREATE PROXYB FROM RATIO OF ILI VISITS TO PREDICTED NON-ILI VISITS.
# Predictions omit the seasonal terms, so as not to suppress ILI
# seasonality.

proxyB <- 100 * ILI$ili_visits / smooth_nonili_visits3

plot (start, proxyB, pch=20); week_lines()
title("ProxyB: 100 x ratio of ILI visits to predicted non-ILI visits")

plot (start, log(proxyB), pch=20, ylim=c(-0.5,2.5))
week_lines()


# CREATE PROXYC FROM ILI VISITS MINUS MODEL PREDICTIONS OF IRRELEVANT PART.
# Uses ili_visits_mod4 to produce predictions for variation in ILI from
# number of providers, from residuals from the model for non-ILI visits,
# and from holidays, which are considered to not relate to actual ILI 
# incidence.  Retains the spline over time, the seasonal component, and
# the model residuals.  The result is divided by a scale factor just to get
# it to numerically match (approximately) the other proxies.

proxyC <- exp (ili_spline_season+residuals(ili_visits_mod4)) / 80

plot (start, proxyC, pch=20); week_lines()
title("ProxyC: Derived from model of ILI visits")

plot (start, log(proxyC), pch=20, ylim=c(-0.5,2.5))
week_lines()


# CREATE PROXYD.  Like proxyC, but without adding in the residuals.

proxyD <- exp (ili_spline_season) / 85

plot (start, proxyD, pch=20); week_lines()
title("ProxyD: Derived from model of ILI visits, without residuals")

plot (start, log(proxyD), pch=20, ylim=c(-0.5,2.5))
week_lines()


# PLOTS COMPARING ILI PROXIES.

par(mfrow=c(4,1))

plot (start, log(ILI$unweighted_pct_ili), pch=20, ylim=c(-0.5,2.5))
week_lines()
title("Unweighted percent ILI")

plot (start, log(ILI$weighted_pct_ili), pch=20, ylim=c(-0.5,2.5))
week_lines()
title("Weighted percent ILI")

plot (start, log(ILI$unweighted_pct_ili/ILI$weighted_pct_ili), pch=20, ylab="")
abline(h=0); week_lines()
title("Log ratio of unweighted % ILI to weighted % ILI")

plot_vs_doy (start, log(ILI$unweighted_pct_ili/ILI$weighted_pct_ili), pch=20)
abline(h=0)
title("Log ratio of unweighted % ILI to weighted % ILI vs DOY")

par(mfrow=c(2,1))

plot_vs_doy (start, log(ILI$unweighted_pct_ili), pch=20)
abline(h=seq(-0.5,2.0,by=0.5),col="gray")
title("Log unweighted percent ILI vs DOY")

plot_vs_doy (start, log(ILI$unweighted_pct_ili), pch=20, 
             xlim=c(34,140), ylim=c(-0.4,0.9))
abline(h=seq(-0.5,1.0,by=0.25),col="gray")

plot_vs_doy (start, log(ILI$weighted_pct_ili), pch=20)
abline(h=seq(-0.5,2.0,by=0.5),col="gray")
title("Log weighted percent ILI vs DOY")

plot_vs_doy (start, log(ILI$weighted_pct_ili), pch=20,
             xlim=c(34,140), ylim=c(-0.4,0.9))
abline(h=seq(-0.5,1.0,by=0.25),col="gray")

plot (log(ILI$weighted_pct_ili), log(ILI$unweighted_pct_ili), 
      pch=20, asp=1, col=1+(week>=40|week<=20))
abline(0,1)
title(
 "Unweighted % versus weighted % (logs, flu season red)")

par(mfrow=c(2,2))

plot (log(ILI$weighted_pct_ili), log(ILI$unweighted_pct_ili), 
      pch=20, asp=1, col=yrcols[year-2013])
abline(0,1)
title("ProxyU versus ProxyW (logs)")

plot (log(ILI$unweighted_pct_ili), log(proxyA), pch=20, asp=1, 
      col=yrcols[year-2013])
abline(0,1)
title("ProxyA versus ProxyU (logs)")

plot (log(ILI$weighted_pct_ili), log(proxyA), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyA versus ProxyW (logs)")

plot (log(proxyA), log(proxyB), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyB versus ProxyA (logs)")

plot (log(proxyA), log(proxyC), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyC versus ProxyA (logs)")

plot (log(proxyA), log(proxyD), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyD versus ProxyA (logs)")

plot (log(proxyC), log(proxyB), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyB versus ProxyC (logs)")

plot (log(proxyC), log(proxyD), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyD versus ProxyC (logs)")

par(mfrow=c(4,1))
plot_two_with_lines (start, 
  log(ILI$weighted_pct_ili), log(ILI$unweighted_pct_ili), 
  pch=19, ylab="Line goes to log(ProxyU)")
week_lines()
title("ProxyU versus ProxyW (logs)")

plot_two_with_lines (start, 
  log(ILI$unweighted_pct_ili), log(proxyA), 
  pch=19, ylab="Line goes to log(ProxyA)")
week_lines()
title("ProxyA versus ProxyU (logs)")

plot_two_with_lines (start, 
  log(ILI$weighted_pct_ili), log(proxyA), 
  pch=19, ylab="Line goes to log(ProxyA)")
week_lines()
title("ProxyA versus ProxyW (logs)")

plot_two_with_lines (start, 
  log(proxyA), log(proxyB), 
  pch=19, ylab="Line goes to log(ProxyB)")
week_lines()
title("ProxyB versus ProxyA (logs)")

plot_two_with_lines (start, 
  log(proxyA), log(proxyC), 
  pch=19, ylab="Line goes to log(ProxyC)")
week_lines()
title("ProxyC versus ProxyA (logs)")

plot_two_with_lines (start, 
  log(proxyA), log(proxyD), 
  pch=19, ylab="Line goes to log(ProxyD)")
week_lines()
title("ProxyD versus ProxyA (logs)")

plot_two_with_lines (start, 
  log(proxyC), log(proxyB), 
  pch=19, ylab="Line goes to log(ProxyB)")
week_lines()
title("ProxyB versus ProxyC (logs)")

plot_two_with_lines (start, 
  log(proxyC), log(proxyD), 
  pch=19, ylab="Line goes to log(ProxyD)")
week_lines()
title("ProxyD versus ProxyC (logs)")

# PLOT ANOMALOUS POINTS IN PROXIES.

par(mfrow=c(2,1))

plot (start, anomalous(log(ILI$weighted_pct_ili)), 
      pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(weighted percent ILI)")

plot (start, anomalous(log(ILI$unweighted_pct_ili)), 
      pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(unweighted percent ILI)")

plot (start, anomalous(log(proxyA)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyA)")

plot (start, anomalous(log(proxyB)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyB)")

plot (start, anomalous(log(proxyC)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyC)")

plot (start, anomalous(log(proxyD)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyD)")


# PRINT STANDARD DEVIATIONS OF LOG PROXIES.

cat("\nStandard deviations of log proxies:\n\n")

print (round (c (W=sd(ILI$weighted_pct), U=sd(ILI$unweighted_pct),
                 A=sd(proxyA), B=sd(proxyB), C=sd(proxyC), D=sd(proxyD)), 4))


# WRITE A FILE WITH THE VARIOUS PROXIES FOR ILI.

ILIproxy <- data.frame (start = as.character(start), year=year, week=week)

ILIproxy$proxyA <- proxyA
ILIproxy$proxyB <- proxyB
ILIproxy$proxyC <- proxyC
ILIproxy$proxyD <- proxyD
ILIproxy$proxyU <- ILI$unweighted_pct
ILIproxy$proxyW <- ILI$weighted_pct

write.table (ILIproxy, "ILIproxy.csv", sep=",",
             quote=FALSE, row.names=FALSE, col.names=TRUE)

# ALL DONE.

dev.off()
