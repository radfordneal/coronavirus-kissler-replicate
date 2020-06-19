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

start_season <- 40 # Week of start of "flu season", according to CDC
end_season <- 20   # Week of end of "flu season" (in year following start)
                   # (except ends week earlier in 2015, since 2014 has 53 weeks)
ILI$season_week <- 
  ifelse (ILI$week>=start_season, ILI$week-start_season+1,
          ILI$week+53-start_season+(ILI$year==2015))

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

yrcont <- (0:(length(start)-1)) / (365.24/7) # Continuous week over whole period

source("../util/util.R")

Christmas_indicator1 <- Christmas_indicator2 <- Christmas_indicator

Christmas_indicator1 [year>2016] <- 0
Christmas_indicator1 [year==2016] <- 0.5 * Christmas_indicator [year==2017]

Christmas_indicator2 [year<2016] <- 0
Christmas_indicator2 [year==2016] <- 0.5 * Christmas_indicator [year==2017]



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
# abline(h=c(13,13.5,14.0))

plot (start, ILI$ili_visits, pch=20); week_lines()
title("ILI physician visits")
plot (start, log(ILI$ili_visits), ylim=c(8.5,11.5), pch=20); week_lines()
title("Log ILI physician visits")
# abline(h=c(8.5,9.0,9.5))

par(mfrow=c(2,1))

plot_vs_doy (start, log(ILI$nonili_visits), pch=20)
abline(h=c(13.0,13.5,14.0,14.5),col="gray")
title("Log non-ILI visits vs day of year from July 1")

plot_vs_doy (start, log(ILI$ili_visits), pch=20)
abline(h=c(8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0),col="gray")
title("Log ILI visits vs day of year from July 1")


# PLOT ILI VISITS OF DIFFERENT AGE GROUPS.

stopifnot (ILI$age_0_4+ILI$age_5_24+ILI$age_25_64+ILI$age_65_up 
            == ILI$ili_visits)  # Check that all ages add up to total

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

par(mfrow=c(2,1))

plot (start, rep(0,length(start)), ylim=c(0,0.5), yaxs="i", type="n",
      ylab="Fraction of ILI visits in each age group")
week_lines()
lines (start, ILI$age_0_4 / ILI$ili_visits, col="red", lwd=3)
lines (start, ILI$age_5_24 / ILI$ili_visits, col="green", lwd=3)
lines (start, ILI$age_25_64 / ILI$ili_visits, col="blue", lwd=3)
lines (start, ILI$age_65_up / ILI$ili_visits, col="black", lwd=3)
title ("ILI visits by age: red:0-4, green:5-24, blue:25-64, black:65+  ")

plot (start, rep(0,length(start)), ylim=c(0.17,0.55), type="n",
      ylab="Fraction of ILI visits in each age group")
week_lines()
lines (start, ILI$age_5_24 / ILI$ili_visits, col="green", lwd=3)
lines (start, (ILI$age_25_64+ILI$age_65_up) / ILI$ili_visits, col="blue", lwd=3)
lines (start, ILI$age_0_4 / ILI$ili_visits, col="red", lwd=3)
title ("ILI visits by age: red:0-4, green:5-24, blue:25+")

par(mfrow=c(4,1))
plot (start, log (ILI$age_65_up / ILI$age_25_64), ylim=c(-2.2,-0.8), pch=20,
      ylab="Log (age 65+ / age 25-64)")
week_lines()
title ("Ratios of ILI by age group")
plot (start, log (ILI$age_25_64 / ILI$age_5_24), ylim=c(-0.8,0.6), pch=20,
      ylab="Log (age 25-64 / age 5-24)")
week_lines()
plot (start, log (ILI$age_25_64 / ILI$age_0_4), ylim=c(-0.4,1.0), pch=20,
      ylab="Log (age 25-64 / age 0-4)")
week_lines()
plot (start, log (ILI$age_5_24 / ILI$age_0_4), ylim=c(-0.4,1.0), pch=20,
      ylab="Log (age 5-24 / age 0-4)")
week_lines()


# PLOT WEIGHTED AND UWEIGHTED ILI PERCENTAGE, AND RATIO.  Also repeat the
# plots of numbers of providers, for reference.

par(mfrow=c(2,1))

proxyW <- ILI$weighted_pct_ili
plot (start, proxyW, pch=20); week_lines()
title("Weighted percent ILI (ProxyW)")
plot (start, log(proxyW), pch=20); week_lines()

proxyU <- ILI$unweighted_pct_ili
plot (start, proxyU, pch=20); week_lines()
title("Unweighted percent ILI (ProxyU)")
plot (start, log(proxyU), pch=20); week_lines()

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

nonili_visits_spline_df <- 9
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
title(paste("Residuals, std. dev. =",round(summary(visits_mod1)$sigma,4)))
abline(h=0)

# spline + providers + holidays

visits_mod2 <- lm (log(nonili_visits) 
                     ~ nonili_visits_spline + log(providers) 
                     + July4th_indicator
                     + Labor_Day_indicator 
                     + Thanksgiving_indicator
                     + Christmas_indicator1 + Christmas_indicator2
                     + New_Year_indicator
                   , data=ILI)
print(summary(visits_mod2))

par(mfrow=c(2,1))

plot (start, log(ILI$nonili_visits), pch=20); week_lines()
lines (start, predict(visits_mod2), col="blue", lwd=2)
title("Non-ILI model, spline+providers+holidays")

plot (start, residuals(visits_mod2), pch=20, ylab="Residuals"); week_lines()
title(paste("Residuals, std. dev. =",round(summary(visits_mod2)$sigma,4)))
abline(h=0)

# spline + providers + holidays + seasonality

nonili_xtra_df_prov <- 1  # df from providers
nonili_xtra_df_ind <- 6   # df from indicators
nonili_xtra_df <- 1 + nonili_xtra_df_prov + nonili_xtra_df_ind

visits_mod3 <- lm (log(nonili_visits) 
                     ~ nonili_visits_spline + log(providers)
                     + July4th_indicator
                     + Labor_Day_indicator 
                     + Thanksgiving_indicator
                     + Christmas_indicator1 + Christmas_indicator2
                     + New_Year_indicator
                     + season_week + I(season_week^2)
                     + sin(1*2*pi*yrcont) + cos(1*2*pi*yrcont)
                     + sin(2*2*pi*yrcont) + cos(2*2*pi*yrcont)
#                    + sin(3*2*pi*yrcont) + cos(3*2*pi*yrcont)
#                    + sin(4*2*pi*yrcont) + cos(4*2*pi*yrcont)
#                    + sin(5*2*pi*yrcont) + cos(5*2*pi*yrcont)
#                    + sin(6*2*pi*yrcont) + cos(6*2*pi*yrcont)
#                    + sin(7*2*pi*yrcont) + cos(7*2*pi*yrcont)
                   , data=ILI, x=TRUE)
print(summary(visits_mod3))

residuals_nonili_visits3 <- residuals(visits_mod3)

smooth_nonili_visits3 <- as.vector (exp (cbind (1,
   nonili_visits_spline, log(ILI$providers) 
  ) %*% coef(visits_mod3) [1:(nonili_visits_spline_df+2)]))

par(mfrow=c(2,1))

plot (start, log(ILI$nonili_visits), pch=20); week_lines()
lines (start, predict(visits_mod3), col="blue", lwd=2)
#lines (start, log(smooth_nonili_visits3), col="orange", lwd=2)
title("Non-ILI model, spline+providers+holidays+season")

plot (start, residuals(visits_mod3), pch=20, ylab="Residuals"); week_lines()
title(paste("Residuals, std. dev. =",round(summary(visits_mod3)$sigma,4)))
abline(h=0)

plot_vs_doy (start, residuals(visits_mod3), pch=20)
abline(h=0)
title("Residuals vs DOY (from July 1)")

plot (start, log(ILI$nonili_visits)-log(smooth_nonili_visits3), pch=20, ylab="")
week_lines(); abline(h=0)
title("Residuals without seasonal and holiday components")

xx <- visits_mod3$x
xx[,1:(nonili_visits_spline_df+nonili_xtra_df)] <- 0

seasonal_visits3 <- as.vector (xx %*% coef(visits_mod3))

plot (start, seasonal_visits3, type="l"); week_lines()
title("Seasonal component of model (log domain)")

xx <- visits_mod3$x
xx[,-(2:(nonili_visits_spline_df+1))] <- 0

trend_visits3 <- as.vector (xx %*% coef(visits_mod3))

plot (start, trend_visits3, type="l"); week_lines()
title("Trend component of model (log domain)")

xx <- visits_mod3$x
xx[,1:(1+nonili_visits_spline_df+nonili_xtra_df_prov)] <- 0
xx[,-(1:(nonili_visits_spline_df+nonili_xtra_df))] <- 0

holiday_visits3 <- as.vector (xx %*% coef(visits_mod3))

plot (start, holiday_visits3, pch=20); week_lines()
title("Holiday component of model (log domain)")


# MODEL THE NUMER OF ILI VISITS.

# Peak in 2014 coincides with Christmas - avoid it all being attributed to that.
Christmas_indicator1_ili <- Christmas_indicator1
Christmas_indicator1_ili[year==2014] <- 0.3*Christmas_indicator1_ili[year==2014]

ili_bound <- range(start)

ili_knots_within_year <- 
  as.Date(c("2000-07-01","2000-08-01","2000-09-01","2000-10-01","2000-11-01",
            "2000-12-01","2001-01-01","2001-02-01","2001-03-01","2001-04-01",
            "2001-05-01","2001-06-01")) -
  as.Date("2000-07-01")

ili_knots <- rep (seq (start[1],start[ILI$syear==2018][1],length=5),
                  each = length(ili_knots_within_year)) +
             ili_knots_within_year

ili_knots <- ili_knots[-1]

ili_xtra_df_prov <- 1  # df from providers
ili_xtra_df_ind <- 6   # df from indicators
ili_xtra_df <- 1 + ili_xtra_df_prov + ili_xtra_df_ind

ili_visits_spline_df <- length(ili_knots) + 3
ili_visits_spline <- bs (start, knots=ili_knots, Bound=ili_bound)

# spline + providers + holidays

ili_visits_mod3 <- lm (log(ili_visits) 
                     ~ log(providers)
                     + July4th_indicator
                     + Labor_Day_indicator 
                     + Thanksgiving_indicator
                     + Christmas_indicator1_ili + Christmas_indicator2
                     + New_Year_indicator
                     + ili_visits_spline
                  , data=ILI, x=TRUE)
print(summary(ili_visits_mod3))

par(mfrow=c(2,1))

plot (start, log(ILI$ili_visits), pch=20); week_lines()
lines (start, predict(ili_visits_mod3), col="blue", lwd=2)
title("ILI model, spline+providers+holidays")

plot (start, residuals(ili_visits_mod3), pch=20, ylab="Residuals")
title(paste("Residuals, std. dev. =",round(summary(ili_visits_mod3)$sigma,4)))
week_lines(); abline(h=0)

plot_vs_doy (start, residuals(ili_visits_mod3), pch=20)
abline(h=0)
title("Residuals vs DOY (from July 1)")

xx <- ili_visits_mod3$x
xx[,2:ili_xtra_df] <- 0
ili_spline3 <- as.vector (xx %*% coef(ili_visits_mod3))

plot(start,ili_spline3,type="l",col="darkgray",lwd=3)
points(start,ili_spline3+residuals(ili_visits_mod3),pch=20,ylab="")
week_lines()
title("Spline, and spline plus residuals")

# spline + providers + holidays + non-ili residuals

ili_visits_mod4 <- lm (log(ili_visits) 
                     ~ log(providers)
                     + July4th_indicator
                     + Labor_Day_indicator 
                     + Thanksgiving_indicator
                     + Christmas_indicator1_ili + Christmas_indicator2
                     + New_Year_indicator
                     + residuals_nonili_visits3
                     + ili_visits_spline
                  , data=ILI, x=TRUE)
print(summary(ili_visits_mod4))

par(mfrow=c(2,1))

plot (start, log(ILI$ili_visits), pch=20); week_lines()
lines (start, predict(ili_visits_mod4), col="blue", lwd=2)
title("ILI model, spline+providers+holidays+noniliresid")

plot (start, residuals(ili_visits_mod4), pch=20, ylab="Residuals")
title(paste("Residuals, std. dev. =",round(summary(ili_visits_mod4)$sigma,4)))
week_lines(); abline(h=0)

plot_vs_doy (start, residuals(ili_visits_mod4), pch=20)
abline(h=0)
title("Residuals vs DOY (from July 1)")

xx <- ili_visits_mod4$x
xx[,2:(ili_xtra_df+1)] <- 0
ili_spline4 <- as.vector (xx %*% coef(ili_visits_mod4))

plot(start,ili_spline4,type="l",col="darkgray",lwd=3)
points(start,ili_spline4+residuals(ili_visits_mod4),pch=20,ylab="")
week_lines()
title("Spline, and spline plus residuals")


# CREATE PROXYA FROM RATIO OF ILI VISITS TO NON-ILI VISITS.

proxyA <- 100 * ILI$ili_visits / ILI$nonili_visits

par(mfrow=c(4,1))

plot (start, proxyA, pch=20); week_lines()
title("ProxyA: 100 x ratio of ILI visits to non-ILI visits")

plot (start, log(proxyA), pch=20, ylim=c(-0.5,2.5))
week_lines()


# CREATE PROXYAX AND PROXYWX BY FUDGING HOLIDAYS IN PROXYA AND PROXYW.

holiday_fudge <- function (proxy)
{
  for (i in which(Thanksgiving_indicator!=0))
  { proxy[i] <- (1/2)*proxy[i-1] + (1/2)*proxy[i+1]
  }

  for (i in which(New_Year_indicator!=0))  # Note: Christmas is not binary
  { proxy[i-1] <- (2/3)*proxy[i-2] + (1/3)*proxy[i+1]
    proxy[i] <- (1/3)*proxy[i-2] + (2/3)*proxy[i+1]
  }

  proxy
}

proxyWX <- holiday_fudge (proxyW)
plot (start, proxyWX, pch=20); week_lines()
week_lines()
title("ProxyWX: ProxyW with holiday fudge")

plot (start, log(proxyWX), pch=20); week_lines()
week_lines()

proxyAX <- holiday_fudge (proxyA)
plot (start, proxyAX, pch=20); week_lines()
week_lines()
title("ProxyAX: ProxyA with holiday fudge")

plot (start, log(proxyAX), pch=20); week_lines()
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
# Uses ili_visits_mod3 to produce predictions, omitting variation in
# ILI due to the number of providers and holidays, which are
# considered to not relate to actual ILI incidence.  Retains the
# spline over time and the model residuals.  The result is divided by
# a scale factor just to get it to numerically match (approximately)
# the other proxies.

proxyC <- exp (ili_spline3+residuals(ili_visits_mod3)) / 80

plot (start, proxyC, pch=20); week_lines()
title("ProxyC: Derived from model of ILI visits")

plot (start, log(proxyC), pch=20, ylim=c(-0.5,2.5))
week_lines()


# CREATE PROXYD FROM ILI VISITS MINUS MODEL PREDICTIONS OF IRRELEVANT PART.  
# Uses ili_visits_mod4 to produce predictions, omitting variation in
# ILI due to the number of providers, the residuals from the model for
# non-ILI visits, and the holidays, since these are considered to not
# relate to actual ILI incidence.  Retains the spline over time, and
# the model residuals.  The result is divided by a scale factor just
# to get it to numerically match (approximately) the other proxies.

proxyD <- exp (ili_spline4+residuals(ili_visits_mod4)) / 85

plot (start, proxyD, pch=20); week_lines()
title("ProxyD: Derived from model of ILI visits, with non-ILI residuals")

plot (start, log(proxyD), pch=20, ylim=c(-0.5,2.5))
week_lines()


# CREATE PROXYE, LIKE PROXYD BUT WITHOUT INCLUDING RESIDUALS.

proxyE <- exp (ili_spline4) / 85

plot (start, proxyE, pch=20); week_lines()
title("ProxyE: ProxyD without ILI model residuals")

plot (start, log(proxyE), pch=20, ylim=c(-0.5,2.5))
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

plot (log(proxyW), log(proxyWX), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyWX versus ProxyW (logs)")

plot (log(proxyA), log(proxyAX), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyAX versus ProxyA (logs)")

plot (log(proxyA), log(proxyB), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyB versus ProxyA (logs)")

plot (log(proxyAX), log(proxyB), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyB versus ProxyAX (logs)")

plot (log(proxyA), log(proxyC), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyC versus ProxyA (logs)")

plot (log(proxyAX), log(proxyC), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyC versus ProxyAX (logs)")

plot (log(proxyA), log(proxyD), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyD versus ProxyA (logs)")

plot (log(proxyAX), log(proxyD), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyD versus ProxyAX (logs)")

plot (log(proxyC), log(proxyB), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyB versus ProxyC (logs)")

plot (log(proxyC), log(proxyD), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyD versus ProxyC (logs)")

plot (log(proxyE), log(proxyD), pch=20, asp=1,
      col=yrcols[year-2013])
abline(0,1)
title("ProxyD versus ProxyE (logs)")

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
  log(proxyW), log(proxyWX), 
  pch=19, ylab="Line goes to log(ProxyWX)")
week_lines()
title("ProxyWX versus ProxyW (logs)")

plot_two_with_lines (start, 
  log(proxyA), log(proxyAX), 
  pch=19, ylab="Line goes to log(ProxyAX)")
week_lines()
title("ProxyAX versus ProxyA (logs)")

plot_two_with_lines (start, 
  log(proxyA), log(proxyB), 
  pch=19, ylab="Line goes to log(ProxyB)")
week_lines()
title("ProxyB versus ProxyA (logs)")

plot_two_with_lines (start, 
  log(proxyAX), log(proxyB), 
  pch=19, ylab="Line goes to log(ProxyB)")
week_lines()
title("ProxyB versus ProxyAX (logs)")

plot_two_with_lines (start, 
  log(proxyA), log(proxyC), 
  pch=19, ylab="Line goes to log(ProxyC)")
week_lines()
title("ProxyC versus ProxyA (logs)")

plot_two_with_lines (start, 
  log(proxyAX), log(proxyC), 
  pch=19, ylab="Line goes to log(ProxyC)")
week_lines()
title("ProxyC versus ProxyAX (logs)")

plot_two_with_lines (start, 
  log(proxyA), log(proxyD), 
  pch=19, ylab="Line goes to log(ProxyD)")
week_lines()
title("ProxyD versus ProxyA (logs)")

plot_two_with_lines (start, 
  log(proxyAX), log(proxyD), 
  pch=19, ylab="Line goes to log(ProxyD)")
week_lines()
title("ProxyD versus ProxyAX (logs)")

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

plot_two_with_lines (start, 
  log(proxyD), log(proxyE), 
  pch=19, ylab="Line goes to log(ProxyE)")
week_lines()
title("ProxyE versus ProxyD (logs)")

plot_two_with_lines (start, 
  log(proxyD), log(proxyW), 
  pch=19, ylab="Line goes to log(ProxyW)")
week_lines()
title("ProxyW versus ProxyD (logs)")

plot_two_with_lines (start, 
  log(proxyD), log(proxyAX), 
  pch=19, ylab="Line goes to log(ProxyAX)")
week_lines()
title("ProxyAX versus ProxyD (logs)")


# PLOT ANOMALOUS POINTS IN PROXIES.

par(mfrow=c(2,1))

plot (start, anomalous(log(proxyW)), 
      pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyW)")

plot (start, anomalous(log(proxyWX)), 
      pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyWX)")

plot (start, anomalous(log(proxyU)), 
      pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyU)")

plot (start, anomalous(log(proxyA)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyA)")

plot (start, anomalous(log(proxyAX)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyAX)")

plot (start, anomalous(log(proxyB)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyB)")

plot (start, anomalous(log(proxyC)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyC)")

plot (start, anomalous(log(proxyD)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyD)")

plot (start, anomalous(log(proxyE)), pch=20, ylim=c(-0.2,0.22))
week_lines(); abline(h=c(0,-0.1,-0.2,0.1,0.2))
title("Anomalous points in log(proxyE)")


# PRINT STANDARD DEVIATIONS OF LOG PROXIES.

cat("\nStandard deviations of log proxies:\n\n")

print (round (c (W=sd(ILI$weighted_pct), U=sd(ILI$unweighted_pct),
                 A=sd(proxyA), B=sd(proxyB), C=sd(proxyC), D=sd(proxyD)), 4))


# WRITE A FILE WITH THE VARIOUS PROXIES FOR ILI.

ILIproxy <- data.frame (start = as.character(start), year=year, week=week)

ILIproxy$proxyW <- proxyW
ILIproxy$proxyWX <- proxyWX
ILIproxy$proxyU <- proxyU
ILIproxy$proxyA <- proxyA
ILIproxy$proxyAX <- proxyAX
ILIproxy$proxyB <- proxyB
ILIproxy$proxyC <- proxyC
ILIproxy$proxyD <- proxyD
ILIproxy$proxyE <- proxyE


write.table (ILIproxy, "ILIproxy.csv", sep=",",
             quote=FALSE, row.names=FALSE, col.names=TRUE)

# ALL DONE.

dev.off()
