# Code to reproduce the estimation of R values for common cold coronaviruses 
# as in the Kissler, et al paper, and to investigate alternative methods.
#
# The coronavirus proxy to use is identified by the R command argument
# after --args (default proxyW, which is the proxy used in the
# Kissler, et al paper).  Use of the daily proxies for coronaviruses
# is specified with the "daily" argument (default is to use weekly proxies,
# as with Kissler, et al).  Filtering of interpolated daily values is 
# specified by "filter" (default no filtering, as with Kissler, et al).  
# Note that "filter" makes no sense in conjunction with "daily".
#
# Produces various plots, written to R-est-<proxy>[-filter].pdf or to
# R-est-<proxy>-daily.pdf.  Weekly estimates for R are written to
# R-est-<proxy>[-filter].csv or to R-est-<proxy>-daily.csv, with columns 
# named <virus>_<est>, with <est> being Rt, Ru, Rt_smoothed, or Ru_smoothed
# (smoothed versions being 21-day averages, rather than 7-day averages).
# The estimate used by Kissler, et al is Ru_smoothed.  The proxy used for 
# each virus is also recorded (as read from the ILI proxy file), as a 
# column named <virus>_<proxy>; for daily proxies, this is the proxy value
# at the start of the week (Sunday).
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


# ESTABLISH WHICH PARAMETERS TO USE, LOOKING AT R'S ARGUMENTS.
# Includes proxy to use and generation interval distribution.

proxy <- "proxyW"       # Defaults as in Kissler, et al
filter <- FALSE
dailyp <- FALSE

args <- commandArgs(trailing=TRUE)

if (any(args=="filter"))
{ args <- args[args!="filter"]
  filter <- TRUE
}

if (any(args=="daily"))
{ args <- args[args!="daily"]
  dailyp <- TRUE
  stopifnot(!filter)
}

if (length(args)>0) 
{ stopifnot(length(args)==1)
  proxy <- args
}

# Generation interval distribution used by Kissler, et al.

SARS_shape <- 2.35
SARS_scale <- 9.48

SARS_gen_interval <- 
  dweibull (1:ceiling(qweibull(0.99, shape=SARS_shape, scale=SARS_scale)),
            shape=SARS_shape, scale=SARS_scale)
print(round(SARS_gen_interval,4))

gen_interval <- SARS_gen_interval


# PLOT SETUP.

pdf (
  paste0 ("R-est-", proxy, if (filter) "-filter", if (dailyp) "-daily", ".pdf"),
  height=8,width=6)
par(mar=c(1.5,2.3,3,0.5),mgp=c(1.4,0.3,0),tcl=-0.22)
yrcols <- c("red","green","blue","orange","darkcyan","darkmagenta")


# READ THE FILE OF CORONAVIRUS INCIDENCE PROXIES.

CoVproxy <- read.csv ("../proxy/CoVproxy.csv", 
                      header=TRUE, stringsAsFactors=FALSE)
CoVproxy$start <- as.Date(CoVproxy$start)

if (dailyp)
{ CoVproxy_daily <- read.csv ("../proxy/CoVproxy-daily.csv",
                              header=TRUE, stringsAsFactors=FALSE)
  CoVproxy_daily$date <- as.Date(CoVproxy_daily$date)
  start_week <- seq(1,nrow(CoVproxy_daily)-6,by=7)
  stopifnot (all (CoVproxy$start == CoVproxy_daily$date[start_week]))
  CoVproxy <- cbind (CoVproxy [, c("start","year","week")],
                     CoVproxy_daily [start_week, ])
  CoVproxy$date <- NULL
} 


# INCLUDE UTILITY FUNCTIONS.  They need "start" and "week" to be defined.

start <- CoVproxy$start  # Start dates of weeks being analysed
year  <- CoVproxy$year   # Year for each week
week  <- CoVproxy$week   # Number of each week in its year

source("../util/util.R")


# INTERPOLATE WEEKLY INCIDENCE TO DAILY INCIDENCE.  Uses the spline method 
# described in te Beest, et al, as is done by Kissler, et al.

interpolate_daily <- function (weekly)
{ 
  cum_weekly <- c(0,cumsum(weekly))
  day_weekly <- seq(0,7*length(weekly),by=7)
  cum_daily <- spline (day_weekly, cum_weekly, method="fmm",
                       xout=seq(0,7*length(weekly))) $ y
  
  daily <- diff(cum_daily)

  if (filter) 
  { daily <- as.vector ( # discard ts class
               filter (daily, c(0.04,0.07,0.10,0.17,0.24,0.17,0.10,0.07,0.04)))
  }

  daily [daily<0] <- 0.0000001

  daily
}


# ESTIMATE DAILY VALUES FOR R, REFERENCING TRANSMISSION TIME.  
# See Wallinga and Lipsitch, equation (4.1).

estimate_Rt <- function (daily,gen_interval)
{ Rt <- rep(NA,length(daily))
  a <- 1:length(gen_interval)
# for (t in (length(gen_interval)+1) : length(daily))
  for (t in (length(gen_interval)+2) : length(daily)) # +2: match Kissler, et al
  { Rt[t] <- daily[t] / sum (daily[t-a]*gen_interval)
  }
  Rt
}


# ESTIMATE DAILY VALUES FOR R, REFERENCING PRIMARY CASE INFECTION TIME. 
# See Wallinga and Lipsitch, equation (4.2).

estimate_Ru <- function (Rt,gen_interval)
{ Ru <- rep(NA,length(Rt))
  a <- 1:length(gen_interval)
  for (u in 1 : (length(Rt)-length(gen_interval)))
  { Ru[u] <- sum (gen_interval * Rt[u+a])
  }
  Ru
}


# CONVERT DAILY VALUES TO WEEKLY VALUES.  Simply takes the geometric
# average of the values for the seven days in the week.  Weekly values
# are NA if any of the daily values for the week are NA.

weekly_from_daily <- function (daily)
{ weekly <- rep(NA,length(daily)/7)
  for (i in seq_along(weekly))
  { weekly[i] <- exp (mean (log (daily [ (7*(i-1)+1) : (7*i) ])))
  }
  weekly
}


# CONVERT DAILY VALUES TO SMOOTHED WEEKLY VALUES.  Takes the geometric
# average of the values for the seven days in the week and the seven
# days before and seven days after the week, with NA values ommited.
# Weekly values are NA only if all 21 daily values averaged are NA.

smoothed_weekly_from_daily <- function (daily)
{ 
  weekly <- rep(NA,length(daily)/7)
  for (i in seq_along(weekly))
  { srt <- max (1, 7*(i-2)+1)
    end <- min (length(daily), 7*(i+1))
    weekly[i] <- ( if (all (is.na (daily[srt:end]))) NA 
                   else exp (mean (log (daily[srt:end]), na.rm=TRUE)))
  }

  weekly[1:3] <- NA                                # Considered unreliable
  weekly[(length(weekly)-2):length(weekly)] <- NA  # by Kissler, et al

  weekly
}


# RUN EVERYTHING FOR EACH VIRUS.

par(mfrow=c(2,1))

R_est <- data.frame (start=start, year=year, week=week)

for (virus in viruses)
{
  virus_proxy <- paste0(virus,"_",proxy)
  weekly <- CoVproxy[,virus_proxy]
  R_est[,paste0(virus,"_proxy")] <- weekly

  daily <- 
    if (dailyp) CoVproxy_daily[,virus_proxy] else interpolate_daily(weekly)
  plot(daily,pch=20,xlab="",ylab="")
  abline(h=0)
  title(paste("Daily proxy for",virus,if(!dailyp)"(interpolated)"))
  plot(log(daily),pch=20,xlab="",ylab="")
  title(paste("Log daily proxy for",virus,if(!dailyp)"(interpolated)"))
  cat("Number of non-missing daily values:",sum(!is.na(daily)),"\n")

  Rt_daily <- estimate_Rt(daily,gen_interval)
  Ru_daily <- estimate_Ru(Rt_daily,gen_interval)
  plot(Rt_daily,pch=20,xlab="",ylab="",ylim=c(0,12.5),yaxs="i")
  abline(h=1)
  title(paste("Daily estimates of Rt for",virus))
  plot(Ru_daily,pch=20,xlab="",ylab="",ylim=c(0,12.5),yaxs="i")
  abline(h=1)
  title(paste("Daily estimates of Ru for",virus))

  Rt_weekly <- weekly_from_daily(Rt_daily)
  Ru_weekly <- weekly_from_daily(Ru_daily)
  Rt_weekly_smoothed <- smoothed_weekly_from_daily(Rt_daily)
  Ru_weekly_smoothed <- smoothed_weekly_from_daily(Ru_daily)
  plot(start,Rt_weekly,pch=20,xlab="",ylab="",ylim=c(0,12.5),yaxs="i")
  week_lines(); abline(h=1); abline(h=seq(2,20,by=2),col="gray")
  title(paste("Weekly estimates of Rt for",virus))
  plot(start,Ru_weekly,pch=20,xlab="",ylab="",ylim=c(0,12.5),yaxs="i")
  week_lines(); abline(h=1); abline(h=seq(2,20,by=2),col="gray")
  title(paste("Weekly estimates of Ru for",virus))
  plot(start,Rt_weekly_smoothed,pch=20,xlab="",ylab="",ylim=c(0,12.5),yaxs="i")
  week_lines(); abline(h=1); abline(h=seq(2,20,by=2),col="gray")
  title(paste("Smoothed weekly estimates of Rt for",virus))
  plot(start,Ru_weekly_smoothed,pch=20,xlab="",ylab="",ylim=c(0,12.5),yaxs="i")
  week_lines(); abline(h=1); abline(h=seq(2,20,by=2),col="gray")
  title(paste("Smoothed weekly estimates of Ru for",virus))

  R_est[,paste0(virus,"_Rt")] <- Rt_weekly
  R_est[,paste0(virus,"_Ru")] <- Ru_weekly
  R_est[,paste0(virus,"_Rt_smoothed")] <- Rt_weekly_smoothed
  R_est[,paste0(virus,"_Ru_smoothed")] <- Ru_weekly_smoothed
}


# WRITE OUT THE ESTIMATES OF R FOR ALL VIRUSES.

R_est$start <- as.character(R_est$start,year=year,week=week)

write.table (R_est, 
  paste0 ("R-est-", proxy, if (filter) "-filter", if (dailyp) "-daily", ".csv"),
  sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

# ALL DONE.

dev.off()
