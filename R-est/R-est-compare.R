# Script to compare R estimates produced from different proxies or
# different methods.
#
# Takes two names of R estimates as argument to the R command, after --args.
# Estimates are read from the files R-est-<est1>.csv and R-est-<est2>.csv.
# Plots are written to R-est-compare-<est1>:<est2>.pdf.
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

args <- commandArgs(trailing=TRUE)

stopifnot (length(args)==2)

arg1 <- args[1]
arg2 <- args[2]

est1 <- read.csv (paste0("R-est-",arg1,".csv"), 
                  header=TRUE, stringsAsFactors=FALSE)
est1$start <- as.Date(est1$start)

est2 <- read.csv (paste0("R-est-",arg2,".csv"),
                  header=TRUE, stringsAsFactors=FALSE)
est2$start <- as.Date(est2$start)

stopifnot (identical (names(est1), names(est2)))
stopifnot (nrow(est1) == nrow(est2))

pdf (paste0 ("R-est-compare-",arg1,":",arg2,".pdf"), height=8, width=6)
par(mar=c(1.5,2.3,3,0.5),mgp=c(1.4,0.3,0),tcl=-0.22)

start <- est1$start
year <- est1$year
week <- est1$week

source("../util/util.R")

minR <- 0.1

par(mfrow=c(2,1))

for (estimate_type in c("Rt","Ru","Ru_smoothed"))
{
  for (virus in viruses)
  { plot_two_with_lines (start, 
      log(pmax(minR,est1[,paste0(virus,"_",estimate_type)])), 
      log(pmax(minR,est2[,paste0(virus,"_",estimate_type)])),
      pch=19, ylab=paste0("Line goes to log(",arg2,")"))
    week_lines()
    title (paste(virus,estimate_type,"-",arg1,":",arg2))
  }
}
