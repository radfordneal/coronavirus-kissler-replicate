# Code simulate data according to a model fit to estimated Rt values, read
# from the R-model directory.
#
# Options are specified by arguments in the R command after --args:
#
#   - The R estimates to use, from corresponding file (required)
#   - Model of "immunity" (required) - i2 (exp decay), i3 (short & long),
#                           i4 (like i3 but with no short-term cross-immunity)
#   - Model of seasonal effect (required) - e2 (sine), e3 (Fourier)
#   - Immune decay constants (defaults in code).
#     For example: decay:NL63=0.9,E229=0.8,OC43=0.7,HKU1=0.6
#     May override just a subst of the values
#   - Long-term immune decay constants for i3/i4 model (defaults in code).
#     For example: ltdecay:NL63=0.91,E229=0.92,OC43=0.93,HKU1=0.94
#     May override just a subst of the values
#
# Produces various plots that are written to the file with name
# reg-sim-<R-estimates>-<R-estimate-type>-<in>-<en>.pdf.  
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

options(warn=1)


# ESTABLISH WHICH PARAMETERS TO USE, LOOKING AT R'S ARGUMENTS.

args <- commandArgs(trailing=TRUE)

cat("ARGUMENTS:",args,"\n\n")

getarg <- function (what)
{ stopifnot (sum (args %in% what) == 1)
  res <- args [args %in% what]
  args <<- args [! (args %in% what)]
  res
}

immune_type <- getarg (c("i2","i3", "i4"))
seffect_type <- getarg (c("e2","e3"))

R_estimates <- args

stopifnot(length(R_estimates)==1)

file_base <- paste0 (R_estimates,"-Rt-s2-",immune_type,"-",seffect_type)


# PLOT SETUP.

pdf (paste0 ("reg-sim-",file_base,".pdf"), height=8, width=6)
par(mar=c(1.5,2.3,3,0.5), mgp=c(1.4,0.3,0), tcl=-0.22)
yrcols <- c("red","green","blue","orange","darkcyan","darkmagenta")


# READ THE FILE OF CORONAVIRUS R ESTIMATES (INCLUDES PROXY).

R_est <- read.csv (paste0("../R-est/R-est-",R_estimates,".csv"),
                   header=TRUE, stringsAsFactors=FALSE)

R_est$start <- as.Date(R_est$start)
R_est$yrs <- (0:(nrow(R_est)-1)) / (365.24/7)


# INCLUDE UTILITY FUNCTIONS.  They need "start" and "week" to be defined.

start <- R_est$start  # Start dates of weeks being analysed
year  <- R_est$year   # Year for each week
week  <- R_est$week   # Number of each week in its year

source("../util/util.R")


# ----- DO EVERYTHING FOR BOTH ALPHACORONAVIRUSES AND BETACORONAVIRUSES -----

for (g in seq_along (virus_groups)) {

  virus_group <- virus_groups[[g]]


# READ THE MODEL.

model <- readRDS (paste0("../R-model/R-model-",file_base,"-",
                         names(virus_groups)[g],".model"))

print(model)

# ----- END OF LOOP OVER THE TWO VIRUS GROUPS -----

}

# ALL DONE.

dev.off()
