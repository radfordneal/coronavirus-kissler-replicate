# Code simulate data according to a model fit to estimated Rt values, read
# from the R-model directory.
#
# Options are specified by arguments in the R command after --args:
#
#   - The R estimates to use, from corresponding file (required)
#   - Model of "immunity" (required) - i2 (exp decay), i3 (short & long),
#                           i4 (like i3 but with no short-term cross-immunity)
#   - Model of seasonal effect (required) - e2 (sine), e3 (Fourier)
#   - Whether heteroskedasticity w.r.t. virus is modelled - het for "yes" 
#     (default "no")
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

het_virus <- sum(args=="het") == 1
if (het_virus) args <- args [! (args %in% "het")]

R_estimates <- args

stopifnot(length(R_estimates)==1)

file_base <- paste0 (R_estimates,"-Rt-s2-",immune_type,"-",seffect_type,
                     if (het_virus) "-het")


# PLOT SETUP.

pdf (paste0 ("reg-sim-",gsub("Rt-s2-","",file_base),".pdf"), height=8, width=6)
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


# INCLUDE FUNCTIONS THAT COMPUTE SEASONAL EFFECTS.

source("../R-model/steffect.R")


# ----- DO EVERYTHING FOR BOTH ALPHACORONAVIRUSES AND BETACORONAVIRUSES -----

yrs <- R_est$yrs
trend_spline <- make_trend_spline (yrs, yrs)

yrsd <- rep(yrs,each=7) + (0:6)/365.24 - 3.5/365.24

for (g in seq_along (virus_groups)) {

  virus_group <- virus_groups[[g]]


# READ THE MODEL.

r <- readRDS (paste0("../R-model/R-model-",file_base,"-",
                     names(virus_groups)[g],".model"))

model <- r$model
imm_decay <- r$imm_decay
ltimm_decay <- r$ltimm_decay

cat("\nMODEL FOR GROUP",names(virus_groups)[g],":",virus_group,"\n\n")
cat("Coefficients:\n\n")
print(as.matrix(coef(model)))
cat("\n")
cat("Immune decay:\n")
print(imm_decay[virus_group])
cat("\n")
cat("Long-term immune decay:\n")
print(ltimm_decay[virus_group])
cat("\n")


# SET GENERATION INTERVAL DISTRIBUTION USED BY KISSLER, ET AL.

SARS_shape <- 2.35
SARS_scale <- 9.48

SARS_gen_interval <-
  dweibull (1:ceiling(qweibull(0.99, shape=SARS_shape, scale=SARS_scale)),
            shape=SARS_shape, scale=SARS_scale)

gen_interval <- SARS_gen_interval # / sum(SARS_gen_interval)
rev_gen_interval <- rev(gen_interval)


# PROXIES FOR VIRUS INCIDENCE.

proxy1 <- R_est[,paste0(virus_group[1],"_proxy")]
proxy2 <- R_est[,paste0(virus_group[2],"_proxy")]


# FUNCTION TO RUN SIMULATIONS.  Returns 'nsims' simulation results, each
# found by a series of five-year simulations with 'warmup' simulations before 
# the final one that is returned.
#
# Simulations are done day-by-day (not week-by-week), using the trend spline
# and seasonal effects, and the long and short term immunity based on previous
# incidence.
#
# Initial values for the history of incidence are taken from the end
# of one of the seasons in the preceding five-year simulation (and set
# somewhat arbitrarily for the first five-year simulation), the same
# one for all simulations done. Some randomness is introduced into the
# Rt values, differently for each simulation.
#
# 'P' is a list of parameter values affecting the simulation.

run_sims <- function (nsims, warmup, 
                      P = list(imm_decay = imm_decay, ltimm_decay = ltimm_decay,
                               mc = coef(model), Rt_noise_sd = 0.05,
                               Rt_offset_sd = 0.05, Rt_offset_alpha = 0.9))
{
  daily_decay <- P$imm_decay ^ (1/7)
  ltdaily_decay <- P$ltimm_decay ^ (1/7)
  mc <- P$mc

  # Pre-compute trend and seasonal effects for all days.
  
  tn <- ncol(trend_spline)
  tseff <-
  ( if (seffect_type=="e2")
      seffect_e2(yrsd,mc)
    else
      seffect_e3(yrsd,mc) + as.vector (predict(trend_spline,yrsd) %*% mc[1:tn])
  )

  # Initial levels for exponentially-decaying past incidence.  Set on the
  # assumption that incidence was low before the start.
  #
  #   t      matrix of short-term exponetial sums for each virus, 2 x nsims
  #   lt     matrix of long-term exponetial sums for each virus, 2 x nsims
  #   past   list of matrices of assumed past incidence values, one matrix
  #          for each virus, dimension length(gen_interval) x nsims

  q <- c (quantile(proxy1,0.1), quantile(proxy2,0.1))
  t <- matrix (q / (1-imm_decay[virus_group]), nrow=2, ncol=nsims)
  tlt <- matrix (q / (1-ltimm_decay[virus_group]), nrow=2, ncol=nsims)
  past <- list (matrix (q[1], nrow=length(gen_interval), ncol=nsims),
                matrix (q[2], nrow=length(gen_interval), ncol=nsims))

  # Space to store simulation results. A list of two matrices, one for each
  # virus, of dimension total number of days x nsims.

  sims <- rep (list (matrix (0, nrow=7*length(start), ncol=nsims)), times=2)

  # Stuff for saving a state for use to initialize the next simulation.
  
  wsave <- sample(rep(1:5,length=warmup+1))
  sv_past <- NULL

  p <- vector("list",2)
  
  for (w in 1:(warmup+1))
  {
    # Do next simulations of a five-year period.
  
    Rt_offset <- rnorm(nsims,0,P$Rt_offset_sd) # initialize AR(1) process that
                                               #   modifies modelled Rt values
  
    for (day in 1:(7*length(start)))
    {
      Rt_offset <- P$Rt_offset_alpha * Rt_offset +
                   sqrt(1-P$Rt_offset_alpha^2) * rnorm(nsims,0,P$Rt_offset_sd)

      for (vi in 1:2)
      { 
        virus <- virus_group[vi]

        # Compute vector of nsims log R values based on trend, seasonality, 
        # and immunity due to past infections.

        log_Rt <- tseff[day] +
                  mc [paste0(virus,"_overall")] +
                  mc [paste0(virus,"_same")] * t [vi,]
        if (immune_type!="i4")
        { log_Rt <- log_Rt + 
                    mc [paste0(virus,"_other")] * t [if (vi==1) 2 else 1,]
        }
        if (immune_type=="i3" || immune_type=="i4")
        { log_Rt <- log_Rt + 
                    mc [paste0(virus,"_samelt")] * tlt [vi,] +
                    mc [paste0(virus,"_otherlt")] * tlt [if (vi==1) 2 else 1,]
        }

        inf <- colSums (past[[vi]] * rev_gen_interval)

        p[[vi]] <- inf * exp (log_Rt + Rt_offset + rnorm(nsims,0,P$Rt_noise_sd))

        past[[vi]] <- rbind (past[[vi]][-1,], p[[vi]])
        sims[[vi]][day,] <- p[[vi]]
  
        t[vi,] <- p[[vi]] + t[vi,]*daily_decay[virus]
        tlt[vi,] <- p[[vi]] + tlt[vi,]*ltdaily_decay[virus]
      }

      # Save end of one of the years to initialize next five-year simulation.

      if ((day+5) %% 365 == 0)
      { wsave[w] <- wsave[w] - 1
        if (wsave[w] == 0)
        { sv_past <- past
          sv_t <- t
          sv_tlt <- tlt
        }
      }
    }
  
    # Set up initial state for next five-year simulations. Randomized a bit,
    # separately for each simulation.
  
    for (vi in 1:2)
    { n <- exp(rnorm(nsims,0,0.2))
      past[[vi]] <- sv_past[[vi]] * rep (n, each=nrow(sv_past[[vi]]))
      t[vi,] <- sv_t[vi,] * n
      tlt[vi,] <- sv_tlt[vi,] * n
    }
  }
  
  # Return results of final simulations, cnverted to weekly values. The
  # returned value is a list of, for each virus, a matrix with dimensions
  # number of weeks x nsims.
  
  w <- seq (4, 7*length(start), by=7)
  wsims <- vector("list",2)
  for (vi in 1:2)
  { wsims[[vi]] <- matrix(NA,length(w),nsims)
    for (s in 1:nsims) wsims[[vi]][,s] <- filter (sims[[vi]][,s], rep(1,7)) [w]
  }

  wsims
}  
  
  
# DO THE SIMULATIONS.

set.seed(1)

warmup <- 10
n_plotted <- 7

wsims <- run_sims (n_plotted, warmup)


# PLOT THE OBSERVED INCIDENCE, THEN SAVED SIMULATIONS.  Plots all use the
# same scales.

par(mfrow=c(4,1))

ylim <-  max (proxy1, proxy2, wsims[[1]], wsims[[2]])

plot (start, rep(0,length(start)),
      ylim=c(0,1.02*ylim), yaxs="i", type="n", ylab="Incidence proxy")

lines (start, proxy1, col="blue")
lines (start, proxy2, col="red")

title (paste
 ("Observed proxies for",virus_group[1],"(blue) and",virus_group[2],"(red)"))

for (s in 1:n_plotted)
{
  plot (start, rep(0,length(start)),
        ylim=c(0,1.02*ylim), yaxs="i", type="n", ylab="Simulated incidence")

  lines (start, wsims[[1]][,s], col="blue")
  lines (start, wsims[[2]][,s], col="red")

  if (s==1) title (paste ("Simulations of",virus_group[1],"and",virus_group[2]))
}

# ----- END OF LOOP OVER THE TWO VIRUS GROUPS -----

}

# ALL DONE.

dev.off()
