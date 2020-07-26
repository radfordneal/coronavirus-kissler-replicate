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
#   - Whether heteroskedasticity w.r.t. virus is modelled - het for "yes" 
#     (default "no")
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


# The simulation is done multiple times for the five seasons, using the
# trend spline and seasonal effects.  Initial values for the history of
# incidence proxy values are taken from the end of one of the five seasons
# in the preceding simulation (and set somewhat arbitrarily for the first
# simulation).  Some randomness is introduced into the initial values,
# and in the daily values generated.
#
# Simulations are done day-by-day (not week-by-week).  So that the
# simulated values are comparable to the weekly incidence proxies,
# they are multiplied by 7 when plotted.

set.seed(1)

nsims_warmup <- 10   # Number of initial simulations to "warm up"
nsims_plotted <- 15  # Number of subsequent simulations that are plotted
nsims <- nsims_warmup + nsims_plotted  # Total number of simulations

mc <- coef(model)

proxy1 <- R_est[,paste0(virus_group[1],"_proxy")]
proxy2 <- R_est[,paste0(virus_group[2],"_proxy")]

daily_decay <- imm_decay ^ (1/7)
ltdaily_decay <- ltimm_decay ^ (1/7)

# Generation interval distribution used by Kissler, et al.

SARS_shape <- 2.35
SARS_scale <- 9.48

SARS_gen_interval <-
  dweibull (1:ceiling(qweibull(0.99, shape=SARS_shape, scale=SARS_scale)),
            shape=SARS_shape, scale=SARS_scale)

gen_interval <- SARS_gen_interval # / sum(SARS_gen_interval)
rev_gen_interval <- rev(gen_interval)

# Pre-compute trend and seasonal effects for all days.

tn <- ncol(trend_spline)
tseff <-
( if (seffect_type=="e2")
    seffect_e2(yrsd)
  else
    seffect_e3(yrsd) + as.vector(predict (trend_spline,yrsd) %*% mc[1:tn])
)

# Initial levels for exponentially-decaying past incidence.

q <- c (quantile(proxy1,0.1), quantile(proxy2,0.1))
t <- q / (1-imm_decay[virus_group])
tlt <- q / (1-ltimm_decay[virus_group])
past <- list (rep(q[1],length(gen_interval)), rep(q[2],length(gen_interval)))

# Do the simulations.

p <- numeric(2)

ylim <-  max (proxy1, proxy2)

sim <- rep (list(numeric(7*length(start))), 2)
sims <- list()

sv_past <- NULL

wsave <- sample(rep(1:5,length=nsims))

for (w in 1:nsims)
{
  # Do one simulation, for all years.

  Rt_offset_sd <- 0.05                     # AR(1) process that modifies
  Rt_offset_alpha <- 0.9                   #   modelled Rt values
  Rt_offset <- rnorm(1,0,Rt_offset_sd)

  Rt_iid_noise_sd <- 0.05                  # Extra iid noise in Rt

  for (i in 1:(7*length(start)))
  {
    for (j in 1:2)
    { virus <- virus_group[j]
      log_Rt <- tseff[i] + mc[paste0(virus,"_same")] * t[j]
      if (immune_type!="i4")
      { log_Rt <- log_Rt + mc[paste0(virus,"_other")] * t [if (j==1) 2 else 1]
      }
      if (immune_type=="i3" || immune_type=="i4")
      { log_Rt <- log_Rt + mc[paste0(virus,"_samelt")] * tlt[j] +
                           mc[paste0(virus,"_otherlt")] * tlt[if(j==1) 2 else 1]
      }
      log_Rt <- log_Rt + mc[paste0(virus,"_overall")]
      inf <- sum (past[[j]]*rev_gen_interval)
      Rt_offset <- Rt_offset_alpha*Rt_offset +
                   sqrt(1-Rt_offset_alpha^2) * rnorm(1,0,Rt_offset_sd)
      p[j] <- exp (log_Rt + Rt_offset + rnorm(1,0,Rt_iid_noise_sd)) * inf
      past[[j]] <- c (past[[j]][-1], p[j])
      sim[[j]][i] <- p[j]

      t[j] <- p[j] + t[j]*daily_decay[virus]
      tlt[j] <- p[j] + tlt[j]*ltdaily_decay[virus]
    }

    if ((i+5) %% 365 == 0)
    { wsave[w] <- wsave[w] - 1
      if (wsave[w] == 0)
      { sv_past <- past
        sv_t <- t
        sv_tlt <- tlt
      }
    }
  }

  # Save simulation for later plotting, if past warm-up.

  if (w > nsims_warmup)
  { sims <- c(sims,list(sim))
    ylim <- max (ylim, 7*sim[[1]], 7*sim[[2]], na.rm=TRUE)
  }

  # Set up initial state for next simulation. Randomized a bit.

  n <- exp(rnorm(2,0,0.2))
  for (j in 1:2)
  { past[[j]] <- sv_past[[j]] * n[j]
  }
  t <- sv_t * n
  tlt <- sv_tlt * n
}

# Plot the observed incidence, then saved simulations.

par(mfrow=c(4,1))

plot (start, rep(0,length(start)),
      ylim=c(0,1.02*ylim), yaxs="i", type="n", ylab="Incidence proxy")

lines (start, proxy1, col="blue")
lines (start, proxy2, col="red")

title (paste
 ("Observed proxies for",virus_group[1],"(blue) and",virus_group[2],"(red)"))

for (k in seq_along(sims))
{
  sim <- sims[[k]]

  plot (start, rep(0,length(start)),
        ylim=c(0,1.02*ylim), yaxs="i", type="n", ylab="Simulated incidence")

  weekly <- filter (sim[[1]],rep(1,7)) [seq(4,length(sim[[1]]),by=7)]
  lines (start, weekly, col="blue")
  weekly <- filter (sim[[2]],rep(1,7)) [seq(4,length(sim[[2]]),by=7)]
  lines (start, weekly, col="red")

  if (k==1) title (paste ("Simulations of",virus_group[1],"and",virus_group[2]))
}

# ----- END OF LOOP OVER THE TWO VIRUS GROUPS -----

}

# ALL DONE.

dev.off()
