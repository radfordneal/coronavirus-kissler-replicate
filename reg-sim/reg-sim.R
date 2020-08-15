# Code to simulate viral incidence according to a model fit to estimated Rt
# values, read from the R-model directory.
#
# Options are specified by arguments in the R command after --args:
#
#   - The R estimates to use, from corresponding file (required)
#   - Model of "immunity" (required) - i2 (exp decay), i3 (short & long),
#     i4 (like i3 but with no short-term cross-immunity)
#   - Model of seasonal effect (required) - e2 (sine), e3 (Fourier)
#   - Whether heteroskedasticity w.r.t. virus is modelled - het for "yes" 
#     (default "no")
#   - Transformation to apply before comparing observed and simulated
#     incidence (required) - identity, sqrt, log
#   - Number of iterations of optimization to do, in the form: opt:<n>
#     Default is to not do optimization
#   - Whether to initialize from regression model (in R-model directory), 
#     which is the default, or a previous optimization run (in this directory),
#     with a given suffix, in the form: init:<suffix>
#   - Suffix for saving plots and parameter estimates, in form: save:<suffix>
#     Default is no suffix
#
# Produces various plots that are written to the file with name
# reg-sim-<R-estimates>-<in>-<en>-<trans>[-<suffix>].pdf.  
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
season_type <- "s2"
seffect_type <- getarg (c("e2","e3"))
itrans_arg <- getarg (c("identity","sqrt","log"))

itrans <- get(itrans_arg)

het_virus <- sum(args=="het") == 1
if (het_virus) args <- args [! (args %in% "het")]

opt_iters <- 0

if (any (substr(args,1,4) == "opt:"))
{ stopifnot (sum (substr(args,1,4) == "opt:") == 1)
  opt_arg <- substr (args [substr(args,1,4) == "opt:"], 5, 100)
  opt_iters <- eval(parse(text=opt_arg))
  args <- args [substr(args,1,4) != "opt:"]
}

init_suffix <- NULL

if (any (substr(args,1,5) == "init:"))
{ stopifnot (sum (substr(args,1,5) == "init:") == 1)
  init_suffix <- substr (args [substr(args,1,5) == "init:"], 6, 100)
  args <- args [substr(args,1,5) != "init:"]
}

save_suffix <- NULL

if (any (substr(args,1,5) == "save:"))
{ stopifnot (sum (substr(args,1,5) == "save:") == 1)
  save_suffix <- substr (args [substr(args,1,5) == "save:"], 6, 100)
  args <- args [substr(args,1,5) != "save:"]
}

R_estimates <- args

stopifnot(length(R_estimates)==1)

file_base <- paste0 (R_estimates,"-Rt-s2-",immune_type,"-",seffect_type,
                     if (het_virus) "-het")
file_base_sim <- paste0("reg-sim-",gsub("Rt-s2-","",file_base),"-",itrans_arg)

if (FALSE)  # Small settings for testing
{ nsims <- 1000         # Number of simulations in full set
  sub <- 30             # Number of simulations in subset
  full_interval <- 10   # Interval for doing full set of simulations
} else      # Settings for serious run
{ nsims <- 100000       # Number of simulations in full set
  sub <- 2500           # Number of simulations in subset
  full_interval <- 30   # Interval for doing full set of simulations
}

n_plotted <- 32         # Number of simulations to plot

Min_inf <- 0.0015       # Minimum infectivity


# PLOT SETUP.

pdf (paste0 (file_base_sim, if (!is.null(save_suffix)) "-", save_suffix,
             ".pdf"), height=8, width=6)
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

if (is.null(init_suffix))
{ P_init <- readRDS (paste0("../R-model/R-model-",file_base,"-",
                            names(virus_groups)[g],".model"))
}
else if (init_suffix=="")
{ P_init <- readRDS (paste0(file_base_sim,"-",names(virus_groups)[g],".model"))
}
else
{ P_init <- readRDS (paste0(file_base_sim,"-",names(virus_groups)[g],
                            "-",init_suffix,".model"))
}

if (is.null(P_init$imm_initial))
{ P_init$imm_initial <- P_init$imm_decay      # to get names
  P_init$imm_initial[] <- 2.0
}
if (is.null(P_init$ltimm_initial))
{ P_init$ltimm_initial <- P_init$ltimm_decay  # to get names
  P_init$ltimm_initial[] <- 5.0
}

if (is.null(P_init$Rt_offset))
{ P_init$Rt_offset <- c (alpha=0.9, sd=0.05)
}

# Transform some parameters to avoid constraints.

P_init$imm_decay <- log (P_init$imm_decay / (1 - P_init$imm_decay))
P_init$ltimm_decay <- log (P_init$ltimm_decay / (1 - P_init$ltimm_decay))
P_init$Rt_offset["alpha"] <- 
  log (P_init$Rt_offset["alpha"] / (1-P_init$Rt_offset["alpha"]))
P_init$Rt_offset["sd"] <- log (abs (P_init$Rt_offset["sd"]))

print_model_parameters <- function (P)
{ if (seffect_type=="e3")
  { cat("Trend:\n")
    print(t(t(P$mc_trend)))
    cat("\n")
  }
  cat("Seasonality:\n")
  print(t(t(P$mc_seasonality)))
  cat("\n")
  cat("Viral coefficients:\n")
  print(t(t(P$mc_viral)))
  cat("\n")
  cat("Immune decay:\n")
  print(t(t(P$imm_decay)))
  cat("\n")
  cat("Long-term immune decay:\n")
  print(t(t(P$ltimm_decay)))
  cat("\n")
  cat("Log initial short-term immunity:\n")
  print(t(t(P$imm_initial)))
  cat("\n")
  cat("Log initial long-term immunity:\n")
  print(t(t(P$ltimm_initial)))
  cat("\n")
  cat("Offset model:\n")
  print (t(t(P$Rt_offset)))
  cat("\n")
}

cat("\n\nMODEL FOR GROUP",names(virus_groups)[g],":",virus_group,"\n\n")

print_model_parameters (P_init)


# SET GENERATION INTERVAL DISTRIBUTION USED BY KISSLER, ET AL.

SARS_shape <- 2.35
SARS_scale <- 9.48

SARS_gen_interval <-
  dweibull (1:ceiling(qweibull(0.99, shape=SARS_shape, scale=SARS_scale)),
            shape=SARS_shape, scale=SARS_scale)

gen_interval <- SARS_gen_interval # / sum(SARS_gen_interval)
rev_gen_interval2 <- rep(rev(gen_interval),2)


# PROXIES FOR VIRUS INCIDENCE.

proxy <- list (R_est[,paste0(virus_group[1],"_proxy")],
               R_est[,paste0(virus_group[2],"_proxy")])

tproxy <- list (itrans(proxy[[1]]), itrans(proxy[[2]]))


# FUNCTION TO RUN SIMULATIONS.  Returns 'nsims' simulation results,
# from 'nsims' five-year simulations.  If 'subset' is non-null, the
# random number generation for the 'nsims' simulations is done as if
# 'full' simulations were being done, with 'nsims' simulations being
# taken from the subset of these defined by the indexes in 'subset'.
#
# 'P' is a list of parameter values affecting the simulation.
#
# Simulations are done day-by-day (not week-by-week), using the trend spline
# and seasonal effects, and the long and short term immunity based on previous
# incidence. However, the results are weekly values, summed over the days in
# each week.
#
# Initial values for the exponential averages of past incidence are 
# generated randomly from vivariate normal distributions found from
# the observed proxies (in conjunction with current parameter values).
# The initial history of past incidence (used with the generation interval 
# distribution) is constant, with value determined from the short-term
# exponential average.
#
# Randomness is introduced into the Rt values, differently for each 
# according to the 'Rt_offset' parameters.
#
# The 'cache' argument may be set to an environment, in which the random
# numbers used are saved, or from which they are taken if they are already
# there.
#
# If the 'info' argument is TRUE, processor time and gc information is 
# printed.
#
# The returned value is a list of, for each virus, a matrix with
# dimensions nsims x number of weeks.

run_sims <- function (nsims, full=nsims, subset=NULL, 
                      cache=NULL, info=FALSE, P=P_init)
{
  stopifnot (nsims == if (is.null(subset)) full else length(subset))

  if (info)
  { start_time <- proc.time()
    cat("GC before simulations:\n")
    print(gc())
  }

  set.seed(seed)
  
  if (is.null(cache))          # just generate random numbers, no cache
  { randn <- function () 
    { if (is.null(subset)) rnorm(full) else rnorm(full)[subset]
    }
  }
  else if (is.null(cache$rn))  # generate random numbers, save them in cache
  { cache$rn <- list()
    next_rn <- 1
    randn <- function () 
    { r <- if (is.null(subset)) rnorm(full) else rnorm(full)[subset]
      cache$rn[[next_rn]] <- r
      next_rn <<- next_rn + 1
      r
    }
  }
  else                         # take random numbers from cache
  { next_rn <- 1
    randn <- function () 
    { r <- cache$rn[[next_rn]]
      next_rn <<- next_rn + 1
      r
    }
  }

  # Pre-compute trend and seasonal effects for all days.
  
  tseff <-
  ( if (seffect_type=="e2")
      seffect_e2(yrsd,P$mc_seasonality)
    else
      seffect_e3(yrsd,P$mc_seasonality,0) +
        as.vector (predict(trend_spline,yrsd) %*% P$mc_trend)
  )

  mc <- P$mc_viral

  alph <- 1 / (1 + exp(-P$Rt_offset["alpha"]))
  sd <- exp(P$Rt_offset["sd"])

  daily_decay <- (1/(1+exp(-P$imm_decay))) ^ (1/7)
  ltdaily_decay <- (1/(1+exp(-P$ltimm_decay))) ^ (1/7)

  # Initial levels for exponentially-decaying averages and past incidence.
  #
  #   t      list of short-term exponetial sums for each virus
  #   tlt    list of long-term exponetial sums for each virus
  #   past   list of matrices of assumed past incidence values, one matrix
  #          for each virus, dimension nsims x length(gen_interval)

  t <- list (exp (P$imm_initial[1] + 0.3*randn()),
             exp (P$imm_initial[2] + 0.3*randn()))
  tlt <- list (exp (P$ltimm_initial[1] + 0.3*randn()),
               exp (P$ltimm_initial[2] + 0.3*randn()))

  sv_t <<- list(t); sv_tlt <<- list(tlt)  # for later plots

  past <- 
    list (matrix (t[[1]]*(1-daily_decay[1]), nsims, length(gen_interval)),
          matrix (t[[2]]*(1-ltdaily_decay[2]), nsims, length(gen_interval)))

  past_next <- rep(1,2)

  # cat("initial past history:\n")
  # print(past)

  # Space to store simulation results. A list of two matrices, one for each
  # virus, of dimension total nsims x number of weeks.

  wsims <- rep (list (matrix (0, nsims, length(start))), times=2)

  Rt_offset <- exp(P$Rt_offset["sd"]) * randn()  # initialize AR(1) process that
                                                 #   modifies modelled Rt values

  # Simulate for all days, for all simulations, adding to weekly results.

  for (day in 1:(7*length(start)))
  {
    Rt_offset <- alph * Rt_offset + sd * sqrt(1-alph^2) * randn()

    for (vi in 1:2)
    { 
      virus <- virus_group[vi]

      # Compute vector of nsims log R values based on trend, seasonality, 
      # and immunity due to past infections.

      log_Rt <- tseff[day] +
                mc [paste0(virus,"_overall")] +
                mc [paste0(virus,"_same")] * t [[vi]]
      if (immune_type!="i4")
      { log_Rt <- log_Rt + 
                 mc [paste0(virus,"_other")] * t [[if (vi==1) 2 else 1]]
      }
      if (immune_type=="i3" || immune_type=="i4")
      { log_Rt <- log_Rt + 
                  mc [paste0(virus,"_samelt")] * tlt [[vi]] +
                  mc [paste0(virus,"_otherlt")] * tlt [[if (vi==1) 2 else 1]]
      }

      rs <- length(gen_interval) + 2 - past_next[vi]
      inf <- as.vector (past[[vi]] %*% 
                        rev_gen_interval2 [rs : (rs+length(gen_interval)-1)])

      p <- pmax(Min_inf,inf) * exp (log_Rt + Rt_offset)
      if (any(is.na(p))) stop("NA in prevalence")

      wk <- ceiling(day/7)
      wsims[[vi]][,wk] <- wsims[[vi]][,wk] + p

      past[[vi]][,past_next[vi]] <- p
      past_next[vi] <- past_next[vi] %% length(gen_interval) + 1

      t[[vi]] <- p + t[[vi]] * daily_decay[virus]
      tlt[[vi]] <- p + tlt[[vi]] * ltdaily_decay[virus]
    }

    if (wk %% 52 == 0 && day %% 7 == 1)
    { sv_t <<- c(sv_t,list(t))
      sv_tlt <<- c(sv_tlt,list(tlt))
    }
  }

  if (info)
  { cat("GC after simulations:\n")
    print(gc())
    cat("Time for simulations:\n")
    print (proc.time()-start_time)
  }

  # Return weekly values.

  wsims
}  


# TRANSFORM SIMULATED VALUES.

itrans_wsims <- function (wsims)
{ twsims <- vector("list",2)
  for (vi in 1:2) twsims[[vi]] <- itrans(wsims[[vi]])
  twsims
}


# FIND ERRORS FOR EACH SIMULATION BASED ON AR(1) MODEL.  The error for
# the first data point for each virus is neglected (effectively
# treating it as if it came from a broad, fixed distribution). Errors
# in other data points are modelled as an AR(1) process, with AR
# parameters err_alpha (one for each virus). The err_sd parameters
# (one for each virus) give the standard deviations of the errors,
# with err_sd^2 * (1-err_alpha^2) being the variance of the
# innovations in the AR(1) process.

sim_errors <- function (twsims, err_alpha, err_sd)
{
  e <- 0
  for (vi in 1:2)
  { tp <- tproxy[[vi]]
    tv <- twsims[[vi]]
    ev <- 0
    res0 <- tp[1] - tv[,1]
    for (j in 2:ncol(tv))
    { res1 <- tp[j] - tv[,j]
      ev <- ev + (res1 - err_alpha[vi]*res0)^2
      res0 <- res1
    }
    ivar <- err_sd[vi]^2 * (1-err_alpha[vi]^2)
    e <- e + ev / ivar
  }
  e
}


# FIND POSTERIOR PROBABILITIES OF ALTERNATIVE INFECTION HISTORIES.

pprob <- function (errors)
{
  p <- exp (-0.5*(errors-min(errors)))
  p / sum(p)
}


# ESTIMATE ERROR ALPHAS AND STANDARD DEVIATIONS.  Assumes that there
# are nsims histories in total.

est_error_model <- function (twsims, init_err_alpha=0.9, init_err_sd=1.0,
                             verbose=FALSE)
{
  if (verbose) cat("\nEstimation of error model\n\n")

  err_alpha <- rep (init_err_alpha, length=2)
  err_sd <- rep (init_err_sd, length=2)

  if (TRUE)  # can set to TRUE to disable estimation of error model
  { if (verbose) cat("Using initial values\n")
    return (list (alpha=err_alpha, sd=err_sd))
  }

  for (i in 0:6)
  { if (verbose) cat ("iter",i,"\n")
    if (i>0)
    { pp <- pprob (sim_errors (twsims, err_alpha, err_sd))
      for (vi in 1:2)
      { tp <- tproxy[[vi]]
        tv <- twsims[[vi]]
        var <- 0
        cov <- 0
        res0 <- tp[1] - tv[,1]
        for (j in 2:ncol(tv))
        { var <- var + res0^2
          res1 <- tp[j] - tv[,j]
          cov <- cov + res1*res0
          res0 <- res1
        }
        var <- sum(pp*var) / (ncol(tv)-1)
        cov <- sum(pp*cov) / (ncol(tv)-1)
        # cat("var",var,"\n")
        # cat("cov",cov,"\n")
        err_alpha[vi] <- cov/var
        esd <- 0
        for (j in 2:ncol(tv))
        { esd <- esd + ((tp[j]-tv[,j]) - err_alpha[vi]*(tp[j-1]-tv[,j-1]))^2
        }
        err_sd[vi] <- sqrt (sum(pp*esd) / (1-err_alpha[vi]^2) / (ncol(tv)-1))
      }
    }
    if (verbose)
    { ll <- log_lik (twsims, err_alpha, err_sd, full=nsims)
      cat ("  err_alpha",round(err_alpha,6),
           ": err_sd",round(err_sd,3),
           ": log likelihood",round(ll,5),
           "\n")
      if (is.na(ll)) stop("log likelihood is NA")
    }
  }

  list (alpha=err_alpha, sd=err_sd)
}


# COMPUTE THE LOG LIKELIHOOD BASED ON A SET OF SIMULATION RESULTS.
# The multiple simulations are taken as a Monte Carlo estimate of the
# distribution of the infection histories for given parameters, with
# the simulations assumed to all those with non-negligible posterior
# probability out of 'full' total.  The likelihood is the average of the
# probability of the observed infection proxies given these histories, 
# using the AR(1) error model. Constant terms involving pi are omitted.

log_lik <- function (twsims, err_alpha, err_sd, errors, full=length(errors))
{
  if (missing(errors))
  { errors <- sim_errors (twsims, err_alpha, err_sd)
  }

  n <- ncol(twsims[[1]])
  mine <- min(errors)

  ( log (sum(exp(-0.5*(errors-mine))) / full) - 0.5*mine 
      - (n-1) * sum(log(err_sd*sqrt(1-err_alpha^2))) )
}


# COMPUTE PROFILE LOG LIKELIHOOD, ESTIMATING ERROR MODEL PARAMETERS.

profile_log_lik <- function (twsims, ...)
{
  est <- est_error_model (twsims, verbose=FALSE)
  log_lik (twsims, est$alpha, est$sd, ...)
}


# FREE MEMORY.

wsims <- twsims <- NULL 
wsims_subset <- twsims_subset <- NULL
wsims_new <- twsims_new <- NULL


# DO SIMULATIONS WITH ORIGINAL PARAMETER ESTIMATES (FROM FILE).

RNGversion("2.15.1")
seed <- 1

# Rprofmemt (nelem=2*nsims+1)

cat ("\nSIMULATIONS WITH ORIGINAL PARAMETER ESTIMATES\n\n")

start_time <- proc.time()

wsims <- run_sims (nsims, info=TRUE)

twsims <- itrans_wsims (wsims)

em <- est_error_model(twsims,verbose=TRUE)
err_alpha <- em$alpha
err_sd <-em$sd

errors <- sim_errors(twsims,err_alpha,err_sd)
pp <- pprob(errors)
# print(errors)
# print(pp)

cat ("\nHighest posterior probabilities:\n")
print (round(sort(pp,decreasing=TRUE)[1:16],6))

ll <- log_lik(twsims,err_alpha,err_sd,errors=errors)
cat ("\nLog likelihood,", length(pp), "simulations:", round(ll,5), "\n")

wmx <- which.max(pp)

cat ("\nRunning subset simulation\n\n")

# print(sort(pp,decreasing=TRUE)[1:sub])
# print(order(pp,decreasing=TRUE)[1:sub])

high <- unique ((order(pp,decreasing=TRUE)[1:sub]-1) %% nsims + 1)
subn <- length(high)  # currently alwas equal to sub, but wasn't before...

# print(high)

cache <- new.env()
wsims_subset <- run_sims (subn, full=nsims, subset=high, cache=cache, 
                          info=TRUE)

if (TRUE)  # enable to check that caching works
{
  cat ("\nRunning subset simulation a second time\n\n")
  
  wsims_subset2 <- 
    run_sims (subn, full=nsims, subset=high, cache=cache, info=TRUE)
  
  stopifnot(identical(wsims_subset,wsims_subset2))
  wsims_subset2 <- NULL
}
  
twsims_subset <- itrans_wsims (wsims_subset)

# print(high)
# print(pp[high])

cat("\n")

errors_subset <- sim_errors(twsims_subset,err_alpha,err_sd)
# print(errors_subset)
# print(pprob(errors_subset))

cat("Lowest probability in top",subn,"is",pp[high[subn]],"\n")
cat("Total probability in is",round(sum(pp[high]),5),"\n")
cat ("Log likelihood based on subset of",subn,"simulations:", 
      round(log_lik(twsims_subset,err_alpha,err_sd,full=nsims),5),
      "\n\n")


# ESTIMATE MODEL PARAMETERS.

if (opt_iters > 0) {

cat("\nESTIMATING MODEL PARAMETERS\n\n")

start_time_est <- proc.time()

if (FALSE)
{ source("estimate-nlm.R")
}
else if (FALSE)
{ source("estimate-nlm-autodiff.R")
}
else if (TRUE)
{ source("estimate-grad-autodiff.R")
}

cat("Time for estimation:\n")
print(proc.time()-start_time_est)

cat("\nInitial and new parameter values:\n\n")
print_model_parameters (mapply (cbind, 
                                Initial=P_init, New=P_new, Change=P_new-P_init))
}

# DO SIMULATIONS WITH NEW PARAMETERS.

if (opt_iters > 0) {

cat ("SIMULATIONS WITH NEW PARAMETER ESTIMATES\n\n")

wsims_new <- run_sims (nsims, P=P_new, info=TRUE)

twsims_new <- itrans_wsims (wsims_new)

em_new <- est_error_model(twsims_new,verbose=TRUE)

errors_new <- sim_errors(twsims_new,em_new$alpha,em_new$sd)
pp_new <- pprob(errors_new)
# print(errors_new)
# print(pp_new)

wmx_new <- which.max(pp_new)

cat ("\nHighest posterior probabilities:\n")
print (round(sort(pp_new,decreasing=TRUE)[1:16],6))

ll_new <- log_lik (twsims_new, em_new$alpha, em_new$sd)

cat ("\nLog likelihood for new parameters,", length(pp_new), "simulations",
      round(ll_new,5), "\n")

cat("GC after post-simulation work:\n")
print(gc())
cat("Total processing time for this group of viruses:\n")
print(proc.time()-start_time)

}

# PLOTS FOR THIS VIRUS GROUP.  Corresponding plots use the same scales.

# PLOT BEST FIT AND OTHER SIMULATIONS.

plot_sims <- function (wsims, twsims, wmx, pp, ll)
{
  for (s in 0:n_plotted)
  {

    plot (start, rep(0,length(start)),
          ylim=c(0,1.02*yupper), yaxs="i", type="n", 
          ylab="Simulated incidence")
  
    ss <- if (s==0) wmx else s
    lines (start, wsims[[1]][ss,], col="blue")
    lines (start, wsims[[2]][ss,], col="red")
  
    if (s==0)  # the best-fit simulation
    { 
      title (paste ("Best fit simulation out of",nsims,
                    "for",virus_group[1],"and",virus_group[2]))
  
      plot (start, rep(0,length(start)),
        ylim=itrans(c(0.98*ylower,1.02*yupper)), yaxs="i", type="n", 
        ylab=paste(if (itrans_arg!="identity") itrans_arg, 
                   "simulated incidence"))
    
      lines (start, twsims[[1]][wmx,], col="blue")
      lines (start, twsims[[2]][wmx,], col="red")
  
      par(mfrow=c(4,2))
  
      plot (sort(pp,decreasing=TRUE)[1:20],pch=20,xlab="",
            ylab="posterior prob")
      title ("Posterior prob. of most-likely histories")
      plot (log10(sort(pp,decreasing=TRUE)[1:20]),pch=20,xlab="",ylim=c(-30,0),
            ylab="log10 of posterior prob")
      title (paste ("Log likelihood based on sum:", round(ll,3)))

      par(mfrow=c(4,1))
    }
    if (s==1) 
    { title (paste ("Other simulations of",virus_group[1],"and",virus_group[2]))
    }
  }
}

yupper <-  max (proxy[[1]], proxy[[2]], 
                wsims[[1]][c(wmx,1:n_plotted),], 
                wsims[[2]][c(wmx,1:n_plotted),])
if (opt_iters > 0) yupper <- max (yupper,
                                  wsims_new[[1]][wmx_new,],
                                  wsims_new[[2]][wmx_new,])
ylower <-  min (proxy[[1]], proxy[[2]], exp(-4),
                wsims[[1]][wmx,],
                wsims[[2]][wmx,])
if (opt_iters > 0) ylower <- min (ylower,
                                  wsims_new[[1]][wmx_new,],
                                  wsims_new[[2]][wmx_new,])

# Observed proxies.

par(mfrow=c(4,1))

plot (start, rep(0,length(start)),
      ylim=c(0,1.02*yupper), yaxs="i", type="n",
      ylab="Incidence proxy")

lines (start, proxy[[1]], col="blue")
lines (start, proxy[[2]], col="red")

title (paste
 ("Observed proxies for",virus_group[1],"(blue) and",virus_group[2],"(red)"))

plot (start, rep(0,length(start)),
      ylim=itrans(c(0.98*ylower,1.02*yupper)), yaxs="i", type="n",
      ylab=paste(if (itrans_arg!="identity") itrans_arg, "incidence proxy"))

lines (start, tproxy[[1]], col="blue")
lines (start, tproxy[[2]], col="red")

# Best fit and other simulations.

plot_sims (wsims, twsims, wmx, pp, ll)

if (opt_iters > 0) {

# Observed proxies, again.

par(mfrow=c(4,1))

par(mfrow=c(4,1))

plot (start, rep(0,length(start)),
      ylim=c(0,1.02*yupper), yaxs="i", type="n",
      ylab="Incidence proxy")

lines (start, proxy[[1]], col="blue")
lines (start, proxy[[2]], col="red")

title (paste
 ("Observed proxies for",virus_group[1],"(blue) and",virus_group[2],"(red)"))

plot (start, rep(0,length(start)),
      ylim=itrans(c(0.98*ylower,1.02*yupper)), yaxs="i", type="n",
      ylab=paste(if (itrans_arg!="identity") itrans_arg, "incidence proxy"))

lines (start, tproxy[[1]], col="blue")
lines (start, tproxy[[2]], col="red")

# Best fit and other simulations with new parameters.

plot_sims (wsims_new, twsims_new, wmx_new, pp_new, ll_new)

# Plot exponential averages at starts of seasons for first 500 simulations.
# Best-fit simulation is larger, in red.  Uses parameters from end of
# optimization.

par(mfrow=c(4,2))

first <- c(1:min(500,nsims),wmx_new)

for (i in 1:6)
{ for (vi in 1:2)
  { plot (pmax(-3,log(sv_t[[i]][[vi]][first])), 
          pmax(3,log(sv_tlt[[i]][[vi]][first])), 
          xlab="log short-term average", ylab="log long-term average", 
          xlim=c(-3,6), ylim=c(3,7),
          pch = ifelse (first==wmx_new, "O", "."), 
          col = 1+(first==wmx_new))
    title (paste (virus_group[vi],"year",i-1))
  }
}

# Plot posterior probabilities of runs before and after parameter change.

par(mfrow=c(1,1))

w <- 1:min(10000,length(pp))
plot (log(pp[w]), log(pp_new[w]), pch=".", asp=1,
      ylab="log new probability", xlab="log probability")
abline(0,1,col="green")
abline(ll-ll_new,1,col="green")
title (paste ("Change in log prob. of simulation runs, log lik. change",
               round(ll_new-ll,2),"  "))

# Plot components of original and new models.

source("../R-model/plot-components.R")

plot_context <- readRDS (paste0 ("../R-model/R-model-",file_base,"-",
                                 names(virus_groups)[g],".context"))

model_x <- plot_context$model_x
model_df <- plot_context$model_df

par(mfcol=c(5,4))
sv <- par (cex.main=3/4, cex.lab=2/3, cex.axis=4/10, mgp=c(0.8,0.18,0))
  
for (virus in virus_group)
{ for (s in unique(model_df$season))
  { plot_components (c(P_init$mc_trend,P_init$mc_seasonality,P_init$mc_viral),
                     model_x, model_df, s, virus, logarithmic=TRUE,
                     title = "Initial")
  }
  for (s in unique(model_df$season))
  { plot_components (c(P_new$mc_trend,P_new$mc_seasonality,P_new$mc_viral),
                     model_x, model_df, s, virus, logarithmic=TRUE,
                     title = "New")
  }
}

par(sv)

}

# Save the new parameter values.

if (opt_iters > 0) {

# Undo transformations done for some parameters.

P_out <- P_new
P_out$imm_decay <- 1 / (1 + exp(-P_out$imm_decay))
P_out$ltimm_decay <- 1 / (1 + exp(-P_out$ltimm_decay))
P_out$Rt_offset["alpha"] <- 1 / (1 + exp(-P_out$Rt_offset["alpha"]))
P_out$Rt_offset["sd"] <- exp(P_out$Rt_offset["sd"])

# Save parameters to file.

if (is.null(save_suffix) || save_suffix=="")
{ saveRDS (P_out, 
           file = paste0(file_base_sim,"-",names(virus_groups)[g],".model"),
           version=2)
}
else
{ saveRDS (P_out,
           file = paste0(file_base_sim,"-",names(virus_groups)[g],
                         "-",save_suffix,".model"), 
           version=2)
}

}


# ----- END OF LOOP OVER THE TWO VIRUS GROUPS -----

}

# ALL DONE.

dev.off()

print(proc.time())
