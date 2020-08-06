# Code to simulate viral incidence according to a model fit to estimated Rt
# values, read from the R-model directory.
#
# Options are specified by arguments in the R command after --args:
#
#   - The R estimates to use, from corresponding file (required)
#   - Model of "immunity" (required) - i2 (exp decay), i3 (short & long),
#                           i4 (like i3 but with no short-term cross-immunity)
#   - Model of seasonal effect (required) - e2 (sine), e3 (Fourier)
#   - Whether heteroskedasticity w.r.t. virus is modelled - het for "yes" 
#     (default "no")
#   - Transformation to apply before comparing observed and simulated
#     incidence (required) - identify, sqrt, log
#
# Produces various plots that are written to the file with name
# reg-sim-<R-estimates>-<in>-<en>-<trans>.pdf.  
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
itrans_arg <- getarg (c("identity","sqrt","log"))

itrans <- get(itrans_arg)

het_virus <- sum(args=="het") == 1
if (het_virus) args <- args [! (args %in% "het")]

R_estimates <- args

stopifnot(length(R_estimates)==1)

file_base <- paste0 (R_estimates,"-Rt-s2-",immune_type,"-",seffect_type,
                     if (het_virus) "-het")

nsims <- 10000
warmup <- 6
keep <- 30
sub <- 10
n_plotted <- 30


# PLOT SETUP.

pdf (paste0 ("reg-sim-",gsub("Rt-s2-","",file_base),"-",itrans_arg,".pdf"),
     height=8, width=6)
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

print_model_parameters <- function (P)
{ cat("Coefficients:\n")
  print(t(t(P$mc)))
  cat("\n")
  cat("Immune decay:\n")
  print(P$imm_decay)
  cat("\n")
  cat("Long-term immune decay:\n")
  print(P$ltimm_decay)
  cat("\n")
  cat("Offset model:\n")
  print (c(Rt_offset_alpha=P$Rt_offset_alpha, Rt_offset_sd=P$Rt_offset_sd))
  cat("\n")
}

model <- r$model
imm_decay <- r$imm_decay[virus_group]
ltimm_decay <- r$ltimm_decay[virus_group]

P_init <- list (mc = coef(model), imm_decay = imm_decay, 
                ltimm_decay = ltimm_decay, Rt_offset_alpha = 0.9,
                Rt_offset_sd = 0.05)

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


# FUNCTION TO RUN SIMULATIONS.  Returns nsims*keep simulation results,
# found by 'nsims' series of five-year simulations, with 'warmup' simulations
# done before 'keep' simulations that are used.  If 'subset' is non-null,
# the random number generation for the 'nsims' simulations is done as if 
# 'full' simulations were being done, with 'nsims' simulations being taken 
# from the subset of these defined by the indexes in 'subset'.
#
# Simulations are done day-by-day (not week-by-week), using the trend spline
# and seasonal effects, and the long and short term immunity based on previous
# incidence. However, the results are weekly values, summed over the days in
# each week.
#
# Initial values for the history of incidence are taken from the end
# of one of the seasons in the preceding five-year simulation (and set
# somewhat arbitrarily for the first five-year simulation), the same
# one for all simulations done. Some randomness is introduced into the
# Rt values, differently for each simulation.
#
# 'P' is a list of parameter values affecting the simulation.
#
# The 'cache' argument may be set to an environment, in which the random
# numbers used are save, or from which they are taken if they are already
# there.
#
# If the 'info' argument is TRUE, processor time and gc information is 
# printed.
#
# The returned value is a list of, for each virus, a matrix with
# dimensions nsims*keep x number of weeks.

run_sims <- function (nsims, warmup, keep, full=nsims, subset=NULL, 
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
  #   t      list of short-term exponetial sums for each virus
  #   tlt    list of long-term exponetial sums for each virus
  #   past   list of matrices of assumed past incidence values, one matrix
  #          for each virus, dimension nsims x length(gen_interval)

  q <- c (quantile(proxy[[1]],0.1), quantile(proxy[[2]],0.1))

  t <- matrix (q / (1-imm_decay[virus_group]), nsims, 2, byrow=TRUE)
  tlt <- matrix (q / (1-ltimm_decay[virus_group]), nsims, 2, byrow=TRUE)

  t <- list (t[,1], t[,2])        # gives faster access
  tlt <- list (tlt[,1], tlt[,2])

  past <- list (matrix (q[1], nsims, length(gen_interval)),
                matrix (q[2], nsims, length(gen_interval)))
  past_next <- rep(1,2)

  # Space to store simulation results. A list of two matrices, one for each
  # virus, of dimension total number of weeks x nsims*keep.

  wsims <- rep (list (matrix (0, nsims*keep, length(start))), times=2)

  # Stuff for saving a state for use to initialize the next simulation.
  
  wsave <- sample(rep(1:5,length=warmup+keep))
  sv_past <- NULL
  
  for (w in 1:(warmup+keep))
  {
    # cat("w =",w,"\n")

    # Do next simulations of a five-year period.

    Rt_offset <- P$Rt_offset_sd * randn()  # initialize AR(1) process that
                                           #   modifies modelled Rt values
                                            
    k <- (nsims * (w-warmup-1) + 1) : (nsims * (w-warmup))

    for (day in 1:(7*length(start)))
    {
      Rt_offset <- P$Rt_offset_alpha * Rt_offset +
                   P$Rt_offset_sd * sqrt(1-P$Rt_offset_alpha^2) * randn()

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

        p <- inf * exp (log_Rt + Rt_offset)

        if (w > warmup)
        { wk <- ceiling(day/7)
          wsims[[vi]][k,wk] <- wsims[[vi]][k,wk] + p
        }

        past[[vi]][,past_next[vi]] <- p
        past_next[vi] <- past_next[vi] %% length(gen_interval) + 1
  
        t[[vi]] <- p + t[[vi]] * daily_decay[virus]
        tlt[[vi]] <- p + tlt[[vi]] * ltdaily_decay[virus]
      }

      # Save end of one of the years to initialize next five-year simulation.

      if ((day+5) %% 365 == 0)
      { wsave[w] <- wsave[w] - 1
        if (wsave[w] == 0)
        { sv_past <- past
          sv_past_next <- past_next
          sv_t <- t
          sv_tlt <- tlt
          # cat("Saved day",day,"\n")
        }
      }
    }
  
    # Set up initial state for next five-year simulations. Randomized a bit,
    # separately for each simulation.
  
    for (vi in 1:2)
    { # cat("start initializing from saved\n")
      n <- exp(0.2*randn())
      past[[vi]] <- sv_past[[vi]] * n
      past_next <- sv_past_next
      t[[vi]] <- sv_t[[vi]] * n
      tlt[[vi]] <- sv_tlt[[vi]] * n
    }

    sv_past <- sv_t <- sv_tlt <- NULL  # free memory
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
# are nsims*keep histories in total.

est_error_model <- function (twsims, init_err_alpha=0, init_err_sd=2,
                             verbose=FALSE)
{
  if (verbose) cat("\nEstimation of error model\n\n")

  err_alpha <- rep (init_err_alpha, length=2)
  err_sd <- rep (init_err_sd, length=2)

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
    { cat ("  err_alpha",round(err_alpha,6),
           ": err_sd",round(err_sd,3),
           ": log likelihood",
                round (log_lik (twsims, err_alpha, err_sd, full=nsims*keep), 3),
           "\n")
    }
  }

  list (err_alpha=err_alpha, err_sd=err_sd)
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
  est <- est_error_model (twsims)
  log_lik (twsims, est$err_alpha, est$err_sd, ...)
}
  

# FREE MEMORY.

wsims <- twsims <- NULL 
wsims_subset <- twsims_subset <- NULL
wsims_new <- twsims_new <- NULL


# DO SIMULATIONS WITH ORIGINAL PARAMETER ESTIMATES (FROM FILE).

RNGversion("2.15.1")
seed <- 1

# Rprofmemt (nelem=2*nsims*keep+1)

cat ("SIMULATIONS WITH ORIGINAL PARAMETER ESTIMATES\n\n")

start_time <- proc.time()

wsims <- run_sims (nsims, warmup, keep, info=TRUE)

twsims <- itrans_wsims (wsims)

em <- est_error_model(twsims,verbose=TRUE)
err_alpha <- em$err_alpha
err_sd <-em$err_sd

errors <- sim_errors(twsims,err_alpha,err_sd)
pp <- pprob(errors)
# print(errors)
# print(pp)

cat ("\nHighest posterior probabilities:\n")
print (round(sort(pp,decreasing=TRUE)[1:16],6))

ll <- log_lik(twsims,err_alpha,err_sd,errors=errors)
cat ("\nLog likelihood,", length(pp), "simulations:", round(ll,3), "\n")

wmx <- which.max(pp)

cat ("\nRunning subset simulation\n\n")

# print(sort(pp,decreasing=TRUE)[1:sub])
# print(order(pp,decreasing=TRUE)[1:sub])

high <- unique ((order(pp,decreasing=TRUE)[1:sub]-1) %% nsims + 1)
subn <- length(high)

# print(high)

cache <- new.env()
wsims_subset <- run_sims (subn, warmup, keep, full=nsims, subset=high, 
                          cache=cache, info=TRUE)

if (FALSE)  # enable to check that caching works
{
  cat ("\nRunning subset simulation a second time\n\n")
  
  wsims_subset2 <- run_sims (subn, warmup, keep, full=nsims, subset=high, 
                             cache=cache, info=TRUE)
  
  stopifnot(identical(wsims_subset,wsims_subset2))
  wsims_subset2 <- NULL
}
  
twsims_subset <- itrans_wsims (wsims_subset)

# print(rep(high,each=keep)+(0:(keep-1))*nsims)
# print(pp[rep(high,each=keep)+(0:(keep-1))*nsims])

cat("\n")

errors_subset <- sim_errors(twsims_subset,err_alpha,err_sd)
# print(errors_subset)
# print(pprob(errors_subset))

cat ("Log likelihood based on subset of",subn,"simulations:", 
      round(log_lik(twsims_subset,err_alpha,err_sd,full=nsims*keep),3),
      "\n\n")

# ESTIMATE MODEL PARAMETERS.

cat("\nESTIMATING MODEL PARAMETERS\n\n")

start_time_est <- proc.time()

N_evals <- 0

opt <- function (x)
{ N_evals <<- N_evals + 1
  P <- P_init
  P$Rt_offset_sd <- exp(x[1])
  P$Rt_offset_alpha <- tanh(x[2])
  -profile_log_lik (full=nsims*keep, itrans_wsims (run_sims (subn, warmup, keep,
    full=nsims*keep, subset=high, P=P, cache=cache, info=FALSE)))
}

P_new <- P_init
estim <- nlm (opt, c (log(P_init$Rt_offset_sd), atanh(P_init$Rt_offset_alpha)),
              fscale=-ll, stepmax=0.1, steptol=1e-30, gradtol=1e-10, iterlim=50,
              print.level=2) $ estimate
P_new$Rt_offset_sd <- exp(estim[1])
P_new$Rt_offset_alpha <- tanh(estim[2])

cat("Number of function evaluations:",N_evals,"\n")
cat("Processing time for estimation:\n")
print(proc.time()-start_time_est)

cat("\nNew parameters:\n\n")
print_model_parameters(P_new)

# DO SIMULATIONS WITH NEW PARAMETERS.

cat ("SIMULATIONS WITH NEW PARAMETER ESTIMATES\n\n")

wsims_new <- run_sims (nsims, warmup, keep, P=P_new, info=TRUE)

twsims_new <- itrans_wsims (wsims_new)

em_new <- est_error_model(twsims_new,verbose=TRUE)

errors_new <- sim_errors(twsims_new,em_new$err_alpha,em_new$err_sd)
pp_new <- pprob(errors_new)
# print(errors_new)
# print(pp_new)

wmx_new <- which.max(pp_new)

cat ("\nHighest posterior probabilities:\n")
print (round(sort(pp_new,decreasing=TRUE)[1:16],6))

ll_new <- log_lik (twsims_new, em_new$err_alpha, em_new$err_sd)

cat ("\nLog likelihood for new parameters,", length(pp_new), "simulations",
      round(ll_new,3), "\n")

cat("GC after post-simulation work:\n")
print(gc())
cat("Total processing time for this group of viruses:\n")
print(proc.time()-start_time)


# PLOTS FOR THIS VIRUS GROUP.  Corresponding plots use the same scales.

par(mfrow=c(4,1))

yupper <-  max (proxy[[1]], proxy[[2]], 
                wsims[[1]][c(wmx,1:n_plotted),], 
                wsims[[2]][c(wmx,1:n_plotted),],
                wsims_new[[1]][wmx_new,],
                wsims_new[[2]][wmx_new,])
ylower <-  min (proxy[[1]], proxy[[2]], exp(-4),
                wsims[[1]][wmx,],
                wsims[[2]][wmx,],
                wsims_new[[1]][wmx_new,],
                wsims_new[[2]][wmx_new,])

# Observed proxies.

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

for (s in 0:n_plotted)
{
  plot (start, rep(0,length(start)),
        ylim=c(0,1.02*yupper), yaxs="i", type="n", 
        ylab="Simulated incidence")

  ss <- if (s==0) wmx else s
  lines (start, wsims[[1]][ss,], col="blue")
  lines (start, wsims[[2]][ss,], col="red")

  if (s==0) 
  { 
    title (paste ("Best fit simulation out of",nsims,"x",keep,
                  "for",virus_group[1],"and",virus_group[2]))

    plot (start, rep(0,length(start)),
      ylim=itrans(c(0.98*ylower,1.02*yupper)), yaxs="i", type="n", 
      ylab=paste(if (itrans_arg!="identity") itrans_arg, "simulated incidence"))
  
    lines (start, twsims[[1]][wmx,], col="blue")
    lines (start, twsims[[2]][wmx,], col="red")

    plot (sort(pp,decreasing=TRUE)[1:20],pch=20,xlab="",
          ylab="posterior prob")
    title (paste (
      "Posterior probabilities of most likely histories, log likelihood", 
       round(ll,3)))
    plot (log10(sort(pp,decreasing=TRUE)[1:20]),pch=20,xlab="",ylim=c(-30,0),
          ylab="log10 of posterior prob")
  }
  if (s==1) 
  { title (paste ("Other simulations of",virus_group[1],"and",virus_group[2]))
  }
}

# Observed proxies, again.

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

# Best fit simulation with new parameters.

plot (start, rep(0,length(start)),
      ylim=c(0,1.02*yupper), yaxs="i", type="n", 
      ylab="Simulated incidence")

lines (start, wsims_new[[1]][wmx_new,], col="blue")
lines (start, wsims_new[[2]][wmx_new,], col="red")

title (paste ("Best fit simulation with new parameters for",
              virus_group[1],"and",virus_group[2]))

plot (start, rep(0,length(start)),
  ylim=itrans(c(0.98*ylower,1.02*yupper)), yaxs="i", type="n", 
  ylab=paste(if (itrans_arg!="identity") itrans_arg, "simulated incidence"))
  
lines (start, twsims_new[[1]][wmx_new,], col="blue")
lines (start, twsims_new[[2]][wmx_new,], col="red")

plot (sort(pp_new,decreasing=TRUE)[1:20],pch=20,xlab="",
      ylab="posterior prob")
title (paste (
  "Posterior probabilities of most likely histories, log likelihood", 
   round(ll_new,3)))
plot (log10(sort(pp_new,decreasing=TRUE)[1:20]),pch=20,xlab="",ylim=c(-30,0),
      ylab="log10 of posterior prob")

# ----- END OF LOOP OVER THE TWO VIRUS GROUPS -----

}

# ALL DONE.

dev.off()

print(proc.time())
