# Code for model parameter estimation using nlm, with autodiff (pqR only).
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

N_evals <- 0

opt <- function (x)
{ 
  N_evals <<- N_evals + 1

  with gradient (x)
  { 
    P <- P_init
    # P$Rt_offset_sd <- exp(x[1])
    # P$Rt_offset_alpha <- tanh(x[2])
  
    # P$mc[] <- x  # keeps names
    # P$mc[18:23] <- x[18:23] / 100
    P$mc[18] <- x / 100
  
    cat("OPT CALLED, with\n  ",round(P$mc,3),"\n")
    cat("diff from P_init is\n  ",P$mc-P_init$mc,"\n")
  
    ws <- run_sims (subn, full=nsims*keep, subset=high, P=P, 
                    cache=cache, info=FALSE)
    tws <- itrans_wsims(ws)
    result <- - profile_log_lik (tws, full=nsims)
    cat("result is",result,"\n")
    cat("gradient of result is",gradient_of(result),"\n")
  
    result
  }
}
options(digits=15)
# x_init <- c (log(P_init$Rt_offset_sd),
#              atanh(P_init$Rt_offset_alpha))

# x_init <- P_init$mc
# x_init[18:23] <- P_init$mc[18:23] * 100
x_init <- P_init$mc[18] * 100

if (TRUE)  # can enable for debugging
{ cat("Initial parameter vector:\n")
  print(x_init)
  opt_value <- opt(x_init)
  cat("Initial opt value:",round(opt_value,3),"\n")
  cat("Gradient of initial opt value:",attr(opt_value,"gradient"),"\n")
  if (TRUE)
  { delta <- 1e-13
    cat("With changes:\n")
    for (i in 1:length(x_init)) 
    { print(opt(x_init+delta*((1:length(x_init))==i)))
    }
    for (i in 1:length(x_init)) 
    { print(opt(x_init-delta*((1:length(x_init))==i)))
    }
  }
  N_evals <- 0
}

P_new <- P_init

if (TRUE) x_new <- x_init  # for debugging
else
{ x_new <- nlm (opt, x_init,
             fscale=-500, stepmax=0.2, steptol=1e-30, gradtol=1e-30, iterlim=50,
             check.analyticals=FALSE, print.level=2) $ estimate
}

# P_new$Rt_offset_sd <- exp(estim[1])
# P_new$Rt_offset_alpha <- tanh(estim[2])

# P_new$mc[] <- x_new  # keeps names
# P_new$mc[18:23] <- x_new[18:23] / 100
P_new$mc[18] <- x_new / 100

cat("Number of function evaluations:",N_evals,"\n")
