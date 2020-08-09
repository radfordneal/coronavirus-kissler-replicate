# Code for model parameter estimation by gradient descent, 
# using autodiff (pqR only).
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
  
    P$mc[18:23] <- x / 100
  
    # cat("OPT CALLED, with\n  ",round(P$mc,3),"\n")
    # cat("diff from P_init is\n  ",P$mc-P_init$mc,"\n")
  
    ws <- run_sims (subn, full=nsims, subset=high, P=P, cache=cache, info=FALSE)
    tws <- itrans_wsims(ws)
    result <- - profile_log_lik (tws, full=nsims)
    # cat("result is",result,"\n")
    # cat("gradient of result is",gradient_of(result),"\n")
  
    result
  }
}

# x_init <- c (log(P_init$Rt_offset_sd),
#              atanh(P_init$Rt_offset_alpha))

x_init <- P_init$mc[18:23] * 100

if (TRUE)  # can enable for debugging
{ cat("Initial parameter vector:\n")
  print(x_init)
  opt_value <- opt(x_init)
  cat("Initial opt value:",round(opt_value,3),"\n")
  cat("Gradient of initial opt value:",attr(opt_value,"gradient"),"\n")
  if (FALSE)
  { delta <- 1e-7
    cat("With",delta,"changes:\n")
    for (i in 1:length(x_init)) 
    { print(opt(x_init+delta*((1:length(x_init))==i)))
    }
    for (i in 1:length(x_init)) 
    { print(opt(x_init-delta*((1:length(x_init))==i)))
    }
  }
  cat("\n")
  N_evals <- 0
}

P_new <- P_init

if (FALSE) x_new <- x_init  # for debugging
else
{ eta <- 4e-5
  alpha <- 0.95
  x_new <- x_init
  p <- 0
  for (iter in 1:100)
  { o <- opt(x_new)
    p <- alpha*p - eta * attr(o,"gradient")
    x_new <- x_new + p
    cat("iteration",iter,"value",round(o,5),"\n")
  }
}

# P_new$Rt_offset_sd <- exp(estim[1])
# P_new$Rt_offset_alpha <- tanh(estim[2])

P_new$mc[18:23] <- x_new / 100

cat("Number of function evaluations:",N_evals,"\n")
