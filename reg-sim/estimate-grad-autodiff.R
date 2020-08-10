# Code for model parameter estimation by gradient descent, 
# using autodiff (pqR only).
#
# Updates parameters in P_init to new values in P_new.
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

cat("Number of parameters:",sum(sapply(P_init,length)),"\n")

if (TRUE)  # can be disabled for debugging
{ 
  eta <- 0*P_init + 5e-4
  eta$mc[18:23] <- 5e-5
  eta$imm_decay <- 1e-6
  eta$ltimm_decay <- 1e-7
  eta$Rt_offset_alpha <- 1e-6
  eta$Rt_offset_sd <- 1e-7

  alpha <- 0.99

  p <- 0*P_init
  P_new <- P_init

  for (iter in 1:100)
  { 
    nll <- with gradient (P_new)
    { ws <- run_sims (subn, full=nsims, subset=high, P=P_new, 
                      cache=cache, info=FALSE)
      tws <- itrans_wsims(ws)
      - profile_log_lik (tws, full=nsims)
    }

    this_eta <- if (iter<15) 0.2*eta else eta
    this_alpha <- if (iter<15) 0 else if (iter<30) alpha^2 else alpha

    g <- attr(nll,"gradient")
    p <- p - this_eta * g
    P_new <- P_new + this_eta * p

    K <- sum(sapply(p^2,sum))/2
    
    cat ("iteration",iter,
         ": log lik",round(-nll,5),
         ": K",round(K,5),
         ": energy",round(nll+K,5))
    if (this_alpha>0)
    { cat(" ->",round(nll+this_alpha^2*K,5))
    }
    cat("\n")

    p <- if (iter<20) 0 else alpha*p
  }

  cat("\nChange in parameter values:\n")
  print(P_new-P_init)
}
