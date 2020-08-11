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

P_new <- P_init
pp_new <- pp

if (TRUE)  # optimization can be disabled for debugging
{ 
  eta <- 0*P_init + 5e-4
  eta$mc_viral[1:6] <- 5e-5
  eta$imm_decay <- 5e-5
  eta$ltimm_decay <- 1e-5
  eta$Rt_offset["alpha"] <- 5e-5
  eta$Rt_offset["sd"] <- 1e-5

  alpha <- 0.994

  p <- 0*P_init

  ws <- run_sims (subn, full=nsims, subset=high, P=P_new, 
                  cache=cache, info=FALSE)
  tws <- itrans_wsims(ws)
  H_prev <- - profile_log_lik (tws, full=nsims)
  p_prev <- p
  P_new_prev <- P_new

  full_rate <- 15
  start_momentum <- 20

  for (iter in 1:n_iter)
  { 
    # At intervals of 'full_interval', simulate a full set of 'nsims'
    # simulations, and reselect the high-probability subset.

    if (iter > 1 && iter%%full_interval == 1)
    {
      ws <- run_sims (subn, full=nsims, subset=high, P=P_new, 
                      cache=cache, info=FALSE)
      tws <- itrans_wsims(ws)
      cat("Log likelihood from current subset of size",subn,":",
           profile_log_lik(tws,full=nsims), "\n")
      cat("Redoing full set of simulations\n")
      wsims_new <- run_sims (nsims, P=P_new, info=FALSE)
      twsims_new <- itrans_wsims (wsims_new)

      em_new <- est_error_model(twsims_new)
      errors_new <- sim_errors(twsims_new,em_new$alpha,em_new$sd)
      pp_new <- pprob(errors_new)
      highL <- logical(nsims); highL[high] <- TRUE
      high <- unique ((order(pp_new,decreasing=TRUE)[1:sub]-1) %% nsims + 1)
      cat("Number of top simulations in common:",sum(highL[high]),"\n")
      highL <- NULL  # free memory
      subn <- length(high) # currently always equal to sub, but wasn't before
      ll_new <- log_lik (twsims_new, em_new$alpha, em_new$sd)

      cat("Lowest probability in top",subn,"is",pp_new[high[subn]],"\n")
      cat("Total probability is",round(sum(pp_new[high]),5),"\n")
      cat("Log likelihood from",nsims,"simulations:",round(ll_new,5),"\n")

      cache <- new.env()
      ws <- run_sims (subn, full=nsims, subset=high, P=P_new, 
                      cache=cache, info=FALSE)
      tws <- itrans_wsims(ws)
      cat("Log likelihood from new subset of size",subn,":",
           profile_log_lik(tws,full=nsims), "\n")
    }

    # Compute new negative log likelihood, with gradient attached.

    nll <- with gradient (P_new)
    { ws <- run_sims (subn, full=nsims, subset=high, P=P_new, 
                      cache=cache, info=FALSE)
      tws <- itrans_wsims(ws)
      - profile_log_lik (tws, full=nsims)
    }

    # Do a gradient descent with momentum update.

    this_eta <- if (iter<full_rate) 0.2*eta else eta
    this_alpha <- if (iter<start_momentum) 0 else alpha

    g <- attr(nll,"gradient")       # Update is done in Hamiltonian
    p <- p - this_eta * g           #   dynamics fashion, so that 'energy'
    P_new <- P_new + this_eta * p   #   would be conserved if exact

    K <- sum(sapply(p^2,sum))/2
    H <- nll + K
    
    cat ("iteration",iter,
         ": log lik",round(-nll,5),
         ": K",round(K,5),
         ": energy",round(H,5))
    if (this_alpha>0)
    { cat(" ->",round(nll+this_alpha^2*K,5))
    }
    cat("\n")

    p <- alpha * p

    # Reduce 'eta' if the energy went up too much.

    if (H > H_prev+0.3)
    { cat("Reducing eta by factor of 0.75 - from",eta,"to",0.8*eta,
          "- and backtracking\n")
      eta <- 0.8 * eta
      P_new <- P_new_prev
      p <- p_prev
    }
    else
    { P_new_prev <- P_new
      p_prev <- p
      H_prev <- H
    }
  }

  cat("\nChange in parameter values:\n\n")
  print_model_parameters(P_new-P_init)
}