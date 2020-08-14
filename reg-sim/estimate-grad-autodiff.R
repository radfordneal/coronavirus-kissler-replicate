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
  eta$mc_trend <- 5e-4
  eta$mc_seasonality <- 5e-4
  eta$mc_viral[1:6] <- 5e-5
  eta$imm_decay <- 5e-4
  eta$ltimm_decay <- 1e-4
  eta$imm_initial <- 8e-4
  eta$ltimm_initial <- 4e-4
  eta$Rt_offset["alpha"] <- 5e-4
  eta$Rt_offset["sd"] <- 2e-4

  cat("Initial value for eta:\n")
  print(eta)

  alpha <- 0.995
  cat("Momentum:",alpha,"\n")

  full_rate <- 15        # When to switch from smaller to full eta
  start_momentum <- 25   # When to swith from zero to small momentum
  full_momentum <- 50    # When to swith from small to full momentum

  p <- 0*P_init

  # Compute negative log likelihood, with gradient attached.

  neg_ll <- function (P) with gradient (P)
  { ws <- run_sims (subn, full=nsims, subset=high, P=P, 
                    cache=cache, info=FALSE)
    tws <- itrans_wsims(ws)
    - profile_log_lik (tws, full=nsims)
  }

  H <- neg_ll (P_new)
  p_prev <- p
  P_new_prev <- P_new

  # A check on the gradient...

  cat("Gradient check:\n")
  delta <- 2e-4

  H2 <- neg_ll(P_new+eta*delta)
  cat ("  Change in H when adding eta *",delta,":",H2-H,"\n")
  gave <- (attr(H,"gradient")+attr(H2,"gradient"))/2
  cat ("  Predicted change from average gradient :",
          sum(sapply(gave*eta*delta,sum)),"\n")

  # Can selectively enable these...

  if (FALSE)
  { eta0 <- 0*eta; 
    eta0$mc_trend <- eta$mc_trend
    H3 <- neg_ll(P_new+eta0*delta)
    cat ("  Change in H when adding eta0 *",delta,"(trend) :",H3-H,"\n")
    ga <- (attr(H,"gradient")+attr(H3,"gradient"))/2
    cat ("  Predicted change from average gradient :",
            sum(sapply(ga*eta0*delta,sum)),"\n")
  }

  if (FALSE)
  { eta0 <- 0*eta; 
    eta0$mc_seasonality <- eta$mc_seasonality
    H3 <- neg_ll(P_new+eta0*delta)
    cat ("  Change in H when adding eta0 *",delta,"(season) :",H3-H,"\n")
    ga <- (attr(H,"gradient")+attr(H3,"gradient"))/2
    cat ("  Predicted change from average gradient :",
            sum(sapply(ga*eta0*delta,sum)),"\n")
  }

  if (FALSE)
  { eta0 <- 0*eta; 
    eta0$mc_viral <- eta$mc_viral
    H3 <- neg_ll(P_new+eta0*delta)
    cat ("  Change in H when adding eta0 *",delta,"(viral) :",H3-H,"\n")
    ga <- (attr(H,"gradient")+attr(H3,"gradient"))/2
    cat ("  Predicted change from average gradient :",
            sum(sapply(ga*eta0*delta,sum)),"\n")
  }

  if (FALSE)
  { eta0 <- 0*eta; 
    eta0$imm_decay <- eta$imm_decay
    H3 <- neg_ll(P_new+eta0*delta)
    cat ("  Change in H when adding eta0 *",delta,"(decay) :",H3-H,"\n")
    ga <- (attr(H,"gradient")+attr(H3,"gradient"))/2
    cat ("  Predicted change from average gradient :",
            sum(sapply(ga*eta0*delta,sum)),"\n")
  }

  if (FALSE)
  { eta0 <- 0*eta; 
    eta0$ltimm_decay <- eta$ltimm_decay
    H3 <- neg_ll(P_new+eta0*delta)
    cat ("  Change in H when adding eta0 *",delta,"(ltdecay) :",H3-H,"\n")
    ga <- (attr(H,"gradient")+attr(H3,"gradient"))/2
    cat ("  Predicted change from average gradient :",
            sum(sapply(ga*eta0*delta,sum)),"\n")
  }

  if (FALSE)
  { eta0 <- 0*eta; 
    eta0$imm_initial <- eta$imm_initial
    H3 <- neg_ll(P_new+eta0*delta)
    cat ("  Change in H when adding eta0 *",delta,"(initial) :",H3-H,"\n")
    ga <- (attr(H,"gradient")+attr(H3,"gradient"))/2
    cat ("  Predicted change from average gradient :",
            sum(sapply(ga*eta0*delta,sum)),"\n")
  }

  if (FALSE)
  { eta0 <- 0*eta; 
    eta0$ltimm_initial <- eta$ltimm_initial
    H3 <- neg_ll(P_new+eta0*delta)
    cat("  Change in H when adding eta0 *",delta,"(ltinitial) :",H3-H,"\n")
    ga <- (attr(H,"gradient")+attr(H3,"gradient"))/2
    cat ("  Predicted change from average gradient :",
            sum(sapply(ga*eta0*delta,sum)),"\n")
  }

  if (FALSE)
  { eta0 <- 0*eta; 
    eta0$Rt_offset <- eta$Rt_offset
    H3 <- neg_ll(P_new+eta0*delta)
    cat ("  Change in H when adding eta0 *",delta,"(offset) :",H3-H,"\n")
    ga <- (attr(H,"gradient")+attr(H3,"gradient"))/2
    cat ("  Predicted change from average gradient :",
            sum(sapply(ga*eta0*delta,sum)),"\n")
  }

  cat ("Average gradient:\n\n")
  print(gave)

  nll <- H

  for (iter in 1:opt_iters)
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

      nll <- neg_ll(P_new)
    }

    # Do a gradient-descent-with-momentum update.  It's done in
    # "leapfrog" style, so that H will be nearly preserved unless 
    # eta is too big, allowing for an adjustment of eta when needed.

    this_eta <- if (iter<full_rate) 0.2*eta else eta
    this_alpha <- ( if (iter<start_momentum) 0 
                    else if (iter<full_momentum) alpha^2 
                    else alpha )

    H0 <- nll + sum(sapply(p^2,sum))/2

    p <- p - (this_eta/2) * attr(nll,"gradient") 
    P_new <- P_new + this_eta * p
    nll <- neg_ll(P_new)
    p <- p - (this_eta/2) * attr(nll,"gradient") 

    K <- sum(sapply(p^2,sum))/2
    H <- nll + K
    
    cat ("iter",sprintf("%4d",iter),
       # ": H0",sprintf("%.5f",H0),
         ": dH",sprintf("%8.5f",H-H0),
         ": ll",sprintf("%.5f",-nll),
         ": K",sprintf("%.5f",K),
         ": H",sprintf("%.5f",H))

    p <- this_alpha * p

    cat (" >",sprintf("%.5f",nll+sum(sapply(p^2,sum))/2), "\n")

    # Reduce 'eta' if H went up too much.

    if (H-H0 > 0.2)
    { eta <- 0.8 * eta
      P_new <- P_new_prev
      p <- p_prev
      nll <- neg_ll(P_new)
      cat("Reducing eta by factor of 0.8 and backtracking - new H is",
           sprintf("%.5f",nll + sum(sapply(p^2,sum))/2),"\n")
    }
    else
    { P_new_prev <- P_new
      p_prev <- p
    }
  }

  cat("\nChange in parameter values:\n\n")
  print_model_parameters(P_new-P_init)
}
