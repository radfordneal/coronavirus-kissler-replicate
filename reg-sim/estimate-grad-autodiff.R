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

source("concat.R")

cat("Number of parameters:",sum(sapply(P_init,length)),"\n")

P_new <- P_init
pp_new <- pp

if (TRUE)  # optimization can be disabled for debugging
{ 
  eta <- 0*P_init + 5e-4
  eta$mc_trend <- 5e-5
  eta$mc_seasonality <- 8e-5
  eta$mc_viral[1:10] <- 2e-5
  eta$mc_viral[11:12] <- 2e-5
  eta$imm_decay <- 7e-3
  eta$ltimm_decay <- 7e-3
  eta$lt2imm_decay <- 7e-3
  eta$imm_initial <- 1e-2
  eta$ltimm_initial <- 1e-2
  eta$lt2imm_initial <- 1e-2
  eta$Rt_offset["alpha"] <- 3e-3
  eta$Rt_offset["sd"] <- 2e-3

  cat("Value for eta:\n")
  print(eta)

  eta_adj <- 2.0
  alpha <- 0.997
  cat("Eta ajustment:",eta_adj," Momentum:",alpha,"\n")

  if (!is.null(momentum))
  { full_rate <- 1
    start_momentum <- 1
    full_momentum <- 1
  }
  else
  { full_rate <- 10        # When to switch from smaller to full eta
    start_momentum <- 20   # When to swith from zero to small momentum
    full_momentum <- 50    # When to swith from small to full momentum
  }

  p <- if (is.null(momentum)) 0*P_init else momentum

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

  rec <- vector("list",opt_iters+1)
  recp <- vector("list",opt_iters+1)

  for (iter in 1:opt_iters)
  { 
    rec[[iter]] <- P_new
    recp[[iter]] <- p

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
    # eta_adj is too big, allowing for change to eta_adj when needed.

    this_eta <- eta_adj * if (iter<full_rate) 0.2*eta else eta
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

    # Reduce 'eta_adj' if H went up too much.

    if (H-H0 > 0.2)
    { eta_adj <- 0.8 * eta_adj
      P_new <- P_new_prev
      p <- p_prev
      nll <- neg_ll(P_new)
      cat("Multiplying eta_adj by 0.8 - now",eta_adj,
          "- and backtracking - new H is",
           sprintf("%.5f",nll + sum(sapply(p^2,sum))/2),"\n")
    }
    else
    { P_new_prev <- P_new
      p_prev <- p
    }
  }

  # Plot changes in parameters, and in corresponding momentum variables.

  par(mfrow=c(4,4))

  rec[[opt_iters+1]] <- P_new
  recp[[opt_iters+1]] <- p
  rec <- concat(rec)
  recp <- concat(recp)

  for (i along(rec))
  { plot(c(0,opt_iters),range(rec[[i]]),type="n",xlab="",ylab="")
    title(names(rec)[[i]])
    for (j in 1..nrow(rec[[i]]))
    { lines(0..opt_iters,rec[[i]][j,],col=1+j)
    }
    plot(c(0,opt_iters),range(recp[[i]]),type="n",xlab="",ylab="")
    abline(h=0)
    for (j in 1..nrow(recp[[i]]))
    { lines(0..opt_iters,recp[[i]][j,],col=1+j)
    }
  }

  momentum <- p
}
