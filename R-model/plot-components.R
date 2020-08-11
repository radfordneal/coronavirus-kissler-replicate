# Code to produce plots of components of a regression model for R.
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


# PLOT COMPONENTS OF MODEL FOR EACH VIRUS AND SEASON.  Similar to plots
# of Kissler, et al. Figure 1.

plot_components <- function (model, model_df, s, virus, logarithmic=FALSE, ...)
{
  trans <- if (logarithmic) log else identity
  itrans <- if (logarithmic) identity else exp
  mc <- coef(model)
  this <- model_df$season==s & model_df$virus==virus
  df <- model_df[this,]
  plot (c(1,season_length), trans(c(1,1)), type="n",
        ylim=trans(c(0.5,2.0)), xlim=c(0,season_length), xaxs="i", xlab="", 
        ylab = if (logarithmic) "Additive effect on log(R)"
               else "Multiplicative effect on R", ...)
  abline(h = if (logarithmic) c(-0.6,-0.4,-0.2,0.2,0.4,0.6) 
             else c(0.5,0.75,1.25,1.5,1.75,2.0),
         col="gray", lwd=0.5)
  if (season_type=="s1") 
  { abline (v = c(5,10,15,20,25,30), col="gray", lwd=0.5)
  }
  abline (v=1, h=trans(1))
  points (trans(model_df[,paste0(virus,"_R")][this]), pch=20, cex=0.75,
          col="pink")
  if (season_type=="s2")  # Show Kissler, et al season as dotted lines
  { abline(v=40+1-start_season,col="gray",lty=3)
    abline(v=40+33-start_season,col="gray",lty=3)
  }
  lines (itrans (predict(model,model_df)[this]), col="red", lwd=2)
  if (seffect_type=="e3")
  { tn <- ncol(trend_spline)
    trend_component <- as.vector (trend_spline[this,] %*% mc[1:tn])
    lines (itrans(trend_component), col="green", lwd=1.5)
  }
  seasonal_component <- switch (seffect_type,
      e1 = as.vector (seasonal_spline[1:season_length,] 
                       %*% mc[1:ncol(seasonal_spline)]),
      e2 = seffect_e2 (df$yrs),
      e3 = seffect_e3 (df$yrs)
    )
  mc0 <- mc
  if (season_type=="s1" && seffect_type!="e1")
  { mc0[paste0(virus,"_overall")] <- 
      mc0[paste0(virus,"_overall")] + seasonal_component[1]
    seasonal_component <- seasonal_component - seasonal_component[1]
  }
  lines (itrans(seasonal_component), col="orange", lwd=1.5)
  same <- df[,paste0(virus,"_same")] * mc[paste0(virus,"_same")]
  lines (itrans(same), col="black", lwd=1.5)
  if (immune_type!="i4")
  { other <- df[,paste0(virus,"_other")] * mc[paste0(virus,"_other")]
    lines (itrans(other), col="gray", lwd=1.5)
  }
  if (immune_type=="i3" || immune_type=="i4")
  { samelt <- df[,paste0(virus,"_samelt")] * mc[paste0(virus,"_samelt")]
    lines (itrans(samelt), col="black", lwd=1.5, lty=2)
    otherlt <- df[,paste0(virus,"_otherlt")] * mc[paste0(virus,"_otherlt")]
    lines (itrans(otherlt), col="gray", lwd=1.5, lty=2)
  }
  if (immune_type=="i1" && seffect_type=="e1")
  { points (1, itrans (mc[paste0(virus,"_season_",s)]
                     - mc[paste0(virus_group[2],"_season_2014")]), 
               pch=20, col="black", cex=0.75)
    points (1, itrans (mc[paste0(virus,"_season_",s)]), 
               pch=20, col="green", cex=0.75)
  }
  else
  { points (1, itrans (mc0[paste0(virus,"_overall")]), 
               pch=20, col="darkgreen", cex=0.75)
  }
  title(paste0(virus," ",s,"-",s+1," ",R_estimates,"   "),line=0.5)
}
