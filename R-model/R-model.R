# Code to reproduce the model for estimated R values of common cold
# coronaviruses of the Kissler, et al paper, and to investigate
# alternative models.  
#
# Options are specified by arguments in the R command after --args, with 
# defaults as in Kissler, et al, as follows:
#
#   - The R estimates to use (default proxyW), from corresponding file
#   - Type of R estimate - Rt, Rt-smoothed, Ru, or Ru-smoothed (default)
#   - Model of seasonal effect - e1 (default, spline), e2 (sine), e3 (Fourier)
#   - Flu season - s1 (default) or s2 (52 weeks, usually the whole year)
#   - Model of "immunity" - i1 (default), i2 (exp decay), i3 (short & long),
#                           i4 (like i3 but with no short-term cross-immunity),
#                           i5 (like i4 but with two-stage immunity)
#   - Whether heteroskedasticity w.r.t. virus is modelled - het for "yes" 
#     (default "no")
#   - Immune decay constants for i2 or i3/i4/i5 model (defaults in code).
#     For example: decay:NL63=0.9,E229=0.8,OC43=0.7,HKU1=0.6
#     May override just a subst of the values
#   - Long-term immune decay constants for i3/i4/i5 model (defaults in code).
#     For example: ltdecay:NL63=0.91,E229=0.92,OC43=0.93,HKU1=0.94
#     May override just a subst of the values
#   - 2nd long-term immune decay constants for i5 model (defaults in code).
#     For example: lt2decay:NL63=0.92,E229=0.93,OC43=0.94,HKU1=0.95
#
# Produces various plots that are written to the file with name
# R-model-<R-estimates>-<R-estimate-type>[-<sn>][-<in>][-<en>][-het].pdf.  
# Information on the model fit is written to standard output. The 
# models themselves are written to files with names of the form
# R-model-<R-estimates>-<R-estimate-type>[-<sn>][-<in>][-<en>][-het]-<g>.model,
# where <g> is the virus group, either "alpha" or "beta".
#
# The e1 seasonal effect model uses a spline, as in the regression model
# of Kissler, et al.  The e2 seasonal effect model uses a single sine 
# wave (expressed as a combination of sin and cos), as in the SEIRS model
# of Kissler, et al.  The e3 seasonal effect model expresses the effect
# as a Fourier series, with sin and cos up to frequency 6; it also has a
# long-term trend component, expressed as a spline with knots at the ends
# of the period and at 1/3 and 2/3 way between.
#
# If the type of R estimate is Rt and the options are s2, i2 or i3, and 
# e2 or e3, simulatons are done using the fitted model parameters, allowing
# the implied dynamics to be seen.  There is randomness in the Rt values
# and in the state at the beginning of the five-year period.
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
library(sandwich)  # needs to be installed with install.packages("sandwich")

options(warn=1)


# ESTABLISH WHICH PARAMETERS TO USE, LOOKING AT R'S ARGUMENTS.

args <- commandArgs(trailing=TRUE)

cat("ARGUMENTS:",args,"\n\n")

getarg <- function (what, default="")
{ if (any (args %in% what))
  { stopifnot (sum (args %in% what) == 1)
    res <- args [args %in% what]
    args <<- args [! (args %in% what)]
    res
  }
  else
    default
}

R_est_type <- getarg (c("Rt","Rt_smoothed","Ru","Ru_smoothed"), "Ru_smoothed")

immune_type <- getarg (c("i1","i2","i3", "i4", "i5"), "i1")
seffect_type <- getarg (c("e1","e2","e3"), "e1")
season_type <- getarg (c("s1","s2"), "s1")

het_virus <- getarg ("het") == "het"

# The default values below for imm_decay with immune_type of "i2" were
# found using the "search" script, with proxyDss-filter proxies.

imm_decay <- ( if (immune_type == "i2")
                 c (NL63=0.9175, E229=0.9850, OC43=0.9500, HKU1=0.9750)
               else
                 c (NL63=0.85, E229=0.85, OC43=0.85, HKU1=0.85) )

ltimm_decay <- ( if (immune_type == "i5")
                   c (NL63=0.98, E229=0.98, OC43=0.98, HKU1=0.98)
                 else
                   c (NL63=0.985, E229=0.985, OC43=0.985, HKU1=0.985)
               )
lt2imm_decay <- c (NL63=0.95, E229=0.95, OC43=0.95, HKU1=0.95)

if (any (substr(args,1,6) == "decay:"))
{ stopifnot (sum (substr(args,1,6) == "decay:") == 1)
  decay_arg <- substr (args [substr(args,1,6) == "decay:"], 7, 100)
  repl_imm_decay <- eval(parse(text=paste0("c(",decay_arg,")")))
  imm_decay[names(repl_imm_decay)] <- repl_imm_decay
  args <- args [substr(args,1,6) != "decay:"]
}

if (any (substr(args,1,8) == "ltdecay:"))
{ stopifnot (sum (substr(args,1,8) == "ltdecay:") == 1)
  ltdecay_arg <- substr (args [substr(args,1,8) == "ltdecay:"], 9, 100)
  repl_ltimm_decay <- eval(parse(text=paste0("c(",ltdecay_arg,")")))
  ltimm_decay[names(repl_ltimm_decay)] <- repl_ltimm_decay
  args <- args [substr(args,1,8) != "ltdecay:"]
}

if (any (substr(args,1,8) == "lt2decay:"))
{ stopifnot (sum (substr(args,1,9) == "lt2decay:") == 1)
  lt2decay_arg <- substr (args [substr(args,1,9) == "lt2decay:"], 10, 100)
  repl_lt2imm_decay <- eval(parse(text=paste0("c(",lt2decay_arg,")")))
  lt2imm_decay[names(repl_lt2imm_decay)] <- repl_lt2imm_decay
  args <- args [substr(args,1,9) != "lt2decay:"]
}

R_estimates <- args

if (length(R_estimates)==0) 
{ R_estimates <- "proxyW"  # As in the Kissler, et al paper
}

stopifnot(length(R_estimates)==1)

cat("imm_decay:\n"); print(imm_decay); cat("\n")
cat("ltimm_decay:\n"); print(ltimm_decay); cat("\n")

file_base <- paste0 ("R-model-",R_estimates,"-",R_est_type,
                     if (season_type!="s1") paste0("-",season_type),
                     if (immune_type!="i1") paste0("-",immune_type),
                     if (seffect_type!="e1") paste0("-",seffect_type),
                     if (het_virus) "-het")

# PLOT SETUP.

pdf (paste0 (file_base,".pdf"), height=8, width=6)
par(mar=c(1.5,2.3,3,0.5), mgp=c(1.4,0.3,0), tcl=-0.22)
yrcols <- c("red","green","blue","orange","darkcyan","darkmagenta")


# READ THE FILE OF CORONAVIRUS R ESTIMATES.

R_est <- read.csv (paste0("../R-est/R-est-",R_estimates,".csv"),
                   header=TRUE, stringsAsFactors=FALSE)

R_est$start <- as.Date(R_est$start)
R_est$yrs <- (0:(nrow(R_est)-1)) / (365.24/7)


# INCLUDE UTILITY FUNCTIONS.  They need "start" and "week" to be defined.

start <- R_est$start  # Start dates of weeks being analysed
year  <- R_est$year   # Year for each week
week  <- R_est$week   # Number of each week in its year

source("../util/util.R")


# INCLUDE FUNCTIONS THAT COMPUTE SEASONAL EFFECTS AND LONG-TERM TREND.

source("steffect.R")


# ----- DO EVERYTHING FOR BOTH ALPHACORONAVIRUSES AND BETACORONAVIRUSES -----

for (g in seq_along(virus_groups)) {

virus_group <- virus_groups[[g]]


# CREATE A DATA FRAME WITH ONLY THE SELECTED R ESTIMATES, FOR FLU SEASON ONLY.
# The new data frame, select_df, has columns start, year, week, season (year
# of start of flu season), season_week (week in the flu season, from 1),
# <virus>_R (for R estimates), <virus>_proxy (for incidence proxy), and
# yrs (year from 0, varying continuously).

# Default, as in Kissler, et al.

start_season <- 40 # Week of start of "flu season"
end_season <- 20   # Week of end of "flu season" (in year following start)
                   # (except ends week earlier in 2015, since 2014 has 53 weeks)

if (season_type=="s2")  # Whole year, except when year is 53 weeks long
{ start_season <- 28
  end_season <- 27
}

R_est_names <- paste0(virus_group,"_",R_est_type)
proxy_names <- paste0(virus_group,"_proxy")

select_df <- R_est [, c("start", "year", "week", "yrs",
                        R_est_names, proxy_names)]
names(select_df) [names(select_df) %in% R_est_names] <- paste0(virus_group,"_R")

select_df$season <- select_df$year - (select_df$week<=end_season)

par(mfrow=c(2,1))

in_season <- select_df$week >= start_season | 
             select_df$week <= end_season-(select_df$year==2015)

for (virus in virus_group)
{ plot (start, select_df[,paste0(virus,"_R")], pch=20, 
        col=ifelse(in_season,"black","gray"), ylab=R_est_type)
  week_lines()
  abline(h=1)
  title (paste("Weekly R estimates for",virus))
  plot (start, log(select_df[,paste0(virus,"_R")]), pch=20, 
        col=ifelse(in_season,"black","gray"), ylab=paste("log",R_est_type))
  week_lines()
  abline(h=0)
}  

select_df <- select_df[in_season,]

select_df$season_week <- 
  ifelse (select_df$week>=start_season, select_df$week-start_season+1,
          select_df$week+53-start_season+(select_df$year==2015))


# COMPUTE CUMULATIVE SEASONAL INCIDENCE FOR EACH CORONAVIRUS IN GROUP.
# These are added to select_df as columns named <virus>_cum.  A second
# set of columns are added named <virus>_cumexp that are exponentially
# weighted sums from the start of the record (not paying attention to
# "seasons"), using the imm_decay values for each virus, wrapping
# around once so that the start has a sum obtained from the average of
# the last four starts-of-years.  Similary, columns named <virus>_cumexplt 
# that are exponentially-weights sums using ltimm_decay values are added,
# and columns named <virus>_cumexplt2 are added if immune_type is "i5".

par(mfrow=c(2,1))

for (virus in virus_group)
{ 
  p <- select_df[,paste0(virus,"_proxy")]
  w <- select_df$week
  c <- rep(0,length(p))
  for (i in seq_along(w)) 
  { if (w[i] == start_season)
    { c[i] = p[i]  
    }
    else
    { c[i] = p[i] + c[i-1]
    }
  }
  select_df[,paste0(virus,"_cum")] <- c
  plot (select_df$start, c, pch=20, ylab="", xlim=range(start))
  abline(h=0)
  abline(h=seq(0,1000,by=50),col="gray",lty=3)
  title (paste ("Cumulative seasonal incidence for",virus))

  p <- R_est[,paste0(virus,"_proxy")]
  c <- rep(0,length(p))
  clt <- rep(0,length(p))
  clt2 <- rep(0,length(p))
  t <- 0
  tlt <- 0
  tlt2 <- 0
  d <- imm_decay[virus]
  dlt <- ltimm_decay[virus]
  dlt2 <- lt2imm_decay[virus]
  for (i in rep(seq_along(p),2))
  { t <- d*t + p[i]
    tlt2 <- dlt2*tlt2 + (1-dlt)*tlt
    tlt <- dlt*tlt + p[i]
    c[i] <- t
    clt[i] <- tlt
    clt2[i] <- tlt2
    if (i==length(p))
    { t <- mean(c[i-52*(0:3)])
      tlt <- mean(clt[i-52*(0:3)])
      tlt2 <- mean(clt2[i-52*(0:3)])
    }
  }
  select_df[,paste0(virus,"_cumexp")] <- c[in_season]
  select_df[,paste0(virus,"_cumexplt")] <- clt[in_season]
  if (immune_type=="i5")
  { select_df[,paste0(virus,"_cumexplt2")] <- clt2[in_season]
  }
  plot (start, c, pch=20, col=ifelse(in_season,"black","gray"),
        ylim=c(0,max(c,clt)), 
        ylab=paste("cum incidence: black st - green lt",
                   if (immune_type=="i5") "- blue lt2"))
  points (start, clt, pch=20, col=ifelse(in_season,"green","lightgreen"))
  if (immune_type=="i5")
  { points (start, clt2, pch=20, col=ifelse(in_season,"blue","lightblue"))
  }
  abline(h=0)
  title (paste ("Exp. cum. incidences for",virus,
                 "- decay",imm_decay[virus],ltimm_decay[virus],
                 if (immune_type=="i5")lt2imm_decay[virus],"   "))
}


# FIT A JOINT MODEL OF LOG(R) FOR BOTH CORONAVIRUSES.  The spline
# modelling the seasonal effect is shared between the two strains.
# The e3 model also has a overall trend shared between the two strains.
# Other aspects are as if the two strains were modelled separately,
# which is accomplished by setting covariate values pertaining to one
# strain to zero when the response is for the other strain.

formula <- "log(R_value) ~ -1"

if (seffect_type=="e3")
{ formula <- paste(formula,"+ trend_spline")
}

switch (seffect_type,
  e1 = formula <- paste(formula,"+ seasonal_spline"),
  e2 = formula <- paste(formula,"+ sin(2*pi*yrs) + cos(2*pi*yrs)"),
  e3 = formula <- paste(formula,"+ sin(1*2*pi*yrs) + cos(1*2*pi*yrs)",
                                "+ sin(2*2*pi*yrs) + cos(2*2*pi*yrs)",
                                "+ sin(3*2*pi*yrs) + cos(3*2*pi*yrs)",
                                "+ sin(4*2*pi*yrs) + cos(4*2*pi*yrs)",
                                "+ sin(5*2*pi*yrs) + cos(5*2*pi*yrs)",
                                "+ sin(6*2*pi*yrs) + cos(6*2*pi*yrs)"
                  )
)

model_df <- NULL
for (virus in virus_group)
{ v_df <- select_df
  row.names(v_df) <- NULL
  v_df$virus <- virus

  v_df [, paste0(virus,"_same")] <- 
    v_df [, paste0 (virus, if (immune_type=="i1") "_cum" else "_cumexp")]
  v_df[,paste0(other_virus_of_type[virus],"_same")] <- 0
  formula <- paste (formula, "+", paste0(virus,"_same"))

  if (immune_type!="i4" && immune_type!="i5") 
  { v_df [, paste0(virus,"_other")] <- 
      v_df [, paste0 (other_virus_of_type[virus], 
                      if (immune_type=="i1") "_cum" else "_cumexp")]
    v_df[,paste0(other_virus_of_type[virus],"_other")] <- 0
    formula <- paste (formula, "+", paste0(virus,"_other"))
  }

  if (immune_type=="i3" || immune_type=="i4" || immune_type=="i5")
  { 
    v_df [, paste0(virus,"_samelt")] <- 
      v_df [, paste0 (virus, "_cumexplt")]
    v_df[,paste0(other_virus_of_type[virus],"_samelt")] <- 0
    formula <- paste (formula, "+", paste0(virus,"_samelt"))

    v_df [, paste0(virus,"_otherlt")] <-  
      v_df [, paste0 (other_virus_of_type[virus],"_cumexplt")]
    v_df[,paste0(other_virus_of_type[virus],"_otherlt")] <- 0
    formula <- paste (formula, "+", paste0(virus,"_otherlt"))
  }

  if (immune_type=="i5")
  { 
    v_df [, paste0(virus,"_samelt2")] <- 
      v_df [, paste0 (virus, "_cumexplt2")]
    v_df[,paste0(other_virus_of_type[virus],"_samelt2")] <- 0
    formula <- paste (formula, "+", paste0(virus,"_samelt2"))

    v_df [, paste0(virus,"_otherlt2")] <-  
      v_df [, paste0 (other_virus_of_type[virus],"_cumexplt2")]
    v_df[,paste0(other_virus_of_type[virus],"_otherlt2")] <- 0
    formula <- paste (formula, "+", paste0(virus,"_otherlt2"))
  }

  model_df <- if (is.null(model_df)) v_df else rbind(model_df,v_df)
}

for (virus in virus_group) 
{ if (immune_type=="i1" && seffect_type=="e1")
  { for (yr in unique(select_df$season))
    { ind <- paste0(virus,"_season_",yr)
      model_df[,ind] <- as.numeric (model_df$virus==virus & model_df$season==yr)
      formula <- paste (formula, "+", ind)
    }
  }
  else
  { ind <- paste0(virus,"_overall")
    model_df[,ind] <- as.numeric (model_df$virus==virus)
    formula <- paste (formula, "+", ind)
  }
}

R_value <- numeric()
for (virus in virus_group) 
{ R <- select_df [ , paste0(virus,"_R")]
  R_value <- c (R_value, R)
}

model_df$R_value <- R_value

season_length <- max(model_df$season_week)
seasonal_spline <- ( if (season_type=="s1")
    bs(model_df$season_week, knots=(1:7)*4+1, Boundary=c(1,season_length))
  else
    bs(model_df$season_week, knots=(1:9)*((season_length-1)/10)+1,
       Boundary=c(1,season_length))
  )

# Spline for slow trend over the years.

trend_spline <- make_trend_spline (
  if (season_type=="s1") model_df$yrs else R_est$yrs, model_df$yrs)

# Fit model, perhaps repeatedly, each time adjusting weights based on previous 
# fit to account for heterskedasticity w.r.t. virus.
  
cat ("\nMODEL FOR", paste(virus_group,collapse=" & "), "\n")
options(digits=9)

var_ratio_2over1 <- 1

for (rpt in if (het_virus) 1:3 else 1)
{
  model <- lm (parse(text=formula)[[1]], data=model_df, x=TRUE, y=TRUE,
               weights=rep(c(var_ratio_2over1,1),each=sum(in_season)))

  # print(summary(model))

  std.errs <- summary(model)$coefficients[,1:2]
  std.errs <- cbind (std.errs, 
    HC3.se = sqrt (diag (vcovHC (model, type="HC3"))))
  std.errs <- cbind (std.errs, 
    NW.se = sqrt (diag (NeweyWest (model, adjust=TRUE, verbose=TRUE))))
  cat("\n")
  std.errs <- cbind (std.errs, 
    NWlag5 = sqrt (diag (NeweyWest (model, adjust=TRUE, lag=5))))
  std.errs <- cbind (std.errs, 
    NWlag10 = sqrt (diag (NeweyWest (model, adjust=TRUE, lag=10))))
  std.errs <- cbind (std.errs, 
    NWlag20 = sqrt (diag (NeweyWest (model, adjust=TRUE, lag=20))))

  print(round(std.errs,5))
  cat("\n")

  # Special computation for standard error in sum of short and long term
  # immunity coefficients for i3/i4/i5 model, using Newey-West covariance
  # matrix of estimates. Also computes standard error for the difference
  # in seasonal effect between week 20 and week 50.

  if (immune_type %in% c("i3","i4","i5"))
  { 
    v <- NeweyWest (model, adjust=TRUE)
    mc <- summary(model)$coefficients[,1]

    for (virus in virus_group)
    { 
      wv <- paste0(virus,c("_same","_samelt", if (immune_type=="i5")"_samelt2"))
      vv <- v[wv,wv]
      cat("Sum of",wv[1],"and",wv[2],
          if (immune_type=="i5") paste("and",wv[3],":") else ":",
          round(sum(mc[wv]),5),
          "NW std err", round(sqrt(sum(vv)),5),
          "\n")
    }

    if (season_type=="s2" && seffect_type=="e3")
    {
      tn <- ncol(trend_spline)
      wv <- (tn+1):(tn+12)
      vv <- v[wv,wv]
      yrs <- 20/52
      sincos20 <- c (sin(1*2*pi*yrs), cos(1*2*pi*yrs),
                     sin(2*2*pi*yrs), cos(2*2*pi*yrs),
                     sin(3*2*pi*yrs), cos(3*2*pi*yrs),
                     sin(4*2*pi*yrs), cos(4*2*pi*yrs),
                     sin(5*2*pi*yrs), cos(5*2*pi*yrs),
                     sin(6*2*pi*yrs), cos(6*2*pi*yrs))
      yrs <- 50/52
      sincos50 <- c (sin(1*2*pi*yrs), cos(1*2*pi*yrs),
                     sin(2*2*pi*yrs), cos(2*2*pi*yrs),
                     sin(3*2*pi*yrs), cos(3*2*pi*yrs),
                     sin(4*2*pi*yrs), cos(4*2*pi*yrs),
                     sin(5*2*pi*yrs), cos(5*2*pi*yrs),
                     sin(6*2*pi*yrs), cos(6*2*pi*yrs))

      cat("\nSeasonality for week50 - week20",":",
           round (as.vector ((sincos50-sincos20) %*% mc[wv]), 5),
           "NW std err", round (sqrt (as.vector ((sincos50-sincos20) %*% vv %*%
                                                 (sincos50-sincos20))), 5), 
           "\n")
    }
  }
  
  resid <- log(R_value) - as.vector(predict(model,model_df))
  
  virus_residuals <- 
    list (resid[1:sum(in_season)], resid[(sum(in_season)+1):(2*sum(in_season))])
  names(virus_residuals) <- virus_group

  cat("\n")  
  for (virus in virus_group)
  { cat ("Residual standard deviation for", paste0(virus,":"),
          round (sd(virus_residuals[[virus]],na.rm=TRUE), 4), "\n")
    cat ("Residual ACF:", 
          round (acf (virus_residuals[[virus]],
                      na.action=na.pass, lag.max=11, plot=FALSE) $ acf, 2), 
         "\n\n")
  }

  var_ratio_2over1 <- mean(virus_residuals[[2]]^2,na.rm=TRUE) /
                      mean(virus_residuals[[1]]^2,na.rm=TRUE)
}


# PLOT MODEL RESIDUALS, AND THEIR ACF.  Residuals are obtained from predictions
# because this allows NAs to be handled correctly.

par(mfrow=c(2,1))

for (virus in virus_group)
{ plot (start[in_season], virus_residuals[[virus]], pch=20, ylab="")
  abline(h=0); week_lines()
  title(paste("Residuals for",virus,"- std. dev.",
               round(sd(virus_residuals[[virus]],na.rm=TRUE),4)))
  acf (virus_residuals[[virus]],na.action=na.pass,lag.max=35,main="")
}


# PLOT COMPONENTS OF MODEL.

source("plot-components.R")

mc <- coef(model)

par(mfcol=c(5,4))
sv <- par (cex.main=3/4, cex.lab=2/3, cex.axis=4/10, mgp=c(0.8,0.18,0))

for (virus in virus_group)
{ for (s in unique(model_df$season))
  { plot_components (mc, model$x, model_df, s, virus,logarithmic=FALSE)
  }
}

for (virus in virus_group)
{ for (s in unique(model_df$season))
  { plot_components (mc, model$x, model_df, s, virus, logarithmic=TRUE)
  }
}

par(sv)

# Again, but bigger, and in a different order.

for (virus in virus_group)
{ par(mfrow=c(3,2))
  for (s in unique(model_df$season))
  { plot_components (mc, model$x, model_df, s, virus, logarithmic=FALSE)
    plot_components (mc, model$x, model_df, s, virus, logarithmic=TRUE)
  }
}


# SAVE THE MODEL COEFFICIENTS AND DECAY VALUES TO A FILE.

saveRDS (list (mc_trend = mc [grepl("trend",names(mc))],
               mc_seasonality = mc [grepl("yrs",names(mc))],
               mc_viral = mc [!grepl("trend",names(mc))
                                & !grepl("yrs",names(mc))],
               imm_decay = imm_decay[virus_group],
               ltimm_decay = ltimm_decay[virus_group]),
         file = paste0(file_base,"-",names(virus_groups)[g],".model"),
         version=2)


# SAVE THE DATA CONTEXT USED FOR PLOTS.

saveRDS (list (model_x = model$x, model_df = model_df),
         file = paste0(file_base,"-",names(virus_groups)[g],".context"),
         version=2)


# ----- END OF LOOP OVER THE TWO VIRUS GROUPS -----

}

# ALL DONE.

dev.off()
