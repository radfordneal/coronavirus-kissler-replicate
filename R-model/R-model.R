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
#   - Flu season - s1 (default) or s2 (almost the whole year)
#   - Model of short-term "immunity" - i1 (default), i2 (exp decay), 
#                                      i3 (two-stage)
#   - Model for long-term immunity - I1 (default, none), I2 (exp decay)
#   - Whether heteroskedasticity w.r.t. virus is modelled - het for "yes" 
#     (default "no")
#   - Immune decay constants for i2 model (default in code, different for I2).
#     For example: decay:NL63=0.96,E229=0.95,OC43=0.97,HKU1=0.98
#   - Whether simulations are run (when possible) - nosim for "no", default
#     is "yes"
#
# Produces various plots that are written to the file with name
# R-model-<R-estimates>-<R-estimate-type>[-<sn>][-<in>][-<en>][-het].pdf.  
# Information on the model fit is written to standard output.
#
# The e1 seasonal effect model uses a spline, as in the regression model
# of Kissler, et al.  The e2 seasonal effect model uses a single sine 
# wave (expressed as a combination of sin and cos), as in the SEIRS model
# of Kissler, et al.  The e3 seasonal effect model expresses the effect
# as a Fourier series, with sin and cos up to frequency 6; it also has a
# long-term trend component, expressed as a spline with knots at the ends
# of the period and at 1/3 and 2/3 way between.
#
# If the type of R estimate is Rt and the options are s2, i2, and either 
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


# ESTABLISH WHICH PARAMETERS TO USE, LOOKING AT R'S ARGUMENTS.

args <- commandArgs(trailing=TRUE)

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

immune_type <- getarg (c("i1","i2","i3"), "i1")
ltimmune_type <- getarg (c("I1","I2"), "I1")
seffect_type <- getarg (c("e1","e2","e3"), "e1")
season_type <- getarg (c("s1","s2"), "s1")

het_virus <- getarg ("het") == "het"

run_sims <- getarg ("nosim") != "nosim"

# The default values below were found using the "search" and
# "search-lt" scripts, with proxyDss-filter proxies.

imm_decay <- 
  ( if (ltimmune_type=="I2") c(NL63=0.90,E229=0.98,OC43=0.96,HKU1=0.80)
    else c(NL63=0.9175,E229=0.9850,OC43=0.9500,HKU1=0.9750) )

ltfac <- 0.00225

if (any (substr(args,1,9) == "ltfactor:"))
{ stopifnot (sum (substr(args,1,9) == "ltfactor:") == 1)
  ltfac_arg <- substr (args [substr(args,1,9) == "ltfactor:"], 10, 100)
  ltfac <- eval(parse(text=ltfac_arg))
  args <- args [substr(args,1,9) != "ltfactor:"]
}

if (any (substr(args,1,6) == "decay:"))
{ stopifnot (sum (substr(args,1,6) == "decay:") == 1)
  decay_arg <- substr (args [substr(args,1,6) == "decay:"], 7, 100)
  imm_decay <- eval(parse(text=paste0("c(",decay_arg,")")))
  args <- args [substr(args,1,6) != "decay:"]
}

stage_decay1 <- 0.94
stage_decay2 <- 0.94

ltimm_decay <- c(NL63=0.98,E229=0.98,OC43=0.98,HKU1=0.98)

if (any (substr(args,1,8) == "ltdecay:"))
{ stopifnot (sum (substr(args,1,8) == "ltdecay:") == 1)
  ltdecay_arg <- substr (args [substr(args,1,8) == "ltdecay:"], 9, 100)
  ltimm_decay <- eval(parse(text=paste0("c(",ltdecay_arg,")")))
  args <- args [substr(args,1,8) != "ltdecay:"]
}

R_estimates <- args

if (length(R_estimates)==0) 
{ R_estimates <- "proxyW"  # As in the Kissler, et al paper
}

stopifnot(length(R_estimates)==1)


# PLOT SETUP.

pdf (paste0 ("R-model-",R_estimates,"-",R_est_type,
             if (season_type=="s2") "-s2",
             if (immune_type=="i2") "-i2" else if (immune_type=="i3") "-i3",
             if (ltimmune_type=="I2") "-I2",
             if (seffect_type=="e2") "-e2" else if (seffect_type=="e3") "-e3",
             if (het_virus) "-het",
             ".pdf"),
     height=8,width=6)
par(mar=c(1.5,2.3,3,0.5),mgp=c(1.4,0.3,0),tcl=-0.22)
yrcols <- c("red","green","blue","orange","darkcyan","darkmagenta")


# READ THE FILE OF CORONAVIRUS R ESTIMATES.

R_est <- read.csv (paste0("../R-est/R-est-",R_estimates,".csv"),
                   header=TRUE, stringsAsFactors=FALSE)

R_est$start <- as.Date(R_est$start)
R_est$yrcont <- (0:(nrow(R_est)-1)) / (365.24/7)


# INCLUDE UTILITY FUNCTIONS.  They need "start" and "week" to be defined.

start <- R_est$start  # Start dates of weeks being analysed
year  <- R_est$year   # Year for each week
week  <- R_est$week   # Number of each week in its year

source("../util/util.R")


# ----- DO EVERYTHING FOR BOTH ALPHACORONAVIRUSES AND BETACORONAVIRUSES -----

for (virus_group in list (alphacoronaviruses, betacoronaviruses)) {


# CREATE A DATA FRAME WITH ONLY THE SELECTED R ESTIMATES, FOR FLU SEASON ONLY.
# The new data frame, select_df, has columns start, year, week, season (year
# of start of flu season), season_week (week in the flu season, from 1),
# <virus>_R (for R estimates), <virus>_proxy (for incidence proxy), and
# yrcont (year from 0, varying continuously).

# Default, as in Kissler, et al.

start_season <- 40 # Week of start of "flu season"
end_season <- 20   # Week of end of "flu season" (in year following start)
                   # (except ends week earlier in 2015, since 2014 has 53 weeks)

if (season_type=="s2") # Alternative
{ start_season <- 28     # Week of start of "flu season"
  end_season <- 26       # Week of end of "flu season" (in year following start)
                         # (except ends week earler in 2015, 2014 has 53 weeks)
}

R_est_names <- paste0(virus_group,"_",R_est_type)
proxy_names <- paste0(virus_group,"_proxy")

select_df <- R_est [, c("start", "year", "week", "yrcont",
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
# set of columns are added named <virus>_cumexp that are
# exponentially-weighted sums from the start of the record (not paying
# attention to "seasons"), using the imm_decay values for each virus,
# wrapping around once so that the start has a sum obtained from the
# average of last four starts-of-years.  Similary, columns named
# <virus>_cumexplt that are exponentially-weights sums using
# ltimm_decay values are added. Two more columns called <virus>_stage1
# and <virus>_stage2 are added for use with the "i3" immunity model.

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
  title (paste ("Cumulative seasonal incidence for",virus))

  p <- R_est[,paste0(virus,"_proxy")]
  c <- rep(0,length(p))
  clt <- rep(0,length(p))
  t <- 0
  tlt <- 0
  d <- imm_decay[virus]
  dlt <- ltimm_decay[virus]
  for (i in rep(seq_along(p),2))
  { t <- d*t + p[i]
    tlt <- dlt*tlt + p[i]
    c[i] <- t
    clt[i] <- tlt
    if (i==length(p))
    { t <- mean(c[i-52*(0:3)])
      tlt <- mean(clt[i-52*(0:3)])
    }
  }
  select_df[,paste0(virus,"_cumexp")] <- c[in_season]
  select_df[,paste0(virus,"_cumexplt")] <- clt[in_season]
  plot (start, c, pch=20, col=ifelse(in_season,"black","gray"),
        ylim=c(0,max(c,clt)), 
        ylab="cum incidence, black short term, green long term")
  points (start, clt, pch=20, col=ifelse(in_season,"green","lightgreen"))
  abline(h=0)
  title (paste ("Exponentially-weighted cumulative incidences for",virus))

  p <- R_est[,paste0(virus,"_proxy")]
  c1 <- c2 <- rep(0,length(p))
  t1 <- t2 <- 0
  for (i in rep(seq_along(p),2))
  { t2 <- stage_decay2*t2 + (1-stage_decay1)*t1
    c2[i] <- t2
    t1 <- stage_decay1*t1 + p[i]
    c1[i] <- t1
    if (i==length(p))
    { t1 <- mean(c1[i-52*(0:3)])
      t2 <- mean(c2[i-52*(0:3)])
    }
  }
  select_df[,paste0(virus,"_stage1")] <- c1[in_season]
  select_df[,paste0(virus,"_stage2")] <- c2[in_season]
  plot (rep(start,2), c(c1,c2), pch=rep(c(19,21),each=length(start)), 
        cex=0.5, ylim=c(0,max(c1+c2)), ylab="", 
        col=rep(ifelse(in_season,"black","gray"),2))
  lines (start, c1+c2)
  lines (start, c, col="gray")
  abline(h=0)
  title (paste ("Occupancy of two stages for", virus))
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
  e2 = formula <- paste(formula,"+ sin(2*pi*yrcont) + cos(2*pi*yrcont)"),
  e3 = formula <- paste(formula,"+ sin(1*2*pi*yrcont) + cos(1*2*pi*yrcont)",
                                "+ sin(2*2*pi*yrcont) + cos(2*2*pi*yrcont)",
                                "+ sin(3*2*pi*yrcont) + cos(3*2*pi*yrcont)",
                                "+ sin(4*2*pi*yrcont) + cos(4*2*pi*yrcont)",
                                "+ sin(5*2*pi*yrcont) + cos(5*2*pi*yrcont)",
                                "+ sin(6*2*pi*yrcont) + cos(6*2*pi*yrcont)"
                  )
)

model_df <- NULL
for (virus in virus_group)
{ v_df <- select_df
  row.names(v_df) <- NULL
  v_df$virus <- virus
  if (immune_type=="i3")
  { v_df [, paste0(virus,"_same")] <- 
      v_df [, paste0(virus,"_stage1")] +
      v_df [, paste0(virus,"_stage2")]
    v_df [, paste0(virus,"_other")] <- 
      v_df [, paste0 (other_virus_of_type[virus], "_stage1")] +
      v_df [, paste0 (other_virus_of_type[virus], "_stage2")]
  }
  else
  { v_df [, paste0(virus,"_same")] <- 
      v_df [, paste0 (virus, 
                      if (immune_type=="i1") "_cum" else "_cumexp")]
    v_df [, paste0(virus,"_other")] <- 
      v_df [, paste0 (other_virus_of_type[virus], 
                      if (immune_type=="i1") "_cum" else "_cumexp")]
  }
  v_df[,paste0(other_virus_of_type[virus],"_same")] <- 0
  v_df[,paste0(other_virus_of_type[virus],"_other")] <- 0
  formula <- paste (formula, "+", paste0(virus,"_same"))
  formula <- paste (formula, "+", paste0(virus,"_other"))
  if (ltimmune_type=="I2")
  { v_df [, paste0(virus,"_samelt")] <- 
      log (1 - ltfac * v_df [, paste0 (virus, "_cumexplt")])
    v_df [, paste0(virus,"_otherlt")] <-  
      log (1 - ltfac * v_df [, paste0 (other_virus_of_type[virus],"_cumexplt")])
    v_df[,paste0(other_virus_of_type[virus],"_samelt")] <- 0
    v_df[,paste0(other_virus_of_type[virus],"_otherlt")] <- 0
    formula <- paste (formula, "+", paste0(virus,"_samelt"))
    formula <- paste (formula, "+", paste0(virus,"_otherlt"))
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

trend_spline <- 
  bs (model_df$yrcont, Boundary=range(model_df$yrcont)+c(-3.5,3.5)/365.24,
      knots = c ((2/3)*min(model_df$yrcont) + (1/3)*max(model_df$yrcont),
                 (1/3)*min(model_df$yrcont) + (2/3)*max(model_df$yrcont)))

# Fit model, perhaps repeatedly, each time adjusting weights based on previous 
# fit to account for heterskedasticity w.r.t. virus.
  
cat ("\nMODEL FOR", paste(virus_group,collapse=" & "), "\n\n")
options(digits=9)

var_ratio_2over1 <- 1

for (rpt in if (het_virus) 1:4 else 1)
{
  model <- lm (parse(text=formula)[[1]], data=model_df, x=TRUE, y=TRUE,
               weights=c(rep(c(var_ratio_2over1,1),each=sum(in_season))))

  print(summary(model))
  
  resid <- log(R_value) - as.vector(predict(model,model_df))
  
  virus_residuals <- 
    list (resid[1:sum(in_season)], resid[(sum(in_season)+1):(2*sum(in_season))])
  names(virus_residuals) <- virus_group
  
  for (virus in virus_group)
  { cat ("Residual standard deviation for", virus, ":",
          round (sd(virus_residuals[[virus]],na.rm=TRUE), 5), "\n")
  }

  var_ratio_2over1 <- mean(virus_residuals[[2]]^2,na.rm=TRUE) /
                      mean(virus_residuals[[1]]^2,na.rm=TRUE)
}


# FUNCTIONS TO COMPUTE VALUE OF SEASONAL EFFECT FOR e2 and e3 MODELS.

seffect_e2 <- function (yrcont)
{ mc <- coef(model)
  sin(2*pi*yrcont)*mc[1] + cos(2*pi*yrcont)*mc[2]
}

seffect_e3 <- function (yrcont)
{ mc <- coef(model)
  tn <- ncol(trend_spline)
  ( sin(1*2*pi*yrcont)*mc[tn+1] + cos(1*2*pi*yrcont)*mc[tn+2]
    + sin(2*2*pi*yrcont)*mc[tn+3] + cos(2*2*pi*yrcont)*mc[tn+4]
    + sin(3*2*pi*yrcont)*mc[tn+5] + cos(3*2*pi*yrcont)*mc[tn+6]
    + sin(4*2*pi*yrcont)*mc[tn+7] + cos(4*2*pi*yrcont)*mc[tn+8]
    + sin(5*2*pi*yrcont)*mc[tn+9] + cos(5*2*pi*yrcont)*mc[tn+10]
    + sin(6*2*pi*yrcont)*mc[tn+11] + cos(6*2*pi*yrcont)*mc[tn+12] )
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


# PLOT COMPONENTS OF MODEL FOR EACH VIRUS AND SEASON.  Similar to plots
# of Kissler, et al. Figure 1.

plot_components <- function (s, virus, logarithmic=FALSE)
{
  trans <- if (logarithmic) log else identity
  itrans <- if (logarithmic) identity else exp
  mc <- coef(model)
  this <- model_df$season==s & model_df$virus==virus
  df <- model_df[this,]
  plot (c(1,season_length), trans(c(1,1)), type="n",
        ylim=trans(c(0.5,2.0)), xlim=c(0,season_length), xaxs="i",
        ylab = if (logarithmic) "Additive effect on log(R)"
               else "Multiplicative effect on R")
  points (trans(model_df[,paste0(virus,"_R")][this]), pch=20, col="pink")
  abline (v=1, h=trans(1))
  abline(h = if (logarithmic) c(-0.6,-0.4,-0.2,0.2,0.4,0.6) else c(0.5,1.5,2.0),
         col="gray", lty=3)
  if (season_type=="s2")  # Show Kissler, et al season as dotted lines
  { abline(v=40+1-start_season,col="gray",lty=3)
    abline(v=40+33-start_season,col="gray",lty=3)
  }
  lines (itrans (predict(model,model_df)[this]), col="red", lwd=2)
  if (seffect_type=="e3")
  { tn <- ncol(trend_spline)
    trend_component <- as.vector (trend_spline[this,] %*% mc[1:tn])
    lines (itrans(trend_component), col="green", lwd=2)
  }
  seasonal_component <- switch (seffect_type,
      e1 = as.vector (seasonal_spline[1:season_length,] 
                       %*% mc[1:ncol(seasonal_spline)]),
      e2 = seffect_e2 (df$yrcont),
      e3 = seffect_e3 (df$yrcont)
    )
  mc0 <- mc
  if (season_type=="s1" && seffect_type!="e1")
  { mc0[paste0(virus,"_overall")] <- 
      mc0[paste0(virus,"_overall")] + seasonal_component[1]
    seasonal_component <- seasonal_component - seasonal_component[1]
  }
  lines (itrans(seasonal_component), col="orange", lwd=2)
  same <- df[,paste0(virus,"_same")] * mc[paste0(virus,"_same")]
  lines (itrans(same), col="black", lwd=2)
  other <- df[,paste0(virus,"_other")] * mc[paste0(virus,"_other")]
  lines (itrans(other), col="gray", lwd=2)
  if (ltimmune_type=="I2")
  { samelt <- df[,paste0(virus,"_samelt")] * mc[paste0(virus,"_samelt")]
    lines (itrans(samelt), col="black", lwd=2, lty=2)
    otherlt <- df[,paste0(virus,"_otherlt")] * mc[paste0(virus,"_otherlt")]
    lines (itrans(otherlt), col="gray", lwd=2, lty=2)
  }
  if (immune_type=="i1" && seffect_type=="e1")
  { points (1, itrans (mc[paste0(virus,"_season_",s)]
                     - mc[paste0(virus_group[2],"_season_2014")]), 
               pch=19, col="gray")
    points (1, itrans (mc[paste0(virus,"_season_",s)]), pch=19)
  }
  else
  { points (1, itrans (mc0[paste0(virus,"_overall")]), pch=19)
  }
  title(paste0(virus," ",s,"-",s+1," ",R_estimates))
}

par(mfrow=c(5,2))

for (s in unique(model_df$season))
{ for (virus in virus_group)
  { plot_components(s,virus,logarithmic=FALSE)
  }
}

for (s in unique(model_df$season))
{ for (virus in virus_group)
  { plot_components(s,virus,logarithmic=TRUE)
  }
}

# Again, but bigger, and in a different order.

for (virus in virus_group)
{ par(mfrow=c(3,2))
  for (s in unique(model_df$season))
  { plot_components(s,virus,logarithmic=FALSE)
    plot_components(s,virus,logarithmic=TRUE)
  }
}


# DO SIMULATIONS FOR SUITABLE MODELS.  Done only if the model is for Rt, and
# season type is s2, immune type is i2 (plus maybe I2), and seasonal effect 
# type is e2 or e3. Can be disabled (for speed) with nosim.
#
# The simulation is done multiple times for the five seasons, using the
# trend spline and seasonal effects.  Initial values for the history of
# incidence proxy values are taken from the end of one of the five seasons
# in the preceding simulation (and set somewhat arbitrarily for the first
# simulation).  Some randomness is introduced into the initial values, 
# and in the daily values generated.
#
# Simulations are done day-by-day (not week-by-week).  So that the
# simulated values are comparable to the weekly incidence proxies,
# they are multiplied by 7 when plotted.

if (!run_sims || R_est_type != "Rt" || immune_type != "i2" 
              || seffect_type != "e2" && seffect_type != "e3")
{ next
}

set.seed(1)

nsims_warmup <- 10   # Number of initial simulations to "warm up"
nsims_plotted <- 15  # Number of subsequent simulations that are plotted
nsims <- nsims_warmup + nsims_plotted  # Total number of simulations

mc <- coef(model)

yrcontd <- rep(R_est$yrcont,each=7) + (0:6)/365.24 - 3.5/365.24

proxy1 <- R_est[,paste0(virus_group[1],"_proxy")]
proxy2 <- R_est[,paste0(virus_group[2],"_proxy")]

daily_decay <- imm_decay ^ (1/7)
ltdaily_decay <- ltimm_decay ^ (1/7)

# Generation interval distribution used by Kissler, et al.

SARS_shape <- 2.35
SARS_scale <- 9.48

SARS_gen_interval <-
  dweibull (1:ceiling(qweibull(0.99, shape=SARS_shape, scale=SARS_scale)),
            shape=SARS_shape, scale=SARS_scale)

gen_interval <- SARS_gen_interval # / sum(SARS_gen_interval)
rev_gen_interval <- rev(gen_interval)

# Pre-compute trend and seasonal effects for all days.

tn <- ncol(trend_spline)
tseff <- 
( if (seffect_type=="e2")
    seffect_e2(yrcontd) 
  else 
    seffect_e3(yrcontd) + as.vector(predict (trend_spline,yrcontd) %*% mc[1:tn])
)

# Initial levels for exponentially-decaying past incidence.

q <- c (quantile(proxy1,0.1), quantile(proxy2,0.1))
t <- q / (1-imm_decay[virus_group])
tlt <- q / (1-ltimm_decay[virus_group])
past <- list (rep(q[1],length(gen_interval)), rep(q[2],length(gen_interval)))

# Do the simulations.

p <- numeric(2)

ylim <-  max (proxy1, proxy2)

sim <- rep (list(numeric(7*length(start))), 2)
sims <- list()

sv_past <- NULL

wsave <- sample(rep(1:5,length=nsims))

for (w in 1:nsims) 
{
  # Do one simulation, for all years.

  for (i in 1:(7*length(start)))
  {
    for (j in 1:2)
    { virus <- virus_group[j]
      log_Rt <- tseff[i] + mc[paste0(virus,"_same")] * t[j] +
                           mc[paste0(virus,"_other")] * t [if (j==1) 2 else 1] +
                           mc[paste0(virus,"_overall")]
      if (ltimmune_type=="I2")
      { log_Rt <- log_Rt + mc[paste0(virus,"_samelt")] * log(1-ltfac*tlt[j]) +
                           mc[paste0(virus,"_otherlt")] * 
                             log (1 - ltfac * tlt[if (j==1) 2 else 1])
      }
      inf <- sum (past[[j]]*rev_gen_interval)
      p[j] <- exp (log_Rt + rnorm(1,0,0.1)) * inf    # Some randomness here
      past[[j]] <- c (past[[j]][-1], p[j])
      sim[[j]][i] <- p[j]

      t[j] <- p[j] + t[j]*daily_decay[virus]
      tlt[j] <- p[j] + tlt[j]*ltdaily_decay[virus]
    }

    if ((i+5) %% 365 == 0)
    { wsave[w] <- wsave[w] - 1
      if (wsave[w] == 0)
      { sv_past <- past
        sv_t <- t
        sv_tlt <- tlt
      }
    }
  }

  # Save simulation for later plotting, if past warm-up.

  if (w > nsims_warmup)
  { sims <- c(sims,list(sim))
    ylim <- max (ylim, 7*sim[[1]], 7*sim[[2]], na.rm=TRUE)
  }

  # Set up initial state for next simulation. Randomized a bit.

  n <- exp(rnorm(2,0,0.2))
  for (j in 1:2) 
  { past[[j]] <- sv_past[[j]] * n[j]
  }
  t <- sv_t * n
  tlt <- sv_tlt * n

}

# Plot the observed incidence, then saved simulations.

par(mfrow=c(4,1))

plot (start, rep(0,length(start)),
      ylim=c(0,1.02*ylim), yaxs="i", type="n", ylab="Incidence proxy")

lines (start, proxy1, col="blue")
lines (start, proxy2, col="red")

title (paste 
 ("Observed proxies for",virus_group[1],"(blue) and",virus_group[2],"(red)"))

for (k in seq_along(sims))
{ 
  sim <- sims[[k]]

  plot (start, rep(0,length(start)),
        ylim=c(0,1.02*ylim), yaxs="i", type="n", ylab="Simulated incidence")

  weekly <- filter (sim[[1]],rep(1,7)) [seq(4,length(sim[[1]]),by=7)]
  lines (start, weekly, col="blue")
  weekly <- filter (sim[[2]],rep(1,7)) [seq(4,length(sim[[2]]),by=7)]
  lines (start, weekly, col="red")

  if (k==1) title (paste ("Simulations of",virus_group[1],"and",virus_group[2]))
}


# ----- END OF LOOP OVER THE TWO VIRUS GROUPS -----

}

# ALL DONE.

dev.off()
