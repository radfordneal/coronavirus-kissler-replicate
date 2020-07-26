# FUNCTIONS FOR SEASONAL EFFECTS AND TRENDS.
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


# Make the spline modelling the slow, long-term trend.

make_trend_spline <- function (yrs, model_yrs)
{   bs (model_yrs, Boundary=range(yrs)+c(-3.5,3.5)/365.24,
        knots = c ((2/3)*min(yrs) + (1/3)*max(yrs),
                   (1/3)*min(yrs) + (2/3)*max(yrs)))
}


# Compute seasonal effects. These functions refer to the global 'model'
# and 'trend_spline' function.

seffect_e2 <- function (yrs)
{ mc <- coef(model)
  sin(2*pi*yrs)*mc[1] + cos(2*pi*yrs)*mc[2]
}

seffect_e3 <- function (yrs)
{ mc <- coef(model)
  tn <- ncol(trend_spline)
  ( sin(1*2*pi*yrs)*mc[tn+1] + cos(1*2*pi*yrs)*mc[tn+2]
    + sin(2*2*pi*yrs)*mc[tn+3] + cos(2*2*pi*yrs)*mc[tn+4]
    + sin(3*2*pi*yrs)*mc[tn+5] + cos(3*2*pi*yrs)*mc[tn+6]
    + sin(4*2*pi*yrs)*mc[tn+7] + cos(4*2*pi*yrs)*mc[tn+8]
    + sin(5*2*pi*yrs)*mc[tn+9] + cos(5*2*pi*yrs)*mc[tn+10]
    + sin(6*2*pi*yrs)*mc[tn+11] + cos(6*2*pi*yrs)*mc[tn+12] )
}
