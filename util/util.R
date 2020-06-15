# Utility functions / variables.  Some of these functions need "start" to be
# set to the start dates of weeks, and "week" to be set to the number of each
# week within its year.
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


# NAMES AND TYPES OF THE COMMON COLD CORONAVIRUS.  229E is called E229
# so that it will be a valid R identifier.

alphacoronaviruses <- c("NL63","E229")
betacoronaviruses  <- c("OC43","HKU1")

viruses <- c (alphacoronaviruses, betacoronaviruses)

other_virus_of_type <- c (NL63="E229", E229="NL63", OC43="HKU1", HKU1="OC43")


# FUNCTION FOR FLAGGING ANOMALOUS VALUES IN A TIME SERIES.  Large
# positive or negative values might indicate anomalies, since the
# correspond to points that differ from the average of both their
# nearest neighbors and their next-nearest neighbors.

anomalous <- function (x) 
{ f1 <- x - filter(x,c(0.5,0,0.5))
  f2 <- x - filter(x,c(0.5,0,0,0,0.5))
  ifelse (abs(f1) < abs(f2), f1, f2)
}


# FUNCTION TO PLOT DATA VERSUS DAY OF YEAR, STARTING AT JULY 1.

plot_vs_doy <- function (date, data, ylab="", ...)
{
  date <- as.Date(date)
  day <- date - as.Date (paste0(substr(as.character(date),1,4),"-07-01"))
  col <- as.numeric(substr(as.character(date),1,4))
  w <- substr(as.character(date),6,7) < "07"
  day[w] <- 184 + date[w] - 
              as.Date (paste0(substr(as.character(date[w]),1,4),"-01-01"))
  col[w] <- col[w]-1

  w <- date < "2014-07-01"
  day[w] <- date[w] - as.Date("2014-07-01")
  col[w] <- 2014

  plot (day, data, col=yrcols[col-2013], xlab="", ylab=ylab, ...)
  abline (v = 183.5)
  abline (v = 183.5 + cumsum(c(31,28,31,30,31,30)), col="gray")
  abline (v = 183.5 - cumsum(c(31,30,31,30,31,31)), col="gray")
}


# FUNCTION FOR DRAWING LINES AT SIGNIFICANT WEEKS.  
#
# A red line is drawn just before week 40   [ delimiting the supposed
# A geen line is drawn just after week 20     cold/flu season ]
#
# Gray lines are drawn just before the weeks of various holidays (only
# some currently enabled).
#
# Dashed / dotted gray lines may be drawn just before the weeks of Easter
# and the First Sedar of Passover (currently disabled).
#
# Blue lines may be drawn before weeks of major hurricanes (currently
# disabled).

New_Year <- 
  c("2015-01-01", "2016-01-01", "2017-01-01", "2018-01-01", "2019-01-01")

MLK_Day <-  # Third Monday of January
  c("2015-01-19", "2016-01-18", "2017-01-16", "2018-01-15", "2019-01-21")

Easter <- 
  c("2015-04-05","2016-03-27","2017-04-16","2018-04-01","2019-04-21")

Passover <- # First Sedar
  c("2015-04-03","2016-04-22","2017-04-10","2018-03-30","2019-04-19")

July4th <- 
  c("2014-07-04","2015-07-04","2016-07-04","2017-07-04","2018-07-04")

Labor_Day <-  # First Monday of September
  c("2014-09-01","2015-09-07","2016-09-05","2017-09-04","2018-09-03")

Thanksgiving <-  # Fourth Thursday of November
  c("2014-11-27","2015-11-26","2016-11-24","2017-11-23","2018-11-22")

Christmas <-
  c("2014-12-25","2015-12-25","2016-12-25","2017-12-25","2018-12-25")

Hurricanes <-
  c( Matthew="2016-10-03", Harvey="2017-08-24", 
     Irma="2017-09-10",    Michael="2018-10-10")

week_lines <- function ()
{ 
  abline (v=start[week==40]-3.5,col="red")
  abline (v=start[week==20]+3.5,col="green")

  for (d in New_Year)
  { abline (v=start[start<=d & d<=start+6]-3.5,col="gray")
  }
# for (d in MLK_Day)
# { abline (v=start[start<=d & d<=start+6]-3.5,col="gray")
# }
# for (d in Passover)
# { abline (v=start[start<=d & d<=start+6]-3.5,col="gray", lty=3)
# }
# for (d in Easter)
# { abline (v=start[start<=d & d<=start+6]-3.5,col="gray", lty=2)
# }
  for (d in July4th)
  { abline (v=start[start<=d & d<=start+6]-3.5,col="gray")
  }
  for (d in Labor_Day)
  { abline (v=start[start+1==d]-3.5,col="gray")
  }
  for (d in Thanksgiving)
  { abline (v=start[start+4==d]-3.5,col="gray")
  }
  for (d in Christmas)
  { abline (v=start[start<=d & d<=start+6]-3.5,col="gray")
  }

# for (d in Hurricanes)
# { abline (v=start[start<=d & d<=start+6]-3.5,col="blue")
# }

}


# PLOT TWO PARALLEL SERIES WITH LINE FROM POINT OF FIRST TO POINT OF SECOND.

plot_two_with_lines <- function (x,y1,y2,...)
{
  plot (c(x,x), c(y1,y2), type="n", ...)
  points (x, y1, pch=19, col="lightpink")
  for (i in seq_along(x))
  { lines (c(x[i],x[i]), c(y1[i],y2[i]), col="black")
  }
}



# INDICATORS FOR USE IN MODELING.

New_Year_indicator <- rep(0,length(start))
for (d in New_Year)
{ New_Year_indicator [start <= d & d < start+7] <- 1
}

MLK_Day_indicator <- rep(0,length(start))
for (d in MLK_Day)
{ MLK_Day_indicator [start <= d & d < start+7] <- 1
}

Easter_indicator <- rep(0,length(start))
for (d in Easter)
{ Easter_indicator [start <= d & d < start+7] <- 1
}

Passover_indicator <- rep(0,length(start))
for (d in Easter)
{ Passover_indicator [start <= d & d < start+7] <- 1
}

July4th_indicator <- rep(0,length(start))
for (d in July4th)
{ July4th_indicator [start <= d & d < start+7] <- 1
}

Labor_Day_indicator  <- as.numeric(as.character((start+1)) %in% Labor_Day)

Thanksgiving_indicator <- as.numeric(as.character((start+4)) %in% Thanksgiving)

Christmas_indicator <- rep(0,length(start))
for (d in Christmas)
{ 
  # Flag week of Christmas.
  Christmas_indicator [start <= d & d < start+7] <- 1

  # Maybe at slightly lower level if Christmas at end of week?
  # Christmas_indicator [start+6 == d] <- 0.9  # Christmas on Saturday
  # Christmas_indicator [start+5 == d] <- 0.8  # Christmas on Friday

  # Flag week before Christmas at lower level if Christmas is early in week
  Christmas_indicator [start+7 == d] <- 0.6  # Christmas on Sunday
  Christmas_indicator [start+8 == d] <- 0.5  # Christmas on Monday
  Christmas_indicator [start+9 == d] <- 0.4  # Christmas on Tuesday
  Christmas_indicator [start+10 == d] <- 0.3  # Christmas on Wednesday
  Christmas_indicator [start+11 == d] <- 0.2  # Christmas on Thursday

  # Note that the week after Christmas is the week of New Year's Day, so
  # nothing is done here for that week.
}
