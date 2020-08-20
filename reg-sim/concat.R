# Concatenate lists.
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

concat <- function (L)
{
  r <- L[[1]]
  for (i in 2..length(L))
  { n <- L[[i]]
    stopifnot (identical (names(n), names(r)))
    for (j along r)
    { r[[j]] <- c(r[[j]],n[[j]])
    }
  }
  for (j along r)
  { dim (r[[j]]) <- c (length(r[[j]])/length(L), length(L))
  }
  r
}
