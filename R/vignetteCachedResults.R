# Copyright 2015-2019 Venelin Mitov
#
# This file is part of POUMM.
#
# POUMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# POUMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with POUMM  If not, see <http://www.gnu.org/licenses/>.


#' Cached objects for the POUMM vignettes and examples
#'
#' A list containing a simulated tree, trait-values and POUMM objects (model 
#' fits). To use these objects in examples you can load them into the global 
#' workspace with the command: 
#' `data(vignetteCachedResults); list2env(vignetteCachedResults, globalenv());`.
#'
#' @format This is a list containing the following named elements:
#' \describe{
#'   \item{g, z, e }{ numeric vectors of simulated genotypic values, phenotypic values and measurement errors. }
#'   \item{tree }{ a simulated phylogenetic tree. }
#'   \item{fitPOUMM, fitPOUMM2, fitH2tMean }{ POUMM fit objects to tree and z. }
#' }
"vignetteCachedResults"
