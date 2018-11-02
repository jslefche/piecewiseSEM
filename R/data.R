#' Data set from Shipley
#'
#' @format A \code{data.frame} with 1900 observations of 9 variables.
#' \describe{
#' \item{site}{Site of observation}
#' \item{tree}{Tree of observation}
#' \item{lat}{Latitude}
#' \item{year}{Year of observation}
#' \item{Date}{Julian date of first bud burst}
#' \item{DD}{Cumulative degree days until first bud burst}
#' \item{Growth}{Increase in stem diameter}
#' \item{Survival}{Plant species richness}
#' \item{Live}{Alive (1) or dead (0)}
#' }
#' @name shipley
#' @docType data
#' @keywords data
"shipley"

#' Data set from Keeley 
#' 
#' @format A \code{data.frame} with 90 observations of 8 variables.
#' \describe{
#' \item{distance}{Distance to coast}
#' \item{elev}{Elevation from sea level}
#' \item{abiotic}{Abiotic favorability}
#' \item{age}{Age of stand before fire}
#' \item{hetero}{Plot heterogeneity}
#' \item{firesev}{Severity of fire}
#' \item{cover}{Cover of plants}
#' \item{rich}{Plant species richness}
#' }
#' @name keeley
#' @docType data
#' @keywords data
"keeley"

#' Data set from Jutila
#'
#' @format A \code{data.frame} with 354 observations of 4 variables.
#' \describe{
#' \item{grazed}{Whether meadows were exposed to grazing: 0 = no, 1 = yes}
#' \item{mass}{Plant biomass in g m[-2]}
#' \item{elev}{Elevation of the plot above mean sea level}
#' \item{rich}{Plant species richness per m[2]}
#' }
#' @name meadows
#' @docType data
#' @keywords data
"meadows"