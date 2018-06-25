##     tas2ws.R Saturation pressure
##
##     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
##
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @title Saturation pressure calculation
#' @description Calculates the saturation pressure using near-surface air temperature and surface pressure data
#' @param tas Near-surface air temperature
#' @param ps Surface pressure
#' @return A climate4R grid of saturation pressure (in Pascals)
#' @author J. Bedia, S. Herrera
#' @keywords internal
#' @import transformeR
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom udunits2 ud.are.convertible
#' @template templateUnits
#' @family derivation
#' @examples \dontrun{
#' data("ps.iberia")
#' data("tas.iberia")
#' ws <- tas2ws(tas.iberia, ps.iberia)
#' library(visualizeR)
#' spatialPlot(climatology(ws), backdrop.theme = "coastline")}

tas2ws <- function(tas, ps) {
    # Consistency checks:
    if (isMultigrid(ps) | isMultigrid(tas)) stop("Multigrids are not an allowed input")
    stopifnot(isGrid(ps) | isGrid(tas))
    # Redim to have members:
    ps %<>%  redim(member = TRUE)
    tas  %<>%  redim(member = TRUE)
    # Check dim
    suppressMessages(checkDim(ps, tas))
    # Check season
    checkSeason(ps, tas)
    # Check units
    u.ps <- getGridUnits(ps)
    u.tas <- getGridUnits(tas)
    if (u.ps != "Pa") {
        if (!ud.are.convertible(u.ps, "Pa")) {
            stop("Non compliant ps units (should be convertible to \'Pascals')")
        }
        message("Converting units ...")
        ps <- udConvertGrid(ps, new.units = "Pa") %>% redim(member = TRUE)
    }
    if (u.tas != "K") {
        if (!ud.are.convertible(u.tas, "K")) {
            stop("Non compliant tas units (should be convertible to \'Kelvin\')")
        }
        message("Converting units ...")
        tas <- udConvertGrid(tas, new.units = "K") %>% redim(member = TRUE)
    }
    message("Calculating saturation pressure ...")
    source(file.path(find.package(package = "convertR"), "constants.R"), local = TRUE)
    coords <- getCoordinates(ps)
    n.mem <- getShape(tas, "member")
    l <- lapply(1:n.mem, function(x) {
        a <- subsetGrid(tas, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        b <- subsetGrid(ps, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        iceMask <- which(a < T0)
        waterMask <- which(a >= T0)
        a[iceMask] <- es0 * exp((6293/T0) - (6293/a[iceMask]) - 0.555 * log(abs(a[iceMask]/T0)))
        a[waterMask] <- es0 * exp((6808/T0) - (6808/a[waterMask]) - 5.09 * log(abs(a[waterMask]/T0)))
        a <- (Rd/Rv) * (a/(b - a))
        ps$Data <- mat2Dto3Darray(a, x = coords$x, y = coords$y)
        a <- b <- NULL
        return(ps)
    })
    tas <- NULL
    wt <- suppressWarnings(bindGrid(l, dimension = "member"))
    wt$Variable$varName <- "ws"
    wt$Variable$level <- NULL
    attr(wt$Variable, "units") <- "Pa"
    attr(wt$Variable, "longname") <- "Saturation_pressure"
    attr(wt$Variable, "description") <- "Estimated saturation pressure"
    message("Done")
    invisible(wt)
}
