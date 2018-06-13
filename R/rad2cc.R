##     rad2cc.R Radiation to cloud cover
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


#' @title Cloud cover estimation from radiation data
#' @description Cloud cover estimation from radiation data
#' @param rsds Shortwave downwelling radiation. Ignored if \code{rtds} is provided
#' @param rlds Longawe downwelling radiation. Ignored if \code{rtds} is provided.
#' @param rtds Total (longwave + shortwave) radiation. If provided, \code{rsds} and \code{rlds} are ignored.
#' @return A climate4R grid of estimated cloud cover
#' @author J. Bedia, S. Herrera
#' @export
#' @importFrom transformeR isGrid redim checkSeason checkDim getGridUnits getRefDates subsetGrid mat2Dto3Darray array3Dto2Dmat bindGrid gridArithmetics
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom udunits2 ud.are.convertible
#' @examples
#' data("rsds.iberia")
#' data("rlds.iberia")
#' cc <- rad2cc(rsds.iberia, rlds.iberia)
#' library(visualizeR)
#' spatialPlot(climatology(cc), rev.colors = TRUE, backdrop.theme = "coastline")
#' # Total radiation can be alternatively used:
#' library(transformeR)
#' rtds <- gridArithmetics(rsds.iberia, rlds.iberia, operator = "+")
#' cc2 <- rad2cc(rtds = rtds)
#' identical(cc2, cc)

rad2cc <- function(rsds = NULL, rlds = NULL, rtds = NULL) {
    # Consistency checks:
    if (is.null(rtds)) {
        stopifnot(!is.null(rsds) & !is.null(rlds))
        stopifnot(isGrid(rsds))
        stopifnot(isGrid(rlds))
        # Redim to have members:
        rsds  %<>%  redim(member = TRUE, var = FALSE)
        rlds %<>%  redim(member = TRUE, var = FALSE)
        # Check dim
        suppressMessages(checkDim(rsds, rlds))
        # Check season
        checkSeason(rsds, rlds)
        # Check units
        u1 <- getGridUnits(rsds)
        u2 <- getGridUnits(rlds)
        if (u1 != "W.m-2") {
            if (!ud.are.convertible(u1, "W.m-2")) {
                stop("Non compliant rsds units (should be convertible to \'W.m-2\')")
            }
            rsds <- udConvertGrid(rsds, new.units = "W.m-2") %>% redim(member = TRUE)
        }
        if (u2 != "W.m-2") {
            if (!ud.are.convertible(u2, "W.m-2")) {
                stop("Non compliant rlds units (should be convertible to \'W.m-2\')")
            }
            rlds <- udConvertGrid(rlds, new.units = "W.m-2") %>% redim(member = TRUE)
        }
        if (!is.null(rtds)) message("NOTE: rtds argument will be ignored, and calculated from rlds and rsds provided")
        rtds <- gridArithmetics(rlds, rsds, operator = "+")
        rsds <- rlds <- NULL
    } else {
        u1 <- getGridUnits(rtds)
        if (u1 != "W.m-2") {
            if (!ud.are.convertible(u1, "W.m-2")) {
                stop("Non compliant rtds units (should be convertible to \'W.m-2\')")
            }
            rtds <- udConvertGrid(rtds, new.units = "W.m-2") %>% redim(member = TRUE)
        }
    }
    rtds %<>% redim(member = TRUE)
    jday <- getRefDates(rtds) %>% as.Date() %>% format("%j") %>% as.integer()
    coords <- getCoordinates(rtds)
    ref.coords <- expand.grid(coords$y, coords$x)[2:1]
    delta <- -23.44 * cos((360 / 365) * (jday + 10) * 2 * pi / 360) %>% matrix(nrow = length(jday), ncol = nrow(ref.coords))
    alpha <- matrix(90 + ref.coords[ ,2], ncol = nrow(ref.coords), nrow = length(jday), byrow = TRUE) - delta
    R0 <- 990 * sin(alpha * 2 * pi / 360) - 30
    n.mem <- getShape(rtds, "member")
    l <- lapply(1:n.mem, function(x) {
        tot.rad <- subsetGrid(rtds, member = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        a <- exp(log((R0 - tot.rad) / (.75 * R0)) / 3.4)
        a[which(R0 - tot.rad < 0)] <- 0
        rtds$Data <- mat2Dto3Darray(a, x = coords$x, y = coords$y)
        return(rtds)
    })
    cc <- suppressWarnings(bindGrid(l, dimension = "member"))
    cc$Variable$varName <- "clt"
    cc$Variable$level <- NULL
    attr(cc$Variable, "units") <- "1"
    attr(cc$Variable, "longname") <- "Cloud_area_fraction"
    attr(cc$Variable, "description") <- "Estimated cloud area fraction from radiation"
    invisible(cc)
}
