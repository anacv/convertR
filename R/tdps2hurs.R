##     tdps2hurs.R Relative humidity from Dew-point temperature 
##
##     Copyright (C) 2020 Santander Meteorology Group (http://www.meteo.unican.es)
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


#' @title Relative humidity from Dew-point temperature
#' @description Estimate the Relative humidity from Dew-point temperature and air temperaure
#' @param tdps Dew-point temperature
#' @param tas Near-surface air temperature
#' @param negatives.to.zero Logical flag indicating if estimated negative values should be truncated to zero \%. Default to TRUE.
#' @return A climate4R CDM grid of estimated relative humidity (in \%)
#' @author J. Bedia
#' @template templateUnits
#' @export
#' @template templateRefPressure
#' @import transformeR
#' @references Lawrence, Mark G., 2005: The relationship between relative humidity and the dewpoint temperature in moist air: A simple conversion and applications. Bull. Amer. Meteor. Soc., 86, 225-233. https://dx.doi.org/10.1175/BAMS-86-2-225 
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom udunits2 ud.are.convertible
#' @importFrom utils packageVersion
#' @seealso hurs2tdps, performing the inverse calculation to derive dew-point temperature from relative humidity and observed temperature
#' @note The formula is a valid aproximation for moist air (RH>50\%), but can yield very inaccurate results otherwise, so use it with caution.
#' @family derivation
#' @family humidity

tdps2hurs <- function(tdps, tas, negatives.to.zero = TRUE) {
    stopifnot(is.logical(negatives.to.zero))
    # Consistency checks:
    if (isMultigrid(tdps) | isMultigrid(tas)) stop("Multigrids are not an allowed input")
    stopifnot(isGrid(tdps) | isGrid(tas))
    # Redim to have members:
    tdps %<>% redim(member = TRUE)
    tas %<>% redim(member = TRUE)
    # Check dim
    suppressMessages(checkDim(tdps, tas))
    # Check season
    checkSeason(tdps, tas)
    # Check units
    tdps.u <- getGridUnits(tdps)
    if (tdps.u != "degC") {
        if (!ud.are.convertible(tdps.u, "degC")) {
            stop("Non compliant tdps units (should be convertible to degrees Celsius)")
        }
        message("[", Sys.time(), "] Converting temperature units ...")
        tdps %<>% udConvertGrid(new.units = "degC") %>% redim(member = TRUE)
    }
    tas.u <- getGridUnits(tas)
    if (tas.u != "degC") {
        if (!ud.are.convertible(tas.u, "degC")) {
            stop("Non compliant tas units (should be convertible to degrees Celsius)")
        }
        message("[", Sys.time(), "] Converting temperature units ...")
        tas %<>% udConvertGrid(new.units = "degC") %>% redim(member = TRUE)
    }
    coords <- getCoordinates(tdps)
    n.mem <- getShape(tdps, "member")
    message("[", Sys.time(), "] Deriving relative humidity ...")
    hurs <- tdps
    l <- lapply(1:n.mem, function(x) {
        dp <- subsetGrid(tdps, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        t <- subsetGrid(tas, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        aux <- 100 - 5*(t - dp)
        if (isTRUE(negatives.to.zero)) aux[which(aux < 0)] <- 0
        hurs$Data <- mat2Dto3Darray(aux, coords$x, coords$y)
        return(hurs)
    })
    hurs <- suppressWarnings(bindGrid(l, dimension = "member"))
    tas <- tdps <- NULL
    hurs$Variable$varName <- "hurs"
    hurs$Variable$level <- NULL
    attr(hurs$Variable, "units") <- "%"
    attr(hurs$Variable, "longname") <- "near-surface relative humidity"
    attr(hurs$Variable, "description") <- "Estimated near-surface relative humidity from dew-point temperature and air temperature"
    attr(hurs, "origin") <- paste0("Calculated with R package 'convertR' v", packageVersion("convertR"))
    attr(hurs, "URL") <- "https://github.com/SantanderMetGroup/convertR"
    message("[", Sys.time(), "] Done.")
    invisible(hurs)
}



#' @title Relative humidity from Dew-point temperature
#' @description Estimate the Relative humidity from Dew-point temperature and air temperaure
#' @param hurs relative humidity 
#' @param tas Near-surface air temperature
#' @return A climate4R CDM grid of estimated dew-point temperature (in degC)
#' @author J. Bedia
#' @template templateUnits
#' @export
#' @template templateRefPressure
#' @import transformeR
#' @references Lawrence, Mark G., 2005: The relationship between relative humidity and the dewpoint temperature in moist air: A simple conversion and applications. Bull. Amer. Meteor. Soc., 86, 225-233. https://dx.doi.org/10.1175/BAMS-86-2-225 
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom udunits2 ud.are.convertible
#' @importFrom utils packageVersion
#' @seealso tdps2hurs, performing the inverse calculation to derive relative humidity from dew-point temperature and observed temperature
#' @note The formula is a valid aproximation for moist air (RH>50\%), but can yield very inaccurate results otherwise, so use it with caution.
#' @family derivation
#' @family humidity

hurs2tdps <- function(hurs, tas) {
    # Consistency checks:
    if (isMultigrid(hurs) | isMultigrid(tas)) stop("Multigrids are not an allowed input")
    stopifnot(isGrid(hurs) | isGrid(tas))
    # Redim to have members:
    hurs %<>% redim(member = TRUE)
    tas %<>% redim(member = TRUE)
    # Check dim
    suppressMessages(checkDim(hurs, tas))
    # Check season
    checkSeason(hurs, tas)
    # Check units
    hurs.u <- getGridUnits(hurs)
    if (hurs.u != "%") {
        if (!ud.are.convertible(hurs.u, "%")) {
            stop("Non compliant hurs units (should be convertible to %)")
        }
        message("[", Sys.time(), "] Converting relative humidity units ...")
        hurs %<>% udConvertGrid(new.units = "%") %>% redim(member = TRUE)
    }
    tas.u <- getGridUnits(tas)
    if (tas.u != "degC") {
        if (!ud.are.convertible(tas.u, "degC")) {
            stop("Non compliant tas units (should be convertible to degrees Celsius)")
        }
        message("[", Sys.time(), "] Converting temperature units ...")
        tas %<>% udConvertGrid(new.units = "degC") %>% redim(member = TRUE)
    }
    coords <- getCoordinates(hurs)
    n.mem <- getShape(hurs, "member")
    message("[", Sys.time(), "] Deriving dew-point temperature ...")
    tdps <- hurs
    l <- lapply(1:n.mem, function(x) {
        rh <- subsetGrid(hurs, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        t <- subsetGrid(tas, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        aux <- t - ((100 - rh) / 5)
        # if (isTRUE(negatives.to.zero)) aux[which(aux < 0)] <- 0
        tdps$Data <- mat2Dto3Darray(aux, coords$x, coords$y)
        return(tdps)
    })
    tdps <- suppressWarnings(bindGrid(l, dimension = "member"))
    tas <- hurs <- NULL
    tdps$Variable$varName <- "tdps"
    tdps$Variable$level <- NULL
    attr(tdps$Variable, "units") <- "degC"
    attr(tdps$Variable, "longname") <- "dewpoint temperature"
    attr(tdps$Variable, "description") <- "Estimated dewpoint temperature from relative humidity and air temperature"
    attr(tdps, "origin") <- paste0("Calculated with R package 'convertR' v", packageVersion("convertR"))
    attr(tdps, "URL") <- "https://github.com/SantanderMetGroup/convertR"
    message("[", Sys.time(), "] Done.")
    invisible(tdps)
}


