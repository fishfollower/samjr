## Reference-point string parser.
##
## Mirrors SAM/.../deterministic_referencepoints.R:11-141 regex-for-regex.

.refpointEnum <- c(None = -99L, FixedF = -1L, StatusQuo = 0L, MSY = 1L,
                   MSYRange = 2L, Max = 3L, xdYPR = 4L, xSPR = 5L,
                   xB0 = 6L, MYPYLdiv = 7L, MYPYL = 8L, MDY = 9L,
                   Crash = 10L, Ext = 11L, Lim = 12L)

##' Parse a reference-point specification string
##'
##' @param x single string, e.g. \code{"MSY"}, \code{"0.35SPR"},
##'   \code{"F=0.2"}, \code{"StatusQuo"}, \code{"0.95MSYRange"}.
##' @return list with \code{rpType} (integer enum) and \code{xVal} (numeric
##'   or \code{numeric(0)}).
##' @examples
##' parseRefpoint("MSY")
##' parseRefpoint("0.35SPR")
##' parseRefpoint("F=0.3")
##' parseRefpoint("0.95MSYRange")
##' parseRefpoint("StatusQuo-1")
##' @export
parseRefpoint <- function(x){
  if(length(x) != 1L) stop("parseRefpoint: one string at a time")
  if(grepl("^StatusQuo$", x))
    x <- "StatusQuo-0"
  else if(grepl("^F=\\.$", x))
    x <- "F=0.0"
  else if(grepl("^F=\\.[[:digit:]]+$", x))
    x <- gsub("=\\.", "=0.", x)
  else if(grepl("^F=[[:digit:]]+\\.$", x))
    x <- gsub("\\.$", ".0", x)
  else if(grepl("^MYPYL$", x))
    x <- "1.0MYPYL"
  else if(grepl("^[[:digit:]]+MYPYL$", x))
    x <- gsub("MYPYL", ".0MYPYL", x)
  else if(grepl("^\\.[[:digit:]]+MYPYL$", x))
    x <- gsub("\\.", "0.", x)
  typePatterns <- list(
    None      = "^$",
    FixedF    = "^F=([[:digit:]]+\\.)?([[:digit:]]+)$",
    StatusQuo = "^StatusQuo(-[[:digit:]]+)?$",
    MSY       = "^MSY$",
    MSYRange  = "^0\\.[[:digit:]]+MSYRange$",
    Max       = "^Max$",
    xdYPR     = "^0\\.[[:digit:]]+dYPR$",
    xSPR      = "^0\\.[[:digit:]]+SPR$",
    xB0       = "^0\\.[[:digit:]]+B0$",
    MYPYLdiv  = "^MYPYLdiv$",
    MYPYL     = "^[[:digit:]]+\\.[[:digit:]]+MYPYL$",
    MDY       = "^MDY$",
    Crash     = "^Crash$",
    Ext       = "^Ext$",
    Lim       = "^Lim$")
  xvalPatterns <- list(
    None      = NA,
    FixedF    = c("(F=)(([[:digit:]]+\\.)?([[:digit:]]+))", "\\2"),
    StatusQuo = c("(StatusQuo)(-?)(([[:digit:]]+)?)", "\\3"),
    MSY       = NA,
    MSYRange  = c("(0\\.[[:digit:]]+)(MSYRange)", "\\1"),
    Max       = NA,
    xdYPR     = c("(0\\.[[:digit:]]+)(dYPR)", "\\1"),
    xSPR      = c("(0\\.[[:digit:]]+)(SPR)", "\\1"),
    xB0       = c("(0\\.[[:digit:]]+)(B0)", "\\1"),
    MYPYLdiv  = NA,
    MYPYL     = c("([[:digit:]]+\\.[[:digit:]]+)(MYPYL)", "\\1"),
    MDY       = NA,
    Crash     = NA,
    Ext       = NA,
    Lim       = NA)
  idx <- which(vapply(typePatterns, function(p) grepl(p, x), logical(1)))
  if(length(idx) == 0L)
    stop("Reference-point specification not recognised: ", x)
  xVal <- numeric(0)
  if(!is.na(xvalPatterns[[idx]][1])){
    pp <- xvalPatterns[[idx]]
    xVal <- as.numeric(gsub(pp[1], pp[2], x))
  }
  list(rpType = .refpointEnum[idx], xVal = xVal, tag = names(.refpointEnum)[idx])
}

##' Merge a new parsed RP into a list of already-parsed RPs.
##'
##' If \code{newRp} has the same \code{rpType} as an existing entry, their
##' \code{xVal} vectors are merged (unique + sorted). Otherwise the new
##' entry is appended.
##' @keywords internal
##' @noRd
mergeRefpoint <- function(rps, newRp){
  if(length(rps) == 0L) return(list(newRp))
  same <- vapply(rps, function(r) isTRUE(r$rpType == newRp$rpType), logical(1))
  if(!any(same)) return(c(rps, list(newRp)))
  rps[same] <- lapply(rps[same], function(r){
    r$xVal <- sort(unique(c(r$xVal, newRp$xVal)))
    r
  })
  rps
}

##' Pretty row labels for a parsed RP entry, one per xVal (or one label for
##' types without xVal).
##' @keywords internal
##' @noRd
refpointName <- function(rp){
  nm <- rp$tag
  if(nm == "None")      return(character(0))
  if(nm == "FixedF")    return(paste0("F=", rp$xVal))
  if(nm == "StatusQuo") return(if(length(rp$xVal) == 0 || all(rp$xVal == 0)) "StatusQuo"
                                else paste0("StatusQuo-", rp$xVal))
  if(nm == "MSY")       return("Fmsy")
  if(nm == "Max")       return("Fmax")
  if(nm == "MYPYLdiv")  return("FMYPYLdiv")
  if(nm == "MDY")       return("Fmdy")
  if(nm == "Crash")     return("Fcrash")
  if(nm == "Ext")       return("Fext")
  if(nm == "MSYRange")  return(c(paste0(rp$xVal, "MSY (low)"),
                                  paste0(rp$xVal, "MSY (high)")))
  if(nm == "xdYPR")     return(paste0("F", rp$xVal, "dYPR"))
  if(nm == "xSPR")      return(paste0("F", 100 * rp$xVal, "%SPR"))
  if(nm == "xB0")       return(paste0("F", 100 * rp$xVal, "%B0"))
  if(nm == "MYPYL")     return(paste0("F", rp$xVal, "MYPYL"))
  if(nm == "Lim")       return("Flim")
  stop("refpointName: unknown tag ", nm)
}
