##' Internal: helper that returns an index sequence with the last value repeated.
##' @param min lower index.
##' @param max upper index.
##' @keywords internal
##' @noRd
setSeq <- function(min, max){
  if(min == max){
    ret <- 1
  }else{
    ret <- c(1:(max - min), max - min)
  }
  ret
}

##' Internal: sequence helper used in \code{defcon}.
##' @param x vector whose length determines the sequence.
##' @keywords internal
##' @noRd
setS <- function(x){
  setSeq(1, length(x))
}

##' Setup a default minimal samjr configuration
##'
##' Generates a basic configuration list from a data object. The dimensions
##' (years, ages, fleets) are inferred from the data. The configuration is
##' intentionally simplistic and is intended as a starting point that the user
##' will modify before calling \code{\link{sam.fit}}.
##'
##' The list contains, among others:
##' \describe{
##' \item{\code{minAge}, \code{maxAge}, \code{maxAgePlusGroup}}{age range and
##' plus-group flag per fleet.}
##' \item{\code{keyLogFsta}}{matrix of integer indices (\code{-1} = no F)
##' linking each (fleet, age) cell to a fishing-mortality state.}
##' \item{\code{corFlag}}{integer vector controlling F-correlation:
##' 0 independent, 1 compound symmetry, 2 AR(1)-like decay. samjr reads only
##' the first element.}
##' \item{\code{keyLogFpar}}{matrix of integer indices linking survey ages
##' to catchability parameters.}
##' \item{\code{keyVarF}, \code{keyVarLogN}, \code{keyVarObs}}{indexing
##' matrices for the F process, N process, and observation variances.}
##' \item{\code{obsCorStruct}}{factor with levels \code{ID}, \code{AR},
##' \code{US} controlling observation cross-age correlation per fleet.}
##' \item{\code{keyCorObs}}{matrix of integer indices for AR or IGAR
##' correlation parameters.}
##' \item{\code{stockRecruitmentModelCode}}{0 = random walk, 1 = Ricker,
##' 2 = Beverton-Holt.}
##' \item{\code{noScaledYears}, \code{keyScaledYears}, \code{keyParScaledYA}}{
##' optional catch-scaling configuration; see \code{\link{sam.fit}}.}
##' \item{\code{fbarRange}}{integer vector of length two giving the age range
##' used to compute average fishing mortality (\eqn{\bar F}).}
##' \item{\code{usePOP}, \code{useHSP}}{switch the close-kin
##' parent-offspring and half-sibling likelihood contributions on (1) or off
##' (0). Both default to 0 and require a \code{ckmr} pair table on the data
##' object; see \code{\link{ckmrData}}.}
##' \item{\code{ckmrPsi}}{exponent in the per-capita reproductive output
##' \eqn{fec = MO \cdot SW^{\psi}}. Fixed at this value unless
##' \code{ckmrEstimatePsi} is 1, in which case it is the starting value.}
##' \item{\code{ckmrEstimatePsi}}{if 1, estimate \eqn{\psi} as
##' \eqn{\exp(\texttt{logPsim1}) + 1} (so \eqn{\psi > 1}) instead of holding
##' it fixed. Only has an effect when close-kin data are in use. Default 0.}
##' \item{\code{ckmrPlusExtra}}{how many years beyond \code{maxAge} a
##' plus-group animal's true age may reach when the half-sibling term
##' marginalises over its unknown birth year. Default 10.}
##' \item{\code{ckmrPlusNodes}}{how many of those candidate birth years are
##' resolved individually; the rest are lumped into one node whose weight is
##' still summed exactly. Default 6, which on a typical plus group carries
##' about 99 percent of the weight in the resolved nodes.}
##' \item{\code{ckmrScale}}{multiplier converting model numbers-at-age to
##' individuals. Close-kin probabilities are inversely proportional to absolute
##' abundance, so this must match the units of the catch data (e.g. 1000 when
##' the catch is in thousands).}
##' }
##'
##' @param dat data list as returned by \code{\link{setup.sam.data}}.
##' @param level either 1 (basic, independent observations) or 2 (with AR
##' correlation on survey fleets).
##' @return a list of configuration entries.
##' @export
defcon <- function(dat, level = 1){
  fleetTypes <- dat$fleetTypes
  ages <- cbind(dat$minAgePerFleet, dat$maxAgePerFleet)
  ages[fleetTypes %in% c(3, 5, 6), ] <- NA
  minAge <- min(ages, na.rm = TRUE)
  maxAge <- max(ages, na.rm = TRUE)
  ages[is.na(ages)] <- minAge
  nAges <- maxAge - minAge + 1
  nFleets <- nrow(ages)
  ret <- list()
  ret$minAge <- minAge
  ret$maxAge <- maxAge
  ret$maxAgePlusGroup <- as.integer(ages[, 2] == max(ages[, 2], na.rm = TRUE))
  x <- matrix(0, nrow = nFleets, ncol = nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i] == 0){
      aa <- ages[i, 1]:ages[i, 2]
      aa <- aa[tapply(dat$logobs[dat$aux[, 2] == i],
                      INDEX = dat$aux[, 3][dat$aux[, 2] == i],
                      function(x) !all(is.na(x)))]
      x[i, aa - minAge + 1] <- setS(aa) + lastMax
      lastMax <- max(x)
    }
  }
  ret$keyLogFsta <- x - 1
  ret$corFlag <- rep(2, length(which(fleetTypes == 0)))

  x <- matrix(0, nrow = nFleets, ncol = nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i] %in% c(1, 2, 3, 4, 6)){
      x[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge + 1)] <-
        setSeq(ages[i, 1], ages[i, 2]) + lastMax
      lastMax <- max(x)
    }
  }
  ret$keyLogFpar <- x - 1
  ret$keyQpow <- matrix(-1, nrow = nFleets, ncol = nAges)

  x <- matrix(0, nrow = nFleets, ncol = nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i] == 0){
      x[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge + 1)] <- lastMax + 1
      lastMax <- max(x)
    }
  }
  ret$keyVarF <- x - 1
  ret$keyVarLogN <- c(1, rep(2, nAges - 1)) - 1

  x <- matrix(0, nrow = nFleets, ncol = nAges)
  lastMax <- 0
  for(i in 1:nrow(x)){
    if(fleetTypes[i] %in% c(0, 1, 2, 3, 4, 6)){
      x[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge + 1)] <- lastMax + 1
      lastMax <- max(x)
    }
  }
  ret$keyVarObs <- x - 1

  ret$obsCorStruct <- factor(rep("ID", nFleets), levels = c("ID", "AR", "US"))
  if(level == 2){
    ret$obsCorStruct[fleetTypes %in% c(2, 4)] <- "AR"
  }
  ret$keyCorObs <- matrix(-1, nrow = nFleets, ncol = nAges - 1)
  colnames(ret$keyCorObs) <- paste(minAge:(maxAge - 1), (minAge + 1):maxAge, sep = "-")
  nextpar <- 0
  for(i in 1:nrow(x)){
    if(ages[i, 1] < ages[i, 2]){
      if((level == 2) & (fleetTypes[i] %in% c(2, 4))){
        ret$keyCorObs[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge)] <- nextpar
        nextpar <- nextpar + 1
      }else{
        ret$keyCorObs[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge)] <- NA
      }
    }
  }

  ret$stockRecruitmentModelCode <- 0
  ret$noScaledYears <- 0
  ret$keyScaledYears <- numeric(0)
  ret$keyParScaledYA <- array(0, c(0, 0))
  ret$fbarRange <- c(ret$minAge, ret$maxAge)
  ret$stockWeightModel <- 0
  ret$keyStockWeightMean <- rep(NA_integer_, nAges)
  ret$keyStockWeightObsVar <- rep(NA_integer_, nAges)
  ret$catchWeightModel <- 0
  ret$keyCatchWeightMean   <- matrix(NA_integer_, nrow = sum(fleetTypes == 0), ncol = nAges)
  ret$keyCatchWeightObsVar <- matrix(NA_integer_, nrow = sum(fleetTypes == 0), ncol = nAges)
  ret$matureModel    <- 0
  ret$keyMatureMean  <- rep(NA_integer_, nAges)
  ret$mortalityModel <- 0
  ret$keyMortalityMean   <- rep(NA_integer_, nAges)
  ret$keyMortalityObsVar <- rep(NA_integer_, nAges)
  ret$predVarObsLink <- matrix(NA_integer_, nrow = nFleets, ncol = nAges)
  for(i in seq_len(nFleets)){
    if(ages[i, 1] < ages[i, 2]){
      ret$predVarObsLink[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge + 1)] <- -1L
    }
  }
  ret$keyXtraSd <- matrix(NA_integer_, nrow = 0, ncol = 4)
  ret$fixVarToWeight <- rep(0L, nFleets)
  ret$usePOP <- 0
  ret$useHSP <- 0
  ret$ckmrPsi <- 1.5
  ret$ckmrEstimatePsi <- 0
  ret$ckmrScale <- 1
  ret$ckmrPlusExtra <- 10
  ret$ckmrPlusNodes <- 6
  ret
}

##' Save a samjr configuration to a SAM-style \code{.cfg} text file
##'
##' Writes a configuration list (as returned by \code{\link{defcon}} or
##' \code{\link{loadConf}}) to a text file in SAM's \code{$name} block
##' format, suitable for round-tripping back through \code{loadConf}.
##'
##' @param x a configuration list.
##' @param file output file path (or \code{""} for stdout).
##' @param overwrite if \code{FALSE} (default), refuse to overwrite an
##' existing file.
##' @return invisibly \code{NULL}.
##' @export
saveConf <- function(x, file = "", overwrite = FALSE){
  writeOne <- function(v, ...){
    if(is.factor(v)){
      cat(" | Possible values are:", paste0('\"', levels(v), '\"'), ...)
      cat("\n", paste0('\"', v, '\"'), "\n", ...)
    }else if(is.matrix(v)){
      if(nrow(v) > 0){
        cat(capture.output(prmatrix(v, rowlab = rep("", nrow(v)),
                                     collab = rep("   ", ncol(v)))),
            sep = "\n", ...)
      }else cat("\n", ...)
    }else{
      cat("\n", v, "\n", ...)
    }
  }
  if(file != "" && file.exists(file) && !overwrite){
    cat("Notice: Did not overwrite existing file\n")
    return(invisible(NULL))
  }
  app <- function() if(file != "") TRUE else FALSE
  cat(paste0("# Configuration saved: ", date()), file = file)
  cat("\n#\n# Where a matrix is specified rows correspond to fleets and columns to ages.\n",
      file = file, append = TRUE)
  cat("# Same number indicates same parameter used\n",
      file = file, append = TRUE)
  cat("# Numbers (integers) start from zero and must be consecutive\n",
      file = file, append = TRUE)
  cat("# Negative numbers indicate that the parameter is not included in the model\n#",
      file = file, append = TRUE)
  for(nm in names(x)){
    cat("\n$", file = file, append = TRUE, sep = "")
    cat(nm, file = file, append = TRUE)
    cat("\n#", file = file, append = TRUE)
    writeOne(x[[nm]], file = file, append = TRUE)
  }
  invisible(NULL)
}

##' Load a samjr configuration from a SAM-style \code{.cfg} text file
##'
##' Reads a configuration file in the format produced by SAM's
##' \code{saveConf} (one block per key, introduced by \code{$name},
##' values on the following lines, \code{#} for comments). The expected
##' types come from the corresponding entries in \code{\link{defcon}(dat)},
##' so the parser knows whether to read a numeric scalar, an integer
##' vector, a matrix or a factor.
##'
##' @param dat data list from \code{\link{setup.sam.data}} (used to derive
##' the expected configuration shape via \code{\link{defcon}}).
##' @param file path to the configuration file.
##' @param patch if \code{TRUE} (default), keys present in
##' \code{defcon(dat)} but missing from the file are filled in with the
##' default values; if \code{FALSE} a missing key is an error.
##' @return a configuration list of the same shape as \code{\link{defcon}(dat)}.
##' @export
loadConf <- function(dat, file, patch = TRUE){
  dconf <- defcon(dat)
  lin <- c(readLines(file), "$end")
  lin <- lin[-grep("^#", lin)]
  keyIdx <- grep("^\\$", lin)
  getIdx <- function(nam){
    idx1 <- grep(paste0("^\\$", nam, "( |$)"), lin) + 1
    idx2 <- min(keyIdx[keyIdx > (idx1 - 1)]) - 1
    if(idx1 <= idx2) idx1:idx2 else NULL
  }
  readByType <- function(template, nam){
    idx <- getIdx(nam)
    if(is.null(idx)) return(template)
    if(is.factor(template)){
      vals <- scan(textConnection(lin[idx]), what = "character", quiet = TRUE)
      return(factor(vals, levels = levels(template)))
    }
    if(is.matrix(template)){
      x <- try(utils::read.table(text = lin[idx], header = FALSE), silent = TRUE)
      if(inherits(x, "try-error")) return(matrix(NA_real_, nrow = 0, ncol = 0))
      return(as.matrix(x))
    }
    scan(textConnection(lin[idx]), quiet = TRUE)
  }
  hasKey <- vapply(names(dconf),
                   function(n) length(grep(paste0("^\\$", n, "( |$)"), lin)) == 1,
                   logical(1))
  if(!all(hasKey) && !patch)
    stop("The configuration file is not compatible with the model version. Use patch=TRUE.")
  conf <- dconf
  for(n in names(dconf)){
    if(hasKey[n]) conf[[n]] <- readByType(dconf[[n]], n)
  }
  conf
}
