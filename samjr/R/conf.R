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
  ret
}

##' Internal: description written above each key by \code{saveConf}.
##' One element per comment line.
##' @keywords internal
##' @noRd
confComments <- function(){
  list(
    minAge = "The minimum age class in the assessment",
    maxAge = "The maximum age class in the assessment",
    maxAgePlusGroup = c(
      "Is last age group considered a plus group for each fleet (1 yes, or 0 no).",
      "If a single value is given it applies to the first fleet, and the other fleets",
      "are treated as plus groups when their last age equals maxAge."),
    keyLogFsta = c(
      "Coupling of the fishing mortality states processes for each age (only the",
      "first row (= fleet) is used).",
      "Sequential numbers indicate that the fishing mortality is estimated individually",
      "for those ages; if the same number is used for two or more ages, F is bound for",
      "those ages (assumed to be the same). Binding fully selected ages will result in a",
      "flat selection pattern for those ages."),
    corFlag = c(
      "Correlation of fishing mortality across ages (0 independent, 1 compound symmetry,",
      "or 2 AR(1)). Only the first value is used.",
      "0: independent means there is no correlation between F across age.",
      "1: compound symmetry means that all ages are equally correlated.",
      "2: AR(1) first order autoregressive - similar ages are more highly correlated than",
      "ages that are further apart, so similar ages have similar F patterns over time.",
      "If the estimated correlation is high, then the F pattern over time for each age",
      "varies in a similar way. E.g. if almost one, then they are parallel (like a",
      "separable model) and if almost zero then they are independent."),
    keyLogFpar = c(
      "Coupling of the survey catchability parameters (normally the first row is",
      "not used, as that is covered by fishing mortality).",
      "For an SSB index fleet only the first column is used."),
    keyQpow = c(
      "Density dependent catchability power parameters (if any).",
      "For coupled ages the index is modelled as Q * N^b (N at the survey time),",
      "where log(b) is estimated. -1 means b = 1."),
    keyVarF = c(
      "Coupling of process variance parameters for log(F)-process (fishing mortality",
      "normally applies to the first (fishing) fleet; therefore only the first row is used)."),
    keyVarLogN = c(
      "Coupling of the recruitment and survival process variance parameters for the",
      "log(N)-process at the different ages.",
      "Note: samjr currently ignores this key and always estimates one variance",
      "parameter for recruitment (first age) and one shared by all older ages."),
    keyVarObs = c(
      "Coupling of the variance parameters for the observations.",
      "First row refers to the coupling of the variance parameters for the catch data",
      "observations by age.",
      "Second and further rows refers to coupling of the variance parameters for the",
      "index data observations by age."),
    obsCorStruct = c(
      "Covariance structure for each fleet (\"ID\" independent, \"AR\" AR(1), or \"US\" for unstructured).",
      "\"AR\" is an irregular AR(1): the correlation between two ages is 0.5^d, where the",
      "distance d is the sum of the estimated steps between neighbouring ages (see keyCorObs)."),
    keyCorObs = c(
      "Coupling of correlation parameters can only be specified if the AR(1) structure is chosen above.",
      "NA's indicate where correlation parameters can be specified (-1 where they cannot).",
      "Columns are pairs of neighbouring ages."),
    stockRecruitmentModelCode = c(
      "Stock recruitment code (0 for plain random walk, 1 for Ricker, and 2 for Beverton-Holt).",
      "Other stock-recruitment functions can be supplied via the logSRfun argument of sam.fit."),
    noScaledYears = "Number of years where catch scaling is applied.",
    keyScaledYears = "A vector of the years where catch scaling is applied.",
    keyParScaledYA = c(
      "A matrix specifying the couplings of scale parameters (nrow = no scaled years,",
      "ncols = no ages)."),
    fbarRange = "First and last age included in the calculation of average fishing mortality (Fbar).",
    stockWeightModel = c(
      "Integer code describing the treatment of stock weights in the model (0 use as known,",
      "1 use as observations to inform stock weight process (GMRF with cohort and within",
      "year correlations), 2 as 1 with extra correlation within the plus group)."),
    keyStockWeightMean = "Coupling of stock-weight process mean parameters (not used if stockWeightModel==0).",
    keyStockWeightObsVar = "Coupling of stock-weight observation variance parameters (not used if stockWeightModel==0).",
    catchWeightModel = c(
      "Integer code describing the treatment of catch weights in the model (0 use as known,",
      "1 use as observations to inform catch weight process (GMRF with cohort and within",
      "year correlations), 2 as 1 with extra correlation within the plus group)."),
    keyCatchWeightMean = c(
      "Coupling of catch-weight process mean parameters (not used if catchWeightModel==0).",
      "Only the first row is used."),
    keyCatchWeightObsVar = c(
      "Coupling of catch-weight observation variance parameters (not used if catchWeightModel==0).",
      "Only the first row is used."),
    matureModel = c(
      "Integer code describing the treatment of proportion mature in the model (0 use as known,",
      "1 use as observations to inform proportion mature process (GMRF with cohort and within",
      "year correlations on logit(proportion mature)), 2 as 1 with extra correlation within",
      "the plus group)."),
    keyMatureMean = "Coupling of mature process mean parameters (not used if matureModel==0).",
    mortalityModel = c(
      "Integer code describing the treatment of natural mortality in the model (0 use as known,",
      "1 use as observations to inform natural mortality process (GMRF with cohort and within",
      "year correlations), 2 as 1 with extra correlation within the plus group)."),
    keyMortalityMean = "Coupling of natural mortality process mean parameters (not used if mortalityModel==0).",
    keyMortalityObsVar = "Coupling of natural mortality observation variance parameters (not used if mortalityModel==0).",
    predVarObsLink = c(
      "Coupling of parameters used in a prediction-variance link for observations.",
      "For linked cells the observation sd is sqrt(log(1 + exp(k))) with",
      "k = log(sd) + (exp(p) - 1) * log(prediction). NA or -1 means no link."),
    keyXtraSd = c(
      "An integer matrix with 4 columns (fleet year age coupling), which allows additional",
      "uncertainty to be estimated for the specified observations (their sd is multiplied",
      "by exp(par))."),
    fixVarToWeight = c(
      "If weight attribute is supplied for observations this option sets the treatment",
      "(0 relative weight, 1 fix variance to weight). One value per fleet (a single value",
      "is recycled).")
  )
}

##' Save a samjr configuration to a SAM-style \code{.cfg} text file
##'
##' Writes a configuration list (as returned by \code{\link{defcon}} or
##' \code{\link{loadConf}}) to a text file in SAM's \code{$name} block
##' format, suitable for round-tripping back through \code{loadConf}.
##' Each block is preceded by comment lines describing the option and its
##' possible values. Matrices and age-indexed vectors get a comment line
##' with the ages (or age pairs, or column names) above the columns, and
##' rows of \code{keyParScaledYA} are labelled with their year.
##'
##' When \code{x} is a fitted \code{sam} object its configuration is saved,
##' and the data are used to list the fleets in the file header and to
##' label each row of the fleet matrices with the fleet name (as a trailing
##' \code{#} comment, which \code{loadConf} ignores).
##'
##' @param x a configuration list, or a fitted \code{sam} object from
##' \code{\link{sam.fit}}.
##' @param file output file path (or \code{""} for stdout).
##' @param overwrite if \code{FALSE} (default), refuse to overwrite an
##' existing file.
##' @return invisibly \code{NULL}.
##' @export
saveConf <- function(x, file = "", overwrite = FALSE){
  if(file != "" && file.exists(file) && !overwrite){
    cat("Notice: Did not overwrite existing file\n")
    return(invisible(NULL))
  }
  data <- NULL
  if(inherits(x, "sam")){
    data <- x$data
    x <- x$conf
  }
  fleetNames <- NULL
  fleetLegend <- character(0)
  if(!is.null(data)){
    fleetTypes <- data$fleetTypes
    nFleets <- length(fleetTypes)
    fleetNames <- attr(data, "fleetNames")
    fleetNames <- if(is.null(fleetNames)) paste("Fleet", seq_len(nFleets)) else trimws(fleetNames)
    aMin <- data$minAgePerFleet
    aMax <- data$maxAgePerFleet
    ageStr <- ifelse(is.na(aMin) | aMin < 0, "no age", paste0("ages ", aMin, "-", aMax))
    typeStr <- vapply(fleetTypes, fleetTypeName, character(1))
    fleetLegend <- c("# Fleets (rows of the fleet matrices):",
                     sprintf("# %3d  %s  %s, %s", seq_len(nFleets),
                             formatC(fleetNames, width = -max(nchar(fleetNames))),
                             typeStr, ageStr),
                     "#")
  }
  txt <- confComments()
  ages <- if(is.null(x$minAge) || is.null(x$maxAge)) integer(0) else x$minAge:x$maxAge
  nAges <- length(ages)
  ageVecKeys <- c("keyVarLogN", "keyStockWeightMean", "keyStockWeightObsVar",
                  "keyMatureMean", "keyMortalityMean", "keyMortalityObsVar")
  colLabels <- function(nm, v){
    n <- if(is.matrix(v)) ncol(v) else length(v)
    if(nm == "keyXtraSd" && n == 4) return(c("fleet", "year", "age", "coupling"))
    if(nm == "keyCorObs" && n == nAges - 1)
      return(paste(ages[-nAges], ages[-1], sep = "-"))
    if(n > 0 && n == nAges && (is.matrix(v) || nm %in% ageVecKeys))
      return(as.character(ages))
    NULL
  }
  rowLabels <- function(nm, v){
    if(!is.matrix(v) || nrow(v) == 0 || nm == "keyXtraSd") return(NULL)
    if(nm == "keyParScaledYA")
      return(if(nrow(v) == length(x$keyScaledYears)) as.character(x$keyScaledYears) else NULL)
    if(is.null(fleetNames)) return(NULL)
    catchNames <- fleetNames[fleetTypes == 0]
    if(nm %in% c("keyCatchWeightMean", "keyCatchWeightObsVar") && nrow(v) == length(catchNames))
      return(catchNames)
    if(nrow(v) == nFleets) return(fleetNames)
    NULL
  }
  formatValues <- function(v, lab, rowLab){
    m <- if(is.matrix(v)) v else matrix(v, nrow = 1)
    if(length(m) == 0) return("")
    s <- matrix(as.character(m), nrow = nrow(m))
    s[is.na(s)] <- "NA"
    if(is.null(lab) && !is.matrix(v)) return(paste0(" ", paste(s, collapse = " ")))
    w <- apply(rbind(nchar(s), if(is.null(lab)) NULL else nchar(lab)), 2, max) + 2
    w <- pmax(w, 4)
    padRow <- function(r) paste(sprintf("%*s", w, r), collapse = "")
    rows <- apply(s, 1, padRow)
    if(!is.null(rowLab)) rows <- paste0(rows, "  # ", rowLab)
    if(is.null(lab)) return(rows)
    c(paste0("#", substring(padRow(lab), 2)), rows)
  }
  out <- c(paste0("# Configuration saved: ", date()),
           "#",
           "# Where a matrix is specified rows correspond to fleets and columns to ages.",
           "# Same number indicates same parameter used",
           "# Numbers (integers) start from zero and must be consecutive",
           "# Negative numbers indicate that the parameter is not included in the model",
           "#",
           fleetLegend)
  for(nm in names(x)){
    v <- x[[nm]]
    out <- c(out, "", paste0("$", nm))
    if(!is.null(txt[[nm]])) out <- c(out, paste("#", txt[[nm]]))
    if(is.factor(v)){
      out <- c(out,
               paste("# Possible values are:", paste0('"', levels(v), '"', collapse = " ")),
               paste0(" ", paste0('"', as.character(v), '"', collapse = " ")))
    }else{
      out <- c(out, formatValues(v, colLabels(nm, v), rowLabels(nm, v)))
    }
  }
  cat(out, sep = "\n", file = file)
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
