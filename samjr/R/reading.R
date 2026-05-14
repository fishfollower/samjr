##' Internal: \code{read.table} wrapper that suppresses the
##' "incomplete final line" warning often triggered by ICES-formatted files.
##' @param ... arguments forwarded to \code{\link[utils]{read.table}}.
##' @return the value returned by \code{read.table}.
##' @keywords internal
##' @noRd
read.table.nowarn <- function(...){
  tryCatch.W.E <- function(expr){
    W <- NULL
    w.handler <- function(w){
      if(!grepl('incomplete final line', w)) W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler), warning = W)
  }
  lis <- tryCatch.W.E(read.table(...))
  if(!is.null(lis$warning)) warning(lis$warning)
  lis$value
}

##' Internal: test whether a number is a whole non-negative number within tolerance.
##' @param x number.
##' @param tol numeric tolerance.
##' @return logical.
##' @keywords internal
##' @noRd
is.whole.positive.number <- function(x, tol = .Machine$double.eps^0.5){
  (abs(x - round(x)) < tol) & (x >= 0)
}

##' Read an ICES survey file
##'
##' Reads a multi-fleet ICES-format survey file and returns a list of matrices,
##' one per survey fleet. Each matrix has years as row names and ages as column
##' names; \code{attr(m, "time")} carries the start/end fraction of the year for
##' the survey.
##'
##' @param filen path to the survey file.
##' @return a named list of survey matrices.
##' @export
read.surveys <- function(filen){
  lin <- readLines(filen, warn = FALSE)[-c(1:2)]
  empty <- which(lapply(lapply(strsplit(lin, split = '[[:space:]]+'),
                               paste, collapse = ''), nchar) == 0)
  if(length(empty) > 0) lin <- lin[-empty]
  lin <- sub("^\\s+", "", lin)
  idx1 <- grep('^[A-Z#]', lin, ignore.case = TRUE)
  idx2 <- c(idx1[-1] - 1, length(lin))
  names <- lin[idx1]
  years <- matrix(as.numeric(unlist(strsplit(lin[idx1 + 1], '[[:space:]]+'))),
                  ncol = 2, byrow = TRUE)
  twofirst <- matrix(as.numeric(unlist(strsplit(lin[idx1 + 2], '[[:space:]]+'))),
                     ncol = 4, byrow = TRUE)[, 1:2, drop = FALSE]
  times <- matrix(as.numeric(unlist(strsplit(lin[idx1 + 2], '[[:space:]]+'))),
                  ncol = 4, byrow = TRUE)[, 3:4, drop = FALSE]
  ages <- matrix(as.numeric(unlist(lapply(strsplit(lin[idx1 + 3], '[[:space:]]+'),
                                          function(x) x[1:2]))),
                 ncol = 2, byrow = TRUE)
  for(i in 1:length(names)){
    if(!is.whole.positive.number(years[i, 1]))
      stop(paste("In file", filen, ": Minimum year is expected to be a positive integer for fleet number", i))
    if(!is.whole.positive.number(years[i, 2]))
      stop(paste("In file", filen, ": Maximum year is expected to be a positive integer for fleet number", i))
    if(years[i, 1] > years[i, 2])
      stop(paste("In file", filen, ": Maximum year must be greater than minimum year for fleet number", i))
    if(ages[i, 1] > ages[i, 2])
      stop(paste("In file", filen, ": Maximum age must be greater than minimum age for fleet number", i))
    if((times[i, 1] < 0) | (times[i, 1] > 1))
      stop(paste("In file", filen, ": Minimum survey time must be within [0,1] for fleet number", i))
    if((times[i, 2] < 0) | (times[i, 2] > 1))
      stop(paste("In file", filen, ": Maximum survey time must be within [0,1] for fleet number", i))
    if(times[i, 2] < times[i, 1])
      stop(paste("In file", filen, ": Maximum survey time must be greater than minimum survey time for fleet number", i))
  }

  as.num <- function(x, na.strings = "NA"){
    stopifnot(is.character(x))
    na <- x %in% na.strings
    x[na] <- 0
    x <- as.numeric(x)
    x[na] <- NA_real_
    x
  }

  onemat <- function(i){
    lin.local <- gsub('^[[:blank:]]*', '', lin[(idx1[i] + 4):idx2[i]])
    nr <- idx2[i] - idx1[i] - 3
    ret <- matrix(as.num(unlist((strsplit(lin.local, '[[:space:]]+')))),
                  nrow = nr, byrow = TRUE)[, , drop = FALSE]
    if(nrow(ret) != (years[i, 2] - years[i, 1] + 1))
      stop(paste("In file", filen, ": Year range does not match number of rows for survey fleet", i))
    if((ncol(ret) - 1) < (ages[i, 2] - ages[i, 1] + 1))
      stop(paste("In file", filen, ": Fewer columns than age range for survey fleet", i))
    if(!is.numeric(ret))
      stop(paste("In file", filen, ": Non numeric values for survey fleet", i))
    ret <- as.matrix(ret[, -1] / ret[, 1])
    rownames(ret) <- years[i, 1]:years[i, 2]
    ret <- ret[, 1:length(ages[i, 1]:ages[i, 2]), drop = FALSE]
    colnames(ret) <- ages[i, 1]:ages[i, 2]
    attr(ret, 'time') <- times[i, ]
    attr(ret, 'twofirst') <- twofirst[i, ]
    ret[ret < 0] <- NA
    ret
  }
  obs <- lapply(1:length(names), onemat)
  names(obs) <- names
  obs
}

##' Read an ICES/CEFAS format data file
##'
##' Reads ICES-formatted catch, weight, mortality, maturity, and survey files
##' and returns a validated matrix (or, for survey files, a list of matrices).
##' The first two lines are ignored and may be used for comments. The third
##' line carries either header numbers (catch-style files: minYear, maxYear,
##' minAge, maxAge, datatype) or a fleet name (survey files).
##'
##' Supported catch-style datatype codes are 1 (full year-by-age matrix),
##' 2 (single row repeated across years), 3 (single scalar repeated) and
##' 5 (single column repeated across ages).
##'
##' @param filen the file to read.
##' @return for catch-style files a numeric matrix with year row names and age
##' column names; for survey files a named list of such matrices (one per fleet).
##' @details
##' First two lines are ignored and can be used for comments. Tests are
##' performed on the format code and on the year/age ranges.
##' @export
read.ices <- function(filen){
  if(grepl("^[0-9]", scan(filen, skip = 2, n = 1, quiet = TRUE, what = ""))){
    head <- scan(filen, skip = 2, n = 5, quiet = TRUE)
    minY <- head[1]; maxY <- head[2]; minA <- head[3]; maxA <- head[4]
    datatype <- head[5]
    if(!is.whole.positive.number(minY))
      stop(paste("In file", filen, ": Minimum year must be a positive integer"))
    if(!is.whole.positive.number(maxY))
      stop(paste("In file", filen, ": Maximum year must be a positive integer"))
    if(!is.whole.positive.number(minA))
      stop(paste("In file", filen, ": Minimum age must be a positive integer"))
    if(!is.whole.positive.number(maxA))
      stop(paste("In file", filen, ": Maximum age must be a positive integer"))
    if(!(datatype %in% c(1, 2, 3, 5)))
      stop(paste("In file", filen, ": Datatype must be one of 1, 2, 3, 5"))
    if(minY > maxY) stop(paste("In file", filen, ": minY > maxY"))
    if(minA > maxA) stop(paste("In file", filen, ": minA > maxA"))

    C <- as.matrix(read.table.nowarn(filen, skip = 5, header = FALSE))

    if(datatype == 1){
      if((maxY - minY + 1) != nrow(C))
        stop(paste("In file", filen, ": Number of rows does not match year range"))
      if((maxA - minA + 1) > ncol(C))
        stop(paste("In file", filen, ": Fewer columns than age range"))
    }
    if(datatype == 2){
      if(1 != nrow(C)) stop(paste("In file", filen, ": For datatype 2 only one row expected"))
      if((maxA - minA + 1) > ncol(C)) stop(paste("In file", filen, ": Fewer columns than age range"))
      C <- C[rep(1, maxY - minY + 1), ]
    }
    if(datatype == 3){
      if(1 != nrow(C)) stop(paste("In file", filen, ": For datatype 3 only one row expected"))
      if(1 != ncol(C)) stop(paste("In file", filen, ": For datatype 3 only one column expected"))
      C <- C[rep(1, maxY - minY + 1), rep(1, maxA - minA + 1)]
    }
    if(datatype == 5){
      if((maxY - minY + 1) != nrow(C))
        stop(paste("In file", filen, ": Number of rows does not match year range"))
      if(1 != ncol(C)) stop(paste("In file", filen, ": For datatype 5 only one column expected"))
      C <- C[, rep(1, maxA - minA + 1)]
    }
    rownames(C) <- minY:maxY
    C <- C[, 1:length(minA:maxA)]
    colnames(C) <- minA:maxA
    if(!is.numeric(C))
      stop(paste("In file", filen, ": Non numeric data values detected"))

    l <- readLines(filen, warn = FALSE)
    lWithAttribName <- grep("^[[:blank:]]*#+'[[:blank:]]*@", l)
    lWithAttribValues <- grep("^[[:blank:]]*#+'[[:blank:]]*[^@]", l)
    lWithAttribValues <- setdiff(lWithAttribValues, lWithAttribName)
    if(length(lWithAttribName) > 0){
      nms <- sapply(l[lWithAttribName],
                    function(v) gsub("[[:blank:]]+$", "",
                                     gsub("(^[[:blank:]]*#+'[[:blank:]]*@)", "", v)))
      valList <- split(l[lWithAttribValues],
                       sapply(lWithAttribValues, function(v) sum(v > lWithAttribName)))
      customAttrib <- lapply(valList, function(l){
        eval(parse(text = paste(gsub("^[[:blank:]]*#+'[[:blank:]]*", "", l), collapse = "\n")))
      })
      names(customAttrib) <- nms
      attributes(C) <- c(attributes(C), customAttrib)
    }
    return(C)
  }else{
    return(read.surveys(filen))
  }
}

##' Read all standard SAM data files from a directory
##'
##' Convenience wrapper that reads the standard ICES-format files in a stock
##' directory and assembles a samjr data list via
##' \code{\link{setup.sam.data}}.
##'
##' Files read:
##' \describe{
##'   \item{\code{cn.dat}}{catch numbers at age (residual fleet).}
##'   \item{\code{cw.dat}}{catch mean weight at age.}
##'   \item{\code{dw.dat}}{discard mean weight at age.}
##'   \item{\code{lw.dat}}{landing mean weight at age.}
##'   \item{\code{mo.dat}}{proportion mature at age.}
##'   \item{\code{nm.dat}}{natural mortality at age.}
##'   \item{\code{pf.dat}}{proportion of F before spawning.}
##'   \item{\code{pm.dat}}{proportion of M before spawning.}
##'   \item{\code{sw.dat}}{stock mean weight at age.}
##'   \item{\code{lf.dat}}{landing fraction at age.}
##'   \item{\code{survey.dat}}{survey indices.}
##' }
##' @param dir directory to read from.
##' @return list as produced by \code{\link{setup.sam.data}}.
##' @export
read.data.files <- function(dir = "."){
  od <- setwd(dir); on.exit(setwd(od))
  cn <- read.ices("cn.dat")
  cw <- read.ices("cw.dat")
  dw <- read.ices("dw.dat")
  lw <- read.ices("lw.dat")
  mo <- read.ices("mo.dat")
  nm <- read.ices("nm.dat")
  pf <- read.ices("pf.dat")
  pm <- read.ices("pm.dat")
  sw <- read.ices("sw.dat")
  lf <- read.ices("lf.dat")
  surveys <- read.ices("survey.dat")
  setup.sam.data(surveys = surveys,
                 residual.fleet = cn,
                 prop.mature = mo,
                 stock.mean.weight = sw,
                 catch.mean.weight = cw,
                 dis.mean.weight = dw,
                 land.mean.weight = lw,
                 prop.f = pf,
                 prop.m = pm,
                 natural.mortality = nm,
                 land.frac = lf)
}

##' Combine the data sources to a samjr data object
##'
##' Builds the data list consumed by \code{\link{sam.fit}}. The interface
##' mirrors the SAM \code{stockassessment::setup.sam.data} but is restricted
##' to a single residual catch fleet and to two-dimensional (year x age)
##' inputs. Multi-fleet 3D-array catch inputs and the
##' \code{recapture}/\code{sum.residual.fleets} machinery are not supported.
##'
##' @param fleets reserved for future use; must be \code{NULL}.
##' @param surveys a single survey matrix or a named list of survey matrices.
##' @param residual.fleet a single residual-catch matrix with year rows and age columns.
##' @param prop.mature year-by-age proportion mature.
##' @param stock.mean.weight year-by-age mean weight in the stock.
##' @param catch.mean.weight year-by-age mean weight in the catch.
##' @param dis.mean.weight year-by-age mean weight of discards.
##' @param land.mean.weight year-by-age mean weight of landings.
##' @param natural.mortality year-by-age natural mortality.
##' @param prop.f year-by-age proportion of F before spawning.
##' @param prop.m year-by-age proportion of M before spawning.
##' @param land.frac year-by-age landing fraction.
##' @return a list of class \code{sam_data} containing fleet metadata
##' (\code{fleetTypes}, \code{sampleTimes}, \code{minAgePerFleet},
##' \code{maxAgePerFleet}), the observation table (\code{aux}, \code{logobs},
##' \code{idx1}, \code{idx2}, \code{nobs}), and the year-by-age biological
##' inputs (\code{propMat}, \code{stockMeanWeight}, \code{catchMeanWeight},
##' \code{natMor}, \code{landFrac}, \code{disMeanWeight},
##' \code{landMeanWeight}, \code{propF}, \code{propM}).
##' @export
setup.sam.data <- function(fleets = NULL, surveys = NULL, residual.fleet = NULL,
                           prop.mature = NULL, stock.mean.weight = NULL,
                           catch.mean.weight = NULL, dis.mean.weight = NULL,
                           land.mean.weight = NULL, natural.mortality = NULL,
                           prop.f = NULL, prop.m = NULL, land.frac = NULL,
                           recapture = NULL){
  if(!is.null(fleets))
    stop("samjr v1: 'fleets' (commercial fleets with effort) not supported")
  if(is.null(residual.fleet))
    stop("samjr v1: a single 'residual.fleet' (matrix) is required")
  if(!(is.matrix(residual.fleet) || is.data.frame(residual.fleet)))
    stop("samjr v1: 'residual.fleet' must be a single matrix; multi-fleet not supported")

  fleet.idx <- 0
  type <- NULL; time <- NULL; name <- NULL
  dat <- data.frame(year = NA_integer_, fleet = NA_integer_,
                    age = NA_integer_, aux = NA_integer_)
  fleetAges <- list()
  weight <- NULL
  corList <- list()
  fleetCovEntries <- list()   # one entry per fleet: list(perObs=, perRow=)

  doone <- function(m){
    yearV <- rownames(m)[row(m)]
    fleet.idx <<- fleet.idx + 1
    fleetV <- rep(fleet.idx, length(yearV))
    ageV <- as.integer(colnames(m)[col(m)])
    fleetAges[[fleet.idx]] <<- as.integer(colnames(m))
    auxV <- as.vector(m)
    dat <<- rbind(dat, data.frame(year = yearV, fleet = fleetV,
                                  age = ageV, aux = auxV))
    A <- attributes(m)
    if(!is.null(A[["cov"]])){
      cv <- A[["cov"]]
      attr(m, "cor") <- lapply(cv, function(x) cov2cor(x))
      wPerRow <- do.call(rbind, lapply(cv, diag))
      weight <<- c(weight, as.vector(wPerRow))
    }else if(!is.null(A[["cov-weight"]])){
      cv <- A[["cov-weight"]]
      attr(m, "cor") <- lapply(cv, function(x) cov2cor(x))
      wPerRow <- do.call(rbind, lapply(cv, diag))
      weight <<- c(weight, 1 / as.vector(wPerRow))
    }else if(!is.null(A[["weight"]])){
      weight <<- c(weight, as.vector(A[["weight"]]))
    }else{
      weight <<- c(weight, rep(NA_real_, length(yearV)))
    }
    A2 <- attributes(m)
    fleetCovEntries[[fleet.idx]] <<- list(cor = A2[["cor"]],
                                          rownames = rownames(m))
  }

  doone(residual.fleet)
  type <- c(type, 0); time <- c(time, 0); name <- c(name, "Residual catch")

  if(!is.null(surveys)){
    if(is.matrix(surveys) || is.data.frame(surveys)){
      doone(surveys)
      thistype <- ifelse(min(as.integer(colnames(surveys))) < (-.5), 3, 2)
      type <- c(type, thistype)
      time <- c(time, mean(attr(surveys, 'time')))
      name <- c(name, "Survey fleet")
    }else{
      lapply(surveys, doone)
      type <- c(type, unlist(lapply(surveys,
                  function(x) ifelse(min(as.integer(colnames(x))) < (-.5), 3, 2))))
      time <- c(time, unlist(lapply(surveys, function(x) mean(attr(x, 'time')))))
      name <- c(name, strtrim(gsub("\\s", "", names(surveys)), 50))
    }
  }

  dat$aux[which(dat$aux <= 0)] <- NA_integer_
  dat <- dat[!is.na(dat$year), ]

  if(!is.null(recapture)){
    fleet.idx <- fleet.idx + 1
    tag <- data.frame(
      year       = recapture$ReleaseY,
      fleet      = fleet.idx,
      age        = recapture$ReleaseY - recapture$Yearclass,
      aux        = exp(recapture$r),
      RecaptureY = recapture$RecaptureY,
      Yearclass  = recapture$Yearclass,
      Nscan      = recapture$Nscan,
      R          = recapture$R,
      Type       = recapture$Type
    )
    extra <- setdiff(names(tag), names(dat))
    for(nm in extra) dat[[nm]] <- NA_integer_
    dat <- rbind(dat, tag)
    weight <- c(weight, rep(NA_real_, nrow(tag)))
    type <- c(type, 5); time <- c(time, 0); name <- c(name, "Recaptures")
    fleetCovEntries[[fleet.idx]] <- list(cor = NULL, rownames = NULL)
  }

  cc <- which(complete.cases(dat[, 1:3]) == TRUE)
  dat <- dat[cc, ]

  o <- order(as.numeric(dat$year), as.numeric(dat$fleet), as.numeric(dat$age))
  dat <- dat[o, ]
  weight <- weight[o]

  newyear <- min(as.numeric(dat$year)):max(as.numeric(dat$year))
  newfleet <- min(as.numeric(dat$fleet)):max(as.numeric(dat$fleet))
  mmfun <- function(f, y, ff){
    idx <- which(dat$year == y & dat$fleet == f)
    if(length(idx) == 0) NA_integer_ else ff(idx) - 1L
  }
  idx1 <- outer(newfleet, newyear, Vectorize(mmfun, c("f", "y")), ff = min)
  idx2 <- outer(newfleet, newyear, Vectorize(mmfun, c("f", "y")), ff = max)

  minAgePerFleet <- tapply(as.integer(dat[, "age"]),
                           INDEX = factor(dat[, "fleet"], seq_len(fleet.idx)),
                           FUN = min)
  maxAgePerFleet <- tapply(as.integer(dat[, "age"]),
                           INDEX = factor(dat[, "fleet"], seq_len(fleet.idx)),
                           FUN = max)

  cutY <- function(x){
    if(is.null(x)) return(NULL)
    rs <- rownames(x)
    ret <- matrix(NA_real_, nrow = length(newyear), ncol = ncol(x))
    rownames(ret) <- as.character(newyear)
    colnames(ret) <- colnames(x)
    mi <- match(as.character(newyear), rs)
    have <- !is.na(mi)
    if(any(have)) ret[have, ] <- x[mi[have], ]
    ret
  }

  if(is.null(prop.m))
    prop.m <- matrix(0, nrow = nrow(residual.fleet), ncol = ncol(residual.fleet))
  if(is.null(land.frac))
    land.frac <- matrix(1, nrow = nrow(residual.fleet), ncol = ncol(residual.fleet))
  if(is.null(prop.f))
    prop.f <- matrix(0, nrow = nrow(residual.fleet), ncol = ncol(residual.fleet))

  nF <- length(type)
  idxCor <- matrix(NA_integer_, nrow = nF, ncol = length(newyear))
  colnames(idxCor) <- as.character(newyear)
  for(fi in seq_len(nF)){
    entry <- fleetCovEntries[[fi]]
    if(is.null(entry$cor)) next
    cl   <- entry$cor
    rnms <- entry$rownames
    okIdx <- which(vapply(cl, function(x) !any(is.na(x)), logical(1)))
    if(length(okIdx) == 0) next
    addBase <- length(corList)
    for(k in seq_along(okIdx)){
      corList[[addBase + k]] <- cl[[okIdx[k]]]
    }
    yrCols <- as.integer(rnms[okIdx])
    matchCols <- match(as.character(yrCols), colnames(idxCor))
    valid <- !is.na(matchCols)
    idxCor[fi, matchCols[valid]] <- (addBase + seq_along(okIdx))[valid]
  }

  ret <- list(
    noFleets = length(type),
    fleetTypes = as.integer(type),
    sampleTimes = time,
    noYears = length(newyear),
    years = newyear,
    minAgePerFleet = minAgePerFleet,
    maxAgePerFleet = maxAgePerFleet,
    nobs = nrow(dat),
    idx1 = idx1,
    idx2 = idx2,
    aux = do.call(cbind, lapply(dat[, -4], as.integer)),
    logobs = log(dat[, 4]),
    weight = as.numeric(weight),
    corList = corList,
    idxCor = idxCor,
    propMat = cutY(prop.mature),
    stockMeanWeight = cutY(stock.mean.weight),
    catchMeanWeight = cutY(catch.mean.weight),
    natMor = cutY(natural.mortality),
    landFrac = cutY(land.frac),
    disMeanWeight = cutY(dis.mean.weight),
    landMeanWeight = cutY(land.mean.weight),
    propF = cutY(prop.f),
    propM = cutY(prop.m)
  )
  attr(ret, "fleetNames") <- name
  class(ret) <- "sam_data"
  ret
}
