##' Log likelihood of a samjr fit
##'
##' Returns the log-likelihood of a fitted samjr model as a \code{logLik}
##' object, with the number of estimated fixed-effect parameters as
##' \code{"df"} attribute and the number of observations as \code{"nobs"}
##' attribute. \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}} dispatch
##' on this method, so \code{AIC(fit)} and \code{BIC(fit)} work without
##' further code.
##'
##' @param object a fitted \code{sam} object from \code{\link{sam.fit}}.
##' @param ... extra arguments (ignored).
##' @return an object of class \code{logLik}.
##' @method logLik sam
##' @importFrom stats logLik
##' @export
logLik.sam <- function(object, ...){
  ret <- -object$opt$objective
  attr(ret, "df")   <- length(object$opt$par)
  attr(ret, "nobs") <- as.integer(object$data$nobs)
  class(ret) <- "logLik"
  ret
}

##' Print a samjr fit
##' @param x a \code{sam} object.
##' @param ... unused.
##' @method print sam
##' @export
print.sam <- function(x, ...){
  cat("samjr model: log likelihood is",
      as.numeric(logLik(x, ...)),
      " Convergence", ifelse(x$opt$convergence == 0, "OK\n", "failed\n"))
}

##' Print a \code{samset}
##' @param x a \code{samset} (list of \code{sam} fits).
##' @param ... unused.
##' @method print samset
##' @export
print.samset <- function(x, ...){
  cat("samset of", length(x), "fits\n")
  if(length(x) == 0) return(invisible(NULL))
  ll  <- vapply(x, function(f) as.numeric(logLik(f)),  numeric(1))
  npar<- vapply(x, function(f) attr(logLik(f), "df"),  integer(1))
  conv<- vapply(x, function(f) f$opt$convergence,      integer(1))
  res <- data.frame(name = if(is.null(names(x))) paste0("M", seq_along(x)) else names(x),
                    `log(L)` = round(ll, 4),
                    `#par`   = npar,
                    AIC      = round(-2 * ll + 2 * npar, 4),
                    conv     = conv,
                    check.names = FALSE)
  print(res, row.names = FALSE)
  invisible(res)
}

##' Summary of a samjr fit
##'
##' Per-year table of recruitment, SSB and Fbar with 95\% confidence
##' intervals (mirrors SAM's \code{summary.sam}).
##' @param object a \code{sam} fit.
##' @param ... unused.
##' @return a numeric matrix with year row names and 9 columns
##' (R/SSB/Fbar each times Estimate/Low/High).
##' @method summary sam
##' @export
summary.sam <- function(object, ...){
  ret <- cbind(round(rectable(object)),
               round(ssbtable(object)),
               round(fbartable(object), 3))
  colnames(ret)[1] <- paste0("R(age ", object$conf$minAge, ")")
  colnames(ret)[4] <- "SSB"
  colnames(ret)[7] <- paste0("Fbar(", object$conf$fbarRange[1], "-",
                              object$conf$fbarRange[2], ")")
  ret
}

##' One-step-ahead (OSA) quantile residuals from a samjr fit
##'
##' Calls \code{\link[RTMB]{oneStepPredict}} on the fitted \code{obj} to
##' compute one-observation-ahead quantile residuals.
##'
##' @param object a fitted \code{sam} object.
##' @param discrete logical; passed to \code{oneStepPredict}.
##' @param subset integer indices of observations to compute residuals
##' for (default: all).
##' @param ... extra arguments forwarded to \code{oneStepPredict}
##' (e.g. \code{method}).
##' @return a data frame of class \code{samres} with the columns from
##' \code{data$aux} (\code{year}, \code{fleet}, \code{age}) bound to the
##' \code{oneStepPredict} output.
##' @method residuals sam
##' @importFrom stats residuals
##' @importFrom RTMB oneStepPredict
##' @export
residuals.sam <- function(object, discrete = FALSE,
                          subset = 1:nobs(object), ...){
  cat("One-observation-ahead residuals. Total number of observations:",
      nobs(object), "\n")
  res <- oneStepPredict(object$obj, observation.name = "logobs",
                        discrete = discrete, subset = subset, ...)
  cat("One-observation-ahead residuals. Done\n")
  ret <- cbind(object$data$aux[subset, ], res)
  attr(ret, "fleetNames") <- attr(object$data, "fleetNames")
  class(ret) <- c("samres", "data.frame")
  ret
}

##' Print a \code{samres} residual data frame
##' @param x a \code{samres} object.
##' @param ... unused.
##' @method print samres
##' @export
print.samres <- function(x, ...){
  attr(x, "fleetNames") <- NULL
  class(x) <- setdiff(class(x), "samres")
  print(x)
}

##' Bubble plot of \code{samres} residuals
##'
##' One panel per fleet; bubble area scales with the absolute value of
##' the residual; positive residuals plotted in blue, negative in red.
##' Bubbles whose \code{|residual|} is closest to 1, 2, 3, ... are
##' labelled in place with that integer in grey, so the figure carries
##' its own size scale without a separate legend.
##' @param x a \code{samres} object.
##' @param bubblescale optional bubble-area scaling factor (default 1).
##' @param ... extra arguments passed to \code{plot}.
##' @method plot samres
##' @importFrom graphics par plot points abline text
##' @importFrom grDevices gray
##' @export
plot.samres <- function(x, bubblescale = 1, ...){
  fleets <- sort(unique(x$fleet))
  fnames <- attr(x, "fleetNames")
  if(is.null(fnames)) fnames <- paste0("Fleet ", fleets)
  nf <- length(fleets)
  oldpar <- par(mfrow = c(nf, 1), mar = c(4, 4, 2, 1))
  on.exit(par(oldpar))
  ok <- which(x$age >= 0 & is.finite(x$residual))
  absR <- abs(x$residual[ok])
  maxAbs <- max(absR)
  scale <- bubblescale * 5 / sqrt(maxAbs)
  ## For each integer k in 1..floor(maxAbs), pick (greedily, without
  ## reuse) the observation whose |residual| is closest to k.
  maxK <- min(floor(maxAbs), length(ok))
  labelRows <- integer(maxK)
  remaining <- absR
  for(k in seq_len(maxK)){
    i <- which.min(abs(remaining - k))
    labelRows[k] <- ok[i]
    remaining[i] <- NA
  }
  labelCex <- sqrt(4) * scale
  labelCol <- gray(0.2, alpha = 0.55)
  for(ff in fleets){
    sel <- x$fleet == ff & x$age >= 0 & is.finite(x$residual)
    if(!any(sel)) next
    yrR <- range(x$year[sel]); agR <- range(x$age[sel])
    plot(x$year[sel], x$age[sel], type = "n",
         xlab = "Year", ylab = "Age",
         xlim = yrR + c(-1, 1),
         ylim = agR + c(-0.7, 0.7),
         main = fnames[ff], ...)
    abline(h = sort(unique(x$age[sel])), lty = "dotted", col = "grey80")
    cols <- ifelse(x$residual[sel] < 0,
                   rgb(1, 0, 0, alpha = 0.5),
                   rgb(0, 0, 1, alpha = 0.5))
    points(x$year[sel], x$age[sel],
           cex = sqrt(abs(x$residual[sel])) * scale,
           pch = 19, col = cols)
    for(k in seq_len(maxK)){
      i <- labelRows[k]
      if(x$fleet[i] == ff){
        text(x$year[i], x$age[i], labels = k,
             cex = labelCex, font = 2, col = labelCol)
      }
    }
  }
  invisible(NULL)
}

##' Combine samjr fits into a \code{samset}
##' @param ... one or more \code{sam} fits.
##' @return a list of fits with class \code{samset}.
##' @method c sam
##' @export
c.sam <- function(...){
  ret <- list(...)
  class(ret) <- "samset"
  ret
}

##' Extract fixed-effect coefficients from a samjr fit
##' @param object a fitted \code{sam} object.
##' @param ... unused.
##' @return a named numeric vector with attributes \code{cov} (the
##' covariance) and \code{sd} (standard errors), class \code{samcoef}.
##' @method coef sam
##' @importFrom stats coef
##' @export
coef.sam <- function(object, ...){
  ret <- object$sdrep$par.fixed
  attr(ret, "cov") <- object$sdrep$cov.fixed
  attr(ret, "sd")  <- sqrt(diag(object$sdrep$cov.fixed))
  class(ret) <- "samcoef"
  ret
}

##' Print method for \code{samcoef}
##' @param x a \code{samcoef} object.
##' @param ... unused.
##' @method print samcoef
##' @export
print.samcoef <- function(x, ...){
  y <- as.vector(x); names(y) <- names(x); print(y)
}

##' Number of observations in a samjr fit
##'
##' @param object a fitted \code{sam} object from \code{\link{sam.fit}}.
##' @param ... extra arguments (ignored).
##' @return an integer.
##' @method nobs sam
##' @importFrom stats nobs
##' @export
nobs.sam <- function(object, ...){
  as.integer(object$data$nobs)
}

##' Concise markdown description of a fitted samjr model
##'
##' Walks \code{fit$data} and \code{fit$conf} to produce a short plain-text
##' (markdown) description of the data, process model, observation model,
##' catch scaling, fbar range, biology processes and basic fit statistics.
##' The text is printed (\code{cat}) and returned invisibly as a single
##' character vector.
##' @param fit a fitted \code{sam} object from \code{\link{sam.fit}}.
##' @return character vector with the markdown text (returned invisibly).
##' @export
modelDescription <- function(fit){
  data <- fit$data; conf <- fit$conf
  L <- character(0)
  push <- function(x) L[length(L) + 1L] <<- x

  fleetTypeName <- function(t) switch(as.character(t),
    "0" = "residual catch (catch-at-age)",
    "1" = "catch-at-age with effort",
    "2" = "survey index (numbers at age)",
    "3" = "biomass / catch index",
    "5" = "tagging",
    "6" = "alternative index",
    "7" = "sum-of-fleets",
    paste0("type ", t))

  ## Group an age vector by its key and produce e.g. "{1-3},{4},{5-6}".
  ## Entries with keyVec < 0 or NA are dropped.
  groupAges <- function(ages, keyVec){
    ok <- !is.na(keyVec) & keyVec >= 0
    if(!any(ok)) return("(none)")
    parts <- tapply(ages[ok], keyVec[ok], function(a){
      a <- sort(unique(a))
      if(length(a) == 1) return(as.character(a))
      if(all(diff(a) == 1)) return(sprintf("%d-%d", min(a), max(a)))
      paste(a, collapse = ",")
    })
    paste0("{", parts, "}", collapse = ", ")
  }
  yearRange <- function(yrs){
    yrs <- sort(unique(yrs))
    if(length(yrs) == 0) return("(no years)")
    if(length(yrs) == (max(yrs) - min(yrs) + 1L))
      sprintf("%d-%d", min(yrs), max(yrs))
    else
      sprintf("%d-%d (%d years, with gaps)",
              min(yrs), max(yrs), length(yrs))
  }

  fnames <- attr(data, "fleetNames")
  nF     <- data$noFleets
  if(is.null(fnames)) fnames <- paste("Fleet", seq_len(nF))
  ages   <- seq(conf$minAge, conf$maxAge)

  ## Header + Data
  push("## SAM model description")
  push("")
  push("### Data")
  push(sprintf("- Years: %d-%d (%d years)",
                min(data$years), max(data$years), length(data$years)))
  push(sprintf("- Age range: %d-%d", conf$minAge, conf$maxAge))
  push(sprintf("- Total observations: %d", length(data$logobs)))
  push(sprintf("- Fleets (%d):", nF))
  for(f in seq_len(nF)){
    aMin <- data$minAgePerFleet[f]; aMax <- data$maxAgePerFleet[f]
    aStr <- if(is.na(aMin) || aMin == -1) "(no age)"
            else sprintf("ages %d-%d", aMin, aMax)
    nobsF <- sum(data$aux[, "fleet"] == f)
    yrs <- data$aux[data$aux[, "fleet"] == f, "year"]
    push(sprintf("  %d. %s - %s, %s, years %s, %d obs",
                  f, fnames[f], fleetTypeName(data$fleetTypes[f]),
                  aStr, yearRange(yrs), nobsF))
  }
  push("")

  ## Process model
  push("### Process model")
  srLab <- c("0" = "random walk on log-recruitment",
             "1" = "Ricker", "2" = "Beverton-Holt")
  srKey <- as.character(conf$stockRecruitmentModelCode)
  push(sprintf("- Recruitment: %s",
                if(srKey %in% names(srLab)) srLab[[srKey]]
                else paste0("code ", srKey)))
  kF <- conf$keyLogFsta[1, ]
  nFstate <- length(unique(kF[kF >= 0]))
  push(sprintf("- F at age: %d state%s, ages coupled as %s",
                nFstate, if(nFstate == 1) "" else "s",
                groupAges(ages, kF)))
  corLab <- c("0" = "independent across ages",
              "1" = "compound symmetry (single shared correlation)",
              "2" = "AR(1)-like decay across ages")
  cKey <- as.character(conf$corFlag[1])
  push(sprintf("- F cross-age correlation: %s",
                if(cKey %in% names(corLab)) corLab[[cKey]]
                else paste0("code ", cKey)))
  vF <- conf$keyVarF[1, ]
  nVF <- length(unique(vF[vF >= 0]))
  push(sprintf("- F process variance: %d sd parameter%s, ages coupled as %s",
                nVF, if(nVF == 1) "" else "s",
                groupAges(ages, vF)))
  push("- N process: 2 sd parameters (recruitment, survival)")
  push("")

  ## Observation model
  push("### Observation model")
  vObs <- conf$keyVarObs
  for(f in seq_len(nF)){
    if(data$fleetTypes[f] == 5) next
    str <- as.character(conf$obsCorStruct[f])
    corStrLab <- switch(str,
      "ID" = "independent across ages",
      "AR" = "AR(1)/IGAR-style correlation across ages",
      "US" = "unstructured correlation",
      str)
    rowSd <- vObs[f, ]
    nSd   <- length(unique(rowSd[rowSd >= 0]))
    push(sprintf("- %s: %s", fnames[f], corStrLab))
    push(sprintf("  - sd: %d parameter%s, ages coupled as %s",
                  nSd, if(nSd == 1) "" else "s",
                  groupAges(ages, rowSd)))
  }
  push("")

  ## Catchability
  push("### Catchability (Q)")
  surveyFleets <- which(data$fleetTypes %in% c(2, 3, 6))
  if(length(surveyFleets) == 0){
    push("- (none - no survey-style fleets)")
  }else{
    for(f in surveyFleets){
      rowQ <- conf$keyLogFpar[f, ]
      nQ <- length(unique(rowQ[rowQ >= 0]))
      push(sprintf("- %s: %d logQ parameter%s, ages coupled as %s",
                    fnames[f], nQ, if(nQ == 1) "" else "s",
                    groupAges(ages, rowQ)))
    }
  }
  push("")

  ## Catch scaling
  if(isTRUE(conf$noScaledYears > 0)){
    push("### Catch scaling")
    yr <- conf$keyScaledYears
    push(sprintf("- %d scaled year(s): %s",
                  conf$noScaledYears, yearRange(yr)))
    keyMat <- conf$keyParScaledYA
    nSc <- max(keyMat[keyMat >= 0], na.rm = TRUE) + 1L
    push(sprintf("- %d distinct scaling parameter(s)", nSc))
    push("")
  }

  ## Fbar range
  push("### Fbar")
  push(sprintf("- Averaged over ages %d-%d",
                conf$fbarRange[1], conf$fbarRange[2]))
  push("")

  ## Biology processes
  bioFlags <- c(stockWeight = conf$stockWeightModel,
                catchWeight = conf$catchWeightModel,
                maturity    = conf$matureModel,
                mortality   = conf$mortalityModel)
  push("### Biology processes")
  for(nm in names(bioFlags)){
    fl <- bioFlags[[nm]]
    lab <- if(is.null(fl) || isTRUE(fl == 0)) "fixed (observed)"
           else sprintf("smoothed (GMRF model %d)", as.integer(fl))
    push(sprintf("- %s: %s", nm, lab))
  }
  push("")

  ## Fit summary
  push("### Fit summary")
  nll  <- if(!is.null(fit$opt$objective)) fit$opt$objective else NA_real_
  conv <- if(!is.null(fit$opt$convergence)) fit$opt$convergence else NA_integer_
  np   <- length(fit$obj$par)
  ll   <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
  aic  <- tryCatch(AIC(fit),                error = function(e) NA_real_)
  push(sprintf("- Negative log-likelihood: %.4f", nll))
  push(sprintf("- log-likelihood: %.4f", ll))
  push(sprintf("- Number of fixed-effect parameters: %d", np))
  push(sprintf("- AIC: %.4f", aic))
  push(sprintf("- Convergence code: %s (%s)", conv,
                if(isTRUE(conv == 0)) "converged" else "see ?nlminb"))
  txt <- paste(L, collapse = "\n")
  cat(txt, "\n", sep = "")
  invisible(txt)
}
