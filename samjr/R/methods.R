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
##' @param x a \code{samres} object.
##' @param bubblescale optional bubble-area scaling factor (default 1).
##' @param ... extra arguments passed to \code{plot}.
##' @method plot samres
##' @importFrom graphics par plot points abline mtext text plot.new plot.window
##' @export
plot.samres <- function(x, bubblescale = 1, ...){
  fleets <- sort(unique(x$fleet))
  fnames <- attr(x, "fleetNames")
  if(is.null(fnames)) fnames <- paste0("Fleet ", fleets)
  oldpar <- par(mfrow = c(length(fleets) + 1, 1), mar = c(4, 4, 2, 1))
  on.exit(par(oldpar))
  maxAbs <- max(abs(x$residual), na.rm = TRUE)
  scale <- bubblescale * 5 / sqrt(maxAbs)
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
  }
  legendVals <- pretty(c(0, maxAbs), n = 4)
  legendVals <- legendVals[legendVals > 0]
  par(mar = c(2, 4, 2, 1))
  plot.new()
  plot.window(xlim = c(0, length(legendVals) + 1), ylim = c(-1, 1))
  xs <- seq_along(legendVals)
  points(xs, rep(0.2, length(xs)),
         cex = sqrt(legendVals) * scale,
         pch = 19, col = rgb(0, 0, 1, alpha = 0.5))
  text(xs, rep(-0.7, length(xs)),
       labels = format(legendVals), cex = 0.9)
  mtext("|residual|", side = 3, line = 0, cex = 0.9)
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
