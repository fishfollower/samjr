## Plot / SR-curve overlay routines for samref.

##' Plot a sam_referencepoints object
##'
##' Base-graphics curves of one equilibrium / per-recruit quantity against
##' F-bar, with vertical lines marking each located reference point.
##'
##' @param x a \code{sam_referencepoints} object.
##' @param type one of \code{"yield"} (default), \code{"ypr"}, \code{"spr"},
##'   \code{"ssb"}, \code{"rec"}, \code{"yearsLost"}, \code{"lifeExpectancy"},
##'   \code{"all"}.
##' @param ... extra arguments passed to \code{plot}.
##' @return invisible \code{NULL}.
##' @method plot sam_referencepoints
##' @export
plot.sam_referencepoints <- function(x,
                                     type = c("yield","ypr","spr","ssb","rec",
                                              "yearsLost","lifeExpectancy","all"),
                                     ...){
  type <- match.arg(type)
  panels <- list(yield = "Yield", ypr = "YieldPerRecruit",
                 spr = "SpawnersPerRecruit", ssb = "Biomass",
                 rec = "Recruitment", yearsLost = "YearsLost",
                 lifeExpectancy = "LifeExpectancy")
  panelLab <- list(yield = "Equilibrium yield", ypr = "Yield per recruit",
                   spr = "Spawners per recruit", ssb = "Equilibrium SSB",
                   rec = "Equilibrium recruitment",
                   yearsLost = "Years lost to fishing",
                   lifeExpectancy = "Life expectancy")
  if(type == "all"){
    keep <- names(panels)
    keep <- keep[vapply(keep, function(k) any(is.finite(x$graphs[[panels[[k]]]])),
                        logical(1))]
    oldpar <- par(mfrow = c(ceiling(length(keep) / 2), 2))
    on.exit(par(oldpar))
    for(k in keep) .plotOnePanel(x, panels[[k]], panelLab[[k]], ...)
    return(invisible(NULL))
  }
  .plotOnePanel(x, panels[[type]], panelLab[[type]], ...)
  invisible(NULL)
}

.plotOnePanel <- function(x, fld, lab, ...){
  yy <- x$graphs[[fld]]
  finite <- is.finite(yy)
  if(!any(finite)){
    plot.new(); title(main = paste(lab, "(not available)")); return(invisible())
  }
  Fs <- x$graphs$F[finite]; yy <- yy[finite]
  plot(Fs, yy, type = "l", xlab = x$fbarlabel, ylab = lab, ...)
  ## RP vertical lines
  if(!is.null(x$tables$F)){
    rpFs <- x$tables$F[, "Estimate"]
    cols <- grDevices::hcl.colors(length(rpFs), palette = "Dark 3")
    for(i in seq_along(rpFs)){
      if(!is.finite(rpFs[i])) next
      abline(v = rpFs[i], col = cols[i], lwd = 2, lty = 2)
    }
    legend("topright", legend = rownames(x$tables$F),
           col = cols, lwd = 2, lty = 2, bty = "n", cex = 0.7)
  }
}

##' Add a stock-recruit curve to the current plot
##'
##' Overlays the fitted Ricker / Beverton-Holt curve on a SSB-vs-R scatter
##' (e.g. one drawn by \code{samjr::srplot}). Uses the delta-method to
##' attach 95\% pointwise CIs.
##'
##' @param fit a samjr \code{sam} fit with SR code 1 or 2.
##' @param ssbRange optional numeric length-2 vector for the x-axis range;
##'   default is taken from the current plotting region.
##' @param n number of grid points.
##' @param col,lwd line aesthetics.
##' @param ciCol fill colour for the CI ribbon.
##' @param ... unused.
##' @return invisible \code{NULL}.
##' @examples
##' \donttest{
##' data(nscodData); data(nscodConf)
##' fit <- samjr::sam.fit(nscodData, nscodConf,
##'                       samjr::defpar(nscodData, nscodConf), silent = TRUE)
##' fit$conf$stockRecruitmentModelCode <- 2L
##' fitBH <- samjr::runwithout(fit)
##' samjr::srplot(fitBH)
##' addRecruitmentCurve(fitBH)
##' }
##' @export
addRecruitmentCurve <- function(fit, ...) UseMethod("addRecruitmentCurve")

##' @rdname addRecruitmentCurve
##' @method addRecruitmentCurve sam
##' @export
addRecruitmentCurve.sam <- function(fit, ssbRange = NULL, n = 200,
                                    col = "black", lwd = 2,
                                    ciCol = grDevices::rgb(0, 0, 0, 0.2),
                                    ...){
  srCode <- fit$conf$stockRecruitmentModelCode
  if(!srCode %in% c(1L, 2L))
    stop("addRecruitmentCurve: SR code ", srCode,
         " has no curve to draw (only Ricker / Beverton-Holt supported)")
  rp <- recPars(fit)
  nm <- if(srCode == 1L) "rickerpar" else "bhpar"
  if(is.null(ssbRange)){
    usr <- par("usr"); ssbRange <- c(usr[1], usr[2])
  }
  ssb <- seq(max(ssbRange[1], 1e-6), ssbRange[2], length = n)
  logR <- predictLogR(srCode, rp, log(ssb))
  cv <- fit$sdrep$cov.fixed
  if(!is.null(cv) && all(nm %in% names(diag(cv)))){
    keep <- names(diag(cv)) %in% nm
    Sigma <- cv[keep, keep, drop = FALSE]
    h <- 1e-4
    J <- matrix(0, n, length(rp))
    for(k in seq_along(rp)){
      dp <- rep(0, length(rp)); dp[k] <- h
      fhi <- predictLogR(srCode, rp + dp, log(ssb))
      flo <- predictLogR(srCode, rp - dp, log(ssb))
      J[, k] <- (fhi - flo) / (2 * h)
    }
    se <- sqrt(pmax(diag(J %*% Sigma %*% t(J)), 0))
    lo <- exp(logR - 1.96 * se); hi <- exp(logR + 1.96 * se)
    polygon(c(ssb, rev(ssb)), c(lo, rev(hi)), col = ciCol, border = NA)
  }
  lines(ssb, exp(logR), col = col, lwd = lwd)
  invisible(NULL)
}
