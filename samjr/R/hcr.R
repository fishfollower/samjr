## Harvest control rules for samref.
##
## Port of SAM's hcr / icesAdviceRule, evaluating a trapezoidal F-vs-SSB
## rule via the self-contained per-year projection in projection.R.

##' Trapezoidal harvest control rule
##'
##' Returns Fbar given SSB:
##' \deqn{F(SSB) = F_{cap} \quad \text{if } SSB < B_{cap}}
##' \deqn{F(SSB) = \min(F_t, \max(F_o,
##'        (SSB - B_o)(F_t - F_o)/(B_{tr} - B_o)))
##'        \quad \text{if } SSB \ge B_{cap}, B_{tr} \ne B_o}
##' \deqn{F(SSB) = F_t \quad \text{if } SSB \ge B_{cap}, B_{tr} = B_o.}
##'
##' Same form as SAM's \code{hcrFun}.
##'
##' @param ssb scalar SSB.
##' @param Ftarget,Btrigger,Forigin,Borigin,Fcap,Bcap rule parameters.
##' @return scalar F.
##' @examples
##' ## ICES-style rule: F = 0 below 100k, ramp to 0.22 by 150k SSB
##' hcrFun( 50000, Ftarget = 0.22, Btrigger = 150000,
##'        Bcap = 100000, Borigin = 100000)   # -> 0
##' hcrFun(120000, Ftarget = 0.22, Btrigger = 150000,
##'        Bcap = 100000, Borigin = 100000)   # -> 0.088
##' hcrFun(200000, Ftarget = 0.22, Btrigger = 150000,
##'        Bcap = 100000, Borigin = 100000)   # -> 0.22
##' @export
hcrFun <- function(ssb, Ftarget, Btrigger, Forigin = 0, Borigin = 0,
                   Fcap = 0, Bcap = 0){
  if(ssb < Bcap) return(Fcap)
  if(Btrigger == Borigin) return(Ftarget)
  Flin <- (ssb - Borigin) * (Ftarget - Forigin) / (Btrigger - Borigin)
  min(Ftarget, max(Forigin, Flin))
}

##' Project a samjr fit forward under a harvest control rule
##'
##' Year-by-year stochastic projection where each year's Fbar is set by
##' \code{\link{hcrFun}(ssb, ...)} applied to that year's simulated SSB.
##' Recruitment is sampled from the SR predictor plus a Gaussian residual
##' with sd taken from the fit (\code{exp(fit$pl$logSdLogN[1])}).
##'
##' @param fit a samjr \code{sam} fit.
##' @param nYears number of forecast years.
##' @param Ftarget,Btrigger,Forigin,Borigin,Fcap,Bcap rule parameters.
##' @param nosim number of Monte-Carlo simulations.
##' @param aveYears years over which to average biology (default last 5).
##' @param selYears years over which to average selectivity (default last
##'   year).
##' @param customSel optional selectivity override.
##' @param seed optional RNG seed.
##' @param ... unused.
##' @return list of class \code{samref_hcr} with components
##'   \code{projection} (a per-year Estimate/Low/High table for SSB, Fbar,
##'   Recruitment, Catch), \code{rule} parameters, and the rule as
##'   \code{hcr} (a closure).
##' @examples
##' \donttest{
##' data(nscodData); data(nscodConf)
##' fit <- samjr::sam.fit(nscodData, nscodConf,
##'                       samjr::defpar(nscodData, nscodConf), silent = TRUE)
##' fit$conf$stockRecruitmentModelCode <- 2L
##' fitBH <- samjr::runwithout(fit)
##' hh <- hcr(fitBH, Ftarget = 0.22, Btrigger = 150000,
##'           Borigin = 100000, Bcap = 100000,
##'           nosim = 100, nYears = 5, seed = 1)
##' hh$projection$SSB
##' }
##' @export
hcr <- function(fit, ...) UseMethod("hcr")

##' @rdname hcr
##' @method hcr sam
##' @export
hcr.sam <- function(fit,
                    nYears = 20,
                    Ftarget,
                    Btrigger,
                    Forigin = 0, Borigin = 0,
                    Fcap = 0, Bcap = 0,
                    nosim = 1000,
                    aveYears  = NULL,
                    selYears  = NULL,
                    customSel = NULL,
                    seed = NULL,
                    ...){
  if(!is.null(seed)) set.seed(seed)
  yrs <- fit$data$years
  if(is.null(aveYears)) aveYears <- yrs[max(1L, length(yrs) - 4L):length(yrs)]
  if(is.null(selYears)) selYears <- yrs[length(yrs)]
  aveIdx <- .aveYearsIndex(fit, aveYears)
  bio <- averageBio(fit, aveIdx)
  sel <- if(!is.null(customSel)) customSel else selectivityFromLogF(fit, selYears)
  srCode <- fit$conf$stockRecruitmentModelCode
  rp <- recPars(fit)
  st <- .initState(fit)
  nAge <- length(st$N)
  maxAgePlusGroup <- isTRUE(fit$conf$maxAgePlusGroup[1] == 1L)
  Nsim <- matrix(rep(st$N, each = nosim), nrow = nosim)
  ssbMat <- matrix(NA_real_, nosim, nYears)
  fMat   <- matrix(NA_real_, nosim, nYears)
  cMat   <- matrix(NA_real_, nosim, nYears)
  rMat   <- matrix(NA_real_, nosim, nYears)
  ruleClosure <- function(ssb)
    hcrFun(ssb, Ftarget, Btrigger, Forigin, Borigin, Fcap, Bcap)
  for(y in seq_len(nYears)){
    for(s in seq_len(nosim)){
      Ns <- Nsim[s, ]
      ssbCur <- sum(Ns * exp(-sel * 0 * bio$propF - bio$natMor * bio$propM) *
                    bio$propMat * bio$stockMeanWeight)
      fval <- ruleClosure(ssbCur)
      step <- .projectOneYear(Ns, fval, bio, sel, srCode, rp,
                              recNoise = stats::rnorm(1, 0, st$sigmaR),
                              maxAgePlusGroup = maxAgePlusGroup)
      Nsim[s, ] <- step$N
      ssbMat[s, y] <- step$ssb
      fMat[s, y]   <- fval
      cMat[s, y]   <- step$catch
      rMat[s, y]   <- step$N[1]
    }
  }
  triplet <- function(M){
    apply(M, 2, function(v)
      c(Estimate = mean(v, na.rm = TRUE),
        Low      = stats::quantile(v, 0.025, na.rm = TRUE, names = FALSE),
        High     = stats::quantile(v, 0.975, na.rm = TRUE, names = FALSE)))
  }
  proj <- list(SSB = triplet(ssbMat), Fbar = triplet(fMat),
               Recruitment = triplet(rMat), Catch = triplet(cMat))
  for(nm in names(proj)){
    rownames(proj[[nm]]) <- c("Estimate","Low","High")
    colnames(proj[[nm]]) <- as.character(st$year0 + seq_len(nYears))
  }
  rule <- list(Ftarget = Ftarget, Btrigger = Btrigger,
               Forigin = Forigin, Borigin = Borigin,
               Fcap = Fcap, Bcap = Bcap)
  ret <- list(projection = proj, rule = rule, hcr = ruleClosure,
              fbarlabel = substitute(bar(F)[X-Y],
                                     list(X = fit$conf$fbarRange[1],
                                          Y = fit$conf$fbarRange[2])))
  attr(ret, "fit") <- fit
  class(ret) <- "samref_hcr"
  ret
}

##' ICES advice rule (F = 0 below Blim, ramp to Fmsy at MSYBtrigger)
##'
##' Thin wrapper around \code{\link{hcr}}.
##' @param fit a samjr fit.
##' @param Fmsy target F.
##' @param MSYBtrigger SSB above which the target is fully applied.
##' @param Blim below this SSB, F is set to zero.
##' @param ... other arguments to \code{\link{hcr.sam}}.
##' @examples
##' \donttest{
##' data(nscodData); data(nscodConf)
##' fit <- samjr::sam.fit(nscodData, nscodConf,
##'                       samjr::defpar(nscodData, nscodConf), silent = TRUE)
##' fit$conf$stockRecruitmentModelCode <- 2L
##' fitBH <- samjr::runwithout(fit)
##' hh <- icesAdviceRule(fitBH, Fmsy = 0.22, MSYBtrigger = 150000,
##'                      Blim = 100000, nosim = 100, nYears = 5, seed = 1)
##' hh$projection$SSB
##' }
##' @export
icesAdviceRule <- function(fit, ...) UseMethod("icesAdviceRule")

##' @rdname icesAdviceRule
##' @method icesAdviceRule sam
##' @export
icesAdviceRule.sam <- function(fit, Fmsy, MSYBtrigger, Blim, ...){
  hcr.sam(fit, Ftarget = Fmsy, Btrigger = MSYBtrigger,
          Forigin = 0, Borigin = Blim,
          Bcap = Blim, Fcap = 0, ...)
}

##' @method print samref_hcr
##' @export
print.samref_hcr <- function(x, ...){
  cat("samref HCR projection\n")
  cat("Rule: ");  print(unlist(x$rule))
  cat("\nSSB (Estimate, 95% CI):\n");   print(round(x$projection$SSB, 1))
  cat("\nFbar:\n");                      print(round(x$projection$Fbar, 4))
  cat("\nRecruitment:\n");               print(round(x$projection$Recruitment, 1))
  cat("\nCatch:\n");                     print(round(x$projection$Catch, 1))
  invisible(x)
}

## -------- HCR plot methods --------

##' @keywords internal
##' @noRd
.plotHcrSeries <- function(mat, ylab, ribbonCol = grDevices::grey(0.86),
                           lineCol = "black", lwd = 2,
                           hlines = numeric(0), hlineLabels = character(0),
                           ...){
  years  <- as.integer(colnames(mat))
  ymin <- min(0, mat["Low", ], na.rm = TRUE)
  ymax <- max(mat["High", ], hlines, na.rm = TRUE)
  plot(years, mat["Estimate", ], type = "n", xlab = "Year", ylab = ylab,
       ylim = c(ymin, ymax * 1.05), ...)
  polygon(c(years, rev(years)),
          c(mat["Low", ], rev(mat["High", ])),
          col = ribbonCol, border = NA)
  lines(years, mat["Estimate", ], col = lineCol, lwd = lwd)
  ## Stagger hline labels vertically when they're within 5% of each other.
  ord <- order(hlines)
  yLab <- hlines
  if(length(yLab) > 1){
    span <- diff(range(c(0, ymax)))
    for(j in 2:length(ord)){
      if(yLab[ord[j]] - yLab[ord[j - 1]] < 0.06 * span)
        yLab[ord[j]] <- yLab[ord[j - 1]] + 0.06 * span
    }
  }
  for(i in seq_along(hlines)){
    abline(h = hlines[i], lty = 2, col = "grey50")
    if(length(hlineLabels) >= i)
      mtext(hlineLabels[i], side = 4, at = yLab[i], line = 0.4,
            col = "grey30", cex = 0.8, las = 1)
  }
  invisible(NULL)
}

##' Plot the SSB / F-bar / Recruitment / Catch trajectory of a samref HCR
##' projection.
##'
##' Each method draws one year-by-year quantity with the 95\% ribbon and
##' (where relevant) dashed reference lines for the rule parameters.
##' Recruitment uses age \code{minAge}.
##'
##' @param fit a \code{samref_hcr} object as returned by \code{\link{hcr}}.
##' @param ... extra graphical parameters passed to \code{plot}.
##' @return invisible \code{NULL}.
##' @name samref_hcr-plots
NULL

##' @rdname samref_hcr-plots
##' @method ssbplot samref_hcr
##' @export
ssbplot.samref_hcr <- function(fit, ...){
  rule <- fit$rule
  .plotHcrSeries(fit$projection$SSB, ylab = "SSB",
                 hlines      = c(rule$Btrigger, rule$Bcap),
                 hlineLabels = c("Btrigger",   "Bcap/Blim"),
                 ...)
}

##' @rdname samref_hcr-plots
##' @method fbarplot samref_hcr
##' @export
fbarplot.samref_hcr <- function(fit, ...){
  rule <- fit$rule
  .plotHcrSeries(fit$projection$Fbar, ylab = fit$fbarlabel,
                 hlines      = c(rule$Ftarget),
                 hlineLabels = c("Ftarget"),
                 ...)
}

##' @rdname samref_hcr-plots
##' @method recplot samref_hcr
##' @export
recplot.samref_hcr <- function(fit, ...){
  fitOrig <- attr(fit, "fit")
  lab <- if(!is.null(fitOrig)) paste0("Recruits (age ", fitOrig$conf$minAge, ")")
         else "Recruitment"
  .plotHcrSeries(fit$projection$Recruitment, ylab = lab, ...)
}

##' @rdname samref_hcr-plots
##' @method catchplot samref_hcr
##' @export
catchplot.samref_hcr <- function(fit, ...){
  .plotHcrSeries(fit$projection$Catch, ylab = "Catch", ...)
}

##' 2x2 overview plot of an \code{samref_hcr} object: SSB, F-bar,
##' Recruitment, Catch.
##' @param x a \code{samref_hcr} object.
##' @param ... extra arguments passed to the panel plots.
##' @return invisible \code{NULL}.
##' @examples
##' \donttest{
##' data(nscodData); data(nscodConf)
##' fit <- samjr::sam.fit(nscodData, nscodConf,
##'                       samjr::defpar(nscodData, nscodConf), silent = TRUE)
##' fit$conf$stockRecruitmentModelCode <- 2L
##' fitBH <- samjr::runwithout(fit)
##' hh <- icesAdviceRule(fitBH, Fmsy = 0.22, MSYBtrigger = 150000,
##'                      Blim = 100000, nosim = 200, nYears = 5, seed = 1)
##' plot(hh)         # 2x2 overview
##' ssbplot(hh)      # SSB only with rule lines
##' fbarplot(hh)
##' }
##' @method plot samref_hcr
##' @export
plot.samref_hcr <- function(x, ...){
  oldpar <- par(mfrow = c(2, 2), mar = c(4.2, 4.3, 1.5, 4.0))
  on.exit(par(oldpar))
  ssbplot.samref_hcr(x, ...)
  fbarplot.samref_hcr(x, ...)
  recplot.samref_hcr(x, ...)
  catchplot.samref_hcr(x, ...)
  invisible(NULL)
}
