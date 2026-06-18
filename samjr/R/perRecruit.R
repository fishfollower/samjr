## Numeric helpers wrapping perRecruitClosure.

##' Per-recruit / equilibrium values at a sequence of F values
##'
##' Convenience wrapper that evaluates \code{\link{perRecruitClosure}} at a
##' vector of Fbar values and returns a data frame, mirroring SAM's
##' \code{.perRecruitR}.
##'
##' @param fit a samjr \code{sam} fit.
##' @param Fsequence numeric F-bar grid (natural scale).
##' @param aveYears,selYears,catchType,customSel passed to
##'   \code{\link{perRecruitClosure}}.
##' @return a data frame with one row per F, columns matching the
##'   \code{perRec} list entries.
##' @examples
##' \donttest{
##' data(nscodData); data(nscodConf)
##' fit <- samjr::sam.fit(nscodData, nscodConf,
##'                       samjr::defpar(nscodData, nscodConf), silent = TRUE)
##' perRecruitTable(fit, Fsequence = seq(0, 1, 0.2), catchType = "landing")
##' }
##' @export
perRecruitTable <- function(fit,
                            Fsequence = seq(0, 2, length = 50),
                            aveYears  = NULL,
                            selYears  = NULL,
                            catchType = c("catch", "landing", "discard"),
                            customSel = NULL){
  catchType <- match.arg(catchType)
  pr <- perRecruitClosure(fit, aveYears = aveYears, selYears = selYears,
                          catchType = catchType, customSel = customSel)
  rp <- recPars(fit)
  fSafe <- pmax(Fsequence, 1e-12)
  rows <- lapply(log(fSafe), function(lf) pr(lf, rp))
  out <- do.call(rbind, lapply(rows, function(r) unlist(r)))
  ret <- as.data.frame(out)
  ret$Fbar <- Fsequence
  ret
}
