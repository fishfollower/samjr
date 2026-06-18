source("../common-fit.R")

## SAM 0.12.0 snapshot of CIs for the same RPs on the same RW fit.
## Captured via one-at-a-time
##   deterministicReferencepoints(fitSAMrw, spec, catchType = "landing")
## on the BB RW fit reconstructed from samjr/testmore/nscod/script.R.
##
## SAM differentiates through the full ADREPORTed tape via sdreport;
## samref uses FD over (rec_pars + logF at selYears) with
## Matrix::solve(jointPrecision) for the marginal covariance. The two
## methods agree on the point estimates byte-for-byte; CI bounds agree to
## ~2% on this fit. The test checks max relative gap stays below 5%.
SAMci <- list(
  "Fmax"    = c(Low = 0.2715, High = 0.3838),
  "F0.1dYPR"= c(Low = 0.1812, High = 0.2701),
  "F35%SPR" = c(Low = 0.2175, High = 0.2767)
)

rp <- referencepoints(fitRW, c("Max", "0.1dYPR", "0.35SPR"),
                       aveYears = ave10, selYears = selLast,
                       catchType = "landing")

tol <- 0.05
out <- character(0)
maxDiff <- 0
for(nm in names(SAMci)){
  for(b in c("Low", "High")){
    s <- SAMci[[nm]][[b]]
    r <- rp$tables$F[nm, b]
    d <- abs(r - s) / abs(s)
    maxDiff <- max(maxDiff, d)
    out <- c(out, sprintf("%-9s %-4s diff = %.4f", nm, b, d))
  }
}
out <- c(out, sprintf("max diff <= %.2f : %s", tol, maxDiff <= tol))
writeLines(out, "res.out")
