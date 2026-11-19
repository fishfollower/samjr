source("../common-fit.R")

## SAM 0.12.0 snapshot of CIs for the same RPs on the same RW fit.
## Captured via one-at-a-time
##   deterministicReferencepoints(fitSAMrw, spec, catchType = "landing")
## on the BB RW fit reconstructed from samjr/testmore/nscod/script.R.
##
## SAM differentiates through the full ADREPORTed tape via sdreport;
## samref uses FD over (rec_pars + logF at selYears) with
## Matrix::solve(jointPrecision) for the marginal covariance. The two
## methods agree on the point estimates; CI bounds agree to within ~0.05%
## on this fit. The test checks each relative gap stays below 5%.
SAMci <- list(
  "Fmax"    = c(Low = 0.2715, High = 0.3838),
  "F0.1dYPR"= c(Low = 0.1812, High = 0.2701),
  "F35%SPR" = c(Low = 0.2175, High = 0.2767)
)

rp <- referencepoints(fitRW, c("Max", "0.1dYPR", "0.35SPR"),
                       aveYears = ave10, selYears = selLast,
                       catchType = "landing")

## Only pass/fail is written: the relative gaps themselves sit at the
## rounding level of the 4-digit snapshot and vary across platforms.
tol <- 0.05
out <- character(0)
allOk <- TRUE
for(nm in names(SAMci)){
  for(b in c("Low", "High")){
    s <- SAMci[[nm]][[b]]
    r <- rp$tables$F[nm, b]
    ok <- is.finite(r) && abs(r - s) / abs(s) <= tol
    allOk <- allOk && ok
    out <- c(out, sprintf("%-9s %-4s within %.2f : %s", nm, b, tol, ok))
  }
}
out <- c(out, sprintf("all within %.2f : %s", tol, allOk))
writeLines(out, "res.out")
