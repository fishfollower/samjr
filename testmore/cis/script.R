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
## ~0.05% on this fit. Only pass/fail against tol goes to res.out: the Fmax
## bounds come from FD over a re-optimised argmax and move at the 1e-4 level
## between LAPACK builds.
SAMci <- list(
  "Fmax"    = c(Low = 0.2715, High = 0.3838),
  "F0.1dYPR"= c(Low = 0.1812, High = 0.2701),
  "F35%SPR" = c(Low = 0.2175, High = 0.2767)
)

rp <- referencepoints(fitRW, c("Max", "0.1dYPR", "0.35SPR"),
                       aveYears = ave10, selYears = selLast,
                       catchType = "landing")

tol <- 0.01
out <- character(0)
for(nm in names(SAMci)){
  for(b in c("Low", "High")){
    s <- SAMci[[nm]][[b]]
    r <- rp$tables$F[nm, b]
    d <- abs(r - s) / abs(s)
    out <- c(out, sprintf("%-9s %-4s within %g%%: %s", nm, b, 100 * tol, d <= tol))
  }
}
writeLines(out, "res.out")
