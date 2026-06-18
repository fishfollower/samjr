source("../common-fit.R")

## Expected values are SAM 0.12.0 outputs:
##   deterministicReferencepoints(fitSAMrw, "Max",     catchType="landing")$tables$F["Max",    "Estimate"]
##   deterministicReferencepoints(fitSAMrw, "0.1dYPR", catchType="landing")$tables$F["0.1dYPR","Estimate"]
##   deterministicReferencepoints(fitSAMrw, "0.35SPR", catchType="landing")$tables$F["0.35SPR","Estimate"]
## Captured once with last-10 averaging + last-year selectivity.

rp <- referencepoints(fitRW, c("Max", "0.1dYPR", "0.35SPR"),
                       aveYears = ave10, selYears = selLast,
                       catchType = "landing")
writeLines(c(
  sprintf("Fmax    = %.4f", rp$tables$F["Fmax",    "Estimate"]),
  sprintf("F0.1    = %.4f", rp$tables$F["F0.1dYPR","Estimate"]),
  sprintf("F35%%SPR = %.4f", rp$tables$F["F35%SPR","Estimate"])
), "res.out")
