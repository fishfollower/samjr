source("../common-fit.R")

## Expected value is the SAM 0.12.0 output:
##   deterministicReferencepoints(fitSAMbh, "MSY",
##                                catchType="landing")$tables$F["MSY","Estimate"]
## (BH fit, last-10 averaging, last-year selectivity).
##
## Note: SAM's MSYRange in 0.12.0 has a known bug (Lower collapses to Fmsy,
## Upper drifts to Y=0). samref locates the actual 0.95*MSY brackets but
## those are NOT a valid SAM comparison on this fit, so this test verifies
## Fmsy only against SAM.

rp <- referencepoints(fitBH, "MSY",
                       aveYears = ave10, selYears = selLast,
                       catchType = "landing")
writeLines(c(
  sprintf("Fmsy = %.4f", rp$tables$F["Fmsy", "Estimate"])
), "res.out")
