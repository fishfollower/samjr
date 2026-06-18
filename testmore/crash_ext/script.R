source("../common-fit.R")

## Expected value is the SAM 0.12.0 output:
##   deterministicReferencepoints(fitSAMbh, "Ext",
##                                catchType="landing")$tables$F["Ext","Estimate"]
## (BH fit, last-10 averaging, last-year selectivity).
##
## Note: SAM rejects Crash for Beverton-Holt (its isCompensatory test marks
## BH as non-compensatory). samref's Crash on BH is well-defined
## mathematically but cannot be cross-checked against SAM here.

rp <- referencepoints(fitBH, "Ext",
                       aveYears = ave10, selYears = selLast,
                       catchType = "landing")
writeLines(c(
  sprintf("Fext = %.4f", rp$tables$F["Fext","Estimate"])
), "res.out")
