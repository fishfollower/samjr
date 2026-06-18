source("../common-fit.R")

## Expected values are SAM 0.12.0 outputs, one RP per call (to dodge SAM's
## multi-RP row-labelling bug):
##   F35%SPR : deterministicReferencepoints(fitSAMbh, "0.35SPR", catchType="landing")
##   F40%SPR : deterministicReferencepoints(fitSAMbh, "0.4SPR",  catchType="landing")
##   F20%B0  : deterministicReferencepoints(fitSAMbh, "0.2B0",   catchType="landing")
## All on the BH fit with default averaging (last 10) and selectivity (last year).

rp <- referencepoints(fitBH, c("0.35SPR", "0.4SPR", "0.2B0"),
                       aveYears = ave10, selYears = selLast,
                       catchType = "landing")
writeLines(c(
  sprintf("F35%%SPR = %.4f", rp$tables$F["F35%SPR","Estimate"]),
  sprintf("F40%%SPR = %.4f", rp$tables$F["F40%SPR","Estimate"]),
  sprintf("F20%%B0  = %.4f", rp$tables$F["F20%B0", "Estimate"])
), "res.out")
