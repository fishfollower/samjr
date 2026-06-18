source("../common-fit.R")

## Expected values are SAM 0.12.0 outputs from
##   deterministicReferencepoints(fitSAMrw, c("F=0.3","F=0.000001"),
##                                catchType="landing")
## (RW fit, last 10 years averaging). Captured once.

pr <- perRecruitClosure(fitRW, aveYears = ave10, selYears = selLast,
                        catchType = "landing")
o0  <- pr(log(1e-6), recPars(fitRW))
o03 <- pr(log(0.30), recPars(fitRW))

writeLines(c(
  sprintf("logYPR(F=0.3) = %.6f", o03$logYPR),
  sprintf("logSPR(F=0.3) = %.6f", o03$logSPR),
  sprintf("logSPR(F~0)   = %.6f", o0$logSPR)
), "res.out")
