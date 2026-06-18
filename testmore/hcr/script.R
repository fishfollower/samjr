source("../common-fit.R")

## Expected values are exact algebraic answers from the trapezoidal rule
## with Ftarget=0.22, Btrigger=150000, Borigin=Bcap=Blim=100000,
## Forigin=Fcap=0. These match SAM 0.12.0's C_hcrR call for the same
## parameters (verified one-off).

hh <- icesAdviceRule(fitBH, Fmsy = 0.22, MSYBtrigger = 150000, Blim = 100000,
                     nosim = 10, nYears = 1, seed = 1,
                     aveYears = ave10, selYears = selLast)
writeLines(c(
  sprintf("F(SSB= 50000) = %.4f", hh$hcr( 50000)),
  sprintf("F(SSB=100000) = %.4f", hh$hcr(100000)),
  sprintf("F(SSB=120000) = %.4f", hh$hcr(120000)),
  sprintf("F(SSB=150000) = %.4f", hh$hcr(150000)),
  sprintf("F(SSB=200000) = %.4f", hh$hcr(200000))
), "res.out")
