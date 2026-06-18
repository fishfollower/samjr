suppressPackageStartupMessages(library(samjr))

## A user-supplied logSRfun/logSRinit must reproduce the built-in SR it
## mirrors. Case study: NEA saithe (1960-2020), where both Beverton-Holt and
## Ricker converge to fully identified SR parameters - so the equivalence is
## checked on the converged fit (nll AND both SR parameters), not just at the
## initial point.
sf   <- readRDS("saithe.rds")
data <- sf$data
conf <- sf$conf

## Fit the built-in SR (code) and the custom logSRfun equivalent; return both.
fitEq <- function(code, srfun, init, builtinName, builtinInit = NULL){
  cc <- conf; cc$stockRecruitmentModelCode <- code
  parB <- defpar(data, cc)
  if(!is.null(builtinInit)) parB[[builtinName]] <- builtinInit
  fitB <- sam.fit(data, cc, parB, silent = TRUE)
  fitC <- sam.fit(data, cc, defpar(data, cc), silent = TRUE,
                  logSRfun = srfun, logSRinit = init)
  list(B = fitB, C = fitC, name = builtinName)
}

verdicts <- function(fit){
  bP <- fit$B$pl[[fit$name]]; cP <- fit$C$pl$logSRpar
  c(conv = fit$B$opt$convergence == 0 && fit$C$opt$convergence == 0,
    nll  = abs(fit$B$opt$objective - fit$C$opt$objective) < 1e-4,
    par  = abs(bP[1] - cP[1]) < 1e-3 && abs(bP[2] - cP[2]) < 1e-3)
}

## Beverton-Holt: default init converges. Ricker: needs a sensible b init
## (b ~ -log(typical SSB)) for both built-in and custom.
bh <- verdicts(fitEq(2L, function(S, pv) pv[1] + log(S) - log(1 + exp(pv[2]) * S),
                     c(1, 1), "bhpar"))
rk <- verdicts(fitEq(1L, function(S, pv) pv[1] + log(S) - exp(pv[2]) * S,
                     c(1, -12), "rickerpar", builtinInit = c(1, -12)))

writeLines(c(
  sprintf("BH     converged:  %s", isTRUE(bh["conv"])),
  sprintf("BH     nll match:  %s", isTRUE(bh["nll"])),
  sprintf("BH     a,b match:  %s", isTRUE(bh["par"])),
  sprintf("Ricker converged:  %s", isTRUE(rk["conv"])),
  sprintf("Ricker nll match:  %s", isTRUE(rk["nll"])),
  sprintf("Ricker a,b match:  %s", isTRUE(rk["par"]))
), "res.out")
