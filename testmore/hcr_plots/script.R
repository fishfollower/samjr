source("../common-fit.R")

## Smoke test: the four HCR plot methods + the overview must all run without
## error on a small icesAdviceRule projection.
hh <- icesAdviceRule(fitBH, Fmsy = 0.22, MSYBtrigger = 150000, Blim = 100000,
                     nosim = 100, nYears = 4, seed = 1,
                     aveYears = ave10, selYears = selLast)

png(tempfile(fileext = ".png"), width = 600, height = 400)
ok1 <- !inherits(try(ssbplot(hh),  silent = TRUE), "try-error"); dev.off()
png(tempfile(fileext = ".png"), width = 600, height = 400)
ok2 <- !inherits(try(fbarplot(hh), silent = TRUE), "try-error"); dev.off()
png(tempfile(fileext = ".png"), width = 600, height = 400)
ok3 <- !inherits(try(recplot(hh),  silent = TRUE), "try-error"); dev.off()
png(tempfile(fileext = ".png"), width = 600, height = 400)
ok4 <- !inherits(try(catchplot(hh),silent = TRUE), "try-error"); dev.off()
png(tempfile(fileext = ".png"), width = 800, height = 600)
ok5 <- !inherits(try(plot(hh),     silent = TRUE), "try-error"); dev.off()

writeLines(c(
  sprintf("ssbplot.samref_hcr:   %s", ok1),
  sprintf("fbarplot.samref_hcr:  %s", ok2),
  sprintf("recplot.samref_hcr:   %s", ok3),
  sprintf("catchplot.samref_hcr: %s", ok4),
  sprintf("plot.samref_hcr:      %s", ok5)
), "res.out")
