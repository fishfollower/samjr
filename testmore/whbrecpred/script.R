library(samjr)

## Defaults for the test run; plotscript.R overrides these before sourcing.
if(!exists("srCodes")) srCodes <- c(0, 62, 102)
if(!exists("nPeel")) nPeel <- 10
if(!exists("resFile")) resFile <- "res.out"
nStep <- 2
stepName <- c("one", "two")
nSim <- 1000
seedBase <- 2026
wscScale <- 1e6

dat <- read.data.files(".")
dat$wsc <- dat$wsc * wscScale
conf <- loadConf(dat, "model.cfg")

fitCode <- function(code){
  cat("Fitting SR code", code, "\n")
  cf <- conf
  cf$stockRecruitmentModelCode <- code
  sam.fit(dat, cf, defpar(dat, cf), silent = TRUE)
}

## Forecast noise is drawn year by year, so the one-step draws are unaffected
## by forecasting further ahead.
predictRec <- function(f, seed){
  lastYear <- max(f$data$years)
  set.seed(seed)
  fc <- forecast(f, fscale = rep(1, nStep + 1), year.base = lastYear, processNoiseF = FALSE,
                 useRproc = TRUE, nosim = nSim)
  t(sapply(seq_len(nStep), function(h){
    q <- quantile(fc[[h + 1]]$rec, c(0.5, 0.025, 0.975), names = FALSE)
    c(step = h, year = lastYear + h, pred = q[1], low = q[2], high = q[3])
  }))
}

fits <- lapply(srCodes, fitCode)
names(fits) <- paste0("sr", srCodes)

recPred <- lapply(names(fits), function(m){
  cat("Running", nPeel, "peels for", m, "\n")
  peels <- retro(fits[[m]], year = nPeel)
  res <- as.data.frame(do.call(rbind, lapply(seq_len(nPeel), function(k) predictRec(peels[[k]], seedBase + k))))
  res <- res[res$year %in% fits[[m]]$data$years, ]
  res <- res[order(res$step, res$year), ]
  full <- rectable(fits[[m]])[as.character(res$year), , drop = FALSE]
  res$full <- full[, "Estimate"]
  res$fullLow <- full[, "Low"]
  res$fullHigh <- full[, "High"]
  res$logErr <- log(res$pred) - log(res$full)
  res
})
names(recPred) <- names(fits)

predYears <- tapply(recPred[[1]]$year, recPred[[1]]$step, range, simplify = FALSE)
rmse <- sapply(recPred, function(r) tapply(r$logErr, r$step, function(e) sqrt(mean(e^2))))
out <- unlist(lapply(seq_len(nStep), function(h)
  sprintf("SR code %-3d RMSE log(R) %s-year-ahead %d-%d: %.3f",
          srCodes, stepName[h], predYears[[h]][1], predYears[[h]][2], rmse[h, ])))
writeLines(out, resFile)
