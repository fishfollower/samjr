srCodes <- c(0, 1, 2, 62, 102)
nPeel <- 15
resFile <- "recpred.out"
source("script.R")

## One page per prediction horizon. Each page has one panel per SR code with
## full run R and the predictions, plus a panel with the log prediction errors.
## Axes are shared across pages so the horizons can be compared directly.
modeCol <- c(sr0 = "#1f5fa8", sr1 = "#eda100", sr2 = "#e87ba4", sr62 = "#4fb0e0", sr102 = "#8c564b")
modeName <- c(sr0 = "SR code 0 (random walk)", sr1 = "SR code 1 (Ricker)",
              sr2 = "SR code 2 (Beverton-Holt)", sr62 = "SR code 62 (AR(1))",
              sr102 = "SR code 102 (AR(1) with WSC in mean)")
modeFormula <- list(
  sr0 = quote(log~R[y] == log~R[y-1] + epsilon[y]),
  sr1 = quote(log~R[y] == a + log~SSB[y-1] - e^b*SSB[y-1] + epsilon[y]),
  sr2 = quote(log~R[y] == a + log~SSB[y-1] - log(1 + e^b*SSB[y-1]) + epsilon[y]),
  sr62 = quote(log~R[y] == mu + rho*(log~R[y-1] - mu) + epsilon[y]),
  sr102 = quote(log~R[y] == (1 - rho)*(mu + alpha~WSC[y-1]) + rho~log~R[y-1] + epsilon[y])
)
models <- names(recPred)
modeCol <- modeCol[models]
modeName <- modeName[models]
errFormula <- quote(epsilon[y] %~% N(0, sigma[R]^2))
stepLabel <- paste0(toupper(substring(stepName, 1, 1)), substring(stepName, 2), "-year-ahead")
toBillion <- function(x) x / 1e6
allYears <- range(unlist(lapply(recPred, function(r) r$year)))
showYears <- (allYears[1] - 1):allYears[2]
fullTab <- lapply(fits, function(f) toBillion(rectable(f)[as.character(showYears), ]))
yLim <- range(unlist(lapply(names(recPred), function(m)
  c(toBillion(recPred[[m]]$low), toBillion(recPred[[m]]$high), fullTab[[m]][, c("Low", "High")])))) * c(1, 3)
yTicks <- c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500)
yTicks <- yTicks[yTicks >= yLim[1] & yTicks <= yLim[2]]
errLim <- range(unlist(lapply(recPred, function(r) r$logErr)), 0)
errLim <- errLim + c(-0.1, 0.9) * diff(errLim)
off <- setNames(seq(-0.32, 0.32, length.out = length(models)), models)

plotStep <- function(h){
  pred <- lapply(recPred, function(r) r[r$step == h, ])
  layout(matrix(1:6, 3, byrow = TRUE))
  par(mar = c(4, 4.5, 3.4, 1), mgp = c(2.6, 0.7, 0), cex = 0.9, oma = c(0, 0, 1.8, 0))
  for(m in names(pred)){
    r <- pred[[m]]
    f <- fullTab[[m]]
    plot(NA, xlim = range(showYears), ylim = yLim, log = "y", las = 1, yaxt = "n",
         xlab = "Year", ylab = "Recruitment (billions)")
    axis(2, at = yTicks, labels = format(yTicks, drop0trailing = TRUE, trim = TRUE), las = 1)
    abline(h = yTicks, col = "grey92")
    polygon(c(showYears, rev(showYears)), c(f[, "Low"], rev(f[, "High"])), col = "grey85", border = NA)
    lines(showYears, f[, "Estimate"], lwd = 2, col = "grey30")
    arrows(r$year, toBillion(r$low), r$year, toBillion(r$high),
           angle = 90, code = 3, length = 0.025, col = modeCol[m], lwd = 1.2)
    points(r$year, toBillion(r$pred), pch = 21, bg = modeCol[m], col = "white", cex = 1.3)
    legend("topleft", c("Full run estimate (+/- 2 sd)", paste(stepLabel[h], "prediction (95%)")),
           lwd = c(2, 1.2), col = c("grey30", modeCol[m]), pch = c(NA, 21), pt.bg = c(NA, modeCol[m]),
           fill = c("grey85", NA), border = NA, bty = "n", cex = 0.8)
    mtext(modeName[m], side = 3, line = 1.7, adj = 0, font = 2)
    mtext(bquote(.(modeFormula[[m]]) * ",  " ~ .(errFormula)), side = 3, line = 0.3, adj = 0, cex = 0.8)
  }
  plot(NA, xlim = allYears + c(-0.5, 0.5), ylim = errLim, las = 1,
       xlab = "Predicted year", ylab = "log(predicted / full run)")
  grid(nx = NA, ny = NULL, col = "grey92", lty = 1)
  abline(h = 0, col = "grey40")
  for(m in names(pred)){
    r <- pred[[m]]
    segments(r$year + off[m], 0, r$year + off[m], r$logErr, col = modeCol[m], lwd = 1.6)
    points(r$year + off[m], r$logErr, pch = 21, bg = modeCol[m], col = "white", cex = 1)
  }
  legend("top", sprintf("%s: RMSE %.2f", modeName, rmse[h, ]), pch = 21, pt.bg = modeCol,
         col = "white", pt.cex = 1.2, bty = "n", cex = 0.85)
  mtext("Prediction error", side = 3, line = 1.7, adj = 0, font = 2)
  mtext(sprintf("%s recruitment predictions (%d-%d)", stepLabel[h], predYears[[h]][1], predYears[[h]][2]),
        side = 3, outer = TRUE, line = 0.3, font = 2, cex = 1.1)
}

pdf("recpred.pdf", width = 11, height = 11)
for(h in seq_len(nStep)) plotStep(h)
dev.off()
