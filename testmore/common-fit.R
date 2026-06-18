## Shared boilerplate: build the nscod fit (RW recruitment) from the
## bundled samjr data + the conf used by samjr/testmore/nscod.
##
## Source from each testmore/<dir>/script.R via source("../common-fit.R").

suppressPackageStartupMessages({
  library(samjr)
})

data(nscodData)
conf <- defcon(nscodData)
conf$keyLogFsta[1, ] <- c(0, 1, 2, 3, 4, 5)
conf$corFlag <- 2
conf$keyLogFpar <- matrix(
  c(-1, -1, -1, -1, -1, -1,
     0,  1,  2,  3,  4, -1,
     5,  6,  7,  8, -1, -1), nrow = 3, byrow = TRUE)
conf$keyVarF[1, ] <- c(0, 1, 1, 1, 1, 1)
conf$keyVarObs <- matrix(
  c( 0,  1,  2,  2,  2,  2,
     3,  4,  4,  4,  4, -1,
     5,  6,  6,  6, -1, -1), nrow = 3, byrow = TRUE)
conf$noScaledYears <- 13
conf$keyScaledYears <- 1993:2005
conf$keyParScaledYA <- row(matrix(NA, nrow = 13, ncol = 6)) - 1
conf$fbarRange <- c(2, 4)

confRW <- conf
parRW  <- defpar(nscodData, confRW)
fitRW  <- sam.fit(nscodData, confRW, parRW, silent = TRUE)

confBH <- conf
confBH$stockRecruitmentModelCode <- 2
parBH  <- defpar(nscodData, confBH)
fitBH  <- sam.fit(nscodData, confBH, parBH, silent = TRUE)

nY <- length(nscodData$years)
aveYears <- nscodData$years[(nY - 14):nY]    # 15-year window
ave10    <- nscodData$years[(nY -  9):nY]    # SAM default (last 10)
selLast  <- nscodData$years[nY]              # SAM default (last year)
