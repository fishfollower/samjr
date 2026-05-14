## Compares samjr retro / mohn / runwithout / leaveout output to a
## STATIC SAM reference (see res.EXP). The reference was generated once
## with stockassessment::retro etc. on the same nscod fit and written
## verbatim into res.EXP using format() with 4 decimals.
##
## Running the test only loads samjr (no live SAM call).

suppressMessages(library(samjr))

cn <- read.ices("cn.dat"); cw <- read.ices("cw.dat"); dw <- read.ices("dw.dat")
lf <- read.ices("lf.dat"); lw <- read.ices("lw.dat"); mo <- read.ices("mo.dat")
nm <- read.ices("nm.dat"); pf <- read.ices("pf.dat"); pm <- read.ices("pm.dat")
sw <- read.ices("sw.dat"); surveys <- read.ices("survey.dat")

dat <- setup.sam.data(surveys = surveys, residual.fleet = cn,
  prop.mature = mo, stock.mean.weight = sw, catch.mean.weight = cw,
  dis.mean.weight = dw, land.mean.weight = lw, prop.f = pf, prop.m = pm,
  natural.mortality = nm, land.frac = lf)
conf <- defcon(dat)
conf$keyLogFsta[1, ] <- c(0, 1, 2, 3, 4, 5)
conf$corFlag        <- 2
conf$keyLogFpar     <- matrix(c(-1, -1, -1, -1, -1, -1,
                                  0,  1,  2,  3,  4, -1,
                                  5,  6,  7,  8, -1, -1), nrow = 3, byrow = TRUE)
conf$keyVarF[1, ]   <- c(0, 1, 1, 1, 1, 1)
conf$keyVarObs      <- matrix(c( 0,  1,  2,  2,  2,  2,
                                  3,  4,  4,  4,  4, -1,
                                  5,  6,  6,  6, -1, -1), nrow = 3, byrow = TRUE)
conf$noScaledYears  <- 13
conf$keyScaledYears <- 1993:2005
conf$keyParScaledYA <- row(matrix(NA, nrow = 13, ncol = 6)) - 1
conf$fbarRange      <- c(2, 4)
par <- defpar(dat, conf)
fit <- sam.fit(dat, conf, par)

ret  <- retro(fit, year = 5)
mret <- mohn(ret)
lout <- leaveout(fit)

con <- file("res.out", open = "w")
cat("# baby logLik         :", format(round(as.numeric(logLik(fit)), 4), nsmall = 4), "\n", file = con)
cat("# baby retro logLik   :",
    paste(format(round(sapply(ret, function(f) as.numeric(logLik(f))), 4), nsmall = 4),
          collapse = " "), "\n", file = con)
cat("# baby retro Mohn     :",
    paste(format(round(mret, 4), nsmall = 4), collapse = " "), "\n", file = con)
cat("# baby leaveout logLik:",
    paste(format(round(sapply(lout, function(f) as.numeric(logLik(f))), 4), nsmall = 4),
          collapse = " "), "\n", file = con)
close(con)
