## Compares samjr procres residuals to a STATIC SAM reference
## (sam-procres.dat). The reference was generated once with
## stockassessment::procres(fit, seed=123456) on the same nscod fit.
##
## samjr and SAM agree on the joint posterior of the standardized
## state-equation innovations to within Hessian-precision (~ 1e-8
## absolute), and the MVN sampling step uses the same algorithm
## (chol(Sigma) + rnorm) and same seed in both packages. With the
## resN/resF state vectors emitted in the same shape and order, the
## seeded sample is therefore bit-close.
##
## Pass criterion: correlation > 0.99999 AND max|diff| < 1e-5,
## matched by year/fleet/age.

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
res <- procres(fit, seed = 123456)

baby <- data.frame(year = res$year, fleet = res$fleet, age = res$age,
                   residual = as.numeric(res$residual))
sam  <- read.table("sam-procres.dat", header = FALSE,
                   col.names = c("year", "fleet", "age", "residual"))
mer  <- merge(baby, sam, by = c("year", "fleet", "age"),
              suffixes = c(".B", ".S"))
ok    <- is.finite(mer$residual.B) & is.finite(mer$residual.S)
rho   <- cor(mer$residual.B[ok], mer$residual.S[ok])
mxabs <- max(abs(mer$residual.B[ok] - mer$residual.S[ok]))

con <- file("res.out", open = "w")
cat("# procres correlation > 0.99999: ", rho   > 0.99999, "\n", sep = "", file = con)
cat("# procres max|diff|  < 1e-5    : ", mxabs < 1e-5,    "\n", sep = "", file = con)
cat("# rows compared                : ", sum(ok),         "\n", sep = "", file = con)
close(con)
