## Compares samjr OSA residuals to a STATIC SAM reference (see res.EXP).
## The reference was generated once with stockassessment::residuals on the
## same nscod fit and written verbatim into res.EXP.
##
## Pass criterion: correlation between the two residual vectors > 0.99
## (matched by year/fleet/age). Running the test only loads samjr.

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
res <- residuals(fit)

baby <- data.frame(year = res$year, fleet = res$fleet, age = res$age,
                   residual = as.numeric(res$residual))
sam  <- read.table("sam-residuals.dat", header = FALSE,
                   col.names = c("year", "fleet", "age", "residual"))
mer  <- merge(baby, sam, by = c("year", "fleet", "age"),
              suffixes = c(".B", ".S"))
ok   <- is.finite(mer$residual.B) & is.finite(mer$residual.S)
rho  <- cor(mer$residual.B[ok], mer$residual.S[ok])

con <- file("res.out", open = "w")
cat("# residual correlation > 0.99:", rho > 0.99, "\n", file = con)
cat("# baby residual correlation : ", format(round(rho, 4), nsmall = 4),
    "\n", file = con)
close(con)
