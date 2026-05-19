## Compares samjr residplot p-values to a STATIC SAM reference
## (sam-residplot.dat). The reference was generated once with
## stockassessment::residplot on the same nscod fit.
##
## residplot's tests (bias t-test, variance chi-square, Ljung-Box in
## age and time, mean-variance F-test, Shapiro-Wilk) are deterministic
## functions of the OSA residuals. samjr and SAM compute OSA residuals
## via RTMB and TMB respectively (both via oneStepPredict) so the
## resulting p-values should agree to roughly Hessian precision.
##
## Pass criterion: correlation > 0.9999 AND max|diff| < 1e-5,
## matched by p-value name.

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
conf$noScaledYears  <- 13; conf$keyScaledYears <- 1993:2005
conf$keyParScaledYA <- row(matrix(NA, nrow = 13, ncol = 6)) - 1
conf$fbarRange      <- c(2, 4)
par <- defpar(dat, conf)
fit <- sam.fit(dat, conf, par)

pdf(tempfile(fileext = ".pdf"))
rp <- residplot(fit)
dev.off()

flat <- function(x, tag){
  if(is.matrix(x)){
    g <- expand.grid(r = seq_len(nrow(x)), c = seq_len(ncol(x)),
                     KEEP.OUT.ATTRS = FALSE)
    data.frame(name = sprintf("%s.%d.%d", tag, g$r, g$c),
               value = as.vector(x))
  }else{
    data.frame(name = sprintf("%s.%d", tag, seq_along(x)),
               value = as.numeric(x))
  }
}
baby <- rbind(
  flat(rp$bias,             "bias"),
  flat(rp$variance,         "variance"),
  flat(rp$correlation.age,  "corage"),
  flat(rp$meanvar,          "meanvar"),
  flat(rp$correlation.time, "cortime"),
  flat(rp$normality,        "normality"))

sam <- read.table("sam-residplot.dat", header = FALSE,
                  col.names = c("name", "value"))
mer <- merge(baby, sam, by = "name", suffixes = c(".B", ".S"))
ok    <- is.finite(mer$value.B) & is.finite(mer$value.S)
rho   <- cor(mer$value.B[ok], mer$value.S[ok])
mxabs <- max(abs(mer$value.B[ok] - mer$value.S[ok]))

con <- file("res.out", open = "w")
cat("# residplot correlation > 0.9999: ", rho   > 0.9999, "\n", sep = "", file = con)
cat("# residplot max|diff|   < 1e-5  : ", mxabs < 1e-5,   "\n", sep = "", file = con)
cat("# rows compared                 : ", sum(ok),        "\n", sep = "", file = con)
close(con)
