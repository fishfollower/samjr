## Compare table outputs of samjr against a snapshot of SAM
## (stockassessment) tables on the canonical nscod fit. The snapshot was
## produced once with stockassessment 0.12.0 and saved as tabSAM.rds; the
## test no longer needs SAM installed at run time. Both packages should
## produce numerically identical Estimate / Low / High columns for the
## standard time-series tables (ssbtable, fbartable, rectable, catchtable,
## tsbtable) and identical raw tables for ntable / faytable / caytable /
## qtable. To regenerate the snapshot run tools/build-tabSAM.R (or
## extractTables(fitSAM) inside an interactive SAM session).

setupAndFit <- function(){
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
  sam.fit(dat, conf, par)
}

extractTables <- function(fit){
  list(
    ssb   = unname(ssbtable(fit)),
    fbar  = unname(fbartable(fit)),
    rec   = unname(rectable(fit)),
    catch = unname(catchtable(fit)),
    tsb   = unname(tsbtable(fit)),
    ntab  = unname(ntable(fit)),
    fay   = unname(faytable(fit)),
    cay   = unname(caytable(fit)),
    qtab  = { q <- qtable(fit); attr(q, "sd") <- NULL; class(q) <- "matrix"; unname(q) },
    logL  = as.numeric(logLik(fit)),
    nobs  = as.integer(nobs(fit)),
    aic   = AIC(fit)
  )
}

## --- 1. samjr ---
suppressMessages(library(samjr))
fitBaby <- setupAndFit()
tabBaby <- extractTables(fitBaby)

## --- 2. SAM stockassessment (loaded from snapshot) ---
tabSAM <- readRDS("tabSAM.rds")

## --- 3. compare ---
tol <- 1e-4
cmp <- function(nm){
  res <- all.equal(tabBaby[[nm]], tabSAM[[nm]],
                   tolerance = tol, check.attributes = FALSE)
  isTRUE(res)
}
cat("logLik\t", cmp("logL"),  "\n", file = "res.out", append = FALSE, sep = "")
cat("nobs\t",   cmp("nobs"),  "\n", file = "res.out", append = TRUE,  sep = "")
cat("AIC\t",    cmp("aic"),   "\n", file = "res.out", append = TRUE,  sep = "")
cat("ssb\t",    cmp("ssb"),   "\n", file = "res.out", append = TRUE,  sep = "")
cat("fbar\t",   cmp("fbar"),  "\n", file = "res.out", append = TRUE,  sep = "")
cat("rec\t",    cmp("rec"),   "\n", file = "res.out", append = TRUE,  sep = "")
cat("catch\t",  cmp("catch"), "\n", file = "res.out", append = TRUE,  sep = "")
cat("tsb\t",    cmp("tsb"),   "\n", file = "res.out", append = TRUE,  sep = "")
cat("ntable\t", cmp("ntab"),  "\n", file = "res.out", append = TRUE,  sep = "")
cat("faytable\t", cmp("fay"), "\n", file = "res.out", append = TRUE,  sep = "")
cat("caytable\t", cmp("cay"), "\n", file = "res.out", append = TRUE,  sep = "")
cat("qtable\t", cmp("qtab"),  "\n", file = "res.out", append = TRUE,  sep = "")
