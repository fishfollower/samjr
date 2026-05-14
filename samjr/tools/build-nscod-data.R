## Rebuild data/nscodData.rda, data/nscodConf.rda, data/nscodParameters.rda
## from the ICES files in SAM/testmore/nscod (run with the samjr package
## already installed).
library(samjr)

scriptDir <- "../testmore/nscod"
if(!dir.exists(scriptDir))
  stop("expected data dir not found: ", normalizePath(scriptDir, mustWork = FALSE))
old <- setwd(scriptDir); on.exit(setwd(old))

cn <- read.ices("cn.dat"); cw <- read.ices("cw.dat"); dw <- read.ices("dw.dat")
lf <- read.ices("lf.dat"); lw <- read.ices("lw.dat"); mo <- read.ices("mo.dat")
nm <- read.ices("nm.dat"); pf <- read.ices("pf.dat"); pm <- read.ices("pm.dat")
sw <- read.ices("sw.dat"); surveys <- read.ices("survey.dat")

nscodData <- setup.sam.data(surveys = surveys, residual.fleet = cn,
  prop.mature = mo, stock.mean.weight = sw, catch.mean.weight = cw,
  dis.mean.weight = dw, land.mean.weight = lw, prop.f = pf, prop.m = pm,
  natural.mortality = nm, land.frac = lf)

nscodConf <- defcon(nscodData)
nscodConf$keyLogFsta[1, ] <- c(0, 1, 2, 3, 4, 5)
nscodConf$corFlag        <- 2
nscodConf$keyLogFpar     <- matrix(
  c(-1, -1, -1, -1, -1, -1,
     0,  1,  2,  3,  4, -1,
     5,  6,  7,  8, -1, -1), nrow = 3, byrow = TRUE)
nscodConf$keyVarF[1, ]   <- c(0, 1, 1, 1, 1, 1)
nscodConf$keyVarObs      <- matrix(
  c( 0,  1,  2,  2,  2,  2,
     3,  4,  4,  4,  4, -1,
     5,  6,  6,  6, -1, -1), nrow = 3, byrow = TRUE)
nscodConf$noScaledYears  <- 13
nscodConf$keyScaledYears <- 1993:2005
nscodConf$keyParScaledYA <- row(matrix(NA, nrow = 13, ncol = 6)) - 1
nscodConf$fbarRange      <- c(2, 4)

nscodParameters <- defpar(nscodData, nscodConf)

setwd(old)
dataDir <- file.path(getwd(), "data")
dir.create(dataDir, showWarnings = FALSE)
save(nscodData,       file = file.path(dataDir, "nscodData.rda"),       compress = "xz")
save(nscodConf,       file = file.path(dataDir, "nscodConf.rda"),       compress = "xz")
save(nscodParameters, file = file.path(dataDir, "nscodParameters.rda"), compress = "xz")
cat("rebuilt:", list.files(dataDir, full.names = FALSE), "\n")
