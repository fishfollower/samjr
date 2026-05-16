##' Conditionally simulate observations from a samjr fit
##'
##' Simulates new observation vectors from a fitted samjr model holding
##' every latent state (\code{logN}, \code{logF}, the missing-observation
##' fills, and the biology processes \code{logSW}/\code{logCW}/
##' \code{logitMO}/\code{logNM}) fixed at its estimated value. Only the
##' observation distributions are redrawn, so RTMB never has to simulate
##' the state-space recurrences.
##'
##' @param object a fitted \code{sam} object from \code{\link{sam.fit}}.
##' @param nsim number of simulated data sets.
##' @param seed optional random seed.
##' @param full.data if \code{TRUE} (default), each element is a full
##' data list with simulated \code{logobs} replacing the original; if
##' \code{FALSE}, only the raw RTMB simulate output is returned.
##' @param ... unused.
##' @return a list of length \code{nsim}.
##' @method simulate sam
##' @importFrom stats simulate
##' @export
simulate.sam <- function(object, nsim = 1, seed = NULL, full.data = TRUE, ...){
  if(!is.null(seed)) set.seed(seed)
  data2 <- object$data
  if(length(object$pl$missing) > 0){
    naIdx <- which(is.na(data2$logobs))
    data2$logobs[naIdx] <- object$pl$missing
  }
  dat <- toBabyDat(data2, object$conf)
  f   <- makeBabyLikelihood(dat)
  parameters <- object$pl
  parameters$missing <- numeric(0)
  objnew <- MakeADFun(f, parameters, random = NULL, silent = TRUE)
  est <- objnew$par
  replicate(nsim, {
    sval <- objnew$simulate(est)
    if(full.data){
      out <- c(data2[names(data2) != "logobs"], sval["logobs"])
      attr(out, "fleetNames") <- attr(object$data, "fleetNames")
      out
    }else{
      sval
    }
  }, simplify = FALSE)
}

##' Conditional simulation study for a samjr fit
##'
##' Draws \code{nsim} simulated data sets from a fitted model with
##' \code{\link{simulate.sam}} (latent states held at their estimates)
##' and re-fits the model to each. Returns a \code{samset} with the
##' original fit attached as the reference.
##'
##' @param fit a fitted \code{sam} object.
##' @param nsim number of simulations.
##' @param ncores number of parallel workers; defaults to
##' \code{parallel::detectCores()}.
##' @return a list of length \code{nsim} with class \code{samset}.
##' @importFrom parallel detectCores makeCluster stopCluster clusterEvalQ clusterExport parLapply
##' @export
simstudy <- function(fit, nsim, ncores = parallel::detectCores()){
  simdata <- simulate(fit, nsim = nsim, full.data = TRUE)
  refit <- function(x){
    par0 <- defpar(x, fit$conf)
    for(nm in names(par0)){
      if(nm %in% names(fit$pl) &&
         length(par0[[nm]]) == length(fit$pl[[nm]])){
        par0[[nm]][] <- fit$pl[[nm]]
      }
    }
    sam.fit(x, fit$conf, par0, map = fit$map, silent = TRUE)
  }
  if(ncores > 1){
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    libVer <- dirname(path.package("samjr"))
    parallel::clusterExport(cl, varlist = "libVer", envir = environment())
    parallel::clusterEvalQ(cl, library(samjr, lib.loc = libVer))
    runs <- parallel::parLapply(cl, simdata, refit)
  }else{
    runs <- lapply(simdata, refit)
  }
  attr(runs, "fit") <- fit
  class(runs) <- "samset"
  runs
}
