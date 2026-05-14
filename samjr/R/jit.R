##' Jitter starting values and re-fit a samjr model
##'
##' Re-fits the model \code{nojit} times from random perturbations of the
##' starting parameter list. The flat numeric parameter vector is
##' obtained via \code{unlist(par)}, jittered with independent
##' \code{rnorm(0, sd)} noise, then put back into the original list shape
##' with \code{utils::relist}. Most parameters are on a log scale, so
##' \code{sd} is roughly a coefficient of variation.
##'
##' Mirrors SAM's \code{stockassessment::jit}, simplified for samjr:
##' the loop is always serial (\code{ncores} is accepted for signature
##' compatibility but ignored).
##'
##' @param fit a fitted \code{sam} object as returned by
##' \code{\link{sam.fit}}.
##' @param nojit number of jittered re-fits.
##' @param par initial parameter list to jitter around. Defaults to
##' \code{defpar(fit$data, fit$conf)}.
##' @param sd standard deviation of the Gaussian jitter applied to the
##' flat parameter vector.
##' @param ncores ignored in samjr (always serial).
##' @return a \code{samset} (list of \code{sam} fits) with the original
##' fit attached as \code{attr(., "fit")}.
##' @export
jit <- function(fit, nojit = 10, par = defpar(fit$data, fit$conf),
                sd = 0.25, ncores = 1){
  parv <- unlist(par)
  pars <- lapply(seq_len(nojit),
                 function(i) utils::relist(parv + rnorm(length(parv), sd = sd), par))
  fits <- lapply(pars, function(p) sam.fit(fit$data, fit$conf, p))
  attr(fits, "fit") <- fit
  attr(fits, "jitflag") <- 1
  class(fits) <- "samset"
  fits
}
