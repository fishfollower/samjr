##' Internal: growth and fecundity grids implied by the configuration.
##'
##' Length at age is log-normal about a von Bertalanffy mean, with the
##' parameters held in \code{conf} rather than estimated. Ages are taken at
##' mid-year, so an animal recorded at age \code{a} has been growing for
##' \code{a + 0.5} years on average.
##'
##' Reproductive output is driven by length, \eqn{\phi(\ell) =
##' (\ell/\ell_{ref})^{\psi}}. Note that \eqn{\psi} here is a \emph{length}
##' exponent: with \eqn{W \propto \ell^3}, a weight exponent of 1.5 in the
##' age-based close-kin model corresponds to 4.5 here.
##'
##' An animal is assumed to hold its length-at-age quantile for life, which is
##' what lets a parent's length be retrojected and what couples an unobserved
##' parent's fecundity across two birth years.
##' @keywords internal
##' @noRd
ckmrlGrid <- function(set, ages){
  aMid  <- ages + 0.5
  Lmean <- set$Linf * (1 - exp(-set$k * (aMid - set$t0)))
  list(Lmean = Lmean, sdL = set$sdLogLength)
}

##' Internal: CKMR-by-length settings, with defaults for configurations
##' predating the feature.
##' @keywords internal
##' @noRd
ckmrlSettings <- function(conf){
  g <- function(nm, def) if(is.null(conf[[nm]])) def else conf[[nm]][1]
  list(usePOPL = as.integer(g("usePOPL", 0L)),
       useHSPL = as.integer(g("useHSPL", 0L)),
       psi     = as.numeric(g("ckmrlPsi", 4.5)),
       estPsi  = as.integer(g("ckmrlEstimatePsi", 0L)),
       scale   = as.numeric(g("ckmrlScale", 1)),
       omega   = as.numeric(g("ckmrlOmega", 1)),
       estOmega = as.integer(g("ckmrlEstimateOmega", 0L)),
       Linf    = as.numeric(g("ckmrlLinf", 41.7)),
       k       = as.numeric(g("ckmrlK", 0.35)),
       t0      = as.numeric(g("ckmrlT0", 0)),
       sdLogLength = as.numeric(g("ckmrlSdLogLength", 0.1)),
       refLength   = as.numeric(g("ckmrlRefLength", 50)))
}

##' Internal: validate a length-based CKMR pair table.
##' @keywords internal
##' @noRd
ckmrlCheck <- function(x){
  need <- c("year1", "len1", "year2", "len2", "nComp")
  if(!is.data.frame(x)) stop("'ckmrl' must be a data frame (see ckmrlData)")
  if(!all(need %in% names(x)))
    stop("'ckmrl' must have columns ", paste(need, collapse = ", "))
  if(is.null(x$nPOP)) x$nPOP <- 0
  if(is.null(x$nHSP)) x$nHSP <- 0
  for(nm in c("year1", "year2")) x[[nm]] <- as.integer(x[[nm]])
  for(nm in c("len1", "len2", "nComp", "nPOP", "nHSP")) x[[nm]] <- as.numeric(x[[nm]])
  if(any(is.na(unlist(x[, c(need, "nPOP", "nHSP")]))))
    stop("'ckmrl' must not contain missing values")
  if(any(x$nComp <= 0)) stop("'ckmrl$nComp' must be positive")
  if(any(x$len1 <= 0 | x$len2 <= 0)) stop("lengths must be positive")
  if(any(x$nPOP < 0 | x$nHSP < 0)) stop("kin counts must not be negative")
  if(any(x$nPOP > x$nComp | x$nHSP > x$nComp))
    stop("more kin pairs than comparisons in at least one row")
  class(x) <- c("ckmrl_data", "data.frame")
  x
}

##' Internal: precompute everything the length-based probabilities need that
##' does not depend on the parameters.
##'
##' The pair table is stored by pair, but the probabilities factorise over
##' \emph{cells}: a row's two members enter only through \eqn{P[a|\ell,y]} for
##' their own (year, length) cell. Ten length classes over three years is
##' thirty cells but four hundred and sixty five pair rows, so working per
##' cell and gathering to rows at the end removes most of the duplication -
##' and lets the sum over the two ages be done as two chained contractions
##' rather than one flat pair-by-age-by-age sweep.
##'
##' With the growth parameters fixed in \code{conf}, the length-at-age density
##' is data and is computed here, once per cell. What remains
##' parameter-dependent is the abundance, the mortality, the maturity when it
##' is a process, the sampled age composition and \eqn{\psi}.
##' @keywords internal
##' @noRd
ckmrlPrep <- function(ckmrl, years, ages, set, propMat = NULL){
  nY <- length(years); nA <- length(ages)
  minY <- min(years); maxY <- max(years)
  minA <- min(ages);  maxA <- max(ages)
  G <- ckmrlGrid(set, ages)
  nPair <- nrow(ckmrl)

  ## ---- cells: the distinct (year, length) combinations the rows refer to ----
  k1 <- paste(ckmrl$year1, ckmrl$len1); k2 <- paste(ckmrl$year2, ckmrl$len2)
  uc <- unique(c(k1, k2))
  cy   <- as.integer(sub(" .*", "", uc))
  clen <- as.numeric(sub(".* ", "", uc))
  ci1 <- match(k1, uc); ci2 <- match(k2, uc); nC <- length(uc)

  troPos <- rep(TRUE, nY)
  if(!is.null(propMat)) troPos <- rowSums(propMat > 0, na.rm = TRUE) > 0
  moPos <- if(is.null(propMat)) matrix(TRUE, nY, nA)
           else (propMat > 0 & !is.na(propMat))
  yI <- function(y) y - minY + 1L

  ## ---- which rows can contribute at all ----
  ## the full (cell, cell, age, age) grid, in plain R and once, purely to find
  ## the rows whose probability is structurally zero
  I  <- rep(seq_len(nC), times = nC*nA*nA)
  J  <- rep(rep(seq_len(nC), each = nC), times = nA*nA)
  AP <- rep(rep(seq_len(nA), each = nC*nC), times = nA)
  AJ <- rep(seq_len(nA), each = nC*nC*nA)
  anyBy <- function(v) matrix(as.vector(matrix(as.numeric(v), nC*nC, nA*nA) %*%
                                        rep(1, nA*nA)) > 0, nC, nC)
  ## parent-offspring, one orientation: I is the parent, J the offspring
  b   <- cy[J] - ages[AJ]
  apb <- ages[AP] - (cy[I] - b)
  okD <- b >= minY & b <= maxY & cy[I] >= b & apb >= minA & apb <= maxA
  okD[is.na(okD)] <- FALSE
  bc  <- pmin(pmax(yI(b), 1L), nY); ac <- pmin(pmax(apb - minA + 1L, 1L), nA)
  okD <- okD & troPos[bc] & moPos[cbind(bc, ac)]
  okD <- anyBy(okD)
  okPOP <- okD | t(okD)
  ## half-sibling: both birth years in range, and some mature shared parent
  b1 <- cy[I] - ages[AP]; b2 <- cy[J] - ages[AJ]
  lo <- pmin(b1, b2); hi <- pmax(b1, b2); d <- hi - lo
  okH <- lo >= minY & hi <= maxY
  okH[is.na(okH)] <- FALSE
  loc <- pmin(pmax(yI(lo), 1L), nY); hic <- pmin(pmax(yI(hi), 1L), nY)
  matePos <- Reduce(`|`, lapply(seq_len(nA), function(a)
    moPos[cbind(loc, rep(a, length(loc)))] &
    moPos[cbind(hic, pmin(a + d, nA))]))
  okHSP <- anyBy(okH & troPos[loc] & troPos[hic] & matePos)

  ## ---- birth years in play, and the gaps between them ----
  bs <- sort(unique(as.vector(outer(cy, ages, "-"))))
  bs <- bs[bs >= minY & bs <= maxY]

  list(nY = nY, nA = nA, nPair = nPair, nC = nC,
       years = years, ages = ages, minY = minY, maxY = maxY,
       minA = minA, maxA = maxA, refLength = set$refLength,
       cy = cy, clen = clen, ci1 = ci1, ci2 = ci2,
       sy = sort(unique(cy)), bset = bs, maxd = max(bs) - min(bs),
       Lmean = G$Lmean, sdLog = G$sdL,
       ## P[length | age] per cell - pure data, since growth is fixed
       PlA = t(vapply(seq_len(nC), function(i)
                stats::dlnorm(clen[i], log(G$Lmean), G$sdL), numeric(nA))),
       okPOP = okPOP, okHSP = okHSP,
       popRow = which(okPOP[cbind(ci1, ci2)]),
       hspRow = which(okHSP[cbind(ci1, ci2)]))
}

##' Internal: length-based close-kin pair probabilities.
##'
##' The single implementation shared by the likelihood and by
##' \code{\link{simulateCKMRL}}, written with \code{offarray} so that every
##' index is the real year, age or length class rather than a hand-built
##' offset into a flat vector.
##'
##' Writing \eqn{\phi(\ell) = (\ell/\ell_{ref})^{\psi}} for fecundity at
##' length and \eqn{A[c,a] = P[a|\ell_c,y_c]} for the age of an animal in cell
##' \eqn{c}, the parent-offspring probability is
##' \deqn{P_{POP}[c_1,c_2] = 2(D[c_1,c_2] + D[c_2,c_1]),\quad
##'   D[c_1,c_2] = \sum_{a_p,a_j} A[c_1,a_p] A[c_2,a_j]
##'   \frac{MO_{a_p',b}\,\phi(\hat\ell)}{TRO_b}}
##' with \eqn{b = y_{c_2} - a_j}, \eqn{a_p'} the parent's age then and
##' \eqn{\hat\ell} its length retrojected to that age. The half-sibling
##' probability sums the birth-year expression over both animals' possible
##' ages, with the unobserved shared parent's length integrated out exactly:
##' it keeps its length-at-age quantile for life, so the required expectations
##' are log-normal moments and have closed forms.
##'
##' The same-cohort term carries a further factor \code{omega}, the
##' lucky-litter effect of Bravington's note: within-season variation in
##' reproductive success that length and age do not explain makes same-cohort
##' pairs share a parent more often than the independent-draw expression
##' gives. It multiplies the \eqn{b_1 = b_2} term alone and is 1 by default.
##'
##' Both sums are done as two chained contractions - over the first age, then
##' the second - rather than over the flat (pair, age, age) grid, because the
##' kernel does not depend on the second animal's length. That keeps the work
##' at \eqn{n_{cell}^2 n_{age}} instead of \eqn{n_{pair} n_{age}^2}.
##'
##' \eqn{P[a|\ell,y]} is built from the length-at-age density and the
##' \emph{population} age composition. Selection - by the fishery, or by
##' whoever picks fish for genotyping - is taken to act on length and year
##' rather than on age, so it cancels out of the conditional distribution once
##' the length is known, and the sample proportions match the population ones.
##' Passing the catch at age instead would apply selectivity twice.
##' @keywords internal
##' @noRd
ckmrlProb <- function(N, Z, MO, Csamp, psi, prep, omega = 1){
  years <- prep$years; ages <- prep$ages; nC <- prep$nC
  minY <- prep$minY; maxY <- prep$maxY; minA <- prep$minA; maxA <- prep$maxA
  lref <- prep$refLength; sy <- prep$sy; bset <- prep$bset
  cc <- seq_len(nC)

  oa <- offarray::offarray
  al <- offarray::autoloop
  Noa <- oa(N,     dimseq = list(y = years, a = ages))
  MOa <- oa(MO,    dimseq = list(y = years, a = ages))
  Za  <- oa(Z,     dimseq = list(y = years, a = ages))
  Cs  <- oa(Csamp, dimseq = list(y = years, a = ages))
  Lm  <- oa(prep$Lmean, dimseq = list(a = ages))
  Pl  <- oa(prep$PlA,   dimseq = list(c = cc, a = ages))
  CY  <- oa(prep$cy,    dimseq = list(c = cc))
  CL  <- oa(prep$clen,  dimseq = list(c = cc))

  ## Fecundity has to be averaged over the length a fish of a given age might
  ## have, and in the half-sibling term over the lengths ONE fish has at two
  ## ages. Both are log-normal moments, so both are exact in closed form. With
  ## l = Lbar_a * d and log d ~ N(0, sd^2),
  ##    E[phi(l_a)]             = phibar(a) * exp(psi^2 sd^2 / 2)
  ##    E[phi(l_a) phi(l_a')]   = phibar(a) phibar(a') * exp(2 psi^2 sd^2)
  ## the second for the same fish at any two ages, a' = a included. Both
  ## corrections are constants, so no quadrature and no tuning parameter.
  ##
  ## This is exact only because: length at age is log-normal; fecundity is a
  ## pure power of length, so the deviation factors out; maturity is indexed by
  ## AGE, so it leaves the expectation untouched; and the spread does not vary
  ## with age. Break any of those - length-based maturity, a saturating
  ## fecundity relation, an age-varying spread - and this is wrong. A
  ## length threshold on maturity still has a closed form, through a truncated
  ## log-normal moment; a general fecundity relation does not, and would want
  ## Gauss-Hermite quadrature rather than the equal-probability grid this
  ## replaced, which converged only as 1/nq.
  kMean <- exp(psi^2 * prep$sdLog^2 / 2)
  kPair <- exp(2 * psi^2 * prep$sdLog^2)
  phiB <- al(indices = list(a = ages), (Lm[a]/lref)^psi)
  mPhi <- al(indices = list(a = ages), phiB[a] * kMean)
  TRO  <- al(indices = list(y = years), SUMOVER = list(a = ages),
             Noa[y,a]*MOa[y,a]*mPhi[a])

  ## age given length and year, per cell, from the sampled composition
  csum <- al(indices = list(y = years), SUMOVER = list(a = ages), Cs[y,a])
  raw  <- al(indices = list(c = cc, a = ages), Pl[c,a]*Cs[CY[c],a]/csum[CY[c]])
  rsum <- al(indices = list(c = cc), SUMOVER = list(a = ages), raw[c,a])
  A    <- al(indices = list(c = cc, a = ages), raw[c,a]/rsum[c])

  ## ---- parent-offspring: sum over the parent's age, then the offspring's ----
  S <- al(indices = list(c = cc, yj = sy, aj = ages), SUMOVER = list(ap = ages), {
    b   <- yj - aj
    apb <- ap - (CY[c] - b)
    ok  <- (b >= minY)*(b <= maxY)*(CY[c] >= b)*(apb >= minA)*(apb <= maxA)
    bc  <- pmin(pmax(b, minY), maxY); ac <- pmin(pmax(apb, minA), maxA)
    ok * A[c,ap] * MOa[bc,ac] * ((Lm[ac]/Lm[ap]*CL[c])/lref)^psi / TRO[bc]
  })
  D   <- al(indices = list(i = cc, j = cc), SUMOVER = list(aj = ages),
            S[i,CY[j],aj]*A[j,aj])
  POP <- al(indices = list(i = cc, j = cc), 2*(D[i,j] + D[j,i]))

  ## ---- half-sibling ----
  ## survival of the shared parent across the gap, one pass per gap
  cum  <- al(indices = list(b = bset, a = ages), 0*Za[minY,minA])
  SURV <- oa(rep(0*N[1], length(bset)*length(ages)*(prep$maxd+1L)),
             dimseq = list(b = bset, a = ages, d = 0:prep$maxd))
  SURV[,,0] <- exp(-cum)
  for(k in seq_len(prep$maxd)){
    cum <- al(indices = list(b = bset, a = ages),
              cum[b,a] + Za[pmin(b+k-1L, maxY), pmin(a+k-1L, maxA)])
    SURV[,,k] <- exp(-cum)
  }
  PBB0 <- al(indices = list(b1 = bset, b2 = bset),
             SUMOVER = list(a = ages), {
    lo <- pmin(b1,b2); hi <- pmax(b1,b2); d <- hi - lo
    a2 <- pmin(a + d, maxA)
    4*Noa[lo,a]*MOa[lo,a]*phiB[a]*SURV[lo,a,d]*MOa[hi,a2]*phiB[a2]*kPair /
      (TRO[lo]*TRO[hi])
  })
  ## the lucky-litter effect: within-season variation in reproductive success
  ## that length and age do not explain makes two animals of the SAME cohort
  ## share a parent more often than the independent-draw expression gives.
  ## omega scales that term and only that term, and is 1 by default.
  PBB <- al(indices = list(b1 = bset, b2 = bset),
            PBB0[b1,b2] * ((b1 == b2)*omega + (b1 != b2)))
  T1 <- al(indices = list(c = cc, yj = sy, a2 = ages), SUMOVER = list(a1 = ages), {
    b1 <- CY[c] - a1; b2 <- yj - a2
    ok <- (b1 >= minY)*(b1 <= maxY)*(b2 >= minY)*(b2 <= maxY)
    ok * A[c,a1] * PBB[pmin(pmax(b1,minY),maxY), pmin(pmax(b2,minY),maxY)]
  })
  HSP <- al(indices = list(i = cc, j = cc), SUMOVER = list(a2 = ages),
            T1[i,CY[j],a2]*A[j,a2])

  ## gather the cell grid back onto the pair rows
  g <- cbind(prep$ci1, prep$ci2)
  list(pPOP = as.array(POP)[g] * as.numeric(prep$okPOP[g]),
       pHSP = as.array(HSP)[g] * as.numeric(prep$okHSP[g]),
       TRO  = as.vector(TRO))
}


##' Build a length-based close-kin mark-recapture pair table
##'
##' The length analogue of \code{\link{ckmrData}}: genotyped individuals are
##' grouped into (year, length-class) cells and one row is produced per pair of
##' distinct cells, with \code{nComp = n1 * n2} comparisons.
##'
##' Unlike the age-based table, pairs drawn from the \emph{same} cell are kept,
##' with \code{n(n-1)/2} comparisons. Two fish of the same length in the same
##' year need not be the same age, so they can be a parent-offspring or a
##' cross-cohort half-sibling pair.
##'
##' @param samples data frame with columns \code{year}, \code{len} (the length
##' class, given by its representative length) and \code{n}.
##' @param pop optional data frame of observed parent-offspring pairs with
##' columns \code{year1}, \code{len1}, \code{year2}, \code{len2}, \code{n};
##' the two cells may be given in either order.
##' @param hsp optional data frame of observed half-sibling pairs, same layout.
##' @return a data frame of class \code{ckmrl_data}.
##' @seealso \code{\link{simulateCKMRL}}, \code{\link{ckmrData}}
##' @export
ckmrlData <- function(samples, pop = NULL, hsp = NULL){
  if(!all(c("year", "len", "n") %in% names(samples)))
    stop("'samples' must have columns year, len and n")
  s <- stats::aggregate(list(n = as.numeric(samples$n)),
                        by = list(year = as.integer(samples$year),
                                  len  = as.numeric(samples$len)), FUN = sum)
  s <- s[s$n > 0, , drop = FALSE]
  if(nrow(s) < 1) stop("no non-empty (year, length) cells")
  s <- s[order(s$year, s$len), , drop = FALSE]
  ij <- if(nrow(s) >= 2) utils::combn(nrow(s), 2) else matrix(integer(0), 2, 0)
  ij <- cbind(ij, rbind(seq_len(nrow(s)), seq_len(nrow(s))))
  same <- ij[1, ] == ij[2, ]
  ret <- data.frame(year1 = s$year[ij[1, ]], len1 = s$len[ij[1, ]],
                    year2 = s$year[ij[2, ]], len2 = s$len[ij[2, ]],
                    nComp = ifelse(same, s$n[ij[1, ]] * (s$n[ij[1, ]] - 1) / 2,
                                   s$n[ij[1, ]] * s$n[ij[2, ]]),
                    nPOP = 0, nHSP = 0)
  ret <- ret[ret$nComp > 0, , drop = FALSE]
  rownames(ret) <- NULL
  k <- ckmrlKey(ret$year1, ret$len1, ret$year2, ret$len2)
  addKin <- function(tab, what){
    if(is.null(tab) || nrow(tab) == 0) return(NULL)
    if(!all(c("year1", "len1", "year2", "len2", "n") %in% names(tab)))
      stop("'", what, "' must have columns year1, len1, year2, len2, n")
    m <- match(ckmrlKey(tab$year1, tab$len1, tab$year2, tab$len2), k)
    if(any(is.na(m)))
      stop(sum(is.na(m)), " ", what, " pair(s) refer to cells not in 'samples'")
    tapply(as.numeric(tab$n), m, sum)
  }
  aP <- addKin(pop, "pop"); if(!is.null(aP)) ret$nPOP[as.integer(names(aP))] <- aP
  aH <- addKin(hsp, "hsp"); if(!is.null(aH)) ret$nHSP[as.integer(names(aH))] <- aH
  class(ret) <- c("ckmrl_data", "data.frame")
  ret
}

##' Internal: order-insensitive key for a pair of (year, length) cells.
##' @keywords internal
##' @noRd
ckmrlKey <- function(y1, l1, y2, l2){
  k1 <- paste(y1, signif(l1, 10), sep = ".")
  k2 <- paste(y2, signif(l2, 10), sep = ".")
  ifelse(k1 <= k2, paste(k1, k2, sep = "|"), paste(k2, k1, sep = "|"))
}

##' Simulate length-based close-kin data from a fitted model
##'
##' Draws a genotyped sample, assigns each fish a length from the
##' length-at-age distribution implied by the configuration, bins the lengths,
##' and simulates the number of parent-offspring and half-sibling pairs found.
##' The probabilities come from the same engine the likelihood uses.
##'
##' Ages are drawn from the model's catch at age, so the sample is
##' age-selective in the way a fishery-dependent sample would be; the
##' probabilities weight by that same composition, which is why they stay
##' consistent.
##'
##' @param fit a fitted \code{sam} object.
##' @param years years in which fish are genotyped.
##' @param n number genotyped per year, recycled over \code{years}.
##' @param lengthBins cut points defining the length classes. \code{NULL}
##' (default) builds a grid spanning the length-at-age distribution.
##' @param nBin number of length classes when \code{lengthBins} is \code{NULL}.
##' @param seed optional random seed.
##' @return a \code{ckmrl_data} pair table.
##' @seealso \code{\link{ckmrlData}}, \code{\link{ckmrltable}}
##' @export
simulateCKMRL <- function(fit, years, n, lengthBins = NULL, nBin = 12,
                          seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  years <- as.integer(years); n <- rep(as.numeric(n), length.out = length(years))
  modYears <- fit$data$years
  ages <- fit$conf$minAge:fit$conf$maxAge
  if(!all(years %in% modYears))
    stop("sampling years outside the model year range: ",
         paste(setdiff(years, modYears), collapse = ", "))
  set <- ckmrlSettings(fit$conf)
  if(set$estPsi == 1L && length(fit$pl$logPsiL) > 0)
    set$psi <- exp(fit$pl$logPsiL[1])
  if(set$estOmega == 1L && length(fit$pl$logOmegaL) > 0)
    set$omega <- exp(fit$pl$logOmegaL[1])
  G <- ckmrlGrid(set, ages)
  if(is.null(lengthBins))
    lengthBins <- seq(exp(stats::qnorm(0.001, log(min(G$Lmean)), G$sdL)),
                      exp(stats::qnorm(0.999, log(max(G$Lmean)), G$sdL)),
                      length = nBin + 1L)
  mids <- lengthBins[-1] - 0.5 * diff(lengthBins)
  ## Ages are drawn from the population, not the catch: the probabilities
  ## assume selection acts on length, so a sample drawn with age-based
  ## selection would not be consistent with them.
  cay <- ntable(fit)[match(years, modYears), , drop = FALSE]
  cay[!is.finite(cay) | cay < 0] <- 0
  samples <- do.call(rbind, lapply(seq_along(years), function(i){
    a <- stats::rmultinom(1, n[i], cay[i, ])[, 1]
    l <- unlist(lapply(seq_along(a), function(j)
      if(a[j] > 0) stats::rlnorm(a[j], log(G$Lmean[j]), G$sdL) else numeric(0)))
    cl <- cut(l, lengthBins, labels = FALSE)
    tb <- table(factor(cl, levels = seq_along(mids)))
    data.frame(year = years[i], len = mids, n = as.numeric(tb))
  }))
  ck <- ckmrlData(samples[samples$n > 0, , drop = FALSE])
  prep <- ckmrlPrep(ck, modYears, ages, set,
                    propMat = if(is.null(fit$conf$matureModel) ||
                                 fit$conf$matureModel == 0) fit$data$propMat else NULL)
  pr <- ckmrlProb(ntable(fit) * set$scale,
                  faytable(fit) + fit$data$natMor,
                  fit$data$propMat, ntable(fit), set$psi, prep, set$omega)
  ck$nPOP <- stats::rpois(nrow(ck), pmax(0, pr$pPOP) * ck$nComp)
  ck$nHSP <- stats::rpois(nrow(ck), pmax(0, pr$pHSP) * ck$nComp)
  attr(ck, "lengthBins") <- lengthBins
  ck
}

##' Observed and expected length-based close-kin pair counts
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return the pair table with \code{ePOP}, \code{rPOP}, \code{eHSP},
##' \code{rHSP} added.
##' @export
ckmrltable <- function(fit, ...) UseMethod("ckmrltable")

##' @rdname ckmrltable
##' @method ckmrltable sam
##' @export
ckmrltable.sam <- function(fit, ...){
  ck <- fit$data$ckmrl
  if(is.null(ck)) stop("this fit carries no length-based CKMR data")
  if(is.null(fit$rep$pPOPL))
    stop("no reported probabilities - was the fit run with usePOPL/useHSPL?")
  ret <- as.data.frame(ck[, c("year1", "len1", "year2", "len2", "nComp",
                              "nPOP", "nHSP")])
  ret$ePOP <- fit$rep$pPOPL * ck$nComp
  ret$rPOP <- (ret$nPOP - ret$ePOP) / sqrt(pmax(ret$ePOP, .Machine$double.eps))
  ret$eHSP <- fit$rep$pHSPL * ck$nComp
  ret$rHSP <- (ret$nHSP - ret$eHSP) / sqrt(pmax(ret$eHSP, .Machine$double.eps))
  ret
}
