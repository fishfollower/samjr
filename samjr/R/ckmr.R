##' Internal: canonical order-insensitive key for a pair of (year, age) cells.
##' @keywords internal
##' @noRd
ckmrKey <- function(y1, a1, y2, a2){
  k1 <- paste(y1, a1, sep = ".")
  k2 <- paste(y2, a2, sep = ".")
  ifelse(k1 <= k2, paste(k1, k2, sep = "|"), paste(k2, k1, sep = "|"))
}

##' Internal: precompute every parameter-independent CKMR index.
##'
##' All of the branching in the close-kin probabilities (which animal can be
##' the parent, birth years, the candidate parent-age range, the cumulative-Z
##' path between two birth years) depends only on the pair table and the model
##' dimensions - never on the parameters. It is therefore computed once, here,
##' and cached in \code{dat}, so that \code{\link{ckmrProb}} reduces to
##' whole-vector gathers on integer indices.
##'
##' Column-major linear indices are used throughout: element \code{[y, a]} of a
##' \code{nY x nA} matrix is at \code{y + (a - 1) * nY}.
##'
##' @param ckmr pair table with columns \code{year1}, \code{age1},
##' \code{year2}, \code{age2}, \code{nComp} (and optionally \code{nPOP},
##' \code{nHSP}).
##' @param years model year vector.
##' @param ages model age vector (\code{minAge:maxAge}).
##' @param plusExtra how many years beyond \code{maxAge} a plus-group animal's
##' true age is allowed to reach when its birth year is marginalised over. The
##' weights fall off geometrically with the plus-group mortality, so the mass
##' beyond a few years is negligible; capping keeps the expansion small on long
##' time series, where the plus group would otherwise nominally admit an animal
##' as old as the series itself.
##' @param plusNodes how many of the leading candidate birth years are resolved
##' one by one. The remainder, out to \code{plusExtra}, are lumped into a single
##' node whose weight is still summed exactly, so no mass is lost - only the
##' birth year within that lump is approximated. The weights fall off
##' geometrically, so a handful of leading nodes carry almost all of it.
##' @param propMat optional observed proportion-mature matrix. When supplied
##' (i.e. maturity is data, not a GMRF process) it is used to drop pairs whose
##' probability is identically zero, which keeps \code{dpois} away from
##' \code{lambda = 0}.
##' @return a list of integer index vectors consumed by \code{\link{ckmrProb}}.
##' @keywords internal
##' @noRd
ckmrPrep <- function(ckmr, years, ages, propMat = NULL, plusExtra = 10L,
                     plusNodes = 6L){
  nY <- length(years); nA <- length(ages)
  minYear <- min(years); minAge <- min(ages)
  y1 <- as.integer(ckmr$year1); a1 <- as.integer(ckmr$age1)
  y2 <- as.integer(ckmr$year2); a2 <- as.integer(ckmr$age2)
  nPair <- length(y1)
  c1 <- y1 - a1; c2 <- y2 - a2

  ## years in which the population has no reproductive output at all would
  ## give TRO == 0; such years cannot carry CKMR information.
  troPos <- rep(TRUE, nY)
  if(!is.null(propMat)) troPos <- rowSums(propMat > 0, na.rm = TRUE) > 0

  ## ---- orientation: slot p = earlier cohort (candidate parent) ----
  isFirst <- c1 < c2
  py <- ifelse(isFirst, y1, y2); pa <- ifelse(isFirst, a1, a2)
  jy <- ifelse(isFirst, y2, y1); ja <- ifelse(isFirst, a2, a1)

  ## ---- parent-offspring ----
  b    <- jy - ja                      # offspring birth year
  lag  <- py - b                       # years from birth to parent sampling
  bI   <- b  - minYear + 1L
  pyI  <- py - minYear + 1L
  paI  <- pa - minAge  + 1L
  ## candidate parent ages at the birth year: a' with min(nA, a' + lag) == paI
  aLo <- ifelse(paI < nA, paI - lag, pmax(1L, nA - lag))
  aHi <- ifelse(paI < nA, paI - lag, nA)
  popOK <- (lag >= 0L) &                                   # lethal sampling
           (bI >= 1L) & (bI <= nY) & (pyI >= 1L) & (pyI <= nY) &
           (paI >= 1L) & (paI <= nA) & (aLo >= 1L) & (aHi >= aLo)
  popOK[is.na(popOK)] <- FALSE
  popOK <- popOK & troPos[pmax(1L, pmin(nY, bI))]
  if(!is.null(propMat)){
    ## drop rows whose candidate parents are all immature (probability is
    ## structurally zero, so they carry no information and would give dpois(0,0))
    cs <- cbind(0, t(apply(matrix(as.numeric(propMat > 0 & !is.na(propMat)),
                                  nY, nA), 1, cumsum)))
    bC  <- pmin(nY, pmax(1L, bI))
    loC <- pmin(nA, pmax(1L, aLo))
    hiC <- pmin(nA, pmax(1L, aHi))
    lo <- ifelse(popOK, cs[cbind(bC, loC)], 1)
    hi <- ifelse(popOK, cs[cbind(bC, hiC + 1L)], 2)
    popOK <- popOK & ((hi - lo) > 0)
  }

  ## ---- half-sibling ----
  ## Both birth years enter this probability directly, and an animal recorded
  ## at the plus-group age has an unknown true age - so its birth year is only
  ## bounded above. Each pair is therefore expanded over the birth years its
  ## members could have had; ckmrProb weights the terms by the plus-group age
  ## composition, which it derives from N and Z. An animal of known age
  ## contributes a single term of weight one, so pairs with no plus-group
  ## member take an unchanged single-term path.
  y1I <- y1 - minYear + 1L; a1I <- a1 - minAge + 1L
  y2I <- y2 - minYear + 1L; a2I <- a2 - minAge + 1L
  plusNodes <- max(1L, as.integer(plusNodes))
  ## Offsets covered by each node for a plus-group animal: the first
  ## plusNodes offsets one by one, then a single lump out to plusExtra whose
  ## representative sits one year past its start - the placement that tracks
  ## the weighted average of a geometrically decaying weight most closely.
  nodeBins <- function(bi){
    maxO <- min(as.integer(plusExtra), bi - 1L)
    if(maxO < 0L) return(list())
    ex <- seq.int(0L, min(plusNodes - 1L, maxO))
    out <- lapply(ex, function(o) list(rep = o, bin = o))
    if(maxO >= plusNodes)
      out <- c(out, list(list(rep = min(plusNodes, maxO) + 1L,
                              bin = seq.int(plusNodes, maxO))))
    out <- lapply(out, function(z){ z$rep <- min(z$rep, maxO); z })
    out[vapply(out, function(z) all(bi - z$bin + nA + minAge - 2L >= 1L), logical(1))]
  }
  ## registry of (sampling year, node) pairs and the birth years each averages
  regKey <- character(0); regY <- integer(0); regBin <- list()
  registerNode <- function(yI, bi, nb){
    k <- paste(yI, bi - nb$rep, sep = ".")
    m <- match(k, regKey)
    if(is.na(m)){
      regKey <<- c(regKey, k); regY <<- c(regY, yI)
      regBin[[length(regKey)]] <<- bi - nb$bin
      m <- length(regKey)
    }
    m
  }
  ## returns representative birth years and, for plus-group animals, node ids
  candBirth <- function(yI, aI){
    bi <- yI - aI - minAge + 1L
    if(is.na(bi) || bi < 1L) return(list(b = integer(0), node = integer(0)))
    if(aI < nA || nA < 2L) return(list(b = bi, node = 0L))
    nb <- nodeBins(bi)
    if(!length(nb)) return(list(b = integer(0), node = integer(0)))
    list(b = bi - vapply(nb, function(z) z$rep, integer(1)),
         node = vapply(nb, function(z) registerNode(yI, bi, z), integer(1)))
  }
  inR <- function(v) !is.na(v) & v >= 1L & v <= nY
  baseOK <- inR(y1I) & inR(y2I) & !is.na(a1I) & !is.na(a2I) &
            a1I >= 1L & a1I <= nA & a2I >= 1L & a2I <= nA
  biRec1 <- y1I - a1I - minAge + 1L
  biRec2 <- y2I - a2I - minAge + 1L
  acc <- vector("list", nPair)
  for(i in seq_len(nPair)){
    if(!baseOK[i]) next
    ## a pair whose recorded cohorts agree is excluded, as before
    if(!is.na(biRec1[i]) && !is.na(biRec2[i]) && biRec1[i] == biRec2[i]) next
    c1 <- candBirth(y1I[i], a1I[i]); c2 <- candBirth(y2I[i], a2I[i])
    k1 <- c1$b; k2 <- c2$b
    if(!length(k1) || !length(k2)) next
    g1 <- rep(k1, times = length(k2)); g2 <- rep(k2, each = length(k1))
    n1 <- rep(c1$node, times = length(k2)); n2 <- rep(c2$node, each = length(k1))
    ## When a plus-group member makes the cohorts ambiguous, the same-cohort
    ## configuration has to stay in the average: such a pair cannot be
    ## recognised as same-cohort and so cannot be excluded from the data
    ## either, and dropping it here would leave those kin observed but never
    ## predicted.
    allowSame <- isTRUE(a1I[i] == nA) || isTRUE(a2I[i] == nA)
    keep <- (allowSame | (g1 != g2)) & troPos[g1] & troPos[g2]
    if(!any(keep)) next
    acc[[i]] <- cbind(i, g1[keep], g2[keep], n1[keep], n2[keep])
  }
  EX <- do.call(rbind, acc)
  if(is.null(EX)) EX <- matrix(integer(0), nrow = 0, ncol = 5)
  hspOK <- logical(nPair)
  if(nrow(EX) > 0) hspOK[unique(EX[, 1])] <- TRUE
  isPlusPair <- (a1I == nA) & !is.na(a1I) | (a2I == nA) & !is.na(a2I)
  isPlusPair[is.na(isPlusPair)] <- FALSE

  popRow <- which(popOK)
  hspRow <- which(hspOK)
  hspSimpleRow <- which(hspOK & !isPlusPair)     # both ages known exactly
  hspMargRow   <- which(hspOK &  isPlusPair)     # at least one plus-group animal
  simple <- popRow[aLo[popRow] == aHi[popRow]]
  plus   <- popRow[aLo[popRow] <  aHi[popRow]]

  ## ---- shared survival grid ----
  ## start years actually needed, and the largest lag over which a cohort has
  ## to be carried forward
  exB1 <- if(nrow(EX)) pmin(EX[, 2], EX[, 3]) else integer(0)
  exB2 <- if(nrow(EX)) pmax(EX[, 2], EX[, 3]) else integer(0)
  ## Many rows and expansions need the same pair of birth years, so the
  ## half-sibling probability is evaluated once per distinct (b1, b2) and then
  ## gathered. Both the single-term and the marginalised paths become lookups.
  sb1 <- pmin(biRec1[hspSimpleRow], biRec2[hspSimpleRow])
  sb2 <- pmax(biRec1[hspSimpleRow], biRec2[hspSimpleRow])
  bbKey <- c(paste(sb1, sb2, sep = "."), paste(exB1, exB2, sep = "."))
  ubb   <- unique(bbKey)
  bbB1  <- as.integer(sub("[.].*", "", ubb))
  bbB2  <- as.integer(sub(".*[.]", "", ubb))
  bbD   <- bbB2 - bbB1
  simpleIdx <- if(length(sb1)) match(paste(sb1, sb2, sep = "."), ubb) else integer(0)
  exIdx     <- if(length(exB1)) match(paste(exB1, exB2, sep = "."), ubb) else integer(0)
  uSet   <- sort(unique(c(bI[plus], bbB1)))
  nU     <- length(uSet)
  maxLag <- max(c(lag[plus], bbD, 0L))
  zLin <- vector("list", maxLag)
  if(maxLag > 0L && nU > 0L){
    Uu <- rep(seq_len(nU), times = nA)     # u varies fastest
    Ua <- rep(seq_len(nA), each  = nU)
    for(k in seq_len(maxLag)){
      zLin[[k]] <- pmin(nY, uSet[Uu] + k - 1L) +
                   (pmin(nA, Ua + k - 1L) - 1L) * nY
    }
  }
  nUA <- nU * nA
  survLin <- function(rowLag, rowU, nRow){
    ## flat index into SURV for every (row, age), row varying fastest
    uPos <- match(rowU, uSet)
    rep(rowLag, times = nA) * nUA +
      (rep(seq_len(nA), each = nRow) - 1L) * nU +
      rep(uPos, times = nA)
  }

  out <- list(nY = nY, nA = nA, nPair = nPair, nU = nU, maxLag = maxLag,
              zLin = zLin, popRow = popRow, hspRow = hspRow)

  ## POP, non-plus-group parent: single candidate age, survival cancels
  out$popSimpleRow    <- simple
  out$popSimpleFecLin <- bI[simple] + (aLo[simple] - 1L) * nY
  out$popSimpleTRO    <- bI[simple]

  ## POP, plus-group parent: average fecundity over the candidate age range
  nP <- length(plus)
  out$popPlusRow <- plus
  if(nP > 0L){
    Ap <- rep(seq_len(nA), each = nP)
    out$popPlusLin   <- rep(bI[plus], times = nA) + (Ap - 1L) * nY
    out$popPlusMask  <- as.numeric(Ap >= rep(aLo[plus], times = nA) &
                                   Ap <= rep(aHi[plus], times = nA))
    out$popPlusSurv  <- survLin(lag[plus], bI[plus], nP)
    out$popPlusNpy   <- pyI[plus] + (paI[plus] - 1L) * nY
    out$popPlusTRO   <- bI[plus]
  }else{
    out$popPlusLin <- out$popPlusMask <- out$popPlusSurv <- integer(0)
    out$popPlusNpy <- out$popPlusTRO <- integer(0)
  }

  ## HSP, pairs whose members both have a known age: one term, weight one
  nBB <- length(ubb)
  out$hspSimpleRow <- hspSimpleRow
  out$hspSimpleIdx <- simpleIdx
  if(nBB > 0L){
    Ab <- rep(seq_len(nA), each = nBB)
    out$hspBBLin1 <- rep(bbB1, times = nA) + (Ab - 1L) * nY
    out$hspBBLin2 <- rep(bbB2, times = nA) +
                     (pmin(nA, Ab + rep(bbD, times = nA)) - 1L) * nY
    out$hspBBSurv <- survLin(bbD, bbB1, nBB)
    out$hspBBTRO1 <- bbB1
    out$hspBBTRO2 <- bbB2
  }else{
    out$hspBBLin1 <- out$hspBBLin2 <- out$hspBBSurv <- integer(0)
    out$hspBBTRO1 <- out$hspBBTRO2 <- integer(0)
  }

  ## HSP, pairs with a plus-group member: expanded over candidate birth years
  ## and padded to a rectangle, so the marginalisation is two matrix products.
  ## Rows are grouped by how many of the two members are in the plus group, so
  ## that pairs needing a handful of terms are not padded out to the square
  ## count the two-plus-group pairs need.
  nPlusMem <- as.integer(!is.na(a1I) & a1I == nA) + as.integer(!is.na(a2I) & a2I == nA)
  mkMarg <- function(rows){
    nM <- length(rows)
    if(nM == 0L) return(list(nM = 0L, row = integer(0)))
    mi <- match(EX[, 1], rows); kp <- !is.na(mi)
    exM <- EX[kp, , drop = FALSE]; mi <- mi[kp]
    o <- order(mi); exM <- exM[o, , drop = FALSE]; mi <- mi[o]
    eNo <- stats::ave(seq_along(mi), mi, FUN = seq_along)
    maxE <- max(eNo); nME <- nM * maxE
    pos <- (eNo - 1L) * nM + mi
    bbi <- rep(1L, nME); mask <- numeric(nME)
    bbi[pos] <- exIdx[kp][o]
    mask[pos] <- 1
    slot <- function(nodeCol){
      nd <- exM[, nodeCol]; isP <- nd > 0L
      list(at = pos[isP], node = nd[isP])
    }
    list(nM = nM, row = rows, maxE = maxE, mask = mask, bbi = bbi,
         w1 = slot(4L), w2 = slot(5L))
  }
  out$hspMarg <- list(mkMarg(which(hspOK & nPlusMem == 1L)),
                      mkMarg(which(hspOK & nPlusMem == 2L)))
  out$hspMargRow <- hspMargRow
  ## Plus-group node weights. Each node averages a set of candidate birth
  ## years; the weight of every (sampling year, birth year) pair is computed
  ## once and summed into its node by a fixed incidence matrix, so no mass is
  ## lost however coarsely the nodes are laid out.
  nNode <- length(regKey)
  out$nPlusNode <- nNode
  if(nNode > 0L){
    allY <- rep(regY, lengths(regBin)); allB <- unlist(regBin)
    key <- paste(allY, allB, sep = ".")
    uk  <- unique(key)
    wYi <- as.integer(sub("[.].*", "", uk)); wBi <- as.integer(sub(".*[.]", "", uk))
    yp  <- wBi + nA + minAge - 1L
    out$wLin <- (yp - 1L) + (nA - 2L) * nY
    out$wCy  <- wYi; out$wCyp <- yp
    out$wDen <- wYi + (nA - 1L) * nY
    M <- matrix(0, nNode, length(uk))
    M[cbind(rep(seq_len(nNode), lengths(regBin)), match(key, uk))] <- 1
    out$binMat <- M
    out$nodeYear <- regY
    ## Sums the nodes belonging to one sampling year, so the weights can be
    ## normalised into a distribution over birth years. Reconstructing the
    ## plus group by the survival recursion and dividing by N[y, A] would
    ## only be a distribution if the plus group accumulated deterministically;
    ## SAM's N carries process noise, and on the mackerel fit the
    ## reconstruction comes to 1.12 times N[y, A] rather than 1. Normalising
    ## also puts the truncated tail back in proportionally, instead of
    ## discarding it.
    out$normMat <- outer(regY, regY, "==") * 1
  }
  ## fixed lower-triangular accumulator for the plus-group column of Z
  out$hspCumL <- outer(seq_len(nY), seq_len(nY), function(i, j) as.numeric(j < i))
  out$hspZcol <- seq_len(nY) + (nA - 1L) * nY
  out
}

##' Internal: close-kin pair probabilities.
##'
##' The single implementation of the parent-offspring and half-sibling
##' probabilities, shared by the likelihood and by
##' \code{\link{simulateCKMR}}. Only \code{*}, \code{+}, \code{exp}, linear
##' indexing and \code{\%*\% rep(1, n)} are used, so the same code runs on
##' plain numeric matrices and on RTMB advectors.
##'
##' Let \eqn{b} be a birth year, \eqn{A} the plus group, \eqn{fec} the
##' per-capita reproductive output and \eqn{TRO_y = \sum_a N_{a,y} fec_{a,y}}.
##' For an adult observed at age \eqn{a_p} in year \eqn{y_p},
##' \deqn{P_{POP} = \frac{2}{TRO_b}
##'   \frac{\sum_{a'} N_{a',b}\, S(a', b \to y_p)\, fec_{a',b}}{N_{a_p, y_p}}}
##' where \eqn{a'} runs over the ages at \eqn{b} that land in the observed
##' class at \eqn{y_p} - a single age when \eqn{a_p < A}, in which case the
##' expression collapses to \eqn{2 fec_{a_p - lag, b} / TRO_b}, and a range
##' when the adult is in the plus group and its age at \eqn{b} is only bounded
##' below. For two juveniles born in \eqn{b_1 < b_2},
##' \deqn{P_{HSP} = \frac{4}{TRO_{b_1} TRO_{b_2}} \sum_a N_{a,b_1}\,
##'   fec_{a,b_1}\, fec_{\min(A, a + d), b_2}\, S(a, b_1 \to b_2).}
##' An animal recorded at the plus-group age has an unknown true age, so its
##' birth year is only bounded above. For such pairs the expression above is
##' averaged over the birth years the members could have had, weighted by the
##' plus-group age composition implied by \eqn{N} and \eqn{Z}.
##'
##' @param N year-by-age abundance, already on the individual scale.
##' @param Z year-by-age total mortality.
##' @param fec year-by-age per-capita reproductive output.
##' @param prep the index list from \code{\link{ckmrPrep}}.
##' @return list with \code{pPOP}, \code{pHSP} (both of length
##' \code{prep$nPair}, zero on rows the model cannot use) and \code{TRO}.
##' @keywords internal
##' @noRd
ckmrProb <- function(N, Z, fec, prep){
  nY <- prep$nY; nA <- prep$nA
  ones <- rep(1, nA)
  TRO  <- as.vector(matrix(N * fec, nY, nA) %*% ones)
  zero <- N[1] * 0
  pPOP <- rep(zero, prep$nPair)
  pHSP <- rep(zero, prep$nPair)

  ## one recursion over lag - not over pairs - builds every survival term.
  ## SURV element (d, u, a) sits at d * nU * nA + (a - 1) * nU + u.
  SURV <- NULL
  if(prep$maxLag > 0L && prep$nU > 0L){
    nUA  <- prep$nU * nA
    SURV <- rep(zero, (prep$maxLag + 1L) * nUA)
    cum  <- rep(zero, nUA)
    SURV[seq_len(nUA)] <- exp(-cum)
    for(k in seq_len(prep$maxLag)){
      cum <- cum + Z[prep$zLin[[k]]]
      SURV[k * nUA + seq_len(nUA)] <- exp(-cum)
    }
  }

  if(length(prep$popSimpleRow) > 0L)
    pPOP[prep$popSimpleRow] <-
      2 * fec[prep$popSimpleFecLin] / TRO[prep$popSimpleTRO]

  nP <- length(prep$popPlusRow)
  if(nP > 0L){
    X <- prep$popPlusMask * N[prep$popPlusLin] * fec[prep$popPlusLin] *
         SURV[prep$popPlusSurv]
    num <- as.vector(matrix(X, nP, nA) %*% ones)
    pPOP[prep$popPlusRow] <-
      2 * num / N[prep$popPlusNpy] / TRO[prep$popPlusTRO]
  }

  nBB <- length(prep$hspBBTRO1)
  if(nBB > 0L){
    X <- N[prep$hspBBLin1] * fec[prep$hspBBLin1] * fec[prep$hspBBLin2] *
         SURV[prep$hspBBSurv]
    num <- as.vector(matrix(X, nBB, nA) %*% ones)
    pBB <- 4 * num / (TRO[prep$hspBBTRO1] * TRO[prep$hspBBTRO2])
    if(length(prep$hspSimpleRow) > 0L)
      pHSP[prep$hspSimpleRow] <- pBB[prep$hspSimpleIdx]
  }

  ## Pairs with a plus-group member: the same expression per candidate birth
  ## year, averaged over the plus-group age composition. That composition
  ## follows from the model's own plus-group accumulation - each year the
  ## class gains those turning maxAge and keeps the survivors of last year's
  ## class - so the number in the plus group at year y born in year b is
  ##   N[y'-1, A-1] exp(-Z[y'-1, A-1]) exp(-sum_{k=y'}^{y-1} Z[k, A]),
  ## with y' the year that cohort reached maxAge. Dividing by N[y, A] turns it
  ## into P(born in b | observed in the plus group at y). The weights need not
  ## sum to one: the shortfall is animals born before the model starts, whose
  ## probability cannot be evaluated, and which contribute zero here just as
  ## they do elsewhere.
  if(any(vapply(prep$hspMarg, function(g) g$nM > 0L, logical(1)))){
    cumZ <- as.vector(prep$hspCumL %*% Z[prep$hspZcol])
    wv <- N[prep$wLin] * exp(-Z[prep$wLin]) *
          exp(-(cumZ[prep$wCy] - cumZ[prep$wCyp]))
    ckmrPlusWeight <- as.vector(prep$binMat %*% wv)
    ckmrPlusWeight <- ckmrPlusWeight /
                      as.vector(prep$normMat %*% ckmrPlusWeight)
    REPORT(ckmrPlusWeight)
    for(g in prep$hspMarg){
      if(g$nM == 0L) next
      nME <- g$nM * g$maxE
      pe <- pBB[g$bbi]
      w <- rep(zero + 1, nME)
      for(sl in list(g$w1, g$w2)){
        if(length(sl$at) == 0L) next
        w[sl$at] <- w[sl$at] * ckmrPlusWeight[sl$node]
      }
      w <- w * g$mask
      pHSP[g$row] <- as.vector(matrix(pe * w, g$nM, g$maxE) %*% rep(1, g$maxE))
    }
  }

  list(pPOP = pPOP, pHSP = pHSP, TRO = TRO)
}

##' Internal: validate a CKMR pair table.
##' @keywords internal
##' @noRd
ckmrCheck <- function(x){
  need <- c("year1", "age1", "year2", "age2", "nComp")
  if(!is.data.frame(x)) stop("'ckmr' must be a data frame (see ckmrData)")
  if(!all(need %in% names(x)))
    stop("'ckmr' must have columns ", paste(need, collapse = ", "))
  if(is.null(x$nPOP)) x$nPOP <- 0
  if(is.null(x$nHSP)) x$nHSP <- 0
  for(nm in c("year1", "age1", "year2", "age2")) x[[nm]] <- as.integer(x[[nm]])
  for(nm in c("nComp", "nPOP", "nHSP")) x[[nm]] <- as.numeric(x[[nm]])
  if(any(is.na(unlist(x[, c(need, "nPOP", "nHSP")]))))
    stop("'ckmr' must not contain missing values")
  if(any(x$nComp <= 0)) stop("'ckmr$nComp' must be positive")
  if(any(x$nPOP < 0 | x$nHSP < 0)) stop("kin counts must not be negative")
  if(any(x$nPOP > x$nComp | x$nHSP > x$nComp))
    stop("more kin pairs than comparisons in at least one row")
  class(x) <- c("ckmr_data", "data.frame")
  x
}

##' Internal: CKMR settings from a configuration, with defaults for
##' configurations predating the CKMR feature.
##' @keywords internal
##' @noRd
ckmrSettings <- function(conf){
  g <- function(nm, def) if(is.null(conf[[nm]])) def else conf[[nm]][1]
  list(usePOP = as.integer(g("usePOP", 0L)),
       useHSP = as.integer(g("useHSP", 0L)),
       psi    = as.numeric(g("ckmrPsi", 1.5)),
       scale  = as.numeric(g("ckmrScale", 1)),
       estPsi = as.integer(g("ckmrEstimatePsi", 0L)),
       plusExtra = as.integer(g("ckmrPlusExtra", 10L)),
       plusNodes = as.integer(g("ckmrPlusNodes", 6L)))
}

##' Build a close-kin mark-recapture pair table
##'
##' Turns what a genotyping study actually produces - who was sampled, and
##' which pairs came back as kin - into the comparison table consumed by
##' \code{\link{setup.sam.data}}.
##'
##' Individuals are grouped into (year, age) cells, and one row is produced per
##' pair of distinct cells with \code{nComp = n1 * n2} comparisons. This is what
##' keeps the table small: ten sampled years by fifteen ages is 150 cells and
##' about 11000 rows, however many thousands of individuals were genotyped.
##'
##' Pairs of animals from the \emph{same} cell are deliberately absent. They are
##' necessarily from the same cohort, which is structurally impossible for a
##' parent-offspring pair and is excluded for half-siblings (same-cohort
##' half-sibs are confounded with full sibs).
##'
##' @param samples a data frame with columns \code{year}, \code{age} and
##' \code{n}: the number of individuals of that age sampled in that year.
##' Duplicated cells are summed.
##' @param pop optional data frame of observed parent-offspring pairs with
##' columns \code{year1}, \code{age1}, \code{year2}, \code{age2}, \code{n}.
##' The two cells may be given in either order.
##' @param hsp optional data frame of observed half-sibling pairs, same layout.
##' @return a data frame of class \code{ckmr_data} with columns \code{year1},
##' \code{age1}, \code{year2}, \code{age2}, \code{nComp}, \code{nPOP},
##' \code{nHSP}.
##' @seealso \code{\link{simulateCKMR}}, \code{\link{setup.sam.data}}
##' @export
ckmrData <- function(samples, pop = NULL, hsp = NULL){
  need <- c("year", "age", "n")
  if(!all(need %in% names(samples)))
    stop("'samples' must have columns year, age and n")
  s <- stats::aggregate(list(n = as.numeric(samples$n)),
                        by = list(year = as.integer(samples$year),
                                  age  = as.integer(samples$age)), FUN = sum)
  s <- s[s$n > 0, , drop = FALSE]
  if(nrow(s) < 2) stop("need at least two non-empty (year, age) cells")
  s <- s[order(s$year, s$age), , drop = FALSE]
  ij <- utils::combn(nrow(s), 2)
  ret <- data.frame(year1 = s$year[ij[1, ]], age1 = s$age[ij[1, ]],
                    year2 = s$year[ij[2, ]], age2 = s$age[ij[2, ]],
                    nComp = s$n[ij[1, ]] * s$n[ij[2, ]],
                    nPOP  = 0, nHSP = 0)
  k <- ckmrKey(ret$year1, ret$age1, ret$year2, ret$age2)
  addKin <- function(tab, what){
    if(is.null(tab) || nrow(tab) == 0) return(invisible(NULL))
    if(!all(c("year1", "age1", "year2", "age2", "n") %in% names(tab)))
      stop("'", what, "' must have columns year1, age1, year2, age2, n")
    kk <- ckmrKey(tab$year1, tab$age1, tab$year2, tab$age2)
    m  <- match(kk, k)
    if(any(is.na(m)))
      stop(sum(is.na(m)), " ", what,
           " pair(s) refer to (year, age) cells that are not in 'samples'")
    tapply(as.numeric(tab$n), m, sum)
  }
  aP <- addKin(pop, "pop"); if(!is.null(aP)) ret$nPOP[as.integer(names(aP))] <- aP
  aH <- addKin(hsp, "hsp"); if(!is.null(aH)) ret$nHSP[as.integer(names(aH))] <- aH
  class(ret) <- c("ckmr_data", "data.frame")
  ret
}

##' Simulate close-kin mark-recapture data from a fitted model
##'
##' Draws a sampling design (how many individuals of each age are genotyped in
##' each year), forms all cross-cell comparisons, and simulates the number of
##' parent-offspring and half-sibling pairs found. The kin probabilities come
##' from exactly the same engine the likelihood uses, so a data set simulated
##' here and then fitted is internally consistent by construction.
##'
##' Counts are drawn from the Poisson pseudo-likelihood itself. That makes the
##' function well suited to study design and to \code{\link{simstudy}}, but it
##' does not reproduce the dependence between overlapping pairs that a full
##' pedigree simulation would.
##'
##' @param fit a fitted \code{sam} object.
##' @param years years in which individuals are genotyped.
##' @param n number genotyped per year (recycled to the length of \code{years}).
##' @param selectivity age composition of the genotyped sample. \code{NULL}
##' (default) uses the model-predicted catch at age (\code{\link{caytable}}) for
##' each year; a vector is used for every year; a matrix is taken as one row per
##' element of \code{years}.
##' @param seed optional random seed.
##' @return a \code{ckmr_data} pair table, ready to pass to
##' \code{\link{setup.sam.data}} as \code{ckmr}.
##' @seealso \code{\link{ckmrData}}, \code{\link{ckmrtable}}
##' @export
simulateCKMR <- function(fit, years, n, selectivity = NULL, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  years <- as.integer(years)
  n     <- rep(as.numeric(n), length.out = length(years))
  modYears <- fit$data$years
  ages     <- fit$conf$minAge:fit$conf$maxAge
  if(!all(years %in% modYears))
    stop("sampling years outside the model year range: ",
         paste(setdiff(years, modYears), collapse = ", "))
  set <- ckmrSettings(fit$conf)
  if(set$estPsi == 1L && length(fit$pl$logPsim1) > 0)
    set$psi <- exp(fit$pl$logPsim1[1]) + 1
  sel <- if(is.null(selectivity)){
           caytable(fit)[match(years, modYears), , drop = FALSE]
         }else if(is.matrix(selectivity)){
           if(nrow(selectivity) != length(years))
             stop("'selectivity' must have one row per sampling year")
           selectivity
         }else{
           matrix(rep(as.numeric(selectivity), each = length(years)),
                  nrow = length(years))
         }
  if(ncol(sel) != length(ages))
    stop("'selectivity' must have one column per model age")
  sel[!is.finite(sel) | sel < 0] <- 0
  samples <- do.call(rbind, lapply(seq_along(years), function(i){
    if(sum(sel[i, ]) <= 0) stop("selectivity sums to zero in year ", years[i])
    data.frame(year = years[i], age = ages,
               n = as.numeric(stats::rmultinom(1, n[i], sel[i, ])[, 1]))
  }))
  ck <- ckmrData(samples)

  prep <- ckmrPrep(ck, modYears, ages,
                   propMat = if(is.null(fit$conf$matureModel) ||
                                fit$conf$matureModel == 0) fit$data$propMat else NULL,
                   plusExtra = set$plusExtra, plusNodes = set$plusNodes)
  N   <- ntable(fit) * set$scale
  Z   <- faytable(fit) + fit$data$natMor
  swPos <- fit$data$stockMeanWeight[is.finite(fit$data$stockMeanWeight) &
                                    fit$data$stockMeanWeight > 0]
  swRef <- if(length(swPos) > 0) exp(mean(log(swPos))) else 1
  fec <- fit$data$propMat * (fit$data$stockMeanWeight / swRef)^set$psi
  pr  <- ckmrProb(N, Z, fec, prep)

  ck$nPOP[prep$popRow] <- stats::rpois(length(prep$popRow),
                                       pr$pPOP[prep$popRow] * ck$nComp[prep$popRow])
  ck$nHSP[prep$hspRow] <- stats::rpois(length(prep$hspRow),
                                       pr$pHSP[prep$hspRow] * ck$nComp[prep$hspRow])
  attr(ck, "pPOP") <- pr$pPOP
  attr(ck, "pHSP") <- pr$pHSP
  ck
}

##' Observed and expected close-kin pair counts
##'
##' Returns the CKMR pair table augmented with the fitted expected number of
##' parent-offspring and half-sibling pairs and the corresponding Pearson
##' residuals. The \code{usePOP} / \code{useHSP} columns flag which rows the
##' likelihood actually used - rows are excluded when the offspring's birth year
##' falls before the first model year, when the candidate parent was sampled
##' before the offspring was born, when no candidate parent age is mature, and
##' (for half-siblings) when the two animals are from the same cohort.
##'
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return a data frame with the pair table plus \code{ePOP}, \code{rPOP},
##' \code{usePOP}, \code{eHSP}, \code{rHSP}, \code{useHSP}.
##' @export
ckmrtable <- function(fit, ...) UseMethod("ckmrtable")

##' @rdname ckmrtable
##' @method ckmrtable sam
##' @export
ckmrtable.sam <- function(fit, ...){
  ck <- fit$data$ckmr
  if(is.null(ck)) stop("this fit carries no CKMR data")
  pP <- fit$rep$pPOP; pH <- fit$rep$pHSP
  if(is.null(pP)) stop("no reported CKMR probabilities - was the fit run with usePOP/useHSP?")
  prep <- fit$dat$ckmrPrep
  ret <- as.data.frame(ck[, c("year1", "age1", "year2", "age2", "nComp",
                              "nPOP", "nHSP")])
  ret$ePOP   <- pP * ck$nComp
  ret$rPOP   <- (ret$nPOP - ret$ePOP) / sqrt(pmax(ret$ePOP, .Machine$double.eps))
  ret$usePOP <- seq_len(nrow(ret)) %in% prep$popRow
  ret$eHSP   <- pH * ck$nComp
  ret$rHSP   <- (ret$nHSP - ret$eHSP) / sqrt(pmax(ret$eHSP, .Machine$double.eps))
  ret$useHSP <- seq_len(nrow(ret)) %in% prep$hspRow
  ret$rPOP[!ret$usePOP] <- NA
  ret$rHSP[!ret$useHSP] <- NA
  ret
}

##' Diagnostic plot of observed against expected close-kin pair counts
##'
##' Individual (year, age) cell pairs carry very few kin, so counts are
##' aggregated before plotting. Observed counts are shown with approximate
##' Poisson intervals against the fitted expectation.
##'
##' @param fit a fitted \code{sam} object.
##' @param by grouping for the aggregation: \code{"gap"} (default) groups by the
##' number of years between the two cohorts, \code{"year"} by the earlier
##' birth year.
##' @param ... passed to \code{\link[graphics]{plot}}.
##' @return invisibly, the aggregated data frame.
##' @export
ckmrplot <- function(fit, by = c("gap", "year"), ...) UseMethod("ckmrplot")

##' @rdname ckmrplot
##' @method ckmrplot sam
##' @export
ckmrplot.sam <- function(fit, by = c("gap", "year"), ...){
  by <- match.arg(by)
  tab <- ckmrtable(fit)
  c1 <- tab$year1 - tab$age1; c2 <- tab$year2 - tab$age2
  grp <- if(by == "gap") abs(c2 - c1) else pmin(c1, c2)
  set <- ckmrSettings(fit$conf)
  panels <- c(if(set$usePOP == 1) "POP", if(set$useHSP == 1) "HSP")
  if(length(panels) == 0) panels <- c("POP", "HSP")
  opar <- par(mfrow = c(1, length(panels))); on.exit(par(opar))
  out <- list()
  for(p in panels){
    use <- tab[[paste0("use", p)]]
    if(!any(use)){ plot.new(); next }
    o <- tapply(tab[[paste0("n", p)]][use], grp[use], sum)
    e <- tapply(tab[[paste0("e", p)]][use], grp[use], sum)
    g <- as.numeric(names(o))
    lo <- pmax(0, o - 2 * sqrt(o)); hi <- o + 2 * sqrt(o)
    plot(g, o, type = "n", ylim = range(0, hi, e),
         xlab = if(by == "gap") "cohort gap (years)" else "birth year",
         ylab = paste("number of", p, "pairs"),
         main = paste(p, "- observed vs expected"), ...)
    grid()
    for(i in seq_along(g)) lines(rep(g[i], 2), c(lo[i], hi[i]), col = "darkgrey")
    points(g, o, pch = 16)
    lines(g, e, col = "darkred", lwd = 2)
    legend("topright", c("observed", "expected"), pch = c(16, NA),
           lty = c(NA, 1), col = c("black", "darkred"), lwd = c(NA, 2), bty = "n")
    out[[p]] <- data.frame(group = g, observed = as.numeric(o),
                           expected = as.numeric(e))
  }
  invisible(out)
}
