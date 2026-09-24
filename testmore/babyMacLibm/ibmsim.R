## Individual-based pedigree simulation of a LENGTH-based close-kin data set.
##
## testmore/babyMibm/ibmsim.R is the starting point; what is new here is that
## every individual carries a length. It is drawn once, at birth, as a
## multiplicative deviation from the von Bertalanffy mean, and the individual
## keeps that deviation for life. That is exactly what the length-based
## probabilities assume, both when they retroject a parent's length back to
## the birth year and when they pair the unobserved shared parent across two
## birth years on one quantile. So this simulation tests whether that
## assumption is implemented correctly - not whether the assumption is true of
## real fish.
##
## Reproductive output follows the individual's own length,
## MO[y, a] * (l / lref)^psi, rather than an age average, so a big fish of a
## given age really does breed more than a small one of the same age.
##
## Nothing in this file evaluates a close-kin probability, so comparing the
## kin it counts against ckmrlProb is an independent check of the formulas -
## something simulateCKMRL cannot give, since it draws from them.

##' Self-consistent abundance: fitted recruitment, deterministic survival.
ibmConsistentN <- function(N, Z){
  nY <- nrow(N); nA <- ncol(N); Nc <- N
  for(y in 2:nY){
    for(a in 2:nA) Nc[y, a] <- Nc[y-1, a-1] * exp(-Z[y-1, a-1])
    Nc[y, nA] <- Nc[y, nA] + Nc[y-1, nA] * exp(-Z[y-1, nA])
  }
  Nc
}

##' Simulate the population, its pedigree and every individual's length.
##'
##' @param Nc year-by-age abundance on the individual scale.
##' @param Z year-by-age total mortality.
##' @param MO year-by-age proportion mature.
##' @param Lmean mean length at each model age.
##' @param sdLog log-normal spread of length at age.
##' @param refLength reference length in the fecundity relation (it cancels).
##' @param psi length exponent in fecundity.
##' @param minAge first modelled age.
##' @return list(mom, dad, coh, dev, aliveBy); 0 marks an unknown parent and
##' dev is the individual's lifelong length-at-age deviation.
ibmPopulationL <- function(Nc, Z, MO, Lmean, sdLog, refLength, psi, minAge){
  nY <- nrow(Nc); nA <- ncol(Nc)
  Ni <- round(Nc); Ni[Ni < 0] <- 0
  total <- sum(Ni[1, ]) + sum(Ni[, 1]) + 10L
  mom <- integer(total); dad <- integer(total); coh <- integer(total)
  sx  <- integer(total); dev <- numeric(total); nxt <- 0L
  newIds <- function(k){
    if(k <= 0) return(integer(0))
    i <- seq.int(nxt + 1L, length.out = k); nxt <<- nxt + k; i
  }
  ## year 1: standing stock, parents unknown (coded 0)
  alive <- vector("list", nA)
  for(a in seq_len(nA)){
    id <- newIds(Ni[1, a])
    if(length(id)){
      coh[id] <- 1L - a
      sx[id]  <- rep(1:2, length.out = length(id))
      dev[id] <- stats::rlnorm(length(id), 0, sdLog)
    }
    alive[[a]] <- id
  }
  pending <- vector("list", nY + minAge + 1L)
  aliveBy <- vector("list", nY)

  for(y in seq_len(nY)){
    aliveBy[[y]] <- alive
    yRec <- y + minAge
    nB <- if(yRec <= nY) Ni[yRec, 1] else 0L
    if(nB > 0){
      pool <- unlist(alive)
      if(length(pool)){
        aPool <- rep(seq_len(nA), lengths(alive))
        ## each individual's own reproductive output: mature for its age, and
        ## fecundity at the length it actually carries
        w <- MO[y, aPool] * (Lmean[aPool] * dev[pool] / refLength)^psi
        w[!is.finite(w)] <- 0
        isF <- sx[pool] == 1L
        okF <- isF & w > 0; okM <- !isF & w > 0
        if(any(okF) && any(okM)){
          kid <- newIds(nB)
          coh[kid] <- y
          sx[kid]  <- rep(1:2, length.out = nB)
          dev[kid] <- stats::rlnorm(nB, 0, sdLog)
          mom[kid] <- sample(pool[okF], nB, replace = TRUE, prob = w[okF])
          dad[kid] <- sample(pool[okM], nB, replace = TRUE, prob = w[okM])
          pending[[yRec]] <- kid
        }
      }
    }
    if(y == nY) break
    ## survive and age: exactly the deterministic numbers, random individuals
    nw <- vector("list", nA)
    for(a in seq_len(nA - 1L)){
      k <- length(alive[[a]]); if(k == 0) next
      ns <- min(k, round(Nc[y, a] * exp(-Z[y, a])))
      if(ns > 0){
        s <- alive[[a]][sample.int(k, ns)]
        tgt <- min(a + 1L, nA)
        nw[[tgt]] <- c(nw[[tgt]], s)
      }
    }
    k <- length(alive[[nA]])                       # plus group survivors
    if(k > 0){
      ns <- min(k, round(Nc[y, nA] * exp(-Z[y, nA])))
      if(ns > 0) nw[[nA]] <- c(nw[[nA]], alive[[nA]][sample.int(k, ns)])
    }
    nw[[1]] <- if(is.null(pending[[y + 1L]])) integer(0) else pending[[y + 1L]]
    alive <- lapply(nw, function(z) if(is.null(z)) integer(0) else z)
  }
  list(mom = mom[seq_len(nxt)], dad = dad[seq_len(nxt)],
       coh = coh[seq_len(nxt)], dev = dev[seq_len(nxt)], aliveBy = aliveBy)
}

##' Draw a genotyped sample and measure it.
##'
##' Ages are drawn from the same composition the likelihood assumes the sample
##' follows - the model's catch at age - because P[age | length, year] is built
##' from it. Sampling is lethal, as in the age-based simulation: an individual
##' taken in one year is removed from the pool, since an adult cannot parent
##' anything born after it was caught.
##'
##' Each sampled fish is measured at its own length, Lmean[age] * dev, and
##' binned. Fish outside the bin range are dropped and counted.
ibmSampleL <- function(pop, years, ages, sampleYears, nPerYear, selectivity,
                       Lmean, lengthBins){
  minYear <- min(years); out <- NULL; taken <- integer(0)
  nPerYear <- rep(nPerYear, length.out = length(sampleYears))
  for(i in seq_along(sampleYears)){
    y <- sampleYears[i] - minYear + 1L
    av <- pop$aliveBy[[y]]
    sel <- if(is.matrix(selectivity)) selectivity[i, ] else selectivity
    want <- as.numeric(stats::rmultinom(1, nPerYear[i], sel)[, 1])
    for(a in seq_along(ages)){
      avail <- av[[a]]
      if(length(taken)) avail <- avail[!(avail %in% taken)]
      k <- min(length(avail), round(want[a]))
      if(k <= 0) next
      pick <- avail[sample.int(length(avail), k)]
      taken <- c(taken, pick)
      len <- Lmean[a] * pop$dev[pick]
      out <- rbind(out, data.frame(id = pick, year = sampleYears[i],
                                   age = ages[a], len = len,
                                   cls = cut(len, lengthBins, labels = FALSE)))
    }
  }
  attr(out, "offGrid") <- sum(is.na(out$cls))
  out[!is.na(out$cls), , drop = FALSE]
}

##' Count the kin actually present among the sampled fish and assemble the
##' length-based comparison table.
##'
##' Half-siblings are pairs sharing exactly one parent: shared-mother plus
##' shared-father pairs, less twice the full siblings. Unlike the age-based
##' table, pairs from the same (year, length) cell are kept - two fish of the
##' same length in one year need not be the same age, so such a pair can be
##' parent-offspring or a cross-cohort half-sibling.
ibmKinL <- function(pop, samp, mids){
  mom <- pop$mom[samp$id]; dad <- pop$dad[samp$id]
  cell <- paste(samp$year, samp$cls, sep = ".")
  uc   <- sort(unique(cell)); ci <- match(cell, uc)

  pos    <- match(c(mom, dad), samp$id)
  hit    <- which(!is.na(pos))
  offIdx <- rep(seq_len(nrow(samp)), 2)[hit]
  popA   <- ci[pos[hit]]; popB <- ci[offIdx]

  sharedPairs <- function(par){
    ok <- par > 0
    if(!any(ok)) return(NULL)
    s <- split(seq_along(par)[ok], par[ok])
    s <- s[lengths(s) >= 2L]
    if(!length(s)) return(NULL)
    do.call(rbind, lapply(s, function(k){
      ij <- utils::combn(length(k), 2L)
      cbind(ci[k[ij[1, ]]], ci[k[ij[2, ]]])
    }))
  }
  key <- function(m) if(is.null(m)) character(0)
                     else paste(pmin(m[, 1], m[, 2]), pmax(m[, 1], m[, 2]), sep = "|")
  tab <- function(k) if(!length(k)) integer(0) else table(k)
  km <- tab(key(sharedPairs(mom)))
  kd <- tab(key(sharedPairs(dad)))
  kf <- tab(key(sharedPairs(ifelse(mom > 0 & dad > 0, mom * 1e9 + dad, 0))))
  nms <- unique(c(names(km), names(kd), names(kf)))
  g <- function(t, n) if(length(t) == 0) rep(0, length(n))
                      else ifelse(is.na(t[n]), 0, as.numeric(t[n]))
  hspCnt <- setNames(g(km, nms) + g(kd, nms) - 2 * g(kf, nms), nms)
  popCnt <- tab(key(cbind(popA, popB)))

  splitKey <- function(nms){
    p  <- do.call(rbind, strsplit(nms, "|", fixed = TRUE))
    c1 <- uc[as.integer(p[, 1])]; c2 <- uc[as.integer(p[, 2])]
    q1 <- do.call(rbind, strsplit(c1, ".", fixed = TRUE))
    q2 <- do.call(rbind, strsplit(c2, ".", fixed = TRUE))
    data.frame(year1 = as.integer(q1[, 1]), len1 = mids[as.integer(q1[, 2])],
               year2 = as.integer(q2[, 1]), len2 = mids[as.integer(q2[, 2])])
  }
  mkTab <- function(cnt){
    if(!length(cnt)) return(NULL)
    d <- splitKey(names(cnt)); d$n <- as.numeric(cnt); d
  }
  cells <- data.frame(year = as.integer(sub("[.].*", "", uc)),
                      len  = mids[as.integer(sub(".*[.]", "", uc))],
                      n    = as.numeric(table(cell)[uc]))
  ckmrlData(cells, pop = mkTab(popCnt), hsp = mkTab(hspCnt[hspCnt > 0]))
}
