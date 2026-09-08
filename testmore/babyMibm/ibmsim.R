## Individual-based pedigree simulation of a close-kin data set.
##
## Builds an explicit population of individuals, each with a recorded mother
## and father, then counts the parent-offspring and half-sibling pairs that
## actually exist among a sample of them. Nothing in this file uses the
## close-kin probability formulas, so fitting - or simply comparing rates to -
## the data it produces is an independent check of those formulas.
##
## The population is built to satisfy the survival recursion exactly:
##   Nc[y+1, a+1] = Nc[y, a] * exp(-Z[y, a]),  plus group accumulating,
## seeded with the assessment's recruitment and first year. SAM's fitted N
## carries process noise and does not satisfy that recursion (realised
## survival ratios run from 0.8 to 5.6, compounding to four orders of
## magnitude over a cohort track), and the close-kin probabilities mix N and
## Z, so an IBM built on the fitted N directly would not be testing them.

##' Self-consistent abundance: fitted recruitment, deterministic survival.
ibmConsistentN <- function(N, Z){
  nY <- nrow(N); nA <- ncol(N); Nc <- N
  for(y in 2:nY){
    for(a in 2:nA) Nc[y, a] <- Nc[y-1, a-1] * exp(-Z[y-1, a-1])
    Nc[y, nA] <- Nc[y, nA] + Nc[y-1, nA] * exp(-Z[y-1, nA])
  }
  Nc
}

##' Simulate the population and its pedigree.
##'
##' @param Nc year-by-age abundance on the individual scale (self-consistent).
##' @param Z   year-by-age total mortality.
##' @param fec year-by-age per-capita reproductive output.
##' @param minAge first modelled age.
##' @return list(mom, dad, coh, aliveBy) with 0 for unknown parents.
ibmPopulation <- function(Nc, Z, fec, minAge){
  nY <- nrow(Nc); nA <- ncol(Nc)
  Ni <- round(Nc); Ni[Ni < 0] <- 0
  total <- sum(Ni[1, ]) + sum(Ni[, 1]) + 10L
  mom <- integer(total); dad <- integer(total); coh <- integer(total)
  sx  <- integer(total); nxt <- 0L
  newIds <- function(k){
    if(k <= 0) return(integer(0))
    i <- seq.int(nxt + 1L, length.out = k); nxt <<- nxt + k; i
  }
  ## year 1: standing stock, parents unknown (coded 0)
  alive <- vector("list", nA)
  for(a in seq_len(nA)){
    id <- newIds(Ni[1, a])
    if(length(id)){ coh[id] <- 1L - a; sx[id] <- rep(1:2, length.out = length(id)) }
    alive[[a]] <- id
  }
  pending <- vector("list", nY + minAge + 1L)   # cohort waiting to recruit
  aliveBy <- vector("list", nY)

  for(y in seq_len(nY)){
    aliveBy[[y]] <- alive
    ## individuals born in row-year y recruit at age index 1 in row y + minAge
    yRec <- y + minAge
    nB <- if(yRec <= nY) Ni[yRec, 1] else 0L
    if(nB > 0){
      pool <- unlist(alive)
      if(length(pool)){
        w <- fec[y, rep(seq_len(nA), lengths(alive))]
        isF <- sx[pool] == 1L
        okF <- isF & w > 0; okM <- !isF & w > 0
        if(any(okF) && any(okM)){
          kid <- newIds(nB)
          coh[kid] <- y
          sx[kid]  <- rep(1:2, length.out = nB)
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
       coh = coh[seq_len(nxt)], aliveBy = aliveBy)
}

##' Draw a genotyped sample, stratified by age as in simulateCKMR's
##' selectivity vector, and return one row per sampled individual.
##'
##' Sampling is lethal: an individual taken in one year is removed from the
##' pool and cannot be taken again. Without that, a fish still alive the next
##' year can be caught twice and appear in the table at two different ages,
##' which contradicts the lethal-sampling assumption the close-kin
##' probabilities are built on (an adult cannot parent anything born after it
##' was sampled) and undercounts kin, since a duplicated parent is only
##' matched at its first occurrence.
ibmSample <- function(pop, years, ages, sampleYears, nPerYear, selectivity){
  minYear <- min(years); out <- NULL; taken <- integer(0)
  nPerYear <- rep(nPerYear, length.out = length(sampleYears))
  for(i in seq_along(sampleYears)){
    y <- sampleYears[i] - minYear + 1L
    av <- pop$aliveBy[[y]]
    want <- as.numeric(stats::rmultinom(1, nPerYear[i], selectivity)[, 1])
    for(a in seq_along(ages)){
      avail <- av[[a]]
      if(length(taken)) avail <- avail[!(avail %in% taken)]
      k <- min(length(avail), round(want[a]))
      if(k <= 0) next
      pick <- avail[sample.int(length(avail), k)]
      taken <- c(taken, pick)
      out <- rbind(out, data.frame(id = pick, year = sampleYears[i],
                                   age = ages[a]))
    }
  }
  out
}

##' Count the parent-offspring and half-sibling pairs actually present among
##' the sampled individuals, and assemble the close-kin comparison table.
##'
##' Half-siblings are pairs sharing exactly one parent: shared-mother pairs
##' plus shared-father pairs, less twice the full siblings, which share both.
##' Pairs whose two members fall in the same (year, age) cell are necessarily
##' from the same cohort - impossible for a parent-offspring pair, and
##' excluded for half-siblings - and are reported separately rather than
##' entered in the table.
ibmKin <- function(pop, samp){
  mom <- pop$mom[samp$id]; dad <- pop$dad[samp$id]
  cell <- paste(samp$year, samp$age, sep = ".")
  uc   <- sort(unique(cell)); ci <- match(cell, uc)

  ## parent-offspring: a sampled individual whose mother or father is also sampled
  pos    <- match(c(mom, dad), samp$id)
  hit    <- which(!is.na(pos))
  offIdx <- rep(seq_len(nrow(samp)), 2)[hit]
  popA   <- ci[pos[hit]]; popB <- ci[offIdx]

  ## pairs of sampled individuals sharing a given parent
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
  kf <- tab(key(sharedPairs(ifelse(mom > 0 & dad > 0,
                                   mom * 1e9 + dad, 0))))
  nms <- unique(c(names(km), names(kd), names(kf)))
  g <- function(t, n) if(length(t) == 0) rep(0, length(n))
                      else ifelse(is.na(t[n]), 0, as.numeric(t[n]))
  hspCnt <- setNames(g(km, nms) + g(kd, nms) - 2 * g(kf, nms), nms)
  popCnt <- tab(key(cbind(popA, popB)))

  ## split the cell keys back into the two cells
  splitKey <- function(nms){
    p <- do.call(rbind, strsplit(nms, "|", fixed = TRUE))
    c1 <- uc[as.integer(p[, 1])]; c2 <- uc[as.integer(p[, 2])]
    q1 <- do.call(rbind, strsplit(c1, ".", fixed = TRUE))
    q2 <- do.call(rbind, strsplit(c2, ".", fixed = TRUE))
    data.frame(year1 = as.integer(q1[, 1]), age1 = as.integer(q1[, 2]),
               year2 = as.integer(q2[, 1]), age2 = as.integer(q2[, 2]),
               same  = c1 == c2)
  }
  mkTab <- function(cnt){
    if(!length(cnt)) return(list(tab = NULL, dropped = 0))
    d <- splitKey(names(cnt)); d$n <- as.numeric(cnt)
    list(tab = d[!d$same, c("year1","age1","year2","age2","n"), drop = FALSE],
         dropped = sum(d$n[d$same]))
  }
  P <- mkTab(popCnt); H <- mkTab(hspCnt[hspCnt > 0])
  cells <- data.frame(year = as.integer(sub("[.].*", "", uc)),
                      age  = as.integer(sub(".*[.]", "", uc)),
                      n    = as.numeric(table(cell)[uc]))
  ck <- ckmrData(cells, pop = P$tab, hsp = H$tab)
  attr(ck, "withinCellPOP") <- P$dropped
  attr(ck, "withinCellHSP") <- H$dropped
  ck
}
