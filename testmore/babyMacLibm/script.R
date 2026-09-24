## As testmore/babyMacL, but the close-kin data comes from an individual-based
## pedigree simulation (ibmsim.R, local to this folder) rather than from
## simulateCKMRL. The IBM builds explicit individuals with recorded mothers,
## fathers and lengths, and counts the kin pairs that actually exist among a
## sample of them. It never evaluates a close-kin probability, so comparing
## its counts against ckmrlProb is an independent check of the length-based
## formulas - which simulateCKMRL cannot give, since it draws from them and is
## self-consistent by construction.
##
## Two things have to be right for the comparison to mean anything, both
## carried over from babyMibm:
##  - the population must satisfy the survival recursion, because the
##    probabilities mix N and Z and SAM's fitted N carries process noise;
##  - sampling must be lethal, or a fish caught twice appears at two lengths
##    and kin are undercounted.
## A third is new here: ages are drawn from the population, because the
## length-based probabilities assume selection acts on length and year, not on
## age, so that conditioning on the measured length absorbs it. Drawing the
## sample with age-based selection would violate that assumption rather than
## test the formula.
library(samjr)

cn<-read.ices("cn.dat"); cw<-read.ices("cw.dat"); dw<-read.ices("dw.dat")
lf<-read.ices("lf.dat"); lw<-read.ices("lw.dat"); mo<-read.ices("mo.dat")
nm<-read.ices("nm.dat"); pf<-read.ices("pf.dat"); pm<-read.ices("pm.dat")
sw<-read.ices("sw.dat"); surveys<-read.ices("survey.dat")
recap<-read.table("tag.dat", header=TRUE)

W<-matrix(NA, nrow=nrow(cn), ncol=ncol(cn))
W[as.numeric(rownames(cn))<2000]<-10
attr(cn,"weight")<-W

dat<-setup.sam.data(surveys=surveys, residual.fleet=cn, prop.mature=mo,
                    stock.mean.weight=sw, catch.mean.weight=cw,
                    dis.mean.weight=dw, land.mean.weight=lw,
                    prop.f=pf, prop.m=pm, natural.mortality=nm,
                    land.frac=lf, recapture=recap)

mkconf<-function(dat){
  conf<-defcon(dat)
  conf$fbarRange <- c(4,8)
  conf$corFlag <- 2
  conf$fixVarToWeight <- 1
  conf$keyLogFsta <- matrix(c(0,1,2,3,4,5,6,7,8,9,9, rep(-1,55)), nrow=6, byrow=TRUE)
  conf$keyLogFpar <- matrix(c(
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,   0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 1, 2, 3, 4, 5, 6, 7, 8, 8,-1,   9,10,11,12,13,14,15,16,17,17,-1,
    18,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1),
    nrow=6, byrow=TRUE)
  conf$keyVarF <- matrix(c(0,0,1,1,1,1,1,1,1,2,2, rep(-1,55)), nrow=6, byrow=TRUE)
  conf$keyVarObs <- matrix(c(
     0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2,   3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 5, 6, 6, 6, 6, 6, 6, 7, 7,-1,   8, 9,10,10,10,10,10,10,10,10,-1,
     4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1),
    nrow=6, byrow=TRUE)
  conf$keyCorObs <- matrix(c(rep(NA,10), rep(-1,10), rep(0,9),-1, rep(1,9),-1,
                             rep(-1,10), rep(-1,10)), nrow=6, byrow=TRUE)
  conf$obsCorStruct <- factor(c("ID","ID","AR","AR","ID","ID"),
                              levels=c("ID","AR","US"))
  conf
}
conf0<-mkconf(dat)
fit0<-sam.fit(dat, conf0, defpar(dat, conf0), silent=TRUE)

source("ibmsim.R")

## ---- individual-based pedigree simulation ----
## The true stock is far too large to simulate one fish at a time, so the
## population is scaled down and ckmrlScale set to match: model numbers-at-age
## times the scale are exactly the simulated individuals.
SCALE <- 0.01
NBIN  <- 12
set   <- samjr:::ckmrlSettings(conf0)   # Linf 41.7, k 0.35, sd 0.1, ref 50, psi 4.5
PSI   <- set$psi
years <- fit0$data$years
ages  <- fit0$conf$minAge:fit0$conf$maxAge
Zmat  <- faytable(fit0) + fit0$data$natMor
MOmat <- fit0$data$propMat
G     <- samjr:::ckmrlGrid(set, ages)
Nc    <- ibmConsistentN(ntable(fit0) * SCALE, Zmat)

lengthBins <- seq(exp(qnorm(0.001, log(min(G$Lmean)), G$sdL)),
                  exp(qnorm(0.999, log(max(G$Lmean)), G$sdL)), length = NBIN + 1L)
mids <- lengthBins[-1] - 0.5*diff(lengthBins)

set.seed(20)
pop  <- ibmPopulationL(Nc, Zmat, MOmat, G$Lmean, set$sdLogLength,
                       set$refLength, PSI, min(ages))
## Ages drawn from the POPULATION, not the catch. The probabilities assume
## selection acts on length and year rather than on age, so that conditioning
## on the measured length absorbs it; a sample drawn with age-based selection
## would violate that assumption instead of testing the formula.
cay  <- Nc; cay[!is.finite(cay) | cay < 0] <- 0
sampleYears <- 2021:2023
samp <- ibmSampleL(pop, years, ages, sampleYears, rep(6000, 3),
                   cay[match(sampleYears, years), , drop=FALSE], G$Lmean, lengthBins)
ckl  <- ibmKinL(pop, samp, mids)

## ---- does the pedigree agree with the probability formulas? ----
prep <- samjr:::ckmrlPrep(ckl, years, ages, set, propMat=MOmat)
prob <- samjr:::ckmrlProb(Nc, Zmat, MOmat, Nc, PSI, prep)
rat <- function(col, rows, p){
  if(!length(rows)) return(NA)
  sum(ckl[[col]][rows]) / sum(p[rows]*ckl$nComp[rows])
}
same  <- ckl$year1==ckl$year2 & ckl$len1==ckl$len2
big   <- pmax(ckl$len1, ckl$len2) >= mids[ceiling(2*NBIN/3)]
sub   <- function(rows, f) rows[f[rows]]
r.pop      <- rat("nPOP", prep$popRow, prob$pPOP)
r.pop.same <- rat("nPOP", sub(prep$popRow,  same), prob$pPOP)
r.pop.big  <- rat("nPOP", sub(prep$popRow,  big),  prob$pPOP)
r.pop.sml  <- rat("nPOP", sub(prep$popRow, !big),  prob$pPOP)
r.hsp      <- rat("nHSP", prep$hspRow, prob$pHSP)
r.hsp.same <- rat("nHSP", sub(prep$hspRow,  same), prob$pHSP)
r.hsp.big  <- rat("nHSP", sub(prep$hspRow,  big),  prob$pHSP)
r.hsp.sml  <- rat("nHSP", sub(prep$hspRow, !big),  prob$pHSP)
## how much age information the length classes actually carry
ageSpread <- mean(tapply(samp$age, samp$cls, function(a) diff(range(a))))

## ---- observed against predicted, grouped by predicted size ----
png("babyMacLibm.png", width=1100, height=520, pointsize=15,
    type=if(capabilities("cairo")) "cairo" else "Xlib")
opar<-par(mfrow=c(1,2), mar=c(4,4,3,1))
for(w in c("POP","HSP")){
  rows <- if(w=="POP") prep$popRow else prep$hspRow
  p    <- if(w=="POP") prob$pPOP   else prob$pHSP
  obs  <- ckl[[paste0("n",w)]][rows]; exp1 <- p[rows]*ckl$nComp[rows]
  grp  <- cut(exp1, quantile(exp1, seq(0,1,length=9)), include.lowest=TRUE)
  o <- tapply(obs, grp, sum); e <- tapply(exp1, grp, sum)
  lo <- pmax(0, o-2*sqrt(o)); hi <- o+2*sqrt(o)
  plot(e, o, xlim=range(0,e,hi), ylim=range(0,e,hi), pch=16,
       xlab="predicted pairs", ylab="observed pairs",
       main=paste(w, "- pedigree vs formula")); grid()
  for(i in seq_along(e)) lines(rep(e[i],2), c(lo[i],hi[i]), col="darkgrey")
  points(e, o, pch=16); abline(0, 1, col="darkred", lwd=2)
  legend("topleft", bty="n", lty=c(NA,1), pch=c(16,NA), col=c("black","darkred"),
         legend=c("observed (95%)","1:1"))
}
par(opar); dev.off()
plotok<-file.exists("babyMacLibm.png") && file.size("babyMacLibm.png") > 10000

writeLines(c(
  paste("base       ", format(fit0$opt$objective, digits=7)),
  paste("ibm.total  ", length(pop$mom)),
  paste("sampled    ", nrow(samp)),
  paste("offgrid    ", attr(samp,"offGrid")),
  paste("npairs     ", nrow(ckl)),
  paste("nPOP       ", sum(ckl$nPOP)),
  paste("nHSP       ", sum(ckl$nHSP)),
  paste("agespread  ", format(ageSpread, digits=3)),
  ## observed kin divided by kin predicted by ckmrlProb; 1 = formula is right
  paste("rat.pop    ", format(r.pop,      digits=4)),
  paste("rat.pop.sml", format(r.pop.sml,  digits=4)),
  paste("rat.pop.big", format(r.pop.big,  digits=4)),
  paste("rat.pop.sam", format(r.pop.same, digits=4)),
  paste("rat.hsp    ", format(r.hsp,      digits=4)),
  paste("rat.hsp.sml", format(r.hsp.sml,  digits=4)),
  paste("rat.hsp.big", format(r.hsp.big,  digits=4)),
  paste("rat.hsp.sam", format(r.hsp.same, digits=4)),
  paste("ok.pop     ", abs(r.pop-1) < 0.10),
  paste("ok.hsp     ", abs(r.hsp-1) < 0.10),
  paste("plot       ", plotok)
), "res.out")
