## As testmore/babyM, but the close-kin data comes from an individual-based
## pedigree simulation (ibmsim.R, local to this folder) rather than from
## simulateCKMR. The IBM builds explicit individuals with recorded mothers and
## fathers and counts the kin pairs that actually exist, so it never touches
## the close-kin probability formulas. Comparing its kin against what ckmrProb
## predicts is therefore an independent check of those formulas, which
## simulateCKMR - drawing from the formulas themselves - cannot provide.
## Deviations from ckmr/babyM.R, documented in ../STATUS.md:
##  - babyM.R also estimated a global natural-mortality multiplier in its two
##    HSP runs; samjr has no such parameter, so those runs differ by that one
##    degree of freedom. psi is estimated, as in babyM.R.
##  - the genotyped sample is uniform over ages, where babyM.R took ages 2-4
##    only. That narrow window leaves every candidate parent on a nearly flat
##    stretch of the fecundity schedule, which confounds psi with abundance and
##    all but removes the POP signal, and it never puts a candidate parent in
##    the plus group.
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
  ## numbers-at-age are in thousands, as in babyM.R's scale=1000
  conf$ckmrScale <- 1000
  conf$ckmrPsi <- 1.5
  conf$ckmrEstimatePsi <- 1
  conf
}

conf0<-mkconf(dat)
fit0<-sam.fit(dat, conf0, defpar(dat, conf0), silent=TRUE)


source("ibmsim.R")

## ---- individual-based pedigree simulation ----
## The true stock is ~8e11 individuals, far too many to simulate one by one,
## so the population is scaled down by 1e5 and conf$ckmrScale set to match:
## model numbers-at-age times ckmrScale are exactly the simulated individuals,
## so the model and the simulation describe the same population.
SCALE <- 0.01
PSI   <- 1.5
years <- fit0$data$years
ages  <- fit0$conf$minAge:fit0$conf$maxAge
maxA  <- max(ages)
Zmat  <- faytable(fit0) + fit0$data$natMor
swPos <- fit0$data$stockMeanWeight[is.finite(fit0$data$stockMeanWeight) &
                                   fit0$data$stockMeanWeight > 0]
fecM  <- fit0$data$propMat * (fit0$data$stockMeanWeight/exp(mean(log(swPos))))^PSI
Nc    <- ibmConsistentN(ntable(fit0) * SCALE, Zmat)

set.seed(20)
pop  <- ibmPopulation(Nc, Zmat, fecM, min(ages))
samp <- ibmSample(pop, years, ages, 2020:2023, rep(6000,4), rep(1,length(ages)))
ckmr <- ibmKin(pop, samp)

## ---- does the pedigree agree with the probability formulas? ----
prep <- samjr:::ckmrPrep(ckmr, years, ages, propMat=fit0$data$propMat,
                         plusExtra=conf0$ckmrPlusExtra,
                         plusNodes=conf0$ckmrPlusNodes)
prob <- samjr:::ckmrProb(Nc, Zmat, fecM, prep)
rat <- function(col, rows, p){
  if(!length(rows)) return(NA)
  sum(ckmr[[col]][rows]) / sum(p[rows]*ckmr$nComp[rows])
}
plusRow <- function(rows) rows[ckmr$age1[rows]==maxA | ckmr$age2[rows]==maxA]
noPlus  <- function(rows) rows[ckmr$age1[rows]< maxA & ckmr$age2[rows]< maxA]
r.pop        <- rat("nPOP", prep$popRow,        prob$pPOP)
r.pop.simple <- rat("nPOP", prep$popSimpleRow,  prob$pPOP)
r.pop.plus   <- rat("nPOP", prep$popPlusRow,    prob$pPOP)
r.hsp        <- rat("nHSP", prep$hspRow,        prob$pHSP)
r.hsp.noplus <- rat("nHSP", noPlus(prep$hspRow),prob$pHSP)
r.hsp.plus   <- rat("nHSP", plusRow(prep$hspRow),prob$pHSP)
## share of sampled plus-group fish whose recorded age is not their real age
misAged <- mean((pop$coh[samp$id] != (samp$year-samp$age-min(years)+1L))[samp$age==maxA])

## ---- the three babyM.R runs, on the pedigree-simulated data ----
dat2<-setup.sam.data(surveys=surveys, residual.fleet=cn, prop.mature=mo,
                     stock.mean.weight=sw, catch.mean.weight=cw,
                     dis.mean.weight=dw, land.mean.weight=lw,
                     prop.f=pf, prop.m=pm, natural.mortality=nm,
                     land.frac=lf, recapture=recap, ckmr=ckmr)
runit<-function(usePOP, useHSP){
  conf<-mkconf(dat2); conf$usePOP<-usePOP; conf$useHSP<-useHSP
  conf$ckmrScale<-SCALE
  par<-defpar(dat2, conf)
  for(nm in names(par))
    if(nm %in% names(fit0$pl) && length(par[[nm]])==length(fit0$pl[[nm]]))
      par[[nm]][]<-fit0$pl[[nm]]
  suppressWarnings(sam.fit(dat2, conf, par, silent=TRUE))
}
fitP <-runit(1,0); fitH <-runit(0,1); fitPH<-runit(1,1)
psiOf<-function(f) exp(f$pl$logPsim1[1])+1

kt<-function(x) exp(x)/1000
panel<-function(f, ttl){
  gv<-function(g,s) kt(g$sdrep$value[names(g$sdrep$value)=="logssb"] +
                       s*2*g$sdrep$sd[names(g$sdrep$value)=="logssb"])
  lo<-min(gv(fit0,-1), gv(f,-1)); hi<-max(gv(fit0,1), gv(f,1))
  plotit(fit0,"logssb",trans=kt,ylab="SSB (1000 t)",ylim=c(min(0,lo),hi),
         col="darkred",cicol=rgb(.55,0,0,.35))
  plotit(f,"logssb",add=TRUE,trans=kt,col="darkblue",cicol=rgb(0,0,.55,.35))
  title(ttl)
  legend("topleft", bty="n", border=NA,
         legend=c("official assessment (95%)","+ close kin (95%)",
                  paste0("psi = ", format(psiOf(f), digits=3))),
         fill=c(rgb(.55,0,0,.35), rgb(0,0,.55,.35), NA))
}
png("babyMibm.png", width=1500, height=520, pointsize=15,
    type=if(capabilities("cairo")) "cairo" else "Xlib")
opar<-par(mfrow=c(1,3), mar=c(4,4,3,1))
panel(fitP,"pedigree data + add pop"); panel(fitH,"pedigree data + add hsp")
panel(fitPH,"pedigree data + add pop & hsp")
par(opar); dev.off()
plotok<-file.exists("babyMibm.png") && file.size("babyMibm.png") > 10000

writeLines(c(
  paste("base       ", format(fit0$opt$objective, digits=7)),
  paste("ibm.alive  ", round(mean(rowSums(Nc)))),
  paste("ibm.total  ", length(pop$mom)),
  paste("sampled    ", nrow(samp)),
  paste("npairs     ", nrow(ckmr)),
  paste("nPOP       ", sum(ckmr$nPOP)),
  paste("nHSP       ", sum(ckmr$nHSP)),
  paste("wcellPOP   ", attr(ckmr,"withinCellPOP")),
  paste("wcellHSP   ", attr(ckmr,"withinCellHSP")),
  ## observed kin divided by kin predicted by ckmrProb; 1 = formula is right
  paste("rat.pop    ", format(r.pop,        digits=4)),
  paste("rat.pop.smp", format(r.pop.simple, digits=4)),
  paste("rat.pop.plu", format(r.pop.plus,   digits=4)),
  paste("rat.hsp    ", format(r.hsp,        digits=4)),
  paste("rat.hsp.nop", format(r.hsp.noplus, digits=4)),
  paste("rat.hsp.plu", format(r.hsp.plus,   digits=4)),
  paste("misaged12  ", format(misAged, digits=3)),
  ## the formulas must reproduce the pedigree wherever ages are known exactly
  paste("ok.pop     ", abs(r.pop.simple-1) < 0.10),
  paste("ok.hsp     ", abs(r.hsp.noplus-1) < 0.10),
  paste("pop        ", format(fitP$opt$objective,  digits=7)),
  paste("hsp        ", format(fitH$opt$objective,  digits=7)),
  paste("pophsp     ", format(fitPH$opt$objective, digits=7)),
  paste("psi.pop    ", format(psiOf(fitP),  digits=5)),
  paste("psi.hsp    ", format(psiOf(fitH),  digits=5)),
  paste("psi.pophsp ", format(psiOf(fitPH), digits=5)),
  paste("plot       ", plotok)
), "res.out")
