## Reproduces the three close-kin runs at the end of ckmr/babyM.R, on the
## mackerel assessment (model_XplConf.08) that script was built around:
##   official data + POP,  official data + HSP,  official data + POP & HSP
## Two deviations from the prototype, both documented in ../STATUS.md:
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

## babyM.R: simDat(years=2020:2023, n=rep(1e6,4)), but sampling uniformly
## over ages rather than babyM.R's selection=c(1,1,1,0,...)
## The seed is picked from a 20-seed scan as the most typical draw on the SSB
## shift and on psi, in all three configurations - not the most favourable one.
## The sampling distribution behind that choice is recorded in ../STATUS.md.
set.seed(20)
sel<-rep(1, 11)
ckmr<-simulateCKMR(fit0, years=2020:2023, n=rep(1000000,4), selectivity=sel)
dat2<-setup.sam.data(surveys=surveys, residual.fleet=cn, prop.mature=mo,
                     stock.mean.weight=sw, catch.mean.weight=cw,
                     dis.mean.weight=dw, land.mean.weight=lw,
                     prop.f=pf, prop.m=pm, natural.mortality=nm,
                     land.frac=lf, recapture=recap, ckmr=ckmr)

## Seed each close-kin run from the baseline assessment rather than from the
## all-zero default. The half-sibling likelihood is bimodal in psi here (babyM.R
## carried a natural-mortality multiplier that absorbed the second mode), and a
## cold start lands in the wrong basin.
runit<-function(usePOP, useHSP){
  conf<-mkconf(dat2); conf$usePOP<-usePOP; conf$useHSP<-useHSP
  par<-defpar(dat2, conf)
  for(nm in names(par))
    if(nm %in% names(fit0$pl) && length(par[[nm]])==length(fit0$pl[[nm]]))
      par[[nm]][]<-fit0$pl[[nm]]
  suppressWarnings(sam.fit(dat2, conf, par, silent=TRUE))
}
fitP <-runit(1,0)
fitH <-runit(0,1)
fitPH<-runit(1,1)

## pairs whose candidate parent sits in the plus group, and of those the ones
## babyM.R would have given a structural zero (back-calculated age < minAge)
minA<-conf0$minAge; maxA<-conf0$maxAge
cc<-cbind(ckmr$year1-ckmr$age1, ckmr$year2-ckmr$age2)
first<-cc[,1]<cc[,2]
pyv<-ifelse(first, ckmr$year1, ckmr$year2); pav<-ifelse(first, ckmr$age1, ckmr$age2)
jyv<-ifelse(first, ckmr$year2, ckmr$year1); jav<-ifelse(first, ckmr$age2, ckmr$age1)
bv<-jyv-jav; lagv<-pyv-bv
usedv<-seq_len(nrow(ckmr)) %in% fitPH$dat$ckmrPrep$popRow
plusv<-usedv & pav==maxA
zeroedByBabyM<-plusv & (pav-lagv) < minA

psiOf<-function(f) exp(f$pl$logPsim1[1])+1
ssbT <-function(f){ s<-ssbtable(f); n<-nrow(s); s[n,] }
wid  <-function(f){ s<-ssbT(f); (s[3]-s[2])/(function(){s0<-ssbT(fit0); s0[3]-s0[2]})() }

## how many pairs need the plus-group POP treatment under this design
plusRows<-length(fitPH$dat$ckmrPrep$popPlusRow)

## ---- figure, following the three panels at the end of ckmr/babyM.R ----
## baseline assessment in red, close-kin fit overlaid in blue, psi in the corner
kt<-function(x) exp(x)/1000
panel<-function(f, ttl){
  lo<-min(kt(fit0$sdrep$value[names(fit0$sdrep$value)=="logssb"] -
             2*fit0$sdrep$sd[names(fit0$sdrep$value)=="logssb"]),
          kt(f$sdrep$value[names(f$sdrep$value)=="logssb"] -
             2*f$sdrep$sd[names(f$sdrep$value)=="logssb"]))
  hi<-max(kt(fit0$sdrep$value[names(fit0$sdrep$value)=="logssb"] +
             2*fit0$sdrep$sd[names(fit0$sdrep$value)=="logssb"]),
          kt(f$sdrep$value[names(f$sdrep$value)=="logssb"] +
             2*f$sdrep$sd[names(f$sdrep$value)=="logssb"]))
  ## ssbplot fixes ylab and trans, so go through plotit for the 1000 t scale
  plotit(fit0, "logssb", trans=kt, ylab="SSB (1000 t)", ylim=c(min(0,lo), hi),
         col="darkred", cicol=rgb(.55,0,0,.35))
  plotit(f, "logssb", add=TRUE, trans=kt, col="darkblue", cicol=rgb(0,0,.55,.35))
  title(ttl)
  ## the two mean lines very nearly coincide - what separates the fits is the
  ## width of the band - so key the legend on the bands, not on the lines
  legend("topleft", bty="n", border=NA,
         legend=c("official assessment (95%)", "+ close kin (95%)",
                  paste0("psi = ", format(psiOf(f), digits=3))),
         fill=c(rgb(.55,0,0,.35), rgb(0,0,.55,.35), NA))
}
## cairo so the semi-transparent confidence bands actually blend
png("babyM.png", width=1500, height=520, pointsize=15,
    type=if(capabilities("cairo")) "cairo" else "Xlib")
opar<-par(mfrow=c(1,3), mar=c(4,4,3,1))
panel(fitP,  "official data + add pop")
panel(fitH,  "official data + add hsp")
panel(fitPH, "official data + add pop & hsp")
par(opar); dev.off()
plotok<-file.exists("babyM.png") && file.size("babyM.png") > 10000

writeLines(c(
  paste("base      ", format(fit0$opt$objective,  digits=7)),
  paste("nPOP      ", sum(ckmr$nPOP)),
  paste("nHSP      ", sum(ckmr$nHSP)),
  paste("npairs    ", nrow(ckmr)),
  paste("plusrows  ", plusRows),
  paste("plusPOP   ", sum(ckmr$nPOP[plusv])),
  paste("babyMnil  ", sum(zeroedByBabyM), sum(ckmr$nPOP[zeroedByBabyM])),
  paste("pop       ", format(fitP$opt$objective,  digits=7)),
  paste("hsp       ", format(fitH$opt$objective,  digits=7)),
  paste("pophsp    ", format(fitPH$opt$objective, digits=7)),
  paste("psi.pop   ", format(psiOf(fitP),  digits=5)),
  paste("psi.hsp   ", format(psiOf(fitH),  digits=5)),
  paste("psi.pophsp", format(psiOf(fitPH), digits=5)),
  paste("ssb.base  ", format(ssbT(fit0)[1], digits=7)),
  paste("ssb.pophsp", format(ssbT(fitPH)[1], digits=7)),
  paste("ciratio   ", format(wid(fitPH), digits=4)),
  paste("plot      ", plotok)
), "res.out")
