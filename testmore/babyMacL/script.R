## The length-based close-kin analogue of babyM, on the same mackerel
## assessment (model_XplConf.08): official data plus POPL, plus HSPL, plus
## both, exactly the three runs babyM does for the age-based terms.
##
## Growth and fecundity are the package defaults, which are the values
## ckmr/babyL.R was built around: von Bertalanffy Linf 41.7 and k 0.35, a
## log-normal spread of 0.1 about it, reference length 50, and psi = 4.5 as a
## length exponent (the equivalent of 1.5 on weight, since W ~ l^3).
## The sampling years and numbers follow babyL.R as well; the length classes
## are far coarser, because the likelihood costs pairs x ages^2 and babyL.R's
## 50 classes would put this case in the hours.
library(samjr)

cn<-read.ices("cn.dat"); cw<-read.ices("cw.dat"); dw<-read.ices("dw.dat")
lf<-read.ices("lf.dat"); lw<-read.ices("lw.dat"); mo<-read.ices("mo.dat")
nm<-read.ices("nm.dat"); pf<-read.ices("pf.dat"); pm<-read.ices("pm.dat")
sw<-read.ices("sw.dat"); surveys<-read.ices("survey.dat")
recap<-read.table("tag.dat", header=TRUE)

W<-matrix(NA, nrow=nrow(cn), ncol=ncol(cn))
W[as.numeric(rownames(cn))<2000]<-10
attr(cn,"weight")<-W

mkdat<-function(ckmrl=NULL)
  setup.sam.data(surveys=surveys, residual.fleet=cn, prop.mature=mo,
                 stock.mean.weight=sw, catch.mean.weight=cw,
                 dis.mean.weight=dw, land.mean.weight=lw,
                 prop.f=pf, prop.m=pm, natural.mortality=nm,
                 land.frac=lf, recapture=recap, ckmrl=ckmrl)

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
  ## numbers-at-age are in thousands, as in babyL.R's scale=1000; the growth
  ## and fecundity settings are left at their defaults on purpose
  conf$ckmrlScale <- 1000
  conf$ckmrlEstimatePsi <- 1
  conf
}

dat<-mkdat()
conf0<-mkconf(dat)
fit0<-sam.fit(dat, conf0, defpar(dat, conf0), silent=TRUE)

## The sample size gives kin counts on the same scale as babyM's 383 POP and
## 331 HSP. Note it is the length-based conf$ckmrlScale that matters here, not
## the age-based conf$ckmrScale: at the wrong one the study finds nothing.
NBIN<-6
NPER<-1000000
set.seed(20240101)
ckmrl<-simulateCKMRL(fit0, years=2021:2023, n=rep(NPER,3), nBin=NBIN)
dat2<-mkdat(ckmrl)

## Seed each close-kin run from the baseline assessment rather than from the
## all-zero default, as babyM does: without a natural-mortality multiplier to
## absorb it the close-kin likelihood is bimodal in psi, and a cold start can
## land in the wrong basin.
runit<-function(usePOPL, useHSPL){
  conf<-mkconf(dat2)
  conf$usePOPL<-usePOPL; conf$useHSPL<-useHSPL
  par<-defpar(dat2, conf)
  for(nm in names(par))
    if(nm %in% names(fit0$pl) && length(par[[nm]])==length(fit0$pl[[nm]]))
      par[[nm]][]<-fit0$pl[[nm]]
  suppressWarnings(sam.fit(dat2, conf, par, silent=TRUE))
}

fitP <-runit(1,0)
fitH <-runit(0,1)
fitPH<-runit(1,1)

psiOf<-function(f) exp(f$pl$logPsiL[1])
ssbT <-function(f){ s<-ssbtable(f); n<-nrow(s); s[n,] }
wid  <-function(f){ s<-ssbT(f); (s[3]-s[2])/(function(){s0<-ssbT(fit0); s0[3]-s0[2]})() }
expOf<-function(f, what){ tt<-ckmrltable(f); sum(tt[[what]]) }

## ---- figure ----
kt<-function(x) exp(x)/1000
panel<-function(f, ttl){
  g<-function(fit, s) kt(fit$sdrep$value[names(fit$sdrep$value)=="logssb"] +
                         s*2*fit$sdrep$sd[names(fit$sdrep$value)=="logssb"])
  lo<-min(g(fit0,-1), g(f,-1)); hi<-max(g(fit0,1), g(f,1))
  plotit(fit0, "logssb", trans=kt, ylab="SSB (1000 t)", ylim=c(min(0,lo), hi),
         col="darkred", cicol=rgb(.55,0,0,.35))
  plotit(f, "logssb", add=TRUE, trans=kt, col="darkblue", cicol=rgb(0,0,.55,.35))
  title(ttl)
  legend("topleft", bty="n", border=NA,
         legend=c("official assessment (95%)", "+ close kin by length (95%)",
                  paste0("psi = ", format(psiOf(f), digits=3))),
         fill=c(rgb(.55,0,0,.35), rgb(0,0,.55,.35), NA))
}
png("babyMacL.png", width=1500, height=520, pointsize=15,
    type=if(capabilities("cairo")) "cairo" else "Xlib")
opar<-par(mfrow=c(1,3), mar=c(4,4,3,1))
panel(fitP,  "official data + add popl")
panel(fitH,  "official data + add hspl")
panel(fitPH, "official data + add popl & hspl")
par(opar); dev.off()
plotok<-file.exists("babyMacL.png") && file.size("babyMacL.png") > 10000

fmt<-function(x, d=7) format(x, digits=d)
writeLines(c(
  paste("base      ", fmt(fit0$opt$objective)),
  paste("npairs    ", nrow(ckmrl)),
  paste("nPOP      ", sum(ckmrl$nPOP)),
  paste("nHSP      ", sum(ckmrl$nHSP)),
  paste("popl      ", fmt(fitP$opt$objective)),
  paste("hspl      ", fmt(fitH$opt$objective)),
  paste("poplhspl  ", fmt(fitPH$opt$objective)),
  paste("psi.popl  ", fmt(psiOf(fitP), 5)),
  paste("psi.hspl  ", fmt(psiOf(fitH), 5)),
  paste("psi.both  ", fmt(psiOf(fitPH), 5)),
  paste("ePOP      ", fmt(expOf(fitPH, "ePOP"), 5)),
  paste("eHSP      ", fmt(expOf(fitPH, "eHSP"), 5)),
  paste("ciratio   ", fmt(wid(fitPH), 4)),
  paste("plot      ", plotok)
), "res.out")
