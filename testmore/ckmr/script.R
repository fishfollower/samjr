library(samjr)
cn<-read.ices("cn.dat")
cw<-read.ices("cw.dat")
dw<-read.ices("dw.dat")
lf<-read.ices("lf.dat")
lw<-read.ices("lw.dat")
mo<-read.ices("mo.dat")
nm<-read.ices("nm.dat")
pf<-read.ices("pf.dat")
pm<-read.ices("pm.dat")
sw<-read.ices("sw.dat")
surveys<-read.ices("survey.dat")

mkconf<-function(dat){
  conf<-defcon(dat)
  conf$keyLogFsta[1,] <- c(0, 1, 2, 3, 4, 5)
  conf$corFlag <- 2
  conf$keyLogFpar <- matrix(
           c(-1, -1, -1, -1, -1, -1,
              0,  1,  2,  3,  4, -1,
              5,  6,  7,  8, -1, -1
             ), nrow=3, byrow=TRUE)
  conf$keyVarF[1,] <- c(0, 1, 1, 1, 1, 1)
  conf$keyVarObs <- matrix(
           c( 0,  1,  2,  2,  2,  2,
              3,  4,  4,  4,  4, -1,
              5,  6,  6,  6, -1, -1
             ), nrow=3, byrow=TRUE)
  conf$noScaledYears <- 13
  conf$keyScaledYears <- 1993:2005
  conf$keyParScaledYA <- row(matrix(NA, nrow=13, ncol=6))-1
  conf$fbarRange <- c(2,4)
  conf$ckmrScale <- 1000   # nscod numbers-at-age are in thousands
  conf$ckmrPsi <- 1.5
  conf
}

base<-setup.sam.data(surveys=surveys, residual.fleet=cn, prop.mature=mo,
                     stock.mean.weight=sw, catch.mean.weight=cw,
                     dis.mean.weight=dw, land.mean.weight=lw,
                     prop.f=pf, prop.m=pm, natural.mortality=nm, land.frac=lf)
conf0<-mkconf(base)
fit0<-sam.fit(base, conf0, defpar(base, conf0), silent=TRUE)

## simulate a close-kin data set from the fitted model
set.seed(20240101)
ckmr<-simulateCKMR(fit0, years=2010:2014, n=rep(20000,5))

dat<-setup.sam.data(surveys=surveys, residual.fleet=cn, prop.mature=mo,
                    stock.mean.weight=sw, catch.mean.weight=cw,
                    dis.mean.weight=dw, land.mean.weight=lw,
                    prop.f=pf, prop.m=pm, natural.mortality=nm, land.frac=lf,
                    ckmr=ckmr)

runit<-function(usePOP, useHSP){
  conf<-mkconf(dat); conf$usePOP<-usePOP; conf$useHSP<-useHSP
  suppressWarnings(sam.fit(dat, conf, defpar(dat, conf), silent=TRUE))
}
fitP <-runit(1,0)
fitH <-runit(0,1)
fitPH<-runit(1,1)

out<-file("res.out")
writeLines(c(
  paste("nPOP      ", sum(ckmr$nPOP)),
  paste("nHSP      ", sum(ckmr$nHSP)),
  paste("npairs    ", nrow(ckmr)),
  paste("base      ", format(fit0$opt$objective,  digits=7)),
  paste("pop       ", format(fitP$opt$objective,  digits=7)),
  paste("hsp       ", format(fitH$opt$objective,  digits=7)),
  paste("pophsp    ", format(fitPH$opt$objective, digits=7)),
  ## CKMR must not add estimated parameters
  paste("npar      ", length(fit0$opt$par), length(fitPH$opt$par)),
  ## expected kin counts must reproduce the observed totals
  paste("expPOP    ", format(sum(ckmrtable(fitPH)$ePOP[ckmrtable(fitPH)$usePOP]), digits=6)),
  paste("expHSP    ", format(sum(ckmrtable(fitPH)$eHSP[ckmrtable(fitPH)$useHSP]), digits=6)),
  ## close-kin data must sharpen terminal SSB
  paste("ciratio   ", format((function(a,b){
     s0<-ssbtable(fit0); s1<-ssbtable(fitPH); n<-nrow(s0)
     (s1[n,3]-s1[n,2])/(s0[n,3]-s0[n,2])})(), digits=4))
), out)
close(out)
