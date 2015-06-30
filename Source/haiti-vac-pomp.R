source("Source/R/leakyvac-pomp-model-inC-seas.R")
palette(brewer.pal(8,"Set1"))
set.seed(243947892)

## 2009 population estimtes but these could be pretty poor for first wave in particular
## due to massive population displacement

pop.portap <- 2.1e6
portap.dat <- get.haiti.data(first.wave.only=F)

## make bsplines
covartab <- make.covartab(-1,nrow(portap.dat)+1,byt=1,degree=6,nbasis=6)

## --------------------------------------------------- ##
## Do calculations to hardcode vaccination into models ##
## --------------------------------------------------- ##

pct.vac <- 0.5
doses <- pop.portap*pct.vac
total.vac.days <- 10
daily.vac.doses <- doses/total.vac.days
inter.dose.time <- 14 # in days
vactime <- 60 #  days after first case, December 20th 2010
one.of.one.start <- vactime
one.of.one.end <- one.of.one.start+total.vac.days
one.of.two.start <- vactime
one.of.two.end <- one.of.two.start+total.vac.days/2
two.of.two.start <- one.of.two.end+inter.dose.time
two.of.two.end <- two.of.two.start + total.vac.days/2

## ---------------------------- ##
## build con pomp model objects ##
## ---------------------------- ##

## note that we can't pass global params to our functions here (i.e.
## vaccination campaign details so we will have to hard code them (in R)
pap.mod.0dose <- build.leaky.model.C(pop=pop.portap,
                                        dat=portap.dat,
                                        my.times="day",
                                        covar=covartab,
                                        my.t0=-1,
                                        model.name="pap_zero_dose")

source("Source/R/leakyvac-pomp-model-inC-1dose-haiti.R")
pap.mod.1dose <- build.leaky.model.C(pop=pop.portap,
                                        dat=portap.dat,
                                        my.times="day",
                                        my.t0=-1,
                                        covar=covartab,
                                        model.name="pap_single_dose")

source("Source/R/leakyvac-pomp-model-inC-2dose-haiti.R")
pap.mod.2dose <- build.leaky.model.C(pop=pop.portap,
                                     dat=portap.dat,
                                     my.times="day",
                                     my.t0=-1,
                                     covar=covartab,
                                     model.name="pap_two_dose")

## ------------------------------------------------------------- ##
## load params from non-vaccination fit and add new ones for the ##
## vaccation model                                               ##
## ------------------------------------------------------------- ##

## load params from no-vac mif
mif.novac <- readRDS("GeneratedData/mif-haiti.rds")
mif.novac.params <- coef(mif.novac)

## start params for vac model
## need to include additional params and states
start.params.pap <- mif.novac.params # get rid of beta 1 since we fit to truncated series past the change point
start.params.pap <- c(start.params.pap,
                      theta1=0.44,
                      theta2=0.77,
                      phi1=0,
                      phi2=0,
                      kappa=0.9,
                      S1.0=0,E1.0=0,I1.0=0,A1.0=0,R1.0=0,
                      S2.0=0,E2.0=0,I2.0=0,A2.0=0,R2.0=0)


## now we will update the state of the system at the time of vaccination
t.vac <- 59 # vaccination start time
est.states <- readRDS("GeneratedData/mif-haiti-states.rds")
time.of.vac.state <- est.states[,t.vac][1:5] / sum(est.states[,t.vac][1:5]) # note this sums to 1 person less than the population since this is the average of people in states across particles at each time
start.params.pap[c('S.0','E.0','I.0','A.0','R.0')] <- time.of.vac.state

## create windowed models
pap.mod.0dose.win <- window(pap.mod.0dose,start=60,end=344)
pap.mod.0dose.win@t0 <- 59
pap.mod.1dose.win <- window(pap.mod.1dose,start=60,end=344)
pap.mod.1dose.win@t0 <- 59
pap.mod.2dose.win <- window(pap.mod.2dose,start=60,end=344)
pap.mod.2dose.win@t0 <- 59

## ---------------------------------- ##
## run simulations with the same seed ##
## ---------------------------------- ##
nsims <- 10000
pap.tm.sim <- simulate(pap.mod.0dose.win,
                       start.params.pap,
                       nsim=nsims,
                       seed=1914679109L,
                       transform=TRUE)

pap.tm.sim.1dose <- simulate(pap.mod.1dose.win,
                             start.params.pap,
                             nsim=nsims,
                             seed=1914679109L,
                             transform=TRUE)

pap.tm.sim.2dose <- simulate(pap.mod.2dose.win,
                             params=start.params.pap,
                             nsim=nsims,
                             seed=1914679109L,
                             transform=TRUE)

 ## ------------------------- ##
 ## Some stats for the papers ##
 ## ------------------------- ##

nodose.quants.pap <- apply(sapply(pap.tm.sim,function(x) x@data[1,]),1,function(y) quantile(y,c(.025,.5,.975)))
onedose.quants.pap <- apply(sapply(pap.tm.sim.1dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
twodose.quants.pap <- apply(sapply(pap.tm.sim.2dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
nodose.mean.pap <- apply(sapply(pap.tm.sim,function(x) x@data[1,]),1,mean)
onedose.mean.pap <- apply(sapply(pap.tm.sim.1dose,function(x) x@data[1,]),1,mean)
twodose.mean.pap <- apply(sapply(pap.tm.sim.2dose,function(x) x@data[1,]),1,mean)

## cases averted calculations
onedose.mat <- sapply(pap.tm.sim.1dose,function(x) x@data[1,])
twodose.mat <- sapply(pap.tm.sim.2dose,function(x) x@data[1,])
nodose.mat <- sapply(pap.tm.sim,function(x) x@data[1,])

## cum sum of cases for each simulation
nodose.sums <- apply(nodose.mat,2,sum)
onedose.sums <- apply(onedose.mat,2,sum)
twodose.sums <- apply(twodose.mat,2,sum)

cases.averted.1dose <- nodose.sums-onedose.sums
mean(cases.averted.1dose)
quantile(cases.averted.1dose,prob=c(.025,.975))

cases.averted.2dose <- nodose.sums-twodose.sums
mean(cases.averted.2dose)
ratio.of.cases.averted <- cases.averted.1dose/cases.averted.2dose

mean(ratio.of.cases.averted)
quantile(ratio.of.cases.averted,prob=c(.025,.975))

#pred.cfr <- c(rep(pred.cfr[1],15),pred.cfr)
pred.cfr <- 0.01
deaths.averted.1dose <- apply((nodose.mat-onedose.mat)*pred.cfr,2,sum)
deaths.averted.2dose <- apply((nodose.mat-twodose.mat)*pred.cfr,2,sum)
ratio.deaths.averted <- deaths.averted.1dose/deaths.averted.2dose
mean(ratio.of.cases.averted)
quantile(ratio.deaths.averted,prob=c(.025,.975))
mean(deaths.averted.1dose-deaths.averted.2dose)
mean(deaths.averted.1dose)
quantile(deaths.averted.1dose,c(0.025,0.975))

## ------------ ##
## The pap plot ##
## ------------ ##
vac.day <- 60
last.day <- 344
pdf("Plots/fig2-pap-50pct-REV.pdf",width=4.5,height=3)
#quartz("",width=4,height=3)
par(mfrow = c(1,1),mar=c(2.5,2.5,.5,.5),oma=c(.5,.5,.5,0.5),mgp=c(1.25,.3,0),tck=-.02)
plot(portap.dat[,2],pch=4,col=AddAlpha("black",.3),cex=.5,axes=FALSE,
     ylab="cases per day",xlab="epidemic day",cex.axis=.8,ylim=c(0,1500))
grid()
#for (i in 1:500) lines(start.week:50,pap.tm.sim[[i]]@data[1,],lty=2,col=AddAlpha("grey",.03))
axis(2)
axis(1)
lines(vac.day:last.day,onedose.mean.pap,col=3,lwd=2,lty=1)
lines(vac.day:last.day,twodose.mean.pap,col=2,lwd=2,lty=1)
lines(vac.day:last.day,nodose.mean.pap,col="darkgrey",lwd=2,lty=2)
legend(x=200,y=1500,
       legend=c("single-dose (mean) epidemic curve","two-dose (mean) epidemic curve","unvaccinated (mean) epidemic curve","data"),
       col=c(3,2,"darkgrey",AddAlpha("black",.5)),
       lwd=c(2,2,2,1),
       pch=c(-1,-1,-1,4),
       lty=c(1,1,2,0),bty="n",cex=.5)
abline(v=60,lty=2,col=AddAlpha("black",.5))
dev.off()

for (i in 1:500) lines(vac.day:last.day,pap.tm.sim[[i]]@data[1,],col=AddAlpha('grey',.05))


## ---------------------------------- ##
## plot with vaccination trajectories ##
## ---------------------------------- ##
## ------------ ##
## The pap plot ##
## ------------ ##
vac.day <- 60
last.day <- 344
pdf("Plots/fig2-pap-50pct-traj-REV.pdf",width=4.5,height=3)
#quartz("",width=4.5,height=3)
par(mfrow = c(1,1),mar=c(2.5,2.5,.5,.5),oma=c(.5,.5,.5,0.5),mgp=c(1.25,.3,0),tck=-.02)
plot(portap.dat[,2],pch=4,col=AddAlpha("black",.3),cex=.5,axes=FALSE,
     ylab="cases per day",xlab="epidemic day",cex.axis=.8,ylim=c(0,1500))
grid()
for (i in 1:150) lines(vac.day:last.day,pap.tm.sim[[i]]@data[1,],lty=2,col=AddAlpha("grey",.05))
for (i in 1:150) lines(vac.day:last.day,pap.tm.sim.1dose[[i]]@data[1,],lty=2,col=AddAlpha(3,.07))
for (i in 1:150) lines(vac.day:last.day,pap.tm.sim.2dose[[i]]@data[1,],lty=2,col=AddAlpha(2,.07))
axis(2)
axis(1)
lines(vac.day:last.day,onedose.mean.pap,col=3,lwd=2,lty=1)
lines(vac.day:last.day,twodose.mean.pap,col=2,lwd=2,lty=1)
lines(vac.day:last.day,nodose.mean.pap,col="darkgrey",lwd=2,lty=2)
legend(x=200,y=1500,
       legend=c("single-dose (mean) epidemic curve","two-dose (mean) epidemic curve","unvaccinated (mean) epidemic curve","data"),
       col=c(3,2,"darkgrey",AddAlpha("black",.5)),
       lwd=c(2,2,2,1),
       pch=c(-1,-1,-1,4),
       lty=c(1,1,2,0),bty="n",cex=.5)
abline(v=60,lty=2,col=AddAlpha("black",.5))
dev.off()

for (i in 1:500) lines(vac.day:last.day,pap.tm.sim[[i]]@data[1,],col=AddAlpha('grey',.05))


###########################
## Exploring other MRSEs ##
###########################
## first 50% greater MRSE
start.params.pap['theta1'] <- 0.44*1.5

pap.tm.sim <- simulate(pap.mod.0dose.win,
                       start.params.pap,
                       nsim=10000,
                       seed=1914679109L,
                       transform=TRUE)

pap.tm.sim.1dose <- simulate(pap.mod.1dose.win,
                             start.params.pap,
                             nsim=10000,
                             seed=1914679109L,
                             transform=TRUE)

pap.tm.sim.2dose <- simulate(pap.mod.2dose.win,
                             params=start.params.pap,
                             nsim=10000,
                             seed=1914679109L,
                             transform=TRUE)

nodose.quants.pap <- apply(sapply(pap.tm.sim,function(x) x@data[1,]),1,function(y) quantile(y,c(.025,.5,.975)))
onedose.quants.pap <- apply(sapply(pap.tm.sim.1dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
twodose.quants.pap <- apply(sapply(pap.tm.sim.2dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
nodose.mean.pap <- apply(sapply(pap.tm.sim,function(x) x@data[1,]),1,mean)
onedose.mean.pap <- apply(sapply(pap.tm.sim.1dose,function(x) x@data[1,]),1,mean)
twodose.mean.pap <- apply(sapply(pap.tm.sim.2dose,function(x) x@data[1,]),1,mean)

## cases averted calculations
onedose.mat <- sapply(pap.tm.sim.1dose,function(x) x@data[1,])
twodose.mat <- sapply(pap.tm.sim.2dose,function(x) x@data[1,])
nodose.mat <- sapply(pap.tm.sim,function(x) x@data[1,])

## cum sum of cases for each simulation
nodose.sums <- apply(nodose.mat,2,sum)
onedose.sums <- apply(onedose.mat,2,sum)
twodose.sums <- apply(twodose.mat,2,sum)

cases.averted.1dose <- nodose.sums-onedose.sums
mean(cases.averted.1dose)
quantile(cases.averted.1dose,prob=c(.025,.975))

cases.averted.2dose <- nodose.sums-twodose.sums
mean(cases.averted.2dose)
ratio.of.cases.averted <- cases.averted.1dose/cases.averted.2dose

mean(ratio.of.cases.averted)
quantile(ratio.of.cases.averted,prob=c(.025,.975))

## then 50% lower MRSE
start.params.pap['theta1'] <- 0.44*0.5

pap.tm.sim <- simulate(pap.mod.0dose.win,
                       start.params.pap,
                       nsim=100000,
                       seed=1914679109L,
                       transform=TRUE)

pap.tm.sim.1dose <- simulate(pap.mod.1dose.win,
                             start.params.pap,
                             nsim=100000,
                             seed=1914679109L,
                             transform=TRUE)

pap.tm.sim.2dose <- simulate(pap.mod.2dose.win,
                             params=start.params.pap,
                             nsim=100000,
                             seed=1914679109L,
                             transform=TRUE)

nodose.quants.pap <- apply(sapply(pap.tm.sim,function(x) x@data[1,]),1,function(y) quantile(y,c(.025,.5,.975)))
onedose.quants.pap <- apply(sapply(pap.tm.sim.1dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
twodose.quants.pap <- apply(sapply(pap.tm.sim.2dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
nodose.mean.pap <- apply(sapply(pap.tm.sim,function(x) x@data[1,]),1,mean)
onedose.mean.pap <- apply(sapply(pap.tm.sim.1dose,function(x) x@data[1,]),1,mean)
twodose.mean.pap <- apply(sapply(pap.tm.sim.2dose,function(x) x@data[1,]),1,mean)

## cases averted calculations
onedose.mat <- sapply(pap.tm.sim.1dose,function(x) x@data[1,])
twodose.mat <- sapply(pap.tm.sim.2dose,function(x) x@data[1,])
nodose.mat <- sapply(pap.tm.sim,function(x) x@data[1,])

## cum sum of cases for each simulation
nodose.sums <- apply(nodose.mat,2,sum)
onedose.sums <- apply(onedose.mat,2,sum)
twodose.sums <- apply(twodose.mat,2,sum)

cases.averted.1dose <- nodose.sums-onedose.sums
mean(cases.averted.1dose)
quantile(cases.averted.1dose,prob=c(.025,.975))

cases.averted.2dose <- nodose.sums-twodose.sums
mean(cases.averted.2dose)
ratio.of.cases.averted <- cases.averted.1dose/cases.averted.2dose

mean(ratio.of.cases.averted)
quantile(ratio.of.cases.averted,prob=c(.025,.975))



### EXPLOREATION OF ALTERNATIVE MODELS

## exploring alternative fits for Haiti
haiti.para <- readRDS(file="GeneratedData/parallel-mif-portaup-6df-fitrho.rds")
logliks <- colMeans(sapply(haiti.para,function(x) x[[2]]))
ordered.liks <- order(logliks,decreasing = T)

###################
## the first one ##
###################

## load params from no-vac mif
ratio.of.cases.averted.list <- list()

for (i in 1:5){
    print(i)
    mif.novac <- haiti.para[ordered.liks[i]][[1]][[1]]
    mif.novac.params <- coef(mif.novac)

    ## start params for vac model
    ## need to include additional params and states
    start.params.pap <- mif.novac.params # get rid of beta 1 since we fit to truncated series past the change point
    start.params.pap <- c(start.params.pap,
                          theta1=0.44,
                          theta2=0.77,
                          phi1=0,
                          phi2=0,
                          kappa=0.9,
                          S1.0=0,E1.0=0,I1.0=0,A1.0=0,R1.0=0,
                          S2.0=0,E2.0=0,I2.0=0,A2.0=0,R2.0=0)

    ## now we will update the state of the system at the time of vaccination
    t.vac <- 59 # vaccination start time

    pf.best <- pfilter(mif.novac,Np=20000,save.states=TRUE)
    est.states <- sapply(pf.best@saved.states,rowMeans)

    time.of.vac.state <- est.states[,t.vac][1:5] / sum(est.states[,t.vac][1:5]) # note this sums to 1 person less than the population since this is the average of people in states across particles at each time
    start.params.pap[c('S.0','E.0','I.0','A.0','R.0')] <- time.of.vac.state

    ## create windowed models
    pap.mod.0dose.win <- window(pap.mod.0dose,start=60,end=344)
    pap.mod.0dose.win@t0 <- 59
    pap.mod.1dose.win <- window(pap.mod.1dose,start=60,end=344)
    pap.mod.1dose.win@t0 <- 59
    pap.mod.2dose.win <- window(pap.mod.2dose,start=60,end=344)
    pap.mod.2dose.win@t0 <- 59

    ## ---------------------------------- ##
    ## run simulations with the same seed ##
    ## ---------------------------------- ##
    nsims <- 25000
    pap.tm.sim <- simulate(pap.mod.0dose.win,
                           start.params.pap,
                           nsim=nsims,
                           seed=1914679109L,
                           transform=TRUE)

    pap.tm.sim.1dose <- simulate(pap.mod.1dose.win,
                                 start.params.pap,
                                 nsim=nsims,
                                 seed=1914679109L,
                                 transform=TRUE)

    pap.tm.sim.2dose <- simulate(pap.mod.2dose.win,
                                 params=start.params.pap,
                                 nsim=nsims,
                                 seed=1914679109L,
                                 transform=TRUE)

    ## ------------------------- ##
    ## Some stats for the papers ##
    ## ------------------------- ##

    nodose.quants.pap <- apply(sapply(pap.tm.sim,function(x) x@data[1,]),1,function(y) quantile(y,c(.025,.5,.975)))
    onedose.quants.pap <- apply(sapply(pap.tm.sim.1dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
    twodose.quants.pap <- apply(sapply(pap.tm.sim.2dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
    nodose.mean.pap <- apply(sapply(pap.tm.sim,function(x) x@data[1,]),1,mean)
    onedose.mean.pap <- apply(sapply(pap.tm.sim.1dose,function(x) x@data[1,]),1,mean)
    twodose.mean.pap <- apply(sapply(pap.tm.sim.2dose,function(x) x@data[1,]),1,mean)

    ## cases averted calculations
    onedose.mat <- sapply(pap.tm.sim.1dose,function(x) x@data[1,])
    twodose.mat <- sapply(pap.tm.sim.2dose,function(x) x@data[1,])
    nodose.mat <- sapply(pap.tm.sim,function(x) x@data[1,])

    ## cum sum of cases for each simulation
    nodose.sums <- apply(nodose.mat,2,sum)
    onedose.sums <- apply(onedose.mat,2,sum)
    twodose.sums <- apply(twodose.mat,2,sum)

    cases.averted.1dose <- nodose.sums-onedose.sums
    mean(cases.averted.1dose)
    quantile(cases.averted.1dose,prob=c(.025,.975))

    cases.averted.2dose <- nodose.sums-twodose.sums
    mean(cases.averted.2dose)
    ratio.of.cases.averted.list[[i]] <- cases.averted.1dose/cases.averted.2dose

    mean(ratio.of.cases.averted)
    quantile(ratio.of.cases.averted,prob=c(.025,.975))
}

## make a table for the supplment
library(xtable)
xtable(
    cbind(t(round(sapply(ratio.of.cases.averted.list,function(x) c(mean(x),quantile(x,prob=c(.025,.975)))),2)),
          t(sapply(1:5,function(x) coef(haiti.para[ordered.liks[x]][[1]][[1]]))))
   ,digits=4)

## let's look at the coefficients for each
sapply(1:5,function(x) coef(haiti.para[ordered.liks[x]][[1]][[1]]))

## let's calculate R0 for each
sapply(1:5,function(x) {
    coefs <- coef(haiti.para[ordered.liks[x]][[1]][[1]])
    sum(covartab[2,-1]*coefs[grep('beta',names(coefs))])/.5
})
                        )


