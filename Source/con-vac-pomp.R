###########################################################
## functions for running Conakry vaccination simulations ##
## for single dose paper                                 ##
###########################################################

source("Source/R/leakyvac-pomp-model-inC.R")
palette(brewer.pal(8,"Set1"))
set.seed(243947892)

## 2010 populatino estimates from Institut National de la Statistique de Guinée
pop.con <- 1656300
## condat
con.dat <- get.conakry.data()

## --------------------------------------------------- ##
## Do calculations to hardcode vaccination into models ##
## --------------------------------------------------- ##

pct.vac <- 0.5
doses <- pop.con*pct.vac
total.vac.days <- 10
daily.vac.doses <- doses/total.vac.days
inter.dose.time <- 14 # in days
vactime <- 60 #  days after first case
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
con.mod.0dose <- build.leaky.model.C(pop=pop.con,
                                     dat=con.dat,
                                     my.times="day",
                                     my.t0=0,
                                     model.name="con_zero_dose")
con.mod.0dose <- window(con.mod.0dose,start=33,end=153)
con.mod.0dose@t0 <- 32

source("Source/R/leakyvac-pomp-model-inC-1dose-con.R")
con.mod.1dose <- build.leaky.model.C(pop=pop.con,
                                     dat=con.dat,
                                     my.times="day",
                                     my.t0=0,
                                     model.name="con_single_dose")
con.mod.1dose <- window(con.mod.1dose,start=33,end=153)
con.mod.1dose@t0 <- 32

source("Source/R/leakyvac-pomp-model-inC-2dose-con.R")
con.mod.2dose <- build.leaky.model.C(pop=pop.con,
                                     dat=con.dat,
                                     my.times="day",
                                     my.t0=0,
                                     model.name="con_two_dose")
con.mod.2dose <- window(con.mod.2dose,start=33,end=153)
con.mod.2dose@t0 <- 32

## ------------------------------------------------------------- ##
## load params from non-vaccination fit and add new ones for the ##
## vaccation model                                               ##
## ------------------------------------------------------------- ##

## load params from no-vac mif
mif.novac <- readRDS("GeneratedData/mif-con.rds")
mif.novac.params <- coef(mif.novac)


## start params for vac model
## need to include additional params and states
start.params.con <- mif.novac.params # get rid of beta 1 since we fit to truncated series past the change point
start.params.con <-start.params.con[-which(names(start.params.con)=='beta2')]
start.params.con <- c(start.params.con,
                      theta1=0.44,
                      theta2=0.77,#0.73,
                      phi1=0,
                      phi2=0,
                      kappa=0.9,
                      S1.0=0,E1.0=0,I1.0=0,A1.0=0,R1.0=0,
                      S2.0=0,E2.0=0,I2.0=0,A2.0=0,R2.0=0)


## now we will update the state of the system at the time of vaccination
t.vac <- 27 # vaccination start time (this is an addition to the t0=32 from the fitting thus vaccination starts at 60)
est.states <- readRDS("GeneratedData/mif-con-states.rds")
time.of.vac.state <- est.states[,t.vac][1:5] / sum(est.states[,t.vac][1:5]) # note this sums to 1 person less than the population since this is the average of people in states across particles at each time
start.params.con[c('S.0','E.0','I.0','A.0','R.0')] <- time.of.vac.state

## create windowed models
con.mod.0dose.win <- window(con.mod.0dose,start=60,end=153)
con.mod.0dose.win@t0 <- 59
con.mod.1dose.win <- window(con.mod.1dose,start=60,end=153)
con.mod.1dose.win@t0 <- 59
con.mod.2dose.win <- window(con.mod.2dose,start=60,end=153)
con.mod.2dose.win@t0 <- 59

## ---------------------------------- ##
## run simulations with the same seed ##
## ---------------------------------- ##
nsims <- 10000#0
con.tm.sim <- simulate(con.mod.0dose.win,
                       start.params.con,
                       nsim=nsims,
                       seed=1914679109L,
                       transform=TRUE)

con.tm.sim.1dose <- simulate(con.mod.1dose.win,
                             start.params.con,
                             nsim=nsims,
                             seed=1914679109L,
                             transform=TRUE)

con.tm.sim.2dose <- simulate(con.mod.2dose.win,
                             params=start.params.con,
                             nsim=nsims,
                             seed=1914679109L,
                             transform=TRUE)


 ## ------------------------- ##
 ## Some stats for the papers ##
 ## ------------------------- ##

nodose.quants.con <- apply(sapply(con.tm.sim,function(x) x@data[1,]),1,function(y) quantile(y,c(.025,.5,.975)))
onedose.quants.con <- apply(sapply(con.tm.sim.1dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
twodose.quants.con <- apply(sapply(con.tm.sim.2dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
nodose.mean.con <- apply(sapply(con.tm.sim,function(x) x@data[1,]),1,mean)
onedose.mean.con <- apply(sapply(con.tm.sim.1dose,function(x) x@data[1,]),1,mean)
twodose.mean.con <- apply(sapply(con.tm.sim.2dose,function(x) x@data[1,]),1,mean)

## cases averted calculations
onedose.mat <- sapply(con.tm.sim.1dose,function(x) x@data[1,])
twodose.mat <- sapply(con.tm.sim.2dose,function(x) x@data[1,])
nodose.mat <- sapply(con.tm.sim,function(x) x@data[1,])

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
pred.cfr <- 0.018
deaths.averted.1dose <- apply((nodose.mat-onedose.mat)*pred.cfr,2,sum)
deaths.averted.2dose <- apply((nodose.mat-twodose.mat)*pred.cfr,2,sum)
ratio.deaths.averted <- deaths.averted.1dose/deaths.averted.2dose
mean(ratio.of.cases.averted)
quantile(ratio.deaths.averted,prob=c(.025,.975))
mean(deaths.averted.1dose-deaths.averted.2dose)
mean(deaths.averted.1dose)
quantile(deaths.averted.1dose,c(0.025,0.975))

## ------------ ##
## The con plot ##
## ------------ ##
vac.day <- 60
last.day <- 153
pdf("Plots/fig2-con-50pct.pdf",width=4.5,height=3)
#quartz("",width=4,height=3)
par(mfrow = c(1,1),mar=c(2.5,2.5,.5,.5),oma=c(.5,.5,.5,0.5),mgp=c(1.25,.3,0),tck=-.02)
plot(con.dat[,2],pch=4,col=AddAlpha("black",.3),cex=.5,axes=FALSE,
     ylab="cases per day",xlab="epidemic day",cex.axis=.8,ylim=c(0,200))
grid()
#for (i in 1:500) lines(start.week:50,con.tm.sim[[i]]@data[1,],lty=2,col=AddAlpha("grey",.03))
axis(2)
axis(1)
lines(vac.day:last.day,onedose.mean.con,col=3,lwd=2,lty=1)
lines(vac.day:last.day,twodose.mean.con,col=2,lwd=2,lty=1)
lines(vac.day:last.day,nodose.mean.con,col="darkgrey",lwd=2,lty=2)
legend(x=2,y=200,
       legend=c("single-dose (mean) epidemic curve","two-dose (mean) epidemic curve","unvaccinated (mean) epidemic curve","data"),
       col=c(3,2,"darkgrey",AddAlpha("black",.5)),
       lwd=c(2,2,2,1),
       pch=c(-1,-1,-1,4),
       lty=c(1,1,2,0),bty="n",cex=.5)
abline(v=60,lty=2,col=AddAlpha("black",.5))
dev.off()


######################
## The con vac plot ##
######################

vac.day <- 60
last.day <- 153
pdf("Plots/fig2-con-50pct-traj.pdf",width=4.5,height=3)
#quartz("",width=4,height=3)
par(mfrow = c(1,1),mar=c(2.5,2.5,.5,.5),oma=c(.5,.5,.5,0.5),mgp=c(1.25,.3,0),tck=-.02)
plot(con.dat[,2],pch=4,col=AddAlpha("black",.3),cex=.5,axes=FALSE,
     ylab="cases per day",xlab="epidemic day",cex.axis=.8,ylim=c(0,200))
grid()
for (i in 1:150) lines(vac.day:last.day,con.tm.sim[[i]]@data[1,],lty=2,col=AddAlpha("grey",.07))
for (i in 1:150) lines(vac.day:last.day,con.tm.sim.1dose[[i]]@data[1,],lty=2,col=AddAlpha(3,.07))
for (i in 1:150) lines(vac.day:last.day,con.tm.sim.2dose[[i]]@data[1,],lty=2,col=AddAlpha(2,.07))
axis(2)
axis(1)
lines(vac.day:last.day,onedose.mean.con,col=3,lwd=2,lty=1)
lines(vac.day:last.day,twodose.mean.con,col=2,lwd=2,lty=1)
lines(vac.day:last.day,nodose.mean.con,col="darkgrey",lwd=2,lty=2)
legend(x=2,y=200,
       legend=c("single-dose (mean) epidemic curve","two-dose (mean) epidemic curve","unvaccinated (mean) epidemic curve","data"),
       col=c(3,2,"darkgrey",AddAlpha("black",.5)),
       lwd=c(2,2,2,1),
       pch=c(-1,-1,-1,4),
       lty=c(1,1,2,0),bty="n",cex=.5)
abline(v=60,lty=2,col=AddAlpha("black",.5))
dev.off()



#########################################
## Exploration of other values of MRSE ##
#########################################

## first for MRSE 50% higher
start.params.con['theta1'] <- 1.5*.44

con.tm.sim <- simulate(con.mod.0dose.win,
                       start.params.con,
                       nsim=10000,
                       seed=1914679109L,
                       transform=TRUE)

con.tm.sim.1dose <- simulate(con.mod.1dose.win,
                             start.params.con,
                             nsim=10000,
                             seed=1914679109L,
                             transform=TRUE)

con.tm.sim.2dose <- simulate(con.mod.2dose.win,
                             params=start.params.con,
                             nsim=10000,
                             seed=1914679109L,
                             transform=TRUE)


nodose.quants.con <- apply(sapply(con.tm.sim,function(x) x@data[1,]),1,function(y) quantile(y,c(.025,.5,.975)))
onedose.quants.con <- apply(sapply(con.tm.sim.1dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
twodose.quants.con <- apply(sapply(con.tm.sim.2dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
nodose.mean.con <- apply(sapply(con.tm.sim,function(x) x@data[1,]),1,mean)
onedose.mean.con <- apply(sapply(con.tm.sim.1dose,function(x) x@data[1,]),1,mean)
twodose.mean.con <- apply(sapply(con.tm.sim.2dose,function(x) x@data[1,]),1,mean)

## cases averted calculations
onedose.mat <- sapply(con.tm.sim.1dose,function(x) x@data[1,])
twodose.mat <- sapply(con.tm.sim.2dose,function(x) x@data[1,])
nodose.mat <- sapply(con.tm.sim,function(x) x@data[1,])

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


## now for it 50% lower
start.params.con['theta1'] <- .5*.44

con.tm.sim <- simulate(con.mod.0dose.win,
                       start.params.con,
                       nsim=10000,
                       seed=1914679109L,
                       transform=TRUE)

con.tm.sim.1dose <- simulate(con.mod.1dose.win,
                             start.params.con,
                             nsim=10000,
                             seed=1914679109L,
                             transform=TRUE)

con.tm.sim.2dose <- simulate(con.mod.2dose.win,
                             params=start.params.con,
                             nsim=10000,
                             seed=1914679109L,
                             transform=TRUE)

nodose.quants.con <- apply(sapply(con.tm.sim,function(x) x@data[1,]),1,function(y) quantile(y,c(.025,.5,.975)))
onedose.quants.con <- apply(sapply(con.tm.sim.1dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
twodose.quants.con <- apply(sapply(con.tm.sim.2dose,function(x) x@data[1,]),1,function(y) quantile(y,c(.25,.5,.75)))
nodose.mean.con <- apply(sapply(con.tm.sim,function(x) x@data[1,]),1,mean)
onedose.mean.con <- apply(sapply(con.tm.sim.1dose,function(x) x@data[1,]),1,mean)
twodose.mean.con <- apply(sapply(con.tm.sim.2dose,function(x) x@data[1,]),1,mean)

## cases averted calculations
onedose.mat <- sapply(con.tm.sim.1dose,function(x) x@data[1,])
twodose.mat <- sapply(con.tm.sim.2dose,function(x) x@data[1,])
nodose.mat <- sapply(con.tm.sim,function(x) x@data[1,])

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
