######################################################################################
## This file contains scripts/functions for producing key papers for the manuscript ##
######################################################################################
source("Source/base-functions.R")
palette(brewer.pal(8,"Dark2"))


## ------------------------------------------------- ##
## Variants of Epidemic Curve with min VE (Figure 1) ##
## ------------------------------------------------- ##

##' 88% VE,
##' get the minimums from runs of min-ves.R, min-leaky.R, and min-aon.R
load("GeneratedData/minVES_betaSAB_VE20p88-1kto500k.rda") #min.ves
load("GeneratedData/minVEleaky_betaSAB_VE20p88-1kto500k.rda") # min.leaky
load("GeneratedData/minVEaon_betaSAB_VE20p88-1kto500k.rda") #min.aon

## a bit of overkill for this plot to import all this but...
load("GeneratedData/vesruns_VE20p88_betaSAB.rda")
uncon.run <- my.biglist[[7]]

make.epicurve.mech.plot(uncon.run=uncon.run,
                        plot.name="epicurvevemech-88",
                        min.ves=min.ves,
                        min.leaky=min.leaky,
                        min.aon=min.aon,
                        vac.starts = seq(0,300,by=1),
                        mrse=TRUE,
                        ve2=0.88,
                        dose.cols = c(1,3),q
                        ve.ylim=c(0.05,.6))


## -------------------------------------------------- ##
## Population vs. Individual-level effect (figure 2?) ##
## -------------------------------------------------- ##

## making heat maps for logIRR for VES runs
## 88% VE, beta for SAB
# load("GeneratedData/vesruns_VE20p88_betaSAB-with500k.rda")
load("GeneratedData/vesruns_VE20p88_betaSAB-12feb15.rda")
tmp <- shape.mlapply.out(my.biglist,
                         vacdoses=c(100000,200000,500000),
                         ve1s=seq(0,0.88,by=0.02),
                         alt.structure=TRUE
                         )
fs.mat.1dose <- tmp$fs.mat.1dose
fs.mat.2dose <- tmp$fs.mat.2dose


make.counterfact.plot2(fs.mat.1dose,
                       fs.mat.2dose,
                       tmp$time.at.risk.1dose,
                       tmp$time.at.risk.2dose,
                       ve1s=seq(0,0.88,by=.02),
                       vac.timings=seq(0,300,by=2),
                       vacdoses=c(100000,200000,500000),
                       add.unvac=F,
                       dosenums = c(1,2,3),
                       mrse=TRUE,
                       pdf=F,
                       irr.all=FALSE,
                       plot.name="vesheatmap_88-500k-REV.pdf")

## 88% VE, SAB 2008 beta
## leaky vaccine (VE_SP)
#load("GeneratedData/leakyruns_VE20p88_betaSAB-with500k.rda")
load("GeneratedData/leakyruns_VE20p88_betaSAB-REV.rda")
#load("GeneratedData/leakyruns_VE20p88_betaSAB.rda")
tmp <- shape.mlapply.out(my.biglist,
                         vacdoses=c(100000,200000,500000),
                         ve1s=seq(0,0.88,by=0.02),
                         alt.structure = TRUE
                         )
fs.mat.1dose <- tmp$fs.mat.1dose
fs.mat.2dose <- tmp$fs.mat.2dose

make.counterfact.plot2(fs.mat.1dose,
                       fs.mat.2dose,
                       tmp$time.at.risk.1dose,tmp$time.at.risk.2dose,
                       ve1s=seq(0,0.88,by=.02),
                       vac.timings=seq(0,300,by=2),
                       vacdoses=c(100000,200000,500000),
                       add.unvac=F,
                       dosenums = c(1:3),
                       mrse=TRUE,
                       pdf=TRUE,plot.name="leaky_heatmap_88-500k-REV.pdf")

#load("GeneratedData/aonruns_VE20p88_betaSAB-with500k.rda")
load("GeneratedData/aonruns_VE20p88_betaSAB-REV.rda")
#load("GeneratedData/aonruns_VE20p88_betaSAB.rda")
tmp <- shape.mlapply.out(my.biglist,
                         vacdoses=c(100000,200000,500000),
                         ve1s=seq(0,0.88,by=0.02),
                         alt.structure=FALSE)

fs.mat.1dose <- tmp$fs.mat.1dose
fs.mat.2dose <- tmp$fs.mat.2dose

make.counterfact.plot2(fs.mat.1dose,fs.mat.2dose,
                       tmp$time.at.risk.1dose,tmp$time.at.risk.2dose,
                       ve1s=seq(0,0.88,by=.02),
                       vac.timings=seq(0,300,by=2),
                       vacdoses=c(100000,200000,500000),
                       add.unvac=F,
                       dosenums = c(1:3),
                       mrse=TRUE,
                       pdf=F,plot.name="aon_heatmap_88-500k-REV.pdf")


### cumulative incidence plots

## inputs
beta.sab <- 0.6538415

my.params.leaky <- list(beta=beta.sab,
                     gamma=1/2,
                     sigma=1/1.4,
                     phi1=0.0, # VE_I for first dose among fully symptomatic
                     phi2=0.0, # VE_I for second dose among fully symptomtic
                     kappa1=0.9, # VE_I for first dose among asym/mildly symptomatic
                     kappa2=0.9, # VE_I for second dose among asym/mildly symptomtic
                     theta0=0,
                     theta1=0.44,
                     theta2=0.80
                     )


## simulate epidemics with specfiic number of doses and timing
vac.timings <- c(4*7,14*7,19*7)
initial.state <- c(.5e6-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
doses <- c(500000*.2,500000*.3)
runs <- expand.grid("time"=vac.timings,"dose"=doses)
times <- seq(0,400,by=1)

uncon.run <- run.single.vac <- run.double.vac <- vector("list",length=nrow(runs))

for (r in 1:nrow(runs)){
    uncon.run[[r]] <- ode(y=initial.state,
                          times=times,
                          func=seir.leaky.dx.dt,
                          parms=my.params.leaky,
                          vac.starts=c(1e6,1e6),
                          vac.ends=c(1e6,1e6),
                          daily.vac=0)

    uncon.run[[r]] <- add.colnames(uncon.run[[r]],"leaky")
    total.vac.days <- 10
    inter.dose.time <- 14 #from last dose

    unvac.states <- "S" #c("S","E","I","R")
    single.vac.states <- "S1" #c("S1","E1","I1","R1")
    double.vac.states <- "S2" #c("S2","E2","I2","R2")

    daily.vac.doses <- runs$dose[r]/total.vac.days

    one.of.one.start <- runs$time[r]
    one.of.one.end <- one.of.one.start+total.vac.days
    one.of.two.start <- runs$time[r]
    one.of.two.end <- one.of.two.start+total.vac.days/2
    two.of.two.start <- one.of.two.end+inter.dose.time
    two.of.two.end <- two.of.two.start + total.vac.days/2

    run.single.vac[[r]] <- ode(y=initial.state,
                               times=times,
                               func=seir.leaky.dx.dt,
                               parms=my.params.leaky,
                               vac.starts=c(one.of.one.start,Inf),
                               vac.ends=c(one.of.one.end,Inf),
                               daily.vac=daily.vac.doses)

    run.single.vac[[r]] <- add.colnames(run.single.vac[[r]],"leaky")

    run.double.vac[[r]] <- ode(y=initial.state,
                               times=times,
                               func=seir.leaky.dx.dt,
                               parms=my.params.leaky,
                               vac.starts=c(one.of.two.start,two.of.two.start),
                               vac.ends=c(one.of.two.end,two.of.two.end),
                               daily.vac=daily.vac.doses)

    run.double.vac[[r]] <- add.colnames(run.double.vac[[r]],"leaky")

}


cols <- c(brewer.pal(3,"Blues"),brewer.pal(3,"Greens"))
par(mfrow=c(1,2))
plot(uncon.run[[1]][,"CI"]/max(uncon.run[[1]][,"CI"]),type="l",col="grey50",xlab="",ylab="")
lines(rowSums(run.single.vac[[1]][,c("CI","CI1")])/max(uncon.run[[1]][,"CI"]),col=cols[6])
lines(rowSums(run.double.vac[[1]][,c("CI","CI1","CI2")])/max(uncon.run[[1]][,"CI"]),col=cols[3])
lines(rowSums(run.single.vac[[2]][,c("CI","CI1")])/max(uncon.run[[1]][,"CI"]),col=cols[6])
lines(rowSums(run.double.vac[[2]][,c("CI","CI1","CI2")])/max(uncon.run[[1]][,"CI"]),col=cols[3])
lines(rowSums(run.single.vac[[3]][,c("CI","CI1")])/max(uncon.run[[1]][,"CI"]),col=cols[6])
lines(rowSums(run.double.vac[[3]][,c("CI","CI1","CI2")])/max(uncon.run[[1]][,"CI"]),col=cols[3])

pdf("Plots/cumcurves_leaky_100kdoses.pdf",width=4,height=2)
par(mar=c(3,3,1,1),mgp=c(1,.1,0),tck=-.02,cex.axis = .8)
plot(uncon.run[[1]][,"CI"]/max(uncon.run[[1]][,"CI"]),type="l",col="grey50",axes=F,
     xlab="epidemic week",
     ylab="cumulative cases")

weeks <- times/7
week.labels <- seq(0,max.weeks,by=8)
week.ats <- approx(x=weeks,y=times,xout=week.labels)$y
axis(1,at=week.ats,labels=week.labels)
axis(2)

grid.labels <- seq(0,max.weeks,by=2)
grid.ats <- approx(x=weeks,y=times,xout=grid.labels)$y
abline(v=grid.ats,col=AddAlpha("grey",.3))
abline(h=seq(0,1,by=.1),col=AddAlpha("grey",.3))


lines(rowSums(run.single.vac[[1]][,c("CI","CI1")])/max(uncon.run[[1]][,"CI"]),col=cols[6])
lines(rowSums(run.double.vac[[1]][,c("CI","CI1","CI2")])/max(uncon.run[[1]][,"CI"]),col=cols[3])
lines(rowSums(run.single.vac[[2]][,c("CI","CI1")])/max(uncon.run[[1]][,"CI"]),col=cols[6])
lines(rowSums(run.double.vac[[2]][,c("CI","CI1","CI2")])/max(uncon.run[[1]][,"CI"]),col=cols[3])
lines(rowSums(run.single.vac[[3]][,c("CI","CI1")])/max(uncon.run[[1]][,"CI"]),col=cols[6])
lines(rowSums(run.double.vac[[3]][,c("CI","CI1","CI2")])/max(uncon.run[[1]][,"CI"]),col=cols[3])
legend("topleft",c("unvaccinated epidemic","1-dose vaccination","2-dose vaccination"),
       col=c("grey50",cols[c(6,3)]),lty=1,
       bty="n",cex=.5)
dev.off()


## ---------------------- ##
## Plots for 2-path model ##
## ---------------------- ##
## NOTE minimization of this objective function was exptrememly figety
## with optim, the bounds were indivually tuned by hand for problematic
## ones. Simple bisection rootfinding should work

pdf("Plots/mrsebypropslow.pdf",width=10,height=6)
## load("GeneratedData/minVES_2path_VE20p80_R1.6-reduced.rda",verbose = T)
## load("GeneratedData/minVES_2path_VE20p88_R1.6-reduced2.rda",verbose = T)
## load("GeneratedData/minVES_2path_VE20p88_R1.6-reduced-0p45max.rda",verbose = T)
## load("GeneratedData/minVES_2path_VE20p88_R1.6-reduced-0p8max.rda",verbose = T)
## load("GeneratedData/minVES_2path_VE20p88_R1.6-reduced-0p75 max.rda",verbose = T)
## load("GeneratedData/minVES_2path_VE20p88_R1.6-reduced-double.rda",verbose = T)
## load("GeneratedData/minVES_2path_VE20p88_R1.6-reduced-double2.rda",verbose = T)
load("GeneratedData/minVES_2path_VE20p88_R1.6-complete.rda",verbose = T)
mins.R1p6 <- mins.betaSAB
nu.seq <- seq(0,1,by=0.1)
cols <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(2,2,2,1))
plot(nu.seq,mins.R1p6[[1]][,1,1]/.88,xlim=c(0,1),ylim=c(0,.75),type="b",col=cols[1],
     xlab="proportion slow transmission",ylab="minimum relative single dose VE")
points(nu.seq,mins.R1p6[[1]][,2,1]/.88,xlim=c(0,1),ylim=c(0,.6),type="b",
       col=cols[1],pch=2)
points(nu.seq,mins.R1p6[[1]][,1,2]/.88,type="b",col=cols[2])
points(nu.seq,mins.R1p6[[1]][,2,2]/.88,type="b",col=cols[2],pch=2)
grid()
legend("bottomright",c("vaccination @ day 0 (1% coverage)",
                       "vaccination @ day 0 (75% coverage)",
                       "vaccination @ peak (1% coverage)",
                       "vaccination @ peak (75% coverage)"),
       col=rep(cols[1:2],each=2),lty=1,pch=rep(1:2,2),bty="n")
#load("GeneratedData/minVES_2path_VE20p88_R3.0-reduced2.rda",verbose = T)
#load("GeneratedData/minVES_2path_VE20p88_R3.0-reduced-0p8 max.rda",verbose = T)
mins.375 <- read.csv("GeneratedData/twopathmins-manual-R3-375k.csv")[,2:3]
load("GeneratedData/minVES_2path_VE20p88_R3.0-reduced-double.rda",verbose = T)
mins.R3p0 <- mins.betaSAB
plot(nu.seq,mins.R3p0[[1]][,1,3]/.88,xlim=c(0,1),ylim=c(0,.75),type="b",col=cols[1],
     xlab="proportion slow transmission",ylab="minimum relative single dose VE",yaxt="n")
axis(4)
points(nu.seq,mins.375[,1]/.88,type="b",col=cols[1],pch=2)
points(nu.seq,mins.R3p0[[1]][,1,4]/.88,type="b",col=cols[2])
points(nu.seq,mins.375[,2]/.88,type="b",col=cols[2],pch=2)
mtext("Proportion Slow Transmission",side=1,outer=T,line=.2)
mtext("Minimum Relative Single Dose VE",side=2,outer=T,line=.2)
mtext("R=1.6",side=3,outer=T,line=-1.5,at=.25)
mtext("R=3.0",side=3,outer=T,line=-1.5,at=.75)
grid()
dev.off()

#######################################
## Let's look at the effect of betas ##
#######################################
source("Source/R/base-functions.R")
## betas to loop over
betas <- seq(.55,2,length=10)
mydoses <- c(1000,500000)
#vac.starts <- c(15,30,60)
#max.time <- 550
palette(brewer.pal(8,"Dark2"))

## ## these were generated from min-aon-varbeta.R etc.
## load("GeneratedData/minVEleaky_varbeta_VE20p80-reduced.rda")
## load("GeneratedData/minVEaon_varbeta_VE20p80-reduced.rda")
## load("GeneratedData/minVES_varbeta_VE20p80-reduced.rda")
## each of these are an array by [time,doses,beta]
load("GeneratedData/minVEaon_varbeta_VE20p88-reduced-2try.rda") #min.aon
load("GeneratedData/minVES_varbeta_VE20p88-reduced2.rda") #min.ves
load("GeneratedData/minVEleaky_varbeta_VE20p88-reduced2.rda") #min.leaky

## now cast in terms of MRSE
min.aon <- min.aon[,,,1]/0.88
## optim was really bad with this so had to do some manually with optim
## with different starts et.
## min.aon[1,2,1] <- 0.3415/.88
## min.aon[2,2,1] <- 0.334/.88
## min.aon[1,2,3] <- 0.4291/.88
## min.aon[1,2,4] <- 0.4398/.88
## min.aon[3,2,7] <-  # min
## min.aon[3,1,8]
## min.aon[3,2,8]
## min.aon[2,1,9]
## min.aon[2,2,9]
## min.aon[3,1,9]
## min.aon[3,2,9]
## min.aon[2,1,10]
## min.aon[2,2,10]
## min.aon[3,1,10]
## min.aon[3,2,10]


min.ves <- min.ves/0.88
min.ves[1,2,3]<- 0.4293/.88
min.ves[2,2,3]<- 0.4058/.88
min.ves[1,2,4]<- 0.4721/.88

## min.leaky[1,2,2] <-  0.3674
## min.leaky[1,2,3]<- 0.3999
## min.leaky[1,2,4] <- 0.3999
## min.leaky[1,1,9] <- 0.3495

min.leaky <- min.leaky/0.88


pdf("Plots/varbetaplot-mrse-500k.pdf",width=15,height = 5)
par(mfrow=c(1,3),mar=c(1,1,1,1),oma=c(4,4,3,1),tck=-0.01,mgp=c(1.5,0.3,0))
plot(-100,-100,xlim=c(1,4),ylim=c(0,.75),xlab="R0",ylab="MRSE")
grid()
lines(betas/.5,min.ves[1,2,],col=1,type="b")
lines(betas/.5,min.ves[2,2,],col=1,type="b")
lines(betas/.5,min.ves[3,2,],col=1,type="b")
text(3.18,.62,"day 0")
text(3.2,.4,"half-way to peak")
text(1.5,.2,"at peak")
#dev.off()

#pdf("Plots/vesleakybyR.pdf",width=5,height = 3)
#par(mar=c(3,3,1,1),tck=-0.01,mgp=c(1.5,0.3,0))
plot(-100,-100,xlim=c(1,4),ylim=c(0,.8),xlab="R0",ylab="MRSE")
grid()
lines(betas/.5,min.leaky[1,2,],col=2,type="b")
lines(betas/.5,min.leaky[2,2,],col=2,type="b")
lines(betas/.5,min.leaky[3,2,],col=2,type="b")
text(3.4,.52,"day 0")
text(2.85,.33,"half-way to peak")
text(1.6,.2,"at peak")
#dev.off()

#pdf("Plots/aonbyR.pdf",width=5,height = 3)
#par(mar=c(3,3,1,1),tck=-0.01,mgp=c(1.5,0.3,0))
plot(-100,-100,xlim=c(1,4),ylim=c(0,.8),xlab="R0",ylab="MRSE")
grid()
## taking out the first point because this in these scenarios so much vaccine
## leads to no epidemic in both situations
lines(betas/.5,min.aon[1,2,],col=3,type="b")
lines(betas/.5,min.aon[2,2,],col=3,type="b")
lines(betas/.5,min.aon[3,2,],col=3,type="b")
text(2.7,.525,"day 0")
text(2.7,.3,"half-way to peak")
text(1.55,.2,"at peak")

mtext("Reproductive Number",side=1,outer=T,line=1)
mtext("MRSE",side=2,outer=T,line=1)
mtext("Susceptibility-Reducing \n Vaccine",side=3,outer=T,line=-1,at=.16)
mtext("Severity-Reducing \n Vaccine",side=3,outer=T,line=-1,at=.5)
mtext("All-or-Nothing \n Vaccine",side=3,outer=T,line=-1,at=.82)
dev.off()

## similar but with only a single coverage level shown ##
pdf("Plots/varbetaplot-mrse.pdf",width=15,height = 5)
par(mfrow=c(1,3),mar=c(1,1,1,1),oma=c(4,4,3,1),tck=-0.01,mgp=c(1.5,0.3,0))
plot(-100,-100,xlim=c(1,4),ylim=c(0,.75),xlab="R0",ylab="MRSE")
grid()
polygon(x=c(betas/.5,rev(betas/.5)),
        y=c(min.ves[1,2,],rev(min.ves[1,1,])),col=AddAlpha(1,.7),border=FALSE)
polygon(x=c(betas/.5,rev(betas/.5)),
        y=c(min.ves[2,2,],rev(min.ves[2,1,])),col=AddAlpha(1,.5),border=FALSE)
polygon(x=c(betas/.5,rev(betas/.5)),
        y=c(min.ves[3,2,],rev(min.ves[3,1,])),col=AddAlpha(1,.3),border=FALSE)
text(3.1,.62,"day 0")
text(3.0,.4,"half-way to peak")
text(1.7,.2,"at peak")
#dev.off()

#pdf("Plots/vesleakybyR.pdf",width=5,height = 3)
#par(mar=c(3,3,1,1),tck=-0.01,mgp=c(1.5,0.3,0))
plot(-100,-100,xlim=c(1,4),ylim=c(0,.8),xlab="R0",ylab="MRSE")
grid()
polygon(x=c(betas/.5,rev(betas/.5)),
        y=c(min.leaky[1,2,],rev(min.leaky[1,1,])),col=AddAlpha(2,.7),border=FALSE)
polygon(x=c(betas/.5,rev(betas/.5)),
        y=c(min.leaky[2,2,],rev(min.leaky[2,1,])),col=AddAlpha(2,.5),border=FALSE)
polygon(x=c(betas/.5,rev(betas/.5)),
        y=c(min.leaky[3,2,],rev(min.leaky[3,1,])),col=AddAlpha(2,.3),border=FALSE)
text(3.4,.41,"day 0")
text(2.6,.33,"half-way to peak")
text(1.6,.2,"at peak")
#dev.off()

#pdf("Plots/aonbyR.pdf",width=5,height = 3)
#par(mar=c(3,3,1,1),tck=-0.01,mgp=c(1.5,0.3,0))
plot(-100,-100,xlim=c(1,4),ylim=c(0,.8),xlab="R0",ylab="MRSE")
grid()
## taking out the first point because this in these scenarios so much vaccine
## leads to no epidemic in both situations
polygon(x=c(c(betas/.5),rev(betas/.5)),
        y=c(min.aon[1,2,],rev(min.aon[1,1,])),col=AddAlpha(3,.7),border=FALSE)
polygon(x=c(c(betas/.5),rev(betas/.5)),
        y=c(min.aon[2,2,],rev(min.aon[2,1,])),col=AddAlpha(3,.5),border=FALSE)
polygon(x=c(c(betas/.5)[],rev(betas/.5)),
        y=c(min.aon[3,2,],rev(min.aon[3,1,])),col=AddAlpha(3,.3),border=FALSE)
text(2.7,.37,"day 0")
text(2.5,.3,"half-way to peak")
text(1.55,.2,"at peak")

mtext("Reproductive Number",side=1,outer=T,line=1)
mtext("MRSE",side=2,outer=T,line=1)
mtext("Susceptibility Reducing \n Vaccine",side=3,outer=T,line=-1,at=.16)
mtext("Severity Reducing \n Vaccine",side=3,outer=T,line=-1,at=.5)
mtext("All-or-Nothing \n Vaccine",side=3,outer=T,line=-1,at=.82)
dev.off()

## ------------------------------------- ##
## Looking at minimum relative efficiacy ##
## ------------------------------------- ##
pdf("Plots/mrse-justification-byvactype-37to88.pdf",
    width=10,height=4)
#quartz("",width=10,height=4)
ve2 <- c(.37,.75,.88)
par(mfrow=c(1,3),oma=c(3,3,3,3),mar=c(2,2,2,2))
palette(brewer.pal(8,"Dark2"))
#load("GeneratedData/minVEaon-para_betaSAB-reduced.rda")
load("GeneratedData/minVEaon-para_betaSAB-reduced-37to88.rda")
vac.times <- seq(0,300,by=20)
plot(vac.times,min.aon[,1,1]/ve2[1],
     col=1,ylim=c(.1,.6),
     ylab="minimum relative efficacy",
     xlab="vaccination start day")
grid()

points(vac.times,min.aon[,1,2]/ve2[2],col=2)
points(vac.times,min.aon[,1,3]/ve2[3],col=3)

points(vac.times,min.aon[,2,1]/ve2[1],pch=2,col=1)
points(vac.times,min.aon[,2,2]/ve2[2],pch=2,col=2)
points(vac.times,min.aon[,2,3]/ve2[3],pch=2,col=3)

points(vac.times,min.aon[,3,1]/ve2[1],pch=3,col=1)
points(vac.times,min.aon[,3,2]/ve2[2],pch=3,col=2)
points(vac.times,min.aon[,3,3]/ve2[3],pch=3,col=3)

## legend("topright",c("ve2=.55, 1k doses","ve2=.75, 1k doses","ve2=.88, 1k doses",
##                     "ve2=.55, 100k doses","ve2=.75, 100k doses","ve2=.88, 100k doses",
##                     "ve2=.55, 500k doses","ve2=.75, 500k doses","ve2=.88, 500k doses"),
##                     ,col=rep(1:3,3),pch=rep(1:3,each=3),bty="n")
## dev.off()

## for leaky
#load("GeneratedData/minVEleaky-para_betaSAB-reduced.rda")
load("GeneratedData/minVEleaky-para_betaSAB-reduced-37to88.rda")
#pdf("Plots/leaky-mre-reduced.pdf")
plot(vac.times,min.leaky[,1,1]/ve2[1],
     col=1,ylim=c(.1,.6),
     ylab="minimum relative efficacy",
     xlab="vaccination start day")
grid()

points(vac.times,min.leaky[,1,2]/ve2[2],col=2)
points(vac.times,min.leaky[,1,3]/ve2[3],col=3)

points(vac.times,min.leaky[,2,1]/ve2[1],pch=2,col=1)
points(vac.times,min.leaky[,2,2]/ve2[2],pch=2,col=2)
points(vac.times,min.leaky[,2,3]/ve2[3],pch=2,col=3)

points(vac.times,min.leaky[,3,1]/ve2[1],pch=3,col=1)
points(vac.times,min.leaky[,3,2]/ve2[2],pch=3,col=2)
points(vac.times,min.leaky[,3,3]/ve2[3],pch=3,col=3)

## legend("topright",c("ve2=.55, 1k doses","ve2=.75, 1k doses","ve2=.88, 1k doses",
##                     "ve2=.55, 100k doses","ve2=.75, 100k doses","ve2=.88, 100k doses",
##                     "ve2=.55, 500k doses","ve2=.75, 500k doses","ve2=.88, 500k doses"),
##                     ,col=rep(1:3,3),pch=rep(1:3,each=3),bty="n")
## dev.off()

load("GeneratedData/minVEs-para_betaSAB-reduced-37to88.rda")
#load("GeneratedData/minVEs-para_betaSAB-reduced.rda")

## pdf("Plots/ves-mre-reduced.pdf")
plot(vac.times,min.ves[,1,1]/ve2[1],
     col=1,ylim=c(.1,.6),
     ylab="minimum relative efficacy",
     xlab="vaccination start day")
grid()
points(vac.times,min.ves[,1,2]/ve2[2],col=2)
points(vac.times,min.ves[,1,3]/ve2[3],col=3)

points(vac.times,min.ves[,2,1]/ve2[1],pch=2,col=1)
points(vac.times,min.ves[,2,2]/ve2[2],pch=2,col=2)
points(vac.times,min.ves[,2,3]/ve2[3],pch=2,col=3)

points(vac.times,min.ves[,3,1]/ve2[1],pch=3,col=1)
points(vac.times,min.ves[,3,2]/ve2[2],pch=3,col=2)
points(vac.times,min.ves[,3,3]/ve2[3],pch=3,col=3)

legend("topright",c("ve2=.37, 1k doses","ve2=.75, 1k doses","ve2=.88, 1k doses",
                    "ve2=.37, 100k doses","ve2=.75, 100k doses","ve2=.88, 100k doses",
                    "ve2=.37, 500k doses","ve2=.75, 500k doses","ve2=.88, 500k doses"),
                    ,col=rep(1:3,3),pch=rep(1:3,each=3),bty="n",cex=.8)

mtext("All-or-nothing",at=1/3/2,outer=T,line=-1)
mtext("VE_SP",at=2/3 - 1/3/2,outer=T,line=-1)
mtext("VE_S",at=1 - 1/3/2,outer=T,line=-1)
mtext("Vaccine Start Day", side=1,at = .5,line=.5,outer=T)
mtext("Minimum Relative Single-Dose Effectiveness", side=2,at = .5,line=.5,outer=T)
dev.off()


## ------------------------------------------ ##
## Exploring the impact of interdose duration ##
## ------------------------------------------ ##
load("GeneratedData/minVEleaky_betaSAB_VE20p88-1kto500k-interdosetimes.rda") #min.leaky
dim(min.leaky) # times, doses, interdose period
times <- seq(0,210,by=2)
cols <- colorRampPalette(palette(brewer.pal(8,"Blues")))(28)

## get uncontrolled epidemic curve from other simulations
load("GeneratedData/vesruns_VE20p88_betaSAB.rda")
uncon.run <- my.biglist[[7]]
uncon.epi <- diff(uncon.run[,"CI"])
pdf("Plots/interdoseduration.pdf",width=10,height=5)
layout(matrix(c(1,2,2,2),nrow=4))
par(oma=c(4,4,4,4),mar=c(1,1,1,1))
plot(0:210,uncon.epi[1:211],type="l",xlab="",ylab="cholera infections",axes=FALSE)
axis(3)
axis(2)
grid()
mtext("cholera infections",side=2,outer=T,at=.85,cex=.75,line=1)
plot(times,min.leaky[,1,1]/.88,ylim=c(0.1,.5),col=cols[1],type="l",axes=FALSE)
sapply(2:28,function(x) lines(times,min.leaky[,1,x]/.88,col=cols[x]))
axis(1)
axis(2)
grid()
mtext("epidemic day",side=1,outer=T,at=.5,cex=.75,line=1.5)
mtext("MRSE",side=2,outer=T,at=.4,cex=.75,line=1.5)
dev.off()
## ------------------------------------------- ##
##  Exploring the impact of time to protection ##
##  on MRSE                                    ##
## ------------------------------------------- ##

load("GeneratedData/minVEleaky_betaSAB_VE20p88-timetoprot-1.0.rda")
mins <- c(min.leaky)
load("GeneratedData/minVEleaky_betaSAB_VE20p88-timetoprot-0.8.rda")
mins <- cbind(mins,min.leaky)
load("GeneratedData/minVEleaky_betaSAB_VE20p88-timetoprot-0.6.rda")
mins <- cbind(mins,min.leaky)
load("GeneratedData/minVEleaky_betaSAB_VE20p88-timetoprot-0.4.rda")
mins <- cbind(mins,min.leaky)
load("GeneratedData/minVEleaky_betaSAB_VE20p88-timetoprot-0.2.rda")
mins <- cbind(mins,min.leaky)

times <- seq(0,300,by=2)
cols <- rev(palette(brewer.pal(8,"Greens"))[3:8])

pdf("Plots/timetoprotection.pdf",width=10,height=5)
layout(matrix(c(1,2,2,2),nrow=4))
par(oma=c(4,4,4,4),mar=c(1,1,1,1))
plot(0:300,uncon.epi[1:301],type="l",xlab="",ylab="cholera infections",axes=FALSE)
axis(3)
axis(2)
grid()
mtext("cholera infections",side=2,outer=T,at=.85,cex=.75,line=1)
plot(times,rep(0.44,length(times))/.88,ylim=c(0.1,.5),col=cols[6],type="l",axes=FALSE)
sapply(1:5,function(x) lines(times,mins[,x]/.88,col=cols[(x+1)]))
axis(1)
axis(2)
grid()
mtext("vaccination start day",side=1,outer=T,at=.5,cex=.75,line=1.5)
mtext("MRSE",side=2,outer=T,at=.4,cex=.75,line=1.5)
legend("bottomleft",paste0(14*seq(1,0,by=-0.2), " days to protection after second dose"),col=cols,lty=1,bty="n")
dev.off()

plot(mins[,1]/.88,type="l",ylim=c(.1,.5))
sapply(2:5,function(x) lines(mins[,x]/.88,col=cols[x]))






require(pomp)
pompExample(gompertz)
coef(gompertz)
mf <- mif(gompertz,
          start=coef(gompertz),
          pars=c('tau','r','X.0'),
          rw.sd=c(tau=0.1,r=.1,sigma=0.1,X.0=0.1),
          ivps="X.0",
          Nmif=5,
          Np=500,
          var.factor=1,
          ic.lag=10,
          cooling.type="hyperbolic",
          cooling.fraction=0.03,
          method="mif2",
          transform=T
    )
 coef(mf)
