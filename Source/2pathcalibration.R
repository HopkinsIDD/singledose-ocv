## ----------------------------------------------- ##
## Script for playing around with 2 path VES model ##
## ----------------------------------------------- ##
source("Source/base-functions.R")

## read in virbio survival data
vib.surv <- read.csv("Data/Li_vibriosurvivial.csv")
plot(vib.surv[,1],vib.surv[,2]/18,xlab="day",ylab="Proportion Culturable")

## reducing sum of squares to emperical survival curve
obj.func <- function(log.pars){
    dat <- vib.surv[,2]/18
    model <- 1-pgamma(vib.surv[,1],rate=exp(log.pars[1]),shape=exp(log.pars[2]))
    lines(vib.surv[,1],model,col=AddAlpha("red",.2))
    sum((dat-model)^2)
}

opt <- optim(par=c(log(2),log(10)),fn=obj.func)
exp(opt$par)
1/exp(opt$par)[1] # 7.5 day duration
round(exp(opt$par))[2] # 3 compartments
curve(1-pgamma(x,rate=exp(opt$par[1]),shape=exp(opt$par[2])),0,60,add=T)

pdf("Plots/LiVibrioSurv.pdf")
plot(vib.surv[,1],vib.surv[,2]/18,xlab="day",ylab="Proportion of Vibrio Isolates Culturable")
curve(1-pgamma(x,rate=exp(opt$par[1]),shape=exp(opt$par[2])),0,60,add=T,lty=2)
abline(h=seq(0,1,.1),col=AddAlpha("grey",.1))
abline(v=seq(0,60,5),col=AddAlpha("grey",.1))
legend("topright",c("data from Li et al., 2014","gamma model fit"),lty=c(0,2),pch=c(1,-10),bty="n")
dev.off()


pdf("Plots/marginalinfectiousperiods.pdf")
curve(dexp(x,rate=1/2),0,80,xlab="days",ylab="density",lty=1)
curve(dgamma(x,rate=exp(opt$par[1]),shape=exp(opt$par[2])),0,80,lty=2,add=T)
legend("topright",c("Direct Transmission Only","Environmental Transmission Only"),
       lty=1:2,bty="n")
title("Marginal Infectious Period Distributions \n in Mixed Mode Cholera Model")
dev.off()

## fast-slow VES model
## set up parameters
ve2 <- .80
params <- list(beta=0.8,
               beta.long =  1.6/(3*7.5),
               gamma=1/2,
               gamma.s=1/7.5,
               sigma=1/1.4,
               phi1=0,
               phi2=0,
               theta1=0.4,
               theta2=ve2,
               nu=.5 #proportion slow
               )

## initial state one seed in each pathway
initial.state <- c(1e5-1,0,0, #S, S1, S2
                   0,0,0,     #E, E1, E2
                   1,0,0,     #I, I1, I2
                   1,0,0,     #I*, I*-1, I*-2
                   0,0,0,     #I1*, I1*-1, I1*-2
                   0,0,0,     #I2*, I2*-1, I2*-2
                   0,0,0,     #R,R1,R2
                   0,0,0,     #V1,V2,W
                   0,0,0      #CI,CI1,CI2
                   )

## simple test run to make sure things are working
test <- ode(initial.state,
            times=seq(0,500,by=1),
            func=ves.2path.dx.dt,
            parms=params,
            vac.starts=c(10,19),
            vac.ends=c(15,23),
            daily.vac=5000,
            n.comps.slow=3)

test <- add.colnames(test,model.type = "2path")
plot(test)

##' now we are going to find the peak day for fixed paramers exepct
##' nu (proportion slow)

nu.seq <- seq(0,1,by=.05)
timings <- array(dim=c(length(nu.seq),2)) # for peak time and end time
for (n in seq_along(nu.seq)){
    params$nu <- nu.seq[n]
    tmp <- ode(initial.state,
               times=seq(0,1000,by=1),
               func=ves.2path.dx.dt,
               parms=params,
               vac.starts=c(1e6,1e6),
               vac.ends=c(1e7,1e7),
               daily.vac=5000,
               n.comps.slow=3)
    tmp <- add.colnames(tmp,model.type = "2path")
    epi.curve <- diff(tmp[,"CI"])
    timings[n,1] <- which.max(epi.curve) #peak
    timings[n,2] <- min(which(epi.curve < .1))
}
colnames(timings) <- c("peak","total")
save(timings,file="GeneratedData/twopath-timings-R1p6.rda")

pdf("Plots/fastslowtimings-R1p6.pdf",width=6,height=4)
par(mfrow=c(1,2))
plot(nu.seq,timings[,1],xlab="Proportion Slow (nu)",ylab="Peak Time (day)")
plot(nu.seq,timings[,2],xlab="Proportion Slow (nu)",ylab="Epidemic Duration (day)")
title("Timings of Uncontrolled Fast/Slow Epidemics (R=1.6)",out=T,line=-2)
dev.off()


cols <- colorRampPalette(brewer.pal(8,"Oranges")[-1])(21)
pdf("Plots/epicurves_proportionslow_R1p6.pdf")
plot(-100,-100,xlim=c(0,500),ylim=c(0,2500),xlab="day",ylab="incident cases")
for (n in seq_along(nu.seq)){
    params$nu <- nu.seq[n]
    tmp <- ode(initial.state,
               times=seq(0,1000,by=1),
               func=ves.2path.dx.dt,
               parms=params,
               vac.starts=c(1e6,1e6),
               vac.ends=c(1e7,1e7),
               daily.vac=5000,
               n.comps.slow=3)
    tmp <- add.colnames(tmp,model.type = "2path")
    epi.curve <- diff(tmp[,"CI"])
    lines(epi.curve,col=cols[n])
    timings[n,1] <- which.max(epi.curve) #peak
    timings[n,2] <- min(which(epi.curve < .1))
}
dev.off()


## now for R=3
beta <- 1.5#
gamma <- 1/2
gamma.s <- 1/7.5
R0 <- beta/gamma
beta.long <- R0/(3*1/gamma.s)

params <- list(beta=beta,
               beta.long = beta.long,
               gamma=gamma,
               gamma.s=gamma.s,
               sigma=1/1.4,
               phi1=0,
               phi2=0,
               theta1=0.4,
               theta2=0.8,
               nu=.5 #proportion slow
               )

# one infected in each path
initial.state <- c(.5e6-1,0,0, #S, S1, S2
                   0,0,0,     #E, E1, E2
                   1,0,0,     #I, I1, I2
                   1,0,0,     #I*, I*-1, I*-2
                   0,0,0,     #I1*, I1*-1, I1*-2
                   0,0,0,     #I2*, I2*-1, I2*-2
                   0,0,0,     #R,R1,R2
                   0,0,0,     #V1,V2,W
                   0,0,0      #CI,CI1,CI2
                   )

nu.seq <- seq(0,1,by=.05)
timings <- array(dim=c(length(nu.seq),2)) # for peak time and end time
for (n in seq_along(nu.seq)){
    params$nu <- nu.seq[n]
    tmp <- ode(initial.state,
               times=seq(0,1000,by=1),
               func=ves.2path.dx.dt,
               parms=params,
               vac.starts=c(1e6,1e6),
               vac.ends=c(1e7,1e7),
               daily.vac=5000,
               n.comps.slow=3)
    tmp <- add.colnames(tmp,model.type = "2path")
    epi.curve <- diff(tmp[,"CI"])
    timings[n,1] <- which.max(epi.curve) #peak
    timings[n,2] <- min(which(epi.curve < .1))
}
colnames(timings) <- c("peak","total")
save(timings,file="GeneratedData/twopath-timings-R3p0.rda")


pdf("Plots/epicurves_proportionslow_R3p0.pdf")
plot(-100,-100,xlim=c(0,200),ylim=c(0,10000),xlab="day",ylab="incident cases")
for (n in seq_along(nu.seq)){
    params$nu <- nu.seq[n]
    tmp <- ode(initial.state,
               times=seq(0,1000,by=1),
               func=ves.2path.dx.dt,
               parms=params,
               vac.starts=c(1e6,1e6),
               vac.ends=c(1e7,1e7),
               daily.vac=5000,
               n.comps.slow=3)
    tmp <- add.colnames(tmp,model.type = "2path")
    epi.curve <- diff(tmp[,"CI"])
    lines(epi.curve,col=cols[n])
}
dev.off()
