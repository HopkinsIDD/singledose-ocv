###################################################################
## An attempt to really nail down unvac fits for zim and conakry ##
###################################################################
source("Source/R/leakyvac-pomp-model-inC-novac.R")
library(pomp)
palette(brewer.pal(8,"Set1"))
set.seed(243947892)

## for some parallel mif runs
library(foreach)
#library(multicore)
library(doMC)
registerDoMC(2)

## ------------------------- ##
## Let's start wiht zimbabwe ##
## ------------------------- ##

## zim popualtion
pop.zim <- 13.4e6

## zimdat
zim.dat <- get.zim.data()[-50,-1]
zim.dat$week <- 1:49
colnames(zim.dat)

## build pomp model object
zim.mod <- build.leaky.model.C(pop=pop.zim,
                               dat = zim.dat,
                               model.name = "zimmod")

## specify starting parameters
## remember these are in units of weeks
E0 <- 10/pop.zim
I0 <- 10/pop.zim
A0 <- 1e-11
R0 <- 0.36
S0 <- 1- R0-I0-E0-A0
guess.params <- c(gamma=1.78,
                  sigma=5,
                  theta=10,
                  beta1=3.30,
                  beta2=6,
                  rho=0.03,
                  theta0=0,
                  S.0=S0,
                  E.0=E0,
                  I.0=I0,
                  A.0=A0,
                  R.0=R0)


#zim.mod.win <- window(zim.mod,start=14,end=49)
#zim.mod.win@t0 <- 13

## first let's start with Trajectory matching
set.seed(19821)
tm.zim <- traj.match(zim.mod,
                     start=guess.params,
                     est=c('beta1','R.0','I.0'),
                     method="subplex",
                     transform=TRUE)
summary(tm.zim)
sim.zim.tm <- simulate(tm.zim,
                       nsim=500,
                       seed=1914679109L,
                       transform=TRUE)

## let's plot
plot(zim.dat[,1],ylim=c(0,20000))

for (i in 1:500) {
    lines(sim.zim.tm[[i]]@data[1,],lty=2,col=AddAlpha(5,.05))
}

## run a pfilter to look at likelihood
zim.pf <- pfilter(tm.zim,
                  Np=1000,
                  save.params=TRUE)
logLik(zim.pf) #-377

## let's mif
set.seed(19822)
mif.zim <- mif(tm.zim,
                start=coef(tm.zim),
                Nmif=50,
                pars=c('beta1','gamma'),
                ivps=c('I.0','S.0'),
                transform=TRUE,
                rw.sd=c(beta1=0.12,
                    I.0=0.12,
                    S.0=0.12,
                    gamma=0.12),
                Np=2000,
                ic.lag=length(zim.mod@data)/2,
                var.factor=1,
                cooling.type="hyperbolic",
                cooling.fraction=0.05,
                method="mif2")


saveRDS(mif.zim,file="GeneratedData/mif-zim-REV.rds")
mif.zim <- readRDS("GeneratedData/mif-zim-REV.rds")
mif.zim <- continue(mif.zim,Nmif=50)

logLik(mif.zim)
logLik(pfilter(mif.zim,Np=5000))

## now let's explore the space around to make sure we aren't stuck in
## a local maximum
estpars <- c('beta1','gamma')
set.seed(19823)
mf.zim <- foreach(i=1:10,
              .inorder=FALSE,
              .options.multicore=list(set.seed=TRUE)
              ) %dopar%
{
    ## let's saple our parameters
    theta.guess <- coef(tm.zim)
    theta.guess[estpars] <- rlnorm(
        n=length(estpars),
        meanlog=log(theta.guess[estpars]),
        sdlog=0.2
        )

    # now sample from I.0
    I.0.count <- runif(1,1,1e4)/pop.zim # people
    theta.guess['S.0'] <- theta.guess['S.0'] - I.0.count
    theta.guess['I.0'] <- I.0.count
    theta.guess['S.0'] <-  theta.guess['S.0']*runif(1,.8,1.2)
    theta.guess['R.0'] <- max(0,1-sum(theta.guess[c('S.0','I.0','E.0','A.0')]))

    m1 <- mif(
        tm.zim,
        Nmif=100,
        start=theta.guess,
        pars=c('beta1','gamma'),
        ivps=c('I.0','S.0'),
        transform=TRUE,
        rw.sd=c(beta1=0.1,gamma=0.1,I.0=0.1,S.0=0.1),
        Np=3000,
        ic.lag=length(zim.mod@data),
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
    )
    ll <- replicate(n=3,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

sapply(mf.zim,function(x) x[[2]])
## compare.mif now depreciated now using mifList object
zim.mflist <- do.call(c,lapply(mf.zim,function(x) x[[1]]))
saveRDS(mf.zim,file="GeneratedData/parallel-mif-zim-REV.rds")
to.pdf(plot(zim.mflist),"Plots/unvac-zimmif-REV.pdf")

## now explore the best mifs a little further
mf.zim <- readRDS("GeneratedData/parallel-mif-zim-REV.rds")
#mf.zim <- readRDS("GeneratedData/parallel-mif-zim.rds")

mult.pomps <- sapply(mf.zim,function(x) x[[1]])
plot(do.call(c,lapply(mf.zim,function(x) x[[1]])))

## compare.mif(mult.pomps)
## let's look at which ones are best
best.pomps <- order(colMeans(sapply(mf.zim,function(x) x[[2]])),decreasing=TRUE)

set.seed(19823)
better.mif.zim <- mif(mult.pomps[[best.pomps[1]]],
                      Nmif=100,
                      pars=c('beta1','gamma'),
                      ivps=c('I.0','S.0'),
                      transform=TRUE,
                      rw.sd=c(beta1=0.12,I.0=0.12,S.0=0.12,gamma=0.12),
                      Np=2000,
                      ic.lag=length(zim.mod@data)/2,
                      var.factor=1,
                      cooling.type="hyperbolic",
                      cooling.fraction=0.05,
                      method="mif2")

logLik(pfilter(better.mif.zim,Np=10000))
to.pdf(plot(better.mif.zim),"Plots/mif-zim-final-diag.pdf")

## and plot
set.seed(4278435)
sim.zim.mif <- simulate(better.mif.zim,
                        nsim=500,
                        transform=TRUE)


## wrapper for plotting simulations and zim data with PIs
make.zim.sim.plot <- function(run,dat,nlines=500){
    zim.mat <- sapply(run,function(x) x@data[1,])
    zim.means <- apply(zim.mat,1,mean)
    zim.ci <- apply(zim.mat,1,function(x) quantile(x,c(.025,.975)))
    plot(dat[,1],ylim=c(0,14000),xlab="epidemic week",ylab="cases per week",pch=4)
    for (i in 1:nlines) {
        lines(sim.zim.mif[[i]]@data[1,],lty=2,col=AddAlpha(4,.05))
    }
    lines(zim.means,col=4)
    lines(zim.ci[1,],col=4,lty=2)
    lines(zim.ci[2,],col=4,lty=2)
    legend("topright",c("simulated epidemic",
                        "mean simulated epidemic",
                        "95% Prediction Interval",
                        "data"),
           col=c(AddAlpha(4,0.1),4,4,"black"),lty=c(1,1,2,-1),pch=c(-1,-1,-1,4),bty="n")
}

## make the pdf
to.pdf(make.zim.sim.plot(sim.zim.mif,zim.dat),"Plots/mif-zim-unvac-REV.pdf")

pdf("Plots/hist-finalsize-uncon-zim-REV.pdf")
hist(colSums(sapply(sim.zim.mif,function(x) x@data[1,])),
     col="grey",border="white",breaks="fd",
     xlab="Final Epidemic Size of Simulation",
     main="Final Size of Zimbabwe Simulations",
     xlim=c(80000,140000))
abline(v=98591,col="orange",lwd=2,lty=2)
text(97000,50,"Reported Epidemic \n Size = 98,591",cex=.9)
dev.off()

## ------------------------------------------------------ ##
## NOTE: we are going to use this fit for our projections ##
## ------------------------------------------------------ ##
saveRDS(better.mif.zim,file="GeneratedData/mif-zim-REV.rds")
better.mif.zim <- readRDS("GeneratedData/mif-zim-REV.rds")
## and save the final states from our particle filter
pf.zim <- pfilter(better.mif.zim,Np=5000,save.states=TRUE)
est.states.zim <- sapply(pf.zim@saved.states,rowMeans)
saveRDS(est.states.zim,file="GeneratedData/mif-zim-states-REV.rds")


## --------------------------------------------------------------- ##
## Let's do a little profiling of our beta and gamma parameters to ##
## see how peaky they look                                         ##
## --------------------------------------------------------------- ##

## first for beta
beta.range <- seq(coef(better.mif.zim)['beta1']*.6,coef(better.mif.zim)['beta1']*1.4,length=15)
mf.zim.beta.prof <- foreach(i=1:length(beta.range),
                            .inorder=FALSE,
                            .options.multicore=list(set.seed=TRUE)
                            ) %dopar%
{

    theta.guess <- coef(better.mif.zim)
    theta.guess['beta1'] <- beta.range[i]

    m1 <- mif(
        better.mif.zim,
        Nmif=50,
        start=theta.guess,
        pars=c('gamma'),
        ivps=c('I.0','S.0'),
        transform=TRUE,
        rw.sd=c(gamma=0.1,I.0=0.1,S.0=0.1),
        Np=2000,
        ic.lag=length(zim.mod@data)/2,
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
        )
    ll <- replicate(n=3,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

beta.logliks <- colMeans(sapply(mf.zim.beta.prof,function(x) x[[2]]))
cis <- max(beta.logliks) - qchisq(.95,1)/2
pdf("Plots/proflik-beta-zim-REV.pdf")
plot(beta.range,beta.logliks,xlab="beta",ylab="log-likelihood",main="Profile Likelihood of Beta (Zimbabwe)")
abline(h=cis,lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(beta.logliks[2:6],beta.range[2:6],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(beta.logliks[8:13],beta.range[8:13],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
text(4,-336,sprintf("95%% CI %.2f-%.2f",
                       approx(beta.logliks[2:6],beta.range[2:6],xout=cis)$y,
                       approx(beta.logliks[8:12],beta.range[8:12],xout=cis)$y))
dev.off()
## approximate 95% CI

# gamma.range <- seq(coef(better.mif.zim)['gamma']*.6,coef(better.mif.zim)['gamma']*2,length=25)
gamma.range <- seq(2.2,5,length=25)

mf.zim.gamma.prof <- foreach(i=1:length(gamma.range),
                             .inorder=FALSE,
                             .options.multicore=list(set.seed=TRUE)
                             ) %dopar%
{

    theta.guess <- coef(better.mif.zim)
    theta.guess['gamma'] <- gamma.range[i]

    m1 <- mif(
        better.mif.zim,
        Nmif=50,
        start=theta.guess,
        pars=c('beta1'),
        ivps=c('I.0','S.0'),
        transform=TRUE,
        rw.sd=c(beta1=0.1,I.0=0.1,S.0=0.1),
        Np=2000,
        ic.lag=length(zim.mod@data)/2,
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
        )
    ll <- replicate(n=3,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

saveRDS(gamma.logliks,file="GeneratedData/gamma_proflik.rds")
plot(gamma.range,colMeans(sapply(mf.zim.gamma.prof,function(x) x[[2]])))
gamma.logliks <- colMeans(sapply(mf.zim.gamma.prof,function(x) x[[2]]))
cis <- max(gamma.logliks) - qchisq(.95,1)/2
pdf("Plots/proflik-gamma-zim-REV.pdf")
plot(gamma.range,gamma.logliks,xlab="gamma",ylab="log-likelihood",main="Profile Likelihood of Gamma (Zimbabwe)")
abline(h=cis,lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(gamma.logliks[5:7],gamma.range[5:7],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(gamma.logliks[20:25],gamma.range[20:25],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
text(3.5,-337.5,sprintf("95%% CI %.2f-%.2f",
                        approx(gamma.logliks[5:7],gamma.range[5:7],xout=cis)$y,
                        approx(gamma.logliks[20:25],gamma.range[20:25],xout=cis)$y))
dev.off()


## NOW Let's calculate the joint profile for beta and gamma to get
## R

beta.seq <- seq(2.5,6,length=30)
gamma.seq <- seq(2,5,length=30)
r.seq <- expand.grid(beta.seq,gamma.seq)

R.prof <- foreach(i=1:nrow(r.seq),
                             .inorder=FALSE,
                             .options.multicore=list(set.seed=TRUE)
                             ) %dopar%
{

    theta.guess <- coef(better.mif.zim)
    theta.guess['gamma'] <- r.seq[i,2]
    theta.guess['beta1'] <- r.seq[i,1]

    m1 <- mif(
        better.mif.zim,
        Nmif=50,
        start=theta.guess,
        ivps=c('I.0','S.0'),
        transform=TRUE,
        rw.sd=c(I.0=0.1,S.0=0.1),
        Np=2000,
        ic.lag=length(zim.mod@data)/2,
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
    )
    ll <- replicate(n=3,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

saveRDS(R.logliks,file="GeneratedData/R_proflik.rds")

## get CI for R
prof.lik <- R.logliks
prof.mat <- matrix(prof.lik[,3],nrow=30)
colnames(prof.mat) <- beta.seq
rownames(prof.mat) <- gamma.seq

ci.lines <- contourLines(prof.mat*2,levels=max(prof.mat)*2 - 3.814/2)
## we will only take the middle
range(approx(seq(0,1,length=30),beta.seq,ci.lines[[2]]$x)$y/approx(seq(0,1,length=30),gamma.seq,ci.lines[[2]]$y)$y)

## now let's look at the initial state of susceptibles
## first for beta
S0.range <- seq(0.89,
                0.93,length=15)

mf.zim.S0.prof <- foreach(i=1:length(S0.range),
                          .inorder=FALSE,
                          .options.multicore=list(set.seed=TRUE)
                          ) %dopar%
{

    theta.guess <- coef(better2.mif.zim)
    theta.guess['S.0'] <- S0.range[i]
    ## need to reduce I.0 so init cond. sum to 1
    theta.guess['R.0'] <- max(0,1- sum(theta.guess[c('E.0','S.0','A.0','I.0')]))

    m1 <- mif(
        better2.mif.zim,
        Nmif=50,
        start=theta.guess,
        pars=c('beta1','gamma'),
        ivps=c('I.0'),
        transform=TRUE,
        rw.sd=c(beta1=0.1,gamma=0.1,I.0=0.1),
        Np=2000,
        ic.lag=length(zim.mod@data),
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
        )
    ll <- replicate(n=3,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

S0.logliks <- colMeans(sapply(mf.zim.S0.prof,function(x) x[[2]]))
cis <- max(S0.logliks) - qchisq(.95,1)/2
pdf("Plots/proflik-S0-zim.pdf")
plot(S0.range,S0.logliks,xlab="S0",ylab="log-likelihood",main="Profile Likelihood of Transmission Parameter (S0)")
abline(h=cis,lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(S0.logliks[5:7],S0.range[5:7],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(S0.logliks[11:13],S0.range[11:13],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
text(3.45,-338,sprintf("95%% CI %.2f-%.2f",
                       approx(S0.logliks[2:6],S0.range[2:6],xout=cis)$y,
                       approx(S0.logliks[7:11],S0.range[7:11],xout=cis)$y))
dev.off()


##################################################################
## For Harare only using data extracted from Fernandez et al.   ##
## Ultiumatley decided not to use this since the data are a bit ##
## suspect and pretty challenging to fit well with an SIR model ##
##################################################################

pop.harare <- 1606000
har.dat <- get.zim.data(harare.only = TRUE)
har.dat$week <- 1:nrow(har.dat)
## build pomp model object
har.mod <- build.leaky.model.C(pop=pop.harare,
                               dat=har.dat,
                               model.name = "harmod")

## specify starting parameters
## remember these are in units of weeks
E0 <- 10/0.04/pop.harare/3
I0 <- 10/0.04/pop.harare/3
A0 <- 100/pop.harare
R0 <- 0.5
S0 <- 1- R0-I0-E0-A0
guess.params <- c(gamma=7/3,
                  sigma=5,
                  theta=10,
                  beta1=3.9,
                  beta2=3.0,
                  rho=0.04,
                  theta0=0.0001,
                  S.0=S0,
                  E.0=E0,
                  I.0=I0,
                  A.0=A0,
                  R.0=R0)


#har.mod.win <- window(har.mod,start=16,end=44)
#har.mod.win@t0 <- 15

## first let's start with Trajector matching
tm.har <- traj.match(har.mod,
                     start=guess.params,
                     est=c('beta1','I.0','E.0','S.0'),
                     method="subplex",
                     maxit=15000,
                     transform=TRUE
                 )

sim.har.tm <- simulate(tm.har,
                       nsim=500,
                       seed=1914679109L,
                       transform=TRUE)

#pdf("Plots/harare-fernandezextract.pdf")
plot(har.dat[,2],ylim=c(0,2000),type="h")
#dev.off()
for (i in 1:500) {
    lines(sim.har.tm[[i]]@data[1,],lty=2,col=AddAlpha(4,.05))
}

mif.har <- mif(tm.har,
               Nmif=100,
               pars=c('beta1'),
               ivps=c('I.0','E.0','S.0'),
               transform=TRUE,
               rw.sd=c(beta1=0.1,S.0=0.1,E.0=0.1,I.0=0.1),
               Np=2000,
               ic.lag=length(har.mod@data),
               var.factor=1,
               cooling.type="hyperbolic",
               cooling.fraction=0.05,
               method="mif2",
               verbose=TRUE)

#mif.zim.mif1 <- mif(mif.zim.mif1,Nmif=50,cooling.fraction=0.80)
#mif3.zim <- mif(mif2.zim,Nmif=50,cooling.fraction=0.80)

sim.zim.mif <- simulate(mif.zim,
                       nsim=500,
                       seed=1914679109L,
                       transform=TRUE)

#dev.off()
for (i in 1:500) {
    lines(sim.zim.mif[[i]]@data[1,],lty=2,col=AddAlpha(2,.05))
}


#####################
## now for conakry ##
#####################

## 2010 populatino estimates from Institut National de la Statistique de Guinée
pop.con <- 1656300
## zimdat
con.dat <- get.conakry.data()

## build pomp model object
con.mod <- build.leaky.model.C(pop=pop.con,
                               dat=con.dat,
                               my.times="day",
                               my.t0=0,
                               model.name="conakrymodel")


## specify starting parameters
## remember these are in units of weeks
E0 <- 1/0.04/pop.con
I0 <- 1/0.04/pop.con
A0 <- 1/pop.con
R0 <- 0.5
S0 <- 1- R0-I0-E0-A0
guess.params.con <- c(gamma=7/3,
                      sigma=5,
                      theta=20,
                      beta1=3.9,
                      beta2=3.0,
                      rho=0.04,
                      theta0=0.0001,
                      S.0=S0,
                      E.0=E0,
                      I.0=I0,
                      A.0=A0,
                      R.0=R0)


# start and stops refer to indices not days (t0=0)
con.mod.win <- window(con.mod,start=33,end=153)
con.mod.win@t0 <- 32

tm.con <- traj.match(con.mod.win,
                     start=guess.params.con,
                     est=c('beta1','E.0','S.0','I.0','gamma'),
                     method="subplex",
                     maxit=15000,
                     transform=TRUE
                 )

sim.con.tm <- simulate(tm.con,
                       nsim=500,
                       seed=1914679109L,
                       transform=TRUE)

plot(con.dat[33:153,2],ylim=c(0,500))

for (i in 1:500) {
    lines(sim.con.tm[[i]]@data[1,],lty=2,col=AddAlpha(5,.05))
}
## not a great fit but it gets up somewhere
logLik(pfilter(tm.con,Np=1000))

mif.con <- mif(tm.con,
               start=coef(tm.con),
               Nmif=100,
               pars=c('beta1','gamma'),
               ivps=c('I.0','S.0'),
               transform=TRUE,
               rw.sd=c(beta1=0.15,gamma=0.15,I.0=0.15,S.0=0.15),
               Np=3000,
               ic.lag=length(con.mod.win@data),
               var.factor=1,
               cooling.type="hyperbolic",
               cooling.fraction=0.05,
               method="mif2",
               verbose=TRUE)

## run a particle filter
pf.con <- pfilter(mif.con,Np=5000,save.states=TRUE)

## get the average state at each time
#est.states.con <- sapply(pf.con@saved.states,rowMeans)
#saveRDS(est.states.con,file="GeneratedData/mif-con-states.rds")

sim.con.mif <- simulate(mif.con,
                        nsim=500,
                        seed=1914679109L,
                        transform=TRUE)

estpars <- c("beta1","gamma")
mf.con <- foreach(i=1:10,
                  .inorder=FALSE,
                  .options.multicore=list(set.seed=TRUE)
                  ) %dopar%
{

    theta.guess <- coef(mif.con)
    theta.guess[estpars] <- rlnorm(
        n=length(estpars),
        meanlog=log(theta.guess[estpars]),
        sdlog=0.1
        )

    ## now sample from I.0
    I.0.count <- runif(1,1,1e4)/pop.con # people
    theta.guess['S.0'] <- theta.guess['S.0'] - I.0.count
    theta.guess['I.0'] <- I.0.count
    theta.guess['S.0'] <-  theta.guess['S.0']*runif(1,.8,1.2)
    theta.guess['R.0'] <- max(0,1-sum(theta.guess[c('S.0','I.0','E.0','A.0')]))

    m1 <- mif(
        tm.con,
        Nmif=100,
        start=theta.guess,
        pars=c('beta1','gamma'),
        ivps=c('I.0','S.0'),
        transform=TRUE,
        rw.sd=c(beta1=0.1,gamma=0.15,I.0=0.15,S.0=0.15),
        Np=2000,
        ic.lag=length(con.mod@data),
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
        )
    ll <- replicate(n=10,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

saveRDS(mf.con,file="GeneratedData/parallel-mif-con-REV.rds")
#mf.con <- readRDS("GeneratedData/parallel-mif-con.rds")
mif.cons <- sapply(mf.con,function(v) v[[1]])
mif.cons.ll <- sapply(mf.con,function(v) v[[2]])
compare.mif(mif.cons)
best.pomps <- order(colMeans(mif.cons.ll),decreasing=TRUE)

better.mif.con <- mif(mf.con[[best.pomps[1]]][[1]],
                      Nmif=200,
                      pars=c('beta1','gamma'),
                      ivps=c('I.0','S.0'),
                      transform=TRUE,
                      rw.sd=c(beta1=0.15,I.0=0.15,S.0=0.15,gamma=0.15),
                      Np=4000,
                      ic.lag=length(con.mod.win@data),
                      var.factor=1,
                      cooling.type="hyperbolic",
                      cooling.fraction=0.05,
                      method="mif2")

logLik(pfilter(better.mif.con,Np=10000))

better.mif.con2 <- mif(better.mif.con,
                       Nmif=100,
                       pars=c('beta1','gamma'),
                       ivps=c('I.0','S.0'),
                       transform=TRUE,
                       rw.sd=c(beta1=0.15,I.0=0.15,S.0=0.15,gamma=0.15),
                       Np=4000,
                       ic.lag=length(con.mod.win@data),
                       var.factor=1,
                       cooling.type="hyperbolic",
                       cooling.fraction=0.05,
                       method="mif2")

logLik(pfilter(better.mif.con2,Np=10000))
## now let's save! this is going to be the object we use in
## future simulations!
#saveRDS(better.mif.con2,file="GeneratedData/mif-con.rds")
better.mif.con2 <- readRDS("GeneratedData/mif-con.rds")

## run a particle filter
pf.con <- pfilter(better.mif.con2,Np=10000,save.states=TRUE)

## get the average state at each time and save it for vac simulations
est.states.con <- sapply(pf.con@saved.states,rowMeans)
saveRDS(est.states.con,file="GeneratedData/mif-con-states.rds")


sim.con.mif <- simulate(better.mif.con2,
                        nsim=500,
                        seed=1914679109L,
                        transform=TRUE)

plot(con.dat[,2],ylim=c(0,250),xlab="epidemic day",ylab="cases per day",pch=4)
for (i in 1:500) {
    lines(33:153,sim.con.mif[[i]]@data[1,],lty=2,col=AddAlpha(3,.05))
}


pdf("Plots/mif-con-unvac.pdf")
con.mat <- sapply(sim.con.mif,function(x) x@data[1,])
con.means <- apply(con.mat,1,mean)
con.ci <- apply(con.mat,1,function(x) quantile(x,c(.025,.975)))
plot(con.dat[,2],ylim=c(0,250),xlab="epidemic day",ylab="cases per day",pch=4)
for (i in 1:500) {
    lines(33:153,sim.con.mif[[i]]@data[1,],lty=2,col=AddAlpha(4,.05))
}
lines(33:153,con.means,col=4)
lines(33:153,con.ci[1,],col=4,lty=2)
lines(33:153,con.ci[2,],col=4,lty=2)
legend("topright",c("simulated epidemic",
                    "mean simulated epidemic",
                    "95% Prediction Interval",
                    "data"),
       col=c(AddAlpha(4,0.1),4,4,"black"),lty=c(1,1,2,-1),pch=c(-1,-1,-1,4),bty="n")
dev.off()


pdf("Plots/hist-finalsize-uncon-con.pdf")
hist(colSums(sapply(sim.con.mif,function(x) x@data[1,])),
     col="grey",border="white",breaks="fd",
     xlab="Final Epidemic Size of Simulation",
     main="Final Size of Conakry Simulations")
abline(v=4566,col="orange",lwd=2,lty=2)
text(4750,50,"Reported Epidemic \n Size = 4,566",cex=.9)
dev.off()


## --------------------------------------------------------------- ##
## Let's do a little profiling of our beta and gamma parameters to ##
## see how peaky they look                                         ##
## --------------------------------------------------------------- ##

## first for beta
beta.range <- seq(coef(better.mif.con2)['beta1']*.6,coef(better.mif.con2)['beta1']*1.4,length=15)

mf.con.beta.prof <- foreach(i=1:length(beta.range),
                            .inorder=FALSE,
                            .options.multicore=list(set.seed=TRUE)
                            ) %dopar%
{

    theta.guess <- coef(better.mif.con2)
    theta.guess['beta1'] <- beta.range[i]

    m1 <- mif(
        better.mif.con2,
        Nmif=70,
        start=theta.guess,
        pars=c('gamma'),
        ivps=c('I.0','S.0'),
        transform=TRUE,
        rw.sd=c(gamma=0.15,I.0=0.15,S.0=0.15),
        Np=2000,
        ic.lag=length(con.mod@data),
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
        )
    ll <- replicate(n=3,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

## a few clearly didn't coverge
theta.guess <- coef(better.mif.con2)
theta.guess['beta1'] <- beta.range[11]

redo.m1 <- mif(
    better.mif.con2,
    Nmif=120,
    start=theta.guess,
    pars=c('gamma'),
    ivps=c('I.0','S.0'),
    transform=TRUE,
    rw.sd=c(gamma=0.15,I.0=0.15,S.0=0.15),
    Np=3000,
    ic.lag=length(con.mod.win@data),
    var.factor=1,
    cooling.type="hyperbolic",
    cooling.fraction=0.05,
    method="mif2"
    )
ll.redo1 <- replicate(n=3,logLik(pfilter(redo.m1,Np=10000)))


beta.logliks <- colMeans(sapply(mf.con.beta.prof,function(x) x[[2]]))
## swap out the mean of the redo
beta.logliks[11] <- mean(ll.redo1)
cis <- max(beta.logliks) - qchisq(.95,1)/2
pdf("Plots/proflik-beta-con.pdf")
plot(beta.range,beta.logliks,
     ylim=c(-500,-400),
     xlim=c(beta.range[1],beta.range[11]),
     xlab="beta",
     ylab="log-likelihood",
     main="Profile Likelihood of Transmission Parameter (Beta)")
abline(h=cis,lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(beta.logliks[5:7],beta.range[5:7],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(beta.logliks[8:10],beta.range[8:10],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
text(2.25,-418,sprintf("95%% CI %.2f-%.2f",approx(beta.logliks[5:7],beta.range[5:7],xout=cis)$y,
                      approx(beta.logliks[8:10],beta.range[8:10],xout=cis)$y))
dev.off()
## approximate 95% CI

gamma.range <- seq(0.58,0.7,length=20)

mf.con.gamma.prof <- foreach(i=1:length(gamma.range),
                             .inorder=FALSE,
                             .options.multicore=list(set.seed=TRUE)
                             ) %dopar%
{

    theta.guess <- coef(better.mif.con2)
    theta.guess['gamma'] <- gamma.range[i]

    m1 <- mif(
        better.mif.con2,
        Nmif=200,
        start=theta.guess,
        pars=c('beta1'),
        ivps=c('I.0','S.0'),
        transform=TRUE,
        rw.sd=c(beta1=0.15,I.0=0.15,S.0=0.15),
        Np=10000,
        ic.lag=length(con.mod@data),
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
        )
    ll <- replicate(n=3,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

plot(gamma.range,colMeans(sapply(mf.con.gamma.prof,function(x) x[[2]])))

gamma.logliks <- colMeans(sapply(mf.con.gamma.prof,function(x) x[[2]]))
plot(gamma.range,gamma.logliks,ylim=c(-425,-410))
cis <- max(gamma.logliks) -qchisq(.95,1)/2
pdf("Plots/proflik-gamma-con.pdf")
plot(gamma.range,gamma.logliks,xlab="gamma",ylab="log-likelihood",main="Profile Likelihood of Gamma")
abline(h=cis,lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(gamma.logliks[5:7],gamma.range[5:7],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(gamma.logliks[14:17],gamma.range[14:17],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
text(2,-337.5,sprintf("95%% CI %.2f-%.2f",approx(gamma.logliks[5:7],gamma.range[5:7],xout=cis)$y,
                      approx(gamma.logliks[14:17],gamma.range[14:17],xout=cis)$y))
dev.off()


beta.seq <- seq(2.5,6,length=30)
gamma.seq <- seq(2,5,length=30)
r.seq <- expand.grid(beta.seq,gamma.seq)

R.prof <- foreach(i=1:nrow(r.seq),
                             .inorder=FALSE,
                             .options.multicore=list(set.seed=TRUE)
                             ) %dopar%
{

    theta.guess <- coef(better.mif.zim)
    theta.guess['gamma'] <- r.seq[i,2]
    theta.guess['beta1'] <- r.seq[i,1]

    m1 <- mif(
        better.mif.zim,
        Nmif=50,
        start=theta.guess,
        ivps=c('I.0','S.0'),
        transform=TRUE,
        rw.sd=c(I.0=0.1,S.0=0.1),
        Np=2000,
        ic.lag=length(zim.mod@data)/2,
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.05,
        method="mif2"
    )
    ll <- replicate(n=3,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

saveRDS(R.logliks,file="GeneratedData/R_proflik.rds")

## get CI for R
prof.lik <- R.logliks
prof.mat <- matrix(prof.lik[,3],nrow=30)
colnames(prof.mat) <- beta.seq
rownames(prof.mat) <- gamma.seq

ci.lines <- contourLines(prof.mat*2,levels=max(prof.mat)*2 - 3.814/2)
## we will only take the middle
range(approx(seq(0,1,length=30),beta.seq,ci.lines[[2]]$x)$y/approx(seq(0,1,length=30),gamma.seq,ci.lines[[2]]$y)$y)



## -------------------------------------- ##
## Now Port au Prince                     ##
## we will start with the first wave only ##
## -------------------------------------- ##

source("Source/R/leakyvac-pomp-model-inC-novac-seasonal.R")

## need to get a better population estimate
## this is likley to be tricky given the IDP
## population at the time
pop.portap <- 2.1e6
## zimdat
portap.dat <- get.haiti.data(first.wave.only=F)

#covartab <- make.covartab(0,nrow(portap.dat)+1,byt=1,degree=3,nbasis=4)
#covartab <- make.covartab(0,nrow(portap.dat)+1,byt=1,degree=5,nbasis=5)
#covartab <- make.covartab(0,nrow(portap.dat)+1,byt=1,degree=4,nbasis=4)
covartab <- make.covartab(0,nrow(portap.dat)+1,byt=1,degree=6,nbasis=6)

## build pomp model object
portap.mod <- build.leaky.model.C.seas(pop=pop.portap,
                                       dat=portap.dat,
                                       my.times="day",
                                       my.t0=0,
                                       covar=covartab,
                                       model.name="papmodel")

#portap.mod.win <- window(portap.mod,start=1,end=297)

## specify starting parameters
## remember these are in units of weeks
E0 <- 10/pop.portap
I0 <- 10/pop.portap
A0 <- 0.0/pop.portap
R0 <- 0.000
S0 <- 1- R0-I0-E0-A0
guess.params.portap <- c(gamma=1/2,
                         sigma=1/1.4,
                         theta=10,
                         beta1=1.1,
                         beta2=.05,
                         beta3=.5,
                         beta4=.2,
                         beta5=.1,
                         beta6=1,
                         iota=1e-10,
                         rho=0.9,#.15
                         theta0=0.0,
                         S.0=S0,
                         E.0=E0,
                         I.0=I0,
                         A.0=A0,
                         R.0=R0)


tm.portap <- traj.match(portap.mod,
                        start=coef(mif.portap.best),#guess.params.portap,
                        est=c('beta1',
                            'beta2',
                            'beta3',
                            'beta4',
                            'beta5',
                            'beta6',
                            'rho',
                            'iota',
                            'I.0',
                            'E.0'),
                        method="Nelder-Mead",
                        maxit=15000,
                        transform=TRUE
                        )
summary(tm.portap)
logLik(pfilter(portap.mod,params=coef(tm.portap),Np=10000))

sim.portap.tm <- simulate(tm.portap,
                          params=coef(tm.portap),
                          nsim=500,
                          seed=1914679109L,
                          transform=TRUE)

plot(portap.dat[,2])

for (i in 1:200) {
    lines(sim.portap.tm[[i]]@data[1,],lty=2,col=AddAlpha(3,.05))
}

mif.portap.6df <- mif(tm.portap,
                      start=coef(mif.portap.best),
                      Nmif=100,
                      ivps = c('E.0','I.0'),
                      transform=TRUE,
                      rw.sd=c(
                          beta1=0.1,
                          beta2=0.1,
                          beta3=0.1,
                          beta4=0.1,
                          beta5=0.1,
                          beta6=0.1,
                          iota=0.1,
                          rho=0.1,
                          E.0=.12,
                          I.0=0.12),
                      Np=5000,
                      ic.lag=length(portap.mod@data),
                      var.factor=1,
                      cooling.type="hyperbolic",
                      cooling.fraction=0.03,
                      method="mif2",
                      verbose=FALSE)

mif.portap.6df.cont <- continue(mif.portap.6df,Nmif=50)
logLik(pfilter(mif.portap.6df.cont,Np=10000))

mif.portap.6df.next2 <- mif(mif.portap.6df.cont,
                           Nmif=50,
                           ivps = c('E.0','I.0'),
                           transform=TRUE,
                           rw.sd=c(
                               beta1=0.1,
                               beta2=0.1,
                               beta3=0.1,
                               beta4=0.1,
                               beta5=0.1,
                               beta6=0.1,
                               iota=0.1,
                               rho=0.1,
                               E.0=.12,
                               I.0=0.12,
                               R.0=0.12),
                           Np=5000,
                           ic.lag=length(portap.mod@data),
                           var.factor=1,
                           cooling.type="hyperbolic",
                           cooling.fraction=0.05,
                           method="mif2",
                           verbose=FALSE)

mif.portap.6df.next3 <- continue(mif.portap.6df.next2,Nmif=50)
mif.portap.6df.next4 <- continue(mif.portap.6df.next3,Nmif=50)

logLik( pfilter(mif.portap.6df,Np=10000))
logLik( pfilter(mif.portap.6df.next4,Np=10000))

#mif.portap.5df.2 <- continue(mif.portap.5df.2,Nmif=50)
#mif.portap <- mif(mif.portap2,Nmif=50)

## run a particle filter
pf.portap <- pfilter(mif.portap.6df,Np=5000,save.states=TRUE)
logLik(pf.portap)

saveRDS(mif.portap.6df.next4,file="GeneratedData/mif-haiti.rds")
mif.portap.6df.next4 <- readRDS("GeneratedData/mif-haiti.rds")
sim.mif.portap <- simulate(mif.portap.6df.next4,
#                           params=tmp,
                           nsim=500,
                           transform=TRUE)

plot(portap.dat[,2],ylim=c(0,2000))

for (i in 1:500) {
    lines(sim.mif.portap[[i]]@data[1,],lty=2,col=AddAlpha(3,.05))
}

par(new=T)
plot((covartab[,2]*coef(mif.portap.6df.cont)["beta1"] +
       covartab[,3]*coef(mif.portap.6df.cont)["beta2"] +
       covartab[,4]*coef(mif.portap.6df.cont)["beta3"] +
       covartab[,5]*coef(mif.portap.6df.cont)["beta4"] +
       covartab[,6]*coef(mif.portap.6df.cont)["beta5"] +
       covartab[,7]*coef(mif.portap.6df.cont)["beta6"])[-c(1:2)]/coef(mif.portap.6df.cont)["gamma"],
      ylab="",type="l",col="red")


estpars <- c("beta1","beta2","beta3","beta4","beta5","beta6","rho","iota")
mf.portap <- foreach(i=1:10,
                  .inorder=FALSE,
                  .options.multicore=list(set.seed=TRUE)
                  ) %dopar%
{

    theta.guess <- coef(mif.portap.6df.cont)
    theta.guess[estpars] <- rlnorm(
        n=length(estpars),
        meanlog=log(theta.guess[estpars]),
        sdlog=0.1
        )

    ## now sample from I.0
    I.0.count <- runif(1,1,100)/pop.portap # people
    E.0.count <- runif(1,1,100)/pop.portap # people
    theta.guess['E.0'] <- E.0.count
    theta.guess['I.0'] <- I.0.count
    theta.guess['S.0'] <- theta.guess['S.0'] - I.0.count - E.0.count
    theta.guess['R.0'] <- max(0,1-sum(theta.guess[c('S.0','I.0','E.0','A.0')]))

    m1 <- mif(
        tm.portap,
        Nmif=100,
        start=theta.guess,
        ivps=c('I.0','E.0'),
        transform=TRUE,
        rw.sd=c(
            beta1=0.1,
            beta2=0.1,
            beta3=0.1,
            beta4=0.1,
            beta5=0.1,
            beta6=0.1,
            iota=0.1,
            rho=0.1,
            I.0=0.1,
            E.0=0.1,
            R.0=0.1),
        Np=5000,
        ic.lag=length(portap.mod@data),
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.03,
        method="mif2"
        )
    ll <- replicate(n=10,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
}

## look at logliks
which.max(
    colMeans(sapply(mf.portap,function(x) x[[2]]))
    )

pf.best <- pfilter(mf.portap[[4]][[1]],Np=20000,save.states=TRUE)
logLik(pf.best)
sim.portap.mif <- simulate(mf.portap[[4]][[1]],
                           nsim=500,
                           seed=1914679109L,
                           transform=TRUE)

mif.portap.best <- mif(test,
                       Nmif=50,
                       ivps = c('E.0','I.0'),
                       transform=TRUE,
                       rw.sd=c(
                           beta1=0.1,
                           beta2=0.1,
                           beta3=0.1,
                           beta4=0.1,
                           beta5=0.1,
                           beta6=0.1,
                           rho=0.1,
                           iota=0.1,
                           E.0=0.12,
                           I.0=0.12),
                       Np=6000,
                       ic.lag=length(portap.mod@data),
                       var.factor=1,
                       cooling.type="hyperbolic",
                       cooling.fraction=0.01,
                       method="mif2",
                       verbose=FALSE)

mif.portap.best <- mif.portap.6df.next4
pf.best <- pfilter(mif.portap.best,Np=20000,save.states=TRUE)
logLik(pf.best)
est.states.portaup <- sapply(pf.best@saved.states,rowMeans)
saveRDS(est.states.portaup,file="GeneratedData/mif-haiti-states.rds")

## ---------------------- ##
## Bring in parallel runs ##
## ---------------------- ##
## mf.portap <- readRDS(file="GeneratedData/parallel-mif-portaup-6df-fitrho.rds")

## mif.portap.best <- mf.portap[[5]][[1]]
##saveRDS(mif.portap.best,file="GeneratedData/mif-haiti-REV.rds")
mf.portap.best <- readRDS(file="GeneratedData/mif-haiti.rds")

sim.portap.mif <- simulate(mif.portap.best,
                           #params=coef(mif.portap.best),
                           nsim=500,
                           seed=1914679109L,
                           transform=TRUE)

pdf("Plots/mif-pap-unvac-6df-best-seas-R0.pdf")
pap.mat <- sapply(sim.portap.mif,function(x) x@data[1,])
pap.means <- apply(pap.mat,1,mean)
pap.ci <- apply(pap.mat,1,function(x) quantile(x,c(.025,.975)))

plot(portap.dat[,2],ylim=c(0,2200),xlab="epidemic day",ylab="cases per day",col=4,pch=4)
for (i in 1:300) {
    lines(sim.portap.mif[[i]]@data[1,],lty=2,col=AddAlpha(4,.02))
}
lines(pap.means,col=4,lwd=2)
lines(pap.ci[1,],col=4,lty=2)
lines(pap.ci[2,],col=4,lty=2)
legend("topright",c("simulated epidemic",
                    "mean simulated epidemic",
                    "95% prediction interval",
                    "seasonal forcing function",
                    "data"),
       col=c(AddAlpha(4,0.1),4,4,3,"black"),lty=c(1,1,2,4,-1),pch=c(-1,-1,-1,-1,4),bty="n")
#dev.off()
par(new=T)
plot(
   # pf.best@states["S",]/colSums(pf.best@states[1:5,])*
     ((covartab[,2]*coef(mif.portap.best)["beta1"] +
       covartab[,3]*coef(mif.portap.best)["beta2"] +
       covartab[,4]*coef(mif.portap.best)["beta3"] +
       covartab[,5]*coef(mif.portap.best)["beta4"] +
       covartab[,6]*coef(mif.portap.best)["beta5"] +
       covartab[,7]*coef(mif.portap.best)["beta6"]
       ))
     [-c(1:2)]/coef(mif.portap.best)["gamma"]
     ,ylab="",axes=F,xlab="",type="l",col=3,lty=4)
axis(4)
dev.off()

pdf("Plots/hist-finalsize-uncon-pap-6df-full.pdf")
hist(colSums(sapply(sim.portap.mif,function(x) x@data[1,1:297])),
     col="grey",border="white",breaks="fd",
     xlab="Final Epidemic Size of Simulation",
     main="Final Size of Port au Prince Simulations")
abline(v=sum(mif.portap.best@data),col="orange",lwd=2,lty=2)
text(129000,70,"Reported Epidemic \n Size  = 119,902",cex=.9)
dev.off()

## compare.mif(
##     sapply(mf.portap,function(x) x[[1]])
##     )


## mif.portap.cont <- mif(mf.portap[[8]][[1]],Nmif=50)
## saveRDS(mf.portap,file="GeneratedData/parallel-mif-portaup-6df-fitrho.rds")
## tmp <- readRDS("GeneratedData/parallel-mif-portaup.rds")


## ## le
#t's do some profileing of params to get CIs
rho.range <- seq(0.9,.95,length=50)

mf.pap.rho.prof <- foreach(i=1:length(rho.range),
                           .inorder=FALSE,
                           .options.multicore=list(set.seed=TRUE)
                           ) %dopar%
{

    theta.guess <- coef(mf.portap.best)
    theta.guess['rho'] <- rho.range[i]


    m1 <- mif(
        mf.portap.best,
        Nmif=50,
        start=theta.guess,
        transform=TRUE,
        pars=c("beta1","beta2","beta3","beta4","beta5","beta6","iota"),
        rw.sd=c(
            beta1=0.1,
            beta2=0.1,
            beta3=0.1,
            beta4=0.1,
            beta5=0.1,
            beta6=0.1,
            iota=0.1),
        Np=5000,
        var.factor=1,
        cooling.type="hyperbolic",
        cooling.fraction=0.03,
        method="mif2",
        verbose=T
    )

    ll <- replicate(n=3,logLik(pfilter(m1,Np=20000)))
    list(mif=m1,ll=ll)
}

#saveRDS(mf.pap.rho.prof,file="GeneratedData/rho_proflik.rds")
mf.pap.rho.prof <- readRDS(file="GeneratedData/rho_proflik3.rds")
rho.logliks <- sapply(mf.pap.rho.prof,function(x) min(x[[2]]))
cis <-max(rho.logliks)- qchisq(.95,1)/2
plot(rho.range,sapply(mf.pap.rho.prof,function(x) max(x[[2]])),ylim=c(-2100,-2000))
points(rho.range,sapply(mf.pap.rho.prof,function(x) min(x[[2]])),pch=3)
,ylim=c(-2100,-2030))
pdf("Plots/proflik-gamma-zim-REV.pdf")
plot(gamma.range,gamma.logliks,xlab="gamma",ylab="log-likelihood",main="Profile Likelihood of Gamma (Zimbabwe)")
abline(h=cis,lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(gamma.logliks[5:7],gamma.range[5:7],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
abline(v=approx(gamma.logliks[20:25],gamma.range[20:25],xout=cis),lty=2,col=AddAlpha("orange",.75),lwd=2)
text(3.5,-337.5,sprintf("95%% CI %.2f-%.2f",
                        approx(gamma.logliks[5:7],gamma.range[5:7],xout=cis)$y,
                        approx(gamma.logliks[20:25],gamma.range[20:25],xout=cis)$y))
dev.off()

## exploring alternative fits for Haiti
haiti.para <- readRDS(file="GeneratedData/parallel-mif-portaup-6df-fitrho.rds")
logliks <- colMeans(sapply(haiti.para,function(x) x[[2]]))
order(logliks,decreasing = T)
