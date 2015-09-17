## -------------------------------------------------------------------------------------------------------------------- ##
## Runs across vaccine timeing and VE1 for VES Model                                                                    ##
## Here we calculate Incidence rates in both vaccinated and unvaccinated                                                ##
## and final size in addtion to the null epidemic                                                                    ## ##
## -------------------------------------------------------------------------------------------------------------------- ##
source("Source/base-functions.R")
library(parallel)
library(dplyr)
library(magrittr)

## key params
my.cores <- 2
ve2 <- 0.88
beta.sab <- 0.6538415

## creating loop variables for mapply statement
vacdoses <- 5e5#c(100000,200000,500000)#seq(5000,355000,by=10000)
ve1s <- seq(0,ve2,by=.02)
vac.timings <- seq(0,300,by=2)
min.inf.thresh <- 0.0001
loop.args <- expand.grid(ve1=ve1s,vactime=vac.timings,doses=vacdoses)

run.ves <- function(vactime,
                    doses,
                    ve1,
                    ve2,
                    my.beta=beta.sab,
                    inf.thresh=min.inf.thresh){

    time.step <- 1
    times <- seq(0,1500,by=time.step)

    seeded.infs <- 1
    initial.state <- c(0.5e6-seeded.infs,0,0,0,0,0,seeded.infs,0,0,0,0,0,0,0,0,0,0,0,0)

    params <- list(beta=my.beta,
                   gamma=1/2,
                   sigma=1/1.4,
                   phi1=0,
                   phi2=0,
                   theta1=ve1,
                   theta2=ve2
                   )

    uncon.run <- ode(y=initial.state,
                     times=times,
                     func=seir.ves.dx.dt,
                     parms=params,
                     vac.starts=c(1e6,1e6),
                     vac.ends=c(1e6,1e6),
                     daily.vac=0)

    uncon.run <- add.colnames(uncon.run)
    total.vac.days <- 10
    inter.dose.time <- 14 #14 #from last dose

    ## used to grab the relavent columns when figuring out the time at risk
    unvac.states <- "S"       # c("S","E","I","R")
    single.vac.states <- "S1" # c("S1","E1","I1","R1")
    double.vac.states <- "S2" # c("S2","E2","I2","R2")


    daily.vac.doses <- doses/total.vac.days

    one.of.one.start <- vactime
    one.of.one.end <- one.of.one.start+total.vac.days
    one.of.two.start <- vactime
    one.of.two.end <- one.of.two.start+total.vac.days/2
    two.of.two.start <- one.of.two.end+inter.dose.time
    two.of.two.end <- two.of.two.start + total.vac.days/2

    run.single.vac <- ode(y=initial.state,
                          times=times,
                          func=seir.ves.dx.dt,
                          parms=params,
                          vac.starts=c(one.of.one.start,Inf),
                          vac.ends=c(one.of.one.end,Inf),
                          daily.vac=daily.vac.doses)

    run.single.vac <- add.colnames(run.single.vac)

    ## ## figure out the time at which the epidemic is effectivley over
    t1 <- which(
        rowSums(run.single.vac[-c(1:50),c("I1","I2","I")]) < inf.thresh)[1] +
          51
    ## get the times at risk for unvac, 1 dose vac, and 2 dose vac
    time.at.risk.1dose <- c(
        time.step*sum(run.single.vac[,unvac.states]),
        time.step*sum(run.single.vac[,single.vac.states]),
        time.step*sum(run.single.vac[,double.vac.states]))

    ## get the cuulative incidence for each vaccination state

    ci.vac.1dose <- run.single.vac %>% data.frame %>% tail(1) %>% summarize(V1 - (S1 + S2)) %>% as.numeric
    fs.mat.1dose <-
        c(tail(run.single.vac[,c("CI","CI1","CI2")],1),ci.vac.1dose)

    run.double.vac <- ode(y=initial.state,
                          times=times,
                          func=seir.ves.dx.dt,
                          parms=params,
                          vac.starts=c(one.of.two.start,two.of.two.start),
                          vac.ends=c(one.of.two.end,two.of.two.end),
                          daily.vac=daily.vac.doses)

    run.double.vac <- add.colnames(run.double.vac)

    t2 <- which(
        rowSums(run.double.vac[-c(1:50),c("I1","I2","I")])
        < inf.thresh)[1]+51
    ## print(paste0("t2=",t2, ":: t1=",t1))
    ## get the times at risk for unvac, 1 dose vac, and 2 dose vac
    time.at.risk.2dose <- c(
        time.step*sum(run.double.vac[,unvac.states]),
        time.step*sum(run.double.vac[,single.vac.states]),
        time.step*sum(run.double.vac[,double.vac.states]))

    ## get the CIs for each vac state in 2 dose campaign
    ci.vac.2dose <- run.double.vac %>%
      data.frame %>%
      tail(1) %>%
      summarize(V1 - (S1 + S2)) %>%
      as.numeric

    fs.mat.2dose <-
        c(tail(run.double.vac[,c("CI","CI1","CI2")],1),ci.vac.2dose)

    rc <- list(fs.mat.1dose=fs.mat.1dose,
               fs.mat.2dose=fs.mat.2dose,
               time.at.risk.1dose = time.at.risk.1dose,
               time.at.risk.2dose = time.at.risk.2dose,
               doses=doses,
               times=times,
               uncon.run=uncon.run,
               ve1=ve1
               )

    return(rc)
}

## my.biglist <- mcmapply(run.ves,
##                        doses=loop.args$doses,
##                        vactime=loop.args$vactime,
##                        ve1=loop.args$ve1,
##                        MoreArgs = list(my.beta=beta.sab,ve2=ve2),
##                        mc.cores=my.cores)

my.biglist <- mapply(run.ves,
                     doses=loop.args$doses,
                     vactime=loop.args$vactime,
                     ve1=loop.args$ve1,
                     MoreArgs = list(my.beta=beta.sab,ve2=ve2))


save(my.biglist,file=sprintf("GeneratedData/vesruns_VE20p%.0f_betaSAB-12feb15.rda",ve2*100))
