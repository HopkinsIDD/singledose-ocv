## ------------------------------------------------------------------ ##
## Here we find the minimum single dose VE for a few different models ##
## by doses and time                                                  ##
## ------------------------------------------------------------------ ##
source("Source/R/base-functions.R")

## two dose eff (wtf)
ve2 <- 0.88

## betas calibrated to 2009 GB epidemic
beta.sab <- 0.6538415

mydoses <- c(5000,355000)#seq(5000,355000,by=10000)
vac.starts <- seq(0,300,by=1)
min.ves  <- array(dim=c(length(vac.starts),length(mydoses)))

params <- list(beta=beta.sab,
               gamma=1/2,
               sigma=1/1.4,
               phi1=0,
               phi2=0,
               theta1=0.325,
               theta2=ve2
               )

initial.state <- c(.5e6-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

for(d in seq_along(mydoses)){
                                        #cat("|")
    for(t in seq_along(vac.starts)){

        ## min VE for a model with vaccine reducing S to disease
        min.ves[t,d] <-
            optim(.2,min.ve.1.obj.func,method="Brent",
                  lower=1e-2,upper=.5,
                  doses=mydoses[d],
                  start.time=vac.starts[t],
                  dx.dt.func=seir.ves.dx.dt,
                  ve2=ve2,
                  params=params,
                  initial.state=initial.state
                  )$par

        cat("start=",vac.starts[t],", d=",mydoses[d],", minVES=",min.ves[t,d],"\n")
    }
}

mins.betaSAB <- list(min.ves=min.ves)
save(mins.betaSAB,file=sprintf("GeneratedData/minVES_betaSAB_VE20p%.0f-reduced.rda",ve2*100))
