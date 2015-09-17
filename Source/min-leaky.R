## ------------------------------------------------------------------ ##
## Here we find the minimum single dose VE for leaky model
## by doses and time                                                  ##
## ------------------------------------------------------------------ ##

source("Source/base-functions.R")

## two dose eff
ve2 <- 0.80

## betas calibrated to 2009 GB epidemic
beta.sab <- 0.6538415
beta.seirb.sab <- 0.09424454
#beta.lusaka <- 0.7000951
#beta.seirb.lusaka <- 0.05305244

mydoses <- c(5000,355000)#seq(5000,355000,by=10000)
vac.starts <- seq(0,300,by=2)
min.ves <- min.aon <- min.seirb <- min.leaky <- array(dim=c(length(vac.starts),length(mydoses)))

params <- list(beta=beta.sab,
               gamma=1/2,
               sigma=1/1.4,
               phi1=0,
               phi2=0,
               theta1=0.325,
               theta2=ve2
               )

params.leaky <- list(beta=beta.sab,
                     gamma=1/2,
                     sigma=1/1.4,
                     phi1=0.0, # VE_I for first dose among fully symptomatic
                     phi2=0.0, # VE_I for second dose among fully symptomtic
                     kappa1=0.9, # VE_I for first dose among asym/mildly symptomatic
                     kappa2=0.9, # VE_I for second dose among asym/mildly symptomtic
                     theta0=0,
                     theta1=0.325,
                     theta2=ve2
                     )


params.seirb <- list(beta=beta.seirb.sab,
                     gamma=1/2,
                     sigma=1/1.4,
                     phi1=0,
                     phi2=0,
                     theta1=0.325,
                     theta2=ve2,
                     kappa = 10^6, #ID 50
                     epsilon=10, #virbrio extretion rate
                     delta=1/1.5
                     )

initial.state <- c(.5e6-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)
initial.state.leaky <- c(.5e6-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
initial.state.seirb <- c(.5e6-1,0,0,0,0,0,1,0,0,0,0,0,1000,0,0,0,0,0,0)

for(d in seq_along(mydoses)){
    #cat("|")
    for(t in seq_along(vac.starts)){


        ## min for both VES and VEI
        min.leaky[t,d] <-
            optim(.2,min.ve.1.obj.func,method="Brent",
                  lower=1e-2,upper=.5,
                  doses=mydoses[d],
                  start.time=vac.starts[t],
                  dx.dt.func=seir.leaky.dx.dt,
                  ve2=ve2,
                  params=params.leaky,
                  initial.state=initial.state.leaky,
                  model.type="leaky"
                  )$par


        cat("start=",vac.starts[t],", d=",mydoses[d],", minLeaky=",min.leaky[t,d],"\n")
    }
}

mins.betaSAB <- list(min.aon=min.aon,min.ves=min.ves,min.seirb=min.seirb,min.leaky=min.leaky)
save(mins.betaSAB,file=sprintf("GeneratedData/minVEleaky_betaSAB_VE20p%.0f-reduced.rda",ve2*100))

