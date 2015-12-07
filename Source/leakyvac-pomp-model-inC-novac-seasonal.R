##################################################
## Model code for leaky vaccine model           ##
##  without vaccine states for fitting purposes ##
##  with seasonal forcing                       ##
##################################################


##' Notes:
##' covert as much as possible to C to speed things up
##'

## -------------------------------- ##
## A few project specific functions ##
## -------------------------------- ##



## gets haiti data
get.haiti.data <- function(first.wave.only=T){
    haiti <- read.csv("Data/portauprince.csv",as.is=T)

    ## first two days have no cases so we will scrub them
    ## we will take the cases.whole column which is data extracted
    ## from mspp.gouv.ht/site/downloads/Rapport journalier MSPP du 30 septembre 2011.pdf
    ## the extraction is not 100% perfect as the final size is off by ~5 cases but this
    ## should be more than fine since we are dealing with many many suspected cases
    if (first.wave.only){
        ## get trough between peaks
        low <- which.min(haiti[100:200,5])
        haiti <- haiti[1:(100+low),]
    }

    ## remove zeros
    haiti <- haiti[-c(1:2),]
    rc <- data.frame(day=1:nrow(haiti),cases=as.numeric(haiti$cases.whole))
    return(rc)
}

## ---------------------------- ##
## Load required libraries etc. ##
## ---------------------------- ##
library(pomp)

## ---------------------------------------- ##
## Step function used in simluating process ##
## ---------------------------------------- ##

seir.leaky.vac.step.C <- '
  double mybeta;
  double rate[5]; // transition rates
  double trans[5]; // transition numbers

   // some population demonitors
  int N = S + E + I + A + R;

  // make seasonal beta term for current time
   mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 +
            beta4*seas4 + beta5*seas5 + beta6*seas6;
  //mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4;

  double foi =  iota + (I + A)*mybeta/N;

   //compute the rates for all the transitions
   rate[0]= foi;  //S -> E
   rate[1]= sigma*(1-theta0);
   rate[2]= sigma*theta0;
   rate[3]= gamma; //from I
   rate[4]= gamma; //from A

  // compute the transition numbers
  reulermultinom(1,S,&rate[0],dt,&trans[0]);
  reulermultinom(2,E,&rate[1],dt,&trans[1]);
  reulermultinom(1,I,&rate[3],dt,&trans[3]);
  reulermultinom(1,A,&rate[4],dt,&trans[4]);

  // balance the equations
  S += -trans[0];
  E += trans[0]-trans[1]-trans[2];
  I += trans[1]-trans[3];
  A += trans[2]-trans[4];
  R += trans[3]+trans[4];
  incid += trans[3]; // incidence is cumulative recoveries
'

## ---------------------------------- ##
## Determinisitic Skeleton of process ##
## ---------------------------------- ##

seir.leaky.vac.skel.C <- '
  double mybeta;
  double rate[5]; // transition rates
  double term[5]; // term numbers

   // some population demonitors
  int N = S + E + I + A+ R;

  mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6;
  //mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4;

 // double foi =  (I+A)*mybeta/N;
  double foi =  iota + (I + A)*mybeta/N;

   //compute the rates for all the transitions
   rate[0]= foi;  //S -> E
   rate[1]= sigma*(1-theta0);
   rate[2]= sigma*theta0;
   rate[3]= gamma; //from I
   rate[4]= gamma; //from A

   // compute transistion terms
   // not totally sure why this is necceary to do sep. but
  // following vignette for now
  term[0] = rate[0]*S;
  term[1] = rate[1]*E;
  term[2] = rate[2]*E;
  term[3] = rate[3]*I;
  term[4] = rate[4]*A;

  DS= -term[0];
  DE= term[0]-term[1]-term[2];
  DI= term[1]-term[3];
  DA= term[2]-term[4];
  DR= term[3]+term[4];
  Dincid = term[3];
'


## ---------------------------------------- ##
## Simulation model for measurement process ##
## ---------------------------------------- ##

seir.leaky.vac.meas.sim.C <- '
  cases = rnbinom_mu(theta,rho*incid);
  //cases = rpois(incid*rho);
'

## -------------------------------------- ##
## Measurement process density evaluation ##
## -------------------------------------- ##


seir.leaky.vac.meas.dens.C <- '
  lik = dnbinom_mu(cases,theta,rho*incid,give_log);
  //lik = dpois(cases, incid*rho, give_log);
'

## ------------------------- ##
## Parameter Transformations ##
## ------------------------- ##
## from estimation scale to natural scale
par.trans.leaky.vac.C <- '
  Tbeta1=exp(beta1);
  Tbeta2=exp(beta2);
  Tbeta3=exp(beta3);
  Tbeta4=exp(beta4);
  Tbeta5=exp(beta5);
  Tbeta6=exp(beta6);
  Tiota = exp(iota);
  Tsigma=exp(sigma);
  Tgamma=exp(gamma);
  Trho=expit(rho);
  Ttheta=exp(theta);
  Ttheta0=expit(theta0);
  from_log_barycentric(&TS_0,&S_0,5);
'

## from estimation scale to natural scale
par.inv.trans.leaky.vac.C <- '
  Tbeta1=log(beta1);
  Tbeta2=log(beta2);
  Tbeta3=log(beta3);
  Tbeta4=log(beta4);
  Tbeta5=log(beta5);
  Tbeta6=log(beta6);
  Tiota = log(iota);
  Tsigma=log(sigma);
  Tgamma=log(gamma);
  Trho=logit(rho);
  Ttheta=log(theta);
  Ttheta0=logit(theta0);
  to_log_barycentric(&TS_0,&S_0,5);
'

## --------------------------------- ##
## construct seasonal B-spline terms ##
## --------------------------------- ##


make.covartab <- function(t0,tmax,byt,nbasis=3,degree=3){

    tbasis <- seq(from=t0,to=tmax,by=byt)

    covartab <-  data.frame(cbind(time=tbasis,
                                  periodic.bspline.basis(
                                      x=tbasis,
                                      nbasis=nbasis,
                                      degree=degree,
                                      period=365.25,
                                      names="seas%d"
                                      )))

    return(covartab)
}

## --------------------- ##
## Build pomp model for  ##
## --------------------- ##

## note that vaccination stuff is currently hard coded into simulater
## and skeleton
##' @param pop - population
##' @param dat - data from get.zim.data() or get.conakry.data()
##' @param my.times - column name for time in dat
##' @param my.t0 - time zero for pomp model
##' @param model.name - name that the C files will take
##' @return pomp model object
##' @author asa
build.leaky.model.C.seas <- function(pop=2.1e6,
                                dat=get.haiti.data(FALSE),
                                my.times="day",
                                my.t0=-1,
                                model.name="leakymod",
                                covar
                                ){


    my.mod <- pompBuilder(
        name=model.name,
        data=dat,
        times=my.times,
        t0=my.t0,
        dmeasure=seir.leaky.vac.meas.dens.C,
        rmeasure=seir.leaky.vac.meas.sim.C,
        step.fn=seir.leaky.vac.step.C,
        step.fn.delta.t=1/2, # 1/2 a day of a week for haiti
        skeleton.type="vectorfield",
        skeleton=seir.leaky.vac.skel.C,
        covar=covar,
        tcovar="time",
        parameter.transform=par.trans.leaky.vac.C,
        parameter.inv.transform=par.inv.trans.leaky.vac.C,
        statenames=c("S","E","I","A","R","incid"),
        paramnames=c("beta1","beta2","beta3","beta4","beta5","beta6",
            "sigma",
            "gamma","rho","iota",
            "theta",
            "theta0",
            "S.0","E.0","I.0","A.0","R.0"),
        ivpnames=c("S.0","E.0","I.0","A.0","R.0"),
        zeronames=c("incid"),
        initializer=function(params, t0, ...) {
            x0 <- setNames(numeric(6),c("S",
                                        "E",
                                        "I",
                                        "A",
                                        "R",
                                        "incid"))
            fracs <- params[c("S.0",
                              "E.0",
                              "I.0",
                              "A.0",
                              "R.0")]
            x0[-length(x0)] <- round(pop*fracs/sum(fracs))
            x0
        }
        ,save=TRUE
        )
    return(my.mod)
}
