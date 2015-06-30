########################################
## Model code for leaky vaccine model ##
########################################

##' Notes:
##' covert as much as possible to C to speed things up
##'

## -------------------------------- ##
## A few project specific functions ##
## -------------------------------- ##
## gets zimbabwe data
get.zim.data <- function(){
    zim <- read.csv("Data/zimepicurve-imputed.csv")
    zim[,1] <- as.Date(zim[,1],format="%m/%d/%y")
    zim[,1] <- as.numeric(zim[,1])
    zim[,2] <- round(zim[,2])
    zim$week <- 0:(nrow(zim)-1)
    #zim <- zim[-nrow(zim),] # for now getting rid of zero
    return(zim)
}

## gets conakry data
get.conakry.data <- function(conakry.only=TRUE){
    gc <- read.csv("Data/cholera\ en\ guinee\ 2012_preliminary.csv",as.is=T)
    ## limit to conakry only
    if(conakry.only)
        gc <- subset(gc,pref=="conakry")

    gc$datee <- as.Date(gc$datee,"%m/%d/%y")
    day.range <- data.frame(day=seq.Date(min(gc$datee),max(gc$datee),by=1))
    ec.con <- tapply(gc[,1],gc$datee,length)

    rc <- data.frame(day=as.Date(names(ec.con)),cases=ec.con)
    rc <- merge(rc,day.range,all.y=TRUE)
    rc <- ifelse(is.na(rc[,2]),0,rc[,2])
    rc <- data.frame(day=0:(nrow(day.range)-1),cases=as.numeric(rc))
    ## ec.full <- tapply(gc.full[,1],gc.full$datee,length)
    ## pdf("Plots/conakryvscountry.pdf")
    ## plot(as.Date(names(ec.full)),as.numeric(ec.full),type="h")
    ## points(as.Date(names(ec.con)),as.numeric(ec.con),type="h",col=2)
    ## dev.off()
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
  double rate[25]; // transition rates
  double trans[25]; // transition numbers
  int time_check_1; // 0 FALSE and 1 TRUE
  int time_check_2; // 0 FALSE and 1 TRUE
  double vac_rate_1=0.0;
  double vac_rate_2=0.0;
  double sec_vac = 0.0;
  double vac_starts_1=16.0;
  double vac_starts_2=18.71;
  double vac_ends_1=16.71;
  double vac_ends_2=19.43; // no second vac
  double daily_vac=4690000.0; // assuming 50% coverage note that this is weekly vac rate here

   // some population demonitors
  int N0 = S + E + I + A +  R;
  int N1 = S1 + E1 + I1 + A1 +  R1;
  int N2 = S2 + E2 + I2 + A2 +  R2;
  int N = N0 + N1 + N2;

  // time checks for vaccinating
   time_check_1 = 0;
   if (t >= vac_starts_1 && t <= vac_ends_1)
       time_check_1 = 1;

   time_check_2 = 0;
   if (t >= vac_starts_2 && t <= vac_ends_2)
       time_check_2 = 1;

   if (time_check_1 == 1)
       vac_rate_1 = daily_vac;

   if (time_check_2 == 1)
       vac_rate_2 = daily_vac;

   double prim_vac = vac_rate_1/N0;

   double foi =  beta1/N*(I+(1-phi1)*I1+(1-phi2)*I2 +
                   A*(1-kappa)+(1-kappa)*A1+(1-kappa)*A2);

   if(time_check_2 == 1)
      sec_vac = vac_rate_2/N1;

   //compute the rates for all the transitions
   rate[0]= foi;  //S -> E
   rate[1]= prim_vac;
   rate[2]= foi; //S1 new infections
   rate[3]= sec_vac; //S1 Vaccinated
   rate[4]= foi;
   rate[5]= sigma*(1-theta0);
   rate[6]= sigma*theta0;
   rate[7]= prim_vac;
   rate[8]= sigma*(1-theta1); //from E1
   rate[9]= sigma*theta1;
   rate[10]= sec_vac;
   rate[11]= sigma*(1-theta2); //from E2
   rate[12]= sigma*theta2;
   rate[13]= gamma; //from I
   rate[14]= prim_vac;
   rate[15]= gamma; //from I1
   rate[16]= sec_vac;
   rate[17]= gamma; // from I2
   rate[18]= gamma; //from A
   rate[19]= prim_vac;
   rate[20]= gamma; //A1
   rate[21]= sec_vac;
   rate[22]= gamma; //A2
   rate[23]= prim_vac; //R
   rate[24]= sec_vac; //R1 (no R2)

  // compute the transition numbers
  reulermultinom(2,S,&rate[0],dt,&trans[0]);
  reulermultinom(2,S1,&rate[2],dt,&trans[2]);
  reulermultinom(1,S2,&rate[4],dt,&trans[4]);
  reulermultinom(3,E,&rate[5],dt,&trans[5]);
  reulermultinom(3,E1,&rate[8],dt,&trans[8]);
  reulermultinom(2,E2,&rate[11],dt,&trans[11]);
  reulermultinom(2,I,&rate[13],dt,&trans[13]);
  reulermultinom(2,I1,&rate[15],dt,&trans[15]);
  reulermultinom(1,I2,&rate[17],dt,&trans[17]);
  reulermultinom(2,A,&rate[18],dt,&trans[18]);
  reulermultinom(2,A1,&rate[20],dt,&trans[20]);
  reulermultinom(1,A2,&rate[22],dt,&trans[22]);
  reulermultinom(1,R,&rate[23],dt,&trans[23]);
  reulermultinom(1,R1,&rate[24],dt,&trans[24]);

  // balance the equations
  S += -trans[0]-trans[1];
  S1 += -trans[2]-trans[3]+trans[1];
  S2 += -trans[4]+trans[3];
  E += trans[0]-trans[5]-trans[6]-trans[7];
  E1 += trans[2]-trans[8]-trans[9]+trans[7]-trans[10];
  E2 += trans[4]-trans[11]-trans[12]+trans[10];
  I += trans[5]-trans[13]-trans[14];
  I1 += trans[8]-trans[15]+trans[14]-trans[16];
  I2 += trans[11]-trans[17]+trans[16];
  A += trans[6]-trans[18]-trans[19];
  A1 += trans[9]-trans[20]+trans[19]-trans[21];
  A2 += trans[12]-trans[22]+trans[21];
  R += trans[13]+trans[18]-trans[23];
  R1 += trans[15]+trans[20]+trans[23]-trans[24];
  R2 += trans[17]+trans[22]+trans[24];
  incid += trans[5]+trans[8]+trans[11]; // incidence is cumulative recoveries
'




## ---------------------------------- ##
## Determinisitic Skeleton of process ##
## ---------------------------------- ##

seir.leaky.vac.skel.C <- '
  double rate[25]; // transition rates
  double term[25]; // term numbers
  int time_check_1; // 0 FALSE and 1 TRUE
  int time_check_2; // 0 FALSE and 1 TRUE
  double vac_rate_1=0.0;
  double vac_rate_2=0.0;
  double sec_vac = 0.0;
  double vac_starts_1=16.0;
  double vac_starts_2=18.71;
  double vac_ends_1=16.71;
  double vac_ends_2=19.43; // no second vac
  double daily_vac=4690000.0; // assuing 40% coverage note that this is weekly vac rate here

   // some population demonitors
  int N0 = S + E + I + A+ R;
  int N1 = S1 + E1 + I1 + A1+ R1;
  int N2 = S2 + E2 + I2 + A2+ R2;
  int N = N0 + N1 + N2;

  // time checks for vaccinating
   time_check_1 = 0;
   if (t >= vac_starts_1 && t <= vac_ends_1)
       time_check_1 = 1;

   time_check_2 = 0;
   if (t >= vac_starts_2 && t <= vac_ends_2)
       time_check_2 = 1;

   if (time_check_1 == 1)
       vac_rate_1 = daily_vac;

   if (time_check_2 == 1)
       vac_rate_2 = daily_vac;

   double prim_vac = vac_rate_1/N0;

   double foi =  beta1/N*(I+(1-phi1)*I1+(1-phi2)*I2 +
                   A*(1-kappa)+(1-kappa)*A1+(1-kappa)*A2);

   if(time_check_2 == 1)
      sec_vac = vac_rate_2/N1;


   // get all the rates of transition
   //compute the rates for all the transitions
   rate[0]= foi;  //S -> E
   rate[1]= prim_vac;
   rate[2]= foi; //S1 new infections
   rate[3]= sec_vac; //S1 Vaccinated
   rate[4]= foi; //S2
   rate[5]= sigma*(1-theta0); //E -> I
   rate[6]= sigma*theta0; //E -> A
   rate[7]= prim_vac;
   rate[8]= sigma*(1-theta1); //from E1
   rate[9]= sigma*theta1;
   rate[10]= sec_vac;
   rate[11]= sigma*(1-theta2); //from E2
   rate[12]= sigma*theta2;
   rate[13]= gamma; //from I
   rate[14]= prim_vac;
   rate[15]= gamma; //from I1
   rate[16]= sec_vac;
   rate[17]= gamma; // from I2
   rate[18]= gamma; //from A
   rate[19]= prim_vac;
   rate[20]= gamma; //A1
   rate[21]= sec_vac;
   rate[22]= gamma; //A2
   rate[23]= prim_vac; //R
   rate[24]= sec_vac; //R1 (no R2)

   // compute transistion terms
   // not totally sure why this is necceary to do sep. but
  // following vignette for now
  term[0] = rate[0]*S;
  term[1] = rate[1]*S;
  term[2] = rate[2]*S1;
  term[3] = rate[3]*S1;
  term[4] = rate[4]*S2;
  term[5] = rate[5]*E;
  term[6] = rate[6]*E;
  term[7] = rate[7]*E;
  term[8] = rate[8]*E1;
  term[9] = rate[9]*E1;
  term[10] = rate[10]*E1;
  term[11] = rate[11]*E2;
  term[12] = rate[12]*E2;
  term[13] = rate[13]*I;
  term[14] = rate[14]*I;
  term[15] = rate[15]*I1;
  term[16] = rate[16]*I1;
  term[17] = rate[17]*I2;
  term[18] = rate[18]*A;
  term[19] = rate[19]*A;
  term[20] = rate[20]*A1;
  term[21] = rate[21]*A1;
  term[22] = rate[22]*A2;
  term[23] = rate[23]*R;
  term[24] = rate[24]*R1;

  DS= -term[0]-term[1];
  DS1=-term[2]-term[3]+term[1];
  DS2=-term[4]+term[3];
  DE=term[0]-term[5]-term[6]-term[7];
  DE1=term[2]-term[8]-term[9]+term[7]-term[10];
  DE2=term[4]-term[11]-term[12]+term[10];
  DI=term[5]-term[13]-term[14];
  DI1=term[8]-term[15]+term[14]-term[16];
  DI2=term[11]-term[17]+term[16];
  DA=term[6]-term[18]-term[19];
  DA1=term[9]-term[20]+term[19]-term[21];
  DA2=term[12]-term[22]+term[21];
  DR=term[13]+term[18]-term[23];
  DR1=term[15]+term[20]+term[23]-term[24];
  DR2=term[17]+term[22]+term[24];
  Dincid = term[5]+term[8]+term[11];
'


## ---------------------------------------- ##
## Simulation model for measurement process ##
## ---------------------------------------- ##


seir.leaky.vac.meas.sim.C <- '
  cases = rnbinom_mu(theta,rho*incid);
'

## -------------------------------------- ##
## Measurement process density evaluation ##
## -------------------------------------- ##


seir.leaky.vac.meas.dens.C <- '
  lik = dnbinom_mu(cases,theta,rho*incid,give_log);
'

## ------------------------- ##
## Parameter Transformations ##
## ------------------------- ##
## from estimation scale to natural scale
par.trans.leaky.vac.C <- '
  Tbeta1=exp(beta1);
  Tsigma=exp(sigma);
  Tgamma=exp(gamma);
  Trho=expit(rho);
  Ttheta=exp(theta);
  Ttheta0=expit(theta0);
  Ttheta1=expit(theta1);
  Ttheta2=expit(theta2);
  Tphi1=expit(phi1);
  Tphi2=expit(phi2);
  Tkappa=expit(kappa);
  from_log_barycentric(&TS_0,&S_0,15);
'

## from estimation scale to natural scale
par.inv.trans.leaky.vac.C <- '
  Tbeta1=log(beta1);
  Tsigma=log(sigma);
  Tgamma=log(gamma);
  Trho=logit(rho);
  Ttheta=log(theta);
  Ttheta0=logit(theta0);
  Ttheta1=logit(theta1);
  Ttheta2=logit(theta2);
  Tphi1=logit(phi1);
  Tphi2=logit(phi2);
  Tkappa=logit(kappa);
  to_log_barycentric(&TS_0,&S_0,15);

'



## --------------------- ##
## Build pomp model for  ##
## --------------------- ##

## note that vaccination stuff is currently hard coded into simulater
## and skeleton
build.leaky.model.C <- function(pop=13.4e6,
                                dat=get.zim.data()[,2:3],
                                my.times="week",
                                my.t0=0,
                                model.name="leakymodWvac"
                                ){


    my.mod <- pompBuilder(
        name=model.name,
        data=dat,
        times=my.times,
        t0=my.t0,
        dmeasure=seir.leaky.vac.meas.dens.C,
        rmeasure=seir.leaky.vac.meas.sim.C,
        step.fn=seir.leaky.vac.step.C,
        step.fn.delta.t=1/10, # 1/20th of a week for zim
        skeleton.type="vectorfield",
        skeleton=seir.leaky.vac.skel.C,
        parameter.transform=par.trans.leaky.vac.C,
        parameter.inv.transform=par.inv.trans.leaky.vac.C,
        statenames=c(
            "S","E","I","A","R",
            "S1","E1","I1","A1","R1",
            "S2","E2","I2","A2","R2",
            "incid"),
        paramnames=c(
            "beta1","sigma","gamma","rho","theta","theta0",
            "theta1","theta2","phi1","phi2",
            "kappa",
            "S.0","E.0","I.0","A.0","R.0",
            "S1.0","E1.0","I1.0","A1.0","R1.0",
            "S2.0","E2.0","I2.0","A2.0","R2.0"),
        zeronames=c("incid"),
                initializer=function(params, t0, ...) {
            x0 <- setNames(numeric(16),c("S","S1","S2",
                                         "E","E1","E2",
                                         "I","I1","I2",
                                         "A","A1","A2",
                                         "R","R1","R2",
                                         "incid"))
            fracs <- params[c("S.0","S1.0","S2.0",
                              "E.0","E1.0","E2.0",
                              "I.0","I1.0","I2.0",
                              "A.0","A1.0","A2.0",
                              "R.0","R1.0","R2.0")]
            x0[-length(x0)] <- round(pop*fracs/sum(fracs))
            x0
        }
        ,save=TRUE
        )
    return(my.mod)
}

