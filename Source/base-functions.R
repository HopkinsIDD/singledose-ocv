
#######################################################################
##  Core functions for Performance of One- versus                    ##
##  Two-Dose Oral Cholera Vaccine Protocols in Response to Outbreaks ##
##  Code by: Andrew Azman                                            ##
##  Last Update: 5-July-2014                                         ##
##                                                                   ##
#######################################################################

## -------------------------- ##
## Section 1: Model Structure ##
## -------------------------- ##
##' Here we have the 'delta' functions for a bunch of differnt models
##' we use in the manuscript

require(deSolve)

##' delta function for a single population
##' all or nothing vaccination
##' @param t - time
##' @param x - state vector
##' @param params  - parameter vector
##' @param vac.starts - vector of vaccination starts, c(1st dose,2nd dose)
##' @param vac.ends - vector of vaccination ends, c(1st dose,2nd dose)
##' @param daily.vac - daily number of doses distributed
##' @return
##' @author Andrew Azman
seir.all.or.nothing.dx.dt <- function(t,
                                      x,
                                      params,
                                      vac.starts,
                                      vac.ends,
                                      daily.vac){

    S0  <- x[1]
    S1 <- x[2]
    S2 <- x[3]
    E0  <- x[4]
    E1  <- x[5]
    E2 <- x[6]
    I0  <- x[7]
    I1 <- x[8]
    I2 <- x[9]
    R0  <- x[10]
    R1 <- x[11]
    R2 <- x[12]
    V1 <- x[13]
    V2 <- x[14]
    W  <- x[15]

    ## note that Vs and W  are simply used for tracking and
    ## aren't real states
    N <- sum(x[1:12])
    N0 <- S0 + E0 + I0 + R0
    N1 <- S1 + E1 + I1 + R1
    N2 <- S2 + E2 + I2 + R2

    ## check that it is a time we want to vaccinate
    time.check <- t >= vac.starts & t <= vac.ends
    vaccinating.check <- any(time.check)

    ## based on the check assign a vaccination rate
    vac.rate <- sapply(1:length(time.check),
                       function(x) ifelse(time.check[x],daily.vac,0))

    ## doing this so we don't have Inf when denominator is 0 before
    ## there are any people with a second dose
    if(time.check[2]){
        sec.vac.S1 <- vac.rate[2]*S1/N1
        sec.vac.E1 <- vac.rate[2]*E1/N1
        sec.vac.I1 <- vac.rate[2]*I1/N1
        sec.vac.R1 <- vac.rate[2]*R1/N1
        wasted.sec <- vac.rate[2]*(1-(S1/N1))
    } else {
        sec.vac.S1<-sec.vac.E1<-sec.vac.I1<-sec.vac.R1<-0
        wasted.sec<-0
    }

    lambda <- params$beta/N*(I0+I1+I2)

    ##print(sprintf("proportion susceptible %s",S/(S+E-E1+I-I1+R-R1)))

    ##################
    ## susceptibles ##
    ##################
    dS0 <- -lambda*S0 - vac.rate[1]*S0/N0
    ## those unprotected after first dose
    dS1 <- -lambda*S1 + (1-params$theta1)*vac.rate[1]*S0/N0 - sec.vac.S1
    ## those still unprotected adter two doses
    dS2 <- -lambda*S2 + (1-params$theta2)/(1-params$theta1)*sec.vac.S1

    ############
    ## latent ##
    ############
    dE0  <- lambda*S0 - params$sigma*E0  - vac.rate[1]*E0/N0
    dE1 <- lambda*S1 - params$sigma*E1 + vac.rate[1]*E0/N0 - sec.vac.E1
    dE2 <- lambda*S2 - params$sigma*E2 + sec.vac.E1

    ################
    ## infectious ##
    ################

    dI0  <- params$sigma*E0 - params$gamma*I0 - vac.rate[1]*I0/N0
    dI1 <- params$sigma*E1 - params$gamma*I1 - sec.vac.I1 + vac.rate[1]*I0/N0
    dI2 <- params$sigma*E2 - params$gamma*I2 + sec.vac.I1

    #############
    ## removed ##
    #############
    dR0  <- params$gamma*I0 - vac.rate[1]*R0/N0
    dR1 <- params$gamma*I1 + params$theta1 * vac.rate[1]*S0/N0 +
        vac.rate[1]*R0/N0 - sec.vac.R1
    dR2 <- params$gamma*I2 + sec.vac.R1 + sec.vac.S1*((1-params$theta1)-(1-params$theta2))/(1-params$theta1)

    ########################################
    ## tracking/redundant state variables ##
    ########################################
    dV1 <- vac.rate[1] # number of 1st dose dispensed
    dV2 <- vac.rate[2] # number of second dose dispensed
    ## we define wasted vaccine as vaccine going to those already infected
    ## recovered or exposed

    ## NOTE: the definiation of W here is only wasted 1st dose
    ## this was changed to help calculated indivudal risk stuff
    dW  <- vac.rate[1]*(1-(S0/N0)) #+ wasted.sec #wasted vaccine
    #paste0("Time ", t,": ",dW)
    dCI  <- params$sigma*E0 #lambda*S #CI in unvaccinated class
    dCI1 <- params$sigma*E1 #lambda*(1-params$theta1)*S1 #CI in single vaccine dose recipients
    dCI2 <- params$sigma*E2 # lambda*(1-params$theta2)*S2 #CI in double dose recipients
    out  <- c(dS0,dS1,dS2,
              dE0,dE1,dE2,
              dI0,dI1,dI2,
              dR0,dR1,dR2,
              dV1,dV2,dW,
              dCI,dCI1,dCI2)
    return(list(out))
}

##' Assumes an imperfect vaccine that acts on
##' susceptibility VEs
##' refered to as susceptiblity reducing vaccine in
##' manuscript
##' @param t - time
##' @param x - state vector
##' @param params - parameter vectors
##' @param vac.starts - vaccination starts vector c(dose1, dose2)
##' @param vac.ends - vaccination ends vector c(dose1,dose2)
##' @param daily.vac - daily number of doses distributed
##' @return
##' @author Andrew Azman
seir.ves.dx.dt <- function(t,
                           x,
                           params,
                           vac.starts,
                           vac.ends,
                           daily.vac){

    ## states
    S  <- x[1]
    S1 <- x[2]
    S2 <- x[3]
    E  <- x[4]
    E1 <- x[5]
    E2 <- x[6]
    I  <- x[7]
    I1 <- x[8]
    I2 <- x[9]
    R  <- x[10]
    R1 <- x[11]
    R2 <- x[12]
    V1 <- x[13]
    V2 <- x[14]
    W1  <- x[15]
    W2  <- x[16]

    ## note that Vs and W  are simply used for tracking and
    ## aren't real states
    N <- sum(x[1:12])
    N0 <- S + E + I + R
    N1 <- S1 + E1 + I1 + R1
    N2 <- S2 + E2 + I2 + R2

    ## check that it is a time we want to vaccinate
    time.check <- t >= vac.starts & t <= vac.ends
    vaccinating.check <- any(time.check)

#    if(t > 145) recover()

    ## based on the check assign a vaccination rate
    vac.rate <- sapply(1:length(time.check),
                       function(x) ifelse(time.check[x],daily.vac,0))

    ## doing this so we don't have Inf when denominator is 0 before
    ## there are any people with a second dose
    if(time.check[2]){
        sec.vac.S <-vac.rate[2]*S1/N1
        sec.vac.E <-vac.rate[2]*E1/N1
        sec.vac.I <-vac.rate[2]*I1/N1
        sec.vac.R <-vac.rate[2]*R1/N1
        wasted.sec <- vac.rate[2]*(1-(S1/N1))
    } else {
        sec.vac.S<-sec.vac.E<-sec.vac.I<-sec.vac.R<-0
        wasted.sec<-0
    }

    lambda <- params$beta/N*(I+(1-params$phi1)*I1 + (1-params$phi2)*I2)
    ##print(sprintf("proportion susceptible %s",S/(S+E-E1+I-I1+R-R1)))

    ##################
    ## susceptibles ##
    ##################
    dS <-
        -lambda*S - vac.rate[1]*S/N0
    ## only those who have not been vaccinated are at risk for getting 1st dose
    dS1 <-
        -lambda*(1-params$theta1)*S1 + vac.rate[1]*S/N0 - sec.vac.S
    ## only those who have received 1st dose are at risk of second dose
    dS2 <-
        -lambda*(1-params$theta2)*S2 + sec.vac.S

    ## if (t < 1)
    ##     print(sprintf("Time: %s --> lambda = %10f,  dS=%5f, dS1=%5f,dS2=%5f",t,lambda*1e6,dS,dS1,dS2))

    ############
    ## latent ##
    ############
    dE  <-
        lambda*S - params$sigma*E - vac.rate[1]*E/N0
    dE1 <-
        lambda*(1-params$theta1)*S1 - params$sigma*E1 + vac.rate[1]*E/N0 - sec.vac.E
    dE2 <- lambda*(1-params$theta2)*S2 - params$sigma*E2 +  sec.vac.E

    ################
    ## infectious ##
    ################
    dI  <-
        params$sigma*E -  params$gamma*I -  vac.rate[1]*I/N0
    dI1 <-
        params$sigma*E1 - params$gamma*I1 + vac.rate[1]*I/N0 - sec.vac.I
    dI2 <- params$sigma*E2 - params$gamma*I2 + sec.vac.I

    #############
    ## removed ##
    #############
    dR  <- params$gamma*I -  vac.rate[1]*R/N0
    dR1 <- params$gamma*I1 + vac.rate[1]*R/N0 - sec.vac.R
    dR2 <- params$gamma*I2 + sec.vac.R

    ########################################
    ## tracking/redundant state variables ##
    ########################################
    dV1 <- vac.rate[1] # number of 1st dose dispensed
    dV2 <- vac.rate[2] # number of second dose dispensed
    dW1 <- vac.rate[1]*(1-(S/N0))
    dW2 <- wasted.sec #wasted vaccine

    ##print(sprintf("time: %s, W: %s",t,dW))
    ##if (dW < -1) recover()
    dCI  <- params$sigma*E #lambda*S  #CI in unvaccinated class
    dCI1 <- params$sigma*E1 #lambda*S1*(1-params$theta1) # #lambda*(1-params$theta1)*S1 #CI in single vaccine dose recipients
    dCI2 <- params$sigma*E2 #lambda*(1-params$theta2)*S2# #lambda*(1-params$theta2)*S2 #CI in double dose recipients
    out  <- c(dS,dS1,dS2,dE,dE1,dE2,
              dI,dI1,dI2,dR,dR1,dR2,
              dV1,dV2,dW1,dW2,
              dCI,dCI1,dCI2)
    return(list(out))
}

##' delta function for SIER model with VE_SP vaccine
##' refered to as severity reducing vaccine in paper (VE_SP)
##' @param t - time
##' @param x - state vector
##' @param params - params vector
##' @param vac.starts - first and second dose vaccinatino start times
##' @param vac.ends - first and second dose vaccination end times
##' @param daily.vac - daily vaccination rate
##' @return
##' @author Andrew Azman
seir.leaky.dx.dt <- function(t,
                             x,
                             params,
                             vac.starts,
                             vac.ends,
                             daily.vac){

    ################################
    ## setting up state variables ##
    ################################

    S <- x[1]
    S1 <- x[2]
    S2 <- x[3]
    E <- x[4]
    E1 <- x[5]
    E2 <- x[6]
    I <- x[7]
    I1 <- x[8]
    I2 <- x[9]
    A <- x[10]
    A1 <- x[11]
    A2 <- x[12]
    R <- x[13]
    R1 <- x[14]
    R2 <- x[15]
    V1 <- x[16]
    V2 <- x[17]
    W <- x[18]
    ## note that Vs and W  are simply used for tracking and
    ## aren't real states
    N <- sum(x[1:15])
    N0 <- S + E + I + A + R
    N1 <- S1 + E1 + I1 + A1 + R1
    N2 <- S2 + E2 + I2 + A2 + R2

    ##   cat("Time: ",t,"\n")
    ##   cat("S=",S,",I=",I,",A=",A,",R=",R,",V=",V,"\n")
    ##   cat("vac rate:",vac.daily.rates*t,"\n")

    #######################################################
    ## getting our vaccination timing and doses straight ##
    #######################################################

    time.check <- t >= vac.starts & t <= vac.ends
    vaccinating.check <- any(time.check)

    vac.rate <- sapply(1:length(time.check),
                       function(x) ifelse(time.check[x],daily.vac,0))

    if(time.check[2]){
        sec.vac.S <-vac.rate[2]*S1/N1
        sec.vac.E <-vac.rate[2]*E1/N1
        sec.vac.I <-vac.rate[2]*I1/N1
        sec.vac.A <-vac.rate[2]*A1/N1
        sec.vac.R <-vac.rate[2]*R1/N1
        wasted.sec <- vac.rate[2]*(1-(S1/N1))
    } else {
        sec.vac.S<-sec.vac.E<-sec.vac.I<-sec.vac.A <- sec.vac.R<-0
        wasted.sec<-0
    }

    #########
    ## FOI ##
    #########

    lambda <- params$beta/N*(I +
                             (1-params$phi1)*I1 +
                             (1-params$phi2)*I2 +
                             A +
                             (1-params$kappa1)*A1 +
                             (1-params$kappa2)*A2 )

    ##################
    ## susceptibles ##
    ##################

    dS<- -lambda*S - vac.rate[1]*S/N0
    ## only those who have not been vaccinated are at risk for getting 1st dose
    dS1<- -lambda*S1 + vac.rate[1]*S/N0 - sec.vac.S
    ## only those who have received 1st dose are at risk of second dose
    dS2<-  -lambda*S2 + sec.vac.S

    ##################
    ## latent folks ##
    ##################

    dE<- lambda*S - params$sigma*E - vac.rate[1]*E/N0
    dE1<- lambda*S1 - params$sigma*E1 + vac.rate[1]*E/N0 - sec.vac.E
    dE2<- lambda*S2 - params$sigma*E2 +  sec.vac.E

    ################
    ## infectious ##
    ################
    dI<- (1-params$theta0)*params$sigma*E - params$gamma*I - vac.rate[1]*I/N0
    dI1<- (1-params$theta1)*params$sigma*E1 - params$gamma*I1 + vac.rate[1]*I/N0 - sec.vac.I
    dI2<- (1-params$theta2)*params$sigma*E2 - params$gamma*I2 + sec.vac.I

    ####################
    ## asymptomatics  ##
    ####################
    dA <-  params$theta0*params$sigma*E - params$gamma*A - vac.rate[1]*A/N0
    dA1 <- params$theta1*params$sigma*E1 - params$gamma*A1 + vac.rate[1]*A/N0 - sec.vac.A
    dA2 <- params$theta2*params$sigma*E2 - params$gamma*A2 + sec.vac.A

    ###################
    ## removed folks ##
    ###################

    dR<- params$gamma*(I+A) - vac.rate[1]*R/N0
    dR1<- params$gamma*(I1+A1) +  vac.rate[1]*R/N0 - sec.vac.R
    dR2<- params$gamma*(I2+A2) + sec.vac.R

    #######################################
    ## redundant/utility state variables ##
    #######################################

    ## tracking total first and second vaccine doses
    dV1<- vac.rate[1]
    dV2<- vac.rate[2]

    ## wasted vaccine is vaccone that went
    ## to anyone except susceptibles
    dW<- vac.rate[1]*(1-(S/N0)) + wasted.sec

    ## tracking incidence of fully symptomatics only
    ## here we track those entering the infectious state
    ## for cholera this is a reasonable assumption
    ## asympotmatics are tracked seperatley
    dCI  <- (1-params$theta0)*params$sigma*E
    dCI1 <- (1-params$theta1)*params$sigma*E1
    dCI2 <- (1-params$theta2)*params$sigma*E2

    out <- c(dS,dS1,dS2,dE,dE1,dE2,dI,dI1,dI2,
             dA,dA1,dA2,dR,dR1,dR2,dV1,dV2,dW,dCI,dCI1,dCI2)
    return(list(out))
}

##' dx function for fast-slow VES model
##' this can accomidate as many slow compartments as
##' needed automatically as long as n.comps.slow is correct
##' @param t - time
##' @param x - state vector
##' @param vac.starts
##' @param vac.ends
##' @param daily.vac
##' @param params - parameter vector
##' @param n.comps.slow - number of slow compartments
##' @return
##' @author asa
ves.2path.dx.dt.generic <- function(t,
                                    x,
                                    vac.starts,
                                    vac.ends,
                                    daily.vac,
                                    params,
                                    n.comps.slow=3){

    ## get the names of each of the slow Infectious states
    I.star.states <- paste0(
        rep(c("I.star",
              paste0("I",1:2,".star")),n.comps.slow),
        rep(1:(n.comps.slow),each=3))

    ## names of all the states
    state.vars.all <- c("S","S1","S2",
                        "E","E1","E2",
                        "I","I1","I2",
                        I.star.states,
                        "R","R1","R2",
                        "V1","V2","W",
                        "CI","CI1","CI2")

    if (length(state.vars.all) != length(x)) stop("input state length not correct")

    dI.star  <- dI1.star <- dI2.star <- numeric(n.comps.slow)

    ## assign our state variables within the environment
    for (y in seq_along(state.vars.all)) assign(state.vars.all[y],x[y])

    n.real.states <- length(x) - 6
    N  <- sum(x[1:(9+n.comps.slow*3+3)])

    N0 <- sum(x[seq(1,n.real.states,by=3)])
    N1 <- sum(x[seq(2,n.real.states,by=3)])
    N2 <- sum(x[seq(3,n.real.states,by=3)])

    I.stars <- x[seq(10,n.real.states-3,by=3)]
    I1.stars <- x[seq(11,n.real.states-3,by=3)]
    I2.stars <- x[seq(12,n.real.states-3,by=3)]

    ## check that it is a time we want to vaccinate
    time.check <- t >= vac.starts & t <= vac.ends
    vaccinating.check <- any(time.check)

    ## based on the check assign a vaccination rate
    vac.rate <- sapply(1:length(time.check),
                       function(x) ifelse(time.check[x],daily.vac,0))

    if(time.check[2]){
        sec.vac.S <-vac.rate[2]*S1/N1
        sec.vac.E <-vac.rate[2]*E1/N1
        sec.vac.I <-vac.rate[2]*I1/N1
        ##  note that this is a vector
        sec.vac.I.star <- vac.rate[2]*I1.stars/N1
        sec.vac.R <-vac.rate[2]*R1/N1
        wasted.sec <- vac.rate[2]*(1-(S1/N1))
    } else {
        sec.vac.S<-sec.vac.E<-sec.vac.I<-sec.vac.R<-0
        sec.vac.I.star <- rep(0,n.comps.slow)
        wasted.sec<-0
    }

    ## force of infection (partitioned)
    lambda.short <- params$beta/N*(I+(1-params$phi1)*I1 + (1-params$phi2)*I2)
    lambda.long <-  params$beta.long/N*
        (sum(I.stars)+(1-params$phi1)*sum(I1.stars) + (1-params$phi2)*sum(I2.stars))

    lambda <- lambda.short + lambda.long

    ##################
    ## susceptibles ##
    ##################
    dS <- -lambda*S - vac.rate[1]*S/N0
    ## only those who have not been vaccinated are at risk for getting 1st dose
    dS1 <- -lambda*(1-params$theta1)*S1 + vac.rate[1]*S/N0 - sec.vac.S
    ## only those who have received 1st dose are at risk of second dose
    dS2 <- -lambda*(1-params$theta2)*S2 + sec.vac.S

    ############
    ## latent ##
    ############
    dE  <- lambda*S - params$sigma*E - vac.rate[1]*E/N0
    dE1 <- lambda*(1-params$theta1)*S1 - params$sigma*E1 + vac.rate[1]*E/N0 - sec.vac.E
    dE2 <- lambda*(1-params$theta2)*S2 - params$sigma*E2 +  sec.vac.E

    #######################
    ## infectious - fast ##
    #######################
    dI  <-  (1-params$nu)*params$sigma*E - params$gamma*I - vac.rate[1]*I/N0
    dI1 <-  (1-params$nu)*params$sigma*E1 - params$gamma*I1 + vac.rate[1]*I/N0 - sec.vac.I
    dI2 <-  (1-params$nu)*params$sigma*E2 - params$gamma*I2 + sec.vac.I

    #######################
    ## infectious - slow ##
    #######################
    ## initials
    dI.star[1]  <-
        params$nu*params$sigma*E - params$gamma.s*I.stars[1] - vac.rate[1]*I.stars[1]/N0
    dI1.star[1] <-
        params$nu*params$sigma*E1 - params$gamma.s*I1.stars[1] + vac.rate[1]*I1.stars[1]/N0 - sec.vac.I.star[1]

    dI2.star[1] <-
         params$nu*params$sigma*E2 - params$gamma.s*I2.stars[1] + sec.vac.I.star[1]

    ## now for subsequent compartments
    if (n.comps.slow > 1){
        for (i in 2:length(I.stars)){
            dI.star[i]  <-
                params$gamma.s*I.stars[i-1]  -  params$gamma.s*I.stars[i] - vac.rate[1]*I.stars[i]/N0
            dI1.star[i] <-
                params$gamma.s*I1.stars[i-1] - params$gamma.s*I1.stars[i] + vac.rate[1]*I1.stars[i]/N0 - sec.vac.I.star[i]
            dI2.star[i] <-
                params$gamma.s*I2.stars[i-1] - params$gamma.s*I2.stars[i] + sec.vac.I.star[i]
        }
    }

    #############
    ## removed ##
    #############
    dR  <- params$gamma*I + params$gamma.s*I.stars[n.comps.slow]  - vac.rate[1]*R/N0
    dR1 <- params$gamma*I1 + params$gamma.s*I1.stars[n.comps.slow]  +  vac.rate[1]*R/N0 - sec.vac.R
    dR2 <- params$gamma*I2 + params$gamma.s*I2.stars[n.comps.slow]  + sec.vac.R

    ########################################
    ## tracking/redundant state variables ##
    ################z########################
    dV1 <- vac.rate[1] # number of 1st dose dispensed
    dV2 <- vac.rate[2] # number of second dose dispensed
    dW  <- vac.rate[1]*(1-(S/N0)) + wasted.sec #wasted vaccine

    ##print(sprintf("time: %s, W: %s",t,dW))
    ##if (dW < -1) recover()
    dCI  <- params$sigma*E #lambda*S #CI in unvaccinated class
    dCI1 <- params$sigma*E1 #lambda*(1-params$theta1)*S1 #CI in single vaccine dose recipients
    dCI2 <- params$sigma*E2 # lambda*(1-params$theta2)*S2 #CI in double dose recipients

    ## totally inefficient but oh well
    I.star.out <- c()
    for (j in 1:n.comps.slow)
        I.star.out <- c(I.star.out,dI.star[j],dI1.star[j],dI2.star[j])

    out  <- c(dS,dS1,dS2,
              dE,dE1,dE2,
              dI,dI1,dI2,
              I.star.out,
              dR,dR1,dR2,
              dV1,dV2,
              dW,
              dCI,dCI1,dCI2)

    return(list(out))
}


##' dx function for fast-slow VES model
##' this one has 3 compartments hardcoded to speed things
##' up (after ftting the number of compartments)
##' @param t - time
##' @param x - state vector
##' @param vac.starts
##' @param vac.ends
##' @param daily.vac
##' @param params - parameter vector
##' @return
##' @author asa
ves.2path.dx.dt <- function(t,
                            x,
                            vac.starts,
                            vac.ends,
                            daily.vac,
                            params){

    dI.star  <- dI1.star <- dI2.star <- numeric(3)

    ## assign our state variables within the environment
    S <- x[1];S1 <- x[2];S2 <- x[3]
    E <- x[4];E1 <- x[5];E2 <- x[6]
    I <- x[7];I1 <- x[8];I2 <- x[9]
    ## I.star1 <- x[10];I1.star1 <- x[11];I2.star1 <- x[12]
    ## I.star2 <- x[11];I1.star2 <- x[12];I2.star2 <- x[13]
    ## I.star3 <- x[12];I1.star3 <- x[13];I2.star3 <- x[15]
    R <- x[16];R1 <- x[17];R2 <- x[18]

    n.real.states <- length(x) - 6
    N  <- sum(x[1:(9+9+3)])

    N0 <- sum(x[seq(1,n.real.states,by=3)])
    N1 <- sum(x[seq(2,n.real.states,by=3)])
    N2 <- sum(x[seq(3,n.real.states,by=3)])

    I.stars <- x[seq(10,n.real.states-3,by=3)]
    I1.stars <- x[seq(11,n.real.states-3,by=3)]
    I2.stars <- x[seq(12,n.real.states-3,by=3)]

    ## check that it is a time we want to vaccinate
    time.check <- t >= vac.starts & t <= vac.ends
    vaccinating.check <- any(time.check)

    ## based on the check assign a vaccination rate
    vac.rate <- sapply(1:length(time.check),
                       function(x) ifelse(time.check[x],daily.vac,0))

    if(time.check[2]){

        sec.vac.S <-vac.rate[2]*S1/N1
        sec.vac.E <-vac.rate[2]*E1/N1
        sec.vac.I <-vac.rate[2]*I1/N1
        ##  note that this is a vector
        sec.vac.I.star <- vac.rate[2]*I1.stars/N1
        sec.vac.R <-vac.rate[2]*R1/N1
        wasted.sec <- vac.rate[2]*(1-(S1/N1))

    } else {

        sec.vac.S<-sec.vac.E<-sec.vac.I<-sec.vac.R<-0
        sec.vac.I.star <- rep(0,3)
        wasted.sec<-0

    }


    lambda.short <- params$beta/N*(I+(1-params$phi1)*I1 + (1-params$phi2)*I2)
    lambda.long <-  params$beta.long/N*
        (sum(I.stars)+(1-params$phi1)*sum(I1.stars) + (1-params$phi2)*sum(I2.stars))

    lambda <- lambda.short + lambda.long

    ##################
    ## susceptibles ##
    ##################
    dS <- -lambda*S - vac.rate[1]*S/N0
    ## only those who have not been vaccinated are at risk for getting 1st dose
    dS1 <- -lambda*(1-params$theta1)*S1 + vac.rate[1]*S/N0 - sec.vac.S
    ## only those who have received 1st dose are at risk of second dose
    dS2 <- -lambda*(1-params$theta2)*S2 + sec.vac.S

    ############
    ## latent ##
    ############
    dE  <- lambda*S - params$sigma*E - vac.rate[1]*E/N0
    dE1 <- lambda*(1-params$theta1)*S1 - params$sigma*E1 + vac.rate[1]*E/N0 - sec.vac.E
    dE2 <- lambda*(1-params$theta2)*S2 - params$sigma*E2 +  sec.vac.E

    #######################
    ## infectious - slow ##
    #######################
    dI  <- (1-params$nu)*params$sigma*E - params$gamma*I - vac.rate[1]*I/N0
    dI1 <-  (1-params$nu)*params$sigma*E1 - params$gamma*I1 + vac.rate[1]*I/N0 - sec.vac.I
    dI2 <-  (1-params$nu)*params$sigma*E2 - params$gamma*I2 + sec.vac.I

    #######################
    ## infectious - fast ##
    #######################
    ## initials
    dI.star[1]  <-
        params$nu*params$sigma*E - params$gamma.s*I.stars[1] - vac.rate[1]*I.stars[1]/N0
    dI1.star[1] <-
        params$nu*params$sigma*E1 - params$gamma.s*I1.stars[1] + vac.rate[1]*I1.stars[1]/N0 - sec.vac.I.star[1]

    dI2.star[1] <-
        params$nu*params$sigma*E2 - params$gamma.s*I2.stars[1] + sec.vac.I.star[1]

    dI.star[2]  <-
        params$gamma.s*I.stars[1]  -  params$gamma.s*I.stars[2] - vac.rate[1]*I.stars[2]/N0
    dI1.star[2] <-
        params$gamma.s*I1.stars[1] - params$gamma.s*I1.stars[2] + vac.rate[1]*I1.stars[2]/N0 - sec.vac.I.star[2]
    dI2.star[2] <-
        params$gamma.s*I2.stars[1] - params$gamma.s*I2.stars[2] + sec.vac.I.star[2]

    dI.star[3]  <-
        params$gamma.s*I.stars[2]  -  params$gamma.s*I.stars[3] - vac.rate[1]*I.stars[3]/N0
    dI1.star[3] <-
        params$gamma.s*I1.stars[2] - params$gamma.s*I1.stars[3] + vac.rate[1]*I1.stars[3]/N0 - sec.vac.I.star[3]
    dI2.star[3] <-
        params$gamma.s*I2.stars[2] - params$gamma.s*I2.stars[3] + sec.vac.I.star[2]


    #############
    ## removed ##
    #############
    dR  <- params$gamma*I + params$gamma.s*I.stars[3]  - vac.rate[1]*R/N0
    dR1 <- params$gamma*I1 + params$gamma.s*I1.stars[3]  +  vac.rate[1]*R/N0 - sec.vac.R
    dR2 <- params$gamma*I2 + params$gamma.s*I2.stars[3]  + sec.vac.R

    ########################################
    ## tracking/redundant state variables ##
    ################z########################
    dV1 <- vac.rate[1] # number of 1st dose dispensed
    dV2 <- vac.rate[2] # number of second dose dispensed
    dW  <- vac.rate[1]*(1-(S/N0)) + wasted.sec #wasted vaccine

    ##print(sprintf("time: %s, W: %s",t,dW))
    ##if (dW < -1) recover()
    dCI  <- params$sigma*E #lambda*S #CI in unvaccinated class
    dCI1 <- params$sigma*E1 #lambda*(1-params$theta1)*S1 #CI in single vaccine dose recipients
    dCI2 <- params$sigma*E2 # lambda*(1-params$theta2)*S2 #CI in double dose recipients

    I.star.out <- c(dI.star[1],dI1.star[1],dI2.star[1],
                    dI.star[2],dI1.star[2],dI2.star[2],
                    dI.star[3],dI1.star[3],dI2.star[3])

    out  <- c(dS,dS1,dS2,
              dE,dE1,dE2,
              dI,dI1,dI2,
              I.star.out,
              dR,dR1,dR2,
              dV1,dV2,
              dW,
              dCI,dCI1,dCI2)

    return(list(out))
}


## ---------------------------------------------------------------- ##
## Section 2: Helper Functions                                      ##
## Note that some of these are not actually used for the manuscript ##
## ---------------------------------------------------------------- ##

##' function used to minizimze when finidng the minimum VEs
##' @param ve1 - 1st-dose VE
##' @param ve2 - 2nd-dose VE
##' @param start.time - start simluation time
##' @param time.step - simulation time step
##' @param max.time - maximum time for simulations
##' @param doses - number of doses available
##' @param N - poulation size
##' @param inter.dose.time - time between doses
##' @param total.vac.days - total vaccination days
##' @param dx.dt.func - delta function to use for simulations
##' @param initial.state - initial state vector
##' @param params - paramter vector
##' @param model.type - string of model type corrspnding to dx.dt.func (redundant but needed)
##' @param hmax - option for maximum step size in numerical integreation
##' @param time.to.prot - vector of time to protection from 1st and second dose
##' @param inc.cutoff - cut off of number of new cases for calling an epidemic over
##' @return absolute difference between single and two-dose campaigns
##' @author asa
min.ve.1.obj.func <- function(ve1=0.35,
                              ve2=0.75,
                              start.time=1,
                              time.step=1,
                              max.time=350,
                              doses=10000,
                              N=5e5,#1e6,
                              inter.dose.time=14,
                              total.vac.days=10,
                              dx.dt.func=seir.ves.dx.dt,
                              initial.state=c(N-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                              params=list(
                                  beta=0.6,
                                  gamma=1/2,
                                  sigma=1/1.4,
                                  phi1=0,
                                  phi2=0,
                                  kappa=10^6,
                                  epsilon=10,
                                  delta=1/1.5,
                                  theta1=ve1,
                                  theta2=ve2
                                  ),
                              model.type="ves",
                              hmax=.2,
                              time.to.prot=c(0,0),
                              inc.cutoff=0.1){

    params$theta1 <- ve1
    params$theta2 <- ve2
    times <- seq(0,max.time,time.step)

    ## fixing the total length of the vaccination campaign

    daily.vac.doses <- doses/total.vac.days

    ## set up timings for both scenarios

    ## first let's look at time.to.prot to see when protection
    ## kicks in. Vaccination kicks in right away in this model, so
    ## we will simulate delayed protection by effectivley starting
    ## the campaigns later

    one.of.one.start <- start.time + time.to.prot[1]
    one.of.one.end   <- one.of.one.start + total.vac.days
    cat(paste0("single dose campaign:",one.of.one.start-1,"-",one.of.one.end-1,"\n"))
    one.of.two.start <- start.time + time.to.prot[1]
    one.of.two.end   <- one.of.two.start+total.vac.days/2
    two.of.two.start <- one.of.two.end-time.to.prot[1]+inter.dose.time+time.to.prot[2]
    two.of.two.end   <- two.of.two.start + total.vac.days/2
    cat(paste0("two-dose campaign, first dose:",one.of.two.start-1,"-",one.of.two.end-1,"\n"))
    cat(paste0("two-dose campaign, second dose:",two.of.two.start-1,"-",two.of.two.end-1,"\n"))

    run.single.vac <- ode(y=initial.state,
                          times=times,
                          func=dx.dt.func,
                          parms=params,
                          vac.starts=c(one.of.one.start,Inf),
                          vac.ends=c(one.of.one.end,Inf),
                          daily.vac=daily.vac.doses,
                          hmax = hmax)

    run.no.vac <- ode(y=initial.state,
                      times=times,
                      func=dx.dt.func,
                      parms=params,
                      vac.starts=c(Inf,Inf),
                      vac.ends=c(Inf,Inf),
                      daily.vac=daily.vac.doses,
                      hmax = hmax)

    run.no.vac <- add.colnames(run.no.vac,model.type)

    run.single.vac <- add.colnames(run.single.vac,model.type)

    ## check here
    long.enough <-
        check.epi.ended(run=run.single.vac,inc.cutoff=inc.cutoff)

    while(!long.enough){
        times <- seq(0,max(times)+50,time.step)
        print(paste0("extending time - single dose (max time =",max(times),")"))

        run.single.vac <- ode(y=initial.state,
                              times=times,
                              func=dx.dt.func,
                              parms=params,
                              vac.starts=c(one.of.one.start,Inf),
                              vac.ends=c(one.of.one.end,Inf),
                              daily.vac=daily.vac.doses,
                              hmax = hmax)

        run.single.vac <- add.colnames(run.single.vac,model.type)
        long.enough <- check.epi.ended(run=run.single.vac,inc.cutoff=inc.cutoff)
    }

    run.double.vac <- ode(y=initial.state,
                          times=times,
                          func=dx.dt.func,
                          parms=params,
                          vac.starts=c(one.of.two.start,two.of.two.start),
                          vac.ends=c(one.of.two.end,two.of.two.end),
                          daily.vac=daily.vac.doses,
                          hmax = hmax)

    run.double.vac <- add.colnames(run.double.vac,model.type)


    ## check here
    long.enough <-
        check.epi.ended(run=run.double.vac,inc.cutoff=inc.cutoff)

    while(!long.enough){
        ##print("extending time - double dose")
        times <- seq(0,max(times)+50,time.step)

        run.double.vac <- ode(y=initial.state,
                              times=times,
                              func=dx.dt.func,
                              parms=params,
                              vac.starts=c(one.of.two.start,two.of.two.start),
                              vac.ends=c(one.of.two.end,two.of.two.end),
                              daily.vac=daily.vac.doses,hmax = hmax)

        run.double.vac <- add.colnames(run.double.vac,model.type)

        long.enough <- check.epi.ended(run=run.double.vac,inc.cutoff=inc.cutoff)
    }

    return(
        abs(
            sum(tail(run.single.vac[,c("CI","CI1","CI2")],1))
            -
            sum(tail(run.double.vac[,c("CI","CI1","CI2")],1))
            )
        )
}

##' adds column names to seir.impvac.dx.dt
##' @param mod.out - output from ode()
##' @return ode data.frame with appropriatley named columns
add.colnames <- function(mod.out,model.type="ves"){
    if (model.type=="ves"){
        colnames(mod.out) <- c("time",
                               unlist(strsplit("S,S1,S2,E,E1,E2,I,I1,I2,R,R1,R2,V1,V2,W1,W2,CI,CI1,CI2",",")))
    } else if(model.type == "leaky"){
        colnames(mod.out) <- c("time",
                               unlist(strsplit("S,S1,S2,E,E1,E2,I,I1,I2,A,A1,A2,R,R1,R2,V1,V2,W,CI,CI1,CI2",",")))
    } else if(model.type == "sirw"){
        colnames(mod.out) <- c("time",
                               unlist(strsplit("S,S1,S2,I,I1,I2,R,R1,R2,C,V1,V2,W,CI,CI1,CI2",",")))
    } else if(model.type == "test"){
        colnames(mod.out) <- c("time",
                               unlist(strsplit("S,S1,S2,I,I1,I2,R,R1,R2,CI,CI1,CI2",",")))

    } else if(model.type == "sirb"){
        colnames(mod.out) <- c("time",
                               unlist(strsplit("S,S1,S2,E,E1,E2,I,I1,I2,R,R1,R2,B,V1,V2,W,CI,CI1,CI2",",")))
    } else if(model.type == "2path"){
        slow.comps <- (dim(mod.out)[2] - 19)/3
        I.star.states <- c("I.star",paste0("I.star",1:(slow.comps-1)))
        I1.star.states <- c("I1.star",paste0("I1.star",1:(slow.comps-1)))
        I2.star.states <- c("I2.star",paste0("I2.star",1:(slow.comps-1)))
        colnames(mod.out) <- c("time",
                               "S","S1","S2",
                               "E","E1","E2",
                               "I","I1","I2",
                               I.star.states,
                               I1.star.states,
                               I2.star.states,
                               "R","R1","R2",
                               "V1","V2","W",
                               "CI","CI1","CI2")


    } else if(model.type =="aon") {
        colnames(mod.out) <- c("time",
                               unlist(strsplit("S0,S1,S2,E0,E1,E2,I0,I1,I2,R0,R1,R2,V1,V2,W,CI,CI1,CI2",",")))

    } else {
        warning("Don't know that model type")
    }

    return(mod.out)
}

##' @title
##' @param pdf
##' @param run.nums
##' @return
##' @author ANDREW AZMAN
make.vecontour.plot <- function(fs.mat.1dose,
                                fs.mat.2dose,
                                uncon.run,
                                vac.timings=seq(0,210,by=7),
                                ve1s=seq(.05,.6,by=.025),
                                vacdoses=seq(50000,900000,by=50000),
                                run.nums=1:18,
                                pdf=FALSE,
                                plot.name="veplot1.pdf"
                                ){

    require(fields)

    if(pdf){
        pdf(paste0("Plots/",plot.name),width=5,height=3)
    }else{
        quartz("",width=5,height=3)
    }
    par(mfrow=c(1,1),mar=c(0.5,3,.5,3),oma=c(1,1,2,3),mgp = c(2,.5,0))
    layout(mat=matrix(c(1,1,2,2,2,2),byrow=T,nrow=3))

                                        # first plot the epidemic curve
    max.time.id <- which(uncon.run[,1]==max(vac.timings))
    plot(c(0,diff(uncon.run[1:max.time.id,c("CI")])),type="l",
         ylab="Incident Cases",
         xlab="",axes=FALSE,
         xlim=c(0,max.time.id))

    abline(v=seq(0,max.time.id,length=length(vac.timings)),col=AddAlpha("grey",.1))
    abline(h=seq(0,max(diff(uncon.run[1:max.time.id,c("CI")])),length=length(ve1s)),col=AddAlpha("grey",.1))

    #recover()
    labs <- seq(0,max(vac.timings)/7,by=1)
    axis(3,
         labels=labs,
         at=seq(0,max.time.id,length=length(labs))
         ,col = "grey",cex.axis=1)

    axis(2,
         col = "grey",
         cex.axis=1)

    for (dosenum in seq_along(run.nums)){
        fs.1dose <- fs.mat.1dose[,,run.nums[dosenum],1] + fs.mat.1dose[,,run.nums[dosenum],2] + fs.mat.1dose[,,run.nums[dosenum],3]
        fs.2dose <- fs.mat.2dose[,,run.nums[dosenum],1] + fs.mat.2dose[,,run.nums[dosenum],2] + fs.mat.2dose[,,run.nums[dosenum],3]
        fs.stat <- 1-(fs.1dose/fs.2dose)
        (zlims <- range(fs.stat))
        ## make pallete with white in the center
        ## lower.col <- colorRampPalette(brewer.pal(8,"Blues"))(50)
        ## upper.col <- colorRampPalette(brewer.pal(8,"Reds"))(50)
        ## my.col <-c(lower.col,"white",upper.col)
        my.col <- colorRampPalette(brewer.pal(8,"Blues"))(18)

                                        #layout(mat=matrix(c(1,2,3,4,5,6),nrow=2))
        if (dosenum == 1){
            par(mar=c(3,3,.5,3))
            plot(-100,-100,
                 xlim=c(0,1),
                 ylim=c(0,1),
                 xlab="Epidemic Week",
             ,axes=FALSE,
                 ylab="Minimum Single Dose VE")
            ## image(t(fs.stat),
            ##       col="white",
            ##       zlim=zlims,
            ##       axes=F)
            abline(v=seq(0,1,length=length(vac.timings)),col=AddAlpha("grey",.1))
            abline(h=seq(0,1,length=length(ve1s)),col=AddAlpha("grey",.1))

            axis(2,labels=
                 round(ve1s*100,1)[seq(1,length(ve1s),by=3)]
                 ,at=seq(0,1,length=length(ve1s))[seq(1,length(ve1s),by=3)]
                 ,col = "grey",cex.axis=1)

            axis(1,
                 labels=labs,
                 ,at=seq(0,1,length=length(labs))
                 ,col = "grey",cex.axis=1)

    }

        contour(t(fs.stat),levels=0,add=T,col=2,
                #my.col[dosenum],
                labels="")
}

    ##     par(oma=c( 0,0,0,1))# reset margin to be much smaller.
    ## set.panel() # reset plotting device

    image.plot(matrix(c(1:length(vacdoses)),nrow=1),legend.only=TRUE,
               col=my.col,smallplot = c(.94,.97,.25,.85),
               axis.args=list(at=seq(1,length(vacdoses),by=3),
                   labels=vacdoses[seq(1,length(vacdoses),by=3)]))
    if(pdf) dev.off()
}



##'
##' @param fs.mat.1dose
##' @param fs.mat.2dose
##' @param time.at.risk.1dose
##' @param time.at.risk.2dose
##' @param ve1s
##' @param vacdoses
##' @param vac.timings
##' @param dosenums
##' @param add.unvac
##' @param pdf
##' @param plot.name
##' @param round.me rounds log.IRR to the nearest .001 to smooth over numerical errors from integration
##' @param mrse - y-axis as relative VE (assumed that ve2 is the max of ve1s vector)
##' @return
##' @author ANDREW AZMAN
make.counterfact.plot <- function(fs.mat.1dose,
                                  fs.mat.2dose,
                                  time.at.risk.1dose,
                                  time.at.risk.2dose,
                                  ve1s=seq(.05,.6,by=.025),
                                  vacdoses=seq(50000,900000,by=50000),
                                  vac.timings = seq(0,210,by=7),
                                  dosenums=c(2,10,18),
                                  add.unvac=TRUE,
                                  pdf=FALSE,
                                  plot.name="ceplot.pdf",
                                  round.me=TRUE,
                                  mrse=FALSE
                                  ){

    require(fields)
    if(pdf){
        pdf(paste0("Plots/",plot.name),width=5,height=3)
    } else {
        quartz("",width=6,height=3)
    }

    log.irr.vac <- list()
    ## calculate the incidencce rate ratios among those vaccinated
    log.irr.vac[[1]] <-
        log((fs.mat.1dose[,,dosenums[1],2]/time.at.risk.1dose[,,dosenums[1],2]) /
            ((fs.mat.2dose[,,dosenums[1],2] + fs.mat.2dose[,,dosenums[1],3])/
             (time.at.risk.2dose[,,dosenums[1],2] + time.at.risk.2dose[,,dosenums[1],3])))

    log.irr.vac[[2]] <-
        log((fs.mat.1dose[,,dosenums[2],2]/time.at.risk.1dose[,,dosenums[2],2]) /
            ((fs.mat.2dose[,,dosenums[2],2] + fs.mat.2dose[,,dosenums[2],3])/
             (time.at.risk.2dose[,,dosenums[2],2] + time.at.risk.2dose[,,dosenums[2],3])))

    log.irr.vac[[3]] <-
        log((fs.mat.1dose[,,dosenums[3],2]/time.at.risk.1dose[,,dosenums[3],2]) /
            ((fs.mat.2dose[,,dosenums[3],2] + fs.mat.2dose[,,dosenums[3],3])/
             (time.at.risk.2dose[,,dosenums[3],2] + time.at.risk.2dose[,,dosenums[3],3])))


    log.irr.unvac <- list()

    log.irr.unvac[[1]] <-
        log((fs.mat.1dose[,,dosenums[1],1]/time.at.risk.1dose[,,dosenums[1],1]) /
            (fs.mat.2dose[,,dosenums[1],1]/time.at.risk.2dose[,,dosenums[1],1]))

    log.irr.unvac[[2]] <-
        log((fs.mat.1dose[,,dosenums[2],1]/time.at.risk.1dose[,,dosenums[2],1]) /
            (fs.mat.2dose[,,dosenums[2],1]/time.at.risk.2dose[,,dosenums[2],1]))

    log.irr.unvac[[3]] <-
        log((fs.mat.1dose[,,dosenums[3],1]/time.at.risk.1dose[,,dosenums[3],1]) /
            (fs.mat.2dose[,,dosenums[3],1]/time.at.risk.2dose[,,dosenums[3],1]))

    if (add.unvac){
        zlims <-
            range(c(log.irr.vac,log.irr.unvac))
        par(mfrow=c(2,3),mar=c(.5,.5,.5,.5),oma=c(3,5,3,7),mgp = c(2,.5,0))
    } else {
        par(mfrow=c(1,3),mar=c(.5,.5,.5,.5),oma=c(3,3,3,5),mgp = c(2,.5,0))
        zlims <- range(log.irr.vac)
    }

    ## make pallete with white in the center
    print(zlims)
    if(zlims[1] > -1e-2) zlims[1] <- -.05
    if(zlims[2] < 1e-2) zlims[2] <- .05

    my.breaks<-c(seq(zlims[1],-1e-2, length.out=31),0,seq(1e-2, zlims[2], length.out=30))
    my.col <- colorRampPalette(c(rev(brewer.pal(5,"Greens")),"white",brewer.pal(5,"Blues")))(61)
    my.contour.levels <- 0#round(min(zlims[1])):round(zlims[2])

    log.irr.vac <- lapply(log.irr.vac,function(x) round(x,3))
    log.irr.unvac <- lapply(log.irr.unvac,function(x) round(x,3))

    image(t(log.irr.vac[[1]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)
    contour(t(log.irr.vac[[1]]),levels=my.contour.levels,add=T)

    axis(2,labels=
         round(ve1s[seq(1,length(ve1s),by=2)]*100,1)
         ,at=seq(0,1,length=length(ve1s))[seq(1,length(ve1s),by=2)]
         ,col = "grey",cex.axis=1)

    if(add.unvac){
        axis(1,labels=rep("",length(seq(1,length(vac.timings),by=2)))
             ,at=seq(0,1,length=length(vac.timings))[seq(1,length(vac.timings),by=2)]
             ,col = "grey",cex.axis=1)
    } else {
        axis(1,labels=
             (round(vac.timings,1)/7)[seq(1,length(vac.timings),by=2)]
             ,at=seq(0,1,length=length(vac.timings))[seq(1,length(vac.timings),by=2)]
             ,col = "grey",cex.axis=1)
    }

    image(t(log.irr.vac[[2]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)
    contour(t(log.irr.vac[[2]]),levels=my.contour.levels,add=T)


    ## axis(2,labels=rep("",length(seq(1,length(ve1s),by=2)))
    ##      ,at=seq(0,1,length=length(ve1s))[seq(1,length(ve1s),by=2)]
    ##      ,col = "grey",cex.axis=1)

    if (add.unvac){
        axis(1,labels=rep("",length(seq(1,length(vac.timings),by=2)))
             ,at=seq(0,1,length=length(vac.timings))[seq(1,length(vac.timings),by=2)]
             ,col = "grey",cex.axis=1)

    } else {
        axis(1,labels=
             round(vac.timings,1)
             ,at=seq(0,1,length=length(vac.timings))
             ,col = "grey",cex.axis=1)
    }

    image(t(log.irr.vac[[3]]),breaks=my.breaks,zlim=zlims, col=my.col,axes=FALSE)
    contour(t(log.irr.vac[[3]]),levels=my.contour.levels,add=T)

    axis(4,labels=
         round(ve1s[seq(1,length(ve1s),by=2)]*100,1)
         ,at=seq(0,1,length=length(ve1s))[seq(1,length(ve1s),by=2)]
         ,col = "grey",cex.axis=1)

    if(add.unvac){
        axis(1,labels=rep("",length(seq(1,length(vac.timings),by=2)))
             ,at=seq(0,1,length=length(vac.timings))[seq(1,length(vac.timings),by=2)]
             ,col = "grey",cex.axis=1)

    } else{
        axis(1,labels=
             round(vac.timings,1)
             ,at=seq(0,1,length=length(vac.timings))
             ,col = "grey",cex.axis=1)
    }

    ##########################
    ## now for unvaccinated ##
    ##########################
    if (add.unvac){
        image(t(log.irr.unvac[[1]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)
        contour(t(log.irr.unvac[[1]]),levels=my.contour.levels,add=T)

        axis(2,labels=
             round(ve1s[seq(1,length(ve1s),by=2)]*100,1)
             ,at=seq(0,1,length=length(ve1s))[seq(1,length(ve1s),by=2)]
             ,col = "grey",cex.axis=1)

        axis(1,labels=
             (round(vac.timings,1)/7)[seq(1,length(vac.timings),by=2)]
             ,at=seq(0,1,length=length(vac.timings))[seq(1,length(vac.timings),by=2)]
             ,col = "grey",cex.axis=1)


        image(t(log.irr.unvac[[2]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)
        contour(t(log.irr.unvac[[2]]),levels=my.contour.levels,add=T)

        ## axis(2,labels=
        ##      rep("",length(seq(1,length(ve1s),by=2)))
        ##      ,at=seq(0,1,length=length(ve1s))[seq(1,length(ve1s),by=2)]
        ##      ,col = "grey",cex.axis=1)

        axis(1,labels=
             (round(vac.timings,1)/7)[seq(1,length(vac.timings),by=2)]
             ,at=seq(0,1,length=length(vac.timings))[seq(1,length(vac.timings),by=2)]
             ,col = "grey",cex.axis=1)

        image(t(log.irr.unvac[[3]]),breaks=my.breaks,zlim=zlims, col=my.col,axes=FALSE)
        contour(t(log.irr.unvac[[3]]),levels=my.contour.levels,add=T)

        axis(4,labels=
             round(ve1s[seq(1,length(ve1s),by=2)]*100,1)
             ,at=seq(0,1,length=length(ve1s))[seq(1,length(ve1s),by=2)]
             ,col = "grey",cex.axis=1)

        axis(1,labels=
             (round(vac.timings,1)/7)[seq(1,length(vac.timings),by=2)]
             ,at=seq(0,1,length=length(vac.timings))[seq(1,length(vac.timings),by=2)]
             ,col = "grey",cex.axis=1)
    }
    ########################
    ## now for the legned ##
    ########################

    par(oma=c( 0,0,0,1))# reset margin to be much smaller.
    set.panel() # reset plotting device

    image.plot(log.irr.vac[[1]],legend.only=TRUE,zlim=zlims,
               breaks = my.breaks,
               col=my.col,
               smallplot = c(.94,.96,.2,.8))

    mtext(paste0(round(vacdoses[dosenums[1]]/.5e4,1),"/",round(vacdoses[dosenums[1]]/.5e4/2,1)),side=3,outer=TRUE,line=-2,at=.25,cex=.8)
    mtext(paste0(round(vacdoses[dosenums[2]]/.5e4,1),"/",round(vacdoses[dosenums[2]]/.5e4/2,1)),side=3,outer=TRUE,line=-2,at=.5,cex=.8)
    mtext(paste0(round(vacdoses[dosenums[3]]/.5e4,1),"/",round(vacdoses[dosenums[3]]/.5e4/2,1)),,side=3,outer=TRUE,line=-2,at=.75,cex=.8)
    mtext("Vaccination Start (Epidemic Week)",side=1,outer=TRUE,line=-1.4,at=.5,cex=.8)
    if(mrse){
        mtext('Relative Single Dose VE',side=2,outer=TRUE,line=-2.6,at=.5,cex=.8)
    } else {
        mtext('Single Dose VE',side=2,outer=TRUE,line=-2.6,at=.5,cex=.8)
    }
    mtext("Vaccinated",side=2,outer=TRUE,line=-1.2,at=.7,cex=.8)
    mtext("Unvaccinated",side=2,outer=TRUE,line=-1.2,at=.3,cex=.8)

    if(pdf) dev.off()
}

##' Returns the R0 for a
##' @param params
##' @param model.type
##' @return
##' @author ANDREW AZMAN
get.R0 <- function(params,model.type){
    if(model.type=="sirw"){
        rc <- (params$beta_C * params$alpha/params$zeta + params$beta_I)/params$gamma
    } else if(model.type == "ves" | model.type == "leaky"){
        rc <- (params[,"beta"])/params[,"gamma"]
    } else {
        stop("unknown model type")
    }
    return(rc)
}


##' This makes Figure 4 from the manuscript
##' @param fs.mat.1dose
##' @param fs.mat.2dose
##' @param time.at.risk.1dose
##' @param time.at.risk.2dose
##' @param ve1s
##' @param vacdoses
##' @param vac.timings
##' @param dosenums
##' @param min.ve1s - array with each column repreenting the min ve by time for a specifc dose
##' @param min.ve1.times - vector of times where min ve was estimates
##' @param add.unvac
##' @param pdf
##' @param plot.name
##' @param round.me
##' @param mrse
##' @param irr.all - do we want the population level threshold to be based on IRR (TRUE) or Final Size (FALSE)
##' @return
##' @author asa
make.counterfact.plot2 <- function(fs.mat.1dose,
                                   fs.mat.2dose,
                                   time.at.risk.1dose,
                                   time.at.risk.2dose,
                                   ve1s=seq(0,0.80,by=.02),
                                   vacdoses=seq(50000,355000,by=10000),
                                   vac.timings = seq(0,300,by=2),
                                   dosenums=c(10,20,30),
                                   min.ve1s,
                                   min.ve1.times,
                                   add.unvac=TRUE,
                                   pdf=FALSE,
                                   plot.name="ceplot.pdf",
                                   round.me=TRUE,
                                   mrse=FALSE,
                                   irr.all=FALSE
                                   ){

    require(fields)
    if(pdf){
        pdf(paste0("Plots/",plot.name),width=9,height=3)
    } else {
        quartz("",width=7.5,height=2.3)
    }

    cex.axis <- 0.9

    ## set up some lists to hold our transformed data
    ## log.irr.vac <- log.irr.unvac <- fs.diff <- log.irr.all <- list()
    arr.vac <- arr.all <- list()

    ## calculate the incidencce rate ratios among those vaccinated
    for (i in 1:length(dosenums)){

        fs.1dose.all <- fs.mat.1dose[,,dosenums[i],1] + fs.mat.1dose[,,dosenums[i],2] +
          fs.mat.1dose[,,dosenums[i],3]
        fs.2dose.all <- fs.mat.2dose[,,dosenums[i],1] + fs.mat.2dose[,,dosenums[i],2] +
          fs.mat.2dose[,,dosenums[i],3]
        arr.all[[i]] <- fs.1dose.all/fs.2dose.all

        fs.1dose.vac <- fs.mat.1dose[,,dosenums[i],4]#fs.mat.1dose[,,dosenums[i],2] + fs.mat.1dose[,,dosenums[i],4]
        fs.2dose.vac <- fs.mat.2dose[,,dosenums[i],4]#fs.mat.2dose[,,dosenums[i],2] + fs.mat.2dose[,,dosenums[i],3] + fs.mat.2dose[,,dosenums[i],4]
        arr.vac[[i]] <- fs.1dose.vac/(2*fs.2dose.vac)

    }

    warning("time at risk for vacinnees is not quite right at the moment")

    ## just to smooth things out a little bit
    n.digits <- 6
    max.t.ind <- 71 # trming to 20 weeks
    log.arr.vac <- lapply(arr.vac,function(x) round(log(x[,1:max.t.ind]),n.digits))
    log.arr.all <- lapply(arr.all,function(x) round(log(x[,1:max.t.ind]),n.digits))

    ## log.irr.unvac <- lapply(log.irr.unvac,function(x) round(x,5))
    ## fs.diff <- lapply(fs.diff,function(x) round(x,5))

    if (add.unvac){
        zlims <- range(c(log.irr.vac,log.irr.unvac))
        par(mfrow=c(2,3),mar=c(.5,.5,.5,.5),oma=c(3,5,3,7),mgp = c(2,.5,0),tck=-.02)
    } else {
        par(mfrow=c(1,3),mar=c(.5,.5,3,.5),oma=c(3,3,0,5),mgp = c(4,.3,0),tck=-.02)
        zlims <- range(c(log(do.call(c,arr.all)),
                         log(do.call(c,arr.vac))))
        ##zlims <- range(log.irr.vac)
    }

    ## make pallete with white in the center
    if(zlims[1] > -1e-2) zlims[1] <- -.05
    if(zlims[2] < 1e-2) zlims[2] <- .05
    print(zlims)

    my.breaks<-c(seq(zlims[1],-1e-2, length.out=31),0,seq(1e-2, zlims[2], length.out=30))
    my.col <- colorRampPalette(c(rev(brewer.pal(5,"Greens")),"white",brewer.pal(5,"Blues")))(61)
    my.contour.levels <- 0#round(min(zlims[1])):round(zlims[2])

    # for x-axis labels

    ## find max weeks
    max.weeks <- max(vac.timings[1:max.t.ind]/7)
    week.labels <- seq(0,max.weeks,by=4)
    week.ats <- seq(0,20/max.weeks,length=6)
    week.grids <- week.ats
    # for y-axis labels

    if (mrse){
        ve1.labels <-seq(0,1,by=.1)
        ve1.at <-seq(0,1,by=.1)
        ve1.grids <- ve1.at
    } else {
        ve1.labels <-seq(0,max(ve1s),by=.05)
        ve1.at <- approx(x=ve1s,y=seq(0,1,length=length(ve1s)),xout=ve1.labels)$y
        ve1.grids <- ve1.at
    }

    ## image(t(log.irr.vac[[1]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)

    image(t(log.arr.vac[[1]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)
    #image(t(arr.vac[[1]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)

    abline(h=ve1.grids,col=AddAlpha("grey",.2))
    abline(v=week.grids,col=AddAlpha("grey",.2))

    #contour(t(log.irr.vac[[1]]),levels=my.contour.levels,add=T)
    ## add contour lines for the vaccinated
    contour(t(log.arr.vac[[1]]),levels=my.contour.levels,add=T)
    #contour(t(arr.vac[[1]]),levels=1,add=T)

    ## add contour lines for the population
    contour(t(log.arr.all[[1]]),levels=my.contour.levels,lty=2,col=2,add=T)

    ## if(irr.all){
    ##     contour(t(log.irr.all[[1]]),levels=0,lty=2,col=2,add=T)
    ## } else {
    ##     contour(t(fs.diff[[1]]),levels=0,lty=2,col=2,add=T)
    ## }

    title(paste0(round(vacdoses[dosenums[1]]/.5e4,1),"/",round(vacdoses[dosenums[1]]/.5e4/2,1),"% Vaccinated")
         ,line=0.5)
    axis(2,labels=ve1.labels,at=ve1.at,col = "grey",cex.axis=cex.axis)

    if(add.unvac){
        axis(1,labels=rep("",length(week.labels))
            ,at=week.ats,col="grey",cex.axis=cex.axis)
    } else {
        axis(1,labels=week.labels,at=week.ats,col="grey",cex.axis=cex.axis)
    }

    ## image(t(log.irr.vac[[2]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)
    image(t(log.arr.vac[[2]]),breaks=my.breaks,zlim=zlims,col=my.col,axes=FALSE)


    abline(h=ve1.grids,col=AddAlpha("grey",.2))
    abline(v=week.grids,col=AddAlpha("grey",.2))
    ## contour(t(log.irr.vac[[2]]),levels=my.contour.levels,add=T)
    contour(t(log.arr.vac[[2]]),levels=my.contour.levels,add=T)
    contour(t(log.arr.all[[2]]),levels=my.contour.levels,lty=2,col=2,add=T)

    ## if(irr.all){
    ##     contour(t(log.irr.all[[2]]),levels=0,lty=2,col=2,add=T)
    ## } else {
    ##     contour(t(fs.diff[[2]]),levels=0,lty=2,col=2,add=T)
    ## }

    title(paste0(round(vacdoses[dosenums[2]]/.5e4,1),"/",round(vacdoses[dosenums[2]]/.5e4/2,1),"% Vaccinated")
         ,line=0.5)


    ## axis(2,labels=rep("",length(seq(1,length(ve1s),by=2)))
    ##      ,at=seq(0,1,length=length(ve1s))[seq(1,length(ve1s),by=2)]
    ##      ,col = "grey",cex.axis=1)

    if (add.unvac){
        axis(1,labels=rep("",length(week.labels))
            ,at=week.ats,col="grey",cex.axis=cex.axis)

    } else {
        axis(1,labels=week.labels,at=week.ats,col="grey",cex.axis=cex.axis)
    }

    axis(2,labels=rep("",length(ve1.labels)),at=ve1.at
        ,col = "grey",cex.axis=cex.axis)
    axis(4,labels=rep("",length(ve1.labels)),at=ve1.at
        ,col = "grey",cex.axis=cex.axis)

    image(t(log.arr.vac[[3]]),breaks=my.breaks,zlim=zlims, col=my.col,axes=FALSE)
    abline(h=ve1.grids,col=AddAlpha("grey",.2))
    abline(v=week.grids,col=AddAlpha("grey",.2))
    contour(t(log.arr.vac[[3]]),levels=my.contour.levels,add=T)
    contour(t(log.arr.all[[3]]),levels=0,lty=2,col=2,add=T)

    ## image(t(log.irr.vac[[3]]),breaks=my.breaks,zlim=zlims, col=my.col,axes=FALSE)
    ## abline(h=ve1.grids,col=AddAlpha("grey",.2))
    ## abline(v=week.grids,col=AddAlpha("grey",.2))
    ## contour(t(log.irr.vac[[3]]),levels=my.contour.levels,add=T)
    ## if(irr.all){
    ##     contour(t(log.irr.all[[3]]),levels=0,lty=2,col=2,add=T)
    ## }
    ## else {
    ##     contour(t(fs.diff[[3]]),levels=0,lty=2,col=2,add=T)
    ## }

    title(paste0(round(vacdoses[dosenums[3]]/.5e4,1),"/",round(vacdoses[dosenums[3]]/.5e4/2,1),"% Vaccinated")
         ,line=0.5)

    axis(2,labels=rep("",length(ve1.labels)),at=ve1.at
        ,col = "grey",cex.axis=cex.axis)

    axis(4,labels=ve1.labels,at=ve1.at,col = "grey",cex.axis=cex.axis)

    if(add.unvac){
        axis(1,labels=rep("",length(week.labels))
            ,at=week.ats,col="grey",cex.axis=0.8)
    } else{
        axis(1,labels=week.labels,at=week.ats,col="grey",cex.axis=cex.axis)
    }


    ########################
    ## now for the legend ##
    ########################

    par(oma=c( 0,0,0,1))# reset margin to be much smaller.
    set.panel() # reset plotting device

   # recover()
image.plot(log.arr.vac[[1]],
           legend.only=TRUE,
           zlim=zlims,
           breaks = my.breaks,
           axis.args=list(at=log(c(1/c(1,2,5,10,100),c(1,2,5,10,100))),labels=c(1/c(1,2,5,10,100),1,2,5,10,100),cex.axis=.5),
           col=my.col,bty='n',
           smallplot = c(.95,.97,.2,.8),
           legend.line=10)

    ## mtext(paste0(round(vacdoses[dosenums[1]]/.5e4,1),"/",round(vacdoses[dosenums[1]]/.5e4/2,1)),side=3,outer=TRUE,line=-2,at=.25,cex=.8)
    ## mtext(paste0(round(vacdoses[dosenums[2]]/.5e4,1),"/",round(vacdoses[dosenums[2]]/.5e4/2,1)),side=3,outer=TRUE,line=-2,at=.5,cex=.8)
    ## mtext(paste0(round(vacdoses[dosenums[3]]/.5e4,1),"/",round(vacdoses[dosenums[3]]/.5e4/2,1)),,side=3,outer=TRUE,line=-2,at=.75,cex=.8)
    mtext("Vaccination Start (Epidemic Week)",side=1,outer=TRUE,line=-1.4,at=.5,cex=.8)
    if(mrse){
        mtext('Relative Single Dose VE',side=2,outer=TRUE,line=-1.6,at=.5,cex=.8)
    }else{
        mtext('Single Dose VE',side=2,outer=TRUE,line=-1.6,at=.5,cex=.8)
    }
    if (add.unvac){
        mtext("Vaccinated",side=2,outer=TRUE,line=-1.2,at=.7,cex=.8)
        mtext("Unvaccinated",side=2,outer=TRUE,line=-1.2,at=.3,cex=.8)
    }

    if(pdf) dev.off()
}


##' Helper function to get final size of an epidemic
##' @param vac
##' @param ci.cols
##' @return
##' @author asa
get.final.size <- function(vac,ci.cols=c("CI","CI1","CI2")){
    tail(rowSums(vac[,ci.cols]),1)
}

##' Using this to calibrate betas to data for different models
##' used to calibrate main model to Guinea Bissau data
##' @param my.beta
##' @param data
##' @param days.to.aggregate
##' @param params
##' @param ode.func
##' @param initial.state
##' @param times
##' @param plot.it
##' @param model.type
##' @return
##' @author asa
obj.func <- function(my.beta,
                     data,
                     days.to.aggregate=1,
                     params,
                     ode.func=seir.ves.dx.dt,
                     initial.state=
                     c(0.5e6-11,0,0,0,0,0,11,0,0,0,0,0,0,0,0,0,0,0),
                     times=seq(0,300,by=1),
                     normalize=TRUE,
                     plot.it=TRUE,
                     model.type="ves"
                     ){

    if (length(my.beta) == 2){
        params$beta <- my.beta[1]
        params$epsilon <- my.beta[2]
    } else{
        params$beta <- my.beta
    }

    run <- ode(y=initial.state,
               times=times,
               func=ode.func,
               parms=params,
               vac.starts=c(1e6,1e6),
               vac.ends=c(1e6,1e6),
               daily.vac=0)


    if (model.type == "sirb"){
        run <- add.colnames(run,"sirb")
    } else if (model.type == "ves") {
        run <- add.colnames(run)
    } else {
        stop("don't know how to add column names to this model")
    }
    sim.curve <- diff(run[,c("CI")],1)

    if(days.to.aggregate > 1){
        zeros <- days.to.aggregate - length(sim.curve) %% days.to.aggregate
        if(zeros>0)
            sim.curve <- c(sim.curve,rep(0,zeros))
        sim.curve <- colSums(matrix(sim.curve,days.to.aggregate))
    }

    if(normalize){
        out.sim.curve <- sim.curve/max(sim.curve)
        out.data <- data/max(data)
    } else{
        out.sim.curve <- sim.curve
        out.data <- data
    }

    ## deal with differences in lenght
    curve.diffs <- length(sim.curve) - length(data)
    if (curve.diffs > 0){
        out.data <- c(out.data,rep(0,curve.diffs))
    } else if (curve.diffs < 0){
        out.sim.curve <- c(out.sim.curve,rep(0,curve.diffs))
    }

    if(plot.it){
        plot(out.data,type="h")#,col="grey")
        lines(out.sim.curve,col=2)
    }

    return(sum((out.data - out.sim.curve)^2))

}


##' gathers data from mclapply runs into usable form
##' @param data
##' @param vacdoses
##' @param ve1s
##' @param vac.timings
##' @param loop.args
##' @param alt.structure - when runs were done in parallel structure is slightly different
##' @return
##' @author asa
shape.mlapply.out <- function(data,
                              vacdoses=seq(5000,355000,by=10000),
                              ve1s= seq(0,0.75,by=.02),
                              vac.timings=seq(0,300,by=2),
                              loop.args=expand.grid(ve1=ve1s,
                                  vactime=vac.timings,
                                  doses=vacdoses),
                              alt.structure=FALSE){

    ## fs.mat.1dose[ve1,vt,dose,]  want array like this

    ## each run returns a list like this:
    ## rc <- list(fs.mat.1dose=fs.mat.1dose, 1
    ##            fs.mat.2dose=fs.mat.2dose, 2
    ##            time.at.risk.1dose = time.at.risk.1dose, 3
    ##            time.at.risk.2dose = time.at.risk.2dose, 4
    ##            doses=doses, 5
    ##            times=times, 6
    ##            uncon.run=uncon.run, 7
    ##            ve1=ve1 8
    ##            )

    ## loop.arg.locations
    ## data.locations

    fs.mat.1dose <- array(dim=c(
                              length(ve1s),length(vac.timings),length(vacdoses),4))
    fs.mat.2dose <- array(dim=c(
                              length(ve1s),length(vac.timings),length(vacdoses),4))

    ## for storing the time at risk for individuals in each vaccination state
    time.at.risk.1dose <- array(dim=c(
                                    length(ve1s),length(vac.timings),length(vacdoses),3))

    time.at.risk.2dose <- array(dim=c(
                                    length(ve1s),length(vac.timings),length(vacdoses),3))

    run.index <- 1
    for(d in seq_along(vacdoses)){
        for(t in seq_along(vac.timings)){
            for(ve in seq_along(ve1s)){
                if (alt.structure){
                    fs.mat.1dose[ve,t,d,] <- data[run.index][[1]]
                    fs.mat.2dose[ve,t,d,] <- data[run.index+1][[1]]
                    time.at.risk.1dose[ve,t,d,] <- data[run.index+2][[1]]
                    time.at.risk.2dose[ve,t,d,] <- data[run.index+3][[1]]
                    run.index <- run.index+8
                } else {
                    fs.mat.1dose[ve,t,d,] <- tail(data[run.index][[1]],1)
                    fs.mat.2dose[ve,t,d,] <- tail(data[run.index+1][[1]],1)
                    time.at.risk.1dose[ve,t,d,] <- tail(data[run.index+2][[1]],1)
                    time.at.risk.2dose[ve,t,d,] <- tail(data[run.index+3][[1]],1)
                    run.index <- run.index+8
                }
            }
        }
    }


    return(list(fs.mat.1dose=fs.mat.1dose,
                fs.mat.2dose=fs.mat.2dose,
                time.at.risk.1dose=time.at.risk.1dose,
                time.at.risk.2dose=time.at.risk.2dose))
}

##' Makes plot with min-ve for different vaccine mechanisms
##' @param uncon.run
##' @param plot.name
##' @param max.time
##' @param vac.starts
##' @param mydoses
##' @param min.ves
##' @param min.aon
##' @param min.leaky
##' @param single.one - if TRUE, then only plots the min.ves
##' @return
##' @author asa
make.epicurve.mech.plot <- function(uncon.run,
                                    plot.name,
                                    max.time = 200,
                                    vac.starts=seq(0,300,by=5),
                                    min.ves,
                                    min.aon,
                                    min.leaky,
                                    single.one=FALSE,
                                    ve.ylim=c(0.05,.5),
                                    epi.plot.type="h",
                                    dose.cols=c(1,2),
                                    mrse=FALSE,
                                    ve2=0.88
                                    ){

    if(!is.null(plot.name)){
        pdf(paste0("Plots/",plot.name,".pdf"),width=6,height=3)
    } else {
        quartz("",width=6,height=3)
    }

    par(mfrow=c(1,1),mar=c(0.5,3,.5,.5),oma=c(3,1,2,.5),mgp = c(2,.5,0))
    layout(mat=matrix(c(1,1,2,2,2,2,2,2),byrow=T,nrow=4))
    plot(diff(uncon.run[,c("CI")]),ylab="Cholera Cases",
         xlab="",axes=FALSE,ylim=c(0,.8),
         xlim=c(0,max.time),type=epi.plot.type,col="darkgrey")

    abline(v=seq(0,max.time,by=7),col=AddAlpha("grey",.1))
    abline(h=seq(0,max(diff(uncon.run[1:max.time,c("CI")])),length=10),
           col=AddAlpha("grey",.1))
    week.ats <- seq(0,max.time,by=7)[seq(1,length(seq(0,max.time,by=7)),by=2)]
    week.labs <- seq(0,max.time,by=7)[seq(1,length(seq(0,max.time,by=7)),by=2)]/7

    axis(3,
         labels=week.labs#seq(0,max.time,by=7)
         ,at=week.ats#seq(0,max.time,length=22)
         ,col = "grey",cex.axis=1)

    axis(2,
         col = "grey",

         cex.axis=1)

    plot(-10,-10,axes=F,xlim=c(0,max.time),ylim=ve.ylim,
         ylab=ifelse(mrse,"MRSE",
             "Minimum Single Dose VE"))

#    recover()

    if(single.one){
        mod.ves <- apply(min.ves,1,range)[,1:length(vac.starts)]/ifelse(mrse,ve2,1)
        polygon(c(vac.starts,rev(vac.starts)),c(mod.ves[dose.cols[1],],rev(mod.ves[dose.cols[2],])),col=AddAlpha(1,.4),border=FALSE)
    } else {

        mod.ves <- min.ves/ifelse(mrse,ve2,1)
        mod.aon <- min.aon/ifelse(mrse,ve2,1)
        mod.leaky <- min.leaky/ifelse(mrse,ve2,1)

        polygon(c(vac.starts,
                  rev(vac.starts)),
                c(mod.ves[,dose.cols[1]],rev(mod.ves[,dose.cols[2]])),
                col=AddAlpha(1,.4),border=FALSE)
        polygon(c(vac.starts,
                  rev(vac.starts)),
                c(mod.leaky[,dose.cols[1]],rev(mod.leaky[,dose.cols[2]])),
                col=AddAlpha(2,.4),border=FALSE)
        polygon(c(vac.starts,
                  rev(vac.starts)),
                c(mod.aon[,dose.cols[1]],rev(mod.aon[,dose.cols[2]])),
                col=AddAlpha(3,.4),border=FALSE)
    }

    axis(2,col="grey")
    axis(1,
         labels=week.labs#seq(0,max.time,by=7)
         ,at=week.ats#seq(0,max.time,length=22)
         ,col = "grey",cex.axis=1)
    abline(v=week.ats,col=AddAlpha("grey",.1))
    abline(h=seq(0,.5,by=0.05),col=AddAlpha("grey",.1))

    abline(h=.428/ifelse(mrse,0.87,1),lty=2,col=AddAlpha("black",.7))
    abline(h=.46/ifelse(mrse,.79,1),lty=2,col=AddAlpha("black",.7))
    text(par("usr")[2]*.8,.435/ifelse(mrse,0.87,1),"Single-Dose VE Estimate, Guinea 2014",cex=.6)
    text(par("usr")[2]*.8,.467/ifelse(mrse,.79,1),"Single-Dose VE Estimate, Zanzibar 2012",cex=.6)

    mtext("Epidemic Week",side=1,outer=T,line=1.5,cex=.8)

    if(!is.null(plot.name)) dev.off()
}


##' Gets latinized samples from random normals for
##' parameters (used in melding routine)
##' Note: this was not used in the paper
##' @param n number of samples
##' @param k numebr of dimenions
##' @param lowers vector of length k with lower uniform bounds on params
##' @param uppers vector of length k with upper uniform bounds on params
##' @param param.names string vector of parameter names (needs to match names in param vector in future steps)
##' @param model.type model type (we may use VES or leaky)
##' @param drop.Rlt1 TRUE if we want to drop all estimates that are
##' @param calc.R.func
get.lhs.samples <-function(n,k,
                           lowers,uppers,
                           param.names=c("beta","gamma"),
                           model.type="ves",
                           drop.Rlt1=TRUE,
                           calc.R.func=NULL){
    require(lhs)
    if (length(lowers) != length(uppers))
        stop("lowers and uppers need to be of same length")
    ## generate n random uniform samples for k dimensions
    x <- randomLHS(n=n, k=k)
    y <- x
    for (i in 1:length(lowers))
        y[,i] <- qunif(x[,i], lowers[i], uppers[i])

    colnames(y) <- param.names

    if (drop.Rlt1){
        Rs <- calc.R.func(y,model.type=model.type)
        if("s0" %in% param.names) Rs <- Rs * y[,"s0"]
        y <- y[-which(Rs < 1),] # drop those with R < 1
    }

    return(y)
}

##' Gets likelihood based resampling weights for melding routine
##' Note: this was not used in the paper
##' @param lhs.samps
##' @param n.samples
##' @param data
##' @param sd.errordist
##' @param runmodel
##' @param params
##' @param N
##' @return
##' @author asa
get.resample.weights <- function(lhs.samps,
                                 data,
                                 sd.errordist=15,
                                 runmodel=run.ves.model,
                                 params=list(
                                     beta=.6,
                                     gamma=1/2,
                                     sigma=1/1.4,
                                     phi1=0,
                                     phi2=0,
                                     theta1=0.4,
                                     theta2=.80,
                                     s0=1,
                                     obs=1
                                     ),
                                 N=0.5e6,
                                 plot.it=FALSE,
                                 aggregate.to="day"
                                 ){

    if(!aggregate.to %in% c("week","day")) stop("don't recognzie this agregation")

    max.time <- ifelse(aggregate.to == "week",max(data[,1])*7,max(data[,1]))
    times <- seq(0,max.time,by=1)

    ## match names of lhs.samps to param names
    samp.to.param.index <- match(colnames(lhs.samps),names(params))
    n <- nrow(lhs.samps)
    ll <- numeric(n)
    ## run model for each set of params
    for (i in 1:n){
        # update params
        if (i %% 500 == 0) cat(".")
        params[samp.to.param.index] <- lhs.samps[i,]

        # assuming a fixed-observed seeded number of infectious
        # this could be relaxed but will be highly correlated
        # with others
        seeded.infs <- data[1,2]/params$obs

        # generate initial state
        initial.state <- c((N*params$s0)-seeded.infs,0,0,
                           0,0,0,
                           seeded.infs,0,0,
                           N*(1-params$s0),0,0,
                           seeded.infs,0,0,
                           0,0,0)

        # run the model
        tmp <- runmodel(params,initial.state,times,
                        return.epi.curve=TRUE,aggregate.to=aggregate.to)

        # calculate the log-likelihood of the params given the residuals
        ll[i] <- sum(dnorm(tmp,mean=data[,2],sd=sd.errordist,log=T))

        if(plot.it){
            plot(data[,2],main=ll[i])
            lines(tmp,col="red")
        }

    }

    # now calculate the weights
    # note that some of these will be quite small
    #epsilon <- 10^-5 # desired precision
    ## See this post: http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
    tmp <- ll - max(ll)
    w <- exp(ll-max(ll))
    #w <- ifelse(tmp < log(epsilon) - log(n),0,w)
    w <- w/sum(w)
    #print(which(w > 0))

    return(w)

}


##' Runs SIER.VES model with an option to output just the epicurve
##' Note: this was not used in the paper
##' @param params
##' @param initial.state
##' @param times
##' @param vac.starts
##' @param vac.ends
##' @param daily.vac
##' @param return.epi.curve
##' @return
##' @author asa
run.my.model <- function(params,
                         initial.state,
                         times,
                         vac.starts=c(1e6,1e6),
                         vac.ends=c(1e6,1e6),
                         daily.vac=0,
                         aggregate.to="week",
                         model.name=seir.ves.dx.dt,
                         model.type="ves",
                         return.epi.curve=T){

    run <- ode(y=initial.state,
               times=times,
               func=model.name,
               parms=params,
               vac.starts=vac.starts,
               vac.ends=vac.ends,
               daily.vac=daily.vac)

    run <- add.colnames(run,model.type=model.type)

    if(return.epi.curve){
        # scale down epi-cirve by observation proportion
        ecurve <- params$obs*diff(rowSums(run[,c("CI","CI1","CI2")]))
        if(aggregate.to=="day"){
            return(ecurve)
        } else if (aggregate.to=="week"){
            ecurve <- tapply(ecurve,rep(1:ceiling(length(ecurve)/7),each=7)[1:length(ecurve)],sum)
            return(ecurve)
        } else {
            stop("Don't know how to aggregate the data based on aggregate.to")
        }
    } else {
        return(run)
        warning("This output is not scaled by the observation proportion like when the epicurve is returned")
    }
}

##' Runs vac sims with parameters from approximate posterior
##' Note: this was not used in the paper
##' @param n.draws
##' @param theta
##' @param w
##' @param runmodel
##' @param params
##' @param times
##' @param N
##' @param data
run.meld.vac.sims <- function(n.draws,
                              theta,
                              w,
                              runmodel,
                              params,
                              times,
                              N=.4e6,
                              data
                              ){

    ## draw the parameters
    samp.pars <- sample(nrow(theta),n.draws,replace = TRUE,prob=w)
    epi.curves <- array(dim=c(n.draws,length(times)-1))

    ## match names of lhs.samps to param names
    samp.to.param.index <- match(colnames(theta),names(params))
    if (any(is.na(samp.to.param.index))) stop("cannot find the correct labels in the param list to coresspond with the those in the the theta array")
    for (i in seq_along(samp.pars)){
        params[samp.to.param.index] <- theta[samp.pars[i],]
        ## print(params[samp.to.param.index])
        seeded.infs <- data[1,2]/params$obs

                                        #recover()
        initial.state <- c((N*params$s0)-seeded.infs,0,0,
                           0,0,0,
                           seeded.infs,0,0,
                           N*(1-params$s0),0,0,
                           seeded.infs,0,0,
                           0,0,0)

        epi.curves[i,] <- runmodel(params,initial.state,times,return.epi.curve=TRUE)
    }

    return(epi.curves)
}


##' Calculates simulated log likelihood
##' using normal errors
##' Note: not used in paper
##' @title
##' @param param.vec
##' @param data
##' @param sd.errordist
##' @param runmodel
##' @param params
##' @param N
##' @param log.params
##' @return
##' @author asa
sim.log.lik <- function(param.vec,
                        data,
                        sd.errordist=15,
                        runmodel=run.ves.model,
                        params=list(
                            beta=.6,
                            gamma=1/2,
                            sigma=1/1.4,
                            phi1=0,
                            phi2=0,
                            theta1=0.4,
                            theta2=.80,
                            s0=1,
                            obs=1
                            ),
                        N=0.5e6,
                        log.params=TRUE
                        ){

    times <- seq(0,max(data[,1]),by=1)
    seeded.infs <- data[1,2]
    if (log.params){
        params$beta <- exp(param.vec[1])
        params$gamma <- exp(param.vec[2])
        params$s0 <- plogis(param.vec[3])
    } else {
        params$beta <- param.vec[1]
        params$gamma <- param.vec[2]
        params$s0 <- param.vec[3]
    }

    initial.state <- c((N*params$s0)-seeded.infs,0,0,
                       0,0,0,
                       seeded.infs,0,0,
                       N*(1-params$s0),0,0,
                       0,0,0,
                       0,0,0)

    tmp <- runmodel(params,initial.state,times,return.epi.curve=TRUE)
    ll <- sum(dnorm(tmp- data[,2],mean=0,sd=sd.errordist,log=T))
    return(ll)
}


##' Helper function to help determine if epidemic has ended
##' @param run named model run
##' @param inc.cutoff threshold value for deciding if epidemic is over
##' @return true of incidence is below cutoff
##' @author asa
check.epi.ended <- function(run,inc.cutoff){
    ifelse(
        tail(diff(rowSums(run[,c("CI","CI1","CI2")])),1) < inc.cutoff,
        TRUE,
        FALSE)
}


##' Simulates vaccination in zimbabwe
##' Note: not used in paper
##' @param n.sims
##' @param sim.start.week
##' @param data
##' @param initial.state
##' @param param.samps
##' @param param.weights
##' @param vac.start
##' @param main.paramsvac.starts
##' @param vac.ends
##' @param daily.vac
##' @param max.time
##' @param N
##' @param dx.dt.func
##' @return
##' @author asa
simulate.zim.vac <- function(n.sims,
                             sim.start.week,
                             data,
                             param.samps,
                             param.weights,
                             main.params,
                             doses=13e6*.2,
                             max.time=350,
                             total.vac.days=10,
                             inter.dose.time=14,
                             N=13e6,
                             dx.dt.func=seir.ves.dx.dt
                             ){

    time.step <- 1
    times <- seq(sim.start.week*7,max.time,time.step)

    daily.vac.doses <- doses/total.vac.days

    ## set up timings for both scenarios
    one.of.one.start <- sim.start.week*7
    one.of.one.end <- one.of.one.start+total.vac.days
    one.of.two.start <- sim.start.week*7
    one.of.two.end <- one.of.two.start+total.vac.days/2
    two.of.two.start <- one.of.two.end+inter.dose.time
    two.of.two.end <- two.of.two.start + total.vac.days/2


    ## adjust initial state to account for reduction in susptibles at the simulatino start time
    n.cum.cases <- sum(data[,2][1:(sim.start.week-1)])/params$obs
    n.infectious.cases <- (data[sim.start.week,2]/7)/params$obs
    initial.state <- c((N*params$s0)-n.cum.cases,0,0,
                       0,0,0,
                       n.infectious.cases,0,0,
                       N-(N*params$s0)+n.cum.cases,0,0,
                       0,0,0,
                       n.cum.cases,0,0)

    # sample paramters
    thetas <- sample(nrow(param.samps),size=n.sims,replace=T,prob=param.weights)
    dvac <- svac <- nvac <- vector("list",length=n.sims)
    for (i in 1:n.sims){
        params$beta <- my.lhs[thetas[i],1]
        params$gamma <- my.lhs[thetas[i],2]

        run.non.vac <- ode(y=initial.state,
                           times=times,
                           func=dx.dt.func,
                           parms=params,
                           vac.starts=c(Inf,Inf),
                           vac.ends=c(Inf,Inf),
                           daily.vac=daily.vac.doses)

        run.non.vac <- add.colnames(run.non.vac,model.type="ves")
        long.enough <- check.epi.ended(run=run.non.vac,inc.cutoff=.1)

        while(!long.enough){
##            print("extending time - non-vac simulation")
            times <- seq(sim.start.week*7,max(times)+50,time.step)

            run.non.vac <- ode(y=initial.state,
                               times=times,
                               func=dx.dt.func,
                               parms=params,
                               vac.starts=c(Inf,Inf),
                               vac.ends=c(Inf,Inf),
                               daily.vac=daily.vac.doses)


            run.non.vac <- add.colnames(run.non.vac,model.type="ves")
            long.enough <- check.epi.ended(run=run.non.vac,inc.cutoff=.1)
        }



        run.single.vac <- ode(y=initial.state,
                              times=times,
                              func=dx.dt.func,
                              parms=params,
                              vac.starts=c(one.of.one.start,Inf),
                              vac.ends=c(one.of.one.end,Inf),
                              daily.vac=daily.vac.doses)

        run.single.vac <- add.colnames(run.single.vac,model.type="ves")
        long.enough <- check.epi.ended(run=run.single.vac,inc.cutoff=.1)

        while(!long.enough){
##            print("extending time - two-dose sim")
            times <- seq(sim.start.week*7,max(times)+50,time.step)

            run.single.vac <- ode(y=initial.state,
                                  times=times,
                                  func=dx.dt.func,
                                  parms=params,
                                  vac.starts=c(one.of.one.start,Inf),
                                  vac.ends=c(one.of.one.end,Inf),
                                  daily.vac=daily.vac.doses)


            run.single.vac <- add.colnames(run.single.vac,model.type="ves")
            long.enough <- check.epi.ended(run=run.single.vac,inc.cutoff=.1)
        }

        run.double.vac <- ode(y=initial.state,
                              times=times,
                              func=dx.dt.func,
                              parms=params,
                              vac.starts=c(one.of.two.start,two.of.two.start),
                              vac.ends=c(one.of.two.end,two.of.two.end),
                              daily.vac=daily.vac.doses)

        run.double.vac <- add.colnames(run.double.vac,model.type="ves")


        ## check here
        long.enough <- check.epi.ended(run=run.double.vac,inc.cutoff=.1)

        while(!long.enough){
##            print("extending time - 2-dose sim")
            times <- seq(sim.start.week*7,max(times)+50,time.step)

            run.double.vac <- ode(y=initial.state,
                                  times=times,
                                  func=dx.dt.func,
                                  parms=params,
                                  vac.starts=c(one.of.two.start,two.of.two.start),
                                  vac.ends=c(one.of.two.end,two.of.two.end),
                                  daily.vac=daily.vac.doses)

            run.double.vac <- add.colnames(run.double.vac,model.type="ves")

            long.enough <- check.epi.ended(run=run.double.vac,inc.cutoff=.1)
        }

    dvac[[i]] <- params$obs*rowSums(run.double.vac[,c("CI","CI1","CI2")])
    svac[[i]] <- params$obs*rowSums(run.single.vac[,c("CI","CI1","CI2")])
    nvac[[i]] <- params$obs*rowSums(run.non.vac[,c("CI","CI1","CI2")])

    }

    return(list(double=dvac,single=svac,unvac=nvac))

}

##' helper function to evaluate(parse(text for a vector
##' Note: not used in paper
##' @param x - vector of strings that correspond to variable names in
##' memory
##' @return
##' @author asa
evpc <- function(x,...){
    sapply(x,function(y) eval(parse(text=y)))
}

##' Adds alpha to a set of colors
##' NOTE: I took this code from somewhere online but don't rememebr where
##' a long time ago
##' @title
##' @param COLORS
##' @param ALPHA
##' @return
AddAlpha <- function(COLORS, ALPHA){
    if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
    RGB <- col2rgb(COLORS, alpha=TRUE)
    RGB[4,] <- round(RGB[4,]*ALPHA)
    NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
    return(NEW.COLORS)
}

