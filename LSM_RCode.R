###############################################################
###############################################################
### ------------------------------------------------------ ####
### ---------------- R Code zur Simulation --------------- ####
### ------------------ und LSM Preisung ------------------ ####
### ------------------------------------------------------ ####


library(dplyr)
library(ggplot2)
library(stargazer)
library(somebm)


simMatrix <- function(S0=100, Nsimulations=100, r=0.01, sigma=0.2, se=123){
  # simMatrix is a function, simulating a Matrix of stock
  # prices for one year. The number of simulations can be
  # set by the parameter Nsimulations. 
  # The function simulates the Stock prices using the 
  # Black-Scholes Model. The Stock price for t>0 is given
  # by: S_t = S_0 * exp{r-sigma^2/2)*t + sigma*W_t}, where W_t is 
  # a Brownian Motion with mu=0 and variance is the time increment.
  # Nsimulations ... Number of simulations (int)
  # S0 ............. Starting value (int)
  # r .............. risk free interest rate (int)
  # sigma .......... the parameter sigma for the Black Scholes
  #                  simulation. 
  # se ............. set seed (int)
  
  set.seed(se)
  M <- matrix(NA, nrow=Nsimulations, ncol=252)
  N <- 252
  for (i in 1:Nsimulations){
    
    # Brownian Motion
    m <- 0
    dW <- rnorm(n=N, mean=m, sd=sqrt(1/252))
    W <- c(0, cumsum(dW))
    
    t <- 1:N
    S <- S0*exp((r-(sigma^2)/2)*(t/252)+sigma*W[t])
    M[i,] <- S
  }
  colnames(M) <- 1:252
  rownames(M) <- paste(1:Nsimulations)
  return(M)
}


M2 <- simMatrix(S0=1, Nsimulations=10, se=1234)
matplot(t(M2), type="l", lty=1, ylab="price underlying", xlab="t")




###############################################################
### ---------------- Martingaleigenschaft ---------------- ####
### ------------------------------------------------------ ####

M3 <- simMatrix(S0=1, Nsimulations=1000, se=1234)
# wir bekommen die mittleren stock preise via
M3means <- colMeans(M3)
# der jährliche zinssatz beträgt 1%
r <- 0.01
# wir berechnen den täglichen, also unterjährigen Abzinsungsfaktor
AZFu <- exp(log(1+0.01)/252)
# Bemerkung: rday=AZFu-1
ru <- AZFu-1
EW <- M3means/((AZFu)^(0:(length(M3means)-1)))
M3means <- data.frame(t=1:252, means=EW)


ggplot(M3means, aes(x=t, y=M3means$means)) +
  geom_point(shape=1) +
  ylim(0.75,1.25) +
  ylab("Arithm. Mittel aller 1000 Trajektorien")


###############################################################
### --------------------- LSM Preisung ------------------- ####
### ------------------------------------------------------ ####


LSM <- function(M, K=1.1, r=0.06){
  # This function performs the Least Squares Monte Carlo
  # algorithm for any simulated matrix M. The output
  # is the mean of the discounted cash flow of returns 
  # where the option was exercised.
  # M ... simulated matrix
  # K ... strike price
  # r ... interest rate
  T <- ncol(M)
  t <- T-1
  N <- nrow(M)
  Ai <- matrix(0, nrow=N, ncol=T)
  lapply(M[,t+1],function(x) max(K-x,0))  %>% unlist -> Ai[,t+1]
  
  for (t in (T-1):2){
    lapply(M[,t],function(x) max(K-x,0))  %>% unlist %>% 
      Filter(f=function(x) x!=0)-> im
    names(im) %>% as.numeric() -> ind
    # discounting
    Y <- rep(0,length(ind))
    for (i in 1:(T-t)){
      Y <- Y+exp(-r*i)*Ai[ind,t+i]
    }
    X <- im
    lm1 <- lm(Y~X+(X^2))
    f1 <- lm1$fitted.values
    # exercise vs keep
    df <- data.frame(beibehalten=f1, ausueben=im)
    ind2 <- apply(df,1, function(x) ifelse(x[2]>x[1], x[2], 0))
    Ai[ind,t] <- ind2
    # set all following to 0
    Ai[Ai[,t]!=0,(t+1):T] <- 0
  }
  
  Ai <- Ai[,2:T]
  Y <- rep(0,N)
  for (i in 1:(T-1)){
    Y <- Y+exp(-r*i)*Ai[,i]
  }
  return(mean(Y))
}


LSM0 <- function(M, K=1.1, r=0.06){
  # This function performs the Least Squares Monte Carlo
  # algorithm for any simulated matrix M. The output
  # is the mean of the discounted cash flow of returns 
  # where the option was exercised.
  # This function uses a mean model for prediction
  # M ... simulated matrix
  # K ... strike price
  # r ... interest rate
  T <- ncol(M)
  t <- T-1
  N <- nrow(M)
  Ai <- matrix(0, nrow=N, ncol=T)
  lapply(M[,t+1],function(x) max(K-x,0))  %>% unlist -> Ai[,t+1]
  
  for (t in (T-1):2){
    lapply(M[,t],function(x) max(K-x,0))  %>% unlist %>% Filter(f=function(x) x!=0)-> im
    names(im) %>% as.numeric() -> ind
    Y <- rep(0,length(ind))
    for (i in 1:(T-t)){
      Y <- Y+exp(-r*i)*Ai[ind,t+i]
    }
    X <- im
    f1 <- rep(mean(Y), length(im))
    # exercise vs keep
    df <- data.frame(beibehalten=f1, ausueben=im)
    ind2 <- apply(df,1, function(x) ifelse(x[2]>x[1], x[2], 0))
    Ai[ind,t] <- ind2
    # set all following to 0
    Ai[Ai[,t]!=0,(t+1):T] <- 0
  }
  
  Ai <- Ai[,2:T]
  Y <- rep(0,N)
  for (i in 1:(T-1)){
    Y <- Y+exp(-r*i)*Ai[,i]
  }
  return(mean(Y))
}


LSM1 <- function(M, K=1.1, r=0.06){
  # This function performs the Least Squares Monte Carlo
  # algorithm for any simulated matrix M. The output
  # is the mean of the discounted cash flow of returns 
  # where the option was exercised.
  # This function uses a single linear model for prediction
  # M ... simulated matrix
  # K ... strike price
  # r ... interest rate
  T <- ncol(M)
  t <- T-1
  N <- nrow(M)
  Ai <- matrix(0, nrow=N, ncol=T)
  lapply(M[,t+1],function(x) max(K-x,0))  %>% unlist -> Ai[,t+1]
  
  for (t in (T-1):2){
    lapply(M[,t],function(x) max(K-x,0))  %>% unlist %>% Filter(f=function(x) x!=0)-> im
    names(im) %>% as.numeric() -> ind
    Y <- rep(0,length(ind))
    for (i in 1:(T-t)){
      Y <- Y+exp(-r*i)*Ai[ind,t+i]
    }
    X <- im
    lm1 <- lm(Y~X)
    f1 <- lm1$fitted.values
    # exercise vs keep
    df <- data.frame(beibehalten=f1, ausueben=im)
    ind2 <- apply(df,1, function(x) ifelse(x[2]>x[1], x[2], 0))
    Ai[ind,t] <- ind2
    # set all following to 0
    Ai[Ai[,t]!=0,(t+1):T] <- 0
  }
  
  Ai <- Ai[,2:T]
  Y <- rep(0,N)
  for (i in 1:(T-1)){
    Y <- Y+exp(-r*i)*Ai[,i]
  }
  return(mean(Y))
}


LSM3 <- function(M, K=1.1, r=0.06){
  # This function performs the Least Squares Monte Carlo
  # algorithm for any simulated matrix M. The output
  # is the mean of the discounted cash flow of returns 
  # where the option was exercised.
  # This function uses a polynom of degree 3 for prediction
  # M ... simulated matrix
  # K ... strike price
  # r ... interest rate
  T <- ncol(M)
  t <- T-1
  N <- nrow(M)
  Ai <- matrix(0, nrow=N, ncol=T)
  lapply(M[,t+1],function(x) max(K-x,0))  %>% unlist -> Ai[,t+1]
  
  for (t in (T-1):2){
    lapply(M[,t],function(x) max(K-x,0))  %>% unlist %>% Filter(f=function(x) x!=0)-> im
    names(im) %>% as.numeric() -> ind
    Y <- rep(0,length(ind))
    for (i in 1:(T-t)){
      Y <- Y+exp(-r*i)*Ai[ind,t+i]
    }
    X <- im
    lm1 <- lm(Y~X+X^2+X^3)
    f1 <- lm1$fitted.values
    # exercise vs keep
    df <- data.frame(beibehalten=f1, ausueben=im)
    ind2 <- apply(df,1, function(x) ifelse(x[2]>x[1], x[2], 0))
    Ai[ind,t] <- ind2
    # set all following to 0
    Ai[Ai[,t]!=0,(t+1):T] <- 0
  }
  
  Ai <- Ai[,2:T]
  Y <- rep(0,N)
  for (i in 1:(T-1)){
    Y <- Y+exp(-r*i)*Ai[,i]
  }
  return(mean(Y))
}


LSM10 <- function(M, K=1.1, r=0.06){
  # This function performs the Least Squares Monte Carlo
  # algorithm for any simulated matrix M. The output
  # is the mean of the discounted cash flow of returns 
  # where the option was exercised.
  # This function uses a polynom of degree 10 for prediction
  # M ... simulated matrix
  # K ... strike price
  # r ... interest rate
  T <- ncol(M)
  t <- T-1
  N <- nrow(M)
  Ai <- matrix(0, nrow=N, ncol=T)
  lapply(M[,t+1],function(x) max(K-x,0))  %>% unlist -> Ai[,t+1]
  
  for (t in (T-1):2){
    lapply(M[,t],function(x) max(K-x,0))  %>% unlist %>% Filter(f=function(x) x!=0)-> im
    names(im) %>% as.numeric() -> ind
    Y <- rep(0,length(ind))
    for (i in 1:(T-t)){
      Y <- Y+exp(-r*i)*Ai[ind,t+i]
    }
    X <- im
    lm1 <- lm(Y~X+X^2+X^3+X^4+X^5+X^6+X^7+X^8+X^9+X^10)
    f1 <- lm1$fitted.values
    # exercise vs keep
    df <- data.frame(beibehalten=f1, ausueben=im)
    ind2 <- apply(df,1, function(x) ifelse(x[2]>x[1], x[2], 0))
    Ai[ind,t] <- ind2
    # set all following to 0
    Ai[Ai[,t]!=0,(t+1):T] <- 0
  }
  
  Ai <- Ai[,2:T]
  Y <- rep(0,N)
  for (i in 1:(T-1)){
    Y <- Y+exp(-r*i)*Ai[,i]
  }
  return(mean(Y))
}

GE0 <- LSM0(M=M3, r=ru)
GE1 <- LSM1(M=M3, r=ru)
GE <- LSM(M=M3, r=ru)
GE3 <- LSM3(M=M3, r=ru)
GE10 <- LSM10(M=M3, r=ru)


DF <- data.frame(Modell=c("Mittelwertmodell","einfach","Polynom 2", "Polynom 3", "Polynom 10"), Optionspreis=c(GE0,GE1,GE,GE3,GE10))

DF
