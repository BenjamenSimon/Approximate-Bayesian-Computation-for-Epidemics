

AbakalikiData <- c(120, 30)

AbakalikiDetectionTimes <- c(0, 13, 20, 22, 25, 25, 25, 26, 30, 35, 38, 40, 40, 42, 42, 47, 50, 51, 55, 55, 56, 57, 58, 60, 60, 61, 66, 66, 71, 76)

AbakalikiDataAlt <- as.matrix(summary(as.factor(AbakalikiDetectionTimes)))
AbakalikiTemporalData <- matrix(nrow = 101, ncol = 4)
colnames(AbakalikiTemporalData) = c('time', 'S', 'I', 'R')
AbakalikiTemporalData[,1] <- seq(0, 100, b=1)
AbakalikiTemporalData[,4] <- 0
AbakalikiTemporalData[(as.integer(rownames(AbakalikiDataAlt))+1), 4] <- as.integer(AbakalikiDataAlt)
AbakalikiTemporalData[,4] <- cumsum(AbakalikiTemporalData[,4])
AbakalikiTemporalData[77,3] <- 0

AbakalikiSummary <- c(2, 6, 3, 7, 8, 4, 0, 76)


#################################################################################################



AbakalikiSIR = function(N, beta, gamma, t.max = 76, t.step = 1){ #t.max = T, t.step = 1 day by default
  S = N-1; I = 1; R = 0 #R=1
  T = seq(0, t.max, b=t.step)
  res = matrix(nrow = length(T), ncol = 4)
  colnames(res) = c('time', 'S', 'I', 'R')
  res[1,] = c(0,S,I,R)
  # If this is the first row of the matrix it means the first R draw
  # was too large.
  
  day = 1
  
  while(day <= length(T)){
    pInfec = 1 - exp(-beta*I*t.step)
    Inew = rbinom(1, size = S, prob = pInfec)
    
    
    pRemove = 1 - exp(-gamma*t.step)
    Rnew = rbinom(1, size = I, prob = pRemove)
    
    S = S - Inew
    I = I + Inew - Rnew
    R = R + Rnew
    
    if(R > 0){
      if( R != AbakalikiTemporalData[(day), 4] ){break}
      #Need it in this order to otherwise puts in last step without checking it
      #So could have complete matrices that don't match the data
      res[day,] = c(T[day], S, I, R)
      
      day = day + 1
    }
  }
  
  res # Outputs a matrix that describes the population of each state at each timestep
}


#################################################################################################


ABCSummaryAbakaliki <- function(res){ #Takes the results of the simulation as an input
  
  Bin013 <- (res[14,4])
  Bin1426 <- (res[27,4] - res[14, 4])
  Bin2739 <- (res[40,4] - res[27, 4])
  Bin4052 <- (res[53,4] - res[40, 4])
  Bin5365 <- (res[66,4] - res[53, 4])
  Bin6678 <- (res[79,4] - res[66, 4])
  Bin79inf <- (res[91,4] - res[79, 4]) # About balancing simulation time and confidence that its finished.
  
  Duration <- min(which(res[,3] == 0)) - 1 # The time at which the last removal occurs
  
  return(c(Bin013, Bin1426, Bin2739, Bin4052, Bin5365, Bin6678, Bin79inf, Duration))
}

ABCSummaryAbakaliki(AbakalikiTemporalData)


#################################################################################################


ABCdistAbakaliki <- function(Tdata, Tsim){ #Takes the summary statistics of the data and simulation as inputs
  
  p <- length(Tdata)
  
  d = rep(NA, p)
  
  for(i in 1:(p-1)){
    d[i] <- dist(c(Tdata[i], Tsim[i]), method = "euclidean")
  }
  
  d[p] = dist(c( (Tdata[p]/50) , (Tsim[p]/50) ), method = "euclidean")
  
  TotalDist <- (sum(d))^(0.5)
  
  return(TotalDist)
}


#################################################################################################


AbakalikiABCSIR = function(N=120, beta, gamma, t.max = 90, t.step = 1){ #t.max = T, t.step = 1 day by default
  S = N-1; I = 1; R = 0 #R=1
  T = seq(0, t.max, b=t.step)
  res = matrix(nrow = length(T), ncol = 4)
  colnames(res) = c('time', 'S', 'I', 'R')
  res[1,] = c(0,S,I,R)
  
  day = 1
  
  while(day <= length(T)){
    pInfec = 1 - exp(-beta*I*t.step)
    Inew = rbinom(1, size = S, prob = pInfec)
    
    
    pRemove = 1 - exp(-gamma*t.step)
    Rnew = rbinom(1, size = I, prob = pRemove)
    
    S = S - Inew
    I = I + Inew - Rnew
    R = R + Rnew
    
    if(R > 0){
      res[day,] = c(T[day], S, I, R)
      
      day = day + 1
    }
  }
  
  res # Outputs a matrix that describes the population of each state at each timestep
}


#################################################################################################


NSim = 1000000 #Number of simulations

PosteriorSamplesABCAbakaliki = matrix(nrow = NSim, ncol = 3) # A matrix to store the results

sin=proc.time()

set.seed(19)

Tdata <- AbakalikiSummary

suppressWarnings(
  
  for(i in 1:NSim)
  {
    betadraw = rexp(1,1)
    gammadraw = rexp(1,1)
    
    stochsim = AbakalikiABCSIR(N=120, beta= betadraw, gamma = gammadraw, t.max = 100, t.step = 1)
    
    Tsim <- ABCSummaryAbakaliki(stochsim)
    
    PosteriorSamplesABCAbakaliki[i, ] <- c(betadraw, gammadraw, ABCdistAbakaliki(Tdata, Tsim))
    
  }
  
)

sout=proc.time()
sout-sin


#################################################################################################


hist(PosteriorSamplesABCAbakaliki[,3], ylim = c(0, 3000), breaks = 100)

TolAbakaliki <- c(5.428) # The tolerance for single distance metric

sum(PosteriorSamplesABCAbakaliki[,3] <= TolAbakaliki[1])


#################################################################################################


Theta_ABC_Abakaliki = PosteriorSamplesABCAbakaliki[which(PosteriorSamplesABCAbakaliki[,3] <= TolAbakaliki[1]), c(1,2)]
summary(Theta_ABC_Abakaliki)
sd(Theta_ABC_Abakaliki[,1])
sd(Theta_ABC_Abakaliki[,2])

hist(Theta_ABC_Abakaliki[,1], main = "Beta") 
summary(Theta_ABC_Abakaliki[,1]) #Beta
boxplot(Theta_ABC_Abakaliki[,1], outline=FALSE, main = "Beta") #Not plotting outliers
sd(Theta_ABC_Abakaliki[,1]) #Beta

hist(Theta_ABC_Abakaliki[,2], main = "Gamma") 
summary(Theta_ABC_Abakaliki[,2]) #Gamma
boxplot(Theta_ABC_Abakaliki[,2], outline=FALSE, main = "Gamma") #Not plotting outliers
sd(Theta_ABC_Abakaliki[,2]) #Gamma


# want $\beta$ to be around $0.00168 +- 0.00047$ and $\gamma$ to be around $0.162 +- 0.050$.




#################################################################################################
#################################################################################################


PostSamplesABC = matrix(nrow = 1000, ncol = 3) # A matrix to store the results

sin=proc.time()

set.seed(109)

Tdata <- AbakalikiSummary

Accepted = 1000

########

sims = 0
a=0

suppressWarnings(
  
  while(a <= Accepted)
  {
    sims = sims + 1
    
    betadraw = rexp(1,1)
    gammadraw = rexp(1,1)
    
    stochsim = AbakalikiABCSIR(N=120, beta= betadraw, gamma = gammadraw, t.max = 100, t.step = 1)
    
    Tsim <- ABCSummaryAbakaliki(stochsim)
    
    dist <- ABCdistAbakaliki(Tdata, Tsim)
    
    if(dist <= 5.428){
      PostSamplesABC[a, ] <- c(betadraw, gammadraw, dist)
      a = a+1
    }else{next}
  }
  
)

sout=proc.time()
sout-sin

summary(PostSamplesABC)
sd(PostSamplesABC[,1])
sd(PostSamplesABC[,2])



R0 <- (120*PostSamplesABC[,1]) / (1-exp(-PostSamplesABC[,2]))
mean(R0)
median(R0)































