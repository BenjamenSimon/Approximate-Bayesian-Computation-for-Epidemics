


#################################################################################################

GastroData <- c(89, 28)

GastroDetectionTimes <- c(0, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7)

GastroDataAlt <- as.matrix(summary(as.factor(GastroDetectionTimes)))
GastroTemporalData <- matrix(nrow = 13, ncol = 4)
colnames(GastroTemporalData) = c('time', 'S', 'I', 'R')
GastroTemporalData[,1] <- seq(0, 12, b=1)
GastroTemporalData[,4] <- 0
GastroTemporalData[(as.integer(rownames(GastroDataAlt))+1), 4] <- as.integer(GastroDataAlt)
GastroTemporalData[,4] <- cumsum(GastroTemporalData[,4])
GastroTemporalData[8,3] <- 0


GastroSummary <- c(1, 4, 2, 3, 3, 10, 5, 0, 7)


#################################################################################################


ABCSummaryGastro <- function(res){ #Takes the results of the simulation as an input
  
  Bin01 <- (res[2,4])
  Bin12 <- (res[3,4] - res[2, 4])
  Bin23 <- (res[4,4] - res[3, 4])
  Bin34 <- (res[5,4] - res[4, 4])
  Bin45 <- (res[6,4] - res[5, 4])
  Bin56 <- (res[7,4] - res[6, 4])
  Bin67 <- (res[8,4] - res[7, 4])
  Bin7inf <- (res[13,4] - res[8, 4]) # About balancing simulation time and confidence that its finished.
  
  Duration <- min(which(res[,3] == 0)) - 1 # The time at which the last removal occurs
  
  return(c(Bin01, Bin12, Bin23, Bin34, Bin45, Bin56, Bin67, Bin7inf, Duration))
}


#################################################################################################


ABCdistGastro <- function(Tdata, Tsim){ #Takes the summary statistics of the data and simulation as inputs
  
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


GastroABCSIR = function(N=89, beta, gamma, t.max = 12, t.step = 1){ #t.max = T, t.step = 1 day by default
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

PosteriorSamplesABCGastro = matrix(nrow = NSim, ncol = 3) # A matrix to store the results

sin=proc.time()

set.seed(19)

Tdata <- GastroSummary

suppressWarnings(
  
  for(i in 1:NSim)
  {
    betadraw = rexp(1,1)
    gammadraw = rexp(1,1)
    
    stochsim = GastroABCSIR(N=89, beta= betadraw, gamma = gammadraw, t.max = 12, t.step = 1)
    
    Tsim <- ABCSummaryGastro(stochsim)
    
    PosteriorSamplesABCGastro[i, ] <- c(betadraw, gammadraw, ABCdistGastro(Tdata, Tsim))
    
  }
  
)

sout=proc.time()
sout-sin # 347.38


#################################################################################################

hist(PosteriorSamplesABCGastro[,3], breaks = 200, xlim = c(3.5,10), ylim = c(0, 200))

TolGastro <- c(8) # The tolerance for single distance metric

sum(PosteriorSamplesABCGastro[,3] <= TolGastro[1])


#################################################################################################


Theta_ABC_Gastro = PosteriorSamplesABCGastro[which(PosteriorSamplesABCGastro[,3] <= TolGastro[1]), c(1,2)]
summary(Theta_ABC_Gastro[,1]) #Beta
summary(Theta_ABC_Gastro[,2]) #Gamma

hist(Theta_ABC_Gastro[,1], main = "Beta")
boxplot(Theta_ABC_Gastro[,1], outline=FALSE, main = "Beta") #Not plotting outliers

hist(Theta_ABC_Gastro[,2], main = "Gamma") 
boxplot(Theta_ABC_Gastro[,2], outline=FALSE, main = "Gamma") #Not plotting outliers


#################################################################################################
#################################################################################################


PostSamplesABC = matrix(nrow = 1000, ncol = 3) # A matrix to store the results

sin=proc.time()

set.seed(109)

Tdata <- GastroSummary

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
    
    stochsim = GastroABCSIR(N=89, beta= betadraw, gamma = gammadraw, t.max = 12, t.step = 1)
    
    Tsim <- ABCSummaryGastro(stochsim)
    
    dist <- ABCdistGastro(Tdata, Tsim)
    
    if(dist <= 4.5){
      PostSamplesABC[a, ] <- c(betadraw, gammadraw, dist)
      a = a+1
    }else{next}
  }
  
)

sout=proc.time()
sout-sin #30.14
sims #86999

summary(PostSamplesABC)
sd(PostSamplesABC[,1])
sd(PostSamplesABC[,2])



R0 <- (89*PostSamplesABC[,1]) / (1-exp(-PostSamplesABC[,2]))
mean(R0)
median(R0)



