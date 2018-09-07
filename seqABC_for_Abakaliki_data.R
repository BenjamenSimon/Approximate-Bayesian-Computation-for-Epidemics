

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


importanceDensity <- function(particles, weights, covar){
  
  thetastar <- sample(nrow(particles), 1, prob = weights)
  thetastar <- particles[thetastar, ]
  
  theta <- rmvn(1, mu = thetastar, sigma = (2 * covar))
  
  return(theta)
  
}


#################################################################################################


WeightFunction <- function(N, theta, particles, prevwt, covar){
  
  sum = 0
  for(i in 1:N){
    sum = sum + prevwt[i] * dmvn(theta, mu = particles[i, ], sigma = (2 * covar))
  }
  sumwt <- sum(prevwt)
  
  weights <- (dexp(theta[1], 1) * dexp(theta[2], 1))/ (sum/sumwt)
  
  return(weights)
  
}


#################################################################################################


library(Hmisc)
library(mvnfast)

AbakalikiSummary <- c(2, 6, 3, 7, 8, 4, 0, 76)

SeqABC <- function(Acc, tol){
  
  PostSample <- matrix(nrow = (Acc * length(tol)), ncol = 4)
  colnames(PostSample) <- c("t", "beta", "gamma", "weight")
  
  sim = 0
  
  t = 1
  
  Tdata <- AbakalikiSummary
  
  while(t <= length(tol)){
    
    print(t)
    
    if(t == 1){
      
      a = 1
      
      while(a <= Acc){
        
        sim = sim + 1
        
        betadraw <- rexp(1)
        gammadraw <- rexp(1)
        
        stochsim <- AbakalikiABCSIR(N=120, beta= betadraw, gamma = gammadraw, t.max = 90, t.step = 1)
        
        suppressWarnings(
          Tsim <- ABCSummaryAbakaliki(stochsim)
        )
        
        if(ABCdistAbakaliki(Tdata, Tsim) > tol[t]){next}else{
          
          wts <- 1
          
          PostSample[(((t-1)*Acc)+a), ] <- c(t, betadraw, gammadraw, wts)
          
          a = a + 1
          
        }
      }
      
      t = t+1
      
    }else{
      
      a = 1
      
      prevtheta <- PostSample[((t-2)*Acc+(1:Acc)), 2:3]
      prevwt <- PostSample[((t-2)*Acc+(1:Acc)), 4]
      covar <- cov.wt(prevtheta, wt = prevwt)$cov
      
      while(a <= Acc){
        
        sim = sim + 1
        
        draws <- importanceDensity(prevtheta, prevwt, covar)
        betadraw <- draws[1]
        gammadraw <- draws[2]
        
        if(dexp(betadraw,1) == 0 | dexp(gammadraw,1) == 0){next}
        
        stochsim <- AbakalikiABCSIR(N=120, beta= betadraw, gamma = gammadraw, t.max = 90, t.step = 1)
        
        suppressWarnings(
          Tsim <- ABCSummaryAbakaliki(stochsim)
        )
        
        if(ABCdistAbakaliki(Tdata, Tsim) > tol[t]){next}else{
          
          wts <- WeightFunction(Acc, draws, prevtheta, prevwt, covar)
          
          PostSample[(((t-1)*Acc)+a), ] <- c(t, betadraw, gammadraw, wts)
          
          a = a + 1
          
        }
      }
      
      t = t+1
      
    }
    
  }
  
  return(list(sim, PostSample))
  
}


#################################################################################################

sin=proc.time()

set.seed(109)

test <- SeqABC(1000, c(12,10,8,6, 5.428, 5, 4.5)) 

sout=proc.time()
sout-sin

test[[1]]

summary(test[[2]][which(test[[2]][,1] == 7), ])
sd(test[[2]][which(test[[2]][,1] == 7), 2])
sd(test[[2]][which(test[[2]][,1] == 7), 3])

hist(test[[2]][which(test[[2]][,1] == 7), 2], breaks = 100)
hist(test[[2]][which(test[[2]][,1] == 7), 3], breaks = 100)


#################################################################################################







































