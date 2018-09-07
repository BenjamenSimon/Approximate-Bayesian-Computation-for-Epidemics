

#################################################################################################


SellkeEpidemic <- function(lambda, N){
  
  Ivec = rexp(N-1,1)
  
  Lvec = 0
  for(i in 1:(N-1)){
    Lvec[i] <- rexp(1,(N - i)/N)
  }
  
  SumIvec = 0
  for(i in 1:(N-1)){
    SumIvec[i] <- lambda * sum(Ivec[1:i])
  }
  
  SumLvec = 0
  for(i in 1:(N-1)){
    SumLvec[i] <- sum(Lvec[1:i])
  }
  
  
  FinalSize = suppressWarnings(which.min(SumIvec > SumLvec))
  if(all(SumIvec > SumLvec)){FinalSize = N}
  
  
  return(FinalSize)
  
} 


#################################################################################################


FinalSize <- 0

for(i in 1:1000){
  FinalSize[i] <- SellkeEpidemic(2, 100)
}

hist(FinalSize, breaks = 100)


#################################################################################################


LambdaBin <- function(N, Final){
  
  # N = Total Population Size
  # Final = The final size required
  
  Ivec = rexp(N-1,1)
  
  Lvec = 0
  for(i in 1:(N-1)){
    Lvec[i] <- rexp(1,(N - i)/N)
  }
  
  
  Bins <- matrix(nrow = N, ncol = 2)
  
  
  SumIvec = 0
  for(i in 1:(N-1)){
    SumIvec[i] <- sum(Ivec[1:i])
  }
  
  SumLvec = 0
  for(i in 1:(N-1)){
    SumLvec[i] <- sum(Lvec[1:i])
  }
  
  Frac = 0
  for(i in 1:(N-1)){
    Frac[i] <- SumLvec[i] / SumIvec[i]
  }
  
  Bins[1,] <- c(0, Frac[1])
  
  for( m in 2:(N-1)){
    
    Lower = max(Frac[1:(m-1)])
    
    Upper = Frac[m]
    
    Bins[m, ] <- c(Lower, Upper)
  }
  
  return(Bins)
  
} 


#################################################################################################


LambdaTest <- LambdaBin(100, 38)

TrueBins <- function(Bins){
  N <- nrow(Bins)
  
  m = which(Bins[,1] < Bins[,2])
  
  TrueBins = matrix(nrow = length(m)+1, ncol = 3)
  TrueBins[,1] <- c(m, N)
  TrueBins[,2:3] <- rbind(Bins[m,], c( Bins[ (m[length(m)]), 2 ], Inf ))
  
  return(TrueBins)
}

LambdaTestBins <- TrueBins(LambdaTest)

plot(LambdaTestBins[,2],LambdaTestBins[,1], type = 's')


#################################################################################################


Kernelfunc <- function(Accepted, Final, tol, method){
  
  k = 0
  
  if(method == "Step"){
    
    for(i in 1:length(Accepted)){
      if(abs(Accepted[i] - Final) <= tol){k[i] = 1}else{k[i] = 0}
    }
    
    return(k)
  }
  
  if(method == "Ep"){
    
    for(i in 1:length(Accepted)){
      if(abs(Accepted[i] - Final) <= 1){k[i] = 3/4 * (1 - (abs(Accepted[i] - Final)/tol))}else{k[i] = 0}
    }
    
    return(k)
  }
}


#################################################################################################


CoupledABC <- function(Acc, N, m, Final, tol){
  
  Results <- matrix(nrow = Acc * (tol*2 +1), ncol = 5)
  b = 1
  i = 1
  sims = 0
  
  while(b <= Acc){
    sims = sims + 1
    
    Lambda <- LambdaBin(N, Final)
    TrueLambda <- TrueBins(Lambda)
    
    #Accepted <- 0
    Accepted <- which(abs(TrueLambda[,1] - Final) <= tol)
    
    if(length(Accepted) == 0){next}
    
    if(length(Accepted) != 0){
      
      k = Kernelfunc(TrueLambda[Accepted, 1], Final, tol, "Ep")
      Res = 0
      
      if(length(k) > 1){
        Res <- cbind(TrueLambda[Accepted, ], k, rep(b, length(k)))
      }else{
        Res <- t(as.matrix(c(TrueLambda[Accepted, ], k, b)))
      }
      
      for(j in i:(i + length(Accepted) - 1)){
        Results[j,1:5] <-  Res[j - (i-1), 1:5]
      }
      
      b = b + 1
      i = i + length(Accepted)
    }
    
  }  
  return(list(Results,sims))
}


#################################################################################################


set.seed(25)

CoupledSim <- CoupledABC(100, 100, 1, 38, 2)
Acceptedc_ABC <- CoupledSim[[2]]
CoupledSim <- CoupledSim[[1]]

(100/Acceptedc_ABC)*100


#################################################################################################


RecoverB <- function(Res){
  
  Acc = max(Res[,5], na.rm = T)
  
  B <- matrix(nrow = Acc, ncol = 4)
  
  for(i in 1:Acc){
    
    Set <- which(Res[,5] == i)
    
    Bmax <- max(Res[Set, 3])
    Bmin <- min(Res[Set,2])
    
    B[i,] <- c(i, Bmin, Bmax, (Bmax - Bmin))
    
  }
  
  return(B)
  
}

B <- RecoverB(CoupledSim)


#################################################################################################


LpiCalc <- function(B, prior, param1, param2){
  
  Acc <- max(B[,1])
  
  # Uniform
  if(prior == "Unif"){
    Lpi <- max(1/B[,4])
  }
  
  
  # Normal
  if(prior == "Normal"){
    mean = param1
    sd = param2
    
    ClosestMean <- 0
    for(i in 1:Acc){
      if(B[i,2] > mean){ClosestMean[i] = B[i,2]}
      if(B[i,3] < mean){ClosestMean[i] = B[i,3]}
      if(B[i,3] > mean & B[i,2] < mean){ClosestMean[i] = mean}
    }
    
    Lpi <- dnorm(ClosestMean[which.min(abs(ClosestMean - mean))], sd)
    
  }
  
  #Exponential
  if(prior == "Exp"){
    min <- min(B[,2])
    Lpi <- dexp(min, param1)
  }
  
  return(Lpi)
}


Lpi <- LpiCalc(B, "Normal", 1.177, 0.211)


#################################################################################################


LamdbaSample <- function(AccSam, B, Res, prior, param1, param2){
  
  a = 1
  
  PostSample <- 0
  draws = 0
  
  while(a <= AccSam){
    draws = draws + 1
    
    Q <- sample(B[,1], size = 1, prob = (B[,4]/sum(B[,4])))
    
    theta <- runif(1, B[Q, 2], B[Q, 3])
    
    accept <- runif(1, 0, 1)
    
    ResQ <- Res[which(Res[,5] == Q), ] 
    
    if(length(ResQ) != ncol(Res)){
      Ai <- max(which(ResQ[,2] < theta))
      k = ResQ[Ai,4]
    }else{
      k = ResQ[4]
    }
    
    if(prior == "Normal"){
      acceptprob <- k * dnorm(theta, param1, param2) / Lpi
    }
    if(prior == "Exp"){
      acceptprob <- k * dexp(theta, param1) / Lpi
    }
    if(prior == "Unif"){
      acceptprob <- 1
    }
    
    
    if(accept < acceptprob){
      PostSample[a] <- theta
      a = a+1
    }else{next}
  }
  return(list(PostSample,draws)) 
}


#################################################################################################


set.seed(19)

sin=proc.time()

CoupledSim <- CoupledABC(10000, 120, 1, 30, 2)
AcceptedStage1 <- CoupledSim[[2]]
CoupledSim <- CoupledSim[[1]]
(10000/AcceptedStage1)*100

B <- RecoverB(CoupledSim)
#Lpi <- LpiCalc(B, "Normal", 1.177, 0.211)
Lpi <- LpiCalc(B, "Exp", 1, 0)

#LambdaPostSample <- LamdbaSample(10000, B, CoupledSim, "Normal", 1.177, 0.211)
LambdaPostSample <- LamdbaSample(10000, B, CoupledSim, "Exp", 1, 0)
AcceptedStage2 <- LambdaPostSample[[2]]
LambdaPostSample <- LambdaPostSample[[1]]

sout=proc.time()
sout-sin # 

(10000/AcceptedStage2)*100

summary(LambdaPostSample)
hist(LambdaPostSample)


#################################################################################################



R0 <- (89*PostSamplesABC[,1]) / (1-exp(-PostSamplesABC[,2]))
mean(R0)
median(R0)






