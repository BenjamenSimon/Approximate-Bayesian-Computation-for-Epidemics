

#################################################################################################


TrueBins <- function(Bins){
  N <- nrow(Bins)
  
  m = which(Bins[,1] < Bins[,2])
  
  TrueBins = matrix(nrow = length(m)+1, ncol = 3)
  TrueBins[,1] <- c(m, N)
  TrueBins[,2:3] <- rbind(Bins[m,], c( Bins[ (m[length(m)]), 2 ], Inf ))
  
  return(TrueBins)
}

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
      if(abs(Accepted[i] - Final) <= tol){k[i] = (1 - (abs(Accepted[i] - Final)/tol))}else{k[i] = 0}
    }
    
    return(k)
  }
}

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


#################################################################################################


MeaslesData <- matrix(nrow = 3, ncol = 3)
MeaslesData[1, ] <- c(0, 1, 2)
MeaslesData[2, ] <- c(18, 11, 6)
MeaslesData[3, ] <- c(79, 189, 149)

m = (18+11+6)
N = (79+189+149)
n1 = 79
n2 = 189
n3 = 149


#################################################################################################


MeaslesBins <- function(tol){
  
  q <- c(1, runif(2, 0, 1))
  
  
  InfectionOrder <- sample( c(rep(0, n1), rep(1, n2), rep(2, n3)), N, replace = FALSE,
                            prob =  c(rep(q[1], n1), rep(q[2], n2), rep(q[3], n3)) )
  datacheck <- c(length(which(InfectionOrder[1:35] == 0)),
                 length(which(InfectionOrder[1:35] == 1)),
                 length(which(InfectionOrder[1:35] == 2)))
  
  if( sum(abs(datacheck - c(18, 11, 6) ) <= tol) != 3){return(0)}
  
  
  s = matrix(ncol = 3, nrow = N+1)
  s[1,] <- c(n1,n2,n3)
  for(i in 1:N){
    if(InfectionOrder[i] == 0){
      s[(i+1), ] <- s[i, ] - c(1,0,0)
    }
    if(InfectionOrder[i] == 1){
      s[(i+1), ] <- s[i, ] - c(0,1,0)
    }
    if(InfectionOrder[i] == 2){
      s[(i+1), ] <- s[i, ] - c(0,0,1)
    }
  }
  
  s = t(s)
  
  alpha = 0
  for(i in 1:N){
    sum = 0
    for(k in 1:3){
      sum = sum + s[k, (i+1)] * q[k]
    }
    
    alpha[i] = 1/N * sum
  }
  
  
  Ivec = rexp(N, 1)
  
  Lvec = 0
  for(i in 1:(N-1)){
    Lvec[i] = rexp(1, alpha[i])
  }
  
  
 #######################################
  
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
  
  for( n in 2:(N-1)){
    
    Lower = max(Frac[1:(n-1)])
    
    Upper = Frac[n]
    
    Bins[n, ] <- c(Lower, Upper)
  }
  
  return(list(q,Bins))
  
  
}


#################################################################################################


SemiCoupledABC <- function(Acc, tol, tolgroup, Final){
  
  Results <- matrix(nrow = Acc * (tol*2 +1), ncol = 7)
  colnames(Results) <- c("Final", "Lower", "Upper", "k", "sim", "q1", "q2")
  b = 1
  i = 1
  sims = 0
  
  while(b <= Acc){
    sims = sims + 1
    
    Lambda <- MeaslesBins(tolgroup)
    if(length(Lambda) != 2){next}
    q <- Lambda[[1]]
    Lambda <- Lambda[[2]]
    TrueLambda <- TrueBins(Lambda)
    
    #Accepted <- 0
    Accepted <- which(abs(TrueLambda[,1] - Final) <= tol)
    
    if(length(Accepted) == 0){next}
    
    if(length(Accepted) != 0){
      
      k = Kernelfunc(TrueLambda[Accepted, 1], Final, tol, "Step")
      Res = 0
      
      if(length(k) > 1){
        Res <- cbind(TrueLambda[Accepted, ], k, rep(b, length(k)), rep(q[2], length(k)), rep(q[3], length(k)))
      }else{
        Res <- t(as.matrix(c(TrueLambda[Accepted, ], k, b, q[2], q[3])))
      }
      
      for(j in i:(i + length(Accepted) - 1)){
        Results[j,1:7] <-  Res[j - (i-1), 1:7]
      }
      
      b = b + 1
      i = i + length(Accepted)
    }
    
  }  
  return(list(Results,sims))
}


#################################################################################################


MeaslesLamdbaSample <- function(AccSam, B, Res, prior, param1, param2){
  
  a = 1
  
  PostSample <- matrix(nrow = AccSam, ncol = 3)
  draws = 0
  
  while(a <= AccSam){
    draws = draws + 1
    
    Q <- sample(B[,1], size = 1, prob = (B[,4]/sum(B[,4])))
    
    theta <- runif(1, B[Q, 2], B[Q, 3])
    
    accept <- runif(1, 0, 1)
    
    ResQ <- Res[which(Res[,5] == Q), ,drop = F] 
    
    #if(length(ResQ) != ncol(Res)){
    Ai <- max(which(ResQ[,2] < theta))
    k = ResQ[Ai,4]
    #}else{
    #  k = ResQ[4]
    #}
    
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
      PostSample[a,1:3] <- c(theta, ResQ[1,6:7])
      a = a+1
    }else{next}
  }
  return(list(PostSample,draws)) 
}


#################################################################################################


set.seed(19)

sin=proc.time()

SemiCoupledSim <- SemiCoupledABC(1000, 2, 0, 35)
AcceptedStage1 <- SemiCoupledSim[[2]]
SemiCoupledSim <- SemiCoupledSim[[1]]
(1000/AcceptedStage1)*100

sout=proc.time()
sout-sin


###############################################################################################


#Lpi <- LpiCalc(B, "Normal", 3, 1.5)

#MeaslesPostSample <- MeaslesLamdbaSample(10000, B, Res = SemiCoupledSim, "Normal", 3, 1.5)

#AcceptedStage2 <- MeaslesPostSample[[2]]
#MeaslesPostSample <- MeaslesPostSample[[1]]
#(10000/AcceptedStage2)*100

#summary(MeaslesPostSample)
#hist(MeaslesPostSample[,1], breaks = 100)
#hist(MeaslesPostSample[,2], breaks = 100)
#hist(MeaslesPostSample[,3], breaks = 100)


#################################################################################################


sin=proc.time()

B <- RecoverB(SemiCoupledSim)

Lpi <- LpiCalc(B, "Exp", 0.1, 0)

MeaslesPostSample <- MeaslesLamdbaSample(10000, B, Res = SemiCoupledSim, "Exp", 0.1, 1)
AcceptedStage2 <- MeaslesPostSample[[2]]
MeaslesPostSample <- MeaslesPostSample[[1]]
(10000/AcceptedStage2)*100

sout=proc.time()
sout-sin

summary(MeaslesPostSample)
sd(MeaslesPostSample[,1])
sd(MeaslesPostSample[,2])
sd(MeaslesPostSample[,3])



hist(MeaslesPostSample[,1], breaks = 100)
hist(MeaslesPostSample[,2], breaks = 100)
hist(MeaslesPostSample[,3], breaks = 100)




#################################################################################################

































































