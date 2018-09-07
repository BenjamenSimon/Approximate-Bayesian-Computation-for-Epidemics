
library(ape)

##############################################################################################

canedata = scan("canedata.txt")

canematrix = matrix(canedata,ncol=8,byrow=TRUE)

colnames(canematrix) <- c("x-coord", "y-coord", "TP1", "TP2", "TP3", "TP4", "TP5", "TP6")

# One row of the matrix per plant
# (x,y) given by the first two columns
# The infectious status (0-susceptible; 1-infectious) at the 6 time points 
# given in columns 3-8. 
#

newcanematrix <- canematrix

newcanematrix <- cbind(seq(1, 1742), newcanematrix)

TimeInf <- rep(99, 1742)

TimeInf[which(canematrix[,8] == 1)] <- 30
TimeInf[which(canematrix[,7] == 1)] <- 23
TimeInf[which(canematrix[,6] == 1)] <- 19
TimeInf[which(canematrix[,5] == 1)] <- 15
TimeInf[which(canematrix[,4] == 1)] <- 11
TimeInf[which(canematrix[,3] == 1)] <- 6

newcanematrix <- cbind(newcanematrix, TimeInf)

summary(as.factor(newcanematrix[,10]))


##############################################################################################


DistanceMatrix <- as.matrix(dist(canematrix[,1:2]))

UniqueDist <- unique(as.vector(DistanceMatrix))

whichUnique <- list()
for(i in 1:length(UniqueDist)){
  whichUnique[[i]] = which(DistanceMatrix == UniqueDist[i])
}

ForceInf <- function(alpha){
  
  UniqueDensity <- exp( (- (UniqueDist)^2 ) / (2 * alpha^2) ) / (sqrt(2 * pi) * alpha)
  
  UniquePairs <- cbind(UniqueDist, UniqueDensity)
  UniquePairs[1,2] <- 0
  
  ForceM <- diag(x = 0, 1742)
  
  for(i in 1:length(UniqueDist)){
    ForceM[whichUnique[[i]]] <- UniquePairs[i,2]
  }
  
  return(ForceM)
  
}


##############################################################################################


RateofInf <- function(Infds, Susbs, Force, Lambda, r){
  
  sum <- sum(rowSums(Force[Susbs, Infds]))
  
  rate = Lambda * ( r*length(Susbs) + sum)
  
  return(rate)
}


##############################################################################################


ProbofInf <- function(Infds, Susbs, Force, Lambda, r, rate){
  
  rowsum <- rowSums(Force[Susbs, Infds])
  
  personalrates <- r*Lambda + Lambda * rowsum
  
  probs = personalrates / rate
  
  return(probs)
}


##############################################################################################


inv.dist.matrix <- 1/DistanceMatrix
diag(inv.dist.matrix) <- 0


SummarySC <- function(res){ #Takes the results of the simulation as an input
  
  newinf <- res[which( (res[,5] <= 11) & (res[,5] > 6) ), , drop = FALSE]
  
  Z <- nrow(newinf)
  if(is.null(Z)){Z = 0}
  
  I <- Moran.I(res[,4], inv.dist.matrix)$observed
  
  return(c(Z,I))
}


SCSum <- SummarySC( newcanematrix[,c(1:3, 4, 10)] ) #change middle one to approprite column 4,5,6,7,8 or 9



##############################################################################################


ContinuousSI <- function(N = 1742, I = 6, alpha, r, lambda){
  
  N = N+9
  S = N - I
  t = 6
  
  time <- c(0, 6, 11)
  snap = 1
  SimSum = rep(NA, 12)
  
  Force <- ForceInf(alpha)
  
  results <- matrix(ncol = 6, nrow = N)
  ########
  colnames(results) <- c("ID", "x-coord", "y-coord", "Inf status", "When inf'd", "Who.inf.under.TP")
  results[,1] <- seq(1, N)
  results[,2] <- rep(seq(1.5, 25.5, 1.5), each = 103)
  results[,3] <- seq(0, 51, 0.5)
  results[,4] <- 0
  results[,6] <- 0
  
  results[620, 2] <- 10.2
  results <- results[-(seq(1734, 1751, 2)), ]
  
  ########
  
  InitInf <- c(295, 452, 1217, 1468, 1597, 1612)
  results[InitInf, 4] <- 1
  results[InitInf, 5] <- 0
  
  Infds <- which(results[,4] == 1)
  
  while(t < 15 & length(Infds) <= (6+18+40)){
    
    Infds <- which(results[,4] == 1)
    
    Susbs <- seq(1,1742)[-Infds]
    
    rate <- RateofInf(Infds, Susbs, Force, lambda, r)
    
    tau = rexp(1, rate)
    
    probs = ProbofInf(Infds, Susbs, Force, lambda, r, rate)
    
    NewI <- sample(Susbs, 1, prob = probs)
    
    t = t + tau
    
    results[NewI, 4] <- 1
    
    results[NewI, 5] <- t
    
  }
  
  results[which(results[,5] < 11) ,6] <- 1
  
  return(results)
  
}

test <- ContinuousSI(N = 1742, I = 6, 0.8, 0.02, 1/30)

SummarySC(test[ ,c(1:3,6,5)])

##############################################################################################


SCABC <- function(Acc, lprior, tolZ, tolI){
  
  Results <- matrix(nrow = Acc, ncol = 5)
  colnames(Results) <- c("alphadraw","rdraw","lambdadraw", "Size.at.TP", "I.at.TP")
  b = 1
  sims = 0
  
  while(b <= Acc){
    sims = sims + 1
    
    alpha.draw = runif(1,0,5)
    r.draw = runif(1,0,1)
    lambda.draw = rexp(1,lprior)
    
    SCSim <- ContinuousSI(N = 1742, I = 6, alpha = alpha.draw, r = r.draw, lambda = lambda.draw)
    
    SumSim <- SummarySC(SCSim)
    
    
    if( abs(SumSim[1] - SCSum[1]) < tolZ & abs(SumSim[2] - SCSum[2]) < tolI ){
      Results[b, 1:3] <- c(alpha.draw,r.draw,lambda.draw)
      Results[b, 4:5] <- SumSim
      
      print(c("Accepted  =", b))
      b = b + 1
    }else{next}
    
  }
  
  return(list(Results,sims))
}


##############################################################################################


sin=proc.time()

set.seed(19)

SCABCtest <- SCABC(10, 30, 2, 0.02)

sout=proc.time()
sout-sin

SCABCtest[[2]]

SCABCtest <- as.data.frame(SCABCtest[[1]])

#write.csv(SCABCtest, file = "~/250SCABCpost")


250/SCABCtest[[2]]*100 #acceptance rate
sout=proc.time()
sout-sin








