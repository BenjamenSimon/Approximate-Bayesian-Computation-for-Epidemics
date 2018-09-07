
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


##############################################################################################

row1 <- which(canematrix[,1] == 1.5)
row2 <- which(canematrix[,1] == 3)
row3 <- which(canematrix[,1] == 4.5)
row4 <- which(canematrix[,1] == 6)
row5 <- which(canematrix[,1] == 7.5)
row6 <- which(canematrix[,1] == 9)
rowoff <- which(canematrix[,1] == 10.2)
row7 <- which(canematrix[,1] == 10.5)
row8 <- which(canematrix[,1] == 12)
row9 <- which(canematrix[,1] == 13.5)
row10 <- which(canematrix[,1] == 15)
row11 <- which(canematrix[,1] == 16.5)
row12 <- which(canematrix[,1] == 18)
row13 <- which(canematrix[,1] == 19.5)
row14 <- which(canematrix[,1] == 21)
row15 <- which(canematrix[,1] == 22.5)
row16 <- which(canematrix[,1] == 24)
row17 <- which(canematrix[,1] == 25.5)

row1inf <- sum(canematrix[row1, 8] == 1)
row2inf <- sum(canematrix[row2, 8] == 1)
row3inf <- sum(canematrix[row3, 8] == 1)
row4inf <- sum(canematrix[row4, 8] == 1)
row5inf <- sum(canematrix[row5, 8] == 1)
row6inf <- sum(canematrix[row6, 8] == 1)
row7inf <- sum(canematrix[row7, 8] == 1)
row8inf <- sum(canematrix[row8, 8] == 1)
row9inf <- sum(canematrix[row9, 8] == 1)
row10inf <- sum(canematrix[row10, 8] == 1)
row11inf <- sum(canematrix[row11, 8] == 1)
row12inf <- sum(canematrix[row12, 8] == 1)
row13inf <- sum(canematrix[row13, 8] == 1)
row14inf <- sum(canematrix[row14, 8] == 1)
row15inf <- sum(canematrix[row15, 8] == 1)
row16inf <- sum(canematrix[row16, 8] == 1)
row17inf <- sum(canematrix[row17, 8] == 1)


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


RateofInf <- function(Infds, Susbs, Force, r){
  
  sum <- sum(rowSums(Force[Susbs, Infds]))
  
  rate = r*length(Susbs) + sum
  
  return(rate)
}


##############################################################################################


ProbofInf <- function(Infds, Susbs, Force, r, rate){
  
  rowsum <- rowSums(Force[Susbs, Infds])
  
  personalrates <- r + rowsum
  
  probs = personalrates / rate
  
  return(probs)
}


##############################################################################################

inv.dist.matrix <- 1/DistanceMatrix
diag(inv.dist.matrix) <- 0

SummarySC <- function(res){ #Takes the results of the simulation as an input
  
  I <- Moran.I(res[,4], inv.dist.matrix)$observed
  
  row1inf <- sum(res[row1, 4] == 1)
  row2inf <- sum(res[row2, 4] == 1)
  row3inf <- sum(res[row3, 4] == 1)
  row4inf <- sum(res[row4, 4] == 1)
  row5inf <- sum(res[row5, 4] == 1)
  row6inf <- sum(res[row6, 4] == 1)
  row7inf <- sum(res[row7, 4] == 1)
  row8inf <- sum(res[row8, 4] == 1)
  row9inf <- sum(res[row9, 4] == 1)
  row10inf <- sum(res[row10, 4] == 1)
  row11inf <- sum(res[row11, 4] == 1)
  row12inf <- sum(res[row12, 4] == 1)
  row13inf <- sum(res[row13, 4] == 1)
  row14inf <- sum(res[row14, 4] == 1)
  row15inf <- sum(res[row15, 4] == 1)
  row16inf <- sum(res[row16, 4] == 1)
  row17inf <- sum(res[row17, 4] == 1)
  
  return(c(I, row1inf, row2inf, row3inf, row4inf, row5inf, row6inf,
              row7inf, row8inf, row9inf, row10inf, row11inf, row12inf,
              row13inf, row14inf, row15inf, row16inf, row17inf))
}


SCSum <- SummarySC( newcanematrix[,c(1:3, 5, 10)] ) #change middle one to approprite column 4,5,6,7,8 or 9


##############################################################################################


SCcoupledSI <- function(alpha, r){
  
  N = 1751
  t = 6
  
  Force <- ForceInf(alpha)
  
  results <- matrix(ncol = 5, nrow = N)
  ########
  colnames(results) <- c("ID", "x-coord", "y-coord", "Inf status", "When inf'd")
  results[,1] <- seq(1, N)
  results[,2] <- rep(seq(1.5, 25.5, 1.5), each = 103)
  results[,3] <- seq(0, 51, 0.5)
  results[,4] <- 0
  
  results[620, 2] <- 10.2
  results <- results[-(seq(1734, 1751, 2)), ]
  
  ########
  
  InitInf <- c(295, 452, 1217, 1468, 1597, 1612)
  results[InitInf, 4] <- 1
  results[InitInf, 5] <- 0
  
  Infds <- which(results[,4] == 1)
  
  while(length(Infds) <= 22){ #2 less than goal because of when we look as who is infected
    
    Infds <- which(results[,4] == 1)
    
    Susbs <- seq(1,1742)[-Infds]
    
    rate <- RateofInf(Infds, Susbs, Force, r)
    
    tau = rexp(1, rate)
    
    probs = ProbofInf(Infds, Susbs, Force, r, rate)
    
    NewI <- sample(Susbs, 1, prob = probs)
    
    t = t + tau
    
    results[NewI, 4] <- 1
    
    results[NewI, 5] <- t
    
  }
  
  return(results)
  
}

test <- SCcoupledSI(0.8, 0.02)

sum(test[,5] >= 6 , na.rm = T)

##############################################################################################


SCcoupledABC <- function(Acc, tolZ, tolI){
  
  Results <- matrix(nrow = Acc, ncol = 5)
  colnames(Results) <- c("alphadraw","rdraw","lambda", "I.at.TP", "Diff.in,row.sums")
  b = 1
  sims = 0
  
  while(b <= Acc){
    sims = sims + 1
    
    alpha.draw = runif(1,0,5)
    r.draw = runif(1,0,1)
    
    SCSim <- SCcoupledSI(alpha = alpha.draw, r = r.draw)
    
    SCSim[,5] <- SCSim[,5]-6
    
    
    #find max time
    last.inf <- max(SCSim[,5], na.rm = T)
    #find lambda to make max time 30
    lambda <- last.inf / 5
    #multiply in lambda
    SCSim[,5] <- SCSim[,5] / lambda
    SCSim[,5] <- SCSim[,5] + 6
    #take summary statistics
    SumSim <- SummarySC(SCSim)

    #compare summary statistics
    diff  <- abs(SCSum - SumSim)
    print(c(diff[1], sum(diff[2:18]), "sim =", sims))

    #accept if...
    
    if(diff[1] <= tolI & sum(diff[2:18]) <= tolZ){
      Results[b, 1:3] <- c(alpha.draw,r.draw,lambda)
      Results[b, 4:5] <- c(SumSim[1], sum(diff[2:18]))
      
      print(b)
      b = b + 1
    }else{next}
    
  }
  
  return(list(Results,sims))
}

sin=proc.time()

set.seed(19)

SCABCtest <- SCcoupledABC(10, 8, 0.02)

SCABCtest[[2]]

SCABCtest <- as.data.frame(SCABCtest[[1]])

sout=proc.time()
sout-sin

#write.csv(SCABCtest, file = "~/250SCscABCpost")

(sout-sin)/60 #mins


##############################################################################################














