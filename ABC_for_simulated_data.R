

#################################################################################################


stochasticSIR = function(N, beta, gamma, t.max, t.step = 1){ #t.max = T, t.step = 1 day by default
  S = N-1; I = 1; R = 0
  T = seq(0, t.max, b=t.step)
  res = matrix(nrow = length(T), ncol = 4)
  colnames(res) = c('time', 'S', 'I', 'R')
  res[1,] = c(0,S,I,R)
  
  for(step in 2:length(T)){
    
    pInfec = 1 - exp(-beta*I*t.step)
    Inew = rbinom(1, size = S, prob = pInfec)
    
    
    pRemove = 1 - exp(-gamma*t.step)
    Rnew = rbinom(1, size = I, prob = pRemove)
    
    S = S - Inew
    I = I + Inew - Rnew
    R = R + Rnew
    
    res[step,] = c(T[step], S, I, R)
  }
  res # Outputs a matrix that describes the population of each state at each timestep
}


#################################################################################################


beta  = 0.0015
gamma = 1/7
N = 200
R_0 = N*beta/gamma
R_0

set.seed(32)
stoch = stochasticSIR(N=N, beta= beta, gamma = gamma, t.max = 100, t.step = 1)


#################################################################################################


ABCSummary <- function(res){ #Takes the results of the simulation as an input
  
  FinalSize <- res[101,4] # The total number of infecteds at the end of the outbreak
  TimeI0 <- ifelse( min( which(res[,3] == 0) ) == Inf,
                    -1000,
                    min( which(res[,3] == 0) ) ) # The time at which the last susceptible becomes infected

  
  return(c(FinalSize,TimeI0))
}


#################################################################################################


ABCdist <- function(Tdata, Tsim){ #Takes the summary statistics of the data and simulation as inputs
  
  p <- length(Tdata)
  
  d = rep(NA, p)
  
  for(i in 1:p){
    d[i] <- dist(c(Tdata[i], Tsim[i]), method = "euclidean")
  } 
  
  return(d)
}


#################################################################################################


NSim = 1000000 #Number of simulations 1,000,000

PosteriorSamplesABC = matrix(nrow = NSim, ncol = 4) # A matrix to store the results

sin=proc.time()

set.seed(109)

Tdata <- ABCSummary(stoch)

########

suppressWarnings(
  
  for(i in 1:NSim)
  {
    betadraw = rexp(1,1)
    gammadraw = rexp(1,1)
    
    stochsim = stochasticSIR(N=200, beta= betadraw, gamma = gammadraw, t.max = 100, t.step = 1)
    
    Tsim <- ABCSummary(stochsim)
    
    PosteriorSamplesABC[i, ] <- c(betadraw, gammadraw, ABCdist(Tdata, Tsim))
    
  }
  
)

sout=proc.time()
sout-sin


#################################################################################################


Tol <- c(50,10) # The tolerance for each parameter

sum(PosteriorSamplesABC[,3] <= Tol[1] & PosteriorSamplesABC[,4] <= Tol[2])

library(ggplot2)

pm3 <- as.data.frame(PosteriorSamplesABC)
colnames(pm3) <- c("beta", "gamma", "difb", "difg")

ggplot(data = pm3, aes(x=difb, y=difg)) + geom_point() + theme_bw()+
  labs(title="Distances",x="|Final Size - 157|", y = "|Duration - 79|") + xlim(0, 100) + ylim(0, 100)

geom_line(col=c(rep("black",7), rep("red", 7*250)), alpha = c(rep(1,7), rep(0.1, 7*3), rep(1, 7), rep(0.1, 7*246)) ) + theme_bw()+
  labs(title="ABC: Newly infected individuals between time points",x="Week", y = "Number infected since last time point")+
  theme(legend.position="none")


#################################################################################################


Theta_ABC = PosteriorSamplesABC[which(PosteriorSamplesABC[,3] <= Tol[1] & PosteriorSamplesABC[,4] <= Tol[2]), c(1,2)]


#hist(Theta_ABC[,1], main = "Beta", breaks = 100) # Should be around 0.00015
#hist(Theta_ABC[,2], main = "Gamma", breaks = 100) # Should be around 0.07

summary(Theta_ABC)



#################################################################################################
#################################################################################################


PostSamplesABC = matrix(nrow = 1000, ncol = 4) # A matrix to store the results

sin=proc.time()

set.seed(109)

Tdata <- ABCSummary(stoch)

Accepted = 1000

########

sims = 0
a=0

suppressWarnings(

  while(a <= Accepted)
  {
    sims = sims + 1
    
    betadraw = rexp(1,1/0.005)
    gammadraw = rexp(1,1/0.25)
    
    stochsim = stochasticSIR(N=200, beta= betadraw, gamma = gammadraw, t.max = 100, t.step = 1)
    
    Tsim <- ABCSummary(stochsim)
    
    dist <- ABCdist(Tdata, Tsim)
    
    if(all(dist <= c(10,10))){
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







