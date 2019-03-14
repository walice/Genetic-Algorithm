# Genetic Algorithm
# Alice Lepissier

## ## ## ## ## ## ## ## ## ## ##
# INDEX                     ####
## ## ## ## ## ## ## ## ## ## ##
# Preamble
# Data
# Matching
# .. Propensity Score Matching
# .. Genetic Matching
# Genetic Algorithm
# .. Asymmetric Double Claw function
# .. Algorithm



## ## ## ## ## ## ## ## ## ## ##
# PREAMBLE                  ####
## ## ## ## ## ## ## ## ## ## ##

setwd("C:/Users/Alice/Box Sync/PhD/Statistics/PSTAT 232")
#install.packages("genalg")
#install.packages("GA")
#devtools::install_github('drizztxx/gatbxr')
#install.packages("Matching")
#install.packages("nor1mix")
#install.packages("rgenoud")
library(genalg)
library(GA)
library(gatbxr)
library(Matching)
library(nor1mix)
library(rgenoud)



## ## ## ## ## ## ## ## ## ## ##
# DATA                      ####
## ## ## ## ## ## ## ## ## ## ##

data("lalonde")
attach(lalonde)
Y <- lalonde$re78
Tr <- lalonde$treat



## ## ## ## ## ## ## ## ## ## ##
# MATCHING                  ####
## ## ## ## ## ## ## ## ## ## ##

# .. Propensity Score Matching ####
glm <- glm(Tr ~ age + educ + black + hisp 
           + married + nodegr + re74 + re75, 
           family = binomial, data = lalonde)
match.ps <- Match(Y = Y, Tr = Tr, X = glm$fitted)
MatchBalance(Tr ~ age + I(age^2) + educ + I(educ^2) 
             + black + hisp + married + nodegr 
             + re74 + I(re74^2) + re75 + I(re75^2) 
             + u74 + u75 + I(re74 * re75)
             + I(age * nodegr) + I(educ * re74)
             + I(educ * re75),
             match.out = match.ps, 
             nboots = 1000, data = lalonde)


# .. Genetic Matching ####

X <- cbind(age, educ, black, hisp, married, 
           nodegr, re74, re75, u74, u75)
BalanceMatrix <- cbind(age, I(age^2), educ, I(educ^2), 
                       black, hisp, married, nodegr, 
                       re74, I(re74^2), re75, I(re75^2), 
                       u74, u75, I(re74 * re75), 
                       I(age * nodegr), I(educ * re74), 
                       I(educ * re75))
GA <- GenMatch(Tr = Tr, X = X, 
               BalanceMatrix = BalanceMatrix, 
               pop.size = 1000)
match.GA <- Match(Y = Y, Tr = Tr, X = X, 
                  Weight.matrix = GA)
MatchBalance(Tr ~ age + I(age^2) + educ + I(educ^2) 
             + black + hisp + married + nodegr 
             + re74 + I(re74^2) + re75 + I(re75^2) 
             + u74 + u75 + I(re74 * re75)
             + I(age * nodegr) + I(educ * re74)
             + I(educ * re75),
             match.out = match.GA, 
             nboots = 1000, data = lalonde)



## ## ## ## ## ## ## ## ## ## ##
# GENETIC ALGORITHM         ####
## ## ## ## ## ## ## ## ## ## ##

# .. Asymmetric Double Claw function ####
claw <- function(xx) {
  x <- xx[1]
  y <- (0.46 * (dnorm(x, -1, 2/3) + dnorm(x, 1, 2/3)) +
            (1/300) * (dnorm(x, -0.5, 0.01) + dnorm(x, -1, 0.01) 
                       + dnorm(x, -1.5, 0.01)) +
            (7/300) * (dnorm(x, 0.5, 0.07) + dnorm(x, 1, 0.07) 
                       + dnorm(x, 1.5, 0.07)))
  return(y)
}
xx <- seq(-3,3, by = 0.01)
yy <- NULL
for (i in 1:length(xx)) {
  yy[i] <- claw(xx[i])
}
plot(xx, yy, type = "l")


# .. Using genoud() function ####
GA1.claw <- genoud(claw, nvars = 1, max = TRUE, pop.size = 3000)
GA1.sol <- GA1.claw$par
GA1.solvalue <- GA1.claw$value
claw(GA1.sol) == GA1.solvalue


# .. Using ga() function ####
GA2.claw <- ga(type = "real-valued", 
               fitness = claw, lower = -3, upper = 3)
summary(GA2.claw)
GA2.pop <- GA2.claw@population
GA2.fitness <- GA2.claw@fitness
plot(GA2.pop)
plot(GA2.fitness)
GA2.sol <- GA2.claw@solution
GA2.solvalue <- GA2.claw@fitnessValue
claw(GA2.sol) == GA2.solvalue

rm(list = ls())


# .. Algorithm ####
GA <- function(fitness.func, pop.size = 500, max.iter = 100,
               lower = NULL, upper = NULL, init.pop = NULL){
  
  # Initialize empty list
  GA.out <- list()
  fitness <- NULL
  prob <- NULL
  
  # Take suggestions for initial population,
  # else sample from uniform
  if (!missing(init.pop)){
    pop <- init.pop
  } else {
    pop <- runif(pop.size, lower, upper)
  }
  
  iter <- 0
  while (iter < max.iter){
    
    iter <- iter + 1
    cat("Iteration: ", iter, "\n")
    
    # Evaluate fitness
    for (i in 1:pop.size) {
      fitness[i] <- fitness.func(pop[i])
    }
    
    # Roulette wheel selection
    previous_prob <- 0
    
    for (i in 1:pop.size){
      prob[i] <- previous_prob + fitness[i]/sum(fitness)
      previous_prob <- prob[i]
    }
    
    newpop <- pop
    for (i in 1:pop.size){
      u <- runif(1)
      if (length(which(prob < u)) == 0){
        newpop[i] <- pop[1]
      } else {
        newpop[i] <- pop[max(which(prob < u))]
      }
    }
    
    # Crossover
    parent_A <- newpop[1:(pop.size/2)]
    parent_B <- newpop[(pop.size/2+1):pop.size]
    p <- runif(1)
    offspring_A <- p*parent_A + (1-p)*parent_B
    offspring_B <- p*parent_B + (1-p)*parent_A
    
    newpop <- c(offspring_A, offspring_B)
    
    # Mutation
    #newpop <- runif(pop.size, min(newpop), max(newpop))
    
    for (i in 1:pop.size) {
      fitness[i] <- fitness.func(newpop[i])
    }
    cat("Fitness value: ", max(fitness), "\n")
    pop <- newpop
  }
  
  fitness.value <- max(fitness)
  solution <- pop[which(fitness == fitness.value)]
  
  GA.out$population <- pop
  GA.out$fitness <- fitness
  GA.out$solution <- solution
  GA.out$fitness.value <- fitness.value
  
  return(GA.out)
}

GA3 <- GA(claw, pop.size = 50, max.iter = 100,
   lower = -3, upper = 3)
plot(GA3$fitness)
