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
# Test Functions
# .. Asymmetric Double Claw function
# .. Rosenbrock function
# .. Rastrigin function
# GA Packages
# .. Using genoud() function
# .. Using ga() function
# Genetic Optimization
# .. Genetic operators
# .. Algorithm



## ## ## ## ## ## ## ## ## ## ##
# PREAMBLE                  ####
## ## ## ## ## ## ## ## ## ## ##

setwd("C:/Users/Alice/Box Sync/PhD/Statistics/PSTAT 232")
#install.packages("genalg")
#install.packages("GA")
#devtools::install_github('drizztxx/gatbxr')
#install.packages("Matching")
#install.packages("MatchIt")
#install.packages("nor1mix")
#install.packages("rgenoud")
library(genalg)
library(GA)
library(gatbxr)
library(Matching)
#library(MatchIt)
library(nor1mix)
library(reshape2)
library(rgenoud)
library(tidyverse)



## ## ## ## ## ## ## ## ## ## ##
# TEST FUNCTIONS            ####
## ## ## ## ## ## ## ## ## ## ##

# .. Asymmetric Double Claw function ####
claw <- function(x) {
  f <- (0.46 * (dnorm(x, -1, 2/3) + dnorm(x, 1, 2/3)) +
            (1/300) * (dnorm(x, -0.5, 0.01) + dnorm(x, -1, 0.01) 
                       + dnorm(x, -1.5, 0.01)) +
            (7/300) * (dnorm(x, 0.5, 0.07) + dnorm(x, 1, 0.07) 
                       + dnorm(x, 1.5, 0.07)))
  return(f)
}
x <- seq(-3,3, by = 0.01)
y <- NULL
for (i in 1:length(x)) {
  y[i] <- claw(x[i])
}
plot(x, y, type = "l")


# .. Rosenbrock function ####
rosenbrock <- function(x1, x2){
  a <- 1
  b <- 100
  f <- (a - x1)^2 + b*(x2 - x1^2)^2
  return(f)
}
x1 <- seq(-2, 2, by = 0.5)
x2 <- seq(-1, 3, by = 0.5)
y <- outer(x1, x2, rosenbrock)
persp3D(x1, x2, y)


# .. Rastrigin function ####
rastrigin <- function(x1, x2){
  f <- 20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
  return(f)
}
x1 <- x2 <- seq(-5.12, 5.12, by = 0.5)
y <- outer(x1, x2, rastrigin)
persp3D(x1, x2, y)



## ## ## ## ## ## ## ## ## ## ##
# GA PACKAGES               ####
## ## ## ## ## ## ## ## ## ## ##

# .. Using genoud() function ####
GA1.claw <- genoud(claw, nvars = 1, max = TRUE, pop.size = 3000,
                   BFGS = FALSE)
GA1.sol <- GA1.claw$par
GA1.solvalue <- GA1.claw$value
claw(GA1.sol) == GA1.solvalue

ros <- function(sol){
  x1 <- sol[1]
  x2 <- sol[2]
  a <- 1
  b <- 100
  f <- (a - x1)^2 + b*(x2 - x1^2)^2
  return(f)
}
GA1.rosenbrock <- genoud(ros, nvars = 2, max = FALSE,
                         Domains = matrix(c(-2, 2, -1, 3), 
                                          nrow = 2, byrow = T),
                         pop.size = 3000, max.generations = 100,
                         BFGS = F)
GA1.sol <- GA1.rosenbrock$par
GA1.solvalue <- GA1.rosenbrock$value
rosenbrock(GA1.sol[1], GA1.sol[2]) == GA1.solvalue


# .. Using ga() function ####
GA2.claw <- ga(type = "real-valued", 
               fitness = claw, lower = -3, upper = 3,
               keepBest = TRUE)
summary(GA2.claw)
GA2.pop <- GA2.claw@population
GA2.fitness <- GA2.claw@fitness
plot(GA2.pop)
plot(GA2.fitness)
GA2.sol <- GA2.claw@solution
GA2.solvalue <- GA2.claw@fitnessValue
claw(GA2.sol) == GA2.solvalue

GA2.rosenbrock <- ga(type = "real-valued",
                    fitness = function(sol) -ros(sol),
                    lower = c(-2, -1), upper = c(2, 3),
                    keepBest = TRUE, seed = 232)
summary(GA2.rosenbrock)
GA2.pop <- GA2.rosenbrock@population
GA2.fitness <- GA2.rosenbrock@fitness
plot(GA2.pop)
plot(GA2.fitness)
GA2.sol <- GA2.rosenbrock@solution
GA2.solvalue <- GA2.rosenbrock@fitnessValue
rosenbrock(GA2.sol[1], GA2.sol[2]) == -GA2.solvalue

GA2.rastrigin <- ga(type = "real-valued",
                    fitness = function(x) -rastrigin(x[1], x[2]),
                    lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
                    seed = 232)
summary(GA2.rastrigin)
GA2.pop <- GA2.rastrigin@population
GA2.fitness <- GA2.rastrigin@fitness
plot(GA2.pop)
plot(GA2.fitness)
GA2.sol <- GA2.rastrigin@solution
GA2.solvalue <- GA2.rastrigin@fitnessValue
rastrigin(GA2.sol[1], GA2.sol[2]) == -GA2.solvalue



## ## ## ## ## ## ## ## ## ## ##
# GENETIC OPTIMIZATION      ####
## ## ## ## ## ## ## ## ## ## ##


# .. Genetic operators ####

# Roulette wheel selection
roulette <- function(pop, fitness){
  pop <- pop
  fitness <- fitness
  pop.size <- nrow(pop)
  prob <- abs(fitness)/sum(abs(fitness))
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  select <- sample(1:pop.size, size = pop.size,
                   prob = prob, replace = T)
  pop <- pop[select, ]
  return(as.matrix(pop))
}

# Selection proportional to fitness with linear scaling
prop.linear.scaling <- function(pop, fitness){
  pop <- pop
  fitness <- fitness
  pop.size <- nrow(pop)
  fitness.min <- min(fitness, na.rm = T)
  if (fmin <- 0){
    fitness <- fitness - fitness.min
    fitness.min <- min(fitness, na.rm = T)
  }
  fitness.mean <- mean(fitness, na.rm = T)
  fitness.max <- max(fitness, na.rm = T)
  
  sfactor <- 2
  eps <- sqrt(.Machine$double.eps)
  # transform f -> f' = a*f + b such that
  if(fitness.min > (sfactor*fitness.mean - fitness.max)/(sfactor-1)){
    delta <- fitness.max - fitness.mean
    a <- (sfactor - 1.0)*fitness.mean/delta
    b <- fitness.mean * (fitness.max - sfactor*fitness.mean)/delta 
  } else {
    delta <- fitness.mean - fitness.min
    a <- fitness.mean/delta
    b <- -1*fitness.min*fitness.mean/delta 
  }
  fscaled <- a*fitness + b
  prob <- abs(fscaled)/sum(abs(fscaled), na.rm = TRUE)
  prob[is.na(prob)] <- eps
  prob <- pmin(pmax(0.0, prob/sum(prob)), 1.0)
  select <- sample(1:pop.size, size = pop.size, 
                   prob = prob, replace = TRUE)
  pop <- as.matrix(pop[select, ])
  return(pop)
}


# .. Algorithm ####
GA <- function(fitness.func, pop.size = 500, max.iter = 100,
               lower = NULL, upper = NULL, init.pop = NULL,
               selection = prop.linear.scaling, prob.crossover = 0.8, 
               prob.mutation = 0.1, percent.elites = 5,
               keep.track = FALSE, seed = NULL){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  # Initialize empty list
  GA.out <- list()
  pop.size <- pop.size
  
  # Take suggestions for initial population,
  # else sample from uniform
  if (!missing(init.pop)){
    pop <- as.matrix(init.pop)
    nvars <- ncol(pop)
  } else {
    nvars <- length(lower)
    pop <- matrix(NA, 
                  nrow = pop.size, ncol = nvars)
    for(j in 1:nvars){
      pop[, j] <- runif(pop.size, lower[j], upper[j]) 
      }
  }
  
  if (keep.track){
    pop.evolution <- NULL
    fitness.evolution <- NULL
    best.individual <- NULL
  }
  
  iter <- 0
  for (iter in 1:max.iter){
    
    # Evaluate fitness
    fitness <- rep(NA, pop.size)
    for(j in 1:nvars){
      for (i in 1:pop.size) {
        fitness[i] <- fitness.func(pop[i,])
      }
    }
    
    # Sort from best individuals to worst
    order.dec <- order(fitness, decreasing = T)
    pop.sorted.d <- as.matrix(pop[order.dec, ])
    fitness.sorted.d <- fitness[order.dec]
    
    # Selection
    pop <- selection(pop, fitness)
    
    # Mating pool
    nmating <- floor(pop.size/2)
    mating.pool <- matrix(sample(1:(2*nmating), 
                                 size = (2*nmating)), 
                          ncol = 2)

    for (i in 1:length(nmating)){
      if(prob.crossover > runif(1)){
        parents.id <- mating.pool[i,]
        parents <- pop[parents.id, ]
        if (is.null(ncol(parents))){
          n <- 1
        } else {
          n <- ncol(parents)
        }
        offspring <- matrix(NA, nrow = 2, ncol = n)
        p <- runif(n)
        parents <- as.matrix(parents)
        offspring[1,] <- p*parents[1,] + (1-p)*parents[2,]
        offspring[2,] <- p*parents[2,] + (1-p)*parents[1,]
        pop[parents.id,] <- offspring
      }
    }
    
    # Mutation
    for (i in 1:length(pop.size)){
      if(prob.mutation > runif(1)){
        mutant <- pop[i, ]
        n <- length(mutant)
        j <- sample(1:n, size = 1)
        mutant[j] <- runif(1, lower[j], upper[j])
        pop[i, ] <- mutant
      }
    }
    
    for(j in 1:nvars){
      for (i in 1:pop.size) {
        fitness[i] <- fitness.func(pop[i,])
      }
    }
    
    # Elites
    elites <- percent.elites/100*pop.size
    order.asc <- order(fitness, na.last = TRUE)
    unique.inds <- which(!duplicated(pop.sorted.d, margin = 1))
    pop[order.asc[1:elites],] <- pop.sorted.d[unique.inds[1:elites],]
    fitness[order.asc[1:elites]] <- fitness.sorted.d[unique.inds[1:elites]]
    
    if(keep.track){
      pop.evolution <- cbind(pop.evolution, pop)
      fitness.evolution <- cbind(fitness.evolution, fitness)
      best.individual <- rbind(best.individual,
                               unique(pop[which(fitness == max(fitness)),]))
    }
    
    cat("Iteration: ", iter, "\n")
    cat("Fitness value: ", max(fitness), "\n")
  }
  
  fitness.value <- max(fitness)
  solution <- unique(pop[which(fitness == fitness.value),])
  
  GA.out$population <- pop
  GA.out$fitness <- fitness
  GA.out$mean.population <- mean(pop)
  GA.out$mean.fitness <- mean(fitness)
  GA.out$solution <- solution
  GA.out$fitness.value <- fitness.value
  GA.out$pop.size <- pop.size
  
  if (keep.track){
    colnames(pop.evolution) <- paste0(rep(seq_len(max.iter), each = nvars), 
                                      "_V", rep(1:nvars))
    colnames(fitness.evolution) <- seq_len(max.iter)
    colnames(best.individual) <- paste0("V", rep(1:nvars))
    GA.out$pop.evolution <- pop.evolution
    GA.out$fitness.evolution <- fitness.evolution
    GA.out$best.individual <- best.individual
  }
    
  return(GA.out)
}


# .. Results ####
# rm(list = setdiff(ls(), c("GA", "claw", "rosenbrock", "rastrigin", 
#                           "roulette", "prop.linear.scaling",
#                           "prob.mutation", "prob.crossover",
#                           "percent.elites", "selection", 
#                           "plot.pop.evol", "plot.fitness.evol")))

GA3.claw <- GA(claw, pop.size = 500, max.iter = 100,
               lower = -10, upper = 10, selection = prop.linear.scaling,
               keep.track = T)
GA3.claw$solution
GA3.claw$fitness.value
plot.pop.evol(GA3.claw, filter = T, 
              gen.sequence = c(1, seq(10, 100, by = 10)))
plot.fitness.evol(GA3.claw, filter = T, 
                  gen.sequence = c(1, seq(10, 100, by = 10)))

GA3.rosenbrock <- GA(fitness.func = function(x) -rosenbrock(x[1], x[2]), 
                     pop.size = 500, max.iter = 100,
                     lower = c(-2, -1), upper = c(2, 3),
                     selection = prop.linear.scaling,
                     keep.track = T)
GA3.rosenbrock$solution
GA3.rosenbrock$fitness.value
plot.pop.evol(GA3.rosenbrock, filter = T, 
              gen.sequence = c(1, seq(10, 100, by = 10)))
plot.fitness.evol(GA3.rosenbrock, filter = T, 
                  gen.sequence = c(1, seq(10, 100, by = 10)))

GA3.rastrigin <- GA(fitness.func = function(x) -rastrigin(x[1], x[2]), 
                    pop.size = 500, max.iter = 100,
                    lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
                    keep.track = T)
GA3.rastrigin$solution
GA3.rastrigin$fitness.value
plot.pop.evol(GA3.rastrigin, filter = T, 
              gen.sequence = c(1, seq(10, 100, by = 10)))
plot.fitness.evol(GA3.rastrigin, filter = T, 
                  gen.sequence = c(1, seq(10, 100, by = 10)))

# x = seq(-3,3, by = 0.01)
# y <- NULL
# for (i in 1:length(x)) {y[i] <- claw(x[i])}
# claw <- data.frame(x = x, y = y)
# 
# best <- GA3.claw$best.individual
# ggplot(claw, aes(x = x, y = y)) +
#   geom_line() +
#   geom_point(x = claw$x[which(claw$y == max(claw$y))],
#              y = max(claw$y),
#              size = 2, col = "red") +
#   geom_point



## ## ## ## ## ## ## ## ## ## ##
# MATCHING                  ####
## ## ## ## ## ## ## ## ## ## ##

# .. Data ####

data("lalonde") # From Matching
attach(lalonde)
Y <- lalonde$re78
D <- lalonde$treat
X <- cbind(age, educ, black, hisp, married, 
           nodegr, re74, re75, u74, u75)


# .. Propensity Score Matching ####
glm <- glm(treat ~ age + educ + black + hisp
           + married + nodegr + re74 + re75, 
           family = binomial, data = lalonde)
match.ps <- Match(Y = Y, Tr = D, X = glm$fitted)
MatchBalance(D ~ age + I(age^2) + educ + I(educ^2) 
             + black + hisp + married + nodegr 
             + re74 + I(re74^2) + re75 + I(re75^2) 
             + u74 + u75 + I(re74 * re75)
             + I(age * nodegr) + I(educ * re74)
             + I(educ * re75),
             match.out = match.ps, 
             nboots = 1000, data = lalonde)


# .. Genetic Matching ####
BalanceMatrix <- cbind(age, I(age^2), educ, I(educ^2), 
                       black, hisp, married, nodegr, 
                       re74, I(re74^2), re75, I(re75^2), 
                       u74, u75, I(re74 * re75), 
                       I(age * nodegr), I(educ * re74), 
                       I(educ * re75))
GA.out <- GenMatch(Tr = D, X = X, 
                   #BalanceMatrix = BalanceMatrix, 
                   pop.size = 100,
                   loss = 2)
match.GA <- Match(Y = Y, Tr = D, X = X, 
                  Weight.matrix = GA.out)
# wmatrix <- GA.out$Weight.matrix
# match.GA <- Match(Y = Y, Tr = D, X = X, 
#                   Weight.matrix = wmatrix)
MatchBalance(D ~ age + I(age^2) + educ + I(educ^2) 
             + black + hisp + married + nodegr 
             + re74 + I(re74^2) + re75 + I(re75^2) 
             + u74 + u75 + I(re74 * re75)
             + I(age * nodegr) + I(educ * re74)
             + I(educ * re75),
             match.out = match.GA, 
             nboots = 1000, data = lalonde)



## ## ## ## ## ## ## ## ## ## ##
# GENETIC MATCHING          ####
## ## ## ## ## ## ## ## ## ## ##

# .. Objective function ####


# .. Results ####
GA.out <- GA(fitness.func = function(x) -matching(x), 
   pop.size = 100, max.iter = 2,
   lower = rep(1, 10), upper = rep(1000, 10),
   keep.track = T)
GA.out$solution
wmatrix <- diag(x = as.numeric(GA.out$solution))

match.GA <- Match(Y = lalonde$re78, Tr = lalonde$treat,
                  X = X, 
                  Weight.matrix = wmatrix)
mine <- MatchBalance(treat ~ age + educ + black + hisp
                     + married + nodegr + re74 + re75,
             match.out = match.GA, 
             nboots = 1000, data = lalonde)

treated.obs <- as.data.frame(match.GA$mdata$X[match.GA$index.treated, ])
control.obs <- as.data.frame(match.GA$mdata$X[match.GA$index.control, ])

qqplot(lalonde$educ[treat == 1], lalonde$educ[treat == 0])
qqplot(treated.obs$educ, control.obs$educ)
