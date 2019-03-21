# Genetic Algorithm with Applications to Matching
# Alice Lepissier
# alice.lepissier@gmail.com


## ## ## ## ## ## ## ## ## ## ##
# INDEX                     ####
## ## ## ## ## ## ## ## ## ## ##
# Preamble
# Test Functions
# .. Asymmetric Double Claw function
# .. Rosenbrock function
# .. Rastrigin function
# GA Packages
# .. Using genoud() function
# .. Using ga() function
# Genetic operators
# .. Roulette wheel selection
# .. Selection proportional to fitness with linear scaling
# .. Simple crossover
# .. Blended crossover
# .. Random uniform mutation
# .. Boundary mutation
# .. Gaussian mutation
# Genetic Algorithm
# Results
# Matching
# .. Data
# .. Propensity Score Matching
# .. Genetic Matching
# Genetic Matching
# .. Objective function
# .. Results



## ## ## ## ## ## ## ## ## ## ##
# PREAMBLE                  ####
## ## ## ## ## ## ## ## ## ## ##

setwd("C:/Users/Alice/Box Sync/PhD/Statistics/PSTAT 232/Final Project") # Laptop
#setwd("C:/boxsync/alepissier/PhD/Statistics/PSTAT 232/Final project") # Bren
#setwd("/home/alice/232/") # Linux server
#install.packages("ebal")
#install.packages("genalg")
#install.packages("GA")
#devtools::install_github('drizztxx/gatbxr')
#install.packages("kableExtra")
#install.packages("Matching")
#install.packages("MatchIt")
#install.packages("nor1mix")
#install.packages("rgenoud")
#install.packages("rgl")
#install.packages("tictoc")
#install.packages("VGAM")
library(ebal)
#library(genalg)
library(GA)
#library(gatbxr)
library(ggridges)
library(kableExtra)
library(imager)
library(Matching)
#library(MatchIt)
#library(nor1mix)
library(reshape2)
library(rgenoud)
library(rgl)
library(tictoc)
library(tidyverse)
library(VGAM) # For Gaussian error function



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
claw.data <- as.data.frame(cbind(x, y))
g <- ggplot(claw.data, aes(x = x, y = y)) +
  geom_line() + 
  geom_point(x = claw.data$x[which(claw.data$y == max(claw.data$y))],
             y = max(claw.data$y),
             col = "red", size = 2) +
  labs(title = "Asymmetric Double Claw function")
ggsave(g, file = "Figures/Asymmetric Double Claw function.png",
       width = 6, height = 4.5, units = "in")


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

# Manipulate viewing angles in interactive window and save snapshot
persp3d(x1, x2, y, col = "yellow")
# zoom <- par3d()$zoom
# userMatrix <- par3d()$userMatrix
# windowRect <- par3d()$windowRect
# open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
# persp3d(x1, x2, y, col = "yellow")
# rgl.snapshot("Figures/Rosenbrock function.png", fmt = "png")


# .. Rastrigin function ####
rastrigin <- function(x1, x2){
  f <- 20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
  return(f)
}
x1 <- x2 <- seq(-5.12, 5.12, by = 0.5)
y <- outer(x1, x2, rastrigin)

# Manipulate viewing angles in interactive window and save snapshot
persp3d(x1, x2, y, col = "green")
# zoom <- par3d()$zoom
# userMatrix <- par3d()$userMatrix
# windowRect <- par3d()$windowRect
# open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
# persp3d(x1, x2, y, col = "green")
# rgl.snapshot("Figures/Rastrigin function.png", fmt = "png")



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
# GENETIC OPERATORS         ####
## ## ## ## ## ## ## ## ## ## ##

# .. Roulette wheel selection ####
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


# .. Selection proportional to fitness with linear scaling ####
# Credit: Luca Scrucca
# I basically implemented his selection operator from the GA package.
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


# .. Simple crossover ####
simple.crossover <- function(parents, params){
  parents <- parents
  nvars <- ncol(parents)
  offspring <- matrix(NA, nrow = 2, ncol = nvars)
  p <- runif(nvars)
  offspring[1,] <- p*parents[1,] + (1-p)*parents[2,]
  offspring[2,] <- p*parents[2,] + (1-p)*parents[1,]
  
  return(offspring)
}


# .. Blended crossover ####
blend.crossover <- function(parents, params){
  parents <- parents
  params <- params
  nvars <- ncol(parents)
  offspring <- matrix(NA, nrow = 2, ncol = nvars)
  alpha <- 0.5
  x.l <- parents[1,] - alpha*(parents[2,] - parents[1,])
  x.u <- parents[2,] + alpha*(parents[2,] - parents[1,])
  for (v in 1:nvars){
    offspring[,v] <- runif(2, min(x.l[v], params$lower[v]),
                           max(x.u[v], params$upper[v]))
  }
  
  return(offspring)
}


# .. Random uniform mutation ####
unif.mutation <- function(mutant, params){
  mutant <- mutant
  params <- params
  n <- length(mutant)
  j <- sample(1:n, size = 1)
  mutant[j] <- runif(1, params$lower[j], params$upper[j])
  
  return(mutant)
}


# .. Boundary mutation ####
boundary.mutation <- function(mutant, params){
  mutant <- mutant
  params <- params
  n <- length(mutant)
  j <- sample(1:n, size = 1)
  p <- runif(1)
  if (p >= 0.5){
    mutant[j] <- params$upper[j]
  } else {
    mutant[j] <- params$lower[j]
  }
  
  return(mutant)
}


# .. Gaussian mutation ####
gaussian.mutation <- function(mutant, params){
  mutant <- mutant
  params <- params
  n <- length(mutant)
  j <- sample(1:n, size = 1)
  
  a <- params$lower[j]
  b <- params$upper[j]
  sigma <- 1/(b - a)
  u.l <- 0.5*(erf((a-mutant[j])/(sqrt(2)*(b-a)*sigma))+1)
  u.u <- 0.5*(erf((b-mutant[j])/(sqrt(2)*(b-a)*sigma))+1)
  p <- runif(1)
  if (p <= 0.5){
    u <- 2*u.l*(1-2*p)
  } else {
    u <- 2*u.u*(2*p-1)
  }
  if (u < -1 | u > 1){
    u <- 0 #erf is not defined ==> erf(0) means no mutation
  }
  mutant[j] <- mutant[j] + sqrt(2)*sigma*(b-a)*erf(u, inverse = T)
  
  return(mutant)
}



## ## ## ## ## ## ## ## ## ## ##
# GENETIC ALGORITHM         ####
## ## ## ## ## ## ## ## ## ## ##

GA <- function(fitness.func, pop.size = 500, max.iter = 100,
               lower = NULL, upper = NULL, init.pop = NULL,
               selection = prop.linear.scaling, 
               crossover = simple.crossover,
               mutation = unif.mutation,
               prob.crossover = 0.8, prob.mutation = 0.1, 
               percent.elites = 5,
               keep.track = FALSE, seed = NULL, verbose = T){
  
  # Initialize
  tic("Run-time")
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  GA.out <- list()
  params <- list()
  params$lower <- lower
  params$upper <- upper
  pop.size <- pop.size
  
  if (keep.track){
    pop.evolution <- NULL
    fitness.evolution <- NULL
    best.individual <- NULL
  }
  
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
  
  iter <- 0
  for (iter in 1:max.iter){
    
    # Evaluate fitness
    fitness <- rep(NA, pop.size)
    for (i in 1:pop.size) {
      fitness[i] <- fitness.func(pop[i,])
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
    
    # Crossover
    for (i in 1:length(nmating)){
      if(prob.crossover > runif(1)){
        parents.id <- mating.pool[i,]
        parents <- pop[parents.id, , drop = FALSE]
        offspring <- crossover(parents, params)
        pop[parents.id,] <- offspring
      }
    }
    
    # Mutation
    for (i in 1:length(pop.size)){
      if(prob.mutation > runif(1)){
        mutant <- pop[i, ]
        mutant <- mutation(mutant, params)
        pop[i, ] <- mutant
      }
    }
    
    # Evaluate fitness
    for (i in 1:pop.size) {
      fitness[i] <- fitness.func(pop[i,])
    }
    
    # Elitism
    # Replace worst individuals with % of best individuals
    # specified in elites
    elites <- percent.elites/100*pop.size
    order.asc <- order(fitness, na.last = TRUE)
    unique.inds <- which(!duplicated(pop.sorted.d, margin = 1))
    pop[order.asc[1:elites],] <- pop.sorted.d[unique.inds[1:elites],]
    fitness[order.asc[1:elites]] <- fitness.sorted.d[unique.inds[1:elites]]
    
    # Keep track of generations
    if(keep.track){
      pop.evolution <- cbind(pop.evolution, pop)
      fitness.evolution <- cbind(fitness.evolution, fitness)
      best.individual <- rbind(best.individual,
                               unique(pop[which(fitness == max(fitness)),]))
      if (verbose){
        cat("Iteration: ", iter, "\n",
            "Fitness value: ", max(fitness), "\n", 
            "Best individual: ", best.individual[iter,], "\n\n")
      }
    } else {
      if (verbose){
        cat("Iteration: ", iter, "\n",
            "Fitness value: ", max(fitness), "\n\n")
      }
    }
  }
  
  # Determine solution and fitness value
  fitness.value <- max(fitness)
  solution <- unique(pop[which(fitness == fitness.value),])
  if (nrow(solution > 1)){
    solution <- apply(solution, 2, mean)
  }
  
  # Spit out results
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
  if (verbose){
    toc()
  }
    
  return(GA.out)
}



## ## ## ## ## ## ## ## ## ## ##
# RESULTS                   ####
## ## ## ## ## ## ## ## ## ## ##

# rm(list = setdiff(ls(), c("GA", "claw", "rosenbrock", "rastrigin",
#                           "roulette", "prop.linear.scaling",
#                           "simple.crossover", "blend.crossover",
#                           "unif.mutation", "boundary.mutation",
#                           "gaussian.mutation",
#                           "prob.mutation", "prob.crossover",
#                           "percent.elites", "selection",
#                           "plot.pop.evol", "plot.fitness.evol")))


# .. Evolution of fitness value and solution ####

GA3.claw <- GA(claw, pop.size = 500, max.iter = 100,
               lower = -3, upper = 3,
               keep.track = T, seed = 232)
GA3.claw$solution
GA3.claw$fitness.value
g <- plot.pop.evol(GA3.claw, filter = T, 
              gen.sequence = c(1, 2, 3, 4, 5, 10, 30, 50, 70, 100))
ggsave(g, file = "Figures/Evolution population claw.pdf",
       width = 6, height = 4, units = "in")
g <- plot.fitness.evol(GA3.claw, filter = T, 
                  gen.sequence = c(1, 2, 3, 4, 5, 10, 30, 50, 70, 100))
ggsave(g, file = "Figures/Evolution fitness claw.pdf",
       width = 6, height = 4, units = "in")

GA3.rosenbrock <- GA(fitness.func = function(x) -rosenbrock(x[1], x[2]), 
                     pop.size = 500, max.iter = 100,
                     lower = c(-2, -1), upper = c(2, 3),
                     keep.track = T, seed = 232)
GA3.rosenbrock$solution
GA3.rosenbrock$fitness.value
g <- plot.pop.evol(GA3.rosenbrock, filter = T, 
              gen.sequence = c(1, seq(10, 100, by = 10)),
              solution = c(1, 1))
ggsave(g, file = "Figures/Evolution population Rosenbrock.pdf",
       width = 6, height = 4, units = "in")
g <- plot.fitness.evol(GA3.rosenbrock, filter = T, 
                  gen.sequence = c(1, seq(10, 100, by = 10)))
ggsave(g, file = "Figures/Evolution fitness Rosenbrock.pdf",
       width = 6, height = 4, units = "in")

GA3.rastrigin <- GA(fitness.func = function(x) -rastrigin(x[1], x[2]), 
                    pop.size = 500, max.iter = 100,
                    lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
                    keep.track = T, seed = 232)
GA3.rastrigin$solution
GA3.rastrigin$fitness.value
g <- plot.pop.evol(GA3.rastrigin, filter = T, 
              gen.sequence = c(1, 2, 3, 4, 5, 10, 30, 50, 70, 100),
              solution = c(0, 0))
ggsave(g, file = "Figures/Evolution population Rastrigin.pdf",
       width = 6, height = 4, units = "in")
g <- plot.fitness.evol(GA3.rastrigin, filter = T, 
                  gen.sequence = c(1, 2, 3, 4, 5, 10, 30, 50, 70, 100))
g <- ggsave(g, file = "Figures/Evolution fitness Rastrigin.pdf",
       width = 6, height = 4, units = "in")


# .. Simulations for different probabilities of mutation ####
nsims <- 10
p.mut <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
store.fitness.pmut <- matrix(NA, nrow = nsims, ncol = length(p.mut))
colnames(store.fitness.pmut) <- p.mut
for (p in 1:length(p.mut)){
  for (s in 1:nsims){
    store.fitness.pmut[s,p] <- GA(fitness.func = function(x) -rosenbrock(x[1], x[2]), 
                                  pop.size = 500, max.iter = 100,
                                  lower = c(-2, -1), upper = c(2, 3),
                                  seed = s, verbose = F,
                                  prob.crossover = 0.8, 
                                  prob.mutation = p.mut[p],
                                  percent.elites = 5,
                                  selection = prop.linear.scaling,
                                  crossover = simple.crossover,
                                  mutation = unif.mutation)$fitness.value
  }
}

store <- melt(as.data.frame(store.fitness.pmut))
g <- ggplot(store,
       aes(x = value, y = variable,
           fill = variable)) +
  geom_density_ridges() +
  labs(x = "Fitness value",
       y = "Probability of mutation",
       title = paste0("Fitness value for Rosenbrock function, ", nsims, " simulations"),
       subtitle = "Selection with fitness scaling\nSimple crossover (p = 0.8); Uniform mutation") +
  guides(fill = FALSE)
ggsave(g, file = "Figures/Simulation probability of mutation.pdf",
       width = 6, height = 4, units = "in")


# .. Simulations for different probabilities of crossover ####
nsims <- 10
p.cross <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
store.fitness.pcross <- matrix(NA, nrow = nsims, ncol = length(p.cross))
colnames(store.fitness.pcross) <- p.cross
for (p in 1:length(p.cross)){
  for (s in 1:nsims){
    store.fitness.pcross[s,p] <- GA(fitness.func = function(x) -rosenbrock(x[1], x[2]), 
                                    pop.size = 500, max.iter = 100,
                                    lower = c(-2, -1), upper = c(2, 3),
                                    seed = s, verbose = F,
                                    prob.crossover = p.cross[p], 
                                    prob.mutation = 0.2,
                                    percent.elites = 5,
                                    selection = prop.linear.scaling,
                                    crossover = simple.crossover,
                                    mutation = unif.mutation)$fitness.value
  }
}

store <- melt(as.data.frame(store.fitness.pcross))
g <- ggplot(store,
            aes(x = value, y = variable,
                fill = variable)) +
  geom_density_ridges() +
  labs(x = "Fitness value",
       y = "Probability of crossover",
       title = paste0("Fitness value for Rosenbrock function, ", nsims, " simulations"),
       subtitle = "Selection with fitness scaling\nSimple crossover; Uniform mutation (p = 0.2)") +
  guides(fill = FALSE)
ggsave(g, file = "Figures/Simulation probability of crossover.pdf",
       width = 6, height = 4, units = "in")


# .. Simulations for different percentages of elites ####
nsims <- 10
p.elites <- c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
store.fitness.pelites <- matrix(NA, nrow = nsims, ncol = length(p.elites))
colnames(store.fitness.pelites) <- p.elites
for (p in 1:length(p.elites)){
  for (s in 1:nsims){
    store.fitness.pelites[s,p] <- GA(fitness.func = function(x) -rosenbrock(x[1], x[2]), 
                                    pop.size = 500, max.iter = 100,
                                    lower = c(-2, -1), upper = c(2, 3),
                                    seed = s, verbose = F,
                                    prob.mutation = 0.2,
                                    prob.crossover = 0.8, 
                                    percent.elites = p.elites[p],
                                    selection = prop.linear.scaling,
                                    crossover = simple.crossover,
                                    mutation = unif.mutation)$fitness.value
  }
}

store <- melt(as.data.frame(store.fitness.pelites))
g <- ggplot(store,
            aes(x = value, y = variable,
                fill = variable)) +
  geom_density_ridges() +
  labs(x = "Fitness value",
       y = "Percentage of elites",
       title = paste0("Fitness value for Rosenbrock function, ", nsims, " simulations"),
       subtitle = "Selection with fitness scaling\nSimple crossover (p=0.8); Uniform mutation (p = 0.2)") +
  guides(fill = FALSE)
ggsave(g, file = "Figures/Simulation percentage of elites.pdf",
       width = 6, height = 4, units = "in")


# .. Simulations for different genetic operators ####
store.genop1 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = prop.linear.scaling,
                   crossover = simple.crossover,
                   mutation = unif.mutation)$fitness.evolution

store.genop2 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = prop.linear.scaling,
                   crossover = simple.crossover,
                   mutation = boundary.mutation)$fitness.evolution

store.genop3 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = prop.linear.scaling,
                   crossover = simple.crossover,
                   mutation = gaussian.mutation)$fitness.evolution

store.genop4 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = prop.linear.scaling,
                   crossover = blend.crossover,
                   mutation = unif.mutation)$fitness.evolution

store.genop5 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = prop.linear.scaling,
                   crossover = blend.crossover,
                   mutation = boundary.mutation)$fitness.evolution

store.genop6 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = prop.linear.scaling,
                   crossover = blend.crossover,
                   mutation = gaussian.mutation)$fitness.evolution

store.genop7 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = roulette,
                   crossover = simple.crossover,
                   mutation = unif.mutation)$fitness.evolution

store.genop8 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = roulette,
                   crossover = simple.crossover,
                   mutation = boundary.mutation)$fitness.evolution

store.genop9 <- GA(fitness.func = claw, 
                   pop.size = 500, max.iter = 100,
                   lower = -3, upper = 3,
                   seed = 232, keep.track = T,
                   prob.mutation = 0.2,
                   prob.crossover = 0.8, 
                   percent.elites = 10,
                   selection = roulette,
                   crossover = simple.crossover,
                   mutation = gaussian.mutation)$fitness.evolution

store.genop10 <- GA(fitness.func = claw, 
                    pop.size = 500, max.iter = 100,
                    lower = -3, upper = 3,
                    seed = 232, keep.track = T,
                    prob.mutation = 0.2,
                    prob.crossover = 0.8, 
                    percent.elites = 10,
                    selection = roulette,
                    crossover = blend.crossover,
                    mutation = unif.mutation)$fitness.evolution

store.genop11 <- GA(fitness.func = claw, 
                    pop.size = 500, max.iter = 100,
                    lower = -3, upper = 3,
                    seed = 232, keep.track = T,
                    prob.mutation = 0.2,
                    prob.crossover = 0.8, 
                    percent.elites = 10,
                    selection = roulette,
                    crossover = blend.crossover,
                    mutation = boundary.mutation)$fitness.evolution

store.genop12 <- GA(fitness.func = claw, 
                    pop.size = 500, max.iter = 100,
                    lower = -3, upper = 3,
                    seed = 232, keep.track = T,
                    prob.mutation = 0.2,
                    prob.crossover = 0.8, 
                    percent.elites = 10,
                    selection = roulette,
                    crossover = blend.crossover,
                    mutation = gaussian.mutation)$fitness.evolution

store.genop1 = apply(store.genop1, 2, mean)
store.genop2 = apply(store.genop2, 2, mean)
store.genop3 = apply(store.genop3, 2, mean)
store.genop4 = apply(store.genop4, 2, mean)
store.genop5 = apply(store.genop5, 2, mean)
store.genop6 = apply(store.genop6, 2, mean)
store.genop7 = apply(store.genop7, 2, mean)
store.genop8 = apply(store.genop8, 2, mean)
store.genop9 = apply(store.genop9, 2, mean)
store.genop10 = apply(store.genop10, 2, mean)
store.genop11 = apply(store.genop11, 2, mean)
store.genop12 = apply(store.genop12, 2, mean)

store.genop <- cbind(store.genop1, store.genop2, store.genop3,
      store.genop4, store.genop5, store.genop6,
      store.genop7, store.genop8, store.genop9,
      store.genop10, store.genop11, store.genop12)
store.genop <- as.data.frame(store.genop) %>% 
  melt %>%
  mutate(Fitness = log(value),
         Generation = rep(1:100, 12))

g <- ggplot(store.genop,
       aes(x = Generation, y = Fitness, col = variable)) +
  geom_line(aes(y = log(0.4113089)), 
            col = "red", linetype = 2) +
  geom_line(size = 1) +
  labs(y = "Mean population fitness (log)",
       title = "Performance of genetic operators on ADC") +
  scale_color_discrete(labels = c("LS + SC + UM", 
                                  "LS + SC + BM", 
                                  "LS + SC + GM", 
                                  "LS + BC + UM", 
                                  "LS + BC + BM", 
                                  "LS + BC + GM", 
                                  "RS + SC + UM", 
                                  "RS + SC + BM", 
                                  "RS + SC + GM", 
                                  "RS + BC + UM", 
                                  "RS + BC + BM", 
                                  "RS + BC + GM")) +
  theme(legend.title = element_blank()) +
  ylim(-1,-0.875)
ggsave(g, file = "Figures/Simulations genetic operators.pdf",
       width = 6, height = 4, units = "in")



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
BalanceMatrix <- cbind(age, I(age^2), educ, I(educ^2), 
                       black, hisp, married, nodegr, 
                       re74, I(re74^2), re75, I(re75^2), 
                       u74, u75, I(re74 * re75), 
                       I(age * nodegr), I(educ * re74), 
                       I(educ * re75))


# .. Balance before matching ####
no.bal <- MatchBalance(D ~ age + I(age^2) + educ + I(educ^2) 
                       + black + hisp + married + nodegr 
                       + re74 + I(re74^2) + re75 + I(re75^2) 
                       + u74 + u75 + I(re74 * re75)
                       + I(age * nodegr) + I(educ * re74)
                       + I(educ * re75),
                       nboots = 1000, data = lalonde)
# varnames <- c("Age", "Age sq.", "Years of education", "Years education sq.",
#               "Black", "Hispanic", "Married", "No degree", "1974 earnings", "1974 earnings sq.",
#               "1975 earnings", "1975 earnings sq.", "Unemployed in 1974", "Unemployed in 1975",
#               "1974 * 1975 earnings", "age * educ", "educ * 1974 earnings", "educ * 1975 earnings")
varnames <- c("age", "age^2", "educ", "educ^2", "black", "hispanic", "married", "no degree", 
              "re1974", "re1974^2", "re1975", "re1975^2", "unemp74", "unemp75", "re1974*re1975", 
              "age*educ", "educ*re1974", "educ*1975")
baltest <- baltest.collect(no.bal,
                           var.names = varnames,
                           after = F)
bal.out.before <- round(baltest[, c("mean.Tr", "mean.Co", "T pval", "KS pval")], 2)
colnames(bal.out) <- c("Mean Treated", "Mean Control", "t p-value", "KS p-value")
kable(bal.out, format = "latex")


# .. Propensity Score Matching ####
glm <- glm(treat ~ age + educ + black + hisp
           + married + nodegr + re74 + re75, 
           family = binomial, data = lalonde)
match.ps <- Match(Y = Y, Tr = D, X = glm$fitted)
ps.bal <- MatchBalance(D ~ age + I(age^2) + educ + I(educ^2) 
             + black + hisp + married + nodegr 
             + re74 + I(re74^2) + re75 + I(re75^2) 
             + u74 + u75 + I(re74 * re75)
             + I(age * nodegr) + I(educ * re74)
             + I(educ * re75),
             match.out = match.ps, 
             nboots = 1000, data = lalonde)
baltest <- baltest.collect(ps.bal,
                           var.names = varnames,
                           after = T)
bal.out.ps <- round(baltest[, c("mean.Tr", "mean.Co", "T pval", "KS pval")], 2)
colnames(bal.out) <- c("Mean Treated", "Mean Control", "t p-value", "KS p-value")
kable(bal.out, format = "latex")


# .. Genetic Matching ####
Genoud.out <- GenMatch(Tr = D, X = X, 
                   BalanceMatrix = BalanceMatrix, 
                   pop.size = 100,
                   loss = 2)
match.Genoud <- Match(Y = Y, Tr = D, X = X, 
                      Weight.matrix = Genoud.out)
genoud.bal <- MatchBalance(D ~ age + I(age^2) + educ + I(educ^2) 
             + black + hisp + married + nodegr 
             + re74 + I(re74^2) + re75 + I(re75^2) 
             + u74 + u75 + I(re74 * re75)
             + I(age * nodegr) + I(educ * re74)
             + I(educ * re75),
             match.out = match.Genoud, 
             nboots = 1000, data = lalonde)
baltest <- baltest.collect(genoud.bal,
                           var.names = varnames,
                           after = T)
bal.out <- round(baltest[, c("mean.Tr", "mean.Co", "T pval", "KS pval")], 2)
colnames(bal.out) <- c("Mean Treated", "Mean Control", "t p-value", "KS p-value")
kable(bal.out, format = "latex")


# All together
bal.out <- cbind(bal.out.before, bal.out.ps, bal.out.GA.pval)
bal.out <- bal.out[,c(3:4,7:8,11:12)]
#colnames(bal.out) <- c("Before matching", "Matching on PS", "Matching using GA")
colnames(bal.out) <- rep(c("t", "KS"), 3)
kable(bal.out, format = "latex")



## ## ## ## ## ## ## ## ## ## ##
# GENETIC MATCHING          ####
## ## ## ## ## ## ## ## ## ## ##

# .. Objective functions ####
minpval <- function(x){
  
  wmatrix <- diag(x, nrow = nvars)
  
  # To counter the error that leading minors are not positive.
  # This happens when the eigenvectors are not positive, and the data is
  # too noisy to estimate a full covariance matrix. So add a penalty
  # term to push away from zero.
  eigenvalues <- eigen(wmatrix, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < sqrt(.Machine$double.eps)){
    wmatrix <- wmatrix + diag(nvars) * sqrt(.Machine$double.eps)
  } 
  
  match.init <- Match(Tr = D, X = X,
                      Weight.matrix = wmatrix)
  
  index.treated <- match.init$index.treated
  index.control <- match.init$index.control
  Tr <- BalanceVars[index.treated,]
  Co <- BalanceVars[index.control,]
  
  storage <- NULL
  for (v in 1:nbalvars){
    storage[v] <- t.test(Tr[,v], Co[,v], paired = T)$p.value
  }
  for (v in (nbalvars+1):(nbalvars*2)){
    storage[v] <- suppressWarnings(ks.test(Tr, Co)$p.value)
  }
  
  loss.func <- min(storage, na.rm = T)
  
  return(loss.func)
}

maxQQdif <- function(x){
  
  wmatrix <- diag(x, nrow = nvars)
  
  # To counter the error that leading minors are not positive.
  # This happens when the eigenvectors are not positive, and the data is
  # too noisy to estimate a full covariance matrix. So add a penalty
  # term to push away from zero.
  eigenvalues <- eigen(wmatrix, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < sqrt(.Machine$double.eps)){
    wmatrix <- wmatrix + diag(nvars) * sqrt(.Machine$double.eps)
  } 
  
  match.init <- Match(Tr = D, X = X,
                      Weight.matrix = wmatrix)
  
  index.treated <- match.init$index.treated
  index.control <- match.init$index.control
  Tr <- BalanceVars[index.treated,]
  Co <- BalanceVars[index.control,]
  
  storage <- NULL
  for (v in 1:nbalvars){
    storage[v] <- qqstats(Tr, Co, standardize = T)$mediandiff
  }
  
  loss.func <- max(storage, na.rm = T)
  
  return(loss.func)
}

# .. Variables to seek balance on ####
PropensityScore <- glm(treat ~ age + educ + black + hisp
                       + married + nodegr + re74 + re75, 
                       family = binomial, data = lalonde)$fitted
BalanceVars <- cbind(BalanceMatrix, PropensityScore)
nvars <- ncol(X)
nbalvars <- ncol(BalanceVars)


# .. Run genetic algorithm to figure out optimal weights ####
GA.out.pval <- GA(fitness.func = minpval, 
                  pop.size = 50, max.iter = 1000,
                  lower = rep(1, nbalvars), upper = rep(1000, nbalvars),
                  percent.elites = 25,
                  prob.mutation = 0.5,
                  keep.track = F, seed = 232)

GA.out.QQ <- GA(fitness.func = function(x) -maxQQdif(x), 
                pop.size = 50, max.iter = 500,
                lower = rep(1, nbalvars), upper = rep(1000, nbalvars),
                mutation = boundary.mutation,
                keep.track = T, seed = 232)

save(GA.out.pval, file = "Results/GA.out.pval.RData")
save(GA.out.QQ, file = "Results/GA.out.QQ.RData")


# .. Perform matching with optimal weights ####
wmatrix.pval <- diag(x = as.numeric(GA.out.pval$solution))
wmatrix.QQ <- diag(x = as.numeric(GA.out.QQ$solution))

match.GA.pval <- Match(Y = Y, Tr = D, X = X, 
                       Weight.matrix = wmatrix.pval,
                       BiasAdjust = T)
match.GA.QQ <- Match(Y = Y, Tr = D, X = X, 
                     Weight.matrix = wmatrix.QQ,
                     BiasAdjust = T)


# .. Evaluate balance ####
GA.pval.bal <- MatchBalance(D ~ age + I(age^2) + educ + I(educ^2) 
             + black + hisp + married + nodegr 
             + re74 + I(re74^2) + re75 + I(re75^2) 
             + u74 + u75 + I(re74 * re75)
             + I(age * nodegr) + I(educ * re74)
             + I(educ * re75),
             match.out = match.GA.pval, 
             nboots = 1000, data = lalonde)
baltest <- baltest.collect(GA.pval.bal,
                           var.names = varnames,
                           after = T)
bal.out.GA.pval <- round(baltest[, c("mean.Tr", "mean.Co", "T pval", "KS pval")], 2)
colnames(bal.out.GA.pval) <- c("Mean Treated", "Mean Control", "t p-value", "KS p-value")
kable(bal.out, format = "latex")

GA.QQ.bal <- MatchBalance(D ~ age + I(age^2) + educ + I(educ^2) 
             + black + hisp + married + nodegr 
             + re74 + I(re74^2) + re75 + I(re75^2) 
             + u74 + u75 + I(re74 * re75)
             + I(age * nodegr) + I(educ * re74)
             + I(educ * re75),
             match.out = match.GA.QQ, 
             nboots = 1000, data = lalonde)
baltest <- baltest.collect(GA.QQ.bal,
                           var.names = varnames,
                           after = T)
bal.out.GA.QQ <- round(baltest[, c("mean.Tr", "mean.Co", "T pval", "KS pval")], 2)
colnames(bal.out.GA.QQ) <- c("Mean Treated", "Mean Control", "t p-value", "KS p-value")
kable(bal.out, format = "latex")


# .. Plot results ####
# treated.obs.pval <- as.data.frame(match.GA.pval$mdata$X[match.GA.pval$index.treated, ])
# control.obs.pval <- as.data.frame(match.GA.pval$mdata$X[match.GA.pval$index.control, ])

treated.obs.pval <- lalonde[match.GA.pval$index.treated,]
control.obs.pval <- lalonde[match.GA.pval$index.control,]

treated.obs.QQ <- as.data.frame(match.GA.QQ$mdata$X[match.GA.QQ$index.treated, ])
control.obs.QQ <- as.data.frame(match.GA.QQ$mdata$X[match.GA.QQ$index.control, ])

for (v in 1:nvars){
  pdf(paste0("Figures/minpval_Empirical QQ plot_", colnames(lalonde)[v], ".pdf"), 
      height = 4, width = 6)
  par(mfrow = c(1,2), oma = c(0, 0, 2, 0))
  qqplot(lalonde[treat == 0, v], lalonde[treat == 1, v],
         main = "Before matching",
         xlab = "Control Observations",
         ylab = "Treated Observations")
  abline(0,1, lty = 2)
  qqplot(control.obs.pval[, v], treated.obs.pval[, v],
         main = "After matching",
         xlab = "Control Observations",
         ylab = "Treated Observations")
  abline(0,1, lty = 2)
  title(paste0("Empirical Q-Q plot of ", colnames(lalonde)[v]), outer = T)
  dev.off()
}

for (v in 1:nvars){
  pdf(paste0("Figures/maxQQdif_Empirical QQ plot_", colnames(lalonde)[v], ".pdf"), 
      height = 4, width = 6)
  par(mfrow = c(1,2), oma = c(0, 0, 2, 0))
  qqplot(lalonde$educ[treat == 0], lalonde$educ[treat == 1],
         main = "Before matching",
         xlab = "Control Observations",
         ylab = "Treated Observations")
  abline(0,1, lty = 2)
  qqplot(control.obs.QQ$educ, treated.obs.QQ$educ,
         main = "After matching",
         xlab = "Control Observations",
         ylab = "Treated Observations")
  abline(0,1, lty = 2)
  title(paste0("Empirical Q-Q plot of ", colnames(lalonde)[v]), outer = T)
  dev.off()
}

