## ## ## ## ## ## ## ## ## ## ##
# DIGIT RECONSTRUCTION      ####
## ## ## ## ## ## ## ## ## ## ##

# .. Data pre-processing ####

# Raw data from https://www.kaggle.com/c/digit-recognizer/data
images.raw <- data.matrix(read.csv("train.csv", header=TRUE))
dim(images.raw)

# Choose digits from 0 to 9
digits <- images.raw[,-1]
dim(digits)

head(which(labels == 0))
digit0 <- matrix(digits[2,], nrow = 1)

head(which(labels == 1))
digit1 <- matrix(digits[3,], nrow = 1)

head(which(labels == 2))
digit2 <- matrix(digits[17,], nrow = 1)

head(which(labels == 3))
digit3 <- matrix(digits[10,], nrow = 1)

head(which(labels == 4))
digit4 <- matrix(digits[40,], nrow = 1)

head(which(labels == 5))
digit5 <- matrix(digits[52,], nrow = 1)

head(which(labels == 6))
digit6 <- matrix(digits[46,], nrow = 1)

head(which(labels == 7))
digit7 <- matrix(digits[30,], nrow = 1)

head(which(labels == 8))
digit8 <- matrix(digits[11,], nrow = 1)

head(which(labels == 9))
digit9 <- matrix(digits[12,], nrow = 1)


# .. Loss function ####
MSE.digit <- function(x){
  x <- matrix(x, nrow = 1)
  MSE <- sum((digit-x)^2)/784
  return(MSE)
}


# .. Run genetic algorithm ####

digit <- digit0
GA.digit0 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit0 <- GA.digit0$solution

digit <- digit1
GA.digit1 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit1 <- GA.digit1$solution

digit <- digit2
GA.digit2 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit2 <- GA.digit2$solution

digit <- digit3
GA.digit3 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit3 <- GA.digit3$solution

digit <- digit4
GA.digit4 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit4 <- GA.digit4$solution

digit <- digit5
GA.digit5 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit5 <- GA.digit5$solution

digit <- digit6
GA.digit6 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit6 <- GA.digit6$solution

digit <- digit7
GA.digit7 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit7 <- GA.digit7$solution

digit <- digit8
GA.digit8 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit8 <- GA.digit8$solution

digit <- digit9
GA.digit9 <- GA(fitness.func = function(x) -MSE.digit(x),
                pop.size = 50, max.iter = 10000,
                lower = rep(0, 784), upper = rep(255, 784),
                crossover = blend.crossover, 
                mutation = unif.mutation,
                prob.mutation = 0.5,
                percent.elites = 50, 
                seed = 232, keep = F, verbose = T)
reconstructed.digit9 <- GA.digit9$solution


# .. Render results ####
digits <- rbind(digit0, digit1, digit2, digit3, digit4,
                digit5, digit6, digit7, digit8, digit9)
recons.dig <- rbind(reconstructed.digit0, reconstructed.digit1,
                    reconstructed.digit2, reconstructed.digit3,
                    reconstructed.digit4, reconstructed.digit5,
                    reconstructed.digit6, reconstructed.digit7,
                    reconstructed.digit8, reconstructed.digit9)

pdf("Figures/Digit reconstruction.pdf", 
    height = 4, width = 6)
par(mfrow = c(2,10), mai = c(0,0,0,0))
for (d in 1:10){
  matrix(digits[d,], nrow=28, byrow=TRUE) %>% 
    apply(2, rev) %>% t %>%
    image(col=gray(1:0), axes=FALSE, asp=1)
}
for (d in 1:10){
  matrix(recons.dig[d,], nrow=28, byrow=TRUE) %>% 
    apply(2, rev) %>% t %>%
    image(col=gray(1:0), axes=FALSE, asp=1)
}
mtext("Real vs. reconstructed digits", at = -4, cex = 1.5)
dev.off()



## ## ## ## ## ## ## ## ## ## ##
# FACE RECONSTRUCTION       ####
## ## ## ## ## ## ## ## ## ## ##

# .. Data pre-processing ####

# Load faces array
load("faces_array.RData")

# Define plot face function
plot_face <- function(image, title = NULL) {
  plot(as.cimg(t(matrix(image, ncol = 100))), 
       axes = FALSE, asp = 1, main = title)
}
plot_facevec <- function(image_vector, title = NULL) {
  plot(as.cimg(t(matrix(image_vector, ncol = 100, byrow = T))), 
       axes = FALSE, asp = 1, main = title)
}

# Vectorize face
vectorize <- function(face){
  facevec <- as.vector(t(face))
  facevec <- round(facevec*255)
  return(facevec)
}


# .. Loss function ####
MSE.face <- function(x){
  MSE <- sum((facevec-x)^2)/2000
  return(MSE)
}


# .. Reconstruct Abraham Lincoln ####

# Choose face
abe <- faces_array[,,53]
plot_face(abe)

facevec <- vectorize(abe)
plot_facevec(facevec)

# Truncate face
facevec1 <- facevec[1:2000]
facevec2 <- facevec[2001:4000]
facevec3 <- facevec[4001:6000]
facevec4 <- facevec[6001:8000]
facevec5 <- facevec[8001:10000]
plot_facevec(facevec1)
plot_facevec(facevec2)
plot_facevec(facevec3)
plot_facevec(facevec4)
plot_facevec(facevec5)

# Reconstruct
reconstructed.face2k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 2000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face2k <- c(reconstructed.face2k, facechunk)
  
}
plot_facevec(reconstructed.face2k)

reconstructed.face10k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 10000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face10k <- c(reconstructed.face10k, facechunk)
  
}
plot_facevec(reconstructed.face10k)

reconstructed.face50k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 50000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face50k <- c(reconstructed.face50k, facechunk)
  
}
plot_facevec(reconstructed.face50k)

pdf("Figures/Abraham Lincoln.pdf", 
    height = 4, width = 6)
par(mfrow = c(1,4), mar = c(0,0,0,0), oma = c(0,0,2,0), mgp = c(0,0,0))
plot_face(abe)
title("Original", line = -6)
plot_facevec(reconstructed.face2k)
title("2K generations", line = -6)
plot_facevec(reconstructed.face10k)
title("10K generations", line = -6)
plot_facevec(reconstructed.face50k)
title("50K generations", line = -6)
title("Abraham Lincoln", outer = T, line = -3, cex.main = 2)
dev.off()


# .. Reconstruct Albert Einstein ####

# Choose face
einstein <- faces_array[,,203]
plot_face(einstein)

facevec <- vectorize(einstein)
plot_facevec(facevec)

# Truncate face
facevec1 <- facevec[1:2000]
facevec2 <- facevec[2001:4000]
facevec3 <- facevec[4001:6000]
facevec4 <- facevec[6001:8000]
facevec5 <- facevec[8001:10000]
plot_facevec(facevec1)
plot_facevec(facevec2)
plot_facevec(facevec3)
plot_facevec(facevec4)
plot_facevec(facevec5)

# Reconstruct
reconstructed.face2k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 1000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face2k <- c(reconstructed.face2k, facechunk)
  
}
plot_facevec(reconstructed.face2k)

reconstructed.face10k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 10000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face10k <- c(reconstructed.face10k, facechunk)
  
}
plot_facevec(reconstructed.face10k)

reconstructed.face50k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 50000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face50k <- c(reconstructed.face50k, facechunk)
  
}
plot_facevec(reconstructed.face50k)

pdf("Figures/Albert Einstein.pdf", 
    height = 4, width = 6)
par(mfrow = c(1,4), mar = c(0,0,0,0), oma = c(0,0,2,0), mgp = c(0,0,0))
plot_face(einstein)
title("Original", line = -6)
plot_facevec(reconstructed.face2k)
title("2K generations", line = -6)
plot_facevec(reconstructed.face10k)
title("10K generations", line = -6)
plot_facevec(reconstructed.face50k)
title("50K generations", line = -6)
title("Albert Einstein", outer = T, line = -3, cex.main = 2)
dev.off()


# .. Reconstruct face with less contrast ####

# Choose face
lowcon <- faces_array[,,222]
plot_face(lowcon)

facevec <- vectorize(lowcon)
plot_facevec(facevec)

# Truncate face
facevec1 <- facevec[1:2000]
facevec2 <- facevec[2001:4000]
facevec3 <- facevec[4001:6000]
facevec4 <- facevec[6001:8000]
facevec5 <- facevec[8001:10000]
plot_facevec(facevec1)
plot_facevec(facevec2)
plot_facevec(facevec3)
plot_facevec(facevec4)
plot_facevec(facevec5)

# Reconstruct
reconstructed.face2k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 1000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face2k <- c(reconstructed.face2k, facechunk)
  
}
plot_facevec(reconstructed.face2k)

reconstructed.face10k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 10000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face10k <- c(reconstructed.face10k, facechunk)
  
}
plot_facevec(reconstructed.face10k)

reconstructed.face50k <- NULL
for (f in 1:5){
  facevec <- eval(as.name(paste0("facevec", f)))
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 50000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$solution
  reconstructed.face50k <- c(reconstructed.face50k, facechunk)
  
}
plot_facevec(reconstructed.face50k)

pdf("Figures/Albio Sires.pdf", 
    height = 4, width = 6)
par(mfrow = c(1,4), mar = c(0,0,0,0), oma = c(0,0,2,0), mgp = c(0,0,0))
plot_face(lowcon)
title("Original", line = -6)
plot_facevec(reconstructed.face2k)
title("2K generations", line = -6)
plot_facevec(reconstructed.face10k)
title("10K generations", line = -6)
plot_facevec(reconstructed.face50k)
title("50K generations", line = -6)
title("Albio Sires", outer = T, line = -3, cex.main = 2)
dev.off()


# .. Simulations ####

facevec <- vectorize(abe)
facevec <- facevec[4001:6000]
pdf("Figures/Lincoln face chunk.pdf",
    height = 4, width = 6)
plot_facevec(facevec)
dev.off()

pm <- seq(0.0, 1, by = 0.1)
MSE.pm <- NULL
for (p in 1:length(pm)){
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 5000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = pm[p],
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$fitness.value
  MSE.pm <- c(MSE.pm, facechunk)
}
save(MSE.pm, file = "Results/MSE.pm.RData")
MSE.pm <- data.frame(PM = pm,
                     MSE = MSE.pm*-1)

g <- ggplot(MSE.pm, aes(x = PM, y = MSE)) + 
  geom_line() +
  scale_x_continuous(breaks = pm) +
  labs(x = "Probability of mutation",
       y = "Mean Squared Error",
       title = "Performance on Lincoln face chunk",
       subtitle = "5,000 generations; Population size 10; Selection with fitness scaling\nSimple crossover (p=0.8); Uniform mutation; 50% elites")
ggsave(g, file = "Figures/MSE probability of mutation.pdf",
       width = 6, height = 4, units = "in")

pop <- seq(10, 1000, by = 50)
MSE.pop <- NULL
for (p in 1:length(pop)){
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = pop[p], max.iter = 5000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = 50,
                  seed = 232, keep = F, verbose = F)$fitness.value
  MSE.pop <- c(MSE.pop, facechunk)
}
save(MSE.pop, file = "Results/MSE.pop.RData")
MSE.pop <- data.frame(Pop = pop,
                      MSE = MSE.pop*-1)

g <- ggplot(MSE.pop, aes(x = Pop, y = MSE)) + 
  geom_line() +
  scale_x_continuous(breaks = c(10, seq(100, 1000, by = 100))) +
  labs(x = "Population size",
       y = "Mean Squared Error",
       title = "Performance on Lincoln face chunk",
       subtitle = "5,000 generations; Selection with fitness scaling\nSimple crossover (p=0.8); Uniform mutation (p=1); 50% elites") +
  theme(axis.text.x = element_text(hjust = 0.7))
ggsave(g, file = "Figures/MSE population size.pdf",
       width = 6, height = 4, units = "in")

pe <- seq(5, 100, by = 5)
MSE.pe <- NULL
for (p in 1:length(pe)){
  facechunk <- GA(fitness.func = function(x) -MSE.face(x),
                  pop.size = 10, max.iter = 5000,
                  lower = rep(0, 2000), upper = rep(255, 2000),
                  crossover = simple.crossover, 
                  mutation = unif.mutation,
                  prob.mutation = 1,
                  percent.elites = pe[p],
                  seed = 232, keep = F, verbose = F)$fitness.value
  MSE.pe <- c(MSE.pe, facechunk)
}
save(MSE.pe, file = "Results/MSE.pe.RData")
MSE.pe <- data.frame(PE = pe,
                     MSE = MSE.pe*-1)

g <- ggplot(MSE.pe, aes(x = PE, y = MSE)) + 
  geom_line() +
  scale_x_continuous(breaks = pe) +
  labs(x = "Percentage elites",
       y = "Mean Squared Error",
       title = "Performance on Lincoln face chunk",
       subtitle = "5,000 generations; Population size 10; Selection with fitness scaling\nSimple crossover (p=0.8); Uniform mutation (p=1)") +
  theme(axis.text.x = element_text(hjust = 0.7))
ggsave(g, file = "Figures/MSE percentage elites.pdf",
       width = 6, height = 4, units = "in")