## ## ## ## ## ## ## ## ## ## ##
# DIGIT RECONSTRUCTION      ####
## ## ## ## ## ## ## ## ## ## ##

# .. Data pre-processing ####

# Raw data from https://www.kaggle.com/c/digit-recognizer/data
images.raw <- data.matrix(read.csv("train.csv", header=TRUE))
dim(images.raw)

# Plotting function
plot_digits = function(img, lab=NULL, figcols = 10, colmap=gray(1:0)){
  
  op = par(no.readonly=TRUE)
  
  figrows = ceiling(nrow(img)/figcols)
  par(mfrow=c(figrows, figcols), mar=rep(.1, 4))
  
  for (i in 1:nrow(img)){ #reverse and transpose each matrix to rotate img
    matrix(img[i,], nrow=28, byrow=TRUE) %>% 
      apply(2, rev) %>% t %>%
      image(col=colmap, axes=FALSE, asp=1)
    box()
    if (!is.null(lab))
      text(0.1, 0.13, col="blue", cex=2, lab[i])
  }
  par(op) #reset the original graphics parameters
}

# Extract true labels
labels <- images.raw[, 1]

# Convert to binary values
images.bin <- 1L*((images.raw[,-1]/255) > 0.5)

# Plot digits
plot_digits(images.bin[1:20,], labels[1:20])
plot_digits(images.raw[1:20,-1], lab = NULL, colmap = rev(gray.colors(255)))

# Choose digits from 0 to 9
digits <- images.raw[,-1]
dim(digits)

head(which(labels == 0))
digit0 <- matrix(digits[2,], nrow = 1)
matrix(digit0, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

head(which(labels == 1))
digit1 <- matrix(digits[3,], nrow = 1)
matrix(digit1, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)

head(which(labels == 2))
digit2 <- matrix(digits[17,], nrow = 1)
matrix(digit2, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)

head(which(labels == 3))
digit3 <- matrix(digits[10,], nrow = 1)
matrix(digit3, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)

head(which(labels == 4))
digit4 <- matrix(digits[40,], nrow = 1)
matrix(digit4, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)

head(which(labels == 5))
digit5 <- matrix(digits[52,], nrow = 1)
matrix(digit5, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)

head(which(labels == 6))
digit6 <- matrix(digits[46,], nrow = 1)
matrix(digit6, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)

head(which(labels == 7))
digit7 <- matrix(digits[30,], nrow = 1)
matrix(digit7, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)

head(which(labels == 8))
digit8 <- matrix(digits[11,], nrow = 1)
matrix(digit8, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)

head(which(labels == 9))
digit9 <- matrix(digits[12,], nrow = 1)
matrix(digit9, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)


# .. Loss function ####
MSE.digit <- function(x){
  x <- matrix(x, nrow = 1)
  MSE <- sum((digit-x)^2)/784
  return(MSE)
}


# .. Run genetic algorithm ####

# Using GA package (to check how mine compares)
digit <- digit0
GA2.digit <- ga(type = "real-valued", 
                fitness = function(x) -MSE.digit(x), 
                lower = rep(0, 784), upper = rep(255, 784),
                seed = 232, popSize = 50, maxiter = 500)
reconstructed.digit2 <- GA2.digit@solution
matrix(reconstructed.digit2, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit2, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

# Using my GA
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
matrix(reconstructed.digit0, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit0, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit1, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit1, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit2, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit2, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit3, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit3, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit4, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit4, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit5, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit5, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit6, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit6, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit7, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit7, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit8, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit8, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

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
matrix(reconstructed.digit9, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
dev.off()
matrix(reconstructed.digit9, nrow=28, byrow=TRUE) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)


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
# This doesn't work. I think it is a combination of the fact that there
# are 10,000 features (an image is 100*100), and the encoding is between
# 0 and 1 rather than on a RGB 0-255 scale.
# I have tried rescaling the images to 0-255 and that doesn't help either.

# Load faces array
load("faces_array.RData")

# Vectorize face pixels
face_mat <- sapply(1:1000, function(i) as.numeric(faces_array[, , i])) %>% t

# Define plot face function
plot_face <- function(image_vector) {
  plot(as.cimg(t(matrix(image_vector, ncol=100))), axes=FALSE, asp=1)
}

# Choose face 222
face <- matrix(face_mat[222,], nrow = 1)
dim(face)
plot_face(face)


# .. Loss function ####
MSE.face <- function(x){
  x <- matrix(x, nrow = 1)
  MSE <- sum((face-x)^2)/10000
  return(MSE)
}


# .. Genetic algorithm using 0-1 encoding ####

# Using GA package
GA.face <- ga(type = "real-valued", 
              fitness = function(x) -MSE.face(x), 
              lower = rep(0, 10000), upper = rep(1, 10000),
              seed = 232, popSize = 50, maxiter = 1000)
reconstructed.face <- GA.face@solution
plot_face(reconstructed.face)

# Using my GA
GA.face <- GA(fitness.func = function(x) -MSE.face(x),
              pop.size = 50, max.iter = 1000,
              lower = rep(0, 10000), upper = rep(1, 10000),
              crossover = blend.crossover, 
              mutation = unif.mutation,
              prob.mutation = 0.7,
              percent.elites = 70,
              seed = 232, keep = F, verbose = T)
reconstructed.face <- GA.face$solution
dim(reconstructed.face)
plot_face(reconstructed.face)


# .. Genetic algorithm using 0-255 encoding ####
face <- face*255
matrix(face, nrow=100, byrow=F) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
matrix(face, nrow=100, byrow=F) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

GA.face <- ga(type = "real-valued", 
               fitness = function(x) -MSE.face(x), 
               lower = rep(0, 10000), upper = rep(255, 10000),
               seed = 232, popSize = 50, maxiter = 1000)
reconstructed.face <- GA.face@solution
matrix(reconstructed.face, nrow=100, byrow=F) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
matrix(reconstructed.face, nrow=100, byrow=F) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)

reconstructed.face <- reconstructed.face/255
plot_face(reconstructed.face)

GA.face <- GA(fitness.func = function(x) -MSE.face(x),
              pop.size = 50, max.iter = 1000,
              lower = rep(0, 10000), upper = rep(255, 10000),
              crossover = blend.crossover, 
              mutation = unif.mutation,
              prob.mutation = 0.7,
              percent.elites = 70,
              seed = 232, keep = F, verbose = T)
reconstructed.face <- GA.face$solution
matrix(reconstructed.face, nrow=100, byrow=F) %>% 
  apply(2, rev) %>% t %>%
  image(col=rev(gray.colors(255)), axes=FALSE, asp=1)
matrix(reconstructed.face, nrow=100, byrow=F) %>% 
  apply(2, rev) %>% t %>%
  image(col=gray(1:0), axes=FALSE, asp=1)