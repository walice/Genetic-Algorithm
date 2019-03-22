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
varnames <- c("age", "age^2", "educ", "educ^2", "black", "hispanic", "married", "no degree", 
              "re1974", "re1974^2", "re1975", "re1975^2", "unemp74", "unemp75", "re1974*re1975", 
              "age*educ", "educ*re1974", "educ*1975")
baltest <- baltest.collect(no.bal,
                           var.names = varnames,
                           after = F)
bal.out.before <- round(baltest[, c("mean.Tr", "mean.Co", "T pval", "KS pval")], 2)
colnames(bal.out.before) <- c("Mean Treated", "Mean Control", "t p-value", "KS p-value")
kable(bal.out.before, format = "latex")


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
colnames(bal.out.ps) <- c("Mean Treated", "Mean Control", "t p-value", "KS p-value")
kable(bal.out.ps, format = "latex")


# .. Genetic Matching using genoud ####
genoud.out <- GenMatch(Tr = D, X = X, 
                       BalanceMatrix = BalanceMatrix, 
                       pop.size = 3000,
                       loss = 2) # To remove lexical optimization
match.genoud <- Match(Y = Y, Tr = D, X = X, 
                      Weight.matrix = genoud.out)
genoud.bal <- MatchBalance(D ~ age + I(age^2) + educ + I(educ^2) 
                           + black + hisp + married + nodegr 
                           + re74 + I(re74^2) + re75 + I(re75^2) 
                           + u74 + u75 + I(re74 * re75)
                           + I(age * nodegr) + I(educ * re74)
                           + I(educ * re75),
                           match.out = match.genoud, 
                           nboots = 1000, data = lalonde)
baltest <- baltest.collect(genoud.bal,
                           var.names = varnames,
                           after = T)
bal.out.genoud <- round(baltest[, c("mean.Tr", "mean.Co", "T pval", "KS pval")], 2)
colnames(bal.out.genoud) <- c("Mean Treated", "Mean Control", "t p-value", "KS p-value")
kable(bal.out.genoud, format = "latex")



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
                  pop.size = 1000, max.iter = 50,
                  lower = rep(1, nbalvars), upper = rep(1000, nbalvars),
                  keep.track = F, seed = 232)

GA.out.QQ <- GA(fitness.func = function(x) -maxQQdif(x), 
                pop.size = 1000, max.iter = 10,
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
kable(bal.out.GA.pval, format = "latex")

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
kable(bal.out.GA.QQ, format = "latex")

# All together
bal.out <- cbind(bal.out.before, bal.out.ps, bal.out.GA.pval)
bal.out <- bal.out[,c(3:4,7:8,11:12)]
#colnames(bal.out) <- c("Before matching", "Matching on PS", "Matching using GA")
colnames(bal.out) <- rep(c("t", "KS"), 3)
kable(bal.out, format = "latex")


# .. Plot results ####
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
