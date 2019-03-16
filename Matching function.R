# D <- as.double(lalonde$treat)
# X <- cbind(lalonde$age, lalonde$educ, lalonde$black, lalonde$hisp, 
#            lalonde$married, lalonde$nodegr, lalonde$re74, lalonde$re75,
#            lalonde$u74, lalonde$u75)
nvars <- ncol(X)
nobs <- nrow(matched.data)

matching <- function(x){
  
  wmatrix <- diag(x, nrow = nvars)
  
  match.init <- Match(Tr = D, X = X,
                      Weight.matrix = wmatrix)
  #matched.data <- match.init$mdata$X
  index.treated <- match.init$index.treated
  index.control <- match.init$index.control
  #weights <- match.init$weights
  #sum(weights) == 185
  
  Tr <- X[index.treated,]
  Co <- X[index.control,]
  
  # paired.t.test <- function(Tr, Co){
  #   # nobs <- length(Tr)
  #   # dif <- Tr - Co
  #   # estimate <- sum(dif * weights)/sum(weights)
  #   # var <- sum(((dif - estimate)^2) * weights)/(sum(weights) * sum(weights))
  #   # 
  #   # if (estimate == 0 && var == 0){
  #   #   return(1)
  #   # }
  #   # 
  #   # statistic <- estimate/sqrt(var)
  #   # p.val <- (1 - pt(abs(statistic), df = nobs-1))*2
  #   p.val <- t.test(Tr, Co, paired = T)$p.value
  #   return(p.val)
  # }
  # 
  # ks.test <- function(Tr, Co){
  #   p.val <- ks.test(Tr, Co)$p.value
  #   return(p.val)
  # }
  # 
  storage <- NULL
  
  for (v in 1:nvars){
    storage[v] <- t.test(Tr, Co, paired = T)$p.value
  }
  for (v in (nvars+1):(nvars*2)){
    storage[v] <- ks.test(Tr, Co)$p.value
  }
  
  loss.func <- min(storage, na.rm = T)
  #S <- cov(X)
  #S.chol <- chol(S)
  #d <- t.test.out
  #distance <- sqrt(t(d) %*% t(S.chol) %*% wmatrix %*% S.chol %*% d)
  
  return(loss.func)
}
