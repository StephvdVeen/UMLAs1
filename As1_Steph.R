if(!require("fastDummies")) install.packages("fastDummies")
if(!require("pracma")) install.packages("pracma")
if(!require("Rculr")) install.packages("RCurl")

load("~/Documents/BDS/Year1/UnsupervisedMachineLearning/Assignment1/FIFA2017_NL.RData")
#summary(fifa)
#create club and position dummies
dummyfeatures <- c('club', 'Position')
XDummies <- data.matrix(dummy_cols(fifa[dummyfeatures],
              select_columns = dummyfeatures, remove_selected_columns = TRUE))
XNumeric <- data.matrix(fifa[,3:35])
#don't include one dummy variable for each category to prevent multicollinearity
mX <- cbind(constant=1, XDummies[,2:21], scale(XNumeric))
#omit observations with NA observation
mXOmitted <- na.omit(mX)


#ordinary PCA self-coded
mv <- svd(mXOmitted)$v #eigenvectors
mLambdas <- as.matrix((svd(mXOmitted)$d)^2) #eigenvalues
#determine the variance explained to decide on number of principal components
mTotVar <- sum(mLambdas)
mVarExplained <- matrix(0L, ncol =1, nrow = nrow(mLambdas))
dPartialVar <- 0.0
for (j in 1: nrow(mLambdas)){
  dPartialVar <- dPartialVar + mLambdas[j,1]
  mVarExplained[j,1] <- dPartialVar / mTotVar
}
plot(mVarExplained)
abline(h = 0.90, col = "red")
abline(h = 0.95, col="blue")

#ordinary PCA package
ordPCA <- prcomp(mXOmitted)
summary(ordPCA)
plot(ordPCA, type = "l")
biplot(ordPCA)


#sparce PCA self-coded(Varimax method)
#updateEigenV updates the eigenvalues for the rank-1 SVD
updateEigenV <- function(a, u, c){
  lambdas <- seq(0, 4, by = 0.1)
  prev_u <- u
  for (i in 1: length(lambdas)){
    i_lambda = lambdas[i]
    dS <- sign(a) * max((abs(a) - i_lambda), 0)
    u <- Norm(dS, p = 2)
    u <- dS / u
    if (Norm(u, p =1) <= c){
      prev_u <- u
    }
    else{
      break
    }
  }
  return(as.matrix(prev_u, ncol =1))
}

#this function performs the penalized rank-1 SVD
#input: data matrix, budget constraints, max number of iterations and convergence
# criterion
# output: singular vectors u and v, and sigma = u'Xv
ranks1SVD <- function(X, c1, c2, imax, eps = 1e-6){
  vV <- as.matrix(svd(X)$v[,1])
  vU <- svd(X)$u[,1]
  i <- 0
  mS <- 1/(nrow(X) - 1) * t(X) %*% X
  #loss_new <- t(vV) %*% mS %*% vV 
  while (i < imax ){
    i <- i + 1
    a <- X %*% vV
    vU <- updateEigenV(a, vU, c1)
    a <- t(X) %*% as.matrix(vU)
    vV <- updateEigenV(a, vV, c2)
    vV <- as.matrix(vV, ncol = 1)
  }
  dSigma <- t(vU) %*% X %*% vV
  return(list("u", vU, "v", vV, "sigma", dSigma))
}

#multifactor matrix decomposition
multifactor <- function(X, K, c1, c2, imax=1000, eps = 1e-6){
  mR <- X
  mVTot <- list()
  mUTot <- list()
  mSigmaTot <- list()
  for (k in 1: K){
    v1SVD <- ranks1SVD(mR, c1, c2, imax)
    vV_k <- as.matrix(v1SVD[[2]], ncol = 1)
    mVTot[[k]] <- vV_k
    vU_k <- as.matrix(v1SVD[[4]], ncol = 1)
    mUTot[[k]] <- vU_k
    dSigma_k <- v1SVD[[6]][1,1]
    mSigmaTot[[k]] <- dSigma_k
    mR <- mR - dSigma_k * vV_k %*% t(vU_k)
  }
  return(list("u_list", mUTot, "v_list", mVTot, "sigma_list", mSigmaTot)) 
}

c1 <- sqrt(nrow(mXOmitted))
c2 <- sqrt(ncol(mXOmitted))
K = 10
sPCACoded <- multifactor(mXOmitted, K, c1, c2, 1000)

#sparce PCA package (Varimax method)
vari <- varimax(ordPCA$rotation[,1:6])$rotmat
SPCA <- ordPCA$rotation[,1:6] %*% vari
plot(SPCA)