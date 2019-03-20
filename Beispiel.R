library(dplyr)

N <- 7
T<- 3
M <- matrix(NA, nrow=N, ncol=T+1)
colnames(M) <- paste0("t", 0:T)
rownames(M) <- 1:N

M[,1] <- 1
r <- rnorm(N*T, sd=0.3)
rM <- matrix(r, nrow=N, ncol=T)
# Cumsums
for (i in 1:3){
  M[,i+1] <- M[,i]+rM[,i]
}

# plot
matplot(t(M), type="l")

# Auszahlungsmatrix
A <- matrix(NA, nrow=N, ncol=T+1)
K <- 1.10
lapply(M[,4],function(x) max(K-x,0))  %>% unlist -> A[,4]

# discounted chashflow 
r <- 0.06
# in the money im
im <- A[,4][A[,4]!=0]
Y <- exp(-r)*im
X <- M[,3][A[,4]!=0]
lm1 <- lm(Y~X+(X^2))
f1 <- lm1$fitted.values

# ausüben vs beibehalten
df <- data.frame(ausueben=f1, beibehalten=im)
A[,3] <- 0
ind <- apply(df,1, function(x) ifelse(x[2]>x[1], x[2], 0))
A[,3][A[,4]!=0] <- ind

#### ZUM ZEITPUNKT t = 2
lapply(M[,3],function(x) max(K-x,0))  %>% unlist -> a

# in the money im
im <- a[a!=0]
Y <- exp(-r)*im
X <- M[,2][a!=0]
lm1 <- lm(Y~X+(X^2))
f1 <- lm1$fitted.values

# ausüben vs beibehalten
df <- data.frame(ausueben=f1, beibehalten=im)
A[,2] <- 0
ind <- apply(df,1, function(x) ifelse(x[2]>x[1], x[2], 0))
A[,2][a!=0] <- ind

colnames(A) <- paste0("t", 0:T)
rownames(A) <- 1:N
