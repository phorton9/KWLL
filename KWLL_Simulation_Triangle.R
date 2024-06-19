library(ggplot2)
library(MASS)

##Set parameters
alpha <- 0.05
beta <- 0.2
rho0 <- 0.8
r <- 0.5
nmax <- 200
n <- r*nmax

##Calculate \hat{\rho}_1 estimated from fixed design
exp1 <- exp(2*(qnorm(1-alpha)+qnorm(1-beta)+sqrt(n-3)*1/2*log((1+rho0)/(1-rho0)))/sqrt(n-3))
rho1 <- (exp1-1)/(exp1+1)


theta1 <- atanh(rho1)
theta0 <- atanh(rho0)
delta <- (theta1 - theta0)/2
thetaS <- (theta1 + theta0)/2
rhoTrue <- (exp(2*thetaS) -1)/(exp(2*thetaS) + 1)



##Covariance matrices
covTrue <- matrix(rhoTrue, 2,2)
diag(covTrue) <- 1
cov1 <- matrix(rho1, 2,2)
diag(cov1) <- 1
cov0 <- matrix(rho0, 2,2)
diag(cov0) <- 1

##Boundary thresholds
A <- .956
B <- .509



##Number of samples
n <- seq(1,nmax)

##Initialize boundaries
boundary1 <- numeric(nmax)
boundary0 <- numeric(nmax)

##Single iteration for visual
data <- mvrnorm(nmax, mu = c(0,0), Sigma = covTrue)
mle <- cor(data[1:4,1], data[1:4,2])
aR <- numeric(nmax)
for (i in 4:nmax){
  mle <- cor(data[1:i,1], data[1:i,2])
  aR[i] <- atanh(mle)*(i-3)
  boundary1[i] <- 2/(delta)*A + (thetaS - delta/2)*(i-3)
  boundary0[i] <- -2/(delta)*B+(thetaS+delta/2)*(i-3)
}

##Find where boundaries cross if there is an intersection
boundLimit <- which(boundary1 < boundary0)[1]


##Calculate minimum crossing point
crosses1 <- which(aR[1:nmax] > boundary1[1:nmax])[1]
crosses0 <- which(aR[1:nmax] < boundary0[1:nmax])[1]
crosses0[is.na(crosses0)] <- 1000
crosses1[is.na(crosses1)] <- 1000
min(crosses1,crosses0)

plotData <- as.data.frame(cbind(n[4:boundLimit],aR[4:boundLimit],boundary1[4:boundLimit],boundary0[4:boundLimit]))
                                                                                                                                                                                                     # panel.background = element_blank(), axis.line = element_line(colour = "black"))+ylim(-3,5)
# Subset for the first 60 observations
plotData_aR <- plotData[1:(min(crosses1,crosses0)-3), ]
len1 <- plotData_aR$V1[length(plotData_aR$V1)]

# Plot using ggplot2
ggplot() +
  geom_line(data = plotData, aes(x = n[4:boundLimit], y = boundary1[4:boundLimit], color = "Reject H0"), size = 1.5) +
  geom_line(data = plotData, aes(x = n[4:boundLimit], y = boundary0[4:boundLimit], color = "Reject H1"), size = 1.5)+
  geom_line(data = plotData_aR, aes(x = n[4:len1], y = aR[4:len1]), linetype = "dashed", size = 0.8) +
  scale_color_manual(values = c("Reject H0" = "blue", "Reject H1" = "red", "aR" = "blue")) +
  xlab('Number of Observations') +
  ylab('Test Statistic') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.position = c(0.2,0.8))

##calculation of Tables 1-2 with non-linear boundary (test statistic - artanh(r))
iters <- 10000

boundary1<- 2/(delta*(n-3))*A + (thetaS - delta/2)
boundary0 <- -2/(delta*(n-3))*B+(thetaS+delta/2)

##Initialize vectors for simulation results
crossed_null <- numeric(iters)
reject_null <- numeric(iters)
accept_null <- numeric(iters)

crossed_alt <- numeric(iters)
reject_alt <- numeric(iters)
accept_alt <- numeric(iters)

crossed_S <- numeric(iters)
reject_S <- numeric(iters)
accept_S <- numeric(iters)

for (x in 1:iters){
  ##Generate data for all three correlations (\rho_0, \rho_1, and \rho^*)
  data_null <- mvrnorm(nmax, mu = c(0,0), Sigma = cov0)
  data_alt <- mvrnorm(nmax, mu = c(0,0), Sigma = cov1)
  data_S <- mvrnorm(nmax, mu = c(0,0), Sigma = covTrue)

  ##Calculate sample correlation and atanh(r)
  aR_null <- numeric(nmax)
  aR_alt <- numeric(nmax)
  aR_S <- numeric(nmax)
  for (i in 4:nmax){
    ##H0
    mle_null <- cor(data_null[1:i,1], data_null[1:i,2])
    aR_null[i] <- atanh(mle_null)
    
    ##H1
    mle_alt <- cor(data_alt[1:i,1], data_alt[1:i,2])
    aR_alt[i] <- atanh(mle_alt)    

    ##ThetaStar
    mle_S <- cor(data_S[1:i,1], data_S[1:i,2])
    aR_S[i] <- atanh(mle_S)    
  }
  
  ##Minimum crossing point and rejection for H0
  crosses1_null <- which(aR_null[4:nmax] >= boundary1[4:nmax])[1]
  crosses0_null <- which(aR_null[4:nmax] <= boundary0[4:nmax])[1]
  crosses0_null[is.na(crosses0_null)] <- 1000
  crosses1_null[is.na(crosses1_null)] <- 1000
  crosses_null <- min(crosses1_null,crosses0_null)
  crossed_null[x] <- crosses_null+3
  reject_null[x] <- crosses_null == crosses1_null
  accept_null[x] <- crosses_null == crosses0_null
  
  ##Minimum crossing point and rejection for H1
  crosses1_alt <- which(aR_alt[4:nmax] >= boundary1[4:nmax])[1]
  crosses0_alt <- which(aR_alt[4:nmax] <= boundary0[4:nmax])[1]
  crosses0_alt[is.na(crosses0_alt)] <- 1000
  crosses1_alt[is.na(crosses1_alt)] <- 1000
  crosses_alt <- min(crosses1_alt,crosses0_alt)
  crossed_alt[x] <- crosses_alt+3
  reject_alt[x] <- crosses_alt == crosses1_alt
  accept_alt[x] <- crosses_alt == crosses0_alt
  
  ##Minimum crossing point and rejection for \rho^*
  crosses1_S <- which(aR_S[4:nmax] >= boundary1[4:nmax])[1]
  crosses0_S <- which(aR_S[4:nmax] <= boundary0[4:nmax])[1]
  crosses0_S[is.na(crosses0_S)] <- 1000
  crosses1_S[is.na(crosses1_S)] <- 1000
  crosses_S <- min(crosses1_S,crosses0_S)
  crossed_S[x] <- crosses_S+3
  reject_S[x] <- crosses_S == crosses1_S
  accept_S[x] <- crosses_S == crosses0_S
}
##Print results for H0
mean(crossed_null) ##Expected sample size
mean(reject_null) ##Type I error
mean(accept_null)
sqrt(var(crossed_null)/iters) ##Standard error of sample size
max(crossed_null) ##Max sample size

##Print results for H1
mean(crossed_alt) ##Expected sample size
mean(reject_alt)
mean(accept_alt) ##Type II error
sqrt(var(crossed_alt)/iters) ##Standard error of sample size
max(crossed_alt) ##Max sample size

##Print results for \rho^*
mean(crossed_S) ##Expected sample size
mean(reject_S)
mean(accept_S)
sqrt(var(crossed_S)/iters) ##Standard error of sample size
max(crossed_S) ##Max sample size


