rm(list=ls())

#load libraries
library(MASS)
library(pscl)
library(coda)
library(fields)
library(MBA)
library(spdep)
library(TruncatedNormal)
library(Metrics)
library(bayestestR)
library(SMFilter)
library(scoringRules)
library(scoringutils)
library(parallel)
library(doParallel)

#load base data
setwd("~/Documents/FSU/Research/Bayesian Analysis with Truncated Distribution/revised2/new codes/")
#setwd("/gpfs/home/ss18eo/LargeSimulation_revised_1/")
#nrep <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#set_num <- 27
load("base_data_DeOliveira.RData")
df_param <- read.csv("df_param.csv")
#set.seed(nrep*set_num*50 + 11483)

#data simulation params
n_obs <- n-N
G <- 200
G0 <- 30
betas <- as.matrix(c(2,5))
a <- 2
b <- 5

#pre-code calculation
Dwinv <- solve(Dw)
A <- Dwinv%*%W
M <- (diag(n)-A)%*%Dwinv
M2 <- M^2
M3 <- M^3
M4 <- M^4

# Calculate X
eigen_obj <- eigen(W)
lambda_1 <- eigen_obj$values[1]
u_1 <- eigen_obj$vectors[,1]
X <- as.matrix(cbind(u_1,runif(n)))

#functions
#function to calculate R^2
func_r2 <- function(actVec,predVec){
  ybar <- mean(actVec)
  return(1 - sum((predVec-actVec)^2)/sum((actVec-ybar)^2))
}

#function to calculate (1/delta)*Dw + (Dw-W) for each rho
Sigma_O_inv <- function(del){
  Sigma_W <- (1/del)*Dw + Dw - W
  return( Sigma_W[c(1:n_obs),c(1:n_obs)] )
}

Delta <- seq(0.1,2,0.1)
#log of det(Sigma^-1) using eigen values for each rho
log_det_Sigma_O_inv <- sapply(1:length(Delta),func <- function(k){
  eigen_vals <- eigen(Sigma_O_inv(Delta[k]))$values
  return( sum(log(eigen_vals)) )
})

#TAR_C
TARC_func <- function(Yn,Xn){
  
  #function to calculate (X'*Sigma^-1*X)^-1 for each rho
  tXSX_inv <- array(0,dim=c(p,p,length(Delta)))
  for(i in 1:length(Delta))
    tXSX_inv[,,i] <- solve(t(Xn)%*%Sigma_O_inv(Delta[i])%*%Xn)
  
  #b of IG distribution for sigma_Y^2 for each rho
  IGbList <- sapply(1:length(Delta),func <- function(k){
    return( t(Yn)%*%Sigma_O_inv(Delta[k])%*%Yn - t(Yn)%*%Sigma_O_inv(Delta[k])%*%Xn%*%tXSX_inv[,,k]%*%t(Xn)%*%Sigma_O_inv(Delta[k])%*%Yn )
  })
  
  #log p(rho,y) with discrete uniform prior on rho, for each rho
  del_log_lik <- sapply(1:length(Delta),func <- function(k){
    return( 0.5*log_det_Sigma_O_inv[k] + 0.5*log(det(tXSX_inv[,,k])) - (a+(n_obs-p)/2)*log(b+0.5*IGbList[k]) )
  })
  
  #calculating probabilities for different rho
  r <- exp( del_log_lik-max(del_log_lik) )
  probVal <- r/sum(r)
  
  #calculating mean of beta for each rho
  beta_mean <- sapply(1:length(Delta),func <- function(k){
    return( tXSX_inv[,,k]%*%t(Xn)%*%Sigma_O_inv(Delta[k])%*%Yn )
  })
  
  #function to sample from beta for each combination of sigma_Y^2 and rho
  beta_dist <- function(x){
    sigma2_Y <- x[1]
    k <- which(Delta==x[2])
    return( mvrnorm(1, beta_mean[,k], sigma2_Y*tXSX_inv[,,k]) )
  }
  
  #sample of rho
  del_sample <- sample(Delta,G,replace = T,prob = probVal)
  
  #sample of sigma_Y^2
  sigma2_Y_sample <- sapply(del_sample,func <- function(del){
    k <- which(Delta==del)
    return( rigamma(1, a+(n_obs-p)/2, b+0.5*IGbList[k]) )
  })
  
  #join sigma_Y^2 and rho
  param_Mat <- as.matrix(rbind(sigma2_Y_sample,del_sample))
  
  #sample of beta
  beta_sample <- apply(param_Mat,2,beta_dist)
  
  return(rbind(beta_sample,param_Mat))
}

TARC_pred_func <- function(Q){
  beta_sample <- Q[1:2,]
  sigma2_Y_sample <- Q[3,]
  del_sample <- Q[4,]
  delunique <- unique(del_sample)
  SigmaMat <- array(0,dim=c(n,n,length(delunique)))
  for(i in 1:length(delunique))
    SigmaMat[,,i] <- delunique[i]*Dwinv - (delunique[i]^2)*M + (delunique[i]^3)*M2 - (delunique[i]^4)*M3 + (delunique[i]^5)*M4
  mu <- X%*%beta_sample
  yhat <- matrix(0,N,G)
  for(g in 1:G)
  {
    k <- which(delunique==del_sample[g])
    Z <- SigmaMat[c((n_obs+1):n),c(1:n_obs),k]%*%Sigma_O_inv(del_sample[g])
    mu_cross <- Z%*%(Y[1:n_obs]-mu[1:n_obs,g])
    sig_cross <- Z%*%SigmaMat[c(1:n_obs),c((n_obs+1):n),k]
    yhat[,g] <- sapply(1:N,func <- function(s0){
      return( rnorm(1,mu[(n_obs+s0),g]+mu_cross[s0],sqrt(sigma2_Y_sample[g]*(SigmaMat[(n_obs+s0),(n_obs+s0),k]-sig_cross[s0,s0]))) )
    })
    #print(paste("Iteration:",g))
  }
  return(yhat)
}

#support of rho
Delta_rho <- c(-0.99,seq(-0.95,0.95,0.05),0.99)[-21]
#function to calculate Dw-rho*W for each rho for CAR
Sigma_rho_inv <- function(rho){
  return( Dw[c(1:n_obs),c(1:n_obs)] - rho*W[c(1:n_obs),c(1:n_obs)] )
}
#log of det(Sigma^-1) using eigen values for each rho
log_det_Sigma_rho_inv <- sapply(1:length(Delta_rho),func <- function(k){
  eigen_vals <- eigen(Sigma_rho_inv(Delta_rho[k]))$values
  return( sum(log(eigen_vals)) )
})

CAR_func <- function(Yn,Xn){
  
  #function to calculate (X'*Sigma^-1*X)^-1 for each rho
  tXSX_inv <- array(0,dim=c(p,p,length(Delta_rho)))
  for(i in 1:length(Delta_rho))
    tXSX_inv[,,i] <- solve(t(Xn)%*%Sigma_rho_inv(Delta_rho[i])%*%Xn)
  
  #b of IG distribution for sigma_Y^2 for each rho
  IGbList <- sapply(1:length(Delta_rho),func <- function(k){
    return( t(Yn)%*%Sigma_rho_inv(Delta_rho[k])%*%Yn - t(Yn)%*%Sigma_rho_inv(Delta_rho[k])%*%Xn%*%tXSX_inv[,,k]%*%t(Xn)%*%Sigma_rho_inv(Delta_rho[k])%*%Yn )
  })
  
  #log p(rho,y) with discrete uniform prior on rho, for each rho
  rho_log_lik <- sapply(1:length(Delta_rho),func <- function(k){
    return( 0.5*log_det_Sigma_rho_inv[k] + 0.5*log(det(tXSX_inv[,,k])) - (a+(n_obs-p)/2)*log(b+0.5*IGbList[k]) )
  })
  
  #calculating probabilities for different rho
  r <- exp( rho_log_lik-max(rho_log_lik) )
  probVal <- r/sum(r)
  
  #calculating mean of beta for each rho
  beta_mean_CAR <- sapply(1:length(Delta_rho),func <- function(k){
    return( tXSX_inv[,,k]%*%t(Xn)%*%Sigma_rho_inv(Delta_rho[k])%*%Yn )
  })
  
  #function to sample from beta for each combination of sigma_Y^2 and rho
  beta_dist <- function(x){
    sigma2_Y <- x[1]
    k <- which(Delta_rho==x[2])
    return( mvrnorm(1, beta_mean_CAR[,k], sigma2_Y*tXSX_inv[,,k]) )
  }
  
  #sample of rho
  rho_sample <- sample(Delta_rho,G,replace = T,prob = probVal)
  
  #sample of sigma_Y^2
  sigma2_Y_sample_CAR <- sapply(rho_sample,func <- function(rho){
    k <- which(Delta_rho==rho)
    return( rigamma(1, a+(n_obs-p)/2, b+0.5*IGbList[k]) )
  })
  
  #join sigma_Y^2 and rho
  param_Mat <- as.matrix(rbind(sigma2_Y_sample_CAR,rho_sample))
  
  
  #sample of beta
  beta_sample_CAR <- apply(param_Mat,2,beta_dist)
  
  return(rbind(beta_sample_CAR,param_Mat))
}

CAR_pred_func <- function(Q){
  beta_sample_CAR <- Q[1:2,]
  sigma2_Y_sample_CAR <- Q[3,]
  rho_sample <- Q[4,]
  yhat_CAR <- matrix(0,n,G0)
  yhat_CAR[c(1:n_obs),] <- X[c(1:n_obs),]%*%beta_sample_CAR[,c(1:G0)]
  for(g in 1:G0)
  {
    SigmaInv <- Dw-rho_sample[g]*W
    Sigma_CAR <- chol2inv(chol(SigmaInv))
    sig <- Sigma_CAR[c((n_obs+1):n),c(1:n_obs)]%*%SigmaInv[c(1:n_obs),c(1:n_obs)]
    yhat_car_mu <- X[(n_obs+1):n,]%*%beta_sample_CAR[,g] + sig%*%(Y[c(1:n_obs)]-yhat_CAR[c(1:n_obs),g])
    yhat_car_sigma <- sigma2_Y_sample_CAR[g]*(Sigma_CAR[c((n_obs+1):n),c((n_obs+1):n)] - sig%*%Sigma_CAR[c(1:n_obs),c((n_obs+1):n)])
    yhat_CAR[(n_obs+1):n,g] <- sapply(1:N, func <- function(k){
      return( rnorm(1,yhat_car_mu[k],sqrt(yhat_car_sigma[k,k])) )
    })
    #print(paste("Iteration:",g))
  }
  return(yhat_CAR[(n_obs+1):n,])
}

df_score <- as.data.frame(matrix(0,2,21))
colnames(df_score) <- c("model","range_true","sigma2_true","beta1","beta1_lb","beta1_ub","beta2","beta2_lb","beta2_ub",
                       "sigma2","sigma2_lb","sigma2_ub","delta","delta_lb","delta_ub",
                       "r2","mae","rmse","crps","int","cvg")

#main code
#data from CAR -> already calculated for each setting
sigma2_val <- 0.2
range_val <- (1/lambda_1)+0.001
Sigma_y_CAR <- sigma2_val*chol2inv(chol(Dw-range_val*W))
#eigen_result <- eigen(Sigma_y_CAR)
#V_t <- sqrt(diag(eigen_result$values))%*%t(eigen_result$vectors) # M=VSV'=(VS^(1/2))%*%(VS^(1/2))'
#V_tX2 <- 3*V_t%*%X[,1] + rnorm(n,mean=0,sd=df_param$mult.sd[set_num])
#X[,2] <- solve(V_t)%*%V_tX2
Y <- mvrnorm(1,mu=X%*%betas,Sigma=Sigma_y_CAR)

#TAR_C
op_TARC <- TARC_func(Y[1:n_obs],X[c(1:n_obs),])
df_score[1,1] <- "TARC"
df_score[1,c(2,3,4,7,10,13)] <- c(range_val,sigma2_val,apply(op_TARC,1,mean))
df_score[1,5:6] <- ci(op_TARC[1,])[2:3]
df_score[1,8:9] <- ci(op_TARC[2,])[2:3]
df_score[1,11:12] <- ci(op_TARC[3,])[2:3]
df_score[1,14:15] <- ci(op_TARC[4,])[2:3]

yhat_TARC <- TARC_pred_func(op_TARC)
null_ind <- which(is.na(apply(yhat_TARC,2,sum))==T)
yhat_mean <- apply(yhat_TARC,1,mean)
yhat_sd <- apply(yhat_TARC,1,sd)
if(length(null_ind)!=0)
{
  yhat_mean <- apply(yhat_TARC[,-null_ind],1,mean)
  yhat_sd <- apply(yhat_TARC[,-null_ind],1,sd)
}
yhat_lb <- qnorm(0.05/2,mean=yhat_mean,sd=yhat_sd)
yhat_ub <- qnorm(1-0.05/2,mean=yhat_mean,sd=yhat_sd)
df_score[1,16:19] <- c(func_r2(Y[(n_obs+1):n],yhat_mean),
                            mae(Y[(n_obs+1):n],yhat_mean),
                            rmse(Y[(n_obs+1):n],yhat_mean),
                            mean(crps.numeric(Y[(n_obs+1):n],family="normal",mean=yhat_mean,sd=yhat_sd)))
df_score[1,20] <- mean(interval_score(Y[(n_obs+1):n],lower=yhat_lb,upper=yhat_ub,
                                           interval_range = rep(95,N),weigh=F))
df_score[1,21] <- mean((yhat_lb <= Y[(n_obs+1):n]) & (Y[(n_obs+1):n] <= yhat_ub))

#CAR
op_CAR <- CAR_func(Y[1:n_obs],X[c(1:n_obs),])
df_score[2,1] <- "CAR"
df_score[2,c(2,3,4,7,10,13)] <- c(range_val,sigma2_val,apply(op_CAR,1,mean))
df_score[2,5:6] <- ci(op_CAR[1,])[2:3]
df_score[2,8:9] <- ci(op_CAR[2,])[2:3]
df_score[2,11:12] <- ci(op_CAR[3,])[2:3]
df_score[2,14:15] <- ci(op_CAR[4,])[2:3]

yhat_CAR <- CAR_pred_func(op_CAR)
null_ind <- which(is.na(apply(yhat_CAR,2,sum))==T)
yhat_mean_CAR <- apply(yhat_CAR,1,mean)
yhat_sd_CAR <- apply(yhat_CAR,1,sd)
if(length(null_ind)!=0)
{
  yhat_mean_CAR <- apply(yhat_CAR[,-null_ind],1,mean)
  yhat_sd_CAR <- apply(yhat_CAR[,-null_ind],1,sd)
}
yhat_lb_CAR <- qnorm(0.05/2,mean=yhat_mean_CAR,sd=yhat_sd_CAR)
yhat_ub_CAR <- qnorm(1-0.05/2,mean=yhat_mean_CAR,sd=yhat_sd_CAR)
df_score[2,16:19] <- c(func_r2(Y[(n_obs+1):n],yhat_mean_CAR),
                           mae(Y[(n_obs+1):n],yhat_mean_CAR),
                           rmse(Y[(n_obs+1):n],yhat_mean_CAR),
                           mean(crps.numeric(Y[(n_obs+1):n],family="normal",mean=yhat_mean_CAR,sd=yhat_sd_CAR)))
df_score[2,20] <- mean(interval_score(Y[(n_obs+1):n],lower=yhat_lb_CAR,upper=yhat_ub_CAR,
                                          interval_range = rep(95,N),weigh=F))
df_score[2,21] <- mean((yhat_lb_CAR <= Y[(n_obs+1):n]) & (Y[(n_obs+1):n] <= yhat_ub_CAR))

fname <- paste0(getwd(),"/DeOliveiraPhiStudyOutput/df_score_3.csv")
write.csv(df_score,fname,row.names = F)

print(paste("range:",range_val,", sigma2:",sigma2_val,", nrep:",nrep,", setting:",set_num))



