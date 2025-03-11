# Variable description:------
#  data: Name of data frame containing continuous, binary and group variables
#  var.end.c: Specify the variable name of the continuous endpoint in double quotation marks
#  var.end.b: Specify the variable name of the binary endpoint in double quotation marks
#  var.group: Specify the variable name of the group in double quotation marks
#  control.group: Specify the category of the group variable that corresponds to the control group in double quotation marks.
#  margin.c: Value of non-inferiority margin for a continuous endpoint
#  margin.b: Value of non-inferiority margin for a binary endpoint
#  alpha: Specify the one-sided significance level

# Results:------
# $IUT Results: The results of the IUT. If it is significant, “Significant” is returned, and if it is not significant, “Not significant” is returned.
# $N_treat: Number of subjects in the treatment group
# $N_control: Number of subjects in the control group
# $statistics_superiority_continuous_ua: Value of the test statistic uA
# $statistics_superiority_continuous_ub: Value of the test statistic uB
# $statistics_noninf_continuous: Value of the test statistic for non-inferiority of the continuous variable
# $statistics_noninf_binary: Value of the test statistic for non-inferiority of the binary variable
#----------------------------


# Library---------------------------
library(MASS)
library(bindata)
library(nleqslv)
library(polycor)
library(mvtnorm)

# Function: mixed.onecontinuous.onebinary------------
mixed.onecontinuous.onebinary <- function(data=df, var.end.c, var.end.b, var.group, control.group, margin.c, margin.b, alpha=0.05){
  
  eval(parse(text=paste("data$",var.group," <- as.character(data$",var.group,")
                        data$d_cont <- data$",var.end.c,"
                        data$d_bin  <- data$",var.end.b, sep="")))
  
  df_t  <- subset(data, Group != control.group)
  df_c  <- subset(data, Group == control.group)
  
  n1 <- nrow(df_t)
  n2 <- nrow(df_c)
  
  f <- function(c) {
    eq <- 0
    for(j in 0:2) {
      eq <- eq + choose(2, j) * (1-pchisq(c, j)) / 4
    }
    return(eq-alpha)
  }
  c <- nleqslv(c(5), f)$x
  
  y11 <- df_t$d_cont
  y12 <- df_t$d_bin
  y21 <- df_c$d_cont
  y22 <- df_c$d_bin
  y1  <- cbind(y11, y12)
  y2  <- cbind(y21, y22)
  
  #Estimate the parameter for Sigma_x----
  y11hat <- mean(y11,na.rm=T)
  y12hat <- mean(y12,na.rm=T)
  y21hat <- mean(y21,na.rm=T)
  y22hat <- mean(y22,na.rm=T)
  
  #Sigma_x----------
  g12 <- qnorm(1-y12hat)
  g22 <- qnorm(1-y22hat)
  
  h1 <- (1/sqrt(2*pi))*exp(-1*(g12^2)/2)
  h2 <- (1/sqrt(2*pi))*exp(-1*(g22^2)/2)
  
  s11 <- sqrt((t(y1[1:n1,1]-y11hat) %*% (y1[1:n1,1]-y11hat))/(n1-1))
  s21 <- sqrt((t(y2[1:n2,1]-y21hat) %*% (y2[1:n2,1]-y21hat))/(n2-1))
  
  polys1 <- polyserial(y1[1:n1,1],y1[1:n1,2])
  rho1   <- (polys1*h1)/sqrt(y12hat*(1-y12hat))
  polys2 <- polyserial(y2[1:n2,1],y2[1:n2,2])
  rho2   <- (polys2*h2)/sqrt(y22hat*(1-y22hat))
  
  sigma11_hat_1 <- (1/n1)*s11^2
  sigma22_hat_1 <- (1/n1)*(y12hat*(1-y12hat))
  sigma12_hat_1 <- rho1*sqrt(sigma11_hat_1)*sqrt(sigma22_hat_1)
  sigma21_hat_1 <- sigma12_hat_1
  
  sigma11_hat_2 <- (1/n2)*s21^2
  sigma22_hat_2 <- (1/n2)*(y22hat*(1-y22hat))
  sigma12_hat_2 <- rho2*sqrt(sigma11_hat_2)*sqrt(sigma22_hat_2)
  sigma21_hat_2 <- sigma12_hat_2
  
  Sigma_hat <- matrix(1, nrow = 2, ncol = 2)
  Sigma_hat[1,1] <- sigma11_hat_1 + sigma11_hat_2
  Sigma_hat[1,2] <- sigma12_hat_1 + sigma12_hat_2
  Sigma_hat[2,1] <- sigma21_hat_1 + sigma21_hat_2
  Sigma_hat[2,2] <- sigma22_hat_1 + sigma22_hat_2
  
  A_hat <- matrix(1, nrow = 2, ncol = 2)
  B_hat <- matrix(1, nrow = 2, ncol = 2)
  
  decom <- function(m) {
    temp <- eigen(m)
    A <- temp$vectors %*% sqrt(solve(diag(temp$values))) %*% t(temp$vectors)
    return(A)
  }
  
  A_hat[c(1,2),c(1,2)] <- decom(Sigma_hat[c(1,2),c(1,2)])
  B_hat[c(1,2),c(1,2)] <- matrix(c(A_hat[1,1],abs(A_hat[1,2]),abs(A_hat[2,1]),A_hat[2,2]), ncol=2)
  
  Y1j_hat <- c(y11hat,y12hat)
  Y2j_hat <- c(y21hat,y22hat)
  X_hat   <- Y1j_hat - Y2j_hat
  
  tr_ua <- t(A_hat[c(1,2),c(1,2)]%*%X_hat)
  tr_ub <- t((det(A_hat[c(1,2),c(1,2)])/det(B_hat[c(1,2),c(1,2)]))*B_hat[c(1,2),c(1,2)]%*%X_hat)
  
  ua_1 <- ifelse(tr_ua[1,1]>0,(tr_ua[1,1])^2,0)
  ua_2 <- ifelse(tr_ua[1,2]>0,(tr_ua[1,2])^2,0)
  stat_ua <- ua_1 + ua_2
  
  ub_1 <- ifelse(tr_ub[1,1]>0,(tr_ub[1,1])^2,0)
  ub_2 <- ifelse(tr_ub[1,2]>0,(tr_ub[1,2])^2,0)
  stat_ub <- ub_1 + ub_2
  
  stat_u <- ifelse(stat_ua>stat_ub,stat_ub,stat_ua)
  
  #Superiority-------
  r_sup <- ifelse(stat_u>c,1,0)
  
  #Non-inferiority-------
  #continuous
  x1 <- y11hat - y21hat
  z1 <- (x1+margin.c)/sqrt(Sigma_hat[1,1])
  r_inf1 <- ifelse(z1>qnorm(1 - alpha),1,0)
  
  #binary
  MLE_noninf <- function(NN1,NN2,p_hat1,p_hat2,margin.b){
    aa <- 1 + (NN2/NN1)
    bb <- -1*(1+(NN2/NN1)+p_hat1+(NN2/NN1)*p_hat2-margin.b*(NN2/NN1+2))
    cc <- margin.b^2 - margin.b*(2*p_hat1 + (NN2/NN1) + 1) + p_hat1 + (NN2/NN1)*p_hat2
    dd <- p_hat1*margin.b*(1-margin.b)
    
    vv <- (bb^3)/(3*aa)^3 - (bb*cc)/(6*aa^2) + dd/(2*aa)
    uu <- sign(vv)*sqrt((bb^2)/(3*aa)^2 - cc/(3*aa))
    ww <- (pi+acos(vv/uu^3))/3
    
    p1_tilde <- 2*uu*cos(ww)-bb/(3*aa)
    p2_tilde <- p1_tilde + margin.b
    c(p1_tilde,p2_tilde)
  }
  
  Y12_tilde <- MLE_noninf(n1,n2,y12hat,y22hat,margin.b)
  
  sigma2_noninf <- (Y12_tilde[1]*(1-Y12_tilde[1]))/n1 + (Y12_tilde[2]*(1-Y12_tilde[2]))/n2
  
  x2 <- y12hat - y22hat
  z2 <- (x2 + margin.b)/sqrt(sigma2_noninf)
  r_inf2 <- ifelse(z2>qnorm(1-alpha),1,0)
  
  result_sup_inf <- ifelse(!is.na(r_sup)&!is.na(r_inf1)&!is.na(r_inf2)&sum(r_sup,r_inf1,r_inf2,na.rm=T)==3,"Significant","Not Significant")
  
  statistics_noninf_continuous <- z1
  statistics_noninf_binary <- z2
  
  statistics_superiority_continuous_ua <- stat_ua
  statistics_superiority_continuous_ub <- stat_ub
  
  Result_test <- result_sup_inf
  
  Result_all <- list(Result_test,
                     n1,n2,
                     statistics_superiority_continuous_ua,
                     statistics_superiority_continuous_ub,
                     statistics_noninf_continuous,
                     statistics_noninf_binary)
  names(Result_all) <- c("IUT Results",
                         "N_treat",
                         "N_control",
                         "statistics_superiority_continuous_ua",
                         "statistics_superiority_continuous_ub",
                         "statistics_noninf_continuous",
                         "statistics_noninf_binary")
  
  return(Result_all)
}

# Example:----

# Generating Data----
set.seed(125)
n <- 100
ID <- 1:n  
Group <- sample(c("A", "B"), n, replace = TRUE)
ContinuousOutcome1 <- rnorm(n, mean = 50, sd = 10)
ContinuousOutcome2 <- rnorm(n, mean = 30, sd = 5)
BinaryOutcome1 <- rbinom(n, size = 1, prob = 0.5)
BinaryOutcome2 <- rbinom(n, size = 1, prob = 0.3)
df <- data.frame(ID, Group, ContinuousOutcome1, ContinuousOutcome2, BinaryOutcome1, BinaryOutcome2)

# Test----
mixed.onecontinuous.onebinary(data=df,
                              var.end.c="ContinuousOutcome1",
                              var.end.b="BinaryOutcome1",
                              var.group="Group",
                              control.group="B",
                              margin.c=1,
                              margin.b=0.05,
                              alpha=0.05)

