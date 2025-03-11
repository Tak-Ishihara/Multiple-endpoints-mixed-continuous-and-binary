# Variable description:----
#  data: Name of data frame containing continuous, binary and group variables
#  var.end.c: Specify the variable name of the continuous endpoint in double quotation marks
#  var.end.b1: Specify the variable name of the binary endpoint1 in double quotation marks
#  var.end.b2: Specify the variable name of the binary endpoint2 in double quotation marks
#  var.group: Specify the variable name of the group in double quotation marks
#  control.group: Specify the category of the group variable that corresponds to the control group in double quotation marks.
#  margin.c: Value of non-inferiority margin for a continuous endpoint
#  margin.b1: Value of non-inferiority margin for a binary endpoint1
#  margin.b2: Value of non-inferiority margin for a binary endpoint2
#  alpha: Specify the one-sided significance level

# Results:----
# $IUT Results: The results of the IUT. If it is significant, “Significant” is returned, and if it is not significant, “Not significant” is returned.
# $N_treat: Number of subjects in the treatment group
# $N_control: Number of subjects in the control group
# $statistics_superiority_continuous_ua: Value of the test statistic uA
# $statistics_superiority_continuous_ub: Value of the test statistic uB
# $statistics_noninf_continuous: Value of the test statistic for non-inferiority of the continuous variable
# $statistics_noninf_binary1: Value of the test statistic for non-inferiority of the binary variable1
# $statistics_noninf_binary2: Value of the test statistic for non-inferiority of the binary variable2
#-------------

#------------#
# Library----
#------------#
library(MASS)
library(bindata)
library(nleqslv)
library(polycor)
library(mvtnorm)

#--------------------------------------------#
# Function: mixed.onecontinuous.twobinary----
#--------------------------------------------#
mixed.onecontinuous.twobinary <- function(data=df,
                                          var.end.c="ContinuousOutcome1",
                                          var.end.b1="BinaryOutcome1",
                                          var.end.b2="BinaryOutcome2",
                                          var.group="Group",
                                          control.group="B",
                                          margin.c=1, 
                                          margin.b1=0.1, 
                                          margin.b2=0.1, 
                                          alpha=0.05){
  eval(parse(text=paste("data$",var.group," <- as.character(data$",var.group,")
                          data$d_cont <- data$",var.end.c,"
                          data$d_bin1 <- data$",var.end.b1,"
                          data$d_bin2 <- data$",var.end.b2, sep="")))
  
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
  y12 <- df_t$d_bin1
  y13 <- df_t$d_bin2
  y21 <- df_c$d_cont
  y22 <- df_c$d_bin1
  y23 <- df_c$d_bin2
  y1 <- cbind(y11, y12, y13)
  y2 <- cbind(y21, y22, y23)
  
  y11hat <- mean(y11,na.rm=T)
  y12hat <- mean(y12,na.rm=T)
  y13hat <- mean(y13,na.rm=T)
  y21hat <- mean(y21,na.rm=T)
  y22hat <- mean(y22,na.rm=T)
  y23hat <- mean(y23,na.rm=T)
  
  s11_1 <- sqrt((t(y1[1:n1,1]-y11hat) %*% (y1[1:n1,1]-y11hat))/(n1-1))
  s11_2 <- sqrt((t(y2[1:n2,1]-y21hat) %*% (y2[1:n2,1]-y21hat))/(n2-1))
  
  #variance
  sigma11_hat_1 <- (1/n1)*s11_1^2
  sigma22_hat_1 <- (1/n1)*(y12hat*(1-y12hat))
  sigma33_hat_1 <- (1/n1)*(y13hat*(1-y13hat))
  sigma11_hat_2 <- (1/n2)*s11_2^2
  sigma22_hat_2 <- (1/n2)*(y22hat*(1-y22hat))
  sigma33_hat_2 <- (1/n2)*(y23hat*(1-y23hat))
  
  #correlation(continuous&binary)
  g12 <- qnorm(1-y12hat)
  g13 <- qnorm(1-y13hat)
  g22 <- qnorm(1-y22hat)
  g23 <- qnorm(1-y23hat)
  
  h12 <- (1/sqrt(2*pi))*exp(-1*(g12^2)/2)
  h13 <- (1/sqrt(2*pi))*exp(-1*(g13^2)/2)
  h22 <- (1/sqrt(2*pi))*exp(-1*(g22^2)/2)
  h23 <- (1/sqrt(2*pi))*exp(-1*(g23^2)/2)
  
  polys12_1 <- polyserial(y1[1:n1,1],y1[1:n1,2])
  rho12_1   <- (polys12_1*h12)/sqrt(y12hat*(1-y12hat))
  polys13_1 <- polyserial(y1[1:n1,1],y1[1:n1,3])
  rho13_1   <- (polys13_1*h13)/sqrt(y13hat*(1-y13hat))
  
  polys12_2 <- polyserial(y2[1:n2,1],y2[1:n2,2])
  rho12_2   <- (polys12_2*h22)/sqrt(y22hat*(1-y22hat))
  polys13_2 <- polyserial(y2[1:n2,1],y2[1:n2,3])
  rho13_2   <- (polys13_2*h23)/sqrt(y23hat*(1-y23hat))
  
  sigma12_hat_1 <- rho12_1*sqrt(sigma11_hat_1)*sqrt(sigma22_hat_1)
  sigma13_hat_1 <- rho13_1*sqrt(sigma11_hat_1)*sqrt(sigma33_hat_1)
  sigma21_hat_1 <- sigma12_hat_1
  sigma31_hat_1 <- sigma13_hat_1
  
  sigma12_hat_2 <- rho12_2*sqrt(sigma11_hat_2)*sqrt(sigma22_hat_2)
  sigma13_hat_2 <- rho13_2*sqrt(sigma11_hat_2)*sqrt(sigma33_hat_2)
  sigma21_hat_2 <- sigma12_hat_2
  sigma31_hat_2 <- sigma13_hat_2
  
  #correlation(binary&binary)
  tetra23_1 <- polychor(y1[1:n1,2],y1[1:n1,3])
  sigma_tetra1 <- matrix(1, nrow = 3, ncol = 3)
  sigma_tetra1[1,2] <- rho12_1
  sigma_tetra1[1,3] <- rho13_1
  sigma_tetra1[2,3] <- tetra23_1
  sigma_tetra1[2,1] <- sigma_tetra1[1,2]   
  sigma_tetra1[3,2] <- sigma_tetra1[2,3]
  sigma_tetra1[3,1] <- sigma_tetra1[1,3]
  p11_re_1 <- pmvnorm(lower=c(-Inf,g12,g13),upper=c(Inf,Inf,Inf),mean=c(0,0,0),corr=sigma_tetra1)[1]
  rho23_1  <- (p11_re_1 - y12hat*y13hat)/(sqrt(y12hat*(1-y12hat))*sqrt(y13hat*(1-y13hat)))
  sigma23_hat_1 <- rho23_1 * sqrt(sigma22_hat_1) * sqrt(sigma33_hat_1)
  sigma32_hat_1 <- sigma23_hat_1
  
  tetra23_2 <- polychor(y2[1:n2,2],y2[1:n2,3])
  sigma_tetra2 <- matrix(1, nrow = 3, ncol = 3)
  sigma_tetra2[1,2] <- rho12_2
  sigma_tetra2[1,3] <- rho13_2
  sigma_tetra2[2,3] <- tetra23_2
  sigma_tetra2[2,1] <- sigma_tetra2[1,2]   
  sigma_tetra2[3,2] <- sigma_tetra2[2,3]
  sigma_tetra2[3,1] <- sigma_tetra2[1,3]
  p11_re_2 <- pmvnorm(lower=c(-Inf,g22,g23),upper=c(Inf,Inf,Inf),mean=c(0,0,0),corr=sigma_tetra2)[1]
  rho23_2  <- (p11_re_2 - y22hat*y23hat)/(sqrt(y22hat*(1-y22hat))*sqrt(y23hat*(1-y23hat)))
  sigma23_hat_2 <- rho23_2 * sqrt(sigma22_hat_2) * sqrt(sigma33_hat_2)
  sigma32_hat_2 <- sigma23_hat_2
  
  Sigma_hat <- matrix(1, nrow = 3, ncol = 3)
  Sigma_hat[1,1] <- sigma11_hat_1 + sigma11_hat_2
  Sigma_hat[1,2] <- sigma12_hat_1 + sigma12_hat_2
  Sigma_hat[1,3] <- sigma13_hat_1 + sigma13_hat_2
  
  Sigma_hat[2,1] <- sigma21_hat_1 + sigma21_hat_2
  Sigma_hat[2,2] <- sigma22_hat_1 + sigma22_hat_2
  Sigma_hat[2,3] <- sigma23_hat_1 + sigma23_hat_2
  
  Sigma_hat[3,1] <- sigma31_hat_1 + sigma31_hat_2
  Sigma_hat[3,2] <- sigma32_hat_1 + sigma32_hat_2
  Sigma_hat[3,3] <- sigma33_hat_1 + sigma33_hat_2
  
  A_hat <- matrix(1, nrow = 3, ncol = 3)
  B_hat <- matrix(1, nrow = 3, ncol = 3)
  
  decom <- function(m) {
    temp <- eigen(m)
    A <- temp$vectors %*% sqrt(solve(diag(temp$values))) %*% t(temp$vectors)
    return(A)
  }
  
  A_hat <- decom(Sigma_hat)
  B_hat <- matrix(c(A_hat[1,1],abs(A_hat[1,2]),abs(A_hat[1,3]),abs(A_hat[2,1]),A_hat[2,2],abs(A_hat[2,3]),abs(A_hat[3,1]),abs(A_hat[3,2]),A_hat[3,3]), ncol=3)
  
  Y1j_hat <- c(y11hat,y12hat,y13hat)
  Y2j_hat <- c(y21hat,y22hat,y23hat)
  X_hat   <- Y1j_hat - Y2j_hat
  
  tr_ua <- t(A_hat%*%X_hat)
  tr_ub <- t((det(A_hat)/det(B_hat))*B_hat%*%X_hat)
  
  ua_1 <- ifelse(tr_ua[1,1]>0,(tr_ua[1,1])^2,0)
  ua_2 <- ifelse(tr_ua[1,2]>0,(tr_ua[1,2])^2,0)
  ua_3 <- ifelse(tr_ua[1,3]>0,(tr_ua[1,3])^2,0)
  stat_ua <- ua_1 + ua_2 + ua_3
  
  ub_1 <- ifelse(tr_ub[1,1]>0,(tr_ub[1,1])^2,0)
  ub_2 <- ifelse(tr_ub[1,2]>0,(tr_ub[1,2])^2,0)
  ub_3 <- ifelse(tr_ub[1,2]>0,(tr_ub[1,2])^2,0)
  stat_ub <- ub_1 + ub_2 + ua_3
  
  stat_u <- ifelse(stat_ua>stat_ub,stat_ub,stat_ua)
  
  #Superiority-------
  r_sup <- ifelse(stat_u>c,1,0)
  
  #Non-inferiority-------
  #continuous
  x1 <- y11hat - y21hat
  z1 <- (x1+margin.c)/sqrt(Sigma_hat[1,1])
  r_inf1 <- ifelse(z1>qnorm(1 - alpha/2),1,0)
  
  #binary
  MLE_noninf <- function(NN1,NN2,p_hat1,p_hat2,epsi){
    aa <- 1 + (NN2/NN1)
    bb <- -1*(1+(NN2/NN1)+p_hat1+(NN2/NN1)*p_hat2-epsi*(NN2/NN1+2))
    cc <- epsi^2 - epsi*(2*p_hat1 + (NN2/NN1) + 1) + p_hat1 + (NN2/NN1)*p_hat2
    dd <- p_hat1*epsi*(1-epsi)
    
    vv <- (bb^3)/(3*aa)^3 - (bb*cc)/(6*aa^2) + dd/(2*aa)
    uu <- sign(vv)*sqrt((bb^2)/(3*aa)^2 - cc/(3*aa))
    ww <- (pi+acos(vv/uu^3))/3
    
    p1_tilde <- 2*uu*cos(ww)-bb/(3*aa)
    p2_tilde <- p1_tilde + epsi
    c(p1_tilde,p2_tilde)
  }
  
  #outcome2
  Y12_tilde <- MLE_noninf(n1,n2,y12hat,y22hat,margin.b1)
  if(is.na(Y12_tilde[1]) | is.na(Y12_tilde[2])){
    error_count_noninf <- error_count_noninf + 1
    next
  }
  sigma2_noninf2 <- (Y12_tilde[1]*(1-Y12_tilde[1]))/n1 + (Y12_tilde[2]*(1-Y12_tilde[2]))/n2
  x2 <- y12hat - y22hat
  z2 <- (x2+margin.b1 )/sqrt(sigma2_noninf2)
  r_inf2 <- ifelse(z2>qnorm(1-alpha/2),1,0)
  
  #outcome3
  Y13_tilde <- MLE_noninf(n1,n2,y13hat,y23hat,margin.b2)
  if(is.na(Y13_tilde[1]) | is.na(Y13_tilde[2])){
    error_count_noninf <- error_count_noninf + 1
    next
  }
  sigma2_noninf3 <- (Y13_tilde[1]*(1-Y13_tilde[1]))/n1 + (Y13_tilde[2]*(1-Y13_tilde[2]))/n2
  x3 <- y13hat - y23hat
  z3 <- (x3+margin.b2 )/sqrt(sigma2_noninf3)
  r_inf3 <- ifelse(z3>qnorm(1-alpha/2),1,0)
  
  result_sup_inf <- ifelse(!is.na(r_sup)&!is.na(r_inf1)&!is.na(r_inf2)&!is.na(r_inf3)&sum(r_sup,r_inf1,r_inf2,r_inf3,na.rm=T)==4,
                           "Significant","Not Significant")
  
  statistics_noninf_continuous <- z1
  statistics_noninf_binary1 <- z2
  statistics_noninf_binary2 <- z3
  
  statistics_superiority_continuous_ua <- stat_ua
  statistics_superiority_continuous_ub <- stat_ub
  
  Result_test <- result_sup_inf
  
  Result_all <- list(Result_test,
                     n1,n2,
                     statistics_superiority_continuous_ua,
                     statistics_superiority_continuous_ub,
                     statistics_noninf_continuous,
                     statistics_noninf_binary1,
                     statistics_noninf_binary2)
  names(Result_all) <- c("IUT Results",
                         "N_treat",
                         "N_control",
                         "statistics_superiority_continuous_ua",
                         "statistics_superiority_continuous_ub",
                         "statistics_noninf_continuous",
                         "statistics_noninf_binary1",
                         "statistics_noninf_binary2")
  
  return(Result_all)
}

#------------#
# Example----
#------------#
# Generating Data----
set.seed(123)
n <- 100
ID <- 1:n  
Group <- sample(c("A", "B"), n, replace = TRUE)
ContinuousOutcome1 <- rnorm(n, mean = 50, sd = 10)
ContinuousOutcome2 <- rnorm(n, mean = 30, sd = 5)
BinaryOutcome1 <- rbinom(n, size = 1, prob = 0.5)
BinaryOutcome2 <- rbinom(n, size = 1, prob = 0.3)
df <- data.frame(ID, Group, ContinuousOutcome1, ContinuousOutcome2, BinaryOutcome1, BinaryOutcome2)

# Test----
mixed.onecontinuous.twobinary(data=df,
                              var.end.c="ContinuousOutcome1",
                              var.end.b1="BinaryOutcome1",
                              var.end.b2="BinaryOutcome2",
                              var.group="Group",
                              control.group="B",
                              margin.c=1, 
                              margin.b1=0.1, 
                              margin.b2=0.1, 
                              alpha=0.05)

#EOF
