a = c(7, 8, 5, 7, 3, 4, 4, 1, 3, 2, 
      6, 2, 9, 7, 3, 4, 2, 1, 4, 7, 
      6, 4, 14, 6, 0, 1, 5, 0, 0, 2, 
      1, 8, 6, 0, 4, 10, 3, 4, 7, 0, 0)

b = c(8, 11, 29, 29, 9, 3, 13, 15, 11, 36,
      6, 5, 12, 14, 22, 7, 8, 30, 24, 36,
      34, 14, 54, 15, 6, 9, 12, 10, 22, 16,
      14, 16, 6, 20, 13, 30, 13, 30, 31, 34, 9)

c = c(11, 8, 4, 4, 0, 4, 13, 13, 7, 12,
      8, 7, 7, 5, 11, 6, 8, 4, 15, 16,
      13, 5, 13, 8, 6, 5, 5, 12, 8, 10,
      7, 15, 7, 5, 2, 12, 2, 5, 15, 34, 0)

d = c(2, 8, 35, 27, 12, 0, 11, 3, 15, 20,
      0, 2, 17, 20, 21, 4, 2, 23, 16, 27,
      8, 34, 61, 13, 0, 10, 10, 2, 16, 11,
      6, 12, 2, 18, 14, 8, 14, 14, 22, 0, 1)

theta = log((a/b)/(c/d))

sd = (1/(a+0.5)+1/(b+0.5)+1/(c+0.5)+1/(d+0.5))^0.5

data = as.data.frame(cbind(a, b, c, d, theta, sd), head = TRUE)

# Objective: whether the probability of occurrence is the same under treatment 
#            control in each individual table

# Tool: intrinsic moment prior 


# Define the functions for the analysis
###---calculate the Bayes factor for the IM prior
Konst <- function(h, a1, b1, a2, b2) {
  const <- rep(0, 2*h+1)
  for (j in 0:(2*h)) {
    const[j+1] <- lchoose(2*h, j) + lbeta(a1+j, b1) + lbeta(a2+2*h-j, b2) -
      lbeta(a1, b1) - lbeta(a2, b2)
    const[j+1] <- exp(const[j+1]) * (-1)^j
  }
  return(sum(const))
}

###---calculate the Bayes factor for the moment prior
bf_10_M <- function(y1, y2, n1, n2, h, a0, b0, a1, b1, a2, b2) {
  
  a0_post <- a0 + y1 + y2
  b0_post <- b0 + n1 + n2 - y1 - y2
  
  a1_post <- a1 + y1
  b1_post <- b1 + n1 - y1
  
  a2_post <- a2 + y2
  b2_post <- b2 + n2 - y2
  # BF10, the last equation
  res = - lbeta(a0_post, b0_post) - lbeta(a1, b1) - lbeta(a2, b2) + lbeta(a0, b0) + 
    lbeta(a1_post, b1_post) + lbeta(a2_post, b2_post)
  
  return(exp(res)*Konst(h, a1_post, b1_post, a2_post, b2_post)/Konst(h, a1, b1, a2, b2))
}

###---calculate Standard Bayes factor with the default prior (last equation on paage 407)
bf_10_standard <- function(y1, y2, n1, n2, h, a0, b0, a1, b1, a2, b2) {
  
  a0_post <- a0 + y1 + y2
  b0_post <- b0 + n1 + n2 - y1 - y2
  
  a1_post <- a1 + y1
  b1_post <- b1 + n1 - y1
  
  a2_post <- a2 + y2
  b2_post <- b2 + n2 - y2
  # BF10, the last equation
  res = - lbeta(a0_post, b0_post) - lbeta(a1, b1) - lbeta(a2, b2) + lbeta(a0, b0) + 
    lbeta(a1_post, b1_post) + lbeta(a2_post, b2_post)
  
  return(exp(res))
}


###---calculate the Bayes factor for the moment prior (matrix form)
BF_10_M <- function(y1, y2, n1, n2, h, t1, t2, a0, b0, a1, b1, a2, b2) {
  
  BF_10 <- matrix(0, nrow=t1+1, ncol=t2+1)
  for (x1 in 0:t1) {
    for (x2 in 0:t2) {
      a1_x <- a1 + x1
      b1_x <- b1 + t1 - x1
      a2_x <- a2 + x2
      b2_x <- b2 + t2 - x2
      
      BF_10[x1+1, x2+1] <- bf_10_M(y1, y2, n1, n2, h, a0, b0, a1_x, b1_x, a2_x, b2_x)
    }
  }
  return(BF_10)
}



m0 <- function(x1, x2, t1, t2, a0, b0) {
  
  a0_post <- a0 + x1 + x2
  b0_post <- b0 + t1 + t2 - x1 - x2
  
  L_res <- lchoose(t1, x1) + lchoose(t2, x2) + lbeta(a0_post, b0_post) -lbeta(a0, b0)
  res <- exp(L_res)
}


###---calculate the marginal likelihood under the null (matrix form)
# Eq 19 sum part
M0 <- function(t1, t2, a0, b0) {
  
  M <- matrix(0, nrow=t1+1, ncol=t2+1)
  for (x1 in 0:t1) {
    for (x2 in 0:t2) {
      M[x1+1, x2+1] <- m0(x1, x2, t1, t2, a0, b0)
    }
  }
  return(M)
}

BF_10_IM <- function(y1, y2, n1, n2, h=1, t1=4, t2=4, 
                     a0=1, b0=1, a1=1, b1=1, a2=1, b2=1) {
  
  return(sum(BF_10_M(y1, y2, n1, n2, h, t1, t2, a0, b0, a1, b1, a2, b2) * M0(t1, t2, a0, b0)))
}

##### Density for the IM prior

density_M <- function (theta1, theta2, h, a1, b1, a2, b2) {
  res <- exp(dbeta(theta1, a1, b1, log=TRUE) + dbeta(theta2, a2, b2, log=TRUE)) * 
    (theta1-theta2)^{2*h}/Konst(h, a1, b1, a2, b2)
  
  return(res)
}


density_IM <- function(theta1, theta2, h, t1, t2, a0, b0,
                       a1, b1, a2, b2) {
  Den <- matrix(0, nrow=t1+1, ncol=t2+1)
  for (x1 in 0:t1) {
    for (x2 in 0:t2) {
      a1_x <- a1 + x1
      b1_x <- b1 + t1 - x1
      a2_x <- a2 + x2
      b2_x <- b2 + t2 - x2
      
      Den[x1+1, x2+1] <- density_M(theta1, theta2, h, a1_x, b1_x, a2_x, b2_x)
    }
  }
  
  M <- M0(t1, t2, a0, b0)
  return( sum(Den * M) )
}



####################### Looking for the optimal t ########################
## Bayes Factor under different h and t
## Choice for h = 0, 1, 2

head(data)

uncont = data[,1]+data[,2]
cont = data[,3]+data[,4]

stor = c()
for (i in (1:length(cont))){
  temp = BF_10_IM(y1=data[,1][i], y2=data[,3][i], n1=uncont[i], n2=cont[i], h=1, t1=4, t2=4, a0=1, b0=1, a1=1, b1=1, a2=1, b2=1)
  stor = c(stor, temp)
}

stor

rankCri = abs(data[,1]/uncont - data[,3]/cont)
stor[order(rankCri)]

plot(stor[order(rankCri)][1:25])

# do an overall analysis
#### Bayes Factor under different h and t
## Choice for h = 0, 1, 2
## (a, b, c, d) = (170,736,352,556).
temp = BF_10_IM(y1=170, y2=352, n1=170+736, n2=352+556, h=1, t1=4, t2=4, a0=1, b0=1, a1=1, b1=1, a2=1, b2=1)


### find the best t for h = 0, 1, 2
# h = 0

thetaCont = data[,1]/(data[,1]+data[,2])
thetaTrt = data[,3]/(data[,3]+data[,4])

#################### BF ####################
# h=0, t=0
stor1 = numeric(41)
for (i in 1:41) {
    temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 0, 
                t1 = 0, t2 = 0, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                b1 = 0.25, a2 = 0.25, b2 = 0.25)
    stor1[i] = temp
}

# h=0, t=8
stor2 = numeric(41)
for (i in 1:41) {
  temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 0, 
                     t1 = 4, t2 = 4, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                     b1 = 0.25, a2 = 0.25, b2 = 0.25)
  stor2[i] = temp
}


# h=0, t=8
stor3 = numeric(41)
for (i in 1:41) {
  temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 0, 
                     t1 = 7, t2 = 7, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                     b1 = 0.25, a2 = 0.25, b2 = 0.25)
  stor3[i] = temp
}
stor3

cbind(stor1, stor2, stor3)

### h=1

# h=1, t=0
stor11 = numeric(41)
for (i in 1:41) {
  temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 1, 
                     t1 = 0, t2 = 0, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                     b1 = 0.25, a2 = 0.25, b2 = 0.25)
  stor11[i] = temp
}

# h=1, t=8
stor12 = numeric(41)
for (i in 1:41) {
  temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 1, 
                     t1 = 4, t2 = 4, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                     b1 = 0.25, a2 = 0.25, b2 = 0.25)
  stor12[i] = temp
}


# h=0, t=8
stor3 = numeric(41)
for (i in 1:41) {
  temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 0, 
                     t1 = 7, t2 = 7, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                     b1 = 0.25, a2 = 0.25, b2 = 0.25)
  stor3[i] = temp
}
stor3

cbind(stor1, stor2, stor3)

# h=0, t=0
stor1 = numeric(41)
for (i in 1:41) {
  temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 0, 
                     t1 = 0, t2 = 0, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                     b1 = 0.25, a2 = 0.25, b2 = 0.25)
  stor1[i] = temp
}

# h=0, t=8
stor2 = numeric(41)
for (i in 1:41) {
  temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 0, 
                     t1 = 4, t2 = 4, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                     b1 = 0.25, a2 = 0.25, b2 = 0.25)
  stor2[i] = temp
}


# h=0, t=8
stor3 = numeric(41)
for (i in 1:41) {
  temp = density_IM (theta1 = thetaCont[i], theta2 = thetaTrt[i], h = 0, 
                     t1 = 7, t2 = 7, a0 = 0.5, b0 = 0.5, a1 = 0.25, 
                     b1 = 0.25, a2 = 0.25, b2 = 0.25)
  stor3[i] = temp
}
stor3

cbind(stor1, stor2, stor3)


