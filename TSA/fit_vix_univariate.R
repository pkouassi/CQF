setwd("H:\\Users\\Ran Zhao\\Research\\Fitting VIX")

# load in data, daily/monthly
vix_indeX_monthly = read.csv("vix_d.csv")
vix_index_daily = read.csv("vix_m.csv")

log.vix.m = log(vix_indeX_monthly[,2])
log.vix.d = log(vix_index_daily[,2])


############################################################################
########################   Empirical TSA   #################################
# stationary check
plot(log.vix.m, type='l', col='blue',main='Monthly log(VIX)')
acf(log.vix.m)   # not stationary

# take the first difference of the log(VIX) process
d.log.vix.m = diff(log.vix.m)
acf(d.log.vix.m)  # stationary
pacf(d.log.vix.m) # preliminary order detection

# fit autoregression model
short.term = 1:6
long.term = 10:15

aic.matrix = matrix(0, nrow = length(short.term), ncol = length(long.term), dimnames = list(c(short.term),c(long.term)))
  
for (i in 1:length(short.term)){
  for (j in 1:length(long.term)){
    fit.para = rep(0,max(long.term)+1)
    fit.para[short.term[i]] = NA   # fit the short term VIX
    fit.para[long.term[j]] = NA   # fit the long term VIX
    #fit.para[max(long.term)+1] = NA  # fit the intercept
    fit.arima = arima(d.log.vix.m, order=c(max(long.term),0,0), fixed=fit.para)
    aic.matrix[i,j] = fit.arima$aic
  }
}
 
# the final fitting model
short.term.index = which(aic.matrix == min(aic.matrix), arr.ind=TRUE)[1]
long.term.index = which(aic.matrix == min(aic.matrix), arr.ind=TRUE)[2]
fit.para.final = rep(0,max(long.term)+1)

fit.para.final[short.term[short.term.index]] = NA   # fit the short term VIX
fit.para.final[long.term[long.term.index]] = NA   # fit the long term VIX
#fit.para.final[max(long.term)+1] = NA  # fit the intercept
fit.arima.final = arima(d.log.vix.m, order=c(max(long.term),0,0), fixed=fit.para.final)
############################################################################

############################################################################
#####################   State Variable Method   ############################
library(dse)

# univariate model 1
output.model1 = optim(c(0.82,0.13,1,0,0.51),vit.fit.dse.likelihood.model1)
k = (1-output.model1$par[1])*12
mu = (output.model1$par[5]/k)*12
sigma_y = output.model1$par[2] * sqrt(12)
k
mu
sigma_y

# univariate model 2
output.model2 = optim(c(0.82,0.13,1,0,0.51),vit.fit.dse.likelihood.model2)
k = (1-output.model2$par[1])*12
mu = (output.model2$par[5]/k)*12
sigma_y = output.model2$par[2] * sqrt(12)
sigma_e = output.model2$par[4]
k
mu
sigma_y
sigma_e

# multivariate model 3
in_arg = c(.2,.6,.05,0.5,0.5)
output.model3 =  optim(in_arg, vit.fit.dse.likelihood.model3)


################
vit.fit.dse = dse::SS(F = matrix(-0.06896509, 1, 1), 
                      Q = matrix(0.1783604, 1, 1),
                      H = matrix(1, 1, 1), 
                      R = matrix(0, 1, 1),               
                      G = matrix(0.4445145, 1, 1),       # k_y * mu_y * dt
                      constants = list(R = matrix(TRUE, 1, 1)),
                      z0 = matrix(0, 1, 1), P0 = matrix(10^5, 1, 1))
vit.m.fit.dse <- estMaxLik(vit.fit.dse, TSdata(input = rep(1,length(log.vix.m)), output = log.vix.m))
vit.m.fit.dse$model
k = (1-vit.m.fit.dse$model$F)*12
mu = (vit.m.fit.dse$model$G/k)*12
sigma_y = vit.m.fit.dse$model$Q * sqrt(12)
k
mu
sigma_y

# two factor model
f <- array(c(.2,0,.8,.6),c(2,2))
h <- array(c(1,0,0,1),c(2,2))
q <- array(c(.1,0,0,.05),c(2,2))
r <- array(c(1,0,0,1),c(2,2))
g <- array(c(0.5,0,0,0.5),c(2,2))
vit.fit.dse.2factor = dse::SS(F = f, Q = q,
                      H = h, R = r, G = g)
vit.m.fit.dse.2factor <- estMaxLik(vit.fit.dse.2factor, TSdata(input = cbind(rep(1,length(log.vix.m)),rep(1,length(log.vix.m))), output = cbind(log.vix.m,log.vix.m)))
vit.m.fit.dse.2factor$model

#vit.fit.dse.likelihood <- function(f=0.5,q=0.5,h=1,r=0,g=0.5){
vit.fit.dse.likelihood.model1 <- function(in_arg){
    f = in_arg[1]
    q = in_arg[2]
    h = in_arg[3]
    r = in_arg[4]
    g = in_arg[5]
    vit.fit.dse.local = dse::SS(F = matrix(f, 1, 1), 
                        Q = matrix(q, 1, 1),
                        H = matrix(h, 1, 1), 
                        R = matrix(r, 1, 1),               
                        G = matrix(g, 1, 1),       # k_y * mu_y * dt
                        constants = list(R = matrix(TRUE, 1, 1),H = matrix(TRUE, 1, 1)),
                        z0 = matrix(0, 1, 1), P0 = matrix(10^5, 1, 1))
  vit.m.fit.dse <- estMaxLik(vit.fit.dse.local, TSdata(input = rep(1,length(log.vix.m)), output = log.vix.m))
  #vit.m.fit.dse$model
  results = vit.m.fit.dse$estimates$results$value
  return(results)
}


vit.fit.dse.likelihood.model2 <- function(in_arg){
  f = in_arg[1]
  q = in_arg[2]
  h = in_arg[3]
  r = in_arg[4]
  g = in_arg[5]
  vit.fit.dse.local = dse::SS(F = matrix(f, 1, 1), 
                              Q = matrix(q, 1, 1),
                              H = matrix(h, 1, 1), 
                              R = matrix(r, 1, 1),
                              G = matrix(g, 1, 1),       # k_y * mu_y * dt
                              #constants = list(R = matrix(TRUE, 1, 1),H = matrix(TRUE, 1, 1)),
                              z0 = matrix(0, 1, 1), P0 = matrix(10^5, 1, 1))
  vit.m.fit.dse <- estMaxLik(vit.fit.dse.local, TSdata(input = rep(1,length(log.vix.m)), output = log.vix.m))
  #vit.m.fit.dse$model
  results = vit.m.fit.dse$estimates$results$value
  return(results)
}


vit.fit.dse.likelihood.model3 <- function(in_arg){
  f <- array(c(1-in_arg[1],0,in_arg[1],in_arg[2]),c(2,2))
  h <- array(c(1,0,0,1),c(2,2))
  q <- array(c(.1,0,0,in_arg[3]),c(2,2))
  r <- array(c(1,0,0,1),c(2,2))
  g <- array(c(in_arg[4],0,0,in_arg[5]),c(2,2))
  vit.fit.dse.2factor.local = dse::SS(F = f,
                                Q = q,
                                H = h, 
                                R = r, 
                                G = g,
                                z0 = matrix(0, 1, 1), P0 = matrix(10^5, 1, 1))
  vit.m.fit.dse.2factor <- estMaxLik(vit.fit.dse.2factor.local, 
                                           TSdata(input = cbind(rep(1,length(log.vix.m)),rep(1,length(log.vix.m))), output = cbind(log.vix.m,log.vix.m)))
  results = vit.m.fit.dse.2factor$estimates$results$value
  return(results)
}

############################################################################






data("Nile", package = "datasets")
m1.dse <- dse::SS(F = matrix(1, 1, 1), Q = matrix(40, 1, 1),
                  H = matrix(1, 1, 1), R = matrix(130, 1, 1),
                  z0 = matrix(0, 1, 1), P0 = matrix(10^5, 1, 1))
m1b.dse <- dse::SS(F = matrix(1, 1, 1), Q = matrix(40, 1, 1),
                   H = matrix(1, 1, 1), R = matrix(130, 1, 1),
                   constants = list(Q = matrix(TRUE, 1, 1), P0 = matrix(TRUE, 1, 1)),
                   z0 = matrix(0, 1, 1), P0 = matrix(10^5, 1, 1))
m2.dse <- dse::SS(F = matrix(0.5, 2, 2), Q = matrix(0, 1, 1),
                  H = matrix(1, 2, 2), R = matrix(130, 1, 1),
                  z0 = matrix(0, 1, 1), P0 = matrix(1, 1, 1))
m2.dse.est <- estMaxLik(m2.dse, TSdata(output = Nile))
m2.dse.est <- estMaxLik(m2.dse, TSdata(output = log.vix.m))

f <- array(c(.5,.3,.2,.4),c(2,2))
h <- array(c(1,0,0,1),c(2,2))
k <- array(c(.5,.3,.2,.4),c(2,2))
ss <- SS(F=f,G=NULL,H=h,K=k)
is.SS(ss)
ss
