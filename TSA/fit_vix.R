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

# univariate model
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
vit.fit.dse.2factor = dse::SS(F = f, Q = q,
                      H = h, R = r)
vit.m.fit.dse.2factor <- estMaxLik(vit.fit.dse.2factor, TSdata(output = cbind(log.vix.m,log.vix.m)))
vit.m.fit.dse.2factor$model

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
