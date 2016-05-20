setwd("H:\\Users\\Ran Zhao\\Research\\Fitting VIX")

# load in data, daily/monthly
vix_indeX_monthly = read.csv("vix_d.csv")
vix_index_daily = read.csv("vix_m.csv")

log.vix.m = log(vix_indeX_monthly[,2])
log.vix.d = log(vix_index_daily[,2])

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

