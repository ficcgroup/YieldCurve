library(readxl)
library(BVAR)
library(xts)
library(tidyquant)

# We aggregate the daily series into monthly averages:
ds = read_excel("C:/Users/USER01/Downloads/ZEROYIELDSCAD.xlsx", na = "na")
ds1 = apply.monthly(as.xts(na.omit(ds[,1:121]), date_col = Date), mean)
dsf = (as.data.frame(ds1)) # since yield values are in decimals, not percent
dsf_adj = dsf*100
# We remove observations for maturities < 1 year and > 10 years:

yields = dsf_adj[,4:40]

#We plot the yields:
ts.plot(ts(yields[c(1,5,9,13,17,21,25,29,33,37)],
           start = c(1991,1),freq=12),
        col = c(1:10),lty = c(1,2,4,5,1,2,4,5,1,2),ylim = c(0,10.6),
        ylab = "Zero-Coupon Yield (%)")
legend("top", 
       legend = c("1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y"), 
       col = c(1:10), cex =0.5, horiz = TRUE, lty = c(1,2,4,5,1,2,4,5,1,2), 
       lwd = c(2,2,2,2), bty = "n")

# We calculate the factor loadings similarly to US yields, but lambda changes.
lambda = 0.106
l1 = rep(1,37)
l2 = c()
l3 = c()

for (t in c(4:40)*3) {
  l2 = c(l2, ((1-exp(-lambda*t))/(lambda*t)))
  l3 = c(l3, (((1-exp(-lambda*t))/(lambda*t)) - exp(-lambda*t)))
}

# We plot the loadings:
ts.plot(ts(l1, freq=4), ts(l2, freq=4), ts(l3, freq=4), col = c(1,2,4), 
        type = "b", lwd = 1.5,ylab = "Loading", xlab = "Maturity (Years)")
axis(1, at = 1:10, labels = 1:10)
legend("right",legend = c("Level", "Slope", "Curvature"), col = c(1,2,4), 
       cex =1,  inset = 0.05, lty = 1, lwd = 2, bty = "n")

# We estimate the Canadian yield factors through OLS:
T1 = length(yields[,1])
beta1 = numeric(T1)
beta2 = numeric(T1)
beta3 = numeric(T1)

for (t in c(1:T1)) {
  lmcoef = as.numeric(lm(as.numeric(yields[t,]) ~ l2 + l3)$coefficients)
  beta1[t] = lmcoef[1]
  beta2[t] = lmcoef[2]
  beta3[t] = lmcoef[3]
}

# we convert the Canadian betas into ts objects with their interpreted names:
level = ts(beta1, start = c(1991,1), freq=12)
slope = ts(beta2, start = c(1991,1), freq=12)
curvature = ts(beta3, start = c(1991,1), freq=12)
yieldfactors = cbind(level, slope, curvature)

# We believe these factors are part of a dependent, autoregressive structure. 
# Therefore, we model them as a VAR process. 
# We create the dataframe object for estimation:
vars = as.data.frame(yieldfactors)

colnames(vars) = colnames(yieldfactors)
rownames(vars) = yearmon(time(yieldfactors))

# We can now set the Minnesota prior:
mn = bv_minnesota(
  lambda = bv_lambda(),
  alpha = bv_alpha(mode = 2),
  psi = bv_psi(),
  b = 0.97
)
# We add the SOC and SUR priors:
priors = bv_priors(
  hyper = c("lambda", "psi"),
  mn = mn
)

# We provide custom settings for the Metropolis-Hastings algorithm:
mh = bv_metropolis(
  scale_hess = 0.05,
  adjust_acc = TRUE,
  adjust_burn = 0.5,
  acc_change = 0.2
)

# We proceed to estimate a BVAR(13):
run1 = bvar(vars, 13, 
            n_draw = 250000, 
            n_burn = 3000, 
            n_thin = 5, 
            priors = priors,
            mh = mh
)

# We visualize the trace and density plots of the Marginal Likelihood and
# lambda: 
plot(run1)

# We forecast 24 months ahead with 90%, 95% credible intervals, and extract
# the forecasts:
pred = predict(run1, 
               conf_bands = c(0.05,0.32), 
               n_thin = 2,
               horizon = 24
)

level_fcast = matrix(nrow = 24, ncol = 5)
slope_fcast = matrix(nrow = 24, ncol = 5)
curvature_fcast = matrix(nrow = 24, ncol = 5)

for (i in c(1:5)) {
  for (j in c(1:24)) {
    level_fcast[j,i] = pred$quants[i,j,1]
    slope_fcast[j,i] = pred$quants[i,j,2]
    curvature_fcast[j,i] = pred$quants[i,j,3]
  }
}

colnames(level_fcast) = c("5%", "32%", "50%", "68%", "95%")
colnames(slope_fcast) = c("5%", "32%", "50%", "68%", "95%")
colnames(curvature_fcast) = c("5%", "32%", "50%", "68%", "95%")

ylds_median = matrix(nrow = 24, ncol = 37)
ylds_32 = matrix(nrow = 24, ncol = 37)
ylds_68 = matrix(nrow = 24, ncol = 37)

for (i in c(1:24)) {
  for (j in c(1:37)) {
    ymed = level_fcast[i,3] + slope_fcast[i,3]*l2[j] + curvature_fcast[i,3]*l3[j]
    y32 = level_fcast[i,2] + slope_fcast[i,2]*l2[j] + curvature_fcast[i,2]*l3[j]
    y68 = level_fcast[i,4] + slope_fcast[i,4]*l2[j] + curvature_fcast[i,4]*l3[j]
    ylds_median[i,j] = ymed
    ylds_32[i,j]= y32
    ylds_68[i,j]= y68
  }
}

lst = c(1,5,9,13,17,21,25,29,33,37)

ts.plot(ts(rep(0,10)), 
        col = c("#8B000033" ,"#CD000033", "#EE0000" ,"#000000"),
        ylab = "Zero-Coupon Yield (%)", xlab = "Maturity (Years)", 
        ylim = c(3.5,5.75), lty = c(1,1,1,1), lwd = c(1.5,1.5,1.5,1.5),
        type = "b")

axis(1, at = 1:10, labels = 1:10)

legend("topright", 
       legend = c("10-2023 (obs.)", "10-2024 (median pred.)",
                  "10-2024 (32% pred.)", 
                  "10-2024 (68% pred.)"), 
       col = c("#000000","#EE0000","#8B000055" ,"#CD000055"), cex =0.8, 
       inset = 0.00, lty = c(1,1,1,1), lwd = c(2,2,2,2), bty = "n")

polygon(c(c(1:10), rev(c(1:10))), c(as.numeric(ylds_32[12,][lst]), 
                                    rev(ylds_68[12,][lst])),
        col = rgb(0.8, 0.8, 0.8, alpha = 0.5), border = NA)
lines(c(1:10), as.numeric(yields[394,])[lst],lty = 1, lwd = 1.5, type = "b",
      col = "#000000")
lines(c(1:10), ylds_median[12,][lst],lty = 1, lwd = 1.5, type = "b", 
      col = "#EE0000")
lines(c(1:10), ylds_68[12,][lst],lty = 1, lwd = 1.5, type = "b", 
      col = "#CD000055")
lines(c(1:10), ylds_32[12,][lst],lty = 1, lwd = 1.5, type = "b", 
      col = "#8B000055")
