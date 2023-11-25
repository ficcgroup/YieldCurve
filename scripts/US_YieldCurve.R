library(readxl)
library(BVAR)
library(xts)
# load data. The data can be found in :
# https://fred.stlouisfed.org/release/tables?rid=354&eid=212994#snid=213004
ds <- read_excel("C:/Users/USER01/Downloads/ZEROYIELDSUS.xls")

# Create ts objects for zero coupon rates
DGS1Z = ts(ds$DGS1Z, start = c(1990,1),freq=12)
DGS2Z = ts(ds$DGS2Z, start = c(1990,1),freq=12)
DGS3Z = ts(ds$DGS3Z, start = c(1990,1),freq=12)
DGS4Z = ts(ds$DGS4Z, start = c(1990,1),freq=12)
DGS5Z = ts(ds$DGS5Z, start = c(1990,1),freq=12)
DGS6Z = ts(ds$DGS6Z, start = c(1990,1),freq=12)
DGS7Z = ts(ds$DGS7Z, start = c(1990,1),freq=12)
DGS8Z = ts(ds$DGS8Z, start = c(1990,1),freq=12)
DGS9Z = ts(ds$DGS9Z, start = c(1990,1),freq=12)
DGS10Z = ts(ds$DGS10Z, start = c(1990,1),freq=12)

yields = cbind(DGS1Z, DGS2Z, DGS3Z, DGS4Z, DGS5Z, DGS6Z, DGS7Z, DGS8Z, DGS9Z,
               DGS10Z)
yields = window(yields, end = c(2023,10))

# plot the yields:
ts.plot(yields, col = c(1:10),ylim = c(0,9.25), lty = c(1,2,4,5,1,2,4,5,1,2),
        ylab = "Zero-Coupon Yield (%)")
legend("top", 
      legend = c("1Y","2Y","3Y","4Y","5Y","6Y","7Y","8Y","9Y","10Y"), 
      col = c(1:10), cex =0.5, horiz = TRUE, lty = c(1,2,4,5,1,2,4,5,1,2), 
      lwd = c(2,2,2,2), bty = "n")

# reminder: Dynamic Nelson Siegel is:
# t : maturity in months
# y(t) = beta1_t + beta2_t*((1-exp(-lambda*t))/(lambda*t)) + 
# beta3_t*(((1-exp(-lambda*t))/(lambda*t)) - exp(-lambda*t))
# set lambda = 0.0609
lambda = 0.0609
lambda = 0.077
# we calculate the factor loadings l1, l2, l3 for t = 1, 2, ..., 10
# loading 1 is just 1 as beta1 is an intercept.
l1 = c(1,1,1,1,1,1,1,1,1,1)
l2 = numeric(10)
l3 = numeric(10)

for (t in c(1:10)*12) {
  l2[t/12] = ((1-exp(-lambda*t))/(lambda*t))
  l3[t/12] = (((1-exp(-lambda*t))/(lambda*t)) - exp(-lambda*t))
}

# visualize the factor loadings:
ts.plot(ts(l1), ts(l2), ts(l3), col = c(1,2,4), 
        type = "b", lwd = 1.5,ylab = "Loading", xlab = "Maturity (Years)")
legend("right",legend = c("Level", "Slope", "Curvature"),col = c(1,2,4),
       cex =1,inset = 0.05, lty = 1, lwd = 2, bty = "n")

# we proceed to step one of two-step DNS: calculate the factors through 
# cross-sectional OLS, by regressing the yields on the factor loadings
T = length(yields)/10

beta1 = numeric(T)
beta2 = numeric(T)
beta3 = numeric(T)

for (t in c(1:T)) {
  lmcoef = as.numeric(lm(as.numeric(yields[t,]) ~ l2 + l3)$coefficients)
  beta1[t] = lmcoef[1]
  beta2[t] = lmcoef[2]
  beta3[t] = lmcoef[3]
}

# we convert the betas into ts objects with their interpreted names:
level = ts(beta1, start = c(1990,1), freq=12)
slope = ts(beta2, start = c(1990,1), freq=12)
curvature = ts(beta3, start = c(1990,1), freq=12)
yieldfactors = cbind(level, slope, curvature)

# We believe these factors are part of a dependent, autoregressive structure. 
# Therefore, we model them as part of a VAR process. 
# We choose to estimate the model through the BVAR methodology originally
# proposed by Litterman (1980), but with the hierarchical procedure in Giannone,
# Lenza & Primiceri (2012).

# we cut the yield factor series to 2023-09:
yieldfactors1 = window(yieldfactors, end = c(2023,10))

# We create the dataframe object for estimation:
vars = as.data.frame((yieldfactors))
colnames(vars) = colnames(yieldfactors)
rownames(vars) = yearmon(time(yields))

# We can now set the Minnesota prior:
mn = bv_minnesota(
  lambda = bv_lambda(),
  alpha = bv_alpha(mode = 2),
  psi = bv_psi(),
  b = 0.97
)

priors = bv_priors(
  hyper = c("lambda","psi"),
  mn = mn
)

# We provide custom settings for the Metropolis-Hastings algorithm:
mh = bv_metropolis(
  scale_hess = 0.025,
  adjust_acc = TRUE,
  adjust_burn = 0.5,
  acc_change = 0.1
)

# We proceed to estimate a BVAR(13):
run1 = bvar(vars, 13, 
            n_draw = 200000, 
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

# We calculate and store the implied yields using the factor loadings:
ylds_median = matrix(nrow = 24, ncol = 10)
ylds_32 = matrix(nrow = 24, ncol = 10)
ylds_68 = matrix(nrow = 24, ncol = 10)

for (i in c(1:24)) {
  for (j in c(1:10)) {
    ymed = level_fcast[i,3] + slope_fcast[i,3]*l2[j] + curvature_fcast[i,3]*l3[j]
    y32 = level_fcast[i,2] + slope_fcast[i,2]*l2[j] + curvature_fcast[i,2]*l3[j]
    y68 = level_fcast[i,4] + slope_fcast[i,4]*l2[j] + curvature_fcast[i,4]*l3[j]
    ylds_median[i,j] = ymed
    ylds_32[i,j]= y32
    ylds_68[i,j]= y68
  }
}

# We plot some yield curve forecasts along with the 09-2023 yield curve:
ts.plot(ts(rep(0,10)), 
        ylab = "Zero-Coupon Yield (%)", xlab = "Maturity (Years)", 
        ylim = c(3.4,5.4))

axis(1, at = 1:10, labels = 1:10)

legend("bottomright", 
       legend = c("10-2023 (obs.)", "10-2024 (median pred.)",
                  "10-2024 (32% pred.)", 
                  "10-2024 (68% pred.)"), 
       col = c("#000000", "#EE0000", "#8B000055" ,"#CD000055"), cex =0.8, 
       inset = 0.01, lty = c(1,1,1,1), lwd = c(2,2,2,2), bty = "n")

polygon(c(c(1:10), rev(c(1:10))), c(as.numeric(ylds_32[12,]), 
                                    rev(ylds_68[12,])),
        col = rgb(0.8, 0.8, 0.8, alpha = 0.5), border = NA)
lines(c(1:10), as.numeric(yields[406,]),lty = 1, lwd = 1.5, type = "b", 
      col = "#000000")
lines(c(1:10), ylds_median[12,],lty = 1, lwd = 1.5, type = "b", 
      col = "#EE0000")
lines(c(1:10), ylds_68[12,],lty = 1, lwd = 1.5, type = "b", 
      col = "#CD000055")
lines(c(1:10), ylds_32[12,],lty = 1, lwd = 1.5, type = "b", 
      col = "#8B000055")