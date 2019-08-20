# instal all necessary packages
list_of_packages <- c("quantmod", "PerformanceAnalytics", "rugarch", "viridis")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new.packages)

library(quantmod)
library(PerformanceAnalytics)
library(rugarch)
library(viridis)

# saving rstudio plot area default format, clear session before running
par.defaults <- par(no.readonly = TRUE)
# run par(par.defaults) line before reproducing graphs

# load data
getSymbols(c("DGS1", "DGS2", "DGS3", "DGS5", "DGS7", "DGS10", "DGS20", "DGS30"), src = "FRED", from = "1990")
fed_yields <- merge(DGS1, DGS2, DGS3, DGS5, DGS7, DGS10, DGS20, DGS30)
# DGS30 omitted due to removal from 2002 - 2006
fed_yields <- na.omit(fed_yields[,1:7])

# potential 2nd data set for analysis, requires Quandl package
# library(Quandl)
# fed_yields <- Quandl("FED/SVENY", type = "xts")

# select beginning of model window
start <- paste0(1995, "/")

# plot instrument yields 
plot.zoo(fed_yields, plot.type = "single", col = viridis(ncol(fed_yields)))
  legend(x = "topleft", legend = colnames(fed_yields), col = viridis(ncol(fed_yields)), cex = 1, lwd = 3)
  title("Fed Yields")

# plot yield changes
fed_yield_changes <- diff.xts(fed_yields)
plot.zoo(fed_yield_changes, ylab = colnames(fed_yield_changes), main = "Fed Yield Changes", col = viridis(ncol(fed_yields)))

# autocorrelation analysis
for (i in 1:ncol(fed_yield_changes)){
  par(mfrow = c(1, 2))
  acf(fed_yield_changes[start,i], ylim=c(-0.2, 1),main =NA, na.action = na.pass)
  acf(abs(fed_yield_changes[start, i]), ylim=c(-0.2, 1),main =NA, na.action = na.pass)
  par(mfrow = c(1, 1))
  title(colnames(fed_yield_changes[,i]))
}

# garch volatility, with a t distribution
spec <- ugarchspec(distribution.model = "sstd")

# fit model to securities
model_fits <- apply(fed_yield_changes[start,], 2, FUN = ugarchfit, spec = spec)

# extract volatilities
fed_vols <- lapply(model_fits, FUN = sigma)
vols <- sapply(fed_vols, FUN = merge)
vols <- as.xts(vols, order.by = index(fed_yield_changes[start,]))

# plot garch volatilities
plot.zoo(vols, plot.type = "single", col = viridis(ncol(fed_yield_changes)), main = "Volatilities")
  legend(x = "topleft", legend = colnames(fed_yields), col = viridis(ncol(fed_yield_changes)), cex = 0.5, lwd = 3)

# plot garch volatilities
plot.zoo(vols, ylab = colnames(fed_yields), ylim = c(0, 0.2), col = viridis(ncol(fed_yield_changes)), main = "Garch Volatilities")

# begin analysis of model fit
# extract scaled residuals to be standardized
residuals <- sapply(model_fits, function(model_fits) scale(residuals(model_fits, standardize = TRUE)))
fed_yield_sd <- sapply(fed_yield_changes[start,], FUN = StdDev)
fed_yield_mean <- colMeans(fed_yield_changes[start,])

# standardize residuals
std_residuals <- residuals * fed_yield_sd + fed_yield_mean
std_residuals <- as.xts(std_residuals, order.by = index(fed_yield_changes[start,]))

# plot standardized residuals
plot.zoo(std_residuals, ylab = colnames(fed_yields), col = viridis(ncol(fed_yield_changes)))

#analyze standardized residual density
for (i in 1:ncol(fed_yield_changes)){
  dens <- density(fed_yield_changes[start, i])
  res_dens <- density(std_residuals[, i])
  
  plot(c(dens, res_dens))
  lines(res_dens, col = "blue")
  title(paste(colnames(fed_yield_changes[,i]), "density"))
  
  norm_dist <- dnorm(seq(-0.4, 0.4, by = .01), mean = fed_yield_mean[i], sd = fed_yield_sd[i])
  lines(seq(-0.4, 0.4, by = .01), 
        norm_dist, 
        col = "green"
  )
  legend("topleft", legend = c("Before GARCH", "After GARCH", "Normal distribution"), 
         col = c("black", "blue", "green"), lty=c(1,1))
}

# Q-Q Plot for quantile distribution of model
for (i in 1:ncol(fed_yield_changes)){
  qqnorm(fed_yield_changes[start, i], ylim = c(-0.5, 0.5), main = NA)
  qqline(fed_yield_changes[start, i], distribution = qnorm)
  
  par(new = TRUE)
  qqnorm(std_residuals[, i], col = "blue", ylim = c(-0.5, 0.5), main = NA)
  qqline(std_residuals[, i], distribution = qnorm, col = "blue")
  title(paste(colnames(fed_yield_changes[,i]), "Normal Q-Q Plot"))
  legend("topleft", c("Before GARCH", "After GARCH"), col = c("black", "blue"), pch=c(1,1))
}
