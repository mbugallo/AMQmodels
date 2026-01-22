
#######################################################
###     READING AND PLOTTING SIMULATION RESULTS     ###
#######################################################

# The code reads Monte Carlo MSE results from multiple replicates
# All results and visualizations are fully reproducible using the provided CSV files
# Results are stored in a list MSEresults for further processing
MSEresults <- list()
index <- 0
for(r in c(10, 25, 50, 100, 200)){
  index <- index + 1
  MSEresults[[index]] <- read.csv(paste0("MSEresultsB/MSEresults_R", r, ".csv"))
}


# Check the Monte Carlo results
# RMSE values are extracted using sapply for each replicate
RMSE.mse.replicatesB <- sapply(MSEresults, function(x){ x$RMSE.mse2 })
par(mar=c(6, 4.5, 4, 2), xpd=F)

# Boxplots are created to visualize RMSE across different numbers of replicates
boxplot(RMSE.mse.replicatesB, names = c("B10", "B25", "B50", "B100", "B200"), ylab = "RMSE",
        cex.axis = 1.5, cex.lab = 1.7, pch=19)

# BIAS values are also extracted and plotted similarly.
BIAS.mse.replicatesB <- sapply(MSEresults, function(x){ x$BIAS.mse2 })
par(mar=c(6, 4.5, 4, 2), xpd=F)
boxplot(BIAS.mse.replicatesB, names = c("B10", "B25", "B50", "B100", "B200"), ylab = "BIAS",
        cex.axis = 1.5, cex.lab = 1.7, pch=19, ylim=c(-1.5, 1.5))
# A horizontal red dashed line at zero highlights unbiasedness in BIAS plots
abline(h=0, col='red', lty=2, lwd=2)

# Comparisons between AMQ analytic, AMQ bootstrap, and Hájek estimators are visualized.
RMSE <- MSEresults[[index]][ 1:3] 
BIAS <- MSEresults[[index]][ 4:6] 

par(mar=c(6, 4.5, 4, 2), xpd=F)
boxplot(RMSE, names = c("AMQ analytic", "AMQ boot", "Hajek"),
        ylab = "RMSE", xlab = "", pch=19, cex.axis = 1.7, cex.lab = 2)

par(mar=c(6, 4.5, 4, 2), xpd=F)
boxplot(BIAS, names = c("AMQ analytic", "AMQ boot", "Hajek"), ylim=c(-5,5),
        ylab = "BIAS", xlab = "", pch=19, cex.axis = 1.7, cex.lab = 2)
abline(h=0, col='red', lty=2, lwd=2)


# Sample sizes are handled separately and incorporated in plots
# RMSE and BIAS are computed for multiple sampling fractions (p=0.05, 0.10, 0.15, 0.25).
results <- read.csv('MSEresultsB/MSEresults_R200.csv')
RMSE.05  <- results[, 1] 
RMSEb.05 <- results[, 2] 
RMSEH.05 <- results[, 3] 

results <- read.csv('MSEresultsB/MSEresults_p010_R200.csv')
RMSE.10  <- results[, 1] 
RMSEb.10 <- results[, 2] 
RMSEH.10 <- results[, 3]

results <- read.csv('MSEresultsB/MSEresults_p015_R200.csv')
RMSE.15  <- results[, 1] 
RMSEb.15 <- results[, 2] 
RMSEH.15 <- results[, 3]

results <- read.csv('MSEresultsB/MSEresults_p025_R200.csv')
RMSE.25  <- results[, 1] 
RMSEb.25 <- results[, 2] 
RMSEH.25 <- results[, 3]

# Par settings ensure consistent margins and label sizes for all plots.
par(mar=c(6, 4.5, 4, 2), xpd=F)
boxplot(cbind(RMSE.05, RMSE.10, RMSE.15, RMSE.25), 
        names = c("p=0.05", "p=0.10", "p=0.15", "p=0.25"), 
        ylim=c(0, 4.5), ylab = "RMSE analytic", xlab = "", 
        pch=19, cex.axis = 1.7, cex.lab = 2)

par(mar=c(6, 4.5, 4, 2), xpd=F)
boxplot(cbind(RMSEb.05, RMSEb.10, RMSEb.15, RMSEb.25),
        names = c("p=0.05", "p=0.10", "p=0.15", "p=0.25"), 
        ylim=c(0, 4.5), ylab = "RMSE bootstrap", xlab = "", 
        pch=19, cex.axis = 1.7, cex.lab = 2)

par(mar=c(6, 4.5, 4, 2), xpd=F)
boxplot(cbind(RMSEH.05, RMSEH.10, RMSEH.15, RMSEH.25), 
        names = c("p=0.05", "p=0.10", "p=0.15", "p=0.25"), 
        ylim=c(0, 4.5), ylab = "RMSE Hajek", xlab = "", 
        pch=19, cex.axis = 1.7, cex.lab = 2)


