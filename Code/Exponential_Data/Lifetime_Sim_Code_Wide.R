library(ggplot2)
library(HDInterval)
library(rjags)
library(R2jags)
library(dplyr)

# Set the true parameters from the mortality curve
true.a  <-  0.001094
true.b <- 0.0488
true.c <- 0.609

# True quadratic
f <- function(x, a = 0.001094, b = 0.0488, c = 0.609){
  a * x^2 - b * x  +  c
}

# Generate 100 data sets and store them in a list
data.list <- list()
set.seed(52)
for(i in 1:100){
  T1 <- rexp(n = 100, f(10))
  T2 <- rexp(n = 100, f(17))
  T3 <- rexp(n = 100, f(21))
  T4 <- rexp(n = 100, f(25))
  T5 <- rexp(n = 100, f(32))
  data.list[[i]] <- data.frame(
    T = c(rep(10, 100), rep(17, 100), rep(21, 100), rep(25, 100), rep(32, 100)), 
    trait = c(T1, T2, T3, T4, T5)
  )
  }

b <- length(data.list)

# Truncation value
# Don't expect lifetime to exceed 100 days
epsilon <- 0.01

# True value vector with a and b in log space
# Spot for sigma and deviance
true_values <- c(log(true.a), log(true.b), true.c, NA, NA)

##########################################Lambda #############################################

# Create list of matrices to store parameter summary statistics
param_list1 <- vector("list", length = 4)
for (i in 1:4) {
# Initialize the inner list
  inner_list1 <- vector("list", length = b)
# Populate the inner list
  for (j in 1:b) {
  # Initialize a 5x4 matrix
    matrix_data1 <- matrix(0, nrow = 5, ncol = 4)
    inner_list1[[j]] <- matrix_data1
    }
  param_list1[[i]] <- inner_list1
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
    }
  # Model
  sink("lambda.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
  mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + 
  ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c)<= epsilon) * epsilon
  trait[i] ~ dexp(mu[i])
  }
      }", file = "lambda.txt")
  # Settings
  parameters <- c("la", "b.l", "c")
  inits <- function(){list(
    # a and b in log space
    la = log(0.001), 
    b.l = log(0.05), 
    c = 0.6
    )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data.raw <- data.list[[i]]
  data <- data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                model.file = "lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list1[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list1[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list1[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list1[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    param_list1[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list1[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list1[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list1[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list1[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list1[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list1[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list1[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for the deviance
    param_list1[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list1[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list1[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list1[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
  }
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new1 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 10, to = 35, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new1 as a function at each of those temperatures
xdf1 <- expand.grid(point = xseq, iteration = 1:nrow(df.new1))

# Apply the quadratic equation for each combination
xdf1$value <- apply(xdf1, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new1[i, 4]) * location^2 - exp(df.new1[i, 1]) * location + df.new1[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf1 <-  xdf1|> mutate(trunc.eval = ifelse(xdf1$value > epsilon, xdf1$value, epsilon), 
                   trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf1 <- xdf1 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                              avg.value = mean(trunc.eval), 
                                              med.value.inv = median(trunc.inv), 
                                              med.value = median(trunc.eval), 
                                              lower.hdi.inv = hdi(trunc.inv)[1], 
                                              upper.hdi.inv = hdi(trunc.inv)[2], 
                                              lower.hdi = hdi(trunc.eval)[1], 
                                              upper.hdi = hdi(trunc.eval)[2])
# Creating columns of the true curve values and the true inverted curve values
true_curve.inv <- function(x) 1/(true.a * x^2 - true.b * x + true.c)
true_curve <- function(x) (true.a * x^2 - true.b * x + true.c)
mean.xdf1$true_curve <- true_curve(xseq)
mean.xdf1$true_curve.inv <- true_curve.inv(xseq)
# Plot the mean mortality rate response
plot1_mean <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = avg.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
labs(title = "Mean Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal()
# Plot the median mortality rate response
plot1_med <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = med.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
labs(title = "Median Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal()
# Plot the mean lifetime response
plot1_mean.inv <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = avg.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
geom_point(aes(x = T, y = trait), data = data.raw) + 
labs(title = "Mean Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Lifetime") + 
theme_minimal()
# Plot the median lifetime response
plot1_med.inv <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = med.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
geom_point(aes(x = T, y = trait), data = data.raw) + 
labs(title = "Median Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Lifetime") + 
theme_minimal()

gridExtra::grid.arrange(plot1_mean, plot1_med, nrow = 1)
gridExtra::grid.arrange(plot1_mean.inv, plot1_med.inv, nrow = 1)


# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
true_values <- c(log(0.001094), log(0.0488), 0.609, NA)
hist_list1 <- vector("list", length = 4)
par(mfrow = c(2, 2))

for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]] 
  # This column is the mean of each param for each chain
  first_column_list1 <- lapply(param_list1[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  combined_vector1 <- unlist(first_column_list1)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(combined_vector1)-50, max(combined_vector1) + 1, length = 1000)
  # Create a histogram and store it in the list
  hist_list1[[i]] <- hist(combined_vector1, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                        xlab = "Values", col = "lightblue", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "red", lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = "green", lty = 2, lwd = 2)
    }
  if(i == 2){
    lines(x, dnorm(x, 0, sqrt(10)), col = "black", lty = 2, lwd = 2)
    }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = "yellow", lty = 2, lwd = 2)
    }
  }


proportions_list_a1 <- lapply(param_list1[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_a1))

proportions_list_b1 <- lapply(param_list1[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_b1))

proportions_list_c1 <- lapply(param_list1[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_c1))

# Create functions to calculate the RMSE for a, b, and c
calculate_rmse_a_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True value
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_b_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True value
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_c_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
rmse_values_a_lambda <- sapply(param_list1[[1]], calculate_rmse_a_lambda)
rmse_values_b_lambda <- sapply(param_list1[[2]], calculate_rmse_b_lambda)
rmse_values_c_lambda <- sapply(param_list1[[3]], calculate_rmse_c_lambda)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
tab.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(tab.lambda) <- c("a", "b", "c")
colnames(tab.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.lambda[1, 1] <- true.a
tab.lambda[2, 1] <- true.b
tab.lambda[3, 1] <- true.c
tab.lambda[1, 2] <- exp(mean(unlist(lapply(param_list1[[1]], function(matrix) matrix[, 1]))))
tab.lambda[2, 2] <- exp(mean(unlist(lapply(param_list1[[2]], function(matrix) matrix[, 1]))))
tab.lambda[3, 2] <- mean(unlist(lapply(param_list1[[3]], function(matrix) matrix[, 1])))
tab.lambda[1, 3] <- exp(hdi(unlist(lapply(param_list1[[1]], function(matrix) matrix[, 1])))[1])
tab.lambda[2, 3] <- exp(hdi(unlist(lapply(param_list1[[2]], function(matrix) matrix[, 1])))[1])
tab.lambda[3, 3] <- hdi(unlist(lapply(param_list1[[3]], function(matrix) matrix[, 1])))[1]
tab.lambda[1, 4] <- exp(hdi(unlist(lapply(param_list1[[1]], function(matrix) matrix[, 1])))[2])
tab.lambda[2, 4] <- exp(hdi(unlist(lapply(param_list1[[2]], function(matrix) matrix[, 1])))[2])
tab.lambda[3, 4] <- hdi(unlist(lapply(param_list1[[3]], function(matrix) matrix[, 1])))[2]
tab.lambda[1, 5] <- mean(unlist(proportions_list_a1))
tab.lambda[2, 5] <- mean(unlist(proportions_list_b1))
tab.lambda[3, 5] <- mean(unlist(proportions_list_c1))
tab.lambda[1, 6] <- mean(rmse_values_a_lambda)
tab.lambda[2, 6] <- mean(rmse_values_b_lambda)
tab.lambda[3, 6] <- mean(rmse_values_c_lambda)

save(param_list1,
     df.new1,
     mean.xdf1,
     hist_list1,
     plot1_mean, 
     plot1_med,
     plot1_mean.inv, 
     plot1_med.inv,
     proportions_list_a1,
     proportions_list_b1,
     proportions_list_c1,
     rmse_values_a_lambda,
     rmse_values_b_lambda,
     rmse_values_c_lambda,
     tab.lambda,
     file="individual_exponential.RData")

##########################################Mean Lambda #############################################

# Create list of matrices to store parameter summary statistics
param_list2 <- vector("list", length = 4)
for (i in 1:4) {
# Initialize the inner list
inner_list2 <- vector("list", length = b)
# Populate the inner list
for (j in 1:b) {
  # Initialize a 5x4 matrix
  matrix_data2 <- matrix(0, nrow = 5, ncol = 4)
  inner_list2[[j]] <- matrix_data2
}
param_list2[[i]] <- inner_list2
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
    }
  data.raw <- data.list[[i]]
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data.list[[i]]$trait, data.list[[i]]$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
    })
  df.mean <- data.frame(
    T = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)), # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
    )
  # Model
  sink("mean_lambda.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
  mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c)<= epsilon) * epsilon
  trait[i] ~ dexp(mu[i])
  }
      }", file = "mean_lambda.txt")
  # Settings
  parameters <- c("la", "b.l", "c")
  inits <- function(){list(
    # a and b in log space
    la = log(0.001), 
    b.l = log(0.05), 
    c = 0.6
    )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data <- df.mean 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                model.file = "mean_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list2[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list2[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list2[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list2[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    param_list2[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list2[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list2[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list2[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list2[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list2[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list2[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list2[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list2[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list2[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list2[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list2[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
  }
  }
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new2 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 10, to = 35, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new2 as a function at each of those temperatures
xdf2 <- expand.grid(point = xseq, iteration = 1:nrow(df.new2))

# Apply the quadratic equation for each combination
xdf2$value <- apply(xdf2, 1, function(row) {
i <- row["iteration"]
location <- row["point"]
exp(df.new2[i, 4]) * location^2 - exp(df.new2[i, 1]) * location + df.new2[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf2 <- xdf2|> mutate(trunc.eval = ifelse(xdf2$value > epsilon, xdf2$value, epsilon), 
                   trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf2 <- xdf2 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                              avg.value = mean(trunc.eval), 
                                              med.value.inv = median(trunc.inv), 
                                              med.value = median(trunc.eval), 
                                              lower.hdi.inv = hdi(trunc.inv)[1], 
                                              upper.hdi.inv = hdi(trunc.inv)[2], 
                                              lower.hdi = hdi(trunc.eval)[1], 
                                              upper.hdi = hdi(trunc.eval)[2])
# Creating columns of the true curve values and the true inverted curve values
true_curve.inv <- function(x) 1/(true.a * x^2 - true.b * x + true.c)
true_curve <- function(x) (true.a * x^2 - true.b * x + true.c)
mean.xdf2$true_curve <- true_curve(xseq)
mean.xdf2$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
plot2_mean <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = avg.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
labs(title = "Mean Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal()
plot2_med <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = med.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
labs(title = "Median Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal()
# Plot mean lifetime response
plot2_mean.inv <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = avg.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
geom_point(aes(x = T, y = trait), data = data.raw) + 
labs(title = "Mean Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Lifetime") + 
theme_minimal()
# Plot median lifetime response
plot2_med.inv <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = med.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
geom_point(aes(x = T, y = trait), data = data.raw) + 
labs(title = "Median Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Lifetime") + 
theme_minimal()

gridExtra::grid.arrange(plot2_mean, plot2_med, nrow = 1)
gridExtra::grid.arrange(plot2_mean.inv, plot2_med.inv, nrow = 1)

# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
true_values <- c(log(true.a), log(true.b), true.c, NA)
hist_list2 <- vector("list", length = 4)
par(mfrow = c(2, 2))
for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]] 
  # This column is the mean of each param for each chain
  first_column_list2 <- lapply(param_list2[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  combined_vector2 <- unlist(first_column_list2)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(combined_vector2)-50, max(combined_vector2) + 1, length = 1000)
  # Create a histogram and store it in the list
  hist_list2[[i]] <- hist(combined_vector2, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                        xlab = "Values", col = "lightblue", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "red", lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = "green", lty = 2, lwd = 2)
    }
  if(i == 2){
    lines(x, dnorm(x, 0, sqrt(10)), col = "black", lty = 2, lwd = 2)
    }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = "yellow", lty = 2, lwd = 2)
    }

  }


proportions_list_a2 <- lapply(param_list2[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_a2))

proportions_list_b2 <- lapply(param_list2[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_b2))

proportions_list_c2 <- lapply(param_list2[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_c2))

calculate_rmse_a_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_b_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1] 
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_c_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
rmse_values_a_lambda_mean <- sapply(param_list2[[1]], calculate_rmse_a_lambda_mean)
rmse_values_b_lambda_mean <- sapply(param_list2[[2]], calculate_rmse_b_lambda_mean)
rmse_values_c_lambda_mean <- sapply(param_list2[[3]], calculate_rmse_c_lambda_mean)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
tab.lambda.mean <- matrix(0, nrow = 3, ncol = 6)
row.names(tab.lambda.mean) <- c("a", "b", "c")
colnames(tab.lambda.mean) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.lambda.mean[1, 1] <- true.a
tab.lambda.mean[2, 1] <- true.b
tab.lambda.mean[3, 1] <- true.c
tab.lambda.mean[1, 2] <- exp(mean(unlist(lapply(param_list2[[1]], function(matrix) matrix[, 1]))))
tab.lambda.mean[2, 2] <- exp(mean(unlist(lapply(param_list2[[2]], function(matrix) matrix[, 1]))))
tab.lambda.mean[3, 2] <- mean(unlist(lapply(param_list2[[3]], function(matrix) matrix[, 1])))
tab.lambda.mean[1, 3] <- exp(hdi(unlist(lapply(param_list2[[1]], function(matrix) matrix[, 1])))[1])
tab.lambda.mean[2, 3] <- exp(hdi(unlist(lapply(param_list2[[2]], function(matrix) matrix[, 1])))[1])
tab.lambda.mean[3, 3] <- hdi(unlist(lapply(param_list2[[3]], function(matrix) matrix[, 1])))[1]
tab.lambda.mean[1, 4] <- exp(hdi(unlist(lapply(param_list2[[1]], function(matrix) matrix[, 1])))[2])
tab.lambda.mean[2, 4] <- exp(hdi(unlist(lapply(param_list2[[2]], function(matrix) matrix[, 1])))[2])
tab.lambda.mean[3, 4] <- hdi(unlist(lapply(param_list2[[3]], function(matrix) matrix[, 1])))[2]
tab.lambda.mean[1, 5] <- mean(unlist(proportions_list_a2))
tab.lambda.mean[2, 5] <- mean(unlist(proportions_list_b2))
tab.lambda.mean[3, 5] <- mean(unlist(proportions_list_c2))
tab.lambda.mean[1, 6] <- mean(rmse_values_a_lambda_mean)
tab.lambda.mean[2, 6] <- mean(rmse_values_b_lambda_mean)
tab.lambda.mean[3, 6] <- mean(rmse_values_c_lambda_mean)

save(param_list2,
     df.new2,
     mean.xdf2,
     hist_list2,
     plot2_mean, 
     plot2_med,
     plot2_mean.inv, 
     plot2_med.inv,
     proportions_list_a2,
     proportions_list_b2,
     proportions_list_c2,
     rmse_values_a_lambda_mean,
     rmse_values_b_lambda_mean,
     rmse_values_c_lambda_mean,
     tab.lambda.mean,
     file="mean_exponential.RData")

##########################################Inverse Lambda #############################################

# Create list of matrices to store parameter summary statistics
param_list3 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  inner_list3 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    matrix_data3 <- matrix(0, nrow = 5, ncol = 4)
    inner_list3[[j]] <- matrix_data3
  }
  param_list3[[i]] <- inner_list3
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("inv_lambda.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  # And nonzero
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(4000)
  sig2 <- sig^2
  tau <- 1/sig2
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c)<= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau) T(0, )
  }
}", file = "inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c", "sig")
  inits <- function(){list(
    # a and b in log space
    la = log(0.001), 
    b.l = log(0.05), 
    c = 0.6, 
    sig = 2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data <- data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "inv_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list3[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list3[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list3[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list3[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    param_list3[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list3[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list3[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list3[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list3[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list3[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list3[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list3[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list3[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list3[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list3[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list3[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    param_list3[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list3[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list3[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list3[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new3 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 10, to = 35, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new3 as a function at each of those temperatures
xdf3 <- expand.grid(point = xseq, iteration = 1:nrow(df.new3))

# Apply the quadratic equation for each combination
xdf3$value <- apply(xdf3, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new3[i, 4]) * location^2 - exp(df.new3[i, 1]) * location + df.new3[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf3 <- xdf3|> mutate(trunc.eval = ifelse(xdf3$value > epsilon, xdf3$value, epsilon), 
                     trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf3 <- xdf3 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                avg.value = mean(trunc.eval), 
                                                med.value.inv = median(trunc.inv), 
                                                med.value = median(trunc.eval), 
                                                lower.hdi.inv = hdi(trunc.inv)[1], 
                                                upper.hdi.inv = hdi(trunc.inv)[2], 
                                                lower.hdi = hdi(trunc.eval)[1], 
                                                upper.hdi = hdi(trunc.eval)[2]) 
# Creating columns of the true curve values and the true inverted curve values
true_curve.inv <- function(x) 1/(true.a * x^2 - true.b * x + true.c)
true_curve <- function(x) (true.a * x^2 - true.b * x + true.c)
mean.xdf3$true_curve <- true_curve(xseq)
mean.xdf3$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
plot3_mean <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot3_med <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot3_mean.inv <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()
# Plot median lifetime response
plot3_med.inv <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()

gridExtra::grid.arrange(plot3_mean, plot3_med, nrow = 1)
gridExtra::grid.arrange(plot3_mean.inv, plot3_med.inv, nrow = 1)


# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
true_values <- c(log(true.a), log(true.b), true.c, NA, NA)
hist_list3 <- vector("list", length = 5)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  first_column_list3 <- lapply(param_list3[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  combined_vector3 <- unlist(first_column_list3)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(combined_vector3)-1, max(combined_vector3) + 1, length = 1000)
  # Create a histogram and store it in the list
  hist_list3[[i]] <- hist(combined_vector3, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "lightblue", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "red", lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = "green", lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 0, sqrt(10)), col = "black", lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = "yellow", lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 4000), col = "purple", lty = 2, lwd = 2)
  }
}


proportions_list_a3 <- lapply(param_list3[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_a3))

proportions_list_b3 <- lapply(param_list3[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_b3))

proportions_list_c3 <- lapply(param_list3[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_c3))

calculate_rmse_a_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_b_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_c_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
rmse_values_a_lambda_inv <- sapply(param_list3[[1]], calculate_rmse_a_lambda_inv)
rmse_values_b_lambda_inv <- sapply(param_list3[[2]], calculate_rmse_b_lambda_inv)
rmse_values_c_lambda_inv <- sapply(param_list3[[3]], calculate_rmse_c_lambda_inv)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
tab.inv.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(tab.inv.lambda) <- c("a", "b", "c")
colnames(tab.inv.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.inv.lambda[1, 1] <- true.a
tab.inv.lambda[2, 1] <- true.b
tab.inv.lambda[3, 1] <- true.c
tab.inv.lambda[1, 2] <- exp(mean(unlist(lapply(param_list3[[1]], function(matrix) matrix[, 1]))))
tab.inv.lambda[2, 2] <- exp(mean(unlist(lapply(param_list3[[2]], function(matrix) matrix[, 1]))))
tab.inv.lambda[3, 2] <- mean(unlist(lapply(param_list3[[3]], function(matrix) matrix[, 1])))
tab.inv.lambda[1, 3] <- exp(hdi(unlist(lapply(param_list3[[1]], function(matrix) matrix[, 1])))[1])
tab.inv.lambda[2, 3] <- exp(hdi(unlist(lapply(param_list3[[2]], function(matrix) matrix[, 1])))[1])
tab.inv.lambda[3, 3] <- hdi(unlist(lapply(param_list3[[3]], function(matrix) matrix[, 1])))[1]
tab.inv.lambda[1, 4] <- exp(hdi(unlist(lapply(param_list3[[1]], function(matrix) matrix[, 1])))[2])
tab.inv.lambda[2, 4] <- exp(hdi(unlist(lapply(param_list3[[2]], function(matrix) matrix[, 1])))[2])
tab.inv.lambda[3, 4] <- hdi(unlist(lapply(param_list3[[3]], function(matrix) matrix[, 1])))[2]
tab.inv.lambda[1, 5] <- mean(unlist(proportions_list_a3))
tab.inv.lambda[2, 5] <- mean(unlist(proportions_list_b3))
tab.inv.lambda[3, 5] <- mean(unlist(proportions_list_c3))
tab.inv.lambda[1, 6] <- mean(rmse_values_a_lambda_inv)
tab.inv.lambda[2, 6] <- mean(rmse_values_b_lambda_inv)
tab.inv.lambda[3, 6] <- mean(rmse_values_c_lambda_inv)

save(param_list3,
     df.new3,
     mean.xdf3,
     hist_list3,
     plot3_mean, 
     plot3_med,
     plot3_mean.inv, 
     plot3_med.inv,
     proportions_list_a3,
     proportions_list_b3,
     proportions_list_c3,
     rmse_values_a_lambda_inv,
     rmse_values_b_lambda_inv,
     rmse_values_c_lambda_inv,
     tab.inv.lambda,
     file="inverse_exponential.RData")

##########################################Mean Inverse Lambda #############################################

# Create list of matrices to store parameter summary statistics
param_list4 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  inner_list4 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    matrix_data4 <- matrix(0, nrow = 5, ncol = 4)
    inner_list4[[j]] <- matrix_data4
  }
  param_list4[[i]] <- inner_list4
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data.raw <- data.list[[i]] 
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data.raw$trait, data.raw$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)), # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("mean_inv_lambda.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(0.5)
  sig2 <- sig^2
  tau <- 1/sig2
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c)<= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau) T(0, )
  }
}", file = "mean_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c", "sig")
  inits <- function(){list(
    # a and b in log space
    la = log(0.001), 
    b.l = log(0.05), 
    c = 0.6, 
    sig = 2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "mean_inv_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list4[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list4[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list4[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list4[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    param_list4[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list4[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list4[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list4[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list4[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list4[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list4[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list4[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list4[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list4[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list4[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list4[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    param_list4[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list4[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list4[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list4[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 10, to = 35, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new4 as a function at each of those temperatures
xdf4.2 <- expand.grid(point = xseq, iteration = 1:nrow(df.new4))

# Apply the quadratic equation for each combination
xdf4.2$value <- apply(xdf4.2, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new4[i, 4]) * location^2 - exp(df.new4[i, 1]) * location + df.new4[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf4.2 <- xdf4.2|> mutate(trunc.eval = ifelse(xdf4.2$value > epsilon, xdf4.2$value, epsilon), 
                          trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf4 <- xdf4.2 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                  avg.value = mean(trunc.eval), 
                                                  med.value.inv = median(trunc.inv), 
                                                  med.value = median(trunc.eval), 
                                                  lower.hdi.inv = hdi(trunc.inv)[1], 
                                                  upper.hdi.inv = hdi(trunc.inv)[2], 
                                                  lower.hdi = hdi(trunc.eval)[1], 
                                                  upper.hdi = hdi(trunc.eval)[2])
# Creating columns of the true curve values and the true inverted curve values
true_curve.inv <- function(x) 1/(true.a * x^2 - true.b * x + true.c)
true_curve <- function(x) (true.a * x^2 - true.b * x + true.c)
mean.xdf4$true_curve <- true_curve(xseq)
mean.xdf4$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
plot4_mean <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot4_med <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot4_mean.inv <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()
# Plot median lifetime response
plot4_med.inv <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()

gridExtra::grid.arrange(plot4_mean, plot4_med, nrow = 1)
gridExtra::grid.arrange(plot4_mean.inv, plot4_med.inv, nrow = 1)

# Create a list to store histograms of the posterior distribution of each parameter
hist_list4 <- vector("list", length = 4)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  first_column_list4 <- lapply(param_list4[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  combined_vector4 <- unlist(first_column_list4)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(combined_vector4)-1, max(combined_vector4) + 1, length = 1000)
  # Create a histogram and store it in the list
  hist_list4[[i]] <- hist(combined_vector4, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "lightblue", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "red", lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = "green", lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 0, sqrt(10)), col = "black", lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = "yellow", lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 0.5), col = "purple", lty = 2, lwd = 2)
  }
}


proportions_list_a4 <- lapply(param_list4[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_a4))

proportions_list_b4 <- lapply(param_list4[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_b4))

proportions_list_c4 <- lapply(param_list4[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_c4))

calculate_rmse_a_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_b_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_c_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
rmse_values_a_lambda_mean_inv <- sapply(param_list4[[1]], calculate_rmse_a_lambda_mean_inv)
rmse_values_b_lambda_mean_inv <- sapply(param_list4[[2]], calculate_rmse_b_lambda_mean_inv)
rmse_values_c_lambda_mean_inv <- sapply(param_list4[[3]], calculate_rmse_c_lambda_mean_inv)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
tab.mean.inv.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(tab.mean.inv.lambda) <- c("a", "b", "c")
colnames(tab.mean.inv.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.mean.inv.lambda[1, 1] <- true.a
tab.mean.inv.lambda[2, 1] <- true.b
tab.mean.inv.lambda[3, 1] <- true.c
tab.mean.inv.lambda[1, 2] <- exp(mean(unlist(lapply(param_list4[[1]], function(matrix) matrix[, 1]))))
tab.mean.inv.lambda[2, 2] <- exp(mean(unlist(lapply(param_list4[[2]], function(matrix) matrix[, 1]))))
tab.mean.inv.lambda[3, 2] <- mean(unlist(lapply(param_list4[[3]], function(matrix) matrix[, 1])))
tab.mean.inv.lambda[1, 3] <- exp(hdi(unlist(lapply(param_list4[[1]], function(matrix) matrix[, 1])))[1])
tab.mean.inv.lambda[2, 3] <- exp(hdi(unlist(lapply(param_list4[[2]], function(matrix) matrix[, 1])))[1])
tab.mean.inv.lambda[3, 3] <- hdi(unlist(lapply(param_list4[[3]], function(matrix) matrix[, 1])))[1]
tab.mean.inv.lambda[1, 4] <- exp(hdi(unlist(lapply(param_list4[[1]], function(matrix) matrix[, 1])))[2])
tab.mean.inv.lambda[2, 4] <- exp(hdi(unlist(lapply(param_list4[[2]], function(matrix) matrix[, 1])))[2])
tab.mean.inv.lambda[3, 4] <- hdi(unlist(lapply(param_list4[[3]], function(matrix) matrix[, 1])))[2]
tab.mean.inv.lambda[1, 5] <- mean(unlist(proportions_list_a4))
tab.mean.inv.lambda[2, 5] <- mean(unlist(proportions_list_b4))
tab.mean.inv.lambda[3, 5] <- mean(unlist(proportions_list_c4))
tab.mean.inv.lambda[1, 6] <- mean(rmse_values_a_lambda_mean_inv)
tab.mean.inv.lambda[2, 6] <- mean(rmse_values_b_lambda_mean_inv)
tab.mean.inv.lambda[3, 6] <- mean(rmse_values_c_lambda_mean_inv)

save(param_list4,
     df.new4,
     mean.xdf4,
     hist_list4,
     plot4_mean, 
     plot4_med,
     plot4_mean.inv, 
     plot4_med.inv,
     proportions_list_a4,
     proportions_list_b4,
     proportions_list_c4,
     rmse_values_a_lambda_mean_inv,
     rmse_values_b_lambda_mean_inv,
     rmse_values_c_lambda_mean_inv,
     tab.mean.inv.lambda,
     file="mean_inverse_exponential.RData")


##########################################Inverse Mean Lambda #############################################

# Create list of matrices to store parameter summary statistics
param_list5 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  inner_list5 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    matrix_data5 <- matrix(0, nrow = 5, ncol = 4)
    inner_list5[[j]] <- matrix_data5
  }
  param_list5[[i]] <- inner_list5
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data.raw <- data.list[[i]] 
  data.raw$inv.trait <- 1/(data.raw$trait)
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data.raw$inv.trait, data.raw$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <-  data.frame(
    T = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)), # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- data$trait
  1/trait
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("inv_mean_lambda.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(400)
  sig2 <- sig^2
  tau <- 1/sig2
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c)<= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau) T(0, )
  }
}", file = "inv_mean_lambda.txt")
  # Settings
  parameters <- c("la", "b.l", "c", "sig")
  inits <- function(){list(
    # a and b in log space
    la = log(0.001), 
    b.l = log(0.05), 
    c = 0.6, 
    sig = 2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "inv_mean_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list5[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list5[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list5[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list5[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    param_list5[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list5[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list5[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list5[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list5[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list5[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list5[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list5[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list5[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list5[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list5[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list5[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    param_list5[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list5[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list5[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list5[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 10, to = 35, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new5 as a function at each of those temperatures
xdf5 <- expand.grid(point = xseq, iteration = 1:nrow(df.new5))

# Apply the quadratic equation for each combination
xdf5$value <- apply(xdf5, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new5[i, 4]) * location^2 - exp(df.new5[i, 1]) * location + df.new5[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf5 <- xdf5|> mutate(trunc.eval = ifelse(xdf5$value > epsilon, xdf5$value, epsilon), 
                      trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf5 <- xdf5 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                avg.value = mean(trunc.eval), 
                                                med.value.inv = median(trunc.inv), 
                                                med.value = median(trunc.eval), 
                                                lower.hdi.inv = hdi(trunc.inv)[1], 
                                                upper.hdi.inv = hdi(trunc.inv)[2], 
                                                lower.hdi = hdi(trunc.eval)[1], 
                                                upper.hdi = hdi(trunc.eval)[2])
# Creating columns of the true curve values and the true inverted curve values
true_curve.inv <- function(x) 1/(true.a * x^2 - true.b * x + true.c)
true_curve <- function(x) (true.a * x^2 - true.b * x + true.c)
mean.xdf5$true_curve <- true_curve(xseq)
mean.xdf5$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
plot5_mean <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot5_med <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot5_mean.inv <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()
# Plot median lifetime response
plot5_med.inv <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()

gridExtra::grid.arrange(plot5_mean, plot5_med, nrow = 1)
gridExtra::grid.arrange(plot5_mean.inv, plot5_med.inv, nrow = 1)

# Create a list to store histograms of the posterior distribution of each parameter
hist_list5 <- vector("list", length = 5)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  first_column_list5 <- lapply(param_list5[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  combined_vector5 <- unlist(first_column_list5)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(combined_vector5)-1, max(combined_vector5) + 1, length = 1000)
  # Create a histogram and store it in the list
  hist_list5[[i]] <- hist(combined_vector5, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "lightblue", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "red", lwd = 2)
  # Add prior based on which parameter 
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = "green", lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 0, sqrt(10)), col = "black", lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = "yellow", lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 400), col = "purple", lty = 2, lwd = 2)
  }
}

proportions_list_a5 <-  lapply(param_list5[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_a5))

proportions_list_b5 <-  lapply(param_list5[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_b5))

proportions_list_c5 <-  lapply(param_list5[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_c5))

calculate_rmse_a_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_b_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_c_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
rmse_values_a_lambda_inv_mean <- sapply(param_list5[[1]], calculate_rmse_a_lambda_inv_mean)
rmse_values_b_lambda_inv_mean <- sapply(param_list5[[2]], calculate_rmse_b_lambda_inv_mean)
rmse_values_c_lambda_inv_mean <- sapply(param_list5[[3]], calculate_rmse_c_lambda_inv_mean)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
tab.inv.mean.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(tab.inv.mean.lambda) <- c("a", "b", "c")
colnames(tab.inv.mean.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.inv.mean.lambda[1, 1] <- true.a
tab.inv.mean.lambda[2, 1] <- true.b
tab.inv.mean.lambda[3, 1] <- true.c
tab.inv.mean.lambda[1, 2] <- exp(mean(unlist(lapply(param_list5[[1]], function(matrix) matrix[, 1]))))
tab.inv.mean.lambda[2, 2] <- exp(mean(unlist(lapply(param_list5[[2]], function(matrix) matrix[, 1]))))
tab.inv.mean.lambda[3, 2] <- mean(unlist(lapply(param_list5[[3]], function(matrix) matrix[, 1])))
tab.inv.mean.lambda[1, 3] <- exp(hdi(unlist(lapply(param_list5[[1]], function(matrix) matrix[, 1])))[1])
tab.inv.mean.lambda[2, 3] <- exp(hdi(unlist(lapply(param_list5[[2]], function(matrix) matrix[, 1])))[1])
tab.inv.mean.lambda[3, 3] <- hdi(unlist(lapply(param_list5[[3]], function(matrix) matrix[, 1])))[1]
tab.inv.mean.lambda[1, 4] <- exp(hdi(unlist(lapply(param_list5[[1]], function(matrix) matrix[, 1])))[2])
tab.inv.mean.lambda[2, 4] <- exp(hdi(unlist(lapply(param_list5[[2]], function(matrix) matrix[, 1])))[2])
tab.inv.mean.lambda[3, 4] <- hdi(unlist(lapply(param_list5[[3]], function(matrix) matrix[, 1])))[2]
tab.inv.mean.lambda[1, 5] <- mean(unlist(proportions_list_a5))
tab.inv.mean.lambda[2, 5] <- mean(unlist(proportions_list_b5))
tab.inv.mean.lambda[3, 5] <- mean(unlist(proportions_list_c5))
tab.inv.mean.lambda[1, 6] <- mean(rmse_values_a_lambda_inv_mean)
tab.inv.mean.lambda[2, 6] <- mean(rmse_values_b_lambda_inv_mean)
tab.inv.mean.lambda[3, 6] <- mean(rmse_values_c_lambda_inv_mean)


save(param_list5,
     df.new5,
     mean.xdf5,
     hist_list5,
     plot5_mean, 
     plot5_med,
     plot5_mean.inv, 
     plot5_med.inv,
     proportions_list_a5,
     proportions_list_b5,
     proportions_list_c5,
     rmse_values_a_lambda_inv_mean,
     rmse_values_b_lambda_inv_mean,
     rmse_values_c_lambda_inv_mean,
     tab.inv.mean.lambda,
     file="inverse_mean_exponential.RData")



##########################################Inverse Lambda NOT TRUNCATED #############################################

# Create list of matrices to store parameter summary statistics
param_list3_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  inner_list3_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    matrix_data3_notrunc <- matrix(0, nrow = 5, ncol = 4)
    inner_list3_notrunc[[j]] <- matrix_data3_notrunc
  }
  param_list3_notrunc[[i]] <- inner_list3_notrunc
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("inv_lambda_nt.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  tau ~ dexp(0.1)
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c)<= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "inv_lambda_nt.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c", "tau")
  inits <- function(){list(
    # a and b in log space
    la = log(0.001), 
    b.l = log(0.05), 
    c = 0.6, 
    tau = 0.2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data.raw <-  data.list[[i]]
  data <- data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "inv_lambda_nt.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list3_notrunc[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list3_notrunc[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list3_notrunc[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list3_notrunc[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    param_list3_notrunc[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list3_notrunc[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list3_notrunc[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list3_notrunc[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list3_notrunc[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list3_notrunc[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list3_notrunc[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list3_notrunc[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list3_notrunc[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list3_notrunc[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list3_notrunc[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list3_notrunc[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    param_list3_notrunc[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list3_notrunc[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list3_notrunc[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list3_notrunc[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new3.NT = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 10, to = 35, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new3.NT as a function at each of those temperatures
xdf3.NT <- expand.grid(point = xseq, iteration = 1:nrow(df.new3.NT))

# Apply the quadratic equation for each combination
xdf3.NT$value <- apply(xdf3.NT, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new3.NT[i, 4]) * location^2 - exp(df.new3.NT[i, 1]) * location + df.new3.NT[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf3.NT <-  xdf3.NT|> mutate(trunc.eval = ifelse(xdf3.NT$value > epsilon, xdf3.NT$value, epsilon), 
                           trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf3.NT <- xdf3.NT |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                      avg.value = mean(trunc.eval), 
                                                      med.value.inv = median(trunc.inv), 
                                                      med.value = median(trunc.eval), 
                                                      lower.hdi.inv = hdi(trunc.inv)[1], 
                                                      upper.hdi.inv = hdi(trunc.inv)[2], 
                                                      lower.hdi = hdi(trunc.eval)[1], 
                                                      upper.hdi = hdi(trunc.eval)[2])
# Creating columns of the true curve values and the true inverted curve values
true_curve.inv <- function(x) 1/(true.a * x^2 - true.b * x + true.c)
true_curve <- function(x) (true.a * x^2 - true.b * x + true.c)
mean.xdf3.NT$true_curve <- true_curve(xseq)
mean.xdf3.NT$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
plot3_mean_NT <- ggplot(mean.xdf3.NT, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot3_med_NT <- ggplot(mean.xdf3.NT, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot3_mean.inv_NT <- ggplot(mean.xdf3.NT, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()
# Plot median lifetime response
plot3_med.inv_NT <- ggplot(mean.xdf3.NT, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()

gridExtra::grid.arrange(plot3_mean_NT, plot3_med_NT, nrow = 1)
gridExtra::grid.arrange(plot3_mean.inv_NT, plot3_med.inv_NT, nrow = 1)

# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
true_values <- c(log(true.a), log(true.b), true.c, NA, NA)
hist_list3_notrunc <- vector("list", length = 5)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  first_column_list3_notrunc <- lapply(param_list3_notrunc[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  combined_vector3_notrunc <- unlist(first_column_list3_notrunc)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(combined_vector3_notrunc)-1, max(combined_vector3_notrunc) + 1, length = 1000)
  # Create a histogram and store it in the list
  hist_list3_notrunc[[i]] <- hist(combined_vector3_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                  xlab = "Values", col = "lightblue", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "red", lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = "green", lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 0, sqrt(10)), col = "black", lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = "yellow", lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 0.1), col = "purple", lty = 2, lwd = 2)
  }
  
}


proportions_list_a3_notrunc <-  lapply(param_list3_notrunc[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_a3_notrunc))

proportions_list_b3_notrunc <-  lapply(param_list3_notrunc[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_b3_notrunc))

proportions_list_c3_notrunc <-  lapply(param_list3_notrunc[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_c3_notrunc))

calculate_rmse_a_lambda_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_b_lambda_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_c_lambda_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
rmse_values_a_lambda_inv_notrunc <- sapply(param_list3_notrunc[[1]], calculate_rmse_a_lambda_inv_notrunc)
rmse_values_b_lambda_inv_notrunc <- sapply(param_list3_notrunc[[2]], calculate_rmse_b_lambda_inv_notrunc)
rmse_values_c_lambda_inv_notrunc <- sapply(param_list3_notrunc[[3]], calculate_rmse_c_lambda_inv_notrunc)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
tab.inv.lambda_notrunc <- matrix(0, nrow = 3, ncol = 6)
row.names(tab.inv.lambda_notrunc) <- c("a", "b", "c")
colnames(tab.inv.lambda_notrunc) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.inv.lambda_notrunc[1, 1] <- true.a
tab.inv.lambda_notrunc[2, 1] <- true.b
tab.inv.lambda_notrunc[3, 1] <- true.c
tab.inv.lambda_notrunc[1, 2] <- exp(mean(unlist(lapply(param_list3_notrunc[[1]], function(matrix) matrix[, 1]))))
tab.inv.lambda_notrunc[2, 2] <- exp(mean(unlist(lapply(param_list3_notrunc[[2]], function(matrix) matrix[, 1]))))
tab.inv.lambda_notrunc[3, 2] <- mean(unlist(lapply(param_list3_notrunc[[3]], function(matrix) matrix[, 1])))
tab.inv.lambda_notrunc[1, 3] <- exp(hdi(unlist(lapply(param_list3_notrunc[[1]], function(matrix) matrix[, 1])))[1])
tab.inv.lambda_notrunc[2, 3] <- exp(hdi(unlist(lapply(param_list3_notrunc[[2]], function(matrix) matrix[, 1])))[1])
tab.inv.lambda_notrunc[3, 3] <- hdi(unlist(lapply(param_list3_notrunc[[3]], function(matrix) matrix[, 1])))[1]
tab.inv.lambda_notrunc[1, 4] <- exp(hdi(unlist(lapply(param_list3_notrunc[[1]], function(matrix) matrix[, 1])))[2])
tab.inv.lambda_notrunc[2, 4] <- exp(hdi(unlist(lapply(param_list3_notrunc[[2]], function(matrix) matrix[, 1])))[2])
tab.inv.lambda_notrunc[3, 4] <- hdi(unlist(lapply(param_list3_notrunc[[3]], function(matrix) matrix[, 1])))[2]
tab.inv.lambda_notrunc[1, 5] <- mean(unlist(proportions_list_a3_notrunc))
tab.inv.lambda_notrunc[2, 5] <- mean(unlist(proportions_list_b3_notrunc))
tab.inv.lambda_notrunc[3, 5] <- mean(unlist(proportions_list_c3_notrunc))
tab.inv.lambda_notrunc[1, 6] <- mean(rmse_values_a_lambda_inv_notrunc)
tab.inv.lambda_notrunc[2, 6] <- mean(rmse_values_b_lambda_inv_notrunc)
tab.inv.lambda_notrunc[3, 6] <- mean(rmse_values_c_lambda_inv_notrunc)

save(param_list3_notrunc,
     df.new3.NT,
     mean.xdf3.NT,
     hist_list3_notrunc,
     plot3_mean_NT, 
     plot3_med_NT,
     plot3_mean.inv_NT, 
     plot3_med.inv_NT,
     proportions_list_a3_notrunc,
     proportions_list_b3_notrunc,
     proportions_list_c3_notrunc,
     rmse_values_a_lambda_inv_notrunc,
     rmse_values_b_lambda_inv_notrunc,
     rmse_values_c_lambda_inv_notrunc,
     tab.inv.lambda_notrunc,
     file="inverse_notrunc.RData")

##########################################Mean Inverse Lambda NOT TRUNCATED #############################################

# Create list of matrices to store parameter summary statistics
param_list4_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  inner_list4_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    matrix_data4_notrunc <- matrix(0, nrow = 5, ncol = 4)
    inner_list4_notrunc[[j]] <- matrix_data4_notrunc
  }
  param_list4_notrunc[[i]] <- inner_list4_notrunc
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data.raw <- data.list[[i]]
  data <- data.list[[i]] 
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <-  data.frame(
    T = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)), # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("mean_inv_lambda_nt.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  tau ~ dexp(0.5)
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c)<= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "mean_inv_lambda_nt.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c", "tau")
  inits <- function(){list(
    # a and b in log space
    la = log(0.001), 
    b.l = log(0.05), 
    c = 0.6, 
    tau = 0.2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "mean_inv_lambda_nt.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list4_notrunc[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list4_notrunc[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list4_notrunc[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list4_notrunc[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    param_list4_notrunc[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list4_notrunc[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list4_notrunc[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list4_notrunc[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list4_notrunc[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list4_notrunc[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list4_notrunc[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list4_notrunc[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list4_notrunc[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list4_notrunc[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list4_notrunc[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list4_notrunc[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    param_list4_notrunc[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list4_notrunc[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list4_notrunc[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list4_notrunc[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4.NT = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 10, to = 35, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new4.NT as a function at each of those temperatures
xdf4.NT <- expand.grid(point = xseq, iteration = 1:nrow(df.new4.NT))

# Apply the quadratic equation for each combination
xdf4.NT$value <- apply(xdf4.NT, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new4.NT[i, 4]) * location^2 - exp(df.new4.NT[i, 1]) * location + df.new4.NT[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf4.NT <-  xdf4.NT|> mutate(trunc.eval = ifelse(xdf4.NT$value > epsilon, xdf4.NT$value, epsilon), 
                           trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf4.NT <- xdf4.NT |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                      avg.value = mean(trunc.eval), 
                                                      med.value.inv = median(trunc.inv), 
                                                      med.value = median(trunc.eval), 
                                                      lower.hdi.inv = hdi(trunc.inv)[1], 
                                                      upper.hdi.inv = hdi(trunc.inv)[2], 
                                                      lower.hdi = hdi(trunc.eval)[1], 
                                                      upper.hdi = hdi(trunc.eval)[2])
# Creating columns of the true curve values and the true inverted curve values
true_curve.inv <- function(x) 1/(true.a * x^2 - true.b * x + true.c)
true_curve <- function(x) (true.a * x^2 - true.b * x + true.c)
mean.xdf4.NT$true_curve <- true_curve(xseq)
mean.xdf4.NT$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
plot4_mean_NT <- ggplot(mean.xdf4.NT, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot4_med_NT <- ggplot(mean.xdf4.NT, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot4_mean.inv_NT <- ggplot(mean.xdf4.NT, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()
# Plot median lifetime response
plot4_med.inv_NT <- ggplot(mean.xdf4.NT, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()

gridExtra::grid.arrange(plot4_mean_NT, plot4_med_NT, nrow = 1)
gridExtra::grid.arrange(plot4_mean.inv_NT, plot4_med.inv_NT, nrow = 1)


# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
hist_list4_notrunc <- vector("list", length = 4)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  first_column_list4_notrunc <- lapply(param_list4_notrunc[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  combined_vector4_notrunc <- unlist(first_column_list4_notrunc)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(combined_vector4_notrunc)-1, max(combined_vector4_notrunc) + 1, length = 1000)
  # Create a histogram and store it in the list
  hist_list4_notrunc[[i]] <- hist(combined_vector4_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                  xlab = "Values", col = "lightblue", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "red", lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = "green", lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 0, sqrt(10)), col = "black", lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = "yellow", lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 0.5), col = "purple", lty = 2, lwd = 2)
  }
}


proportions_list_a4_notrunc <-  lapply(param_list4_notrunc[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_a4_notrunc))

proportions_list_b4_notrunc <-  lapply(param_list4_notrunc[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_b4_notrunc))

proportions_list_c4_notrunc <-  lapply(param_list4_notrunc[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_c4_notrunc))

calculate_rmse_a_lambda_mean_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_b_lambda_mean_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_c_lambda_mean_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
rmse_values_a_lambda_mean_inv_notrunc <- sapply(param_list4_notrunc[[1]], calculate_rmse_a_lambda_mean_inv_notrunc)
rmse_values_b_lambda_mean_inv_notrunc <- sapply(param_list4_notrunc[[2]], calculate_rmse_b_lambda_mean_inv_notrunc)
rmse_values_c_lambda_mean_inv_notrunc <- sapply(param_list4_notrunc[[3]], calculate_rmse_c_lambda_mean_inv_notrunc)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
tab.mean.inv.lambda_notrunc <- matrix(0, nrow = 3, ncol = 6)
row.names(tab.mean.inv.lambda_notrunc) <- c("a", "b", "c")
colnames(tab.mean.inv.lambda_notrunc) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.mean.inv.lambda_notrunc[1, 1] <- true.a
tab.mean.inv.lambda_notrunc[2, 1] <- true.b
tab.mean.inv.lambda_notrunc[3, 1] <- true.c
tab.mean.inv.lambda_notrunc[1, 2] <- exp(mean(unlist(lapply(param_list4_notrunc[[1]], function(matrix) matrix[, 1]))))
tab.mean.inv.lambda_notrunc[2, 2] <- exp(mean(unlist(lapply(param_list4_notrunc[[2]], function(matrix) matrix[, 1]))))
tab.mean.inv.lambda_notrunc[3, 2] <- mean(unlist(lapply(param_list4_notrunc[[3]], function(matrix) matrix[, 1])))
tab.mean.inv.lambda_notrunc[1, 3] <- exp(hdi(unlist(lapply(param_list4_notrunc[[1]], function(matrix) matrix[, 1])))[1])
tab.mean.inv.lambda_notrunc[2, 3] <- exp(hdi(unlist(lapply(param_list4_notrunc[[2]], function(matrix) matrix[, 1])))[1])
tab.mean.inv.lambda_notrunc[3, 3] <- hdi(unlist(lapply(param_list4_notrunc[[3]], function(matrix) matrix[, 1])))[1]
tab.mean.inv.lambda_notrunc[1, 4] <- exp(hdi(unlist(lapply(param_list4_notrunc[[1]], function(matrix) matrix[, 1])))[2])
tab.mean.inv.lambda_notrunc[2, 4] <- exp(hdi(unlist(lapply(param_list4_notrunc[[2]], function(matrix) matrix[, 1])))[2])
tab.mean.inv.lambda_notrunc[3, 4] <- hdi(unlist(lapply(param_list4_notrunc[[3]], function(matrix) matrix[, 1])))[2]
tab.mean.inv.lambda_notrunc[1, 5] <- mean(unlist(proportions_list_a4_notrunc))
tab.mean.inv.lambda_notrunc[2, 5] <- mean(unlist(proportions_list_b4_notrunc))
tab.mean.inv.lambda_notrunc[3, 5] <- mean(unlist(proportions_list_c4_notrunc))
tab.mean.inv.lambda_notrunc[1, 6] <- mean(rmse_values_a_lambda_mean_inv_notrunc)
tab.mean.inv.lambda_notrunc[2, 6] <- mean(rmse_values_b_lambda_mean_inv_notrunc)
tab.mean.inv.lambda_notrunc[3, 6] <- mean(rmse_values_c_lambda_mean_inv_notrunc)

save(param_list4_notrunc,
     df.new4.NT,
     mean.xdf4.NT,
     hist_list4_notrunc,
     plot4_mean_NT, 
     plot4_med_NT,
     plot4_mean.inv_NT, 
     plot4_med.inv_NT,
     proportions_list_a4_notrunc,
     proportions_list_b4_notrunc,
     proportions_list_c4_notrunc,
     rmse_values_a_lambda_mean_inv_notrunc,
     rmse_values_b_lambda_mean_inv_notrunc,
     rmse_values_c_lambda_mean_inv_notrunc,
     tab.mean.inv.lambda_notrunc,
     file="mean_inverse_notrunc.RData")


##########################################Inverse Mean Lambda NOT TRUNCATED #############################################

# Create list of matrices to store parameter summary statistics
param_list5_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  inner_list5_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    matrix_data5_notrunc <- matrix(0, nrow = 5, ncol = 4)
    inner_list5_notrunc[[j]] <- matrix_data5_notrunc
  }
  param_list5_notrunc[[i]] <- inner_list5_notrunc
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data.raw <- data.list[[i]]
  data <- data.list[[i]] 
  data$inv.trait <- 1/(data$trait)
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$inv.trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <-  data.frame(
    T = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- data$trait
  data.raw <-  data.list[[i]]
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("inv_mean_lambda_nt.txt")
  cat("model{
  # Priors
  # If mu<0,  make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  tau ~ dexp(0.1)
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c)<= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "inv_mean_lambda_nt.txt")
  # Settings
  parameters <- c("la", "b.l", "c", "tau")
  inits <- function(){list(
    # a and b in log space
    la = log(0.001), 
    b.l = log(0.05), 
    c = 0.6, 
    tau = 0.2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "inv_mean_lambda_nt.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list5_notrunc[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list5_notrunc[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list5_notrunc[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list5_notrunc[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    param_list5_notrunc[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list5_notrunc[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list5_notrunc[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list5_notrunc[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list5_notrunc[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list5_notrunc[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list5_notrunc[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list5_notrunc[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list5_notrunc[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list5_notrunc[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list5_notrunc[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list5_notrunc[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    param_list5_notrunc[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list5_notrunc[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list5_notrunc[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list5_notrunc[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5.NT = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 10, to = 35, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new5.NT as a function at each of those temperatures
xdf5.NT <- expand.grid(point = xseq, iteration = 1:nrow(df.new5.NT))

# Apply the quadratic equation for each combination
xdf5.NT$value <- apply(xdf5.NT, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new5.NT[i, 4]) * location^2 - exp(df.new5.NT[i, 1]) * location + df.new5.NT[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf5.NT <-  xdf5.NT|> mutate(trunc.eval = ifelse(xdf5.NT$value > epsilon, xdf5.NT$value, epsilon), 
                           trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf5.NT <- xdf5.NT |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                      avg.value = mean(trunc.eval), 
                                                      med.value.inv = median(trunc.inv), 
                                                      med.value = median(trunc.eval), 
                                                      lower.hdi.inv = hdi(trunc.inv)[1], 
                                                      upper.hdi.inv = hdi(trunc.inv)[2], 
                                                      lower.hdi = hdi(trunc.eval)[1], 
                                                      upper.hdi = hdi(trunc.eval)[2])
# Creating columns of the true curve values and the true inverted curve values
true_curve.inv <- function(x) 1/(true.a * x^2 - true.b * x + true.c)
true_curve <- function(x) (true.a * x^2 - true.b * x + true.c)
mean.xdf5.NT$true_curve <- true_curve(xseq)
mean.xdf5.NT$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
plot5_mean_NT <- ggplot(mean.xdf5.NT, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot5_med_NT <- ggplot(mean.xdf5.NT, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = "deepskyblue2", linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot5_mean.inv_NT <- ggplot(mean.xdf5.NT, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Response") + 
  theme_minimal()
# Plot median lifetime response
plot5_med.inv_NT <- ggplot(mean.xdf5.NT, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "coral1", alpha = 0.3) + 
  geom_line(aes(y = true_curve.inv), color = "deepskyblue2", linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Response") + 
  theme_minimal()

gridExtra::grid.arrange(plot5_mean_NT, plot5_med_NT, nrow = 1)
gridExtra::grid.arrange(plot5_mean.inv_NT, plot5_med.inv_NT, nrow = 1)

# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
hist_list5_notrunc <- vector("list", length = 5)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  first_column_list5_notrunc <- lapply(param_list5_notrunc[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  combined_vector5_notrunc <- unlist(first_column_list5_notrunc)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(combined_vector5_notrunc)-1, max(combined_vector5_notrunc) + 1, length = 3500)
  # Create a histogram and store it in the list
  hist_list5_notrunc[[i]] <- hist(combined_vector5_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                  xlab = "Values", col = "lightblue", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "red", lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = "green", lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 0, sqrt(10)), col = "black", lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = "yellow", lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 0.1), col = "purple", lty = 2, lwd = 2)
  }
}


proportions_list_a5_notrunc <-  lapply(param_list5_notrunc[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_a5_notrunc))

proportions_list_b5_notrunc <-  lapply(param_list5_notrunc[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_b5_notrunc))

proportions_list_c5_notrunc <-  lapply(param_list5_notrunc[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(proportions_list_c5_notrunc))

calculate_rmse_a_lambda_inv_mean_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_b_lambda_inv_mean_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

calculate_rmse_c_lambda_inv_mean_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
rmse_values_a_lambda_inv_mean_notrunc <- sapply(param_list5_notrunc[[1]], calculate_rmse_a_lambda_inv_mean_notrunc)
rmse_values_b_lambda_inv_mean_notrunc <- sapply(param_list5_notrunc[[2]], calculate_rmse_b_lambda_inv_mean_notrunc)
rmse_values_c_lambda_inv_mean_notrunc <- sapply(param_list5_notrunc[[3]], calculate_rmse_c_lambda_inv_mean_notrunc)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
tab.inv.mean.lambda_notrunc <- matrix(0, nrow = 3, ncol = 6)
row.names(tab.inv.mean.lambda_notrunc) <- c("a", "b", "c")
colnames(tab.inv.mean.lambda_notrunc) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.inv.mean.lambda_notrunc[1, 1] <- true.a
tab.inv.mean.lambda_notrunc[2, 1] <- true.b
tab.inv.mean.lambda_notrunc[3, 1] <- true.c
tab.inv.mean.lambda_notrunc[1, 2] <- exp(mean(unlist(lapply(param_list5_notrunc[[1]], function(matrix) matrix[, 1]))))
tab.inv.mean.lambda_notrunc[2, 2] <- exp(mean(unlist(lapply(param_list5_notrunc[[2]], function(matrix) matrix[, 1]))))
tab.inv.mean.lambda_notrunc[3, 2] <- mean(unlist(lapply(param_list5_notrunc[[3]], function(matrix) matrix[, 1])))
tab.inv.mean.lambda_notrunc[1, 3] <- exp(hdi(unlist(lapply(param_list5_notrunc[[1]], function(matrix) matrix[, 1])))[1])
tab.inv.mean.lambda_notrunc[2, 3] <- exp(hdi(unlist(lapply(param_list5_notrunc[[2]], function(matrix) matrix[, 1])))[1])
tab.inv.mean.lambda_notrunc[3, 3] <- hdi(unlist(lapply(param_list5_notrunc[[3]], function(matrix) matrix[, 1])))[1]
tab.inv.mean.lambda_notrunc[1, 4] <- exp(hdi(unlist(lapply(param_list5_notrunc[[1]], function(matrix) matrix[, 1])))[2])
tab.inv.mean.lambda_notrunc[2, 4] <- exp(hdi(unlist(lapply(param_list5_notrunc[[2]], function(matrix) matrix[, 1])))[2])
tab.inv.mean.lambda_notrunc[3, 4] <- hdi(unlist(lapply(param_list5_notrunc[[3]], function(matrix) matrix[, 1])))[2]
tab.inv.mean.lambda_notrunc[1, 5] <- mean(unlist(proportions_list_a5_notrunc))
tab.inv.mean.lambda_notrunc[2, 5] <- mean(unlist(proportions_list_b5_notrunc))
tab.inv.mean.lambda_notrunc[3, 5] <- mean(unlist(proportions_list_c5_notrunc))
tab.inv.mean.lambda_notrunc[1, 6] <- mean(rmse_values_a_lambda_inv_mean_notrunc)
tab.inv.mean.lambda_notrunc[2, 6] <- mean(rmse_values_b_lambda_inv_mean_notrunc)
tab.inv.mean.lambda_notrunc[3, 6] <- mean(rmse_values_c_lambda_inv_mean_notrunc)


save(param_list5_notrunc,
     df.new5.NT,
     mean.xdf5.NT,
     hist_list5_notrunc,
     plot5_mean_NT, 
     plot5_med_NT,
     plot5_mean.inv_NT, 
     plot5_med.inv_NT,
     proportions_list_a5_notrunc,
     proportions_list_b5_notrunc,
     proportions_list_c5_notrunc,
     rmse_values_a_lambda_inv_mean_notrunc,
     rmse_values_b_lambda_inv_mean_notrunc,
     rmse_values_c_lambda_inv_mean_notrunc,
     tab.inv.mean.lambda_notrunc,
     file="inverse_mean_notrunc.RData")


















