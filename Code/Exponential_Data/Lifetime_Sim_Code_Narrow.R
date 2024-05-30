library(ggplot2)
library(HDInterval)
library(rjags)
library(R2jags)
library(dplyr)

# Set the true parameters from the mortality curve
true.a <- 0.001094
true.b <- 0.0488
true.c <- 0.609

# True quadratic
f <- function(x, a = 0.001094, b = 0.0488, c = 0.609){
  a * x^2 - b * x + c
}

# Generate 100 data sets and store them in a list
# Narrow temperature interval
N_data.list <- list()
set.seed(12)
for(i in 1:100){
  T1 <- rexp(n = 100, f(13))
  T2 <- rexp(n = 100, f(17))
  T3 <- rexp(n = 100, f(21))
  T4 <- rexp(n = 100, f(25))
  T5 <- rexp(n = 100, f(29))
  N_data.list[[i]] <- data.frame(
    T = c(rep(13, 100), rep(17, 100), rep(21, 100), rep(25, 100), rep(29, 100)), 
    trait = c(T1, T2, T3, T4, T5)
  )
}

# Truncation value
# Don't expect lifetime to exceed 100 days
epsilon <- 0.01
b <- length(N_data.list)

# True value vector with a and b in log space
# Spot for sigma and deviance
true_values <- c(log(true.a), log(true.b), true.c, NA, NA)

########################################## Lambda #############################################

# Create list of matrices to store parameter summary statistics
N_param_list1 <- vector("list", length = 4)
for (i in 1:4) {
  # Initialize the inner list
  N_inner_list1 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    N_matrix_data1 <- matrix(0, nrow = 5, ncol = 4)
    N_inner_list1[[j]] <- N_matrix_data1
  }
  N_param_list1[[i]] <- N_inner_list1
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("n_lambda.txt")
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
    ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) <= epsilon) * epsilon
    trait[i] ~ dexp(mu[i])
  }
  }", file = "n_lambda.txt")
  
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
  N_data.raw <- N_data.list[[i]]
  data <- N_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "n_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc fit
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    N_param_list1[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list1[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list1[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list1[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    N_param_list1[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list1[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list1[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list1[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list1[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list1[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list1[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list1[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list1[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list1[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list1[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list1[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
  }
  
  
}

# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new1_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 13, to = 29, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new1_N as a function at each of those temperatures
xdf1_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new1_N))

# Apply the quadratic equation for each combination
xdf1_N$value <- apply(xdf1_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new1_N[i, 4]) * location^2 - exp(df.new1_N[i, 1]) * location + df.new1_N[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf1_N <- xdf1_N|> mutate(trunc.eval = ifelse(xdf1_N$value > epsilon, xdf1_N$value, epsilon), 
                         trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf1 <- xdf1_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
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
N_mean.xdf1$true_curve <- true_curve(xseq)
N_mean.xdf1$true_curve.inv <- true_curve.inv(xseq)
# Plot the mean mortality rate response
N_plot1_mean <- ggplot(N_mean.xdf1, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot the median mortality rate response
N_plot1_med <- ggplot(N_mean.xdf1, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot the mean lifetime response
N_plot1_mean.inv <- ggplot(N_mean.xdf1, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot the median lifetime response
N_plot1_med.inv <- ggplot(N_mean.xdf1, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(N_plot1_mean, N_plot1_med, nrow = 1)
gridExtra::grid.arrange(N_plot1_mean.inv, N_plot1_med.inv, nrow = 1)

# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
true_values <- c(log(0.001094), log(0.0488), 0.609, NA)
N_hist_list1 <- vector("list", length = 4)
par(mfrow = c(2, 2))
for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  N_first_column_list1 <- lapply(N_param_list1[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  N_combined_vector1 <- unlist(N_first_column_list1)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(N_combined_vector1) - 1, max(N_combined_vector1) + 1, length = 1000)
  # Create a histogram and store it in the list
  N_hist_list1[[i]] <- hist(N_combined_vector1, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "lavender", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "forestgreen", lwd = 2)
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


N_proportions_list_a1 <- lapply(N_param_list1[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_a1))

N_proportions_list_b1 <- lapply(N_param_list1[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_b1))

N_proportions_list_c1 <- lapply(N_param_list1[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_c1))

# Create functions to calculate the RMSE for a, b, and c
N_calculate_rmse_a_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_b_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_c_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
N_rmse_values_a_lambda <- sapply(N_param_list1[[1]], N_calculate_rmse_a_lambda)
N_rmse_values_b_lambda <- sapply(N_param_list1[[2]], N_calculate_rmse_b_lambda)
N_rmse_values_c_lambda <- sapply(N_param_list1[[3]], N_calculate_rmse_c_lambda)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
N_tab.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(N_tab.lambda) <- c("a", "b", "c")
colnames(N_tab.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
N_tab.lambda[1, 1] <- true.a
N_tab.lambda[2, 1] <- true.b
N_tab.lambda[3, 1] <- true.c
N_tab.lambda[1, 2] <- exp(mean(unlist(lapply(N_param_list1[[1]], function(matrix) matrix[, 1]))))
N_tab.lambda[2, 2] <- exp(mean(unlist(lapply(N_param_list1[[2]], function(matrix) matrix[, 1]))))
N_tab.lambda[3, 2] <- mean(unlist(lapply(N_param_list1[[3]], function(matrix) matrix[, 1])))
N_tab.lambda[1, 3] <- exp(hdi(unlist(lapply(N_param_list1[[1]], function(matrix) matrix[, 1])))[1])
N_tab.lambda[2, 3] <- exp(hdi(unlist(lapply(N_param_list1[[2]], function(matrix) matrix[, 1])))[1])
N_tab.lambda[3, 3] <- hdi(unlist(lapply(N_param_list1[[3]], function(matrix) matrix[, 1])))[1]
N_tab.lambda[1, 4] <- exp(hdi(unlist(lapply(N_param_list1[[1]], function(matrix) matrix[, 1])))[2])
N_tab.lambda[2, 4] <- exp(hdi(unlist(lapply(N_param_list1[[2]], function(matrix) matrix[, 1])))[2])
N_tab.lambda[3, 4] <- hdi(unlist(lapply(N_param_list1[[3]], function(matrix) matrix[, 1])))[2]
N_tab.lambda[1, 5] <- mean(unlist(N_proportions_list_a1))
N_tab.lambda[2, 5] <- mean(unlist(N_proportions_list_b1))
N_tab.lambda[3, 5] <- mean(unlist(N_proportions_list_c1))
N_tab.lambda[1, 6] <- mean(N_rmse_values_a_lambda)
N_tab.lambda[2, 6] <- mean(N_rmse_values_b_lambda)
N_tab.lambda[3, 6] <- mean(N_rmse_values_c_lambda)

save(N_param_list1,
     xdf1_N,
     df.new1_N,
     N_mean.xdf1,
     N_hist_list1,
     N_plot1_mean, 
     N_plot1_med,
     N_plot1_mean.inv, 
     N_plot1_med.inv,
     N_proportions_list_a1,
     N_proportions_list_b1,
     N_proportions_list_c1,
     N_rmse_values_a_lambda,
     N_rmse_values_b_lambda,
     N_rmse_values_c_lambda,
     N_tab.lambda,
     file="individual_exp_narrow.RData")


########################################## Mean Lambda #############################################

# Create list of matrices to store parameter summary statistics
N_param_list2 <- vector("list", length = 4)
for (i in 1:4) {
  # Initialize the inner list
  N_inner_list2 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    N_matrix_data2 <- matrix(0, nrow = 5, ncol = 4)
    N_inner_list2[[j]] <- N_matrix_data2
  }
  N_param_list2[[i]] <- N_inner_list2
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  N_data.raw <- N_data.list[[i]]
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(N_data.list[[i]]$trait, N_data.list[[i]]$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(13, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(29, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  # Model
  sink("n_mean_lambda.txt")
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
    ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) <= epsilon) * epsilon
    trait[i] ~ dexp(mu[i])
  }
  }", file = "n_mean_lambda.txt")
  
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
                  model.file = "n_mean_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into an mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into an mcmc list
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    N_param_list2[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list2[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list2[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list2[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    N_param_list2[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list2[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list2[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list2[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list2[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list2[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list2[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list2[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list2[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list2[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list2[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list2[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
  }
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new2_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 13, to = 29, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new2_N as a function at each of those temperatures
xdf2_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new2_N))

# Apply the quadratic equation for each combination
xdf2_N$value <- apply(xdf2_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new2_N[i, 4]) * location^2 - exp(df.new2_N[i, 1]) * location + df.new2_N[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf2_N <- xdf2_N|> mutate(trunc.eval = ifelse(xdf2_N$value > epsilon, xdf2_N$value, epsilon), 
                         trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf2 <- xdf2_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
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
N_mean.xdf2$true_curve <- true_curve(xseq)
N_mean.xdf2$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
N_plot2_mean <- ggplot(N_mean.xdf2, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median mortality rate response
N_plot2_med <- ggplot(N_mean.xdf2, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean lifetime response
N_plot2_mean.inv <- ggplot(N_mean.xdf2, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median lifetime response
N_plot2_med.inv <- ggplot(N_mean.xdf2, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(N_plot2_mean, N_plot2_med, nrow = 1)
gridExtra::grid.arrange(N_plot2_mean.inv, N_plot2_med.inv, nrow = 1)

# Create a list to store histograms of the posterior distribution of each parameter
N_hist_list2 <- vector("list", length = 4)
par(mfrow = c(2, 2))
for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  N_first_column_list2 <- lapply(N_param_list2[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  N_combined_vector2 <- unlist(N_first_column_list2)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(N_combined_vector2) - 1, max(N_combined_vector2) + 1, length = 1000)
  # Create a histogram and store it in the list
  N_hist_list2[[i]] <- hist(N_combined_vector2, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "lavender", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "forestgreen", lwd = 2)
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


N_proportions_list_a2 <- lapply(N_param_list2[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_a2))

N_proportions_list_b2 <- lapply(N_param_list2[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_b2))

N_proportions_list_c2 <- lapply(N_param_list2[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_c2))

N_calculate_rmse_a_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_b_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_c_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1] 
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
N_rmse_values_a_lambda_mean <- sapply(N_param_list2[[1]], N_calculate_rmse_a_lambda_mean)
N_rmse_values_b_lambda_mean <- sapply(N_param_list2[[2]], N_calculate_rmse_b_lambda_mean)
N_rmse_values_c_lambda_mean <- sapply(N_param_list2[[3]], N_calculate_rmse_c_lambda_mean)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
N_tab.lambda.mean <- matrix(0, nrow = 3, ncol = 6)
row.names(N_tab.lambda.mean) <- c("a", "b", "c")
colnames(N_tab.lambda.mean) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
N_tab.lambda.mean[1, 1] <- true.a
N_tab.lambda.mean[2, 1] <- true.b
N_tab.lambda.mean[3, 1] <- true.c
N_tab.lambda.mean[1, 2] <- exp(mean(unlist(lapply(N_param_list2[[1]], function(matrix) matrix[, 1]))))
N_tab.lambda.mean[2, 2] <- exp(mean(unlist(lapply(N_param_list2[[2]], function(matrix) matrix[, 1]))))
N_tab.lambda.mean[3, 2] <- mean(unlist(lapply(N_param_list2[[3]], function(matrix) matrix[, 1])))
N_tab.lambda.mean[1, 3] <- exp(hdi(unlist(lapply(N_param_list2[[1]], function(matrix) matrix[, 1])))[1])
N_tab.lambda.mean[2, 3] <- exp(hdi(unlist(lapply(N_param_list2[[2]], function(matrix) matrix[, 1])))[1])
N_tab.lambda.mean[3, 3] <- hdi(unlist(lapply(N_param_list2[[3]], function(matrix) matrix[, 1])))[1]
N_tab.lambda.mean[1, 4] <- exp(hdi(unlist(lapply(N_param_list2[[1]], function(matrix) matrix[, 1])))[2])
N_tab.lambda.mean[2, 4] <- exp(hdi(unlist(lapply(N_param_list2[[2]], function(matrix) matrix[, 1])))[2])
N_tab.lambda.mean[3, 4] <- hdi(unlist(lapply(N_param_list2[[3]], function(matrix) matrix[, 1])))[2]
N_tab.lambda.mean[1, 5] <- mean(unlist(N_proportions_list_a2))
N_tab.lambda.mean[2, 5] <- mean(unlist(N_proportions_list_b2))
N_tab.lambda.mean[3, 5] <- mean(unlist(N_proportions_list_c2))
N_tab.lambda.mean[1, 6] <- mean(N_rmse_values_a_lambda_mean)
N_tab.lambda.mean[2, 6] <- mean(N_rmse_values_b_lambda_mean)
N_tab.lambda.mean[3, 6] <- mean(N_rmse_values_c_lambda_mean)

save(N_param_list2,
     xdf2_N,
     df.new2_N,
     N_mean.xdf2,
     N_hist_list2,
     N_plot2_mean, 
     N_plot2_med,
     N_plot2_mean.inv, 
     N_plot2_med.inv,
     N_proportions_list_a2,
     N_proportions_list_b2,
     N_proportions_list_c2,
     N_rmse_values_a_lambda_mean,
     N_rmse_values_b_lambda_mean,
     N_rmse_values_c_lambda_mean,
     N_tab.lambda.mean,
     file="mean_exp_narrow.RData")

########################################## Inverse Lambda #############################################

# Create list of matrices to store parameter summary statistics
N_param_list3 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  N_inner_list3 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    N_matrix_data3 <- matrix(0, nrow = 5, ncol = 4)
    N_inner_list3[[j]] <- N_matrix_data3
  }
  N_param_list3[[i]] <- N_inner_list3
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("n_inv_lambda.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(15000)
  sig2 <- sig^2
  tau <- 1/sig2
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + 
    ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) <= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau) T(0, )
  }
}", file = "n_inv_lambda.txt")
  
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
  N_data.raw <- N_data.list[[i]]
  
  data <- N_data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "n_inv_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, DIC = T, working.directory = getwd())
  # Turn into mcmc object
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  # Turn into mcmc fit
  samp <- as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    N_param_list3[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list3[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list3[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list3[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    N_param_list3[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list3[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list3[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list3[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list3[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list3[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list3[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list3[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list3[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list3[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list3[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list3[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    N_param_list3[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list3[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list3[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list3[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new3_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 13, to = 29, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new3_N as a function at each of those temperatures
xdf3_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new3_N))

# Apply the quadratic equation for each combination
xdf3_N$value <- apply(xdf3_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new3_N[i, 4]) * location^2 - exp(df.new3_N[i, 1]) * location + df.new3_N[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf3_N <- xdf3_N|> mutate(trunc.eval = ifelse(xdf3_N$value > epsilon, xdf3_N$value, epsilon), 
                         trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf3 <- xdf3_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
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
N_mean.xdf3$true_curve <- true_curve(xseq)
N_mean.xdf3$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
N_plot3_mean <- ggplot(N_mean.xdf3, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median mortality rate response
N_plot3_med <- ggplot(N_mean.xdf3, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean lifetime response
N_plot3_mean.inv <- ggplot(N_mean.xdf3, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median lifetime response
N_plot3_med.inv <- ggplot(N_mean.xdf3, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(N_plot3_mean, N_plot3_med, nrow = 1)
gridExtra::grid.arrange(N_plot3_mean.inv, N_plot3_med.inv, nrow = 1)

# Create a list to store histograms of the posterior distribution of each parameter
N_hist_list3 <- vector("list", length = 5)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  N_first_column_list3 <- lapply(N_param_list3[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  N_combined_vector3 <- unlist(N_first_column_list3)
  x <- seq(min(N_combined_vector3) - 1, max(N_combined_vector3) + 1, length = 1000)
  # Create a histogram and store it in the list
  N_hist_list3[[i]] <- hist(N_combined_vector3, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "lavender", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "forestgreen", lwd = 2)
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
    lines(x, dexp(x, 5000), col = "purple", lty = 2, lwd = 2)
  }
}


N_proportions_list_a3 <- lapply(N_param_list3[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_a3))

N_proportions_list_b3 <- lapply(N_param_list3[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_b3))

N_proportions_list_c3 <- lapply(N_param_list3[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_c3))

N_calculate_rmse_a_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_b_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_c_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
N_rmse_values_a_lambda_inv <- sapply(N_param_list3[[1]], N_calculate_rmse_a_lambda_inv)
N_rmse_values_b_lambda_inv <- sapply(N_param_list3[[2]], N_calculate_rmse_b_lambda_inv)
N_rmse_values_c_lambda_inv <- sapply(N_param_list3[[3]], N_calculate_rmse_c_lambda_inv)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
N_tab.inv.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(N_tab.inv.lambda) <- c("a", "b", "c")
colnames(N_tab.inv.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
N_tab.inv.lambda[1, 1] <- true.a
N_tab.inv.lambda[2, 1] <- true.b
N_tab.inv.lambda[3, 1] <- true.c
N_tab.inv.lambda[1, 2] <- exp(mean(unlist(lapply(N_param_list3[[1]], function(matrix) matrix[, 1]))))
N_tab.inv.lambda[2, 2] <- exp(mean(unlist(lapply(N_param_list3[[2]], function(matrix) matrix[, 1]))))
N_tab.inv.lambda[3, 2] <- mean(unlist(lapply(N_param_list3[[3]], function(matrix) matrix[, 1])))
N_tab.inv.lambda[1, 3] <- exp(hdi(unlist(lapply(N_param_list3[[1]], function(matrix) matrix[, 1])))[1])
N_tab.inv.lambda[2, 3] <- exp(hdi(unlist(lapply(N_param_list3[[2]], function(matrix) matrix[, 1])))[1])
N_tab.inv.lambda[3, 3] <- hdi(unlist(lapply(N_param_list3[[3]], function(matrix) matrix[, 1])))[1]
N_tab.inv.lambda[1, 4] <- exp(hdi(unlist(lapply(N_param_list3[[1]], function(matrix) matrix[, 1])))[2])
N_tab.inv.lambda[2, 4] <- exp(hdi(unlist(lapply(N_param_list3[[2]], function(matrix) matrix[, 1])))[2])
N_tab.inv.lambda[3, 4] <- hdi(unlist(lapply(N_param_list3[[3]], function(matrix) matrix[, 1])))[2]
N_tab.inv.lambda[1, 5] <- mean(unlist(N_proportions_list_a3))
N_tab.inv.lambda[2, 5] <- mean(unlist(N_proportions_list_b3))
N_tab.inv.lambda[3, 5] <- mean(unlist(N_proportions_list_c3))
N_tab.inv.lambda[1, 6] <- mean(N_rmse_values_a_lambda_inv)
N_tab.inv.lambda[2, 6] <- mean(N_rmse_values_b_lambda_inv)
N_tab.inv.lambda[3, 6] <- mean(N_rmse_values_c_lambda_inv)

save(N_param_list3,
     xdf3_N,
     df.new3_N,
     N_mean.xdf3,
     N_hist_list3,
     N_plot3_mean, 
     N_plot3_med,
     N_plot3_mean.inv, 
     N_plot3_med.inv,
     N_proportions_list_a3,
     N_proportions_list_b3,
     N_proportions_list_c3,
     N_rmse_values_a_lambda_inv,
     N_rmse_values_b_lambda_inv,
     N_rmse_values_c_lambda_inv,
     N_tab.inv.lambda,
     file="inverse_exp_narrow.RData")

########################################## Mean Inverse Lambda #############################################

# Create list of matrices to store parameter summary statistics
N_param_list4 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  N_inner_list4 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    N_matrix_data4 <- matrix(0, nrow = 5, ncol = 4)
    N_inner_list4[[j]] <- N_matrix_data4
  }
  N_param_list4[[i]] <- N_inner_list4
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  N_data.raw <- N_data.list[[i]]
  data <- N_data.list[[i]] 
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(13, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(29, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("n_mean_inv_lambda.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(1000)
  sig2 <- sig^2
  tau <- 1/sig2
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + 
    ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) <= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau) T(0, )
  }
}", file = "n_mean_inv_lambda.txt")
  
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
                  model.file = "n_mean_inv_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    N_param_list4[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list4[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list4[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list4[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    N_param_list4[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list4[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list4[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list4[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list4[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list4[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list4[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list4[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list4[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list4[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list4[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list4[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    N_param_list4[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list4[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list4[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list4[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
  
}

# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 13, to = 29, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new4_N as a function at each of those temperatures
xdf4_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new4_N))

# Apply the quadratic equation for each combination
xdf4_N$value <- apply(xdf4_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new4_N[i, 4]) * location^2 - exp(df.new4_N[i, 1]) * location + df.new4_N[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf4_N <- xdf4_N|> mutate(trunc.eval = ifelse(xdf4_N$value > epsilon, xdf4_N$value, epsilon), 
                         trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf4 <- xdf4_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
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
N_mean.xdf4$true_curve <- true_curve(xseq)
N_mean.xdf4$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
N_plot4_mean <- ggplot(N_mean.xdf4, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median mortality rate response
N_plot4_med <- ggplot(N_mean.xdf4, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean lifetime response
N_plot4_mean.inv <- ggplot(N_mean.xdf4, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median lifetime response
N_plot4_med.inv <- ggplot(N_mean.xdf4, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
gridExtra::grid.arrange(N_plot4_mean, N_plot4_med, nrow = 1)
gridExtra::grid.arrange(N_plot4_mean.inv, N_plot4_med.inv, nrow = 1)

# Create a list to store histograms of the posterior distribution of each parameter
N_hist_list4 <- vector("list", length = 4)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  N_first_column_list4 <- lapply(N_param_list4[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  N_combined_vector4 <- unlist(N_first_column_list4)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(N_combined_vector4) - 1, max(N_combined_vector4) + 1, length = 1000)
  # Create a histogram and store it in the list
  N_hist_list4[[i]] <- hist(N_combined_vector4, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "lavender", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "forestgreen", lwd = 2)
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
    lines(x, dexp(x, 1000), col = "purple", lty = 2, lwd = 2)
  }
}


N_proportions_list_a4 <- lapply(N_param_list4[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_a4))

N_proportions_list_b4 <- lapply(N_param_list4[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_b4))

N_proportions_list_c4 <- lapply(N_param_list4[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_c4))

N_calculate_rmse_a_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1] 
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_b_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_c_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
N_rmse_values_a_lambda_mean_inv <- sapply(N_param_list4[[1]], N_calculate_rmse_a_lambda_mean_inv)
N_rmse_values_b_lambda_mean_inv <- sapply(N_param_list4[[2]], N_calculate_rmse_b_lambda_mean_inv)
N_rmse_values_c_lambda_mean_inv <- sapply(N_param_list4[[3]], N_calculate_rmse_c_lambda_mean_inv)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
N_tab.mean.inv.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(N_tab.mean.inv.lambda) <- c("a", "b", "c")
colnames(N_tab.mean.inv.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
N_tab.mean.inv.lambda[1, 1] <- true.a
N_tab.mean.inv.lambda[2, 1] <- true.b
N_tab.mean.inv.lambda[3, 1] <- true.c
N_tab.mean.inv.lambda[1, 2] <- exp(mean(unlist(lapply(N_param_list4[[1]], function(matrix) matrix[, 1]))))
N_tab.mean.inv.lambda[2, 2] <- exp(mean(unlist(lapply(N_param_list4[[2]], function(matrix) matrix[, 1]))))
N_tab.mean.inv.lambda[3, 2] <- mean(unlist(lapply(N_param_list4[[3]], function(matrix) matrix[, 1])))
N_tab.mean.inv.lambda[1, 3] <- exp(hdi(unlist(lapply(N_param_list4[[1]], function(matrix) matrix[, 1])))[1])
N_tab.mean.inv.lambda[2, 3] <- exp(hdi(unlist(lapply(N_param_list4[[2]], function(matrix) matrix[, 1])))[1])
N_tab.mean.inv.lambda[3, 3] <- hdi(unlist(lapply(N_param_list4[[3]], function(matrix) matrix[, 1])))[1]
N_tab.mean.inv.lambda[1, 4] <- exp(hdi(unlist(lapply(N_param_list4[[1]], function(matrix) matrix[, 1])))[2])
N_tab.mean.inv.lambda[2, 4] <- exp(hdi(unlist(lapply(N_param_list4[[2]], function(matrix) matrix[, 1])))[2])
N_tab.mean.inv.lambda[3, 4] <- hdi(unlist(lapply(N_param_list4[[3]], function(matrix) matrix[, 1])))[2]
N_tab.mean.inv.lambda[1, 5] <- mean(unlist(N_proportions_list_a4))
N_tab.mean.inv.lambda[2, 5] <- mean(unlist(N_proportions_list_b4))
N_tab.mean.inv.lambda[3, 5] <- mean(unlist(N_proportions_list_c4))
N_tab.mean.inv.lambda[1, 6] <- mean(N_rmse_values_a_lambda_mean_inv)
N_tab.mean.inv.lambda[2, 6] <- mean(N_rmse_values_b_lambda_mean_inv)
N_tab.mean.inv.lambda[3, 6] <- mean(N_rmse_values_c_lambda_mean_inv)

save(N_param_list4,
     xdf4_N,
     df.new4_N,
     N_mean.xdf4,
     N_hist_list4,
     N_plot4_mean, 
     N_plot4_med,
     N_plot4_mean.inv, 
     N_plot4_med.inv,
     N_proportions_list_a4,
     N_proportions_list_b4,
     N_proportions_list_c4,
     N_rmse_values_a_lambda_mean_inv,
     N_rmse_values_b_lambda_mean_inv,
     N_rmse_values_c_lambda_mean_inv,
     N_tab.mean.inv.lambda,
     file="mean_inverse_exp_narrow.RData")

########################################## Inverse Mean Lambda #############################################

# Create list of matrices to store parameter summary statistics
N_param_list5 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  N_inner_list5 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    N_matrix_data5 <- matrix(0, nrow = 5, ncol = 4)
    N_inner_list5[[j]] <- N_matrix_data5
  }
  N_param_list5[[i]] <- N_inner_list5
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data.raw <- N_data.list[[i]]
  data <- N_data.list[[i]] 
  data$inv.trait <- 1/(data$trait)
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$inv.trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(13, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(29, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- data$trait
  
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("n_inv_mean_lambda.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(7000)
  sig2 <- sig^2
  tau <- 1/sig2
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + 
    ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) <= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau) T(0, )
  }
}", file = "n_inv_mean_lambda.txt")
  
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
  
  N_data.raw <- N_data.list[[i]]
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "n_inv_mean_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    N_param_list5[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list5[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list5[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list5[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    N_param_list5[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list5[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list5[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list5[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list5[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list5[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list5[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list5[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list5[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list5[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list5[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list5[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    N_param_list5[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list5[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list5[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list5[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 13, to = 29, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new5_N as a function at each of those temperatures
xdf5_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new5_N))

# Apply the quadratic equation for each combination
xdf5_N$value <- apply(xdf5_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new5_N[i, 4]) * location^2 - exp(df.new5_N[i, 1]) * location + df.new5_N[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf5_N <- xdf5_N|> mutate(trunc.eval = ifelse(xdf5_N$value > epsilon, xdf5_N$value, epsilon), 
                         trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf5 <- xdf5_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
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
N_mean.xdf5$true_curve <- true_curve(xseq)
N_mean.xdf5$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
N_plot5_mean <- ggplot(N_mean.xdf5, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median mortality rate response
N_plot5_med <- ggplot(N_mean.xdf5, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean lifetime response
N_plot5_mean.inv <- ggplot(N_mean.xdf5, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median lifetime response
N_plot5_med.inv <- ggplot(N_mean.xdf5, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(N_plot5_mean, N_plot5_med, nrow = 1)
gridExtra::grid.arrange(N_plot5_mean.inv, N_plot5_med.inv, nrow = 1)

# Create a list to store histograms of the posterior distribution of each parameter
N_hist_list5 <- vector("list", length = 5)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  N_first_column_list5 <- lapply(N_param_list5[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  N_combined_vector5 <- unlist(N_first_column_list5)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(N_combined_vector5) - 1, max(N_combined_vector5) + 1, length = 1000)
  # Create a histogram and store it in the list
  N_hist_list5[[i]] <- hist(N_combined_vector5, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "lavender", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "forestgreen", lwd = 2)
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
    lines(x, dexp(x, 7000), col = "purple", lty = 2, lwd = 2)
  }
}

N_proportions_list_a5 <- lapply(N_param_list5[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_a5))

N_proportions_list_b5 <- lapply(N_param_list5[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_b5))

N_proportions_list_c5 <- lapply(N_param_list5[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_c5))

N_calculate_rmse_a_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_b_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_c_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
N_rmse_values_a_lambda_inv_mean <- sapply(N_param_list5[[1]], N_calculate_rmse_a_lambda_inv_mean)
N_rmse_values_b_lambda_inv_mean <- sapply(N_param_list5[[2]], N_calculate_rmse_b_lambda_inv_mean)
N_rmse_values_c_lambda_inv_mean <- sapply(N_param_list5[[3]], N_calculate_rmse_c_lambda_inv_mean)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
N_tab.inv.mean.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(N_tab.inv.mean.lambda) <- c("a", "b", "c")
colnames(N_tab.inv.mean.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
N_tab.inv.mean.lambda[1, 1] <- true.a
N_tab.inv.mean.lambda[2, 1] <- true.b
N_tab.inv.mean.lambda[3, 1] <- true.c
N_tab.inv.mean.lambda[1, 2] <- exp(mean(unlist(lapply(N_param_list5[[1]], function(matrix) matrix[, 1]))))
N_tab.inv.mean.lambda[2, 2] <- exp(mean(unlist(lapply(N_param_list5[[2]], function(matrix) matrix[, 1]))))
N_tab.inv.mean.lambda[3, 2] <- mean(unlist(lapply(N_param_list5[[3]], function(matrix) matrix[, 1])))
N_tab.inv.mean.lambda[1, 3] <- exp(hdi(unlist(lapply(N_param_list5[[1]], function(matrix) matrix[, 1])))[1])
N_tab.inv.mean.lambda[2, 3] <- exp(hdi(unlist(lapply(N_param_list5[[2]], function(matrix) matrix[, 1])))[1])
N_tab.inv.mean.lambda[3, 3] <- hdi(unlist(lapply(N_param_list5[[3]], function(matrix) matrix[, 1])))[1]
N_tab.inv.mean.lambda[1, 4] <- exp(hdi(unlist(lapply(N_param_list5[[1]], function(matrix) matrix[, 1])))[2])
N_tab.inv.mean.lambda[2, 4] <- exp(hdi(unlist(lapply(N_param_list5[[2]], function(matrix) matrix[, 1])))[2])
N_tab.inv.mean.lambda[3, 4] <- hdi(unlist(lapply(N_param_list5[[3]], function(matrix) matrix[, 1])))[2]
N_tab.inv.mean.lambda[1, 5] <- mean(unlist(N_proportions_list_a5))
N_tab.inv.mean.lambda[2, 5] <- mean(unlist(N_proportions_list_b5))
N_tab.inv.mean.lambda[3, 5] <- mean(unlist(N_proportions_list_c5))
N_tab.inv.mean.lambda[1, 6] <- mean(N_rmse_values_a_lambda_inv_mean)
N_tab.inv.mean.lambda[2, 6] <- mean(N_rmse_values_b_lambda_inv_mean)
N_tab.inv.mean.lambda[3, 6] <- mean(N_rmse_values_c_lambda_inv_mean)

save(N_param_list5,
     xdf5_N,
     df.new5_N,
     N_mean.xdf5,
     N_hist_list5,
     N_plot5_mean, 
     N_plot5_med,
     N_plot5_mean.inv, 
     N_plot5_med.inv,
     N_proportions_list_a5,
     N_proportions_list_b5,
     N_proportions_list_c5,
     N_rmse_values_a_lambda_inv_mean,
     N_rmse_values_b_lambda_inv_mean,
     N_rmse_values_c_lambda_inv_mean,
     N_tab.inv.mean.lambda,
     file="inverse_mean_exp_narrow.RData")





########################################## Inverse Lambda NOT TRUNCATED #############################################

# Create list of matrices to store parameter summary statistics
N_param_list3_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  N_inner_list3_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    N_matrix_data3_notrunc <- matrix(0, nrow = 5, ncol = 4)
    N_inner_list3_notrunc[[j]] <- N_matrix_data3_notrunc
  }
  N_param_list3_notrunc[[i]] <- N_inner_list3_notrunc
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("n_inv_lambda_nt.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10)
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(7000)
  sig2 <- sig^2
  tau <- 1/sig2
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + 
    ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) <= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "n_inv_lambda_nt.txt")
  
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
  data <- N_data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  N_data.raw <- N_data.list[[i]]
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "n_inv_lambda_nt.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    N_param_list3_notrunc[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list3_notrunc[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list3_notrunc[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list3_notrunc[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    N_param_list3_notrunc[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list3_notrunc[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list3_notrunc[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list3_notrunc[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list3_notrunc[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list3_notrunc[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list3_notrunc[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list3_notrunc[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list3_notrunc[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list3_notrunc[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list3_notrunc[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list3_notrunc[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    N_param_list3_notrunc[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list3_notrunc[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list3_notrunc[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list3_notrunc[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new3_N.NT = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 13, to = 29, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new3_N.NT as a function at each of those temperatures
xdf3_N.NT <- expand.grid(point = xseq, iteration = 1:nrow(df.new3_N.NT))

# Apply the quadratic equation for each combination
xdf3_N.NT$value <- apply(xdf3_N.NT, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new3_N.NT[i, 4]) * location^2 - exp(df.new3_N.NT[i, 1]) * location + df.new3_N.NT[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf3_N.NT <- xdf3_N.NT|> mutate(trunc.eval = ifelse(xdf3_N.NT$value > epsilon, xdf3_N.NT$value, epsilon), 
                               trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf3_NT <- xdf3_N.NT |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
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
N_mean.xdf3_NT$true_curve <- true_curve(xseq)
N_mean.xdf3_NT$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
N_plot3_mean_NT <- ggplot(N_mean.xdf3_NT, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median mortality rate response
N_plot3_med_NT <- ggplot(N_mean.xdf3_NT, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean lifetime response
N_plot3_mean.inv_NT <- ggplot(N_mean.xdf3_NT, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median lifetime response
N_plot3_med.inv_NT <- ggplot(N_mean.xdf3_NT, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(N_plot3_mean_NT, N_plot3_med_NT, nrow = 1)
gridExtra::grid.arrange(N_plot3_mean.inv_NT, N_plot3_med.inv_NT, nrow = 1)

# Create a list to store histograms of the posterior distribution of each parameter
N_hist_list3_notrunc <- vector("list", length = 5)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  N_first_column_list3_notrunc <- lapply(N_param_list3_notrunc[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  N_combined_vector3_notrunc <- unlist(N_first_column_list3_notrunc)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(N_combined_vector3_notrunc) - 1, max(N_combined_vector3_notrunc) + 1, length = 1000)
  # Create a histogram and store it in the list
  N_hist_list3_notrunc[[i]] <- hist(N_combined_vector3_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                  xlab = "Values", col = "lavender", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "forestgreen", lwd = 2)
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
    lines(x, dexp(x, 7000), col = "purple", lty = 2, lwd = 2)
  }
}


N_proportions_list_a3_notrunc <- lapply(N_param_list3_notrunc[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_a3_notrunc))

N_proportions_list_b3_notrunc <- lapply(N_param_list3_notrunc[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_b3_notrunc))

N_proportions_list_c3_notrunc <- lapply(N_param_list3_notrunc[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_c3_notrunc))

N_calculate_rmse_a_lambda_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_b_lambda_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_c_lambda_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
N_rmse_values_a_lambda_inv_notrunc <- sapply(N_param_list3_notrunc[[1]], N_calculate_rmse_a_lambda_inv_notrunc)
N_rmse_values_b_lambda_inv_notrunc <- sapply(N_param_list3_notrunc[[2]], N_calculate_rmse_b_lambda_inv_notrunc)
N_rmse_values_c_lambda_inv_notrunc <- sapply(N_param_list3_notrunc[[3]], N_calculate_rmse_c_lambda_inv_notrunc)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
N_tab.inv.lambda_notrunc <- matrix(0, nrow = 3, ncol = 6)
row.names(N_tab.inv.lambda_notrunc) <- c("a", "b", "c")
colnames(N_tab.inv.lambda_notrunc) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
N_tab.inv.lambda_notrunc[1, 1] <- true.a
N_tab.inv.lambda_notrunc[2, 1] <- true.b
N_tab.inv.lambda_notrunc[3, 1] <- true.c
N_tab.inv.lambda_notrunc[1, 2] <- exp(mean(unlist(lapply(N_param_list3_notrunc[[1]], function(matrix) matrix[, 1]))))
N_tab.inv.lambda_notrunc[2, 2] <- exp(mean(unlist(lapply(N_param_list3_notrunc[[2]], function(matrix) matrix[, 1]))))
N_tab.inv.lambda_notrunc[3, 2] <- mean(unlist(lapply(N_param_list3_notrunc[[3]], function(matrix) matrix[, 1])))
N_tab.inv.lambda_notrunc[1, 3] <- exp(hdi(unlist(lapply(N_param_list3_notrunc[[1]], function(matrix) matrix[, 1])))[1])
N_tab.inv.lambda_notrunc[2, 3] <- exp(hdi(unlist(lapply(N_param_list3_notrunc[[2]], function(matrix) matrix[, 1])))[1])
N_tab.inv.lambda_notrunc[3, 3] <- hdi(unlist(lapply(N_param_list3_notrunc[[3]], function(matrix) matrix[, 1])))[1]
N_tab.inv.lambda_notrunc[1, 4] <- exp(hdi(unlist(lapply(N_param_list3_notrunc[[1]], function(matrix) matrix[, 1])))[2])
N_tab.inv.lambda_notrunc[2, 4] <- exp(hdi(unlist(lapply(N_param_list3_notrunc[[2]], function(matrix) matrix[, 1])))[2])
N_tab.inv.lambda_notrunc[3, 4] <- hdi(unlist(lapply(N_param_list3_notrunc[[3]], function(matrix) matrix[, 1])))[2]
N_tab.inv.lambda_notrunc[1, 5] <- mean(unlist(N_proportions_list_a3_notrunc))
N_tab.inv.lambda_notrunc[2, 5] <- mean(unlist(N_proportions_list_b3_notrunc))
N_tab.inv.lambda_notrunc[3, 5] <- mean(unlist(N_proportions_list_c3_notrunc))
N_tab.inv.lambda_notrunc[1, 6] <- mean(N_rmse_values_a_lambda_inv_notrunc)
N_tab.inv.lambda_notrunc[2, 6] <- mean(N_rmse_values_b_lambda_inv_notrunc)
N_tab.inv.lambda_notrunc[3, 6] <- mean(N_rmse_values_c_lambda_inv_notrunc)


save(N_param_list3_notrunc,
     xdf3_N.NT,
     df.new3_N.NT,
     N_mean.xdf3_NT,
     N_hist_list3_notrunc,
     N_plot3_mean_NT, 
     N_plot3_med_NT,
     N_plot3_mean_NT, 
     N_plot3_med.inv_NT,
     N_proportions_list_a3_notrunc,
     N_proportions_list_b3_notrunc,
     N_proportions_list_c3_notrunc,
     N_rmse_values_a_lambda_inv_notrunc,
     N_rmse_values_b_lambda_inv_notrunc,
     N_rmse_values_c_lambda_inv_notrunc,
     N_tab.inv.lambda_notrunc,
     file="inverse_exp_narrow_NT.RData")

########################################## Mean Inverse Lambda NOT TRUNCATED #############################################

# Create list of matrices to store parameter summary statistics
N_param_list4_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  N_inner_list4_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    N_matrix_data4_notrunc <- matrix(0, nrow = 5, ncol = 4)
    N_inner_list4_notrunc[[j]] <- N_matrix_data4_notrunc
  }
  N_param_list4_notrunc[[i]] <- N_inner_list4_notrunc
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- N_data.list[[i]] 
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(13, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(29, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  N_data.raw <- N_data.list[[i]]
  
  data <- df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("n_mean_inv_lambda_nt.txt")
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
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + 
    ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) <= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "n_mean_inv_lambda_nt.txt")
  
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
                  model.file = "n_mean_inv_lambda_nt.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    N_param_list4_notrunc[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list4_notrunc[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list4_notrunc[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list4_notrunc[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    N_param_list4_notrunc[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list4_notrunc[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list4_notrunc[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list4_notrunc[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list4_notrunc[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list4_notrunc[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list4_notrunc[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list4_notrunc[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list4_notrunc[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list4_notrunc[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list4_notrunc[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list4_notrunc[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    N_param_list4_notrunc[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list4_notrunc[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list4_notrunc[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list4_notrunc[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4_N.NT = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 13, to = 29, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new4_N.NT as a function at each of those temperatures
xdf4_N.NT <- expand.grid(point = xseq, iteration = 1:nrow(df.new4_N.NT))

# Apply the quadratic equation for each combination
xdf4_N.NT$value <- apply(xdf4_N.NT, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new4_N.NT[i, 4]) * location^2 - exp(df.new4_N.NT[i, 1]) * location + df.new4_N.NT[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf4_N.NT <- xdf4_N.NT|> mutate(trunc.eval = ifelse(xdf4_N.NT$value > epsilon, xdf4_N.NT$value, epsilon), 
                               trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf4_NT <- xdf4_N.NT |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
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
N_mean.xdf4_NT$true_curve <- true_curve(xseq)
N_mean.xdf4_NT$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
N_plot4_mean_NT <- ggplot(N_mean.xdf4_NT, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median mortality rate response
N_plot4_med_NT <- ggplot(N_mean.xdf4_NT, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean lifetime response
N_plot4_mean.inv_NT <- ggplot(N_mean.xdf4_NT, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median lifetime response
N_plot4_med.inv_NT <- ggplot(N_mean.xdf4_NT, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(N_plot4_mean_NT, N_plot4_med_NT, nrow = 1)
gridExtra::grid.arrange(N_plot4_mean.inv_NT, N_plot4_med.inv_NT, nrow = 1)


# Create a list to store histograms of the posterior distribution of each parameter
N_hist_list4_notrunc <- vector("list", length = 4)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  N_first_column_list4_notrunc <- lapply(N_param_list4_notrunc[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  N_combined_vector4_notrunc <- unlist(N_first_column_list4_notrunc)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(N_combined_vector4_notrunc)-1, max(N_combined_vector4_notrunc) +1, length = 1000)
  # Create a histogram and store it in the list
  N_hist_list4_notrunc[[i]] <- hist(N_combined_vector4_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                    xlab = "Values", col = "lavender", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "forestgreen", lwd = 2)
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


N_proportions_list_a4_notrunc <- lapply(N_param_list4_notrunc[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_a4_notrunc))

N_proportions_list_b4_notrunc <- lapply(N_param_list4_notrunc[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_b4_notrunc))

N_proportions_list_c4_notrunc <- lapply(N_param_list4_notrunc[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_c4_notrunc))

N_calculate_rmse_a_lambda_mean_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_b_lambda_mean_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_c_lambda_mean_inv_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1] 
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
N_rmse_values_a_lambda_mean_inv_notrunc <- sapply(N_param_list4_notrunc[[1]], N_calculate_rmse_a_lambda_mean_inv_notrunc)
N_rmse_values_b_lambda_mean_inv_notrunc <- sapply(N_param_list4_notrunc[[2]], N_calculate_rmse_b_lambda_mean_inv_notrunc)
N_rmse_values_c_lambda_mean_inv_notrunc <- sapply(N_param_list4_notrunc[[3]], N_calculate_rmse_c_lambda_mean_inv_notrunc)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
N_tab.mean.inv.lambda_notrunc <- matrix(0, nrow = 3, ncol = 6)
row.names(N_tab.mean.inv.lambda_notrunc) <- c("a", "b", "c")
colnames(N_tab.mean.inv.lambda_notrunc) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
N_tab.mean.inv.lambda_notrunc[1, 1] <- true.a
N_tab.mean.inv.lambda_notrunc[2, 1] <- true.b
N_tab.mean.inv.lambda_notrunc[3, 1] <- true.c
N_tab.mean.inv.lambda_notrunc[1, 2] <- exp(mean(unlist(lapply(N_param_list4_notrunc[[1]], function(matrix) matrix[, 1]))))
N_tab.mean.inv.lambda_notrunc[2, 2] <- exp(mean(unlist(lapply(N_param_list4_notrunc[[2]], function(matrix) matrix[, 1]))))
N_tab.mean.inv.lambda_notrunc[3, 2] <- mean(unlist(lapply(N_param_list4_notrunc[[3]], function(matrix) matrix[, 1])))
N_tab.mean.inv.lambda_notrunc[1, 3] <- exp(hdi(unlist(lapply(N_param_list4_notrunc[[1]], function(matrix) matrix[, 1])))[1])
N_tab.mean.inv.lambda_notrunc[2, 3] <- exp(hdi(unlist(lapply(N_param_list4_notrunc[[2]], function(matrix) matrix[, 1])))[1])
N_tab.mean.inv.lambda_notrunc[3, 3] <- hdi(unlist(lapply(N_param_list4_notrunc[[3]], function(matrix) matrix[, 1])))[1]
N_tab.mean.inv.lambda_notrunc[1, 4] <- exp(hdi(unlist(lapply(N_param_list4_notrunc[[1]], function(matrix) matrix[, 1])))[2])
N_tab.mean.inv.lambda_notrunc[2, 4] <- exp(hdi(unlist(lapply(N_param_list4_notrunc[[2]], function(matrix) matrix[, 1])))[2])
N_tab.mean.inv.lambda_notrunc[3, 4] <- hdi(unlist(lapply(N_param_list4_notrunc[[3]], function(matrix) matrix[, 1])))[2]
N_tab.mean.inv.lambda_notrunc[1, 5] <- mean(unlist(N_proportions_list_a4_notrunc))
N_tab.mean.inv.lambda_notrunc[2, 5] <- mean(unlist(N_proportions_list_b4_notrunc))
N_tab.mean.inv.lambda_notrunc[3, 5] <- mean(unlist(N_proportions_list_c4_notrunc))
N_tab.mean.inv.lambda_notrunc[1, 6] <- mean(N_rmse_values_a_lambda_mean_inv_notrunc)
N_tab.mean.inv.lambda_notrunc[2, 6] <- mean(N_rmse_values_b_lambda_mean_inv_notrunc)
N_tab.mean.inv.lambda_notrunc[3, 6] <- mean(N_rmse_values_c_lambda_mean_inv_notrunc)


save(N_param_list4_notrunc,
     xdf4_N.NT,
     df.new4_N.NT,
     N_mean.xdf4_NT,
     N_hist_list4_notrunc,
     N_plot4_mean_NT, 
     N_plot4_med_NT,
     N_plot4_mean_NT, 
     N_plot4_med.inv_NT,
     N_proportions_list_a4_notrunc,
     N_proportions_list_b4_notrunc,
     N_proportions_list_c4_notrunc,
     N_rmse_values_a_lambda_mean_inv_notrunc,
     N_rmse_values_b_lambda_mean_inv_notrunc,
     N_rmse_values_c_lambda_mean_inv_notrunc,
     N_tab.mean.inv.lambda_notrunc,
     file="mean_inverse_exp_narrow_NT.RData")




########################################## Inverse Mean Lambda NOT TRUNCATED #############################################

# Create list of matrices to store parameter summary statistics
N_param_list5_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  N_inner_list5_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    N_matrix_data5_notrunc <- matrix(0, nrow = 5, ncol = 4)
    N_inner_list5_notrunc[[j]] <- N_matrix_data5_notrunc
  }
  N_param_list5_notrunc[[i]] <- N_inner_list5_notrunc
}
# Create JAGS model for each of the 100 datasets
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  N_data.raw <- N_data.list[[i]]
  data <- N_data.list[[i]] 
  data$inv.trait <- 1/(data$trait)
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$inv.trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(13, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(29, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- data$trait
  
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("n_inv_mean_lambda_nt.txt")
  cat("model{
  # Priors
  # If mu<0, make it small number
  la ~ dnorm(0, 1/10) 
  b.l ~ dnorm(0, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(7000)
  sig2 <- sig^2
  tau <- 1/sig2
  # Set epsilon to avoid negative or small values
  epsilon <- 0.01
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) * ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) > epsilon) + 
    ((exp(la) * temp[i]^2 - exp(b.l) * temp[i] + c) <= epsilon) * epsilon
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "n_inv_mean_lambda_nt.txt")
  
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
                  model.file = "n_inv_mean_lambda_nt.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    N_param_list5_notrunc[[1]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list5_notrunc[[1]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list5_notrunc[[1]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list5_notrunc[[1]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for b in log space
    N_param_list5_notrunc[[2]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list5_notrunc[[2]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list5_notrunc[[2]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list5_notrunc[[2]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list5_notrunc[[3]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list5_notrunc[[3]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list5_notrunc[[3]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list5_notrunc[[3]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list5_notrunc[[4]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list5_notrunc[[4]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list5_notrunc[[4]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list5_notrunc[[4]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for tau
    N_param_list5_notrunc[[5]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list5_notrunc[[5]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list5_notrunc[[5]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list5_notrunc[[5]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
  }
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5_N.NT = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 13, to = 29, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new5_N.NT as a function at each of those temperatures
xdf5_N.NT <- expand.grid(point = xseq, iteration = 1:nrow(df.new5_N.NT))

# Apply the quadratic equation for each combination
xdf5_N.NT$value <- apply(xdf5_N.NT, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new5_N.NT[i, 4]) * location^2 - exp(df.new5_N.NT[i, 1]) * location + df.new5_N.NT[i, 2]
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf5_N.NT <- xdf5_N.NT|> mutate(trunc.eval = ifelse(xdf5_N.NT$value > epsilon, xdf5_N.NT$value, epsilon), 
                               trunc.inv = 1/trunc.eval)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf5_NT <- xdf5_N.NT |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
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
N_mean.xdf5_NT$true_curve <- true_curve(xseq)
N_mean.xdf5_NT$true_curve.inv <- true_curve.inv(xseq)
# Plot mean mortality rate response
N_plot5_mean_NT <- ggplot(N_mean.xdf5_NT, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median mortality rate response
N_plot5_med_NT <- ggplot(N_mean.xdf5_NT, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve), color = "violetred1", linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean lifetime response
N_plot5_mean.inv_NT <- ggplot(N_mean.xdf5_NT, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median lifetime response
N_plot5_med.inv_NT <- ggplot(N_mean.xdf5_NT, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkslategray2", alpha = 0.3) +
  geom_line(aes(y = true_curve.inv), color = "violetred1", linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = N_data.raw) +
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(N_plot5_mean_NT, N_plot5_med_NT, nrow = 1)
gridExtra::grid.arrange(N_plot5_mean.inv_NT, N_plot5_med.inv_NT, nrow = 1)

# Create a list to store histograms of the posterior distribution of each parameter
N_hist_list5_notrunc <- vector("list", length = 5)
par(mfrow = c(2, 3))
# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  N_first_column_list5_notrunc <- lapply(N_param_list5_notrunc[[i]], function(matrix) matrix[, 1])
  # Combine the vectors into a single vector
  N_combined_vector5_notrunc <- unlist(N_first_column_list5_notrunc)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(N_combined_vector5_notrunc) - 1, max(N_combined_vector5_notrunc) + 1, length = 1000)
  # Create a histogram and store it in the list
  N_hist_list5_notrunc[[i]] <- hist(N_combined_vector5_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                    xlab = "Values", col = "lavender", border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = "forestgreen", lwd = 2)
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
    lines(x, dexp(x, 7000), col = "purple", lty = 2, lwd = 2)
  }
}

N_proportions_list_a5_notrunc <- lapply(N_param_list5_notrunc[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_a5_notrunc))

N_proportions_list_b5_notrunc <- lapply(N_param_list5_notrunc[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.b) > col3 & log(true.b) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_b5_notrunc))

N_proportions_list_c5_notrunc <- lapply(N_param_list5_notrunc[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(N_proportions_list_c5_notrunc))

N_calculate_rmse_a_lambda_inv_mean_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_b_lambda_inv_mean_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.b), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

N_calculate_rmse_c_lambda_inv_mean_notrunc <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
N_rmse_values_a_lambda_inv_mean_notrunc <- sapply(N_param_list5_notrunc[[1]], N_calculate_rmse_a_lambda_inv_mean_notrunc)
N_rmse_values_b_lambda_inv_mean_notrunc <- sapply(N_param_list5_notrunc[[2]], N_calculate_rmse_b_lambda_inv_mean_notrunc)
N_rmse_values_c_lambda_inv_mean_notrunc <- sapply(N_param_list5_notrunc[[3]], N_calculate_rmse_c_lambda_inv_mean_notrunc)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
N_tab.inv.mean.lambda_notrunc <- matrix(0, nrow = 3, ncol = 6)
row.names(N_tab.inv.mean.lambda_notrunc) <- c("a", "b", "c")
colnames(N_tab.inv.mean.lambda_notrunc) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
N_tab.inv.mean.lambda_notrunc[1, 1] <- true.a
N_tab.inv.mean.lambda_notrunc[2, 1] <- true.b
N_tab.inv.mean.lambda_notrunc[3, 1] <- true.c
N_tab.inv.mean.lambda_notrunc[1, 2] <- exp(mean(unlist(lapply(N_param_list5_notrunc[[1]], function(matrix) matrix[, 1]))))
N_tab.inv.mean.lambda_notrunc[2, 2] <- exp(mean(unlist(lapply(N_param_list5_notrunc[[2]], function(matrix) matrix[, 1]))))
N_tab.inv.mean.lambda_notrunc[3, 2] <- mean(unlist(lapply(N_param_list5_notrunc[[3]], function(matrix) matrix[, 1])))
N_tab.inv.mean.lambda_notrunc[1, 3] <- exp(hdi(unlist(lapply(N_param_list5_notrunc[[1]], function(matrix) matrix[, 1])))[1])
N_tab.inv.mean.lambda_notrunc[2, 3] <- exp(hdi(unlist(lapply(N_param_list5_notrunc[[2]], function(matrix) matrix[, 1])))[1])
N_tab.inv.mean.lambda_notrunc[3, 3] <- hdi(unlist(lapply(N_param_list5_notrunc[[3]], function(matrix) matrix[, 1])))[1]
N_tab.inv.mean.lambda_notrunc[1, 4] <- exp(hdi(unlist(lapply(N_param_list5_notrunc[[1]], function(matrix) matrix[, 1])))[2])
N_tab.inv.mean.lambda_notrunc[2, 4] <- exp(hdi(unlist(lapply(N_param_list5_notrunc[[2]], function(matrix) matrix[, 1])))[2])
N_tab.inv.mean.lambda_notrunc[3, 4] <- hdi(unlist(lapply(N_param_list5_notrunc[[3]], function(matrix) matrix[, 1])))[2]
N_tab.inv.mean.lambda_notrunc[1, 5] <- mean(unlist(N_proportions_list_a5_notrunc))
N_tab.inv.mean.lambda_notrunc[2, 5] <- mean(unlist(N_proportions_list_b5_notrunc))
N_tab.inv.mean.lambda_notrunc[3, 5] <- mean(unlist(N_proportions_list_c5_notrunc))
N_tab.inv.mean.lambda_notrunc[1, 6] <- mean(N_rmse_values_a_lambda_inv_mean_notrunc)
N_tab.inv.mean.lambda_notrunc[2, 6] <- mean(N_rmse_values_b_lambda_inv_mean_notrunc)
N_tab.inv.mean.lambda_notrunc[3, 6] <- mean(N_rmse_values_c_lambda_inv_mean_notrunc)

save(N_param_list5_notrunc,
     xdf5_N.NT,
     df.new5_N.NT,
     N_mean.xdf5_NT,
     N_hist_list5_notrunc,
     N_plot5_mean_NT, 
     N_plot5_med_NT,
     N_plot5_mean_NT, 
     N_plot5_med.inv_NT,
     N_proportions_list_a5_notrunc,
     N_proportions_list_b5_notrunc,
     N_proportions_list_c5_notrunc,
     N_rmse_values_a_lambda_inv_mean_notrunc,
     N_rmse_values_b_lambda_inv_mean_notrunc,
     N_rmse_values_c_lambda_inv_mean_notrunc,
     N_tab.inv.mean.lambda_notrunc,
     file="inverse_mean_exp_narrow_NT.RData")


