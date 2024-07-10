library(ggplot2)
library(HDInterval)
library(rjags)
library(R2jags)
library(dplyr)

# Set the true parameters from the mortality curve
true.a  <-  0.001094
true.Topt <- 22.3
true.c <- 0.066

# True quadratic
f <- function(x, a = 0.001094, mu = 22.3, c = 0.066){
  a * (x - mu)^2 + c
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

# True value vector with a and b in log space
# Spot for sigma and deviance
true_values <- c(log(true.a), true.Topt, true.c, NA, NA)

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
  la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dexp(mu[i])
  }
      }", file = "lambda.txt")
  # Settings
  parameters <- c("la", "Topt", "c")
  inits <- function(){list(
    # a and b in log space
    la = log(0.1), 
    Topt = 20, 
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
    param_list1[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list1[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list1[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list1[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    param_list1[[2]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list1[[2]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list1[[2]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list1[[2]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list1[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list1[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list1[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list1[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for the deviance
    param_list1[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list1[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list1[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list1[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
  }
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new1 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new1 as a function at each of those temperatures
xdf1 <- expand.grid(point = xseq, iteration = 1:nrow(df.new1))

# Apply the quadratic equation for each combination
xdf1$value <- apply(xdf1, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new1[i, 3]) * (location - df.new1[i, 4])^2  + df.new1[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf1 <-  xdf1|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf1 <- xdf1 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                              avg.value = mean(value), 
                                              med.value.inv = median(trunc.inv), 
                                              med.value = median(value), 
                                              lower.hdi.inv = hdi(trunc.inv)[1], 
                                              upper.hdi.inv = hdi(trunc.inv)[2], 
                                              lower.hdi = hdi(value)[1], 
                                              upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf1$true_curve <- f(xseq)
# Plot the mean mortality rate response
plot1_mean <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = avg.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
labs(title = "Mean Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal()
# Plot the median mortality rate response
plot1_med <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = med.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
labs(title = "Median Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal()
# Plot the mean lifetime response
plot1_mean.inv <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = avg.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
geom_point(aes(x = T, y = trait), data = data.raw) + 
labs(title = "Mean Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Lifetime") + 
theme_minimal()
# Plot the median lifetime response
plot1_med.inv <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = med.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
geom_point(aes(x = T, y = trait), data = data.raw) + 
labs(title = "Median Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Lifetime") + 
theme_minimal()

gridExtra::grid.arrange(plot1_mean, plot1_med, nrow = 1)
gridExtra::grid.arrange(plot1_mean.inv, plot1_med.inv, nrow = 1)
ggsave("mortality_hdi.png", plot = plot1_mean)
ggsave("lifetime_hdi.png", plot = plot1_mean.inv)

# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
true_values <- c(log(true.a), true.Topt, true.c, NA)
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
                          xlab = "Values", col = rocket(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = rocket(12)[1], lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 22, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = rocket(12)[5], lty = 2, lwd = 2)
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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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
  predicted <- rep(true.Topt, length(observed))  
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
row.names(tab.lambda) <- c("a", "Topt", "c")
colnames(tab.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.lambda[1, 1] <- true.a
tab.lambda[2, 1] <- true.Topt
tab.lambda[3, 1] <- true.c
tab.lambda[1, 2] <- exp(mean(unlist(lapply(param_list1[[1]], function(matrix) matrix[, 1]))))
tab.lambda[2, 2] <- mean(unlist(lapply(param_list1[[2]], function(matrix) matrix[, 1])))
tab.lambda[3, 2] <- mean(unlist(lapply(param_list1[[3]], function(matrix) matrix[, 1])))
tab.lambda[1, 3] <- exp(hdi(unlist(lapply(param_list1[[1]], function(matrix) matrix[, 1])))[1])
tab.lambda[2, 3] <- hdi(unlist(lapply(param_list1[[2]], function(matrix) matrix[, 1])))[1]
tab.lambda[3, 3] <- hdi(unlist(lapply(param_list1[[3]], function(matrix) matrix[, 1])))[1]
tab.lambda[1, 4] <- exp(hdi(unlist(lapply(param_list1[[1]], function(matrix) matrix[, 1])))[2])
tab.lambda[2, 4] <- hdi(unlist(lapply(param_list1[[2]], function(matrix) matrix[, 1])))[2]
tab.lambda[3, 4] <- hdi(unlist(lapply(param_list1[[3]], function(matrix) matrix[, 1])))[2]
tab.lambda[1, 5] <- mean(unlist(proportions_list_a1))
tab.lambda[2, 5] <- mean(unlist(proportions_list_b1))
tab.lambda[3, 5] <- mean(unlist(proportions_list_c1))
tab.lambda[1, 6] <- mean(rmse_values_a_lambda)
tab.lambda[2, 6] <- mean(rmse_values_b_lambda)
tab.lambda[3, 6] <- mean(rmse_values_c_lambda)

save(param_list1,
     xdf1,
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
  la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10)
  c ~ dexp(0.5) # Has to be positive
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dexp(mu[i])
  }
      }", file = "mean_lambda.txt")
  # Settings
  parameters <- c("la", "Topt", "c")
  inits <- function(){list(
    # a and b in log space
    la = log(0.1), 
    Topt = 20, 
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
    param_list2[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list2[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list2[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list2[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt 
    param_list2[[2]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list2[[2]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list2[[2]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list2[[2]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list2[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list2[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list2[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list2[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list2[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list2[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list2[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list2[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
  }
  }
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new2 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new2 as a function at each of those temperatures
xdf2 <- expand.grid(point = xseq, iteration = 1:nrow(df.new2))

# Apply the quadratic equation for each combination
xdf2$value <- apply(xdf2, 1, function(row) {
i <- row["iteration"]
location <- row["point"]
(exp(df.new2[i, 3]) * (location - df.new2[i, 4])^2  + df.new2[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf2 <- xdf2|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf2 <- xdf2 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                              avg.value = mean(value), 
                                              med.value.inv = median(trunc.inv), 
                                              med.value = median(value), 
                                              lower.hdi.inv = hdi(trunc.inv)[1], 
                                              upper.hdi.inv = hdi(trunc.inv)[2], 
                                              lower.hdi = hdi(value)[1], 
                                              upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf2$true_curve <- f(xseq)
# Plot mean mortality rate response
plot2_mean <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = avg.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
labs(title = "Mean Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal()
plot2_med <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = med.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
labs(title = "Median Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal()
# Plot mean lifetime response
plot2_mean.inv <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = avg.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
geom_point(aes(x = T, y = trait), data = data.raw) + 
labs(title = "Mean Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Lifetime") + 
theme_minimal()
# Plot median lifetime response
plot2_med.inv <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = med.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
geom_point(aes(x = T, y = trait), data = data.raw) + 
labs(title = "Median Curve with Interval Bands and True Curve", 
     x = "Temperature", 
     y = "Lifetime") + 
theme_minimal()

gridExtra::grid.arrange(plot2_mean, plot2_med, nrow = 1)
gridExtra::grid.arrange(plot2_mean.inv, plot2_med.inv, nrow = 1)

# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
true_values <- c(log(true.a), true.Topt, true.c, NA)
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
                        xlab = "Values", col = rocket(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = rocket(12)[1], lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
    }
  if(i == 2){
    lines(x, dnorm(x, 22, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
    }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = rocket(12)[5], lty = 2, lwd = 2)
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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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
  predicted <- rep(true.Topt, length(observed))  
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
row.names(tab.lambda.mean) <- c("a", "Topt", "c")
colnames(tab.lambda.mean) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.lambda.mean[1, 1] <- true.a
tab.lambda.mean[2, 1] <- true.Topt
tab.lambda.mean[3, 1] <- true.c
tab.lambda.mean[1, 2] <- exp(mean(unlist(lapply(param_list2[[1]], function(matrix) matrix[, 1]))))
tab.lambda.mean[2, 2] <- mean(unlist(lapply(param_list2[[2]], function(matrix) matrix[, 1])))
tab.lambda.mean[3, 2] <- mean(unlist(lapply(param_list2[[3]], function(matrix) matrix[, 1])))
tab.lambda.mean[1, 3] <- exp(hdi(unlist(lapply(param_list2[[1]], function(matrix) matrix[, 1])))[1])
tab.lambda.mean[2, 3] <- hdi(unlist(lapply(param_list2[[2]], function(matrix) matrix[, 1])))[1]
tab.lambda.mean[3, 3] <- hdi(unlist(lapply(param_list2[[3]], function(matrix) matrix[, 1])))[1]
tab.lambda.mean[1, 4] <- exp(hdi(unlist(lapply(param_list2[[1]], function(matrix) matrix[, 1])))[2])
tab.lambda.mean[2, 4] <- hdi(unlist(lapply(param_list2[[2]], function(matrix) matrix[, 1])))[2]
tab.lambda.mean[3, 4] <- hdi(unlist(lapply(param_list2[[3]], function(matrix) matrix[, 1])))[2]
tab.lambda.mean[1, 5] <- mean(unlist(proportions_list_a2))
tab.lambda.mean[2, 5] <- mean(unlist(proportions_list_b2))
tab.lambda.mean[3, 5] <- mean(unlist(proportions_list_c2))
tab.lambda.mean[1, 6] <- mean(rmse_values_a_lambda_mean)
tab.lambda.mean[2, 6] <- mean(rmse_values_b_lambda_mean)
tab.lambda.mean[3, 6] <- mean(rmse_values_c_lambda_mean)

save(param_list2,
     xdf2,
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
  la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(0.5)
  sig2 <- sig^2
  tau <- 1/sig2
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <-  exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau) T(0,)
  }
}", file = "inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "Topt", "c", "sig")
  inits <- function(){list(
    # a and b in log space
    la = log(0.1), 
    Topt = 20, 
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
    param_list3[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list3[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list3[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list3[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    param_list3[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list3[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list3[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list3[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list3[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list3[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list3[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list3[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list3[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list3[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list3[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list3[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    param_list3[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list3[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list3[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list3[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new3 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new3 as a function at each of those temperatures
xdf3 <- expand.grid(point = xseq, iteration = 1:nrow(df.new3))

# Apply the quadratic equation for each combination
xdf3$value <- apply(xdf3, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new3[i, 3]) * (location - df.new3[i, 5])^2  + df.new3[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf3 <- xdf3|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf3 <- xdf3 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                avg.value = mean(value), 
                                                med.value.inv = median(trunc.inv), 
                                                med.value = median(value), 
                                                lower.hdi.inv = hdi(trunc.inv)[1], 
                                                upper.hdi.inv = hdi(trunc.inv)[2], 
                                                lower.hdi = hdi(value)[1], 
                                                upper.hdi = hdi(value)[2]) 
# Creating columns of the true curve values and the true inverted curve values
mean.xdf3$true_curve <- f(xseq)
# Plot mean mortality rate response
plot3_mean <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot3_med <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot3_mean.inv <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()
# Plot median lifetime response
plot3_med.inv <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()

gridExtra::grid.arrange(plot3_mean, plot3_med, nrow = 1)
gridExtra::grid.arrange(plot3_mean.inv, plot3_med.inv, nrow = 1)
ggsave("inv_mortality_hdi.png", plot = plot3_mean)
ggsave("inv_lifetime_hdi.png", plot = plot3_mean.inv)

# Assign true values 
# Create a list to store histograms of the posterior distribution of each parameter
true_values <- c(log(true.a), true.Topt, true.c, NA, NA)
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
                          xlab = "Values", col = rocket(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = rocket(12)[1], lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 22, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 0.5), col = rocket(12)[5], lty = 2, lwd = 2)
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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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
  predicted <- rep(true.Topt, length(observed))  
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
row.names(tab.inv.lambda) <- c("a", "Topt", "c")
colnames(tab.inv.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.inv.lambda[1, 1] <- true.a
tab.inv.lambda[2, 1] <- true.Topt
tab.inv.lambda[3, 1] <- true.c
tab.inv.lambda[1, 2] <- exp(mean(unlist(lapply(param_list3[[1]], function(matrix) matrix[, 1]))))
tab.inv.lambda[2, 2] <- mean(unlist(lapply(param_list3[[2]], function(matrix) matrix[, 1])))
tab.inv.lambda[3, 2] <- mean(unlist(lapply(param_list3[[3]], function(matrix) matrix[, 1])))
tab.inv.lambda[1, 3] <- exp(hdi(unlist(lapply(param_list3[[1]], function(matrix) matrix[, 1])))[1])
tab.inv.lambda[2, 3] <- hdi(unlist(lapply(param_list3[[2]], function(matrix) matrix[, 1])))[1]
tab.inv.lambda[3, 3] <- hdi(unlist(lapply(param_list3[[3]], function(matrix) matrix[, 1])))[1]
tab.inv.lambda[1, 4] <- exp(hdi(unlist(lapply(param_list3[[1]], function(matrix) matrix[, 1])))[2])
tab.inv.lambda[2, 4] <- hdi(unlist(lapply(param_list3[[2]], function(matrix) matrix[, 1])))[2]
tab.inv.lambda[3, 4] <- hdi(unlist(lapply(param_list3[[3]], function(matrix) matrix[, 1])))[2]
tab.inv.lambda[1, 5] <- mean(unlist(proportions_list_a3))
tab.inv.lambda[2, 5] <- mean(unlist(proportions_list_b3))
tab.inv.lambda[3, 5] <- mean(unlist(proportions_list_c3))
tab.inv.lambda[1, 6] <- mean(rmse_values_a_lambda_inv)
tab.inv.lambda[2, 6] <- mean(rmse_values_b_lambda_inv)
tab.inv.lambda[3, 6] <- mean(rmse_values_c_lambda_inv)

save(param_list3,
     xdf3,
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
  group_means <- tapply(data.list[[i]]$trait, data.list[[i]]$T, function(x) {
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
  la ~ dnorm(0, 1/10)
  Topt ~ dnorm(22, 1/10)
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(0.5)
  sig2 <- sig^2
  tau <- 1/sig2
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <-  exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau) T(0,)
  }
}", file = "mean_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "Topt", "c", "sig")
  inits <- function(){list(
    # a and b in log space
    la = log(0.1), 
    Topt = 20, 
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
    param_list4[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list4[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list4[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list4[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt 
    param_list4[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list4[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list4[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list4[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list4[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list4[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list4[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list4[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list4[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list4[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list4[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list4[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    param_list4[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list4[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list4[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list4[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new4 as a function at each of those temperatures
xdf4 <- expand.grid(point = xseq, iteration = 1:nrow(df.new4))

# Apply the quadratic equation for each combination
xdf4$value <- apply(xdf4, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new4[i, 3]) * (location - df.new4[i, 5])^2  + df.new4[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf4 <- xdf4|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf4 <- xdf4 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                  avg.value = mean(value), 
                                                  med.value.inv = median(trunc.inv), 
                                                  med.value = median(value), 
                                                  lower.hdi.inv = hdi(trunc.inv)[1], 
                                                  upper.hdi.inv = hdi(trunc.inv)[2], 
                                                  lower.hdi = hdi(value)[1], 
                                                  upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf4$true_curve <- f(xseq)
# Plot mean mortality rate response
plot4_mean <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot4_med <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot4_mean.inv <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()
# Plot median lifetime response
plot4_med.inv <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
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
                          xlab = "Values", col = rocket(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = rocket(12)[1], lwd = 2)
  # Add prior based on which parameter
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 22, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 0.5), col = rocket(12)[5], lty = 2, lwd = 2)
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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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
  predicted <- rep(true.Topt, length(observed))  
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
row.names(tab.mean.inv.lambda) <- c("a", "Topt", "c")
colnames(tab.mean.inv.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.mean.inv.lambda[1, 1] <- true.a
tab.mean.inv.lambda[2, 1] <- true.Topt
tab.mean.inv.lambda[3, 1] <- true.c
tab.mean.inv.lambda[1, 2] <- exp(mean(unlist(lapply(param_list4[[1]], function(matrix) matrix[, 1]))))
tab.mean.inv.lambda[2, 2] <- mean(unlist(lapply(param_list4[[2]], function(matrix) matrix[, 1])))
tab.mean.inv.lambda[3, 2] <- mean(unlist(lapply(param_list4[[3]], function(matrix) matrix[, 1])))
tab.mean.inv.lambda[1, 3] <- exp(hdi(unlist(lapply(param_list4[[1]], function(matrix) matrix[, 1])))[1])
tab.mean.inv.lambda[2, 3] <- hdi(unlist(lapply(param_list4[[2]], function(matrix) matrix[, 1])))[1]
tab.mean.inv.lambda[3, 3] <- hdi(unlist(lapply(param_list4[[3]], function(matrix) matrix[, 1])))[1]
tab.mean.inv.lambda[1, 4] <- exp(hdi(unlist(lapply(param_list4[[1]], function(matrix) matrix[, 1])))[2])
tab.mean.inv.lambda[2, 4] <- hdi(unlist(lapply(param_list4[[2]], function(matrix) matrix[, 1])))[2]
tab.mean.inv.lambda[3, 4] <- hdi(unlist(lapply(param_list4[[3]], function(matrix) matrix[, 1])))[2]
tab.mean.inv.lambda[1, 5] <- mean(unlist(proportions_list_a4))
tab.mean.inv.lambda[2, 5] <- mean(unlist(proportions_list_b4))
tab.mean.inv.lambda[3, 5] <- mean(unlist(proportions_list_c4))
tab.mean.inv.lambda[1, 6] <- mean(rmse_values_a_lambda_mean_inv)
tab.mean.inv.lambda[2, 6] <- mean(rmse_values_b_lambda_mean_inv)
tab.mean.inv.lambda[3, 6] <- mean(rmse_values_c_lambda_mean_inv)

save(param_list4,
     xdf4,
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
   la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10)
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(0.5)
  sig2 <- sig^2
  tau <- 1/sig2
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <-  exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau) T(0,)
  }
}", file = "inv_mean_lambda.txt")
  # Settings
  parameters <- c("la", "Topt", "c", "sig")
  inits <- function(){list(
    # a and b in log space
    la = log(0.1), 
    Topt = 20, 
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
    param_list5[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list5[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list5[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list5[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    param_list5[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list5[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list5[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list5[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    param_list5[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    param_list5[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    param_list5[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    param_list5[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    param_list5[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    param_list5[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    param_list5[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    param_list5[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    param_list5[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list5[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list5[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list5[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new5 as a function at each of those temperatures
xdf5 <- expand.grid(point = xseq, iteration = 1:nrow(df.new5))

# Apply the quadratic equation for each combination
xdf5$value <- apply(xdf5, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new5[i, 3]) * (location - df.new5[i, 5])^2  + df.new5[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf5 <- xdf5|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf5 <- xdf5 |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                avg.value = mean(value), 
                                                med.value.inv = median(trunc.inv), 
                                                med.value = median(value), 
                                                lower.hdi.inv = hdi(trunc.inv)[1], 
                                                upper.hdi.inv = hdi(trunc.inv)[2], 
                                                lower.hdi = hdi(value)[1], 
                                                upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf5$true_curve <- true_curve(xseq)
# Plot mean mortality rate response
plot5_mean <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot median mortality rate response
plot5_med <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(title = "Median Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal()
# Plot mean lifetime response
plot5_mean.inv <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = T, y = trait), data = data.raw) + 
  labs(title = "Mean Curve with Interval Bands and True Curve", 
       x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal()
# Plot median lifetime response
plot5_med.inv <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
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
                          xlab = "Values", col = rocket(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = true_values[i], col = rocket(12)[1], lwd = 2)
  # Add prior based on which parameter 
  if(i == 1){
    lines(x, dnorm(x, 0, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 2){
    lines(x, dnorm(x, 22, sqrt(10)), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 3){
    lines(x, dexp(x, 0.5), col = rocket(12)[5], lty = 2, lwd = 2)
  }
  if(i == 5){
    lines(x, dexp(x, 0.5), col = rocket(12)[5], lty = 2, lwd = 2)
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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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
  predicted <- rep(true.Topt, length(observed))  
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
row.names(tab.inv.mean.lambda) <- c("a", "Topt", "c")
colnames(tab.inv.mean.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
tab.inv.mean.lambda[1, 1] <- true.a
tab.inv.mean.lambda[2, 1] <- true.Topt
tab.inv.mean.lambda[3, 1] <- true.c
tab.inv.mean.lambda[1, 2] <- exp(mean(unlist(lapply(param_list5[[1]], function(matrix) matrix[, 1]))))
tab.inv.mean.lambda[2, 2] <- mean(unlist(lapply(param_list5[[2]], function(matrix) matrix[, 1])))
tab.inv.mean.lambda[3, 2] <- mean(unlist(lapply(param_list5[[3]], function(matrix) matrix[, 1])))
tab.inv.mean.lambda[1, 3] <- exp(hdi(unlist(lapply(param_list5[[1]], function(matrix) matrix[, 1])))[1])
tab.inv.mean.lambda[2, 3] <- hdi(unlist(lapply(param_list5[[2]], function(matrix) matrix[, 1])))[1]
tab.inv.mean.lambda[3, 3] <- hdi(unlist(lapply(param_list5[[3]], function(matrix) matrix[, 1])))[1]
tab.inv.mean.lambda[1, 4] <- exp(hdi(unlist(lapply(param_list5[[1]], function(matrix) matrix[, 1])))[2])
tab.inv.mean.lambda[2, 4] <- hdi(unlist(lapply(param_list5[[2]], function(matrix) matrix[, 1])))[2]
tab.inv.mean.lambda[3, 4] <- hdi(unlist(lapply(param_list5[[3]], function(matrix) matrix[, 1])))[2]
tab.inv.mean.lambda[1, 5] <- mean(unlist(proportions_list_a5))
tab.inv.mean.lambda[2, 5] <- mean(unlist(proportions_list_b5))
tab.inv.mean.lambda[3, 5] <- mean(unlist(proportions_list_c5))
tab.inv.mean.lambda[1, 6] <- mean(rmse_values_a_lambda_inv_mean)
tab.inv.mean.lambda[2, 6] <- mean(rmse_values_b_lambda_inv_mean)
tab.inv.mean.lambda[3, 6] <- mean(rmse_values_c_lambda_inv_mean)


save(param_list5,
     xdf5,
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


