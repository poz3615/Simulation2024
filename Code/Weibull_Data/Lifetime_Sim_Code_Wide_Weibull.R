library(rjags)
library(dplyr)
library(ggplot2)
library(R2jags)
library(HDInterval)

# Set true mortality curve parameters
# Parameters chosen to resemble curve with exponential data
true.a  <-  0.001094
true.Topt <- 22.3
true.c <- 0.066

# Making median lifetime function (concave down)
# Exponentiated to avoid negative values
my.med <- function(x, a = 0.001094, mu = 22.3, c = 0.066){
  1 / (a * (x - mu)^2 + c)
}

# Set reasonable shape parameter
sh.est <- 3
# Lambda function
my.lambda <- function(x, a = 0.001094, mu = 22.3, c = 0.066, sh.est = 3){ 
  my.med(x, a, mu, c)/(log(2))^(1/sh.est) # # fx0 is the median here
}


# Generate 100 datasets and store them in a list
Weibull_data.list <- list()
set.seed(52)
for(i in 1:100){
  T1W <- rweibull(n = 100, scale = my.lambda(10), shape = sh.est)
  T2W <- rweibull(n = 100, scale = my.lambda(17), shape = sh.est)
  T3W <- rweibull(n = 100, scale = my.lambda(21), shape = sh.est)
  T4W <- rweibull(n = 100, scale = my.lambda(25), shape = sh.est)
  T5W <- rweibull(n = 100, scale = my.lambda(32), shape = sh.est)
  Weibull_data.list[[i]] <- data.frame(
    T = c(rep(10, 100), rep(17, 100), rep(21, 100), rep(25, 100), rep(32, 100)), 
    trait = c(T1W, T2W, T3W, T4W, T5W)
  )
}

w <- length(Weibull_data.list)

# True value vector with a in log space
Weibull_true_values <- c(log(true.a), true.Topt, true.c, NA, sh.est)

# Raw data for HDI plot
data.raw.w <- Weibull_data.list[[100]]

##########################################Lambda #############################################

# Vector for RMSEs for each sample of each model
rmse.range.ind.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list1 <- vector("list", length = 4)
for (i in 1:4) {
  # Initialize the inner list
  Weibull_inner_list1 <- vector("list", length = w)
  # Populate the inner list
  for (j in 1:w) {
    # Initialize a 5x4 matrix
    Weibull_matrix_data1 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_inner_list1[[j]] <- Weibull_matrix_data1
  }
  Weibull_param_list1[[i]] <- Weibull_inner_list1
}
# Create JAGS model for each of the 100 datasets
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("Weibull_lambda.txt")
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
  }", file = "Weibull_lambda.txt")
  
  # Settings
  parameters <- c("la", "Topt", "c")

  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  inits = vector('list', nc)
  GenInits = function() {
    la <- rnorm(1, 0, 10)
    c <- rexp(1, 0.5)
    Topt <- rnorm(1, 22, 10)
    list(
      la = la,             
      c = c,
      Topt = Topt
    )
  }
  for(i in 1:nc){
    inits[[i]] = GenInits()
  }
  data.raw.w <- Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "Weibull_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    Weibull_param_list1[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    Weibull_param_list1[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    Weibull_param_list1[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    Weibull_param_list1[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    Weibull_param_list1[[2]][[i]][j, 1] <- mean(samp[[j]][, 4])
    Weibull_param_list1[[2]][[i]][j, 2] <- median(samp[[j]][, 4])
    Weibull_param_list1[[2]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    Weibull_param_list1[[2]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for c
    Weibull_param_list1[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    Weibull_param_list1[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    Weibull_param_list1[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    Weibull_param_list1[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    Weibull_param_list1[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    Weibull_param_list1[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    Weibull_param_list1[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    Weibull_param_list1[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
  }
  df.rmse1.wb <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df1.wb <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse1.wb))
  # Apply the quadratic equation for each combination
  rmse.df1.wb$value <- apply(rmse.df1.wb, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse1.wb[k, 3]) * (location - df.rmse1.wb[k, 4])^2 + df.rmse1.wb[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse1.wb), function(k) {
    evaluated_curve <- rmse.df1.wb$value[rmse.df1.wb$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.ind.wb <- c(rmse.range.ind.wb, rmse_values)  
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new1_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a sequence of temperatures from 10 to 35
xdf1_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new1_W))
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new1_W as a function at each of those temperatures
xdf1_W$value <- apply(xdf1_W, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new1_W[i, 3]) *  (location - df.new1_W[i, 4])^2 + df.new1_W[i, 1]
})
# To keep consistent with exponential data
xdf1_W <- xdf1_W|> mutate(trunc.inv = log(2)/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf1_W <- xdf1_W |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf1_W$true_curve <- f(xseq)
# Plot the mean lifetime response
plot1_mean_W <- ggplot(mean.xdf1_W, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot the median lifetime response
plot1_med_W <- ggplot(mean.xdf1_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot the mean mortality rate response
plot1_mean.inv_W <- ggplot(mean.xdf1_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  geom_line(aes(y = log(2)/true_curve), color = plasma(10)[2], linetype = "dashed") +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot the median mortality rate response
plot1_med.inv_W <- ggplot(mean.xdf1_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = log(2)/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(plot1_mean_W, plot1_med_W, nrow = 1)
gridExtra::grid.arrange(plot1_mean.inv_W, plot1_med.inv_W, nrow = 1)
ggsave("ind.mort.wb.png", plot = plot1_med_W, width = 5, height = 5)
ggsave("ind.life.wb.png", plot = plot1_med.inv_W, width = 5, height = 5)

# Create a list to store histograms of the posterior distribution of each parameter
Weibull_hist_list1 <- vector("list", length = 4)
parameter.names <- c("a", "Topt", "c", "Deviance")
png("hist.ind.wb.png", width = 8, height = 6, units = "in", res = 250)
par(mfrow = c(2, 2))
for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  Weibull_first_column_list1 <- lapply(Weibull_param_list1[[i]], function(matrix) matrix[, 2])
  # Combine the vectors into a single vector
  Weibull_combined_vector1 <- unlist(Weibull_first_column_list1)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(Weibull_combined_vector1) - 50, max(Weibull_combined_vector1) + 1, length = 1000)
  # Create a histogram and store it in the list
  Weibull_hist_list1[[i]] <- hist(Weibull_combined_vector1, main = paste("Posterior Medians of Parameter", parameter.names[i]), 
                          xlab = "Posterior Median", col = plasma(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = Weibull_true_values[i], col = plasma(12)[1], lwd = 2)

}
dev.off()


Weibull_proportions_list_a1 <- lapply(Weibull_param_list1[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_a1))

Weibull_proportions_list_b1 <- lapply(Weibull_param_list1[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_b1))

Weibull_proportions_list_c1 <- lapply(Weibull_param_list1[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_c1))

Weibull_calculate_rmse_a_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1] 
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_b_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.Topt, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_c_lambda <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda <- sapply(Weibull_param_list1[[1]], Weibull_calculate_rmse_a_lambda)
Weibull_rmse_values_b_lambda <- sapply(Weibull_param_list1[[2]], Weibull_calculate_rmse_b_lambda)
Weibull_rmse_values_c_lambda <- sapply(Weibull_param_list1[[3]], Weibull_calculate_rmse_c_lambda)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
Weibull_tab.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(Weibull_tab.lambda) <- c("a", "Topt", "c")
colnames(Weibull_tab.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
Weibull_tab.lambda[1, 1] <- true.a
Weibull_tab.lambda[2, 1] <- true.Topt
Weibull_tab.lambda[3, 1] <- true.c
Weibull_tab.lambda[1, 2] <- exp(mean(unlist(lapply(Weibull_param_list1[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.lambda[2, 2] <- (mean(unlist(lapply(Weibull_param_list1[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.lambda[3, 2] <- mean(unlist(lapply(Weibull_param_list1[[3]], function(matrix) matrix[, 1])))
Weibull_tab.lambda[1, 3] <- exp(hdi(unlist(lapply(Weibull_param_list1[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.lambda[2, 3] <- (hdi(unlist(lapply(Weibull_param_list1[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.lambda[3, 3] <- hdi(unlist(lapply(Weibull_param_list1[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.lambda[1, 4] <- exp(hdi(unlist(lapply(Weibull_param_list1[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.lambda[2, 4] <- (hdi(unlist(lapply(Weibull_param_list1[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.lambda[3, 4] <- hdi(unlist(lapply(Weibull_param_list1[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.lambda[1, 5] <- mean(unlist(Weibull_proportions_list_a1))
Weibull_tab.lambda[2, 5] <- mean(unlist(Weibull_proportions_list_b1))
Weibull_tab.lambda[3, 5] <- mean(unlist(Weibull_proportions_list_c1))
Weibull_tab.lambda[1, 6] <- mean(Weibull_rmse_values_a_lambda)
Weibull_tab.lambda[2, 6] <- mean(Weibull_rmse_values_b_lambda)
Weibull_tab.lambda[3, 6] <- mean(Weibull_rmse_values_c_lambda)

save(Weibull_param_list1,
     xdf1_W,
     df.new1_W,
     mean.xdf1_W,
     Weibull_hist_list1,
     plot1_mean_W, 
     plot1_med_W,
     plot1_mean.inv_W, 
     plot1_med.inv_W,
     Weibull_proportions_list_a1,
     Weibull_proportions_list_b1,
     Weibull_proportions_list_c1,
     Weibull_rmse_values_a_lambda,
     Weibull_rmse_values_b_lambda,
     Weibull_rmse_values_c_lambda,
     Weibull_tab.lambda,
     rmse.range.ind.wb,
     file="Weibull_individual.RData")


##########################################Mean Lambda #############################################

# Vector for RMSEs for each sample of each model
rmse.range.mean.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list2 <- vector("list", length = 4)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_inner_list2 <- vector("list", length = w)
  # Populate the inner list
  for (j in 1:w) {
    # Initialize a 5x4 matrix
    Weibull_matrix_data2 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_inner_list2[[j]] <- Weibull_matrix_data2
  }
  Weibull_param_list2[[i]] <- Weibull_inner_list2
}
# Create JAGS model for each of the 100 datasets
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(Weibull_data.list[[i]]$trait, Weibull_data.list[[i]]$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)),  
    trait = unlist(group_means)  
  )
  # Model
  sink("Weibull_mean_lambda.txt")
  cat("model{
  # Priors
   la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10)
  c ~ dexp(10) # Has to be positive
  sig ~ dexp(0.5)
  sig2 <- sig^2
  tau <- 1/sig2
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <-  exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau) T(0,)
  }
}", file = "Weibull_mean_lambda.txt")
  # Settings
  parameters <- c("la", "Topt", "c", "sig")

  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  inits = vector('list', nc)
  GenInits = function() {
    la <- rnorm(1, 0, 10)
    c <- rexp(1, 10)
    Topt <- rnorm(1, 22, 10)
    sig <- rexp(1, 0.5)
    list(
      la = la,             
      c = c,
      Topt = Topt,
      sig = sig
    )
  }
  for(i in 1:nc){
    inits[[i]] = GenInits()
  }
  
  data.raw.w <- Weibull_data.list[[i]]
  data <- df.mean 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "Weibull_mean_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    Weibull_param_list2[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    Weibull_param_list2[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    Weibull_param_list2[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    Weibull_param_list2[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    Weibull_param_list2[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    Weibull_param_list2[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    Weibull_param_list2[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    Weibull_param_list2[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    Weibull_param_list2[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    Weibull_param_list2[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    Weibull_param_list2[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    Weibull_param_list2[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    Weibull_param_list2[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    Weibull_param_list2[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    Weibull_param_list2[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    Weibull_param_list2[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    Weibull_param_list2[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    Weibull_param_list2[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    Weibull_param_list2[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    Weibull_param_list2[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse2.wb <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df2.wb <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse2.wb))
  # Apply the quadratic equation for each combination
  rmse.df2.wb$value <- apply(rmse.df2.wb, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse2.wb[k, 3]) * (location - df.rmse2.wb[k, 5])^2 + df.rmse2.wb[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse2.wb), function(k) {
    evaluated_curve <- rmse.df2.wb$value[rmse.df2.wb$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.mean.wb <- c(rmse.range.mean.wb, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new2_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new2_W as a function at each of those temperatures
xdf2_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new2_W))

# Apply the quadratic equation for each combination
xdf2_W$value <- apply(xdf2_W, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new2_W[i, 3]) *  (location - df.new2_W[i, 5])^2 + df.new2_W[i, 1]
})
# To keep consistent with exponential data
xdf2_W <- xdf2_W|> mutate(trunc.inv = log(2)/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf2_W <- xdf2_W |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf2_W$true_curve <- f(xseq)
# Plot mean lifetime response
plot2_mean_W <- ggplot(mean.xdf2_W, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  #geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median lifetime response
plot2_med_W <- ggplot(mean.xdf2_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  #geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean mortality rate response
plot2_mean.inv_W <- ggplot(mean.xdf2_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = log(2)/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median mortality rate response
plot2_med.inv_W <- ggplot(mean.xdf2_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = log(2)/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(plot2_mean_W, plot2_med_W, nrow = 1)
gridExtra::grid.arrange(plot2_mean.inv_W, plot2_med.inv_W, nrow = 1)
ggsave("mean.mort.wb.png", plot = plot2_med_W, width = 5, height = 5)
ggsave("mean.life.wb.png", plot = plot2_med.inv_W, width = 5, height = 5)

# Create a list to store histograms of the posterior distribution of each parameter
Weibull_hist_list2 <- vector("list", length = 4)
parameter.names.sig <- c("a", "Topt", "c", "Deviance", "Sigma")
png("hist.mean.wb.png", width = 8, height = 6, units = "in", res = 250)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  Weibull_first_column_list2 <- lapply(Weibull_param_list2[[i]], function(matrix) matrix[, 2])
  # Combine the vectors into a single vector
  Weibull_combined_vector2 <- unlist(Weibull_first_column_list2)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(Weibull_combined_vector2) - 50, max(Weibull_combined_vector2) + 1, length = 1000)
  # Create a histogram and store it in the list
  Weibull_hist_list2[[i]] <- hist(Weibull_combined_vector2, main = paste("Posterior Medians of Parameter", parameter.names.sig[i]), 
                          xlab = "Posterior Median", col = plasma(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = Weibull_true_values[i], col = plasma(12)[1], lwd = 2)
 
}
dev.off()

Weibull_proportions_list_a2 <- lapply(Weibull_param_list2[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_a2))

Weibull_proportions_list_b2 <- lapply(Weibull_param_list2[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_b2))

Weibull_proportions_list_c2 <- lapply(Weibull_param_list2[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_c2))

Weibull_calculate_rmse_a_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_b_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.Topt, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_c_lambda_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda_mean <- sapply(Weibull_param_list2[[1]], Weibull_calculate_rmse_a_lambda_mean)
Weibull_rmse_values_b_lambda_mean <- sapply(Weibull_param_list2[[2]], Weibull_calculate_rmse_b_lambda_mean)
Weibull_rmse_values_c_lambda_mean <- sapply(Weibull_param_list2[[3]], Weibull_calculate_rmse_c_lambda_mean)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
Weibull_tab.lambda.mean <- matrix(0, nrow = 3, ncol = 6)
row.names(Weibull_tab.lambda.mean) <- c("a", "Topt", "c")
colnames(Weibull_tab.lambda.mean) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
Weibull_tab.lambda.mean[1, 1] <- true.a
Weibull_tab.lambda.mean[2, 1] <- true.Topt
Weibull_tab.lambda.mean[3, 1] <- true.c
Weibull_tab.lambda.mean[1, 2] <- exp(mean(unlist(lapply(Weibull_param_list2[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.lambda.mean[2, 2] <- (mean(unlist(lapply(Weibull_param_list2[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.lambda.mean[3, 2] <- mean(unlist(lapply(Weibull_param_list2[[3]], function(matrix) matrix[, 1])))
Weibull_tab.lambda.mean[1, 3] <- exp(hdi(unlist(lapply(Weibull_param_list2[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.lambda.mean[2, 3] <- (hdi(unlist(lapply(Weibull_param_list2[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.lambda.mean[3, 3] <- hdi(unlist(lapply(Weibull_param_list2[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.lambda.mean[1, 4] <- exp(hdi(unlist(lapply(Weibull_param_list2[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.lambda.mean[2, 4] <- (hdi(unlist(lapply(Weibull_param_list2[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.lambda.mean[3, 4] <- hdi(unlist(lapply(Weibull_param_list2[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.lambda.mean[1, 5] <- mean(unlist(Weibull_proportions_list_a2))
Weibull_tab.lambda.mean[2, 5] <- mean(unlist(Weibull_proportions_list_b2))
Weibull_tab.lambda.mean[3, 5] <- mean(unlist(Weibull_proportions_list_c2))
Weibull_tab.lambda.mean[1, 6] <- mean(Weibull_rmse_values_a_lambda_mean)
Weibull_tab.lambda.mean[2, 6] <- mean(Weibull_rmse_values_b_lambda_mean)
Weibull_tab.lambda.mean[3, 6] <- mean(Weibull_rmse_values_c_lambda_mean)

save(Weibull_param_list2,
     xdf2_W,
     df.new2_W,
     mean.xdf2_W,
     Weibull_hist_list2,
     plot2_mean_W, 
     plot2_med_W,
     plot2_mean.inv_W, 
     plot2_med.inv_W,
     Weibull_proportions_list_a2,
     Weibull_proportions_list_b2,
     Weibull_proportions_list_c2,
     Weibull_rmse_values_a_lambda_mean,
     Weibull_rmse_values_b_lambda_mean,
     Weibull_rmse_values_c_lambda_mean,
     Weibull_tab.lambda.mean,
     rmse.range.mean.wb,
     file="Weibull_mean.RData")

##########################################Inverse Lambda #############################################

# Vector for RMSEs for each sample of each model
rmse.range.inv.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list3 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_inner_list3 <- vector("list", length = w)
  # Populate the inner list
  for (j in 1:w) {
    # Initialize a 5x4 matrix
    Weibull_matrix_data3 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_inner_list3[[j]] <- Weibull_matrix_data3
  }
  Weibull_param_list3[[i]] <- Weibull_inner_list3
}
# Create JAGS model for each of the 100 datasets
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("Weibull_inv_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10)
  c ~ dexp(0.5) # Has to be positive
  sig ~ dexp(0.5)
  sig2 <- sig^2
  tau <- 1/sig2
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "Weibull_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "Topt", "c", "sig")

  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  inits = vector('list', nc)
  GenInits = function() {
    la <- rnorm(1, 0, 10)
    c <- rexp(1, 0.5)
    Topt <- rnorm(1, 22, 10)
    sig <- rexp(1, 0.5)
    list(
      la = la,             
      c = c,
      Topt = Topt,
      sig = sig
    )
  }
  for(i in 1:nc){
    inits[[i]] = GenInits()
  }
  
  data.raw.w <- Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "Weibull_inv_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    Weibull_param_list3[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    Weibull_param_list3[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    Weibull_param_list3[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    Weibull_param_list3[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    Weibull_param_list3[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    Weibull_param_list3[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    Weibull_param_list3[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    Weibull_param_list3[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    Weibull_param_list3[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    Weibull_param_list3[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    Weibull_param_list3[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    Weibull_param_list3[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    Weibull_param_list3[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    Weibull_param_list3[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    Weibull_param_list3[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    Weibull_param_list3[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    Weibull_param_list3[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    Weibull_param_list3[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    Weibull_param_list3[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    Weibull_param_list3[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse3.wb <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df3.wb <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse3.wb))
  # Apply the quadratic equation for each combination
  rmse.df3.wb$value <- apply(rmse.df3.wb, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse3.wb[k, 3]) * (location - df.rmse3.wb[k, 5])^2 + df.rmse3.wb[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse3.wb), function(k) {
    evaluated_curve <- rmse.df3.wb$value[rmse.df3.wb$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.inv.wb <- c(rmse.range.inv.wb, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new3_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new3_W as a function at each of those temperatures
xdf3_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new3_W))

# Apply the quadratic equation for each combination
xdf3_W$value <- apply(xdf3_W, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new3_W[i, 3]) * (location - df.new3_W[i, 5])^2 + df.new3_W[i, 1]
})
# To keep consistent with exponential data
xdf3_W <- xdf3_W|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf3_W <- xdf3_W |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf3_W$true_curve <- f(xseq)
# Plot mean lifetime response
plot3_mean_W <- ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median lifetime response
plot3_med_W <- ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean mortality rate response
plot3_mean.inv_W <- ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median mortality rate response
plot3_med.inv_W <- ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()


gridExtra::grid.arrange(plot3_mean_W, plot3_med_W, nrow = 1)
gridExtra::grid.arrange(plot3_mean.inv_W, plot3_med.inv_W, nrow = 1)
ggsave("inv.mort.wb.png", plot = plot3_med_W, width = 5, height = 5)
ggsave("inv.life.wb.png", plot = plot3_med.inv_W, width = 5, height = 5)


# Create a list to store histograms of the posterior distribution of each parameter
Weibull_hist_list3 <- vector("list", length = 5)
parameter.names.sig <- c("a", "Topt", "c", "Deviance", "Sigma")
png("hist.inv.wb.png", width = 10, height = 6, units = "in", res = 250)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  Weibull_first_column_list3 <- lapply(Weibull_param_list3[[i]], function(matrix) matrix[, 2])
  # Combine the vectors into a single vector
  Weibull_combined_vector3 <- unlist(Weibull_first_column_list3)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(Weibull_combined_vector3) - 1, max(Weibull_combined_vector3) + 1, length = 1000)
  # Create a histogram and store it in the list
  Weibull_hist_list3[[i]] <- hist(Weibull_combined_vector3, main = paste("Posterior Medians of Parameter", parameter.names.sig[i]), 
                          xlab = "Posterior Median", col = plasma(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = Weibull_true_values[i], col = plasma(12)[1], lwd = 2)
}
dev.off()


Weibull_proportions_list_a3 <- lapply(Weibull_param_list3[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_a3))

Weibull_proportions_list_b3 <- lapply(Weibull_param_list3[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_b3))

Weibull_proportions_list_c3 <- lapply(Weibull_param_list3[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_c3))

Weibull_calculate_rmse_a_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_b_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1] 
  # True values
  predicted <- rep(true.Topt, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_c_lambda_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1] 
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda_inv <- sapply(Weibull_param_list3[[1]], Weibull_calculate_rmse_a_lambda_inv)
Weibull_rmse_values_b_lambda_inv <- sapply(Weibull_param_list3[[2]], Weibull_calculate_rmse_b_lambda_inv)
Weibull_rmse_values_c_lambda_inv <- sapply(Weibull_param_list3[[3]], Weibull_calculate_rmse_c_lambda_inv)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
Weibull_tab.inv.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(Weibull_tab.inv.lambda) <- c("a", "Topt", "c")
colnames(Weibull_tab.inv.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
Weibull_tab.inv.lambda[1, 1] <- true.a
Weibull_tab.inv.lambda[2, 1] <- true.Topt
Weibull_tab.inv.lambda[3, 1] <- true.c
Weibull_tab.inv.lambda[1, 2] <- exp(mean(unlist(lapply(Weibull_param_list3[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.inv.lambda[2, 2] <- (mean(unlist(lapply(Weibull_param_list3[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.inv.lambda[3, 2] <- mean(unlist(lapply(Weibull_param_list3[[3]], function(matrix) matrix[, 1])))
Weibull_tab.inv.lambda[1, 3] <- exp(hdi(unlist(lapply(Weibull_param_list3[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.inv.lambda[2, 3] <- (hdi(unlist(lapply(Weibull_param_list3[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.inv.lambda[3, 3] <- hdi(unlist(lapply(Weibull_param_list3[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.inv.lambda[1, 4] <- exp(hdi(unlist(lapply(Weibull_param_list3[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.inv.lambda[2, 4] <- (hdi(unlist(lapply(Weibull_param_list3[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.inv.lambda[3, 4] <- hdi(unlist(lapply(Weibull_param_list3[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.inv.lambda[1, 5] <- mean(unlist(Weibull_proportions_list_a3))
Weibull_tab.inv.lambda[2, 5] <- mean(unlist(Weibull_proportions_list_b3))
Weibull_tab.inv.lambda[3, 5] <- mean(unlist(Weibull_proportions_list_c3))
Weibull_tab.inv.lambda[1, 6] <- mean(Weibull_rmse_values_a_lambda_inv)
Weibull_tab.inv.lambda[2, 6] <- mean(Weibull_rmse_values_b_lambda_inv)
Weibull_tab.inv.lambda[3, 6] <- mean(Weibull_rmse_values_c_lambda_inv)

save(Weibull_param_list3,
     xdf3_W,
     df.new3_W,
     mean.xdf3_W,
     Weibull_hist_list3,
     plot3_mean_W, 
     plot3_med_W,
     plot3_mean.inv_W, 
     plot3_med.inv_W,
     Weibull_proportions_list_a3,
     Weibull_proportions_list_b3,
     Weibull_proportions_list_c3,
     Weibull_rmse_values_a_lambda_inv,
     Weibull_rmse_values_b_lambda_inv,
     Weibull_rmse_values_c_lambda_inv,
     Weibull_tab.inv.lambda,
     rmse.range.inv.wb,
     file="Weibull_inverse.RData")


##########################################Mean Inverse Lambda #############################################


# Vector for RMSEs for each sample of each model
rmse.range.mean.inv.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list4 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_inner_list4 <- vector("list", length = w)
  # Populate the inner list
  for (j in 1:w) {
    # Initialize a 5x4 matrix
    Weibull_matrix_data4 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_inner_list4[[j]] <- Weibull_matrix_data4
  }
  Weibull_param_list4[[i]] <- Weibull_inner_list4
}
# Create JAGS model for each of the 100 datasets
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- Weibull_data.list[[i]] 
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data.raw.w <- Weibull_data.list[[i]]
  data <- df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("Weibull_mean_inv_lambda.txt")
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
    mu[i] <- exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "Weibull_mean_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "Topt", "c", "sig")

  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  inits = vector('list', nc)
  GenInits = function() {
    la <- rnorm(1, 0, 10)
    c <- rexp(1, 0.5)
    Topt <- rnorm(1, 22, 10)
    sig <- rexp(1, 0.5)
    list(
      la = la,             
      c = c,
      Topt = Topt,
      sig = sig
    )
  }
  for(i in 1:nc){
    inits[[i]] = GenInits()
  }
  
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "Weibull_mean_inv_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    Weibull_param_list4[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    Weibull_param_list4[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    Weibull_param_list4[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    Weibull_param_list4[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    Weibull_param_list4[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    Weibull_param_list4[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    Weibull_param_list4[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    Weibull_param_list4[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    Weibull_param_list4[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    Weibull_param_list4[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    Weibull_param_list4[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    Weibull_param_list4[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    Weibull_param_list4[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    Weibull_param_list4[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    Weibull_param_list4[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    Weibull_param_list4[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    Weibull_param_list4[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    Weibull_param_list4[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    Weibull_param_list4[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    Weibull_param_list4[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse4.wb <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df4.wb <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse4.wb))
  # Apply the quadratic equation for each combination
  rmse.df4.wb$value <- apply(rmse.df4.wb, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse4.wb[k, 3]) * (location - df.rmse4.wb[k, 5])^2 + df.rmse4.wb[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse4.wb), function(k) {
    evaluated_curve <- rmse.df4.wb$value[rmse.df4.wb$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.mean.inv.wb <- c(rmse.range.mean.inv.wb, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new4_W as a function at each of those temperatures
xdf4_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new4_W))

# Apply the quadratic equation for each combination
xdf4_W$value <- apply(xdf4_W, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new4_W[i, 3]) *  (location - df.new4_W[i, 5])^2 + df.new4_W[i, 1]
})
# To keep consistent with exponential data
xdf4_W <- xdf4_W|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf4_W <- xdf4_W |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf4_W$true_curve <- f(xseq)
# Plot mean lifetime response
plot4_mean_W <- ggplot(mean.xdf4_W, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median lifetime response
plot4_med_W <- ggplot(mean.xdf4_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean mortality rate response
plot4_mean.inv_W <- ggplot(mean.xdf4_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median mortality rate response
plot4_med.inv_W <- ggplot(mean.xdf4_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(plot4_mean_W, plot4_med_W, nrow = 1)
gridExtra::grid.arrange(plot4_mean.inv_W, plot4_med.inv_W, nrow = 1)
ggsave("mean.inv.mort.wb.png", plot = plot4_med_W, width = 5, height = 5)
ggsave("mean.inv.life.wb.png", plot = plot4_med.inv_W, width = 5, height = 5)


# Create a list to store histograms of the posterior distribution of each parameter
Weibull_hist_list4 <- vector("list", length = 4)
parameter.names.sig <- c("a", "Topt", "c", "Deviance", "Sigma")
png("hist.mean.inv.wb.png", width = 10, height = 6, units = "in", res = 250)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  Weibull_first_column_list4 <- lapply(Weibull_param_list4[[i]], function(matrix) matrix[, 2])
  # Combine the vectors into a single vector
  Weibull_combined_vector4 <- unlist(Weibull_first_column_list4)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(Weibull_combined_vector4) - 1, max(Weibull_combined_vector4) + 1, length = 1000)
  # Create a histogram and store it in the list
  Weibull_hist_list4[[i]] <- hist(Weibull_combined_vector4, main = paste("Posterior Medians of Parameter", parameter.names.sig[i]), 
                          xlab = "Posterior Median", col = plasma(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = Weibull_true_values[i], col = plasma(12)[1], lwd = 2)
}
dev.off()

Weibull_proportions_list_a4 <- lapply(Weibull_param_list4[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_a4))

Weibull_proportions_list_b4 <- lapply(Weibull_param_list4[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_b4))

Weibull_proportions_list_c4 <- lapply(Weibull_param_list4[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_c4))

Weibull_calculate_rmse_a_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_b_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.Topt, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_c_lambda_mean_inv <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda_mean_inv <- sapply(Weibull_param_list4[[1]], Weibull_calculate_rmse_a_lambda_mean_inv)
Weibull_rmse_values_b_lambda_mean_inv <- sapply(Weibull_param_list4[[2]], Weibull_calculate_rmse_b_lambda_mean_inv)
Weibull_rmse_values_c_lambda_mean_inv <- sapply(Weibull_param_list4[[3]], Weibull_calculate_rmse_c_lambda_mean_inv)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
Weibull_tab.mean.inv.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(Weibull_tab.mean.inv.lambda) <- c("a", "Topt", "c")
colnames(Weibull_tab.mean.inv.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
Weibull_tab.mean.inv.lambda[1, 1] <- true.a
Weibull_tab.mean.inv.lambda[2, 1] <- true.Topt
Weibull_tab.mean.inv.lambda[3, 1] <- true.c
Weibull_tab.mean.inv.lambda[1, 2] <- exp(mean(unlist(lapply(Weibull_param_list4[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.mean.inv.lambda[2, 2] <- (mean(unlist(lapply(Weibull_param_list4[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.mean.inv.lambda[3, 2] <- mean(unlist(lapply(Weibull_param_list4[[3]], function(matrix) matrix[, 1])))
Weibull_tab.mean.inv.lambda[1, 3] <- exp(hdi(unlist(lapply(Weibull_param_list4[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.mean.inv.lambda[2, 3] <- (hdi(unlist(lapply(Weibull_param_list4[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.mean.inv.lambda[3, 3] <- hdi(unlist(lapply(Weibull_param_list4[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.mean.inv.lambda[1, 4] <- exp(hdi(unlist(lapply(Weibull_param_list4[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.mean.inv.lambda[2, 4] <- (hdi(unlist(lapply(Weibull_param_list4[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.mean.inv.lambda[3, 4] <- hdi(unlist(lapply(Weibull_param_list4[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.mean.inv.lambda[1, 5] <- mean(unlist(Weibull_proportions_list_a4))
Weibull_tab.mean.inv.lambda[2, 5] <- mean(unlist(Weibull_proportions_list_b4))
Weibull_tab.mean.inv.lambda[3, 5] <- mean(unlist(Weibull_proportions_list_c4))
Weibull_tab.mean.inv.lambda[1, 6] <- mean(Weibull_rmse_values_a_lambda_mean_inv)
Weibull_tab.mean.inv.lambda[2, 6] <- mean(Weibull_rmse_values_b_lambda_mean_inv)
Weibull_tab.mean.inv.lambda[3, 6] <- mean(Weibull_rmse_values_c_lambda_mean_inv)

save(Weibull_param_list4,
     xdf4_W,
     df.new4_W,
     mean.xdf4_W,
     Weibull_hist_list4,
     plot4_mean_W, 
     plot4_med_W,
     plot4_mean.inv_W, 
     plot4_med.inv_W,
     Weibull_proportions_list_a4,
     Weibull_proportions_list_b4,
     Weibull_proportions_list_c4,
     Weibull_rmse_values_a_lambda_mean_inv,
     Weibull_rmse_values_b_lambda_mean_inv,
     Weibull_rmse_values_c_lambda_mean_inv,
     Weibull_tab.mean.inv.lambda,
     rmse.range.mean.inv.wb,
     file="Weibull_mean_inverse.RData")

##########################################Inverse Mean Lambda #############################################

# Vector for RMSEs for each sample of each model
rmse.range.inv.mean.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list5 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_inner_list5 <- vector("list", length = w)
  # Populate the inner list
  for (j in 1:w) {
    # Initialize a 5x4 matrix
    Weibull_matrix_data5 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_inner_list5[[j]] <- Weibull_matrix_data5
  }
  Weibull_param_list5[[i]] <- Weibull_inner_list5
}
# Create JAGS model for each of the 100 datasets
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- Weibull_data.list[[i]] 
  data$inv.trait <- 1/(data$trait)
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$inv.trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    T = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data.raw.w <- Weibull_data.list[[i]]
  data <- df.mean
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  # Model
  sink("Weibull_inv_mean_lambda.txt")
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
    mu[i] <- exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau)
  }
}", file = "Weibull_inv_mean_lambda.txt")
  # Settings
  parameters <- c("la", "Topt", "c", "sig")

  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  inits = vector('list', nc)
  GenInits = function() {
    la <- rnorm(1, 0, 10)
    c <- rexp(1, 0.5)
    Topt <- rnorm(1, 22, 10)
    sig <- rexp(1, 0.5)
    list(
      la = la,             
      c = c,
      Topt = Topt,
      sig = sig
    )
  }
  for(i in 1:nc){
    inits[[i]] = GenInits()
  }
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "Weibull_inv_mean_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    Weibull_param_list5[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    Weibull_param_list5[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    Weibull_param_list5[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    Weibull_param_list5[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    Weibull_param_list5[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    Weibull_param_list5[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    Weibull_param_list5[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    Weibull_param_list5[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    Weibull_param_list5[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    Weibull_param_list5[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    Weibull_param_list5[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    Weibull_param_list5[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    Weibull_param_list5[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    Weibull_param_list5[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    Weibull_param_list5[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    Weibull_param_list5[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    Weibull_param_list5[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    Weibull_param_list5[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    Weibull_param_list5[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    Weibull_param_list5[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse5.wb <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df5.wb <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse5.wb))
  # Apply the quadratic equation for each combination
  rmse.df5.wb$value <- apply(rmse.df5.wb, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse5.wb[k, 3]) * (location - df.rmse5.wb[k, 5])^2 + df.rmse5.wb[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse5.wb), function(k) {
    evaluated_curve <- rmse.df5.wb$value[rmse.df5.wb$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.inv.mean.wb <- c(rmse.range.inv.mean.wb, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new5_W as a function at each of those temperatures
xdf5_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new5_W))

# Apply the quadratic equation for each combination
xdf5_W$value <- apply(xdf5_W, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new5_W[i, 3]) *  (location - df.new5_W[i, 5])^2 + df.new5_W[i, 1]
})
# To keep consistent with exponential data
xdf5_W <- xdf5_W|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf5_W <- xdf5_W |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
mean.xdf5_W$true_curve <- f(xseq)
# Plot mean lifetime response
plot5_mean_W <- ggplot(mean.xdf5_W, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median lifetime response
plot5_med_W <- ggplot(mean.xdf5_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean mortality rate response
plot5_mean.inv_W <- ggplot(mean.xdf5_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median mortality rate response
plot5_med.inv_W <- ggplot(mean.xdf5_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(plot5_mean_W, plot5_med_W, nrow = 1)
gridExtra::grid.arrange(plot5_mean.inv_W, plot5_med.inv_W, nrow = 1)
ggsave("inv.mean.mort.wb.png", plot = plot5_med_W, width = 5, height = 5)
ggsave("inv.mean.life.wb.png", plot = plot5_med.inv_W, width = 5, height = 5)

# Create a list to store histograms of the posterior distribution of each parameter
Weibull_hist_list5 <- vector("list", length = 5)
parameter.names.sig <- c("a", "Topt", "c", "Deviance", "Sigma")
png("hist.inv.mean.wb.png", width = 10, height = 6, units = "in", res = 250)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  Weibull_first_column_list5 <- lapply(Weibull_param_list5[[i]], function(matrix) matrix[, 2])
  # Combine the vectors into a single vector
  Weibull_combined_vector5 <- unlist(Weibull_first_column_list5)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(Weibull_combined_vector5) - 1, max(Weibull_combined_vector5) + 1, length = 1000)
  # Create a histogram and store it in the list
  Weibull_hist_list5[[i]] <- hist(Weibull_combined_vector5, main = paste("Posterior Medians of Parameter", parameter.names.sig[i]), 
                          xlab = "Posterior Median", col = plasma(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = Weibull_true_values[i], col = plasma(12)[1], lwd = 2)
}
dev.off()

Weibull_proportions_list_a5 <- lapply(Weibull_param_list5[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_a5))

Weibull_proportions_list_b5 <- lapply(Weibull_param_list5[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_b5))

Weibull_proportions_list_c5 <- lapply(Weibull_param_list5[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_c5))

Weibull_calculate_rmse_a_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_b_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.Topt, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_c_lambda_inv_mean <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda_inv_mean <- sapply(Weibull_param_list5[[1]], Weibull_calculate_rmse_a_lambda_inv_mean)
Weibull_rmse_values_b_lambda_inv_mean <- sapply(Weibull_param_list5[[2]], Weibull_calculate_rmse_b_lambda_inv_mean)
Weibull_rmse_values_c_lambda_inv_mean <- sapply(Weibull_param_list5[[3]], Weibull_calculate_rmse_c_lambda_inv_mean)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
Weibull_tab.inv.mean.lambda <- matrix(0, nrow = 3, ncol = 6)
row.names(Weibull_tab.inv.mean.lambda) <- c("a", "Topt", "c")
colnames(Weibull_tab.inv.mean.lambda) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
Weibull_tab.inv.mean.lambda[1, 1] <- true.a
Weibull_tab.inv.mean.lambda[2, 1] <- true.Topt
Weibull_tab.inv.mean.lambda[3, 1] <- true.c
Weibull_tab.inv.mean.lambda[1, 2] <- exp(mean(unlist(lapply(Weibull_param_list5[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.inv.mean.lambda[2, 2] <- (mean(unlist(lapply(Weibull_param_list5[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.inv.mean.lambda[3, 2] <- mean(unlist(lapply(Weibull_param_list5[[3]], function(matrix) matrix[, 1])))
Weibull_tab.inv.mean.lambda[1, 3] <- exp(hdi(unlist(lapply(Weibull_param_list5[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.inv.mean.lambda[2, 3] <- (hdi(unlist(lapply(Weibull_param_list5[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.inv.mean.lambda[3, 3] <- hdi(unlist(lapply(Weibull_param_list5[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.inv.mean.lambda[1, 4] <- exp(hdi(unlist(lapply(Weibull_param_list5[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.inv.mean.lambda[2, 4] <- (hdi(unlist(lapply(Weibull_param_list5[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.inv.mean.lambda[3, 4] <- hdi(unlist(lapply(Weibull_param_list5[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.inv.mean.lambda[1, 5] <- mean(unlist(Weibull_proportions_list_a5))
Weibull_tab.inv.mean.lambda[2, 5] <- mean(unlist(Weibull_proportions_list_b5))
Weibull_tab.inv.mean.lambda[3, 5] <- mean(unlist(Weibull_proportions_list_c5))
Weibull_tab.inv.mean.lambda[1, 6] <- mean(Weibull_rmse_values_a_lambda_inv_mean)
Weibull_tab.inv.mean.lambda[2, 6] <- mean(Weibull_rmse_values_b_lambda_inv_mean)
Weibull_tab.inv.mean.lambda[3, 6] <- mean(Weibull_rmse_values_c_lambda_inv_mean)


save(Weibull_param_list5,
     xdf5_W,
     df.new5_W,
     mean.xdf5_W,
     Weibull_hist_list5,
     plot5_mean_W, 
     plot5_med_W,
     plot5_mean.inv_W, 
     plot5_med.inv_W,
     Weibull_proportions_list_a5,
     Weibull_proportions_list_b5,
     Weibull_proportions_list_c5,
     Weibull_rmse_values_a_lambda_inv_mean,
     Weibull_rmse_values_b_lambda_inv_mean,
     Weibull_rmse_values_c_lambda_inv_mean,
     Weibull_tab.inv.mean.lambda,
     rmse.range.inv.mean.wb,
     file="Weibull_inverse_mean.RData")



########################################## Weibull ###########################################

# Vector for RMSEs for each sample of each model
rmse.range.ind.wb.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list6 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_inner_list6 <- vector("list", length = w)
  # Populate the inner list
  for (j in 1:w) {
    # Initialize a 5x4 matrix
    Weibull_matrix_data6 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_inner_list6[[j]] <- Weibull_matrix_data6
  }
  Weibull_param_list6[[i]] <- Weibull_inner_list6
}
# Create JAGS model for each of the 100 datasets
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("Weibull_weibull.txt")
  cat("model{
  #Priors
  la~dnorm(0, 1/10) 
  Topt~dnorm(22, 1/10) 
  c~dexp(0.5) #Has to be positive
  shape~dexp(0.001)
  #Likelihood
  for (i in 1:N.obs){
    trait[i] ~ dgen.gamma(1, 1/lambda[i], shape) #
    lambda[i] <- (med[i])/(log(2))^(1/shape)
    med[i] <- 1/(exp(la) * (temp[i] - Topt)^2 + c)
  }
  }", file = "Weibull_weibull.txt")
  
  # Settings
  parameters <- c("la", "Topt", "c", "shape")

  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  inits = vector('list', nc)
  GenInits = function() {
    la <- rnorm(1, 0, 10)
    c <- rexp(1, 0.5)
    Topt <- rnorm(1, 22, 10)
    shape <- rexp(1, 0.001)
    list(
      la = la,             
      c = c,
      Topt = Topt,
      shape = shape
    )
  }
  for(i in 1:nc){
    inits[[i]] = GenInits()
  }
  
  data.raw.w <- Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data = jag.data, inits = inits, parameters.to.save = parameters, 
                  model.file = "Weibull_weibull.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    Weibull_param_list6[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    Weibull_param_list6[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    Weibull_param_list6[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    Weibull_param_list6[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    Weibull_param_list6[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    Weibull_param_list6[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    Weibull_param_list6[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    Weibull_param_list6[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    Weibull_param_list6[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    Weibull_param_list6[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    Weibull_param_list6[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    Weibull_param_list6[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    Weibull_param_list6[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    Weibull_param_list6[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    Weibull_param_list6[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    Weibull_param_list6[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for shape
    Weibull_param_list6[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    Weibull_param_list6[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    Weibull_param_list6[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    Weibull_param_list6[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse6.wb <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df6.wb <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse6.wb))
  # Apply the quadratic equation for each combination
  rmse.df6.wb$value <- apply(rmse.df6.wb, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse6.wb[k, 3]) * (location - df.rmse6.wb[k, 5])^2 + df.rmse6.wb[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq) 
  rmse_values <- sapply(1:nrow(df.rmse6.wb), function(k) {
    evaluated_curve <- rmse.df6.wb$value[rmse.df6.wb$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.ind.wb.wb <- c(rmse.range.ind.wb.wb, rmse_values)

}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new6_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new6_W as a function at each of those temperatures
xdf6_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new6_W))

xdf6_W$value <- apply(xdf6_W, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(df.new6_W[i, 3]) * (location - df.new6_W[i, 5])^2  + df.new6_W[i, 1]
})
# To keep consistent with exponential data
xdf6_W <- xdf6_W|> mutate(trunc.inv = log(2)/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
mean.xdf6_W <- xdf6_W |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values 
mean.xdf6_W$true_curve <- f(xseq)
# Plot mean lifetime response
plot6_mean_W <- ggplot(mean.xdf6_W, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot median lifetime response
plot6_med_W <- ggplot(mean.xdf6_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal()
# Plot mean mortality rate response
plot6_mean.inv_W <- ggplot(mean.xdf6_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = log(2)/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()
# Plot median mortality rate response
plot6_med.inv_W <- ggplot(mean.xdf6_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = log(2)/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal()

gridExtra::grid.arrange(plot6_mean_W, plot6_med_W, nrow = 1)
gridExtra::grid.arrange(plot6_mean.inv_W, plot6_med.inv_W, nrow = 1)
ggsave("ind.wb.mort.wb.png", plot = plot6_med_W, width = 5, height = 5)
ggsave("ind.wb.life.wb.png", plot = plot6_med.inv_W, width = 5, height = 5)

# Create a list to store histograms of the posterior distribution of each parameter
Weibull_hist_list6 <- vector("list", length = 5)
parameter.names.sh <- c("a", "Topt", "c", "Deviance", "Shape")
png("hist.wb.ind.wb.png", width = 10, height = 6, units = "in", res = 250)
par(mfrow = c(2, 3))
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  # This column is the mean of each param for each chain
  Weibull_first_column_list6 <- lapply(Weibull_param_list6[[i]], function(matrix) matrix[, 2])
  # Combine the vectors into a single vector
  Weibull_combined_vector6 <- unlist(Weibull_first_column_list6)
  # Sequence of x values to overlay prior distribution onto histograms
  x <- seq(min(Weibull_combined_vector6) - 1, max(Weibull_combined_vector6) + 1, length = 1000)
  # Create a histogram and store it in the list
  Weibull_hist_list6[[i]] <- hist(Weibull_combined_vector6, main = paste("Posterior Medians of Parameter", parameter.names.sh[i]), 
                                    xlab = "Posterior Median", col = plasma(12)[11], border = "black", breaks = 10, freq = FALSE)
  abline(v = Weibull_true_values[i], col = plasma(12)[1], lwd = 2)
}
dev.off()

Weibull_proportions_list_a6 <- lapply(Weibull_param_list6[[1]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for a
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(log(true.a) > col3 & log(true.a) < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_a6))

Weibull_proportions_list_b6 <- lapply(Weibull_param_list6[[2]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for b
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_b6))

Weibull_proportions_list_c6 <- lapply(Weibull_param_list6[[3]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for c
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(true.c > col3 & true.c < col4)
  return(proportion)
  })
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_c6))

Weibull_proportions_list_shape6 <- lapply(Weibull_param_list6[[5]], function(matrix) {
  # Extract columns for mean, lower hdi, and upper hdi of posterior samples for shape
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  # Calculate average proportion of coverage for each model
  proportion <- mean(sh.est > col3 & sh.est < col4)
  return(proportion)
})
# Calculate average proportion of coverage for all 100 models
mean(unlist(Weibull_proportions_list_shape6))

Weibull_calculate_rmse_a_wb <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(log(true.a), length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_b_wb <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]
  # True values
  predicted <- rep(true.Topt, length(observed))  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
  }

Weibull_calculate_rmse_c_wb <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(true.c, length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}

Weibull_calculate_rmse_shape_wb <- function(mat) {
  # Get the first column, which is posterior mean of the chain
  observed <- mat[, 1]  
  # True values
  predicted <- rep(sh.est, length(observed)) 
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_rmse_values_a_wb <- sapply(Weibull_param_list6[[1]], Weibull_calculate_rmse_a_wb)
Weibull_rmse_values_b_wb <- sapply(Weibull_param_list6[[2]], Weibull_calculate_rmse_b_wb)
Weibull_rmse_values_c_wb <- sapply(Weibull_param_list6[[3]], Weibull_calculate_rmse_c_wb)
Weibull_rmse_values_shape_wb <- sapply(Weibull_param_list6[[5]], Weibull_calculate_rmse_shape_wb)

# Create a table with true value, posterior mean,
# lower and upper hdi of posterior mean, coverage, and RMSE
Weibull_tab.wb <- matrix(0, nrow = 4, ncol = 6)
row.names(Weibull_tab.wb) <- c("a", "Topt", "c", "shape")
colnames(Weibull_tab.wb) <- c("True Value", "Mean", "Lower HDI of Mean", "Upper HDI of Mean", "Coverage", "RMSE")
Weibull_tab.wb[1, 1] <- true.a
Weibull_tab.wb[2, 1] <- true.Topt
Weibull_tab.wb[3, 1] <- true.c
Weibull_tab.wb[4, 1] <- sh.est
Weibull_tab.wb[1, 2] <- exp(mean(unlist(lapply(Weibull_param_list6[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.wb[2, 2] <- (mean(unlist(lapply(Weibull_param_list6[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.wb[3, 2] <- mean(unlist(lapply(Weibull_param_list6[[3]], function(matrix) matrix[, 1])))
Weibull_tab.wb[4, 2] <- mean(unlist(lapply(Weibull_param_list6[[5]], function(matrix) matrix[, 1])))
Weibull_tab.wb[1, 3] <- exp(hdi(unlist(lapply(Weibull_param_list6[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.wb[2, 3] <- (hdi(unlist(lapply(Weibull_param_list6[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.wb[3, 3] <- hdi(unlist(lapply(Weibull_param_list6[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.wb[4, 3] <- hdi(unlist(lapply(Weibull_param_list6[[5]], function(matrix) matrix[, 1])))[1]
Weibull_tab.wb[1, 4] <- exp(hdi(unlist(lapply(Weibull_param_list6[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.wb[2, 4] <- (hdi(unlist(lapply(Weibull_param_list6[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.wb[3, 4] <- hdi(unlist(lapply(Weibull_param_list6[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.wb[4, 4] <- hdi(unlist(lapply(Weibull_param_list6[[5]], function(matrix) matrix[, 1])))[2]
Weibull_tab.wb[1, 5] <- mean(unlist(Weibull_proportions_list_a6))
Weibull_tab.wb[2, 5] <- mean(unlist(Weibull_proportions_list_b6))
Weibull_tab.wb[3, 5] <- mean(unlist(Weibull_proportions_list_c6))
Weibull_tab.wb[4, 5] <- mean(unlist(Weibull_proportions_list_shape6))
Weibull_tab.wb[1, 6] <- mean(Weibull_rmse_values_a_wb)
Weibull_tab.wb[2, 6] <- mean(Weibull_rmse_values_b_wb)
Weibull_tab.wb[3, 6] <- mean(Weibull_rmse_values_c_wb)
Weibull_tab.wb[4, 6] <- mean(Weibull_rmse_values_shape_wb)

save(Weibull_param_list6,
     xdf6_W,
     df.new6_W,
     mean.xdf6_W,
     Weibull_hist_list6,
     plot6_mean_W, 
     plot6_med_W,
     plot6_mean.inv_W, 
     plot6_med.inv_W,
     Weibull_proportions_list_a6,
     Weibull_proportions_list_b6,
     Weibull_proportions_list_c6,
     Weibull_rmse_values_a_wb,
     Weibull_rmse_values_b_wb,
     Weibull_rmse_values_c_wb,
     Weibull_tab.wb,
     rmse.range.ind.wb.wb,
     file="Weibull_weibull_data.RData")

