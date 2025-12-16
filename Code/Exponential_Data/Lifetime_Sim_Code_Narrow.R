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
# Narrow temperature interval
N_data.list <- list()
set.seed(5212)
for(i in 1:100){
  T1 <- rexp(n = 100, f(13))
  T2 <- rexp(n = 100, f(17))
  T3 <- rexp(n = 100, f(21))
  T4 <- rexp(n = 100, f(25))
  T5 <- rexp(n = 100, f(29))
  N_data.list[[i]] <- data.frame(
    Temp = c(rep(13, 100), rep(17, 100), rep(21, 100), rep(25, 100), rep(29, 100)), 
    trait = c(T1, T2, T3, T4, T5)
  )
}
for(i in 1:100){
  N_data.list[[i]]$trait <- ifelse(N_data.list[[i]]$trait < 1, 1, N_data.list[[i]]$trait)
}

b <- length(N_data.list)

# True value vector with a and b in log space
# Spot for sigma and deviance
true_values <- c(log(true.a), true.Topt, true.c, NA, NA)
xseq <- seq(from = 5, to = 37, length.out = 250)

########################################## Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data1.n <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.ind.n <- c()
# Create list of matrices to store parameter summary statistics
N_param_list1 <- vector("list", length = 4)
# Create dataframe to store uncertainty results
n.exp.ind.uncertainty <- data.frame()
n.exp.ind.whisk <- data.frame()
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
  la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10) 
  c ~ dexp(0.5) # Has to be positive
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dexp(mu[i])
  }
  }", file = "n_lambda.txt")
  
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
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  
  N_data.raw <- N_data.list[[i]]
  data <- N_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$Temp
  
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
  xseq <- seq(from = 5, to = 37, length.out = 250)
  for (j in 1:nc) {
    # Compute HDI bounds for parameters of interest (Topt and inv of max lifetime)
    hdi_bounds.Topt <- hdi(samp[[j]][, "Topt"])
    hdi_bounds.life <- hdi(samp[[j]][, "c"])
    # Build a row
    new_row <- data.frame(
      Method = "Individual Data",
      Model = i,
      Chain = j,
      # Compute lower and upper HDI bounds of samples, as well as length of interval
      Lower.Topt = hdi_bounds.Topt[1],
      Upper.Topt = hdi_bounds.Topt[2],
      Length.Topt = hdi_bounds.Topt[2] - hdi_bounds.Topt[1],
      Lower.life = hdi_bounds.life[1],
      Upper.life = hdi_bounds.life[2],
      Length.life = hdi_bounds.life[2] - hdi_bounds.life[1]
    )
    # Append to results dataframe
    n.exp.ind.uncertainty <- rbind(n.exp.ind.uncertainty, new_row)
  }
  # Put chains together for extreme temps
  all_chains <- do.call(rbind, samp)
  df <- as.data.frame(all_chains)
  x_values <- numeric()
  for (row in 1:nrow(df)) {
    a <- exp(df[row, "la"])         
    b <- -2 * df[row, "Topt"] * a
    c <- df[row, "c"] + a * (df[row, "Topt"])^2
    
    discrim <- b^2 - 4 * a * (c - 1)
    if (discrim >= 0) {
      x1 <- (-b + sqrt(discrim)) / (2 * a)
      x2 <- (-b - sqrt(discrim)) / (2 * a)
      x_values <- c(x_values, x1, x2)
    }
  }
  
  # Separate into Low and High critical temperatures
  temp_df <- data.frame(
    Low = x_values[c(FALSE, TRUE)],
    High = x_values[c(TRUE, FALSE)]
  )
  
  # Compute boxplot whisker stats across all datasets
  q1Low <- quantile(temp_df$Low, 0.25)
  q3Low <- quantile(temp_df$Low, 0.75)
  iqrLow <- q3Low - q1Low
  minLow <- q1Low - 1.5 * iqrLow
  maxLow <- q3Low + 1.5 * iqrLow
  rangeLow <- maxLow - minLow
  
  q1High <- quantile(temp_df$High, 0.25)
  q3High <- quantile(temp_df$High, 0.75)
  iqrHigh <- q3High - q1High
  minHigh <- q1High - 1.5 * iqrHigh
  maxHigh <- q3High + 1.5 * iqrHigh
  rangeHigh <- maxHigh - minHigh
  
  # Save results for this model
  new_row <- data.frame(
    Method = "Individual Data", 
    Model = i,
    Range.Low = rangeLow,
    Range.High = rangeHigh
  )
  
  n.exp.ind.whisk <- rbind(n.exp.ind.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data1.n) {
  #   pdf(paste0("traceplot.exp.ind.n.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    N_param_list1[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list1[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list1[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list1[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    N_param_list1[[2]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list1[[2]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list1[[2]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list1[[2]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list1[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list1[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list1[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list1[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list1[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list1[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list1[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list1[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
  }
  df.rmse1.n <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df1.n <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse1.n))
  # Apply the quadratic equation for each combination
  rmse.df1.n$value <- apply(rmse.df1.n, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse1.n[k, 3]) * (location - df.rmse1.n[k, 4])^2 + df.rmse1.n[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse1.n), function(k) {
    evaluated_curve <- rmse.df1.n$value[rmse.df1.n$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.ind.n <- c(rmse.range.ind.n, rmse_values)
  
}

# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new1_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = -10, to = 60, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new1_N as a function at each of those temperatures
xdf1_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new1_N))

# Apply the quadratic equation for each combination
xdf1_N$value <- apply(xdf1_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new1_N[i, 3]) * (location - df.new1_N[i, 4])^2  + df.new1_N[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf1_N <- xdf1_N|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf1 <- xdf1_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
N_mean.xdf1$true_curve <- f(xseq)
# Plot the mean mortality rate response
N_plot1_mean <- ggplot(N_mean.xdf1, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the median mortality rate response
N_plot1_med <- ggplot(N_mean.xdf1, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the mean lifetime response
N_plot1_mean.inv <- ggplot(N_mean.xdf1, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the median lifetime response
N_plot1_med.inv <- ggplot(N_mean.xdf1, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

N_plot1_ex.temp <- ggplot(N_mean.xdf1, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") + 
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") + 
  #coord_cartesian(ylim = c(0, 40)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("n.ind.ex.plot.png", plot = N_plot1_ex.temp, width = 5, height = 5)


gridExtra::grid.arrange(N_plot1_mean, N_plot1_med, nrow = 1)
gridExtra::grid.arrange(N_plot1_mean.inv, N_plot1_med.inv, nrow = 1)
ggsave("n.ind.mort.exp.png", plot = N_plot1_med, width = 5, height = 5)
ggsave("n.ind.life.exp.png", plot = N_plot1_med.inv, width = 5, height = 5)


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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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

save(N_param_list1,
     xdf1_N,
     df.new1_N,
     N_mean.xdf1,
     N_plot1_mean, 
     N_plot1_med,
     N_plot1_mean.inv, 
     N_plot1_med.inv,
     N_proportions_list_a1,
     N_proportions_list_b1,
     N_proportions_list_c1,
     rmse.range.ind.n,
     file="individual_exp_narrow.RData")


########################################## Mean Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data2.n <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.mean.n <- c()
# Create list of matrices to store parameter summary statistics
N_param_list2 <- vector("list", length = 4)
# Create dataframe to store uncertainty results
n.exp.mean.uncertainty <- data.frame()
n.exp.mean.whisk <- data.frame()
for (i in 1:5) {
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
  group_means <- tapply(N_data.list[[i]]$trait, N_data.list[[i]]$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    Temp = c(rep(13, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(29, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  # Model
  sink("n_mean_lambda.txt")
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
    trait[i] ~ dnorm(1/mu[i], tau) T(0,)
  }
}", file = "n_mean_lambda.txt")
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
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  data <- df.mean 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$Temp
  
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
  xseq <- seq(from = 5, to = 37, length.out = 250)
  for (j in 1:nc) {
    # Compute HDI bounds for parameters of interest (Topt and inv of max lifetime)
    hdi_bounds.Topt <- hdi(samp[[j]][, "Topt"])
    hdi_bounds.life <- hdi(samp[[j]][, "c"])
    # Build a row
    new_row <- data.frame(
      Method = "Mean Data",
      Model = i,
      Chain = j,
      # Compute lower and upper HDI bounds of samples, as well as length of interval
      Lower.Topt = hdi_bounds.Topt[1],
      Upper.Topt = hdi_bounds.Topt[2],
      Length.Topt = hdi_bounds.Topt[2] - hdi_bounds.Topt[1],
      Lower.life = hdi_bounds.life[1],
      Upper.life = hdi_bounds.life[2],
      Length.life = hdi_bounds.life[2] - hdi_bounds.life[1]
    )
    # Append to results dataframe
    n.exp.mean.uncertainty <- rbind(n.exp.mean.uncertainty, new_row)
  }
  # Put chains together for extreme temps
  all_chains <- do.call(rbind, samp)
  df <- as.data.frame(all_chains)
  x_values <- numeric()
  for (row in 1:nrow(df)) {
    a <- exp(df[row, "la"])         
    b <- -2 * df[row, "Topt"] * a
    c <- df[row, "c"] + a * (df[row, "Topt"])^2
    
    discrim <- b^2 - 4 * a * (c - 1)
    if (discrim >= 0) {
      x1 <- (-b + sqrt(discrim)) / (2 * a)
      x2 <- (-b - sqrt(discrim)) / (2 * a)
      x_values <- c(x_values, x1, x2)
    }
  }
  
  # Separate into Low and High critical temperatures
  temp_df <- data.frame(
    Low = x_values[c(FALSE, TRUE)],
    High = x_values[c(TRUE, FALSE)]
  )
  
  # Compute boxplot whisker stats across all datasets
  q1Low <- quantile(temp_df$Low, 0.25)
  q3Low <- quantile(temp_df$Low, 0.75)
  iqrLow <- q3Low - q1Low
  minLow <- q1Low - 1.5 * iqrLow
  maxLow <- q3Low + 1.5 * iqrLow
  rangeLow <- maxLow - minLow
  
  q1High <- quantile(temp_df$High, 0.25)
  q3High <- quantile(temp_df$High, 0.75)
  iqrHigh <- q3High - q1High
  minHigh <- q1High - 1.5 * iqrHigh
  maxHigh <- q3High + 1.5 * iqrHigh
  rangeHigh <- maxHigh - minHigh
  
  # Save results for this model
  new_row <- data.frame(
    Method = "Mean Data", 
    Model = i,
    Range.Low = rangeLow,
    Range.High = rangeHigh
  )
  
  n.exp.mean.whisk <- rbind(n.exp.mean.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data2.n) {
  #   pdf(paste0("traceplot.exp.mean.n.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    N_param_list2[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list2[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list2[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list2[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    N_param_list2[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list2[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list2[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list2[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list2[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list2[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list2[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list2[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list2[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list2[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list2[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list2[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    N_param_list2[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list2[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list2[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list2[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse2.n <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df2.n <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse2.n))
  # Apply the quadratic equation for each combination
  rmse.df2.n$value <- apply(rmse.df2.n, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse2.n[k, 3]) * (location - df.rmse2.n[k, 5])^2 + df.rmse2.n[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse2.n), function(k) {
    evaluated_curve <- rmse.df2.n$value[rmse.df2.n$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.mean.n <- c(rmse.range.mean.n, rmse_values)
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new2_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = 5, to = 37, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new2_N as a function at each of those temperatures
xdf2_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new2_N))

# Apply the quadratic equation for each combination
xdf2_N$value <- apply(xdf2_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new2_N[i, 3]) * (location - df.new2_N[i, 5])^2  + df.new2_N[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf2_N <- xdf2_N|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf2 <- xdf2_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
N_mean.xdf2$true_curve <- f(xseq)
# Plot mean mortality rate response
N_plot2_mean <- ggplot(N_mean.xdf2, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
N_plot2_med <- ggplot(N_mean.xdf2, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean lifetime response
N_plot2_mean.inv <- ggplot(N_mean.xdf2, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
N_plot2_med.inv <- ggplot(N_mean.xdf2, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

N_plot2_ex.temp <- ggplot(N_mean.xdf2, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") + 
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") + 
  #coord_cartesian(ylim = c(0, 40)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("n.mean.ex.plot.png", plot = N_plot2_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(N_plot2_mean, N_plot2_med, nrow = 1)
gridExtra::grid.arrange(N_plot2_mean.inv, N_plot2_med.inv, nrow = 1)
ggsave("n.mean.mort.exp.png", plot = N_plot2_med, width = 5, height = 5)
ggsave("n.mean.life.exp.png", plot = N_plot2_med.inv, width = 5, height = 5)


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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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

save(N_param_list2,
     xdf2_N,
     df.new2_N,
     N_mean.xdf2,
     N_plot2_mean, 
     N_plot2_med,
     N_plot2_mean.inv, 
     N_plot2_med.inv,
     N_proportions_list_a2,
     N_proportions_list_b2,
     N_proportions_list_c2,
     rmse.range.mean.n,
     file="mean_exp_narrow.RData")

########################################## Inverse Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data3.n <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.inv.n <- c()
# Create list of matrices to store parameter summary statistics
N_param_list3 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
n.exp.inv.uncertainty <- data.frame()
n.exp.inv.whisk <- data.frame()
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
  la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10)
  c ~ dgamma(10,1) # Has to be positive
  sig ~ dexp(0.5)
  sig2 <- sig^2
  tau <- 1/sig2
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <-  exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau) T(0,)
  }
}", file = "n_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "Topt", "c", "sig")

  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  inits = vector('list', nc)
  GenInits = function() {
    la <- rnorm(1, 0, 10)
    c <- rgamma(1, 10, 1)
    Topt <- rnorm(1, 22, 10)
    sig <- rexp(1, 0.5)
    list(
      la = la,             
      c = c,
      Topt = Topt,
      sig = sig
    )
  }
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  
  N_data.raw <- N_data.list[[i]]
  data <- N_data.list[[i]] 
  trait <- ifelse(1/(data$trait) > 1, 1, 1/(data$trait))
  N.obs <- length(trait)
  temp <- data$Temp
  
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
  xseq <- seq(from = 5, to = 37, length.out = 250)
  for (j in 1:nc) {
    # Compute HDI bounds for parameters of interest (Topt and inv of max lifetime)
    hdi_bounds.Topt <- hdi(samp[[j]][, "Topt"])
    hdi_bounds.life <- hdi(samp[[j]][, "c"])
    # Build a row
    new_row <- data.frame(
      Method = "Individual Inverse Data",
      Model = i,
      Chain = j,
      # Compute lower and upper HDI bounds of samples, as well as length of interval
      Lower.Topt = hdi_bounds.Topt[1],
      Upper.Topt = hdi_bounds.Topt[2],
      Length.Topt = hdi_bounds.Topt[2] - hdi_bounds.Topt[1],
      Lower.life = hdi_bounds.life[1],
      Upper.life = hdi_bounds.life[2],
      Length.life = hdi_bounds.life[2] - hdi_bounds.life[1]
    )
    # Append to results dataframe
    n.exp.inv.uncertainty <- rbind(n.exp.inv.uncertainty, new_row)
  }
  # Put chains together for extreme temps
  all_chains <- do.call(rbind, samp)
  df <- as.data.frame(all_chains)
  x_values <- numeric()
  for (row in 1:nrow(df)) {
    a <- exp(df[row, "la"])         
    b <- -2 * df[row, "Topt"] * a
    c <- df[row, "c"] + a * (df[row, "Topt"])^2
    
    discrim <- b^2 - 4 * a * (c - 1)
    if (discrim >= 0) {
      x1 <- (-b + sqrt(discrim)) / (2 * a)
      x2 <- (-b - sqrt(discrim)) / (2 * a)
      x_values <- c(x_values, x1, x2)
    }
  }
  
  # Separate into Low and High critical temperatures
  temp_df <- data.frame(
    Low = x_values[c(FALSE, TRUE)],
    High = x_values[c(TRUE, FALSE)]
  )
  
  # Compute boxplot whisker stats across all datasets
  q1Low <- quantile(temp_df$Low, 0.25)
  q3Low <- quantile(temp_df$Low, 0.75)
  iqrLow <- q3Low - q1Low
  minLow <- q1Low - 1.5 * iqrLow
  maxLow <- q3Low + 1.5 * iqrLow
  rangeLow <- maxLow - minLow
  
  q1High <- quantile(temp_df$High, 0.25)
  q3High <- quantile(temp_df$High, 0.75)
  iqrHigh <- q3High - q1High
  minHigh <- q1High - 1.5 * iqrHigh
  maxHigh <- q3High + 1.5 * iqrHigh
  rangeHigh <- maxHigh - minHigh
  
  # Save results for this model
  new_row <- data.frame(
    Method = "Individual Inverse Data", 
    Model = i,
    Range.Low = rangeLow,
    Range.High = rangeHigh
  )
  
  n.exp.inv.whisk <- rbind(n.exp.inv.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data3.n) {
  #   pdf(paste0("traceplot.exp.inv.n.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    N_param_list3[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list3[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list3[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list3[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    N_param_list3[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list3[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list3[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list3[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list3[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list3[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list3[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list3[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list3[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list3[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list3[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list3[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    N_param_list3[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list3[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list3[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list3[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse3.n <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df3.n <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse3.n))
  # Apply the quadratic equation for each combination
  rmse.df3.n$value <- apply(rmse.df3.n, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse3.n[k, 3]) * (location - df.rmse3.n[k, 5])^2 + df.rmse3.n[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse3.n), function(k) {
    evaluated_curve <- rmse.df3.n$value[rmse.df3.n$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.inv.n <- c(rmse.range.inv.n, rmse_values)
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new3_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = -10, to = 60, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new3_N as a function at each of those temperatures
xdf3_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new3_N))

# Apply the quadratic equation for each combination
xdf3_N$value <- apply(xdf3_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new3_N[i, 3]) * (location - df.new3_N[i, 5])^2  + df.new3_N[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf3_N <- xdf3_N|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf3 <- xdf3_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
N_mean.xdf3$true_curve <- f(xseq)
# Plot mean mortality rate response
N_plot3_mean <- ggplot(N_mean.xdf3, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
N_plot3_med <- ggplot(N_mean.xdf3, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = 1/trait), data = data) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean lifetime response
N_plot3_mean.inv <- ggplot(N_mean.xdf3, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
N_plot3_med.inv <- ggplot(N_mean.xdf3, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

N_plot3_ex.temp <- ggplot(N_mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") + 
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") + 
  #coord_cartesian(ylim = c(0, 40)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("n.inv.ex.plot.png", plot = N_plot3_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(N_plot3_mean, N_plot3_med, nrow = 1)
gridExtra::grid.arrange(N_plot3_mean.inv, N_plot3_med.inv, nrow = 1)
ggsave("n.inv.mort.exp.png", plot = N_plot3_med, width = 5, height = 5)
ggsave("n.inv.life.exp.png", plot = N_plot3_med.inv, width = 5, height = 5)


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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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

save(N_param_list3,
     xdf3_N,
     df.new3_N,
     N_mean.xdf3,
     N_plot3_mean, 
     N_plot3_med,
     N_plot3_mean.inv, 
     N_plot3_med.inv,
     N_proportions_list_a3,
     N_proportions_list_b3,
     N_proportions_list_c3,
     rmse.range.inv.n,
     file="inverse_exp_narrow.RData")

########################################## Inverse of Mean Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data4.n <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.inv.mean.n <- c()
# Create list of matrices to store parameter summary statistics
N_param_list4 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
n.exp.iom.uncertainty <- data.frame()
n.exp.iom.whisk <- data.frame()
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
  group_means <- tapply(data$trait, data$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    Temp = c(rep(13, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(29, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$Temp
  # Model
  sink("n_inv_mean_lambda.txt")
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
}", file = "n_inv_mean_lambda.txt")
  
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
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  
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
  xseq <- seq(from = 5, to = 37, length.out = 250)
  for (j in 1:nc) {
    # Compute HDI bounds for parameters of interest (Topt and inv of max lifetime)
    hdi_bounds.Topt <- hdi(samp[[j]][, "Topt"])
    hdi_bounds.life <- hdi(samp[[j]][, "c"])
    # Build a row
    new_row <- data.frame(
      Method = "Inverse of Mean Data",
      Model = i,
      Chain = j,
      # Compute lower and upper HDI bounds of samples, as well as length of interval
      Lower.Topt = hdi_bounds.Topt[1],
      Upper.Topt = hdi_bounds.Topt[2],
      Length.Topt = hdi_bounds.Topt[2] - hdi_bounds.Topt[1],
      Lower.life = hdi_bounds.life[1],
      Upper.life = hdi_bounds.life[2],
      Length.life = hdi_bounds.life[2] - hdi_bounds.life[1]
    )
    # Append to results dataframe
    n.exp.iom.uncertainty <- rbind(n.exp.iom.uncertainty, new_row)
  }
  # Put chains together for extreme temps
  all_chains <- do.call(rbind, samp)
  df <- as.data.frame(all_chains)
  x_values <- numeric()
  for (row in 1:nrow(df)) {
    a <- exp(df[row, "la"])         
    b <- -2 * df[row, "Topt"] * a
    c <- df[row, "c"] + a * (df[row, "Topt"])^2
    
    discrim <- b^2 - 4 * a * (c - 1)
    if (discrim >= 0) {
      x1 <- (-b + sqrt(discrim)) / (2 * a)
      x2 <- (-b - sqrt(discrim)) / (2 * a)
      x_values <- c(x_values, x1, x2)
    }
  }
  
  # Separate into Low and High critical temperatures
  temp_df <- data.frame(
    Low = x_values[c(FALSE, TRUE)],
    High = x_values[c(TRUE, FALSE)]
  )
  
  # Compute boxplot whisker stats across all datasets
  q1Low <- quantile(temp_df$Low, 0.25)
  q3Low <- quantile(temp_df$Low, 0.75)
  iqrLow <- q3Low - q1Low
  minLow <- q1Low - 1.5 * iqrLow
  maxLow <- q3Low + 1.5 * iqrLow
  rangeLow <- maxLow - minLow
  
  q1High <- quantile(temp_df$High, 0.25)
  q3High <- quantile(temp_df$High, 0.75)
  iqrHigh <- q3High - q1High
  minHigh <- q1High - 1.5 * iqrHigh
  maxHigh <- q3High + 1.5 * iqrHigh
  rangeHigh <- maxHigh - minHigh
  
  # Save results for this model
  new_row <- data.frame(
    Method = "Inverse of Mean Data",  
    Model = i,
    Range.Low = rangeLow,
    Range.High = rangeHigh
  )
  
  n.exp.iom.whisk <- rbind(n.exp.iom.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data4.n) {
  #   pdf(paste0("traceplot.exp.inv.mean.n.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    N_param_list4[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list4[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list4[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list4[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    N_param_list4[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list4[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list4[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list4[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list4[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list4[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list4[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list4[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list4[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list4[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list4[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list4[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    N_param_list4[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list4[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list4[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list4[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse4.n <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df4.n <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse4.n))
  # Apply the quadratic equation for each combination
  rmse.df4.n$value <- apply(rmse.df4.n, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse4.n[k, 3]) * (location - df.rmse4.n[k, 5])^2 + df.rmse4.n[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse4.n), function(k) {
    evaluated_curve <- rmse.df4.n$value[rmse.df4.n$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.inv.mean.n <- c(rmse.range.inv.mean.n, rmse_values)
  
}

# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = -10, to = 60, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new4_N as a function at each of those temperatures
xdf4_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new4_N))

# Apply the quadratic equation for each combination
xdf4_N$value <- apply(xdf4_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new4_N[i, 3]) * (location - df.new4_N[i, 5])^2  + df.new4_N[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf4_N <- xdf4_N|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf4 <- xdf4_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
N_mean.xdf4$true_curve <- f(xseq)
# Plot mean mortality rate response
N_plot4_mean <- ggplot(N_mean.xdf4, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
N_plot4_med <- ggplot(N_mean.xdf4, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = 1/trait), data = data) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean lifetime response
N_plot4_mean.inv <- ggplot(N_mean.xdf4, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
N_plot4_med.inv <- ggplot(N_mean.xdf4, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )


N_plot4_ex.temp <- ggplot(N_mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") + 
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") + 
  coord_cartesian(ylim = c(0, 2), xlim = c(51, 53)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# CHANGED TO ZOOM
ggsave("n.inv.mean.ex.plot.zoom.png", plot = N_plot4_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(N_plot4_mean, N_plot4_med, nrow = 1)
gridExtra::grid.arrange(N_plot4_mean.inv, N_plot4_med.inv, nrow = 1)
ggsave("n.inv.mean.mort.exp.png", plot = N_plot4_med, width = 5, height = 5)
ggsave("n.inv.mean.life.exp.png", plot = N_plot4_med.inv, width = 5, height = 5)

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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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

save(N_param_list4,
     xdf4_N,
     df.new4_N,
     N_mean.xdf4,
     N_plot4_mean, 
     N_plot4_med,
     N_plot4_inv.mean, 
     N_plot4_med.inv,
     N_proportions_list_a4,
     N_proportions_list_b4,
     N_proportions_list_c4,
     rmse.range.inv.mean.n,
     file="inverse_mean_exp_narrow.RData")

########################################## Mean of Inverse Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data5.n <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.mean.inv.n <- c()
# Create list of matrices to store parameter summary statistics
N_param_list5 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
n.exp.moi.uncertainty <- data.frame()
n.exp.moi.whisk <- data.frame()
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
  data$t.trait <- ifelse(data$trait < 1, 1, data$trait)
  data$inv.trait <- 1/(data$t.trait)
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data$inv.trait, data$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    Temp = c(rep(13, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(29, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- ifelse(data$trait > 1, 1, data$trait)
  N.obs <- length(trait)
  temp <- data$Temp
  # Model
  sink("n_mean_inv_lambda.txt")
  cat("model{
  # Priors
   la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10) 
  c ~ dgamma(10, 1) # Has to be positive
  sig ~ dexp(0.5)
  sig2 <- sig^2
  tau <- 1/sig2
  # Likelihood
  for (i in 1:N.obs){
    mu[i] <-  exp(la) * (temp[i] - Topt)^2 + c
    trait[i] ~ dnorm(mu[i], tau) T(0,)
  }
}", file = "n_mean_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "Topt", "c", "sig")

  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  inits = vector('list', nc)
  GenInits = function() {
    la <- rnorm(1, 0, 10)
    c <- rgamma(1, 10, 1)
    Topt <- rnorm(1, 22, 10)
    sig <- rexp(1, 0.5)
    list(
      la = la,             
      c = c,
      Topt = Topt,
      sig = sig
    )
  }
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  
  N_data.raw <- N_data.list[[i]]
  
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
  xseq <- seq(from = 5, to = 37, length.out = 250)
  for (j in 1:nc) {
    # Compute HDI bounds for parameters of interest (Topt and inv of max lifetime)
    hdi_bounds.Topt <- hdi(samp[[j]][, "Topt"])
    hdi_bounds.life <- hdi(samp[[j]][, "c"])
    # Build a row
    new_row <- data.frame(
      Method = "Mean of Inverse Data",
      Model = i,
      Chain = j,
      # Compute lower and upper HDI bounds of samples, as well as length of interval
      Lower.Topt = hdi_bounds.Topt[1],
      Upper.Topt = hdi_bounds.Topt[2],
      Length.Topt = hdi_bounds.Topt[2] - hdi_bounds.Topt[1],
      Lower.life = hdi_bounds.life[1],
      Upper.life = hdi_bounds.life[2],
      Length.life = hdi_bounds.life[2] - hdi_bounds.life[1]
    )
    # Append to results dataframe
    n.exp.moi.uncertainty <- rbind(n.exp.moi.uncertainty, new_row)
  }
  # Put chains together for extreme temps
  all_chains <- do.call(rbind, samp)
  df <- as.data.frame(all_chains)
  x_values <- numeric()
  for (row in 1:nrow(df)) {
    a <- exp(df[row, "la"])         
    b <- -2 * df[row, "Topt"] * a
    c <- df[row, "c"] + a * (df[row, "Topt"])^2
    
    discrim <- b^2 - 4 * a * (c - 1)
    if (discrim >= 0) {
      x1 <- (-b + sqrt(discrim)) / (2 * a)
      x2 <- (-b - sqrt(discrim)) / (2 * a)
      x_values <- c(x_values, x1, x2)
    }
  }
  
  # Separate into Low and High critical temperatures
  temp_df <- data.frame(
    Low = x_values[c(FALSE, TRUE)],
    High = x_values[c(TRUE, FALSE)]
  )
  
  # Compute boxplot whisker stats across all datasets
  q1Low <- quantile(temp_df$Low, 0.25)
  q3Low <- quantile(temp_df$Low, 0.75)
  iqrLow <- q3Low - q1Low
  minLow <- q1Low - 1.5 * iqrLow
  maxLow <- q3Low + 1.5 * iqrLow
  rangeLow <- maxLow - minLow
  
  q1High <- quantile(temp_df$High, 0.25)
  q3High <- quantile(temp_df$High, 0.75)
  iqrHigh <- q3High - q1High
  minHigh <- q1High - 1.5 * iqrHigh
  maxHigh <- q3High + 1.5 * iqrHigh
  rangeHigh <- maxHigh - minHigh
  
  # Save results for this model
  new_row <- data.frame(
    Method = "Mean of Inverse Data", 
    Model = i,
    Range.Low = rangeLow,
    Range.High = rangeHigh
  )
  
  n.exp.moi.whisk <- rbind(n.exp.moi.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data5.n) {
  #   pdf(paste0("traceplot.exp.mean.inv.n.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    N_param_list5[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    N_param_list5[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    N_param_list5[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    N_param_list5[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt
    N_param_list5[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    N_param_list5[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    N_param_list5[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    N_param_list5[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
    # For each chain, store mean, median, and hdi bounds for c
    N_param_list5[[3]][[i]][j, 1] <- mean(samp[[j]][, 1])
    N_param_list5[[3]][[i]][j, 2] <- median(samp[[j]][, 1])
    N_param_list5[[3]][[i]][j, 3] <- hdi(samp[[j]][, 1])[1]
    N_param_list5[[3]][[i]][j, 4] <- hdi(samp[[j]][, 1])[2]
    # For each chain, store mean, median, and hdi bounds for deviance
    N_param_list5[[4]][[i]][j, 1] <- mean(samp[[j]][, 2])
    N_param_list5[[4]][[i]][j, 2] <- median(samp[[j]][, 2])
    N_param_list5[[4]][[i]][j, 3] <- hdi(samp[[j]][, 2])[1]
    N_param_list5[[4]][[i]][j, 4] <- hdi(samp[[j]][, 2])[2]
    # For each chain, store mean, median, and hdi bounds for sig
    N_param_list5[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    N_param_list5[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    N_param_list5[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    N_param_list5[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse5.n <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df5.n <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse5.n))
  # Apply the quadratic equation for each combination
  rmse.df5.n$value <- apply(rmse.df5.n, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse5.n[k, 3]) * (location - df.rmse5.n[k, 5])^2 + df.rmse5.n[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse5.n), function(k) {
    evaluated_curve <- rmse.df5.n$value[rmse.df5.n$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.mean.inv.n <- c(rmse.range.mean.inv.n, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5_N = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 13 to 29
xseq <- seq(from = -10, to = 60, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new5_N as a function at each of those temperatures
xdf5_N <- expand.grid(point = xseq, iteration = 1:nrow(df.new5_N))

# Apply the quadratic equation for each combination
xdf5_N$value <- apply(xdf5_N, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  (exp(df.new5_N[i, 3]) * (location - df.new5_N[i, 5])^2  + df.new5_N[i, 1])
})
# Apply truncation so that if the evaluated function gives a value smaller than epsilon,
# it gets truncated to epsilon
xdf5_N <- xdf5_N|> mutate(trunc.inv = 1/value)
# At each temperature point, summarize to get the mean, median, and hdi bounds
# of the evaluated function
N_mean.xdf5 <- xdf5_N |> group_by(point) |> summarize(avg.value.inv = mean(trunc.inv), 
                                                    avg.value = mean(value), 
                                                    med.value.inv = median(trunc.inv), 
                                                    med.value = median(value), 
                                                    lower.hdi.inv = hdi(trunc.inv)[1], 
                                                    upper.hdi.inv = hdi(trunc.inv)[2], 
                                                    lower.hdi = hdi(value)[1], 
                                                    upper.hdi = hdi(value)[2])
# Creating columns of the true curve values and the true inverted curve values
N_mean.xdf5$true_curve <- f(xseq)
# Plot mean mortality rate response
N_plot5_mean <- ggplot(N_mean.xdf5, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  #geom_point(aes(x = T, y = 1/trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
N_plot5_med <- ggplot(N_mean.xdf5, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean lifetime response
N_plot5_mean.inv <- ggplot(N_mean.xdf5, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
N_plot5_med.inv <- ggplot(N_mean.xdf5, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = N_data.raw) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )


N_plot5_ex.temp <- ggplot(N_mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = cividis(10)[1], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = cividis(10)[10], linetype = "dashed") + 
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") + 
  #coord_cartesian(ylim = c(0, 40)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("n.mean.inv.ex.plot.png", plot = N_plot5_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(N_plot5_mean, N_plot5_med, nrow = 1)
gridExtra::grid.arrange(N_plot5_mean.inv, N_plot5_med.inv, nrow = 1)
ggsave("n.mean.inv.mort.exp.png", plot = N_plot5_med, width = 5, height = 5)
ggsave("n.mean.inv.life.exp.png", plot = N_plot5_med.inv, width = 5, height = 5)

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
  proportion <- mean(true.Topt > col3 & true.Topt < col4)
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

save(N_param_list5,
     xdf5_N,
     df.new5_N,
     N_mean.xdf5,
     N_plot5_mean, 
     N_plot5_med,
     N_plot5_mean.inv, 
     N_plot5_med.inv,
     N_proportions_list_a5,
     N_proportions_list_b5,
     N_proportions_list_c5,
     rmse.range.mean.inv.n,
     file="mean_inverse_exp_narrow.RData")





