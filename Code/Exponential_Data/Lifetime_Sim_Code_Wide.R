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
set.seed(1252)
for(i in 1:100){
  T1 <- rexp(n = 100, f(10))
  T2 <- rexp(n = 100, f(17))
  T3 <- rexp(n = 100, f(21))
  T4 <- rexp(n = 100, f(25))
  T5 <- rexp(n = 100, f(32))
  data.list[[i]] <- data.frame(
    Temp = c(rep(10, 100), rep(17, 100), rep(21, 100), rep(25, 100), rep(32, 100)), 
    trait = c(T1, T2, T3, T4, T5)
  )
}
for(i in 1:100){
  data.list[[i]]$trait <- ifelse(data.list[[i]]$trait < 1, 1, data.list[[i]]$trait)
}

b <- length(data.list)

# True value vector with a in log space
# Spot for sigma and deviance
true_values <- c(log(true.a), true.Topt, true.c, NA, NA)

xseq <- seq(from = 5, to = 37, length.out = 250)

########################################## Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
# traceplot.data1 <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.ind <- c()
# Create list of matrices to store parameter summary statistics
param_list1 <- vector("list", length = 4)
# Create dataframe to store uncertainty results
exp.ind.uncertainty <- data.frame()
exp.ind.whisk <- data.frame()
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
  
  data.raw <- data.list[[i]]
  data <- data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$Temp
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
    exp.ind.uncertainty <- rbind(exp.ind.uncertainty, new_row)
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
  
  exp.ind.whisk <- rbind(exp.ind.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data1) {
  #   pdf(paste0("traceplot.exp.ind.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
  df.rmse1 <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df1 <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse1))
  # Apply the quadratic equation for each combination
  rmse.df1$value <- apply(rmse.df1, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse1[k, 3]) * (location - df.rmse1[k, 4])^2 + df.rmse1[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse1), function(k) {
    evaluated_curve <- rmse.df1$value[rmse.df1$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.ind <- c(rmse.range.ind, rmse_values)
}

# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new1 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = -10, to = 60, length.out = 250)
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
labs(x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the median mortality rate response
plot1_med <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = med.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
labs(x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the mean lifetime response
plot1_mean.inv <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = avg.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
geom_point(aes(x = Temp, y = trait), data = data.list[[100]]) + 
labs(x = "Temperature", 
     y = "Lifetime") + 
theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the median lifetime response
plot1_med.inv <- ggplot(mean.xdf1, aes(x = point)) + 
geom_line(aes(y = med.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
geom_point(aes(x = Temp, y = trait), data = data.list[[100]]) + 
labs(x = "Temperature", 
     y = "Lifetime") + 
theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )


plot1_ex.temp <- ggplot(mean.xdf1, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
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
ggsave("ind.ex.plot.png", plot = plot1_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(plot1_mean, plot1_med, nrow = 1)
gridExtra::grid.arrange(plot1_mean.inv, plot1_med.inv, nrow = 1)
ggsave("ind.mort.exp.png", plot = plot1_med, width = 5, height = 5)
ggsave("ind.life.exp.png", plot = plot1_med.inv, width = 5, height = 5)

# Assign true values 
true_values <- c(log(true.a), true.Topt, true.c, NA)


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

save(param_list1,
     xdf1,
     df.new1,
     mean.xdf1,
     plot1_mean, 
     plot1_med,
     plot1_mean.inv, 
     plot1_med.inv,
     proportions_list_a1,
     proportions_list_b1,
     proportions_list_c1,
     rmse.range.ind,
     file="individual_exponential.RData")

########################################## Mean Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
# traceplot.data2 <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.mean <- c()
# Create list of matrices to store parameter summary statistics
param_list2 <- vector("list", length = 4)
# Create dataframe to store uncertainty results
exp.mean.uncertainty <- data.frame()
exp.mean.whisk <- data.frame()
for (i in 1:5) {
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
  group_means <- tapply(data.list[[i]]$trait, data.list[[i]]$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
    })
  df.mean <- data.frame(
    Temp = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)), # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
    )
  # Model
  sink("mean_lambda.txt")
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
}", file = "mean_lambda.txt")
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
                model.file = "mean_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    exp.mean.uncertainty <- rbind(exp.mean.uncertainty, new_row)
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
  
  exp.mean.whisk <- rbind(exp.mean.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data2) {
  #   pdf(paste0("traceplot.exp.mean.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
  for(j in 1:nc){
    # For each chain, store mean, median, and hdi bounds for a in log space
    # Statistics are all in a different column, each row is a chain
    # 100 matrices of 5x4 for each parameter
    param_list2[[1]][[i]][j, 1] <- mean(samp[[j]][, 3])
    param_list2[[1]][[i]][j, 2] <- median(samp[[j]][, 3])
    param_list2[[1]][[i]][j, 3] <- hdi(samp[[j]][, 3])[1]
    param_list2[[1]][[i]][j, 4] <- hdi(samp[[j]][, 3])[2]
    # For each chain, store mean, median, and hdi bounds for Topt 
    param_list2[[2]][[i]][j, 1] <- mean(samp[[j]][, 5])
    param_list2[[2]][[i]][j, 2] <- median(samp[[j]][, 5])
    param_list2[[2]][[i]][j, 3] <- hdi(samp[[j]][, 5])[1]
    param_list2[[2]][[i]][j, 4] <- hdi(samp[[j]][, 5])[2]
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
    # For each chain, store mean, median, and hdi bounds for sigma
    param_list2[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list2[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list2[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list2[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse2 <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df2 <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse2))
  # Apply the quadratic equation for each combination
  rmse.df2$value <- apply(rmse.df2, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse2[k, 3]) * (location - df.rmse2[k, 5])^2 + df.rmse2[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse2), function(k) {
    evaluated_curve <- rmse.df2$value[rmse.df2$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.mean <- c(rmse.range.mean, rmse_values)
}

# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new2 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = -10, to = 60, length.out = 250)
# Create a grid with point as each temperature, and iteration being
# the number of rows we selected from samp
# The idea is to evaluate each row of df.new2 as a function at each of those temperatures
xdf2 <- expand.grid(point = xseq, iteration = 1:nrow(df.new2))

# Apply the quadratic equation for each combination
xdf2$value <- apply(xdf2, 1, function(row) {
i <- row["iteration"]
location <- row["point"]
(exp(df.new2[i, 3]) * (location - df.new2[i, 5])^2  + df.new2[i, 1])
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
labs(x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
plot2_med <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = med.value), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
labs(x = "Temperature", 
     y = "Mortality Rate") + 
theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean lifetime response
plot2_mean.inv <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = avg.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
geom_point(aes(x = Temp, y = trait), data = data) + 
labs(x = "Temperature", 
     y = "Lifetime") + 
theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot2_med.inv <- ggplot(mean.xdf2, aes(x = point)) + 
geom_line(aes(y = med.value.inv), color = "black") + 
geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
geom_point(aes(x = Temp, y = trait), data = data) + 
  coord_cartesian(ylim = c(0, 100)) +
labs(x = "Temperature", 
     y = "Lifetime") + 
theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )


plot2_ex.temp <- ggplot(mean.xdf2, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
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
ggsave("mean.ex.plot.png", plot = plot2_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(plot2_mean, plot2_med, nrow = 1)
gridExtra::grid.arrange(plot2_mean.inv, plot2_med.inv, nrow = 1)
ggsave("mean.mort.exp.png", plot = plot2_med, width = 5, height = 5)
ggsave("mean.life.exp.png", plot = plot2_med.inv, width = 5, height = 5)

# Assign true values 
true_values <- c(log(true.a), true.Topt, true.c, NA)

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


save(param_list2,
     xdf2,
     df.new2,
     mean.xdf2,
     plot2_mean, 
     plot2_med,
     plot2_mean.inv, 
     plot2_med.inv,
     proportions_list_a2,
     proportions_list_b2,
     proportions_list_c2,
     rmse.range.mean,
     file="mean_exponential.RData")

########################################## Inverse Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data3 <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.inv <- c()
# Create list of matrices to store parameter summary statistics
param_list3 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
exp.inv.uncertainty <- data.frame()
exp.inv.whisk <- data.frame()
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
  c ~ dgamma(10,1) # Has to be positive
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
  data <- data.list[[i]] 
  trait <- ifelse(1/(data$trait) > 1, 1, 1/(data$trait))
  N.obs <- length(trait)
  temp <- data$Temp
  
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
    exp.inv.uncertainty <- rbind(exp.inv.uncertainty, new_row)
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
  
  exp.inv.whisk <- rbind(exp.inv.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data3) {
  #   pdf(paste0("traceplot.exp.inv.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
    # For each chain, store mean, median, and hdi bounds for sig2
    param_list3[[5]][[i]][j, 1] <- mean(samp[[j]][, 4])
    param_list3[[5]][[i]][j, 2] <- median(samp[[j]][, 4])
    param_list3[[5]][[i]][j, 3] <- hdi(samp[[j]][, 4])[1]
    param_list3[[5]][[i]][j, 4] <- hdi(samp[[j]][, 4])[2]
  }
  df.rmse3 <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df3 <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse3))
  # Apply the quadratic equation for each combination
  rmse.df3$value <- apply(rmse.df3, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse3[k, 3]) * (location - df.rmse3[k, 5])^2 + df.rmse3[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse3), function(k) {
    evaluated_curve <- rmse.df3$value[rmse.df3$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.inv <- c(rmse.range.inv, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new3 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = -10, to = 60, length.out = 250)
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
  labs(x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
plot3_med <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = 1/trait), data = data) + 
  labs(x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean lifetime response
plot3_mean.inv <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = trait), data = data.list[[100]]) + 
  labs(x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot3_med.inv <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = trait), data = data.list[[100]]) + 
  labs(x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

plot3_ex.temp <- ggplot(mean.xdf3, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
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
ggsave("inv.ex.plot.png", plot = plot3_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(plot3_mean, plot3_med, nrow = 1)
gridExtra::grid.arrange(plot3_mean.inv, plot3_med.inv, nrow = 1)
ggsave("inv.mort.exp.png", plot = plot3_med, width = 5, height = 5)
ggsave("inv.life.exp.png", plot = plot3_med.inv, width = 5, height = 5)

# Assign true values 
true_values <- c(log(true.a), true.Topt, true.c, NA, NA)


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

save(param_list3,
     xdf3,
     df.new3,
     mean.xdf3,
     plot3_mean, 
     plot3_med,
     plot3_mean.inv, 
     plot3_med.inv,
     proportions_list_a3,
     proportions_list_b3,
     proportions_list_c3,
     rmse.range.inv,
     file="inverse_exponential.RData")

########################################## Inverse of Mean Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data4 <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.inv.mean <- c()
# Create list of matrices to store parameter summary statistics
param_list4 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
exp.iom.uncertainty <- data.frame()
exp.iom.whisk <- data.frame()
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
  group_means <- tapply(data.list[[i]]$trait, data.list[[i]]$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    Temp = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)), # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- 1/(data$trait)
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
                  model.file = "inv_mean_lambda.txt", n.thin = nt, n.chains = nc, n.burnin = nb, 
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
    exp.iom.uncertainty <- rbind(exp.iom.uncertainty, new_row)
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
  
  exp.iom.whisk <- rbind(exp.iom.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data4) {
  #   pdf(paste0("traceplot.exp.inv.mean.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
  df.rmse4 <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df4 <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse4))
  # Apply the quadratic equation for each combination
  rmse.df4$value <- apply(rmse.df4, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse4[k, 3]) * (location - df.rmse4[k, 5])^2 + df.rmse4[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse4), function(k) {
    evaluated_curve <- rmse.df4$value[rmse.df4$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.inv.mean <- c(rmse.range.inv.mean, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = -10, to = 60, length.out = 250)
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
  labs(x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
plot4_med <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = 1/trait), data = data) + 
  labs(x = "Temperature", 
       y = "Mortality Rate") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean lifetime response
plot4_mean.inv <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = trait), data = data.list[[100]]) + 
  labs(x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot4_med.inv <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = trait), data = data.list[[100]]) + 
  labs(x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

plot4_ex.temp <- ggplot(mean.xdf4, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") + 
  coord_cartesian(ylim = c(0, 2), xlim = c(48, 52)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# CHANGED TO ZOOM IN
ggsave("inv.mean.ex.plot.zoom.png", plot = plot4_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(plot4_mean, plot4_med, nrow = 1)
gridExtra::grid.arrange(plot4_mean.inv, plot4_med.inv, nrow = 1)
ggsave("inv.mean.mort.exp.png", plot = plot4_med, width = 5, height = 5)
ggsave("inv.mean.life.exp.png", plot = plot4_med.inv, width = 5, height = 5)


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


save(param_list4,
     xdf4,
     df.new4,
     mean.xdf4,
     plot4_mean, 
     plot4_med,
     plot4_mean.inv, 
     plot4_med.inv,
     proportions_list_a4,
     proportions_list_b4,
     proportions_list_c4,
     rmse.range.inv.mean,
     file="inverse_mean_exponential.RData")


########################################## Mean of Inverse Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data5 <- sort(sample(1:b, 5))
# Vector for RMSEs for each sample of each model
rmse.range.mean.inv <- c()
# Create list of matrices to store parameter summary statistics
param_list5 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
exp.moi.uncertainty <- data.frame()
exp.moi.whisk <- data.frame()
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
  data.raw$t.trait <- ifelse(data.raw$trait < 1, 1, data.raw$trait)
  data.raw$inv.trait <- 1/(data.raw$t.trait)
  # Get 10 groups of 10, and calculate the means
  group_means <- tapply(data.raw$inv.trait, data.raw$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <-  data.frame(
    Temp = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)), # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data <- df.mean
  trait <- ifelse(data$trait > 1, 1, data$trait)
  N.obs <- length(trait)
  temp <- data$Temp
  # Model
  sink("mean_inv_lambda.txt")
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
}", file = "mean_inv_lambda.txt")
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
    exp.moi.uncertainty <- rbind(exp.moi.uncertainty, new_row)
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
  
  exp.moi.whisk <- rbind(exp.moi.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data5) {
  #   pdf(paste0("traceplot.exp.mean.inv.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
  df.rmse5 <- samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  # Create a grid for temperatures and iterations
  rmse.df5 <- expand.grid(point = xseq, iteration = 1:nrow(df.rmse5))
  # Apply the quadratic equation for each combination
  rmse.df5$value <- apply(rmse.df5, 1, function(row) {
    k <- row["iteration"]
    location <- row["point"]
    (exp(df.rmse5[k, 3]) * (location - df.rmse5[k, 5])^2 + df.rmse5[k, 1])
  })
  # Calculate RMSE between the evaluated curve and the true curve
  true_curve <- f(xseq)   
  rmse_values <- sapply(1:nrow(df.rmse5), function(k) {
    evaluated_curve <- rmse.df5$value[rmse.df5$iteration == k]
    sqrt(mean((evaluated_curve - true_curve)^2))
  })
  # Append the RMSE values for this model to the main vector
  rmse.range.mean.inv <- c(rmse.range.mean.inv, rmse_values)
}

# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5 = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = -10, to = 60, length.out = 250)
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
mean.xdf5$true_curve <- f(xseq)
# Plot mean mortality rate response
plot5_mean <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = avg.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  # geom_point(aes(x = T, y = 1/trait), data = data.raw) + 
  labs(x = "Temperature", 
       y = "Mortality Rate") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
plot5_med <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = med.value), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = trait), data = data) + 
  labs(x = "Temperature", 
       y = "Mortality Rate") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean lifetime response
plot5_mean.inv <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = avg.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = trait), data = data.list[[100]]) + 
  labs(x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot5_med.inv <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
  geom_point(aes(x = Temp, y = trait), data = data.list[[100]]) + 
  labs(x = "Temperature", 
       y = "Lifetime") + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

plot5_ex.temp <- ggplot(mean.xdf5, aes(x = point)) + 
  geom_line(aes(y = med.value.inv), color = "black") + 
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = rocket(10)[5], alpha = 0.3) + 
  geom_line(aes(y = 1/true_curve), color = rocket(10)[1], linetype = "dashed") + 
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
ggsave("mean.inv.ex.plot.png", plot = plot5_ex.temp, width = 5, height = 5)

gridExtra::grid.arrange(plot5_mean, plot5_med, nrow = 1)
gridExtra::grid.arrange(plot5_mean.inv, plot5_med.inv, nrow = 1)
ggsave("mean.inv.mort.exp.png", plot = plot5_med, width = 5, height = 5)
ggsave("mean.inv.life.exp.png", plot = plot5_med.inv, width = 5, height = 5)

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


save(param_list5,
     xdf5,
     df.new5,
     mean.xdf5,
     plot5_mean, 
     plot5_med,
     plot5_mean.inv, 
     plot5_med.inv,
     proportions_list_a5,
     proportions_list_b5,
     proportions_list_c5,
     rmse.range.mean.inv,
     file="mean_inverse_exponential.RData")




