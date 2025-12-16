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

# True quadratic
f <- function(x, a = 0.001094, mu = 22.3, c = 0.066){
  a * (x - mu)^2 + c
}

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
    Temp = c(rep(10, 100), rep(17, 100), rep(21, 100), rep(25, 100), rep(32, 100)), 
    trait = c(T1W, T2W, T3W, T4W, T5W)
  )
}
for(i in 1:100){
  Weibull_data.list[[i]]$trait <- ifelse(Weibull_data.list[[i]]$trait < 1, 1, Weibull_data.list[[i]]$trait)
}

w <- length(Weibull_data.list)

# True value vector with a in log space
Weibull_true_values <- c(log(true.a), true.Topt, true.c, NA, sh.est)

# Raw data for HDI plot
data.raw.w <- Weibull_data.list[[100]]

##########################################Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data.wb1 <- sort(sample(1:w, 5))
# Vector for RMSEs for each sample of each model
rmse.range.ind.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list1 <- vector("list", length = 4)
# Create dataframe to store uncertainty results
wb.ind.uncertainty <- data.frame()
wb.ind.whisk <- data.frame()
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
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  data.raw.w <- Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$Temp
  
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
  xseq <- seq(from = 5, to = 37, length.out = 250)
  for (j in 1:nc) {
    # Compute HDI bounds for parameters of interest (Topt and inv of max lifetime)
    hdi_bounds.Topt <- hdi(samp[[j]][, "Topt"])
    hdi_bounds.life <- hdi(samp[[j]][, "c"])
    # Build a row
    new_row <- data.frame(
      Method = "Individual Data, Exp",
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
    wb.ind.uncertainty <- rbind(wb.ind.uncertainty, new_row)
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
    Method = "Individual Data, Exp",
    Model = i,
    Range.Low = rangeLow,
    Range.High = rangeHigh
  )
  
  wb.ind.whisk <- rbind(wb.ind.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data.wb1) {
  #   pdf(paste0("traceplot.wb.ind.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
xdf1_W <- xdf1_W|> mutate(trunc.inv = 1/value)
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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the median lifetime response
plot1_med_W <- ggplot(mean.xdf1_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the mean mortality rate response
plot1_mean.inv_W <- ggplot(mean.xdf1_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_point(aes(x = Temp, y = trait), data = data.raw.w) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot the median mortality rate response
plot1_med.inv_W <- ggplot(mean.xdf1_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 30)) +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )


plot1_ex.temp.wb <- ggplot(mean.xdf1_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  #coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("ind.ex.plot.wb.png", plot = plot1_ex.temp.wb, width = 5, height = 5)

gridExtra::grid.arrange(plot1_mean_W, plot1_med_W, nrow = 1)
gridExtra::grid.arrange(plot1_mean.inv_W, plot1_med.inv_W, nrow = 1)
ggsave("ind.mort.wb.png", plot = plot1_med_W, width = 5, height = 5)
ggsave("ind.life.wb.png", plot = plot1_med.inv_W, width = 5, height = 5)

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


save(Weibull_param_list1,
     xdf1_W,
     df.new1_W,
     mean.xdf1_W,
     plot1_mean_W, 
     plot1_med_W,
     plot1_mean.inv_W, 
     plot1_med.inv_W,
     Weibull_proportions_list_a1,
     Weibull_proportions_list_b1,
     Weibull_proportions_list_c1,
     rmse.range.ind.wb,
     file="Weibull_individual.RData")


##########################################Mean Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data.wb2 <- sort(sample(1:w, 5))
# Vector for RMSEs for each sample of each model
rmse.range.mean.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list2 <- vector("list", length = 4)
# Create dataframe to store uncertainty results
wb.mean.uncertainty <- data.frame()
wb.mean.whisk <- data.frame()
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
  group_means <- tapply(Weibull_data.list[[i]]$trait, Weibull_data.list[[i]]$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    Temp = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)),  
    trait = unlist(group_means)  
  )
  # Model
  sink("Weibull_mean_lambda.txt")
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
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  
  data.raw.w <- Weibull_data.list[[i]]
  data <- df.mean 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$Temp
  
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
    wb.mean.uncertainty <- rbind(wb.mean.uncertainty, new_row)
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
  
  wb.mean.whisk <- rbind(wb.mean.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data.wb2) {
  #   pdf(paste0("traceplot.wb.mean.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
xseq <- seq(from = -10, to = 60, length.out = 250)
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
xdf2_W <- xdf2_W|> mutate(trunc.inv = 1/value)
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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot2_med_W <- ggplot(mean.xdf2_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  #geom_point(aes(x = Temp, y = trait), data = data) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean mortality rate response
plot2_mean.inv_W <- ggplot(mean.xdf2_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
plot2_med.inv_W <- ggplot(mean.xdf2_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  coord_cartesian(ylim = c(0,30)) +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

plot2_ex.temp.wb <- ggplot(mean.xdf2_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  #coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("mean.ex.plot.wb.png", plot = plot2_ex.temp.wb, width = 5, height = 5)

gridExtra::grid.arrange(plot2_mean_W, plot2_med_W, nrow = 1)
gridExtra::grid.arrange(plot2_mean.inv_W, plot2_med.inv_W, nrow = 1)
ggsave("mean.mort.wb.png", plot = plot2_med_W, width = 5, height = 5)
ggsave("mean.life.wb.png", plot = plot2_med.inv_W, width = 5, height = 5)


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

save(Weibull_param_list2,
     xdf2_W,
     df.new2_W,
     mean.xdf2_W,
     plot2_mean_W, 
     plot2_med_W,
     plot2_mean.inv_W, 
     plot2_med.inv_W,
     Weibull_proportions_list_a2,
     Weibull_proportions_list_b2,
     Weibull_proportions_list_c2,
     rmse.range.mean.wb,
     file="Weibull_mean.RData")

##########################################Inverse Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data.wb3 <- sort(sample(1:w, 5))
# Vector for RMSEs for each sample of each model
rmse.range.inv.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list3 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
wb.inv.uncertainty <- data.frame()
wb.inv.whisk <- data.frame()
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
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  
  data.raw.w <- Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$Temp
  
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
    wb.inv.uncertainty <- rbind(wb.inv.uncertainty, new_row)
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
  
  wb.inv.whisk <- rbind(wb.inv.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data.wb3) {
  #   pdf(paste0("traceplot.wb.inv.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
xseq <- seq(from = -10, to = 60, length.out = 250)
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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot3_med_W <- ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = 1/trait), data = data) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean mortality rate response
plot3_mean.inv_W <- ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  #geom_point(aes(x = Temp, y = 1/trait), data = data) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
plot3_med.inv_W <- ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  #geom_point(aes(x = Temp, y = 1/trait), data = data) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

plot3_ex.temp.wb <- ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  #coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("inv.ex.plot.wb.png", plot = plot3_ex.temp.wb, width = 5, height = 5)


gridExtra::grid.arrange(plot3_mean_W, plot3_med_W, nrow = 1)
gridExtra::grid.arrange(plot3_mean.inv_W, plot3_med.inv_W, nrow = 1)
ggsave("inv.mort.wb.png", plot = plot3_med_W, width = 5, height = 5)
ggsave("inv.life.wb.png", plot = plot3_med.inv_W, width = 5, height = 5)



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

save(Weibull_param_list3,
     xdf3_W,
     df.new3_W,
     mean.xdf3_W,
     plot3_mean_W, 
     plot3_med_W,
     plot3_mean.inv_W, 
     plot3_med.inv_W,
     Weibull_proportions_list_a3,
     Weibull_proportions_list_b3,
     Weibull_proportions_list_c3,
     rmse.range.inv.wb,
     file="Weibull_inverse.RData")


##########################################Inverse of Mean Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data.wb4 <- sort(sample(1:w, 5))
# Vector for RMSEs for each sample of each model
rmse.range.inv.mean.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list4 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
wb.iom.uncertainty <- data.frame()
wb.iom.whisk <- data.frame()
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
  group_means <- tapply(data$trait, data$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    Temp = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data.raw.w <- Weibull_data.list[[i]]
  data <- df.mean
  trait <- 1/(data$trait)
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
  for(k in 1:nc){
    inits[[k]] = GenInits()
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
    wb.iom.uncertainty <- rbind(wb.iom.uncertainty, new_row)
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
  
  wb.iom.whisk <- rbind(wb.iom.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data.wb4) {
  #   pdf(paste0("traceplot.wb.inv.mean.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
  rmse.range.inv.mean.wb <- c(rmse.range.inv.mean.wb, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new4_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = -10, to = 60, length.out = 250)
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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot4_med_W <- ggplot(mean.xdf4_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = 1/trait), data = data) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean mortality rate response
plot4_mean.inv_W <- ggplot(mean.xdf4_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
plot4_med.inv_W <- ggplot(mean.xdf4_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

plot4_ex.temp.wb <- ggplot(mean.xdf4_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  coord_cartesian(ylim = c(0, 2), xlim = c(49, 52)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# CHANGED TO ZOOM IN
ggsave("inv.mean.ex.plot.wb.zoom.png", plot = plot4_ex.temp.wb, width = 5, height = 5)

gridExtra::grid.arrange(plot4_mean_W, plot4_med_W, nrow = 1)
gridExtra::grid.arrange(plot4_mean.inv_W, plot4_med.inv_W, nrow = 1)
ggsave("inv.mean.mort.wb.png", plot = plot4_med_W, width = 5, height = 5)
ggsave("inv.mean.life.wb.png", plot = plot4_med.inv_W, width = 5, height = 5)


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

save(Weibull_param_list4,
     xdf4_W,
     df.new4_W,
     mean.xdf4_W,
     plot4_mean_W, 
     plot4_med_W,
     plot4_inv.mean_W, 
     plot4_med.inv_W,
     Weibull_proportions_list_a4,
     Weibull_proportions_list_b4,
     Weibull_proportions_list_c4,
     rmse.range.inv.mean.wb,
     file="Weibull_inverse_mean.RData")

##########################################Mean of Inverse Lambda #############################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data.wb5 <- sort(sample(1:w, 5))
# Vector for RMSEs for each sample of each model
rmse.range.mean.inv.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list5 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
wb.moi.uncertainty <- data.frame()
wb.moi.whisk <- data.frame()
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
  group_means <- tapply(data$inv.trait, data$Temp, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean <- data.frame(
    Temp = c(rep(10, 10), rep(17, 10), rep(21, 10), rep(25, 10), rep(32, 10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data.raw.w <- Weibull_data.list[[i]]
  data <- df.mean
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$Temp
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
  for(k in 1:nc){
    inits[[k]] = GenInits()
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
    wb.moi.uncertainty <- rbind(wb.moi.uncertainty, new_row)
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
  
  wb.moi.whisk <- rbind(wb.moi.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data.wb5) {
  #   pdf(paste0("traceplot.wb.mean.inv.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
  rmse.range.mean.inv.wb <- c(rmse.range.mean.inv.wb, rmse_values)
  
}
# Take samp, the mcmc list, chain 1 and save every fourth row
# Saved every fourth row to avoid autocorrelation
df.new5_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
# Create a sequence of temperatures from 10 to 35
xseq <- seq(from = -10, to = 60, length.out = 250)
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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot5_med_W <- ggplot(mean.xdf5_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean mortality rate response
plot5_mean.inv_W <- ggplot(mean.xdf5_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
plot5_med.inv_W <- ggplot(mean.xdf5_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

plot5_ex.temp.wb <- ggplot(mean.xdf5_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  #coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("mean.inv.ex.plot.wb.png", plot = plot5_ex.temp.wb, width = 5, height = 5)

gridExtra::grid.arrange(plot5_mean_W, plot5_med_W, nrow = 1)
gridExtra::grid.arrange(plot5_mean.inv_W, plot5_med.inv_W, nrow = 1)
ggsave("mean.inv.mort.wb.png", plot = plot5_med_W, width = 5, height = 5)
ggsave("mean.inv.life.wb.png", plot = plot5_med.inv_W, width = 5, height = 5)


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

save(Weibull_param_list5,
     xdf5_W,
     df.new5_W,
     mean.xdf5_W,
     plot5_mean_W, 
     plot5_med_W,
     plot5_mean.inv_W, 
     plot5_med.inv_W,
     Weibull_proportions_list_a5,
     Weibull_proportions_list_b5,
     Weibull_proportions_list_c5,
     rmse.range.mean.inv.wb,
     file="Weibull_mean_inverse.RData")



########################################## Weibull ###########################################
# Choose 5 (or any number of choosing) random datasets for convergence checking
# Uncomment to check convergence
traceplot.data.wb6 <- sort(sample(1:w, 5))
# Vector for RMSEs for each sample of each model
rmse.range.ind.wb.wb <- c()
# Create list of matrices to store parameter summary statistics
Weibull_param_list6 <- vector("list", length = 5)
# Create dataframe to store uncertainty results
wb.ind.wb.uncertainty <- data.frame()
wb.ind.wb.whisk <- data.frame()
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
  la ~ dnorm(0, 1/10) 
  Topt ~ dnorm(22, 1/10) 
  c ~ dexp(0.5) #Has to be positive
  shape ~ dexp(0.5)
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
    shape <- rexp(1, 0.5)
    list(
      la = la,             
      c = c,
      Topt = Topt,
      shape = shape
    )
  }
  for(k in 1:nc){
    inits[[k]] = GenInits()
  }
  
  data.raw.w <- Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$Temp
  
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
  xseq <- seq(from = 5, to = 37, length.out = 250)
  for (j in 1:nc) {
    # Compute HDI bounds for parameters of interest (Topt and inv of max lifetime)
    hdi_bounds.Topt <- hdi(samp[[j]][, "Topt"])
    hdi_bounds.life <- hdi(samp[[j]][, "c"])
    # Build a row
    new_row <- data.frame(
      Method = "Individual Data, WB",
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
    wb.ind.wb.uncertainty <- rbind(wb.ind.wb.uncertainty, new_row)
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
    Method = "Individual Data, WB", 
    Model = i,
    Range.Low = rangeLow,
    Range.High = rangeHigh
  )
  
  wb.ind.wb.whisk <- rbind(wb.ind.wb.whisk, new_row)
  # Check traceplots for 5 randomly sampled models
  # Uncomment to check convergence
  # if (i %in% traceplot.data.wb6) {
  #   pdf(paste0("traceplot.wb.wb.ind.", i, ".pdf"))
  #   plot(mod.fit.mcmc)
  #   dev.off()
  # }
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
xseq <- seq(from = -10, to = 60, length.out = 250)
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
xdf6_W <- xdf6_W|> mutate(trunc.inv = 1/value)
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
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median lifetime response
plot6_med_W <- ggplot(mean.xdf6_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = true_curve), color = plasma(10)[2], linetype = "dashed") +
  # geom_point(aes(x = T, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Mortality Rate") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot mean mortality rate response
plot6_mean.inv_W <- ggplot(mean.xdf6_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data.raw.w) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
# Plot median mortality rate response
plot6_med.inv_W <- ggplot(mean.xdf6_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_point(aes(x = Temp, y = trait), data = data.raw.w) +
  coord_cartesian(ylim = c(0, 30)) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )

plot6_ex.temp.wb <- ggplot(mean.xdf6_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = plasma(10)[8], alpha = 0.3) +
  geom_line(aes(y = 1/true_curve), color = plasma(10)[2], linetype = "dashed") +
  geom_hline(yintercept = 1, linewidth = 1, color = "green") +
  geom_point(data = data.frame(x = c(-6.935, 51.542), y = c(1, 1)), 
             aes(x = x, y = y), 
             shape = 1, 
             size = 10) +
  labs(x = "Temperature", 
       y = "Lifetime") +
  #coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis title size
  )
ggsave("ind.wb.ex.plot.wb.png", plot = plot6_ex.temp.wb, width = 5, height = 5)

gridExtra::grid.arrange(plot6_mean_W, plot6_med_W, nrow = 1)
gridExtra::grid.arrange(plot6_mean.inv_W, plot6_med.inv_W, nrow = 1)
ggsave("ind.wb.mort.wb.png", plot = plot6_med_W, width = 5, height = 5)
ggsave("ind.wb.life.wb.png", plot = plot6_med.inv_W, width = 5, height = 5)


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

save(Weibull_param_list6,
     xdf6_W,
     df.new6_W,
     mean.xdf6_W,
     plot6_mean_W, 
     plot6_med_W,
     plot6_mean.inv_W, 
     plot6_med.inv_W,
     Weibull_proportions_list_a6,
     Weibull_proportions_list_b6,
     Weibull_proportions_list_c6,
     rmse.range.ind.wb.wb,
     file="Weibull_weibull_data.RData")

