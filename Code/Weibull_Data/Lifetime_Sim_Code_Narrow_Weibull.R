library(rjags)
library(dplyr)

true.exp.a<- 0.00826
true.exp.b<- 0.37
true.exp.c<- 1.406
f.exp<-function(x,a=0.00826,b=0.37,c=1.406){
  exp(a*x^2-b*x+c)
}
sh.est<-3
my.inv.lambda<-function(x,a=0.00826,b=0.37,c=1.406,sh=3){ 
  (log(2))^(1/sh)/f.exp(x,a,b,c) ## fx0 is the median here
}

Weibull_N_data.list<-list()
set.seed(52)
for(i in 1:100){
  T1W<-rweibull(n=100,scale=my.inv.lambda(13),shape=sh.est)
  T2W<-rweibull(n=100,scale=my.inv.lambda(17),shape=sh.est)
  T3W<-rweibull(n=100,scale=my.inv.lambda(21),shape=sh.est)
  T4W<-rweibull(n=100,scale=my.inv.lambda(25),shape=sh.est)
  T5W<-rweibull(n=100,scale=my.inv.lambda(29),shape=sh.est)
  Weibull_N_data.list[[i]]<-data.frame(
    T=c(rep(13,100),rep(17,100),rep(21,100),rep(25,100),rep(29,100)),
    trait=c(T1W,T2W,T3W,T4W,T5W)
  )
}
b<-length(Weibull_N_data.list)

Weibull_true_values<-c(true.exp.a, true.exp.b, true.exp.c, NA,NA)

########################################## Lambda #############################################

Weibull_N_param_list1 <- vector("list", length = 4)
for (i in 1:4) {
  # Initialize the inner list
  Weibull_N_inner_list1 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data1 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list1[[j]] <- Weibull_N_matrix_data1
  }
  Weibull_N_param_list1[[i]] <- Weibull_N_inner_list1
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("Weibull_n_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - b.l * temp[i] + c)
    trait[i] ~ dexp(mu[i])
  }
  }",file="Weibull_n_lambda.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c","shape")
  inits<-function(){list(
    la = log(0.001),
    b.l = 0.2,
    c = 0.6,
    shape=sh.est
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data <- Weibull_N_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="Weibull_n_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list1[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list1[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list1[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list1[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list1[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list1[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list1[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list1[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list1[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list1[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list1[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list1[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list1[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list1[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list1[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list1[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
  }
  
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf1<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                   med.value=median(value),
                                                   lower.hdi=hdi(value)[1],
                                                   upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf1$true_curve <- true_curve(xseq)
  Weibull_N_plot1_mean<-ggplot(Weibull_N_mean.xdf1, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot1_med<-ggplot(Weibull_N_mean.xdf1, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot1_mean,Weibull_N_plot1_med,nrow=1)

Weibull_N_hist_list1 <- vector("list", length = 4)

# Loop through indices 1 to 4
for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list1 <- lapply(Weibull_N_param_list1[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector1 <- unlist(Weibull_N_first_column_list1)
  x<-seq(min(Weibull_N_combined_vector1)-1,max(Weibull_N_combined_vector1)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list1[[i]] <- hist(Weibull_N_combined_vector1, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                            xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 0.1), col="purple", lty=2, lwd=2)
  }
}


Weibull_N_proportions_list_a1<- lapply(Weibull_N_param_list1[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a1))

Weibull_N_proportions_list_b1<- lapply(Weibull_N_param_list1[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b1))

Weibull_N_proportions_list_c1<- lapply(Weibull_N_param_list1[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c1))

Weibull_N_calculate_rmse_a_lambda <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_lambda <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_lambda <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_lambda <- sapply(Weibull_N_param_list1[[1]], Weibull_N_calculate_rmse_a_lambda)
Weibull_N_rmse_values_b_lambda <- sapply(Weibull_N_param_list1[[2]], Weibull_N_calculate_rmse_b_lambda)
Weibull_N_rmse_values_c_lambda <- sapply(Weibull_N_param_list1[[3]], Weibull_N_calculate_rmse_c_lambda)

Weibull_N_tab.lambda<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.lambda)<-c("a","b","c")
colnames(Weibull_N_tab.lambda)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_N_tab.lambda[1,1]<-true.exp.a
Weibull_N_tab.lambda[2,1]<-true.exp.b
Weibull_N_tab.lambda[3,1]<-true.exp.c
Weibull_N_tab.lambda[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list1[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.lambda[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list1[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.lambda[3,2]<-mean(unlist(lapply(Weibull_N_param_list1[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.lambda[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list1[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.lambda[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list1[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.lambda[3,3]<-hdi(unlist(lapply(Weibull_N_param_list1[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.lambda[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list1[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.lambda[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list1[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.lambda[3,4]<-hdi(unlist(lapply(Weibull_N_param_list1[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.lambda[1,5]<-mean(unlist(Weibull_N_proportions_list_a1))
Weibull_N_tab.lambda[2,5]<-mean(unlist(Weibull_N_proportions_list_b1))
Weibull_N_tab.lambda[3,5]<-mean(unlist(Weibull_N_proportions_list_c1))
Weibull_N_tab.lambda[1,6]<-mean(Weibull_N_rmse_values_a_lambda)
Weibull_N_tab.lambda[2,6]<-mean(Weibull_N_rmse_values_b_lambda)
Weibull_N_tab.lambda[3,6]<-mean(Weibull_N_rmse_values_c_lambda)


########################################## Mean Lambda #############################################

Weibull_N_param_list2 <- vector("list", length = 4)
for (i in 1:4) {
  # Initialize the inner list
  Weibull_N_inner_list2 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data2 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list2[[j]] <- Weibull_N_matrix_data2
  }
  Weibull_N_param_list2[[i]] <- Weibull_N_inner_list2
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  group_means <- tapply(Weibull_N_data.list[[i]]$trait, Weibull_N_data.list[[i]]$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean<- data.frame(
    T=c(rep(13,10),rep(17,10),rep(21,10),rep(25,10),rep(29,10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  # Model
  sink("Weibull_n_mean_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - b.l * temp[i] + c)
    trait[i] ~ dexp(mu[i])
  }
  }",file="Weibull_n_mean_lambda.txt")
  
  # Settins
  parameters <- c("la", "b.l", "c")
  inits<-function(){list(
    la = log(0.05),
    b.l = 0.2,
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
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="Weibull_n_mean_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list2[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list2[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list2[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list2[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list2[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list2[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list2[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list2[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list2[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list2[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list2[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list2[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list2[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list2[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list2[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list2[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
  }
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf2<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                   med.value=median(value),
                                                   lower.hdi=hdi(value)[1],
                                                   upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf2$true_curve <- true_curve(xseq)
  Weibull_N_plot2_mean<-ggplot(Weibull_N_mean.xdf2, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot2_med<-ggplot(Weibull_N_mean.xdf2, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot2_mean,Weibull_N_plot2_med,nrow=1)

Weibull_N_hist_list2 <- vector("list", length = 4)

# Loop through indices 1 to 4
for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list2 <- lapply(Weibull_N_param_list2[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector2 <- unlist(Weibull_N_first_column_list2)
  x<-seq(min(Weibull_N_combined_vector2)-1,max(Weibull_N_combined_vector2)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list2[[i]] <- hist(Weibull_N_combined_vector2, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                            xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 0.1), col="purple", lty=2, lwd=2)
  }
}


Weibull_N_proportions_list_a2<- lapply(Weibull_N_param_list2[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a2))

Weibull_N_proportions_list_b2<- lapply(Weibull_N_param_list2[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b2))

Weibull_N_proportions_list_c2<- lapply(Weibull_N_param_list2[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c2))

Weibull_N_calculate_rmse_a_lambda_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_lambda_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_lambda_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_lambda_mean <- sapply(Weibull_N_param_list2[[1]], Weibull_N_calculate_rmse_a_lambda_mean)
Weibull_N_rmse_values_b_lambda_mean <- sapply(Weibull_N_param_list2[[2]], Weibull_N_calculate_rmse_b_lambda_mean)
Weibull_N_rmse_values_c_lambda_mean <- sapply(Weibull_N_param_list2[[3]], Weibull_N_calculate_rmse_c_lambda_mean)


Weibull_N_tab.lambda.mean<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.lambda.mean)<-c("a","b","c")
colnames(Weibull_N_tab.lambda.mean)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_N_tab.lambda.mean[1,1]<-true.exp.a
Weibull_N_tab.lambda.mean[2,1]<-true.exp.b
Weibull_N_tab.lambda.mean[3,1]<-true.exp.c
Weibull_N_tab.lambda.mean[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list2[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.lambda.mean[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list2[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.lambda.mean[3,2]<-mean(unlist(lapply(Weibull_N_param_list2[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.lambda.mean[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list2[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.lambda.mean[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list2[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.lambda.mean[3,3]<-hdi(unlist(lapply(Weibull_N_param_list2[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.lambda.mean[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list2[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.lambda.mean[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list2[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.lambda.mean[3,4]<-hdi(unlist(lapply(Weibull_N_param_list2[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.lambda.mean[1,5]<-mean(unlist(Weibull_N_proportions_list_a2))
Weibull_N_tab.lambda.mean[2,5]<-mean(unlist(Weibull_N_proportions_list_b2))
Weibull_N_tab.lambda.mean[3,5]<-mean(unlist(Weibull_N_proportions_list_c2))
Weibull_N_tab.lambda.mean[1,6]<-mean(Weibull_N_rmse_values_a_lambda_mean)
Weibull_N_tab.lambda.mean[2,6]<-mean(Weibull_N_rmse_values_b_lambda_mean)
Weibull_N_tab.lambda.mean[3,6]<-mean(Weibull_N_rmse_values_c_lambda_mean)



########################################## Inverse Lambda #############################################

Weibull_N_param_list3 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_N_inner_list3 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data3 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list3[[j]] <- Weibull_N_matrix_data3
  }
  Weibull_N_param_list3[[i]] <- Weibull_N_inner_list3
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  #Model
  sink("Weibull_n_inv_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  sig2~dexp(7000)
  tau<-1/sig2
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - b.l * temp[i] + c)
    trait[i] ~ dnorm(mu[i],tau) T(0,)
  }
}",file="Weibull_n_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c","sig2")
  inits<-function(){list(
    la = log(0.001),
    b.l = 0.2,
    c = 0.6,
    sig2=2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data <- Weibull_N_data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="Weibull_n_inv_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list3[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list3[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list3[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list3[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list3[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list3[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list3[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list3[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list3[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list3[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list3[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list3[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list3[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list3[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list3[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list3[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_N_param_list3[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_N_param_list3[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_N_param_list3[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_N_param_list3[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf3<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                   med.value=median(value),
                                                   lower.hdi=hdi(value)[1],
                                                   upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf3$true_curve <- true_curve(xseq)
  Weibull_N_plot3_mean<-ggplot(Weibull_N_mean.xdf3, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot3_med<-ggplot(Weibull_N_mean.xdf3, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot3_mean,Weibull_N_plot3_med,nrow=1)

Weibull_N_hist_list3 <- vector("list", length = 5)

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list3 <- lapply(Weibull_N_param_list3[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector3 <- unlist(Weibull_N_first_column_list3)
  x<-seq(min(Weibull_N_combined_vector3)-1,max(Weibull_N_combined_vector3)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list3[[i]] <- hist(Weibull_N_combined_vector3, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                            xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 7000), col="purple", lty=2, lwd=2)
  }
}


Weibull_N_proportions_list_a3<- lapply(Weibull_N_param_list3[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a3))

Weibull_N_proportions_list_b3<- lapply(Weibull_N_param_list3[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b3))

Weibull_N_proportions_list_c3<- lapply(Weibull_N_param_list3[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c3))

Weibull_N_calculate_rmse_a_lambda_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_lambda_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_lambda_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_lambda_inv <- sapply(Weibull_N_param_list3[[1]], Weibull_N_calculate_rmse_a_lambda_inv)
Weibull_N_rmse_values_b_lambda_inv <- sapply(Weibull_N_param_list3[[2]], Weibull_N_calculate_rmse_b_lambda_inv)
Weibull_N_rmse_values_c_lambda_inv <- sapply(Weibull_N_param_list3[[3]], Weibull_N_calculate_rmse_c_lambda_inv)

Weibull_N_tab.inv.lambda<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.inv.lambda)<-c("a","b","c")
colnames(Weibull_N_tab.inv.lambda)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_N_tab.inv.lambda[1,1]<-true.exp.a
Weibull_N_tab.inv.lambda[2,1]<-true.exp.b
Weibull_N_tab.inv.lambda[3,1]<-true.exp.c
Weibull_N_tab.inv.lambda[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list3[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.inv.lambda[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list3[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.inv.lambda[3,2]<-mean(unlist(lapply(Weibull_N_param_list3[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.inv.lambda[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list3[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.inv.lambda[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list3[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.inv.lambda[3,3]<-hdi(unlist(lapply(Weibull_N_param_list3[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.inv.lambda[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list3[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.inv.lambda[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list3[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.inv.lambda[3,4]<-hdi(unlist(lapply(Weibull_N_param_list3[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.inv.lambda[1,5]<-mean(unlist(Weibull_N_proportions_list_a3))
Weibull_N_tab.inv.lambda[2,5]<-mean(unlist(Weibull_N_proportions_list_b3))
Weibull_N_tab.inv.lambda[3,5]<-mean(unlist(Weibull_N_proportions_list_c3))
Weibull_N_tab.inv.lambda[1,6]<-mean(Weibull_N_rmse_values_a_lambda_inv)
Weibull_N_tab.inv.lambda[2,6]<-mean(Weibull_N_rmse_values_b_lambda_inv)
Weibull_N_tab.inv.lambda[3,6]<-mean(Weibull_N_rmse_values_c_lambda_inv)



########################################## Mean Inverse Lambda #############################################

Weibull_N_param_list4 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_N_inner_list4 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data4 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list4[[j]] <- Weibull_N_matrix_data4
  }
  Weibull_N_param_list4[[i]] <- Weibull_N_inner_list4
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- Weibull_N_data.list[[i]] 
  group_means <- tapply(data$trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean<- data.frame(
    T=c(rep(13,10),rep(17,10),rep(21,10),rep(25,10),rep(29,10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data<-df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  #Model
  sink("weibull_n_mean_inv_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  sig2~dexp(250)
  tau<-1/sig2
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - b.l * temp[i] + c)
    trait[i] ~ dnorm(mu[i],tau) T(0,)
  }
}",file="weibull_n_mean_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c","sig2")
  inits<-function(){list(
    la = log(0.001),
    b.l = 0.2,
    c = 0.6,
    sig2=2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="weibull_n_mean_inv_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list4[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list4[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list4[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list4[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list4[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list4[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list4[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list4[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list4[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list4[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list4[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list4[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list4[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list4[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list4[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list4[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_N_param_list4[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_N_param_list4[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_N_param_list4[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_N_param_list4[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf4<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                   med.value=median(value),
                                                   lower.hdi=hdi(value)[1],
                                                   upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf4$true_curve <- true_curve(xseq)
  Weibull_N_plot4_mean<-ggplot(Weibull_N_mean.xdf4, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot4_med<-ggplot(Weibull_N_mean.xdf4, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot4_mean,Weibull_N_plot4_med,nrow=1)

Weibull_N_hist_list4 <- vector("list", length = 4)

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list4 <- lapply(Weibull_N_param_list4[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector4 <- unlist(Weibull_N_first_column_list4)
  x<-seq(min(Weibull_N_combined_vector4)-1,max(Weibull_N_combined_vector4)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list4[[i]] <- hist(Weibull_N_combined_vector4, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                            xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 250), col="purple", lty=2, lwd=2)
  }
}


Weibull_N_proportions_list_a4<- lapply(Weibull_N_param_list4[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a4))

Weibull_N_proportions_list_b4<- lapply(Weibull_N_param_list4[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b4))

Weibull_N_proportions_list_c4<- lapply(Weibull_N_param_list4[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c4))

Weibull_N_calculate_rmse_a_lambda_mean_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_lambda_mean_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_lambda_mean_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_lambda_mean_inv <- sapply(Weibull_N_param_list4[[1]], Weibull_N_calculate_rmse_a_lambda_mean_inv)
Weibull_N_rmse_values_b_lambda_mean_inv <- sapply(Weibull_N_param_list4[[2]], Weibull_N_calculate_rmse_b_lambda_mean_inv)
Weibull_N_rmse_values_c_lambda_mean_inv <- sapply(Weibull_N_param_list4[[3]], Weibull_N_calculate_rmse_c_lambda_mean_inv)

Weibull_N_tab.mean.inv.lambda<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.mean.inv.lambda)<-c("a","b","c")
colnames(Weibull_N_tab.mean.inv.lambda)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage", "RMSE")
Weibull_N_tab.mean.inv.lambda[1,1]<-true.exp.a
Weibull_N_tab.mean.inv.lambda[2,1]<-true.exp.b
Weibull_N_tab.mean.inv.lambda[3,1]<-true.exp.c
Weibull_N_tab.mean.inv.lambda[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list4[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.mean.inv.lambda[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list4[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.mean.inv.lambda[3,2]<-mean(unlist(lapply(Weibull_N_param_list4[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.mean.inv.lambda[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list4[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.mean.inv.lambda[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list4[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.mean.inv.lambda[3,3]<-hdi(unlist(lapply(Weibull_N_param_list4[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.mean.inv.lambda[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list4[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.mean.inv.lambda[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list4[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.mean.inv.lambda[3,4]<-hdi(unlist(lapply(Weibull_N_param_list4[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.mean.inv.lambda[1,5]<-mean(unlist(Weibull_N_proportions_list_a4))
Weibull_N_tab.mean.inv.lambda[2,5]<-mean(unlist(Weibull_N_proportions_list_b4))
Weibull_N_tab.mean.inv.lambda[3,5]<-mean(unlist(Weibull_N_proportions_list_c4))
Weibull_N_tab.mean.inv.lambda[1,6]<-mean(Weibull_N_rmse_values_a_lambda_mean_inv)
Weibull_N_tab.mean.inv.lambda[2,6]<-mean(Weibull_N_rmse_values_b_lambda_mean_inv)
Weibull_N_tab.mean.inv.lambda[3,6]<-mean(Weibull_N_rmse_values_c_lambda_mean_inv)


########################################## Inverse Mean Lambda #############################################

Weibull_N_param_list5 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_N_inner_list5 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data5 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list5[[j]] <- Weibull_N_matrix_data5
  }
  Weibull_N_param_list5[[i]] <- Weibull_N_inner_list5
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- Weibull_N_data.list[[i]] 
  data$inv.trait <- 1/(data$trait)
  group_means <- tapply(data$inv.trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean<- data.frame(
    T=c(rep(13,10),rep(17,10),rep(21,10),rep(25,10),rep(29,10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data<-df.mean
  trait <- data$trait
  
  N.obs <- length(trait)
  temp <- data$T
  #Model
  sink("weibull_n_inv_mean_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  sig2~dexp(325)
  tau<-1/sig2
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - b.l * temp[i] + c)
    trait[i] ~ dnorm(mu[i],tau) T(0,)
  }
}",file="weibull_n_inv_mean_lambda.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c","sig2")
  inits<-function(){list(
    la = log(0.001),
    b.l = 0.2,
    c = 0.6,
    sig2=2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="weibull_n_inv_mean_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list5[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list5[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list5[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list5[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list5[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list5[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list5[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list5[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list5[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list5[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list5[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list5[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list5[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list5[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list5[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list5[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_N_param_list5[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_N_param_list5[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_N_param_list5[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_N_param_list5[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf5<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                   med.value=median(value),
                                                   lower.hdi=hdi(value)[1],
                                                   upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf5$true_curve <- true_curve(xseq)
  Weibull_N_plot5_mean<-ggplot(Weibull_N_mean.xdf5, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot5_med<-ggplot(Weibull_N_mean.xdf5, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot5_mean,Weibull_N_plot5_med,nrow=1)

Weibull_N_hist_list5 <- vector("list", length = 5)

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list5 <- lapply(Weibull_N_param_list5[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector5 <- unlist(Weibull_N_first_column_list5)
  x<-seq(min(Weibull_N_combined_vector5)-1,max(Weibull_N_combined_vector5)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list5[[i]] <- hist(Weibull_N_combined_vector5, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                            xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 325), col="purple", lty=2, lwd=2)
  }
}

Weibull_N_proportions_list_a5<- lapply(Weibull_N_param_list5[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a5))

Weibull_N_proportions_list_b5<- lapply(Weibull_N_param_list5[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b5))

Weibull_N_proportions_list_c5<- lapply(Weibull_N_param_list5[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c5))

Weibull_N_calculate_rmse_a_lambda_inv_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_lambda_inv_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_lambda_inv_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_lambda_inv_mean <- sapply(Weibull_N_param_list5[[1]], Weibull_N_calculate_rmse_a_lambda_inv_mean)
Weibull_N_rmse_values_b_lambda_inv_mean <- sapply(Weibull_N_param_list5[[2]], Weibull_N_calculate_rmse_b_lambda_inv_mean)
Weibull_N_rmse_values_c_lambda_inv_mean <- sapply(Weibull_N_param_list5[[3]], Weibull_N_calculate_rmse_c_lambda_inv_mean)

Weibull_N_tab.inv.mean.lambda<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.inv.mean.lambda)<-c("a","b","c")
colnames(Weibull_N_tab.inv.mean.lambda)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_N_tab.inv.mean.lambda[1,1]<-true.exp.a
Weibull_N_tab.inv.mean.lambda[2,1]<-true.exp.b
Weibull_N_tab.inv.mean.lambda[3,1]<-true.exp.c
Weibull_N_tab.inv.mean.lambda[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list5[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.inv.mean.lambda[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list5[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.inv.mean.lambda[3,2]<-mean(unlist(lapply(Weibull_N_param_list5[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.inv.mean.lambda[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list5[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.inv.mean.lambda[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list5[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.inv.mean.lambda[3,3]<-hdi(unlist(lapply(Weibull_N_param_list5[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.inv.mean.lambda[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list5[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.inv.mean.lambda[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list5[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.inv.mean.lambda[3,4]<-hdi(unlist(lapply(Weibull_N_param_list5[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.inv.mean.lambda[1,5]<-mean(unlist(Weibull_N_proportions_list_a5))
Weibull_N_tab.inv.mean.lambda[2,5]<-mean(unlist(Weibull_N_proportions_list_b5))
Weibull_N_tab.inv.mean.lambda[3,5]<-mean(unlist(Weibull_N_proportions_list_c5))
Weibull_N_tab.inv.mean.lambda[1,6]<-mean(Weibull_N_rmse_values_a_lambda_inv_mean)
Weibull_N_tab.inv.mean.lambda[2,6]<-mean(Weibull_N_rmse_values_b_lambda_inv_mean)
Weibull_N_tab.inv.mean.lambda[3,6]<-mean(Weibull_N_rmse_values_c_lambda_inv_mean)






########################################## Weibull ###########################################
Weibull_N_param_list6 <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_N_inner_list6 <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data6 <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list6[[j]] <- Weibull_N_matrix_data6
  }
  Weibull_N_param_list6[[i]] <- Weibull_N_inner_list6
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("Weibull_n_weibull.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  shape~dexp(0.001)
  #Likelihood
  for (i in 1:N.obs){
    trait[i] ~ dgen.gamma(1, lambda[i], shape) #
    lambda[i] <- (log(2))^(1/shape)/(med[i])
    med[i] <- exp(exp(la) * temp[i]^2 - b.l * temp[i] + c)
  }
  }",file="Weibull_n_weibull.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c","shape")
  inits<-function(){list(
    la = log(0.001),
    b.l = 0.2,
    c = 0.6,
    shape=3
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data <- Weibull_N_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="Weibull_n_weibull.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list6[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list6[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list6[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list6[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list6[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list6[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list6[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list6[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list6[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list6[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list6[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list6[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list6[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list6[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list6[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list6[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # shape
    Weibull_N_param_list6[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_N_param_list6[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_N_param_list6[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_N_param_list6[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf6<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                           med.value=median(value),
                                                           lower.hdi=hdi(value)[1],
                                                           upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf6$true_curve <- true_curve(xseq)
  Weibull_N_plot6_mean<-ggplot(Weibull_N_mean.xdf6, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot6_med<-ggplot(Weibull_N_mean.xdf6, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot6_mean,Weibull_N_plot6_med,nrow=1)

Weibull_N_hist_list6 <- vector("list", length = 5)

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list6 <- lapply(Weibull_N_param_list6[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector6 <- unlist(Weibull_N_first_column_list6)
  x<-seq(min(Weibull_N_combined_vector6)-1,max(Weibull_N_combined_vector6)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list6[[i]] <- hist(Weibull_N_combined_vector6, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                    xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 0.001), col="purple", lty=2, lwd=2)
  }
}


Weibull_N_proportions_list_a6<- lapply(Weibull_N_param_list6[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a6))

Weibull_N_proportions_list_b6<- lapply(Weibull_N_param_list6[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b6))

Weibull_N_proportions_list_c6<- lapply(Weibull_N_param_list6[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c6))

Weibull_N_calculate_rmse_a_wb <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_wb <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_wb <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_wb <- sapply(Weibull_N_param_list6[[1]], Weibull_N_calculate_rmse_a_wb)
Weibull_N_rmse_values_b_wb <- sapply(Weibull_N_param_list6[[2]], Weibull_N_calculate_rmse_b_wb)
Weibull_N_rmse_values_c_wb <- sapply(Weibull_N_param_list6[[3]], Weibull_N_calculate_rmse_c_wb)

Weibull_N_tab.wb<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.wb)<-c("a","b","c")
colnames(Weibull_N_tab.wb)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_N_tab.wb[1,1]<-true.exp.a
Weibull_N_tab.wb[2,1]<-true.exp.b
Weibull_N_tab.wb[3,1]<-true.exp.c
Weibull_N_tab.wb[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list6[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.wb[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list6[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.wb[3,2]<-mean(unlist(lapply(Weibull_N_param_list6[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.wb[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list6[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.wb[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list6[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.wb[3,3]<-hdi(unlist(lapply(Weibull_N_param_list6[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.wb[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list6[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.wb[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list6[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.wb[3,4]<-hdi(unlist(lapply(Weibull_N_param_list6[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.wb[1,5]<-mean(unlist(Weibull_N_proportions_list_a6))
Weibull_N_tab.wb[2,5]<-mean(unlist(Weibull_N_proportions_list_b6))
Weibull_N_tab.wb[3,5]<-mean(unlist(Weibull_N_proportions_list_c6))
Weibull_N_tab.wb[1,6]<-mean(Weibull_N_rmse_values_a_wb)
Weibull_N_tab.wb[2,6]<-mean(Weibull_N_rmse_values_b_wb)
Weibull_N_tab.wb[3,6]<-mean(Weibull_N_rmse_values_c_wb)

########################################## Inverse Lambda NOT TRUNCATED #############################################

Weibull_N_param_list3_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_N_inner_list3_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data3_notrunc <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list3_notrunc[[j]] <- Weibull_N_matrix_data3_notrunc
  }
  Weibull_N_param_list3_notrunc[[i]] <- Weibull_N_inner_list3_notrunc
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  #Model
  sink("weibull_n_inv_lambda_nt.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  tau~dexp(0.1)
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - b.l * temp[i] + c)
    trait[i] ~ dnorm(mu[i],tau)
  }
}",file="weibull_n_inv_lambda_nt.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c","tau")
  inits<-function(){list(
    la = log(0.001),
    b.l = 0.2,
    c = 0.6,
    tau=0.2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data <- N_data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="weibull_n_inv_lambda_nt.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list3_notrunc[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list3_notrunc[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list3_notrunc[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list3_notrunc[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list3_notrunc[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list3_notrunc[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list3_notrunc[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list3_notrunc[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list3_notrunc[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list3_notrunc[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list3_notrunc[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list3_notrunc[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list3_notrunc[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list3_notrunc[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list3_notrunc[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list3_notrunc[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_N_param_list3_notrunc[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_N_param_list3_notrunc[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_N_param_list3_notrunc[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_N_param_list3_notrunc[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf3_NT<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                      med.value=median(value),
                                                      lower.hdi=hdi(value)[1],
                                                      upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf3_NT$true_curve <- true_curve(xseq)
  Weibull_N_plot3_mean_NT<-ggplot(Weibull_N_mean.xdf3_NT, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot3_med_NT<-ggplot(Weibull_N_mean.xdf3_NT, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot3_mean_NT,Weibull_N_plot3_med_NT,nrow=1)

Weibull_N_hist_list3_notrunc <- vector("list", length = 5)

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list3_notrunc <- lapply(Weibull_N_param_list3_notrunc[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector3_notrunc <- unlist(Weibull_N_first_column_list3_notrunc)
  x<-seq(min(Weibull_N_combined_vector3_notrunc)-1,max(Weibull_N_combined_vector3_notrunc)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list3_notrunc[[i]] <- hist(Weibull_N_combined_vector3_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                    xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 0.1), col="purple", lty=2, lwd=2)
  }
}


Weibull_N_proportions_list_a3_notrunc<- lapply(Weibull_N_param_list3_notrunc[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a3_notrunc))

Weibull_N_proportions_list_b3_notrunc<- lapply(Weibull_N_param_list3_notrunc[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b3_notrunc))

Weibull_N_proportions_list_c3_notrunc<- lapply(Weibull_N_param_list3_notrunc[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c3_notrunc))

Weibull_N_calculate_rmse_a_lambda_inv_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_lambda_inv_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_lambda_inv_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_lambda_inv_notrunc <- sapply(Weibull_N_param_list3_notrunc[[1]], Weibull_N_calculate_rmse_a_lambda_inv_notrunc)
Weibull_N_rmse_values_b_lambda_inv_notrunc <- sapply(Weibull_N_param_list3_notrunc[[2]], Weibull_N_calculate_rmse_b_lambda_inv_notrunc)
Weibull_N_rmse_values_c_lambda_inv_notrunc <- sapply(Weibull_N_param_list3_notrunc[[3]], Weibull_N_calculate_rmse_c_lambda_inv_notrunc)

Weibull_N_tab.inv.lambda_notrunc<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.inv.lambda_notrunc)<-c("a","b","c")
colnames(Weibull_N_tab.inv.lambda_notrunc)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_N_tab.inv.lambda_notrunc[1,1]<-true.exp.a
Weibull_N_tab.inv.lambda_notrunc[2,1]<-true.exp.b
Weibull_N_tab.inv.lambda_notrunc[3,1]<-true.exp.c
Weibull_N_tab.inv.lambda_notrunc[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list3_notrunc[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.inv.lambda_notrunc[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list3_notrunc[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.inv.lambda_notrunc[3,2]<-mean(unlist(lapply(Weibull_N_param_list3_notrunc[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.inv.lambda_notrunc[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list3_notrunc[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.inv.lambda_notrunc[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list3_notrunc[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.inv.lambda_notrunc[3,3]<-hdi(unlist(lapply(Weibull_N_param_list3_notrunc[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.inv.lambda_notrunc[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list3_notrunc[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.inv.lambda_notrunc[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list3_notrunc[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.inv.lambda_notrunc[3,4]<-hdi(unlist(lapply(Weibull_N_param_list3_notrunc[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.inv.lambda_notrunc[1,5]<-mean(unlist(Weibull_N_proportions_list_a3_notrunc))
Weibull_N_tab.inv.lambda_notrunc[2,5]<-mean(unlist(Weibull_N_proportions_list_b3_notrunc))
Weibull_N_tab.inv.lambda_notrunc[3,5]<-mean(unlist(Weibull_N_proportions_list_c3_notrunc))
Weibull_N_tab.inv.lambda_notrunc[1,6]<-mean(Weibull_N_rmse_values_a_lambda_inv_notrunc)
Weibull_N_tab.inv.lambda_notrunc[2,6]<-mean(Weibull_N_rmse_values_b_lambda_inv_notrunc)
Weibull_N_tab.inv.lambda_notrunc[3,6]<-mean(Weibull_N_rmse_values_c_lambda_inv_notrunc)



########################################## Inverse Mean Lambda NOT TRUNCATED #############################################

Weibull_N_param_list5_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_N_inner_list5_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data5_notrunc <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list5_notrunc[[j]] <- Weibull_N_matrix_data5_notrunc
  }
  Weibull_N_param_list5_notrunc[[i]] <- Weibull_N_inner_list5_notrunc
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- Weibull_N_data.list[[i]] 
  data$inv.trait <- 1/(data$trait)
  group_means <- tapply(data$inv.trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean<- data.frame(
    T=c(rep(13,10),rep(17,10),rep(21,10),rep(25,10),rep(29,10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data<-df.mean
  trait <- data$trait
  
  N.obs <- length(trait)
  temp <- data$T
  #Model
  sink("weibull_n_inv_mean_lambda_nt.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  tau~dexp(0.1)
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - b.l * temp[i] + c)
    trait[i] ~ dnorm(mu[i],tau)
  }
}",file="weibull_n_inv_mean_lambda_nt.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c","tau")
  inits<-function(){list(
    la = log(0.001),
    b.l = 0.2,
    c = 0.6,
    tau=0.2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="weibull_n_inv_mean_lambda_nt.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list5_notrunc[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list5_notrunc[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list5_notrunc[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list5_notrunc[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list5_notrunc[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list5_notrunc[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list5_notrunc[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list5_notrunc[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list5_notrunc[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list5_notrunc[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list5_notrunc[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list5_notrunc[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list5_notrunc[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list5_notrunc[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list5_notrunc[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list5_notrunc[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_N_param_list5_notrunc[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_N_param_list5_notrunc[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_N_param_list5_notrunc[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_N_param_list5_notrunc[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf5_NT<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                      med.value=median(value),
                                                      lower.hdi=hdi(value)[1],
                                                      upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf5_NT$true_curve <- true_curve(xseq)
  Weibull_N_plot5_mean_NT<-ggplot(Weibull_N_mean.xdf5_NT, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot5_med_NT<-ggplot(Weibull_N_mean.xdf5_NT, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot5_mean_NT,Weibull_N_plot5_med_NT,nrow=1)

Weibull_N_hist_list5_notrunc <- vector("list", length = 5)

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list5_notrunc <- lapply(Weibull_N_param_list5_notrunc[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector5_notrunc <- unlist(Weibull_N_first_column_list5_notrunc)
  x<-seq(min(Weibull_N_combined_vector5_notrunc)-1,max(Weibull_N_combined_vector5_notrunc)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list5_notrunc[[i]] <- hist(Weibull_N_combined_vector5_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                    xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 0.1), col="purple", lty=2, lwd=2)
  }
}

Weibull_N_proportions_list_a5_notrunc<- lapply(Weibull_N_param_list5_notrunc[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a5_notrunc))

Weibull_N_proportions_list_b5_notrunc<- lapply(Weibull_N_param_list5_notrunc[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b5_notrunc))

Weibull_N_proportions_list_c5_notrunc<- lapply(Weibull_N_param_list5_notrunc[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c5_notrunc))

Weibull_N_calculate_rmse_a_lambda_inv_mean_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_lambda_inv_mean_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_lambda_inv_mean_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_lambda_inv_mean_notrunc <- sapply(Weibull_N_param_list5_notrunc[[1]], Weibull_N_calculate_rmse_a_lambda_inv_mean_notrunc)
Weibull_N_rmse_values_b_lambda_inv_mean_notrunc <- sapply(Weibull_N_param_list5_notrunc[[2]], Weibull_N_calculate_rmse_b_lambda_inv_mean_notrunc)
Weibull_N_rmse_values_c_lambda_inv_mean_notrunc <- sapply(Weibull_N_param_list5_notrunc[[3]], Weibull_N_calculate_rmse_c_lambda_inv_mean_notrunc)

Weibull_N_tab.inv.mean.lambda_notrunc<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.inv.mean.lambda_notrunc)<-c("a","b","c")
colnames(Weibull_N_tab.inv.mean.lambda_notrunc)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_N_tab.inv.mean.lambda_notrunc[1,1]<-true.exp.a
Weibull_N_tab.inv.mean.lambda_notrunc[2,1]<-true.exp.b
Weibull_N_tab.inv.mean.lambda_notrunc[3,1]<-true.exp.c
Weibull_N_tab.inv.mean.lambda_notrunc[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list5_notrunc[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.inv.mean.lambda_notrunc[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list5_notrunc[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.inv.mean.lambda_notrunc[3,2]<-mean(unlist(lapply(Weibull_N_param_list5_notrunc[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.inv.mean.lambda_notrunc[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list5_notrunc[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.inv.mean.lambda_notrunc[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list5_notrunc[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.inv.mean.lambda_notrunc[3,3]<-hdi(unlist(lapply(Weibull_N_param_list5_notrunc[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.inv.mean.lambda_notrunc[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list5_notrunc[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.inv.mean.lambda_notrunc[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list5_notrunc[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.inv.mean.lambda_notrunc[3,4]<-hdi(unlist(lapply(Weibull_N_param_list5_notrunc[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.inv.mean.lambda_notrunc[1,5]<-mean(unlist(Weibull_N_proportions_list_a5_notrunc))
Weibull_N_tab.inv.mean.lambda_notrunc[2,5]<-mean(unlist(Weibull_N_proportions_list_b5_notrunc))
Weibull_N_tab.inv.mean.lambda_notrunc[3,5]<-mean(unlist(Weibull_N_proportions_list_c5_notrunc))
Weibull_N_tab.inv.mean.lambda_notrunc[1,6]<-mean(Weibull_N_rmse_values_a_lambda_inv_mean_notrunc)
Weibull_N_tab.inv.mean.lambda_notrunc[2,6]<-mean(Weibull_N_rmse_values_b_lambda_inv_mean_notrunc)
Weibull_N_tab.inv.mean.lambda_notrunc[3,6]<-mean(Weibull_N_rmse_values_c_lambda_inv_mean_notrunc)


########################################## Mean Inverse Lambda NOT TRUNCATED #############################################

Weibull_N_param_list4_notrunc <- vector("list", length = 5)
for (i in 1:5) {
  # Initialize the inner list
  Weibull_N_inner_list4_notrunc <- vector("list", length = b)
  # Populate the inner list
  for (j in 1:b) {
    # Initialize a 5x4 matrix
    Weibull_N_matrix_data4_notrunc <- matrix(0, nrow = 5, ncol = 4)
    Weibull_N_inner_list4_notrunc[[j]] <- Weibull_N_matrix_data4_notrunc
  }
  Weibull_N_param_list4_notrunc[[i]] <- Weibull_N_inner_list4_notrunc
}
for(i in 1:b){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- Weibull_N_data.list[[i]] 
  group_means <- tapply(data$trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean<- data.frame(
    T=c(rep(13,10),rep(17,10),rep(21,10),rep(25,10),rep(29,10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data<-df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  #Model
  sink("weibull_n_mean_inv_lambda_nt.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b.l~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  tau~dexp(0.1)
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- (exp(la) * temp[i]^2 - b.l * temp[i] + c)
    trait[i] ~ dnorm(mu[i],tau)
  }
}",file="weibull_n_mean_inv_lambda_nt.txt")
  
  # Settings
  parameters <- c("la", "b.l", "c","tau")
  inits<-function(){list(
    la = log(0.001),
    b.l = 0.2,
    c = 0.6,
    tau=0.2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="weibull_n_mean_inv_lambda_nt.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_N_param_list4_notrunc[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_N_param_list4_notrunc[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_N_param_list4_notrunc[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_N_param_list4_notrunc[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_N_param_list4_notrunc[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_N_param_list4_notrunc[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_N_param_list4_notrunc[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_N_param_list4_notrunc[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_N_param_list4_notrunc[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_N_param_list4_notrunc[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_N_param_list4_notrunc[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_N_param_list4_notrunc[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_N_param_list4_notrunc[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_N_param_list4_notrunc[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_N_param_list4_notrunc[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_N_param_list4_notrunc[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_N_param_list4_notrunc[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_N_param_list4_notrunc[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_N_param_list4_notrunc[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_N_param_list4_notrunc[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  df.new = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf <- expand.grid(point = xseq, iteration = 1:nrow(df.new))
  
  # Apply the quadratic equation for each combination
  xdf$value <- apply(xdf, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    1/(exp(exp(df.new[i, 4]) * location^2 - df.new[i, 1] * location + df.new[i, 2]))
  })
  Weibull_N_mean.xdf4_NT<-xdf |> group_by(point) |> summarize(avg.value=mean(value),
                                                      med.value=median(value),
                                                      lower.hdi=hdi(value)[1],
                                                      upper.hdi=hdi(value)[2])
  true_curve <- function(x) 1/(exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  Weibull_N_mean.xdf4_NT$true_curve <- true_curve(xseq)
  Weibull_N_plot4_mean_NT<-ggplot(Weibull_N_mean.xdf4_NT, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
  Weibull_N_plot4_med_NT<-ggplot(Weibull_N_mean.xdf4_NT, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "greenyellow", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="blue2",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Response")+
    theme_minimal()
}

gridExtra::grid.arrange(Weibull_N_plot4_mean_NT,Weibull_N_plot4_med_NT,nrow=1)

Weibull_N_hist_list4_notrunc <- vector("list", length = 4)

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_N_first_column_list4_notrunc <- lapply(Weibull_N_param_list4_notrunc[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_N_combined_vector4_notrunc <- unlist(Weibull_N_first_column_list4_notrunc)
  x<-seq(min(Weibull_N_combined_vector4_notrunc)-1,max(Weibull_N_combined_vector4_notrunc)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_N_hist_list4_notrunc[[i]] <- hist(Weibull_N_combined_vector4_notrunc, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                    xlab = "Values", col = "lightpink2", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "red", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x,0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 0.1), col="purple", lty=2, lwd=2)
  }
}


Weibull_N_proportions_list_a4_notrunc<- lapply(Weibull_N_param_list4_notrunc[[1]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(log(true.exp.a) > col3 & log(true.exp.a) < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_a4_notrunc))

Weibull_N_proportions_list_b4_notrunc<- lapply(Weibull_N_param_list4_notrunc[[2]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.b > col3 & true.exp.b < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_b4_notrunc))

Weibull_N_proportions_list_c4_notrunc<- lapply(Weibull_N_param_list4_notrunc[[3]], function(matrix) {
  # Extract columns
  col2 <- matrix[, 1]
  col3 <- matrix[, 3]
  col4 <- matrix[, 4]
  
  # Calculate proportion
  proportion <- mean(true.exp.c > col3 & true.exp.c < col4)
  return(proportion)
})

# Display the proportions
#print(proportions_list_median)
mean(unlist(Weibull_N_proportions_list_c4_notrunc))

Weibull_N_calculate_rmse_a_lambda_mean_inv_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_b_lambda_mean_inv_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_N_calculate_rmse_c_lambda_mean_inv_notrunc <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_N_rmse_values_a_lambda_mean_inv_notrunc <- sapply(Weibull_N_param_list4_notrunc[[1]], Weibull_N_calculate_rmse_a_lambda_mean_inv_notrunc)
Weibull_N_rmse_values_b_lambda_mean_inv_notrunc <- sapply(Weibull_N_param_list4_notrunc[[2]], Weibull_N_calculate_rmse_b_lambda_mean_inv_notrunc)
Weibull_N_rmse_values_c_lambda_mean_inv_notrunc <- sapply(Weibull_N_param_list4_notrunc[[3]], Weibull_N_calculate_rmse_c_lambda_mean_inv_notrunc)

Weibull_N_tab.mean.inv.lambda_notrunc<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_N_tab.mean.inv.lambda_notrunc)<-c("a","b","c")
colnames(Weibull_N_tab.mean.inv.lambda_notrunc)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage", "RMSE")
Weibull_N_tab.mean.inv.lambda_notrunc[1,1]<-true.exp.a
Weibull_N_tab.mean.inv.lambda_notrunc[2,1]<-true.exp.b
Weibull_N_tab.mean.inv.lambda_notrunc[3,1]<-true.exp.c
Weibull_N_tab.mean.inv.lambda_notrunc[1,2]<-exp(mean(unlist(lapply(Weibull_N_param_list4_notrunc[[1]], function(matrix) matrix[, 1]))))
Weibull_N_tab.mean.inv.lambda_notrunc[2,2]<-exp(mean(unlist(lapply(Weibull_N_param_list4_notrunc[[2]], function(matrix) matrix[, 1]))))
Weibull_N_tab.mean.inv.lambda_notrunc[3,2]<-mean(unlist(lapply(Weibull_N_param_list4_notrunc[[3]], function(matrix) matrix[, 1])))
Weibull_N_tab.mean.inv.lambda_notrunc[1,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list4_notrunc[[1]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.mean.inv.lambda_notrunc[2,3]<-exp(hdi(unlist(lapply(Weibull_N_param_list4_notrunc[[2]], function(matrix) matrix[, 1])))[1])
Weibull_N_tab.mean.inv.lambda_notrunc[3,3]<-hdi(unlist(lapply(Weibull_N_param_list4_notrunc[[3]], function(matrix) matrix[, 1])))[1]
Weibull_N_tab.mean.inv.lambda_notrunc[1,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list4_notrunc[[1]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.mean.inv.lambda_notrunc[2,4]<-exp(hdi(unlist(lapply(Weibull_N_param_list4_notrunc[[2]], function(matrix) matrix[, 1])))[2])
Weibull_N_tab.mean.inv.lambda_notrunc[3,4]<-hdi(unlist(lapply(Weibull_N_param_list4_notrunc[[3]], function(matrix) matrix[, 1])))[2]
Weibull_N_tab.mean.inv.lambda_notrunc[1,5]<-mean(unlist(Weibull_N_proportions_list_a4_notrunc))
Weibull_N_tab.mean.inv.lambda_notrunc[2,5]<-mean(unlist(Weibull_N_proportions_list_b4_notrunc))
Weibull_N_tab.mean.inv.lambda_notrunc[3,5]<-mean(unlist(Weibull_N_proportions_list_c4_notrunc))
Weibull_N_tab.mean.inv.lambda_notrunc[1,6]<-mean(Weibull_N_rmse_values_a_lambda_mean_inv_notrunc)
Weibull_N_tab.mean.inv.lambda_notrunc[2,6]<-mean(Weibull_N_rmse_values_b_lambda_mean_inv_notrunc)
Weibull_N_tab.mean.inv.lambda_notrunc[3,6]<-mean(Weibull_N_rmse_values_c_lambda_mean_inv_notrunc)

