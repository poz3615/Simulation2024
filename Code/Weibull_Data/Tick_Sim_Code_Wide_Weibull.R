library(rjags)
library(dplyr)
library(ggplot2)

true.exp.a<- 0.00826
true.exp.b<- 0.37
true.exp.c<- 1.406
# Making quadratic function
my.med<-function(x,a=0.00826,b=0.37,c=1.406){
  exp(-(a*x^2-b*x+c))
}
# Set reasonable shape parameter
sh.est<-3
# Inverse lambda function
my.lambda<-function(x,a=0.00826,b=0.37,c=1.406,sh.est=3){ 
  my.med(x,a,b,c)/(log(2))^(1/sh.est) ## fx0 is the median here
}

Weibull_data.list<-list()
set.seed(52)
for(i in 1:100){
  T1W<-rweibull(n=100,scale=my.inv.lambda(10),shape=sh.est)
  T2W<-rweibull(n=100,scale=my.inv.lambda(17),shape=sh.est)
  T3W<-rweibull(n=100,scale=my.inv.lambda(21),shape=sh.est)
  T4W<-rweibull(n=100,scale=my.inv.lambda(25),shape=sh.est)
  T5W<-rweibull(n=100,scale=my.inv.lambda(32),shape=sh.est)
  Weibull_data.list[[i]]<-data.frame(
    T=c(rep(10,100),rep(17,100),rep(21,100),rep(25,100),rep(32,100)),
    trait=c(T1W,T2W,T3W,T4W,T5W)
  )
}

w<-length(Weibull_data.list)

epsilon<-0.01

Weibull_true_values<-c(log(true.exp.a), true.exp.b, true.exp.c, NA,NA)

########################################## Lambda #############################################

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
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("Weibull_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp(exp(la) * temp[i]^2 - b * temp[i] + c)
    trait[i] ~ dexp(mu[i])
  }
  }",file="Weibull_lambda.txt")
  
  # Settings
  parameters <- c("la", "b", "c")
  inits<-function(){list(
    la = log(0.005),
    b = 0.4,
    c = 1.5
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data.raw.w<-Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="Weibull_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_param_list1[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_param_list1[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_param_list1[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_param_list1[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_param_list1[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_param_list1[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_param_list1[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_param_list1[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_param_list1[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_param_list1[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_param_list1[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_param_list1[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_param_list1[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_param_list1[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_param_list1[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_param_list1[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
  }
  
  df.new1_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf1_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new1_W))
  
  xdf1_W$value <- apply(xdf1_W, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    exp(-(exp(df.new1_W[i, 4]) * location^2 - df.new1_W[i, 1] * location + df.new1_W[i, 2]))
  })
  xdf1_W<- xdf1_W|> mutate(trunc.eval=ifelse(xdf1_W$value>0,xdf1_W$value,epsilon),
                       trunc.inv=log(2)/trunc.eval)
  mean.xdf1_W<-xdf1_W |> group_by(point) |> summarize(avg.value.inv=mean(trunc.inv),
                                                  avg.value=mean(trunc.eval),
                                                  med.value.inv=median(trunc.inv),
                                                  med.value=median(trunc.eval),
                                                  lower.hdi.inv=hdi(trunc.inv)[1],
                                                  upper.hdi.inv=hdi(trunc.inv)[2],
                                                  lower.hdi=hdi(trunc.eval)[1],
                                                  upper.hdi=hdi(trunc.eval)[2])
  true_curve.inv <- function(x) log(2)/exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  true_curve <- function(x) exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  mean.xdf1_W$true_curve <- true_curve(xseq)
  mean.xdf1_W$true_curve.inv <- true_curve.inv(xseq)
  plot1_mean_W<-ggplot(mean.xdf1_W, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot1_med_W<-ggplot(mean.xdf1_W, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot1_mean.inv_W<-ggplot(mean.xdf1_W, aes(x = point)) +
    geom_line(aes(y = avg.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortlity Rate")+
    theme_minimal()
  plot1_med.inv_W<-ggplot(mean.xdf1_W, aes(x = point)) +
    geom_line(aes(y = med.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
}

gridExtra::grid.arrange(plot1_mean_W,plot1_med_W,nrow=1)
gridExtra::grid.arrange(plot1_mean.inv_W,plot1_med.inv_W,nrow=1)

Weibull_hist_list1 <- vector("list", length = 4)
par(mfrow=c(2,2))

# Loop through indices 1 to 4
for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_first_column_list1 <- lapply(Weibull_param_list1[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_combined_vector1 <- unlist(Weibull_first_column_list1)
  x<-seq(min(Weibull_combined_vector1)-50,max(Weibull_combined_vector1)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_hist_list1[[i]] <- hist(Weibull_combined_vector1, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "mintcream", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "darkmagenta", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x, 0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 0.1), col="purple", lty=2, lwd=2)
  }
}


Weibull_proportions_list_a1<- lapply(Weibull_param_list1[[1]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_a1))

Weibull_proportions_list_b1<- lapply(Weibull_param_list1[[2]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_b1))

Weibull_proportions_list_c1<- lapply(Weibull_param_list1[[3]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_c1))

Weibull_calculate_rmse_a_lambda <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_b_lambda <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_c_lambda <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda <- sapply(Weibull_param_list1[[1]], Weibull_calculate_rmse_a_lambda)
Weibull_rmse_values_b_lambda <- sapply(Weibull_param_list1[[2]], Weibull_calculate_rmse_b_lambda)
Weibull_rmse_values_c_lambda <- sapply(Weibull_param_list1[[3]], Weibull_calculate_rmse_c_lambda)

Weibull_tab.lambda<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_tab.lambda)<-c("a","b","c")
colnames(Weibull_tab.lambda)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_tab.lambda[1,1]<-true.exp.a
Weibull_tab.lambda[2,1]<-true.exp.b
Weibull_tab.lambda[3,1]<-true.exp.c
Weibull_tab.lambda[1,2]<-exp(mean(unlist(lapply(Weibull_param_list1[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.lambda[2,2]<-(mean(unlist(lapply(Weibull_param_list1[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.lambda[3,2]<-mean(unlist(lapply(Weibull_param_list1[[3]], function(matrix) matrix[, 1])))
Weibull_tab.lambda[1,3]<-exp(hdi(unlist(lapply(Weibull_param_list1[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.lambda[2,3]<-(hdi(unlist(lapply(Weibull_param_list1[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.lambda[3,3]<-hdi(unlist(lapply(Weibull_param_list1[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.lambda[1,4]<-exp(hdi(unlist(lapply(Weibull_param_list1[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.lambda[2,4]<-(hdi(unlist(lapply(Weibull_param_list1[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.lambda[3,4]<-hdi(unlist(lapply(Weibull_param_list1[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.lambda[1,5]<-mean(unlist(Weibull_proportions_list_a1))
Weibull_tab.lambda[2,5]<-mean(unlist(Weibull_proportions_list_b1))
Weibull_tab.lambda[3,5]<-mean(unlist(Weibull_proportions_list_c1))
Weibull_tab.lambda[1,6]<-mean(Weibull_rmse_values_a_lambda)
Weibull_tab.lambda[2,6]<-mean(Weibull_rmse_values_b_lambda)
Weibull_tab.lambda[3,6]<-mean(Weibull_rmse_values_c_lambda)


########################################## Mean Lambda #############################################

Weibull_param_list2 <- vector("list", length = 4)
for (i in 1:4) {
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
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  group_means <- tapply(Weibull_data.list[[i]]$trait, Weibull_data.list[[i]]$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean<- data.frame(
    T=c(rep(10,10),rep(17,10),rep(21,10),rep(25,10),rep(32,10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  # Model
  sink("Weibull_mean_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp(exp(la) * temp[i]^2 - b * temp[i] + c)
    trait[i] ~ dexp(mu[i])
  }
  }",file="Weibull_mean_lambda.txt")
  
  # Settins
  parameters <- c("la", "b", "c")
  inits<-function(){list(
    la = log(0.005),
    b = 0.4,
    c = 1.5
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data.raw.w<-Weibull_data.list[[i]]
  data <- df.mean 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="Weibull_mean_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_param_list2[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_param_list2[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_param_list2[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_param_list2[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_param_list2[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_param_list2[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_param_list2[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_param_list2[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_param_list2[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_param_list2[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_param_list2[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_param_list2[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_param_list2[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_param_list2[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_param_list2[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_param_list2[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
  }
  
  df.new2_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf2_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new2_W))
  
  # Apply the quadratic equation for each combination
  xdf2_W$value <- apply(xdf2_W, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    exp(-(exp(df.new2_W[i, 4]) * location^2 - df.new2_W[i, 1] * location + df.new2_W[i, 2]))
  })
  xdf2_W<- xdf2_W|> mutate(trunc.eval=ifelse(xdf2_W$value>0,xdf2_W$value,epsilon),
                           trunc.inv=log(2)/trunc.eval)
  mean.xdf2_W<-xdf2_W |> group_by(point) |> summarize(avg.value.inv=mean(trunc.inv),
                                                      avg.value=mean(trunc.eval),
                                                      med.value.inv=median(trunc.inv),
                                                      med.value=median(trunc.eval),
                                                      lower.hdi.inv=hdi(trunc.inv)[1],
                                                      upper.hdi.inv=hdi(trunc.inv)[2],
                                                      lower.hdi=hdi(trunc.eval)[1],
                                                      upper.hdi=hdi(trunc.eval)[2])
  true_curve.inv <- function(x) log(2)/exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  true_curve <- function(x) exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  mean.xdf2_W$true_curve <- true_curve(xseq)
  mean.xdf2_W$true_curve.inv <- true_curve.inv(xseq)
  plot2_mean_W<-ggplot(mean.xdf2_W, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot2_med_W<-ggplot(mean.xdf2_W, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot2_mean.inv_W<-ggplot(mean.xdf2_W, aes(x = point)) +
    geom_line(aes(y = avg.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
  plot2_med.inv_W<-ggplot(mean.xdf2_W, aes(x = point)) +
    geom_line(aes(y = med.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
}

gridExtra::grid.arrange(plot2_mean_W,plot2_med_W,nrow=1)
gridExtra::grid.arrange(plot2_mean.inv_W,plot2_med.inv_W,nrow=1)


Weibull_hist_list2 <- vector("list", length = 4)
par(mfrow=c(2,2))
# Loop through indices 1 to 4
for (i in 1:4) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_first_column_list2 <- lapply(Weibull_param_list2[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_combined_vector2 <- unlist(Weibull_first_column_list2)
  x<-seq(min(Weibull_combined_vector2)-50,max(Weibull_combined_vector2)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_hist_list2[[i]] <- hist(Weibull_combined_vector2, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "mintcream", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "darkmagenta", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x, 0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 0.1), col="purple", lty=2, lwd=2)
  }
}


Weibull_proportions_list_a2<- lapply(Weibull_param_list2[[1]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_a2))

Weibull_proportions_list_b2<- lapply(Weibull_param_list2[[2]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_b2))

Weibull_proportions_list_c2<- lapply(Weibull_param_list2[[3]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_c2))

Weibull_calculate_rmse_a_lambda_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_b_lambda_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_c_lambda_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda_mean <- sapply(Weibull_param_list2[[1]], Weibull_calculate_rmse_a_lambda_mean)
Weibull_rmse_values_b_lambda_mean <- sapply(Weibull_param_list2[[2]], Weibull_calculate_rmse_b_lambda_mean)
Weibull_rmse_values_c_lambda_mean <- sapply(Weibull_param_list2[[3]], Weibull_calculate_rmse_c_lambda_mean)


Weibull_tab.lambda.mean<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_tab.lambda.mean)<-c("a","b","c")
colnames(Weibull_tab.lambda.mean)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_tab.lambda.mean[1,1]<-true.exp.a
Weibull_tab.lambda.mean[2,1]<-true.exp.b
Weibull_tab.lambda.mean[3,1]<-true.exp.c
Weibull_tab.lambda.mean[1,2]<-exp(mean(unlist(lapply(Weibull_param_list2[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.lambda.mean[2,2]<-(mean(unlist(lapply(Weibull_param_list2[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.lambda.mean[3,2]<-mean(unlist(lapply(Weibull_param_list2[[3]], function(matrix) matrix[, 1])))
Weibull_tab.lambda.mean[1,3]<-exp(hdi(unlist(lapply(Weibull_param_list2[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.lambda.mean[2,3]<-(hdi(unlist(lapply(Weibull_param_list2[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.lambda.mean[3,3]<-hdi(unlist(lapply(Weibull_param_list2[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.lambda.mean[1,4]<-exp(hdi(unlist(lapply(Weibull_param_list2[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.lambda.mean[2,4]<-(hdi(unlist(lapply(Weibull_param_list2[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.lambda.mean[3,4]<-hdi(unlist(lapply(Weibull_param_list2[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.lambda.mean[1,5]<-mean(unlist(Weibull_proportions_list_a2))
Weibull_tab.lambda.mean[2,5]<-mean(unlist(Weibull_proportions_list_b2))
Weibull_tab.lambda.mean[3,5]<-mean(unlist(Weibull_proportions_list_c2))
Weibull_tab.lambda.mean[1,6]<-mean(Weibull_rmse_values_a_lambda_mean)
Weibull_tab.lambda.mean[2,6]<-mean(Weibull_rmse_values_b_lambda_mean)
Weibull_tab.lambda.mean[3,6]<-mean(Weibull_rmse_values_c_lambda_mean)



########################################## Inverse Lambda #############################################

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
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  #Model
  sink("Weibull_inv_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  sig~dexp(15000)
  sig2<-sig^2
  tau<-1/sig2
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp((exp(la) * temp[i]^2 - b * temp[i] + c))
    trait[i] ~ dnorm(mu[i],tau) T(0,)
  }
}",file="Weibull_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "b", "c","sig")
  inits<-function(){list(
    la = log(0.005),
    b = 0.4,
    c = 1.5,
    sig=2
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data.raw.w<-Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="Weibull_inv_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_param_list3[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_param_list3[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_param_list3[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_param_list3[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_param_list3[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_param_list3[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_param_list3[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_param_list3[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_param_list3[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_param_list3[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_param_list3[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_param_list3[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_param_list3[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_param_list3[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_param_list3[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_param_list3[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_param_list3[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_param_list3[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_param_list3[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_param_list3[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  
  
}

df.new3_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
xseq<-seq(from=10,to=35,length.out=250)
xdf3_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new3_W))

# Apply the quadratic equation for each combination
xdf3_W$value.inv <- apply(xdf3_W, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  1/exp(exp(df.new3_W[i, 4]) * location^2 - df.new3_W[i, 1] * location + df.new3_W[i, 2])
})
xdf3_W$value <- apply(xdf3_W, 1, function(row) {
  i <- row["iteration"]
  location <- row["point"]
  exp(exp(df.new3_W[i, 4]) * location^2 - df.new3_W[i, 1] * location + df.new3_W[i, 2])
})
xdf3_W<- xdf3_W|> mutate(trunc.eval=ifelse(xdf3_W$value>0,xdf3_W$value,epsilon),
                         trunc.inv=1/trunc.eval)
mean.xdf3_W<-xdf3_W |> group_by(point) |> summarize(avg.value.inv=mean(trunc.inv),
                                                    avg.value=mean(trunc.eval),
                                                    med.value.inv=median(trunc.inv),
                                                    med.value=median(trunc.eval),
                                                    lower.hdi.inv=hdi(trunc.inv)[1],
                                                    upper.hdi.inv=hdi(trunc.inv)[2],
                                                    lower.hdi=hdi(trunc.eval)[1],
                                                    upper.hdi=hdi(trunc.eval)[2])
true_curve.inv <- function(x) 1/exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c)
true_curve <- function(x) exp(true.exp.a * x^2 - true.exp.b * x + true.exp.c)
mean.xdf3_W$true_curve <- true_curve(xseq)
mean.xdf3_W$true_curve.inv <- true_curve.inv(xseq)
plot3_mean_W<-ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = avg.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
  geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
  geom_point(aes(x=T,y=1/trait),data=data.raw.w)+
  labs(title = "Mean Curve with Interval Bands and True Curve",
       x = "Temperature",
       y = "Mortality Rate")+
  theme_minimal()
plot3_med_W<-ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = med.value), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
  geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
  geom_point(aes(x=T,y=1/trait),data=data.raw.w)+
  labs(title = "Median Curve with Interval Bands and True Curve",
       x = "Temperature",
       y = "Mortality Rate")+
  theme_minimal()
plot3_mean.inv_W<-ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = avg.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
  geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
  geom_point(aes(x=T,y=trait),data=data.raw.w)+
  labs(title = "Mean Curve with Interval Bands and True Curve",
       x = "Temperature",
       y = "Lifetime")+
  theme_minimal()
plot3_med.inv_W<-ggplot(mean.xdf3_W, aes(x = point)) +
  geom_line(aes(y = med.value.inv), color = "black") +
  geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
  geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
  geom_point(aes(x=T,y=trait),data=data.raw.w)+
  labs(title = "Median Curve with Interval Bands and True Curve",
       x = "Temperature",
       y = "Lifetime")+
  theme_minimal()


gridExtra::grid.arrange(plot3_mean_W,plot3_med_W,nrow=1)
gridExtra::grid.arrange(plot3_mean.inv_W,plot3_med.inv_W,nrow=1)



Weibull_hist_list3 <- vector("list", length = 5)
par(mfrow=c(2,3))

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_first_column_list3 <- lapply(Weibull_param_list3[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_combined_vector3 <- unlist(Weibull_first_column_list3)
  x<-seq(min(Weibull_combined_vector3)-1,max(Weibull_combined_vector3)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_hist_list3[[i]] <- hist(Weibull_combined_vector3, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "mintcream", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "darkmagenta", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x, 0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 3500), col="purple", lty=2, lwd=2)
  }
}


Weibull_proportions_list_a3<- lapply(Weibull_param_list3[[1]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_a3))

Weibull_proportions_list_b3<- lapply(Weibull_param_list3[[2]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_b3))

Weibull_proportions_list_c3<- lapply(Weibull_param_list3[[3]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_c3))

Weibull_calculate_rmse_a_lambda_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_b_lambda_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_c_lambda_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda_inv <- sapply(Weibull_param_list3[[1]], Weibull_calculate_rmse_a_lambda_inv)
Weibull_rmse_values_b_lambda_inv <- sapply(Weibull_param_list3[[2]], Weibull_calculate_rmse_b_lambda_inv)
Weibull_rmse_values_c_lambda_inv <- sapply(Weibull_param_list3[[3]], Weibull_calculate_rmse_c_lambda_inv)

Weibull_tab.inv.lambda<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_tab.inv.lambda)<-c("a","b","c")
colnames(Weibull_tab.inv.lambda)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_tab.inv.lambda[1,1]<-true.exp.a
Weibull_tab.inv.lambda[2,1]<-true.exp.b
Weibull_tab.inv.lambda[3,1]<-true.exp.c
Weibull_tab.inv.lambda[1,2]<-exp(mean(unlist(lapply(Weibull_param_list3[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.inv.lambda[2,2]<-(mean(unlist(lapply(Weibull_param_list3[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.inv.lambda[3,2]<-mean(unlist(lapply(Weibull_param_list3[[3]], function(matrix) matrix[, 1])))
Weibull_tab.inv.lambda[1,3]<-exp(hdi(unlist(lapply(Weibull_param_list3[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.inv.lambda[2,3]<-(hdi(unlist(lapply(Weibull_param_list3[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.inv.lambda[3,3]<-hdi(unlist(lapply(Weibull_param_list3[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.inv.lambda[1,4]<-exp(hdi(unlist(lapply(Weibull_param_list3[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.inv.lambda[2,4]<-(hdi(unlist(lapply(Weibull_param_list3[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.inv.lambda[3,4]<-hdi(unlist(lapply(Weibull_param_list3[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.inv.lambda[1,5]<-mean(unlist(Weibull_proportions_list_a3))
Weibull_tab.inv.lambda[2,5]<-mean(unlist(Weibull_proportions_list_b3))
Weibull_tab.inv.lambda[3,5]<-mean(unlist(Weibull_proportions_list_c3))
Weibull_tab.inv.lambda[1,6]<-mean(Weibull_rmse_values_a_lambda_inv)
Weibull_tab.inv.lambda[2,6]<-mean(Weibull_rmse_values_b_lambda_inv)
Weibull_tab.inv.lambda[3,6]<-mean(Weibull_rmse_values_c_lambda_inv)


########################################## Mean Inverse Lambda #############################################

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
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- Weibull_data.list[[i]] 
  group_means <- tapply(data$trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean<- data.frame(
    T=c(rep(10,10),rep(17,10),rep(21,10),rep(25,10),rep(32,10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data.raw.w<-Weibull_data.list[[i]]
  data<-df.mean
  trait <- 1/(data$trait)
  N.obs <- length(trait)
  temp <- data$T
  #Model
  sink("Weibull_mean_inv_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  sig2~dexp(250)
  tau<-1/sig2
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp(exp(la) * temp[i]^2 - b * temp[i] + c)
    trait[i] ~ dnorm(mu[i],tau) T(0,)
  }
}",file="Weibull_mean_inv_lambda.txt")
  
  # Settings
  parameters <- c("la", "b", "c","sig2")
  inits<-function(){list(
    la = log(0.005),
    b = 0.4,
    c = 1.5,
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
                  model.file="Weibull_mean_inv_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_param_list4[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_param_list4[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_param_list4[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_param_list4[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_param_list4[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_param_list4[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_param_list4[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_param_list4[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_param_list4[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_param_list4[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_param_list4[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_param_list4[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_param_list4[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_param_list4[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_param_list4[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_param_list4[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_param_list4[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_param_list4[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_param_list4[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_param_list4[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  
  df.new4_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf4_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new4_W))
  
  # Apply the quadratic equation for each combination
  xdf4_W$value <- apply(xdf4_W, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    exp(-(exp(df.new4_W[i, 4]) * location^2 - df.new4_W[i, 1] * location + df.new4_W[i, 2]))
  })
  mean.xdf4_W<-xdf4_W |> group_by(point) |> summarize(avg.value.inv=mean(1/value),
                                                      avg.value=mean(value),
                                                      med.value.inv=median(1/value),
                                                      med.value=median(value),
                                                      lower.hdi.inv=hdi(1/value)[1],
                                                      upper.hdi.inv=hdi(1/value)[2],
                                                      lower.hdi=hdi(value)[1],
                                                      upper.hdi=hdi(value)[2])
  true_curve.inv <- function(x) 1/exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  true_curve <- function(x) exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  mean.xdf4_W$true_curve <- true_curve(xseq)
  mean.xdf4_W$true_curve.inv <- true_curve.inv(xseq)
  plot4_mean_W<-ggplot(mean.xdf4_W, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot4_med_W<-ggplot(mean.xdf4_W, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot4_mean.inv_W<-ggplot(mean.xdf4_W, aes(x = point)) +
    geom_line(aes(y = avg.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=1/trait),data=data.raw.w)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
  plot4_med.inv_W<-ggplot(mean.xdf4_W, aes(x = point)) +
    geom_line(aes(y = med.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=1/trait),data=data.raw.w)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
}

gridExtra::grid.arrange(plot4_mean_W,plot4_med_W,nrow=1)
gridExtra::grid.arrange(plot4_mean.inv_W,plot4_med.inv_W,nrow=1)



Weibull_hist_list4 <- vector("list", length = 4)
par(mfrow=c(2,3))

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_first_column_list4 <- lapply(Weibull_param_list4[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_combined_vector4 <- unlist(Weibull_first_column_list4)
  x<-seq(min(Weibull_combined_vector4)-1,max(Weibull_combined_vector4)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_hist_list4[[i]] <- hist(Weibull_combined_vector4, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "mintcream", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "darkmagenta", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x, 0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 250), col="purple", lty=2, lwd=2)
  }
}


Weibull_proportions_list_a4<- lapply(Weibull_param_list4[[1]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_a4))

Weibull_proportions_list_b4<- lapply(Weibull_param_list4[[2]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_b4))

Weibull_proportions_list_c4<- lapply(Weibull_param_list4[[3]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_c4))

Weibull_calculate_rmse_a_lambda_mean_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_b_lambda_mean_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_c_lambda_mean_inv <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda_mean_inv <- sapply(Weibull_param_list4[[1]], Weibull_calculate_rmse_a_lambda_mean_inv)
Weibull_rmse_values_b_lambda_mean_inv <- sapply(Weibull_param_list4[[2]], Weibull_calculate_rmse_b_lambda_mean_inv)
Weibull_rmse_values_c_lambda_mean_inv <- sapply(Weibull_param_list4[[3]], Weibull_calculate_rmse_c_lambda_mean_inv)

Weibull_tab.mean.inv.lambda<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_tab.mean.inv.lambda)<-c("a","b","c")
colnames(Weibull_tab.mean.inv.lambda)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage", "RMSE")
Weibull_tab.mean.inv.lambda[1,1]<-true.exp.a
Weibull_tab.mean.inv.lambda[2,1]<-true.exp.b
Weibull_tab.mean.inv.lambda[3,1]<-true.exp.c
Weibull_tab.mean.inv.lambda[1,2]<-exp(mean(unlist(lapply(Weibull_param_list4[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.mean.inv.lambda[2,2]<-(mean(unlist(lapply(Weibull_param_list4[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.mean.inv.lambda[3,2]<-mean(unlist(lapply(Weibull_param_list4[[3]], function(matrix) matrix[, 1])))
Weibull_tab.mean.inv.lambda[1,3]<-exp(hdi(unlist(lapply(Weibull_param_list4[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.mean.inv.lambda[2,3]<-(hdi(unlist(lapply(Weibull_param_list4[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.mean.inv.lambda[3,3]<-hdi(unlist(lapply(Weibull_param_list4[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.mean.inv.lambda[1,4]<-exp(hdi(unlist(lapply(Weibull_param_list4[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.mean.inv.lambda[2,4]<-(hdi(unlist(lapply(Weibull_param_list4[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.mean.inv.lambda[3,4]<-hdi(unlist(lapply(Weibull_param_list4[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.mean.inv.lambda[1,5]<-mean(unlist(Weibull_proportions_list_a4))
Weibull_tab.mean.inv.lambda[2,5]<-mean(unlist(Weibull_proportions_list_b4))
Weibull_tab.mean.inv.lambda[3,5]<-mean(unlist(Weibull_proportions_list_c4))
Weibull_tab.mean.inv.lambda[1,6]<-mean(Weibull_rmse_values_a_lambda_mean_inv)
Weibull_tab.mean.inv.lambda[2,6]<-mean(Weibull_rmse_values_b_lambda_mean_inv)
Weibull_tab.mean.inv.lambda[3,6]<-mean(Weibull_rmse_values_c_lambda_mean_inv)


########################################## Inverse Mean Lambda #############################################

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
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  data <- Weibull_data.list[[i]] 
  data$inv.trait <- 1/(data$trait)
  group_means <- tapply(data$inv.trait, data$T, function(x) {
    split_values <- split(x, rep(1:10, each = 10))
    means <- sapply(split_values, mean)
    return(means)
  })
  df.mean<- data.frame(
    T=c(rep(10,10),rep(17,10),rep(21,10),rep(25,10),rep(32,10)),  # Extract T values as numeric
    trait = unlist(group_means)  # Extract means and unlist the result
  )
  data.raw.w<-Weibull_data.list[[i]]
  data<-df.mean
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  #Model
  sink("Weibull_inv_mean_lambda.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  sig2~dexp(150)
  tau<-1/sig2
  #Likelihood
  for (i in 1:N.obs){
    mu[i] <- exp(exp(la) * temp[i]^2 - b * temp[i] + c)
    trait[i] ~ dnorm(mu[i],tau) T(0,)
  }
}",file="Weibull_inv_mean_lambda.txt")
  # Settings
  parameters <- c("la", "b", "c","sig2")
  inits<-function(){list(
    la = log(0.005),
    b = 0.4,
    c = 1.5,
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
                  model.file="Weibull_inv_mean_lambda.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_param_list5[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_param_list5[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_param_list5[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_param_list5[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_param_list5[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_param_list5[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_param_list5[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_param_list5[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_param_list5[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_param_list5[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_param_list5[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_param_list5[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_param_list5[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_param_list5[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_param_list5[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_param_list5[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # tau
    Weibull_param_list5[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_param_list5[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_param_list5[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_param_list5[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  
  df.new5_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf5_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new5_W))
  
  # Apply the quadratic equation for each combination
  xdf5_W$value <- apply(xdf5_W, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    exp(-(exp(df.new5_W[i, 4]) * location^2 - df.new5_W[i, 1] * location + df.new5_W[i, 2]))
  })
  xdf5_W<- xdf5_W|> mutate(trunc.eval=ifelse(xdf5_W$value>0,xdf5_W$value,epsilon),
                           trunc.inv=1/trunc.eval)
  mean.xdf5_W<-xdf5_W |> group_by(point) |> summarize(avg.value.inv=mean(trunc.inv),
                                                      avg.value=mean(trunc.eval),
                                                      med.value.inv=median(trunc.inv),
                                                      med.value=median(trunc.eval),
                                                      lower.hdi.inv=hdi(trunc.inv)[1],
                                                      upper.hdi.inv=hdi(trunc.inv)[2],
                                                      lower.hdi=hdi(trunc.eval)[1],
                                                      upper.hdi=hdi(trunc.eval)[2])
  true_curve.inv <- function(x) 1/exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  true_curve <- function(x) exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  mean.xdf5_W$true_curve <- true_curve(xseq)
  mean.xdf5_W$true_curve.inv <- true_curve.inv(xseq)
  plot5_mean_W<-ggplot(mean.xdf5_W, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot5_med_W<-ggplot(mean.xdf5_W, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot5_mean.inv_W<-ggplot(mean.xdf5_W, aes(x = point)) +
    geom_line(aes(y = avg.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    #geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
  plot5_med.inv_W<-ggplot(mean.xdf5_W, aes(x = point)) +
    geom_line(aes(y = med.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    #geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
}

gridExtra::grid.arrange(plot5_mean_W,plot5_med_W,nrow=1)
gridExtra::grid.arrange(plot5_mean.inv_W,plot5_med.inv_W,nrow=1)

Weibull_hist_list5 <- vector("list", length = 5)
par(mfrow=c(2,3))
# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_first_column_list5 <- lapply(Weibull_param_list5[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_combined_vector5 <- unlist(Weibull_first_column_list5)
  x<-seq(min(Weibull_combined_vector5)-1,max(Weibull_combined_vector5)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_hist_list5[[i]] <- hist(Weibull_combined_vector5, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                          xlab = "Values", col = "mintcream", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "darkmagenta", lwd = 2)
  if(i==1){
    lines(x,dnorm(x,0, sqrt(10)), col="green", lty=2, lwd=2)
  }
  if(i==2){
    lines(x, dexp(x, 0.5), col="black", lty=2, lwd=2)
  }
  if(i==3){
    lines(x, dexp(x, 0.5), col="yellow", lty=2, lwd=2)
  }
  if(i==5){
    lines(x, dexp(x, 12), col="purple", lty=2, lwd=2)
  }
}

Weibull_proportions_list_a5<- lapply(Weibull_param_list5[[1]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_a5))

Weibull_proportions_list_b5<- lapply(Weibull_param_list5[[2]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_b5))

Weibull_proportions_list_c5<- lapply(Weibull_param_list5[[3]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_c5))

Weibull_calculate_rmse_a_lambda_inv_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_b_lambda_inv_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_c_lambda_inv_mean <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_rmse_values_a_lambda_inv_mean <- sapply(Weibull_param_list5[[1]], Weibull_calculate_rmse_a_lambda_inv_mean)
Weibull_rmse_values_b_lambda_inv_mean <- sapply(Weibull_param_list5[[2]], Weibull_calculate_rmse_b_lambda_inv_mean)
Weibull_rmse_values_c_lambda_inv_mean <- sapply(Weibull_param_list5[[3]], Weibull_calculate_rmse_c_lambda_inv_mean)

Weibull_tab.inv.mean.lambda<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_tab.inv.mean.lambda)<-c("a","b","c")
colnames(Weibull_tab.inv.mean.lambda)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_tab.inv.mean.lambda[1,1]<-true.exp.a
Weibull_tab.inv.mean.lambda[2,1]<-true.exp.b
Weibull_tab.inv.mean.lambda[3,1]<-true.exp.c
Weibull_tab.inv.mean.lambda[1,2]<-exp(mean(unlist(lapply(Weibull_param_list5[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.inv.mean.lambda[2,2]<-(mean(unlist(lapply(Weibull_param_list5[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.inv.mean.lambda[3,2]<-mean(unlist(lapply(Weibull_param_list5[[3]], function(matrix) matrix[, 1])))
Weibull_tab.inv.mean.lambda[1,3]<-exp(hdi(unlist(lapply(Weibull_param_list5[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.inv.mean.lambda[2,3]<-(hdi(unlist(lapply(Weibull_param_list5[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.inv.mean.lambda[3,3]<-hdi(unlist(lapply(Weibull_param_list5[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.inv.mean.lambda[1,4]<-exp(hdi(unlist(lapply(Weibull_param_list5[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.inv.mean.lambda[2,4]<-(hdi(unlist(lapply(Weibull_param_list5[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.inv.mean.lambda[3,4]<-hdi(unlist(lapply(Weibull_param_list5[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.inv.mean.lambda[1,5]<-mean(unlist(Weibull_proportions_list_a5))
Weibull_tab.inv.mean.lambda[2,5]<-mean(unlist(Weibull_proportions_list_b5))
Weibull_tab.inv.mean.lambda[3,5]<-mean(unlist(Weibull_proportions_list_c5))
Weibull_tab.inv.mean.lambda[1,6]<-mean(Weibull_rmse_values_a_lambda_inv_mean)
Weibull_tab.inv.mean.lambda[2,6]<-mean(Weibull_rmse_values_b_lambda_inv_mean)
Weibull_tab.inv.mean.lambda[3,6]<-mean(Weibull_rmse_values_c_lambda_inv_mean)






########################################## Weibull ###########################################

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
for(i in 1:w){
  if (i %% 10 == 0) {
    print(i)
  }
  # Model
  sink("Weibull_weibull.txt")
  cat("model{
  #Priors
  #If mu<0, make it small number
  #And nonzero
  la~dnorm(0, 1/10) #has to be positive exp(1)?
  b~dexp(0.5) #Has to be positive, because function has neg sign
  c~dexp(0.5) #Has to be positive
  shape~dexp(0.001)
  #Likelihood
  for (i in 1:N.obs){
    trait[i] ~ dgen.gamma(1, 1/lambda[i], shape) #
    lambda[i] <- (med[i])/(log(2))^(1/shape)
    med[i] <- exp(-(exp(la) * temp[i]^2 - b * temp[i] + c))
  }
  }",file="Weibull_weibull.txt")
  
  # Settings
  parameters <- c("la", "b", "c","shape")
  inits<-function(){list(
    la = log(0.005),
    b = 0.4,
    c = 1.5,
    shape=sh.est
  )}
  ni <- 60000 
  nb <- 30000 
  nt <- 8 
  nc <- 5 
  data.raw.w<-Weibull_data.list[[i]]
  data <- Weibull_data.list[[i]] 
  trait <- data$trait
  N.obs <- length(trait)
  temp <- data$T
  
  # List data
  jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
  
  # Run model  
  mod.fit <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                  model.file="Weibull_weibull.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                  n.iter=ni, DIC=T, working.directory=getwd())
  mod.fit.mcmc <- as.mcmc(mod.fit) 
  closeAllConnections()
  samp<-as.mcmc.list(mod.fit.mcmc)
  for(j in 1:nc){
    #la
    Weibull_param_list6[[1]][[i]][j,1]<-mean(samp[[j]][,4])
    Weibull_param_list6[[1]][[i]][j,2]<-median(samp[[j]][,4])
    Weibull_param_list6[[1]][[i]][j,3]<-hdi(samp[[j]][,4])[1]
    Weibull_param_list6[[1]][[i]][j,4]<-hdi(samp[[j]][,4])[2]
    # b.l
    Weibull_param_list6[[2]][[i]][j,1]<-mean(samp[[j]][,1])
    Weibull_param_list6[[2]][[i]][j,2]<-median(samp[[j]][,1])
    Weibull_param_list6[[2]][[i]][j,3]<-hdi(samp[[j]][,1])[1]
    Weibull_param_list6[[2]][[i]][j,4]<-hdi(samp[[j]][,1])[2]
    # c
    Weibull_param_list6[[3]][[i]][j,1]<-mean(samp[[j]][,2])
    Weibull_param_list6[[3]][[i]][j,2]<-median(samp[[j]][,2])
    Weibull_param_list6[[3]][[i]][j,3]<-hdi(samp[[j]][,2])[1]
    Weibull_param_list6[[3]][[i]][j,4]<-hdi(samp[[j]][,2])[2]
    # dev
    Weibull_param_list6[[4]][[i]][j,1]<-mean(samp[[j]][,3])
    Weibull_param_list6[[4]][[i]][j,2]<-median(samp[[j]][,3])
    Weibull_param_list6[[4]][[i]][j,3]<-hdi(samp[[j]][,3])[1]
    Weibull_param_list6[[4]][[i]][j,4]<-hdi(samp[[j]][,3])[2]
    # shape
    Weibull_param_list6[[5]][[i]][j,1]<-mean(samp[[j]][,5])
    Weibull_param_list6[[5]][[i]][j,2]<-median(samp[[j]][,5])
    Weibull_param_list6[[5]][[i]][j,3]<-hdi(samp[[j]][,5])[1]
    Weibull_param_list6[[5]][[i]][j,4]<-hdi(samp[[j]][,5])[2]
  }
  
  df.new6_W = samp[[1]][seq(1, nrow(samp[[1]]), 4), ]
  xseq<-seq(from=10,to=35,length.out=250)
  xdf6_W <- expand.grid(point = xseq, iteration = 1:nrow(df.new6_W))
  
  xdf6_W$value <- apply(xdf6_W, 1, function(row) {
    i <- row["iteration"]
    location <- row["point"]
    exp(-(exp(df.new6_W[i, 4]) * location^2 - df.new6_W[i, 1] * location + df.new6_W[i, 2]))
  })
  xdf6_W<- xdf6_W|> mutate(trunc.eval=ifelse(xdf6_W$value>0,xdf6_W$value,epsilon),
                           trunc.inv=log(2)/trunc.eval)
  mean.xdf6_W<-xdf6_W |> group_by(point) |> summarize(avg.value.inv=mean(trunc.inv),
                                                      avg.value=mean(trunc.eval),
                                                      med.value.inv=median(trunc.inv),
                                                      med.value=median(trunc.eval),
                                                      lower.hdi.inv=hdi(trunc.inv)[1],
                                                      upper.hdi.inv=hdi(trunc.inv)[2],
                                                      lower.hdi=hdi(trunc.eval)[1],
                                                      upper.hdi=hdi(trunc.eval)[2])
  true_curve.inv <- function(x) log(2)/exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  true_curve <- function(x) exp(-(true.exp.a * x^2 - true.exp.b * x + true.exp.c))
  mean.xdf6_W$true_curve <- true_curve(xseq)
  mean.xdf6_W$true_curve.inv <- true_curve.inv(xseq)
  plot6_mean_W<-ggplot(mean.xdf6_W, aes(x = point)) +
    geom_line(aes(y = avg.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot6_med_W<-ggplot(mean.xdf6_W, aes(x = point)) +
    geom_line(aes(y = med.value), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve),color="hotpink",linetype="dashed")+
    geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Lifetime")+
    theme_minimal()
  plot6_mean.inv_W<-ggplot(mean.xdf6_W, aes(x = point)) +
    geom_line(aes(y = avg.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    #geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Mean Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
  plot6_med.inv_W<-ggplot(mean.xdf6_W, aes(x = point)) +
    geom_line(aes(y = med.value.inv), color = "black") +
    geom_ribbon(aes(ymin = lower.hdi.inv, ymax = upper.hdi.inv), fill = "darkorange1", alpha = 0.3)+
    geom_line(aes(y=true_curve.inv),color="hotpink",linetype="dashed")+
    #geom_point(aes(x=T,y=trait),data=data.raw.w)+
    labs(title = "Median Curve with Interval Bands and True Curve",
         x = "Temperature",
         y = "Mortality Rate")+
    theme_minimal()
}

gridExtra::grid.arrange(plot6_mean_W,plot6_med_W,nrow=1)
gridExtra::grid.arrange(plot6_mean.inv_W,plot6_med.inv_W,nrow=1)


Weibull_hist_list6 <- vector("list", length = 5)
par(mfrow=c(2,3))

# Loop through indices 1 to 4
for (i in 1:5) {
  # Extract the first column from each matrix in param_list[[i]]
  Weibull_first_column_list6 <- lapply(Weibull_param_list6[[i]], function(matrix) matrix[, 1])
  
  # Combine the vectors into a single vector
  Weibull_combined_vector6 <- unlist(Weibull_first_column_list6)
  x<-seq(min(Weibull_combined_vector6)-1,max(Weibull_combined_vector6)+1, length=1000)
  # Create a histogram and store it in the list
  Weibull_hist_list6[[i]] <- hist(Weibull_combined_vector6, main = paste("Combined Histogram Mean (Parameter", i, ")"), 
                                    xlab = "Values", col = "mintcream", border = "black", breaks = 10,freq=FALSE)
  abline(v = Weibull_true_values[i], col = "darkmagenta", lwd = 2)
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


Weibull_proportions_list_a6<- lapply(Weibull_param_list6[[1]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_a6))

Weibull_proportions_list_b6<- lapply(Weibull_param_list6[[2]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_b6))

Weibull_proportions_list_c6<- lapply(Weibull_param_list6[[3]], function(matrix) {
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
mean(unlist(Weibull_proportions_list_c6))

Weibull_calculate_rmse_a_wb <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(log(true.exp.a), length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_b_wb <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.b, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}
Weibull_calculate_rmse_c_wb <- function(mat) {
  observed <- mat[, 1]  # Assuming the first column is the observed values
  predicted <- rep(true.exp.c, length(observed))  # Replace this with your predicted values
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  return(rmse)
}


# Apply the function to each element in the list
Weibull_rmse_values_a_wb <- sapply(Weibull_param_list6[[1]], Weibull_calculate_rmse_a_wb)
Weibull_rmse_values_b_wb <- sapply(Weibull_param_list6[[2]], Weibull_calculate_rmse_b_wb)
Weibull_rmse_values_c_wb <- sapply(Weibull_param_list6[[3]], Weibull_calculate_rmse_c_wb)

Weibull_tab.wb<-matrix(0,nrow=3,ncol=6)
row.names(Weibull_tab.wb)<-c("a","b","c")
colnames(Weibull_tab.wb)<-c("True Value","Mean","Lower HDI of Mean","Upper HDI of Mean","Coverage","RMSE")
Weibull_tab.wb[1,1]<-true.exp.a
Weibull_tab.wb[2,1]<-true.exp.b
Weibull_tab.wb[3,1]<-true.exp.c
Weibull_tab.wb[1,2]<-exp(mean(unlist(lapply(Weibull_param_list6[[1]], function(matrix) matrix[, 1]))))
Weibull_tab.wb[2,2]<-(mean(unlist(lapply(Weibull_param_list6[[2]], function(matrix) matrix[, 1]))))
Weibull_tab.wb[3,2]<-mean(unlist(lapply(Weibull_param_list6[[3]], function(matrix) matrix[, 1])))
Weibull_tab.wb[1,3]<-exp(hdi(unlist(lapply(Weibull_param_list6[[1]], function(matrix) matrix[, 1])))[1])
Weibull_tab.wb[2,3]<-(hdi(unlist(lapply(Weibull_param_list6[[2]], function(matrix) matrix[, 1])))[1])
Weibull_tab.wb[3,3]<-hdi(unlist(lapply(Weibull_param_list6[[3]], function(matrix) matrix[, 1])))[1]
Weibull_tab.wb[1,4]<-exp(hdi(unlist(lapply(Weibull_param_list6[[1]], function(matrix) matrix[, 1])))[2])
Weibull_tab.wb[2,4]<-(hdi(unlist(lapply(Weibull_param_list6[[2]], function(matrix) matrix[, 1])))[2])
Weibull_tab.wb[3,4]<-hdi(unlist(lapply(Weibull_param_list6[[3]], function(matrix) matrix[, 1])))[2]
Weibull_tab.wb[1,5]<-mean(unlist(Weibull_proportions_list_a6))
Weibull_tab.wb[2,5]<-mean(unlist(Weibull_proportions_list_b6))
Weibull_tab.wb[3,5]<-mean(unlist(Weibull_proportions_list_c6))
Weibull_tab.wb[1,6]<-mean(Weibull_rmse_values_a_wb)
Weibull_tab.wb[2,6]<-mean(Weibull_rmse_values_b_wb)
Weibull_tab.wb[3,6]<-mean(Weibull_rmse_values_c_wb)



