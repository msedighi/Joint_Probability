library(d3heatmap)
library(geigen)

min_d <- 5;
max_d <- 5;

Num_d <- max_d - min_d + 1

Num_Samples <- 2000

# Total_FreeEnergy_0 <- matrix(NA,Num_d,Num_d)
# Total_FreeEnergy <- matrix(NA,Num_d,Num_d)
# TotalIndependent_FreeEnergy <- matrix(NA,Num_d,Num_d)
# Partition_Function_0 <- matrix(NA,Num_d,Num_d)
# Partition_Function <- matrix(NA,Num_d,Num_d)
# 
# Wave_Energy <- array(NA,dim=c(Num_Samples,Num_d,Num_d))
# Max_Wave_Energy <- array(NA,dim=c(Num_Samples,Num_d,Num_d))
# Min_Wave_Energy <- array(NA,dim=c(Num_Samples,Num_d,Num_d))
# Median_Wave_Energy <- array(NA,dim=c(Num_Samples,Num_d,Num_d))

Laplacian_Energy <- array(NA,dim=c(Num_Samples,Num_d))

# Hamiltonian_Spectrum <- array(0,dim=c(Num_Samples,max_d,Num_d,Num_d))
Laplacian_Spectrum <- array(0,dim=c(Num_Samples,max_d,Num_d))

for (d1 in seq(min_d, max_d)) 
{
  V_1 <- matrix(1,d1,1)
  i1 <- d1 - min_d + 1; # Index 
  
  # for (d2 in seq(min_d, max_d)) 
  for (d2 in d1) 
  {
    
    V_2 <- matrix(1,d2,1)
    i2 <- d2 - min_d + 1; # Index 
    
    for (cnt in seq(1,Num_Samples))
    {
      
      N_12 <- matrix(runif(d1*d2,0,100),d1,d2)
      N_1_vector <- N_12 %*% V_2 # Sum the Rows
      N_2_vector <- t(V_1) %*% N_12  # Sum the Columns
      N_1 <-diag(as.vector(N_1_vector)) # Diagonal Metric (Number Operator)  for parameter 1
      N_2 <-diag(as.vector(N_2_vector)) # Diagonal Metric (Number Operator) for parameter 2
      
      N <- as.numeric(t(V_1) %*% N_12 %*% V_2)
      N_12_Independent <- (N_1_vector %*% N_2_vector)/N
      
      # Validation!
      #     N_1_vector - N_12_Independent %*% V_2
      #     N_1_vector - N_12 %*% V_2
      #     N_2_vector - t(V_1) %*% N_12_Independent
      #     N_2_vector - t(V_1) %*% N_12
      
      Pure_Correlation <- N_12/N_12_Independent
      FreeEnergy <- -log(Pure_Correlation)
      
      #       Total_FreeEnergy_0[i1,i2] <- t(V_1) %*% FreeEnergy %*% V_2
      #       Total_FreeEnergy[i1,i2] <- t(V_1) %*% (N_12*FreeEnergy) %*% V_2
      #       TotalIndependent_FreeEnergy[i1,i2] <- t(V_1) %*% (N_12_Independent*FreeEnergy) %*% V_2
      #   
      #       Partition_Function_0[i1,i2] <- t(V_1) %*% Pure_Correlation %*% V_2
      #       Partition_Function[i1,i2] <- t(N_1_vector) %*% Pure_Correlation %*% t(N_2_vector) # This is equal to N^2
      
      Laplacian <- Pure_Correlation - matrix(1,d1,d2)
      Laplacian_SVD <- svd(Laplacian)
      
      # Validation!
      # Laplacian - Laplacian_SVD$u %*% diag(Laplacian_SVD$d) %*% t(Laplacian_SVD$v)
      # Laplacian_SVD$d - diag(t(Laplacian_SVD$u) %*% Laplacian %*% Laplacian_SVD$v)
      #**Zero Eigenvector**   N_2_vector - Laplacian_SVD$v[,5]/sum(Laplacian_SVD$v[,5])*N
      #**Zero Eigenvector**   N_1_vector - Laplacian_SVD$u[,5]/sum(Laplacian_SVD$u[,5])*N
      
      # This is Important!!
      # Laplacian_SVD$u %*% diag(Laplacian_SVD$d) = Laplacian %*% Laplacian_SVD$v
      
      Laplacian_Spectrum[cnt,1:d2,i1] <- Laplacian_SVD$d; 
      Laplacian_Energy[cnt,i1] <- sum(Laplacian_SVD$d);
      
      # Hamiltonian <- t(Pure_Correlation) %*% N_1 %*% Pure_Correlation - matrix(N, d2, d2)
      # Wave_Energy[cnt,i1,i2] <- sum(diag(Hamiltonian %*% N_2))/N^2 # This is the sum of all Hamiltonian energies i.e. Sum of all Hamiltonian eigenvalues
      #       Hamiltonian_Spectrum[cnt,d2:1,i1,i2] <- geigen(Hamiltonian,solve(N_2),symmetric = TRUE, only.values = TRUE)$values/N^2 # Generalized eigenvalue problem i.e. Hamiltonian*psi_2 = lambda*N_2*psi_2
      #       Max_Wave_Energy[cnt,i1,i2] <- max(Hamiltonian_Spectrum[cnt,d2:1,i1,i2])
      #       Min_Wave_Energy[cnt,i1,i2] <- min(Hamiltonian_Spectrum[cnt,d2:1,i1,i2])
      #       Median_Wave_Energy[cnt,i1,i2] <- median(Hamiltonian_Spectrum[cnt,d2:1,i1,i2])
    }
    
  }
}

# Mean_Laplacian_Energy <- diag(apply(Laplacian_Energy,c(2,3),mean);
# Mean_eigz_Laplacian_Energy <- apply(apply(Laplacian_Spectrum, c(2,3,4), mean),1,diag);
Mean_Laplacian_Energy <- colMeans(Laplacian_Energy);
Mean_eigz_Laplacian_Energy <- apply(Laplacian_Spectrum, c(2,3), mean);
                              
x11()
matplot(Mean_Laplacian_Energy, type = "l", xlab = NULL, ylab = NULL)
x11()
matplot(t(Mean_eigz_Laplacian_Energy), type = "l", xlab = NULL, ylab = NULL)

# Total_FreeEnergy_0
# Total_FreeEnergy
# TotalIndependent_FreeEnergy
# 
# Partition_Function

# There are the following Questions to investigate:
#   1. What happens if the number of "system parameters" (features of the system) change? Namely, solving for N_123 (N_12..k) i.e. adding a N_3, V_3, etc
#     This helps to theoretically understand how every parameter(feature) contributes to the total Free-Energy.
#
#   2. A lot more sampling of random matrices should be collected for a significant result.
#
#   3. What happens when "dimension(resolution) of each feature" (d1, d2, etc) becomes larger? (Partially Done!) 

stop()

High_Quantile_eigz_Wave_Energy <- apply(apply(Hamiltonian_Spectrum, c(2,3,4), quantile, probs=c(.95)),1,diag)
Low_Quantile_eigz_Wave_Energy <- apply(apply(Hamiltonian_Spectrum, c(2,3,4), quantile, probs=c(.05)),1,diag)

Mean_eigz_Wave_Energy <- apply(apply(Hamiltonian_Spectrum, c(2,3,4), mean),1,diag)
Std_eigz_Wave_Energy <- apply(apply(Hamiltonian_Spectrum, c(2,3,4), sd),1,diag)

Mean_Sum_Wave_Energy <- apply(Wave_Energy,c(2,3),mean)
High_Quantile_Sum_Wave_Energy <- apply(Wave_Energy,c(2,3),quantile, probs=c(c(.95)))
Low_Quantile_Sum_Wave_Energy <- apply(Wave_Energy,c(2,3),quantile, probs=c(c(.05)))

# Mean_Max_Wave_Energy <- apply(Max_Wave_Energy,c(2,3),mean)
# Std_Max_Wave_Energy <- apply(Max_Wave_Energy,c(2,3),sd)
# Mean_Min_Wave_Energy <- apply(Min_Wave_Energy,c(2,3),mean)
# Std_Min_Wave_Energy <- apply(Min_Wave_Energy,c(2,3),sd)
# Mean_Median_Wave_Energy <- apply(Median_Wave_Energy,c(2,3),mean)
# Std_Median_Wave_Energy <- apply(Median_Wave_Energy,c(2,3),sd)

# diag(Mean_Max_Wave_Energy)
# diag(Std_Max_Wave_Energy)

# Max_Wave_Energy_Range <- rbind(diag(Mean_Max_Wave_Energy+2*Std_Max_Wave_Energy),diag(Mean_Max_Wave_Energy),diag(Mean_Max_Wave_Energy-2*Std_Max_Wave_Energy))
# Min_Wave_Energy_Range <- rbind(diag(Mean_Min_Wave_Energy+2*Std_Min_Wave_Energy),diag(Mean_Min_Wave_Energy),diag(Mean_Min_Wave_Energy-2*Std_Min_Wave_Energy))
# Median_Wave_Energy_Range <- rbind(diag(Mean_Median_Wave_Energy+2*Std_Median_Wave_Energy),diag(Mean_Median_Wave_Energy),diag(Mean_Median_Wave_Energy-2*Std_Median_Wave_Energy))
Sum_Wave_Energy_Range <- rbind(diag(High_Quantile_Sum_Wave_Energy),diag(Mean_Sum_Wave_Energy),diag(Low_Quantile_Sum_Wave_Energy))
# Wave_Energy_Values <- t(rbind(Max_Wave_Energy_Range,Median_Wave_Energy_Range,Min_Wave_Energy_Range,Sum_Wave_Energy_Range))


x11()
matplot(t(Sum_Wave_Energy_Range), type = "l", xlab = NULL, ylab = NULL)
# matplot(Wave_Energy_Values, type = "l", col=c(3,3,3,4,4,4,2,2,2,1,1,1),xlab = NULL, ylab = NULL)
title(main = "Hamiltonian Energy Spectrum for Random Sampling of 2-D PDFs", xlab = "Dimension")

scale_factor <- matrix(rep(min_d:max_d,max_d),max_d-min_d+1,max_d)

x11()
# matplot(scale_factor*Mean_eigz_Wave_Energy, type = "l", col=rep(4,max_d), xlab = NULL, ylab = NULL, ylim=range(c(0,2)))
# par(new = TRUE)
matplot(t(scale_factor*High_Quantile_eigz_Wave_Energy), type = "l", col=rep(3,max_d), xlab = NULL, ylab = NULL, ylim=range(c(0,2)))
par(new = TRUE)
matplot(t(scale_factor*Low_Quantile_eigz_Wave_Energy), type = "l", col=rep(2,max_d), xlab = NULL, ylab = NULL, ylim=range(c(0,2)))
title(main = "Hamiltonian Energy Spectrum for Random Sampling of 2-D PDFs", xlab = "Dimension")

# x11()
# plot(diag(Std_Max_Wave_Energy))
# x11()
# plot(diag(Std_Min_Wave_Energy))

# x11()
# d3heatmap(Mean_Max_Wave_Energy, dendrogram = "none")
# x11()
# d3heatmap(Std_Max_Wave_Energy, dendrogram = "none")
