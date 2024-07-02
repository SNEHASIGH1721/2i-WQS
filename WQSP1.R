#WQS
#WQS regression model has been widely used in the 
#field of environmental epidemiology to address the impact of multiple environmental chemical exposure on certain outcome. WQS is established to address highly correlated and collinear data. We 
#utilized WQS regression to estimate joint associations of all five OPE 
#metabolites with CVD and calculated the empirical weight of each 
#element. Each weighted index represented the contribution of each OPE 
#metabolite substance included in the mixture to the positive association. 
#To begin with, concentrations of all OPE metabolites were ranked in 
#quartiles. Next, all data were randomly split into training set and validation set). The weighted index is constrained to a 
#range of 0–1, and to add up to 1 (Yu et al., 2022a). Based on previous 
#published studies, 1000 bootstrap replicates were performed to estimate 
#the magnitude of the effect and corresponding 95% confidence intervals. 
#Data were randomly split into 60% validation set and 40% test set.
#researchers might be interested in taking into account the relationship between the exposures and the outcome while summarizing the complex exposure to the mixture of interest. The weighted quantile sum (WQS), developed specifically for the context of environmental mixtures analysis, is an increasingly common approach that allows evaluating a mixture-outcome association by 
#creating a summary score of the mixture in a supervised fashion
#https://bookdown.org/andreabellavia/mixtures/the-weighted-quantile-sum-wqs-and-its-extensions.html use this for additional information and help 
getwd()
load("Datafiltered31NEHA.RData")
ls()
# converting the metabolites into the quarters
#1 DPHP_cr_In
filtered_data31 <- filtered_data31 %>%
  mutate(DPHP_quartiled = cut(DPHP_cr_ln, 
                              breaks = quantile(DPHP_cr_ln, probs = 0:4/4, na.rm = TRUE),
                              include.lowest = TRUE, 
                              labels = c("Q1", "Q2", "Q3", "Q4")))
#2 BDCPP_cr_ln
filtered_data31 <- filtered_data31 %>%
  mutate(BDCPP_quartiled = cut(BDCPP_cr_ln, 
                               breaks = quantile(BDCPP_cr_ln, probs = 0:4/4, na.rm = TRUE),
                               include.lowest = TRUE, 
                               labels = c("Q1", "Q2", "Q3", "Q4")))
#3 BCPP_cr_ln
filtered_data31 <- filtered_data31 %>%
  mutate(BCPP_quartiled = cut(BCPP_cr_ln, 
                              breaks = quantile(BCPP_cr_ln, probs = 0:4/4, na.rm = TRUE),
                              include.lowest = TRUE, 
                              labels = c("Q1", "Q2", "Q3", "Q4")))
#4 BCEP_cr_ln
filtered_data31 <- filtered_data31 %>%
  mutate(BCEP_quartiled = cut(BCEP_cr_ln, 
                              breaks = quantile(BCEP_cr_ln, probs = 0:4/4, na.rm = TRUE),
                              include.lowest = TRUE, 
                              labels = c("Q1", "Q2", "Q3", "Q4")))
#5 DBUP_cr_ln
filtered_data31 <- filtered_data31 %>%
  mutate(DBUP_quartiled = cut(DBUP_cr_ln, 
                              breaks = quantile(DBUP_cr_ln, probs = 0:4/4, na.rm = TRUE),
                              include.lowest = TRUE, 
                              labels = c("Q1", "Q2", "Q3", "Q4")))
names(filtered_data31)
str(filtered_data31)
library(bkmr)
#install.packages("gWQS")#This is the package for WQS
library(gWQS)
# Convert factor quartiles to integers
filtered_data31$DPHP_quartiled <- as.integer(filtered_data31$DPHP_quartiled)
filtered_data31$BDCPP_quartiled <- as.integer(filtered_data31$BDCPP_quartiled)
filtered_data31$BCPP_quartiled <- as.integer(filtered_data31$BCPP_quartiled)
filtered_data31$BCEP_quartiled <- as.integer(filtered_data31$BCEP_quartiled)
filtered_data31$DBUP_quartiled <- as.integer(filtered_data31$DBUP_quartiled)


#install.packages("coda")# For MCMX diagnostics
filtered_data31$DPHP_cr_ln<-NULL
filtered_data31$BDCPP_cr_ln<-NULL
filtered_data31$BCPP_cr_ln<-NULL
filtered_data31$BCEP_cr_ln<-NULL
filtered_data31$DBUP_cr_ln<-NULL
filtered_data31$ID<-NULL
filtered_data31$DPHP<-NULL
filtered_data31$BDCPP<-NULL
filtered_data31$BCPP<-NULL
filtered_data31$BCEP<-NULL
filtered_data31$DBUP<-NULL
filtered_data31$W.ALL<-NULL
filtered_data31$BCPP<-NULL
filtered_data31$BCEP<-NULL
filtered_data31$DBUP<-NULL
filtered_data31$W.MEC<-NULL
filtered_data31$PSU<-NULL
filtered_data31$STRATA<-NULL
filtered_data31$RA_status_binary <- ifelse(filtered_data31$RA_status == "RA", 1, 0)
filtered_data31$RA_status<-NULL
filtered_data31$Age<-NULL
filtered_data31$DPHP_cr<-NULL
filtered_data31$BDCPP_cr<-NULL
filtered_data31$BCPP_cr<-NULL
filtered_data31$BCEP_cr<-NULL
filtered_data31$DBUP_cr<-NULL
filtered_data31$logDPHP<-NULL
filtered_data31$logBDCPP<-NULL
filtered_data31$logBCPP<-NULL
filtered_data31$logBCEP<-NULL
filtered_data31$logDBUP<-NULL


names(filtered_data31)
head(filtered_data31)
library(dplyr)


exposure<-names(filtered_data31[, 31:35])
# Z-standardize the exposure variables

##
# Assuming exposure is correctly defined positive association,
results1 <- gwqs(RA_status_binary ~ wqs, mix_name = exposure, data = filtered_data31, q = 4, 
                 validation = 0.6, b = 100, b1_pos = TRUE, b_constr = FALSE, 
                 family = "binomial", seed = 123)
gwqs_barplot(results1, tau=NULL)
# fitted values vs rediduals scatterplot
gwqs_fitted_vs_resid(results1)
gwqs_weights_tab(results1)

####
# Create a data frame for the weights
weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.4190, 0.2500, 0.1590, 0.1060, 0.0661)
)

# Calculate the total weight
total_weight <- sum(weights_df$mean_weight)

# Calculate the percentage for each weight
weights_df$percentage <- (weights_df$mean_weight / total_weight) * 100

# Calculate the mean weight for the red line
mean_weight <- mean(weights_df$mean_weight)

# Load the ggplot2 library
library(ggplot2)
# Create the bar plot (this one will fill the bar with gray color and it for positive association)
p <- ggplot(weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  geom_text(aes(label = paste0(sprintf("%.2f", percentage), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.3) +
  geom_vline(xintercept = mean_weight, color = "red", linetype = "dashed") +
  coord_flip() + # Flip axes to make it horizontal
  labs(y = "Weights", x = "Exposure", title = "Weighted Quantile Sum Regression Weights") +
  theme_minimal() # Use a minimal theme for a clean look

# Print the plot
print(p)


###To adjust the positive WQS regression model for confounders we can add them in the model as presented here
# Assuming exposure is correctly defined this code snippet is not negative association
results1_adjusted <- gwqs(RA_status_binary ~ wqs, mix_name = exposure, data = filtered_data31, q = 4, 
                 validation = 0.6, b = 10, b1_pos = FALSE, b_constr = FALSE, 
                 family = "gaussian", seed = 123)
gwqs_barplot(results1_adjusted, tau=NULL)
gwqs_weights_tab(results1_adjusted)


# Create a data frame for the weights
weights_df <- data.frame(
  mix_name = c("DPHP_quartiled", "DBUP_quartiled", "BCEP_quartiled", "BCPP_quartiled", "BDCPP_quartiled"),
  mean_weight = c(0.4190, 0.2500, 0.1590, 0.1060, 0.0661)
)

# Calculate the total weight
total_weight <- sum(weights_df$mean_weight)

# Calculate the percentage for each weight
weights_df$percentage <- (weights_df$mean_weight / total_weight) * 100

# Calculate the mean weight for the red line
mean_weight <- mean(weights_df$mean_weight)

# Load the ggplot2 library
library(ggplot2)

# Create the bar plot with red line and percentage labels
p <- ggplot(weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f%%", percentage), y = mean_weight), vjust = 1.5, color = "black") +
  geom_vline(xintercept = mean_weight, color = "red", linetype = "dashed", size = 3) +
  #geom_text(aes(label = mix_name, y = -0.02), angle = 90, hjust = 1, vjust = 0.5, color = "black")
  #coord_flip() + # Flip axes to make it horizontal
  labs(x = "Exposure", y = "WQS index weights", title = "Weighted Quantile Sum Regression Weights") +
  theme_minimal() +# Use a minimal theme for a clean look
  theme(legend.position = "right")
  # Position the legend at the bottom
  #theme(legend.position = "none") # Remove legend



# Print the plot
print(p)
#the code end here
library(ggplot2)

# ... [your existing code to create the weights_df data frame] ...

p <- ggplot(weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f%%", percentage), y = mean_weight), vjust = 1.5, color = "black") +
  geom_vline(xintercept = mean_weight, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Exposure", y = "WQS index weights", title = "Weighted Quantile Sum Regression Weights") +
  theme_minimal() +
  theme(legend.position = "right", legend.background = element_rect(fill = "grey90"), legend.key = element_blank()) +
  scale_fill_manual(values = c("BCEP_quartiled" = "green", "BCPP_quartiled" = "blue", 
                               "BDCPP_quartiled" = "yellow", "DBUP_quartiled" = "orange", 
                               "DPHP_quartiled" = "purple"))

# Print the plot
print(p)

####
#Positive Index
###To adjust the positive WQS regression model for confounders we can add them in the model as presented here
# Assuming exposure is correctly defined this code snippet is not negative association
results1_adjusted <- gwqs(RA_status_binary ~ wqs, mix_name = exposure, data = filtered_data31, q = 4, 
                          validation = 0.6, b = 10, b1_pos = TRUE, b_constr = FALSE, 
                          family = "gaussian", seed = 123)
gwqs_barplot(results1_adjusted, tau=NULL)
gwqs_weights_tab(results1_adjusted)


# Create a data frame for the weights
weights_df <- data.frame(
  mix_name = c("DPHP_quartiled", "DBUP_quartiled", "BCEP_quartiled", "BCPP_quartiled", "BDCPP_quartiled"),
  mean_weight = c(0.284, 0.236, 0.188, 0.171, 0.121)
)

# Calculate the total weight
total_weight <- sum(weights_df$mean_weight)

# Calculate the percentage for each weight
weights_df$percentage <- (weights_df$mean_weight / total_weight) * 100

# Calculate the mean weight for the red line
mean_weight <- mean(weights_df$mean_weight)

# Load the ggplot2 library
library(ggplot2)

# Create the bar plot with red line and percentage labels
p <- ggplot(weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f%%", percentage), y = mean_weight), vjust = 1.5, color = "black") +
  geom_vline(xintercept = mean_weight, color = "red", linetype = "dashed", size = 3) +
  #coord_flip() + # Flip axes to make it horizontal
  labs(x = "Weights", y = "Exposure", title = "Weighted Quantile Sum Regression Weights-Ng") +
  theme_minimal() + # Use a minimal theme for a clean look
  theme(legend.position = "none") # Remove legend

# Print the plot
print(p)




#This is the  OR ratio plot in box
#the code end here
library(ggplot2)

# Create a data frame for the Rheumatoid arthritis indices and their Odds Ratios (ORs)
RA_data <- data.frame(
  Index = c("Rheumatoid arthritis positive index", "Rheumatoid arthritis negative index"),
  OR = c(2, 0.5),
  lower_CI = c(1, 0.4),
  upper_CI = c(4, 0.7)
)

# Calculate the ORs and their lower and upper error values for error bars
RA_data$lower_error <- RA_data$OR - RA_data$lower_CI
RA_data$upper_error <- RA_data$upper_CI - RA_data$OR

# Create the plot
p <- ggplot(RA_data, aes(y = Index, x = OR)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower_CI, xmax = upper_CI), width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(x = "Odds Ratio (OR)", title = "Association of Rheumatoid arthritis  with OPEs Exposure Indices") +
  #coord_flip() + # Flip axes to display horizontal plot
  theme_minimal()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) # Add a box around the plot

# Print the plot
print(p)










