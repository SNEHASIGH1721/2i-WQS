#this code snippet is for the fully adjusted model 
getwd()
load("Datafiltered31NEHA.RData")
ls()
install.packages("pacman")
install.packages("devtools")

pacman::p_load(future, future.apply, dplyr, reshape2, MASS, qgcomp, gWQS)
# converting the metabolites into the quarters

names(filtered_data31)
str(filtered_data31)
library(bkmr)
library(dplyr)
#install.packages("gWQS")#This is the package for WQS
library(gWQS)
filtered_data31<- filtered_data31%>%
  ("logBDCPP","logBCPP","logBCEP","logDBUP","logDPHP","DPHP_cr","BDCPP_cr","BCPP_cr","BCEP_cr","DBUP_cr","ID","W.ALL","W.MEC","PSU","STRATA","DPHP","BDCPP","BCPP","BCEP","DBUP")

columns_to_delete <- c("logBDCPP", "logBCPP", "logBCEP", "logDBUP", "logDPHP", "DPHP_cr", "BDCPP_cr", "BCPP_cr", "BCEP_cr", "DBUP_cr", "ID", "W.ALL", "W.MEC", "PSU", "STRATA", "DPHP", "BDCPP", "BCPP", "BCEP", "DBUP")

filtered_data31 <- filtered_data31[, !(names(filtered_data31) %in% columns_to_delete)]

#install.packages("coda")# For MCMX diagnostics
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
set.seed(123)  # For reproducibility

# Define the proportion for training
train_proportion <- 0.7

# Create a random index for splitting
train_index <- sample(seq_len(nrow(filtered_data31)), size = train_proportion * nrow(filtered_data31))

# Split the data into training and validation sets
train_data <- filtered_data31[train_index, ]
validation_data <- filtered_data31[-train_index, ]



exposure<-names(filtered_data31[, 16:20])
# Z-standardize the exposure variables
# Fit the WQS model on training data with both positive and negative indices
results_train2 <- gwqs(
  RA_status_binary ~ pwqs + nwqs+Gender+Ethnicity+Marital_status+Education+Citizenship+Urinary_creatinine+Alcohol+BMI+Diabetes+Smoking+Age.cat+Activity_level+PIR, 
  mix_name = exposure, 
  data = train_data, 
  q = 4, 
  rh = 100, 
  validation = 0.6, 
  b = 100, 
  b_constr = FALSE, 
  family = binomial(), 
  seed = 123,
  lambda = 100,plan_strategy = "multisession", 
  solve_dir_issue = "inverse"  
)
# Summarize the training results
summary(results_train2)
results_train2$final_weights
##
# Manually calculate the WQS indices for the validation data
final_weights <- results_train2$final_weights
##
gwqs_barplot(results_train2, tau=NULL)
gwqs_weights_tab(results_train2)

# Create a matrix of exposure variables in the same order as the final weights
exposure_matrix <- as.matrix(validation_data[, exposure])

# Calculate positive and negative WQS indices for the validation data
validation_data$pwqs <- exposure_matrix %*% final_weights$`Estimate pos`
validation_data$nwqs <- exposure_matrix %*% final_weights$`Estimate neg`

# Create a logistic regression model to make predictions
validation_model <- glm(RA_status_binary ~ pwqs + nwqs, data = validation_data, family = binomial())

# Make predictions on the validation data
validation_data$predicted <- predict(validation_model, newdata = validation_data, type = "response")

# Ensure the factor levels are consistent
validation_predictions <- factor(validation_predictions, levels = c(0, 1))
validation_actual <- factor(validation_data$RA_status_binary, levels = c(0, 1))

conf_matrix <- confusionMatrix(validation_predictions, validation_actual)

# Print the confusion matrix
print(conf_matrix)
###
###
###
###
library(ggplot2)
library(tidyr)
library(dplyr)

# Define the proportion for training
train_proportion <- 0.7

# Create a random index for splitting
train_index <- sample(seq_len(nrow(filtered_data31)), size = train_proportion * nrow(filtered_data31))

# Split the data into training and validation sets
train_data <- filtered_data31[train_index, ]
validation_data <- filtered_data31[-train_index, ]

# Define exposure variables
exposure <- names(filtered_data31)[31:35]

# Fit the WQS model on training data with both positive and negative indices
results_train <- gwqs(
  RA_status_binary ~ pwqs + nwqs, 
  mix_name = exposure, 
  data = train_data, 
  q = 4, 
  rh = 3, 
  validation = 0.6, 
  b = 100, 
  b_constr = FALSE, 
  family = binomial(), 
  seed = 123,
  lambda = 100, 
  plan_strategy = "multisession", 
  solve_dir_issue = "inverse"
)

# Summarize the training results
summary(results_train)
final_weights <- results_train$final_weights
final_weights

# Define the number of mixture components
num_components <- 5

# Calculate the prespecified cutoff
cutoff <- 1 / num_components

# Identify significant weights for positive and negative indices
significant_weights_pos <- final_weights[final_weights$`Estimate pos` > cutoff, ]
significant_weights_neg <- final_weights[final_weights$`Estimate neg` > cutoff, ]

# Print the significant weights
print("Significant weights for positive index:")
print(significant_weights_pos)

print("Significant weights for negative index:")
print(significant_weights_neg)

# Create the dataframe for plotting
final_weights_long <- final_weights %>%
  gather(key = "Direction", value = "Weight", -mix_name)

# Plot the weights with the cutoff line
ggplot(final_weights_long, aes(x = mix_name, y = Weight, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "red") +
  facet_wrap(~ Direction) +
  labs(title = "Box-plot of Weights Associated with Positive and Negative Indices",
       x = "Nutrient", y = "Weight (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####

# Create the dataframe for plotting
final_weights_long <- final_weights %>%
  pivot_longer(cols = c(`Estimate pos`, `Estimate neg`), names_to = "Direction", values_to = "Weight") %>%
  mutate(Direction = ifelse(Direction == "Estimate pos", "Positive", "Negative"))

# Plot the weights with the cutoff line
ggplot(final_weights_long, aes(x = mix_name, y = Weight, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "red") +
  facet_wrap(~ Direction, scales = "free_y") +
  labs(title = "Box-plot of Weights Associated with Positive and Negative Indices",
       x = "Nutrient", y = "Weight (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####
# Create a data frame for the weights
weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.381, 0.018, 0.040, 0.339, 0.219)
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
postiveindex <- ggplot(weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  geom_text(aes(label = paste0(sprintf("%.2f", percentage), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.3) +
  geom_hline(yintercept = cutoff, color = "red", linetype = "dashed") +
  coord_flip() + # Flip axes to make it horizontal
  labs(y = "Weights", x = "Exposure", title = "Weighted Quantile Sum Regression Weights") +
  theme_minimal() # Use a minimal theme for a clean look

# Print the plot
print(postiveindex)

####
# Create a data frame for the weights
weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.019, 0.614, 0.288, 0.052, 0.026)
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
negativeindex <- ggplot(weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  geom_text(aes(label = paste0(sprintf("%.2f", percentage), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.3) +
  geom_hline(yintercept = cutoff, color = "red", linetype = "dashed") +
  coord_flip() + # Flip axes to make it horizontal
  labs(y = "Weights", x = "Exposure", title = "Weighted Quantile Sum Regression Weights") +
  theme_minimal() # Use a minimal theme for a clean look

# Print the plot
print(negativeindex)




# Assuming exposure is correctly defined positive association,
#results1 <- gwqs(RA_status_binary ~ wqs, mix_name = exposure, data = filtered_data31, q = 4, 
#                 validation = 0.6, b = 100, b1_pos = TRUE, b_constr = FALSE, 
#                 family = "binomial", seed = 123)
###
###
###
library(ggplot2)

# Define a consistent color palette
colors <- c("DPHP" = "palegreen4", "BDCPP" = "mediumpurple3", "BCPP" = "orange1", "DBUP" = "salmon1", "BCEP" = "rosybrown1")

# Data for positive index
positive_weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.381, 0.018, 0.040, 0.339, 0.219)
)

# Calculate the total weight
total_weight_positive <- sum(positive_weights_df$mean_weight)

# Calculate the percentage for each weight
positive_weights_df$percentage <- (positive_weights_df$mean_weight / total_weight_positive) * 100

# Calculate the prespecified cutoff
cutoff <- 1 / length(positive_weights_df$mix_name)

# Create the bar plot for the positive index
positive_index_plot <- ggplot(positive_weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = paste0(sprintf("%.2f", percentage), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.3) +
  geom_hline(yintercept = cutoff, color = "red", linetype = "dashed")+ 
  coord_flip()+  # Flip axes to make it horizontal
  labs(y = "WQS index weights", x = "OPEs", title = "Rheumatoid arthritis Positive Index") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) # Remove legend and adjust axis text

# Data for negative index
negative_weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.019, 0.614, 0.288, 0.052, 0.026)
)

# Calculate the total weight
total_weight_negative <- sum(negative_weights_df$mean_weight)

# Calculate the percentage for each weight
negative_weights_df$percentage <- (negative_weights_df$mean_weight / total_weight_negative) * 100

# Create the bar plot for the negative index
negative_index_plot <- ggplot(negative_weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = paste0(sprintf("%.2f", percentage), "%")), 
            position = position_dodge(width = 0.9), vjust = -0.3) +
  geom_hline(yintercept = cutoff, color = "red", linetype = "dashed")+ 
  coord_flip()+  # Flip axes to make it horizontal
  labs(y = "WQS index weights", x = "OPEs", title = "Rheumatoid arthritis Positive Index") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) # Remove legend and adjust axis text

# Print the plots
print(positive_index_plot)
print(negative_index_plot)
library(gridExtra)
# Print both plots in one figure
grid.arrange(positive_index_plot, negative_index_plot, ncol = 1)
####
###
###
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)

# Define a consistent color palette
colors <- c("DPHP" = "palegreen4", "BDCPP" = "mediumpurple3", "BCPP" = "orange1", "DBUP" = "salmon1", "BCEP" = "rosybrown1")

# Data for positive index
positive_weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.381, 0.018, 0.040, 0.339, 0.219),
  index_type = "RA positive Index"
)

# Data for negative index
negative_weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.019, 0.614, 0.288, 0.052, 0.026),
  index_type = "RA negative Index"
)

# Combine data
combined_weights_df <- rbind(positive_weights_df, negative_weights_df)

# Calculate the prespecified cutoff
cutoff <- 1 / length(unique(combined_weights_df$mix_name))

# Calculate the percentage for each weight
combined_weights_df <- combined_weights_df %>%
  group_by(index_type) %>%
  mutate(total_weight = sum(mean_weight)) %>%
  ungroup() %>%
  mutate(percentage = (mean_weight / total_weight) * 100)

# Create the bar plot
index_plot <- ggplot(combined_weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = paste0(sprintf("%.2f", percentage), "%")), 
            position = position_dodge(width = 0.9), hjust = -0.1) +
  geom_hline(yintercept = cutoff, color = "red", linetype = "dashed") +
  coord_flip() +  # Flip axes to make it horizontal
  labs(y = "WQS index weights", x = "OPEs") +
  facet_grid(index_type ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "grey80", color = "grey50"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 1))

# Print the plot
print(index_plot)
####
####
####
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)

# Define a consistent color palette
colors <- c("DPHP" = "palegreen4", "BDCPP" = "mediumpurple3", "BCPP" = "orange1", "DBUP" = "salmon1", "BCEP" = "rosybrown1")

# Data for positive index
positive_weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.381, 0.018, 0.040, 0.339, 0.219),
  index_type = "RA positive Index"
)

# Data for negative index
negative_weights_df <- data.frame(
  mix_name = c("DPHP", "BDCPP", "BCPP", "DBUP", "BCEP"),
  mean_weight = c(0.019, 0.614, 0.288, 0.052, 0.026),
  index_type = "RA negative Index"
)

# Combine data
combined_weights_df <- rbind(positive_weights_df, negative_weights_df)

# Calculate the prespecified cutoff
cutoff <- 1 / length(unique(combined_weights_df$mix_name))

# Calculate the percentage for each weight
combined_weights_df <- combined_weights_df %>%
  group_by(index_type) %>%
  mutate(total_weight = sum(mean_weight)) %>%
  ungroup() %>%
  mutate(percentage = (mean_weight / total_weight) * 100)

# Ensure the order of index_type
combined_weights_df$index_type <- factor(combined_weights_df$index_type, levels = c("RA positive Index", "RA negative Index"))

# Create the bar plot
index_plot <- ggplot(combined_weights_df, aes(x = reorder(mix_name, mean_weight), y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = paste0(sprintf("%.2f", percentage), "%")), 
            position = position_dodge(width = 0.9), hjust = -0.1) +
  geom_hline(yintercept = cutoff, color = "red", linetype = "dashed") +
  coord_flip() +  # Flip axes to make it horizontal
  labs(y = "WQS index weights", x = "OPEs") +
  facet_grid(index_type ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "grey80", color = "grey50"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 1))

# Print the plot
print(index_plot)
####
library(car)
vif_values <- vif(lm(RA_status_binary ~ Gender + Ethnicity + Marital_status + Education + Citizenship + Urinary_creatinine + Alcohol + BMI + Diabetes + Smoking + Age.cat + Activity_level + PIR, data = train_data))
print(vif_values)
