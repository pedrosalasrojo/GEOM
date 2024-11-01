
#         Author: Pedro Salas Rojo, Paolo Brunori, and Pedro Torres-Lopez
#         Name of project: Estimate Ex-ante IOp with Trees and Random Forests

rm(list = ls(all.names = TRUE)) 
library(tidyverse)
library(haven)
library(dineq)
library(Metrics)

# Get ex ante functions
source("/Exante_IOp/functions_exante.R")

# Get data ----
data <- read.csv("/Exante_IOp/Exante_IOp/data.csv")

# Store the factored circumstances with the labels you want to show in the plotted tree
data$Sex <- factor(data$sex)
data$Ethnicity <- factor(data$eth)
data$F_Occup <- factor(data$focc)
data$F_Educ <- factor(data$fedu)
data$M_Occup <- factor(data$mocc)
data$M_Educ <- factor(data$medu)

# Write the model including all circumstances of interest
model <- income ~ Sex + Ethnicity + F_Occup + F_Educ + M_Occup + M_Educ

# CTree ----
# Set seed for replicability
set.seed(2)

# Set cross-validation method and number of folds. See package "caret" for
# details
cv5 <- trainControl(method = "cv", number = 5, 
                    verboseIter = FALSE)        

# Define grid of (1-alpha) used to tune the algorithm. See package "caret" for
# details
grid_tr <- expand.grid(mincriterion = seq(0.8, 0.99, 0.01))

# Tune the algorithm with the cross validation and grid defined above

tune <- tune_tree(model = model, 
                  cv = cv5, 
                  grid = grid_tr,
                  data = data)

tune_mincrit <- tune[["mincriterion"]]
tune_rmse <- tune[["RMSE"]]
tune_results <- tune[["results"]]

# Run tree
# Compute tree with the mincriterion selected in the tuning. 
# We need to detach the package "party".

tree <- get_tree(model = model, 
                 data = data,
                 mincri = tune_mincrit)

# Predict types (terminal nodes)
data$types <- predict(tree, type="node")

# Plot tree.
plot_tree(tree, data = data, dep = "income", 
          wts = "weights", norm = TRUE, font = 12)

# Results ctree
# Use terminal nodes to estimate the weighted average by types
data <- data %>%
  group_by(types) %>%
  mutate(y_tilde = weighted.mean(income, weights)) %>%
  ungroup()

# Obtain absolute IOp with Gini and MLD, as well as the number of types and the RMSE
iop <- data %>%
  summarise(name = "Ctree",
            types = length(unique(types)),
            gini = gini.wtd(y_tilde, weights),
            mld = mld.wtd(y_tilde, weights),
            RMSE = as.numeric(tune_rmse))


# CForest ----

# Set seed for replicability
set.seed(3)

# Set cross-validation method and number of folds
cv3 <- trainControl(method = "cv", number = 3, 
                    verboseIter = FALSE)      

# Define grid of "mtry" used to tune the algorithm
grid_tr <- expand.grid(mtry = seq(1, 4, 1))

# Tune the algorithm with the cross validation and grid defined above.
# A better tuning implies a higher grid and number of trees, but it takes more time!
tune <- tune_forest(model = model,
                    data = data,
                    cv = cv3,
                    grid = grid_tr,
                    ntree = 5,
                    mincri = 0)

tune_mtry <- tune[["mtry"]]
tune_rmse <- tune[["RMSE"]]
tune_results <- tune[["results"]]

# Run random forest
# Compute the random forest with the mtry selected in the tuning.
forest <- get_forest(model = model,
                     data = data,
                     ntree = 5,
                     mtry = tune[["mtry"]],
                     mincri = 0)

# Predict income, use indexes to make it faster.
data$index<-sample(1:10, dim(data)[1], replace = TRUE)

data$y_tilde_rf_exante <- NA

for (ind in 1:10){
  data$y_tilde_rf_exante[data$index==ind]<-predict(forest, newdata=data[data$index==ind,])
  print(paste0("Random Forest prediction: ", ind," of ", 10, " folds."))
}

# Results from Random Forest. Note that in Random Forest we do not generate types
iop2 <- data %>%
  summarise(name = "Random Forest",
            types = "No Types",
            gini = gini.wtd(y_tilde_rf_exante, weights),
            mld = mld.wtd(y_tilde_rf_exante, weights),
            RMSE = as.numeric(tune_rmse))

results <- rbind(iop, iop2)

# IOp results from Ctree and Random Forest
results

# Mocking tune for other parameters ----

# Imagine that you do not want to tune the alpha/mincriterion as in the former
# example, but the maximum depth of the tree.

# First, set the grid of values for different depths. In this case, from 1 to 5
depth_try <- seq(1, 5, 1)
set.seed(1)

# Create folds to create train and test sample (see below)
data$flag = sample(1:5, nrow(data), replace = TRUE)
rmse_test <- NA

# Loop to perform k-fold cross validation
for (kf in seq(1, 5, 1)) {
  
  # Create train and test sample
  train <- data[data$flag!=kf,]   
  test  <- data[data$flag==kf,]  
  
  # Loop for each depth you want to try
  for (depth in depth_try) {

  # Get the tree. Note that the data comes from the "train" sample. All other
  # parameters are fixed. If you want to tune several parameters at the same
  # time, it is only a matter of "nesting" one tuning loop inside the other.
  # Note that depending on the grid of values, this can be time consuming.
    
    tree <- get_tree(model = model, 
                     data = train,
                     mincri = 0.99,
                     minbu = 50,
                     maxd = depth)
    
    # Predict, now using newdate = test. This is the "out of sample" prediction
    test$pred <- predict(tree, type="response", newdata = test)
    
    # Estimate RMSE or any other loss function in your test sample, between the real
    # income value and the prediction.
    rmse <- rmse(test$income, test$pred)
    
    # Store and continue with the next fold/value
    rmse <- cbind(kf, depth, rmse)
    rmse_test <- as.data.frame(na.omit(rbind(rmse_test, rmse)))
    
  }
}

# Now, rmse_test contains the RMSE obtained for each possible depth and
# fold. To tune it is only a matter of averaging across depths.

rmse_test <- rmse_test %>%
  dplyr::select(depth, rmse) %>%
  group_by(depth) %>%
  summarise_all(mean) %>%
  ungroup()

# Get the tuned value of depth as the one corresponding to the smallest 
# out of sample RMSE.

tuned_depth <- rmse_test$depth[which.min(rmse_test$rmse)]

# Plot depth values (x) and y (RMSE) to visually check the depth value associated to the smallest
# RMSE. Select that in your tree

plot(x = rmse_test$depth, y = rmse_test$rmse)
abline(v = tuned_depth, col = "red")

tree <- get_tree(model = model, 
                 data = train,
                 mincri = 0.99,
                 minbu = 50,
                 maxd = depth)

# Continue with your analysis.
# Note that this k-fold procedure can be used to tune any other parameter, also
# in the random forests and in most ML algorithms. 
