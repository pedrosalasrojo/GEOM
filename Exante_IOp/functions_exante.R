
#         Author: Pedro Salas Rojo, Paolo Brunori, and Pedro Torres-Lopez
#         Name of project: Ex-Ante IOp Functions

library(caret)
library(partykit)
library(glmnet)
library(grDevices)
library(stringr)
library(dineq)
library(gtools)

# Function to tune the tree ----

tune_tree <- function(data, model, cv, grid, minbu = 50) {
  
  print(paste0("This function temporarily loads the package 'party'."))
  
  library(party)
  
  set.seed(1)
  
  tr_train <- train(model,
                    data = data, 
                    method = "ctree", 
                    trControl = cv,  
                    tuneGrid = grid,
                    controls = ctree_control(testtype = "Bonferroni", 
                                             minbucket = minbu))
  
  results <- tr_train[["results"]]
  mincri <- tr_train[["bestTune"]][["mincriterion"]]    
  rmse <- round(mean(tr_train[["resample"]][["RMSE"]]), 2)
  
  detach("package:party", unload = TRUE)
  
  print(paste0("The selected alpha is: ", (1-mincri)))
  print(paste0("The selected mincriterion is: ", (mincri)))
  print(paste0("The selected alpha delivers an RMSE of: ",rmse))
  
  return(list(`mincriterion` = mincri, `RMSE` = rmse, `results` = results)) 
}

# Function to get the tree ----

get_tree <- function(data, model, mincri = 0.99, minbu = 50, maxd = Inf) {
  
  set.seed(1)
  
  tree <- partykit::ctree(model,
                          data = data, 
                          control = ctree_control(testtype = "Bonferroni", 
                                                  mincriterion = mincri,
                                                  minbucket = minbu,
                                                  maxdepth = maxd))
  return(tree)
}

# Function to plot trees ----

plot_tree <- function(data, tree, dep, wts, 
                      norm = FALSE, font = 6) {
  
  ct_node <- as.list(tree$node)
  data$types <- predict(tree, type = "node")
  
  data$dep <- data[[dep]]
  data$wts <- data[[wts]]
  
  pred <- data %>%
    group_by(types) %>%
    mutate(x = stats::weighted.mean(x = dep, 
                                    w = wts)) %>%
    dplyr::select(types, x)  %>%
    summarise_all(funs(mean), na.rm = TRUE) %>%
    ungroup()
  
  a <- data %>%
    mutate(m = stats::weighted.mean(x = dep, wts))
  
  mean_pop <- round(mean(a$m),3)
  
  pred <- as.data.frame(pred)
  qi <- pred
  
  for (t in 1:length(qi[,1])){
    typ<-as.numeric(names(table(data$types)[t]))
    qi[t,2]<-length(data$types[data$types==typ])/length(data$types) 
  }
  
  if (norm == TRUE) {
    print(paste0("1 corresponds to the weighted mean: ", mean_pop))
    pred$x <- pred$x/mean_pop
    dig <- 3
  } else {
    pred$x  <-pred$x
    dig <- 0
  }
  
  for(t in 1:nrow(pred)) {
    ct_node[[pred[t,1]]]$info$prediction <- as.numeric(paste(format(round(pred[t, -1], 
                                                                          digits = dig), nsmall = 2)))
    ct_node[[pred[t,1]]]$info$nobs       <- as.numeric(paste(format(round(100*qi[t, -1]  , 
                                                                          digits = 2), nsmall = 2)))
  }
  
  tree$node <- as.partynode(ct_node)
  
  plot((tree), terminal_panel=node_terminal, gp = gpar(fontsize = font),
       tp_args = list(FUN = function(node) 
         c("Rel. Type Mean",node$prediction, "Pop. Share (%)", node$nobs)))
}

# Function to tune the forest ----

tune_forest <- function(data, model, cv, grid, 
                        ntree = 50, mincri = 0, minbu = 10) {
  
  set.seed(1)
  
  print(paste0("This function temporarily loads the package 'party'."))
  library(party)
  
  rf_train <- caret::train(model,
                           data,
                           method = "cforest",
                           trControl = cv,
                           tuneGrid = grid,
                           controls =  cforest_control(ntree = ntree,
                                                       mincriterion = mincri,
                                                       testtype = "Bonferroni",
                                                       minbucket = minbu))
  
  results <- rf_train[["results"]]
  mtry  <- rf_train[["bestTune"]][["mtry"]]        
  rmse <- round(mean(rf_train[["resample"]][["RMSE"]]), 2)  
  detach("package:party", unload = TRUE)
  
  print(paste0("The selected mtry is: ", mtry))
  print(paste0("The selected mtry delivers an RMSE of: ",rmse))
  return(list(`mtry` = mtry, `RMSE` = rmse, `results` = results))
}

# Function to get the forest ----

get_forest <- function(data, model, ntree = 50, 
                       mincri = 0, minbu = 10,
                       mtry = "default") {
  
  set.seed(1)
  
  if(mtry == "default"){
    mvar <- ceiling((stringi::stri_count(as.vector(as.character(mod[3])), fixed = "+"))^0.5)
  } else {
    mvar <- mtry
  }
  
  forest <- partykit::cforest(model,
                              data = data,
                              ntree = ntree,
                              mtry = mvar,
                              trace = TRUE,
                              control = ctree_control(testtype = "Bonferroni",
                                                      mincriterion = mincri,
                                                      minbucket = minbu))
  
  return(forest)
}
