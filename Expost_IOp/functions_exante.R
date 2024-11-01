
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

tune_tree <- function(model, data, cv, grid, minbu = 100, plot = TRUE) {
  
  print(paste0("This function temporarily loads the package 'party'."))
  library(party)
  
  tr_train <- caret::train(model,
                          data = data, 
                          method = "ctree", 
                          trControl = cv,  
                          tuneGrid = grid,
                          controls = ctree_control(testtype = "Bonferroni", 
                                                   teststat = "quad", 
                                                   minbucket = minbu))
  
  if(plot == TRUE) {
    p <- print(plot(tr_train, main = "RMSE by 1 - alpha", 
                    xlab = "Mincriterion (1-alpha)", lwd=3.0))      
    print(p)      
  } else {
    p <- NULL
  }
  
  results <- tr_train[["results"]]
  mincri <- tr_train[["bestTune"]][["mincriterion"]]    
  rmse <- round(mean(tr_train[["resample"]][["RMSE"]]), 2)
  
  detach("package:party", unload = TRUE)
  
  print(paste0("The selected alpha is: ", (1-mincri)))
  print(paste0("The selected mincriterion is: ", (mincri)))
  print(paste0("The selected alpha delivers an RMSE of: ",rmse))
  
  return(list(`mincriterion` = mincri, `RMSE` = rmse, 
              `plot` = p, `results` = results)) 
}

# Function to get the tree ----

get_tree <- function(model, data, mincri = 0.95, minbu = 100, maxd = Inf) {
  
  tree <- partykit::ctree(model,
                          data = data, 
                          control = ctree_control(testtype = "Bonferroni", 
                                                  teststat = "quad", 
                                                  mincriterion = mincri,
                                                  minbucket = minbu,
                                                  maxdepth = maxd))
  return(tree)
}

# Function to plot trees ----

plot_tree <- function(tree, data, dep, wts=NA, norm = FALSE) {
  
  ct_node <- as.list(tree$node)
  data$types <- predict(tree, type = "node")

  data <- data %>%
    mutate(weights_sam = ifelse(is.na({{wts}}), 1, {{wts}}))
  
  pred <- data %>%
    group_by(types) %>%
    mutate(x = stats::weighted.mean(x = {{dep}}, 
                                           w = weights_sam)) %>%
    summarise_all(funs(mean), na.rm = TRUE)    %>%
    ungroup() %>%
    dplyr::select(types, x) 
  
 a <- data %>%
   mutate(m = stats::weighted.mean(x = {{dep}}, weights_sam))
 
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
  
  plot(tree,  terminal_panel=node_terminal, 
       tp_args = list(FUN = function(node) 
         c("Exp. outcome",node$prediction, "Pop. Share (%)", node$nobs)))
  
}

# Function to tune the forest ----

tune_forest <- function(model, data, cv, grid, ntree = 50, 
                        mincri = 0, minbu = 100,
                        plot = TRUE) {
  
  print(paste0("This function temporarily loads the package 'party'."))
  library(party)
  
  rf_train <- caret::train(model,
                           data,
                           method = "cforest",
                           trControl = cv,
                           tuneGrid = grid,
                           controls =  cforest_control(ntree = ntree,
                                                       mincriterion = mincri,
                                                       teststat = "quad",
                                                       testtype = "Bonferroni",
                                                       minbucket = minbu))
  
  if(plot == TRUE) {
    p <- plot(rf_train, main = "RMSE by mtry", lwd=3.0)   
  } else {
    p <- NULL   
  }
  
  results <- rf_train[["results"]]
  mtry  <- rf_train[["bestTune"]][["mtry"]]        
  rmse <- round(mean(rf_train[["resample"]][["RMSE"]]), 2)  
  detach("package:party", unload = TRUE)
  
  print(paste0("The selected mtry is: ", mtry))
  print(paste0("The selected mtry delivers an RMSE of: ",rmse))
  return(list(`mtry` = mtry, `RMSE` = rmse, `plot` = p, `results` = results))
}

# Function to get table with circumstance shares ----

table_circ <- function(data, dep, wts, circum, types){
  
  data$dep <- data[[dep]]
  data$wts <- data[[wts]]
  data$types <- data[[types]]
  
  lst = list()
  
  n = 1
  for(i in circum){
    
    data$circ_var_loop <- data[[i]]
    
    meanpop <- weighted.mean(data$dep, data$wts)
    allpop <- nrow(data)
    
    # Income by circumstances
    output_circ <- data%>%
      group_by(circ_var_loop) %>%
      summarise(mean = weighted.mean(dep, wts)/meanpop) %>%
      dplyr::select(mean)
    
    output_circ <- round(as.vector(t(output_circ)), 2)
    
    # Income by types
    output_type <- data %>%
      group_by(types) %>%
      summarise(mean = weighted.mean(dep, wts)/meanpop) %>%
      dplyr::select(mean)
    
    output_type <- round(as.vector(t(output_type)), 2)
    output_type <- c("Output", output_type)
    
    # Population by circumstances
    table_popcirc <- round(100*prop.table(table(data$circ_var_loop)), 2)
    
    # Population by types
    table_poptype <- round(100*prop.table(table(data$types_tree)), 2)
    
    # Table
    table_prop <- round(100*prop.table(table(data$circ_var_loop, data$types)), 2)
    table_prop <- rbind(table_prop, colSums(table_prop))
    table_prop <- cbind(table_prop, rowSums(table_prop))
    table_prop <- t(cbind(c(output_circ,"Type Share"), as.matrix(table_prop)))
    table_prop <- (cbind(c(output_type,"Circ Share"), as.matrix(table_prop)))
    
    # Order by increasing expected outcome
    row_1 <- table_prop[1,]
    row_last <- table_prop[nrow(table_prop),]
    
    tab <- table_prop[-1,]
    tab <- as.data.frame(tab[-nrow(tab),])
    name <- names(tab)[1]
    
    tab <- tab %>%
      arrange(get(name))
    
    table_prop <- rbind(row_1, as.matrix(tab), row_last)
    
    lst[[n]] = assign(paste0("Table_Circ_",i), table_prop)
    n = n+1
    
  }
  
  return(lst)
  
}

# Function to get the forest ----

get_forest <- function(model, data, ntree = 50, mtry = "default", 
                       mincri = 0, minbu = 100) {
  
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
                                                      teststat = "quad",
                                                      mincriterion = mincri,
                                                      minbucket = minbu))
  
  return(forest)
}

# Function to get relative importance from forest ----

rel_imp <- function(forest) {
  imp <- partykit::varimp(forest)
  relimp <- round(100*imp/max(imp), 2)
  relimp <- relimp[order(-relimp)]
  print(relimp)
  return(`relimp` = relimp)
}

# Function to get Van Kerm's correction ----

vk_exante <- function(data, dep, types, wts = NA, reps = 10){
  
  gini<-1:reps
  mld<-1:reps
  
  data$inc <- data[[dep]]
  data$types <- data[[types]]
  
  if(is.na(wts)){
    data$weights <- 1
  } else {
    data$weights <- data[[wts]]
  }

  for (b in 1:reps){
    set.seed(b)
    datvk<-data[c("weights", "types", "inc")]
    datvk$inc<-gtools::permute(datvk$inc)
    datvk <- datvk %>%
      group_by(types) %>%
      mutate(mean_yc = weighted.mean(inc, weights)) %>%
      ungroup()
    gini[b] = round(gini.wtd(datvk$mean_yc, datvk$weights),4)
    mld[b] = round(mld.wtd(datvk$mean_yc, datvk$weights),4)
  }
  vkrn_correction<-c(mean(gini),mean(mld))
  return(`vkrn_correction` = vkrn_correction)
}

# Function to tune the LASSO ----

tune_lasso <- function(model, data, cv, wts = NULL, 
                       grid, plot = TRUE, fam = "gaussian") {
  
  data <- data
  
  if(is.null(wts)){
    data$weights <- rep(1, times = nrow(data))
  } else {
    data$weights <- data[[wts]]
  }
  
  lasso_tr <- caret::train(model,
                           data = data,
                           method = "glmnet",
                           weights = weights,
                           family = fam,
                           trControl = cv,
                           tuneGrid = grid)
  
  if(plot == TRUE) {
    print(plot(y = lasso_tr$results$RMSE, x = lasso_tr$results$lambda,
               main = "RMSE by lambda value", xlab = "Lambda", ylab = "RMSE"))
    print(abline(v = lasso_tr[["bestTune"]][["lambda"]]))  
    plot <- recordPlot()
  } else {
    plot <- NULL
  }
  
  results <- lasso_tr[["results"]]
  lambda <- lasso_tr[["bestTune"]][["lambda"]]         
  rmse <- round(mean(lasso_tr[["resample"]][["RMSE"]]), 2)  
  
  print(paste0("The selected lambda is: ", lambda))
  print(paste0("The selected lambda delivers an RMSE of: ",rmse))
  return(list(`lambda` = lambda, `RMSE` = rmse, 
              `plot` = plot, `results` = results))
}

# Function to get the LASSO ----

get_lasso <- function(dep, indep, data, lambda, wts = NULL, 
                      plot = TRUE, fam = "gaussian") {
  
  data <- data
  
  if(is.null(wts)){
    data$weights <- rep(1, times = nrow(data))
  } else {
    data$weights <- data[[wts]]
  }
  
  lasso_mod <- glmnet(vec, dep, alpha=1, 
                      weights = data$weights,
                      family = fam)
  
  if(plot == TRUE) {
    plot(lasso_mod, xvar = "lambda")
    abline(v=log(tune[["lambda"]]), lty="dashed", col="black")
    plot <- recordPlot()
    print(plot)
  } else {
    plot <- NA
  }
  
  lasso <- glmnet(vec, dep, alpha=1, lambda = lambda, 
                  weights = data$weights)
  coeff <- lasso$beta
  
  return(list(`lasso` = lasso, `coeff` = coeff, `plot` = plot))
  
}

