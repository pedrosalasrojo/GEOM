
#         Author: Pedro Salas Rojo and Paolo Brunori
#         Date: 10/2023
#         Name of project: Shapley Functions

library(resample)
library(utils)
library(gtools)
library(data.table)
library(dineq)
library(stringr)
library(stats)


shapley <- function(data, model, vars, type, ntree = 1, wts = NA, 
                    mincri = 0, minbu = 100, resample = 0.632,
                    depname = NA, centiles = 99, order = 4, rel.ineq = TRUE){
  
  time_1 <- Sys.time()
  total <- 2^length(vars)
  num <- 0
  mod_1 <- model
  mod_2 <- update(mod_1, . ~ 1)
  g_results <- NA
  
  if (is.na(wts)) {
    data$weights <- 1
  } else {
    data$weights <- data[[wts]]
  }
  
  data$income <- data[[depname]]
  
  # Use tree to estimate ineq_base ----
  
  ineq_base <- NA
  
  for(iter in seq(1, ntree, 1)) {
    
  #  print(paste0("Baseline Inequality, Iter : ", iter, " of ", ntree, " trees"))
    
    set.seed(iter)
    perm <- NULL
    perm <- data[sample(1:nrow(data), nrow(data)*resample, replace = FALSE),]  # Bag (Necessary for boostrapping)
    
    if(type=="ctree"){
      
    tree <- get_tree(model = mod_1, 
                     data = perm,
                     mincri = as.numeric(mincri), 
                     minbu = as.numeric(minbu))
    
    perm$types <- predict(tree, type="node")
    
    perm <- perm %>%
      group_by(types) %>%
      mutate(y_tilde = weighted.mean(income, weights)) %>%
      ungroup()
    
    if (rel.ineq==TRUE){
      res_ineq <- gini.wtd(perm$y_tilde, perm$weights)
    } else {
      res_ineq <- modi::weighted.var(perm$y_tilde, perm$weights)
    }

    } else if (type=="ols") {
      
      model_ols = lm(model, weights=weights, data=perm)
      perm$y_tilde_ols=fitted.values(model_ols)        
      
      
      if (rel.ineq==TRUE){
        res_ineq <- gini.wtd(perm$y_tilde_ols, perm$weights)
      } else {
        res_ineq <- modi::weighted.var(perm$y_tilde_ols, perm$weights) 
      }
      
      
    }  else if (type=="trafotree") {
    
      trtree <- get_trtree(model = mod_1, 
                              dep = depname,
                              data = perm, 
                              order = order, 
                              mincri = as.numeric(mincri), 
                              minbu = as.numeric(minbu), 
                              centiles = centiles,
                              rel.ineq = rel.ineq,
                              lenv = TRUE, share_lenv = 0.1)
      
      iop <- trtree[["trafodata"]]
      
      if (rel.ineq==TRUE){
        res_ineq <- gini.wtd(iop$y_tilde, iop$weights)
      } else {
        res_ineq <- modi::weighted.var(iop$y_tilde, iop$weights)} 
      
    } else {
      stop("The method in 'type' does not exist in the Shapley Function")
    }
    
    ineq_base <- rbind(ineq_base, res_ineq)
    
  }
  
  ineq_base <- na.omit(as.data.frame(ineq_base))
  
  g_results <- na.omit(rbind(g_results, cbind(ineq_base, "baseline_inequality")))
  names(g_results) <- c("ineq", "comb")
  
  ineq_base <- ineq_base %>%                                    
    summarise_all(funs(mean), na.rm = TRUE) 
  
  ineq_base <- as.numeric(ineq_base)
  
  # Loop to get values for each permutation ----
  
  for(value_comb in seq(1, length(circum) , 1)) {           # Values of combinations
    
    comb <- combn(x = vars,                      # Variables to combine (circumstances in this case)
                  m = value_comb,                # Number of combinations. Iterate this to get all possible combinations
                  simplify = TRUE)               # Eliminate combinations with the same meaning (ab = ba)
    
    mat <- matrix(NA, nrow = ncol(comb), ncol = 2)  # Define matrix of results. First column: name of the combination, 
    # in the second goes the ineq value.
    
    for(roww in seq(1, ncol(comb), 1)) {         # For each one of the combinations in "comb"
      
      gap <- NA
      
      for(iter in seq(1, ntree, 1)) {
        
        if(iter == 1){
          num <- num + 1
        } else {
          num <- num 
        }
        
        print(paste0("Combination: ", num," out of: ", total, ", iter : ", iter, " of ", ntree, " trees"))
        
        set.seed(iter)
        nom <- ""
        perm <- NULL
        perm <- data[sample(1:nrow(data), nrow(data)*resample, replace = FALSE),]  # Bag (Necessary for boostrapping)
        
        permu_vars <- NA
        
        for (name_var in seq (1, nrow(comb), 1)) {
          
          nom <- paste0(nom,"", comb[name_var, roww])  #Generate name with all variables involved in the combination
          
          circ <- comb[name_var, roww]                 #Select variable that is permuted
          
          perm[[circ]] <- 0                         #Substitute values by 0, so they are not used. Loop by rows if more than one circumstance has to be permuted.
          
          permu_vars <- na.omit(rbind(permu_vars, comb[name_var, roww]))
          
        }
        
        # Use tree to estimate the especific contributions to IOp after permutations
        
        if(type=="ctree"){
          
          tree <- get_tree(model = mod_1, 
                           data = perm,
                           mincri = as.numeric(mincri), 
                           minbu = as.numeric(minbu))
          
          perm$types <- predict(tree, type="node")
          
          perm <- perm %>%
            group_by(types) %>%
            mutate(y_tilde = weighted.mean(income, weights)) %>%
            ungroup()
          
          if (rel.ineq==TRUE){
            
            if(value_comb != length(circum)){
              ineq_perm <- round(gini.wtd(perm$y_tilde, perm$weights), 4)
            } else {
              ineq_perm <- 0
            }

          } else {
            
            if(value_comb != length(circum)){
              ineq_perm <- round(modi::weighted.var(perm$y_tilde, perm$weights) , 4)
            } else {
              ineq_perm <- 0
            }
            
          }
          
          } else if (type=="trafotree") {
          
          trtree <- get_trtree(model = mod_1, 
                               dep = depname,
                               data = perm, 
                               order = order, 
                               mincri = as.numeric(mincri), 
                               minbu = as.numeric(minbu), 
                               centiles = centiles,
                               rel.ineq = rel.ineq,
                               lenv = TRUE, share_lenv = 0.1)
          
          iop <- trtree[["trafodata"]]
          
          if (rel.ineq==TRUE){
            
            if(value_comb != length(circum)){
              ineq_perm <- round(gini.wtd(iop$y_tilde, iop$weights), 4)
            } else {
              ineq_perm <- 0
            }
            
          } else {
            
            if(value_comb != length(circum)){
              ineq_perm <- round(modi::weighted.var(iop$y_tilde, iop$weights) , 4)
            } else {
              ineq_perm <- 0
            }
            
          }  
          
          } else if (type=="ols") {
            
            model_ols = lm(model, weights=weights, data=perm)
            perm$y_tilde_ols=fitted.values(model_ols)        

            if (rel.ineq==TRUE){
              
              if(value_comb != length(circum)){
                ineq_perm <- round(gini.wtd(perm$y_tilde_ols, perm$weights), 4)
              } else {
                ineq_perm <- 0
              }
              
            } else {
              
              if(value_comb != length(circum)){
                ineq_perm <- round(modi::weighted.var(perm$y_tilde_ols, perm$weights) , 4)
              } else {
                ineq_perm <- 0
              }
            }
              
          
        } else{
          stop("The method in 'type' does not exist in the Shapley Function")
        }        

        gap <- rbind(gap, ineq_perm)
        
      }
      
      gap <- as.data.frame(gap)
      g_gap <- na.omit(cbind(gap, paste0(permu_vars, collapse = "_")))
      names(g_gap) <- c("ineq", "comb")
      g_results <- rbind(g_results, g_gap)

      gap <- gap %>%                                    
        summarise_all(funs(mean), na.rm = TRUE) 
      
      mat[roww, 1] <- nom                              
      mat[roww, 2] <- round(as.numeric(gap), 3)        
      
      assign(paste0("matrix_", value_comb), mat)
      
    }
  }
  
  # Get Contributions and Shapley Value  ----
  long <- length(circum) - 1                    
  dem <- factorial(length(circum))
  
  # Shapley weights
  shap_we <- matrix(NA, nrow = long, ncol = 1)   
  for(i in seq(1, long, 1)) {
    num <- factorial(i)*factorial(length(circum)-i-1) 
    shap_we[i, 1] <- round(num/dem, 4)               
  }
  
  # Define matrix of results
  
  all_mat <- NA
  for(k in seq(1,length(circum),1)){
    nam <- (paste0("matrix_",k))
    all_mat <- cbind(all_mat, nam)
  }
  all_mat <- all_mat[,-1]
  shapval <- matrix(NA, ncol = 2, nrow = length(circum))
  
  for(i in circum) {
    contributions <- NA
    
    # Contribution baseline - ineq when only 1 is iterated
    vec <- ifelse(grepl(i, matrix_1[,1]), TRUE, FALSE )   
    num <- factorial(0)*factorial(length(circum)-0-1) 
    shap_we <- num/dem
    cont_1 <- shap_we*(ineq_base - as.numeric(matrix_1[match(TRUE,vec), 2])) 
    contributions <- rbind(contributions, cont_1)
    
    # Get remaining contributions (this is, not the baseline)
    for(k in seq(2, length(circum), 1)){
      
      orig = get(paste0("matrix_", k))                #Get the matrix with the objective (origin) matrix
      prev = get(paste0("matrix_", k-1))              #Get the matrix over which we compare inequality
      
      if(k<=long) {
        l <- k-1                                     # k selects the number of circumstances, but the coalition is always one less!
        num <- factorial(l)*factorial(length(circum)-l-1) #Numerator of weights in Shapley depending on the coalition
        shap_we <- num/dem
      } else {
        num <- factorial(long)*factorial(length(circum)-long-1) # Last contribution, when only the final coalition (all permuted)
        shap_we <- num/dem                               # is considered
      }
      
      vec <- ifelse(grepl(i, orig[,1]), TRUE, FALSE)    #Which value in the objective matrix (2, 3, ..., m-1) when (2, 3, ..., m-1) 
      #variables are permuted, include the circumstance in the first column (names)
      
      for(j in seq(1, length(vec), 1)) {                    #For all values combinations (that we have stored in vec)
        
        if(vec[j]==TRUE){                                 # If the name of the circumstance is TRUE in vec
          y <- ifelse(vec[j] == TRUE, orig[j,2], NA)        # Store the correspondent value of inequality
          
          h <- ifelse(vec[j] == TRUE, orig[j,1], NA)        # Store the name of the combination
          h <- stringr::str_remove(h, i)                                 # Remove from the name in h, the letters corresponding to the circumstance
          
          vec2 <- ifelse(grepl(h, prev[,1]), TRUE, FALSE )  # Search in the matrix including one permutation less, the combination 
          # corresponding to the name stored in h
          cont <- shap_we*(as.numeric(prev[match(TRUE,vec2), 2]) - as.numeric(y)) # Estimate inequality as the inequality in the matrix
          # with one permutation less - inequality in the objective matrix. Multiply by weight.
          
          # assign(paste0("cont_", k,"_",j), cont)                  # Store the contribution. 
          contributions <- rbind(contributions, cont)
          
        } else {                                              # If vec does not contain the objective circumstance, ignore.
          NULL
        }
      }
    }
    
    contributions <- na.omit(contributions)
    contributions <- sum(contributions)
    
    row_mat <- match(i,circum)
    
    shapval[row_mat, 1] <- i                           
    shapval[row_mat, 2] <- contributions   
    
  }
  
  time_2 <- Sys.time()
  
  # Show how much time has passed between time_1 and time_2, to see how long does the loop takes
  print(round(time_2 - time_1, 2))
  
  # Check that the residual is zero (complete decomposition)
  (residual <- sum(as.numeric(shapval[,2])) + as.numeric(orig[,2]) - ineq_base)
  
  print(paste0("Residual: ", residual))
  
  # Get Marginal contribution
  (shapval)
  
  # Get Relative contribution
  rel_shapval <- shapval
  
  for(r in seq(1, nrow(shapval), 1)) {
    rel_shapval[r,2] <- 100*as.numeric(shapval[r,2])/ineq_base
  }
  
  (rel_shapval)
  
  # Check that the sum is 100
  print(paste0("Relative sum rel_shapval: ", sum(as.numeric(rel_shapval[,2]))))
  
  print(paste0("Relative sum resid: ", sum(100*as.numeric(orig[,2])/ineq_base)))
  
  print(paste0("Relative sum both: ", sum(as.numeric(rel_shapval[,2])) + 100*as.numeric(orig[,2])/ineq_base))
  
  # Get maximum importance = 100 and index accordingly the other variables
  
  maximp <- max(as.numeric(shapval[,2]))
  rel_shap_max <- shapval
  
  for(r in seq(1, nrow(shapval), 1)) {
    rel_shap_max[r,2] <- 100*as.numeric(shapval[r,2])/as.numeric(maximp)
  }
  
  (rel_shap_max)
  
  return(list(`shapval` = shapval, `rel_shapval` = rel_shapval, `rel_shap_max` = rel_shap_max,
              `all_results` = g_results)) 
  
}

shapley_eop <- function(data, model, vars, ntree = 1, wts = NA, 
                        mincri = 0, minbu = 100, resample = 0.632,
                        depname = NA, centiles = 99, order = 4, 
                        rel.ineq = TRUE, lenv = TRUE, share_lenv = 0.1){
  
  time_1 <- Sys.time()
  total <- 2^length(vars)
  num <- 1
  mod_1 <- model
  mod_2 <- update(mod_1, . ~ 1)
  max_eop<-weighted.mean(data[[depname]])

  if (is.na(wts)) {
    data$weights <- 1
  } else {
    data$weights <- data[[wts]]
  }
  
  data$income <- data[[depname]]
  
  # Use tree to estimate ineq_base ----
  
  ineq_base   <- NA
  eop_base    <- NA
  eopsh_base  <- NA
  
  for(iter in seq(1, ntree, 1)) {
    
    print(paste0("Baseline Inequality, Iter : ", iter, " of ", ntree, " trees"))
    
    set.seed(iter)
    perm <- NULL
    perm <- data[sample(1:nrow(data), nrow(data)*resample, replace = FALSE),]  # Bag (Necessary for boostrapping)
    
    trtree <- get_trtree(model = mod_1, 
                         dep = depname,
                         data = perm, 
                         order = order, 
                         mincri = as.numeric(mincri), 
                         minbu = as.numeric(minbu), 
                         centiles = centiles, 
                         rel.ineq=rel.ineq,
                         lenv = lenv,
                         share_lenv = share_lenv)
    
    iop <- trtree[["trafodata"]]
    
    if (rel.ineq==TRUE){
      res_ineq <- gini.wtd(iop$y_tilde, iop$weights)
    } else {
      res_ineq <- modi::weighted.var(iop$y_tilde, iop$weights)
    } 
    
    if (lenv == TRUE){
      eop <- trtree[["eop"]]
    } else {
      eop <- NA
    }
    
    sh_eop <- trtree[["eopx"]]
    
    
    ineq_base <- rbind(ineq_base, res_ineq)
    eop_base <- rbind(eop_base, eop)
    eopsh_base <- rbind(eopsh_base, sh_eop)
    
  }
  
  ineq_base <- na.omit(as.data.frame(ineq_base))
  eop_base <- na.omit(as.data.frame(eop_base))
  eopsh_base <- na.omit(as.data.frame(eopsh_base))
  
  ineq_base <- ineq_base %>%                                    
    summarise_all(funs(mean), na.rm = TRUE) 
  
  eop_base<-eop_base %>%                                    
    summarise_all(funs(mean), na.rm = TRUE) 
  
  eopsh_base<-eopsh_base %>%                                    
    summarise_all(funs(mean), na.rm = TRUE) 
  
  ineq_base <- as.numeric(ineq_base)
  eop_base <- as.numeric(eop_base)
  eopsh_base <- as.numeric(eopsh_base)
  
  # Loop to get values for each permutation ----
  
  for(value_comb in seq(1, length(circum) , 1)) {           # Values of combinations
    
    comb <- combn(x = vars,                      # Variables to combine (circumstances in this case)
                  m = value_comb,                # Number of combinations. Iterate this to get all possible combinations
                  simplify = TRUE)               # Eliminate combinations with the same meaning (ab = ba)
    
    mat <- matrix(NA, nrow = ncol(comb), ncol = 2)  # Define matrix of results. First column: name of the combination, 
    # in the second goes the ineq value.
    mat_eop  <- mat # same for lower envelop
    mat_eopx <- mat
    
    for(roww in seq(1, ncol(comb), 1)) {         # For each one of the combinations in "comb"
      
      gap <- NA
      gap_eop <- NA
      gap_eopx <- NA
      
      for(iter in seq(1, ntree, 1)) {
        
        if(iter == 1){
          num <- num + 1
        } else {
          num <- num 
        }
        
        print(paste0("Combination: ", num," out of: ", total, ", iter : ", iter, " of ", ntree, " trees"))
        
        set.seed(iter)
        nom <- ""
        perm <- NULL
        perm <- data[sample(1:nrow(data), nrow(data)*resample, replace = FALSE),]  # Bag (Necessary for boostrapping)
        
        for (name_var in seq (1, nrow(comb), 1)) {
          
          nom <- paste0(nom,"", comb[name_var, roww])  #Generate name with all variables involved in the combination
          
          circ <- comb[name_var, roww]                 #Select variable that is permuted
          
          perm[[circ]] <- 0                         #Substitute values by 0, so they are not used. Loop by rows if more than one circumstance has to be permuted.
          
        }
        
        # Use tree to estimate the specific contributions to IOp after permutations
        
        trtree <- get_trtree(model = mod_1, 
                             dep = depname,
                             data = perm, 
                             order = order, 
                             mincri = as.numeric(mincri), 
                             minbu = as.numeric(minbu), 
                             centiles = centiles,
                             rel.ineq= rel.ineq,
                             lenv=lenv,
                             share_lenv = share_lenv)
        
        iop <- trtree[["trafodata"]]
        
        if(value_comb != length(circum)){
          sh_eop <- trtree[["eopx"]]
          #sh_eop <- (1 - (sh_eop-eopsh_base)/(max_eop-eopsh_base))
        } else {
          sh_eop <- max_eop
        }
        
        if (lenv==TRUE){ 
          
          if(value_comb != length(circum)){
            eop_perm <- trtree[["eop"]]
            # eop_perm <- (1 - (eop_perm-eop_base)/(max_eop-eop_base))
          } else {
            eop_perm <- max_eop
          }
          
        } else {
          eop_perm <- NA
        }
        
        if (rel.ineq==TRUE){
          
          if(value_comb != length(circum)){
            ineq_perm <- round(gini.wtd(iop$y_tilde, iop$weights), 4)
          } else {
            ineq_perm <- 0
          }
          
        } else {
          
          if(value_comb != length(circum)){
            ineq_perm <- round(modi::weighted.var(iop$y_tilde, iop$weights) , 4)
          } else {
            ineq_perm <- 0
          }
          
        }  
        
        gap <- rbind(gap, ineq_perm)
        gap_eop <- rbind(gap_eop, eop_perm)
        gap_eopx <- rbind(gap_eopx, sh_eop)
        
      }
      
      # iop
      
      gap <- as.data.frame(gap)
      
      gap <- gap %>%                                    
        summarise_all(funs(mean), na.rm = TRUE) 
      
      mat[roww, 1] <- nom                              
      mat[roww, 2] <- round(as.numeric(gap), 3)        
      
      assign(paste0("matrix_", value_comb), mat)
      
      # eop
      
      gap_eop <- na.omit(as.data.frame(gap_eop))
      gap_eopx <- na.omit(as.data.frame(gap_eopx))
      
      gap_eop <- gap_eop %>%                                    
        summarise_all(funs(mean), na.rm = TRUE) 
      
      gap_eopx <- gap_eopx %>%                                    
        summarise_all(funs(mean), na.rm = TRUE) 
      
      mat_eop[roww, 1] <- nom                              
      mat_eop[roww, 2] <- round(as.numeric(gap_eop), 3)        
      
      mat_eopx[roww, 1] <- nom                              
      mat_eopx[roww, 2] <- round(as.numeric(gap_eopx), 3)        
      
      assign(paste0("matrix_eop_", value_comb), mat_eop)
      assign(paste0("matrix_eopx_", value_comb), mat_eopx)
      
    }
  }
  
  # Get Contributions and Shapley Value  ----
  long <- length(circum) - 1                    
  dem <- factorial(length(circum))
  
  # Shapley weights
  shap_we <- matrix(NA, nrow = long, ncol = 1)   
  for(i in seq(1, long, 1)) {
    num <- factorial(i)*factorial(length(circum)-i-1) 
    shap_we[i, 1] <- round(num/dem, 4)               
  }
  
  # Define matrix of results
  
  shapval <- matrix(NA, ncol = 2, nrow = length(circum))
  shapval_eop <- shapval
  shapval_eopx <- shapval
  
  for(i in circum) {
    
    contributions <- NA
    contributions_eop <- NA
    contributions_eopx <- NA
    
    # Contribution baseline - ineq when only 1 is iterated
    vec <- ifelse(grepl(i, matrix_1[,1]), TRUE, FALSE )   
    vec_eop <- ifelse(grepl(i,matrix_eop_1[,1]), TRUE, FALSE )  
    vec_eopx <- ifelse(grepl(i,matrix_eopx_1[,1]), TRUE, FALSE )  
    
    num <- factorial(0)*factorial(length(circum)-0-1) 
    shap_we <- num/dem
    cont_1 <- shap_we*(ineq_base - as.numeric(matrix_1[match(TRUE,vec), 2])) 
    cont_1_eop <- shap_we*(as.numeric(matrix_eop_1[match(TRUE,vec_eop), 2]) - eop_base) 
    cont_1_eopx <- shap_we*(as.numeric(matrix_eopx_1[match(TRUE,vec_eopx), 2]) - eopsh_base) 
    contributions <- rbind(contributions, cont_1)
    contributions_eop <- rbind(contributions_eop, cont_1_eop)
    contributions_eopx <- rbind(contributions_eopx, cont_1_eopx)
    
    # Get remaining contributions (this is, not the baseline)
    for(k in seq(2, length(circum), 1)){
      
      orig = get(paste0("matrix_", k))                #Get the matrix with the objective (origin) matrix
      prev = get(paste0("matrix_", k-1))              #Get the matrix over which we compare inequality
      orig_eop = get(paste0("matrix_eop_", k))                #Get the matrix with the objective (origin) matrix
      prev_eop = get(paste0("matrix_eop_", k-1))              #Get the matrix over which we compare inequality
      orig_eopx = get(paste0("matrix_eopx_", k))                #Get the matrix with the objective (origin) matrix
      prev_eopx = get(paste0("matrix_eopx_", k-1))              #Get the matrix over which we compare inequality
      
      if(k<=long) {
        l <- k-1                                     # k selects the number of circumstances, but the coalition is always one less!
        num <- factorial(l)*factorial(length(circum)-l-1) #Numerator of weights in Shapley depending on the coalition
        shap_we <- num/dem
      } else {
        num <- factorial(long)*factorial(length(circum)-long-1) # Last contribution, when only the final coalition (all permuted)
        shap_we <- num/dem                               # is considered
      }
      
      vec <- ifelse(grepl(i, orig[,1]), TRUE, FALSE)    #Which value in the objective matrix (2, 3, ..., m-1) when (2, 3, ..., m-1) 
      #variables are permuted, include the circumstance in the first column (names)
      
      vec_eop <- ifelse(grepl(i, orig_eop[,1]), TRUE, FALSE)  
      vec_eopx <- ifelse(grepl(i, orig_eopx[,1]), TRUE, FALSE)  
      
      for(j in seq(1, length(vec), 1)) {                    #For all values combinations (that we have stored in vec)
        
        if(vec[j]==TRUE){                                 # If the name of the circumstance is TRUE in vec
          y <- ifelse(vec[j] == TRUE, orig[j,2], NA)        # Store the correspondent value of inequality
          
          h <- ifelse(vec[j] == TRUE, orig[j,1], NA)        # Store the name of the combination
          h <- stringr::str_remove(h, i)                                 # Remove from the name in h, the letters corresponding to the circumstance
          
          vec2 <- ifelse(grepl(h, prev[,1]), TRUE, FALSE )  # Search in the matrix including one permutation less, the combination 
          # corresponding to the name stored in h
          cont <- shap_we*(as.numeric(prev[match(TRUE,vec2), 2]) - as.numeric(y)) # Estimate inequality as the inequality in the matrix
          # with one permutation less - inequality in the objective matrix. Multiply by weight.
          
          # assign(paste0("cont_", k,"_",j), cont)                  # Store the contribution. 
          contributions <- rbind(contributions, cont)
          
        } else {                                              # If vec does not contain the objective circumstance, ignore.
          NULL
        }
        
        if(vec_eop[j]==TRUE){                                 # If the name of the circumstance is TRUE in vec
          y_eop <- ifelse(vec_eop[j] == TRUE, orig_eop[j,2], NA)        # Store the correspondent value of inequality
          
          h_eop <- ifelse(vec_eop[j] == TRUE, orig_eop[j,1], NA)        # Store the name of the combination
          h_eop <- stringr::str_remove(h_eop, i)                                 # Remove from the name in h, the letters corresponding to the circumstance
          
          vec2_eop <- ifelse(grepl(h_eop, prev_eop[,1]), TRUE, FALSE )  # Search in the matrix including one permutation less, the combination 
          # corresponding to the name stored in h
          cont_eop <- shap_we*(as.numeric(y_eop) - as.numeric(prev_eop[match(TRUE,vec2_eop), 2])) # Estimate inequality as the inequality in the matrix
          # with one permutation less - inequality in the objective matrix. Multiply by weight.
          
          # assign(paste0("cont_", k,"_",j), cont)                  # Store the contribution. 
          contributions_eop <- rbind(contributions_eop, cont_eop)
          
        } else {                                              # If vec does not contain the objective circumstance, ignore.
          NULL
        }
        
        if(vec_eopx[j]==TRUE){                                 # If the name of the circumstance is TRUE in vec
          y_eopx <- ifelse(vec_eopx[j] == TRUE, orig_eopx[j,2], NA)        # Store the correspondent value of inequality
          
          h_eopx <- ifelse(vec_eopx[j] == TRUE, orig_eopx[j,1], NA)        # Store the name of the combination
          h_eopx <- stringr::str_remove(h_eopx, i)                                 # Remove from the name in h, the letters corresponding to the circumstance
          
          vec2_eopx <- ifelse(grepl(h_eopx, prev_eopx[,1]), TRUE, FALSE )  # Search in the matrix including one permutation less, the combination 
          # corresponding to the name stored in h
          cont_eopx <- shap_we*(as.numeric(as.numeric(y_eopx) - as.numeric(prev_eopx[match(TRUE,vec2_eopx), 2]))) # Estimate inequality as the inequality in the matrix
          # with one permutation less - inequality in the objective matrix. Multiply by weight.
          
          # assign(paste0("cont_", k,"_",j), cont)                  # Store the contribution. 
          contributions_eopx <- rbind(contributions_eopx, cont_eopx)
          
        } else {                                              # If vec does not contain the objective circumstance, ignore.
          NULL
        }
      }
    }
    
    contributions <- na.omit(contributions)
    contributions <- sum(contributions)
    
    contributions_eop <- na.omit(contributions_eop)
    contributions_eop <- sum(contributions_eop)
    
    contributions_eopx <- na.omit(contributions_eopx)
    contributions_eopx <- sum(contributions_eopx)
    
    row_mat <- match(i,circum)
    
    shapval[row_mat, 1] <- i                           
    shapval[row_mat, 2] <- contributions   
    
    shapval_eop[row_mat, 1] <- i                           
    shapval_eop[row_mat, 2] <- contributions_eop   
    
    shapval_eopx[row_mat, 1] <- i                           
    shapval_eopx[row_mat, 2] <- contributions_eopx  
    
  }
  
  time_2 <- Sys.time()
  
  # Show how much time has passed between time_1 and time_2, to see how long does the loop takes
  print(round(time_2 - time_1, 2))
  
  # Check that the residual is zero (complete decomposition)
  (residual <- sum(as.numeric(shapval[,2])) + as.numeric(orig[,2]) - ineq_base)
  
  print(paste0("Residual of Trafo Shapley: ", residual))
  
  # Get Marginal contribution
  (shapval)
  (shapval_eop)
  (shapval_eopx)
  
  # Get Relative contribution
  rel_shapval <- shapval
  
  for(r in seq(1, nrow(shapval), 1)) {
    rel_shapval[r,2] <- 100*as.numeric(shapval[r,2])/ineq_base
  }
  
  (rel_shapval)
  
  # Check that the sum is 100
  print(paste0("Relative sum rel_shapval: ", sum(as.numeric(rel_shapval[,2]))))
  
  print(paste0("Relative sum resid: ", sum(100*as.numeric(orig[,2])/ineq_base)))
  
  print(paste0("Relative sum both: ", sum(as.numeric(rel_shapval[,2])) + 100*as.numeric(orig[,2])/ineq_base))
  
  # Get maximum importance = 100 and index accordingly the other variables
  
  maximp <- max(as.numeric(shapval[,2]))
  rel_shap_max <- shapval
  
  maximp_eop <- max(as.numeric(shapval_eop[,2]))
  rel_shap_max_eop <- shapval_eop
  
  maximp_eopx <- max(as.numeric(shapval_eopx[,2]))
  rel_shap_max_eopx <- shapval_eopx
  
  for(r in seq(1, nrow(shapval), 1)) {
    rel_shap_max[r,2] <- round(100*as.numeric(shapval[r,2])/as.numeric(maximp), 3)
    rel_shap_max_eop[r,2] <- round(100*as.numeric(shapval_eop[r,2])/as.numeric(maximp_eop), 3)
    rel_shap_max_eopx[r,2] <- round(100*as.numeric(shapval_eopx[r,2])/as.numeric(maximp_eopx), 3)
  }
  
  (rel_shap_max)
  (rel_shap_max_eop)  
  return(list(`shapval` = shapval, `rel_shapval` = rel_shapval, `rel_shap_max` = rel_shap_max, 
              `shapval_eop` = shapval_eop, `rel_shap_max_eop` = rel_shap_max_eop,
              `shapval_eopx` = shapval_eopx, `rel_shap_max_eopx` = rel_shap_max_eopx,
              `all_results` = g_results)) 
  
}


