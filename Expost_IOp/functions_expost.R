
#         Author: Pedro Salas Rojo, Paolo Brunori, and Pedro Torres-Lopez
#         Name of project: Store Ex-Post functions

library(trtf)
library(tram)
library(dineq)
library(lattice)
library(latticeExtra)
library(multcomp)
library(memisc)
library(Matrix)
options(Matrix.warn = FALSE)
library(colorspace)
library(grid)
library(libcoin)
library(inum)
library(partykit)
library(ATR)
library(trtf)
library(mlt)
library(grid)
library(modi)
library(resample)
library(scales)
library(gtools)
library(RColorBrewer)

# Set training function ----

tune_trafotree <- function(model, data, folds = 5, minorder = 2, 
                           maxorder = 10, 
                           mindiff = 0.1, plot = TRUE) {
  
  data_tune <- data %>%
    mutate(flag = sample(1:folds, nrow(data), replace = TRUE))
  
  loss_poli_test <- NA
  loss_type <- NA
  mod_1 <- model
  mod_2 <- update(mod_1, . ~ 1)
  
  for (kf in seq(1, folds, 1)) {
    
    train <- data_tune[data_tune$flag!=kf,]   
    test  <- data_tune[data_tune$flag==kf,]  
    
    tab2_poli_test <- NA

    for (ord in seq(minorder, maxorder, 1)) {
      
      #print(paste0("fold: ",kf,", order: ",ord))
      
      basic_param_train <- tram::BoxCox(mod_2, data = train, order = ord) 
      
      test$llik_poli_test <- as.numeric(logLik(basic_param_train, newdata = test)) 
      
      tab_poli_test <- test %>%
        summarise(val_llik_poli = round(mean(llik_poli_test), 2),
                  ord = ord) 
      
      names(tab_poli_test)[1] <- paste0("LL_",kf)
      tab2_poli_test <- rbind(tab_poli_test, tab2_poli_test)
      
    }
    
    loss_poli_test <- merge(tab2_poli_test, loss_poli_test, all = TRUE) 
  }
  
  loss_poli_test <- loss_poli_test %>%
    dplyr::select(-y)
  loss_poli_test <- na.omit(loss_poli_test)
  loss_poli_test$mean_poli <- round(rowMeans(subset(loss_poli_test, select = names(loss_poli_test[-1])), 
                                             na.rm = TRUE), 2)

  results <- as.data.frame(seq(minorder, maxorder, 1))
  names(results)[1] <- "ord"
  results <- merge(results, loss_poli_test , by = "ord")

  results <- results %>%
    dplyr::mutate(opt = (mean_poli - min(mean_poli))/min(mean_poli),
                  diff = opt - lag(opt),
                  max = ifelse(abs(diff)<=mindiff & opt < 0, TRUE, FALSE))
  
  res <- results %>%
    dplyr::select(ord, mean_poli, opt, diff, max)
  
  if(isTRUE(any(res$max==TRUE, na.rm = TRUE))){

  order <- min(results$ord[results$max==TRUE & !is.na(results$max)])
  loglik  <- min(results$mean_poli[results$max==TRUE & !is.na(results$max)])
  
  if(plot == TRUE) {
    plot <- (ggplot(data=results, mapping=aes(x=ord), show.legend = TRUE)+
               geom_line(mapping=aes(y = mean_poli), size = 3, col = "gray15") + 
               geom_vline(xintercept = min(results$ord[results$max==TRUE & !is.na(results$max)]), 
                          lty = 4, color = "tomato", size = 2) +
               guides(fill = guide_legend(keywidth = 1, keyheight = 1),
                      linetype=guide_legend(keywidth = 3, keyheight = 1),
                      colour=guide_legend(keywidth = 3, keyheight = 1)) +
               labs(x="Order of Bernstein Polynomial", y="Out of Sample LogLikelihood") +
               ggtitle("Out of Sample Log-Likelihood by orders of Bernstein Polynomial") +  
               ylim(min(results$mean_poli)-50, max(results$mean_poli)+50) +
               theme_bw(base_size = 20))
    
    print(plot)  
  } else {
    NULL
  }
  
  print(paste0("The selected order is: ", order))
  print(paste0("The selected order delivers an loglik of: ",loglik))
  
  return(list(`order` = order, `loglik` = loglik,  `res` = res,  `plot` = plot)) 
  
  
  } else {
    
  cat(paste0("WARNING: For the range of order values provided: (",minorder,",",
               maxorder,") there is no value in which the loglikelihood 
               improvement is smaller than ",mindiff,". As a default, the order
               selected is that especified in maxorder:",maxorder,". 
               Explore the 'res' object to get more information"))
  
  plot <- (ggplot(data=results, mapping=aes(x=ord), show.legend = TRUE)+
               geom_line(mapping=aes(y = mean_poli), size = 3, col = "gray15") + 
               geom_vline(xintercept = as.numeric(maxorder), 
                          lty = 4, color = "tomato", size = 2) +
               guides(fill = guide_legend(keywidth = 1, keyheight = 1),
                      linetype=guide_legend(keywidth = 3, keyheight = 1),
                      colour=guide_legend(keywidth = 3, keyheight = 1)) +
               labs(x="Order of Bernstein Polynomial", y="Out of Sample LogLikelihood") +
               ggtitle("Out of Sample Log-Likelihood by orders of Bernstein Polynomial") +  
               ylim(min(results$mean_poli)-50, max(results$mean_poli)+50) +
               theme_bw(base_size = 20))
  print(plot)  
  
  order <- maxorder
  
  return(list(`order` = order, `loglik` = NA,  `res` = res,  `plot` = plot)) 
    
  }
  
}

# Function to get transformation tree ----

get_trtree <- function(model, dep, data, order = 5,
                       mincri = 0.99, minbu = 100, centiles = 99, maxd = Inf,
                       rel.ineq = TRUE, lenv = TRUE, share_lenv = 0.1) {
  
  
  if(nrow(data)<=3000){
    warning("Estimations with less than 3000 observations may deliver unreliable results.", call. = FALSE)
  }
  
  mod_1 <- model
  mod_2 <- update(mod_1, . ~ 1)
  
  basic_param <- tram::BoxCox(mod_2, data = data, order = order,
                              bounds = c(min(data[[dep]]), max(data[[dep]])))
  
  
  trtree <- trafotree(basic_param, formula = mod_1,
                      data = data,
                      teststat = "quad",          
                      testtype = "Bonferroni",    
                      mincriterion = mincri,           
                      minbucket = minbu,
                      maxdepth = maxd)                
  
  data$types <- predict(trtree, newdata = data,  type ="node") 
  
  nti <- length(table(data$types))     
  typ <- unique(data$types)
  
  qi <- as.numeric(table(data$types))/sum(table(data$types))
  
  TR <- matrix(0, nti+1, nrow=centiles)     
  qtl <- round(seq(1-centiles/(centiles+1),
                   centiles/(centiles+1), 
                   1-centiles/(centiles+1)), nchar(centiles)) 
  TR[,nti+1]<-qtl 
  
  QI<-TR
  for (cents in 1:centiles){
    QI[cents,1:nti]<-qi
  }
  
  ni <- length(qtl)   
  spond<-data.frame(table(data$types))
  data$tranche<-NA
  data$yhat<-NA
  
  func_a <- function(x, i){
    ctm_w <- ctm(response = x, todistr = "MaxExtrVal")  
    mlt_w <- mlt(ctm_w, data = i) 
    return(mlt_w)
  }
  
  func_b <- function(x, i){
    ctm_w <- ctm(response = x, todistr = "Normal")  
    mlt_w <- mlt(ctm_w, data = i) 
    return(mlt_w)
  }
  
  func_c <- function(x, i){
    ctm_w <- ctm(response = x, todistr = "Logistic")  
    mlt_w <- mlt(ctm_w, data = i) 
    return(mlt_w)
  }
  
  for (ttp in 1:nti) {
    
    yy <- data[[dep]][data$types==names(table(data$types))[ttp]]  
    inc = as.data.frame(yy)                         
    var_w <- numeric_var("yy", bounds = c(min(yy), max(yy)), support = c(min(yy), max(yy)), 
                         add = c(0, 0), )    
    c(sapply(nd_w <- mkgrid(var_w, ni), range))    
    B_w <- Bernstein_basis(var_w, order = order, ui = "increasing",
                           extrapolate = TRUE) 
    
    mlt <- tryCatch(func_a(x = B_w, i = inc),
                    error = function(e){
                      tryCatch(func_b(x = B_w, i = inc),
                               error = function(e){
                                 func_c(x = B_w, i = inc)
                               }, warning = function(w){
                                 func_c(x = B_w, i = inc)
                               })}, warning = function(w){
                                 tryCatch(func_b(x = B_w, i = inc),
                                          error = function(e){
                                            func_c(x = B_w, i = inc)
                                          }, warning = function(w){
                                            func_c(x = B_w, i = inc)
                                          })})
    
    nd_w$d<-predict(mlt, type = "quantile", prob=qtl) 
    ECDF<-ecdf(inc$yy)        
    inc$tranche<-(round(ECDF(inc$yy), nchar(centiles)))       
    inc$tranche[inc$tranche==0]<-1-centiles/(centiles+1)    
    inc$tranche[inc$tranche==1]<-centiles/(centiles+1)
    TR[,ttp]<-predict(mlt, type = "quantile", prob=qtl)     
    data$tranche[data$types==names(table(data$types))[ttp]]<-round(inc$tranche, nchar(centiles))  
    data$yhat[data$types==names(table(data$types))[ttp]]<-predict(mlt, type = "quantile", prob=inc$tranche)
    spond$Freq[ttp]<-sum(data$types==names(table(data$types))[ttp])
    
  }
  
  spond$Freq<-spond$Freq/sum(spond$Freq)
  pctls<-as.numeric(names(table(data$tranche)))
  data$mj<- NA
  
  for (p in 1:length(pctls)){
    data$mj[data$tranche==pctls[p]]<-weighted.mean(TR[,1:nti][TR[,nti+1]==pctls[p]], spond$Freq)
  }
  
  if (rel.ineq==TRUE){
    data$y_tilde<-data$yhat/data$mj  
  } else {   
    data$y_tilde<-data$yhat - data$mj   
  }
  
  if (lenv==TRUE){
    
    EOP.x<-1:centiles
    EO.p<-1:centiles
    TRE<-TR[,1:nti]
    
    if (dim(as.matrix(TRE))[2]>1){
      
      for (qt in 1:centiles){
        
        TRE[qt,1:nti]<-TR[qt,order(TR[qt,1:nti])] # types ordered according to their outcome in each quantile 
        QI[qt,1:nti]<-QI[qt,order(TR[qt,1:nti])]
        upto<-as.numeric(table(cumsum(QI[qt,1:nti])>share_lenv)[1])
        
        if (upto==nti){
          upto<-upto-1
        }
        
        EOP.x[qt]<-sum(QI[qt,1:(upto+1)]*TRE[qt,1:(upto+1)])/sum(QI[qt,1:(upto+1)]) # weighted average 
        EO.p[qt]<-min(TR[qt,order(TR[qt,1:nti])])
        
        if (QI[qt,1]> share_lenv){ 
          EOP.x[qt] <- EO.p[qt] 
        } # in case the lowest type 
        
      }
      
      EOP<-sum(apply(TRE, 1, min))*(1/(centiles+1))
      EOPX<-sum(EOP.x)*(1/(centiles+1))
      
    } else {
      
      EOP<-mean(as.numeric(TRE))  
      EOPX<-mean(as.numeric(TRE))  
      
    }
    
  } else {
    
    EOP <- NA
    EOPX<- NA
    
  }
  
  data$eop_thres <- NA
  data$eopx_thres <- NA
  
  for (p in 1:length(pctls)){
    data$eopx_thres[data$tranche==pctls[p]] <- EOP.x[p] 
    data$eop_thres[data$tranche==pctls[p]] <- EO.p[p] 
    
  }

  data$eop_belong <- ifelse(data$yhat <= data$eop_thres, 1, 0)
  data$eopx_belong <- ifelse(data$yhat <= data$eopx_thres, 1, 0)
  
  return(list(`trafodata` = data, `tree` = trtree, `qtl` = qtl, `tr` = TR, `eop` = EOP, `eopx` = EOPX )) 
  
}

# Function to get nice log nodes ----

node_dense<- function(obj,
                      log = T,
                      width = 0.5,
                      yscale = NULL,
                      ylines = 0,
                      cex = 0.5,
                      id = TRUE,
                      mainlab = NULL, 
                      gp = gpar(),
                      colordens = "#E0112B",
                      font_t = 20) {
  if (log == T){
    y <- log(obj$fitted[["(response)"]]) # Log-transform values
    stopifnot(is.numeric(y))
    
    ymin <- min(y)
    ymax <- max(y)
  } else {
    y <- obj$fitted[["(response)"]] # Normal values
    stopifnot(is.numeric(y))
    
    ymin <- min(y)
    ymax <- max(y)
  }
  
  total <- length(y)
  
  # Weights 
  data <<- data # Data is a global variable: it should be the data frame with which we are working
  y_total <- obj$fitted[["(response)"]]
  obj$fitted[["(weights)"]] <- data$weights 
  mean_y <- weighted.mean(y_total, w = obj$fitted[["(weights)"]])
  
  types <- length(unique(obj$fitted[["(fitted)"]]))
  first <- min(unique(obj$fitted[["(fitted)"]]))
  
  print
  
  if (is.null(yscale)) 
    yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
  
  # panel function for density of log income
  rval <- function(node) {
    
    ## extract data
    nid <- id_node(node)
    dat <- data_party(obj, nid)
    
    if (log == T){
      yn <- log(dat[["(response)"]]) # Log-transform values in node
    } else{
      yn <- dat[["(response)"]]
    }
    
    nobs <- length(yn)
    
    wn  <- dat[["(weights)"]]
    if(is.null(wn)) wn <- rep(1, length(yn))
    
    dens <- density(yn)
    
    top_vp <- viewport(layout = grid.layout(nrow   = 2,
                                            ncol   = 1,
                                            widths = grid::unit(c(1), 
                                                                c("null")),  
                                            heights = grid::unit(c(1, 1),
                                                                 c("lines", "null")
                                            )
    ),
    width  = grid::unit(0.8, "npc"), 
    height = grid::unit(1, "npc") - grid::unit(0.5, "lines"),
    name   = paste("node_boxplot", nid, sep = ""),
    gp     = gp)
    
    pushViewport(top_vp)
    
    grid.rect(gp = gpar(fill = NA,
                        col  = 0))
    
    plot <- viewport(layout.pos.col = 1,
                     layout.pos.row = 2,
                     name           = paste0("node_boxplot", nid, "plot"),
                     clip           = T)
    pushViewport(plot)
    
    if (nid == first) {
      x_axis_text <- element_text(size = font_t)
      x_axis_line <- element_line()
      x_axis_tick <- element_line()
    } else {
      x_axis_text <- element_blank()
      x_axis_line <- element_blank()
      x_axis_tick <- element_blank()
    }
    
    a <- ggplot(as.data.frame(yn)) + aes(x = yn) +
      geom_density(adjust = 2,
                   fill   = colordens) +
      geom_vline(xintercept = log(mean_y), size = 2) +
      theme_classic() +
      xlab("") + 
      ylab("") +
      xlim(max(c(ymin - 0.5, 1)),  ymax + 0.5) +
      theme(
        axis.title.y     = element_text(
          size = font_t
        ),
        axis.ticks.y     = element_blank(),
        axis.text.y      = element_blank(),
        axis.line.y      = element_blank(),
        axis.text.x      = element_text(size=font_t*3),
        axis.line.x      = x_axis_line,
        axis.ticks.x     = x_axis_tick, 
        panel.background = element_rect(fill = 'transparent'),
        plot.background  = element_rect(fill = 'transparent', color=NA),
        plot.margin      = margin(t    = 0, 
                                  r    = 5,
                                  b    = 0,
                                  l    = 0,
                                  unit = "pt")
      ) +
      annotate(
        "text",
        x     = max(c(ymin - 0.5, 1)),
        y     = max(dens$y) - sd(dens$y)*1.5,
        label = paste("Type ", nid, "\n",
                      "Pop. = ", round(nobs/total*100, 2), "%",
                      "\ny = ", round(weighted.mean(dat[["(response)"]], w = wn) / mean_y, 2), 
                      sep = ""),
        size  = font_t,
        hjust = 0
      )
    
    print(a, newpage=FALSE)
    
    upViewport(2)
  }
  
  return(rval)
  
}
class(node_dense) <- "grapcon_generator"

# Arrange colors for plots ----

colplot <- function(data, dep, types, grouping_var){
  
  set.seed(1)
  
  data$dep <- data[[dep]]
  data$tps <- data[[types]]
  data$gro <- data[[grouping_var]]
  
  data <- data %>%
    group_by(tps) %>%
    arrange(gro, decreasing = FALSE) %>%
    mutate(color_types = paste0(unique(gro[order(gro, decreasing = TRUE)]), 
                                collapse = "/"),
           mu = mean(dep)) 
  
  n_types <- sort(names(table(data$color_types)))
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  rownam <- c("Dark2", "Paired")
  listcol <- c("#E04F1D", "#7AE028", "#E0112B", "#28E0B7", "#BD00DB",
               "#6212E0", "#529415", "#000000", "#00B9F0", "#F0B05F", 
               "#EF8F8B", "#E5A5FF", "#DEE18E", "#945900", "#E5C61E")
  
  if(length(unique(n_types))<=length(listcol)){
    
    col=sample(listcol, length(unique(n_types)), replace = FALSE)
    
  } else {
    
    color <- listcol
    col=sample(listcol, length(unique(n_types)) - length(listcol), replace = TRUE)
    col <- c(col, color)
    
  }
  
  return(list(`plot_data` = data, `col` = col, `n_types` = n_types))
  
}
# Plot Mountains ----

plot_mountain <- function(data, dep, group, fill, col, 
                         limit_x, x_lab = "", y_lab = "",  labs = "", 
                         legend_pos = "bottom", title = ""){
  
  print(paste0("Note: Density plots use 'logs' as default."))
  
  data$dep <- data[[dep]]
  data$group <- data[[group]]
  data$fill <- data[[fill]]
  
  den <- density(log(data$dep))
  
  
  maxden <- max(den$y)
  
  myplot <- ggplot(data = data, aes(x=log(data$dep), y=..count.., 
                                    group = data$group, fill = factor(data$fill))) +
    geom_density(alpha = 0.9, bw = 0.25, position = "stack", trim = FALSE)
  
  pg <- ggplot_build(myplot)
  
  max_val <- pg[["layout"]][["panel_scales_y"]][[1]][["range"]][["range"]][2]
  
  myplot <- myplot +
    scale_fill_manual(values = col) +
    xlim(limit_x) +
    scale_y_continuous(labels= function(x)round(x*maxden/max_val, 2))+
    xlab(x_lab) + ylab(y_lab) +
    labs(fill = labs, color = "Type") +
    ggtitle(title) +
    theme_bw(base_size =19)
  myplot <- myplot + 
    theme(legend.position=legend_pos)  
  
  return(myplot)
  
}

# Plot Ponytail ----

plot_ponytail <- function(data, dep, quart_var, tr, 
                          nti, group, fill, col, limit_x, 
                          x_lab = "", y_lab = "", labs = "", 
                          legend_pos = "", title = ""){
  
  quart <- as.data.frame(quart_var)
  
  depen <- data[[dep]]

  myplot <- ggplot() +
    stat_ecdf(data = data, aes(depen, group = .data[[group]], 
                               color = factor(.data[[fill]])), size = 1) +
    scale_color_manual(values = col) +
    xlab(x_lab) + ylab(y_lab) +
    labs(fill = labs) +
    ggtitle(title) +
    theme_bw(base_size =19)
  
  for (ttp in 1:nti){
    quart$VV<-(tr[,ttp])

    myplot <- myplot+geom_line(data=quart, aes(x=VV,y=quart_var), color = "black",
                   linetype = "dashed") 
  }
  
  myplot <- myplot + theme(legend.position = legend_pos) + labs( color = "Types")+
    coord_cartesian(xlim=limit_x)

  
  return(myplot)
  
}

# Get log interpolations ----

get_tr <- function(dep, types, data, order = 5, centiles = 99, 
                   rel.ineq = TRUE, lenv = TRUE, share_lenv = 0.1) {
  
  nti <- length(table(data[[types]]))     
  TR <- matrix(0, nti+1, nrow=centiles)     
  qtl <- round(seq(1-centiles/(centiles+1),
                   centiles/(centiles+1), 
                   1-centiles/(centiles+1)), nchar(centiles)) 
  TR[,nti+1]<-qtl 
  ni <- length(qtl)   
  spond<-data.frame(table(data[[types]]))
  data$tranche<-NA
  data$yhat<-NA
  
  qi <- as.numeric(table(data[[types]]))/sum(table(data[[types]]))
  
  QI<-TR
  for (cents in 1:centiles){
    QI[cents,1:nti]<-qi
  }
  
  func_a <- function(x, i){
    ctm_w <- ctm(response = x, todistr = "MaxExtrVal")  
    mlt_w <- mlt(ctm_w, data = i) 
    return(mlt_w)
  }
  
  func_b <- function(x, i){
    ctm_w <- ctm(response = x, todistr = "Normal")  
    mlt_w <- mlt(ctm_w, data = i) 
    return(mlt_w)
  }
  
  func_c <- function(x, i){
    ctm_w <- ctm(response = x, todistr = "Logistic")  
    mlt_w <- mlt(ctm_w, data = i) 
    return(mlt_w)
  }
  
  for (ttp in 1:nti) {
    
    yy <- data[[dep]][data[[types]]==names(table(data[[types]]))[ttp]]  
    inc = as.data.frame(yy)                         
    var_w <- numeric_var("yy", bounds = c(min(yy), max(yy)), 
                         support = c(min(yy), max(yy)), 
                         add = c(0, 0), )    
    c(sapply(nd_w <- mkgrid(var_w, ni), range))    
    B_w <- Bernstein_basis(var_w, order = order, ui = "increasing",
                           extrapolate = TRUE) 
    
    mlt <- tryCatch(func_a(x = B_w, i = inc),
                    error = function(e){
                      tryCatch(func_b(x = B_w, i = inc),
                               error = function(e){
                                 func_c(x = B_w, i = inc)
                               }, warning = function(w){
                                 func_c(x = B_w, i = inc)
                               })}, warning = function(w){
                                 tryCatch(func_b(x = B_w, i = inc),
                                          error = function(e){
                                            func_c(x = B_w, i = inc)
                                          }, warning = function(w){
                                            func_c(x = B_w, i = inc)
                                          })})
    
    nd_w$d<-predict(mlt, type = "quantile", prob=qtl) 
    ECDF<-ecdf(inc$yy)        
    inc$tranche<-(round(ECDF(inc$yy), nchar(centiles)))       
    inc$tranche[inc$tranche==0]<-1-centiles/(centiles+1)    
    inc$tranche[inc$tranche==1]<-centiles/(centiles+1)
    TR[,ttp]<-predict(mlt, type = "quantile", prob=qtl)     
    data$tranche[data[[types]]==names(table(data[[types]]))[ttp]]<-round(inc$tranche, nchar(centiles))  
    data$yhat[data[[types]]==names(table(data[[types]]))[ttp]]<-predict(mlt, type = "quantile", prob=inc$tranche)
    spond$Freq[ttp]<-sum(data[[types]]==names(table(data[[types]]))[ttp])
    
  }
  
  spond$Freq<-spond$Freq/sum(spond$Freq)
  pctls<-as.numeric(names(table(data$tranche)))
  data$mj<- NA
  
  for (p in 1: length(pctls)){
    data$mj[data$tranche==pctls[p]]<-weighted.mean(TR[,1:nti][TR[,nti+1]==pctls[p]], spond$Freq)
  }
  
  if (rel.ineq==TRUE){
    data$y_tilde<-data$yhat/data$mj  
  } else {   
    data$y_tilde<-data$yhat - data$mj   
  }
  
  if (lenv==TRUE){
    
    EOP.x<-1:centiles
    EO.p<-1:centiles
    TRE<-TR[,1:nti]
    
    if (dim(as.matrix(TRE))[2]>1){
      
      for (qt in 1:centiles){
        
        TRE[qt,1:nti]<-TR[qt,order(TR[qt,1:nti])] # types ordered according to their outcome in each quantile 
        QI[qt,1:nti]<-QI[qt,order(TR[qt,1:nti])]
        upto<-as.numeric(table(cumsum(QI[qt,1:nti])>share_lenv)[1])
        
        if (upto==nti){
          upto<-upto-1
        }
        
        EOP.x[qt]<-sum(QI[qt,1:(upto+1)]*TRE[qt,1:(upto+1)])/sum(QI[qt,1:(upto+1)]) # weighted average 
        EO.p[qt]<-min(TR[qt,order(TR[qt,1:nti])])
        
        if (QI[qt,1]> share_lenv){ 
          EOP.x[qt] <- EO.p[qt] 
        } # in case the lowest type 
        
      }
      
      EOP<-sum(apply(TRE, 1, min))*(1/(centiles+1))
      EOPX<-sum(EOP.x)*(1/(centiles+1))
      
    } else {
      
      EOP<-mean(as.numeric(TRE))  
      EOPX<-mean(as.numeric(TRE))  
      
    }
    
  } else {
    
    EOP <- NA
    EOPX<- NA
    
  }
  
  data$eop_thres <- NA
  data$eopx_thres <- NA
  
  for (p in 1:length(pctls)){
    data$eopx_thres[data$tranche==pctls[p]] <- EOP.x[p] 
    data$eop_thres[data$tranche==pctls[p]] <- EO.p[p] 
    
  }
  
  data$eop_belong <- ifelse(data$yhat <= data$eop_thres, 1, 0)
  data$eopx_belong <- ifelse(data$yhat <= data$eopx_thres, 1, 0)
  
  return(list(`trafodata` = data, `qtl` = qtl, `tr` = TR, `eop` = EOP, `eopx` = EOPX )) 
  
}
