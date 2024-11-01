
#         Author: Pedro Salas Rojo and Paolo Brunori
#         Date: 10/2023
#         Name of project: Estimate Ex-post IOp with several methods

rm(list = ls(all.names = TRUE)) 
library(tidyverse)
library(haven)
library(dineq)

# Get ex ante functions
source("C:/Users/SALASROJ/Dropbox/BRUNORI&SALAS_/database_github_etc/TO SHARE/Expost_IOp/functions_expost.R")

# Get data ----
data <- read.csv("C:/Users/SALASROJ/Dropbox/BRUNORI&SALAS_/database_github_etc/TO SHARE/Expost_IOp/data.csv")

# Store the factored circumstances with the labels you want to show in the plotted tree

data$Sex <- factor(data$sex)
data$Ethnicity <- factor(data$eth)
data$F_Occup <- factor(data$focc)
data$F_Educ <- factor(data$fedu)
data$M_Occup <- factor(data$mocc)
data$M_Educ <- factor(data$medu)

model <- income ~ Sex + Ethnicity + F_Occup + F_Educ + M_Occup + M_Educ

# Tune Trafotree ----

# Set seed for replicability
set.seed(1)

# Launch tuning functions
tune <- tune_trafotree(data = data, model = model, folds = 5, minorder = 2,
               maxorder = 10, mindiff = 0.01, plot = TRUE)

# Get tuning results.
order <- tune[["order"]]
loglik <- tune[["loglik"]]
results <- tune[["res"]]
plot <- tune[["plot"]]
plot(plot)

# Run Trafotree ----

# Launch trafotree, the order is the one previously selected
trafotree <- get_trtree(data = data, model = model, dep = "income",
                        order = order, mincri = 0.99, minbu = 100, 
                        centiles = 99, lenv = TRUE, share_lenv = 0.1)

# Get trafotree output
tr_data <- trafotree[["trafodata"]]
tree <- trafotree[["tree"]]
qtl <- trafotree[["qtl"]]
tr <- trafotree[["tr"]]
eop <- trafotree[["eop"]]
eopx <- trafotree[["eopx"]]

# Results Trafotree ----

iop <- tr_data %>%
  summarise(name = "TrTree",
            types = length(unique(types)),
            gini = gini.wtd(y_tilde, weights),
            mld = mld.wtd(y_tilde, weights),
            loglik = as.numeric(loglik),
            eop = eop,
            eopx = eopx)

# Plot Trafotree ----

# BEFORE PLOTTING THE TRAFOTREE, MAKE SURE THAT YOU HAVE A DATAFRAME NAMED "data"
# WHERE data$income is your outcome, and data$weights are your weights. 
# If you have no weights, then generate data$weights = 1. Otherwise the function
# does not work!

# data <- tr_data %>%
#  dplyr::select(income, weights)

plot <- plot(ATR::rotate(tree), 
             terminal_panel = node_dense(tree, font_t = 3),
             gp = gpar(fontsize = 8))

# Other Plots ----

# Get log interpolation for plot

tr_data$loginc <- log(tr_data$income)

logtr <- get_tr(data = tr_data, dep = "loginc", types = "types",
            order = 20, centiles = 99)

logdata <- logtr[["trafodata"]]
qtl_logdata <- logtr[["qtl"]]
tr_logdata <- logtr[["tr"]]

# Get colors for plots

col <- colplot(data = logdata, dep = "loginc", types = "types",
               grouping_var = "types")

data_plot <- col[["plot_data"]]
colors <- col[["col"]]
nti <- length(unique(data_plot$types))

# Get ECDF (ponytail) plot

ecdf <- plot_ponytail(data = data_plot, dep = "loginc", group = "types", fill = "color_types",
                      quart_var = qtl_logdata, 
                      tr = tr_logdata, nti = nti, 
                      limit_x = c(min(data_plot$loginc), max(data_plot$loginc)), 
                      col = colors, 
                      x_lab = "Log Income", y_lab = "ECDF", 
                      legend_pos = "bottom")

plot(ecdf)

# Get Density (mountain) plot

mountain <- plot_mountain(data = data_plot, dep = "income", 
                          col = colors,   limit_x = c(5, 12), 
                         group = "types", fill = "color_types", x_lab = "Log Income", 
                         y_lab = "Density", labs = "Types: ")

plot(mountain)

# Get lenv

plot <- plot_ponytail(data_plot, dep = "loginc", quart_var = qtl, 
                  tr = tr_logdata, nti = nti, 
                  limit_x = c(min(data_plot$loginc), max(data_plot$loginc)), 
                  group = "types", fill = "color_types", 
                  col = rep("grey", length(colors)), 
                  x_lab = "Log Income", y_lab = "ECDF", 
                  legend_pos = "bottom")

plot <- plot +
  stat_ecdf(data = logdata, aes(x=eopx_thres), col = "#E0112B", size = 1)

plot <- plot + theme(legend.position="")
plot <- plot + coord_flip(xlim = c(min(data_plot$loginc)-0.5, max(data_plot$loginc)+0.2)) 
plot(plot)
