
#         Author: Pedro Salas Rojo and Paolo Brunori
#         Date: 10/2023
#         Name of project:Apply Shapley-value decomposition to IOp

rm(list = ls(all.names = TRUE)) 
library(dplyr)

source("C:/Users/user/Dropbox/BRUNORI&SALAS_/database_github_etc/TO SHARE/Expost_IOp/Expost_IOp/functions_exante.R")
source("C:/Users/user/Dropbox/BRUNORI&SALAS_/database_github_etc/TO SHARE/Expost_IOp/Expost_IOp/functions_expost.R")
source("C:/Users/user/Dropbox/BRUNORI&SALAS_/database_github_etc/TO SHARE/Expost_IOp/Expost_IOp/functions_shapley.R")

# Get data, define model, define circumstances

data <- read.csv("C:/Users/user/Dropbox/BRUNORI&SALAS_/database_github_etc/TO SHARE/Expost_IOp/Expost_IOp/data.csv")

circum <- c("eth", "focc", "fedu", "sex")

mod <- income ~ factor(eth) + factor(focc) + factor(fedu) + factor(sex)

# Examples

# Ctree Shapley

 ctreeshap <- shapley(model = mod, data = data, depname = "income", wts = "weights", 
                      vars = circum, type="ctree", ntree = 5, mincri = 0.5, minbu = 100,
                      resample = 0.7)
 
 ctreeshap[["rel_shapval"]]
 ctreeshap[["all_results"]]
 
 # Ctree Trafotree Shapley
 
 trafoshap <- shapley(model = mod, data = data, depname = "income", wts = "weights", 
                      vars = circum, type="trafotree", ntree = 2, mincri = 0.5, minbu = 100,
                      resample = 0.7, centiles = 99, order = 4)
 
 trafoshap[["rel_shapval"]]

 # EOp Shapley
 
eopshap <- shapley_eop(model = mod, data = data, depname = "income", wts = "weights", 
                       vars = circum, ntree = 2, mincri = 0.5, minbu = 100,
                       resample = 0.7, centiles = 99, order = 4, share_lenv = 0.1, lenv = TRUE)

eopshap[["rel_shap_max_eop"]]
eopshap[["rel_shap_max_eopx"]]


