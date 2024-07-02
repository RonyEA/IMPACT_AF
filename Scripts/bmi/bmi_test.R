library(fst)
library(igraph)
library(data.table)

source("/home/rony/projects/IMPACTaf/Rpackage/IMPACTaf_model_pkg/R/Design_class.R")


design <- Design$new("/home/rony/projects/IMPACTaf/japan_backup/inputs_jp/sim_design.yaml")
dt <- read_fst("./japan_backup/inputs_jp/pop_estimates/observed_population_japan.fst", as.data.table = TRUE)


# dt <- dt[age > 20& age < 30]
# dt <- dt[year == 2020]
# dt[, prbl := pops / sum(pops)][, `:=`(reg = NULL, pops = NULL)]
# dtinit <- dt[sample(.N, 10, TRUE, prbl)]
# dtinit 
# dt <- dt[age == 21]
# siz <- dtinit[age == 21, .N]
# dtfut <- dt[sample(.N, siz*5 , TRUE, prbl)]
# dtfut[, age := age - rep(1:5, siz)]
# dt <- rbind(dtfut, dtinit)
# dt[, prbl := NULL]
# dt


dt <- dt[age %in% c(design$sim_prm$ageL:design$sim_prm$ageH)] # TODO: Needs fix? Check old version
dt[, prbl := pops / sum(pops)][, `:=`(reg = NULL, pops = NULL)]

dtinit <- dt[sample(.N, design$sim_prm$n, TRUE, prbl)]
dtinit 

dt <- dt[age == design$sim_prm$ageL]

siz <- dtinit[age == design$sim_prm$ageL, .N]

dtfut <- dt[sample(.N, siz * design$sim_prm$sim_horizon_max, TRUE, prbl)]
# dtfut[, .N, by=.(age, sex, year)]
# dtinit[age==30, .N, by=.(age, sex, year)]
dtfut[, age := age - rep(1:design$sim_prm$sim_horizon_max, siz)]

dt <- rbind(dtfut, dtinit)
dt[, prbl := NULL]

dt[, `:=`(pid = .I)]