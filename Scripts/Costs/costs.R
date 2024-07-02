library(data.table)
age <- 65:90 # age > 65
sex <- c("men", "women")
death <- c(0,1)
clusters <- c(1, 2, 3, 6, 7)

dt <- CJ(age = age, sex = sex, MMCluster = clusters, death = death)
dt <- dt[, MMCluster := as.factor(MMCluster)]
dt <- dt[, death := as.factor(death)]
df <- data.frame(dt)

tt <- readRDS("/mnt/storage_slow/AFFIRMO/Padova_cost_lms/lm_total_cost_MMcluster.rda")

dt[, predicted_total_cost := exp(predict(tt, newdata = df))]

dt[, wghts_age := weighted.mean(predicted_total_cost, w = predicted_total_cost/sum(predicted_total_cost)), by=age]
dt[, wghts_mmcluster := weighted.mean(predicted_total_cost, w = predicted_total_cost/sum(predicted_total_cost)), by=MMCluster]
dt[, wghts_death := weighted.mean(predicted_total_cost, w = predicted_total_cost/sum(predicted_total_cost)), by=death]

dt
