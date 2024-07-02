## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncd_Engl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

setwd("/home/rony/projects/IMPACTaf/")
# setwd("/home/ckyprid/My_Models/IMPACTncd_Engl/")
# # For ages 20 to 90
# univariate_analysis <- FALSE
# diagnostics         <- FALSE
# plots               <- TRUE
seed                <- 30L


# if (!require(CKutils)) {
#   if (!require(remotes)) install.packages("remotes")
#   remotes::install_github("ChristK/CKutils")
# }
# dependencies(c("qs", "fst", "MASS", "splines", "matrixStats", "reldist", "future", "future.apply", "data.table"))
# options(future.fork.enable = TRUE)
# plan(multiprocess)

# if (file.exists("./preparatory_work/HSE_ts.fst")) {
#   HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
# } else {
#   source("./preparatory_work/preprocess_HSE.R", local = TRUE)
# }
# source("./preparatory_work/aux_functions.R", local = TRUE)

library(fst)
library(data.table)
library(MASS)
library(splines)
library(matrixStats)
library(future.apply)
library(qs)



dt <- read_fst("/mnt/rony/UoL/AFFIRMO/EHIS_microdata/inputs/bmi_IT.fst", as.data.table = TRUE)


dt[, smk := ordered(smk)]

set.seed(seed)

if (univariate_analysis) {
  # age
  #age_scaled <- scale(20:90, 49.6, 17)

  dt[, .(smk_mean = mean(as.integer(smk))), keyby = .(age)
     ][, scatter.smooth(age, smk_mean, ylim = c(0, 3))]


  #dt[, .(active_days_mean = wtd.mean(as.integer(active_days), weight = wt_int)), keyby = .(age)
  #   ][, scatter.smooth(age, active_days_mean, ylim = c(0, 7))]

  m_age1 <- polr(
    smk ~ ns(age, 2),
    #weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )

  tt <- predict(m_age1, type = "p", newdata = data.frame(age = 30:85))
  lines(30:85, apply(tt, 1, function(x) mean(sample(3, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- polr(
    smk ~ poly(age, 2),
    #weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )

  tt <- predict(m_age2, type = "p", newdata = data.frame(age = 30:85))
  lines(30:85, apply(tt, 1, function(x) mean(sample(3, 1e5, TRUE, x))), col = "red1")

  #tt <- predict(m_age2, type = "p", newdata = data.frame(age = age_scaled))
  #lines(age_scaled, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "red1")

  m_age3 <- polr(
    smk ~ bs(age, 3),
    #weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )

  tt <- predict(m_age3, type = "p", newdata = data.frame(age = 30:85))
  lines(30:85, apply(tt, 1, function(x) mean(sample(3, 1e5, TRUE, x))), col = "green1")

  m_age4 <- polr(
    smk ~ log(age),
    #weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )

  tt <- predict(m_age3, type = "p", newdata = data.frame(age = 30:85))
  lines(30:85, apply(tt, 1, function(x) mean(sample(3, 1e5, TRUE, x))), col = "green1")

  #tt <- predict(m_age3, type = "p", newdata = data.frame(age = age_scaled))
  #lines(age_scaled, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "green1")

  setDT(BIC(m_age1, m_age2, m_age3, m_age4), keep.rownames = TRUE, key = "BIC")[] # m_age3

  # year (This is likely meaningless as we project quintiles)

  dt[, .(smk_mean = mean(as.integer(smk))), keyby = .(year)
     ][, plot(year, smk_mean, xlim = c(14, 20), ylim = c(0, 3))]

  m_year1 <- polr(
    smk ~ year,
    #weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )

  tt <- predict(m_age2, type = "p", newdata = data.frame(age = 30:85))
  lines(30:85, apply(tt, 1, function(x) mean(sample(3, 1e5, TRUE, x))), col = "blue1")

  #tt <- predict(m_year1 , type = "p", newdata = data.frame(year = 14:20))
  #lines(14:20, apply(tt, 1, function(x) mean(sample(3, 1e5, TRUE, x))), col = "blue1")

  m_year2 <- polr(
    smk ~ year,
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year2 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "red1")

  setDT(BIC(m_year1, m_year2), keep.rownames = TRUE, key = "BIC")[] # equal I will accept the linear

}


smk_model <- polr(
  smk ~ year * poly(age, 2) *  sex,
  #weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = TRUE
) 




mod_min <- polr(
  smk ~ year + poly(age, 2) + sex,
  #weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = T
)

smk_model2 <- stepAIC(mod_min,
                             scope = list(
                              lower = ~ year + poly(age, 2) + sex,
                               upper = ~ (
                                 year + poly(age, 2) + sex
                               ) ^ 3
                             ),
                             direction = "both",
                             #k = log(nrow(dt))
                             k = 2
)


smk_model3 <- polr(
  smk ~ year + poly(age, 3) + sex + poly(age, 3):sex + 
    year:poly(age, 3) + year:sex + year:poly(age, 3):sex,
  #weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = T
)

smk_model4 <- polr(
  smk ~ year + ns(age, 3) + sex + ns(age, 3):sex + 
    year:ns(age, 3) + year:sex + year:ns(age, 3):sex,
  #weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = T
)


  setDT(BIC(smk_model2, smk_model3, smk_model4), keep.rownames = TRUE, key = "BIC")[] # equal I will accept the linear

smk_model2$data <- copy(dt)

qsave(smk_model2, "./models/lifecourse_models/smk_IT_model.qs", preset = "high")

smk_model2 <- qread("/models/lifecourse_models/smk_IT_model.qs")

print("Model saved")

trms <- all.vars(formula(smk_model2))[-1] # -1 excludes dependent var
trms

newdata <- CJ(year = 1:50, age = 30:90, sex = unique(dt$sex)) #
newdata

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("smk", 1:3)) := data.table(rowCumsums(predict(smk_model2, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))


newdata <- rbindlist(newdata)

# newdata[, age := age_int]
# newdata[, c("age_int", "pa7") := NULL]
# newdata[, "smk3" := NULL]

write_fst(newdata, "./inputs/exposure_distributions/smk_IT_table.fst", 100L)



if (plots) {
  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  xlab_nam <- expression(bold(Smoking ~ Status))

  tbl <- read_fst("./inputs/exposure_distributions/smk_IT_table.fst", as.data.table = TRUE)
  smk_model <- qread("./models/lifecourse_models/smk_IT_model.qs")

  zz <- clone_dt(smk_model$data, 10)
  zz[, smk := NULL]
  # zz[, age := round(age*17+49.6)]
  zz[, rank_smk := runif(.N)]

  nam <- intersect(names(zz), names(tbl))

  zz[tbl, smk := (rank_smk > smk1) + (rank_smk > smk2) + 1L, on = nam]

  zz[, `:=` (
    type = "Modelled",
    smk = factor(
      smk,
      levels = 1:3,
      labels = 1:3,
      ordered = TRUE
    ),
    .id = NULL
  )]
  zz[, rank_smk := NULL]
  zz <- rbind(zz, smk_model$data[, type := "Observed"])
  zz[, smk := as.integer(as.character(smk))]

  zz <- zz[!is.na(smk)]


  future({
    dir.create("./validation/synthpop_models", FALSE)
    zz[, weight := 1, by = type]
    png(
      "./validation/synthpop_models/smk_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", smk],
                        zz[type == "Modelled", smk],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        discrete = TRUE)
    dev.off()
  })

  future(plot_synthpop_val(zz, active_days, "agegrp10", "wt_int", "Active days by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "year", "wt_int", "Active days by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "sex", "wt_int", "Active days by sex", xlab_nam, FALSE, FALSE))

  future(plot_synthpop_val(zz, active_days, "qimd", "wt_int", "Active days by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "sha", "wt_int", "Active days by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "ethnicity", "wt_int", "Active days by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, c("year", "agegrp10"), "wt_int", "Active days by year and agegroup", xlab_nam, FALSE, FALSE))
}



