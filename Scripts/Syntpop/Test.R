library(IMPACTaf)
library(igraph)
library(fst)
library(data.table)
library(dqrng)
library(CKutils)
library(gamlss)

source("./global.R")

source("./Rpackage/IMPACTaf_model_pkg/R/Design_class.R")

design_ <- Design$new("./inputs/sim_design.yaml")


# file <- "../Data/Final_Pop_Data/pop_proj_eu.fst"
file <- "./inputs/pop_observed/pop_observed_eu.fst"

dt <- read_fst(file, as.data.table = TRUE)
dt_meta <- metadata_fst(file)
stopifnot(
  "Population size file need to be keyed by year" =
    identical("year", dt_meta$keys[1])
)

countries <- design_$sim_prm$country

# countrywise synth pop.
if (("country" %in% colnames(dt)) & (length(unique(dt[, country])) > 1)) {

    dt <- dt[country %in% countries]

    dt <- dt[, .(pops = sum(pops)), by = .(year, age, sex)]

}

dt <- dt[year == 2000L + design_$sim_prm$init_year]
dt <- dt[age %in% c(design_$sim_prm$ageL:design_$sim_prm$ageH)] # TODO: Needs fix? Check old version

dt <- dt[, prbl := pops / sum(pops)][, `:=`(pops = NULL)] 


#dt <- dt[, prbl := pops / sum(pops)][, `:=`(reg = NULL, pops = NULL)] 

# I do not explicitly set.seed because I do so in the gen_synthpop()
dtinit <- dt[sample(.N, design_$sim_prm$n, TRUE, prbl)] 
setkey(dtinit, year, age)

if (design_$sim_prm$logs) {
  message("Generate the cohorts of ", design_$sim_prm$ageL, " year old")
}

dt <- dt[age == design_$sim_prm$ageL]
siz <- dtinit[age == design_$sim_prm$ageL, .N]

dtfut <- dt[sample(.N, siz * design_$sim_prm$sim_horizon_max, TRUE, prbl)]
dtfut[, age := age - rep(1:design_$sim_prm$sim_horizon_max, siz)]

dt <- rbind(dtfut, dtinit)
dt[, prbl := NULL]

dt[, `:=`(pid = .I)]
new_n <- nrow(dt)


# NOTE!! from now on year in the short form i.e. 13 not 2013
dt[, `:=`(pid = .I)]
new_n <- nrow(dt)


# Generate correlated ranks for the individuals ----
if (design_$sim_prm$logs) {
  message("Generate correlated ranks for the individuals")
}

cm_mean <- as.matrix(
  read_fst(
    "./inputs/exposure_distributions/exposure_corr_mean.fst", # Change-for-IMPACT-NCD-Japan
    as.data.table = TRUE
  ),
  rownames = "rn"
)
          
diag(cm_mean) <- 1

# ??generate_corr_unifs(new_n, cm_mean)
# Change-for-IMPACT-NCD-Japan
rank_mtx <- generate_corr_unifs(new_n, cm_mean)
if (design_$sim_prm$logs) message("generate correlated uniforms")

rank_mtx <- data.table(rank_mtx)

# sum((cor(rank_mtx) - cm_mean) ^ 2)
if (design_$sim_prm$logs) message("correlated ranks matrix to data.table")

# Adding Exposures : RONY
dt[, c(
  "rank_bmi",
  "rankstat_smk",
  "rankstat_mm_cluster"
) := rank_mtx[, list(
  bmi_r,
  smk_r,
  mm_cluster_r
), ]]

rm(rank_mtx)

# Project forward for simulation and back project for lags  ----
if (design_$sim_prm$logs) message("Project forward and back project")

dt <-
  clone_dt(
    dt,
    design_$sim_prm$sim_horizon_max +
      design_$sim_prm$maxlag + 1L
  )

dt[.id <= design_$sim_prm$maxlag, `:=`(
  age = age - .id,
  year = year - .id
)]
dt[.id > design_$sim_prm$maxlag, `:=`(
  age  = age + .id - design_$sim_prm$maxlag - 1L,
  year = year + .id - design_$sim_prm$maxlag - 1L
)]
# dt <-
#   dt[between(age, design_$sim_prm$ageL - design_$sim_prm$maxlag, design_$sim_prm$ageH)]
# delete unnecessary ages
del_dt_rows(
  dt,
  !between(
    dt$age,
    design_$sim_prm$ageL - design_$sim_prm$maxlag,
    design_$sim_prm$ageH
  ),
  environment()
)

dt[, `:=`(.id = NULL)]

# Change-for-IMPACT-NCD-Japan
if (max(dt$age) >= 85L) {
  dt[, age100 := age]
  dt[age >= 85L, age := 85L]
}

# to_agegrp(dt, 20L, 85L, "age", "agegrp20", to_factor = TRUE)
# to_agegrp(dt, 10L, 85L, "age", "agegrp10", to_factor = TRUE)
# to_agegrp(dt,  5L, 85L, "age", "agegrp5" , to_factor = TRUE)

# Simulate exposures -----

# Random walk for ranks ----
if (design_$sim_prm$logs) message("Random walk for ranks")

setkeyv(dt, c("pid", "year"))
setindexv(dt, c("year", "age", "sex")) # STRATA

dt[, pid_mrk := mk_new_simulant_markers(pid)]

dt[, lapply(
  .SD,
  fscramble_trajectories,
  pid_mrk,
  design_$sim_prm$jumpiness
),
.SDcols = patterns("^rank_")
]
# ggplot2::qplot(year, rank_ssb, data = dt[pid %in% sample(1e1, 1)], ylim = c(0,1))

# Change-for-IMPACT-NCD-Japan
# Set limit age ranges
# Temp <- read_fst("/home/rstudio/IMPACT_NCD_data/NHNS_data/Output_data_organized/GAMLSS_created/HSE_ts.fst", as.data.table = TRUE)[between(Age, 20L, max(dt$age))]
# limit_age <- Temp[, .(min = min(Age), max = max(Age))]
# rm(Temp)
limit_age <- data.table(min = min(dt$age), max = max(dt$age))

dt <- dt[year >= 2003]
dt[, sex:=factor(sex)]

# Generate smk ----
# Prediction for smoking.

if (design_$sim_prm$logs) message("Generaet Smoking Status")

tbl <- read_fst("./inputs/exposure_distributions/smk_IT_table.fst", as.data.table = TRUE)

tbl[, year:= year+2000L]
tbl <- tbl[year >= 2003]

col_nam <- setdiff(names(tbl), intersect(names(dt), names(tbl)))
col_nam

absorb_dt(dt, tbl)

dt[, smk := factor(
    (rankstat_smk > smk1) + 
    (rankstat_smk > smk2) +
    (rankstat_smk > smk3) +
    1,
    levels = 1:3, labels = 1:3, ordered = TRUE
    )]
# dt[, rankstat_smk := NULL]

if (!design_$sim_prm$keep_simulants_rn) col_nam <- c(col_nam, "rankstat_smk")

dt[, c(col_nam) := NULL]

# Generate mm_clusters ----
tbl <- 
  read_fst("./inputs/exposure_distributions/af_mm_cluster_table.fst",
  as.data.table = TRUE
  )

col_nam <-
  setdiff(names(tbl), intersect(names(dt), names(tbl)))

absorb_dt(dt, tbl)
# print(str(dt))

dt[, mm_cluster := factor(
  (rankstat_mm_cluster > no_morbidity) + 
  (rankstat_mm_cluster > unspecific) +
  (rankstat_mm_cluster > neuropsychiatric) +
  (rankstat_mm_cluster > complex) +
  (rankstat_mm_cluster > eye) +
  (rankstat_mm_cluster > musculoskeletal) +
  (rankstat_mm_cluster > metabolic) +
  (rankstat_mm_cluster > cardiovascluar)
)]

if (!design_$sim_prm$keep_simulants_rn) col_nam <- c(col_nam, "rankstat_mm_cluster")
dt[, c(col_nam) := NULL]


# Generate bmi ----
# countrywise bmi_dt has to be loaded : RONY
if (design_$sim_prm$logs) message("Generate BMI")

# the path comes from design_yml
tbl <-
  read_fst("./inputs/exposure_distributions/bmi_IT_table.fst",
    as.data.table = TRUE
  )


tbl[, year := year+2000]
tbl <- tbl[year >= 2003]


print(str(dt))
print(str(tbl))

col_nam <-
  setdiff(names(tbl), intersect(names(dt), names(tbl)))

absorb_dt(dt, tbl)

dt[, smk := as.integer(smk)]
          

print(str(dt))

dt[, bmi := qGB2(rank_bmi, mu, sigma, nu, tau), ] # , n_cpu = design_$sim_prm$n_cpu)]
dt[bmi < 10, bmi := 10] # Truncate BMI predictions to avoid unrealistic values.
dt[bmi > 70, bmi := 70] # Truncate BMI predictions to avoid unrealistic values.

if (!design_$sim_prm$keep_simulants_rn) col_nam <- c(col_nam, "rank_bmi")
dt[, c(col_nam) := NULL]

dt_trial <- copy(dt)

tt <-
    read_fst("./inputs/pop_projections/pop_combined_eu.fst", as.data.table = TRUE) # Change-for-IMPACT-NCD-Japan

tt <- tt[country == "Italy"]

tt[age == "100+", age := 100]
tt[, age := as.integer(age)]

tt_trial <- copy(tt)

tt_trial <- tt_trial[
    between(age, min(dt_trial$age), max(dt_trial$age)) &
        between(year, min(dt_trial$year), max(dt_trial$year)),
    .(pops = sum(pops)),
    keyby = .(year, age, sex)
]


dt_trial[, wt_immrtl := .N, by = .(year, age, sex)]

absorb_dt(dt_trial, tt_trial)

dt_trial[, wt_immrtl := pops / (wt_immrtl * design_$sim_prm$n_synthpop_aggregation)]

dt_trial[, pops := NULL]



dt[, `:=` (
  age_group  = fcase(
    age >= 35 & age <= 44, "35-44",
    age >= 45 & age <= 54, "45-54",
    age >= 55 & age <= 64, "55-64",
    age >= 65 & age <= 74, "65-74",
    age >= 75, "75+"
  )
)]


strata <- c("mc", design_$sim_prm$strata_for_output)

dt[, c(lapply(.SD, function(x, wt) sum(x*wt), wt)),.SDcols = "utility_score", keyby = strata]
