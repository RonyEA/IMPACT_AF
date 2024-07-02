library(fst)
library(data.table)

source("./global.R")

source("./Rpackage/IMPACTaf_model_pkg/R/SynthPop_class.R")
source("./Rpackage/IMPACTaf_model_pkg/R/Design_class.R")


# load dt
file <- "../Data/Final_Pop_Data/pop_proj_eu.fst"

dt <- read_fst(file, as.data.table = TRUE)
dt_meta <- metadata_fst(file)
stopifnot(
  "Population size file need to be keyed by year" =
    identical("year", dt_meta$keys[1])
)

countries <- design$sim_prm$country

# countrywise synth pop.
if (("country" %in% colnames(dt)) & (length(unique(dt[, country])) > 1)) {

    dt <- dt[country %in% countries]

    dt <- dt[, .(pops = sum(pops)), by = .(country, year, age, sex)]

}

dt <- dt[year == 2000L + design$sim_prm$init_year]

dt <- dt[age == "100+", age:=100]
dt <- dt[, age := as.numeric(age)]

dt <- dt[age %in% c(design$sim_prm$ageL:design$sim_prm$ageH)] # TODO: Needs fix? Check old version
dt <- dt[, prbl := pops / sum(pops)][, `:=`(pops = NULL)] 



# I do not explicitly set.seed because I do so in the gen_synthpop()
dtinit <- dt[sample(.N, design$sim_prm$n, TRUE, prbl)] 
setkey(dtinit, year, age)


if (design$sim_prm$logs) {
  message("Generate the cohorts of ", design$sim_prm$ageL, " year old")
}

dt <- dt[age == design$sim_prm$ageL]

siz <- dtinit[age == design$sim_prm$ageL, .N]

dtfut <- dt[sample(.N, siz * design$sim_prm$sim_horizon_max, TRUE, prbl)]
dtfut[, age := age - rep(1:design$sim_prm$sim_horizon_max, siz)]


dt <- rbind(dtfut, dtinit)
dt[, prbl := NULL]
      


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
  "rankstat_smk"
) := rank_mtx[, list(
  bmi_r,
  smk_r
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

limit_age <- data.table(min = min(dt$age), max = max(dt$age))


# Generate smk ----
# Prediction for smoking.
dt[, smk := as.factor(1L)] # never smoker
dt[, rankstat_smk := NULL]

# Generate bmi ----
# countrywise bmi_dt has to be loaded : RONY
if (design_$sim_prm$logs) message("Generate BMI")

# the path comes from design_yml
tbl <-
  read_fst("./inputs/exposure_distributions/bmi_IT_table.fst",
    as.data.table = TRUE
  )

tbl[, year := year+2000]

tbl[sex == "men", sex := "Male"]
tbl[sex == "women", sex := "Female"]

dt <- dt[year >= 2003]

print(str(dt))

print(str(tbl))

col_nam <-
  setdiff(names(tbl), intersect(names(dt), names(tbl)))
# if (.Platform$OS.type == "unix") {
#  lookup_dt(dt, tbl, check_lookup_tbl_validity = FALSE) #TODO: Lookup_dt
# } else {
absorb_dt(dt, tbl)
# }

dt[, smk := as.integer(smk)]


print(str(dt))

dt[, bmi := qGB2(rank_bmi, mu, sigma, nu, tau), ] # , n_cpu = design_$sim_prm$n_cpu)]
dt[bmi < 10, bmi := 10] # Truncate BMI predictions to avoid unrealistic values.
dt[bmi > 70, bmi := 70] # Truncate BMI predictions to avoid unrealistic values.

if (!design_$sim_prm$keep_simulants_rn) col_nam <- c(col_nam, "rank_bmi")
dt[, c(col_nam) := NULL]

## --------------------------------------------------
dt[, `:=`(
  pid_mrk = NULL
  # to be recreated when loading synthpop
)]

# ????? 20230206  # all exposure names  we do not need rank_ rankstat_
xps_tolag <- c(
      "bmi",
      "smk"
)

xps_nam <- paste0(xps_tolag, "_curr_xps")
setnames(dt, xps_tolag, xps_nam)

if ("age100" %in% names(dt)) {
  dt[, age := NULL]
  setnames(dt, "age100", "age")
}

dt[, sex := factor(sex)]
dt[, year := as.integer(year)]
#dt[, smk_curr_xps := as.integer(smk_curr_xps)]

setkey(dt, pid, year) # Just in case
setcolorder(dt, c("pid", "year", "age", "sex")) # STRATA
setindexv(dt, c("year", "age", "sex")) # STRATA

