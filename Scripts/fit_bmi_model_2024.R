library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk

# Set some variables we will use later in automation
univariable_analysis <- FALSE
diagnostics <- FALSE
plots <- FALSE

Sys.setenv("DISPLAY"=":0")

#ds <- fread("./dt.csv") # Load the dataset
#View(ds)

bmi_df <- read_fst("./inputs/bmi_ehis.fst", as.data.table = TRUE)

# use middle of age grp.
# Change col names. Make similar to project dataset.

head(bmi_df)

cols <- c("AGE", "YEAR", "SEX", "COUNTRY", "ALC", "SMK", "BMI")

bmi_df <- bmi_df[, ..cols]

head(bmi_df)


bmi_df[, `:=`(
  SEX = factor(SEX),
  SMK = factor(SMK),
  ALC = factor(ALC)
)]


# First identify what distribution is best fit for SBP (the left hand side
# variable)
plot(density(bmi_df$BMI))

summary(ds$sbp)
# This is a continuous distribution that cannot be 0

# Let's see what distribution work best for this. We will use the function
# fitDist from gamlss. Please have a read of the documentation for this function
marg_distr <- fitDist(sbp,
  k = 2, # Default is for AIC. use log(length(ds$sbp)) for BIC
  type = "realAll", # Other options "realline", "realplus", "real0to1", "counts", "binom"
  try.gamlss = TRUE, # Better but slow
  extra = NULL,
  data = ds,
  trace = TRUE
)

# Check the best 10 distr
head(marg_distr$fits, 10)

# lets plot some of them
m1 <- gamlss(sbp ~ 1,
  family = names(marg_distr$fits[1]),
  # weights = ds$wt_nurse,
  # method = mixed(20, 20)
  data = ds
)

# Let's create a helper function
sample_from_gamlss <- function(m, sample_size = 1e4) {
  # m is a gamlss model
  stopifnot(inherits(m, "gamlss")) # Stop if model is not gamlss class
  stopifnot(unique(sapply(
    lapply(
      m$parameters,
      function(x) {
        as.numeric(coef(m1, x))
      }
    ), length
  )) == 1L) # stop if predictors are present
  parm <- lapply(m$parameters, function(x) {
    fitted(m, x)[1]
  })
  names(parm) <- m$parameters
  y <-
    do.call(paste0("r", m$family[1]), c("n" = sample_size, parm)) # Sample from the model
}

plot(density(ds$sbp), lwd = 3, ylim = c(0, 0.03))
lines(density(sample_from_gamlss(m1)), col = alpha("red", 0.4))

m2 <- gamlss(sbp ~ 1, family = names(marg_distr$fits[4]), data = ds)
lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))

m3 <- gamlss(sbp ~ 1, family = names(marg_distr$fits[20]), data = ds)
lines(density(sample_from_gamlss(m3)), col = alpha("green", 0.4))

# An alternative plot
reldist(sample_from_gamlss(m1), ds$sbp, method = "bgk", bar = TRUE, show = "effect")

# What would the plot looked like if the distribution was off
reldist(sample_from_gamlss(m1) + 1, ds$sbp, method = "bgk", bar = TRUE, show = "effect")
reldist(sample_from_gamlss(m1) - 1, ds$sbp, method = "bgk", bar = TRUE, show = "effect")


distr_nam <- names(marg_distr$fits[4]) # I will pick GG
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3

if (univariable_analysis) {
  # Age
  ds[, .(sbp_median = wtd.quantile(sbp, weight = wt_nurse)), keyby = .(age)][, scatter.smooth(age, sbp_median)]

  m_age0 <- gamlss(
    sbp ~ age,
    family = distr_nam,
    # weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age0, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "purple")

  m_age1 <- gamlss(
    sbp ~ pb(age),
    family = distr_nam,
    # weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")

  m_age2 <- gamlss(
    sbp ~ poly(age, 3),
    family = distr_nam,
    # weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")

  m_age3 <- gamlss(
    sbp ~ cs(age),
    family = distr_nam,
    # weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")
  GAIC.table(m_age0, m_age1, m_age2, m_age3)
  centiles(m_age1, xvar = ds$age)
  centiles(m_age2, xvar = ds$age)

  ds[, .(sbp_median = wtd.quantile(sbp, weight = wt_nurse)), keyby = .(year)][, scatter.smooth(year, sbp_median, xlim = c(3, 40), ylim = c(110, 130))]

  # Year #TODO
  m_year1 <- gamlss(
    sbp ~ year,
    family = distr_nam,
    # weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "blue1")

  m_year2 <- gamlss(
    sbp ~ log(year),
    family = distr_nam,
    # weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "red1")

  m_year3 <- gamlss(
    sbp ~ log(year + 100),
    family = distr_nam,
    # weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")

  GAIC.table(m_year1, m_year2, m_year3)
  centiles(m_year1, xvar = ds$year)
  centiles(m_year2, xvar = ds$year)

  m_sm1 <- gamlss(
    sbp ~ smok_status,
    family = distr_nam,
    # weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  summary(m_sm1)
}

# Theory driven model selection
sbp_model <- gamlss(
  sbp ~ log(year + 100) + pb(age) + pcat(qimd) + sex + pcat(sha) + pcat(smok_status),
  ~ log(year + 100) + pb(age),
  ~ log(year + 100) + pb(age),
  family = distr_nam,
  # weights = ds$wt_nurse,
  data = ds,
  method = mixed(20, 400),
  control = gamlss.control(c.crit = 1e-2)
)
sbp_model <- update(sbp_model, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy
sbp_model$mydata <- copy(ds) # Optional and potential security breach

qsave(sbp_model, "./sbp_model.qs", preset = "high")


# Data driven model selection
mod_min <- gamlss(
  sbp ~ log(year + 100) + pb(age),
  family = distr_nam,
  # weights = ds_trn$wt_nurse,
  data = ds,
  method = mixed(5, 100),
  control = con1 # for speed but looses accuracy if not 1e-3
)

# Option 1
sbp_modelA <- stepGAICAll.A(
  mod_min,
  scope = list(
    lower = ~ log(year + 100) + pb(age),
    upper = ~ ( log(year + 100) + pb(age) + pcat(qimd) + sex + pcat(sha) + pcat(smok_status))^3
  ),
  sigma.scope = list(
    lower = ~1,
    upper = ~ log(year + 100) + pb(age) + pcat(qimd) + sex + pcat(sha) + pcat(smok_status)
  ),
  nu.scope = list(
    lower = ~1,
    upper = ~ log(year + 100) + pb(age) + pcat(qimd) + sex + pcat(sha) + pcat(smok_status)
  ),
  parallel = "multicore",
  ncpus = 12L
)
sbp_modelA
sbp_modelA <- update(sbp_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy
qsave(sbp_modelA, "./sbp_modelA.qs", preset = "high")


# Option 2
sbp_modelB <- stepGAICAll.B(
  mod_min,
  scope = list(
    lower = ~ log(year + 100) + pb(age),
    upper = ~ log(year + 100) + pb(age) + pcat(qimd) + sex + pcat(sha) + pcat(smok_status)
  ),
  parallel = "multicore",
  ncpus = 12L
)

sbp_modelB
sbp_modelB <- update(sbp_modelB, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy
qsave(sbp_modelB, "./sbp_modelB.qs", preset = "high")


GAIC.table(sbp_modelA, sbp_modelB, sbp_model)



tt <- chooseDist(sbp_modelA,
  type = "realplus",
  trace = TRUE, data = ds,
  parallel = "multicore", ncpus = 15L
) # Double check that the distribution is still a good one
qsave(tt, "./new_distr.qs", preset = "high")
# tt <- qread("./new_distr.qs")



which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC
which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC

sbp_modelbest <- update(sbp_modelA, family = "GB2") # Note that GB2 has 4 parameters and we only use predictors for the 1st 3

GAIC.table(sbp_modelbest, sbp_modelA, sbp_modelB, sbp_model)
qsave(sbp_modelbest, "./sbp_modelbest.qs", preset = "high")


# Code to create a table with predictions
trms <- all.vars(formula(sbp_modelbest))[-1] # -1 excludes dependent var
newdata <- CJ(
  year = 3:50, age = 20:90, sex = unique(ds$sex), qimd = unique(ds$qimd),
  sha = unique(ds$sha), smok_status = unique(ds$smok_status)
)
# newdata holds all the possible combinations of predictors

# This is to be able to parallelise
newdata <- split(newdata, by = "year")
newdata <- # assignment necessary! Copies of data.tables are happening
  lapply(
    newdata,
    function(x) {
      x[, (sbp_modelbest$parameters) := predictAll(sbp_modelbest, .SD, data = ds), .SDcols = trms]
    }
  )
newdata <- rbindlist(newdata) # bind the chunks back
View(head(newdata, 1e3))
 # copulus
write_fst(newdata, path = "./sbp_table.fst", compress = 100L) # This is what my models use as an input

print("Table saved")

if (diagnostics) {
  sbp_model <- qread("./sbp_modelbest.qs")

  plot(sbp_model)

  wp(sbp_model) # detrended QQ-plot
  wp(resid = resid(sbp_model)) # equivalen to the one above
  wp(resid = resid(sbp_model), ylim.all = 80 * sqrt(1 / length(resid(sbp_model))))
  tt <- wp(sbp_model, xvar = ds$age, n.inter = 6)
  wp(sbp_model, xvar = ds$year, n.inter = 10)
  wp(resid = resid(sbp_model), xvar = ~ ds$qimd)

  dtop(sbp_model, xvar = ds$age)

  rqres.plot(sbp_model)
  rqres.plot(sbp_model, type = "QQ")
}

if (plots) {
  sbp_model <- qread("./sbp_modelbest.qs")
  dir.create("./validation/synthpop_models", FALSE)

  source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())


  xlab_nam <- expression(bold(SBP ~ (mmHg)))
  sbp_model_tbl <- read_fst("./sbp_table.fst", as.data.table = TRUE)
  zz <-
    validate_gamlss_tbl(ds[age > 19], sbp_model_tbl, 50, "sbp", paste0("q", sbp_model$family[1]))[
      between(sbp, quantile(sbp, 0.01), quantile(sbp, 0.99))
    ]
  zz[, weight := wt_nurse / sum(wt_nurse), by = "type"]

  png(
    "./validation/synthpop_models/SMK_rel_dist.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", sbp],
    zz[type == "Modelled", sbp],
    zz[type == "Observed", weight],
    zz[type == "Modelled", weight],
    main = xlab_nam,
    100
  )
  dev.off()

  plot_synthpop_val(zz, sbp, "agegrp10", "wt_nurse", "SBP by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "year", "wt_nurse", "SBP by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "qimd", "wt_nurse", "SBP by QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "sha", "wt_nurse", "SBP by SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "smok_status", "wt_nurse", "SBP by smoking status", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, c("year", "agegrp10"), "wt_nurse", "SBP by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, c("year", "qimd"), "wt_nurse", "SBP by year and QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, c("year", "sha"), "wt_nurse", "SBP by year and SHA", xlab_nam, FALSE, FALSE)
}

