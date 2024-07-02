library(data.table)
library(ggplot2)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk

# Set some variables we will use later in automation
univariable_analysis <- FALSE
diagnostics <- FALSE
plots <- FALSE

ds <- read_fst("./inputs/bmi_IT.fst", as.data.table = TRUE)

ds[, `:=`(sex = factor(sex), smk = factor(smk))]

head(ds)

plot(density(ds$bmi))

summary(ds$bmi)

bmi_density <- ggplot(ds, aes(x = bmi)) +
    geom_density() +
    labs(x = "BMI", y = "Density", title = "Density plot of BMI - Italy") +
    theme_gray() + # Apply a minimal theme
    theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        axis.line = element_line(color = "black"),
    )
bmi_density


bmi_age <- ggplot(ds, aes(x = age, y = bmi)) +
    geom_point(color = "blue", alpha = 0.7, size = 3) +
    labs(x = "Age", y = "BMI", title = "Scatterplot of Age VS BMI - Italy") +
    theme_minimal() + # Apply a minimal theme
    theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 10),
        axis.line = element_line(color = "black"),
    )

print(bmi_age)

ggsave("./outputs/plots/IT/bmi_age.jpg", plot = bmi_age, width = 10, height = 5, dpi = 600)
ggsave("./outputs/plots/IT/bmi_density.jpg", plot = bmi_density, width = 10, height = 5, dpi = 600)


# Checking which distribution fits best based on AIC (k=2)
marg_distr <- fitDist(bmi, k = 2, type = "realAll", try.gamlss = TRUE, extra = NULL, data = ds, trace = TRUE)
qsave(marg_distr, "./outputs/models/marg_distr.qs", preset = "high")

marg_distr <- qread("./outputs/models/marg_distr.qs")

head(marg_distr$fits, 10)

# distribution -> SEP4
m1 <- gamlss(bmi ~ 1, family = names(marg_distr$fits[1]), data = ds)
qsave(m1, "./outputs/models/m1_model.qs", preset = "high")

m1 <- qread("./outputs/models/m1_model.qs")


# distribution -> JSUo
m2 <- gamlss(bmi ~ 1, family = names(marg_distr$fits[3]), data = ds)
qsave(m2, "./outputs/models/m2_model.qs", preset = "high")

m2 <- qread("./outputs/models/m1_model.qs")


# distribution -> BCPEo
m3 <- gamlss(bmi ~ 1, family = names(marg_distr$fits[5]), data = ds)
qsave(m3, "./outputs/models/m3_model.qs", preset = "high")

m3 <- qread("./outputs/models/m3_model.qs")


# distribution -> IGAMMA
m4 <- gamlss(bmi ~ 1, family = names(marg_distr$fits[6]), data = ds)
qsave(m4, "./outputs/models/m4_model.qs", preset = "high")

m4 <- qread("./outputs/models/m4_model.qs")


# distribution -> GG
m5 <- gamlss(bmi ~ 1, family = names(marg_distr$fits[7]), data = ds)
qsave(m5, "./outputs/models/m5_model.qs", preset = "high")

m5 <- qread("./outputs/models/m5_model.qs")



# Sample from distribution
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


plot(density(ds$bmi))
lines(density(sample_from_gamlss(m1)), col = alpha("red", 0.4))

bmi_observed_sampled <- ggplot() +
    geom_density(aes(x = ds$bmi, color = "observed"), alpha = 0.5) +
    # geom_density(aes(x = sample_from_gamlss(m1), color = "sampled_SPE4"), alpha = 0.5) +
    # geom_density(aes(x = sample_from_gamlss(m2), color = "sampled_JSUo"), alpha = 0.5) +
    # geom_density(aes(x = sample_from_gamlss(m3), color = "sampled_BCPEo"), alpha = 0.5) +
    # geom_density(aes(x = sample_from_gamlss(m4), color = "sampled_IGAMMA"), alpha = 0.5) +
    geom_density(aes(x = sample_from_gamlss(m5), color = "sampled_GG"), alpha = 0.1) +
    scale_color_manual(values = c(
        "observed" = "black",
        # "sampled_SPE4" = "green",
        # "sampled_JSUo" = "blue",
        # "sampled_BCPEo" = "purple",
        # "sampled_IGAMMA" = "yellow",
        "sampled_GG" = "red"
    )) +
    labs(
        title = paste0("Density plot of observed BMI and Sampled BMI (Distribution: ", names(marg_distr$fits[1]), ")"),
        x = "BMI",
        y = "Density"
    ) +
    theme_minimal()

bmi_observed_sampled

ggsave("./outputs/plots/IT/bmi_observed_sampled.jpg", plot = bmi_observed_sampled, width = 10, height = 5, dpi = 600)


distr_nam <- names(marg_distr$fits[7]) # I will pick GG
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3


# Age for predicting BMI
bmi_median <- ds[, .(bmi_median = wtd.quantile(bmi)), keyby = .(age)]
bmi_median


m_age0 <- gamlss(bmi ~ age, family = distr_nam, data = ds, method = mixed(20, 20))
m_age1 <- gamlss(bmi ~ pb(age), family = distr_nam, data = ds, method = mixed(20, 20))
m_age2 <- gamlss(bmi ~ poly(age, 3), family = distr_nam, data = ds, method = mixed(20, 20))
m_age3 <- gamlss(bmi ~ cs(age), family = distr_nam, data = ds, method = mixed(20, 20))




age0_pred <- centiles.pred(m_age0, xname = "age", xvalues = 20:85, cent = 50, data = ds)
colnames(age0_pred) <- c("age", "bmi_median")

age1_pred <- centiles.pred(m_age1, xname = "age", xvalues = 20:85, cent = 50, data = ds)
colnames(age1_pred) <- c("age", "bmi_median")

age2_pred <- centiles.pred(m_age2, xname = "age", xvalues = 20:85, cent = 50, data = ds)
colnames(age2_pred) <- c("age", "bmi_median")

age3_pred <- centiles.pred(m_age3, xname = "age", xvalues = 20:85, cent = 50, data = ds)
colnames(age3_pred) <- c("age", "bmi_median")


bmi_median_plot <- ggplot() +
    geom_point(data = bmi_median, aes(x = age, y = bmi_median)) +
    geom_line(data = age0_pred, aes(x = age, y = bmi_median, color = "age")) +
    geom_line(data = age1_pred, aes(x = age, y = bmi_median, color = "pb(age)")) +
    geom_line(data = age2_pred, aes(x = age, y = bmi_median, color = "poly(age)")) +
    geom_line(data = age3_pred, aes(x = age, y = bmi_median, color = "cs(age)")) +
    scale_color_manual(name = "Age Groups", values = c("age" = "black", "pb(age)" = "green", "poly(age)" = "steelblue", "cs(age)" = "red"))

bmi_median_plot

ggsave("./outputs/plots/IT/bmi_median_vs_age.jpg", plot = bmi_median_plot, width = 10, height = 5, dpi = 600)


GAIC.table(m_age0, m_age1, m_age2, m_age3)

centiles_age1_plot <- centiles(m_age1, xvar = ds$age)
ggsave("./outputs/plots/IT/centiles_age1_plot.jpg", plot = centiles_age1_plot, width = 10, height = 5, dpi = 600)



# Year for predicting BMI
bmi_median <- ds[, .(bmi_median = wtd.quantile(bmi)), keyby = .(year)]
bmi_median

m_year0 <- gamlss(bmi ~ year, family = distr_nam, data = ds, method = mixed(20, 20))
m_year1 <- gamlss(bmi ~ log(year), family = distr_nam, data = ds, method = mixed(20, 20))
m_year2 <- gamlss(bmi ~ log(year + 100), family = distr_nam, data = ds, method = mixed(20, 20))


year0_pred <- centiles.pred(m_year0, xname = "year", xvalues = 3:40, cent = 50, data = ds)
colnames(year0_pred) <- c("year", "bmi_median")

year1_pred <- centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds)
colnames(year1_pred) <- c("year", "bmi_median")

year2_pred <- centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds)
colnames(year2_pred) <- c("year", "bmi_median")


year_median_plot <- ggplot() +
    geom_point(data = bmi_median, aes(x = year, y = bmi_median)) +
    geom_line(data = year0_pred, aes(x = year, y = bmi_median, color = "years")) +
    geom_line(data = year1_pred, aes(x = year, y = bmi_median, color = "log(year)")) +
    geom_line(data = year2_pred, aes(x = year, y = bmi_median, color = "log(year + 100)")) +
    scale_color_manual(name = "Year", values = c("years" = "green", "log(year)" = "red", "log(year + 100)" = "steelblue"))


ggsave("./outputs/plots/IT/bmi_median_vs_year.jpg", plot = year_median_plot, width = 10, height = 5, dpi = 600)

GAIC.table(m_year0, m_year1, m_year2)

centiles_age1_plot <- centiles(m_year2, xvar = ds$year)

# ggsave("./outputs/plots/IT/centiles_age1_plot.jpg", plot = centiles_age1_plot, width = 10, height = 5, dpi = 600)

con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3

# Data driven model selection
mod_min <- gamlss(
    # bmi ~ log(year) + pb(age),
    bmi ~ 1,
    family = distr_nam,
    # weights = ds_trn$wt_nurse,
    data = ds,
    method = mixed(5, 100),
    control = con1 # for speed but looses accuracy if not 1e-3
)

qsave(mod_min, "./outputs/models/mod_min.qs", preset = "high")

# Option 1
bmi_modelA <- stepGAICgll.A(
    mod_min,
    scope = list(
        lower = ~1,
        upper = ~ (log(year) + pb(age) + sex + pcat(smk))^2
    ),
    sigma.scope = list(
        lower = ~1,
        upper = ~ log(year) + pb(age) + sex + pcat(smk)
    ),
    nu.scope = list(
        lower = ~1,
        # upper = ~ log(year) + pb(age)
        upper = ~ log(year) + pb(age) + sex + pcat(smk)
    ),
    parallel = "multicore",
    ncpus = 12L
)

bmi_modelA

bmi_modelA <- update(bmi_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy

qsave(sbp_modelA, "./outputs/models/bmi_modelA.qs", preset = "high")

bmi_modelA <- qread("./outputs/models/bmi_modelA.qs")


tt <- chooseDist(bmi_modelA,
    type = "realplus",
    trace = TRUE, data = ds,
    parallel = "multicore", ncpus = 15L
) # Double check that the distribution is still a good one

qsave(tt, "./outputs/models/new_distr.qs", preset = "high")



bmi_modelbest <- update(bmi_modelA, family = "GB2") # Note that GB2 has 4 parameters and we only use predictors for the 1st 3


GAIC.table(bmi_modelA, bmi_modelbest)

qsave(bmi_modelbest, "./outputs/models/bmi_modelbest.qs", preset = "high")


# Code to create a table with predictions
trms <- all.vars(formula(bmi_modelbest))[-1] # -1 excludes dependent var
newdata <- CJ(
  year = 3:50, age = 20:90, sex = unique(ds$sex), smk = unique(ds$smk)
)
# newdata holds all the possible combinations of predictors

# This is to be able to parallelise
newdata <- split(newdata, by = "year")
newdata <- # assignment necessary! Copies of data.tables are happening
  lapply(
    newdata,
    function(x) {
      x[, (bmi_modelbest$parameters) := predictAll(bmi_modelbest, .SD, data = ds), .SDcols = trms]
    }
  )

newdata <- rbindlist(newdata) # bind the chunks back
View(head(newdata, 1e3))

write_fst(newdata, path = "./inputs/bmi_IT_table.fst", compress = 100L) # This is what my models use as an input

print("Table saved")

diagnostics = TRUE


bmi_model <- qread("./outputs/models/bmi_modelbest.qs")

plot(bmi_model)

wp(bmi_model) # detrended QQ-plot

wp(resid = resid(bmi_model)) # equivalen to the one above

wp(resid = resid(bmi_model), ylim.all = 80 * sqrt(1 / length(resid(bmi_model))))

tt <- wp(bmi_model, xvar = ds$age, n.inter = 6)

wp(bmi_model, xvar = ds$year, n.inter = 10)


dtop(bmi_model, xvar = ds$age)

rqres.plot(bmi_model)

rqres.plot(bmi_model, type = "QQ")



bmi_model <- qread("./outputs/models/bmi_modelbest.qs")

dir.create("./validation/synthpop_models", FALSE)

source("misc_functions.R")
library(dqrng)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


xlab_nam <- expression(bold(SBP ~ (mmHg)))

bmi_model_tbl <- read_fst("./inputs/bmi_IT_table.fst", as.data.table = TRUE)

zz <-
  validate_gamlss_tbl(ds[age > 19], bmi_model_tbl, 50, "bmi", paste0("q", bmi_model$family[1]))[
    between(bmi, quantile(bmi, 0.01), quantile(bmi, 0.99))
  ]

zz[, weight := 1]

png(
  "./validation/synthpop_models/BMI_IT_rel_dist.png",
  3840,
  2160,
  pointsize = 48
)

reldist_diagnostics(zz[type == "Observed", bmi],
  zz[type == "Modelled", bmi],
  zz[type == "Observed", weight],
  zz[type == "Modelled", weight],
  main = xlab_nam,
  100
)
dev.off()

library(CKutils)

to_agegrp(zz, 10, agegrp_colname = "agegrp10")

plot_synthpop_val(zz, bmi, "agegrp10", "weight", "BMI by agegroup", xlab_nam, FALSE, TRUE)

plot_synthpop_val(zz, bmi, "year", "weight", "BMI by year", xlab_nam, FALSE, FALSE)

plot_synthpop_val(zz, bmi, "smk", "weight", "BMI by smoking status", xlab_nam, FALSE, FALSE)

plot_synthpop_val(zz, bmi, c("year", "agegrp10"), "weight", "BMI by year and agegroup", xlab_nam, FALSE, FALSE)


