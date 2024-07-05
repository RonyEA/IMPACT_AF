source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design_calibration.yaml")

# IMPACTncd$
#   del_logs()$
#   del_outputs()


self <- IMPACTncd$.__enclos_env__$self
private <- IMPACTncd$.__enclos_env__$private

replace <- FALSE
mc <- 1:10

self$reconstruct_large_files()

country <- self$design$sim_prm$country

export_xps <- self$design$sim_prm$export_xps # save the original value to be restored later

self$design$sim_prm$export_xps <- FALSE # turn off export_xps to speed up the calibration


private$create_empty_calibration_prms_file(replace = replace)


clbr <- fread("./simulation/calibration_prms.csv", 
                colClasses = list(numeric = c("af_incd_clbr_fctr",
                                               "nonmodelled_ftlt_clbr_fctr")))

memedian <- function(x) {
  out_med <- median(x)
  out_mean <- mean(x)

  if (out_med == 0 & out_mean == 0) {
    out <- 1
  } else if (out_med == 0) {
     out <- out_mean
  } else {
    out <- out_med
  }

  return(out)
}

if (replace) {
  age_start <- self$design$sim_prm$ageL
} else { # if replace == FALSE
  # if all ages exist skip calibration
  if (dim(clbr[af_incd_clbr_fctr == 1 | nonmodelled_ftlt_clbr_fctr == 1])[1] == 0) {
    message("All ages have been calibrated. Skipping calibration.")
    return(invisible(self))
  }
  age_start <- clbr[af_incd_clbr_fctr == 1 | nonmodelled_ftlt_clbr_fctr == 1, min(age)] # Unsafe but rarely
  message(paste0("Starting calibration from age ", age_start, "."))
}

age_ <- age_start + 0

# Run the simulation from min to max age
for (age_ in age_start:self$design$sim_prm$ageH) {
  # Run the simulation and export summaries. TODO restrict ages for efficiency.
  self$
    del_logs()$
    del_outputs()$
    run(mc, multicore = TRUE, "sc0")$
    export_summaries(multicore = TRUE, type = c("incd", "prvl", "dis_mrtl"), single_year_of_age = TRUE) # 
  # Incidence calibration
  # load the uncalibrated results
  unclbr <- fread(file.path(self$design$sim_prm$output_dir, "summaries", "incd_scaled_up.csv.gz"), 
                  select = c("year", "age", "sex", "mc", "popsize", "af_incd"))


########### Test #########################
  trial <- unclbr[age == age_, .(af_incd = as.integer(af_incd/popsize)), keyby = .(age, sex, year, mc)]

  trial[, .(af_incd = mean(af_incd)), keyby = .(age, sex, year)]
########### Test #########################


  unclbr <- unclbr[age == age_, .(af_incd = af_incd/popsize), keyby = .(age, sex, year, mc)
    ][, .(af_incd = memedian(af_incd)), keyby = .(age, sex, year)]

  # for CHD
  # fit a log-log linear model to the uncalibrated results and store the coefficients
  unclbr[af_incd > 0, c("intercept_unclbr", "trend_unclbr") := as.list(coef(lm(log(af_incd)~log(year)))), by = sex]
  unclbr[, intercept_unclbr := nafill(intercept_unclbr, "const", max(intercept_unclbr, na.rm = TRUE)), by = sex] # NOTE I use max just to return a value. It doesn't matter what value it is.
  unclbr[, trend_unclbr := nafill(trend_unclbr, "const", max(trend_unclbr, na.rm = TRUE)), by = sex] # NOTE I use max just to return a value. It doesn't matter what value it is.
  # load benchmark 
  benchmark <- read_fst(file.path("./inputs/disease_burden", country, "af_incd.fst"), columns = c("age", "sex", "year", "mu") , as.data.table = TRUE)[age == age_,] 
  # fit a log-log linear model to the benchmark incidence and store the coefficients 
  benchmark[year >= self$design$sim_prm$init_year_long, c("intercept_bnchmrk", "trend_bnchmrk") := as.list(coef(lm(log(mu)~log(year)))), by = sex]
  # calculate the calibration factors that the uncalibrated log-log model
  # need to be multiplied with so it can match the benchmark log-log model
  unclbr[benchmark[year == max(year)], af_incd_clbr_fctr := exp(intercept_bnchmrk + trend_bnchmrk * log(year)) / exp(intercept_unclbr + trend_unclbr * log(year)), on = c("age", "sex")] # Do not join on year!
  unclbr[, c("intercept_unclbr", "trend_unclbr") := NULL]
  # keep only year, age, sex, and calibration factors (to be multiplied
  # with p0)
  unclbr[, `:=` (af_prvl_correction = af_incd * (af_incd_clbr_fctr - 1),
                af_incd = NULL, stroke_incd = NULL, intercept_unclbr = NULL, trend_unclbr = NULL)]
  clbr[unclbr, on = c("year", "age", "sex"), `:=` (
    af_incd_clbr_fctr = i.af_incd_clbr_fctr
  )]

  # Case fatality calibration
  # Because we do incd and case fatality correction in the same step, we
  # need to estimate the expected changes on prvl because of the incd
  # calibration, before we proceed with the case fatality calibration.
  # Note that the calibration factor (multiplier) is 1/prvl as we
  # currently have mortality rates in the ftlt files.
  prvl <- fread(file.path(self$design$sim_prm$output_dir, "summaries", "prvl_scaled_up.csv.gz"), 
                  select = c("year", "age", "sex", "mc", "popsize", "af_prvl"))[age == age_,]
  # prvl <- prvl[, `:=` (
  #   chd_ftlt_clbr_fctr = (chd_prvl - chd_prvl*((stroke_mrtl + nonmodelled_mrtl)/popsize) + chd_prvl_correction * popsize)/chd_mrtl,
  #   stroke_ftlt_clbr_fctr = (stroke_prvl - stroke_prvl*((chd_mrtl + nonmodelled_mrtl)/popsize) + stroke_prvl_correction * popsize)/stroke_mrtl,
  #   nonmodelled_ftlt_clbr_fctr = popsize/(popsize - chd_mrtl - stroke_mrtl)
  #   )][, .(chd_ftlt_clbr_fctr = mean(chd_ftlt_clbr_fctr),
  #                 stroke_ftlt_clbr_fctr = mean(stroke_ftlt_clbr_fctr), 
  #                 nonmodelled_ftlt_clbr_fctr = mean(nonmodelled_ftlt_clbr_fctr)),
  #                  keyby = .(age, sex, year)]

  prvl <- prvl[, .(
    af_prvl = af_prvl/popsize,
    #stroke_prvl = stroke_prvl/popsize,
    popsize, age, sex, year, mc)
    ][, .(af_prvl = mean(af_prvl),
          popsize = mean(popsize)), keyby = .(age, sex, year)]
  prvl[unclbr, on = c("year", "age", "sex"), `:=` (
    af_prvl_correction = i.af_prvl_correction # Note corrections for prvl are rates
  )]
  #benchmark <- read_fst(file.path("./inputs/disease_burden", "chd_ftlt.fst"), columns = c("age", "sex", "year", "mu2") , as.data.table = TRUE)[age == age_,]
  #prvl[benchmark, on = c("age", "sex", "year"), chd_mrtl := mu2]
  #benchmark <- read_fst(file.path("./inputs/disease_burden", "stroke_ftlt.fst"), columns = c("age", "sex", "year", "mu2") , as.data.table = TRUE)[age == age_,]
  #prvl[benchmark, on = c("age", "sex", "year"), stroke_mrtl := mu2]
  benchmark <- read_fst(file.path("./inputs/disease_burden", country,  "nonmodelled_ftlt.fst"), columns = c("age", "sex", "year", "mu2") , as.data.table = TRUE)[age == age_,]
  prvl[benchmark, on = c("age", "sex", "year"), nonmodelled_mrtl := mu2]
  prvl[, `:=` (
    nonmodelled_ftlt_clbr_fctr = 1)]
  #prvl[, `:=` (
  #  af_ftlt_clbr_fctr = 1 / (chd_prvl + chd_prvl_correction), #  - stroke_mrtl - nonmodelled_mrtl
  #  stroke_ftlt_clbr_fctr = 1 / (stroke_prvl + stroke_prvl_correction), #  - chd_mrtl - nonmodelled_mrtl
  #  nonmodelled_ftlt_clbr_fctr = 1/(1 - chd_mrtl - stroke_mrtl))]
  # Fix the calibration factors for the ages that have been calibrated
  if (age_ > age_start) {
  # NOTE here age is age_1L
  mrtl <- fread(file.path(self$design$sim_prm$output_dir, "summaries", "dis_mrtl_scaled_up.csv.gz"), 
                  select = c("year", "age", "sex", "mc", "popsize", "af_deaths", "nonmodelled_deaths"))[age == age_ - 1L,]
  mrtl <- mrtl[, .(
    #af_mrtl = af_deaths/popsize,
    #stroke_mrtl = stroke_deaths/popsize,
    nonmodelled_mrtl = nonmodelled_deaths/popsize,
    popsize, age, sex, year, mc)
    ][, .(nonmodelled_mrtl = mean(nonmodelled_mrtl),
          popsize = mean(popsize)), keyby = .(age, sex, year)]
  #benchmark <- read_fst(file.path("./inputs/disease_burden", "chd_ftlt.fst"), columns = c("age", "sex", "year", "mu2") , as.data.table = TRUE)[age == age_ - 1L,]
  #mrtl[benchmark, on = c("age", "sex", "year"), chd_ftlt_clbr_fctr := mu2/chd_mrtl]
  #benchmark <- read_fst(file.path("./inputs/disease_burden", "stroke_ftlt.fst"), columns = c("age", "sex", "year", "mu2") , as.data.table = TRUE)[age == age_ - 1L,]
  #mrtl[benchmark, on = c("age", "sex", "year"), stroke_ftlt_clbr_fctr := mu2/stroke_mrtl]
  benchmark <- read_fst(file.path("./inputs/disease_burden", country, "nonmodelled_ftlt.fst"), columns = c("age", "sex", "year", "mu2") , as.data.table = TRUE)[age == age_ - 1L,]
  mrtl[benchmark, on = c("age", "sex", "year"), nonmodelled_ftlt_clbr_fctr := mu2/nonmodelled_mrtl]
  #mrtl[chd_ftlt_clbr_fctr == Inf, chd_ftlt_clbr_fctr := 1] # to avoid Inf through division by 0
  #mrtl[stroke_ftlt_clbr_fctr == Inf, stroke_ftlt_clbr_fctr := 1] # to avoid Inf through division by 0
  mrtl[nonmodelled_ftlt_clbr_fctr == Inf, nonmodelled_ftlt_clbr_fctr := 1] # to avoid Inf through division by 0
  clbr[mrtl, on = c("year", "age", "sex"), `:=` (
    #chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr * chd_ftlt_clbr_fctr,
    #stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr * stroke_ftlt_clbr_fctr,
    nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr * nonmodelled_ftlt_clbr_fctr
  )] 
  }
  if (age_ == self$design$sim_prm$ageH) {
    # shortcut for age == 99 hopefully with tiny bias
    mrtl[, age := age + 1L]
    prvl[mrtl, on = c("year", "age", "sex"), `:=` (
    #chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr * chd_ftlt_clbr_fctr,
    #stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr * stroke_ftlt_clbr_fctr,
    nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr * nonmodelled_ftlt_clbr_fctr
  )] 
  }
  clbr[prvl, on = c("year", "age", "sex"), `:=` (
    #chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr,
    #stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr,
    nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr
  )] 
  fwrite(clbr, "./simulation/calibration_prms.csv") # NOTE this needs to be inside the loop so it influences the simulation during the loop over ages
} # end loop over ages

self$design$sim_prm$export_xps <- export_xps # restore the original value
invisible(self)


















































sp <- SynthPop$new(mc_, self$design)

lapply(self$diseases, function(x) {
  if (self$design$sim_prm$logs) print(x$name)
  x$
   gen_parf(sp, self$design, self$diseases)$
   set_init_prvl(sp = sp, design_ = self$design)
})

private$primary_prevention_scn(sp) # apply primary pevention scenario

lapply(self$diseases, function(x) {
  x$set_rr(sp, self$design)$
    set_incd_prb(sp, self$design)$
    set_dgns_prb(sp, self$design)$
    set_mrtl_prb(sp, self$design)
})

private$secondary_prevention_scn(sp) # apply secondary pevention scenario

mc_ = 1L
scenario_nam = "sc0"

l <- private$mk_scenario_init(sp, scenario_nam)
simcpp(sp$pop, l, sp$mc)


sp$pop[, tmp := sum(wt_immrtl), keyby = .(year, age, sex)]

set(sp$pop, NULL, "wt", 0)

sp$pop[!is.na(all_cause_mrtl), wt := wt_immrtl * tmp / sum(wt_immrtl),
    by = .(year, age, sex)
]

sp$pop <- sp$pop[all_cause_mrtl >= 0L &
         year >= self$design$sim_prm$init_year_long &
         between(age, self$design$sim_prm$ageL, self$design$sim_prm$ageH), ]

setkey(sp$pop, pid, year)

sp$pop[, pid_mrk := mk_new_simulant_markers(pid)]

# apply ESP weights
to_agegrp(sp$pop, 5, 99)

absorb_dt(sp$pop, private$esp_weights)

sp$pop[, wt_esp := wt_esp * unique(wt_esp) / sum(wt_esp),
       by = .(year, agegrp, sex)] # NOTE keyby changes the key