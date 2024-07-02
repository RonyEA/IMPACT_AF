source("./global.R")
IMPACTncd <- Simulation$new("./inputs/sim_design_calibration.yaml")

self <- IMPACTncd$.__enclos_env__$self
private <- IMPACTncd$.__enclos_env__$private

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