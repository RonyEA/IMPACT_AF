# change L6 to L10 and remove "#" from L54 to L84 for IMPACT-NCD-JAPAN
source("./global.R")
design <- Design$new("./inputs/sim_design.yaml")
# RR ----
# Create a named list of Exposure objects for the files in ./inputs/RR
fl <- list.files(path = "./inputs/RR", pattern = ".csvy$", full.names = TRUE)
RR <- future_lapply(fl, Exposure$new, design, future.seed = 950480304L)
# RR <- vector("list", length(fl))
# for (i in seq_along(fl)) {
#     print(fl[i])
#     RR[[i]] <- Exposure$new(fl[i], design)
# }
names(RR) <- sapply(RR, function(x) x$get_name())
invisible(future_lapply(RR, function(x) {
    x$gen_stochastic_effect(design, overwrite = TRUE, smooth = FALSE)
},
future.seed = 627524136L))
# NOTE smooth cannot be exported to Design for now, because the first time
# this parameter changes we need logic to overwrite unsmoothed files
rm(fl)
#

# Generate diseases ----
diseases <- lapply(design$sim_prm$diseases, function(x) {
    x[["design_"]] <- design
    x[["RR"]] <- RR
    do.call(Disease$new, x)
})
names(diseases) <- sapply(design$sim_prm$diseases, `[[`, "name")

mk_scenario_init2 <- function(scenario_name, diseases_, sp, design_) {
    # scenario_suffix_for_pop <- paste0("_", scenario_name) # TODO get suffix from design
    scenario_suffix_for_pop <- scenario_name
    list(
        "exposures"          = design_$sim_prm$exposures,
        "scenarios"          = design_$sim_prm$scenarios, # to be generated programmatically
        "scenario"           = scenario_name,
        "kismet"             = design_$sim_prm$kismet, # If TRUE random numbers are the same for each scenario.
        "init_year"          = design_$sim_prm$init_year,
        "pids"               = "pid",
        "years"              = "year",
        "ages"               = "age",
        "ageL"               = design_$sim_prm$ageL,
        "all_cause_mrtl"     = paste0("all_cause_mrtl", scenario_suffix_for_pop),
        "cms_score"          = paste0("cms_score", scenario_suffix_for_pop),
        "cms_count"          = paste0("cms_count", scenario_suffix_for_pop),
        "strata_for_outputs" = c("pid", "year", "age", "sex"),
        "diseases"           = lapply(diseases_, function(x) x$to_cpp(sp, design_))
    )
}

# sim <- SynthPop$new(0L, design)
# sim$write_synthpop(1:500)
# sim$delete_synthpop(NULL)
# ll <- sim$gen_synthpop_demog(design)
sp  <- SynthPop$new(1L, design)

self <- diseases$af$.__enclos_env__$self
private <- diseases$af$.__enclos_env__$private

design_ <- design
diseases_ <- diseases

# lapply(diseases, function(x) x$harmonise_epi_tables(sp))

# ---- gen_parf_files
popsize = 100
check = design_$sim_prm$logs
keep_intermediate_file = TRUE
bUpdateExistingDiseaseSnapshot = TRUE




parf_dt <- read_fst(private$parf_filenam, as.data.table = TRUE)

yrs <- design_$sim_prm$init_year_long

tt <- self$get_incd(yrs, sp$mc_aggr, design_ = design_)[, year := NULL]
if ("country" %in% names(tt)) tt[, country := NULL]

absorb_dt(parf_dt, tt) # now mu is incd


clbr <- fread("./simulation/calibration_prms.csv",  
            select = c("year", "age", "sex", paste0(self$name, "_incd_clbr_fctr"))
            )


setnafill(parf_dt, "c", fill = 0, cols = "mu") # fix for prostate and breast cancer
  
parf_dt <- absorb_dt(fread("./simulation/calibration_prms.csv",
  select =
    c("year", "age", "sex", paste0(self$name, "_incd_clbr_fctr"))
), parf_dt)



colnam <-
    setdiff(names(parf_dt), intersect(names(sp$pop), names(parf_dt)))



private$parf <- parf_dt[sp$pop, on = .NATURAL, ..colnam]

if ("p0" %in% names(private$parf)) {
    setnafill(private$parf, "c", fill = 0, cols = "p0") # fix for prostate and breast cancer
}

if ("mu" %in% names(private$parf)) {
    setnafill(private$parf, "c", fill = 0, cols = "mu") # fix for prostate and breast cancer
}






































lapply(diseases, function(x) {
    print(x)
    x$set_init_prvl(sp, design)
})

lapply(diseases, function(x) {
    print(x)
    x$set_rr(sp, design)
})

lapply(diseases, function(x) {
    print(x)
    x$set_incd_prb(sp, design)
})
lapply(diseases, function(x) {
    print(x)
    x$set_dgns_prb(sp, design)
})
lapply(diseases, function(x) {
    print(x)
    x$set_mrtl_prb(sp, design)
})


# old <- read_fst("/home/ckyprid/My_Models/IMPACTncd_Japan/backup_inputs_sep23/disease_burden/nonmodelled_ftlt.fst", as.data.table = T)
# new <- read_fst("/home/ckyprid/My_Models/IMPACTncd_Japan/inputs/disease_burden/nonmodelled_ftlt.fst", as.data.table = T)
# # new[, c("mu", "mu_lower", "mu_upper") := list(mu/1e5, mu_lower/1e5, mu_upper/1e5)]
# setnames(new, c("Rate", "Rate_lower", "Rate_upper"), c("mu2", "mu_lower", "mu_upper"))
# write_fst(new, "/home/ckyprid/My_Models/IMPACTncd_Japan/inputs/disease_burden/nonmodelled_ftlt.fst")

l <- mk_scenario_init2("", diseases, sp, design)
simcpp(sp$pop, l, sp$mc)

sp$update_pop_weights()
sp$pop[, mc := sp$mc_aggr]

# self <- IMPACTncd$.__enclos_env__$self
# private <- IMPACTncd$.__enclos_env__$private

self <- diseases$nonmodelled$.__enclos_env__$self
private <- diseases$nonmodelled$.__enclos_env__$private
design_ <- design
diseases_ <- diseases
popsize <- 100
check <- design_$sim_prm$logs
keep_intermediate_file <- TRUE
bUpdateExistingDiseaseSnapshot <- TRUE
mc_iter <- mc_ <- 1


