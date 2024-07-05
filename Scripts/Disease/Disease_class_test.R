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

# lapply(diseases, function(x) x$harmonise_epi_tables(sp))

lapply(diseases, function(x) {
    print(x)
    x$gen_parf_files(design)
})
lapply(diseases, function(x) {
    print(x)
    x$gen_parf(sp, design, diseases)
})

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


