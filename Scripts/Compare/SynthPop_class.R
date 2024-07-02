## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by
## Chris Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngl is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the Free
## Software Foundation; either version 3 of the License, or (at your option) any
## later version. This program is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
## Public License for more details. You should have received a copy of the GNU
## General Public License along with this program; if not, see
## <http://www.gnu.org/licenses/> or write to the Free Software Foundation,
## Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



# From
# https://stackoverflow.com/questions/33424233/how-do-i-tell-an-r6-class-what-to-do-with-square-brackets
# Allows data.table syntax to the R6class object directly. Assumes it has a
# field 'pop' that is a data.table
#' @export
`[.SynthPop` <- function(x, ...) x$pop[...]

#' R6 Class representing a synthetic population
#'
#' @description
#' A synthpop has a `pop` field that contains the life course of simulants in a
#' `data.table`.
#'
#' @details
#' To be completed...
#'
#' @export
SynthPop <-
    R6::R6Class(
        classname = "SynthPop",

        # public ------------------------------------------------------------------
        public = list(
            #' @field mc The Monte Carlo iteration of the synthetic population fragment. Every
            #'   integer generates a unique synthetic population fragment.
            mc = NA,

            #' @field mc_aggr The Monte Carlo iteration of the synthetic population to
            #'   be used when multiple synthetic population fragments getting
            #'   aggregated. For instance if the synthpop consists of 2 fragments,
            #'   mc_aggr will be the same for both, but mc will differ. It ensures
            #'   correct seeds for the RNGs during the simulation for the RRs and the
            #'   lags.
            mc_aggr = NA,

            #' @field metadata Metadata of the synthpop.
            metadata = NA,

            #' @field pop The data.table that contains the life-course of simulants.
            #'   If the file exists, it is loaded from disk. If it doesn't, it is
            #'   first generated, then saved to disk, and then loaded from disk.
            pop = NA,

            # initialize ----
            #' @description Create a new SynthPop object.
            #' If a synthpop file in \code{\link[fst]{fst-package}} format already
            #' exists, then the synthpop is loaded from there. Otherwise it is
            #' generated from scratch and then saved as `filename` in
            #' \code{\link[fst]{fst-package}} format. Two additional files are saved
            #' for each 'synthpop'. A metadata file, and an index file.
            #' @param mc_ The Monte Carlo iteration of the synthetic population. Each
            #'   integer generates a unique synthetic population. If `mc = 0` an
            #'   object with an empty synthpop is initiated.
            #' @param design_ A \code{\link[IMPACTaf]{Design}} object.
            #' @param synthpop_dir_ The directory where 'SynthPop' objects are stored.
            #'   The synthpop file in \code{\link[fst]{fst-package}} format. If
            #'   `filename` already exists, then the synthpop is loaded from there.
            #'   Otherwise it is generated from scratch and then saved as `filename`
            #'   in \code{\link[fst]{fst-package}} format. Two additional files are
            #'   saved for each 'synthpop'. A metadata file, and an index file.
            #' @return A new `SynthPop` object.
            #' @examples
            #' design <- Design$new("./validation/design_for_trends_validation.yaml")
            #' POP$write_synthpop(1:6, design)
            #' POP <- SynthPop$new(4L, design)
            #' POP$print()
            #' POP$count_synthpop()
            #'
            #' POP$delete_synthpop(1L)
            #' POP$delete_synthpop(5:6)
            #' POP$get_filename()
            initialize = function(mc_, design_) {
                stopifnot(length(mc_) == 1L, is.numeric(mc_), ceiling(mc_) >= 0L)
                stopifnot("Design" %in% class(design_))

                mc_ <- as.integer(ceiling(mc_))
                # Create synthpop_dir if it doesn't exists
                # NOTE code below is duplicated in Simulation class. This is intentional
                if (!dir.exists(design_$sim_prm$synthpop_dir)) {
                    dir.create(design_$sim_prm$synthpop_dir, recursive = TRUE)
                    message(paste0(
                        "Folder ", design_$sim_prm$synthpop_dir,
                        " was created"
                    ))
                }

                # get unique lsoas
                # lsoas <- private$get_unique_LSOAs(design_)

                private$checksum <- private$gen_checksum(design_)

                self$mc <- mc_
                self$mc_aggr <-
                    as.integer(ceiling(mc_ / design_$sim_prm$n_synthpop_aggregation))

                private$design <- design_
                private$synthpop_dir <- design_$sim_prm$synthpop_dir

                if (mc_ > 0) {
                    # Logic to reuse a synthpop with larger age range if it exists (an expansion could be used for sim horizon)) (WIP)
                    #   if (design_$sim_prm$ageH == 99L) {
                    #     private$filename <- private$gen_synthpop_filename(mc_, private$checksum, design_)
                    #   } else {
                    #     for (age_ in design_$sim_prm$ageH:99)
                    #     original_ageH <- design_$sim_prm$ageH
                    #     design_$sim_prm$ageH <- age_
                    #     new_checksum <- private$gen_checksum(design_)
                    #     potential_filename <- private$gen_synthpop_filename(mc_, new_checksum, design_)
                    #     design_$sim_prm$ageH <- original_ageH

                    #     if (all(sapply(private$filename, file.exists))) {
                    #       private$filename <- potential_filename
                    #       break
                    #     }
                    # }
                    private$filename <- private$gen_synthpop_filename(mc_, private$checksum, design_)
                    # logic for the synthpop load
                    files_exist <- sapply(private$filename, file.exists)
                    if (all(!files_exist)) {
                        # No files exist. Create the synthpop and store the file on disk (no
                        # parallelism)
                        private$gen_synthpop(
                            mc_,
                            private$filename,
                            design_
                        )
                    } else if (file.exists(private$filename$metafile) &&
                        !all(files_exist)) {
                        # Metafile exists but not all three files. It means that most likely
                        # a generate_synthpop() is still running. So the function waits
                        # until the file is created before it proceeds to load it. Note that
                        # if this is not the case then the loop is infinite!!!
                        while (!all(sapply(private$filename, file.exists))) {
                            Sys.sleep(5)
                            if (design_$sim_prm$logs) {
                                message("Metafile exists without a synthpop file. Check for synthpop after 5 sec.")
                            }
                        }
                        # Ensure the file write is complete (size stable)
                        if (design_$sim_prm$logs) {
                            message("Synthpop file found.")
                        }

                        sz1 <- file.size(private$filename$synthpop)
                        Sys.sleep(3)
                        sz2 <- file.size(private$filename$synthpop)
                        while (sz1 != sz2) {
                            if (design_$sim_prm$logs) {
                                message("Synthpop file size increases.")
                            }
                            sz1 <- file.size(private$filename$synthpop)
                            Sys.sleep(3)
                            sz2 <- file.size(private$filename$synthpop)
                        }
                        if (design_$sim_prm$logs) {
                            message("Synthpop file stabilised.")
                        }
                    } else if (!file.exists(private$filename$metafile) &&
                        !all(files_exist)) {
                        # Metafile doesn't exist but some other files exist. In this case
                        # delete everything and start from scratch
                        self$delete_incomplete_synthpop()
                        private$gen_synthpop(
                            mc_,
                            private$filename,
                            design_
                        )
                    }
                    # No need to provision for case when all file present. The following
                    # lines handle this case anyway

                    if (design_$sim_prm$load_simulants_rn) {
                        exclude_cols_ <- c()
                    } else { # if not load_simulants_rn = TRUE
                        exclude_cols_ <- c(
                            "rank_Fruit_vege",
                            "rankstat_Smoking_act",
                            "rankstat_Smoking_ex",
                            "rankstat_Med_HT",
                            "rankstat_Med_HL",
                            "rankstat_Med_DM",
                            "rank_PA_days",
                            "rank_BMI",
                            "rank_HbA1c",
                            "rank_LDLc",
                            "rank_SBP",
                            "rankstat_Smoking_number"
                        )
                    }
                    self$pop <- private$get_synthpop(exclude_cols = exclude_cols_)
                    self$metadata <- yaml::read_yaml(private$filename$metafile)

                    if (design_$sim_prm$logs) self$print()
                }
                invisible(self)
            },

            #  update_design ----
            #' @description
            #' Updates the Design object that is stored in the SynthPop object.
            #' @param design_ A design object with the simulation parameters.
            #' @return The invisible self for chaining.

            update_design = function(design_ = design) {
                if (!inherits(design, "Design")) {
                    stop("Argument design_ needs to be a Design object.")
                }

                private$design <- design
                invisible(self)
            },

            # update_pop_weights ----
            #' @description
            #' Updates the wt_immrtl to account for mortality in baseline scenario.
            #' @param scenario_nam The scenario name. Logic is different if "sc0".
            #' @return The invisible self for chaining.
            update_pop_weights = function(scenario_nam = "sc0") {
                if (scenario_nam == "sc0" && !"wt" %in% names(self$pop)) { # baseline
                    self$pop[, tmp := sum(wt_immrtl), keyby = .(year, age, sex)]
                    set(self$pop, NULL, "wt", 0)
                    self$pop[!is.na(all_cause_mrtl), wt := wt_immrtl * tmp / sum(wt_immrtl),
                        by = .(year, age, sex)
                    ]

                    self$pop[, tmp := NULL]
                } else if (scenario_nam != "sc0" && !"wt" %in% names(self$pop)) {
                    # For policy scenarios.

                    fnam <- file.path(private$design$sim_prm$output_dir, paste0("lifecourse/", self$mc_aggr, "_lifecourse.csv.gz"))

                    x <- file.path(private$design$sim_prm$output_dir, paste0("lifecourse/", self$mc_aggr, "_lifecourse.csv.gz"))

                    t0 <- fread(fnam,
                        select = list(integer = c("pid", "year"), character = "scenario", numeric = "wt"),
                        key = c("scenario", "pid", "year")
                    )[scenario == "sc0", ] # wt for sc0

                    # For some reason pid and year get read incorrectly as character sometimes
                    t0[, pid := as.integer(pid)]
                    self$pop[, pid := as.integer(pid)]
                    t0[, year := as.integer(year)]
                    self$pop[, year := as.integer(year)]

                    self$pop[t0, on = c("pid", "year"), wt := i.wt]
                    self$pop[is.na(all_cause_mrtl), wt := 0]
                    self$pop[is.na(wt), wt := wt_immrtl]
                } else {
                    stop("The baseline scenario need to be named 'sc0' and simulated first, before any policy scenarios.") # TODO more informative message
                }

                invisible(self)
            },

            # delete_synthpop ----
            #' @description
            #' Delete (all) synthpop files in the synthpop directory.
            #' @param mc_ If `mc_ = NULL`, delete all files in the synthpop directory.
            #'   If `mc_` is an integer vector delete the specific synthpop files
            #'   including the metadata and index files.
            #' @param check_checksum If  `TRUE` only delete files with the same
            #'   checksum as the synthpop. Only relevant when `mc_ = NULL`.
            #' @param invert If `TRUE` (default is `FALSE`) keeps files with the same
            #'   checksum as the synthpop and deletes all other synthpops. Only
            #'   relevant when `mc_ = NULL` and `check_checksum = TRUE`.
            #' @return The invisible `SynthPop` object.
            delete_synthpop = function(mc_, check_checksum = TRUE, invert = FALSE) {
                if (missing(mc_)) stop("Use mc_ = NULL if you want to delete all synthpop files.")
                if (is.null(mc_)) {
                    if (check_checksum) {
                        fl <- list.files(
                            private$synthpop_dir,
                            pattern = paste0("^synthpop_", private$checksum),
                            full.names = TRUE,
                            recursive = TRUE
                        )
                        if (invert) {
                            fl2 <- list.files(
                                private$synthpop_dir,
                                pattern = "^synthpop_",
                                full.names = TRUE,
                                recursive = TRUE
                            )
                            fl <- setdiff(fl2, fl)
                        }
                    } else {
                        fl <- list.files(
                            private$synthpop_dir,
                            pattern = "^synthpop_",
                            full.names = TRUE,
                            recursive = TRUE
                        )
                    }
                    file.remove(fl)
                } else if (length(mc_) == 1L &&
                    is.numeric(mc_) && ceiling(mc_) > 0L) {
                    fl <- unlist(
                        private$gen_synthpop_filename(
                            mc_,
                            private$checksum,
                            private$design
                        )
                    )
                    file.remove(fl)
                } else if (length(mc_) > 1L &&
                    all(is.numeric(mc_)) && all(ceiling(mc_) > 0L)) {
                    fl <-
                        lapply(
                            mc_,
                            private$gen_synthpop_filename,
                            private$checksum,
                            private$design
                        )
                    fl <- unlist(fl)
                    file.remove(fl)
                } else {
                    message("mc_ need to be NULL or numeric. Nothing was deleted.")
                }

                return(invisible(self))
            },

            # delete_incomplete_synthpop ----
            #' @description
            #' Check that every synthpop file has a metafile and an index file. Delete
            #' any orphan files.
            #' @param check_checksum If  `TRUE` only delete incomplete group files
            #'   with the same checksum as the synthpop.
            #' @return The invisible `SynthPop` object.
            delete_incomplete_synthpop =
                function(check_checksum = TRUE) {
                    if (check_checksum) {
                        f1 <- paste0("^synthpop_", private$checksum, ".*\\.fst$")
                        f2 <- paste0("^synthpop_", private$checksum, ".*_meta\\.yaml$")
                    } else {
                        f1 <- "^synthpop_.*\\.fst$"
                        f2 <- "^synthpop_.*_meta\\.yaml$"
                    }

                    files <-
                        list.files(private$synthpop_dir, f1)
                    # remove indx files
                    files <- sub("\\.fst$", "", files)
                    metafiles <-
                        list.files(private$synthpop_dir, f2)
                    metafiles <- sub("_meta\\.yaml$", "", metafiles)

                    to_remove <- setdiff(metafiles, files)
                    if (length(to_remove) > 0) {
                        to_remove <- paste0(to_remove, "_meta.yaml")
                        file.remove(file.path(private$synthpop_dir, to_remove))
                    }

                    to_remove <- setdiff(files, metafiles)
                    if (length(to_remove) > 0) {
                        to_remove2 <- paste0(to_remove, ".fst")
                        file.remove(file.path(private$synthpop_dir, to_remove2))
                    }

                    return(invisible(self))
                },

            # check_integridy ----
            #' @description
            #' Check the integrity of (and optionally delete) .fst files by checking
            #' their metadata are readable.
            #' @param remove_malformed If `TRUE`, delete all malformed .fst files and
            #'   their associated files.
            #' @param check_checksum If  `TRUE` only check files with the same
            #'   checksum as the synthpop.
            #' @return The invisible `SynthPop` object.
            check_integridy =
                function(remove_malformed = FALSE,
                         check_checksum = TRUE) {
                    if (check_checksum) {
                        pat <- paste0("^synthpop_", private$checksum, ".*\\.fst$")
                    } else {
                        pat <- "^synthpop_.*\\.fst$"
                    }

                    files <-
                        list.files(private$synthpop_dir,
                            pat,
                            full.names = TRUE
                        )
                    if (length(files) > 0L) {
                        malformed <- sapply(files, function(x) {
                            out <- try(metadata_fst(x), silent = TRUE)
                            out <- inherits(out, "try-error")
                            out
                        }, USE.NAMES = FALSE)


                        des <- sum(malformed)

                        if (remove_malformed) {
                            if (des == 0L) {
                                message(paste0(des, " malformed fst file(s)"))
                            } else {
                                # des != 0L
                                message(paste0(des, " malformed fst file(s)..."))
                                to_remove <- files[malformed]

                                # then remove other files
                                to_remove <- gsub(".fst$", "", to_remove)

                                # _meta.yaml
                                tr <- paste0(to_remove, "_meta.yaml")
                                file.remove(tr[file.exists(tr)])
                                # .fst
                                tr <- paste0(to_remove, ".fst")
                                file.remove(tr[file.exists(tr)])

                                message("...now deleted!")
                            }
                        } else {
                            # remove_malformed = FALSE
                            message(paste0(des, " malformed fst file(s)"))
                        }
                        file_indx <- read_fst(file, as.data.table = TRUE, columns = "year")[, .(from = min(.I), to = max(.I)), keyby = "year"][year == 2000L + design_$sim_prm$init_year]
                    } else {
                        # if length(files) == 0
                        message("no .fst files found.")
                    }
                    return(invisible(self))
                },


            # count_synthpop ----
            #' @description
            #' Count the synthpop files in a directory. It includes files without
            #' metafiles and index files.
            #' @return The invisible `SynthPop` object.
            count_synthpop =
                function() {
                    out <- list()
                    # folder size
                    files <-
                        list.files(private$synthpop_dir, full.names = TRUE)
                    if (length(files) > 0L) {
                        vect_size <- sapply(files, file.size)
                        out$`synthpop folder size (Gb)` <-
                            signif(sum(vect_size) / (1024^3), 4) # Gb

                        # synthpops with same checksum
                        files <- list.files(
                            private$synthpop_dir,
                            paste0("^synthpop_", private$checksum, ".*\\.fst$")
                        )

                        out$`synthpop meta files with same checksum` <-
                            length(list.files(
                                private$synthpop_dir,
                                paste0("^synthpop_", private$checksum, ".*_meta\\.yaml$")
                            ))

                        # synthpops with any checksum
                        files <-
                            list.files(private$synthpop_dir, "^synthpop_.*\\.fst$")

                        out$`synthpop meta files with any checksum` <-
                            length(list.files(private$synthpop_dir, "^synthpop_.*_meta\\.yaml$"))


                        cat(paste0(names(out), ": ", out, "\n"))
                    } else {
                        # if length(files) == 0L
                        cat("no files found.")
                    }
                    return(invisible(self))
                },

            # get_checksum ----
            #' @description
            #' Get the synthpop checksum.
            #' @param x One of "all", "synthpop" or "metafile". Can be abbreviated.
            #' @return The invisible `SynthPop` object.
            get_checksum = function() {
                out <- private$checksum
                names(out) <- "Checksum"
                cat(paste0(names(out), ": ", out))
                invisible(self)
            },

            # get_filename ----
            #' @description
            #' Get the synthpop file paths.
            #' @param x One of "all", "synthpop" or "metafile". Can be abbreviated.
            #' @return The invisible `SynthPop` object.
            get_filename = function(x = c("all", "synthpop", "metafile")) {
                if (self$mc == 0L) {
                    print("Not relevant because mc = 0L")
                } else {
                    x <- match.arg(x)
                    switch(x,
                        all      = print(private$filename),
                        synthpop = print(private$filename[["synthpop"]]),
                        metafile = print(private$filename[["metafile"]])
                    )
                }
                invisible(self)
            },


            # get_design ----
            #' @description
            #' Get the synthpop design.
            #' @return The invisible `SynthPop` object.
            get_design = function() {
                # print(private$design)
                # invisible(self)
                private$design
            },

            # get_dir ----
            #' @description
            #' Get the synthpop dir.
            #' @return The invisible `SynthPop` object.
            get_dir = function() {
                print(private$synthpop_dir)
                invisible(self)
            },

            # gen_synthpop_demog ----
            #' @description
            #' Generate synthpop sociodemographics, random sample of the population.
            #' @param design_ A Design object,
            #' @param month April or July are accepted. Use July for mid-year
            #'   population estimates.
            #' @return An invisible `data.table` with sociodemographic information.
            # Change-for-IMPACT-NCD-Japan, we use population in October because the official population estimate was baed on population in October
            gen_synthpop_demog =
                function(design_) {
                    # load dt
                    # file <- "../Data/Final_Pop_Data/pop_proj_eu.fst"
                    file <- "./inputs/pop_observed/pop_observed_eu.fst"

                    dt <- read_fst(file, as.data.table = TRUE)
                    dt_meta <- metadata_fst(file)
                    stopifnot(
                        "Population size file need to be keyed by year" =
                            identical("year", dt_meta$keys[1])
                    )

                    countries <- design_$sim_prm$country
                    print(countries)

                    # countrywise synth pop.
                    if (("country" %in% colnames(dt)) & (length(unique(dt[, country])) > 1)) {
                        dt <- dt[country %in% countries]

                        dt <- dt[, .(pops = sum(pops)), by = .(year, age, sex)]
                    }

                    dt <- dt[year == 2000L + design_$sim_prm$init_year]

                    dt <- dt[age == "100+", age := 100]
                    dt <- dt[, age := as.numeric(age)]

                    dt <- dt[age %in% c(design_$sim_prm$ageL:design_$sim_prm$ageH)] # TODO: Needs fix? Check old version

                    
                    dt <- dt[, prbl := pops / sum(pops)][, `:=`(pops = NULL)]
                    # dt <- dt[, prbl := pops / sum(pops)][, `:=`(reg = NULL, pops = NULL)]

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

                    return(invisible(dt))
                },

            # write_synthpop ----
            #' @description
            #' Generate synthpop files in parallel, using foreach, and writes them to
            #' disk. It skips files that are already on disk.
            #' Note: the backend for foreach needs to be initialised before calling
            #' the function.
            #' @param mc_ An integer vector for the Monte Carlo iteration of the
            #'   synthetic population. Each integer generates a unique synthetic
            #'   population.
            #' @return The invisible `SynthPop` object.
            write_synthpop = function(mc_) {
                stopifnot(all(is.numeric(mc_)), all(ceiling(mc_) > 0L))
                on.exit(self$delete_incomplete_synthpop(), add = TRUE)
                mc_ <- as.integer(ceiling(mc_))

                if (.Platform$OS.type == "windows") {
                    # TODO update to make compatible with windows
                    cl <-
                        makeCluster(private$design$sim_prm$clusternumber) # used for clustering. Windows compatible
                    registerDoParallel(cl)
                } else {
                    registerDoParallel(private$design$sim_prm$clusternumber) # used for forking. Only Linux/OSX compatible
                }

                foreach(
                    mc_iter = mc_,
                    .inorder = FALSE,
                    .verbose = private$design$sim_prm$logs,
                    .packages = c(
                        "R6",
                        "gamlss.dist",
                        # For distr in prevalence.R
                        "dqrng",
                        "qs",
                        "fst",
                        "CKutils",
                        "IMPACTncdEngl",
                        "data.table"
                    ),
                    .export = NULL,
                    .noexport = NULL # c("time_mark")
                ) %dopar% {
                    data.table::setDTthreads(private$design$sim_prm$n_cpus)
                    fst::threads_fst(private$design$sim_prm$n_cpus)
                    filename <-
                        private$gen_synthpop_filename(
                            mc_iter,
                            private$checksum,
                            private$design
                        )

                    # logic for the synthpop load
                    files_exist <- sapply(filename, file.exists)
                    if (all(!files_exist)) {
                        # No files exist. Create the synthpop and store
                        # the file on disk
                        private$gen_synthpop(
                            mc_iter,
                            filename,
                            private$design
                        )
                    } else if (file.exists(filename$metafile) &&
                        !all(files_exist)) {
                        # Metafile exists but not all three files. It means
                        # that most likely a generate_synthpop() is still running. So the
                        # function waits until the file is created before it proceeds to
                        # load it. Note that if this is not the case then the loop is
                        # infinite!!!
                        while (!all(sapply(filename, file.exists))) {
                            Sys.sleep(5)
                        }

                        # Ensure the file write is complete (size stable)
                        sz1 <- file.size(filename$synthpop)
                        Sys.sleep(3)
                        sz2 <- file.size(filename$synthpop)
                        while (sz1 != sz2) {
                            sz1 <- file.size(filename$synthpop)
                            Sys.sleep(3)
                            sz2 <- file.size(filename$synthpop)
                        }
                    } else if (!file.exists(filename$metafile) &&
                        !all(files_exist)) {
                        # Metafile doesn't exist but some other files exist. In this case
                        # delete everything and start from scratch
                        self$delete_incomplete_synthpop()
                        private$gen_synthpop(
                            mc_iter,
                            filename,
                            private$design
                        )
                    }
                    # No need to provision for case when all files present.

                    return(NULL)
                }
                if (exists("cl")) {
                    stopCluster(cl)
                }

                invisible(self)
            },

            # get_risks ----
            #' @description Get the risks for all individuals in a synthetic
            #'   population for a disease.
            #' @param disease_nam The disease that the risks will be returned.
            #' @return A data.table with columns for pid, year, and all associated
            #'   risks if disease_nam is specified. Else a list of data.tables for all
            #'   diseases.
            get_risks = function(disease_nam) {
                if (missing(disease_nam)) {
                    return(private$risks)
                } else {
                    stopifnot(is.character(disease_nam))
                    return(private$risks[[disease_nam]])
                }
            },

            # store_risks ----
            #' @description Stores the disease risks for all individuals in a synthetic
            #'   population in a private list.
            #' @param disease_nam The disease that the risks will be stored.
            #' @return The invisible self for chaining.
            store_risks = function(disease_nam) {
                stopifnot(is.character(disease_nam))

                nam <- grep("_rr$", names(self$pop), value = TRUE)

                private$risks[[disease_nam]] <-
                    self$pop[, .SD, .SDcols = c("pid", "year", nam)]

                self$pop[, (nam) := NULL]
                invisible(self)
            },

            # print ----
            #' @description
            #' Prints the synthpop object metadata.
            #' @return The invisible `SynthPop` object.
            print = function() {
                print(c(
                    "path" = ifelse(self$mc == 0L,
                        "Not relevant because mc = 0L",
                        private$filename$synthpop
                    ),
                    "checksum" = private$checksum,
                    "mc" = self$mc,
                    self$metadata
                ))
                invisible(self)
            }
        ),



        # private -----------------------------------------------------------------
        private = list(
            filename = NA,
            checksum = NA,
            # The design object with the simulation parameters.
            design = NA,
            synthpop_dir = NA,
            risks = list(), # holds the risks for all individuals

            # Special deep copy for data.table. Use POP$clone(deep = TRUE) to
            # dispatch. Otherwise a reference is created
            deep_clone = function(name, value) {
                if ("data.table" %in% class(value)) {
                    data.table::copy(value)
                } else if ("R6" %in% class(value)) {
                    value$clone()
                } else {
                    # For everything else, just return it. This results in a shallow
                    # copy of s3.
                    value
                }
            },

            # get a smaller design list only with characteristics that are important
            # for synthpop creation and define the uniqueness of the object. I.e. if
            # these parameters are different the synthpop has to have different
            # filename and vice-versa
            get_unique_characteristics = function(design_) {
                design_$sim_prm[c(
                    "n",
                    "sim_horizon_max",
                    "init_year_long",
                    "maxlag",
                    "ageL",
                    "ageH",
                    "jumpiness"
                )]
            },

            # gen_checksum ----
            # gen synthpop unique checksum for the given set of inputs
            gen_checksum =
                function(design_) {
                    # get a md5 checksum based on function arguments
                    # First get function call arguments
                    fcall <- private$get_unique_characteristics(design_)

                    years_age_id <-
                        digest(paste(fcall, sep = ",", collapse = ","),
                            serialize = FALSE
                        )
                    return(years_age_id)
                },

            # gen synthpop filename for the given set of inputs
            gen_synthpop_filename =
                function(mc_, checksum_, design_) {
                    return(list(
                        "synthpop" = normalizePath(
                            paste0(
                                design_$sim_prm$synthpop_dir, "/synthpop_", checksum_, "_",
                                design_$sim_prm$country,
                                "_",
                                mc_,
                                ".fst"
                            ),
                            mustWork = FALSE
                        ),
                        "metafile" = normalizePath(
                            paste0(
                                design_$sim_prm$synthpop_dir,
                                "/synthpop_",
                                checksum_,
                                "_",
                                design_$sim_prm$country,
                                "_",
                                mc_,
                                "_meta.yaml"
                            ),
                            mustWork = FALSE
                        )
                    ))
                },
            del_incomplete = function(filename_) {
                if (file.exists(filename_$metafile) &&
                    (!file.exists(filename_$synthpop)
                    )) {
                    suppressWarnings(sapply(filename_, file.remove))
                }
            },

            # gen_synthpop ----
            gen_synthpop = # returns NULL. Writes synthpop on disk
                function(mc_,
                         filename_,
                         design_) {
                    # Save synthpop metadata
                    if (!file.exists(filename_$metafile)) {
                        yaml::write_yaml(
                            private$get_unique_characteristics(design_),
                            filename_$metafile
                        )
                    }
                    # NOTE In shiny app if 2 users click the  button at the same time, 2
                    # functions will run almost concurrently with potential race condition

                    # To avoid edge cases when the function stopped prematurely and a
                    # metafile was created while the file was not. On.exit ensures that
                    # either both files exist or none.

                    on.exit(private$del_incomplete(filename_), add = TRUE)

                    dqRNGkind("pcg64")
                    SEED <-
                        2121870L # sample(1e7, 1) # Hard-coded for reproducibility
                    set.seed(SEED + mc_)
                    dqset.seed(SEED, mc_)

                    # Generate synthpops with sociodemographic and exposures information.

                    dt <- self$gen_synthpop_demog(design_)

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
                    dt <- dt[, sex := factor(sex)]

                    # Generate smk ----
                    # Prediction for smoking.

                    if (design_$sim_prm$logs) message("Generate Smoking Status")

                    tbl <- read_fst("./inputs/exposure_distributions/smk_IT_table.fst", as.data.table = TRUE)

                    tbl[, year := year + 2000L]
                    tbl <- tbl[year >= 2003]

                    col_nam <- setdiff(names(tbl), intersect(names(dt), names(tbl)))

                    absorb_dt(dt, tbl)

                    dt[, smk := factor(
                        (rankstat_smk > smk1) +
                            (rankstat_smk > smk2) +
                            (rankstat_smk > smk3) +
                            1,
                        levels = 1:3, labels = 1:3, ordered = FALSE
                    )]

                    # dt[, smk := (rankstat_smk > smk1) +
                    #     (rankstat_smk > smk2) +
                    #     (rankstat_smk > smk3) +
                    #     1]

                    dt[, smok_cig := fcase(
                        smk == 1, 0,
                        smk == 2, 1,
                        smk == 3, 15
                    )]

                    if (!design_$sim_prm$keep_simulants_rn) col_nam <- c(col_nam, "rankstat_smk")
                    dt[, c(col_nam) := NULL]

                    # Generate mm_clusters ----
                    if (design_$sim_prm$logs) message("Generate Multimorbidity Clusters")
                    tbl <-
                        read_fst("./inputs/exposure_distributions/af_mm_cluster_table.fst",
                            as.data.table = TRUE
                        )

                    col_nam <-
                        setdiff(names(tbl), intersect(names(dt), names(tbl)))

                    absorb_dt(dt, tbl)

                    dt[, mm_cluster := 
                        (rankstat_mm_cluster > no_morbidity) +
                            (rankstat_mm_cluster > unspecific) +
                            (rankstat_mm_cluster > neuropsychiatric) +
                            (rankstat_mm_cluster > complex) +
                            (rankstat_mm_cluster > eye) +
                            (rankstat_mm_cluster > musculoskeletal) +
                            (rankstat_mm_cluster > metabolic) +
                            (rankstat_mm_cluster > cardiovascluar)
                    ]

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

                    # tbl[,smk := as.integer(smk)]

                    tbl[, year := year + 2000]
                    tbl <- tbl[year >= 2003]

                    col_nam <-
                        setdiff(names(tbl), intersect(names(dt), names(tbl)))

                    absorb_dt(dt, tbl)

                    if (design_$sim_prm$logs) message("Generating BMI from qGB2")

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
                        "smk",
                        "smok_cig",
                        "mm_cluster"
                    )

                    xps_nam <- paste0(xps_tolag, "_curr_xps")
                    setnames(dt, xps_tolag, xps_nam)

                    if ("age100" %in% names(dt)) {
                        dt[, age := NULL]
                        setnames(dt, "age100", "age")
                    }

                    dt[, sex := factor(sex)]
                    dt[, year := as.integer(year)]
                    # dt[, smk_curr_xps := as.integer(smk_curr_xps)]

                    setkey(dt, pid, year) # Just in case
                    setcolorder(dt, c("pid", "year", "age", "sex")) # STRATA
                    setindexv(dt, c("year", "age", "sex")) # STRATA

                    if (design_$sim_prm$logs) message("Writing synthpop to disk")

                    write_fst(
                        dt,
                        filename_$synthpop,
                        90
                    ) # 100 is too slow
                    return(invisible(NULL))
                },

            # get_synthpop ----
            # Load a synthpop file from disk in full or in chunks.
            get_synthpop =
                function(exclude_cols = c()) {
                    mm_synthpop <- metadata_fst(private$filename$synthpop)
                    mm_synthpop <- setdiff(mm_synthpop$columnNames, exclude_cols)

                    # Read synthpop

                    dt <- read_fst(private$filename$synthpop,
                        columns = mm_synthpop,
                        as.data.table = TRUE
                    )
                    dt <- dt[between(
                        year - 2000L,
                        private$design$sim_prm$init_year - private$design$sim_prm$maxlag,
                        private$design$sim_prm$init_year + private$design$sim_prm$sim_horizon_fromGUI
                    ) &
                        between(
                            age,
                            private$design$sim_prm$ageL - private$design$sim_prm$maxlag,
                            private$design$sim_prm$ageH
                        )]

                    # Ensure pid does not overlap for files from different mc
                    new_n <-
                        it <- as.integer(ceiling(self$mc %% private$design$sim_prm$n_synthpop_aggregation))
                    if ((max(dt$pid) + it * 1e8) >= .Machine$integer.max) stop("pid larger than int32 limit.")
                    dt[, pid := as.integer(pid + it * 1e8)]

                    dt[, pid_mrk := mk_new_simulant_markers(pid)] # TODO Do I need this?

                    # dt[, pid_mrk := mk_new_simulant_markers(pid)] # TODO Do I need this?

                    # Ensure pid does not overlap for files from different mc
                    # new_n <- uniqueN(dt$pid)
                    # it <- as.integer(ceiling(self$mc %% private$design$sim_prm$n_synthpop_aggregation))
                    # it[it == 0L] <- private$design$sim_prm$n_synthpop_aggregation
                    # it <- it - 1L
                    # if (max(dt$pid + (private$design$sim_prm$n_synthpop_aggregation - 1) * new_n) < .Machine$integer.max) {
                    #  dt[, pid := as.integer(pid + it * new_n)]
                    # } else stop("pid larger than int32 limit.")

                    # generate population weights
                    private$gen_pop_weights(dt, private$design)

                    set(dt, NULL, "all_cause_mrtl", 0L)
                    set(dt, NULL, "cms_score", 0) # CMS score of diagnosed conditions
                    set(dt, NULL, "cms_count", 0L) # Count of diagnosed CMS conditions

                    setkey(dt, pid, year)
                    # dt[, dead := identify_longdead(all_cause_mrtl, pid_mrk)]
                    # dt[, ncc := clamp(
                    #   ncc - (chd_prvl > 0) - (stroke_prvl > 0) -
                    #     (poststroke_dementia_prvl > 0) -
                    #     (htn_prvl > 0) - (t2dm_prvl > 0) - (af_prvl > 0) -
                    #     (copd_prvl > 0) - (lung_ca_prvl > 0) -
                    #     (colon_ca_prvl > 0) -
                    #     (breast_ca_prvl > 0),
                    #   0L,
                    #   10L
                    # )]
                    # to be added back in the qaly fn. Otherwise when I prevent disease
                    # the ncc does not decrease.


                    invisible(dt)
                },

            # gen_pop_weights ----
            # Calculate weights so that their sum is the population of the area based
            # on ONS. It takes into account synthpop aggregation. So you need to sum
            # all the synthpops belong to the same aggregation to reach the total pop.
            # NOTE Wt are still incomplete because they assume everyone remains alive.
            # So baseline population underestimated as clearly some die
            gen_pop_weights = function(dt, design) {
                tt <-
                    read_fst("./inputs/pop_projections/pop_combined_eu.fst", as.data.table = TRUE) # Change-for-IMPACT-NCD-Japan
                
                tt <- tt[country == design$sim_prm$country]
                
                tt[age == "100+", age := 100]
                tt[, age := as.integer(age)]
                
                tt <- tt[
                    between(age, min(dt$age), max(dt$age)) &
                        between(year, min(dt$year), max(dt$year)),
                    .(pops = sum(pops)),
                    keyby = .(year, age, sex)
                ]
                
                
                dt[, wt_immrtl := .N, by = .(year, age, sex)]
                
                absorb_dt(dt, tt)
                
                dt[, wt_immrtl := pops / (wt_immrtl * design$sim_prm$n_synthpop_aggregation)]
                dt[, pops := NULL]

                invisible(dt)
            }
        )
    )
