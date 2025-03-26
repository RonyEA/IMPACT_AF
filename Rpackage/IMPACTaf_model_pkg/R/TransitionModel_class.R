#' TransitionModel Class
#'
#' This class encapsulates the functionality for running a Transition model simulation
#' with patient transitions, initial states, and a customizable transition matrix.

#'
#' @import data.table
#' @import matrixStats
#' @import yaml
#' 
#' @export
TransitionModel <- R6::R6Class(
  "TransitionModel",
  public = list(
    # TODO: Model path can be retrieved from a yaml?
    # TODO: n_patients, n_cycles and initial_transition.csv path from yaml?
    config = NULL,
    results = NULL,
    init_trans_mat = NULL,
    mc = NULL,
    mm_cluster_rename = NULL,


    #' @description
    #' Initializes the TransitionModel by loading configuration settings from a fixed YAML file path,
    #' reading the initial transition matrix, and preparing internal structures needed for
    #' simulating patient transitions.
    #'
    #' @details
    #' This method sets up the model by:
    #' \itemize{
    #'   \item Loading model configuration from `./inputs/transitions_config.yaml`
    #'   \item Renaming multimorbidity clusters for consistency
    #'   \item Loading and preparing the initial transition matrix
    #'   \item Counting available Monte Carlo (MC) iterations based on input files
    #'   \item Preloading transition probabilities from saved model coefficients
    #' }
    #'
    #' @return None. Side-effects include populating internal model fields and printing an initialization message.
    #'
    #' @examples
    #' tm <- TransitionModel$new()
    #'
    #' @keywords internal
    initialize = function() {
      self$config <- read_yaml("./inputs/transitions_config.yaml")

      self$mm_cluster_rename <- c(
        "Cardiovascular" = 1,
        "Complex" = 2,
        "Eye" = 3,
        "MSK" = 4,
        "Metabolic" = 5,
        "Neuropsychiatric" = 6,
        "Unspecified" = 7
      )

      self$init_trans_mat <- self$load_init_trans_mat()

      self$results <- list()

      self$mc <- length(list.files(path = "./outputs/lifecourse/", pattern = "_lifecourse\\.csv\\.gz$", full.name = TRUE))

      self$load_transition_probs()


      # self$get_transition_probs()
      message("Model Initialisation Completed")
    },


    #' Load and Filter Cohort Data
    #'
    #' @description
    #' Loads the life course data for a specific Monte Carlo iteration and filters the dataset
    #' to include only individuals aged 65 and above. Returns a subset of relevant variables.
    #'
    #' @param mc_counter Integer. The Monte Carlo iteration index used to identify the appropriate lifecourse file.
    #'
    #' @return A `data.table` with filtered cohort data containing the following columns:
    #' \itemize{
    #'   \item \code{pid} - Patient ID
    #'   \item \code{year} - Simulation year
    #'   \item \code{sex} - Sex of the individual
    #'   \item \code{agegrp} - Age group
    #'   \item \code{mm_cluster} - Multimorbidity cluster ID
    #'   \item \code{af_prvl} - Atrial fibrillation prevalence
    #' }
    #'
    #' @examples
    #' tm <- TransitionModel$new("path/to/config.yaml")
    #' cohort_data <- tm$get_cohort(1)
    get_cohort = function(mc_counter) {
      lifecourse_file <- sprintf("/home/ron/Projects/IMPACT_AF/outputs/lifecourse/%s_lifecourse.csv.gz", mc_counter)
      lifecourse <- fread(lifecourse_file)
      lifecourse <- lifecourse[age >= 65]
      lifecourse <- lifecourse[, .(pid, year, sex, agegrp, mm_cluster, af_prvl)]

      message("Cohort filtered: ", nrow(self$results), " patients above age 65.")

      return(lifecourse)
    },
    

    #' Load Initial Transition Matrix
    #'
    #' @description
    #' Loads the initial transition probability matrix for each patient subgroup defined by age group,
    #' sex, and multimorbidity cluster. The function reads from a file path defined in the configuration
    #' YAML and computes cumulative probabilities for the transition to each health state.
    #'
    #' @details
    #' This matrix is used to assign initial health states to patients during simulation initialization.
    #' The cumulative sum across health state columns enables probabilistic sampling.
    #'
    #' @return A `data.table` with cumulative probabilities for initial transitions, indexed by age group,
    #' sex, and multimorbidity cluster.
    #'
    #' @examples
    #' tm <- TransitionModel$new("path/to/config.yaml")
    #' init_matrix <- tm$load_init_trans_mat()
    #'
    #' @keywords internal
    load_init_trans_mat = function() {
      init_trans_mat <- fread(self$config$transition_matrix_path)
      init_trans_mat[, mm_cluster := self$mm_cluster_rename[mm_cluster]]

      # TODO: is it a robust way to use indexing using nums?
      init_trans_mat[,
        (names(init_trans_mat)[4:18]) := as.data.table(t(apply(.SD, 1, cumsum))),
        .SDcols = names(init_trans_mat)[4:18]
      ]

      return(init_trans_mat)
    },


    #' Initialize Patient States
    #'
    #' @description
    #' Initializes the health states for patients in the cohort based on random sampling
    #' using probabilities from the initial transition matrix. This is typically the first
    #' step in running a simulation for a given Monte Carlo iteration.
    #'
    #' @param mc_counter Integer. The Monte Carlo iteration index used to retrieve the relevant
    #' life course file and determine initial states.
    #'
    #' @details
    #' The method:
    #' \itemize{
    #'   \item Loads and filters the cohort for patients aged 65 and over.
    #'   \item Merges the cohort with the initial transition probabilities using age group, sex, and multimorbidity cluster.
    #'   \item Uses a random rank to probabilistically assign each patient to a health state based on cumulative probabilities.
    #' }
    #'
    #' @return A `data.table` representing the life course of each patient, now including the initialized `CurrentEvent` state.
    #'
    #' @examples
    #' tm <- TransitionModel$new("config/transition_config.yaml")
    #' initialized_lifecourse <- tm$initialize_states(1)
    #'
    #' @keywords internal
    initialize_states = function(mc_counter) {
      lifecourse <- self$get_cohort(mc_counter)

      p_c1 <- lifecourse[af_prvl > 0][, .(agegrp, year = min(year), sex, mm_cluster, af_prvl), by = pid]
      setorder(p_c1, pid, year, af_prvl)

      p_c1 <- p_c1[!duplicated(p_c1, by = c("pid", "year"))]
      p_c1[, rank := runif(.N)]

      setcolorder(p_c1, c("pid", "year", "agegrp", "sex", "mm_cluster", "af_prvl", "rank"))

      # merging initial transitions
      p_c1 <- p_c1[self$init_trans_mat, on = .(agegrp, sex, mm_cluster)]

      health_states <- setdiff(names(p_c1), c("pid", "year", "agegrp", "sex", "mm_cluster", "rank", "af_prvl", "mc"))

      # Melting the data to long format as colwise operations are more friendly to cumsum of probs
      p_c1_long <- melt(p_c1,
        id.vars = c("pid", "agegrp", "sex", "mm_cluster", "year", "rank"),
        measure.vars = health_states, variable.name = "state", value.name = "probs"
      )
      # Sorting the data
      setorder(p_c1_long, pid, probs)

      # Find the state where rank is greater than or equal to cumprob
      initial_states <- p_c1_long[, CurrentEvent := state[findInterval(rank[1], probs) + 1], by = .(pid, agegrp, sex, mm_cluster)][, .SD[1], by = pid]
      initial_states[, c("rank", "state", "probs") := NULL]

      cols <- intersect(names(lifecourse), names(initial_states))
      lifecourse <- merge(lifecourse, initial_states, by = cols, all.x = TRUE)

      return(lifecourse)
    },


    #' Load and Simulate Transition Probabilities
    #'
    #' @description
    #' Simulates transition probabilities between health states using coefficients from
    #' a multinomial logistic regression model. The method uses Monte Carlo sampling
    #' to reflect uncertainty in model parameters.
    #'
    #' @details
    #' This method:
    #' \itemize{
    #'   \item Loads the coefficient estimates and variance-covariance matrix of a transition model.
    #'   \item Uses multivariate normal sampling to generate simulated coefficients for each Monte Carlo iteration.
    #'   \item Applies the coefficients to a design matrix to calculate transition probabilities via softmax.
    #'   \item Assigns contextual variables (e.g., age group, sex, mm_cluster, CurrentEvent).
    #'   \item Saves the resulting simulated probabilities to disk for each Monte Carlo iteration.
    #' }
    #'
    #' @return None. Saves simulated transition probabilities as `.csv` files in the `./outputs/transitions/` directory.
    #'
    #' @examples
    #' tm <- TransitionModel$new("config/transition_config.yaml")
    #' tm$load_transition_probs()
    #'
    #' @keywords internal
    load_transition_probs = function() {
      # Extract coefficients and variance-covariance matrix
      coef_model <- readRDS("./inputs/transitions/Present_coef_model.rds") # Coefficients
      vcov_model <- readRDS("./inputs/transitions/Present_vcov_model.rds") # Variance-covariance matrix

      simulated_coefs <- MASS::mvrnorm(self$mc, mu = as.vector(t(coef_model)), Sigma = vcov_model)

      design_matrix <- readRDS("./inputs/transitions/Present_design_matrix.rds")
      rename_design_matrix <- c(
        "mypredclass7Cardiovascular" = "Cardiovascular",
        "mypredclass7Complex" = "Complex",
        "mypredclass7Eye" = "Eye",
        "mypredclass7MSK" = "MSK",
        "mypredclass7Metabolic" = "Metabolic",
        "mypredclass7Neuropsychiatric" = "Neuropsychiatric",
        "mypredclass7Unspecified" = "Unspecified"
      )
      # Rename only the matching columns
      colnames(design_matrix)[colnames(design_matrix) %in% names(rename_design_matrix)] <-
        rename_design_matrix[colnames(design_matrix)[colnames(design_matrix) %in% names(rename_design_matrix)]]


      num_categories <- length(rownames(coef_model)) + 1 # Number of outcome categories``
      simulated_probs <- array(NA, dim = c(self$mc, nrow(design_matrix), num_categories))

      # Efficient matrix operations - avoid loop wherever possible
      beta_sim_list <- lapply(1:self$mc, function(i) {
        beta_sim <- matrix(simulated_coefs[i, ], nrow = num_categories - 1, byrow = TRUE)
        beta_sim <- rbind(0, beta_sim) # Add baseline row

        # Compute probabilities using softmax
        linear_preds <- design_matrix %*% t(beta_sim)
        exp_preds <- exp(linear_preds)
        probs <- exp_preds / rowSums(exp_preds)

        as.data.table(cbind(design_matrix, probs))
      })


      final_sim_probs_dt <- rbindlist(beta_sim_list, idcol = "mc")

      setnames(final_sim_probs_dt, paste0("V", 29:47), c("NoEvent", rownames(coef_model)))

      print(colnames(final_sim_probs_dt))

      final_sim_probs_dt[, CurrentEvent := fifelse(
        `CurrentEventAB` == 1, "AB",
        fifelse(
          `CurrentEventAHS` == 1, "AHS",
          fifelse(
            `CurrentEventAHS_AB` == 1, "AHS_AB",
            fifelse(
              `CurrentEventAIS` == 1, "AIS",
              fifelse(
                `CurrentEventAIS_AB` == 1, "AIS_AB",
                fifelse(
                  `CurrentEventPB` == 1, "PB",
                  fifelse(
                    `CurrentEventPB_AHS` == 1, "PB_AHS",
                    fifelse(
                      `CurrentEventPB_AIS` == 1, "PB_AIS",
                      fifelse(
                        `CurrentEventPHS` == 1, "PHS",
                        fifelse(
                          `CurrentEventPIS` == 1, "PIS",
                          fifelse(
                            `CurrentEventPS_AB` == 1, "PS_AB",
                            fifelse(
                              `CurrentEventPS_PB` == 1, "PS_PB",
                              "NoEvent"
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )]

      # Assign 'sex'
      final_sim_probs_dt[, sex := fifelse(sexM == 1, "men", "women")]

      # Assign 'age_group'
      final_sim_probs_dt[, age_group := fifelse(
        `age_group70-74` == 1, "70-74",
        fifelse(
          `age_group75-79` == 1, "75-79",
          fifelse(
            `age_group80-84` == 1, "80-84",
            fifelse(
              `age_group85-89` == 1, "85-89",
              fifelse(
                `age_group90-94` == 1, "90-94",
                fifelse(
                  `age_group95-99` == 1, "95-99",
                  fifelse(
                    `age_group100-104` == 1, "100-104",
                    fifelse(
                      `age_group105-109` == 1, "105-109",
                      "65-69"
                    )
                  )
                )
              )
            )
          )
        )
      )]

      # Assign mm_cluster
      final_sim_probs_dt[, mm_cluster := fifelse(
        Neuropsychiatric == 1, "Neuropsychiatric",
        fifelse(
          Complex == 1, "Complex",
          fifelse(
            Eye == 1, "Eye",
            fifelse(
              MSK == 1, "MSK",
              fifelse(
                Metabolic == 1, "Metabolic",
                fifelse(
                  Cardiovascular == 1, "Cardiovascular",
                 "Unspecified"
                )
              )
            )
          )
        )
      )]


      final_sim_probs_dt[, colnames(design_matrix) := NULL]

      setnames(final_sim_probs_dt, "age_group", "agegrp")
      final_sim_probs_dt[, mm_cluster := self$mm_cluster_rename[mm_cluster]]

      for (mc_counter in 1:self$mc) {
        data.table::fwrite(final_sim_probs_dt[mc == mc_counter], sprintf("./outputs/transitions/basecase/%s_simulated_probs.csv", mc_counter))
      }

      message(paste("Transition Probabilities Saved to Disk @", "./outputs/transitions/basecase/"))
    },


    #' Run Markov Transition Simulation
    #'
    #' @description
    #' Executes the Markov simulation model for all patients across all defined simulation years
    #' and Monte Carlo iterations. Simulates annual transitions between health states.
    #'
    #' @details
    #' For each Monte Carlo iteration:
    #' \itemize{
    #'   \item Initializes patient states using the initial transition matrix.
    #'   \item Loads the simulated transition probabilities for that iteration.
    #'   \item Iteratively updates each patient's state for each year based on the current state and transition probabilities.
    #'   \item Uses random sampling to select the next state based on cumulative transition probabilities.
    #'   \item Filters out terminal states (e.g., death) and handles special state corrections.
    #'   \item Writes the simulated life course with transitions to disk for that iteration.
    #' }
    #'
    #' @return None. Results are written to `./outputs/lifecourse/transitions/` as compressed `.csv.gz` files.
    #'
    #' @examples
    #' tm <- TransitionModel$new("config/transition_config.yaml")
    #' tm$run_simulation()
    #'
    #' @keywords internal
    run_simulation = function() {
      for (mc_counter in 1:self$mc) {
        lifecourse <- self$initialize_states(mc_counter)
        transitions <- fread(sprintf("./outputs/transitions/basecase/%s_simulated_probs.csv", mc_counter))

        message(paste0("Running Simulation for mc: ", mc_counter))

        for (present_year in (min(lifecourse[, year])):(max(lifecourse[, year]))) {
          # message(present_year)

          temp <- lifecourse[year == present_year]

          temp[CurrentEvent == "PIS_PHS", CurrentEvent := "PIS"]

          temp <- temp[!is.na(CurrentEvent) & !(CurrentEvent %in% c(
            "", "DTH", "AB_Death",
            "AHS_AB_Death", "AHS_Death",
            "AIS_AB_Death", "AIS_Death",
            "CurrentDeath"
          ))]
          temp[, rank := runif(.N)]

          temp <- merge(temp, transitions, by = c("agegrp", "sex", "mm_cluster", "CurrentEvent"), all.x = TRUE)

          temp[, (names(temp)[10:28]) := as.data.table(t(apply(.SD, 1, cumsum))), .SDcols = names(temp)[10:28]]

          next_states <- setdiff(colnames(temp), c("pid", "agegrp", "age", "sex", "mm_cluster", "year", "af_prvl", "rank", "CurrentEvent", "mc"))

          temp_long <- melt(temp,
            id.vars = c("pid", "agegrp", "sex", "mm_cluster", "year", "rank", "mc"),
            measure.vars = next_states, variable.name = "state", value.name = "probs"
          )

          # Sorting the data
          setorder(temp_long, mc, pid, probs)

          n_states <- temp_long[, NextState := state[findInterval(rank[1], probs) + 1], by = .(pid, agegrp, sex, mm_cluster, mc)][, .SD[1], by = .(pid, mc)]
          n_states[, year := year + 1]

          lifecourse <- merge(lifecourse, n_states[, .(pid, year, NextState)], by = c("pid", "year"), all.x = TRUE)
          lifecourse[!is.na(NextState), CurrentEvent := NextState]
          lifecourse[, NextState := NULL] # Remove extra column after merging
        }
        message(paste0("Uploading results for mc ", mc_counter, " @ ",  sprintf("./outputs/lifecourse/transitions/%s_lifecourse.csv.gz", mc_counter)))
        data.table::fwrite(lifecourse, sprintf("./outputs/lifecourse/transitions/basecase/%s_lifecourse.csv.gz", mc_counter))
      }
    },

    #' Run Scenario-Based Transition Simulation
    #'
    #' Applies a set of probability adjustments to the base transition probabilities 
    #' and simulates patient health state transitions for each Monte Carlo iteration.
    #' This function allows testing of counterfactual or intervention scenarios 
    #' by modifying specific transition probabilities (e.g., reducing AIS transitions by 10%),
    #' while ensuring that each row of transition probabilities remains normalized (sums to 1).
    #'
    #' @param prob_changes A named list of transition probability multipliers. 
    #'   The names correspond to health states (e.g., "AIS", "AB") and the values 
    #'   are numeric scaling factors (e.g., list(AIS = 0.9) reduces AIS transitions by 10%).
    #' @param scenario_name A character string indicating the name of the scenario. 
    #'   Used to label output files (e.g., "reduced_AIS_10pct").
    #'
    #' @details
    #' For each Monte Carlo iteration:
    #' \itemize{
    #'   \item Loads the basecase transition probability table.
    #'   \item Applies the specified scaling factors to the selected states.
    #'   \item Re-normalizes each row to ensure all transition probabilities sum to 1.
    #'   \item Saves the modified transition probabilities.
    #'   \item Initializes the cohort and simulates transitions year by year.
    #'   \item Uses cumulative probabilities and uniform sampling to assign next states.
    #'   \item Saves the simulated lifecourse data for the scenario.
    #' }
    #'
    #' Output files are saved to:
    #' \itemize{
    #'   \item `./outputs/transitions/scenarios/` for adjusted probabilities
    #'   \item `./outputs/lifecourse/transitions/scenarios/` for simulated lifecourses
    #' }
    #'
    #' This function does not modify the original basecase outputs.
    #'
    #' @return No value is returned. Output is written to disk.
    run_scenario = function(prob_changes, scenario_name) {
      for (mc_counter in 1:self$mc) {
        trans_probs_dt <- fread(sprintf("./outputs/transitions/basecase/%s_simulated_probs.csv", mc_counter, scenario_name))

        state_cols <- setdiff(
          names(trans_probs_dt),
          c("mc", "agegrp", "sex", "mm_cluster", "CurrentEvent") # adapt if needed
        )

        for (state in names(prob_changes)) {
          factor <- prob_changes[[state]]

          # Reduce or increase the specified state's probability
          trans_probs_dt[[state]] <- trans_probs_dt[[state]] * factor
        }

        # Re-normalize so each row sums to 1
        row_sums <- rowSums(trans_probs_dt[, ..state_cols])
        for (state in state_cols) {
          trans_probs_dt[[state]] <- trans_probs_dt[[state]] / row_sums
        }

        data.table::fwrite(trans_probs_dt, sprintf("./outputs/transitions/scenarios/%s_%s_simulated_probs.csv", mc_counter, scenario_name))

        lifecourse <- self$initialize_states(mc_counter)
        transitions <- fread(sprintf("./outputs/transitions/scenarios/%s_%s_simulated_probs.csv", mc_counter, scenario_name))

        message(paste0("Running Simulation for mc: ", mc_counter))

        for (present_year in (min(lifecourse[, year])):(max(lifecourse[, year]))) {
          # message(present_year)

          temp <- lifecourse[year == present_year]

          temp[CurrentEvent == "PIS_PHS", CurrentEvent := "PIS"]

          temp <- temp[!is.na(CurrentEvent) & !(CurrentEvent %in% c(
            "", "DTH", "AB_Death",
            "AHS_AB_Death", "AHS_Death",
            "AIS_AB_Death", "AIS_Death",
            "CurrentDeath"
          ))]
          temp[, rank := runif(.N)]

          temp <- merge(temp, transitions, by = c("agegrp", "sex", "mm_cluster", "CurrentEvent"), all.x = TRUE)

          temp[, (names(temp)[10:28]) := as.data.table(t(apply(.SD, 1, cumsum))), .SDcols = names(temp)[10:28]]

          next_states <- setdiff(colnames(temp), c("pid", "agegrp", "age", "sex", "mm_cluster", "year", "af_prvl", "rank", "CurrentEvent", "mc"))

          temp_long <- melt(temp,
            id.vars = c("pid", "agegrp", "sex", "mm_cluster", "year", "rank", "mc"),
            measure.vars = next_states, variable.name = "state", value.name = "probs"
          )

          # Sorting the data
          setorder(temp_long, mc, pid, probs)

          n_states <- temp_long[, NextState := state[findInterval(rank[1], probs) + 1], by = .(pid, agegrp, sex, mm_cluster, mc)][, .SD[1], by = .(pid, mc)]
          n_states[, year := year + 1]

          lifecourse <- merge(lifecourse, n_states[, .(pid, year, NextState)], by = c("pid", "year"), all.x = TRUE)
          lifecourse[!is.na(NextState), CurrentEvent := NextState]
          lifecourse[, NextState := NULL] # Remove extra column after merging
        }
        message(paste0("Uploading results for mc ", mc_counter, " @ ",  sprintf("./outputs/lifecourse/transitions/scenarios/%s_%s_lifecourse.csv.gz", mc_counter, scenario_name)))
        data.table::fwrite(lifecourse, sprintf("./outputs/lifecourse/transitions/scenarios/%s_%s_lifecourse.csv.gz", mc_counter, scenario_name))
      }
    }
  )
)
