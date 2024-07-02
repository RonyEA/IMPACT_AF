library(data.table)
library(ggplot2)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(R6)
library(parallel)


ModelBMI <- R6::R6Class(
    public = list(

        main_dir = NULL,
        plots_folder = NULL,
        data_folder = NULL,
        univar_models_folder = NULL,
        country = NULL,
        bmi_data = NULL,
        bmi_distribution_time = NULL,
        fit_10best_distr_bmi_time = NULL,

        initialize = function(country) {
            self$country <- country
            self$set_paths()
            self$get_bmi_data()
            
        },

        set_paths = function() {
            self$main_dir <- "/home/rony/projects/IMPACTaf/"
            self$plots_folder <- paste0(self$main_dir, "outputs/plots/bmi/")
            self$data_folder <- paste0(self$main_dir, "inputs/bmi_datastore/")
            self$univar_models_folder <- paste0(self$main_dir, "outputs/models/univariate_bmi_dist/", self$country, "/")
        },

        get_bmi_data = function() {
            self$bmi_data <- private$read_data()
        },

        get_density_plot = function() {
            path = paste0(self$plots_folder, self$country, "_bmi_density.qs") 
            if (file.exists(path)) {
                density_plot <- qread(path)
            } else {
                private$bmi_density_plot()
                density_plot <- qread(path)
            }
            return(density_plot)
        },

        get_marg_distr = function() {
            path = paste0(self$univar_models_folder, self$country, "_bmi_marg_distr.qs")
            if (file.exists(path)) {
                marg_distr <- qread(path)
            } else {
                private$bmi_distribution()
                paste0("Time taken by get_marg_distr: ", self$bmi_distribution_time)
                marg_distr <- qread(path)
            }
            return(marg_distr)
        },

        get_best_bmi_marg_distr_fitted = function() {
            path = paste0(self$univar_models_folder, self$country, "_best_bmi_marg_distr_fitted.qs")
            if (file.exists(path)) {
                models <- qread(path)
            } else {
                private$best_bmi_marg_distr_fit()
                paste0("Time taken to fit 10 best marg distributions: ", self$fit_10best_distr_bmi_time)
                models <- qread(path)
            }
            return(models)
        },


    ),

    private = list(

        read_data = function() {
            path <- paste0(self$data_folder,"bmi_", self$country, ".fst")
            data <- read_fst(path, as.data.table = TRUE)

            data <- data[, `:=` (sex = factor(sex), smk = factor(smk))]

            return(data)
        },

        bmi_density_plot = function() {
            bmi_density <- ggplot(ds, aes(x = bmi)) +
                geom_density() +
                labs(x = "BMI", y = "Density", title = "Density plot of BMI - Italy") +
                theme_minimal() + # Apply a minimal theme
                theme(
                    plot.title = element_text(face = "bold", size = 16),
                    axis.title = element_text(face = "bold", size = 14),
                    axis.text = element_text(size = 12),
                    axis.line = element_line(color = "#0e0909"),
                )
            qsave(bmi_density, paste0(self$plots_folder, self$country, "_bmi_density.qs"))
        },

        bmi_distribution = function() {
            start_time <- Sys.time()
            marg_distr <- fitDist(bmi, k = 2, type = "realAll", try.gamlss = TRUE, extra = NULL, data = self$bmi_data, trace = TRUE)
            qsave(marg_distr, paste0(self$univar_models_folder, self$country, "_bmi_marg_distr.qs"), preset = "high")
            end_time <- Sys.time()
            self$bmi_distribution_time <- end_time - start_time
        },

        fit_marg_distr_model = function(distr) {
            print("Model Function is Called")
            model <- gamlss(bmi~1, family = distr, data = self$bmi_data)
            return(model)
        },

        best_bmi_marg_distr_fit = function() {
            start_time <- Sys.time()
            models <- mclapply(self$distr, fit_marg_distr_model)
            qsave(models, paste0(self$univar_models_folder, self$country, "_best_bmi_marg_distr_fitted.qs"))
            end_time <- Sys.time()
            self$fit_marg_distr_model <- end_time - start_time
        }
    )
)