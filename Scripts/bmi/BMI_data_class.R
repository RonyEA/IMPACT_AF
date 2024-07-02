library(R6)
library(countrycode)
library(fst)
library(data.table)

GetBMIData <- R6::R6Class(
    public = list(
        
        w1_countries = NULL,
        w2_countries = NULL,
        w3_countries = NULL,

        ehis_countries = NULL,

        data = NULL,

        initialize = function() {
           self$data <- private$combine_all_waves() 
           self$ehis_countries <- private$update_countries()
        },

        get_w1_bmi = function() {
            return(private$prepare_wave1_bmi_data())
        },

        get_w2_bmi = function() {
            return(private$prepare_wave2_bmi_data())
        },

        get_w3_bmi = function() {
            return(private$prepare_wave3_bmi_data())
        },

        update_datastore = function(path) {
            dt <- self$data

            for (c in self$ehis_countries) {
                temp_dt <- dt[country == c]
                temp_dt <- temp_dt[order(year, age, sex)]
                file <- paste0(path, "bmi_", c, ".fst")
                write_fst(temp_dt, file)

                print(paste0(file, "-> Created"))
            }
        },

        show_missing = function(dt, country, is_short = FALSE) {

            if(is_short){
                country <- countrycode(country, origin = "country.name", destination = "iso2c")
            }

            cols <- c("age", "year", "sex", "smk", "ht", "wt") 
            percentage_missing <- dt[, lapply(.SD, function(x) round(mean(is.na(x)) * 100, 2)), .SDcols = cols, by = country]
            print(percentage_missing)
        },

        replace_negatives_with_NA = function(col) {
            col[col < 0] <- NA
            return(col)
        }

    ),

    private = list(
        prepare_wave1_bmi_data = function() {
            path_w1 <- "/mnt/rony/UoL/AFFIRMO/EHIS_microdata/EHIS wave 1/Data EHIS"
            dt_w1 <- fread(paste0(path_w1, "/EHIS1.csv"))

            #w1_cols <- c("AGEGROUP", "YEAR", "SEX", "IP01", "HS04C", "HS04D", "HS04E", "HS04F", "HS04K", "AL01", "SK04", "BMI01", "BMI02")
            w1_cols <- c("AGEGROUP", "YEAR", "SEX", "IP01", "SK04", "BMI01", "BMI02")

            #w1_new_cols <- c("age", "year", "sex", "country", "mi", "chd", "hpb", "strk", "diab", "alc" , "smk", "ht", "wt") 
            w1_new_cols <- c("age", "year", "sex", "country", "smk", "ht", "wt") 

            dt_w1 <- dt_w1[, ..w1_cols]
            colnames(dt_w1) <- w1_new_cols

            dt_w1 <- dt_w1[year > 2006 & year <= 2009]
            dt_w1 <- dt_w1[, names(dt_w1) := lapply(.SD, self$replace_negatives_with_NA), .SDcols = names(dt_w1)]
            dt_w1 <- na.omit(dt_w1)

            self$w1_countries <- unique(dt_w1[, country])
            return(dt_w1)
        },

        prepare_wave2_bmi_data = function() {
            path <- "/mnt/rony/UoL/AFFIRMO/EHIS_microdata/EHIS wave 2/"

            files <- list.files(path = path, pattern = ".csv")
            files <- paste(path, files, sep="")
            # length(files)

            data_list <- lapply(files, function(file) {
                data <- fread(file, select = NULL, fill = TRUE)
                setnames(data, toupper(names(data)))
                return(data)
            })

            dt_w2 <- rbindlist(data_list, fill=TRUE)

            #w2_cols <- c("AGE", "REFYEAR", "SEX", "COUNTRY", "CD1C", "CD1D", "CD1E", "CD1F", "CD1J", "AL1", "SK1", "BM1", "BM2")
            w2_cols <- c("AGE", "REFYEAR", "SEX", "COUNTRY", "SK1", "BM1", "BM2")
            #w2_new_cols <- c("age", "year", "sex", "country", "mi", "chd", "hpb", "strk", "diab", "alc" , "smk", "ht", "wt") 
            w2_new_cols <- c("age", "year", "sex", "country", "smk", "ht", "wt") 

            dt_w2 <- dt_w2[, ..w2_cols]
            colnames(dt_w2) <- w2_new_cols

            dt_w2 <- dt_w2[, names(dt_w2) := lapply(.SD, self$replace_negatives_with_NA), .SDcols = names(dt_w2)]

            dt_w2 <- dt_w2[year >= 2013 & year <= 2015]
            dt_w2 <- na.omit(dt_w2)

            #dt_w2 <- dt_w2[,SMK := 4-SMK]
            #dt_w2 <- dt_w2[SMK > 1,SMK := 2]

            self$w2_countries <- unique(dt_w2[, country])
            return(dt_w2)
        },

        prepare_wave3_bmi_data = function() {
            path = "/mnt/rony/UoL/AFFIRMO/EHIS_microdata/EHIS wave 3/"

            files <- list.files(path = path, pattern = ".csv")
            files <- paste(path, files, sep="")
            # length(files)

            data_list <- lapply(files, function(file) {
                data <- fread(file, select = NULL, fill = TRUE)
                setnames(data, toupper(names(data)))
                return(data)
            })

            dt_w3 <- rbindlist(data_list, fill=TRUE)

            #w3_cols <- c("AGE", "REFDATE", "SEX", "COUNTRY", "CD1C", "CD1D", "CD1E", "CD1F", "CD1J", "AL1", "SK1", "BM1", "BM2")
            w3_cols <- c("AGE", "REFDATE", "SEX", "COUNTRY", "SK1", "BM1", "BM2")
            #w3_new_cols <- c("age", "year", "sex", "country", "mi", "chd", "hpb", "strk", "diab", "alc" , "smk", "ht", "wt") 
            w3_new_cols <- c("age", "year", "sex", "country", "smk", "ht", "wt") 

            dt_w3 <- dt_w3[, ..w3_cols]
            colnames(dt_w3) <- w3_new_cols

            dt_w3 <- dt_w3[, names(dt_w3) := lapply(.SD, self$replace_negatives_with_NA), .SDcols = names(dt_w3)]

            dt_w3 <- na.omit(dt_w3)

            dt_w3 <- dt_w3[, year := year %/% 100]
            dt_w3 <- dt_w3[country == "IT", year := year %/% 100]

            dt_w3 <- dt_w3[year >= 2018 & year <= 2020]

            self$w3_countries <- unique(dt_w3[, country])
            return(dt_w3)
        },

        combine_all_waves = function() {

            fin_df <- rbindlist(list(self$get_w1_bmi(), self$get_w2_bmi(), self$get_w3_bmi()))

            fin_df <- fin_df[age != "75+" & age != "80+"]
            fin_df <- fin_df[, ht := ht/100]
            fin_df <- fin_df[, bmi := wt/(ht^2)]

            fin_df_85 <- fin_df[age=="85+"]
            fin_df_lessthan85 <- fin_df[age != "85+"]

            fin_df_lessthan85[, c("lower_bound", "upper_bound") := tstrsplit(age, "-", fixed = TRUE)]
            fin_df_lessthan85[, c("lower_bound", "upper_bound") := lapply(.SD, as.numeric), .SDcols = c("lower_bound", "upper_bound")]
            fin_df_lessthan85[, age := (lower_bound + upper_bound) / 2]

            fin_df_lessthan85[, c("lower_bound", "upper_bound") := NULL]

            fin_df_85[, age:= 85]

            fin_df <- rbindlist(list(fin_df_lessthan85, fin_df_85))

            fin_df[, age := as.numeric(age)]

            fin_df[, sex := ifelse(sex==1, "men", "women")]
            fin_df[, c("ht", "wt") := NULL]

            fin_df[, year := year - 2000]
            fin_df[, smk := 4-smk]

            #col_order <- c("country", "year", "age", "sex", "mi", "chd", "hpb", "strk", "diab", "alc" , "smk", "bmi") 
            col_order <- c("country", "year", "age", "sex", "smk", "bmi") 
            setcolorder(fin_df, col_order)

            return(fin_df[])
        },

        update_countries = function() {
            countries <- Reduce(union, list(self$w1_countries, self$w2_countries, self$w3_countries))
            return(countries)
        }
    )

)

