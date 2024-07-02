library(ungroup)
library(data.table)
library(fst)


dt <- fread("./GBD/stroke_prvl_incd.csv")

dt[location == "United Kingdom of Great Britain and Northern Ireland", location := "United Kingdom"]
dt <- dt[age != "<1 year"]
dt <- dt[metric=="Rate"]
dt

dt[, Lower_age := fcase(
    # age == "<1 year", 0,
    age != "95+ years", as.numeric(sub("-.*", "", age)),
    age == "95+ years", 95
)]

dt[sex=="Male", sex:= "men"]
dt[sex=="Female", sex:= "women"]

# dt[, c("val", "lower", "upper") := lapply(.SD, function(x) x / 100000), .SDcols = c("val", "lower", "upper")]
dt[, c("val", "lower", "upper") := lapply(.SD, function(x) ifelse(x == 0, 1e-10, x)), .SDcols = c("val", "lower", "upper")]
setkey(dt, location, sex, Lower_age)

x <- dt[year == 2019 & location == "Italy" & sex == "Male" & measure == "Prevalence",][,Lower_age]
y <- dt[year == 2019 & location == "Italy" & sex == "Male" & measure == "Prevalence",][,val]
nlast <- 35

years <- 1991:2021
loc <- unique(dt[, location])
sex <- c("men", "women")
measure <- c("Prevalence", "Incidence")

result_list <- list()

combinations <- expand.grid(years = years, loc = loc, sex = sex, measure = measure)
combination_list <- Map(list, combinations$years, combinations$loc, combinations$sex, combinations$measure)


age_model <- function(params) {
    x <- dt[year == params[[1]] & location == params[[2]] & sex == params[[3]] & measure == params[[4]],][,Lower_age]
    y <- dt[year == params[[1]] & location == params[[2]] & sex == params[[3]] & measure == params[[4]],][,val]
    y_upper <- dt[year == params[[1]] & location == params[[2]] & sex == params[[3]] & measure == params[[4]],][,upper]
    y_lower <- dt[year == params[[1]] & location == params[[2]] & sex == params[[3]] & measure == params[[4]],][,lower]
    nlast <- 25
    m <- pclm(x, y, nlast)
    m_lower <- pclm(x, y_lower, nlast)
    m_upper <- pclm(x, y_upper, nlast)

    return(data.table(
        year = params[[1]],
        country = params[[2]],
        sex = params[[3]],
        age = c(1:length(fitted(m))),
        measure_name = params[[4]],
        mu = fitted(m),
        mu_lower = fitted(m_lower),
        mu_upper = fitted(m_upper)
    ))
}


# dt_list <- lapply(combination_list[1:3], age_model)
dt_list <- lapply(combination_list, age_model)

dt_fitted <- rbindlist(dt_list)

dt_fitted <- dt_fitted[age>19 & age<100]

dt_incd <- dt_fitted[measure_name=="Incidence"]
dt_prvl <- dt_fitted[measure_name=="Prevalence"]


for (cntry in unique(dt_fitted[, country])) {
    if (!dir.exists(paste0("./inputs/disease_burden/", cntry,"/"))){
        dir.create(paste0("./inputs/disease_burden/", cntry,"/"))
    }
    if (!file.exists(paste0("./inputs/disease_burden/", cntry,"/stroke_incd.fst"))) {
        temp_dt <- dt_incd[country==cntry]
        setkey(temp_dt, year)

        write_fst(temp_dt, paste0("./inputs/disease_burden/", cntry,"/stroke_incd.fst"))
        write_fst(temp_dt[, .(from=.I[1], to=.I[.N]), by=year], paste0("./inputs/disease_burden/", cntry, "/stroke_incd_indx.fst"))
    }
    if (!file.exists(paste0("./inputs/disease_burden/", cntry,"/stroke_prvl.fst"))) {
        temp_dt <- dt_prvl[country==cntry]
        setkey(temp_dt, year)

        write_fst(temp_dt, paste0("./inputs/disease_burden/", cntry,"/stroke_prvl.fst"))
        write_fst(temp_dt[, .(from=.I[1], to=.I[.N]), by=year], paste0("./inputs/disease_burden/", cntry, "/stroke_prvl_indx.fst"))
    }
}



for (cntry in unique(dt_fitted[, country])) {
    if (!dir.exists(paste0("./inputs/disease_burden/", cntry,"/"))){
        dir.create(paste0("./inputs/disease_burden/", cntry,"/"))
        if (!file.exists(paste0("./inputs/disease_burden/", cntry,"/stroke_incd.fst"))) {
            temp_dt <- dt_incd[country==cntry]
            setkey(temp_dt, year)

            write_fst(temp_dt, paste0("./inputs/disease_burden/", cntry,"/stroke_incd.fst"))
            write_fst(temp_dt[, .(from=.I[1], to=.I[.N]), by=year], paste0("./inputs/disease_burden/", cntry, "/stroke_incd_indx.fst"))
        }
        if (!file.exists(paste0("./inputs/disease_burden/", cntry,"/stroke_prvl.fst"))) {
            temp_dt <- dt_prvl[country==cntry]
            setkey(temp_dt, year)

            write_fst(temp_dt, paste0("./inputs/disease_burden/", cntry,"/stroke_prvl.fst"))
            write_fst(temp_dt[, .(from=.I[1], to=.I[.N]), by=year], paste0("./inputs/disease_burden/", cntry, "/stroke_prvl_indx.fst"))
        }
    }
}







