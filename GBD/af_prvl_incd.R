library(ungroup)
library(data.table)
library(fst)


dt <- fread("./GBD/af_prvl_incd.csv")
dt[location == "United Kingdom of Great Britain and Northern Ireland", location := "United Kingdom"]

dt <- dt[age != "<1 year"]
dt

dt[, age := sub(" years", "", age)]
dt[age=="95+", age := "95-100"]

dt[, measure_name := measure]
dt[, measure:=NULL]


dt[, cause_name := "af"]
dt[, cause := NULL]

dt[, age_group := age]
dt[, age:=NULL]

# dt[, c("lower_age", "upper_age") := tstrsplit(age_group, "-", fixed = TRUE)]

dt[sex=="Male", sex:= "men"]
dt[sex=="Female", sex:= "women"]

years <- 1990:2021
loc <- unique(dt[, location])
sex <- c("men", "women")
age <- 2:100
measure <- c("Prevalence", "Incidence")

result_list <- list()

combinations <- data.table(expand.grid(year = years, location = loc, sex = sex, age = age, measure_name = measure))
combinations[, age_group := ""]

for (age_grp in unique(dt[, age_group])) {
    lower_age <- as.integer(strsplit(age_grp, "-", fixed = T)[[1]][1])
    upper_age <- as.integer(strsplit(age_grp, "-", fixed = T)[[1]][2])

    combinations[age >= lower_age & age <= upper_age, age_group := age_grp]
}

dt <- merge(combinations, dt, by = c("measure_name", "location", "sex", "year", "age_group"))









# dt[, Lower_age := fcase(
#     # age == "<1 year", 0,
#     age != "95+ years", as.numeric(sub("-.*", "", age)),
#     age == "95+ years", 95
# )]

# dt[, Upper_age := fcase(
#     # age == "<1 year", 0,
#     age != "95+ years", as.numeric(sub("-.*", "", age)),
#     age == "95+ years", 95
# )]


# dt[, c("val", "lower", "upper") := lapply(.SD, function(x) x / 100000), .SDcols = c("val", "lower", "upper")]
dt[, c("val", "lower", "upper") := lapply(.SD, function(x) ifelse(x == 0, 1e-10, x)), .SDcols = c("val", "lower", "upper")]
setkey(dt, location, sex, age)

# x <- dt[year == 2019 & location == "Italy" & sex == "Male" & measure == "Prevalence",][,Lower_age]
# y <- dt[year == 2019 & location == "Italy" & sex == "Male" & measure == "Prevalence",][,val]
# nlast <- 35

# years <- 1990:2021
# loc <- unique(dt[, location])
# sex <- c("men", "women")
# age <- 2:100
# measure <- c("Prevalence", "Incidence")

# result_list <- list()

# combinations <- expand.grid(years = years, loc = loc, sex = sex, age = age, measure = measure)

# combination_list <- Map(list, combinations$years, combinations$loc, combinations$sex, combinations$measure)


# gbd[year==2021 & sex == "men", sum(val)/100000]
# af_incd[year==2021 & sex == "men", sum(mu)/100000]

# af_incd <- read_fst("/home/rony/projects/IMPACT_AF/inputs/disease_burden/Italy/af_incd.fst", as.data.table = TRUE)
# head(af_incd)

# af_incd[year==2021 & sex == "men"]

# age_model <- function(params) {
#     x <- dt[year == params[[1]] & location == params[[2]] & sex == params[[3]] & measure == params[[4]],][,Lower_age]
#     y <- dt[year == params[[1]] & location == params[[2]] & sex == params[[3]] & measure == params[[4]],][,val]

#     y_upper <- dt[year == params[[1]] & location == params[[2]] & sex == params[[3]] & measure == params[[4]],][,upper]
#     y_lower <- dt[year == params[[1]] & location == params[[2]] & sex == params[[3]] & measure == params[[4]],][,lower]

#     nlast <- 5

#     m <- pclm(x, y, nlast)
#     m_lower <- pclm(x, y_lower, nlast)
#     m_upper <- pclm(x, y_upper, nlast)

#     return(data.table(
#         year = params[[1]],
#         country = params[[2]],
#         sex = params[[3]],
#         age = c(1:length(fitted(m))),
#         measure_name = params[[4]],
#         mu = fitted(m),
#         mu_lower = fitted(m_lower),
#         mu_upper = fitted(m_upper)
#     ))
# }


# # dt_list <- lapply(combination_list[1:3], age_model)
# dt_list <- lapply(combination_list, age_model)

# dt_fitted <- rbindlist(dt_list)

dt <- dt[age>19 & age<100]

dt[, country := location]
dt[, location := NULL]

dt[, mu := val]
dt[, mu_lower := lower]
dt[, mu_upper := upper]

dt[, val := NULL]
dt[, lower := NULL]
dt[, upper := NULL]

dt[, age_group := NULL]
dt[, metric := NULL]

dt_incd <- dt[measure_name=="Incidence"]
dt_prvl <- dt[measure_name=="Prevalence"]




cntry <- unique(dt[, country])

# for (cntry in unique(dt_fitted[, country])) {
#     if (!dir.exists(paste0("./inputs/disease_burden/", cntry,"/"))){
#         dir.create(paste0("./inputs/disease_burden/", cntry,"/"))
#         if (file.exists(paste0("./inputs/disease_burden/", cntry,"/af_incd.fst"))) {
#             temp_dt <- dt_incd[country==cntry]
#             setkey(temp_dt, year, sex, age)

#             write_fst(temp_dt, paste0("./inputs/disease_burden/", cntry,"/af_incd.fst"))
#             write_fst(temp_dt[, .(from=.I[1], to=.I[.N]), by=year], paste0("./inputs/disease_burden/", cntry, "/af_incd_indx.fst"))
#         }
#         if (file.exists(paste0("./inputs/disease_burden/", cntry,"/af_prvl.fst"))) {
#             temp_dt <- dt_incd[country==cntry]
#             setkey(temp_dt, year, sex, age)

#             write_fst(temp_dt, paste0("./inputs/disease_burden/", cntry,"/af_prvl.fst"))
#             write_fst(temp_dt[, .(from=.I[1], to=.I[.N]), by=year], paste0("./inputs/disease_burden/", cntry, "/af_prvl_indx.fst"))
#         }
#     }
# }



for (cntry in unique(dt[, country])) {
    temp_dt <- dt_incd[country==cntry]
    setkey(temp_dt, year, age, sex)

    write_fst(temp_dt, paste0("./inputs/disease_burden/", cntry,"/af_incd.fst"))
    write_fst(temp_dt[, .(from=.I[1], to=.I[.N]), by=year], paste0("./inputs/disease_burden/", cntry, "/af_incd_indx.fst"))
    
    temp_dt <- dt_prvl[country==cntry]
    setkey(temp_dt, year, age, sex)

    write_fst(temp_dt, paste0("./inputs/disease_burden/", cntry,"/af_prvl.fst"))
    write_fst(temp_dt[, .(from=.I[1], to=.I[.N]), by=year], paste0("./inputs/disease_burden/", cntry, "/af_prvl_indx.fst"))
}




