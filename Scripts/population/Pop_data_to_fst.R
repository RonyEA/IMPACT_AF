library(readxl)
library(eurostat)
library(data.table)
library(fst)

countries <- eu_countries$name
countries <- c(countries, "United Kingdom")

###################################### Population Estimates ##################################################

## Population Data for Male
pop_m <- read_excel("../Data/Population_Data/Pop_Male.xlsx", sheet = 1, skip = 16)
pop_m <- data.table(pop_m)

cols <- colnames(pop_m)
cols_to_sel <- cols[c(3, 11:112)]

pop_m <- pop_m[, ..cols_to_sel]
setnames(pop_m, cols[3], "Country")


# Subsetting for european countries
pop_m_eu <- pop_m[Country %in% countries]
pop_m_eu <- melt(pop_m_eu, id.vars= c("Country", "Year"), variable.name = "Age", value.name = "Population")
pop_m_eu

## Population Data for Female
pop_f <- read_excel("../Data/Population_Data/Pop_Female.xlsx", sheet = 1, skip = 16)
# pop_f <- read_excel("./Pop_Female.xlsx", sheet = 1, skip = 16)
pop_f <- data.table(pop_f)
pop_f <- pop_f[, ..cols_to_sel]
setnames(pop_f, cols[3], "Country")

# Subsetting for european countries
pop_f_eu <- pop_f[Country %in% countries]
pop_f_eu <- melt(pop_f_eu, id.vars= c("Country", "Year"), variable.name = "Age", value.name = "Population")

# Combine Male and Female Data
pop_m_eu <- pop_m_eu[, Population := as.numeric(Population) * 1000]
pop_m_eu <- pop_m_eu[, Sex := "men"]

pop_f_eu <- pop_f_eu[, Population := as.numeric(Population) * 1000]
pop_f_eu <- pop_f_eu[, Sex := "women"]

pop_eu <- rbindlist(list(pop_m_eu, pop_f_eu))
colnames(pop_eu) <- c("country", "year", "age", "pops", "sex")


pop_eu[, age:=as.integer(age)]
pop_eu[, age:=age-1]

setkey(pop_eu, year)

write_fst(pop_eu, "./inputs/pop_observed/pop_observed_eu.fst")


###################################### Population Projections ##################################################

## Population Data for Male
pop_m <- read_excel("../Data/Population_Data/Pop_proj_Male.xlsx", sheet = 3, skip = 16)
pop_m <- data.table(pop_m)

cols <- colnames(pop_m)
cols_to_sel <- cols[c(3, 11:112)]

pop_m <- pop_m[, ..cols_to_sel]
setnames(pop_m, cols[3], "Country")


# Subsetting for european countries
pop_m_eu <- pop_m[Country %in% countries]
pop_m_eu <- melt(pop_m_eu, id.vars= c("Country", "Year"), variable.name = "Age", value.name = "Population")
pop_m_eu

## Population Data for Female
pop_f <- read_excel("../Data/Population_Data/Pop_proj_Female.xlsx", sheet = 1, skip = 16)
# pop_f <- read_excel("./Pop_Female.xlsx", sheet = 1, skip = 16)
pop_f <- data.table(pop_f)
pop_f <- pop_f[, ..cols_to_sel]
setnames(pop_f, cols[3], "Country")

# Subsetting for european countries
pop_f_eu <- pop_f[Country %in% countries]
pop_f_eu <- melt(pop_f_eu, id.vars= c("Country", "Year"), variable.name = "Age", value.name = "Population")

# Combine Male and Female Data
pop_m_eu <- pop_m_eu[, Population := as.numeric(Population) * 1000]
pop_m_eu <- pop_m_eu[, Sex := "men"]

pop_f_eu <- pop_f_eu[, Population := as.numeric(Population) * 1000]
pop_f_eu <- pop_f_eu[, Sex := "women"]

pop_eu <- rbindlist(list(pop_m_eu, pop_f_eu))
colnames(pop_eu) <- c("country", "year", "age", "pops", "sex")


pop_eu[, age:=as.integer(age)]
pop_eu[, age:=age-1]

setkey(pop_eu, year)

write_fst(pop_eu, "./inputs/pop_projections/pop_projection_eu.fst")

pop_obs <- read_fst("./inputs/pop_observed/pop_observed_eu.fst", as.data.table = TRUE)
pop_proj <- read_fst("./inputs/pop_projections/pop_projection_eu.fst", as.data.table = TRUE)

pop_combined <- rbindlist(list(pop_obs, pop_proj))
setkey(pop_combined, year)

write_fst(pop_combined, "./inputs/pop_projections/pop_combined_eu.fst")