source("./global.R")
source("./Rpackage/IMPACTaf_model_pkg/R/Design_class.R")

library(fst)

# Check data used.
pop_eu <- read_fst("../Data/Final_Pop_Data/pop_proj_eu.fst", as.data.table = TRUE)
head(pop_eu)

pop_jp <- read_fst("./inputs_jp/pop_projections/projected_population_japan.fst",as.data.table = TRUE)
head(pop_jp)




design <- Design$new("./inputs/sim_design.yaml")

sim <- SynthPop$new(0L, design)




sim$write_synthpop(1:500)
sim$delete_synthpop(NULL)
ll <- sim$gen_synthpop_demog(design)
sp <- SynthPop$new(1L, design)
