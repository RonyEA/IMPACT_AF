library(data.table)
library(fst)

mort_df <- fread("/home/rony/projects/IMPACTaf/inputs/mortalityIT.txt")

mort_df[Age=="110+", Age:="110"]
mort_df[, Age:=as.numeric(Age)]
mort_df[, Female:=as.numeric(Female)]
mort_df[, Male:=as.numeric(Male)]
mort_df[, Total:=as.numeric(Total)]

mort_df <- mort_df[Age <= 100]
write_fst(mort_df, "/home/rony/projects/IMPACTaf/inputs/mortality/IT_mort.fst")
