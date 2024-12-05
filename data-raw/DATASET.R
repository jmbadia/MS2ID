##CODE to prepare `DATASET` for calculate adducts
if(F){
#A. Code to obtain 'adducts' dataset------------------
#load adducts table (enviPat)
adducts <- readRDS("data/adductsTable.rds")
# change positive or negative polarity with usual terms (1 or 0)
adducts[adducts$Ion_mode=="positive","Ion_mode"] <- 1
adducts[adducts$Ion_mode=="negative","Ion_mode"] <- 0

#desar data nomes per us intern del package (e,g, taula adducts)
# usethis::use_data(adducts, overwrite = TRUE, internal = TRUE)

#save as rda, default format loaded automatically with the package
save(adducts, file = "data/adducts.rda")
}
