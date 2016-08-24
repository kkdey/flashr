

######  Data Preparation  #############

flash_deng_2 <- get(load("../data/flash_deng_voom_Aug152016.rda"));
saveRDS(flash_deng_2, file="../data/flash_deng.rds")
flash_deng_ex <- readRDS("../data/flash_deng.rds")
save(flash_deng_ex, file="../data/flash_deng_ex.rda")
data("flash_deng_ex")
devtools::use_data(flash_deng_ex, overwrite=TRUE)
