setwd("E:/DISC/reproducibility")
utilities_path = "./source/utilities.r"
source(utilities_path)
DeSCI_imputation = readh5_imputation("./data/MELANOMA/DISC.hdf5", with_outliers = T)
save_h5("./data/MELANOMA/DISC.loom", t(DeSCI_imputation))
DeSCI_imputation = readh5_imputation("./data/MELANOMA/ds_0.5/DISC.hdf5", with_outliers = T)
save_h5("./data/MELANOMA/ds_0.5/DISC.loom", t(DeSCI_imputation))
DeSCI_imputation = readh5_imputation("./data/PBMC/ds_0.5/DISC.hdf5", with_outliers = T)
save_h5("./data/PBMC/ds_0.5/DISC.loom", t(DeSCI_imputation))


