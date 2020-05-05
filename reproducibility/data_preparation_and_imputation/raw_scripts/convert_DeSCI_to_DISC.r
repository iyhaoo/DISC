setwd("E:/DISC/reproducibility")
utilities_path = "./source/utilities.r"
source(utilities_path)
DeSCI_imputation = readh5_imputation("./data/MELANOMA/DISC.hdf5", with_outliers = T)
save_h5("./data/MELANOMA/DISC.loom", t(DeSCI_imputation))
DeSCI_imputation = readh5_imputation("./data/MELANOMA/ds_0.5/DISC.hdf5", with_outliers = T)
save_h5("./data/MELANOMA/ds_0.5/DISC.loom", t(DeSCI_imputation))
DeSCI_imputation = readh5_imputation("./data/PBMC/ds_0.5/DISC.hdf5", with_outliers = T)
save_h5("./data/PBMC/ds_0.5/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_2/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_2/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/DISC.loom", t(DeSCI_imputation))

DeSCI_imputation = readh5_imputation("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/imputation.hdf5", with_outliers = T)
save_h5("/home/yuanhao/data/fn/melanoma/ds/downsampling_first_repeat_1/DeSCI_2.7.4.33/ds_0.5/results/epoch_1170/DISC.loom", t(DeSCI_imputation))
