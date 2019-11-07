source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/generic_functions.r")
method_names = c("Raw", "DeSCI", "SAVER", "MAGIC", "DCA", "scScope", "scVI")
replicates = paste0("downsampling_first_repeat_", seq(5))
#ds_mode = c("ds_0.3", "ds_0.5", "ds_0.7")
ds_mode = c("ds_0.5")
raw_path = "/home/yuanhao/data/fn/pbmc3k/pbmc3k_filtered.loom"
output_dir = paste0(stringi::stri_trim_right(raw_path, "[/]"), "MAE_ds")
dir.create(output_dir, showWarnings = F, recursive = T)
compare_gene = get_loom_gene(raw_path)
for(this_mode in ds_mode){
  for(ii in replicates){
    ds_path = paste0("/home/yuanhao/data/fn/pbmc3k/ds/", ii, "/imputation/pbmc3k_filtered_", this_mode, "_mc_10_mce_1.loom")
    compare_gene = intersect(compare_gene, get_loom_gene(ds_path))
  }
}
cat("Compare gene number: ", length(compare_gene), "\n")
raw_mat = get_gene_bc_mat(raw_path)[compare_gene, ]
mae_gt0_list = list()
mae_eq0_list = list()
mae_all_list = list()
for(this_mode in ds_mode){
  switch(this_mode, 
         ds_0.3 = {scale_factor = 1 / 0.3},
         ds_0.5 = {scale_factor = 1 / 0.5},
         ds_0.7 = {scale_factor = 1 / 0.7}
         )
  mae_gt0_list[[this_mode]] = matrix(nrow = length(method_names), ncol = length(replicates), dimnames = list(method_names, replicates))
  mae_eq0_list[[this_mode]] = matrix(nrow = length(method_names), ncol = length(replicates), dimnames = list(method_names, replicates))
  mae_all_list[[this_mode]] = matrix(nrow = length(method_names), ncol = length(replicates), dimnames = list(method_names, replicates))
  for(ii in replicates){
    ds_path = paste0("/home/yuanhao/data/fn/pbmc3k/ds/", ii, "/imputation/pbmc3k_filtered_", this_mode, "_mc_10_mce_1.loom")
    ds_mat = readh5_loom(ds_path)[compare_gene, ]
    for(method in method_names){
      switch(method,
             Raw={
               impute_result = ds_mat
             },
             "DeSCI"={
               impute_result = readh5_imputation(get_optimal_point1(paste0("/home/yuanhao/data/fn/pbmc3k/ds/", ii, "/DeSCI_2.7.4.33/", this_mode, "/log.txt")), with_outliers = T)[compare_gene, ]
             },{
               impute_result = readh5_imputation(paste0("/home/yuanhao/data/fn/pbmc3k/ds/", ii, "/imputation/pbmc3k_filtered_", this_mode, "_", method, "_mc_10_mce_1.hdf5"))[compare_gene, ]
             })
      eq0 = sum(raw_mat > 0 & ds_mat == 0)
      gt0 = sum(raw_mat > 0 & ds_mat > 0)
      sae_gt0 = sum(sapply(rownames(impute_result), function(x){
        expressed_mask = raw_mat[x, ] > 0 & ds_mat[x, ] > 0
        return(sum(abs(raw_mat[x, expressed_mask] - (impute_result[x, expressed_mask] * scale_factor))))
      }))
      sae_eq0 = sum(sapply(rownames(impute_result), function(x){
        expressed_mask = raw_mat[x, ] > 0 & ds_mat[x, ] == 0
        return(sum(abs(raw_mat[x, expressed_mask] - (impute_result[x, expressed_mask] * scale_factor))))
      }))
      mae_gt0 = sae_gt0 / gt0
      mae_eq0 = sae_eq0 / eq0
      mae_all = (sae_gt0 + sae_eq0) / (gt0 + eq0)
      mae_gt0_list[[this_mode]][method, ii] = mae_gt0
      mae_eq0_list[[this_mode]][method, ii] = mae_eq0
      mae_all_list[[this_mode]][method, ii] = mae_all
    }
  }
}
save(mae_gt0_list, mae_eq0_list, mae_all_list, file = paste0(output_dir, "/mae_tmp.rdata"))
mae_gt0_mat = matrix(nrow = length(method_names), ncol = length(ds_mode), dimnames = list(method_names, ds_mode))
mae_gt0_ste_mat = matrix(nrow = length(method_names), ncol = length(ds_mode), dimnames = list(method_names, paste0(ds_mode, "_ste")))
for(this_mode in ds_mode){
  mae_gt0_mat[, this_mode] = Matrix::rowMeans(mae_gt0_list[[this_mode]])
  mae_gt0_ste_mat[, paste0(this_mode, "_ste")] = apply(mae_gt0_list[[this_mode]], 1, function(x) sqrt(var(x)/length(x)))
}
write.table(cbind(mae_gt0_mat, mae_gt0_ste_mat), paste0(output_dir, "/mae_gt0_with_outliers.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)


mae_eq0_mat = matrix(nrow = length(method_names), ncol = length(ds_mode), dimnames = list(method_names, ds_mode))
mae_eq0_ste_mat = matrix(nrow = length(method_names), ncol = length(ds_mode), dimnames = list(method_names, paste0(ds_mode, "_ste")))
for(this_mode in ds_mode){
  mae_eq0_mat[, this_mode] = Matrix::rowMeans(mae_eq0_list[[this_mode]])
  mae_eq0_ste_mat[, paste0(this_mode, "_ste")] = apply(mae_eq0_list[[this_mode]], 1, function(x) sqrt(var(x)/length(x)))
}
write.table(cbind(mae_eq0_mat, mae_eq0_ste_mat), paste0(output_dir, "/mae_eq0_with_outliers.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)



mae_all_mat = matrix(nrow = length(method_names), ncol = length(ds_mode), dimnames = list(method_names, ds_mode))
mae_all_ste_mat = matrix(nrow = length(method_names), ncol = length(ds_mode), dimnames = list(method_names, paste0(ds_mode, "_ste")))
for(this_mode in ds_mode){
  mae_all_mat[, this_mode] = Matrix::rowMeans(mae_all_list[[this_mode]])
  mae_all_ste_mat[, paste0(this_mode, "_ste")] = apply(mae_all_list[[this_mode]], 1, function(x) sqrt(var(x)/length(x)))
}
write.table(cbind(mae_all_mat, mae_all_ste_mat), paste0(output_dir, "/mae_all_with_outliers.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)


