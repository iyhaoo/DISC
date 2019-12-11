source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/generic_functions.r")
method_names = c("Observed", "DeSCI", "SAVER", "MAGIC", "DCA", "scScope", "scVI")
replicates = paste0("downsampling_first_repeat_", seq(5))
#ds_mode = c("ds_0.3", "ds_0.5", "ds_0.7")
ds_mode = c("ds_0.3")
raw_path = "/home/yuanhao/data/fn/pbmc3k/pbmc3k_filtered.loom"
output_dir = paste0(stringi::stri_trim_right(raw_path, "[/]"), "CMD_vst_ds")
dir.create(output_dir, showWarnings = F, recursive = T)
gene_bc_mat = readh5_loom(raw_path)
vst_file = paste0(output_dir, "/vst_gene.tsv")
if(file.exists(vst_file)){
  hvg_info = read.table(vst_file)
  print("load vst_file")
}else{
  hvg_info = FindVariableFeatures_vst_by_genes(gene_bc_mat)
  hvg_info = hvg_info[order(hvg_info$variance.standardized, decreasing = T), ]
  write.table(hvg_info, paste0(output_dir, "/vst_gene.tsv"), sep = "\t", quote = F, row.names = T, col.names = T)
}
for(this_mode in ds_mode){
  cmd_mat = matrix(nrow = length(replicates), ncol = length(method_names), dimnames = list(replicates, method_names))
  raw_mat = gene_bc_mat
  compare_gene = rownames(raw_mat)
  for(ii in replicates){
    ds_path = paste0("/home/yuanhao/data/fn/pbmc3k/ds/", ii, "/imputation/pbmc3k_filtered_", this_mode, "_mc_10_mce_1.loom")
    compare_gene = intersect(compare_gene, get_loom_gene(ds_path))
  }
  use_genes = rownames(hvg_info[rownames(hvg_info) %in% compare_gene, ])[1:300]
  cat("use_genes: ", length(use_genes), "\n")
  raw_mat = raw_mat[use_genes, ]
  for(ii in replicates){
    cor_all = list()
    for(method in c("Raw", method_names)){
      cor_all[[method]] = list()
      cor_all[[method]] = matrix(nrow = length(use_genes), ncol = length(use_genes), dimnames = list(use_genes, use_genes))
    }
    for(method in c("Raw", method_names)){
      switch(method,
             Raw = {
               impute_result = raw_mat
             },
             Observed={
               impute_result = readh5_loom(paste0("/home/yuanhao/data/fn/pbmc3k/ds/", ii, "/imputation/pbmc3k_filtered_", this_mode, "_mc_10_mce_1.loom"))[use_genes, ]
             },
             "DeSCI"={
               impute_result = readh5_imputation(get_optimal_point33(paste0("/home/yuanhao/data/fn/pbmc3k/ds/", ii, "/DeSCI_2.7.4.33/", this_mode, "/log.txt")), with_outliers = T)[use_genes, ]
             },{
               impute_result = readh5_imputation(paste0("/home/yuanhao/data/fn/pbmc3k/ds/", ii, "/imputation/pbmc3k_filtered_", this_mode, "_", method, "_mc_10_mce_1.hdf5"))[use_genes, ]
             })
      this_mat = delect_lt0.5(impute_result)
      no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      clusterExport(cl, varlist = c("this_mat", "method"))
      return_list = parLapply(cl, seq(nrow(this_mat) - 1), function(x){
        return_vector = rep(NA, nrow(this_mat) - x)
        ii_express = this_mat[x, ]
        ii_mask = ii_express > 0
        for(jj in (x + 1):nrow(this_mat)){
          jj_express = this_mat[jj, ]
          if(method %in% c("")){
            express_mask = rep(T, length(ii_mask))
          }else{
            express_mask = ii_mask | jj_express > 0
          }
          if(sum(express_mask) > 0){
            this_corr = cor(ii_express[express_mask], jj_express[express_mask], use = "pairwise.complete.obs")
            if(is.na(this_corr)){
              return_vector[jj - x] = 0
            }else{
              return_vector[jj - x] = this_corr
            }
          }else{
            return_vector[jj - x] = 0
          }
        }
        return(list("return_vector" = return_vector))
      })
      stopCluster(cl)
      cor_all[[method]][1, 1] = 1
      for(jj in 1:length(return_list)){
        cor_all[[method]][jj, (jj + 1): nrow(this_mat)] = return_list[[jj]][["return_vector"]]
        cor_all[[method]][(jj + 1): nrow(this_mat), jj] = return_list[[jj]][["return_vector"]]
        cor_all[[method]][(jj + 1), (jj + 1)] = 1
      }
      print(method)
    }
    saveRDS(cor_all, paste0(output_dir, "/", ds_mode, "_", ii, "cor_all.rds"))
    for(method in method_names){
      cmd_mat[ii, method] = calc_cmd(cor_all[["Raw"]], cor_all[[method]])
    }
  }
  print(cmd_mat)
  saveRDS(cmd_mat, paste0(output_dir, "/", ds_mode, "_cmd.rds"))
  saveRDS(colMeans(cmd_mat), paste0(output_dir, "/", ds_mode, "_cmd_result.rds"))
}






