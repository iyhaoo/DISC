source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/generic_functions.r")
loom_path = "/home/yuanhao/data/fn/neuron1.3m/1M_neurons_filt_gene_500_5000.loom"
gene_bc_mat = t(h5read(loom_path, "matrix"))
colnames(gene_bc_mat) = h5read(loom_path, "col_attrs/CellID")
row.names(gene_bc_mat) = h5read(loom_path, "row_attrs/Gene")
hvf.info = FindVariableFeatures_vst_by_genes(gene_bc_mat)
hvg = rownames(hvf.info)[order(hvf.info$variance.standardized, decreasing = T)]
write.table(hvf.info[order(hvf.info$variance.standardized, decreasing = T), ], "/home/yuanhao/data/fn/neuron1.3m/1M_neurons_filt_gene_500_5000_unique_rename_hvg_1.tsv", sep = "\t", quote = F, row.names = T, col.names = T)





