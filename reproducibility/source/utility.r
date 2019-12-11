library(Matrix)
library(ggplot2)
library(rhdf5)
library(parallel)
library(reldist)
library(Seurat)
library(future)
library(stringi)
#  read data
cal_RMSD = function(pd_array, window_size){
  apply(sapply(1:(length(pd_array) - window_size) + window_size - 1, function(x){
    pd_array[(x - window_size + 1):x]
  }),2 , function(y){
    sqrt(sum((y - mean(y))^2) / window_size)
  })
}
get_optimal_point = function(log_path){
  training_log = readLines(log_path)
  cell_number = as.numeric(stri_extract(training_log[5], regex = "[0-9]+"))
  feature_loss = as.numeric(stri_extract(training_log, regex = "(?<=feature_loss: )[0-9.]+"))
  feature_loss = feature_loss[!is.na(feature_loss)]
  window_size = max(c(ceiling(250000 / cell_number), 5))
  feature_RMSD = cal_RMSD(feature_loss, window_size)
  names(feature_RMSD) = seq(length(feature_RMSD)) + window_size - 1
  feature_convergence_estimate = feature_RMSD * sqrt(seq(length(feature_RMSD)) + window_size - 1)
  estimation = feature_RMSD
  optimal_point = round(as.numeric(names(which.min(feature_convergence_estimate))) / window_size) * window_size
  if(optimal_point > (length(feature_loss) - 1)){
    optimal_point = (round(optimal_point / window_size) - 1) * window_size
  }
  return(paste0(stri_trim_right(log_path, pattern = "[/]"), "results/epoch_", optimal_point, "/imputation.hdf5"))
}
get_optimal_point1 = function(log_path){
  training_log = readLines(log_path)
  cell_number = as.numeric(stri_extract(training_log[5], regex = "[0-9]+"))
  cell_window = 10000
  batch_window = 100
  warm_up = 500000
  feature_loss = as.numeric(stri_extract(training_log, regex = "(?<=feature_loss: )[0-9.]+"))
  feature_loss = feature_loss[!is.na(feature_loss)]
  eval_interval = ceiling(250000 / cell_number)
  summary_table = read.table(paste0(stri_trim_right(log_path, pattern = "[/]"), "summary.tsv"), sep = "\t", header = T)
  feature_convergence_estimate = summary_table["feature_std"] * sqrt(nrow(summary_table) * cell_number / 3000000 + 1)
  estimation = feature_convergence_estimate
  optimal_point = round((which.min(feature_convergence_estimate[!is.na(feature_convergence_estimate)]) + sum(is.na(feature_convergence_estimate))) / eval_interval) * eval_interval
  if(optimal_point > (length(feature_loss) - 1)){
    optimal_point = (round(optimal_point / eval_interval) - 1) * eval_interval
  }
  return(paste0(stri_trim_right(log_path, pattern = "[/]"), "results/epoch_", optimal_point, "/imputation.hdf5"))
}

get_optimal_point33 = function(log_path){
  training_log = readLines(log_path)
  cell_number = as.numeric(stri_extract(training_log[5], regex = "[0-9]+"))
  cell_window = 10000
  warm_up = 500000
  feature_loss = as.numeric(stri_extract(training_log, regex = "(?<=feature_loss: )[0-9.]+"))
  feature_loss = feature_loss[!is.na(feature_loss)]
  eval_interval = ceiling(250000 / cell_number)
  summary_table = read.table(paste0(stri_trim_right(log_path, pattern = "[/]"), "summary.tsv"), sep = "\t", header = T)
  feature_convergence_estimate = summary_table["feature_std"] * sqrt(nrow(summary_table) * cell_number / 5000000 + 1)
  estimation = feature_convergence_estimate
  optimal_point = round((which.min(feature_convergence_estimate[!is.na(feature_convergence_estimate)]) + sum(is.na(feature_convergence_estimate))) / eval_interval) * eval_interval
  if(optimal_point > (length(feature_loss) - 1)){
    optimal_point = (round(optimal_point / eval_interval) - 1) * eval_interval
  }
  return(paste0(stri_trim_right(log_path, pattern = "[/]"), "results/epoch_", optimal_point, "/imputation.hdf5"))
}

get_optimal_point37 = function(log_path){
  training_log = readLines(log_path)
  cell_number = as.numeric(stri_extract(training_log[5], regex = "[0-9]+"))
  cell_window = 10000
  warm_up = 500000
  feature_loss = as.numeric(stri_extract(training_log, regex = "(?<=feature_loss: )[0-9.]+"))
  feature_loss = feature_loss[!is.na(feature_loss)]
  eval_interval = ceiling(250000 / cell_number)
  summary_table = read.table(paste0(stri_trim_right(log_path, pattern = "[/]"), "summary.tsv"), sep = "\t", header = T)
  feature_convergence_estimate = summary_table["feature_std"] * sqrt(nrow(summary_table) * cell_number / 5000000 + 2.5)
  estimation = feature_convergence_estimate
  optimal_point = round((which.min(feature_convergence_estimate[!is.na(feature_convergence_estimate)]) + sum(is.na(feature_convergence_estimate))) / eval_interval) * eval_interval
  if(optimal_point > (length(feature_loss) - 1)){
    optimal_point = (round(optimal_point / eval_interval) - 1) * eval_interval
  }
  return(paste0(stri_trim_right(log_path, pattern = "[/]"), "results/epoch_", optimal_point, "/imputation.hdf5"))
}


readh5_imputation = function(h5_path, use_genes=NULL, used_cells=NULL, with_outliers=F){
  print(h5_path)
  gene_name = h5read(h5_path, "gene_name")
  cell_id = h5read(h5_path, "cell_id")
  if(!is.null(used_cells)){
    if(is.character(used_cells)){
      cell_grasp_index = which(cell_id %in% used_cells)
    }else{
      cell_grasp_index = used_cells
    }
  }else{
    cell_grasp_index = c(1: length(cell_id))
  }
  if(!is.null(use_genes)){
    if(is.character(use_genes)){
      gene_grasp_index = which(gene_name %in% use_genes)
    }else{
      gene_grasp_index = use_genes
    }
  }else{
    gene_grasp_index = c(1: length(gene_name))
  }
  gene_bc_mat = h5read(h5_path, "imputation", index = list(gene_grasp_index, cell_grasp_index))
  if(with_outliers){
    gene_bc_mat = gene_bc_mat + h5read(h5_path, "outliers", index = list(gene_grasp_index, cell_grasp_index))
  }
  gene_bc_mat[gene_bc_mat == -1] = NA
  rownames(gene_bc_mat) = gene_name[gene_grasp_index]
  colnames(gene_bc_mat) = cell_id[cell_grasp_index]
  print(dim(gene_bc_mat))
  h5closeAll()
  return(gene_bc_mat)
}
readh5_loom = function(loom_path, use_genes=NULL, used_cells=NULL, is_feature=FALSE){
  print(loom_path)
  gene_name = h5read(loom_path, "row_attrs/Gene")
  cell_id = h5read(loom_path, "col_attrs/CellID")
  if(!is.null(used_cells)){
    if(is.character(used_cells)){
      cell_grasp_index = which(cell_id %in% used_cells)
    }else{
      cell_grasp_index = used_cells
    }
  }else{
    cell_grasp_index = c(1: length(cell_id))
  }
  if(!is.null(use_genes)){
    if(is.character(use_genes)){
      gene_grasp_index = which(gene_name %in% use_genes)
    }else{
      gene_grasp_index = use_genes
    }
  }else{
    gene_grasp_index = c(1: length(gene_name))
  }
  gene_bc_mat = t(h5read(loom_path, "matrix", index = list(cell_grasp_index, gene_grasp_index)))
  if(!is_feature){
    gene_bc_mat[gene_bc_mat == -1] = NA
  }
  row.names(gene_bc_mat) = gene_name[gene_grasp_index]
  colnames(gene_bc_mat) = cell_id[cell_grasp_index]
  print(dim(gene_bc_mat))
  return(gene_bc_mat)
}

get_gene_bc_mat = function(gene_bc_mat_path, use_genes=NULL, used_cells=NULL, with_outliers=F){
  file_format = get_last_element(unlist(strsplit(gene_bc_mat_path, ".", fix = T)))
  if(file_format == "loom"){
    gene_bc_mat = readh5_loom(gene_bc_mat_path, use_genes = use_genes, used_cells = used_cells)
  }else{
    gene_bc_mat = readh5_imputation(gene_bc_mat_path, use_genes = use_genes, used_cells = used_cells, with_outliers = with_outliers)
  }
  return(gene_bc_mat)
}

readh5_feature = function(h5_path){
  print(h5_path)
  cell_id = h5read(h5_path, "cell_id")
  gene_bc_mat = h5read(h5_path, "feature")
  colnames(gene_bc_mat) = cell_id
  rownames(gene_bc_mat) = paste0("feature_", seq(nrow(gene_bc_mat)))
  print(dim(gene_bc_mat))
  h5closeAll()
  return(gene_bc_mat)
}

get_loom_gene = function(loom_path){
  return(h5read(loom_path, "row_attrs/Gene"))
}

get_imputation_gene = function(running_info){
  return(h5read(running_info, "imputation_gene_name"))
}

get_imputation_gene1 = function(running_info){
  return(h5read(running_info, "target_gene"))
}

readh5_cell_type = function(h5_path){
  print(h5_path)
  cell_type = h5read(h5_path, "cell_type")
  names(cell_type) = h5read(h5_path, "cell_id")
  h5closeAll()
  return(cell_type)
}

readh5_sample_id = function(h5_path){
  print(h5_path)
  sample_id = h5read(h5_path, "col_attrs/SampleID")
  names(sample_id) = h5read(h5_path, "col_attrs/CellID")
  h5closeAll()
  return(sample_id)
}

save_h5 = function(output_path, bc_gene_mat){
  h5createFile(output_path)
  h5createDataset(file = output_path,
                  dataset = "matrix",
                  dims = dim(bc_gene_mat),
                  storage.mode = "double",
                  chunk=c(1, ncol(bc_gene_mat)))
  h5write(bc_gene_mat, file=output_path, name="matrix")
  h5createGroup(output_path, "col_attrs")
  h5write(rownames(bc_gene_mat), output_path,"col_attrs/CellID")
  h5createGroup(output_path, "row_attrs")
  h5write(colnames(bc_gene_mat), output_path,"row_attrs/Gene")
  h5createGroup(output_path, "layers")
  h5createGroup(output_path, "col_graphs")
  h5createGroup(output_path, "row_graphs")
}
# downsampling
log.pi.res <- function(z){
  exp(z)/(1+exp(z))
}

sample.for.dropout <- function(z){
  dropout.list <- NULL
  for(i in 1:length(z)){
    dropout.list <- c(dropout.list, sample(0:1, 1, prob=c(z[i], 1-z[i])))
  }
  return(dropout.list)
}

change.rate <- function(x1, x2, coef){
  log.x1 <- log2(x1 + 0.1)
  log.x2 <- log2(x2 + 0.1)
  log.x1.hat <- coef[1] + log.x1*coef[2]
  log.x2.hat <- coef[1] + log.x2*coef[2]
  rate.x1 <- log.pi.res(log.x1.hat)
  rate.x2 <- log.pi.res(log.x2.hat)
  dropout <- (rate.x2 - rate.x1)/(1 - rate.x1)
  dropout[dropout<0] <- 0.0001
  dropout.label <- sample.for.dropout(dropout)
  return(x2*dropout.label)
}

sample_reads = function(x, le){
  reads <- sum(x)
  new.reads <- round(reads*le)
  newsample <- sample(rep(names(x), x), size = new.reads, replace = F)
  newsample.tab <- table(newsample)
  newsample_vector = rep(0, length(x))
  names(newsample_vector) = names(x)
  newsample_vector[names(newsample.tab)] = newsample.tab
  return(newsample_vector)
}

downsampling_cell = function(le, gene_bc_mat){
  gene_bc_mat[Matrix::rowSums(gene_bc_mat) > 0, ] = apply(gene_bc_mat[Matrix::rowSums(gene_bc_mat) > 0, ], 2, function(x){
    return(sample_reads(x, le))
  })
  return(gene_bc_mat)
}

downsampling_pool = function(le, gene_bc_mat){
  dim_names = dimnames(gene_bc_mat)
  gene_bc_mat = as(gene_bc_mat, "dgTMatrix")
  transcript_i = c()
  transcript_j = c()
  transcript_table = table(gene_bc_mat@x)
  expression_values_c = names(transcript_table)
  if(sum(as.integer(expression_values_c) != as.numeric(expression_values_c)) != 0){
    stop("float exist!")
  }
  for(ii in expression_values_c){
    for(jj in seq(transcript_table[ii])){
      transcript_num = as.integer(ii)
      this_mask = gene_bc_mat@x == transcript_num
      for(kk in seq(transcript_num)){
        transcript_i = c(transcript_i, gene_bc_mat@i[this_mask])
        transcript_j = c(transcript_j, gene_bc_mat@j[this_mask])
      }
    }
  }
  reads <- length(transcript_i)
  new_reads <- round(reads * le)
  use_index = sample(seq(reads), size = new_reads, replace = F)
  transcript_i = transcript_i[use_index]
  transcript_j = transcript_j[use_index]
  ds_mat = Matrix::sparseMatrix(i = transcript_i + 1, j = transcript_j + 1, x = rep(1, length(transcript_i)), dimnames = dim_names)
  return(as.matrix(ds_mat))
}

downsampling_viper_dependent = function(le, gene_bc_mat){
  gene_mask = rowSums(gene_bc_mat) > 0
  logxx = log2(gene_bc_mat[gene_mask, ] + 0.1)
  zero.rate = round(rowSums(gene_bc_mat[gene_mask, ] == 0) / ncol(gene_bc_mat), 2)
  pi_ratio = log(zero.rate / (1 - zero.rate))
  pi_ratio[is.infinite(pi_ratio)] = NA
  coeff = summary(lm(pi_ratio ~ rowMeans(logxx)))$coefficients[, 1]
  gene_bc_mat[gene_mask, ] = apply(gene_bc_mat[gene_mask, ], 2, function(x){
    ds_cell = sample_reads(x, le)
    na_mask = !is.na(ds_cell)
    ds_cell[na_mask] = change.rate(x[na_mask], ds_cell[na_mask], coeff)
    return(ds_cell)
  })
  return(gene_bc_mat)
}

downsampling_viper_independent = function(le, gene_bc_mat){
  gene_mask = rowSums(gene_bc_mat) > 0
  #keep_rate = 0.95 ^ ((1 - le) / 0.05) #  raw = 0.98
  keep_rate = 0.5
  gene_bc_mat[gene_mask, ] = apply(gene_bc_mat[gene_mask, ], 2, function(x){
    ds_cell = sample_reads(x, le)
    na_mask = !is.na(ds_cell)
    ds_cell[na_mask] = ds_cell[na_mask] * sample(c(0, 1), sum(na_mask), replace = TRUE, prob = c(1 - keep_rate, keep_rate))
    return(ds_cell)
  })
  return(gene_bc_mat)
}

downsampling_for_splatter = function(gene_bc_mat, d_min=0.01, d_max=0.05){
  gene_bc_mat[Matrix::rowSums(gene_bc_mat) > 0, ] = t(apply(gene_bc_mat[Matrix::rowSums(gene_bc_mat) > 0, ], 1, function(x){
    return(sample_reads(x, runif(1, d_min, d_max)))
  }))
  return(gene_bc_mat)
}

#  strings
delete_last_element <- function(x){
  return(x[1: (length(x) - 1)])
}
get_last_element <- function(x){
  return(x[length(x)])
}
get_first_element <- function(x){
  return(x[1])
}
firstup <- function(x) {
  x = tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

delete_lt0.5 = function(gene_bc_mat){
  gene_bc_mat[gene_bc_mat < 0.5] = 0
  return(gene_bc_mat)
}

rmse = function(x, y){
  return(sqrt(mean((x - y)^2)))
}

rsq = function (x, y){
  return(cor(x, y) ^ 2)
}

jaccard <- function(x, y){
  length(intersect(x, y))/length(union(x, y))
}

quadct = function(x, y, xx, yy){
  n = length(xx)
  ix1 = xx <= x
  ix2 = yy <= y
  a = sum(ix1 & ix2) / n
  b = sum(ix1 & !ix2) / n
  c = sum(!ix1 & ix2) / n
  d = 1 - a - b - c
  return(list(a = a, b = b, c = c, d = d))
}

maxdist = function(x1, y1, x2, y2){
  n1 = length(x1)
  D1 = t(sapply(seq(n1), function(x){
    list1 = quadct(x1[x], y1[x], x1, y1)
    list2 = quadct(x1[x], y1[x], x2, y2)
    return(c(list1$a - list2$a, list1$b - list2$b, list1$c - list2$c, list1$d - list2$d))
  }))
  # re-assign the point to maximize difference,
  # the discrepancy is significant for N < ~50
  D1[, 1] = D1[, 1] - (1 / n1)
  dmin = -min(D1)
  dmax = max(D1) + 1 / n1
  return(max(c(dmin, dmax)))
}

avgmaxdist = function(x1, y1, x2, y2){
  D1 = maxdist(x1, y1, x2, y2)
  D2 = maxdist(x2, y2, x1, y1)
  return((D1 + D2) / 2)
}
  
ks2d2s = function(x1, y1, x2, y2){
  #Two-dimensional Kolmogorov-Smirnov test on two samples.
  #Parameters
  #----------
  #x1, y1 : ndarray, shape (n1, )
  #Data of sample 1.
  #x2, y2 : ndarray, shape (n2, )
  #Data of sample 2. Size of two samples can be different.
  #Returns
  #-------
  #D : float
  #KS statistic.
  #References
  #----------
  #Peacock, J.A. 1983, Two-Dimensional Goodness-of-Fit Testing in Astronomy, Monthly Notices of the Royal Astronomical Society, vol. 202, pp. 615-627
  #Fasano, G. and Franceschini, A. 1987, A Multidimensional Version of the Kolmogorov-Smirnov Test, Monthly Notices of the Royal Astronomical Society, vol. 225, pp. 155-170
  #Press, W.H. et al. 2007, Numerical Recipes, section 14.8
  if(!(length(x1) == length(y1) & length(x2) == length(y2))){
    stop()
  }
  D = avgmaxdist(x1, y1, x2, y2)
  return(D)
}
  

cal_2D_KL_divergence = function(p_raw, q_raw){
  p1 = as.vector(p_raw)[as.vector(p_raw) > 0]
  q1 = as.vector(q_raw)[as.vector(p_raw) > 0]
  q1[q1 == 0] = 0.5
  KL1 = scipy$stats$entropy(np_array(p1), np_array(q1))
  p2 = as.vector(p_raw)[as.vector(q_raw) > 0]
  q2 = as.vector(q_raw)[as.vector(q_raw) > 0]
  p2[p2 == 0] = 0.5
  KL2 = scipy$stats$entropy(np_array(q2), np_array(p2))
  return(KL1 + KL2)
}

cal_KL_divergence = function(p, q){
  p = p[!is.na(p)]
  q = q[!is.na(q)]
  p_mat = as.matrix(aggregate(rep(1 / length(p), length(p)), by = list(round(p)), FUN = "sum"))
  q_mat = as.matrix(aggregate(rep(1 / length(q), length(q)), by = list(round(q)), FUN = "sum"))
  max_value = max(c(p_mat[, 1], q_mat[, 1]))
  p_vec = rep(0, max_value + 1)
  names(p_vec) = c(0:max_value)
  q_vec = rep(0, max_value + 1)
  names(q_vec) = c(0:max_value)
  p_vec[as.character(p_mat[, 1])] = p_mat[, 2]
  q_vec[as.character(q_mat[, 1])] = q_mat[, 2]
  return(cal_2D_KL_divergence(p_vec, q_vec))
}

gamma_result = function(saver_result, num_of_obs=5){
  saver_est = saver_result$estimate
  saver_se = saver_result$se
  saver_sf = saver_result$info$size.factor
  result_mat = t(sapply(1:nrow(saver_est), FUN = function(ii){
    rgamma(n = num_of_obs * ncol(saver_est),
           shape = saver_est[ii, ]^2/saver_se[ii, ]^2,
           rate = saver_est[ii, ]/saver_se[ii, ]^2/saver_sf)
  }))
  rownames(result_mat) = rownames(saver_est)
  colnames(result_mat) = colnames(saver_est)
  return(result_mat)
}

saver_filter_fun = function(gene_bc_mat){
  gene_bc_mat[, gene_bc_mat["GAPDH", ] < quantile(gene_bc_mat["GAPDH", ], 0.9, na.rm = TRUE) &
                gene_bc_mat["GAPDH", ] > quantile(gene_bc_mat["GAPDH", ], 0.1, na.rm = TRUE)]
}

saver_norm_fun = function(gene_bc_mat){
  sweep(gene_bc_mat, 2, gene_bc_mat["GAPDH", ] / mean(gene_bc_mat["GAPDH", ]), "/")
}

binom_result = function(saver_result, fish_gene_bc_mat){
  saver_est = saver_result$estimate
  num_cells = ncol(saver_est)
  intersect_genes = intersect(row.names(saver_est), row.names(fish_gene_bc_mat))
  fish_filt_as_saver <- saver_filter_fun(fish_gene_bc_mat)
  filt_as_saver <- saver_filter_fun(saver_est[intersect_genes, ])
  mult_factor <- sapply(matrix(intersect_genes), function(i)
    mean(fish_filt_as_saver[i, ], na.rm = TRUE)/mean(filt_as_saver[i, ]))
  t(sapply(matrix(intersect_genes), function(i)
    rnbinom(num_cells,
            mu = saver_est[i, ] * mult_factor[i],
            size = saver_est[i, ]^2/saver_result$se[i, ]^2)))
}

binom_result_mean = function(saver_result, fish_gene_bc_mat){
  saver_est = saver_result$estimate
  num_cells = ncol(saver_est)
  intersect_genes = intersect(row.names(saver_est), row.names(fish_gene_bc_mat))
  mult_factor <- sapply(matrix(intersect_genes), function(i)
    mean(fish_gene_bc_mat[i, ], na.rm = TRUE)/mean(saver_est[i, ]))
  t(sapply(matrix(intersect_genes), function(i)
    rnbinom(num_cells,
            mu = saver_est[i, ] * mult_factor[i],
            size = saver_est[i, ]^2/saver_result$se[i, ]^2)))
}

saver_density_norm_method = function(impute_result, fish_gene_bc_mat){
  fish_filt_as_saver <- saver_filter_fun(fish_gene_bc_mat)
  intersect_genes = intersect(rownames(fish_gene_bc_mat), rownames(impute_result))
  filt_as_saver <- saver_filter_fun(impute_result[intersect_genes, ])
  norm_as_saver <- saver_norm_fun(filt_as_saver)
  mult_factor <- sapply(matrix(intersect_genes), function(x) mean(fish_filt_as_saver[x, ], na.rm = TRUE)/mean(filt_as_saver[x, ]))
  return(sweep(norm_as_saver, 1, mult_factor, "*"))
}

mean_norm_fun = function(gene_bc_mat, fish_gene_bc_mat, return_factor = F){
  intersect_genes = intersect(rownames(fish_gene_bc_mat), rownames(gene_bc_mat))
  if(return_factor){
    return(sapply(intersect_genes, function(x){
      mean(fish_gene_bc_mat[x, ], na.rm = TRUE) / mean(gene_bc_mat[x, ])
    }))
  }else{
    return(t(sapply(intersect_genes, function(x){
      gene_bc_mat[x, ] * mean(fish_gene_bc_mat[x, ], na.rm = TRUE) / mean(gene_bc_mat[x, ])
    })))
  }
}

grid_switch = function(target_mat, input_mat, output_all = F){
  x_max = round(max(c(target_mat[, 1], input_mat[, 1]))) + 1
  y_max = round(max(c(target_mat[, 2], input_mat[, 2]))) + 1
  grid_output = list()
  if(output_all){
    target_mat_grid = Matrix::spMatrix(nrow = x_max, ncol = y_max)
    target_tmp1 = cbind(apply(target_mat, 1, paste, collapse = "_"), 1 / nrow(target_mat))
    target_tmp2 = aggregate(as.numeric(target_tmp1[, 2]), by = list(target_tmp1[, 1]), FUN = "sum")
    target_tmp3 = matrix(unlist(strsplit(target_tmp2[, 1], split = "_")), nrow = 2)
    target_mat_grid@i = as.integer(target_tmp3[1, ])
    target_mat_grid@j = as.integer(target_tmp3[2, ])
    target_mat_grid@x = target_tmp2[, 2]
    grid_output$target_grid = target_mat_grid
  }
  ### input ###
  input_mat_grid = Matrix::spMatrix(nrow = x_max, ncol = y_max)
  input_tmp1 = cbind(apply(round(input_mat), 1, paste, collapse = "_"), 1 / nrow(input_mat))
  input_tmp2 = aggregate(as.numeric(input_tmp1[, 2]), by = list(input_tmp1[, 1]), FUN = "sum")
  input_tmp3 = matrix(unlist(strsplit(input_tmp2[, 1], split = "_")), nrow = 2)
  input_mat_grid@i = as.integer(input_tmp3[1, ])
  input_mat_grid@j = as.integer(input_tmp3[2, ])
  input_mat_grid@x = input_tmp2[, 2]
  ### output ###
  grid_output$input_grid = input_mat_grid
  return(grid_output)
}



grid_compare = function(target_mat, input_mat){
  x_max = round(max(c(target_mat[, 1], input_mat[, 1]))) + 1
  y_max = round(max(c(target_mat[, 2], input_mat[, 2]))) + 1
  target_mat_grid = Matrix::spMatrix(nrow = x_max, ncol = y_max)
  input_mat_grid = Matrix::spMatrix(nrow = x_max, ncol = y_max)
  target_tmp1 = cbind(apply(target_mat, 1, paste, collapse = "_"), 1 / nrow(target_mat))
  target_tmp2 = aggregate(as.numeric(target_tmp1[, 2]), by = list(target_tmp1[, 1]), FUN = "sum")
  target_tmp3 = matrix(unlist(strsplit(target_tmp2[, 1], split = "_")), nrow = 2)
  target_mat_grid@i = as.integer(target_tmp3[1, ])
  target_mat_grid@j = as.integer(target_tmp3[2, ])
  target_mat_grid@x = target_tmp2[, 2]
  ### input ###
  input_tmp1 = cbind(apply(round(input_mat), 1, paste, collapse = "_"), 1 / nrow(input_mat))
  input_tmp2 = aggregate(as.numeric(input_tmp1[, 2]), by = list(input_tmp1[, 1]), FUN = "sum")
  input_tmp3 = matrix(unlist(strsplit(input_tmp2[, 1], split = "_")), nrow = 2)
  input_mat_grid@i = as.integer(input_tmp3[1, ])
  input_mat_grid@j = as.integer(input_tmp3[2, ])
  input_mat_grid@x = input_tmp2[, 2]
  ### output ###
  grid_compare = list()
  grid_compare$result = sum(abs(target_mat_grid - input_mat_grid)) / 2 ### 0 - 1: Not overlap degree
  grid_compare$target_grid = target_mat_grid
  grid_compare$input_grid = input_mat_grid
  return(grid_compare)
}

compute_use_color = function(grid_mat, max_points){
  grid_mat = as.matrix(grid_mat)
  grid_hue = sort(unique(c(-log(grid_mat))), decreasing = T) / log(max_points) * 0.72 - 0.1
  grid_hue = c(0.72, grid_hue[is.finite(grid_hue)])
  grid_hue[grid_hue < 0] = grid_hue[grid_hue < 0] + 1
  return(hsv(grid_hue, 0.85, 1))
}

heatmap_plot = function(grid_mat, use_colors, xlab = "", ylab = "", main = "", col.main = "black"){
  grid_mat = as.matrix(grid_mat)
  plot(0, 0, type = "n",
       xlim = c(0, nrow(grid_mat)),
       ylim = c(0, ncol(grid_mat)),
       axes = FALSE, xlab = xlab, ylab = ylab, main = main, col.main=col.main, cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
  xleft = rep(c(0:(nrow(grid_mat) - 1)), ncol(grid_mat))
  xright = rep(c(1:nrow(grid_mat)), ncol(grid_mat))
  ybottom = rep(c(0:(ncol(grid_mat) - 1)), each = nrow(grid_mat))
  ytop = rep(c(1:ncol(grid_mat)), each = nrow(grid_mat))
  grid_factor = as.numeric(as.factor(grid_mat))
  colors = use_colors[grid_factor]
  rect(xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, col=colors, border = NA)
}



barplot_usage = function(data_vector, main, bar_color, text_color=NULL, use_data_order=F, decreasing=F, standard_error=NULL, show_number=F, ylim=NULL, cex.main=2, use_log1p=F, use_border=T, ...){
  data_order = c(1, order(data_vector[-1], decreasing = decreasing) + 1)
  if(use_data_order){
    data_vector = data_vector[data_order]
    bar_color = bar_color[data_order]
    if(!is.null(text_color)){
      text_color = text_color[data_order]
    }
    if(!is.null(standard_error)){
      standard_error = standard_error[data_order]
    }
  }
  if(use_log1p){
    plot_data = log1p(data_vector)
    ylim = log1p(ylim)
    if(!is.null(standard_error)){
      plot_data_up = log1p(data_vector + standard_error)
      plot_data_down = log1p(data_vector - standard_error)
    }
  }else{
    plot_data = data_vector
    if(!is.null(standard_error)){
      plot_data_up = data_vector + standard_error
      plot_data_down = data_vector - standard_error
    }
  }
  if(is.null(ylim)){
    if(is.null(standard_error)){
      ylim = c(NA, max(plot_data) * 1.01)
    }else{
      ylim = c(NA, max(plot_data + plot_data_up) * 1.01)
    }
    if(!is.null(text_color)){
      ylim[1] = -0.1 * ylim[2]
    }else{
      ylim[1] = 0
    }
  }
  if(use_border){
    border = par("fg")
  }else{
    border = NA
  }
  bp = barplot(plot_data, main = main, las=1, names.arg="", col = bar_color, ylim = ylim, cex.axis = 1.2, border = border, cex.main = cex.main, ...)
  if(!is.null(standard_error)){
    arrows(bp, plot_data_down, bp, plot_data_up, length=0.05, angle=90, code=3)
  }
  if(!is.null(text_color)){
    text_position = -0.015 * max(data_vector)
    text(bp, sapply(data_vector,
                    function(x){
                      if(x >= 0){
                        return(text_position)
                      }
                      else{
                        return(x - text_position)
                      }
                    }),
         srt = 90, adj= 1, xpd = TRUE, labels = names(data_vector), cex=1.5, col = text_color)
    if(show_number){
      text(bp, sapply(data_vector,
                      function(x){
                        if(x >= 0){
                          return(x)
                        }
                        else{
                          return(0)
                        }}),
           srt = 90, adj= 0, xpd = TRUE, labels = round(data_vector, digits = 3) , cex=1.5, col = text_color)
    }
  }
}


barplot_usage_new = function(data_vector, main, method_color, use_log1p=F, use_data_order=F, decreasing=F, standard_error=NULL, cex.main = 2, ylim = NULL, ...){
  data_order = c(1, order(data_vector[-1], decreasing = decreasing) + 1)
  if(use_data_order){
    data_vector = data_vector[data_order]
    method_color = method_color[data_order]
    if(!is.null(standard_error)){
      standard_error = standard_error[data_order]
    }
  }
  if(is.null(ylim)){
    if(is.null(standard_error)){
      ylim = c(0, max(data_vector) * 1.01)
    }else{
      ylim = c(0, max(data_vector + standard_error) * 1.01)
    }
  }
  if(use_log1p){
    plot_data = log1p(data_vector)
    ylim = log1p(ylim)
    if(!is.null(standard_error)){
      plot_data_up = log1p(data_vector + standard_error)
      plot_data_down = log1p(data_vector - standard_error)
    }
  }else{
    plot_data = data_vector
    if(!is.null(standard_error)){
      plot_data_up = data_vector + standard_error
      plot_data_down = data_vector - standard_error
    }
  }
  bp = barplot(plot_data, main = main, las=1, names.arg="", col = method_color, cex.axis = 1.2, cex.main = cex.main, border = NA, ylim = ylim, ...)
  if(!is.null(standard_error)){
    arrows(bp, plot_data_down, bp, plot_data_up, length=0.05, angle=90, code=3)
  }
}

calc_cmd <- function(R1, R2) {
  na_mask = is.na(R1) | is.na(R2)
  R1[na_mask] = 0
  R2[na_mask] = 0
  traceR1R2 <- sum(R1 * t(R2))
  R1.norm <- norm(R1, type = "F")
  R2.norm <- norm(R2, type = "F")
  return(1-traceR1R2/(R1.norm*R2.norm))
}

correlogram_plot = function(correlation_matrix, xlab = "", ylab = "", main = "", col.main = "black", plot.cex = 0.8, cex.axis = 0.9, diag = T, use_label = F, plot_color_legend = F, plot_NA_legend = F){
  correlation_matrix = as.matrix(correlation_matrix)
  if(!use_label){
    par(mar = c(0, 0, 2.5, 0))
  }
  plot(0, 0, type = "n",
       xlim = c(0, nrow(correlation_matrix)),
       ylim = rev(c(0, ncol(correlation_matrix))),
       axes = FALSE, xlab = xlab, ylab = ylab, main = main, col.main=col.main, cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
  xleft = rep(c(0:(nrow(correlation_matrix) - 1)), ncol(correlation_matrix))
  xright = rep(c(1:nrow(correlation_matrix)), ncol(correlation_matrix))
  ybottom = rep(c(0:(ncol(correlation_matrix) - 1)), each = nrow(correlation_matrix))
  ytop = rep(c(1:ncol(correlation_matrix)), each = nrow(correlation_matrix))
  use_colors = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061"))(200)
  plot_color = c(use_colors[1:100], "#FFFFFF", use_colors[101:200])
  display_mask = matrix(T, nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
  if(diag){
    for(ii in 1:(nrow(correlation_matrix)-1)){
      display_mask[(ii+1):nrow(correlation_matrix), ii] = F
    }
  }
  display_vector = as.vector(display_mask)
  na_mask = is.na(correlation_matrix)
  na_vector = as.vector(na_mask)
  correlation_matrix[na_mask] = 0
  correlation_matrix = round((correlation_matrix + 1) * 100) + 1
  colors = plot_color[as.vector(correlation_matrix)]
  colors[na_vector] = "#FFFFFF"
  #rect(xleft = xleft[display_vector], ybottom = ybottom[display_vector], xright = xright[display_vector], ytop = ytop[display_vector], col=colors[display_vector], border = NA)
  #if(sum(na_vector & display_vector) > 0){
  #  text(xleft[na_vector & display_vector] + 0.5, ybottom[na_vector & display_vector] + 0.5, "NA", cex = plot.cex)
  #}
  rect(xleft = xleft[na_vector & display_vector], ybottom = ybottom[na_vector & display_vector], xright = xright[na_vector & display_vector], ytop = ytop[na_vector & display_vector], col="white", border = NA)
  if(sum(na_vector & display_vector) > 0){
    for(ii in which(na_vector & display_vector)){
      lines(x = c(xleft[ii], xright[ii]), y = c(ytop[ii], ybottom[ii]), col="gray80", lwd = 2)
      lines(x = c(xleft[ii], xright[ii]), y = c(ybottom[ii], ytop[ii]), col="gray80", lwd = 2)
    }
    has_NA = T
  }else{
    has_NA = F
  }
  rect(xleft = xleft[!na_vector & display_vector], ybottom = ybottom[!na_vector & display_vector], xright = xright[!na_vector & display_vector], ytop = ytop[!na_vector & display_vector], col=colors[!na_vector & display_vector], border = NA)
  if(use_label){
    axis(1, at = seq(nrow(correlation_matrix)) - 0.5, labels = rownames(correlation_matrix), cex.axis = cex.axis, tick = F, las = 2, hadj = 1)
    axis(2, at = seq(ncol(correlation_matrix)) - 0.5, labels = colnames(correlation_matrix), cex.axis = cex.axis, tick = F, las = 1, hadj = 1)
  }
  if(plot_color_legend){
    par(mar = c(3, 0, 3, 0))
    xright = 30
    plot(c(0,length(plot_color) - 1), c(-1, 1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
    ticks.at = seq(from=-1, to=1, length=3)
    axis(side = 4, pos=xright, at=ticks.at, labels=ticks.at, las=1, cex.axis=1.5)
    if(plot_NA_legend){
      mtext(side = 1, "X", line=1.5,at=par("usr")[1]+0.1*diff(par("usr")[1:2]),cex=1.2, col = "gray80")
      mtext(side = 1, ": NA", line=1.5, at=par("usr")[1]+0.5*diff(par("usr")[1:2]), cex=1.2, col = "black")
    }
    this_scale = length(plot_color) / 2
    for (ii in 1:length(plot_color)) {
      y = (ii - 1) / this_scale - 1
      rect(xleft = 0, ybottom = y, xright = xright, ytop = y + 1 / this_scale, col=plot_color[ii], border=NA)
    }
  }
  return(has_NA)
}

layout_correlogram_plot = function(cor_all, lab_cex = 3, this_xlab = "Gene2", this_ylab = "Gene1", plot_width = 3.5, plot_height = 4, plot_row = 2, use_order=NULL, ...){
  plot_region = t(matrix(seq(length(cor_all)), ncol = plot_row))
  this_height = rep(plot_height, nrow(plot_region))
  this_width = rep(plot_width, ncol(plot_region))
  this_index = max(plot_region)
  layout_mat = plot_region
  plot_col = rep("p", ncol(plot_region))
  use_legend = T
  if(use_legend){
    this_index = this_index + 1
    this_width = c(this_width, 1)
    legend_region = matrix(c(rep(0, plot_row - 1), this_index), nrow = plot_row)
    layout_mat = cbind(layout_mat, legend_region)
    plot_col = c(plot_col, 0)
  }
  if(!is.na(this_ylab)){
    this_index = this_index + 1
    this_width = c(1, this_width)
    ylab_region = matrix(rep(this_index, plot_row), nrow = plot_row)
    layout_mat = cbind(ylab_region, layout_mat)
    plot_col = c(0, plot_col)
  }
  if(!is.na(this_xlab)){
    this_index = this_index + 1
    this_height = c(this_height, 1)
    xlab_region = matrix(plot_col, nrow = 1)
    xlab_region[xlab_region == "p"] = this_index
    layout_mat = rbind(layout_mat, xlab_region)
  }
  layout(mat = layout_mat, heights = this_height, widths = this_width)
  cmd_vector = c()
  plot_NA_legend = F
  for(ii in seq(length(cor_all))){
    if(ii == length(cor_all)){
      plot_color_legend = T
    }else{
      plot_color_legend = F
    }
    this_title = names(cor_all)[ii]
    if(!is.null(use_order)){
      plot_cor_mat = cor_all[[this_title]][use_order, use_order]
    }else{
      plot_cor_mat = cor_all[[this_title]]
    }
    if(this_title == "DeSCI"){
      plot_NA_legend = plot_NA_legend + correlogram_plot(plot_cor_mat, main = this_title, col.main = "red", cex.axis = 0.7, plot_color_legend = plot_color_legend, plot_NA_legend = plot_NA_legend, ...)
    }else{
      plot_NA_legend = plot_NA_legend + correlogram_plot(plot_cor_mat, main = this_title, col.main = "black", cex.axis = 0.7, plot_color_legend = plot_color_legend, plot_NA_legend = plot_NA_legend, ...)
    }
    cmd_vector = c(cmd_vector, calc_cmd(cor_all[["FISH"]], cor_all[[this_title]]))
  }
  if(!is.na(this_ylab)){
    par(mar = rep(0, 4))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, this_ylab, srt = 90, cex = lab_cex)
  }
  if(!is.na(this_xlab)){
    par(mar = rep(0, 4))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, this_xlab, cex = lab_cex)
  }
  return(cmd_vector)
}

layout_scatter = function(result_list, method_names, mask_or_index=NULL, color_point=NULL, lab_cex = 3, plot_width = 3.5, plot_height = 4, plot_row = 2, this_xlab=NULL, this_ylab=NULL, ...){
  plot_number = ceiling(length(method_names) / plot_row) * plot_row
  plot_region = t(matrix(seq(plot_number), ncol = plot_row))
  plot_region[plot_region > length(method_names)] = 0
  this_height = rep(plot_height, nrow(plot_region))
  this_width = rep(plot_width, ncol(plot_region))
  this_index = max(plot_region)
  layout_mat = plot_region
  plot_col = rep("p", ncol(plot_region))
  if(!is.null(this_ylab)){
    this_index = this_index + 1
    this_width = c(1, this_width)
    ylab_region = matrix(rep(this_index, plot_row), nrow = plot_row)
    layout_mat = cbind(ylab_region, layout_mat)
    plot_col = c(0, plot_col)
  }
  if(!is.null(this_xlab)){
    this_index = this_index + 1
    this_height = c(this_height, 1)
    xlab_region = matrix(plot_col, nrow = 1)
    xlab_region[xlab_region == "p"] = this_index
    layout_mat = rbind(layout_mat, xlab_region)
  }
  layout(mat = layout_mat, heights = this_height, widths = this_width)
  fish_names = names(result_list[["FISH"]])
  fish_values = as.numeric(result_list[["FISH"]])
  names(fish_values) = fish_names
  if(!is.null(mask_or_index)){
    fish_values = fish_values[mask_or_index]
  }
  par(mar = c(2, 2, 3, 1))
  all_rmse = c()
  cex = 4.25
  cex.text = 1.5
  pch = 21
  for(ii in method_names){
    if(ii == "DeSCI"){
      col.main = "red"
    }else{
      col.main = "black"
    }
    this_names = names(result_list[[ii]])
    this_values = as.numeric(result_list[[ii]])
    names(this_values) = this_names
    if(!is.null(mask_or_index)){
      this_values = this_values[mask_or_index]
    }
    this_rmse = rmse(fish_values, this_values)
    all_rmse = c(all_rmse, this_rmse)
    if(is.null(color_point)){
      normal_genes_mask = rep(TRUE, length(mask_or_index))
    }else{
      special_genes_mask = mask_or_index %in% names(color_point)
      normal_genes_mask = !special_genes_mask
    }
    plot(fish_values[normal_genes_mask], this_values[normal_genes_mask], pch = pch, bg = "#bdd7e7", cex.main = cex.text, cex.axis = cex.text, cex = cex, ylab = "", main = ii, col.main = col.main, font.main = 1, bty="n", ...)
    if(!is.null(color_point)){
      special_index = which(special_genes_mask)
      for(jj in seq(sum(special_genes_mask))){
        points(fish_values[special_index[jj]], this_values[special_index[jj]], pch = pch, bg = color_point[jj], cex = cex)
      }
    }
    legend("bottomright", paste0("RMSE = ", round(rmse(fish_values, this_values), 4)), bty = "n", cex = cex.text)
    abline(0, 1, col = "gray", lty = 2)
  }
  names(all_rmse) = method_names
  if(!is.null(this_ylab)){
    par(mar = rep(0, 4))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, this_ylab, srt = 90, cex = lab_cex)
  }
  if(!is.null(this_xlab)){
    par(mar = rep(0, 4))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, this_xlab, cex = lab_cex)
  }
  return(all_rmse)
}

adjustedRandIndex = function (x, y){
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1)))
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

FindVariableFeatures_vst_by_genes = function(gene_bc_mat){
  clip.max <- sqrt(ncol(gene_bc_mat))
  loess.span = 0.3
  hvf.info <- data.frame(mean = rowMeans(gene_bc_mat))
  hvf.info$variance <- apply(gene_bc_mat, 1, var)
  hvf.info$variance.expected <- 0
  hvf.info$variance.standardized <- 0
  not.const <- hvf.info$variance > 0
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, ],
    span = loess.span
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  # use c function to get variance after feature standardization
  hvf.info$variance.standardized = sapply(seq(nrow(gene_bc_mat)), function(x){
    standardized_value = (gene_bc_mat[x, ] - hvf.info$mean[x]) / sqrt(hvf.info$variance.expected[x])
    standardized_value[standardized_value > clip.max] = clip.max
    return(var(standardized_value))
  })
  return(hvf.info)
}

smoothScatter1 = function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white",
                                                                                            blues9)), nrpoints = 100, ret.selection = FALSE, pch = ".",
                           cex = 1, col = "black", transformation = function(x) x^0.25,
                           postPlotHook = box, xlab = NULL, ylab = NULL, xlim, ylim,
                           xaxs = par("xaxs"), yaxs = par("yaxs"), ...)
{
  if (!is.numeric(nrpoints) || nrpoints < 0 || length(nrpoints) !=
      1)
    stop("'nrpoints' should be numeric scalar with value >= 0.")
  nrpoints <- round(nrpoints)
  ret.selection <- ret.selection && nrpoints > 0
  xlabel <- if (!missing(x))
    deparse(substitute(x))
  ylabel <- if (!missing(y))
    deparse(substitute(y))
  xy <- xy.coords(x, y, xlabel, ylabel)
  xlab <- if (is.null(xlab))
    xy$xlab
  else xlab
  ylab <- if (is.null(ylab))
    xy$ylab
  else ylab
  x <- cbind(xy$x, xy$y)[I <- is.finite(xy$x) & is.finite(xy$y),
                         , drop = FALSE]
  if (ret.selection)
    iS <- which(I)
  if (!missing(xlim)) {
    stopifnot(is.numeric(xlim), length(xlim) == 2, is.finite(xlim))
    x <- x[I <- min(xlim) <= x[, 1] & x[, 1] <= max(xlim),
           , drop = FALSE]
    if (ret.selection)
      iS <- iS[I]
  }
  else {
    xlim <- range(x[, 1])
  }
  if (!missing(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim) == 2, is.finite(ylim))
    x <- x[I <- min(ylim) <= x[, 2] & x[, 2] <= max(ylim),
           , drop = FALSE]
    if (ret.selection)
      iS <- iS[I]
  }
  else {
    ylim <- range(x[, 2])
  }
  map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
  xm <- map$x1
  ym <- map$x2
  dens <- map$fhat
  dens[] <- transformation(dens)
  image(xm, ym, z = dens, col = colramp(16), xlab = xlab,
        ylab = ylab, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs,
        ...)
  if (!is.null(postPlotHook))
    postPlotHook()
  if (nrpoints > 0) {
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    stopifnot((nx <- length(xm)) == nrow(dens), (ny <- length(ym)) ==
                ncol(dens))
    ixm <- 1L + as.integer((nx - 1) * (x[, 1] - xm[1])/(xm[nx] -
                                                          xm[1]))
    iym <- 1L + as.integer((ny - 1) * (x[, 2] - ym[1])/(ym[ny] -
                                                          ym[1]))
    sel <- order(dens[cbind(ixm, iym)])[seq_len(nrpoints)]
    x <- x[sel, , drop = FALSE]
    points(x, pch = pch, cex = cex, col = col)
    if (ret.selection)
      iS[sel]
  }
}

dropout_rate_plot = function(gene_bc_mat, main){
  dropout_rate = rowSums(gene_bc_mat == 0) / ncol(gene_bc_mat)
  mean_gene_expression = log1p(rowSums(gene_bc_mat) / rowSums(gene_bc_mat > 0))
  plot(mean_gene_expression, dropout_rate,
       main = main, xlim = c(0, max(mean_gene_expression, na.rm=T)), ylim = c(0, 1),
       pch = 16, cex = 1, col = alpha("#4752cc", 0.3),
       xlab = "log(mean expression + 1)", ylab = "Dropout rate for each gene")
}

###cell_type_mapping###
cluster_evaluation_pbmc = function(this_markers, this_metadata, prior_cell_type=NULL){
  if(!is.null(prior_cell_type)){
    use_cell <-intersect(rownames(this_metadata), names(prior_cell_type))
    this_metadata = this_metadata[use_cell, ]
    this_metadata$cell_type = prior_cell_type[use_cell]
    cell_number = length(use_cell)
  }else{
    this_metadata$cell_type
    cell_number = nrow(this_metadata)
  }
  this_markers$cluster = as.numeric(as.character(this_markers$cluster))
  ind1 <-which(this_markers[,7]=="IL7R")
  a <- this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="CD14")
  b <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="LYZ")
  c <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="MS4A1")
  d <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="CD8A")
  e <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="FCGR3A")
  f <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="MS4A7")
  g <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="GNLY")
  h <- this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="NKG7")
  i <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="FCER1A")
  j <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="CST3")
  k <-this_markers[ind1,][,6]
  ind1 <-which(this_markers[,7]=="PPBP")
  l <-this_markers[ind1,][,6]

  CD4 <-setdiff(a,e)
  CD14 <-intersect(b,c)
  B <-d
  CD8 <-e
  Mono <-intersect(f,g)
  NK <-setdiff(intersect(h,i),e)
  DC <-intersect(j,k)
  Platelet <-l

  union <-c(CD4,CD14,B,CD8,Mono,NK,DC,Platelet)
  dup <- unique(union[duplicated(union)])
  CD41 <-setdiff(CD4,dup)
  CD141 <-setdiff(CD14,dup)
  B1 <-setdiff(B,dup)
  CD81 <-setdiff(CD8,dup)
  Mono1 <-setdiff(Mono,dup)
  NK1 <-setdiff(NK,dup)
  DC1 <-setdiff(DC,dup)
  Platelet1 <-setdiff(Platelet,dup)

  assigned <-unique(c(CD41,CD141,B1,CD81,Mono1,NK1,DC1,Platelet1))
  A <-unique(this_markers[,6])
  unassignment = setdiff(A,assigned)

  Others = unique(c(unassignment, dup))

  return_list = list()
  return_list[["CD4 T"]] = CD41
  return_list[["CD14+ Mono"]] = CD141
  return_list[["B"]] = B1
  return_list[["CD8 T"]] = CD81
  return_list[["FCGR3A+ Mono"]] = Mono1
  return_list[["NK"]] = NK1
  return_list[["DC"]] = DC1
  return_list[["Platelet"]] = Platelet1
  return_list[["Unresolved"]] = Others
  return_list[["Unassignment"]] = unassignment
  return_list[["Duplication"]] = dup
  this_metadata$assignment = "cluster_outliers"
  for(ii in names(return_list)){
    return_list[[ii]] = sort(return_list[[ii]])
    this_metadata$assignment[this_metadata$seurat_clusters %in% return_list[[ii]]] = ii
  }
  return_list[["Unassignment Rate"]] = sum(this_metadata$seurat_clusters %in% unassignment) / cell_number
  return_list[["Cluster outlier Rate"]] = sum(this_metadata$assignment == "cluster_outliers") / cell_number
  return_list[["Mixed Rate"]] = sum(this_metadata$seurat_clusters %in% dup) / cell_number
  cluster_number = length(unique(this_metadata$seurat_clusters))
  unknown_cluster = c(unassignment, dup)
  if(!is.null(prior_cell_type)){
    purity_mask = this_metadata$assignment == this_metadata$cell_type
    purity = sum(purity_mask) / cell_number
    ARI <-adjustedRandIndex(this_metadata$assignment, this_metadata$cell_type)
    recall_type = unique(this_metadata$cell_type[purity_mask])
    return_list[["recall_type"]] = recall_type
    return_list[["lost_type"]] = setdiff(unique(this_metadata$cell_type), recall_type)
    recall_number = length(recall_type)
    result_summary <-c(purity, ARI, cluster_number, recall_number, length(unknown_cluster))
    names(result_summary) = c("Accuracy", "ARI", "Clusters", "Types", "Unknown") #  Accuracy means purity
  }else{
    result_summary <-c(cluster_number, length(unknown_cluster))
    names(result_summary) = c("Clusters", "Unknown") #  Accuracy means purity
    assignment = this_metadata$assignment
    names(assignment) = rownames(this_metadata)
    return_list[["assignment"]] = assignment
  }
  return_list[["result"]] = result_summary
  return(return_list)
}



seurat_clustering = function(expression_path, feature_path = NULL, cell_type_rds = NULL, gene_selection_rds = NULL, output_dir = NULL, pca_dim = 50, res = 1.4, min_pct = 0.1, show = FALSE){
  if(is.null(feature_path)){
    output_dir_candidate = paste(c(delete_last_element(unlist(strsplit(expression_path, "/", fixed = T))), "cluster_evaluation", get_last_element(delete_last_element(unlist(strsplit(expression_path, "[/\\.]", perl = T)))), "pca"), collapse = "/")
  }else{
    cell.embeddings = t(readh5_loom(feature_path, is_feature = TRUE))
    output_dir_candidate = paste(c(delete_last_element(unlist(strsplit(expression_path, "/", fixed = T))), "cluster_evaluation", get_last_element(delete_last_element(unlist(strsplit(expression_path, "[/\\.]", perl = T)))), "feature"), collapse = "/")
  }
  if(is.null(output_dir)){
    output_dir = output_dir_candidate
  }
  gene_bc_mat = get_gene_bc_mat(expression_path)
  if(!is.null(gene_selection_rds)){
    gene_bc_mat = gene_bc_mat[gsub("_", "-", rownames(gene_bc_mat)) %in% readRDS(gene_selection_rds), ] 
  }
  seurat_obj = CreateSeuratObject(as.data.frame(gene_bc_mat))
  seurat_obj <- NormalizeData(object = seurat_obj)
  dir.create(output_dir, showWarnings = F, recursive = T)
  if(is.null(feature_path)){
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(object = seurat_obj, verbose = !show)
    pca_dim_use = min(c(pca_dim * 2, nrow(gene_bc_mat)))
    seurat_obj <- RunPCA(object = seurat_obj, npcs = pca_dim_use)
    tmp_plot = ElbowPlot(seurat_obj, ndims = pca_dim_use)
    if(show){
      print(tmp_plot)
    }
    ggsave(plot = tmp_plot, filename = paste0(output_dir, "/elbow.pdf"), height = 8, width = 11)
  }else{
    seurat_obj@reductions = list()
    seurat_obj@reductions[["pca"]] = Seurat:::DimReduc(cell.embeddings = cell.embeddings, assay.used = "RNA", key = "feature_")
    pca_dim = ncol(cell.embeddings)
  }
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:pca_dim, tsne.method = "Rtsne", check_duplicates = FALSE)
  seurat_obj = RunUMAP(seurat_obj, dims = 1:pca_dim)
  if(!is.null(cell_type_rds)){
    all_cell_type = readRDS(cell_type_rds)
    cell_type_overlap_mask = colnames(gene_bc_mat) %in% cell_type
    cell_type = all_cell_type[cell_type_overlap_mask]
    plot_seurat_obj = subset(seurat_obj, cells = names(cell_type))
    plot_seurat_obj@active.ident = factor(cell_type, levels = sort(unique(cell_type)))
    tmp_plot = TSNEPlot(plot_seurat_obj)
    if(show){
      print(tmp_plot)
    }
    ggsave(plot = tmp_plot, filename = paste0(output_dir, "/tsne_cell_type.pdf"), height = 8, width = 11)
    tmp_plot = UMAPPlot(plot_seurat_obj)
    if(show){
      print(tmp_plot)
    }
    ggsave(plot = tmp_plot, filename = paste0(output_dir, "/umap_cell_type.pdf"), height = 8, width = 11)
  }
  seurat_obj <- FindNeighbors(object = seurat_obj, dims = 1:pca_dim, force.recalc = TRUE)
  seurat_obj <- FindClusters(object = seurat_obj, resolution = res)
  tmp_plot = TSNEPlot(seurat_obj)
  if(show){
    print(tmp_plot)
  }
  ggsave(plot = tmp_plot, filename = paste0(output_dir, "/tsne_cluster.pdf"), height = 8, width = 11)
  tmp_plot = UMAPPlot(seurat_obj)
  if(show){
    print(tmp_plot)
  }
  ggsave(plot = tmp_plot, filename = paste0(output_dir, "/umap_cluster.pdf"), height = 8, width = 11)
  this_metadata = seurat_obj@meta.data
  write.table(this_metadata, paste0(output_dir, "/metadata.txt"), row.names = T, col.names = T, quote = F)
  this_markers <- FindAllMarkers(object = seurat_obj, only.pos = TRUE, min.pct = min_pct, logfc.threshold = 0.25, verbose = !show)
  write.table(this_markers, paste0(output_dir, "/markers.txt"), row.names = T, col.names = T, quote = F)
  write.table(seurat_obj@active.ident, paste0(output_dir, "/cluster_cell_type.txt"), row.names = T, col.names = F, quote = F)
  cluster_cell_type = seurat_obj@active.ident
  save(this_markers, this_metadata, cluster_cell_type, file = paste0(output_dir, "/tmp.rdata"))
}



