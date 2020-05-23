utilities_path = "./source/utilities.r"
source(utilities_path)
### Trajectory evaluation
#  The trajectory performance using monocle following these scripts(https://github.com/Winnie09/imputationBenchmark/blob/93f27e890a86fdc732257a4036bf38a52faf9f33/trajectory/code/hca/monocle2/01_get_score.R, https://github.com/Winnie09/imputationBenchmark/blob/93f27e890a86fdc732257a4036bf38a52faf9f33/trajectory/code/hca/tscan/01_get_score.R).
###  monocle2
library(monocle)
get_cds_monocle2 = function(gene_bc_mat){
  #  Make a new CDS and use DDRTree for dimension reduction.
  pd = new("AnnotatedDataFrame", data = data.frame(row.names = colnames(gene_bc_mat), cell = colnames(gene_bc_mat)))
  fd = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(gene_bc_mat), gene_short_name = rownames(gene_bc_mat)))
  cds = newCellDataSet(gene_bc_mat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  print("Reducing dimension...")
  cds = reduceDimension(cds)
  return(orderCells(cds))
}
get_score_monocle2 = function(cds, cell_type, correct_order, wrong_order, output_dir = NULL, type_level = NULL, show_interactively = FALSE){
  print("Looking for the root state...")
  used_cells = as.character(pData(cds)$cell)
  if(!is.null(output_dir) | show_interactively){
    if(!is.null(output_dir)){
      dir.create(output_dir, recursive = T, showWarnings = F)
    }
    p = plot_cell_trajectory(cds, color_by = "State")
    if(!is.null(output_dir)){
      ggsave(paste0(output_dir, "/state.pdf"), p)
    }
    if(show_interactively){
      print(p)
    }
    pData(cds)$CellType = cell_type[used_cells]
    p = plot_cell_trajectory(cds, color_by = "CellType")
    if(!is.null(output_dir)){
      ggsave(paste0(output_dir, "/celltype.pdf"), p)
    }
    if(show_interactively){
      print(p)
    }
    if(!is.null(type_level)){
      pData(cds)$Level = type_level[cell_type[used_cells]]
      p = plot_cell_trajectory(cds, color_by = "Level")
      if(!is.null(output_dir)){
        ggsave(paste0(output_dir, "/level.pdf"), p)
      }
      if(show_interactively){
        print(p)
      }
    }
  }
  cell_states = as.numeric(as.character(pData(cds)$State))
  names(cell_states) = used_cells
  unique_states = unique(cell_states)
  checkroot = sapply(unique_states, function(x){
    cds = orderCells(cds, root_state = x)
    return(length(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell))
  })
  candidate_root_states = sort(unique_states[checkroot > 0])
  result_mat = matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("pair_number", "acc")))
  for(ii in candidate_root_states){
    cds = orderCells(cds, root_state = ii)
    this_output_dir = paste0(output_dir, "/rootstate_", ii)
    dir.create(this_output_dir, recursive = T, showWarnings = F)
    if(!is.null(output_dir) | show_interactively){
      p = plot_cell_trajectory(cds, color_by = "Pseudotime")
      if(!is.null(output_dir)){
        ggsave(paste0(this_output_dir, "/pseudotime.pdf"), p)
      }
      if(show_interactively){
        print(p)
      }
    }
    all_branch_points = cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    if(length(all_branch_points) > 0){
      for(jj in seq(length(all_branch_points))){
        cds_tmp = cds
        tryCatch({
          cds_reduced = buildBranchCellDataSet(cds_tmp, branch_point = jj)
          df = data.frame(pData(cds_reduced),stringsAsFactors = F)[used_cells, ]
          if(!is.null(output_dir) | show_interactively){
            pData(cds_tmp)$Pseudotime = df[, "Pseudotime"]
            pData(cds_tmp)$Branch = df[, "Branch"]
            pData(cds_tmp)$State = df[, "State"]
            p = plot_cell_trajectory(cds_tmp, color_by = "Pseudotime")
            if(!is.null(output_dir)){
              ggsave(paste0(this_output_dir, "/branchpoint_", jj, "_pseudotime.pdf"), p)
            }
            if(show_interactively){
              print(p)
            }
            p = plot_cell_trajectory(cds_tmp, color_by = "Branch")
            if(!is.null(output_dir)){
              ggsave(paste0(this_output_dir, "/branchpoint_", jj, "_branch.pdf"), p)
            }
            if(show_interactively){
              print(p)
            }
            p = plot_cell_trajectory(cds_tmp, color_by = "State")
            if(!is.null(output_dir)){
              ggsave(paste0(this_output_dir, "/branchpoint_", jj, "_state.pdf"), p)
            }
            if(show_interactively){
              print(p)
            }
            pData(cds_tmp)$CellType = df[, "CellType"]
            p = plot_cell_trajectory(cds_tmp, color_by = "CellType")
            if(!is.null(output_dir)){
              ggsave(paste0(this_output_dir, "/branchpoint_", jj, "_celltype.pdf"), p)
            }
            if(show_interactively){
              print(p)
            }
            if(!is.null(type_level)){
              pData(cds_tmp)$Level = df[, "Level"]
              p = plot_cell_trajectory(cds_tmp, color_by = "Level")
              if(!is.null(output_dir)){
                ggsave(paste0(this_output_dir, "/branchpoint_", jj, "_level.pdf"), p)
              }
              if(show_interactively){
                print(p)
              }
            }
          }
          df = df[order(df$Pseudotime), ]
          score = rowSums(sapply(unique(df$Branch),function(x){
            branch_cell = as.character(df[df$Branch == x, 1])
            branch_celltype = cell_type[branch_cell]
            index_pair = combn(length(branch_cell), 2)
            if(min(index_pair[2,] - index_pair[1,]) < 0){
              stop("index_pair error")
            }
            branch_cellorder = sprintf("%s_%s",branch_celltype[index_pair[1, ]], branch_celltype[index_pair[2, ]])
            return(c(sum(branch_cellorder %in% correct_order), sum(branch_cellorder %in% wrong_order)))
          }))
          pair_number = sum(score)
          acc = score[1] / pair_number
          this_branch_point_results = matrix(c(pair_number, acc), nrow = 1)
          rownames(this_branch_point_results) = paste0("RS_", ii, "_BP_", jj)
          result_mat = rbind(result_mat, this_branch_point_results)
        }, error = function(e){
          cat(ii, " - ", jj, "\n")
          print(e)
        })
      }
    }
  }
  return(result_mat)
}

###
cell_level_df = data.frame(level = c(1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5), 
                           immunepath = c(1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0), 
                           monopath = c(1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0), 
                           erypath = c(1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1), 
                           stringsAsFactors = F)
rownames(cell_level_df) = c("HSC", "MPP", "LMPP", "CMP", "CLP", "GMP", "MEP", "Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "Mono", "Ery")
path_name = c("immunepath", "monopath", "erypath")
order_list = list(correct = list(), wrong = list(), cell_type = list())
for(ii in path_name){
  path_celltype = rownames(cell_level_df)[cell_level_df[, ii] == 1]
  order_list[["cell_type"]][[ii]] = path_celltype
  cell_type_pair = as.matrix(apply(expand.grid(path_celltype, path_celltype), 2, as.character))
  cell_type_pair = cell_type_pair[cell_type_pair[, 1] != cell_type_pair[, 2], ]
  correct_mat = cell_type_pair[cell_level_df[cell_type_pair[, 1], "level"] < cell_level_df[cell_type_pair[, 2], "level"], ]
  wrong_mat = cell_type_pair[cell_level_df[cell_type_pair[, 1], "level"] > cell_level_df[cell_type_pair[, 2], "level"], ]
  order_list[["correct"]][[ii]] = apply(correct_mat, 1, paste, collapse = "_")
  order_list[["wrong"]][[ii]] = apply(wrong_mat, 1, paste, collapse = "_")
}
type_level = as.character(cell_level_df[, "level"])
names(type_level) = rownames(cell_level_df)
correct_order_all = unique(unlist(order_list[["correct"]]))
wrong_order_all = unique(unlist(order_list[["wrong"]]))
