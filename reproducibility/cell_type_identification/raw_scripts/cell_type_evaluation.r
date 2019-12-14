args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3){
  stop("R --slave < this_code.r --args <clustering result rdata> <cell type rds> <mapping function>")
}
source("/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/generic_functions.r")
output_dir = paste(delete_last_element(unlist(strsplit(args[1], "/", fixed = T))), collapse = "/")
load(args[1])
switch(args[3],
       "pbmc"={
         mapping_function = cluster_evaluation_pbmc
       },{
         stop("no such mapping")
       })
evaluate_result = cluster_evaluation_pbmc(this_markers, this_metadata, readRDS(args[2]))
print(evaluate_result)
this_index = evaluate_result[["result"]]
print(this_index)
write.table(t(this_index), paste0(output_dir, "/cluster_evaluation_result.txt"), sep = "\t", row.names = F, col.names = T)





