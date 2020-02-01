utilities_path = "/home/yuanhao/single_cell/scripts/evaluation_pipeline/evaluation/utilities.r"
source(utilities_path)
result_list = list()
method_names = c("Observed", "DISC", "SAVER", "scImpute", "VIPER", "MAGIC", "DCA", "deepImpute", "scScope", "scVI")
for(ii in method_names){
  result_list[[ii]] = list()
  for(jj in c(1, 2, 3, 4, 5)){}
    if(ii %in% c("Observed", "DISC", "SAVER", "MAGIC", "DCA", "scScope", "scVI")){
      if(ii %in% c("SAVER", "MAGIC", "DCA", "scScope", "scVI")){
        use_name = ii
      }else if(ii == "Observed"){
        use_name = "RAW"
      }else if(ii == "DISC"){
        use_name = "DeSCI"
      }else{
        stop("Invaild use_name!")
      }
      cell_type_result = as.matrix(read.table(paste0("/home/wucheng/imputation/pbmc/ds/repeat", jj, "/", use_name, "/Recall1.txt"), sep = "\t"))
      overall_result = as.matrix(read.table(paste0("/home/wucheng/imputation/pbmc/ds/repeat", jj, "/", use_name, "/Accuracy.txt"), sep = "\t", header = T))
      this_result = list(Jaccard = )
    }
}



















