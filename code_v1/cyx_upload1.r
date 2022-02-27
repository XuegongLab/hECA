suppressMessages(suppressWarnings(library(Seurat)))

args <- commandArgs(trailingOnly=TRUE)
query_dir <- args[1]
out_dir <- args[2]
# delete \r
query_dir <- gsub('\r','',query_dir)
out_dir <- gsub('\r','',out_dir)

print(query_dir)
print(out_dir)

as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# load object
dataobj <- readRDS(query_dir)
print("load finished")

# build matrix
if (nrow(dataobj@meta.data) <=40000){
    data.matrix <- as.data.frame(as.matrix(dataobj@assays$RNA@data))
}else{
    data.matrix <- as.data.frame(as_matrix(dataobj@assays$RNA@data))
}
annotation.matrix <- dataobj@meta.data

# write files
write.csv(data.matrix, file=paste0(out_dir,"/UniformedExpression.csv"))
write.csv(annotation.matrix, file=paste0(out_dir,"/Annotation.csv"))