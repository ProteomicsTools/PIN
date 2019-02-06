#' @importFrom data.table fread


#' @export
import_openswath <- function(input_search_results, remove_prefixInFileName = FALSE) {
    
  message("Message: It starts reading raw OpenSWATH search results...")
  
  raw <- fread(input=input_search_results, head=T)

  names(raw)[1] <- "PeptideIon"

  raw$PeptideIon <- paste(sapply(strsplit(raw$PeptideIon, "_"), "[[", 2), "_", sapply(strsplit(raw$PeptideIon, "_"), "[[", 3), sep="")

  #wenguang: in the search results from euler portal, the filename was in the format as "/scratch/71239421.tmpdir/xuep_J180621_SW_3.mzXML.gz". this will remove this unnecessary prefix...
  if(remove_prefixInFileName == TRUE) {
    raw$filename <- sapply(strsplit(raw$filename, "/"), "[[", 4)
  }

  data <- raw[-which(grepl("DECOY", raw$ProteinName)), ]

  summarize_data(data)

#data$filename <- sapply(strsplit(data$filename, "/"), "[[", 4)

  return(data)

} 




#' @export
read_sample_annotation <- function(input_file) {
  
  anno <- data.frame(read.table(file = input_file, sep = "\t", head = T, stringsAsFactors = F))
  
  if (!"sample_id" %in% names(anno) || !"filename" %in% names(anno)) {
    STOP("Please choose a valid sample annotation table that contains at least sample_id and filename columns.")
  }
  
  return(anno)
  
}
