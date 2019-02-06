
#' @export
generate_peptide_table <- function(search_results, sptxt, sample_annotation, remove_prefixInFileName = FALSE, normalize_intensity = TRUE) {

  lib_info <- generate_lib_info(input_sptxt=sptxt)

  data <- import_openswath(input_search_results=search_results, remove_prefixInFileName = remove_prefixInFileName)

  anno <- read_sample_annotation(input_file=sample_annotation)

  all_peptides <- long2wide(data, level="PeptideIon")

  proteotypic_peptides <- keep_only_proteotypic(all_peptides)

  peptides_extend <- merge(proteotypic_peptides, lib_info, by.x="PeptideIon", by.y="PeptideIon")


#Wenguang: Tested severeal times. This parameter would not change the PIN based results a lot. By turning off the normalization, it will increase the sensitivity of degradation detection (e.g. the range of PINs increases) as expected, which is not the usual case. So by default the intensity would be normalized...
if(normalize_intensity == TRUE) {

  index_intensity <- which( names(peptides_extend) %in% unique(data$filename) )
  peptides_extend_normalized <- normalize_data(peptides_extend, index_intensity)

  cons_peptides <- merge_replicates(peptides_extend_normalized, anno)

} else {

  cons_peptides <- merge_replicates(peptides_extend, anno)

}

  cons_peptides_long <- wide2long(cons_peptides)

  return(cons_peptides_long)

#Wenguang: this will generate a peptide table in "long" format.

}



#' @export
perform_PIN_analysis <- function(peptide_table, intensity_normalization = "sqrt", remove_zero_rows = FALSE) {

  proteins_long <- pept2prot(peptide_table, intensity_normalization = intensity_normalization)
  cons_proteins <- long2wide(proteins_long, level="Protein")

  semi_cons_peptides_long <- keep_only_semi(peptide_table)

  semi_proteins_long <- pept2prot(semi_cons_peptides_long, intensity_normalization = intensity_normalization)
  semi_cons_proteins <- long2wide(semi_proteins_long, level="Protein")

  iPIS <- generate_iPIS_matrix(cons_proteins, semi_cons_proteins)
  write_tsv(iPIS)

  PIN <- calculate_PIN_new(iPIS, remove_zero_rows)
  #PIN <- calculate_PIN(iPIS, remove_zero_rows)
  write_tsv(PIN)

  message("Message: PIN analysis is done.")

}
