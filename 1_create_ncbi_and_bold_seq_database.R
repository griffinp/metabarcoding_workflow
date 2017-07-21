library(metabarcodedb)
library(stringr)

gene <- "COI"
# some Genbank entries are spelt differently!
gene_alternate_name_1 <- "CO1"
gene_alternate_name_2 <- "COX1"
co1_gene_names <- c(gene, gene_alternate_name_1, gene_alternate_name_2)

#chunk_size <- 50
#block_size <- 2000
#max_records <- 200000

################

# Import and tidy vector of taxa for which we want to get sequences for all species

output_dir <- "~/Documents/Research_Projects/metabarcoding/metaR/ref_files/taxon_list1"
input_dir <- "~/Documents/Research_Projects/metabarcoding/metaR/"

taxon_list1 <- read.table(file=paste(input_dir, "/Taxon_list1.txt", sep=""), sep="\t", stringsAsFactors = FALSE, header=TRUE)[,"mixed_levels"] %>% str_trim()
taxon_list2 <- read.table(file=paste(input_dir, "/Taxon_list2.txt", sep=""), sep="\t", stringsAsFactors = FALSE, header=TRUE)[,"mixed_levels"] %>% str_trim()

for(taxon in taxon_list2){
  output_dir <- "~/Documents/Research_Projects/metabarcoding/metaR/ref_files/taxon_list2"
  position <- which(taxon_list2==taxon)
  all_blocks <- get_ncbi_search_blocks(taxon=taxon, gene_names=co1_gene_names, record_block_size = 600)
  all_ids <- extract_ids_from_ncbi_search_blocks(all_blocks)
  if(length(all_ids[[1]])>0){
    check_and_do_mapping(taxon=taxon, mapping_file_name=make_mapping_name(taxon=taxon, directory=output_dir),
                         id_list=all_ids)
    check_and_do_classification(taxon=taxon, classification_file_name=make_classification_name(taxon=taxon, directory=output_dir),
                                id_list = all_ids, chunk_size=60, block_size=400)
    check_and_get_seqs(taxon=taxon, raw_fasta_file_name = make_raw_fasta_name(taxon=taxon, directory=output_dir, suffix="_COI.fasta"),
                       mapping_file_name = make_mapping_name(taxon=taxon, directory=output_dir),
                       ncbi_record_blocks = all_blocks)

  } else{
  }
}

for(taxon in taxon_list2){
  output_dir <- "~/Documents/Research_Projects/metabarcoding/metaR/ref_files/taxon_list2"
  raw_fasta_file_name <- make_raw_fasta_name(taxon=taxon, directory=output_dir, suffix="_COI.fasta")
  dada2_fasta_file_name <- make_dada2_fasta_name(taxon=taxon, directory=output_dir)
  if(file.exists(raw_fasta_file_name) & file.exists(dada2_fasta_file_name)==FALSE){
    make_dada2_ref_database(taxon=taxon, raw_fasta_file_name=raw_fasta_file_name,
                            mapping_file_name = make_mapping_name(taxon=taxon, directory=output_dir),
                            classification_file_name = make_classification_name(taxon=taxon, directory=output_dir),
                            dada2_fasta_file_name = dada2_fasta_file_name, include.identifier="uid")
  }
}



###############
# Trying BOLD #
###############

output_dir <- "~/Documents/Research_Projects/metabarcoding/metaR/ref_files/taxon_list1"

for(i in 1:length(taxon_list1)){
  taxon <- taxon_list1[i]
  marker <- "COI-5P|COI-3P"
  mapping_file_name <- make_mapping_name(taxon=taxon, directory=output_dir)
  dada2_file_name <- paste(output_dir, "/", taxon, "_extra_bold_seq_dada2.fasta", sep="")
  
  check_and_get_bold_records(taxon=taxon, marker=marker, mapping_file_name=mapping_file_name,
                             dada2_file_name=dada2_file_name)
}



