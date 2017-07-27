##################################################################################
# This step requires prior installation of the BLAST+ software on your system.   #
# It can be obtained from                                                        #
# https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download #
##################################################################################

library(microclass)

database_directory <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/ref_db"
dada2_database_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/ref_db/dada2_ref_db_2017-07-25.fasta"
dada2_database_path_gz <- paste(dada2_database_path, ".gz", sep="")
blast_database_name <- "custom_blast_db_2017-07-25"


if(!file.exists(dada2_database_path) & !file.exists(dada2_database_path_gz)){ 
  message(paste("Cannot find FASTA file", dada2_database_path, "to convert to BLAST database"))
} else{
    if(file.exists(dada2_database_path_gz) & !file.exists(dada2_database_path)){
    system(command=paste("gunzip", dada2_database_path))
    }
  refdb <- readFasta(dada2_database_path)
  message(paste("Making reference BLAST database", blast_database_name))
  ref_blastdb <- blastDbase16S(name=blast_database_name, sequence=refdb$Sequence, taxon=refdb$Header)
  saveRDS(ref_blastdb, file=paste(database_directory, blast_database_name, ".rds", sep=""))
  system(command=paste("mv ", blast_database_name, "* ", database_directory, sep=""))
  message(paste("Finished making", blast_database_name))
}
