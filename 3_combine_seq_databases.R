#################################
# Concatenating all dada2 files #
#################################

# concatenate all files containing 'dada2' into one
# nb BOLD dada2 files have " " for missing fields whereas others have ""
# (not sure why - tried to fix this)
# not sure if this will be a problem for dada2?
# no, doesn't appear to cause any problems

base_dir <- "~/Documents/Research_Projects/metabarcoding/metaR/ref_files"
search_string_1 <- paste(base_dir, "/taxon_list1/*dada2*", sep="")
search_string_2 <- paste(base_dir, "/taxon_list2/*dada2*", sep="")
search_string_3 <- paste(base_dir, "/private_ref_sequences/*dada2*", sep="")
search_string <- paste(c(search_string_1, search_string_2, search_string_3), collapse=" ")
dada2_matches <- system(command=paste('ls', search_string), intern=TRUE)
collapsed <- paste(dada2_matches, collapse=" ")
system(command=paste("cat", collapsed, "> ~/Documents/Research_Projects/metabarcoding/metaR/ref_files/dada2_ref_db_2017-07-25.fasta"))
system(command=paste("gzip", "~/Documents/Research_Projects/metabarcoding/metaR/ref_files/dada2_ref_db_2017-07-25.fasta"))
