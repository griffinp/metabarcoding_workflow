######################################################################
# This section aims to obtain taxonomic info and format sequences    #
# that Mel Carew obtained for local Australian aquatic invertebrates #
# but that are not yet in GenBank or BOLD                            #
######################################################################

library(rentrez)
library(ape)
library(taxize)
library(metabarcodedb)

#############
# FUNCTIONS #
#############

replace_sp_value <- function(classif, new_sp_value){
  classif[classif$rank=="species","name"] <- new_sp_value
  return(classif)
}

# add_identifier_line <- function(classif, identifier_value){
#   # nb this is now a package function
#   new_classif <- rbind(classif, c(name=identifier_value, rank="identifier", id=""))
#   return(new_classif)
# }

#########################

# Mel's sequences are in the format:
# >Paracalliopiidae(MRD16Amp1)
# CACCCTTTATTTTATTTTAGCCGCATGAGCTAGTATAGTAGGCACTTCTCTAAGAGTTATTATCCGAACAGAATTAAGGGCCCCAGGCAACTTAATCGGAGACGATCAGATTTACAATACCGTAGTAACAGCTCACGCATTCGTTATAATTTTTTTTTTAGTCATACCAGCCATAATTGGTGGTTTTGGTAATTGACTAGTCCCTCTTATACTAGGAAGACCAGATATAGCTTTCCCCCGAATAAACAACATAAGATTCTGACTTTTACCCCCTTCTCTGACCCTATTACTAATAAGAGGCCTAGTAGAAAGAGGGGTGGGCACAGGCTGAACCGTCTACCCGCCTCTTGCTGGCAACATCGCCCATAGTGGAGCCTCAGTAGATCTGGCTATTTTTTCTCTTCATTTAGCAGGAGCCTCCTCTATTTTAGGCGCTATCAATTTTATCTCCACAGTAATTAACATGCGGGCCCCGGCCATGCCAATAGACCAAATCCCATTATTTGTATGGTCTGTATTTATCACAGCCATCTTACTATTATTATCTTTACCTGTTTTAGCCGGTGCAATTACCATACTACTTACAGACCGTAATTTAAACACATCTTTTTTTGACCCATCTGGAGGAGGGGACCCAATTTTATACCAGCATTTATTT
# >Physa_acuta(3aPhy1)
# AACATTGTATTTAATTTTTGGGATTTGGTGTGGCTTGGTCGGTACAGGCTTAAGCTTGTTAATTCGTTTGGAATTAGGAACATCTCTGGTACTATTGGATGAACATTTTTATAATGTAATCGTTACAGCACATGCTTTTGTAATGATTTTTTTTATAGTTATACCTATAATAATTGGAGGCTTTGGGAATTGAATAGTACCTATATTAATTGGTGCTCCCGATATAAGCTTTCCTCGAATAAATAATATAAGATTTTGACTTTTACCCCCTTCATTTATCTTATTATTATGTAGGTCTATAGTTGAGGGTGGGGTAGGAACTGGGTGAACTGTTTATCCTCCACTATCAGGACCTGTAGCTCACTCTGGTTCATCAGTAGATCTTGCTATTTTCTCATTACACTTAGCTGGGTTATCATCTATTTTAGGTGCTATTAATTTTATTACTACAATTTTTAATATACGTTCTCCTGGTATTACACTGGAACGAATAAGTTTATTTGTTTGATCAGTGTTAATTACTGCATTTTTATTATTATTGTCATTGCCTGTTTTAGCGGGGGCTATTACTATACTATTAACTGATCGAAATTTTAATACTAGGTTCTTTGATCCAAGAGGGGGGGGAGACCCTATTCTATATCAACATCTATTT
# >Ablabesmyia_sp1(S04Tp11)
# AACTTTATATTTTATTTTTGGGGCCTGAGCTGGAATAGTGGGTACTTCCCTTAGTATCCTTATTCGAACAGAATTAGGACACCCAGGAGCTTTAATCGGAGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCATTTGTTATAATTTTTTTTATAGTAATACCTATTTTAATTGGTGGATTTGGAAATTGACTAGTACCCCTTATACTAGGTGCCCCAGATATAGCATTTCCACGAATAAATAATATAAGATTTTGACTACTTCCCCCCTCATTAACTTTATTGTTATCTAGTTCTATTGTAGAAAATGGAGCAGGAACTGGTTGAACCGTTTACCCCCCTTTAGCCTCAGGTATTGCCCATGCCGGAGCTTCTGTAGATTTAGCTATTTTCTCTCTTCATTTAGCAGGAATCTCTTCAATTTTAGGAGCTGTAAATTTTATTACTACAGTTATTAATATACGATCTACAGGAATTACATTAGACCGAATACCCTTATTTGTTTGGTCTGTTGTAATTACTGCTATTTTATTGCTTTTATCCCTACCAGTTTTAGCTGGGGCTATTACTATATTATTGACAGATCGAAATTTAAATACTTCTTTCTTTGACCCCGCAGGAGGGGGAGACCCCATTTTATACCAGCACTTATTT

#######################################
# Importing Mel's sequences
#######################################


output_dir <- "~/Documents/Research_Projects/metabarcoding/metaR/ref_files/private_ref_sequences"

private_seq_1 <- ape::read.dna(file="~/Documents/Research_Projects/metabarcoding/mel_data/Private_DNA_barcodes_11_16.fasta", format = "fasta")
original_names_1 <- names(private_seq_1)

# Make some name edits as per instructions in the document 'DNAbarcades1234_adjustments.txt'

name_changes_1 <- data.frame(current_names=c("Ceratopogonidae_EPAsp15(cer15103yar3421)", "Ceratopogonidae_EPAsp15(cer15103yar3582)",
                                           "Ceratopogonidae_EPAsp17(cer17s1bdg1)", "Ceratopogonidae_EPAsp17(cer17s1dan1)",
                                           "Ceratopogonidae_EPAsp17(cer17s1dpw1)", "Ceratopogonidae_EPAsp17(cer17s1ler0150)",
                                           "Ceratopogonidae_EPAsp17(cer17s1n7d1)", "Ceratopogonidae_EPAsp3(cer3s1cdc003)",
                                           "Ceratopogonidae_EPAsp48(cer48s1ple1)", "Ceratopogonidae_EPAsp50(cer50pleghd1)",
                                           "Ceratopogonidae_EPAsp50(cer50s1bdg1)", "Dinotoperia_thwaitesi(SC15ASGrip1)",
                                           "Hellyethira_simplex(NSHydp1)", "Nousia_AV2(2aLept)", "Telephlebiidae_sp1(MAR2_Tele1)"),
                           replacement_names=c("Ceratopogonidae_sp15(cer15103yar3421)", "Ceratopogonidae_sp5(cer15103yar3582)", 
                                               "Ceratopogonidae_sp4(cer17s1bdg1)", "Ceratopogonidae_sp2(cer17s1dan1)", 
                                               "Ceratopogonidae_sp9(cer17s1dpw1)", "Ceratopogonidae_sp3(cer17s1ler0150)", 
                                               "Ceratopogonidae_sp2(cer17s1n7d1)", "Ceratopogonidae_sp16(cer3s1cdc003)", 
                                               "Ceratopogonidae_sp8(cer48s1ple1)", "Ceratopogonidae_sp16(cer50pleghd1)", 
                                               "Ceratopogonidae_sp10(cer50s1bdg1)", "Plecoptera_sp1(SC15ASGrip1)", 
                                               "Hellyethira_spACC4813 (NSHydp1)", "Nousia_spAAX5671 (2aLept)", 
                                               "Austroaeschna_pinheyi(MAR2_Tele1)"), stringsAsFactors=FALSE)

for(i in 1:nrow(name_changes_1)){
  current_name <- name_changes_1[i, "current_names"]
  replacement_name <- name_changes_1[i, "replacement_names"]
  names(private_seq_1)[which(names(private_seq_1)==current_name)] <- replacement_name
  #print(c(current_name, replacement_name))
}

private_seq_2 <- ape::read.dna(file="~/Dropbox/mel_metabarcoding/DNAbarcodes1234_282\ Sequences.fasta", format = "fasta")

private_seq <- c(private_seq_1, private_seq_2)

# look for duplicate sequence names and remove duplicates

duplicates_to_remove <- which(duplicated(names(private_seq)))
private_seq <- private_seq[-duplicates_to_remove]

#######################################
# making some name edits on the raw   #
# file as per Mel's suggestions, and  #
# extracting the lowest possible      #
# taxonomic level                     #
# for use in taxonomy searching       #
#######################################

private_seq_taxon_name <- str_split(names(private_seq), pattern="[()]", simplify=TRUE)[,1]
private_seq_codes <- str_split(names(private_seq), pattern="[()]", simplify=TRUE)[,2]

private_seq_taxon_name_split <- str_split(private_seq_taxon_name, pattern="_", simplify=TRUE)
private_seq_taxon_name_split <- apply(private_seq_taxon_name_split, MARGIN = 2, FUN=str_trim)

# Make some spelling corrections based on Mel's input in the document 
# 'no_taxonomy_hits_list_170227_corrections.txt' and issues identified downstream

private_seq_taxon_name_split[8,1] <- "Ablabesmyia"
private_seq_taxon_name_split[10,1] <- "Botryocladius"
private_seq_taxon_name_split[12,1] <- "Chironomus"
private_seq_taxon_name_split[78:80,2] <- "februarius"
private_seq_taxon_name_split[115:118,2] <- "pseudoppositus"
private_seq_taxon_name_split[168:169,2] <- "bilinearis"
private_seq_taxon_name_split[198:202,1] <- "Corynoneura"
private_seq_taxon_name_split[287,1] <- "Cricotopus"
private_seq_taxon_name_split[c(288, 321),1] <- "Cryptochironomus"
private_seq_taxon_name_split[367,2] <- "spK4"
private_seq_taxon_name_split[419,2] <- "spS04"
private_seq_taxon_name_split[434,1] <- "Parachironomus"
private_seq_taxon_name_split[444:447,1] <- "Parakiefferiella"
private_seq_taxon_name_split[470:471,2] <- "grimmii"
private_seq_taxon_name_split[514,2] <- "kathleenae"
private_seq_taxon_name_split[563:564,2] <- "griseoguttatum"
private_seq_taxon_name_split[569:570,2] <- "spM1"
private_seq_taxon_name_split[c(583:595,604:611),2] <- "spC"
private_seq_taxon_name_split[596:603,2] <- "spE"
private_seq_taxon_name_split[685,2] <- "villosimanus"
private_seq_taxon_name_split[732,2] <- "fuscithorax"
private_seq_taxon_name_split[786:787,1] <- "Cloeon"
private_seq_taxon_name_split[c(788:790, 848:849, 924:927),1] <- "Paracalliopidae"
private_seq_taxon_name_split[c(800, 803),1] <- "Austrochiltonia"
private_seq_taxon_name_split[821:824,1] <- "Dinotoperla"
private_seq_taxon_name_split[821,2] <- "thwaitesi"
private_seq_taxon_name_split[830:832,1] <- "Ischnura"
private_seq_taxon_name_split[841:842,1] <- "Necterosoma"
private_seq_taxon_name_split[843,2] <- ""
private_seq_taxon_name_split[861:862,1] <- "Tamasia"
private_seq_taxon_name_split[841:842,2] <- "penicillatum"
private_seq_taxon_name_split[867,1] <- "Veliidae"
private_seq_taxon_name_split[885,1] <- "Simulium"
private_seq_taxon_name_split[886:895,1] <- "Scirtidae"
private_seq_taxon_name_split[896,1] <- "Scirtes"
private_seq_taxon_name_split[954,1] <- "Notalina"
private_seq_taxon_name_split[957,1] <- "Tubificidae"
private_seq_taxon_name_split[981,1] <- "Hydrachnidae"
private_seq_taxon_name_split[995,1] <- "Eusiridae"
private_seq_taxon_name_split[1037,1] <- "Ceratopogonidae"


new_spp <- grep(private_seq_taxon_name_split[,2], pattern="sp[0-9A-Z]")
more_new_spp <- grep(private_seq_taxon_name_split[,2], pattern="sp$")

corrected_taxon_name <- apply(private_seq_taxon_name_split, MARGIN = 1, FUN = paste, collapse=" ") %>% str_trim()
species_level <- grep(pattern="\\s",corrected_taxon_name)


private_seq_taxon_search_string <- rep(NA, times=nrow(private_seq_taxon_name_split))

for(i in 1:nrow(private_seq_taxon_name_split)){
  if(i%in%c(new_spp, more_new_spp)){
    private_seq_taxon_search_string[i] <- private_seq_taxon_name_split[i,1]%>%str_trim()
  } else{
    private_seq_taxon_search_string[i] <- paste(private_seq_taxon_name_split[i,1]%>%str_trim(), private_seq_taxon_name_split[i,2]%>%str_trim(), sep=" ")
  }
}


########################################################
# Use parsed taxon name to search for taxonomic ID in  #
# NCBI taxonomy database                               #
########################################################

no_taxonomy_hits <- c()
taxon_id_info <- data.frame(taxon_id=rep(NA, times=length(private_seq_taxon_search_string)),
                            species_string=rep(NA, times=length(private_seq_taxon_search_string)))
taxonomy_info <- list()
for(i in 1:length(private_seq_taxon_search_string)){
  tax_test <- entrez_search(db="taxonomy", term=private_seq_taxon_search_string[i])
  message(paste("Getting taxon id for", private_seq_taxon_search_string[i], ": record", i, "of", length(private_seq_taxon_search_string)))
  if(length(tax_test$ids)>0){
    taxon_id_info$taxon_id[i] <- tax_test$ids
  } else{
    # Try just genus for species-level names returning no hits
    genus_species <- str_split(private_seq_taxon_search_string[i], pattern = " ", simplify = TRUE)
    if(length(genus_species) > 1){
      tax_test <- entrez_search(db="taxonomy", term=genus_species[1])
      taxon_id_info$species_string[i] <- private_seq_taxon_search_string[i]
      if(length(tax_test$ids)>0){
        taxon_id_info$taxon_id[i] <- tax_test$ids
      }
      else{
        no_taxonomy_hits <- c(no_taxonomy_hits, paste("Record #", i, private_seq_taxon_search_string[i]))
      }
    }
    else{
      no_taxonomy_hits <- c(no_taxonomy_hits, paste("Record #", i, private_seq_taxon_search_string[i]))
    }
  }
}


obtained_ids <- taxon_id_info$taxon_id[!is.na(taxon_id_info$taxon_id)]

##############################################
# Obtain full classification info            #
# from NCBI Taxonomy database                #
# for all taxa for which a TaxonID was found #
##############################################

classifications <- list()
for(i in 1:nrow(taxon_id_info)){
  message(paste("Getting full classification for", private_seq_taxon_search_string[i], ": record", i, "of", nrow(taxon_id_info)))
  if(is.na(taxon_id_info$taxon_id[i])==FALSE){
    classifications[[i]] <- classification(taxon_id_info$taxon_id[i], db="ncbi")
  } else{
    classifications[[i]] <- NA
  }
}

# Manually make classifications for the taxa in no_taxonomy_hits:
# Byrrocryptus and Megacentron erebeum
# check vector locations for these...
no_taxonomy_hits
private_seq_taxon_search_string[which(is.na(taxon_id_info$taxon_id))]
which(is.na(taxon_id_info$taxon_id))

classifications[[804]] <- list(data.frame(name=c("cellular organisms",
                                                 "Eukaryota","Opisthokonta", "Metazoa",
                                                 "Eumetazoa", "Bilateria",
                                                 "Protostomia", "Ecdysozoa",
                                                 "Panarthropoda", "Arthropoda",
                                                 "Insecta", "Coleoptera", "Polyphaga",
                                                 "Byrrhoidea", "Ptilodactylidae", "Byrrocryptus", ""),
                                          rank=c("no rank", "superkingdom", "no rank", 
                                                 "kingdom", "no rank", "no rank", "no rank",
                                                 "no rank", "no rank", "phylum", "class", "order",
                                                 "suborder", "superfamily", "family", "genus", "species"),
                                          id=as.character(c(131567,2759,33154,33208,6072,33213,
                                                            33317,1206794,88770,6656,50557,7041,41084,107942,186984,NA,NA)), stringsAsFactors=FALSE))
classifications[[805]] <- classifications[[804]]
classifications[[806]] <- classifications[[804]]
classifications[[1038]] <- classifications[[804]]
classifications[[1039]] <- classifications[[804]]
classifications[[1040]] <- classifications[[804]]
classifications[[1041]] <- classifications[[804]]
classifications[[1042]] <- classifications[[804]]
classifications[[1043]] <- classifications[[804]]
classifications[[1044]] <- classifications[[804]]
classifications[[1045]] <- classifications[[804]]
classifications[[837]] <- list(data.frame(name=c("cellular organisms",
                                                 "Eukaryota","Opisthokonta", "Metazoa",
                                                 "Eumetazoa", "Bilateria",
                                                 "Protostomia", "Ecdysozoa",
                                                 "Panarthropoda", "Arthropoda",
                                                 "Insecta", "Diptera", "Culicomorpha",
                                                 "Chironomoidea", "Chironomidae", "Megacentron", "Megacentron erebeum"),
                                          rank=c("no rank", "superkingdom", "no rank", 
                                                 "kingdom", "no rank", "no rank", "no rank",
                                                 "no rank", "no rank", "phylum", "class", "order", "infraorder",
                                                 "superfamily", "family", "genus", "species"),
                                          id=as.character(c(131567,2759,33154,33208,6072,33213,
                                                            33317,1206794,88770,6656,50557,7147,43786,41828,7149,NA,NA)), stringsAsFactors=FALSE))

############################################
# Edit classifications with functions from #
# metabarcodedb package to retain only     #
# desired taxonomic levels                 #
############################################

library(metabarcodedb)
classifs <- lapply(classifications, "[[", 1)
reduced_classifications <- lapply(classifs, reduce_taxonomic_levels, taxonomic_levels=c("kingdom", "phylum", "class", "order",
                                                                                        "family", "genus", "species"))
# replace species value with species string from earlier, but only for species-level IDs
for(i in species_level){
  reduced_classifications[[i]] <- replace_sp_value(reduced_classifications[[i]], new_sp_value=corrected_taxon_name[i])
}

# add identifier code as key-value pair

for(i in 1:length(reduced_classifications)){
  reduced_classifications[[i]] <- add_identifier_line_to_classification(classif=reduced_classifications[[i]],
                                                      identifier_value=paste("mc", private_seq_codes[i], sep=":"))
}

formatted_classifications <- lapply(reduced_classifications, format_taxonomic_levels)
new_fasta <- private_seq
names(new_fasta) <- formatted_classifications
dada2_fasta_file_name <- paste(output_dir, "private_dna_seq_2017-07-20_dada2.fasta", sep="/")
ape::write.dna(new_fasta, file = dada2_fasta_file_name, format = "fasta", nbcol=1, blocksep=0, colw=1e06)

# replace any blank lines...

system(command=paste("awk 'NF > 0' <", dada2_fasta_file_name, ">", paste(dada2_fasta_file_name, "_2", sep="")))
system(command=paste("mv", paste(dada2_fasta_file_name, "_2", sep=""), dada2_fasta_file_name))

# replace any weird characters

system(command=paste("sed -E 's/(^>.*$)/#\\1#/'", dada2_fasta_file_name,
                     "| tr -d '\r' | tr -d '\n' | sed -e 's/$/#/' | tr '#' '\n' | sed -e '/^$/d' >", paste(dada2_fasta_file_name, "_2", sep="")))
system(command=paste("mv", paste(dada2_fasta_file_name, "_2", sep=""), dada2_fasta_file_name))
system(command=paste("sed '/^>/ s/$/;/'", dada2_fasta_file_name, "| sed 's/; ;/;;/g' | sed 's/ /_/g' >", paste(dada2_fasta_file_name, "_2", sep="")))
system(command=paste("mv", paste(dada2_fasta_file_name, "_2", sep=""), dada2_fasta_file_name))
system(command=paste("awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }'", dada2_fasta_file_name, ">", paste(dada2_fasta_file_name, "_2", sep="")))
system(command=paste("mv", paste(dada2_fasta_file_name, "_2", sep=""), dada2_fasta_file_name))

