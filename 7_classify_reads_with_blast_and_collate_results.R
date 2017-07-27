#install.packages("RCurl", lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")

#source("https://bioconductor.org/biocLite.R")

#biocLite("RCurl", lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")
#biocLite("dada2",
#         lib.loc="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2",
#         lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")

library(dada2)
library(microclass)

#############
# FUNCTIONS #
#############

blastClassify <- function (sequence, bdb) 
  # Custom function for classifying sequences using BLAST
  # adapted from microclass::blastClassify16S to return more than
  # just the top hit
{
  n <- length(sequence)
  tags <- paste("Query", 1:n, sep = "_")
  fdta <- data.frame(Header = tags, Sequence = sequence, stringsAsFactors = F)
  writeFasta(fdta, out.file = "query.fasta")
  cmd <- paste("blastn -task megablast -query query.fasta -db ", bdb, " -num_alignments 400", 
               " -out bres.txt -outfmt \"6 qseqid qlen sseqid length pident bitscore\"", 
               sep = "")
  system(cmd)
  # btab <- read.table("bres.txt", sep = "\t", header = F, stringsAsFactors = F)
  # file.remove("query.fasta")
  # btab <- btab[order(btab[, 5], decreasing = T), ]
  # btab <- btab[which(!duplicated(btab[, 1])), ]
  # tax.hat <- gsub("_[0-9]+$", "", btab[, 3])
  # pident_max <- max(btab[, 5])
  # pident_max_subset <- btab[,which(btab[,5]>(pident_max-0.001))]
  # pident_max_sseqids <- pident_max_subset[seq(3, ncol(pident_max_subset), by=6),]
  #idty <- (btab[, 5]/100) * btab[, 4]/btab[, 2] + pmax(0, btab[, 
  #                                                             2] - btab[, 4]) * 0.25
  # taxon.hat <- rep("unclassified", n)
  # identity <- rep(0, n)
  # idx <- match(btab[, 1], tags)
  # taxon.hat[idx] <- tax.hat
  # #identity[idx] <- idty
  # return(data.frame(Taxon = taxon.hat, Identity = identity, 
  #                   stringsAsFactors = F))
}

##################
# SETTING UP     #
##################


path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/trimmed/A_E"
filt_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/trimmed/A_E/filtered"

ref_blastdb <- readRDS("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/dada2_blast_ref_db.rds")

blast_classif <- blastClassify(sequence=colnames(seqtab.nochim), bdb=ref_blastdb)

saveRDS(blast_classif, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_blast_classification_custom.rds")


#classif <- readRDS("A_E_blast_classification.rds")

classif <- read.csv("bres.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)

classif$seqnumber <- as.numeric(unlist(lapply(strsplit(classif[,1], split="_"), "[[", 2)))

extract_top_pident_hits <- function(blast_result){
  pident_max <- max(blast_result[,5])
  pident_max_subset <- blast_result[which(blast_result[,5] > (pident_max-0.001)),c(1,3,5)]
  return(pident_max_subset)
}

condense_taxa <- function(taxa_vector){
  if(length(taxa_vector)==1){
    output_levels <- taxa_vector
    taxa_split <- strsplit(taxa_vector, split=";")[[1]]
    finest_level <- taxa_split[7]
    if(finest_level%in%c("", "_")){
      finest_level <- taxa_split[max(which(taxa_split[1:7]%in%c("", "_")==FALSE))]
    }
  } else{
    taxa_table <- t(as.data.frame(strsplit(taxa_vector, split=";"), stringsAsFactors = FALSE))[,1:8]
    row.names(taxa_table) <- 1:nrow(taxa_table)
    taxa_table <- unique(taxa_table)
    if(nrow(taxa_table)==1){
      output_levels <- paste(taxa_table, collapse=";")
      finest_level <- taxa_table[,7]
      if(finest_level%in%c("", "_")){
        finest_level <- taxa_table[1,max(which(taxa_table[1,1:7]!=""))]
      }
    } else{
      number_taxon_levels <- apply(taxa_table, MARGIN = 2, FUN=function(x){length(levels(as.factor(x)))})
      lowest_nonvarying_column <- which(number_taxon_levels>1)[1]-1
      nonvarying_levels <- paste(taxa_table[1,1:lowest_nonvarying_column], collapse=";")
      finest_level <- taxa_table[1,lowest_nonvarying_column]
      if(finest_level%in%c("", "_")){
        finest_level <- taxa_table[1,max(which(taxa_table[1,1:lowest_nonvarying_column]%in%c("", "_")==FALSE))]
      }
      if(lowest_nonvarying_column<7){
        varying_levels_pre <- apply(taxa_table[,(lowest_nonvarying_column+1):ncol(taxa_table)], MARGIN=1, FUN=function(x){paste(x, collapse=";")})
        varying_levels <- paste(varying_levels_pre, collapse="|")
        nonassigned_level_length <- length(number_taxon_levels>1)-1
        nonassigned_levels <- rep("", times=nonassigned_level_length)
        output_levels <- paste(nonvarying_levels, varying_levels, nonassigned_levels, sep=";")
      }
      if(lowest_nonvarying_column==7){
        varying_levels <- paste(taxa_table[,(lowest_nonvarying_column+1)], collapse="|")
        output_levels <- paste(nonvarying_levels, varying_levels, sep=";")
      }
    }
  }
  return(c(finest_level, output_levels))
}

info_per_seq <- data.frame(query_seq_number=numeric(), top_pident=numeric(), lowest_unique_match=character(), matching_taxa=character(),
                           stringsAsFactors=FALSE)
for(i in 1:max(classif$seqnumber)){
  classif_subset <- classif[classif$seqnumber==i,]
  if(nrow(classif_subset)==0){
    info_per_seq[i,] <- c(query_seq_number=i, top_pident=0, lowest_unique_match="unidentified", matching_taxa="unidentified")
  } else{
    top_hits <- extract_top_pident_hits(classif_subset)
    message(paste(nrow(top_hits), "hits with identity score of", top_hits[1,3], "for sequence #", i))
    top_taxa <- condense_taxa(top_hits[,2])
    info_per_seq[i,] <- c(query_seq_number=i, top_pident=top_hits[1,3], 
                          lowest_unique_match=top_taxa[1], matching_taxa=top_taxa[2])
  }
}

info_per_seq$sequence <- colnames(seqtab.nochim)
write.csv(info_per_seq, file="A_E_blast_summary.csv", row.names = FALSE, quote=FALSE)

seqtab_assigned <- seqtab.nochim
colnames(seqtab_assigned) <- info_per_seq$lowest_unique_match
seqtab_assigned <- seqtab_assigned[,order(info_per_seq$lowest_unique_match)]
pident_sorted <- info_per_seq[order(info_per_seq$lowest_unique_match), "top_pident"]
seq_sorted <- colnames(seqtab.nochim)[order(info_per_seq$lowest_unique_match)]
matching_taxa_sorted <- info_per_seq[order(info_per_seq$lowest_unique_match), "matching_taxa"]

seqtab_out <- rbind(seqtab_assigned, pident_sorted, seq_sorted, matching_taxa_sorted)
row.names(seqtab_out)[((nrow(seqtab_out)-2):nrow(seqtab_out))] <- c("top_pident", "sequence", "all_top_match_taxa")

write.csv(seqtab_out, file="Seq_abundance_A_E_with_lowest_unique_assignment.csv", 
          quote=FALSE)





# Collating at species level (this ignores the percent identity values though)

collate_at_species_level <- function(info_per_seq, seqtab_assigned){
  #if(exists("species_collated")){rm(species_collated)}
  unique_taxa_matches <- unique(info_per_seq$lowest_unique_match)
  for(i in 1:length(unique_taxa_matches)){
    unique_taxon <- unique_taxa_matches[i]
    cols <- which(colnames(seqtab_assigned)==unique_taxon)
    tempsubset <- seqtab_assigned[,cols]
    if(length(cols)==1|("_"%in%unlist(strsplit(unique_taxon, split=""))==FALSE)){
      if(i==1){
        species_collated <- as.data.frame(tempsubset)
        colnames(species_collated) <- unique_taxon
      } else{
        existing_colnames <- colnames(species_collated)
        species_collated <- cbind(species_collated,tempsubset)
        #print(colnames(tempsubset))
        if(length(cols)==1){
          colnames(species_collated) <- c(existing_colnames, unique_taxon)
        } else{
          new_colnames <- paste(rep(unique_taxon, times=ncol(as.data.frame(tempsubset))),
                                1:ncol(as.data.frame(tempsubset)), sep="_")
          colnames(species_collated) <- c(existing_colnames, new_colnames)
        }
      }
    } else{
      if(i==1){
        species_collated <- as.data.frame(rowSums(tempsubset))
        colnames(species_collated) <- unique_taxon
      } else{
        existing_colnames <- colnames(species_collated)
        species_collated <- cbind(species_collated, rowSums(tempsubset))
        colnames(species_collated) <- c(existing_colnames, unique_taxon)
      }
    }
  }
return(species_collated)
}


matches_above_97.5 <- as.numeric(info_per_seq$top_pident)>97.499
matches_below_97.5 <- as.numeric(info_per_seq$top_pident)<=97.499

species_collated <- collate_at_species_level(info_per_seq=info_per_seq, seqtab_assigned=seqtab_assigned)

species_collated <- species_collated[,order(colnames(species_collated))]
write.csv(species_collated, file="Seq_abundance_A_E_lowest_unique_assignment_grouped_by_species.csv", quote=FALSE)

presabs <- as.matrix((species_collated > 0) + 0)
rowSums(presabs)

pdf("Species_collated_presence_heatmap.pdf", width=7, height=14)
image(presabs, col=c("white", "red"), xaxt = "n", yaxt="n")
axis(side = 1, labels = rownames(presabs),
     at = seq(0, by = 1/(nrow(presabs)-1), length.out = nrow(presabs)), las=3,
     cex.axis=0.3)
axis(side = 2, labels=colnames(presabs),
     at=seq(0, by=1/(ncol(presabs)-1), length.out=ncol(presabs)), las=2,
     cex.axis=0.3)
dev.off()

#taxa_assignment <- assignTaxonomy(seqtab.nochim, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/dada2_fixed_uppercase_nolong.fasta.gz", 
#                       multithread=TRUE, taxLevels = c("Kingdom", "Phylum", "Class",
#                                                       "Order", "Family", "Genus", "Species", "ID"))

#taxa <- taxa_assignment["taxa"]
#bootstraps <- taxa_assignment["boot"]

#write.csv(unname(taxa), file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/taxa_table_A_E.csv")
#saveRDS(taxa, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_taxa.rds")

#write.csv(bootstraps, file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/bootstrap_table_A_E.csv")
#saveRDS(bootstraps, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_bootstraps.rds")


# seqtab <- readRDS("~/Documents/metabarcoding/test_run1_data/A_E_seqtab_nochim.rds")
# 
# presabs <- as.matrix((seqtab > 0) + 0)
# rowSums(presabs)
# 
# image(presabs, col=c("white", "red"), xaxt = "n", yaxt="n")
# axis(side = 1, labels = rownames(presabs),
#      at = seq(0, by = 1/(nrow(presabs)-1), length.out = nrow(presabs)), las=3,
#      cex.axis=0.3)
# 
# 
# temptab <- seqtab
# colnames(temptab) <- rep("", ncol(temptab))
# tempdf <- as.data.frame(temptab[12:13,1:200])


