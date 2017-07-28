library(microclass)

#############
# FUNCTIONS #
#############

blastClassify <- function (sequence, bdb, bdb_path, output_file) 
  # Custom function for classifying sequences using BLAST
  # adapted from microclass::blastClassify16S to return more than
  # just the top hit
{
  setwd(bdb_path)
  n <- length(sequence)
  tags <- paste("Query", 1:n, sep = "_")
  fdta <- data.frame(Header = tags, Sequence = sequence, stringsAsFactors = F)
  writeFasta(fdta, out.file = "query.fasta")
  cmd <- paste("blastn -task megablast -query query.fasta -db ", bdb, " -num_alignments 400", 
               " -out ", output_file, " -outfmt \"6 qseqid qlen sseqid length pident bitscore\"", 
               sep = "")
  system(cmd)
  if(file.exists("query.fasta")){
    file.remove("query.fasta")
  }
}

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
      lowest_nonvarying_column <- tail(which(number_taxon_levels>1), 1)-1
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

collate_at_species_level <- function(info_per_seq, seqtab_assigned){
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

##################
# SETTING UP     #
##################

input_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/trimmed"
results_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/results"
ref_db_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/ref_db"
adapter_table <- read.csv("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/scripts/primer_sequences_for_R_script_run2.csv", sep=",", 
                          header=TRUE, stringsAsFactors=FALSE)
adapter_pairs <- adapter_table$Pair_name
ref_blastdb <- readRDS("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/ref_db/custom_blast_db_2017-07-25.rds")

#################################
# Loop through all primer pairs #
#################################

for(pp in 1:length(adapter_pairs)){
  primer_pair <- adapter_pairs[pp] 
  
  output_file <- paste(results_path, "/", primer_pair, "_blast_classification.txt", sep="")
  input_file <- paste(results_path, "/", primer_pair, "_seqtab_nochim.rds", sep="")
  input_tab <- readRDS(input_file)
  classification_summary_path <- paste(results_path, "/", primer_pair, "_blast_summary.csv", sep="")
  lowest_unique_classification_path <- paste(results_path, "/", primer_pair, "_lowest_unique_classification.csv", sep="")
  classification_collated_by_species_path <- paste(results_path, "/", primer_pair, "_classification_collated_by_species.csv", sep="")
  classification_collated_by_species_plot_path <- paste(results_path, "/", primer_pair, "_classification_collated_by_species_presabs.pdf", sep="")
  
  #############################
  # BLAST classification step #
  #############################
  
  
  if(file.exists(output_file)){
    classif <- read.csv(output_file, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  } else{
    message(paste("BLAST classification in progress for", input_file))
    blastClassify(sequence=colnames(input_tab), bdb=ref_blastdb, bdb_path=ref_db_path, output_file=output_file)
    classif <- read.csv(output_file, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  }
  
  #################################
  # Extract some info from        #
  # the BLAST classification file #
  #################################
  
  classif$seqnumber <- as.numeric(unlist(lapply(strsplit(classif[,1], split="_"), "[[", 2)))
  
  if(!file.exists(classification_summary_path)){
    message(paste("Summarising BLAST classification results for", input_file))
    classification_summary <- data.frame(query_seq_number=numeric(), top_pident=numeric(), lowest_unique_match=character(), matching_taxa=character(),
                                         stringsAsFactors=FALSE)
    for(i in 1:max(classif$seqnumber)){
      classif_subset <- classif[classif$seqnumber==i,]
      if(nrow(classif_subset)==0){
        classification_summary[i,] <- c(query_seq_number=i, top_pident=0, lowest_unique_match="unidentified", matching_taxa="unidentified")
      } else{
        top_hits <- extract_top_pident_hits(classif_subset)
        message(paste(nrow(top_hits), "hits with identity score of", top_hits[1,3], "for sequence #", i))
        top_taxa <- condense_taxa(top_hits[,2])
        classification_summary[i,] <- c(query_seq_number=i, top_pident=top_hits[1,3], 
                                        lowest_unique_match=top_taxa[1], matching_taxa=top_taxa[2])
      }
    }
    classification_summary$sequence <- colnames(input_tab)
    write.csv(classification_summary, file=classification_summary_path, row.names = FALSE, quote=FALSE)
  } else{
    classification_summary <- read.csv(classification_summary_path, stringsAsFactors=FALSE)
  }
  
  ###################################
  # Create and output table of      #
  # sequence abundance at           #
  # lowest unique taxon assignment, #
  # or import if already exists     #
  ###################################
  
  if(!file.exists(lowest_unique_classification_path)){
    message(paste("Collating BLAST classification to lowest unique taxon for", input_file))
    seqtab_assigned <- input_tab
    colnames(seqtab_assigned) <- classification_summary$lowest_unique_match
    seqtab_assigned <- seqtab_assigned[,order(classification_summary$lowest_unique_match)]
    pident_sorted <- classification_summary[order(classification_summary$lowest_unique_match), "top_pident"]
    seq_sorted <- colnames(input_tab)[order(classification_summary$lowest_unique_match)]
    matching_taxa_sorted <- classification_summary[order(classification_summary$lowest_unique_match), "matching_taxa"]
    
    seqtab_out <- rbind(seqtab_assigned, pident_sorted, seq_sorted, matching_taxa_sorted)
    row.names(seqtab_out)[((nrow(seqtab_out)-2):nrow(seqtab_out))] <- c("top_pident", "sequence", "all_top_match_taxa")
    
    write.csv(seqtab_out, file=lowest_unique_classification_path, 
              quote=FALSE)
  } else{
    seqtab_out <- read.csv(lowest_unique_classification_path, stringsAsFactors = FALSE,
                           header=TRUE)
  }
  
  #############################
  # Collate at species level  #
  # (this ignores the percent #
  # identity values though)   #
  #############################
  
  message(paste("Collating at species level for", input_file))
  
  matches_above_97.5 <- as.numeric(classification_summary$top_pident)>97.499
  matches_below_97.5 <- as.numeric(classification_summary$top_pident)<=97.499
  
  seqtab_assigned <- input_tab
  colnames(seqtab_assigned) <- classification_summary$lowest_unique_match
  seqtab_assigned <- seqtab_assigned[,order(classification_summary$lowest_unique_match)]
  
  species_collated <- collate_at_species_level(info_per_seq=classification_summary, seqtab_assigned=seqtab_assigned)
  
  species_collated <- species_collated[,order(colnames(species_collated))]
  write.csv(species_collated, file=classification_collated_by_species_path, quote=FALSE)
  
  presabs <- as.matrix((species_collated > 0) + 0)
  rowSums(presabs)
  
  pdf(classification_collated_by_species_plot_path, width=7, height=14)
  image(presabs, col=c("white", "red"), xaxt = "n", yaxt="n")
  axis(side = 1, labels = rownames(presabs),
       at = seq(0, by = 1/(nrow(presabs)-1), length.out = nrow(presabs)), las=3,
       cex.axis=0.15)
  axis(side = 2, labels=colnames(presabs),
       at=seq(0, by=1/(ncol(presabs)-1), length.out=ncol(presabs)), las=2,
       cex.axis=0.15)
  dev.off()
}
