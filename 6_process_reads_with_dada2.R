##### Custom function for classifying sequences using BLAST
##### adapted from microclass::blastClassify16S to return more than the top hit

blastClassify <- function (sequence, bdb) 
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

#path <- "~/Documents/metabarcoding/test_run1_data"
path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/trimmed/A_E"
filt_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/trimmed/A_E/filtered"

dir.create(filt_path)

fnFs <- sort(list.files(path, pattern="_R1.fastq"))
fnRs <- sort(list.files(path, pattern="_R2.fastq"))

if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

sample.names <- sapply(strsplit(fnFs, "_step"), `[`, 1)

fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#plotQualityProfile(fnFs[1])

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#out <- filterAndTrim(fwd=fnFs, filt=filtFs, rev=fnRs, filt.rev=filtRs,
#                     maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
#                     compress=TRUE, multithread=TRUE, verbose=TRUE)
#out

# errF <- learnErrors(filtFs, nread=2e6, multithread=TRUE)
# errR <- learnErrors(filtRs, nread=2e6, multithread=TRUE)
# 
# saveRDS(errF, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_errF.rds")
# saveRDS(errR, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_errR.rds")

#errF <- readRDS("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_errF.rds")
#errR <- readRDS("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_errR.rds")


# Unclear whether this approach assumes all files are for the 
# same amplicon (which may cause probs for Mel's design - need to split files 
# and run each sub-amplicon separately)

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

#dds <- vector("list", length(sample.names))
#names(dds) <- sample.names

# mergers <- vector("list", length(sample.names))
# names(mergers) <- sample.names
# print(sample.names)
# for(sam in sample.names) {
#   cat("Processing:", sam, "\n")
#   derepF <- derepFastq(filtFs[[sam]])
#   ddF <- dada(derepF, err=errF, multithread=TRUE)
#   derepR <- derepFastq(filtRs[[sam]])
#   ddR <- dada(derepR, err=errR, multithread=TRUE)
#   merger <- mergePairs(ddF, derepF, ddR, derepR)
#   mergers[[sam]] <- merger
# }
# rm(derepF); rm(derepR)
# 
# seqtab <- makeSequenceTable(mergers)
# table(nchar(getSequences(seqtab)))
# 
# saveRDS(mergers, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_mergers.rds")
# saveRDS(seqtab, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_seqtab.rds")

#mergers <- readRDS(file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_mergers.rds")
#mergers <- readRDS(file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_mergers.rds")
#seqtab <- readRDS(file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_seqtab.rds")

#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(310,403)]

# seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# sum(seqtab.nochim)/sum(seqtab2)
# 
# saveRDS(seqtab.nochim, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_seqtab_nochim.rds")
seqtab.nochim <- readRDS("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/A_E_seqtab_nochim.rds")

#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
#colnames(track) <- c("merged", "tabled", "nonchim")
#rownames(track) <- sample.names
#print(track)

# taxa1 <- assignTaxonomy(seqtab.nochim, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/subset1.fasta.gz", 
#                         multithread=TRUE, taxLevels = c("Kingdom", "Phylum", "Class",
#                                                         "Order", "Family", "Genus", "Species", "ID"))
# 
# write.csv(unname(taxa1), file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/taxa_table_A_E_subset1.csv")
# 
# taxa2 <- assignTaxonomy(seqtab.nochim, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/subset2.fasta.gz", 
#                         multithread=TRUE, taxLevels = c("Kingdom", "Phylum", "Class",
#                                                         "Order", "Family", "Genus", "Species", "ID"))
# 
# write.csv(unname(taxa2), file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/taxa_table_A_E_subset2.csv")
# 
# taxa3 <- assignTaxonomy(seqtab.nochim, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/subset3.fasta.gz", 
#                         multithread=TRUE, taxLevels = c("Kingdom", "Phylum", "Class",
#                                                         "Order", "Family", "Genus", "Species", "ID"))
# 
# write.csv(unname(taxa3), file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/taxa_table_A_E_subset3.csv")
# 
# taxa4 <- assignTaxonomy(seqtab.nochim, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/subset4.fasta.gz", 
#                         multithread=TRUE, taxLevels = c("Kingdom", "Phylum", "Class",
#                                                         "Order", "Family", "Genus", "Species", "ID"))
# 
# write.csv(unname(taxa4), file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/taxa_table_A_E_subset4.csv")
# 
# taxa5 <- assignTaxonomy(seqtab.nochim, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/subset5.fasta.gz", 
#                         multithread=TRUE, taxLevels = c("Kingdom", "Phylum", "Class",
#                                                         "Order", "Family", "Genus", "Species", "ID"))
# 
# write.csv(unname(taxa5), file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/taxa_table_A_E_subset5.csv")
# 
# taxa6 <- assignTaxonomy(seqtab.nochim, "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/subset6.fasta.gz", 
#                         multithread=TRUE, taxLevels = c("Kingdom", "Phylum", "Class",
#                                                         "Order", "Family", "Genus", "Species", "ID"))
# 
# write.csv(unname(taxa6), file="/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/results/taxa_table_A_E_subset6.csv")
