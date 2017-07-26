library(dada2)

input_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/trimmed"
results_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/results"
adapter_table <- read.csv("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/scripts/primer_sequences_for_R_script_run2.csv", sep=",", 
                          header=TRUE, stringsAsFactors=FALSE)
adapter_pairs <- adapter_table$Pair_name
seq_lengths_to_keep <- data.frame(Pair_name=adapter_pairs, min_length=rep(100, times=length(adapter_pairs)), 
                                  max_length=rep(600, times=length(adapter_pairs)))
### add code to adjust seq_lengths_to_keep for specific primer pairs here if needed ###

###############################
# Move cutadapt-trimmed files #
# (step 5 output only) into a #
# separate directory for easy #
# access                      #
###############################

trimmed_paths <- c()
for(i in adapter_pairs){
  trimmed_path <- paste(input_path, i, sep="/")
  dir.create(filt_path)
  system(command=paste("mv ", input_path, "/*", i, "_step5* ", trimmed_path, sep=""))
  trimmed_paths <- c(trimmed_paths, trimmed_path)
}

#################################
# tackle the first primer pair: #
# sort the files                #
# and filter them with          #
# dada2::filterAndTrim          #
#################################

primer_pair <- adapter_pairs[1]
trimmed_path <- trimmed_paths[1]

fnFs <- sort(list.files(trimmed_path, pattern="_R1.fastq"))
fnRs <- sort(list.files(trimmed_path, pattern="_R2.fastq"))

if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

sample.names <- sapply(strsplit(fnFs, "_step"), `[`, 1)

fnFs <- file.path(trimmed_path, fnFs)
fnRs <- file.path(trimmed_path, fnRs)

#plotQualityProfile(fnFs[1])

filt_path <- file.path(trimmed_path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fwd=fnFs, filt=filtFs, rev=fnRs, filt.rev=filtRs,
                     maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose=TRUE)
out

#############################
# Do error estimation,      #
# or import previous        #
# error estimates if they   #
# were generated previously #
#############################

errF_path <- paste(results_path, "/", primer_pair, "_errF.rds", sep="")
errR_path <- paste(results_path, "/", primer_pair, "_errR.rds", sep="")

if(exists(errF_path) & exists(errR_path)){
  errF <- readRDS(errF_path)
  errR <- readRDS(errR_path)
} else{
  errF <- learnErrors(filtFs, nread=2e6, multithread=TRUE)
  errR <- learnErrors(filtRs, nread=2e6, multithread=TRUE)
  saveRDS(errF, errF_path)
  saveRDS(errR, errR_path)
}

# Output error plots
pdf(file=paste(results_path, "/", primer_pair, "_error_plots.pdf", sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

dds <- vector("list", length(sample.names))
names(dds) <- sample.names

####################################
# Merge and de-replicate sequences #
# or import previously-made files  #
# if they already exist            #
####################################

merger_path <- paste(results_path, "/", primer_pair, "_mergers.rds", sep="")
seqtab_path <- paste(results_path, "/", primer_pair, "_seqtab.rds", sep="")

if(file.exists(merger_path) & file.exists(seqtab_path)){
  mergers <- readRDS(merger_path)
  seqtab <- readRDS(seqtab_path)
} else{
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  print(sample.names)
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
  }
  rm(derepF); rm(derepR)
  seqtab <- makeSequenceTable(mergers)
  saveRDS(mergers, merger_path)
  saveRDS(seqtab, seqtab_path)
}

seqlength_table <- table(nchar(getSequences(seqtab)))
write.table(seqlength_table, file=paste(results_path, "/", primer_pair, "seqlength_table.txt", sep=""))

##############################################
# Exclude very-short or very-long sequences; #
# this requires manual examination of the    #
# seqlength_table object to decide what are  #
# sensible sequence lengths. Set these at    #
# the top of this script...                  #
##############################################

min_length <- seq_lengths_to_keep[Pair_name==primer_pair, "min_length"]
max_length <- seq_lengths_to_keep[Pair_name==primer_pair, "max_length"]
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(min_length, max_length)]

###################################
# Remove chimeric sequences       #
# using dada2::removeBimeraDenovo #
# command; or import instead if   #
# already exists                  #
###################################

seqtab.nochim_path <- paste(results_path, "/", primer_pair, "_seqtab_nochim.rds")

if(file.exists(seqtab.nochim_path)){
  seqtab.nochim <- readRDS(seqtab.nochim_path)
} else{
  seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
  print(paste("After chimera removal, sequence table has dimensions", dim(seqtab.nochim)))
  print(paste(sum(seqtab.nochim)/sum(seqtab2), "sequences remaining after chimera removal"))
  saveRDS(seqtab.nochim, seqtab.nochim_path)
}

############################
# Make a summary table to  #
# display how many reads   #
# were retained after each #
# processing step          #
############################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("merged", "tabled", "nonchim")
rownames(track) <- sample.names
write.table(track, file=paste(results_path, "/", primer_pair, "_read_number_summary.txt", sep=""))

