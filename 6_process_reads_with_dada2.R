#install.packages("RCurl", lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")

#source("https://bioconductor.org/biocLite.R")

#biocLite("RCurl", lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")
#biocLite("dada2",
#         lib.loc="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2",
#         lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")


library(dada2)

input_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/trimmed"
results_path <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/results"
adapter_table <- read.csv("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/scripts/primer_sequences_for_R_script_run2.csv", sep=",", 
                          header=TRUE, stringsAsFactors=FALSE)
adapter_pairs <- adapter_table$Pair_name
# an example default value for sequence lengths to retain
#seq_lengths_to_keep <- data.frame(Pair_name=adapter_pairs, min_length=rep(100, times=nrow(adapter_table)), 
#                                  max_length=rep(600, times=nrow(adapter_table)))
seq_lengths_to_keep <- data.frame(Pair_name=adapter_pairs, min_length=c(220, 100, 244), 
                                  max_length=c(322, 600, 249))
### added code to adjust seq_lengths_to_keep for specific primer pairs,
### based on the expected sequence length and on the distribution of observed lengths
### after trimming - which can be seen in the files ***.seqlength_table.txt 

###############################
# Move cutadapt-trimmed files #
# (step 5 output only) into a #
# separate directory for easy #
# access                      #
###############################

trimmed_paths <- c()
for(i in adapter_pairs){
  trimmed_path <- paste(input_path, i, sep="/")
  dir.create(trimmed_path)
  system(command=paste("mv ", input_path, "/*", i, "_step5* ", trimmed_path, sep=""))
  trimmed_paths <- c(trimmed_paths, trimmed_path)
}

####################################
# tackle each primer pair in turn: #
# sort the files                   #
# and filter them with             #
# dada2::filterAndTrim             #
####################################

for(j in 1:length(adapter_pairs)){

  primer_pair <- adapter_pairs[j]
  trimmed_path <- trimmed_paths[j]
  
  message(paste("Processing primer pair", primer_pair))
  
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
  filtered_output_object_path <- paste(filt_path, "/", "out.rds", sep="")
  
  if(!file.exists(filtFs[1])){
    out <- filterAndTrim(fwd=fnFs, filt=filtFs, rev=fnRs, filt.rev=filtRs,
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose=TRUE)
    out
    saveRDS(out, filtered_output_object_path)
  } else{
    message("Some filtered files already exist; skipping filterAndTrim step")
    out <- readRDS(filtered_output_object_path)
  }
  
  #############################
  # Do error estimation,      #
  # or import previous        #
  # error estimates if they   #
  # were generated previously #
  #############################
  
  errF_path <- paste(results_path, "/", primer_pair, "_errF.rds", sep="")
  errR_path <- paste(results_path, "/", primer_pair, "_errR.rds", sep="")
  
  if(file.exists(errF_path) & file.exists(errR_path)){
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
  write.table(seqlength_table, quote=FALSE, row.names=FALSE,
              file=paste(results_path, "/", primer_pair, "_seqlength_table.txt", sep=""))
  
  ##############################################
  # Exclude very-short or very-long sequences; #
  # this requires manual examination of the    #
  # seqlength_table object to decide what are  #
  # sensible sequence lengths. Set these at    #
  # the top of this script...                  #
  ##############################################
  
  min_length <- seq_lengths_to_keep[which(seq_lengths_to_keep$Pair_name==primer_pair), "min_length"]
  max_length <- seq_lengths_to_keep[which(seq_lengths_to_keep$Pair_name==primer_pair), "max_length"]
  seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(min_length, max_length)]
  
  ###################################
  # Remove chimeric sequences       #
  # using dada2::removeBimeraDenovo #
  # command; or import instead if   #
  # already exists                  #
  ###################################
  
  seqtab.nochim_path <- paste(results_path, "/", primer_pair, "_seqtab_nochim.rds", sep="")
  
  if(file.exists(seqtab.nochim_path)){
    seqtab.nochim <- readRDS(seqtab.nochim_path)
  } else{
    seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
    print(paste(round(sum(seqtab.nochim)/sum(seqtab2), digits = 3), "of sequences remaining after chimera removal"))
    saveRDS(seqtab.nochim, seqtab.nochim_path)
  }
  
  ############################
  # Make a summary table to  #
  # display how many reads   #
  # were retained after each #
  # processing step          #
  ############################
  
  getN <- function(x){sum(getUniques(x))}
  track <- cbind(out, sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
  colnames(track) <- c("trimmed", "filtered", "merged", "tabled", "nonchim")
  rownames(track) <- sample.names
  write.table(track, file=paste(results_path, "/", primer_pair, "_read_number_summary.txt", sep=""),
              quote=FALSE)
}
