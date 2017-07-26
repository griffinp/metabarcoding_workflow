# R script to search for and remove primers from mel's metabarcoding round 2 sequencing

raw_dir <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/raw/AGRF_CAGRF13385_AWP8B/"
trimmed_dir <- "/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/trimmed/"
pc_trimmed_threshold <- 0.5
adapter_table <- read.csv("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run2/scripts/primer_sequences_for_R_script_run2.csv", sep=",", 
                    header=TRUE, stringsAsFactors=FALSE)
raw_files <- list.files(path=raw_dir, pattern="L001_R1.fastq$")

for(each_file in raw_files){
  
  basename <- strsplit(each_file, split="L001_R1.fastq")[[1]][1]
  # start with a raw R1 and R2 file
  
  input_R1<-paste(raw_dir, basename, "L001_R1.fastq", sep="")
  input_R2<-paste(raw_dir, basename, "L001_R2.fastq", sep="")
  trimmed_R1<-paste(trimmed_dir, basename, "_qtrimmed_R1.fastq", sep="")
  trimmed_R2<-paste(trimmed_dir, basename, "_qtrimmed_R2.fastq", sep="")
  
  message(paste("Input files are", input_R1, "and", input_R2))
  
  # basic quality-trim (STEP 1)
  step1_command <- paste("cutadapt -q 20 -a XXX -A XXX -o", trimmed_R1, "--paired-output", trimmed_R2,
                         input_R1, input_R2)
  message("Quality trimming, step 1")
  system(command=step1_command)

  # loop through possible adapter pairs
  
  
  for(i in 1:nrow(adapter_table)){
    
    temprow <- adapter_table[i,]
    temp_pair <- temprow$Pair_name
    tempP1_seq <- temprow$P1_seq
    tempP2_seq <- temprow$P2_seq
    tempP1_rc_seq <- temprow$P1_rc_seq
    tempP2_rc_seq <- temprow$P2_rc_seq
    # run cutadapt to remove anchored 5' P1f matches at start of R1 only (STEP 2)
    
    stepcode <- paste("_", temp_pair, "_step2", sep="")
    adapter_seq <- tempP1_seq
    input_R1<-paste(trimmed_dir, basename, "_qtrimmed_R1.fastq", sep="")
    input_R2<-paste(trimmed_dir, basename, "_qtrimmed_R2.fastq", sep="")
    trimmed_R1<-paste(trimmed_dir, basename, stepcode, "_trimmed_R1.fastq", sep="")
    trimmed_R2<-paste(trimmed_dir, basename, stepcode, "_trimmed_R2.fastq", sep="")
    untrimmed_R1<-paste(trimmed_dir, basename, stepcode, "_untrimmed_R1.fastq", sep="")
    untrimmed_R2<-paste(trimmed_dir, basename, stepcode, "_untrimmed_R2.fastq", sep="")
    
    step2_command <- paste("cutadapt -g ^", adapter_seq, " -e 0.1 --untrimmed-output ",
                           untrimmed_R1, " --untrimmed-paired-output ", untrimmed_R2, " -o ", 
                           trimmed_R1, " --paired-output ", trimmed_R2, " ", input_R1, " ", input_R2, sep="")
    message(paste("Adapter trimming, step 2 for", temp_pair))
    #message(paste("Command is:", step2_command))
    system(step2_command)
    
    number_raw_R1 <- as.numeric(system(command=paste("wc -l", input_R1, "| awk {'print $1'}"), intern = TRUE))/4
    number_step2_trimmed_R1 <- as.numeric(system(command=paste("wc -l", trimmed_R1, "| awk {'print $1'}"), intern = TRUE))/4
    #print(paste("str(number_raw_R1) is", str(number_raw_R1)))
    pc_trimmed <- (number_step2_trimmed_R1*100)/number_raw_R1
    
    message(paste(pc_trimmed, "% of reads were trimmed in step 2 for", temp_pair))
    
    # if less than threshold% of reads were trimmed,
    
    if(pc_trimmed < pc_trimmed_threshold){
      # delete trimmed and untrimmed output files
      system(command=paste("rm", trimmed_R1, trimmed_R2, untrimmed_R1, untrimmed_R2))

      input_R1 <- paste(trimmed_dir, basename, "_qtrimmed_R1.fastq", sep="")
      input_R2 <- paste(trimmed_dir, basename, "_qtrimmed_R2.fastq", sep="")
      message("Moving to next primer pair")
      
    } else{ # if more than threshold% of reads were trimmed,
      # run cutadapt on trimmed output files to remove anchored 3' P2 matches (should remove them at end of R2 only) (STEP 3)
      # for now, using a workaround to avoid a bug in cutadapt: swapping the P1 and P2 files.
      stepcode <- paste("_", temp_pair, "_step3", sep="")
      adapter_seq <- tempP2_seq
      step3_stepcode <- stepcode
      input_R1<-trimmed_R1
      input_R2<-trimmed_R2
      trimmed_R1<-paste(trimmed_dir, basename, stepcode, "_trimmed_R1.fastq", sep="")
      trimmed_R2<-paste(trimmed_dir, basename, stepcode, "_trimmed_R2.fastq", sep="")
      untrimmed_R1<-paste(trimmed_dir, basename, stepcode, "_untrimmed_R1.fastq", sep="")
      untrimmed_R2<-paste(trimmed_dir, basename, stepcode, "_untrimmed_R2.fastq", sep="")
      
      step3_command <- paste("cutadapt -g ^", adapter_seq, " -e 0.1 --untrimmed-output ",
                             untrimmed_R2, " --untrimmed-paired-output ", untrimmed_R1, " -o ", 
                             trimmed_R2, " --paired-output ", trimmed_R1, " ", input_R2, 
                             " ", input_R1, sep="")
      message(paste("Adapter trimming, step 3 for", temp_pair))
      #message(paste("Command is:", step3_command))
      system(step3_command)
      
      # run cutadapt on new trimmed output files to remove non-anchored P1f matches in forward orientation (STEP 4)
      
      stepcode <- paste("_", temp_pair, "_step4", sep="")
      adapter_seq <- tempP2_rc_seq
      input_R1<-trimmed_R1
      input_R2<-trimmed_R2
      trimmed_R1<-paste(trimmed_dir, basename, stepcode, "_trimmed_R1.fastq", sep="")
      trimmed_R2<-paste(trimmed_dir, basename, stepcode, "_trimmed_R2.fastq", sep="")
      # no longer producing 'untrimmed' output files in this step, since all reads are already
      # identified as belonging to this primer pair
      untrimmed_R1<-paste(trimmed_dir, basename, stepcode, "_untrimmed_R1.fastq", sep="")
      untrimmed_R2<-paste(trimmed_dir, basename, stepcode, "_untrimmed_R2.fastq", sep="")
      
      #step4_command <- paste("cutadapt -g ", adapter_seq, " -e 0.1 --untrimmed-output ",
      #                       untrimmed_R1, " --untrimmed-paired-output ", untrimmed_R2, " -o ", 
      #                       trimmed_R1, " --paired-output ", trimmed_R2, " ", input_R1, 
      #                       " ", input_R2, sep="")
      step4_command <- paste("cutadapt -a ", adapter_seq, " -e 0.05 --minimum-length 20 -o ", 
                             trimmed_R1, " --paired-output ", trimmed_R2, " ", input_R1, 
                             " ", input_R2, sep="")
      message(paste("Adapter trimming, step 4 for", temp_pair))
      #message(paste("Command is:", step4_command))
      system(step4_command)
      
      # run cutadapt on new trimmed output files to remove non-anchored P1r matches in reverse orientation
      # for now, using a workaround to avoid a bug in cutadapt: swapping the P1 and P2 files
      
      stepcode <- paste("_", temp_pair, "_step5", sep="")
      adapter_seq <- tempP1_rc_seq
      input_R1<-trimmed_R1
      input_R2<-trimmed_R2
      trimmed_R1<-paste(trimmed_dir, basename, stepcode, "_trimmed_R1.fastq", sep="")
      trimmed_R2<-paste(trimmed_dir, basename, stepcode, "_trimmed_R2.fastq", sep="")
      # no longer producing 'untrimmed' output files in this step, since all reads are already
      # identified as belonging to this primer pair
      untrimmed_R1<-paste(trimmed_dir, basename, stepcode, "_untrimmed_R1.fastq", sep="")
      untrimmed_R2<-paste(trimmed_dir, basename, stepcode, "_untrimmed_R2.fastq", sep="")
      
      #step5_command <- paste("cutadapt -G ", adapter_seq, " -e 0.1 --untrimmed-output ",
      #                       untrimmed_R1, " --untrimmed-paired-output ", untrimmed_R2, " -o ", 
      #                       trimmed_R1, " --paired-output ", trimmed_R2, " ", input_R1, 
      #                       " ", input_R2, sep="")
      step5_command <- paste("cutadapt -a ", adapter_seq, " -e 0.05 --minimum-length 20 -o ", 
                             trimmed_R2, " --paired-output ", trimmed_R1, " ", input_R2, 
                             " ", input_R1, sep="")
      message(paste("Adapter trimming, step 5 for", temp_pair))
      #message(paste("Command is:", step5_command))
      system(step5_command)
      
    }
  }
}
