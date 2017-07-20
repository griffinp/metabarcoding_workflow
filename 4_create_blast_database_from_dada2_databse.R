# Testing the DADA2 pipeline on A_E amplicons

#install.packages("RCurl", lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")

#source("https://bioconductor.org/biocLite.R")

#biocLite("RCurl", lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")
#biocLite("dada2",
#         lib.loc="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2",
#         lib="/vlsci/VR0267/pgriffin/R/x86_64-unknown-linux-gnu-library/3.3.2")

library(dada2)
packageVersion("dada2")
library(microclass)


if(file.exists("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/dada2_fixed_uppercase_nolong.fasta.gz")){
  system(command="gunzip /vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/dada2_fixed_uppercase_nolong.fasta.gz")
}
refdb <- readFasta("/vlsci/VR0267/pgriffin/hsm/mel_metabarcoding_run1/ref_db/dada2_fixed_uppercase_nolong.fasta")
message("Making ref_blastdb")
ref_blastdb <- blastDbase16S(name="dada2_fixed_uppercase_nolong_blastdb", sequence=refdb$Sequence, taxon=refdb$Header)
message("Finished making ref_blastdb")
