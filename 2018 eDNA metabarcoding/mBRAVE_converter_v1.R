#181108 convert fasta to fastq
setwd("~/Desktop/Tutorial_eDNA")


data <- readLines("~/Desktop/Tutorial_eDNA/J_U_cluster_otus/_data/5_subset/5_OTU_sub_0.01.fasta")

ID <- data[seq(1, length(data), 2)] # get sequence names
ID <- sub(">", "@", ID) # change > to @
seq <- data[seq(2, length(data), 2)] # get sequences
sequL <- nchar(seq) # get sequence length


for (i in 1:length(ID)){

temp <- paste(ID[i], seq[i], "+", paste(rep("K", sequL[i]), collapse=""), sep="\n")
cat(temp, file="fasta_exp.fastq", append=T, sep="\n")
}


i<- 1






