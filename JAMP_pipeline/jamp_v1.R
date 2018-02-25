# 170201 jamp pipeline
setwd("~/Documents/UNI_und_VORLESUNGEN/GitHub/JAMP/")

#install.packages(c("bold", "XML", "seqinr"), dependencies=T)
#install.packages("JAMP", repos = NULL, type="source")


library("JAMP") #v0.34
library("bold")
library("XML")



# base directory
setwd("~/Documents/UNI_und_VORLESUNGEN/14 Guelph/1 TEACHING/2018 metabarcoding course/1 JAMP_empty")

Demultiplexing_shifted(file1="../16_S10_L001_R1_001_run1.fastq", file2="../16_S10_L001_R2_001_run1.fastq", tags="../_converter/indexe_1.csv", combinations="../_converter/combos_1.csv")

# Paired end merging
U_merge_PE(fastq_pctid=75)

# check for PhiX
system2("usearch", "-usearch_global A_Demultiplexing_shifted/_data/N_debris_r1.txt -db PhiX.fasta -id 0.9 -strand both -blast6out PhiX_table.txt")



# trimm primers
Cutadapt(forward="GGWACWGGWTGAACWGTWTAYCCYCC", # mlCOIintF
reverse="TANACYTCNGGRTGNCCRAARAAYCA") # jgHCO, I (inosin) replaced with N

# trim reads in different orrientation
PEmerged <- list.files("~/Documents/UNI_und_VORLESUNGEN/14 Guelph/1 TEACHING/2018 metabarcoding course/1 JAMP/B_U_merge_PE/_data", full.names=T)

Cutadapt(files= PEmerged, forward="TANACYTCNGGRTGNCCRAARAAYCA", # jgHCO, I (inosin) replaced with N
reverse="GGWACWGGWTGAACWGTWTAYCCYCC") # mlCOIintF


U_revcomp(RC=T)


# merge forward and now reverse complement reverse reads!

dir.create("F_merge/_data", recursive=T)

FW <- list.files("C_Cutadapt/_data", full.names=T)
RC <- list.files("E_U_revcomp/_data", full.names=T)

i <- 1

for (i in 1:length(FW)){
system2("cat", paste(FW[i], " ", RC[i], " > F_merge/_data/", sub("E_U_revcomp/_data/(.*_).*", "\\1", RC[i]), "merged.txt", sep=""))
}

cat(file="log.txt", append=T, c("\nPROCESSING MODULE:", "F_merge"), sep="\n")


# discard with non target length
Minmax(min=(313-10), max=(313+10))

# discard reads above 1 expected error
U_max_ee(max_ee=1)

# subsample to lowest sample size, should be done if samples are widely different in sequencing depth (as one starts with)
U_subset(sample_size=60000)


#cluster OTUs
U_cluster_otus(filter=0.01)
file.rename("J_U_cluster_otus", "J_U_cluster_otus - 60k")

#cluster OTUs (without subsetting)
no_subset <- list.files("~/Documents/UNI_und_VORLESUNGEN/14 Guelph/1 TEACHING/2018 metabarcoding course/1 JAMP/H_U_max_ee/_data", full.names=T)

U_cluster_otus(files= no_subset, filter=0.01)


# assign taxonomy to OTUs without sub setting! K_U_cluster_otus
Bold_web_hack(file="K_BOLD_TAX.txt")


# setting read counts below 0.01% to 0
data <- read.csv("K_U_cluster_otus/5_OTU_table_0.01.csv", stringsAsFactors=F)
colSum <- colSums(data[,3:9])
i <- 3
for(i in 3:9){
keep <- data[-nrow(data),i]/colSum[i-2]*100>=0.01 # find values blow 0.01
data[nrow(data),i] <- data[nrow(data),i]+sum(data[,i][!keep]) # below 0.01
data[,i][!keep] <- 0 # set lavues to 0
}
write.csv(data, "K_U_cluster_otus/data_cleaned.csv", row.names=F)

data2 <- data
#relative abundance
for(i in 3:9){
data2[,i] <- round(data2[,i]/colSum[i-2]*100, 4) #converte data into relative abundance
}
write.csv(data2, "K_U_cluster_otus/data_cleaned_rel.csv", row.names=F)




