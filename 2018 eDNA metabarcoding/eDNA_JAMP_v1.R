# 181025 jamp eDNA
# @VascoElbrecht

# See wiki for full installation guide!
# https://github.com/VascoElbrecht/Metabarcoding_Workshop/wiki/1)-Software-to-install

# Installing dependencies needed fro JAMP
install.packages(c("bold", "XML", "seqinr", "devtools", "fastqcr"), dependencies=T)
# Load devtools and install package directly from GitHub
library("devtools")
install_github("VascoElbrecht/PrimerMiner", subdir="PrimerMiner")
install_github("VascoElbrecht/JAMP", subdir="JAMP")




setwd("~/Desktop/Tutorial_eDNA")
list.files()

library("JAMP")

# JAMP generates a new folder for each processing step! Should your files alreay be demultiplexed or you want to start somwhere else in the pipeline with preprocessed reads, you can generate an empty folder and place your files to be processed in "_data"
Empty_folder()

#To delete the last generated folder, run
Remove_last_folder() # this is yout "undo" function. Do not delete folders "randomly" as this breaks the log.txt


# let's read the SRA IDs!
data <- read.csv("SRA_list_v1.csv", stringsAsFactors=F)

head(data) #you can also use str(data)

# Downloading sequences:
SRA(ID=data$SRA.accession.number, rename=data$JAMP_Name)

# Quick read quality check!
FastQC()



# Paired end merging
U_merge_PE(fastq_pctid=75)

data <- read.csv("B_U_merge_PE/_stats/B_sequ_length_abund.csv", stringsAsFactors=F)

data

# check for PhiX
#system2("usearch", "-usearch_global A_Demultiplexing_shifted/_data/N_debris_r1.txt -db PhiX.fasta -id 0.9 -strand both -blast6out PhiX_table.txt")


# mlCOIintF+jgHCO2198
# trimm primers
Cutadapt(forward="GGWACWGGWTGAACWGTWTAYCCYCC", # mlCOIintF
reverse="TAIACYTCIGGRTGICCRAARAAYCA", LDist=T) # jgHCO2198


# discard with non target length
Minmax(min=(313-10), max=(313+10))

# discard reads above 1 expected error
U_max_ee(max_ee=1)

# subsample to lowest sample size, should be done if samples are widely different in sequencing depth, but usually not needed

U_subset(sample_size=25000) # will produce some errors because files contain less than 25000 sequences (would need to move negative controls over manually)
Remove_last_folder()


#cluster OTUs
U_cluster_otus(filter=0.01)

# assign taxonomy to OTUs without sub setting! K_U_cluster_otus
Bold_web_hack(file="K_bold_results.txt")




# second primer set
Next_primer_set <- list.files("B_U_merge_PE/_data", full.names=T)

Cutadapt(Next_primer_set, forward="GGTCAACAAATCATAAAGATATTGG", #LCO1490
reverse="GGIGGRTAIACIGTTCAICC", LDist=T) # Ill_C_R

#renamefiles!
names <- list.files("G_Cutadapt/_data", full.name=T)
names2 <- sub("_data/", "-data/ill_", names)

file.rename(names, names2)

Minmax(min=(325-10), max=(325+10))

U_max_ee(max_ee=1)

U_cluster_otus(filter=0.01)








# haplotyping
# from merged data:
no_subset <- list.files("D_Minmax/_data", full.names=T)

# Keep only sequences of 313 bp length
Minmax(file=no_subset, min=313, max=313)

# Stricter EE filtering
U_max_ee(max_ee=0.2)

# Extracting haplotypes from the metabarcoding data
Denoise(minsize=5, minrelsize=0.001, OTUmin=0.01, minHaploPresence=1)





