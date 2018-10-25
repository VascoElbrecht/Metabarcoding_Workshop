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



SRA(ID=data$SRA.accession.number, rename=data$JAMP_Name)

FastQC()



# Paired end merging
U_merge_PE(fastq_pctid=75)

# check for PhiX
system2("usearch", "-usearch_global A_Demultiplexing_shifted/_data/N_debris_r1.txt -db PhiX.fasta -id 0.9 -strand both -blast6out PhiX_table.txt")



# trimm primers
Cutadapt(forward="GGWACWGGWTGAACWGTWTAYCCYCC", # mlCOIintF
reverse="TANACYTCNGGRTGNCCRAARAAYCA", bothsides=T) # jgHCO, I (inosin) replaced with N
#by using "bothsides=T", forward or reverse primers are detected on both ends. This is not nessesary for fusion primers.

# discard with non target length
Minmax(min=(313-10), max=(313+10))

# discard reads above 1 expected error
U_max_ee(max_ee=1)

# subsample to lowest sample size, should be done if samples are widely different in sequencing depth (as one starts with)
U_subset(sample_size=60000)




#cluster OTUs
U_cluster_otus(filter=0.01)
file.rename("G_U_cluster_otus", "G_U_cluster_otus - 60k")

#cluster OTUs (without subsetting)
no_subset <- list.files("E_U_max_ee/_data", full.names=T)

U_cluster_otus(files= no_subset, filter=0.01)

# assign taxonomy to OTUs without sub setting! K_U_cluster_otus
Bold_web_hack(file="K_bold_results.txt")





# haplotyping
# from merged data:
no_subset <- list.files("D_Minmax/_data", full.names=T)

# Keep only sequences of 313 bp length
Minmax(file=no_subset, min=313, max=313)

# Stricter EE filtering
U_max_ee(max_ee=0.2)

# Extracting haplotypes from the metabarcoding data
Denoise(minsize=5, minrelsize=0.001, OTUmin=0.01, minHaploPresence=1)





