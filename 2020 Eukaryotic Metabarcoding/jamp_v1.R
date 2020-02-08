# 190204 jamp pipeline


setwd("~/Documents/University/2019 Bonn/physalia2020")


# Connect to SSH
sudo ssh -i ~/Documents/metaberlin.pem ubuntu@52.43.173.99


# change working directory
cd 

# list files 
ls -l

mkdir metabarcoding

chmod 777 metabarcoding

.libPaths("/usr/local/lib/R/site-library")

# Download the sequencing data
wget https://ndownloader.figshare.com/files/6503952 -O 16_S10_L001_R1_001_run1.fastq.gz
wget https://ndownloader.figshare.com/files/6503991 -O 16_S10_L001_R2_001_run1.fastq.gz

gunzip -k 16_S10_L001_R1_001_run1.fastq.gz
gunzip -k 16_S10_L001_R2_001_run1.fastq.gz



#get demultiplexing files (tags, tag combinations)
wget https://raw.githubusercontent.com/VascoElbrecht/JAMP/master/Tutorial/_converter/combos_1.csv
wget https://raw.githubusercontent.com/VascoElbrecht/JAMP/master/Tutorial/_converter/indexe_1.csv


# look at PhiX

wget https://raw.githubusercontent.com/VascoElbrecht/Metabarcoding_Workshop/master/2019%20Eukaryotic%20Metabarcoding/PhiX.fasta



# Upload PhiX.fasta - AWS -> FileZilla
system2("usearch", "-usearch_global 16_S10_L001_R1_001_run1.fastq -db PhiX.fasta -id 0.9 -strand both -blast6out PhiX_table.txt")

usearch -usearch_global 16_S10_L001_R1_001_run1.fastq -db PhiX.fasta -id 0.9 -strand both -blast6out PhiX_table.txt



# Start processing with R
sudo R


# base directory
setwd("~/Documents/GitHub/Metabarcoding_Workshop/2019 Eukaryotic Metabarcoding")


# load JAMP
library("JAMP")



#Create folders (to drop data into)
Empty_folder()

#Delete folders
Remove_last_folder()



# Demultiplex by sample
Demultiplexing_shifted(file1="16_S10_L001_R1_001_run1.fastq", file2="16_S10_L001_R2_001_run1.fastq", tags="indexe_1.csv", combinations="combos_1.csv")

file1="16_S10_L001_R1_001_run1.fastq"
file2="16_S10_L001_R2_001_run1.fastq"
tags="indexe_1.csv"
combinations="combos_1.csv"



# Analyse sequence quality
FastQC(exe="/usr/local/bin/FastQC/fastqc")


# Merge paired end reads
Merge_PE(LDist=T, exe="usearch")


# trimm primers
Cutadapt(forward="GGWACWGGWTGAACWGTWTAYCCYCC", reverse="TANACYTCNGGRTGNCCRAARAAYCA", bothsides=T, LDist=T) 


# discard with non target length
Minmax(min=(313-10), max=(313+10), LDist=T)

# discard reads above 1 expected error
U_max_ee(max_ee=0.5)



# subsample to lowest sample size, should be done if samples are widely different in sequencing depth (as one starts with)
U_subset(sample_size=50000)


#cluster OTUs
U_cluster_otus(filter=0.01)
file.rename("G_U_cluster_otus", "G_U_cluster_otus - 60k")

#cluster OTUs (without subsetting)
no_subset <- list.files("~/Documents/University/2019 Bonn/physalia2020/E_U_max_ee/_data", full.names=T)

U_cluster_otus(files= no_subset, filter=0.01)


# assign taxonomy to OTUs without sub setting! K_U_cluster_otus
Bold_web_hack(file="K_BOLD_TAX.txt")





# haplotyping

# temporary fix, old version used!
#source("https://raw.githubusercontent.com/VascoElbrecht/JAMP/9d4404d5a31c8635a58cf8ffd7628869831dae8f/JAMP/R/Denoise.R")



# from merged data:
no_subset <- list.files("D_Minmax/_data", full.names=T)

U_revcomp(files= no_subset, RC=T)
#Remove_last_folder()


U_truncate(trunclen=310)

#Keep only sequences of 313 bp length
#Minmax(min=310, max=310)

U_revcomp(RC=T)

# stricter EE filtering
U_max_ee(max_ee=0.2)


ee <- list.files("L_U_max_ee/_data", full.names=T)


# Denoising ESVs
Denoise(files=ee, minsize=5, minrelsize=0.001, OTUmin=0.01, minHaploPresence=2, poolsample=F, renameSamples="(.*)_PE_RC_RC_cut.*")

Remove_last_folder()


# remove all intermediate data (do only run when doen with all processing)
delete_data()

