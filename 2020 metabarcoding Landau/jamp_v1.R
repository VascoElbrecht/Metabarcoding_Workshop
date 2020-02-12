# 20200209 jamp pipeline
# study - https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12789



setwd("~/Desktop/meta_landau")


# Connect to SSH
sudo ssh -i ~/Documents/metaberlin.pem ubuntu@52.43.173.99


# change working directory
cd 

# list files 
ls -l

mkdir metabarcoding

chmod 777 metabarcoding

cd metabarcoding

#.libPaths("/usr/local/lib/R/site-library")


# Start processing with R
sudo R # running things as sudo can be dangerous! 



# load JAMP
library("JAMP")



#Create folders (to drop data into)
Empty_folder()

#Delete folders
Remove_last_folder()

# Download data from NCBI
SRA(c("SRX5996960", "SRX5996959", "SRX5996934", "SRX5996933", "SRX5996864", "SRX5996863"), c("T01_A", "T01_B", "K01_A", "K01_B", "Sa01_A", "Sa01_B"))


# If download too slow, cancel with [cmd]+[c]


Empty_folder()
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226345/T01_71_BF23_BR22_r1.fastq.1 -O A_Empty_Folder/_data/T01_A_R1.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226345/T01_71_BF23_BR22_r2.fastq.1 -O A_Empty_Folder/_data/T01_A_R2.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226374/T01_101_BR20_BF20_r1.fastq.1 -O A_Empty_Folder/_data/T01_B_R1.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226374/T01_101_BR20_BF20_r2.fastq.1 -O A_Empty_Folder/_data/T01_B_R2.fastq")

system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226399/K01_116_BR2B_BF24_r1.fastq.1 -O A_Empty_Folder/_data/K01_A_R1.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226399/K01_116_BR2B_BF24_r2.fastq.1 -O A_Empty_Folder/_data/K01_A_R2.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226400/K01_86_BF22_BR21_r1.fastq.1 -O A_Empty_Folder/_data/K01_B_R1.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226400/K01_86_BF22_BR21_r2.fastq.1 -O A_Empty_Folder/_data/K01_B_R2.fastq")

system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226469/Sa01_61_BF22_BR20_r1.fastq.1 -O A_Empty_Folder/_data/Sa01_A_R1.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226469/Sa01_61_BF22_BR20_r2.fastq.1 -O A_Empty_Folder/_data/Sa01_A_R2.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226470/Sa01_91_BR21_BF21_r1.fastq.1 -O A_Empty_Folder/_data/Sa01_B_R1.fastq")
system2("wget", "https://sra-pub-src-1.s3.amazonaws.com/SRR9226470/Sa01_91_BR21_BF21_r2.fastq.1 -O A_Empty_Folder/_data/Sa01_B_R2.fastq")






# Analyse sequence quality
FastQC(exe="/usr/local/bin/fastqc/")



# Merge paired end reads
# also works with Vsearch (to avoid 32 bit limitation of usearch)
Merge_PE(LDist=T, exe="usearch")



cbind(list.files("B_merge_PE/_data"), c(T, F, F, T, F, T))

# reverse complement in reverse orientatioon
U_revcomp(RC=c(T, F, F, T, F, T))


# trimm primers
Cutadapt(forward="GCHCCHGAYATRGCHTTYCC", reverse="TCDGGRTGNCCRAARAAYCA", bothsides=F, LDist=T) 


# discard with non target length
Minmax(min=(421-10), max=(421+10), LDist=T)

# discard reads above 1 expected error
U_max_ee(max_ee=1)



# subsample to lowest sample size, should be done if samples are widely different in sequencing depth (as one starts with)
#U_subset(sample_size=50000)


#cluster OTUs
U_cluster_otus(filter=0.01)


# assign taxonomy to OTUs without sub setting! K_U_cluster_otus
Bold_web_hack(file="K_BOLD_TAX.txt")





# haplotyping

# temporary fix, old version used!
#source("https://raw.githubusercontent.com/VascoElbrecht/JAMP/9d4404d5a31c8635a58cf8ffd7628869831dae8f/JAMP/R/Denoise.R")



# from merged data:
no_subset <- list.files("E_Minmax/_data", full.names=T)

# Exact length
Minmax(file=no_subset, min=421, max=421)


# stricter EE filtering
U_max_ee(max_ee=0.5)


# Denoising ESVs
Denoise(minsize=5, minrelsize=0.001, OTUmin=0.01, minHaploPresence=2, poolsample=F, renameSamples="(.*)_PE_RC_cut.*")

#Remove_last_folder()


# remove all intermediate data (do only run when doen with all processing)
delete_data()

