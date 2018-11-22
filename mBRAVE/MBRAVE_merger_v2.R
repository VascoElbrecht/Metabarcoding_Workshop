#mBRAVE merger
setwd("~/Documents/GitHub/Metabarcoding_Workshop/mBRAVE")

files <- list.files(path = "tables",full.names=T)


ID <- sub(".*/.*\\.(.*)\\..*\\..*", "\\1", files)

list2 <- read.csv(files[1], sep="\t", stringsAsFactors=F)
list <- list2[c(1,9)]

list2 <- list2[1:8]


for (g in 2:length(files)){
data <- read.csv(files[g], sep="\t", stringsAsFactors=F)
list2 <- rbind(list2, data[1:8])

list <- merge(list,data[c(1,9)],by.x="BIN.Taxon.ID",by.y="BIN.Taxon.ID",all=T)
list[is.na(list)] <- 0
}
names(list) <- c("BIN", ID)



list3 <- list2[!duplicated(list2[1]),]

list4 <- merge(list, list3, by.x="BIN", by.y="BIN.Taxon.ID")

write.csv(list4,file="merged.csv")







