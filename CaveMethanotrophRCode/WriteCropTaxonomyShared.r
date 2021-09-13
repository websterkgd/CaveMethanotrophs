#This file crops the Mothur Taxonomy file based on 
#the FASTA file accepted by NCBI 

rm(list=ls())

#get the citation for the R package
citation(package = "base", lib.loc = NULL)

#import packages
require("seqinr")

#importing base FASTA file
fo <- seqinr::read.fasta("CvMcrb201605.final.0.03.pick.fasta")

#importing accepted FASTA File  
fa <- seqinr::read.fasta("CvMcrb20190213.final.0.03.pick.filt.fasta")

#creating Otu labels

L.na <- rep('Otu', length(fo))

for (i in 1:length(L.na)){
  L.na[i] <- paste('Otu',
                   formatC(i,width=6,format="d", 
                           flag="0"),
                   sep = "")
}

#creating crop vector
#first removing all characters after the tab ##I could have used apply here
n.t <- names(fo)
for (i in 1:length(n.t)){
  n.t[i] <- gsub("\t.*","",n.t[i])
}

cr <- which(n.t %in% names(fa))

c.otu <- L.na[c(cr)]

#importing Taxonomy
require(phyloseq)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file = 'CvMcrb03.bac.final.0.03.taxonomy', parseFunction = parse_taxonomy_default)

c.tax <- TAX@.Data[which(row.names(TAX@.Data)
                         %in% c.otu),]  

write.table(c.tax, file = "CvMcrb03_20190225.taxonomy",row.names=FALSE,
            sep=";",quote=FALSE) #temp taxonomy file

# ###Writing Final NCBI Taxonomy file
# npc.Tax <- read.table(file='CvMcrb03_20190225.taxonomy',sep ="",
#                         header=FALSE, row.names = NULL)
# 
# empty <- as.data.frame(rep(NA, length(c.tax[,1])))
# names <- as.data.frame(rownames(c.tax))
# n.Tax<-cbind(names,empty,npc.Tax[2:(length(npc.Tax[,1])),1])
# colnames(n.Tax) <-c('OTU','Size','Taxonomy')
#  
# # #Writing to table
# write.table(n.Tax, file = "CvMcrb03_20190225.taxonomy",row.names=FALSE,
#             sep="\t",quote=FALSE)

#importing the shared file
OTU <- import_mothur(mothur_shared_file = 'CvMcrb03.bac.final.shared', parseFunction = parse_taxonomy_default)

#transposing OTU_table
t.OTU <- as.data.frame(t(OTU@.Data))
#cropping the OTU table
#OTUs in columns
c.t.OTU <- t.OTU[,which(colnames(t.OTU) %in% c.otu)]
#Add similarity label
label <- rep(0.03, length(c.t.OTU[,1]))
#Add OTU sums
numOtus <-rowSums(c.t.OTU)

#create rownames vector rownames
Group <- row.names(c.t.OTU)
#putting data into .shared format
s.c.t.OTU <- cbind(label,Group,numOtus,c.t.OTU)

# #Writing OTU Table as shared file
# write.table(s.c.t.OTU, file = "CvMcrb03_20190225.shared",
#             row.names=F,sep="\t",quote=FALSE)
