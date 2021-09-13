#This file links the FASTA names from USEARCH
# to Otu names for use in other analysis

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

#importing USEARCH file
m.u <- read.table(file='mMethanotrophs98_20190224.txt',sep="\t",
                  dec =".")  

length(unique(m.u[,2])) # 22 OTUS

#creating matching dataframe
m.df <- cbind(c.otu,names(fa))

#mathcing USEARCH to names
mm.u <- m.df[c(which(names(fa) %in% m.u[,2])),]

#call phyloseq to identify methanotrophs identified by mothur
require(phyloseq)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file='CvMcrb03_20190225.taxonomy', 
                     parseFunction=parse_taxonomy_default)

#pulling in Methanotrophs from the methylocystaceae
Fcys.Gcys <- subset_taxa(TAX, Rank6=="Methylocystis") #Genus
Fcys.Gsin <- subset_taxa(TAX, Rank6=="Methylosinus") #Genus

#pulling in methylocella
Gcel <- subset_taxa(TAX, Rank6=="Methylocella")

#pulling in methanotrophs from the gammaproteobacteria
Ococ.Fcoc <- subset_taxa(TAX, Rank5=="Methylococcaceae") #Family
Ococ.Gcrn<- subset_taxa(TAX, Rank6=="Crenothrix") #Genus

#merging phyloseq
mtrp.m <- merge_phyloseq(Fcys.Gcys,Fcys.Gsin,Gcel,
                         Ococ.Fcoc,Ococ.Gcrn)

#creating vector of unique names
m.un <- unique(c(mm.u[,1],rownames(mtrp.m@.Data)))

#matching OTU names to Fasta
m.n <- m.df[c(which(c.otu %in% m.un)),]

#creating methanotroph FASTA
m.fa <- fa[c(which(names(fa) %in% m.n[,2]))]

names(m.fa) <- m.n[,1]

# 
# seqinr::write.fasta(sequences = m.fa, names = m.n[,1], nbchar = 6000,
#                     file.out ="OtuMethanotrophName.fasta")
