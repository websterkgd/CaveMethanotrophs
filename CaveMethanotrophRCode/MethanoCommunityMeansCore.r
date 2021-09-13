#This r code explores basic descriptors of the cave methanotroph community

rm(list=ls())

#get the citation for the R package
citation(package = "base", lib.loc = NULL)

# loading the phyloseq package for microbial analysis
library("phyloseq")
packageVersion("phyloseq")

# loading the ggplot2 package
library("ggplot2")
packageVersion("ggplot2")
theme_set(theme_bw())

#calling vegan
require('vegan')

#importing the shared file 
OTU <- import_mothur(mothur_shared_file = 'CvMcrb03_20190225.shared', parseFunction = parse_taxonomy_default)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file = 'CvMcrb03_20190225.taxonomy', parseFunction = parse_taxonomy_default)

#create a physeq table for analysis
physeq = phyloseq(OTU, TAX)

#Creating transformations of the data 
f.physeq = transform_sample_counts(physeq, function(x) x / sum(x)) # fractional abundance transformation of the data

# import meta data
md = read.table('sample_meta_data+d.txt', header = TRUE, sep = "", dec = ".", row.names = 1) 
  
#Let Phyloseq see the meta data
sample_meta_data <- sample_data(md, errorIfNULL = TRUE)

sample_meta_data$Location <- as.factor(sample_meta_data$Location)
sample_meta_data$Bin.Loc <- as.factor(sample_meta_data$Bin.Loc)

#merge meta data with other data
AllCaves <- merge_phyloseq(f.physeq, sample_meta_data)

#Importing Identified Methanotrophs
Id.m <- read.table('NCBIPutMntrphIDCave_20190228.txt',
                   header = T, sep ="\t")

Id.m <- Id.m[-c(which(is.na(Id.m[,2]) == T)),]

#creating methanotroph OTU table
m.otu <-AllCaves@otu_table@.Data[
  c(which(rownames(AllCaves@otu_table@.Data) %in%
                Id.m[,1])),]

m.otu <- otu_table(m.otu,taxa_are_rows = T)

#creating methanotroph Tax table
m.tax <-AllCaves@tax_table@.Data[
  c(which(rownames(AllCaves@tax_table@.Data) %in%
            Id.m[,1])),]

m.tax[,7] <- paste(as.factor(Id.m[1:22,2]))

m.tax <- tax_table(m.tax)

C.Mtf <- merge_phyloseq(m.otu,m.tax,sample_meta_data)

#pulling in specific methanotrophs
gUSc <- subset_taxa(C.Mtf, Rank7=='USCg')
gUSa <- subset_taxa(C.Mtf, Rank7=='USCa')
gsin <- subset_taxa(C.Mtf, Rank7=='Methylosinus')
gbac <- subset_taxa(C.Mtf, Rank7=='Methylobacter')
gcoc <- subset_taxa(C.Mtf, Rank7=='Methylococcus')
gcys <- subset_taxa(C.Mtf, Rank7=='Methylocystis')
gmon <- subset_taxa(C.Mtf, Rank7=='Methylomonas')
gcel <- subset_taxa(C.Mtf, Rank7=='Methylocella')
gglo <- subset_taxa(C.Mtf, Rank7=='Methyloglobulus')
gsom <- subset_taxa(C.Mtf, Rank7=='Methylosoma')
gvol <- subset_taxa(C.Mtf, Rank7=='Methylovolum')
gcap <- subset_taxa(C.Mtf, Rank7=='Methylocapsa')
gmic <- subset_taxa(C.Mtf, Rank7=='Methylomicrobium')

#subsetting by type
TUSC <- merge_phyloseq(gUSc,gUSa)
T.II <- merge_phyloseq(gcel,gcys,gsin,gcap)
Tp.I <- merge_phyloseq(gbac,gcoc,gmon,gglo,
                       gsom,gvol,gmic)

#sample_sums() shows how many individuals were observed in the data set
#sample_sums(data_set)
sample_sums(C.Mtf) 
sample_sums(TUSC)
sample_sums(T.II)
sample_sums(Tp.I)

((length(sample_sums(C.Mtf)))-
         length(which(sample_sums(C.Mtf)==0)))/
    length(sample_sums(C.Mtf)) #97.6

#range methanotrophs
range(sample_sums(C.Mtf)) #0 0.06
median(sample_sums(C.Mtf)) #0.0052

#median proportion of different clade
median(sample_sums(TUSC)/sample_sums(C.Mtf),na.rm = T) #0.993
median(sample_sums(gUSc)/sample_sums(C.Mtf),na.rm = T) #0.7008
median(sample_sums(gUSa)/sample_sums(C.Mtf),na.rm = T) #0.253
median(sample_sums(T.II)/sample_sums(C.Mtf),na.rm = T) #0.0012
median(sample_sums(Tp.I)/sample_sums(C.Mtf),na.rm = T) #0.0

###Core Methanotroph Microbiome 
#transform to sorenson
mD <- as.data.frame(m.otu@.Data)
for (i in 1:length(mD[,1])){
  for (j in 1:length(mD[1,])){
    if ((mD[i,j] > 0)==T){
      mD[i,j] <- 1
      }
  }
}

C.Mtf.core.22 = rownames(mD[which(rowSums(mD) >= 26),]) # 2 OTUs

length(which(colSums(mD[1:2,]) == 2)) #36 
                      

###How many methanotroph species are present at each location?

smD <- colSums(mD)

m_smD <- mean(smD) #3

hist(smD,12) #plots rough histogram of distribution

sd_smD <- sd(smD) #2