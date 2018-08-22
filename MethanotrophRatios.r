#This script analyses ratios of methylotrophs and methanotrophs

rm(list=ls())
getwd()
setwd('C:/Users/websterkgd/GitHub/CaveMethanotrophs')

#get the citation for the R package
citation(package = "base", lib.loc = NULL)

# loading the phyloseq package for microbial analysis
library("phyloseq")
packageVersion("phyloseq")

# loading the ggplot2 package
library("ggplot2")
packageVersion("ggplot2")
theme_set(theme_bw())

# calling vegan
#install.packages("vegan")
require('vegan')

#importing the shared file 
OTU <- import_mothur(mothur_shared_file = 'CvMcrb03.bac.final.shared', parseFunction = parse_taxonomy_default)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file = 'CvMcrb03.bac.final.0.03.taxonomy', parseFunction = parse_taxonomy_default)

#create a physeq table for analysis
physeq = phyloseq(OTU, TAX)

#Creating transformations of the data 
f.physeq = transform_sample_counts(physeq, function(x) x / sum(x)) # fractional abundance transformation of the data

# import meta data
md.d = read.table('sample_meta_data+d.txt', header = TRUE, sep = "", dec = ".", row.names = 1) 
md = read.table('sample_meta_data-d.txt', header = TRUE, sep = "", dec = ".", row.names = 1) 
  
#Let Phyloseq see the meta data
sample_meta_data <- sample_data(md, errorIfNULL = TRUE)

#Pruning NAs in the "design matrix" for analysis of distance on microbial communities
sample_meta_data.d <- sample_data(md.d, errorIfNULL = TRUE)

sample_meta_data$Location <- as.factor(sample_meta_data$Location)
sample_meta_data$Bin.Loc <- as.factor(sample_meta_data$Bin.Loc)

sample_meta_data.d$Location <- as.factor(sample_meta_data.d$Location)
sample_meta_data.d$Bin.Loc <- as.factor(sample_meta_data.d$Bin.Loc)

#remove na from meta data
p.sample_meta_data.d <- na.omit(sample_meta_data.d) 

#merge meta data with other data
AllCaves <- merge_phyloseq(f.physeq, sample_meta_data)

#Rank1 =K, Rank2 =P, Rank3 =C, Rank4 =O, Rank5 =F, Rank6 =G, Rank7 =S
#Select the taxa of interest
#Nazaries et al., 2013 helped me look for methanotrophs. 

Fcys <- subset_taxa(AllCaves, Rank5=="Methylocystaceae") #Family Fcys
Fcys.Gcys <- subset_taxa(Fcys, Rank6=="Methylocystis") #Genus
Fcys.Gsin <- subset_taxa(Fcys, Rank6=="Methylosinus") #Genus
Fcys.Guct <- subset_taxa(Fcys, Rank6=="uncultured") #Genus 
Fcys.Gucs <- subset_taxa(Fcys, Rank6=="unclassified") #Genus

#Proteobacteria Betaproteobacteria Methylophilales
Fphl <- subset_taxa(AllCaves, Rank5=="Methylophilaceae") # = Fphl, Family
Fphl.Gphl <-subset_taxa(Fphl, Rank6=="Methylophilus") #Genus
Fphl.Gbac <-subset_taxa(Fphl, Rank6=="Methylobacillus") #Genus
Fphl.GD28 <-subset_taxa(Fphl, Rank6=="LD28") #Genus
Fphl.Guct <-subset_taxa(Fphl, Rank6=="uncultured") #Genus
Fphl.Gucs <-subset_taxa(Fphl, Rank6=="unclassified") #Genus

#Proteobacteria Gammaproteobacteria
Ococ <- subset_taxa(AllCaves, Rank4=="Methylococcales") #Order Ococ
Ococ.Fcoc <- subset_taxa(Ococ, Rank5=="Methylococcaceae") #Family
Ococ.Fucs <- subset_taxa(Ococ, Rank5=="unclassified") #Family
Ococ.Fcoc.Gbac <- subset_taxa(Ococ.Fcoc, Rank6=="Methylobacter") #Genus
Ococ.Fcoc.Gsar <- subset_taxa(Ococ.Fcoc, Rank6=="Methylosarcina") #Genus
Ococ.Fcoc.Gcoc <- subset_taxa(Ococ.Fcoc, Rank6=="Methylococcus") #Genus
Ococ.Fcoc.Gcal <- subset_taxa(Ococ.Fcoc, Rank6=="Methylocaldum") #Genus
Ococ.Fcoc.Gmon <- subset_taxa(Ococ.Fcoc, Rank6=="Methylomonas") #Genus
Ococ.Fcoc.Gmic <- subset_taxa(Ococ.Fcoc, Rank6=="Methylomicrobium") #Genus
Ococ.Fcoc.Gsom <- subset_taxa(Ococ.Fcoc, Rank6=="Methylosoma") #Genus
Ococ.Fcoc.Gucs <- subset_taxa(Ococ.Fcoc, Rank6=="unclassified") #Genus
# family = Crenotrichaceae
Ococ.C.Gcrn<- subset_taxa(Ococ, Rank6=="Crenothrix") #Genus

#Proteobacteria Betaproteobacteria Burkholderiales Comamdaceae 
Glib <- subset_taxa(AllCaves, Rank6=="Methylibium") #Genus Glib

#Proteobacteria Alphaproteobacteria Rhizobiales
Fbac <- subset_taxa(AllCaves, Rank5=="Methylobacteriaceae") #Family Fbac
Fbac.Gbac <- subset_taxa(Fbac, Rank6=="Methylobacterium") #Genus
Fbac.Guct <- subset_taxa(Fbac, Rank6=="uncultured") #Genus
Fbac.Gucs <- subset_taxa(Fbac, Rank6=="unclassified") #Genus

#Verrucomicrobia
V.Caci <- subset_taxa(AllCaves, Rank3=="Acidimethylosilex") #V.Caci Class

#Proteobacteria Alphaproteobacteria Rhizobiales  Beijerinckiaceae
Gcel <- subset_taxa(AllCaves, Rank6=="Methylocella") # Genus Gcel

#Proteobacteria Alphaproteobacteria Rhizobiales
Fbei <- subset_taxa(AllCaves, Rank5=="Beijerinckiaceae") #Family

Fbei.Gucs <- subset_taxa(Fbei, Rank6=="unclassified")
Fbei.Guct <- subset_taxa(Fbei, Rank6=="uncultured")

#total methanotrophs
Tman <- merge_phyloseq(Fcys, Ococ, Gcel,V.Caci,Gcel) # pulling the total methanotrophic community

#total methylotrophs
Tlot <- merge_phyloseq(Fphl, Glib, Fbac) 

#ratios of methanotrophs to whole community

MtW <- sample_sums(Tman)/sample_sums(f.physeq)
MtW <- na.omit(MtW)
is.na(MtW) <- sapply(MtW, is.infinite)
MtW <- na.omit(MtW)
MtW[MtW == 0] <- NA
MtW <- na.omit(MtW)
range(MtW) #1.1E-05 2.2E-02
median(MtW) # 0.0023
quantile(MtW, seq(0, 1, by=0.05)) #Q1 = 5.4E-04, Q3 = 4.7E-03
boxplot(MtW)
title(ylab = "Fractional Methanotroph Abundance")

#Testing MtW for normality
plot(density(MtW))
shapiro.test(MtW) #0.0019
qqnorm(MtW)
qqline(MtW, col = 2) 

