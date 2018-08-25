#This r code plots the abundance of the methylocystaceae and methylococcales
#against methane concentration

rm(list=ls())
getwd()
setwd("C:/Users/websterkgd/GitHub/CaveMethanotrophs")

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
md = read.table('sample_meta_data+d.txt', header = TRUE, sep = "", dec = ".", row.names = 1) 
  
#Let Phyloseq see the meta data
sample_meta_data <- sample_data(md, errorIfNULL = TRUE)

sample_meta_data$Location <- as.factor(sample_meta_data$Location)
sample_meta_data$Bin.Loc <- as.factor(sample_meta_data$Bin.Loc)

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

#Verrucomicrobia
V.Caci <- subset_taxa(AllCaves, Rank3=="Acidimethylosilex") #V.Caci Class

#Proteobacteria Alphaproteobacteria Rhizobiales  Beijerinckiaceae
Gcel <- subset_taxa(AllCaves, Rank6=="Methylocella") # Genus Gcel

#Total methanotrophic community
Tman <- merge_phyloseq(Fcys, Ococ, Gcel) # pulling the total methanotrophic community

## Metahnotroph subgroups vs grainsize proportion
# Ococ vs Gravel.
OcocVGrav <- cor.test(x=sample_meta_data$X.Gravel,
                      y=log10(sample_sums(Ococ)), method = "spearman")
OcocVGrav #rho = 0.06 p =0.7

#Fcys vs Gravel
FcysVGrav <- cor.test(x=sample_meta_data$X.Gravel,
                      y=log10(sample_sums(Fcys)), method = "spearman")
FcysVGrav #rho = 0.39 p =0.01 s = 10231

#Ococ vs Sand
OcocVSand <- cor.test(x=sample_meta_data$X.Sand,
                      y=log10(sample_sums(Ococ)), method = "spearman")
OcocVSand #rho = 0.04 p =0.81

#Fcys vs Sand
FcysVSand <- cor.test(x=sample_meta_data$X.Sand,
                      y=log10(sample_sums(Fcys)), method = "spearman")
FcysVSand #rho = 0.12 p =0.12

#Ococ vs Silt
OcocVSilt <- cor.test(x=sample_meta_data$X.Silt,
                      y=log10(sample_sums(Ococ)), method = "spearman")
OcocVSilt #rho = -0.05 p =0.76

#Fcys vs Silt
FcysVSilt <- cor.test(x=sample_meta_data$X.Silt,
                      y=log10(sample_sums(Fcys)), method = "spearman")
FcysVSilt #rho = -0.14 p =0.36

#Ococ vs Clay
OcocVClay <- cor.test(x=sample_meta_data$X.Clay,
                      y=log10(sample_sums(Ococ)), method = "spearman")
OcocVClay #rho = -0.02 p =0.90

#Fcys vs Clay
FcysVClay <- cor.test(x=sample_meta_data$X.Clay,
                      y=log10(sample_sums(Fcys)), method = "spearman")
FcysVClay #rho = 0.08 p =0.61

#Tman vs Gravel
TmanVGrav <- cor.test(x=sample_meta_data$X.Gravel,
                      y=log10(sample_sums(Tman)), method = "spearman")
TmanVGrav #rho = 0.37 p =0.02 s = 6680

#Plot of Methylocystaceae against Gravel 
#export as 4.5 in X 4.5 in
plot(log10(sample_sums(Fcys)) ~ sample_meta_data$X.Gravel,
     yaxt = "n", ylab = "",
     ylim =c(-5.5,-1),xlim=c(0,7),
     xlab = expression('Gravel (%)'), pch = 16, cex = 1.5, col="gray") 
points(sample_meta_data$X.Gravel,log10(sample_sums(Fcys)),
       pch = 1, cex = 1.5)
title(ylab=expression("Relative Abundance"),
      mgp=c(2.8,0,1), cex.lab=1)
axis(side = 2, lwd.ticks = 1, cex.axis = 1, las = 1, tick =T,
     labels=expression(10^-1,10^-2,10^-3,10^-4,10^-5),
     at=c(-1,-2,-3,-4,-5))
text(5.15,-4.5, expression(rho), srt =-10) 
text(5.9,-4.5, as.expression("= 0.39")) 
text(5.9,-5, as.expression(italic(S)~"="~10231)) 
