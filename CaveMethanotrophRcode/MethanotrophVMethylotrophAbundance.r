#This script analyses the relative abundance of methylotrophs to methanotrophs
# and of methylotrophs to methane concenta

rm(list=ls())
getwd()
setwd(getwd())

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
require('stringi')

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

#Pruning NAs in the "design matrix" for analysis of distance on microbial communities
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

#Total methanotrophic community
Tman <- merge_phyloseq(Fcys, Ococ, Gcel) # pulling the total methanotrophic community

#Total methylotrophic community
Tlot <- merge_phyloseq(Fphl, Glib, Fbac) # pulling the total methylotrophic community

#Methylotroph abundance vs Methanotroph Abundance
#Vector of the relative abundance of methylotrophs and methanotrophs
Tlot.ss <- sample_sums(Tlot)
Tman.ss <- sample_sums(Tman)

#using a spearman's test to determine if the function is monotonic
corr <- cor.test(x=Tman.ss, y=Tlot.ss, method = 'spearman')
#S = 3849.1 rho = 0.71 p =1e-7; function is monotonic

###Creating plot of methanotroph abundance against methylotroph abundance
#Export as 5in x 5in
plot(Tlot.ss ~ Tman.ss,
     xlab ="Relative Abundance Methanotrophs",
     ylab ="Relative Abundance Methylotrophs",
     pch=16, col="gray", cex=1.5)
points(Tman.ss, Tlot.ss,
       pch = 1, cex = 1.5, col = "black")
text(0.0166,0.00395, expression(rho), srt =-10) 
text(0.0189,0.004, as.expression("= 0.71")) 
text(0.0183,0.0034, as.expression(italic(S)~"="~3849)) 

#Methylotrophs against CH4 concentration
Tlot.ss.fz <- Tlot.ss #create a new vector to filter zeros
Tlot.ss.fz <- as.matrix(Tlot.ss.fz)
Tlot.ss.fz[c(10,11,21,22,27,30,32,43)] <- NA
CH4.ml.f <-md$CH4_conc.ppm.
CH4.ml.f <- as.matrix(CH4.ml.f)
CH4.ml.f[c(10,11,21,22,27,30,32,43)] <- NA
mod.Tlot <- summary(lm(log10(Tlot.ss.fz)~CH4.ml.f)) #p0.001 r2 =0.27
corr.mlvCH4 <- cor.test(x=CH4.ml.f, y=Tlot.ss.fz, method = 'spearman')
#S = 8913.8, p = 0.15, rho = -0.24
