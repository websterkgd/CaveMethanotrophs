#This tests the relationships between the microbial community in the sampled caves to their 
#relationship with distance from other caves

rm(list=ls())
getwd()
setwd('C:/Users/websterkgd/GitHub/CaveMethanotrophs')

#get the citation for the R package
citation(package = "base", lib.loc = NULL)

# loading the phyloseq package for microbial analysis
require('phyloseq')
# loading the ggplot2 package
require("ggplot2")
packageVersion("ggplot2")

# calling vegan
#install.packages("vegan")
require('vegan')
require('ape')
require('seqinr')
require('fossil')
require('simba')
require('picante')

#importing the shared file 
OTU <- import_mothur(mothur_shared_file = 'CvMcrb03.bac.final.shared', parseFunction = parse_taxonomy_default)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file = 'CvMcrb03.bac.final.0.03.taxonomy', parseFunction = parse_taxonomy_default)

#create a physeq table for analysis
physeq = phyloseq(OTU, TAX)

#Creating transformations of the data 
f.physeq = transform_sample_counts(physeq, function(x) x / sum(x)) # fractional abundance transformation of the data

# import meta data
md = read.table('sample_meta_data-d.txt', header = TRUE, sep = "", dec = ".", row.names = 1) 
  
#Let Phyloseq see the meta data
sample_meta_data <- sample_data(md, errorIfNULL = TRUE)

#Turns location into factor data

sample_meta_data$Location <- as.factor(sample_meta_data$Location)
sample_meta_data$Bin.Loc <- as.factor(sample_meta_data$Bin.Loc)

#merge meta data with other data
AllCaves <- merge_phyloseq(f.physeq, sample_meta_data)

#Rank1 =K, Rank2 =P, Rank3 =C, Rank4 =O, Rank5 =F, Rank6 =G, Rank7 =S
#Select the taxa of interest

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

#Proteobacteria Alphaproteobacteria Rhizobiales  Beijernickaceae
Gcel <- subset_taxa(AllCaves, Rank6=="Methylocella") # Genus Gcel

#Pruning absences from matrices to for decay analyses

#pruning Fcys for future community similarity calculations

prune.Fcys <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2c-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5d-Skyline", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "7c-Grand", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BC-II", "BiC-SR", "BiC-TBBV-I", "CDVL-BS", "CDVL-Cas-Peq", "CDVL-Mid-Spr", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Sky-Hi-12m", "TVC-BP", "Woodward-3D")

kpNZ.Fcys <-prune_samples(prune.Fcys, Fcys)

#pruning Fphl for future community similarity calculations

prune.Fphl <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2c-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "6d-Shenandoah", "7c-Grand", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BiC-SR", "BiC-TBBV-I", "BiC-TBBV-II", "CDVL-BS", "CDVL-Mid-Spr", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Sky-Hi-12m", "TVC-BP", "Woodward-3D")

kpNZ.Fphl <-prune_samples(prune.Fphl, Fphl)

#pruning Ococ for future similarity calculations

prune.Ococ <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "6d-Shenandoah", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BiC-SR", "CDVL-BS", "CDVL-Mid-Spr", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "TVC-BP", "Woodward-3D")

kpNZ.Ococ <-prune_samples(prune.Ococ, Ococ)

#pruning Ococ.Fucs for future community similarity calculations

prune.Ococ.Fucs <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5e-Skyline", "6b-Shenandoah", "7f-Grand-II", "7f-Grand-III", "BiC-SR", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "TVC-BP")

kpNZ.Ococ.Fucs <-prune_samples(prune.Ococ.Fucs, Ococ.Fucs)

#pruning Ococ.Fcoc.Gcal for future community similarity calculations

prune.Ococ.Fcoc.Gcal <- c("1c-Lincoln", "6c-Shenandoah","6d-Shenandoah", "BiC-SR", "CDVL-Mid-Spr", "Mnt-Spr-Cv-12m", "SC-LDER", "SC-SLUMR", "ShiC-IBDSP", "ShiC-SP")

kpNZ.Ococ.Fcoc.Gcal <-prune_samples(prune.Ococ.Fcoc.Gcal, Ococ.Fcoc.Gcal)

#pruning Ococ.Fcoc.Gucs for future community similarity calculations

prune.Ococ.Fcoc.Gucs <- c("10d-Bristol-I", "11d-Appalachian", "1b-Lincoln", "1c-Lincoln","1d-Indian", "4c-Lost-River", "5e-Skyline", "6c-Shenandoah", "6d-Shenandoah", "7f-Grand-II", "8d-Bridge", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "TVC-BP", "Woodward-3D")

kpNZ.Ococ.Fcoc.Gucs <-prune_samples(prune.Ococ.Fcoc.Gucs, Ococ.Fcoc.Gucs)

#pruning of Glib for future community similarity calculations

prune.Glib <- c("1c-Lincoln", "9b-Dixie", "Mnt-Spr-Cv-12m", "SBC-DV", "SC-SLUMR", "ShiC-IBDSP", "ShiC-SP", "TVC-BP")

kpNZ.Glib <-prune_samples(prune.Glib, Glib)

#pruning Fbac for future community similarity calculations

prune.Fbac <- c("10d-Bristol-I","11d-Appalachian", "1c-Lincoln", "1d-Indian", "2c-Indian","2e-Indian","4b-Lost-River", "4c-Lost-River", "6c-Shenandoah","7c-Grand","7f-Grand-II","7f-Grand-III","8d-Bridge", "9b-Dixie","BiC-SR", "Mnt-Spr-Cv-12m","SBC-WF","SC-LDER", "SC-SLUMR","ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Woodward-3D")

kpNZ.Fbac <-prune_samples(prune.Fbac, Fbac)

#Total Methylotrophs and methanotrophs

Tmet <- merge_phyloseq(Fcys, Fbac, Ococ, Fphl, Gcel, Glib, V.Caci)

prune.Tmet <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2c-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5d-skyline", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "6d-Shenandoah", "7c-Grand", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BC-II", "BiC-SR", "BiC-TBBV-I", "BiC-TBBV-II", "CDVL-BS", "CDVL-cas-Paq", "CVDL-Mid-Spr", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Sky-Hi-12m", "TVC-BP", "Woodward-3D")

kpNZ.Tmet <-prune_samples(prune.Tmet, Tmet)

#Filtering out CVL samples

#pruning Fcys for future community similarity calculations

prune.Fcys.cvl <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2c-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5d-Skyline", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "7c-Grand", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BC-II", "BiC-SR", "BiC-TBBV-I", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Sky-Hi-12m", "TVC-BP", "Woodward-3D")

kpNZ.Fcys.cvl <-prune_samples(prune.Fcys.cvl, Fcys)

#pruning Fphl for future community similarity calculations

prune.Fphl.cvl <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2c-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "6d-Shenandoah", "7c-Grand", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BiC-SR", "BiC-TBBV-I", "BiC-TBBV-II", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Sky-Hi-12m", "TVC-BP", "Woodward-3D")

kpNZ.Fphl <-prune_samples(prune.Fphl.cvl, Fphl)

#pruning Ococ for future similarity calculations

prune.Ococ.cvl <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "6d-Shenandoah", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BiC-SR", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "TVC-BP", "Woodward-3D")

kpNZ.Ococ.cvl <-prune_samples(prune.Ococ.cvl, Ococ)

#pruning Ococ.Fucs for future community similarity calculations

prune.Ococ.Fucs.cvl <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5e-Skyline", "6b-Shenandoah", "7f-Grand-II", "7f-Grand-III", "BiC-SR", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "TVC-BP")

kpNZ.Ococ.Fucs.cvl <-prune_samples(prune.Ococ.Fucs.cvl, Ococ.Fucs)

#pruning Ococ.Fcoc.Gcal for future community similarity calculations

prune.Ococ.Fcoc.Gcal.cvl <- c("1c-Lincoln", "6c-Shenandoah","6d-Shenandoah", "BiC-SR", "CDVL-Mid-Spr", "Mnt-Spr-Cv-12m", "SC-LDER", "SC-SLUMR", "ShiC-IBDSP", "ShiC-SP")

kpNZ.Ococ.Fcoc.Gcal <-prune_samples(prune.Ococ.Fcoc.Gcal.cvl, Ococ.Fcoc.Gcal)

#pruning Ococ.Fcoc.Gucs for future community similarity calculations

prune.Ococ.Fcoc.Gucs.cvl <- c("10d-Bristol-I", "11d-Appalachian", "1b-Lincoln", "1c-Lincoln","1d-Indian", "4c-Lost-River", "5e-Skyline", "6c-Shenandoah", "6d-Shenandoah", "7f-Grand-II", "8d-Bridge", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "TVC-BP", "Woodward-3D")

kpNZ.Ococ.Fcoc.Gucs.cvl <-prune_samples(prune.Ococ.Fcoc.Gucs, Ococ.Fcoc.Gucs)

#pruning of Glib for future community similarity calculations

prune.Glib.cvl <- c("1c-Lincoln", "9b-Dixie", "Mnt-Spr-Cv-12m", "SBC-DV", "SC-SLUMR", "ShiC-IBDSP", "ShiC-SP", "TVC-BP")

kpNZ.Glib.cvl <-prune_samples(prune.Glib.cvl, Glib)

#pruning Fbac for future community similarity calculations

prune.Fbac.cvl <- c("10d-Bristol-I","11d-Appalachian", "1c-Lincoln", "1d-Indian", "2c-Indian","2e-Indian","4b-Lost-River", "4c-Lost-River", "6c-Shenandoah","7c-Grand","7f-Grand-II","7f-Grand-III","8d-Bridge", "9b-Dixie","BiC-SR", "Mnt-Spr-Cv-12m","SBC-WF","SC-LDER", "SC-SLUMR","ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Woodward-3D")

kpNZ.Fbac.cvl <-prune_samples(prune.Fbac.cvl, Fbac)

#Total Methylotrophs and methanotrophs

Tmet <- merge_phyloseq(Fcys, Fbac, Ococ, Fphl, Gcel, Glib, V.Caci)

prune.Tmet.cvl <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2c-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5d-skyline", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "6d-Shenandoah", "7c-Grand", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BC-II", "BiC-SR", "BiC-TBBV-I", "BiC-TBBV-II", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Sky-Hi-12m", "TVC-BP", "Woodward-3D")

kpNZ.Tmet.cvl <-prune_samples(prune.Tmet.cvl, Tmet)

# pulling the total methanotrophic community
Tman <- merge_phyloseq(Fcys, Ococ, Gcel) 

prune.Tman.cvl <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", "1d-Indian", "2c-Indian", "2e-Indian", "4b-Lost-River", "4c-Lost-River", "5d-skyline", "5e-Skyline", "6b-Shenandoah", "6c-Shenandoah", "6d-Shenandoah", "7c-Grand", "7f-Grand-II", "7f-Grand-III", "8d-Bridge", "9b-Dixie", "BC-II", "BiC-SR", "BiC-TBBV-I", "JC-ER-WS", "Mnt-Spr-Cv-12m", "SBC-DV", "SBC-WF", "SC-LDER", "SC-SLUMR", "ShiC-BDWF", "ShiC-IBDSP", "ShiC-SP", "Sky-Hi-12m", "TVC-BP", "Woodward-3D")

kpNZ.Tman.cvl <-prune_samples(prune.Tman.cvl, Tman)

####### Distance decay for methanotrophs
Tman.lat.cvl <- as.numeric(get_variable(kpNZ.Tman.cvl, "Lattitude"))
Tman.lon.cvl <- as.numeric(get_variable(kpNZ.Tman.cvl, "Longitude"))

Tman.struc.dist.cvl <- phyloseq::distance(kpNZ.Tman.cvl, "bray")
Tman.coord.dist.cvl <- dist(as.matrix(Tman.lat.cvl, 
                                      Tman.lon.cvl)) #geographical distance
Tman.coord.dist.cvl <- Tman.coord.dist.cvl*111.325 #degrees to km

#transform environmental variables into numeric data
Tman.CH4.cvl <- as.numeric(get_variable(kpNZ.Tman.cvl, "CH4_conc.ppm."))
Tman.CO2.cvl <- as.numeric(get_variable(kpNZ.Tman.cvl, "CO2_conc.ppm."))
Tman.loc.cvl <- as.numeric(get_variable(kpNZ.Tman.cvl, "Location"))

#calculate the euclidean distance between the samples in enivronmental variables

Tman.env.dist.cvl <- 1 - vegdist(cbind(Tman.CH4.cvl, Tman.CO2.cvl, 
                                      Tman.loc.cvl))

#transform distance matrices into lists
Tman.struc.dist.ls.cvl <- liste(Tman.struc.dist.cvl, entry="Tman.struc.cvl")
Tman.env.dist.ls.cvl <- liste(Tman.env.dist.cvl, entry ="Tman.env.cvl")
Tman.coord.dist.ls.cvl <- liste(Tman.coord.dist.cvl, entry = "Tman.dist.cvl")

#create a dataframe containing the similarity of the environment and the community
Tman.df.cvl <- data.frame(Tman.coord.dist.ls.cvl, 
                         Tman.env.dist.ls.cvl[3], Tman.struc.dist.ls.cvl[3])
names(Tman.df.cvl)[4:5] <- c("Tman.env.cvl", "Tman.struc.cvl")
attach(Tman.df.cvl)

#plot the distance decay for environmental pattern vs distance

par(mfrow = c(1,1))
plot(Tman.dist.cvl, Tman.env.cvl, xlab = "Geographic Distance", ylab = "Environmental Distance", main = "Distance-Decay for the Environment", col="SteelBlue")

Tman.DD_CA_e.cvl <- lm(Tman.env.cvl ~ Tman.dist.cvl)
Tman.DD_CA_e.cvl

abline(Tman.DD_CA_e.cvl)

par(mfrow = c(1,1))
plot(Tman.dist.cvl, Tman.struc.cvl, xlab = "Relative Geographic Distance", 
     ylab = "Community Distance", 
     main = "Distance-Decay for Methanotrophs", col="SteelBlue")

Tman.DD_CA_s.cvl <- lm(Tman.struc.cvl ~ Tman.dist.cvl)
summary(Tman.DD_CA_s.cvl)
abline(Tman.DD_CA_s.cvl)

mantel(Tman.coord.dist.cvl,Tman.struc.dist.cvl, permutations =999)
#r 0.08 p 0.14

diffslope(Tman.dist.cvl, Tman.env.cvl, Tman.dist.cvl, Tman.struc.cvl) # not significant
