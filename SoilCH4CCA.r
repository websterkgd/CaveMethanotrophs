#CCA analysis of methane+Co2_soil vs mehtanotroph community
#JC_Wr and SC-SLUMR is Filtered from all analyses because of NAs in soil

rm(list=ls())
setwd('C:/Users/websterkgd/GitHub/CaveMethanotrophs/')

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
require('vegan')
require('car')

#importing the shared file 
OTU <- import_mothur(mothur_shared_file = 'CvMcrb03.bac.final.shared', parseFunction = parse_taxonomy_default)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file = 'CvMcrb03.bac.final.0.03.taxonomy', parseFunction = parse_taxonomy_default)

#create a physeq table for analysis
physeq = phyloseq(OTU, TAX)

#Creating transformations of the data
f.physeq = transform_sample_counts(physeq, function(x) x / sum(x)) 
# fractional abundance transformation of the data

# import mineral data
md = read.table('sample_meta_data+d.txt', header = TRUE, sep = "", dec = ".", row.names = 1) 

#Let Phyloseq see the meta data
sample_meta_data <- sample_data(md, errorIfNULL = TRUE)

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

#pruning Fcys for future community similarity calculations
prune.Fcys <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln", 
                "1d-Indian","2c-Indian","2e-Indian","4b-Lost-River","4c-Lost-River", 
                "5d-Skyline","5e-Skyline","6b-Shenandoah","6c-Shenandoah","7c-Grand", 
                "7f-Grand-II","7f-Grand-III","8d-Bridge","9b-Dixie","BC-II","BiC-SR", 
                "BiC-TBBV-I","CDVL-BS","CDVL-Cas-Peq","CDVL-Mid-Spr","JC-ER-WS",
                "Mnt-Spr-Cv-12m","SBC-DV","SBC-WF","SC-LDER","ShiC-BDWF",
                "ShiC-IBDSP","ShiC-SP","Sky-Hi-12m","TVC-BP", "Woodward-3D")

kpNZ.Fcys <-prune_samples(prune.Fcys, Fcys)

#pruning Fphl for future community similarity calculations

prune.Fphl <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln","1d-Indian", 
                "2c-Indian","2e-Indian","4b-Lost-River","4c-Lost-River","5e-Skyline", 
                "6b-Shenandoah","6c-Shenandoah","6d-Shenandoah","7c-Grand","7f-Grand-II", 
                "7f-Grand-III","8d-Bridge","9b-Dixie","BiC-SR","BiC-TBBV-I","BiC-TBBV-II", 
                "CDVL-BS","CDVL-Mid-Spr","JC-ER-WS","Mnt-Spr-Cv-12m","SBC-DV","SBC-WF", 
                "SC-LDER","ShiC-BDWF","ShiC-IBDSP","ShiC-SP","Sky-Hi-12m","TVC-BP", 
                "Woodward-3D")

kpNZ.Fphl <-prune_samples(prune.Fphl, Fphl)

#pruning Ococ for future similarity calculations

prune.Ococ <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln","1d-Indian", 
                "2e-Indian","4b-Lost-River","4c-Lost-River","5e-Skyline","6b-Shenandoah", 
                "6c-Shenandoah","6d-Shenandoah","7f-Grand-II","7f-Grand-III","8d-Bridge", 
                "9b-Dixie","BiC-SR","CDVL-BS","CDVL-Mid-Spr","JC-ER-WS","Mnt-Spr-Cv-12m", 
                "SBC-DV","SBC-WF","SC-LDER","ShiC-BDWF","ShiC-IBDSP","ShiC-SP","TVC-BP", 
                "Woodward-3D")

kpNZ.Ococ <-prune_samples(prune.Ococ, Ococ)

#pruning Ococ.Fucs for future community similarity calculations

prune.Ococ.Fucs <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln","1d-Indian", 
                     "2e-Indian","4b-Lost-River","4c-Lost-River","5e-Skyline","6b-Shenandoah", 
                     "7f-Grand-II","7f-Grand-III","BiC-SR","JC-ER-WS","Mnt-Spr-Cv-12m","SBC-WF", 
                     "SC-LDER","ShiC-BDWF","TVC-BP")

kpNZ.Ococ.Fucs <-prune_samples(prune.Ococ.Fucs, Ococ.Fucs)

#pruning Ococ.Fcoc.Gcal for future community similarity calculations

prune.Ococ.Fcoc.Gcal <- c("1c-Lincoln","6c-Shenandoah","6d-Shenandoah","BiC-SR","CDVL-Mid-Spr", 
                          "Mnt-Spr-Cv-12m","SC-LDER","ShiC-IBDSP","ShiC-SP")

kpNZ.Ococ.Fcoc.Gcal <-prune_samples(prune.Ococ.Fcoc.Gcal, Ococ.Fcoc.Gcal)

#pruning Ococ.Fcoc.Gucs for future community similarity calculations

prune.Ococ.Fcoc.Gucs <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln","1d-Indian", 
                          "4c-Lost-River","5e-Skyline","6c-Shenandoah","6d-Shenandoah","7f-Grand-II", 
                          "8d-Bridge","Mnt-Spr-Cv-12m","SBC-DV","SBC-WF","SC-LDER","ShiC-BDWF", 
                          "ShiC-IBDSP","ShiC-SP","TVC-BP","Woodward-3D")

kpNZ.Ococ.Fcoc.Gucs <-prune_samples(prune.Ococ.Fcoc.Gucs, Ococ.Fcoc.Gucs)

#pruning of Glib for future community similarity calculations

prune.Glib <- c("1c-Lincoln","9b-Dixie","Mnt-Spr-Cv-12m","SBC-DV","ShiC-IBDSP","ShiC-SP","TVC-BP")

kpNZ.Glib <-prune_samples(prune.Glib, Glib)

#pruning Fbac for future community similarity calculations

prune.Fbac <- c("10d-Bristol-I","11d-Appalachian","1c-Lincoln","1d-Indian","2c-Indian","2e-Indian",
                "4b-Lost-River","4c-Lost-River","6c-Shenandoah","7c-Grand","7f-Grand-II","7f-Grand-III",
                "8d-Bridge","9b-Dixie","BiC-SR","Mnt-Spr-Cv-12m","SBC-WF","ShiC-BDWF", 
                "ShiC-IBDSP","ShiC-SP","Woodward-3D")

kpNZ.Fbac <-prune_samples(prune.Fbac, Fbac)

# pruning the Beijerinckiaceae

prune.Fbei <-c("10d-Bristol-I","2c-Indian","2e-Indian","6c-Shenandoah","7f-Grand-II","7f-Grand-III","9b-Dixie", 
               "BC-II","BiC-SR","CDVL-Mid-Spr","Mnt-Spr-Cv-12m","SBC-DV","SBC-WF","SC-LDER","ShiC-IBDSP", 
               "ShiC-SP","Woodward-3D")

prune.Fbei.Guct <-c("10d-Bristol-I","2c-Indian","2e-Indian","6c-Shenandoah","7f-Grand-II","7f-Grand-III", 
                    "9b-Dixie","BC-II","BiC-SR","CDVL-Mid-Spr","Mnt-Spr-Cv-12m","SBC-DV","SBC-WF","SC-LDER", 
                    "ShiC-IBDSP", "ShiC-SP","Woodward-3D")

kpNZ.Fbei <-prune_samples(prune.Fbei, Fbei)

Tmet <- merge_phyloseq(Fcys, Fbac, Ococ, Fphl, Gcel, Glib, V.Caci)

prune.Tmet <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln","1d-Indian","2c-Indian","2e-Indian", 
                "4b-Lost-River","4c-Lost-River","5d-skyline","5e-Skyline","6b-Shenandoah","6c-Shenandoah", 
                "6d-Shenandoah","7c-Grand","7f-Grand-II","7f-Grand-III","8d-Bridge","9b-Dixie","BC-II","BiC-SR", 
                "BiC-TBBV-I","BiC-TBBV-II","CDVL-BS","CDVL-cas-Paq","CVDL-Mid-Spr","JC-ER-WS","Mnt-Spr-Cv-12m", 
                "SBC-DV","SBC-WF","SC-LDER","ShiC-BDWF","ShiC-IBDSP","ShiC-SP","Sky-Hi-12m","TVC-BP", 
                "Woodward-3D")

kpNZ.Tmet <-prune_samples(prune.Tmet, Tmet)

Tman <- merge_phyloseq(Fcys, Ococ, Gcel) # pulling the total methanotrophic community

prune.Tman <- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln","1d-Indian","2c-Indian", 
                "2e-Indian","4b-Lost-River","4c-Lost-River","5d-skyline","5e-Skyline","6b-Shenandoah", 
                "6c-Shenandoah","6d-Shenandoah","7c-Grand","7f-Grand-II","7f-Grand-III","8d-Bridge", 
                "9b-Dixie","BC-II","BiC-SR","BiC-TBBV-I","CDVL-BS","CDVL-Cas-Peq","CDVL-Mid-Spr","JC-ER-WS", 
                "Mnt-Spr-Cv-12m","SBC-DV","SBC-WF","SC-LDER","ShiC-BDWF","ShiC-IBDSP","ShiC-SP","Sky-Hi-12m", 
                "TVC-BP","Woodward-3D")

kpNZ.Tman <-prune_samples(prune.Tman, Tman)

Tloc <- merge_phyloseq(Fbac, Fphl, Glib) # total methylotrophs

prune.Tloc<- c("10d-Bristol-I","11d-Appalachian","1b-Lincoln","1c-Lincoln","1d-Indian","2c-Indian", 
               "2e-Indian","4b-Lost-River","4c-Lost-River","5e-Skyline","6b-Shenandoah","6c-Shenandoah", 
               "6d-Shenandoah","7c-Grand","7f-Grand-II","7f-Grand-III","8d-Bridge","9b-Dixie","BiC-SR", 
               "BiC-TBBV-I","BiC-TBBV-II","CDVL-BS","CVDL-Mid-Spr","JC-ER-WS","Mnt-Spr-Cv-12m","SBC-DV", 
               "SBC-WF","SC-LDER","ShiC-BDWF","ShiC-IBDSP","ShiC-SP","Sky-Hi-12m","TVC-BP", 
               "Woodward-3D")

kpNZ.Tloc <-prune_samples(prune.Tloc, Tloc)

#####CCA
##Methylococcales OTU
kpNZ.Ococ.env <- as.data.frame(cbind(get_variable(kpNZ.Ococ, "CH4_conc.ppm."), 
                                     get_variable(kpNZ.Ococ, "CO2_conc.ppm."),
                                     get_variable(kpNZ.Ococ, "X.Gravel"),
                                     get_variable(kpNZ.Ococ, "X.Sand"),
                                     get_variable(kpNZ.Ococ, "X.Silt"),
                                     get_variable(kpNZ.Ococ, "X.Clay")))
names(kpNZ.Ococ.env) <- c("CH4","CO2","Gravel","Sand","Silt","Clay")
kpNZ.Ococ.env <-as.matrix(kpNZ.Ococ.env)

Ococ.cca <- vegan::cca(t(kpNZ.Ococ@otu_table)~kpNZ.Ococ.env) # 
anova(Ococ.cca, by ="axis") #p0.06 
Ococ.cca.fit <- envfit(Ococ.cca, kpNZ.Ococ.env, perm=999) 
Ococ.cca.fit #p0.001 # Methanococcales related to CH4 and to clay

#calculate explained variance
Ococ.cca.explainvar1 <- round(Ococ.cca$CCA$eig[1]/
                                sum(c(Ococ.cca$CCA$eig, Ococ.cca$CA$eig)), 3) * 100
Ococ.cca.explainvar2 <- round(Ococ.cca$CCA$eig[2]/
                                sum(c(Ococ.cca$CCA$eig, Ococ.cca$CA$eig)), 3) * 100

#define plot parameters 

par(mar = c(5,5,4,4)+0.1)

#intiate plot
plot(scores(Ococ.cca, display = "wa"), xlim = c(-4,3), ylim=c(-4, 4),
     xlab = paste("CCA 1(", Ococ.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2(", Ococ.cca.explainvar2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

#add axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

#add points & labels
points(scores(Ococ.cca, display = "wa"),
       pch = 19, cex = 3, bg = "gray", col = "gray")

#Add environmental vectors
vectors <- scores(Ococ.cca, display = "bp")
row.names(vectors) <- c("CH4","CO2","Gravel","Sand","Silt","Clay")
arrows(0,0, vectors[,1]*2, vectors[,2]*2,
       lwd=2, lty=1, length = 0.2, col = "red")
text(vectors[,1]*2,vectors[,2]*2, pos = 3, labels = row.names(vectors), col ="red")
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, col ="red", lwd = 2.2,
     at = pretty(range(vectors[,1])) *2, labels = pretty(range(vectors[,1])))
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, col ="red", lwd = 2.2,
     at = pretty(range(vectors[,1])) *2, labels = pretty(range(vectors[,1])))

##Methylobacteriaceae OTU
kpNZ.Fbac.env <- as.data.frame(cbind(get_variable(kpNZ.Fbac, "CH4_conc.ppm."), 
                                     get_variable(kpNZ.Fbac, "CO2_conc.ppm."),
                                     get_variable(kpNZ.Fbac, "X.Gravel"),
                                     get_variable(kpNZ.Fbac, "X.Sand"),
                                     get_variable(kpNZ.Fbac, "X.Silt"),
                                     get_variable(kpNZ.Fbac, "X.Clay")))
names(kpNZ.Fbac.env) <- c("CH4","CO2","Gravel","Sand","Silt","Clay")
kpNZ.Fbac.env <-as.matrix(kpNZ.Fbac.env)

Fbac.cca <- vegan::cca(t(kpNZ.Fbac@otu_table)~kpNZ.Fbac.env) #! It ran! 
anova(Fbac.cca, by ="axis") #CCA1 is significant 0.02 very similar to CH4
Fbac.cca.fit <- envfit(Fbac.cca, kpNZ.Fbac.env, perm=999) 
Fbac.cca.fit # Methylbacetiaceae related to CH4 p0.02 and silt p0.04

#calculate explained variance
Fbac.cca.explainvar1 <- round(Fbac.cca$CCA$eig[1]/
                                sum(c(Fbac.cca$CCA$eig, Fbac.cca$CA$eig)), 3) * 100
Fbac.cca.explainvar2 <- round(Fbac.cca$CCA$eig[2]/
                                sum(c(Fbac.cca$CCA$eig, Fbac.cca$CA$eig)), 3) * 100

#define plot parameters 

par(mar = c(5,5,4,4)+0.1)

#intiate plot
plot(scores(Fbac.cca, display = "wa"), xlim = c(-4,3), ylim=c(-4, 4),
     xlab = paste("CCA 1(", Fbac.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2(", Fbac.cca.explainvar2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

#add axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

#add points & labels
points(scores(Fbac.cca, display = "wa"),
       pch = 19, cex = 2, bg = "gray", col = "gray")

#Add environmental vectors
vectors <- scores(Fbac.cca, display = "bp")
row.names(vectors) <- c("CH4","CO2","Gravel","Sand","Silt","Clay")
arrows(0,0, vectors[,1]*2, vectors[,2]*2,
       lwd=2, lty=1, length = 0.2, col = "red")
text(vectors[,1]*2,vectors[,2]*2, pos = 3, labels = row.names(vectors), col ="red")
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, col ="red", lwd = 2.2,
     at = pretty(range(vectors[,1])) *2, labels = pretty(range(vectors[,1])))
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, col ="red", lwd = 2.2,
     at = pretty(range(vectors[,1])) *2, labels = pretty(range(vectors[,1])))


##Methylophilales OTU
kpNZ.Fphl.env <- as.data.frame(cbind(get_variable(kpNZ.Fphl, "CH4_conc.ppm."), 
                                     get_variable(kpNZ.Fphl, "CO2_conc.ppm."),
                                     get_variable(kpNZ.Fphl, "X.Gravel"),
                                     get_variable(kpNZ.Fphl, "X.Sand"),
                                     get_variable(kpNZ.Fphl, "X.Silt"),
                                     get_variable(kpNZ.Fphl, "X.Clay")))
names(kpNZ.Fphl.env) <- c("CH4","CO2","Gravel","Sand","Silt","Clay")
kpNZ.Fphl.env <-as.matrix(kpNZ.Fphl.env)

Fphl.cca <- vegan::cca(t(kpNZ.Fphl@otu_table)~kpNZ.Fphl.env) #! It ran! 
anova(Fphl.cca, by ="axis") #no CCA axis is sig
Fphl.cca.fit <- envfit(Fphl.cca, kpNZ.Fphl.env, perm=999) 
Fphl.cca.fit #p0.03 clay

#calculate explained variance
Fphl.cca.explainvar1 <- round(Fphl.cca$CCA$eig[1]/
                                sum(c(Fphl.cca$CCA$eig, Fphl.cca$CA$eig)), 3) * 100
Fphl.cca.explainvar2 <- round(Fphl.cca$CCA$eig[2]/
                                sum(c(Fphl.cca$CCA$eig, Fphl.cca$CA$eig)), 3) * 100

#define plot parameters 

par(mar = c(5,5,4,4)+0.1)

#intiate plot
plot(scores(Fphl.cca, display = "wa"), xlim = c(-4,3), ylim=c(-4, 4),
     xlab = paste("CCA 1(", Fphl.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2(", Fphl.cca.explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

#add axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

#add points & labels
points(scores(Fphl.cca, display = "wa"),
       pch = 19, cex = 2, bg = "gray", col = "gray")

#Add environmental vectors
vectors <- scores(Fphl.cca, display = "bp")
row.names(vectors) <- c("CH4","CO2","Gravel","Sand","Silt","Clay")
arrows(0,0, vectors[,1]*2, vectors[,2]*2,
       lwd=2, lty=1, length = 0.2, col = "red")
text(vectors[,1]*2,vectors[,2]*2, pos = 3, labels = row.names(vectors), col ="red")
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, col ="red", lwd = 2.2,
     at = pretty(range(vectors[,1])) *2, labels = pretty(range(vectors[,1])))
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, col ="red", lwd = 2.2,
     at = pretty(range(vectors[,1])) *2, labels = pretty(range(vectors[,1])))


############## plot #export as 5 x 5
plot(scores(Ococ.cca, display = "wa"), xlim =c(-8, 2), ylim=c(-4, 8),
     xlab = paste("CCA 1 (", Ococ.cca.explainvar1, " %)", sep = ""), 
     cex.main = 2, ylab = "",
     pch = 16, cex = 1.4, type = "n", cex.lab = 1.4, cex.axis = 1.4, axes = FALSE)
title(ylab=paste("CCA 2 (", Ococ.cca.explainvar2, " %)", sep = ""), cex.lab = 1.4, line=3)

#add axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

#add points & labels
points(scores(Ococ.cca, display = "wa"),
       pch = 16, cex = 1.5, bg = "gray", col = "gray")

#Add environmental vectors
vectors <- scores(Ococ.cca, display = "bp")
arrows(0,0, vectors[c(1,6),1]*4, vectors[c(1,6),2]*4,
       lwd=2, lty=1, length = 0.2, col = "black")
text(vectors[c(1,6),1]*5.1,vectors[c(1,6),2]*5, 
     labels=expression("CH"[4],"Clay"), col ="black", cex = 1.5)

