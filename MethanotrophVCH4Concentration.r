#This r code plots the abundance of the methylocystaceae and methylococcales
#against methane concentration

rm(list=ls())
getwd()
setwd("C:/Users/websterkgd/GitHub/cave-mob/data_analyses")

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
md = read.table('sample_meta_data-d.txt', header = TRUE, sep = "", dec = ".", row.names = 1) 
  
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



## Numeric plots

ss.Fcys <-sample_sums(Fcys)
ss.Fcys.f <- c(ss.Fcys[1:9],ss.Fcys[11:20],ss.Fcys[22:24],ss.Fcys[26:29],
               ss.Fcys[31],ss.Fcys[33:42])

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:9], 
               sample_meta_data$CH4_conc.ppm.[11:20],
               sample_meta_data$CH4_conc.ppm.[22:24],
               sample_meta_data$CH4_conc.ppm.[26:29], 
               sample_meta_data$CH4_conc.ppm.[31], 
               sample_meta_data$CH4_conc.ppm.[33:42])

ss.Ococ <-sample_sums(Ococ)
ss.Ococ.f <- c(ss.Ococ[1:5],ss.Ococ[7:9],ss.Ococ[12:15],ss.Ococ[17:20],ss.Ococ[23],
               ss.Ococ[26],ss.Ococ[28:29],ss.Ococ[31],ss.Ococ[33:39],
               ss.Ococ[41:42])

smd.O.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:5], 
                 sample_meta_data$CH4_conc.ppm.[7:9],
                 sample_meta_data$CH4_conc.ppm.[12:15],
                 sample_meta_data$CH4_conc.ppm.[17:20], 
                 sample_meta_data$CH4_conc.ppm.[23], 
                 sample_meta_data$CH4_conc.ppm.[26],
                 sample_meta_data$CH4_conc.ppm.[28:29],
                 sample_meta_data$CH4_conc.ppm.[31],
                 sample_meta_data$CH4_conc.ppm.[33:39],
                 sample_meta_data$CH4_conc.ppm.[41:42])
 
plot(log10(ss.Ococ.f)~smd.O.CH4.F) # re filter. 
plot(log10(ss.Fcys.f)~smd.CH4.F)

#export as 5.5in X 5.5in
plot(smd.O.CH4.F,log10(ss.Ococ.f), yaxt = "n", ylab = "",
     ylim =c(-6,-1),xlim=c(0,4),
     xlab = expression('CH'[4]*' Concentration (ppmv)'),
     col="gray", pch=16,cex =1.3) 
points(smd.O.CH4.F,log10(ss.Ococ.f), col='black',pch=1, cex=1.3) 
points(smd.CH4.F,log10(ss.Fcys.f), col='black',pch=16, cex=1.3) 
title(ylab=expression("Relative Abundance"),
      mgp=c(2.8,0,1), cex.lab=1)
axis(side = 2, lwd.ticks = 1, cex.axis = 1, las = 1, tick =T,
     labels=expression(10^-1,10^-2,10^-3,10^-4,10^-5,10^-6),
     at=c(-1,-2,-3,-4,-5,-6))
legend(2, -0.5, legend=c("Methylocystaceae", "Methylococcales"),
       pch=c(16,16), col=c("black","gray"), cex=1, bty="n")
legend(2, -0.5, legend=c("Methylocystaceae", "Methylococcales"),
       pch=c(16,1), cex=1, bty="n")
abline(lm(log10(ss.Ococ.f)~smd.O.CH4.F),lty=2)
abline(lm(log10(ss.Fcys.f)~smd.CH4.F),lty=2)

###Ancova Analysis for slopes

mv.u.Fcys <-sample_sums(Fcys)
mv.u.Fcys.f <- c(mv.u.Fcys[1:31],mv.u.Fcys[33:43])
l.mv.u.Fcys.f <-log1p(mv.u.Fcys.f)

mv.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:31],
              sample_meta_data$CH4_conc.ppm.[33:43])

mv.u.Ococ <-sample_sums(Ococ)
mv.u.Ococ.f <- c(mv.u.Ococ[1:31],mv.u.Ococ[33:43])
l.mv.u.Ococ.f <-log1p(mv.u.Ococ.f)

mv.O.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:31],
                sample_meta_data$CH4_conc.ppm.[33:43])

#rel abund ~ CH4 + type + (CH4 x type)

RA <- c(as.numeric(l.mv.u.Ococ.f),as.numeric(l.mv.u.Fcys.f))
mn <- c(as.numeric(mv.O.CH4.F),as.numeric(mv.CH4.F))
tp <- c(1:84)
tp[1:42] <-1
tp[43:84] <-2

#ancova for slopes
mod1 <- aov(RA~mn*tp)
mod2 <- aov(RA~mn+tp)
anova(mod1,mod2) # f =1.6 p =0.21

#########