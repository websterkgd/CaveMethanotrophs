#CCA analysis of methanotrophs and the mineral assemblage

rm(list=ls())
getwd()

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

#get the citation for the R package
citation(package = "base", lib.loc = NULL)

# loading the phyloseq package for microbial analysis
library("phyloseq")
packageVersion("phyloseq")

#importing the shared file 
OTU <- import_mothur(mothur_shared_file = 'CvMcrb03_20190225.shared', 
                     parseFunction = parse_taxonomy_default)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file = 'CvMcrb03_20190225.taxonomy',
                     parseFunction = parse_taxonomy_default)

#create a physeq table for analysis
physeq = phyloseq(OTU, TAX)

#Creating transformations of the data 
# fractional abundance transformation of the data
f.physeq = transform_sample_counts(physeq, function(x) x / sum(x)) 

# import meta data #just mineral da
md = read.table('MinAbnMatr.txt', header = TRUE, sep = "", 
                dec = ".", row.names = 1) 

#restricting mineral data to numeric data
md.n <- md[,-1]

#Let Phyloseq see the meta data
sample_meta_data <- sample_data(md.n, errorIfNULL = TRUE)

#merge meta data with other data
AllCaves <- merge_phyloseq(f.physeq, sample_meta_data)

#Importing Identified Methanotrophs
Id.m <- read.table('NCBIPutMntrphIDCave_20220616.txt',
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

m.tax[,7] <- paste(as.factor(Id.m[c(1:20,23:25,27:36,38:58,60:64,66:71,73:82),2]))

m.tax <- tax_table(m.tax)

#removing samples with no methanotrophs
m.otu@.Data <-m.otu@.Data[,-c(9)]

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
LA.M <- merge_phyloseq(T.II, Tp.I)

###CCA Total Methanotrophs
tm.d <- as.data.frame(C.Mtf@sam_data@.Data)
rownames(tm.d) <- colnames(C.Mtf@otu_table@.Data)
colnames(tm.d) <- colnames(md.n)
tm.d <- as.matrix(tm.d)

C.Mtf.cca <- vegan::cca(t(C.Mtf@otu_table)[,-c(9,11:22)]~
                          tm.d[,-10]) 
anova(C.Mtf.cca, by ="axis") # nothing sig
C.Mtf.cca.fit <- envfit(C.Mtf.cca, tm.d[,-10], perm=9999) 
C.Mtf.cca.fit 
#clinochlorite 0.008
#muscovite 0.03
#microcline 0.007

summary(C.Mtf.cca)

require('CCA')

#calculate explained variance
C.Mtf.cca.explainvar1 <- 
  round(C.Mtf.cca$CCA$eig[1]/
          sum(c(C.Mtf.cca$CCA$eig, C.Mtf.cca$CA$eig)), 3) * 100
#34.2

C.Mtf.cca.explainvar2 <- 
  round(C.Mtf.cca$CCA$eig[2]/
          sum(c(C.Mtf.cca$CCA$eig, C.Mtf.cca$CA$eig)), 3) * 100
#20.8

#define plot parameters 
par(mar = c(5,5,4,4)+0.1)

#intiate plot
plot(scores(C.Mtf.cca, display = "wa"), xlim = c(-1.5,3), ylim=c(-2, 8),
     xlab = paste("CCA 1 (", C.Mtf.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2 (", C.Mtf.cca.explainvar2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, 
     axes = FALSE)

#plot points
points(scores(C.Mtf.cca, display = "wa"),
       pch = 19, cex = 2, bg = "gray", col = "gray")

#add axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

#Add environmental vectors
vectors <- scores(C.Mtf.cca, display = "bp")
arrows(0,0, vectors[c(2,3,5,7),1]*1.2, vectors[c(2,3,5,7),2]*1.2,
       lwd=2, lty=1, length = 0.2, col = "red")
text(vectors[c(2,3,5,7),1]*1.3,vectors[c(2,3,5,7),2]*1.3, pos = 3, 
     labels = c("Muscovite","Clinochlor",
                "Orthoclase","Microcline"), col ="red")

#Low affinity
#filter samples with 0 T2s
LA.M <- subset_samples(LA.M, sample_names(LA.M) != "10d-Bristol-I")
LA.M <- subset_samples(LA.M, sample_names(LA.M) != "5b-Skyline")

#get mineral data
la.d <- as.data.frame(LA.M@sam_data@.Data)
rownames(la.d) <- colnames(LA.M@otu_table@.Data)
colnames(la.d) <- colnames(md.n)
la.d <- as.matrix(la.d)

LA.M.cca <- vegan::cca(t(LA.M@otu_table)~
                          la.d) #
anova(LA.M.cca, by ="axis") # CCA1 0.001
LA.M.cca.fit <- envfit(LA.M.cca, la.d, perm=9999)
LA.M.cca.fit
#Quartz at ~ p = 0.05

#calculate explained variance
LA.M.cca.explainvar1 <- round(LA.M.cca$CCA$eig[1]/
                                 sum(c(LA.M.cca$CCA$eig, LA.M.cca$CA$eig)), 3) * 100
#42.9

LA.M.cca.explainvar2 <- round(LA.M.cca$CCA$eig[2]/
                                 sum(c(LA.M.cca$CCA$eig, LA.M.cca$CA$eig)), 3) * 100
#8.3



#Type USC

#get mineral data
tu.d <- as.data.frame(TUSC@sam_data@.Data)
rownames(tu.d) <- colnames(TUSC@otu_table@.Data)
colnames(tu.d) <- colnames(md.n)
tu.d <- as.matrix(tu.d)

TUSC.cca <- vegan::cca(t(TUSC@otu_table)~
                         tu.d) #
anova(TUSC.cca, by ="axis") #
TUSC.cca.fit <- envfit(TUSC.cca, tu.d, perm=9999)
TUSC.cca.fit
#Clinochlorite p = 0.03
#Microcline p = 0.02


##### Plot
##### Export as 5 X 5
# define plot parameters 
par(mar = c(5,5,4,4)+0.1)

#intiate plot
plot(scores(C.Mtf.cca, display = "wa"), xlim = c(-2,3), ylim=c(-2, 8),
     xlab = paste("CCA 1 (", C.Mtf.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2 (", C.Mtf.cca.explainvar2, "%)", sep = ""), 
     pch = 16, cex = 2.0, 
     type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

#plot points
points(scores(C.Mtf.cca, display = "wa"),
       pch = 19, cex = 2, bg = "gray", 
       col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.6))

#add axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

#Add environmental vectors
vectors <- scores(C.Mtf.cca, display = "bp")
arrows(0,0, vectors[c(2,3,7),1]*1.2, vectors[c(2,3,7),2]*1.2,
       lwd=2, lty=1, length = 0.2, col = "red")
text(-1.5,0, pos = 3, 
     labels = c("Muscovite"), col ="red")
text(-0.9,0.5, pos = 3, 
     labels = "Clinochlorite", col ="red")
text(-1.0,-1.4, pos = 3, 
     labels = "Microcline", col ="red")
