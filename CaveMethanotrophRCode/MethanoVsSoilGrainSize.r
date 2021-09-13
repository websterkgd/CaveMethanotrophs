#This r code plots the abundance of the methylocystaceae and methylococcales
#against soil grainsize

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
#install.packages("vegan")
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
sample_sums(TUSC) 

## Numeric plots
ss.C.Mtf <-sample_sums(C.Mtf)
ss.C.Mtf.f <- c(ss.C.Mtf[1:31],ss.C.Mtf[33:43])

smd.Grav <- c(sample_meta_data$X.Gravel[1:31], 
              sample_meta_data$X.Gravel[33:43])

smd.Sand <- c(sample_meta_data$X.Sand[1:31], 
              sample_meta_data$X.Sand[33:43])

smd.Silt <- c(sample_meta_data$X.Silt[1:31], 
              sample_meta_data$X.Silt[33:43])

smd.Clay <- c(sample_meta_data$X.Clay[1:31], 
              sample_meta_data$X.Clay[33:43])

plot(ss.C.Mtf.f ~ smd.Grav)
plot(ss.C.Mtf.f ~ smd.Sand)
plot(ss.C.Mtf.f ~ smd.Silt)
plot(ss.C.Mtf.f ~ smd.Clay)

cor.test(smd.Grav,ss.C.Mtf.f, method = "spearman")
#rho = -0.14, p = 0.39
cor.test(smd.Sand,ss.C.Mtf.f, method = "spearman")
#rho = -0.11, p = 0.48
cor.test(smd.Silt,ss.C.Mtf.f, method = "spearman")
#rho = -0.17, p = 0.29
cor.test(smd.Clay,ss.C.Mtf.f, method = "spearman")
#rho = -0.10, p = 0.55

#Filtering CVL
ss.C.Mtf.cl <- c(ss.C.Mtf[1:25],ss.C.Mtf[29:31],ss.C.Mtf[33:43]) # CVL filtered

smd.Grav.cl <- c(sample_meta_data$X.Grav[1:25],
               sample_meta_data$X.Grav[29:31],
               sample_meta_data$X.Grav[33:43]) # CVL filtered

smd.Sand.cl <- c(sample_meta_data$X.Sand[1:25],
                 sample_meta_data$X.Sand[29:31],
                 sample_meta_data$X.Sand[33:43]) # CVL filtered

smd.Silt.cl <- c(sample_meta_data$X.Silt[1:25],
                 sample_meta_data$X.Silt[29:31],
                 sample_meta_data$X.Silt[33:43]) # CVL filtered

smd.Clay.cl <- c(sample_meta_data$X.Clay[1:25],
                 sample_meta_data$X.Clay[29:31],
                 sample_meta_data$X.Clay[33:43]) # CVL filtered

plot(ss.C.Mtf.cl ~ smd.Grav.cl)
plot(ss.C.Mtf.cl ~ smd.Sand.cl)
plot(ss.C.Mtf.cl ~ smd.Silt.cl)
plot(ss.C.Mtf.cl ~ smd.Clay.cl)

cor.test(smd.Grav.cl,ss.C.Mtf.cl, method = "spearman")
#rho = -0.2, p = 0.2
cor.test(smd.Sand.cl,ss.C.Mtf.cl, method = "spearman")
#rho = -0.07, p = 0.7
cor.test(smd.Silt.cl,ss.C.Mtf.cl, method = "spearman")
#rho = -0.17, p = 0.3
cor.test(smd.Clay.cl,ss.C.Mtf.cl, method = "spearman")
#rho = -0.08, p = 0.6
 
#Type I
ss.Tp.I <-sample_sums(Tp.I)
ss.Tp.I.f <- c(ss.Tp.I[1:31],ss.Tp.I[33:43])

smd.Grav <- c(sample_meta_data$X.Gravel[1:31], 
              sample_meta_data$X.Gravel[33:43])

smd.Sand <- c(sample_meta_data$X.Sand[1:31], 
              sample_meta_data$X.Sand[33:43])

smd.Silt <- c(sample_meta_data$X.Silt[1:31], 
              sample_meta_data$X.Silt[33:43])

smd.Clay <- c(sample_meta_data$X.Clay[1:31], 
              sample_meta_data$X.Clay[33:43])

plot(ss.Tp.I.f ~ smd.Grav)
plot(ss.Tp.I.f ~ smd.Sand)
plot(ss.Tp.I.f ~ smd.Silt)
plot(ss.Tp.I.f ~ smd.Clay)

cor.test(smd.Grav,ss.Tp.I.f, method = "spearman")
#rho = -0.07, p = 0.65
cor.test(smd.Sand,ss.Tp.I.f, method = "spearman")
#rho = -0.13, p = 0.39
cor.test(smd.Silt,ss.Tp.I.f, method = "spearman")
#rho = -0.03, p = 0.83
cor.test(smd.Clay,ss.Tp.I.f, method = "spearman")
#rho = 0.10, p = 0.53

#Filtering CVL
ss.Tp.I.cl <- c(ss.Tp.I[1:25],ss.Tp.I[29:31],ss.Tp.I[33:43])

smd.Grav.cl <- c(sample_meta_data$X.Grav[1:25],
                 sample_meta_data$X.Grav[29:31],
                 sample_meta_data$X.Grav[33:43]) # CVL filtered

smd.Sand.cl <- c(sample_meta_data$X.Sand[1:25],
                 sample_meta_data$X.Sand[29:31],
                 sample_meta_data$X.Sand[33:43]) # CVL filtered

smd.Silt.cl <- c(sample_meta_data$X.Silt[1:25],
                 sample_meta_data$X.Silt[29:31],
                 sample_meta_data$X.Silt[33:43]) # CVL filtered

smd.Clay.cl <- c(sample_meta_data$X.Clay[1:25],
                 sample_meta_data$X.Clay[29:31],
                 sample_meta_data$X.Clay[33:43]) # CVL filtered

plot(ss.Tp.I.cl ~ smd.Grav.cl)
plot(ss.Tp.I.cl ~ smd.Sand.cl)
plot(ss.Tp.I.cl ~ smd.Silt.cl)
plot(ss.Tp.I.cl ~ smd.Clay.cl)

cor.test(smd.Grav.cl,ss.Tp.I.cl, method = "spearman")
#rho = -0.03, p = 0.8
cor.test(smd.Sand.cl,ss.Tp.I.cl, method = "spearman")
#rho = -0.13, p = 0.4
cor.test(smd.Silt.cl,ss.Tp.I.cl, method = "spearman")
#rho = 0.03, p = 0.9
cor.test(smd.Clay.cl,ss.Tp.I.cl, method = "spearman")
#rho = 0.08, p = 0.7

#Type II
ss.T.II <-sample_sums(T.II)
ss.T.II.f <- c(ss.T.II[1:31],ss.T.II[33:43])

smd.Grav <- c(sample_meta_data$X.Gravel[1:31], 
              sample_meta_data$X.Gravel[33:43])

smd.Sand <- c(sample_meta_data$X.Sand[1:31], 
              sample_meta_data$X.Sand[33:43])

smd.Silt <- c(sample_meta_data$X.Silt[1:31], 
              sample_meta_data$X.Silt[33:43])

smd.Clay <- c(sample_meta_data$X.Clay[1:31], 
              sample_meta_data$X.Clay[33:43])

plot(ss.T.II.f ~ smd.Grav)
plot(ss.T.II.f ~ smd.Sand)
plot(ss.T.II.f ~ smd.Silt)
plot(ss.T.II.f ~ smd.Clay)

cor.test(smd.Grav,ss.T.II.f, method = "spearman")
#rho = 0.15, p = 0.34
cor.test(smd.Sand,ss.T.II.f, method = "spearman")
#rho = 0.15, p = 0.33
cor.test(smd.Silt,ss.T.II.f, method = "spearman")
#rho = -0.21, p = 0.19
cor.test(smd.Clay,ss.T.II.f, method = "spearman")
#rho = 0.10, p = 0.51

#Filtering CVL
ss.T.II.cl <- c(ss.T.II[1:25],ss.T.II[29:31],ss.T.II[33:43]) 

smd.Grav.cl <- c(sample_meta_data$X.Grav[1:25],
                 sample_meta_data$X.Grav[29:31],
                 sample_meta_data$X.Grav[33:43]) # CVL filtered

smd.Sand.cl <- c(sample_meta_data$X.Sand[1:25],
                 sample_meta_data$X.Sand[29:31],
                 sample_meta_data$X.Sand[33:43]) # CVL filtered

smd.Silt.cl <- c(sample_meta_data$X.Silt[1:25],
                 sample_meta_data$X.Silt[29:31],
                 sample_meta_data$X.Silt[33:43]) # CVL filtered

smd.Clay.cl <- c(sample_meta_data$X.Clay[1:25],
                 sample_meta_data$X.Clay[29:31],
                 sample_meta_data$X.Clay[33:43]) # CVL filtered

plot(ss.T.II.cl ~ smd.Grav.cl)
plot(ss.T.II.cl ~ smd.Sand.cl)
plot(ss.T.II.cl ~ smd.Silt.cl)
plot(ss.T.II.cl ~ smd.Clay.cl)

cor.test(smd.Grav.cl,ss.T.II.cl, method = "spearman")
#rho = 0.09, p = 0.6
cor.test(smd.Sand.cl,ss.T.II.cl, method = "spearman")
#rho = 0.19, p = 0.2
cor.test(smd.Silt.cl,ss.T.II.cl, method = "spearman")
#rho = -0.24, p = 0.14
cor.test(smd.Clay.cl,ss.T.II.cl, method = "spearman")
#rho = 0.08, p = 0.62
