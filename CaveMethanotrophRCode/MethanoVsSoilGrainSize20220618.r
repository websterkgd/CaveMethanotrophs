#This r code plots the abundance of the methanotrophs against soil grainsize

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
#rho = -0.13, p = 0.42
cor.test(smd.Sand,ss.C.Mtf.f, method = "spearman")
#rho = -0.09, p = 0.56
cor.test(smd.Silt,ss.C.Mtf.f, method = "spearman")
#rho = 0.17, p = 0.30
cor.test(smd.Clay,ss.C.Mtf.f, method = "spearman")
#rho = 0.00, p = 0.99

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
#rho = -0.19, p = 0.26
cor.test(smd.Sand.cl,ss.C.Mtf.cl, method = "spearman")
#rho = -0.05, p = 0.77
cor.test(smd.Silt.cl,ss.C.Mtf.cl, method = "spearman")
#rho = 0.17, p = 0.30
cor.test(smd.Clay.cl,ss.C.Mtf.cl, method = "spearman")
#rho = -0.20, p = 0.23
 
# Low Affinity
ss.LA.M <-sample_sums(LA.M)
ss.LA.M.f <- c(ss.LA.M[1:31],ss.LA.M[33:43])

smd.Grav <- c(sample_meta_data$X.Gravel[1:31], 
              sample_meta_data$X.Gravel[33:43])

smd.Sand <- c(sample_meta_data$X.Sand[1:31], 
              sample_meta_data$X.Sand[33:43])

smd.Silt <- c(sample_meta_data$X.Silt[1:31], 
              sample_meta_data$X.Silt[33:43])

smd.Clay <- c(sample_meta_data$X.Clay[1:31], 
              sample_meta_data$X.Clay[33:43])

plot(ss.LA.M.f ~ smd.Grav)
plot(ss.LA.M.f ~ smd.Sand)
plot(ss.LA.M.f ~ smd.Silt)
plot(ss.LA.M.f ~ smd.Clay)

cor.test(smd.Grav,ss.LA.M.f, method = "spearman")
#rho = 0.00, p = 0.95
cor.test(smd.Sand,ss.LA.M.f, method = "spearman")
#rho = -0.05, p = 0.74
cor.test(smd.Silt,ss.LA.M.f, method = "spearman")
#rho = -0.05, p = 0.75
cor.test(smd.Clay,ss.LA.M.f, method = "spearman")
#rho = 0.19, p = 0.23

#Filtering CVL
ss.LA.M.cl <- c(ss.LA.M[1:25],ss.LA.M[29:31],ss.LA.M[33:43])

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

plot(ss.LA.M.cl ~ smd.Grav.cl)
plot(ss.LA.M.cl ~ smd.Sand.cl)
plot(ss.LA.M.cl ~ smd.Silt.cl)
plot(ss.LA.M.cl ~ smd.Clay.cl)

cor.test(smd.Grav.cl,ss.LA.M.cl, method = "spearman")
#rho = -0.08, p = 0.62
cor.test(smd.Sand.cl,ss.LA.M.cl, method = "spearman")
#rho = -0.02, p = 0.90
cor.test(smd.Silt.cl,ss.LA.M.cl, method = "spearman")
#rho = 0.08, p = 0.66
cor.test(smd.Clay.cl,ss.LA.M.cl, method = "spearman")
#rho = 0.13, p = 0.43

# High Affinity
ss.TUSC <-sample_sums(TUSC)
ss.TUSC.f <- c(ss.TUSC[1:31],ss.TUSC[33:43])

smd.Grav <- c(sample_meta_data$X.Gravel[1:31], 
              sample_meta_data$X.Gravel[33:43])

smd.Sand <- c(sample_meta_data$X.Sand[1:31], 
              sample_meta_data$X.Sand[33:43])

smd.Silt <- c(sample_meta_data$X.Silt[1:31], 
              sample_meta_data$X.Silt[33:43])

smd.Clay <- c(sample_meta_data$X.Clay[1:31], 
              sample_meta_data$X.Clay[33:43])

plot(ss.TUSC.f ~ smd.Grav)
plot(ss.TUSC.f ~ smd.Sand)
plot(ss.TUSC.f ~ smd.Silt)
plot(ss.TUSC.f ~ smd.Clay)

cor.test(smd.Grav,ss.TUSC.f, method = "spearman")
#rho = -0.14, p = 0.38
cor.test(smd.Sand,ss.TUSC.f, method = "spearman")
#rho = -0.09, p = 0.57
cor.test(smd.Silt,ss.TUSC.f, method = "spearman")
#rho = 0.16, p = 0.31
cor.test(smd.Clay,ss.TUSC.f, method = "spearman")
#rho = 0.04, p = 0.98

#Filtering CVL
ss.TUSC.cl <- c(ss.TUSC[1:25],ss.TUSC[29:31],ss.TUSC[33:43])

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

plot(ss.TUSC.cl ~ smd.Grav.cl)
plot(ss.TUSC.cl ~ smd.Sand.cl)
plot(ss.TUSC.cl ~ smd.Silt.cl)
plot(ss.TUSC.cl ~ smd.Clay.cl)

cor.test(smd.Grav.cl,ss.TUSC.cl, method = "spearman")
#rho = -0.21, p = 0.21
cor.test(smd.Sand.cl,ss.TUSC.cl, method = "spearman")
#rho = -0.04, p = 0.79
cor.test(smd.Silt.cl,ss.TUSC.cl, method = "spearman")
#rho = 0.17, p = 0.32
cor.test(smd.Clay.cl,ss.TUSC.cl, method = "spearman")
#rho = -0.19, p = 0.25

