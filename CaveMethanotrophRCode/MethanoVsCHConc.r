#This r code plots the abundance of methanotrophs against methane concentration

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

#calling vegan
require('vegan')

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

# import meta data
md = read.table('sample_meta_data+d.txt', header = TRUE, sep = "", dec = ".", 
                row.names = 1) 
  
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

## Numeric plots
ss.C.Mtf <-sample_sums(C.Mtf)
ss.C.Mtf.f <- c(ss.C.Mtf[1:31],ss.C.Mtf[33:43])

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:31], 
               sample_meta_data$CH4_conc.ppm.[33:43])

ss.C.Mtf.f <- c(ss.C.Mtf[1:25],ss.C.Mtf[29:31],ss.C.Mtf[33:43]) # CVL filtered

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:25], 
               sample_meta_data$CH4_conc.ppm.[29:31],
               sample_meta_data$CH4_conc.ppm.[33:43]) # CVL filtered

plot(smd.CH4.F,ss.C.Mtf.f)

cor.test(smd.CH4.F, ss.C.Mtf.f, method = "spearman") # I can filter CVL here
#p = 0.01, rho = 0.40 #CVL filtered, S = 5920

summary(lm(smd.CH4.F~ss.C.Mtf.f))

#Now bring in varying groups of methanotrophs
ss.TUSC <-sample_sums(TUSC)
ss.TUSC.f <- c(ss.TUSC[1:31],ss.TUSC[33:43])

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:31], 
               sample_meta_data$CH4_conc.ppm.[33:43])

ss.TUSC.f <- c(ss.TUSC[1:25],ss.TUSC[29:31],ss.TUSC[33:43]) # CVL filtered

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:25], 
               sample_meta_data$CH4_conc.ppm.[29:31],
               sample_meta_data$CH4_conc.ppm.[33:43]) # CVL filtered

plot(smd.CH4.F,ss.TUSC.f)

cor.test(smd.CH4.F, ss.TUSC.f, method = "spearman") # I can filter CVL here
#p = 0.01 rho = 0.40 #S 5887 #CVL filtered #basically same pattern as whole data set
 
ss.T.II <-sample_sums(T.II)
ss.T.II.f <- c(ss.T.II[1:31],ss.T.II[33:43])

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:31], 
               sample_meta_data$CH4_conc.ppm.[33:43])

ss.T.II.f <- c(ss.T.II[1:25],ss.T.II[29:31],ss.T.II[33:43]) # CVL filtered

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:25], 
               sample_meta_data$CH4_conc.ppm.[29:31],
               sample_meta_data$CH4_conc.ppm.[33:43]) # CVL filtered

plot(smd.CH4.F,ss.T.II.f)

cor.test(smd.CH4.F, ss.T.II.f, method = "spearman") # I can filter CVL here
#p = 0.79 rho = 0.04 #CVL filtered s = 9462

ss.Tp.I <-sample_sums(Tp.I)
ss.Tp.I.f <- c(ss.Tp.I[1:31],ss.Tp.I[33:43])

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:31], 
               sample_meta_data$CH4_conc.ppm.[33:43])

ss.Tp.I.f <- c(ss.Tp.I[1:25],ss.Tp.I[29:31],ss.Tp.I[33:43]) # CVL filtered

smd.CH4.F <- c(sample_meta_data$CH4_conc.ppm.[1:25], 
               sample_meta_data$CH4_conc.ppm.[29:31],
               sample_meta_data$CH4_conc.ppm.[33:43]) # CVL filtered

plot(smd.CH4.F,ss.Tp.I.f)

cor.test(smd.CH4.F, ss.Tp.I.f, method = "spearman") # I can filter CVL here
#p = 0.94 rho = -0.01 #CVL filtered s = 10002

plot(smd.CH4.F,(ss.Tp.I.f+ss.T.II.f+ss.TUSC.f))

#####Plot for publication  ####Change to percent i.e. multiply by 100
#####Export as 5 X 5
plot(smd.CH4.F,ss.C.Mtf.f,
     xlab = paste("Methane Concentration (ppmv)"),
     ylab = paste("Mehanotroph Relative Abundance (%)"),
     pch = 16, 
     col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.6), cex = 2,las=1,
     xlim = c(0,3.5))
text(2.44,0.04, expression(rho), srt =-10) 
text(2.76,0.04, as.expression("= 0.40")) 
text(2.69,0.037, as.expression(italic(S)~"="~5920)) 

