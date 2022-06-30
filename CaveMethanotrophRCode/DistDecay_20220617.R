#Examining distance decay relationships among
#the Distance decay analysis of methanotrophs

rm(list=ls())
getwd()

#loading data into phyloseq
require(phyloseq)
require(biomformat)
require(vegan)
require(ggplot2)
require(ecodist)

#importing the shared file 
OTU <- import_mothur(mothur_shared_file = 'CvMcrb03_20190225.shared', parseFunction = parse_taxonomy_default)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file = 'CvMcrb03_20190225.taxonomy', parseFunction = parse_taxonomy_default)

#create a physeq table for analysis
physeq = phyloseq(OTU, TAX)

#Creating transformations of the data 
f.physeq = transform_sample_counts(physeq, function(x) x / sum(x)) # fractional abundance transformation of the data

# import meta data 
md = read.table('sample_meta_data+d.txt', header = TRUE, sep = "", 
                dec = ".", row.names = 1)

#restricting meta data to CH4, CO2, grain size, lat and long
md.n <- md[,c(2,3,7,8,13:16)]

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

#removing samples with no methanotrophs and Mock-Com and no grain size
m.otu@.Data <-m.otu@.Data[,-c(21,30,32,36)]

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
T.LA <- merge_phyloseq(T.II,Tp.I)

###basic descriptors
plot_bar(C.Mtf, "Abundance")
ntaxa(C.Mtf) # 75 species
M.s <- rowSums(C.Mtf@otu_table@.Data)
M.s
M.c.s <- colSums(C.Mtf@otu_table@.Data)
M.c.s

#universal environmental variables
coord.dist <- dist(as.matrix(C.Mtf@sam_data@.Data[[3]],
                             C.Mtf@sam_data@.Data[[4]]))
coord.dist <- coord.dist/max(coord.dist) #normalizing

#log transforming appropriate variables in e
e <- cbind(log10(C.Mtf@sam_data@.Data[[1]]),
           log10(C.Mtf@sam_data@.Data[[2]]),
           log1p(C.Mtf@sam_data@.Data[[5]]),
           log10(C.Mtf@sam_data@.Data[[6]]),
           log10(C.Mtf@sam_data@.Data[[7]]),
           log10(C.Mtf@sam_data@.Data[[8]]))

require(stats)
env.dist <- vegdist(e,"euclidean")
env.dist <- 1 - env.dist/(max(env.dist))

#distancedecay total methanotrophs
M.struc.dist <- 1 - phyloseq::distance(C.Mtf, "bray") #bray-curtis similarity

#transforming into list
M.e.f <- cbind(coord.dist, env.dist, M.struc.dist)
M.e.fr <- M.e.f[which(M.struc.dist !=0),]

#regression analysis
M.mno.dd <- lm(M.e.fr[,3] ~ M.e.fr[,1]) # Community structure vs physical distance
plot(M.e.fr[,3] ~ M.e.fr[,1]) #huge gap
abline(M.mno.dd,col="red")
summary(M.mno.dd) #p < 2*10^(-16) m = -0.35, r2=0.19

#env
M.env.dd <- lm(M.e.fr[,2] ~ M.e.fr[,1]) #environmental structure vs physical distance
plot(M.e.fr[,2] ~ M.e.fr[,1])
abline(M.env.dd,col="red")
summary(M.env.dd) #p < 6.6*10^(-5) m = -0.09, r2=0.02

#com vs env
i.e.M.e <- (1-M.e.fr[,2])
M.e_c.dd <- lm(M.e.fr[,3] ~ i.e.M.e) #community vs environment
plot(M.e.fr[,3] ~ i.e.M.e)
abline(M.e_c.dd,col="red")
summary(M.e_c.dd) #p = 0.008 m = -0.12, r2=0.009

x.tm <- slopetest2(M.e.fr[,3],M.e.fr[,1],M.e.fr[,2],M.e.fr[,1])

#Filtering CVL for total 
#removing CVL
CM.cvl <- prune_samples(colnames(
  C.Mtf@otu_table@.Data[,which(sample_names(C.Mtf)!=
                                 c("CDVL-BS","CDVL-Cas-Peq","CDVL-Mid-Spr"))]),
              C.Mtf)

#universal environmental variables
cvlf.dist <- dist(as.matrix(CM.cvl@sam_data@.Data[[3]],
                             CM.cvl@sam_data@.Data[[4]]))
cvlf.dist <- cvlf.dist/max(cvlf.dist) #normalizing

#log transforming appropriate variables in e
cvlf.e <- cbind(log10(CM.cvl@sam_data@.Data[[1]]),
           log10(CM.cvl@sam_data@.Data[[2]]),
           log1p(CM.cvl@sam_data@.Data[[5]]),
           log10(CM.cvl@sam_data@.Data[[6]]),
           log10(CM.cvl@sam_data@.Data[[7]]),
           log10(CM.cvl@sam_data@.Data[[8]]))

cf.e.dist <- vegdist(cvlf.e,"euclidean")
cf.e.dist <- 1 - cf.e.dist/(max(cf.e.dist))

#distancedecay total methanotrophs cvl filt
CM.f.struc.dist <- 
  1 - phyloseq::distance(CM.cvl, "bray") #bray-curtis similarity

#combining for code below
CM.f.e.f <- cbind(cvlf.dist, cf.e.dist, CM.f.struc.dist)
CM.f.e.fr <- subset(CM.f.e.f, CM.f.e.f[,3] !=0)

CM.f.mno.dd <- lm(CM.f.e.fr[,3] ~ CM.f.e.fr[,1])
plot(CM.f.e.fr[,3] ~ CM.f.e.fr[,1])
abline(CM.f.mno.dd,col="red")
summary(CM.f.mno.dd) #p < 0.013 m = -0.09, r2=0.01

#env
CM.f.env.dd <- lm(CM.f.e.fr[,2] ~ CM.f.e.fr[,1])
plot(CM.f.e.fr[,2] ~ CM.f.e.fr[,1])
abline(CM.f.env.dd,col="red")
summary(CM.f.env.dd) #p = 0.12 m = -0.05, r2=0.004

z.tm.cvl <- slopetest2(CM.f.e.fr[,3],CM.f.e.fr[,1],
                       CM.f.e.fr[,2],CM.f.e.fr[,1])
#z = -1.27 p < 0.05; x.tm = -39
# methanotrophs change faster than the environment

#distancedecay USC methanotrophs cvl filt
#removing CVL
TH.cvl <- prune_samples(colnames(
  TUSC@otu_table@.Data[,which(sample_names(TUSC)!=
                                c("CDVL-BS","CDVL-Cas-Peq","CDVL-Mid-Spr"))]),
  TUSC)

# USC environmental variables
TH.cvlf.dist <- dist(as.matrix(TH.cvl@sam_data@.Data[[3]],
                               TH.cvl@sam_data@.Data[[4]]))
TH.cvlf.dist <- TH.cvlf.dist/max(TH.cvlf.dist) #normalizing

#log transforming appropriate variables in e
TH.cvlf.e <- cbind(log10(TH.cvl@sam_data@.Data[[1]]),
                   log10(TH.cvl@sam_data@.Data[[2]]),
                   log1p(TH.cvl@sam_data@.Data[[5]]),
                   log10(TH.cvl@sam_data@.Data[[6]]),
                   log10(TH.cvl@sam_data@.Data[[7]]),
                   log10(TH.cvl@sam_data@.Data[[8]]))

TH.cf.e.dist <- vegdist(TH.cvlf.e,"euclidean")
TH.cf.e.dist <- 1 - TH.cf.e.dist/(max(TH.cf.e.dist))

TH.f.struc.dist <-
  1 - phyloseq::distance(TH.cvl, "bray") #bray-curtis similarity

#transforming into matrix
TH.f.e.f <- cbind(TH.cvlf.dist,TH.cf.e.dist, TH.f.struc.dist)
TH.f.e.fr <- subset(TH.f.e.f, TH.f.e.f[,3] !=0)

TH.f.mno.dd <- lm(TH.f.e.fr[,3] ~ TH.f.e.fr[,1])
plot(TH.f.e.fr[,3] ~ TH.f.e.fr[,1])
abline(TH.f.mno.dd,col="red")
summary(TH.f.mno.dd) #p < 0.01 m = -0.09, r2=0.01

#env
TH.f.env.dd <- lm(TH.f.e.fr[,2] ~ TH.f.e.fr[,1])
plot(TH.f.e.fr[,2] ~ TH.f.e.fr[,1])
abline(TH.f.env.dd,col="red")
summary(TH.f.env.dd) #p < 0.13 m = -0.05, r2= 0.004

z.TH.cvl <- slopetest2(TH.f.e.fr[,3],TH.f.e.fr[,1],
                       TH.f.e.fr[,2],TH.f.e.fr[,1])
#z = -43 p = 0.01
# methanotrophs change faster than env

### Low affinity combination
#distancedecay LA methanotrophs cvl filt
#removing CVL
LA.cvl <- prune_samples(colnames(
  T.LA@otu_table@.Data[,which(sample_names(T.LA)!=
                                c("CDVL-BS","CDVL-Cas-Peq","CDVL-Mid-Spr"))]),
  T.LA)

# USC environmental variables
LA.cvlf.dist <- dist(as.matrix(LA.cvl@sam_data@.Data[[3]],
                               LA.cvl@sam_data@.Data[[4]]))
LA.cvlf.dist <- LA.cvlf.dist/max(LA.cvlf.dist) #normalizing

#log transforming appropriate variables in e
LA.cvlf.e <- cbind(log10(LA.cvl@sam_data@.Data[[1]]),
                   log10(LA.cvl@sam_data@.Data[[2]]),
                   log1p(LA.cvl@sam_data@.Data[[5]]),
                   log10(LA.cvl@sam_data@.Data[[6]]),
                   log10(LA.cvl@sam_data@.Data[[7]]),
                   log10(LA.cvl@sam_data@.Data[[8]]))

LA.cf.e.dist <- vegdist(LA.cvlf.e,"euclidean")
LA.cf.e.dist <- 1 - LA.cf.e.dist/(max(LA.cf.e.dist))

LA.f.struc.dist <-
  1 - phyloseq::distance(LA.cvl, "bray") #bray-curtis similarity

#transforming into matrix
LA.f.e.f <- cbind(LA.cvlf.dist,LA.cf.e.dist, 
                   LA.f.struc.dist)
LA.f.e.fr <- subset(LA.f.e.f, LA.f.e.f[,3] !=0)

LA.f.mno.dd <- lm(LA.f.e.fr[,3] ~ LA.f.e.fr[,1])
plot(LA.f.e.fr[,3] ~ LA.f.e.fr[,1])
abline(LA.f.mno.dd,col="red")
summary(LA.f.mno.dd) #p = 0.09 m = 0.10, r2=0.009

#env
LA.f.env.dd <- lm(LA.f.e.fr[,2] ~ LA.f.e.fr[,1])
plot(LA.f.e.fr[,2] ~ LA.f.e.fr[,1])
abline(LA.f.env.dd,col="red")
summary(LA.f.env.dd) #p = 0.42 m = -0.04, r2= 0.002

z.LA.cvl <- slopetest2(LA.f.e.fr[,3],LA.f.e.fr[,1],
                       LA.f.e.fr[,2],LA.f.e.fr[,1])
#z = 4 p < 0.05 x = 47
# methanotrophs change slower than env


#######creating 1 x 2 plot
par(mfrow=c(1,2), mai = c(0.7, 0.8, 0.3, 0.5),
    mgp=c(2.2,1,0)) # export as 4 X 8

#USC geog dist
plot(TH.f.e.fr[,3] ~ TH.f.e.fr[,1],
     pch = 16, 
     col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3), 
     cex = 1.5,las=1,
     xlab="Relative Geographic Distance", 
     ylab ="High Affinty Community Similarity")
abline(TH.f.mno.dd,col="red")

#USC env vs dist
plot(TH.f.e.fr[,2] ~ TH.f.e.fr[,1],
     pch = 16, 
     col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3), 
     cex = 1.5,las=1,
     xlab="Relative Geographic Distance", 
     ylab ="Relative Environmental Similarity")
abline(TH.f.env.dd,col="red")




