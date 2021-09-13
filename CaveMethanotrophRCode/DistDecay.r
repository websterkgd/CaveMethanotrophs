#Examining distance decay relationships among
#the taxa

rm(list=ls())
getwd()

#loading data into phyloseq
require(phyloseq)
require(biomformat)
require(vegan)
require(ggplot2)
require(simba)
require(ecodist)

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
md = read.table('sample_meta_data+d.txt', header = TRUE, sep = "", 
                dec = ".", row.names = 1)

#restricting meta data to CH4, CO2, grain size, lat and long
md.n <- md[,c(2,3,7,8,13:16)]

#Let Phyloseq see the meta data
sample_meta_data <- sample_data(md.n, errorIfNULL = TRUE)

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
ntaxa(C.Mtf) # 22 species
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

env.dist <- vegdist(e,"euclidean")
env.dist <- 1 - env.dist/(max(env.dist))
env.dist.ls <- liste(env.dist, entry="env")
coord.dist.ls <- liste(coord.dist, entry="dist")

#distancedecay total methanotrophs
M.struc.dist <- 1 - phyloseq::distance(C.Mtf, "bray") #bray-curtis similarity

#transforming into list
M.struc.dist.ls <- liste(M.struc.dist, entry="struc")

M.e.fr <- cbind(coord.dist.ls,env.dist.ls[,3], M.struc.dist.ls[,3])
names(M.e.fr)[4:5] <- c("env","struc")
attach(M.e.fr)
M.e.fr <- subset(M.e.fr, struc !=0)

M.mno.dd <- lm(M.e.fr$struc ~ M.e.fr$dist)
plot(M.e.fr$struc ~ M.e.fr$dist)
abline(M.mno.dd,col="red")
summary(M.mno.dd) #p < 2*10^(-16) m = -0.36, r2=0.16

#env
M.env.dd <- lm(M.e.fr$env ~ M.e.fr$dist)
plot(M.e.fr$env ~ M.e.fr$dist)
abline(M.env.dd,col="red")
summary(M.env.dd) #p < 2*10^(-16) m = -0.09, r2=0.02

#com vs env
i.e.M.e <- (1-M.e.fr$env)
M.e_c.dd <- lm(M.e.fr$struc ~ i.e.M.e)
plot(M.e.fr$struc ~ i.e.M.e)
abline(M.e_c.dd,col="red")
summary(M.e_c.dd) #p < 0.003 m = -0.15, r2=0.01

d.tm <- diffslope2(M.e.fr$dist,M.e.fr$struc,M.e.fr$dist,M.e.fr$env) 
#dif = -0.27
# methanotrophs change faster than the environment

dI.tm <- diffslope2(M.e.fr$dist,M.e.fr$struc,i.e.M.e,M.e.fr$struc)
#dif = -0.21
# Com v E is steeper closer to species sorting in mass effects

###distancedecay Type I methanotrophs
#Remove samples without type 1s
Tp.I <- prune_samples(names(which(sample_sums(Tp.I)!=0)),Tp.I)

TypI.struc.dist <- 1 - phyloseq::distance(Tp.I, "bray") #bray-curtis similarity

#universal environmental variables
T1.C.dist <- dist(as.matrix(Tp.I@sam_data@.Data[[3]],
                             Tp.I@sam_data@.Data[[4]]))
T1.C.dist <- T1.C.dist/max(T1.C.dist) #normalizing

#log transforming appropriate variables in e
T1.e <- cbind(log10(Tp.I@sam_data@.Data[[1]]),
           log10(Tp.I@sam_data@.Data[[2]]),
           log1p(Tp.I@sam_data@.Data[[5]]),
           log10(Tp.I@sam_data@.Data[[6]]),
           log10(Tp.I@sam_data@.Data[[7]]),
           log10(Tp.I@sam_data@.Data[[8]]))

T1.e.dist <- vegdist(T1.e,"euclidean")
T1.e.dist <- 1 - T1.e.dist/(max(T1.e.dist))
T1.e.dist.ls <- liste(T1.e.dist, entry="env")
T1.C.dist.ls <- liste(T1.C.dist, entry="dist")

#transforming into list
TypI.struc.dist.ls <- liste(TypI.struc.dist, entry="struc")

TypI.e.fr <- cbind(T1.C.dist.ls,T1.e.dist.ls[,3], TypI.struc.dist.ls[,3])
names(TypI.e.fr)[4:5] <- c("env","struc")
attach(TypI.e.fr)
TypI.e.fr <- subset(TypI.e.fr, struc !=0)

TypI.mno.dd <- lm(TypI.e.fr$struc ~ TypI.e.fr$dist)
plot(TypI.e.fr$struc ~ TypI.e.fr$dist)
abline(TypI.mno.dd,col="red")
summary(TypI.mno.dd) #p 0.29 m = 0.09, r2=0.001

TypI.mno.ed <- lm(TypI.e.fr$env ~ TypI.e.fr$dist)
plot(TypI.e.fr$env ~ TypI.e.fr$dist)
abline(TypI.mno.ed,col="red")
summary(TypI.mno.ed) #p < 0.78 m = 0.02, r2=0.001

#com vs env
i.e.TypI <- (1-TypI.e.fr$env)
TypI.e_c.dd <- lm(TypI.e.fr$struc ~ i.e.TypI)
plot(TypI.e.fr$struc ~ i.e.TypI)
abline(TypI.e_c.dd,col="red")
summary(TypI.e_c.dd) #p < 0.27 m = -0.14, r2=0.01

d.TpI<- diffslope2(TypI.e.fr$dist, TypI.e.fr$env, TypI.e.fr$dist, TypI.e.fr$struc)
#dif = -0.07p = 0.19, result is accurate
# No difference

dI.TypI <- diffslope2(TypI.e.fr$dist,TypI.e.fr$struc,i.e.TypI,TypI.e.fr$struc) 
#dif = 0.23 p = 0.006 
# Com v dist is steeper closer to nuetral

###distancedecay Type II methanotrophs
#remove samples without Type IIs
T.II <- prune_samples(names(which(sample_sums(T.II)!=0)),T.II)

TpII.struc.dist <- 1 - phyloseq::distance(T.II, "bray") #bray-curtis similarity

#universal environmental variables
T2.C.dist <- dist(as.matrix(T.II@sam_data@.Data[[3]],
                            T.II@sam_data@.Data[[4]]))
T2.C.dist <- T2.C.dist/max(T2.C.dist) #normalizing

#log transforming appropriate variables in e
T2.e <- cbind(log10(T.II@sam_data@.Data[[1]]),
              log10(T.II@sam_data@.Data[[2]]),
              log1p(T.II@sam_data@.Data[[5]]),
              log10(T.II@sam_data@.Data[[6]]),
              log10(T.II@sam_data@.Data[[7]]),
              log10(T.II@sam_data@.Data[[8]]))

T2.e.dist <- vegdist(T2.e,"euclidean")
T2.e.dist <- 1 - T2.e.dist/(max(T2.e.dist))
T2.e.dist.ls <- liste(T2.e.dist, entry="env")
T2.C.dist.ls <- liste(T2.C.dist, entry="dist")

#transforming into list
TpII.struc.dist.ls <- liste(TpII.struc.dist, entry="struc")

TpII.e.fr <- cbind(T2.C.dist.ls,T2.e.dist.ls[,3], TpII.struc.dist.ls[,3])
names(TpII.e.fr)[4:5] <- c("env","struc")
attach(TpII.e.fr)
TpII.e.fr <- subset(TpII.e.fr, struc !=0)

TpII.mno.dd <- lm(TpII.e.fr$struc ~ TpII.e.fr$dist)
plot(TpII.e.fr$struc ~ TpII.e.fr$dist)
abline(TpII.mno.dd,col="red")
summary(TpII.mno.dd) #p 0.29 m = 0.09, r2=0.005

TpII.mno.ed <- lm(TpII.e.fr$env ~ TpII.e.fr$dist)
plot(TpII.e.fr$env ~ TpII.e.fr$dist)
abline(TpII.mno.ed,col="red")
summary(TpII.mno.ed) #p < 0.15 m = -0.08, r2=0.005

#com vs env
i.e.TpII <- (1-TpII.e.fr$env)
TpII.e_c.dd <- lm(TpII.e.fr$struc ~ i.e.TpII)
plot(TpII.e.fr$struc ~ i.e.TpII)
abline(TpII.e_c.dd,col="red")
summary(TpII.e_c.dd) #p < 0.60 m = -0.05, r2=0.001

d.TpII<- diffslope2(TpII.e.fr$dist, TpII.e.fr$env, TpII.e.fr$dist, TpII.e.fr$struc)
#dif = -0.17 p = 0.004, result is accurate
# methanos change slower than env

dI.TpII <- diffslope2(TpII.e.fr$dist,TpII.e.fr$struc,i.e.TpII,TpII.e.fr$struc) 
#dif = 0.14 p = 0.05 
# Com v dist is steeper closer to nuetral

# ###distancedecay USC methanotrophs
# #Remove samples without type 1s
TUSC <- prune_samples(names(which(sample_sums(TUSC)!=0)),TUSC)

T.HA.struc.dist <- 1 - phyloseq::distance(TUSC, "bray") #bray-curtis similarity

# #universal environmental variables
TH.C.dist <- dist(as.matrix(TUSC@sam_data@.Data[[3]],
                            TUSC@sam_data@.Data[[4]]))
TH.C.dist <- TH.C.dist/max(TH.C.dist) #normalizing

#log transforming appropriate variables in e
TH.e <- cbind(log10(TUSC@sam_data@.Data[[1]]),
              log10(TUSC@sam_data@.Data[[2]]),
              log1p(TUSC@sam_data@.Data[[5]]),
              log10(TUSC@sam_data@.Data[[6]]),
              log10(TUSC@sam_data@.Data[[7]]),
              log10(TUSC@sam_data@.Data[[8]]))

TH.e.dist <- vegdist(TH.e,"euclidean")
TH.e.dist <- 1 - TH.e.dist/(max(TH.e.dist))
TH.e.dist.ls <- liste(TH.e.dist, entry="env")
TH.C.dist.ls <- liste(TH.C.dist, entry="dist")

#transforming into list
T.HA.struc.dist.ls <- liste(T.HA.struc.dist, entry="struc")

T.HA.e.fr <- cbind(TH.C.dist.ls,TH.e.dist.ls[,3], T.HA.struc.dist.ls[,3])
names(T.HA.e.fr)[4:5] <- c("env","struc")
attach(T.HA.e.fr)
T.HA.e.fr <- subset(T.HA.e.fr, struc !=0)

T.HA.mno.dd <- lm(T.HA.e.fr$struc ~ T.HA.e.fr$dist)
plot(T.HA.e.fr$struc ~ T.HA.e.fr$dist)
abline(T.HA.mno.dd,col="red")
summary(T.HA.mno.dd) #p<0.01 m = -0.37, r2=0.17

T.HA.mno.ed <- lm(T.HA.e.fr$env ~ T.HA.e.fr$dist)
plot(T.HA.e.fr$env ~ T.HA.e.fr$dist)
abline(T.HA.mno.ed,col="red")
summary(T.HA.mno.ed) #p < 0.01 m = -0.09, r2=0.02

#com vs env
i.e.T.HA <- (1-T.HA.e.fr$env)
T.HA.e_c.dd <- lm(T.HA.e.fr$struc ~ i.e.T.HA)
plot(T.HA.e.fr$struc ~ i.e.T.HA)
abline(T.HA.e_c.dd,col="red")
summary(T.HA.e_c.dd) #p < 0.01 m = -0.16, r2=0.01

d.T.HA<- diffslope2(T.HA.e.fr$dist, T.HA.e.fr$env, T.HA.e.fr$dist, T.HA.e.fr$struc)
#dif = 0.27 p = 0.001, result is accurate
# Methanos change fast than env

dI.T.HA <- diffslope2(T.HA.e.fr$dist,T.HA.e.fr$struc,i.e.T.HA,T.HA.e.fr$struc)
#dif = 0.21 p = 0.001
# Com v dist is steeper closer to nuetral


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
cf.e.dist.ls <- liste(cf.e.dist, entry="env")
cvlf.dist.ls <- liste(cvlf.dist, entry="dist")

#distancedecay total methanotrophs cvl filt
CM.f.struc.dist <- 
  1 - phyloseq::distance(CM.cvl, "bray") #bray-curtis similarity

#transforming into list
CM.f.struc.dist.ls <- liste(CM.f.struc.dist, entry="struc")

CM.f.e.fr <- cbind(cvlf.dist.ls,cf.e.dist.ls[,3], CM.f.struc.dist.ls[,3])
names(CM.f.e.fr)[4:5] <- c("env","struc")
attach(CM.f.e.fr)
CM.f.e.fr <- subset(CM.f.e.fr, struc !=0)

CM.f.mno.dd <- lm(CM.f.e.fr$struc ~ CM.f.e.fr$dist)
plot(CM.f.e.fr$struc ~ CM.f.e.fr$dist)
abline(CM.f.mno.dd,col="red")
summary(CM.f.mno.dd) #p < 0.008 m = -0.14, r2=0.02

#env
CM.f.env.dd <- lm(CM.f.e.fr$env ~ CM.f.e.fr$dist)
plot(CM.f.e.fr$env ~ CM.f.e.fr$dist)
abline(CM.f.env.dd,col="red")
summary(CM.f.env.dd) #p < 0.13 m = -0.05, r2=0.003

#com vs env
i.e.CM.f.e <- (1-CM.f.e.fr$env)
CM.f.e_c.dd <- lm(CM.f.e.fr$struc ~ i.e.CM.f.e)
plot(CM.f.e.fr$struc ~ i.e.CM.f.e)
abline(CM.f.e_c.dd,col="red")
summary(CM.f.e_c.dd) #p < 0.32 m = -0.05, r2=0.02

d.tm.cvl <- diffslope2(CM.f.e.fr$dist,CM.f.e.fr$struc,
                       CM.f.e.fr$dist,CM.f.e.fr$env)
#dif = -0.09 p = 0.004
# methanotrophs change faster than the environment

dI.tm.cvl <- diffslope2(CM.f.e.fr$dist,CM.f.e.fr$struc,
                        i.e.CM.f.e,CM.f.e.fr$struc)
#dif = -0.08
# Com v E is steeper closer to species sorting in mass effects

#distancedecay type 1 methanotrophs cvl filt
#removing CVL
T1.cvl <- prune_samples(colnames(
  Tp.I@otu_table@.Data[,which(sample_names(Tp.I)!=
                                 c("CDVL-BS","CDVL-Cas-Peq","CDVL-Mid-Spr"))]),
  Tp.I)

#type 1 environmental variables
T1.cvlf.dist <- dist(as.matrix(T1.cvl@sam_data@.Data[[3]],
                            T1.cvl@sam_data@.Data[[4]]))
T1.cvlf.dist <- T1.cvlf.dist/max(T1.cvlf.dist) #normalizing

#log transforming appropriate variables in e
T1.cvlf.e <- cbind(log10(T1.cvl@sam_data@.Data[[1]]),
                log10(T1.cvl@sam_data@.Data[[2]]),
                log1p(T1.cvl@sam_data@.Data[[5]]),
                log10(T1.cvl@sam_data@.Data[[6]]),
                log10(T1.cvl@sam_data@.Data[[7]]),
                log10(T1.cvl@sam_data@.Data[[8]]))

T1.cf.e.dist <- vegdist(T1.cvlf.e,"euclidean")
T1.cf.e.dist <- 1 - T1.cf.e.dist/(max(T1.cf.e.dist))
T1.cf.e.dist.ls <- liste(T1.cf.e.dist, entry="env")
T1.cvlf.dist.ls <- liste(T1.cvlf.dist, entry="dist")

T1.f.struc.dist <-
  1 - phyloseq::distance(T1.cvl, "bray") #bray-curtis similarity

#transforming into list
T1.f.struc.dist.ls <- liste(T1.f.struc.dist, entry="struc")

T1.f.e.fr <- cbind(T1.cvlf.dist.ls,T1.cf.e.dist.ls[,3], T1.f.struc.dist.ls[,3])
names(T1.f.e.fr)[4:5] <- c("env","struc")
attach(T1.f.e.fr)
T1.f.e.fr <- subset(T1.f.e.fr, struc !=0)

T1.f.mno.dd <- lm(T1.f.e.fr$struc ~ T1.f.e.fr$dist)
plot(T1.f.e.fr$struc ~ T1.f.e.fr$dist)
abline(T1.f.mno.dd,col="red")
summary(T1.f.mno.dd) #p < 0.29 m = 0.10, r2=0.001

#env
T1.f.env.dd <- lm(T1.f.e.fr$env ~ T1.f.e.fr$dist)
plot(T1.f.e.fr$env ~ T1.f.e.fr$dist)
abline(T1.f.env.dd,col="red")
summary(T1.f.env.dd) #p < 0.78 m = 0.02, r2= -0.01

#com vs env
i.e.T1.f.e <- (1-T1.f.e.fr$env)
T1.f.e_c.dd <- lm(T1.f.e.fr$struc ~ i.e.T1.f.e)
plot(T1.f.e.fr$struc ~ i.e.T1.f.e)
abline(T1.f.e_c.dd,col="red")
summary(T1.f.e_c.dd) #p < 0.27 m = -0.14, r2=0.003

d.T1.cvl <- diffslope2(T1.f.e.fr$dist,T1.f.e.fr$struc,
                       T1.f.e.fr$dist,T1.f.e.fr$env)
#dif = 0.08 p = 0.182
# env changes faster than methanotrophs no dif

dI.T1.cvl <- diffslope2(T1.f.e.fr$dist,T1.f.e.fr$struc,
                        i.e.T1.f.e,T1.f.e.fr$struc)
#dif = 0.24 p = 0.006
# Com v E is steeper closer to neutral

#distancedecay type 2 methanotrophs cvl filt
#removing CVL
T2.cvl <- prune_samples(colnames(
  T.II@otu_table@.Data[,which(sample_names(T.II)!=
                                c("CDVL-BS","CDVL-Cas-Peq","CDVL-Mid-Spr"))]),
  T.II)

#type 2 environmental variables
T2.cvlf.dist <- dist(as.matrix(T2.cvl@sam_data@.Data[[3]],
                               T2.cvl@sam_data@.Data[[4]]))
T2.cvlf.dist <- T2.cvlf.dist/max(T2.cvlf.dist) #normalizing

#log transforming appropriate variables in e
T2.cvlf.e <- cbind(log10(T2.cvl@sam_data@.Data[[1]]),
                   log10(T2.cvl@sam_data@.Data[[2]]),
                   log1p(T2.cvl@sam_data@.Data[[5]]),
                   log10(T2.cvl@sam_data@.Data[[6]]),
                   log10(T2.cvl@sam_data@.Data[[7]]),
                   log10(T2.cvl@sam_data@.Data[[8]]))

T2.cf.e.dist <- vegdist(T2.cvlf.e,"euclidean")
T2.cf.e.dist <- 1 - T2.cf.e.dist/(max(T2.cf.e.dist))
T2.cf.e.dist.ls <- liste(T2.cf.e.dist, entry="env")
T2.cvlf.dist.ls <- liste(T2.cvlf.dist, entry="dist")

T2.f.struc.dist <-
  1 - phyloseq::distance(T2.cvl, "bray") #bray-curtis similarity

#transforming into list
T2.f.struc.dist.ls <- liste(T2.f.struc.dist, entry="struc")

T2.f.e.fr <- cbind(T2.cvlf.dist.ls,T2.cf.e.dist.ls[,3], T2.f.struc.dist.ls[,3])
names(T2.f.e.fr)[4:5] <- c("env","struc")
attach(T2.f.e.fr)
T2.f.e.fr <- subset(T2.f.e.fr, struc !=0)

T2.f.mno.dd <- lm(T2.f.e.fr$struc ~ T2.f.e.fr$dist)
plot(T2.f.e.fr$struc ~ T2.f.e.fr$dist)
abline(T2.f.mno.dd,col="red")
summary(T2.f.mno.dd) #p < 0.29 m = 0.09, r2=0.0005

#env
T2.f.env.dd <- lm(T2.f.e.fr$env ~ T2.f.e.fr$dist)
plot(T2.f.e.fr$env ~ T2.f.e.fr$dist)
abline(T2.f.env.dd,col="red")
summary(T2.f.env.dd) #p < 0.15 m = -0.08, r2= 0.005

#com vs env
i.e.T2.f.e <- (1-T2.f.e.fr$env)
T2.f.e_c.dd <- lm(T2.f.e.fr$struc ~ i.e.T2.f.e)
plot(T2.f.e.fr$struc ~ i.e.T2.f.e)
abline(T2.f.e_c.dd,col="red")
summary(T2.f.e_c.dd) #p < 0.60 m = -0.05, r2=-0.003

d.T2.cvl <- diffslope2(T2.f.e.fr$dist,T2.f.e.fr$struc,
                       T2.f.e.fr$dist,T2.f.e.fr$env)
#dif = 0.17 p = 0.009
# env changes faster than methanotrophs 

dI.T2.cvl <- diffslope2(T2.f.e.fr$dist,T2.f.e.fr$struc,
                        i.e.T2.f.e,T2.f.e.fr$struc)
#dif = 0.14 p = 0.03
# closer to neutral

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
TH.cf.e.dist.ls <- liste(TH.cf.e.dist, entry="env")
TH.cvlf.dist.ls <- liste(TH.cvlf.dist, entry="dist")

TH.f.struc.dist <-
  1 - phyloseq::distance(TH.cvl, "bray") #bray-curtis similarity

#transforming into list
TH.f.struc.dist.ls <- liste(TH.f.struc.dist, entry="struc")

TH.f.e.fr <- cbind(TH.cvlf.dist.ls,TH.cf.e.dist.ls[,3], TH.f.struc.dist.ls[,3])
names(TH.f.e.fr)[4:5] <- c("env","struc")
attach(TH.f.e.fr)
TH.f.e.fr <- subset(TH.f.e.fr, struc !=0)

TH.f.mno.dd <- lm(TH.f.e.fr$struc ~ TH.f.e.fr$dist)
plot(TH.f.e.fr$struc ~ TH.f.e.fr$dist)
abline(TH.f.mno.dd,col="red")
summary(TH.f.mno.dd) #p < 0.01 m = -0.14, r2=0.02

#env
TH.f.env.dd <- lm(TH.f.e.fr$env ~ TH.f.e.fr$dist)
plot(TH.f.e.fr$env ~ TH.f.e.fr$dist)
abline(TH.f.env.dd,col="red")
summary(TH.f.env.dd) #p < 0.14 m = -0.05, r2= 0.002

#com vs env
i.e.TH.f.e <- (1-TH.f.e.fr$env)
TH.f.e_c.dd <- lm(TH.f.e.fr$struc ~ i.e.TH.f.e)
plot(TH.f.e.fr$struc ~ i.e.TH.f.e)
abline(TH.f.e_c.dd,col="red")
summary(TH.f.e_c.dd) #p < 0.32 m = -0.05, r2= 0

d.TH.cvl <- diffslope2(TH.f.e.fr$dist,TH.f.e.fr$struc,
                       TH.f.e.fr$dist,TH.f.e.fr$env)
#dif = -0.09 p = 0.01
# methanotrophs change faster than env

dI.TH.cvl <- diffslope2(TH.f.e.fr$dist,TH.f.e.fr$struc,
                        i.e.TH.f.e,TH.f.e.fr$struc)
#dif = -0.08 p = 0.03
# Com v dist is steeper 

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
LA.cf.e.dist.ls <- liste(LA.cf.e.dist, entry="env")
LA.cvlf.dist.ls <- liste(LA.cvlf.dist, entry="dist")

LA.f.struc.dist <-
  1 - phyloseq::distance(LA.cvl, "bray") #bray-curtis similarity

#transforming into list
LA.f.struc.dist.ls <- liste(LA.f.struc.dist, entry="struc")

LA.f.e.fr <- cbind(LA.cvlf.dist.ls,LA.cf.e.dist.ls[,3], 
                   LA.f.struc.dist.ls[,3])
names(LA.f.e.fr)[4:5] <- c("env","struc")
attach(LA.f.e.fr)
LA.f.e.fr <- subset(LA.f.e.fr, struc !=0)

LA.f.mno.dd <- lm(LA.f.e.fr$struc ~ LA.f.e.fr$dist)
plot(LA.f.e.fr$struc ~ LA.f.e.fr$dist)
abline(LA.f.mno.dd,col="red")
summary(LA.f.mno.dd) #p = 0.14 m = 0.10, r2=0.004

#env
LA.f.env.dd <- lm(LA.f.e.fr$env ~ LA.f.e.fr$dist)
plot(LA.f.e.fr$env ~ LA.f.e.fr$dist)
abline(LA.f.env.dd,col="red")
summary(LA.f.env.dd) #p = 0.18 m = -0.07, r2= 0.003

#com vs env
i.e.LA.f.e <- (1-LA.f.e.fr$env)
LA.f.e_c.dd <- lm(LA.f.e.fr$struc ~ i.e.LA.f.e)
plot(LA.f.e.fr$struc ~ i.e.LA.f.e)
abline(LA.f.e_c.dd,col="red")
summary(LA.f.e_c.dd) #p = 0.45 m = -0.06, r2= -0.002

d.LA.cvl <- diffslope2(LA.f.e.fr$dist,LA.f.e.fr$struc,
                       LA.f.e.fr$dist,LA.f.e.fr$env)
#dif = 0.17 p = 0.002
# methanotrophs change slower than env

dI.LA.cvl <- diffslope2(LA.f.e.fr$dist,LA.f.e.fr$struc,
                        i.e.LA.f.e,LA.f.e.fr$struc)

#dif = 0.16 p = 0.004
# methanotrophs change slower than env

#######creating 1 x 2 plot
par(mfrow=c(1,2), mai = c(0.7, 0.8, 0.3, 0.5),
    mgp=c(2.2,1,0)) # export as 8 X 4

#USC geog dist
plot(TH.f.e.fr$struc ~ TH.f.e.fr$dist,
     pch = 16, 
     col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3), 
     cex = 1.5,las=1,
     xlab="Relative Geographic Distance", 
     ylab ="High Affinty Community Similarity")
abline(TH.f.mno.dd,col="red")

#USC env vs dist
plot(TH.f.e.fr$struc ~ TH.f.e.fr$env,
     pch = 16, 
     col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3), 
     cex = 1.5,las=1,
     xlab="Relative Environmental Distance", 
     ylab ="High Affinty Community Similarity")
abline(TH.f.env.dd,col="red")


### Plot
#creating 1 x 2 plot
par(mfrow=c(1,2), mai = c(0.7, 0.8, 0.3, 0.5),
    mgp=c(2.2,1,0)) # export as 8 X 4

#USC geog dist
plot(TH.f.e.fr$struc ~ TH.f.e.fr$dist,
     pch = 16, 
     col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3), 
     cex = 1.5,las=1,
     xlab="Relative Geographic Distance", 
     ylab ="High Affinty Community Similarity")
abline(TH.f.mno.dd,col="red")

#USC env vs dist
plot(TH.f.e.fr$env ~ TH.f.e.fr$dist,
     pch = 16, 
     col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3), 
     cex = 1.5,las=1,
     xlab="Relative Geographic Distance", 
     ylab ="Relative Environmental Similarity")
abline(TH.f.env.dd,col="red")
