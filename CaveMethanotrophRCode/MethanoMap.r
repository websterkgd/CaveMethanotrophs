#Creating maps of methanotrophs

rm(list=ls())
getwd()

#Attempting to load data into phyloseq
require(phyloseq)
require(biomformat)
require(vegan)
require(ggplot2)
require(simba)
require(devtools)
require(dplyr)
require(stringr)
require(mapdata)
require(ggmap)

#importing the shared file 
OTU <- import_mothur(mothur_shared_file = 'CvMcrb03_20190225.shared', parseFunction = parse_taxonomy_default)

# Importing the taxonomy file 
TAX <- import_mothur(mothur_constaxonomy_file = 'CvMcrb03_20190225.taxonomy', parseFunction = parse_taxonomy_default)

#create a physeq table for analysis
physeq = phyloseq(OTU, TAX)

#Creating transformations of the data 
f.physeq = transform_sample_counts(physeq, function(x) x / sum(x)) # fractional abundance transformation of the data

# import meta data
md <- read.table('sample_meta_data+d.txt', header = TRUE, sep = "", dec = ".", row.names = 1) 

#filtering CVL for plotting
md <- rbind(as.data.frame(md[1:25,]),
        as.data.frame(md[29:31,]),
        as.data.frame(md[33:43,])) # CVL filtered

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

#building 3 panel figure
E_usa <- map_data("state", 
                  region = c("North Carolina","Virginia","Washington D.C.",
                             "Maryland", "Tennessee", "Kentucky","Ohio",
                             "West Virginia","Indiana","Pennsylvania",
                             "Delaware","New Jersey"))

#Total High Affinty
#define color scale
TS.c <- (log1p(sample_sums(TUSC)*(10^2)))

#define size
TS.s <- log1p(sample_sums(TUSC))
TS.s <- replace(TS.s, TS.s==0, NA)

p.TS <- ggplot() + 
  geom_polygon(data = E_usa, aes(x=long, y = lat, group = group), 
               fill = NA, color = "dark gray") + 
  coord_fixed(1.3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle(expression(paste("USC Methanotrophs")))+
  geom_point(aes(x = -1*md[,8], y = md[,7]), 
             color = rgb(red = 0.9*TS.c/max(TS.c), 
                         green =0.3*TS.c/max(TS.c), 
                         blue = 0.05*TS.c/max(TS.c), alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*1.2*10^3*TS.s))+
  geom_point(aes(x = -90, y = 41), 
             color = rgb(red = 0.9*(log1p(0.001*(10^2))/max(TS.c)), 
                         green = 0.3*(log1p(0.001*(10^2))/max(TS.c)), 
                         blue = 0.05*(log1p(0.001*(10^2))/max(TS.c)), alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*1.2*10^3*log1p(0.001)))+
  geom_point(aes(x = -90, y = 40), 
             color = rgb(red = 0.9*(log1p(0.01*(10^2))/max(TS.c)), 
                         green =0.3*(log1p(0.01*(10^2))/max(TS.c)), 
                         blue = 0.05*(log1p(0.01*(10^2))/max(TS.c)), alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*1.2*10^3*log1p(0.01)))+
  geom_point(aes(x = -90, y = 38.5), 
             color = rgb(red = 0.9*(log1p(0.05*(10^2))/max(TS.c)), 
                         green = 0.3*(log1p(0.05*(10^2))/max(TS.c)),  
                         blue = 0.05*(log1p(0.05*(10^2))/max(TS.c)), alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*1.2*10^3*log1p(0.05)))+
  annotate("text", x=-88.5, y=41, label="0.1*'%'", parse=TRUE)+
  annotate("text", x=-88.5, y=40, label="1*'%'", parse=TRUE)+
  annotate("text", x=-88.5, y=38.5, label="5*'%'", parse=TRUE)

#T.II methanotrophs
#define color scale
T2.c <- (log1p(sample_sums(T.II)*(10^6)))

#define size
T2.s <- log1p(sample_sums(T.II))
T2.s <- replace(T2.s, T2.s==0, NA)

p.T2 <- ggplot() + 
  geom_polygon(data = E_usa, aes(x=long, y = lat, group = group), 
               fill = NA, color = "dark gray") + 
  coord_fixed(1.3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle(expression(paste("Type 2 Methanotrophs")))+
  geom_point(aes(x = -1*md[,8], y = md[,7]), 
             color = rgb(red = 0.7, 
                         green = 0.4*(T2.c/max(T2.c)), 
                         blue = 0.2, alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*10^5*T2.s))+
  geom_point(aes(x = -90.5, y = 41), 
             color = rgb(red = 0.7, 
                         green = 0.4*((log1p(0.00001*(10^6)))/max(T2.c)), 
                         blue = 0.2, alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*10^5*log1p(0.00001)))+
  geom_point(aes(x = -90.5, y = 40), 
             color = rgb(red = 0.7, 
                         green = 0.4*((log1p(0.0001*(10^6)))/max(T2.c)), 
                         blue = 0.2, alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*10^5*log1p(0.0001)))+
  geom_point(aes(x = -90.5, y = 38.5),
             color = rgb(red = 0.7,
                         green = 0.4*((log1p(0.0006*(10^6)))/max(T2.c)), 
                         blue = 0.2, alpha = 0.5),
             pch = 16,
             cex = sqrt((4/pi)*10^5*log1p(0.0006)))+
  annotate("text", x=-88.8, y=41, label="0.001*'%'", parse=TRUE)+
  annotate("text", x=-88.8, y=40, label="0.01*'%'", parse=TRUE)+
  annotate("text", x=-88.8, y=38.5, label="0.06*'%'", parse=TRUE)

#Tp.I methanotrophs
#define color scale
T1.c <- (log1p(sample_sums(Tp.I)*(10^6)))

#define size
T1.s <- log1p(sample_sums(Tp.I))
T1.s <- replace(T1.s, T1.s==0, NA)

p.T1 <- ggplot() + 
  geom_polygon(data = E_usa, aes(x=long, y = lat, group = group), 
               fill = NA, color = "dark gray") + 
  coord_fixed(1.3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle(expression(paste("Type 1 Methanotrophs")))+
  geom_point(aes(x = -1*md[,8], y = md[,7]), 
             color = rgb(red = 0.2, 
                         green = 0.8*(T1.c/max(T1.c)), 
                         blue = 0.6, alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*2*10^5*T1.s))+
  geom_point(aes(x = -90.5, y = 41), 
             color = rgb(red = 0.2, 
                         green = 0.8*((log1p(0.00001*(10^6)))/max(T1.c)), 
                         blue = 0.6, alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*2*10^5*log1p(0.00001)))+
  geom_point(aes(x = -90.5, y = 40), 
             color = rgb(red = 0.2, 
                         green = 0.8*((log1p(0.0001*(10^6)))/max(T1.c)), 
                         blue = 0.6, alpha = 0.5),
             pch = 16, 
             cex = sqrt((4/pi)*2*10^5*log1p(0.0001)))+
  geom_point(aes(x = -90.5, y = 38.5),
             color = rgb(red = 0.2,
                         green = 0.8*((log1p(0.0002*(10^6)))/max(T1.c)), 
                         blue = 0.6, alpha = 0.5),
             pch = 16,
             cex = sqrt((4/pi)*2*10^5*log1p(0.0002)))+
  annotate("text", x=-88.8, y=41, label="0.001*'%'", parse=TRUE)+
  annotate("text", x=-88.8, y=40, label="0.01*'%'", parse=TRUE)+
  annotate("text", x=-88.8, y=38.5, label="0.02*'%'", parse=TRUE)


#4 plot fix up.

require("cowplot")

#export as 8.5 x 10.5
plot_grid(p.TS, p.T2, p.T1,
          labels = c("A","B","C"),
          ncol = 2)

###################

