install.packages("installr")
library(installr)
updateR()

library(dada2); packageVersion("dada2")
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.20") # change the ref argument to get other versions
install.packages("dada2")
library(dada2); packageVersion("dada2")

install.packages("BiocManager")
library("BiocManager")
BiocManager::install("dada2", version = "3.11")
library(dada2); packageVersion("dada2")
install.packages("latticeExtra")
library(latticeExtra)
BiocManager::install("dada2", version = "3.15")
library(dada2); packageVersion("dada2")

devtools::install_version("Matrix", version = "1.3.2", repos = "http://cran.us.r-project.org")
library("Matrix")
remove.packages("Matrix")
packageVersion("Matrix")

library(dada2); packageVersion("dada2")
library("devtools")
library(dada2); packageVersion("dada2")
library("BiocManager")
library(dada2); packageVersion("dada2")
library("Matrix")


install.packages('dada2')


path <- "~/Desktop/NGS_run/Social_isolation_fastq" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
fnFs

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")
length(fnFs)
length(fnRs)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=11, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
out
write.csv(out, "~/Desktop/out.csv")

out#Learning rate
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, minOverlap = 8, filtRs, verbose=TRUE)
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

#Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy using Greengenes
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/NGS_run/silva_nr_v132_train_set.fa", multithread=TRUE)
#Assigning taxonomy using Silva
#taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/NGS_run/paired_end_social_isolation_fecal_samples_04132022/silva_nr_v132_train_set.fa", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print

df1 <- data.frame(taxa)
write.csv(df1, "~/Desktop/test100.csv")

#Phyloseq
install.packages("phyloseq")
install.packages("Biostrings")
BiocManager::install("Biostrings")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
samdf <- data.frame(Subject=subject)
rownames(samdf) <- samples.out

sample_data(ps)
#gender <- substr(subject,1,1)
#subject <- substr(subject,2,999)
#day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
#samdf$When <- "Early"
#samdf$When[samdf$Day>100] <- "Late"

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
gp <- tax_glom(ps, taxrank="Genus", NArm= TRUE) #Remove NA
ps.proportion <- transform_sample_counts(gp, function(OTU) OTU/sum(OTU))
#ps.top20 <- prune_taxa(top20, ps.top20)

ps.proportion
otu_table(ps.proportion)

#Bar plot
library(RColorBrewer)
nb.cols <- 16 #baesd on the number of taxons
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

#colours <- brewer.pal(64, "Set3") 
#colours
df2 <- psmelt(ps.proportion)
df2
write.csv(df2,"~/Desktop/test111.csv")
df2
p = ggplot(df2, aes(x=Subject, y=Abundance, fill=Family))
p = p + geom_bar(stat="identity", position="fill") + 
  scale_fill_manual(values=mycolors) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

p

#Generate OTU table 
otu <- otu_table(ps.proportion)

df_otu <- as.data.frame(otu)

df_otu

identical(colnames(df_otu), taxa_names(ps.proportion))
colnames(df_otu) <- as.data.frame(tax_table(ps.proportion))$Phylum
df_otu

write.csv(df_otu,"~/Desktop/test112.csv")

colnames(df_otu)
rownames(df_otu)
#######################################################################################################
#Construct phylogenetic tree
BiocManager::install("DECIPHER")
install.packages("RSQLite")
install.packages("RSQLite")
install.packages("phangorn")
remotes::install_github("KlausVigo/phangorn")
BiocManager::install("phangorn")

library(RSQLite)
library(DECIPHER)
library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
install.packages('phangorn')

install.packages("installr")
library(installr)
require(installr)
updateR()

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

dm <- dist.ml(phang.align)

treeNJ <- NJ(dm) # Note, tip order != sequence order

fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

#merge into phyloseq
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),phy_tree(fit))

ps


#######################################################################################################
#Dada2 alpha diversity

install.packages("phyloseq")
BiocManager::install("phyloseq")
remotes::install_github("vmikk/metagMisc")
install.packages("vegan")
install.packages('PhyloMeasures')
BiocManager::install("PhyloMeasures")
BiocManager::install("metagMisc") #faith PD calculation
library(metagMisc)
library(biomeUtils)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())
library(metagMisc)
library(vegan)
library(PhyloMeasures)
install.packages("~/Downloads/PhyloMeasures_2.1.tar", repos = NULL, type='source')

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
#gender <- substr(subject,1,1)
#subject <- substr(subject,2,999)
#day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject)
#samdf$When <- "Early"
#samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa), phy_tree(fitGTR$tree))

#Rooting tree
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))

ps
ps <- prune_samples(sample_names(ps) != "snd-44", ps) # Remove mock sample

View(tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
sample_data(ps)



phyloseq_phylo_div(ps, measures=(c("PD")))

source("http://bioconductor.org/biocLite.R")
BiocManager::install("biomeUtils")


gp <- tax_glom(ps, taxrank="Genus", NArm= TRUE) #Remove NA
estimate_richness(gp, measures=c("InvSimpson", "Shannon"))
plot_richness(gp, x="Subject", measures=c("Shannon", "InvSimpson"))

calculatePD(ps, justDF = FALSE, include_root = TRUE)

ad <- estimate_richness(gp)
ad$Shannon
dim(ad) 
dim(sg)
ad

ad <- phyloseq_phylo_div(ps,
  measures = c("PD"))
  

write.csv(ad, "~/Desktop/Shapira_Lab_Projects/SocialIsolation_gutmicrobiome/alphadiversity_PD.csv")
write.csv(ad, "~/Desktop/alphadiversity_PD.csv")

pd_snd <- read.csv("~/Desktop/Shapira_Lab_Projects/SocialIsolation_gutmicrobiome/Phylogenetic_diversity_snd.csv")
pd_snd

ggplot(pd_snd, aes(x=subject, y=PD)) + geom_point(size = 3) +
  theme(axis.text.x=element_text(angle=90,hjust=1))

sample_names(ps)
#Creating vector for sample type
single <- rep(c("1_single"),times=21)
single_F <- rep(c("single_F"),times=9)
single_M <- rep(c("single_M"),times=11)
single_F_2 <- rep(c("single_F"),times=1)

group_F <- rep(c("group_F"),times=11)
group_M <- rep(c("group_M"),times=12)


group <- rep(c("2_group"), times= 23)

#Creting vector for female only
single_F <- rep(c("single_F"),times=10)
group_F <- rep(c("group_F"),times=11)

#Creating vector for grouped female, grouped male, and singles
single <- rep(c("1_single"),times=21)
group_F <- rep(c("group_F"),times=11)
group_M <- rep(c("group_M"),times=12)

sg <- c(single, group_F, group_M)
sample_data(ps)$group <- factor(sg)
sample_data(ps)



sg_vector <- cbind(single, group)
sg <- c(single_M, group_M)
sg
ad_sampletype <- cbind(ad, sg)
ad_sampletype

sg <- c(single_F, single_M, single_F_2, group_F, group_M)
sg <- c(single_F, group_F)

sample_data(ps)$group_sex <- factor(sg)
sample_data(ps)$group <- factor(sg)
sample_data(ps)
#Creating vector for more defined groups
singleF <- rep(c("singleF"),times=9)
singleM <- rep(c("singleM"),times=11)
singleF_2 <- c("singleF")
group1 <- rep(c("group1"), times=3)
group2 <- rep(c("group2"), times=4)
group3 <- rep(c("group3"), times=4)
group4 <- rep(c("group4"), times=4)
group5 <- rep(c("group5"), times=4)
group6 <- rep(c("group6"), times=4)

#Creating vector for stress level groups



grouping <- c(singleF, singleM, singleF_2,group1,group2,group3,group4,group5,group6)
stress_level <- c('high','low', 'high','low','high','low','low','high', 'high', 'low', 'low', 'low','low', 'low' , 'high', 'high','low', 'high', 'high','low', 'high' ,'high' ,'high','low', 'high', 'high', 'low', 'low', 'high', 'low', 'high', 'low', 'low' , 'low','low', 'low', 'low', 'low', 'low', 'high', 'low','low', 'low', 'low')

sample_data(ps)$stress_level <- factor(stress_level)

sample_data(ps)


##########
#take average of shannon and sd for shannon groupby sg
##take average simpson and sd for shannon groupby sg

######################################################################################################
#PCA
# Transform data to proportions as appropriate for Bray-Curtis distances


ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ps.prop

ord.pcoa.wuni <- ordinate(ps.prop, method="PCoA", distance="bray")
ord.pcoa.wuni

ps.prop

sample_data(ps.prop)

ord.pcoa.wuni <- ordinate(ps.prop, method="PCoA", distance="wunifrac")

theme_set(theme_bw()) + theme_set(plot.background = element_blank(), panel.grid.major = element_blank())

theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

allGroupsColors <- c("grey0", "red", "deepskyblue")
p1 <- plot_ordination(ps.prop, ord.pcoa.wuni, color = "group", title="PCoA") + scale_color_manual(values = allGroupsColors) + geom_point(size = 5)
p1 + theme(panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank())

#p1 + geom_mark_ellipse(aes(fill=group),show.legend = FALSE)

geom_mark_ellipse()

library(ggplot2)
library(ggplot)
install.packages("ggforce")
library(ggforce)

BiocManager::install(mark_ellipse)


p1 + facet_wrap(~Phylum, 3)

library(vegan)

ordihull(ord.pcoa.wuni, groups, display = "sites", draw = c("lines","polygon", "none"),
         col = NULL, alpha = 127, show.groups, label = FALSE,  ...)




help("distance")

otu_table(ps.prop)

sample_data(ps)

######################################################################################################
#PERMANOVA
# samples x species as input
library(vegan)


otu <- otu_table(ps)
meta <- as(sample_data(ps), "data.frame")
meta

permanova <- adonis2(otu ~ group, data =meta, permutations=999, method="bray")

meta


print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])

permanova

# P-value
print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])

install.packages("pvclust")


library(vegan)
metadata <- as(sample_data(ps), "data.frame")
distance(ps)

adonis2(otu_table(ps))


?distance
permanova <- adonis2(distance(ps, 'unifrac') ~ group, data = metadata, permutations=999)

print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])

permanova

otu_table(ps)

distance(ps, 'wunifrac')
distance(ps, 'bray')

ps
wunifr <- UniFrac(ps, weighted=TRUE, normalized= TRUE, parallel = FALSE)
uunifr <- UniFrac(ps, weighted=FALSE, normalized= TRUE, parallel = FALSE)

permanova
permanova
print(as.data.frame(permanova$aov.tab)["stress_level", "Pr(>F)"])
sample_data(ps)

adonis2


######################################################################################################
#Weighted unifrac PCoA
library(ggplot2)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa), phy_tree(fitGTR$tree))

#Rooting tree
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))

ps
ps <- prune_samples(sample_names(ps) != "snd-44", ps) # Remove mock sample

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.p.wuni <- ordinate(ps, method="PCoA", distance="unifrac")

p1 <- plot_ordination(ps.prop, ord.p.wuni, color = "group", title="PCoA")
p1
p1 + facet_wrap(~more_grouping, ncol= 3)

ps.prop




######################################################################################################
#Plotting

library(RColorBrewer)
coul <- brewer.pal(5, "Set2") 
barplot(height=data$value, names=data$name, col=coul )

plot_bar(ps.proportion,fill="Genus")

###Customizing color using ggplot
phylums <- c('Proteobacteria','Bacteroidetes','Firmicutes')

df$Phylum[!df$Phylum %in% phylums] <- "Others"
df$Genus[!df$Phylum %in% phylums] <- "Others"

df$Genus[df$Phylum=="Proteobacteria" & 
            !df$Genus %in% c('Alcaligenaceae','Enterobacteriaceae')] <- "Other Protobacteria"

df$Genus[df$Phylum=="Bacteroidetes" &
            !df$Genus %in% c('Bacteroidaceae','Rikenellaceae','Porphyromonadaceae')] <- "Other Bacteroidetes"

df$Genus[df$Phylum=="Firmicutes" & 
            !df$Genus %in% c('Lactobacillaceae','Clostridiaceae','Ruminococcaceae','Lachnospiraceae')] <- "Other Firmicutes"

library(forcats)
library(dplyr)
df2 <- select(df, Sample, Phylum, Genus) %>%
  mutate(Phylum=factor(Phylum, levels=c(phylums, "Others")),
         Genus=fct_reorder(Genus, 10*as.integer(Phylum) + grepl("Others", Genus))) %>%
  
  group_by(Genus) %>%  # For this dataset only
  sample_n(100)         # Otherwise, unnecessary

colours <- brewer.pal(11, "Set3") 

gp.ch = subset_taxa(ps.proportion, Phylum == "Chlamydiae")
df2 <- psmelt(ps.proportion)
library(ggplot2)
ggplot(df2, aes(x=Subject, y = Abundance,fill = Genus)) + 
  geom_bar(position="fill", colour = "grey") +  # Stacked 100% barplot
  scale_fill_manual("", values=colours) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +  # Vertical x-axis tick labels
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y="Relative abundance")

colnames(df2)
head(df2)
df2
colours <- brewer.pal(11, "Set3") 
df2 <- psmelt(ps.proportion)

p = ggplot(df2, aes(x=Subject, y=Abundance, fill=Family))
p = p + geom_bar(color="black", stat="identity", position="stack") + 
  scale_fill_manual("", values=colours) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))
print(p)
