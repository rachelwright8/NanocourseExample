# Introduction to RNA-Seq Nanocourse
## Overall objective: Demonstrate gene expression analysis methods
### 1. Sample clustering
### 2. Functional enrichment
### 3. Graphical representation

# SET WORKING DIRECTORY
setwd("~/Desktop/Nanocourse/")

# INSTALL PACKAGES (ONLY DO THIS ONCE)
# source("https://bioconductor.org/biocLite.R")
# biocLite("airway")
# biocLite("DESeq2")

# LOAD DATA -----
library(airway)
data(airway) #loads RangedSummarizedExperiment called "airway" that contains the read counts per gene in treated and untreated smooth muscle cells.

# see where the files were downloaded on your computer
dir <- system.file("extdata", package="airway", mustWork=TRUE)

# look at the files
list.files(dir) 

# Extract the sample table that contains sample metadata
sampleTable <- read.csv(file.path(dir,"sample_table.csv"),row.names=1)
sampleTable

# GET COUNT DATA -----
library("Rsamtools") # this packages calculates gene counts from bam files

# extract the bam files
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000) # process 2mil reads at a time

seqinfo(bamfiles[1])

# define gene models using the human gtf file in the airway package
library("GenomicFeatures")
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
# produce GRangesList of exons grouped by gene
(ebg <- exonsBy(txdb, by="gene"))

# read counting
library("GenomicAlignments")
library("BiocParallel")

register(SerialParam()) # use multiple cores

# CREATE SUMMARIZED EXPERIMENT ------
# create a summarized experiment that counts reads mapped to a genome
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )


# various ways of looking at this object
se
dim(se)
head(assay(se),3)
str(metadata(rowRanges(se)))

# where are the metadata?
colData(se) # this is empty right now
(colData(se) <- DataFrame(sampleTable)) # add your sample data to the colData (metadata) "slot" in the summarized experiment object

# DIFFERENTIAL GENE EXPRESSION ANALYSIS USING DESEQ2
library(DESeq2) # for differential gene expression analysis
library(arrayQualityMetrics) # for outlier analysis

summary(as.data.frame(colData(airway)))

# create a DESeq data set
dds = DESeqDataSet(airway, design= ~ 1 + dex)

# Find sample outliers
rld = rlogTransformation(dds, blind=T)
e = ExpressionSet(assay(rld), AnnotatedDataFrame(as.data.frame(colData(rld))))
arrayQualityMetrics(e, intgroup = c("cell"), force=T)
# no outliers

# DESeq pipeline
deds <- DESeq(dds)

# results by treatment
res = results(deds, contrast=c("dex", "trt", "untrt"))
summary(res)
# LFC > 0 (up)     : 1877, 5.6% 
#LFC < 0 (down)   : 1505, 4.5% 
head(res)
# log2 fold change (MAP): dex trt vs untrt 
# Wald test p-value: dex trt vs untrt 

# write results
write.csv( as.data.frame(res), file = "resultsTable.csv")

# Get RLD  ----------------------------------------------------------------
rld = rlogTransformation(deds)
head(assay(rld))
colnames(rld) = paste(coldataOut$gen, coldataOut$vlp, coldataOut$heat,sep="")
head(assay(rld))

# rename samples to something meaningful
colnames(rld)
head(sampleTable)
colnames(rld) = paste(sampleTable$cell, sampleTable$dex, sep="_")


# LOAD/SAVE (start here once you've done the above steps once)---------
save(airway, deds, res, sampleTable, rld, file = "DESeq2data.Rdata")
load("DESeq2data.Rdata")

# Explore with plots -------------------------------------------------
library(pheatmap) # for "pretty" heatmaps

quartz() # or windows() if you're using a PC
pheatmap(cor(assay(rld)),border_color=NA, main="Sample Heatmap")

# Diagnostics ------------------------------------------------------

#####-------------Dispersions plot
quartz()
plotDispEsts(deds, main="Dispersion Plot Baseline")

####-----------MA plot
maplot <- function (res, thresh=0.05, thresh2=1, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh, log2FoldChange>thresh2), points(baseMean, log2FoldChange, col="red", pch=20, cex=0.75))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=NA, cex=textcx, col=2))
  }
}
quartz()
maplot(res, main="MA Plot")

# Write results for making heatmaps ---------------------------------------

#Get pvals
head(res)
dim(res)
#24475     6
rldd = assay(rld)
head(rldd)
rlddFilt = rldd[row.names(rldd) %in% row.names(res),]
head(rlddFilt)
dim(rlddFilt)
#64102     8

vals = cbind("pval"=res$pvalue,"padj"=res$padj)
head(vals)
length(vals[,1])
table(complete.cases(vals))
# FALSE  TRUE 
# 45536 18566 

#Make rlogdata and pvals table
head(rlddFilt)
rldpvals = as.data.frame(cbind(rlddFilt,vals))
head(rldpvals)
dim(rldpvals)
# 64102    10

# subset for not NA
rldpvals = rldpvals[!is.na(rldpvals$padj),]
dim(rldpvals)
# 18566    10
table(complete.cases(rldpvals))

write.csv(rldpvals, "RLDandPVALS.csv", quote= F)

# Write results for GO/KOG analysis by negative log pvalue-----------------------------------

head(res)
logs = data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
logs$logP = as.numeric(as.character(logs$logP))
sign = rep(1,nrow(logs))
sign[res$log2FoldChange<0]=-1  ##change to correct model
table(sign)
#-1     1 
#17001 47101 
logs$logP = logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_logP.csv",sep=",")

# Get LFC for GO/KOG analysis ---------------------------------------------------------

head(res)
lfc = as.data.frame(res[2])
head(lfc)
write.csv(lfc, "LFC.csv", quote=F)

# Get GO/gene annotations ------------------------------------------------------------------
#biocLite("biomaRt") # only run once
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
biomartResults <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "go_id", "description"),
                        filters = "ensembl_gene_id", values = c(row.names(rld)),
                        mart = mart)
head(biomartResults)
biogo = biomartResults[c(1,3)]
head(biogo)
write.table(biogo, row.names=F,file = "biomartGO.tab", sep="\t", quote=F)
# NRify in perl for goMWU analysis

biogene = biomartResults[c(1,4)]
head(biogene)
write.table(biogene, row.names=F,file = "biomartGENE.tab", sep="\t", quote=F)

# GO ----------------------------------------------------------------------
source("gomwu.functions.R") # make sure this is in your working directory

# Edit these to match your data file names: 
input="GO_logP.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="biomartGONR.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", 
           largest=0.1,  
           smallest=5,   
           clusterCutHeight=0.25)

# Plotting results
quartz()
gomwuPlot(input,goAnnotations,goDivision,
          absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.5 if you are doing Fisher's exact test for standard GO enrichment.
          level1=0.01, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
          level2=0.005, # FDR cutoff to print in regular (not italic) font.
          level3=0.001, # FDR cutoff to print in large bold font.
          txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
          treeHeight=0.3, # height of the hierarchical clustering tree
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.01,level2=0.005,level3=0.001.  

# PCoA --------------------------------------------------------------------
# principal coordinate analysis 
library(vegan)
library(rgl)
library(ape)

rldd = assay(rld)
head(rldd) #expression data to use
head(sampleTable) #conditions data to use


# principal coordinate calculation
dd.veg = vegdist(t(rldd), "manhattan")
div.dd.veg = dd.veg/1000
head(div.dd.veg)

dd.pcoa = pcoa(div.dd.veg) 
head(dd.pcoa)
scores = dd.pcoa$vectors

# PERMANOVA  
adonis(t(rldd)~cell,data=sampleTable,method="manhattan") #0.168
adonis(t(rldd)~dex,data=sampleTable,method="manhattan") #0.028

pco1=scores[,1]
TukeyHSD(aov(pco1~sampleTable$cell))# all different
TukeyHSD(aov(pco1~sampleTable$dex))# 5.5e-06

pco2=scores[,2]
TukeyHSD(aov(pco2~sampleTable$cell))# look at N080611-N061011 and N61311-N080611
TukeyHSD(aov(pco2~sampleTable$dex))#0.8512763

pco3=scores[,3]
TukeyHSD(aov(pco3~sampleTable$cell))
TukeyHSD(aov(pco3~sampleTable$dex))

# make colors
cellCol = ifelse(sampleTable$cell == "N052611", "red", ifelse(sampleTable$cell == "N061011", "orange", ifelse(sampleTable$cell == "N080611", "black", "blue")))
cellCol

dexShape = ifelse(sampleTable$dex == "trt", 17, 16)
dexShape

# Plot (ordispider by cell type)
quartz()
plot(scores[,1], scores[,2], col = cellCol, xlab="PCo1 (27.7)", pch = dexShape, ylab="PCo2 (17.6)", main="PCoA of Smooth Muscle Cells by Line and Treatment")
ordispider(scores,sampleTable$cell,label=F)
legend(locator(1), col=c("red", "orange", "black", "blue"), legend=c("N052611", "N061011", "N080611", "N61311"), pch=16, lwd=3, cex=0.8, bty="n")
legend(locator(1), col="black", legend=c("treated", "untreated"), pch=c(16,17), lwd=1, cex=0.8, bty="n")

# Plot (ordispider by treatment )
quartz()
plot(scores[,1], scores[,2], col = cellCol, xlab="PCo1 (27.7)", pch = dexShape, ylab="PCo2 (17.6)", main="PCoA of Smooth Muscle Cells by Line and Treatment")
ordispider(scores,sampleTable$dex,label=F)
legend(locator(1), col=c("red", "orange", "black", "blue"), legend=c("N052611", "N061011", "N080611", "N61311"), pch=16, lwd=3, cex=0.8, bty="n")
legend(locator(1), col="black", legend=c("treated", "untreated"), pch=c(16,17), lwd=1, cex=0.8, bty="n")

# Creating gene expression heatmaps -------------


# heatmap -----------------------------
library(pheatmap)

# read in annotations file
gg = read.delim("biomartGENE.tab")
head(gg)

# read in RLD and pvals
d0 = read.csv("RLDandPVALS.csv")
head(d0)
row.names(d0) = d0$X
d0$X = NULL
d = d0[c(1:8)]
head(d)

pvals = d0[c(9:10)]
head(pvals)

# subset by pvals
cutoff = 1e-20
sig = d[pvals$padj<cutoff & !is.na(pvals$padj),]
length(sig[,1])

conds = d[(row.names(d) %in% row.names(sig)),]
explc=conds

# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(explc)){
  s=subset(gg,ensembl_gene_id==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$description[1],i,sep="."))
    expg=rbind(expg,explc[i,])
  } 
}
row.names(expg)=gnames
head(expg)
dim(expg) 
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them
head(explc)
names(explc)


heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.6)(100)

quartz()
ph=pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,border_color=NA,clustering_distance_rows="correlation", cluster_cols=T)
