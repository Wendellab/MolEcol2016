###### Start analysis on Biocrunch
# ssh biocrunch.ent.iastate.edu
# cd ~/jfw-lab/Papers/GallagherMolecularEcology2015/ProfilinAnalysis/
# R


############### Step 1. Basic data processing and cleaning  ############### 
################

getwd()
library(WGCNA);
library(RColorBrewer)
library(flashClust);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

################## load AD1 expression data
data = read.table("AD1.fiber.counts",header=TRUE, sep="\t");
# Take a quick look at what is in the data set:
dim(data);  #37223   25
names(data);
names(data)<-gsub("AD1_|.sort.bam","",names(data))
names(data)
#  Make each row corresponds to a gene and column to a sample or auxiliary information.
datExprT = as.data.frame(t(data[, -1]));
names(datExprT) = data[,1];
pdf("s1.AD.boxplot.pdf")
boxplot(log2(t(datExprT) ), las=2)
dev.off()
## note maxxa_20dpa doesn't look good; anyway, proceed
ADdata<-datExprT
ADdpa<-as.numeric(gsub(".*_|dpa","",rownames(ADdata)))

# Check outliers
# sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(t(ADdata),type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor1=ifelse(Z.k<thresholdZ.k,"red","black")
thresholdZ.k= -2.5
outlierColor2=ifelse(Z.k<thresholdZ.k ,"red","black")

pdf("s1.AD1.sample_dendrogram_and_trait_heatmap.pdf")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation:
# where red indicates high values
dpaColors=numbers2colors(ADdpa,signed=FALSE)
datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2, dpa=dpaColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
colors=datColors,main="Sample dendrogram and trait heatmap")
dev.off()
## Samples clubstering show no outliers, but neither make much sense


# repeat above code for Adata
#################### load A1 expression data
data = read.table("AD1.fiber.counts",header=TRUE, sep="\t");
# Take a quick look at what is in the data set:
dim(data);  #37223   25
names(data);
names(data)<-gsub(".trim.fq.gz.sort.bam","",names(data))
names(data)
#  Make each row corresponds to a gene and column to a sample or auxiliary information.
datExprT = as.data.frame(t(data[, -1]));
names(datExprT) = data[,1];
pdf("s1.Adip.boxplot.pdf")
boxplot(log2(t(datExprT) ), las=2)
dev.off()
## note maxxa_20dpa doesn't look good; anyway, proceed first
Adata<-datExprT
Adpa<-as.numeric(gsub(".*_|dpa","",rownames(Adata)))

# Check outliers
# sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(t(Adata),type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor1=ifelse(Z.k<thresholdZ.k,"red","black")
thresholdZ.k= -2.5
outlierColor2=ifelse(Z.k<thresholdZ.k ,"red","black")

pdf("s1.Adip.sample_dendrogram_and_trait_heatmap.pdf")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation:
# where red indicates high values
dpaColors=numbers2colors(Adpa,signed=FALSE)
datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2, dpa=dpaColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
colors=datColors,main="Sample dendrogram and trait heatmap")
dev.off()
## Samples clubstering show JMS_10dpa as outliers
## samples are actually clustered by dpa better 


################
info<-read.table(file="accesionInfo.txt",header=TRUE,sep="\t")
library(ggplot2)

# plot PCA, check clustering and outlier again
pdf("s1.pca.pdf")

# AD1
data = prcomp(ADdata)
summary(data)
coldata<-data.frame( sample = rownames(ADdata), accession = gsub("_.*","", rownames(ADdata)), dpa = gsub(".*_","", rownames(ADdata)) )
coldata<-merge(coldata,info, all.x=TRUE,by="accession")
qplot(PC1, PC2, main="AD1",color=coldata$domestication, shape=coldata$dpa, data=as.data.frame(data$x)) +
xlab(paste0("PC1: ",summary(data)[[6]][2,1]*100, "% variance")) +
ylab(paste0("PC1: ",summary(data)[[6]][2,2]*100, "% variance"))

# A1
data = prcomp(Adata)
summary(data)
coldata<-data.frame( sample = rownames(Adata), accession = gsub("_.*","", rownames(Adata)), dpa = gsub(".*_","", rownames(Adata)) )
coldata<-merge(coldata,info, all.x=TRUE,by="accession")
qplot(PC1, PC2, main="A1",color=coldata$domestication, shape=coldata$dpa, data=as.data.frame(data$x)) +
xlab(paste0("PC1: ",summary(data)[[6]][2,1]*100, "% variance")) +
ylab(paste0("PC1: ",summary(data)[[6]][2,2]*100, "% variance"))

dev.off()


############### Step 2 Use DESeq2 rlog normalization  ###############
################
getwd()
library(WGCNA);
library(RColorBrewer)
library(flashClust);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

info<-read.table(file="accesionInfo.txt",header=TRUE,sep="\t")
library(DESeq2)


################## load AD1 expression data
data = read.table("AD1.fiber.counts",header=TRUE, sep="\t");
# Take a quick look at what is in the data set:
dim(data);  #37223   25
names(data);
names(data)<-gsub("AD1_|.sort.bam","",names(data))
names(data)
#  Make each row corresponds to a gene and column to a sample or auxiliary information.
datExprT = as.data.frame(t(data[, -1]));
names(datExprT) = data[,1];
ADdata<-datExprT
ADdpa<-as.numeric(gsub(".*_|dpa","",rownames(ADdata)))

# normalization using DESeq2 rlog
count<-as.data.frame( t(datExprT) )
coldata<-data.frame( sample = rownames(ADdata), accession = gsub("_.*","", rownames(ADdata)), dpa = gsub(".*_","", rownames(ADdata)) )
coldata<-merge(coldata,info, all.x=TRUE,by="accession")
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
# rlog transformation, note that with default blind=TRUE, design is actually ~1
rld <- rlog(dds)
head(assay(rld))
# write transformed data
expr<-as.data.frame(assay(rld) )
names(expr)<-names(count)
write.table(expr,"s2.AD1.rlog.txt", sep="\t")
ADdata_rld<-t(expr)


# plots
pdf("s2.AD1.rlog.pdf")

#PCA plot
#Another way to visualize sample-to-sample distances is a principal-components analysis (PCA). In this ordination method, the data points (i.e., here, the samples) are projected onto the 2D plane such that they spread out in the two directions which explain most of the differences in the data. The x-axis is the direction (or principal component) which separates the data points the most. The amount of the total variance which is contained in the direction is printed in the axis label.
plotPCA(rld, intgroup = c("domestication", "dpa"))
# Here, we have used the function plotPCA which comes with DESeq2. The two terms specified by intgroup are the interesting groups for labeling the samples; they tell the function to use them to choose colors.

# We can also build the PCA plot from scratch using ggplot2. This is done by asking the plotPCA function to return the data used for plotting rather than building the plot. See the ggplot2 documentation for more details on using ggplot.
library(genefilter)
sumPCA<-
function (x, intgroup = "condition", ntop = 500)
{
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]),
    1, paste, collapse = " : "))
    return(pca)
}
data<-sumPCA(rld, intgroup = c("domestication", "dpa"))
summary(data)
library(ggplot2)
qplot(PC1, PC2, main="top 500", color=coldata$domestication, shape=coldata$dpa, data=as.data.frame(data$x)) +
xlab(paste0("PC1: ",summary(data)[[6]][2,1]*100, "% variance")) +
ylab(paste0("PC2: ",summary(data)[[6]][2,2]*100, "% variance"))
# From both PCA visualizations, we see that the differences between genomes are not stronger than the differences due to dpa.

# all genes
data<-sumPCA(rld, intgroup = c("genome", "dpa"), ntop=37223)
qplot(PC1, PC2, main="all 37223", color=coldata$domestication, shape=coldata$dpa, data=as.data.frame(data$x)) +
xlab(paste0("PC1: ",summary(data)[[6]][2,1]*100, "% variance")) +
ylab(paste0("PC2: ",summary(data)[[6]][2,2]*100, "% variance"))
# From both PCA visualizations, we see that the differences between genomes are not stronger than the differences due to dpa.

# Check outliers
# sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(expr,type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor1=ifelse(Z.k<thresholdZ.k,"red","black")
thresholdZ.k= -2.5
outlierColor2=ifelse(Z.k<thresholdZ.k ,"red","black")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation:
# where red indicates high values
dpaColors=numbers2colors(ADdpa,signed=FALSE)
datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2, dpa=dpaColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
colors=datColors,main="Sample dendrogram and trait heatmap")

dev.off()
# clustering looks a lot better, grouping by dpa then domestication


################## load A1 expression data
data = read.table("Adip.fiber.counts",header=TRUE, sep="\t");
# Take a quick look at what is in the data set:
dim(data);  #37223   25
names(data);
names(data)<-gsub(".trim.fq.gz.sort.bam","",names(data))
names(data)
#  Make each row corresponds to a gene and column to a sample or auxiliary information.
datExprT = as.data.frame(t(data[, -1]));
names(datExprT) = data[,1];

Adata<-datExprT
Adpa<-as.numeric(gsub(".*_|dpa","",rownames(Adata)))

# normalization using DESeq2 rlog
count<-as.data.frame( t(datExprT) )
coldata<-data.frame( sample = rownames(Adata), accession = gsub("_.*","", rownames(Adata)), dpa = gsub(".*_","", rownames(Adata)) )
coldata<-merge(coldata,info, all.x=TRUE,by="accession")
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
# rlog transformation, note that with default blind=TRUE, design is actually ~1
rld <- rlog(dds)
head(assay(rld))
# write transformed data
expr<-as.data.frame(assay(rld) )
names(expr)<-names(count)
write.table(expr,"s2.Adip.rlog.txt", sep="\t")
Adata_rld<-t(expr)

# plots
pdf("s2.Adip.rlog.pdf")

#PCA plot
#Another way to visualize sample-to-sample distances is a principal-components analysis (PCA). In this ordination method, the data points (i.e., here, the samples) are projected onto the 2D plane such that they spread out in the two directions which explain most of the differences in the data. The x-axis is the direction (or principal component) which separates the data points the most. The amount of the total variance which is contained in the direction is printed in the axis label.
plotPCA(rld, intgroup = c("domestication", "dpa"))
# Here, we have used the function plotPCA which comes with DESeq2. The two terms specified by intgroup are the interesting groups for labeling the samples; they tell the function to use them to choose colors.

# We can also build the PCA plot from scratch using ggplot2. This is done by asking the plotPCA function to return the data used for plotting rather than building the plot. See the ggplot2 documentation for more details on using ggplot.
library(genefilter)
sumPCA<-
function (x, intgroup = "condition", ntop = 500)
{
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]),
    1, paste, collapse = " : "))
    return(pca)
}
data<-sumPCA(rld, intgroup = c("domestication", "dpa"))
summary(data)
library(ggplot2)
qplot(PC1, PC2, main="top 500", color=coldata$domestication, shape=coldata$dpa, data=as.data.frame(data$x)) +
xlab(paste0("PC1: ",summary(data)[[6]][2,1]*100, "% variance")) +
ylab(paste0("PC2: ",summary(data)[[6]][2,2]*100, "% variance"))
# From both PCA visualizations, we see that the differences between genomes are not stronger than the differences due to dpa.

# all genes
data<-sumPCA(rld, intgroup = c("genome", "dpa"), ntop=37223)
qplot(PC1, PC2, main="all 37223", color=coldata$domestication, shape=coldata$dpa, data=as.data.frame(data$x)) +
xlab(paste0("PC1: ",summary(data)[[6]][2,1]*100, "% variance")) +
ylab(paste0("PC2: ",summary(data)[[6]][2,2]*100, "% variance"))
# From both PCA visualizations, we see that the differences between genomes are not stronger than the differences due to dpa.

# Check outliers
# sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(expr,type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor1=ifelse(Z.k<thresholdZ.k,"red","black")
thresholdZ.k= -2.5
outlierColor2=ifelse(Z.k<thresholdZ.k ,"red","black")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation:
# where red indicates high values
dpaColors=numbers2colors(Adpa,signed=FALSE)
datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2, dpa=dpaColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
colors=datColors,main="Sample dendrogram and trait heatmap")

dev.off()
# clustering looks a lot better, grouping by dpa then domestication

############################
# Next we make a multi-set data, considering 3 different combinations
nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("AD1 wild and domesticated genomes", "A1 wild and domesticated genomes")
shortLabels = c("AD1", "A1")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = ADdata_rld);
multiExpr[[2]] = list(data = Adata_rld);

# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)
# $nSets 3
# $nGenes 37223
# $nSamples 24 24
# $structureOK  TRUE

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding ZERO genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
# Update exprSize
checkSets(multiExpr)
# $nSets 3
# $nGenes 35580
# $nSamples 24 24
# $structureOK  TRUE

# check if profilin genes were kept
profilin<- c( "Gorai.009G028500", "Gorai.009G028400", "Gorai.004G282400", "Gorai.013G229500", "Gorai.010G078400", "Gorai.001G025300", "Gorai.003G061200")
match(profilin, colnames(multiExpr[[1]]$data))
# 21811 21810  9655 35157 26729   249  5743

pdf(file = "s2.SampleClusteringS.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;

save(multiExpr,nSets, setLabels, shortLabels, file = "R-02-dataInput.RData")



############### Step 3.  Choosing the soft-thresholding power: analysis of network topology  ###############
################

library(WGCNA)
library(RColorBrewer)
library(ggplot2);
options(stringsAsFactors = FALSE)

lnames = load(file = "R-02-dataInput.RData")
lnames
nSets
nGenes<-checkSets(multiExpr)$nGenes #35580
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels
shortLabels


# Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
types<-c("unsigned", "signed", "signed hybrid")

for (type in types)
{
    powers = c(c(1:10), seq(from = 12, to=40, by=2))
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type)[[2]])      }
    collectGarbage()
    
    # Plot the results:
    colors=brewer.pal(5,"Set1")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
        for (col in 1:length(plotCols))
        {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        }
    }
    
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    # sizeGrWindow(8, 6)
    pdf(paste("s3.ChooseSoftThresholdPower_",gsub(".* ","", type), ".pdf", sep="") )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
        if (set==1)
        {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
            addGrid()
        }
    if (col==1)
    {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
        labels=powers,cex=cex1,col=colors[set]);
    } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
    if (col==1)
    {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
    } else
    legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
    
    dev.off()
    assign(paste("powerTables.",gsub(".* ","", type),sep=""),powerTables)
# repeat above for powerTables.signed, and powerTables.hybrid
}

save(powerTables.unsigned,powerTables.signed,powerTables.hybrid , file = "R-03-choosePower.RData")


# Inspect "s3.ChooseSoftThresholdPower_signed.pdf". Power=12 is good.
# Basically, the default power is 6 for unsigned network, 12 for signed network; if the observed power choice is smaller than default, I will choose the observed power, otherwise use default.





############### Step 4.  Network Construction  ###############
################
library(WGCNA)
options(stringsAsFactors = FALSE);
library(RColorBrewer)
library(ggplot2);

lnames = load(file = "R-02-dataInput.RData")
lnames     #  "multiExpr"   "nSets"       "setLabels"   "shortLabels" "datTraits"
nSets      # 2
setLabels  #  "AD1 wild and domesticated genomes" "A1 wild and domesticated genomes"
shortLabels#  "AD1" "A1"
nGenes<-checkSets(multiExpr)$nGenes
nGenes     #  35580

powers<-c(12)

# work with individual genome, then work with different soft threshold
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    for (j in powers )
    {
        softPower = j
        print(paste("Start building network for ",shortLabels[set]," using soft threshold ",j,"......",sep=""))
        # Network construction
        
        net = blockwiseModules(
             # Input data
             subDat,
             # Data checking options
             checkMissingData = TRUE,
             
             # Options for splitting data into blocks
             blocks = NULL,
             randomSeed = 12345,
             maxBlockSize = nGenes,  # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
             
             # Network construction arguments: correlation options, use bicor instead of default pearson
             corType = "pearson",
             # Adjacency and topology overlap function options
             power = j, networkType = "signed", TOMType = "signed",
             
             # Saving or returning TOM
             saveTOMs = TRUE,
             saveTOMFileBase = paste(shortLabels[set],"_power",j,"_TOM",sep=""),
             
             # Basic tree cut options
             deepSplit = 2,  #default, known to reasonable
             minModuleSize = min(30, ncol(subDat)/2 ), #default 20, use 30 for transcriptome
             pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
             
             # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
             mergeCutHeight = 0.25,
             
             # others
             reassignThreshold = 0,
             numericLabels = TRUE,
             verbose = 3)
             
        assign(paste(shortLabels[set],"net",j,sep=""), net)
        }
}
save(list=grep(".+net.+",ls(), value=TRUE), file = "R-04-buildNetwork.RData")

# OK this is pretty tricky here. The use of "bicor" correlation is supposed to be more robust than "pearson" correlation, because it is good at removing outliers from correlation calculation. BUT we have a different problem here: "pearson" network results in over 30 modules, while "bicor" gives only 2 modules; "pearson" resulted modules are much better correlated with gene significance than "bicor" modules. This is because some genome&dpa specific expressions (represented by only 3 samples) were treated as outliers, and their effects were removed. Below code and pdf output illustrated this point. SO we should stick to "Pearson" correlation.
# below code is optional, and I will skip it this time.
##########################
pdf("s3.CompareBicor&Pearson.pdf")
netB = blockwiseModules( multiExpr[[1]]$data, power = 12, maxBlockSize = nGenes, corType = "bicor", networkType = "signed", TOMType = "signed")
netP = blockwiseModules( multiExpr[[1]]$data, power = 12, maxBlockSize = nGenes, corType = "pearson", networkType = "signed", TOMType = "signed")
# Use seed weight to define a gene significance variable
GS.weight=as.numeric(cor(multiExpr[[1]]$data,datTraits$seed_weight,use="p"))
# This translates the numeric values into colors
GS.weightColor=numbers2colors(GS.weight,signed=T)
# Use oil to define a gene significance variable
GS.oil=as.numeric(cor(multiExpr[[1]]$data,datTraits$oil_content,use="p"))
# This translates the numeric values into colors
GS.oilColor=numbers2colors(GS.oil,signed=T)
# Use dpa to define a gene significance variable
GS.dpa=as.numeric(cor(multiExpr[[1]]$data,datTraits$dpa,use="p"))
# This translates the numeric values into colors
GS.dpaColor=numbers2colors(GS.dpa,signed=T)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netB$dendrograms[[1]],colors=data.frame(labels2colors(netB$colors),labels2colors(netP$colors), GS.weightColor, GS.oilColor, GS.dpaColor), groupLabels=c("Bicor","Pearson", "Seed_weight", "oil_content", "dpa"),dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05, main="Cluster dendrogram, bicor")
plotDendroAndColors(netP$dendrograms[[1]],colors=data.frame(labels2colors(netP$colors),labels2colors(netB$colors), GS.weightColor, GS.oilColor, GS.dpaColor), groupLabels=c("Pearson","Bicor", "Seed_weight", "oil_content", "dpa"),dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05, main="Cluster dendrogram, pearson")
dev.off()
##########################


############### Step 5.  About profilin  ###############
################
library(WGCNA)
options(stringsAsFactors = FALSE);

lnames = load(file = "R-02-dataInput.RData")
lnames     #  "multiExpr"   "nSets"       "setLabels"   "shortLabels" "datTraits"
nSets      # 2
setLabels  #  "AD1 wild and domesticated genomes" "A1 wild and domesticated genomes"
shortLabels#  "AD1" "A1"
nGenes<-checkSets(multiExpr)$nGenes
nGenes     #  35580

lnames = load(file = "R-04-buildNetwork.RData")

# 8 profilin genes
profilin<- c( "Gorai.009G028500", "Gorai.009G028400", "Gorai.004G282400", "Gorai.013G229500", "Gorai.010G078400", "Gorai.001G025300", "Gorai.003G061200", "Gorai.010G078500")
PRF<-c("PRF1","PRF2","PRF3","PRF4","PRF5","profilin","profilin","PRF6")
prf<-match(profilin, colnames(multiExpr[[1]]$data))
# 21811 21810  9655 35157 26729   249  5743

# Plot profilin gene expression profiles in AD1 and A1
multiExpr[[1]]$data[,"Gorai.009G028500"]
info<-read.table(file="accesionInfo.txt",header=TRUE,sep="\t")
# AD1 and A1 coldata
coldata<-list()
for(i in 1:nSets)
{
    coldata[[i]]<-data.frame( sample = rownames(multiExpr[[i]]$data), accession = gsub("_.*","", rownames(multiExpr[[i]]$data)), dpa = gsub(".*_","", rownames(multiExpr[[i]]$data)) )
    coldata[[i]]<-merge(coldata[[i]],info, all.x=TRUE,by="accession")
    coldata[[i]]$dpa <- factor(gsub("dpa","",coldata[[i]]$dpa) , levels=c("5", "10", "15","20"))

}
# do plots
# NOTE: failed on biocrunch, I had to do these plots locally
source("summarySE.r")
library(ggplot2)
pdf("s4.profilinExpressionProfiles.pdf")
for(n in 1:8) {
    plots <- list()
    for(i in 1:2) #loop for AD1 and A1
    {
    ss<-as.factor(coldata[[i]]$sample)
    # barplot
    # par(mfrow=c(2,1))
    # barplot(multiExpr[[i]]$data[,"Gorai.009G028500"],  main="", cex.main=2, ylab="RNA-seq with rlog transformation",xlab="fiber development (dpa)", names.arg=ss, las=2)
    
    # ggplot line, anova
    df<-data.frame(expr=multiExpr[[i]]$data[,profilin[n] ], coldata[[i]] )
    fit<-aov(expr~sample,df)
    dfc<-summarySE(df, measurevar="expr", groupvars=c("dpa", "domestication"))
    
    plots[[i]]<-ggplot(dfc, aes(x=dpa, y=expr, group = domestication, colour=domestication)) +
    geom_line() +  geom_point() +
    geom_errorbar(aes(ymin=expr-se, ymax=expr+se), width=.3,) +
    ggtitle(paste(shortLabels[[i]],"-",  PRF[n], "-",profilin[n] ))+
    theme_bw() +
    theme(plot.title=element_text( size=11),legend.position = "bottom")
    }
    multiplot(plotlist = plots,  layout=matrix(1:4, nrow=2, byrow=TRUE) )
}
dev.off()
# examine pdf file: PRF1, 2, 3

# Are all profilin genes assigned to the same module?
AD1net12$colors[prf]
#  5  5 12  2  1  1 12 3
A1net12$colors[prf]
#  2 10  7  7  6  1  6 2
# well, guess not.


# Look at module eigengenes
for (set in 1:2 ) # AD1, A1
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    genome <-  shortLabels[set]
    net<-get(paste(genome,"net12",sep="") )
    gg<-net$goodGenes    #get the good genes descision
    #    adjacency = adjacency(subDat, power = 12, type = "signed")
    Nmodules= dim(net$MEs)[2]
    
    #   load(net$TOMFiles)
    print(net$TOMFiles)
    print(paste("Number of modules in ",genome," network is ",Nmodules,sep=""))
    colorsa<- labels2colors(net$colors)
    
    pdf(paste("s4.",genome,"_modules.pdf",sep=""))
    
    # Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
    # Displaying module heatmap and the eigengene
    # sizeGrWindow(8,7);
    MEs<-net$MEs
    plots <- list()  # new empty list
    # dpa=as.factor(rep(seq(10,40,by=10),each=3))
    ss<-as.factor(coldata[[set]]$sample)
    for(me in 0:(Nmodules-1)) {
        which.module=paste("ME",me,sep="")
        module.color=labels2colors(me)
        #heatmap
        par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
        plotMat(t(scale(subDat[,net$colors==me ]) ),
        nrgcols=30,rlabels=T,rcols=module.color,
        main=paste(which.module, module.color, sep=": "), cex.main=2)
        #barplot
        par(mar=c(5, 4.2, 0, 0.7))
        barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
        ylab="eigengene expression", names.arg=paste(coldata[[set]]$domestication,coldata[[set]]$dpa), las=2)
        #line, anova
        df<-data.frame(ME=MEs[,which.module], ss, module = which.module, coldata[[set]] )
        fit<-aov(ME~domestication,df)  # basically looking at differences between wild and domesticated on all stages
        dfc<-summarySE(df, measurevar="ME", groupvars=c("dpa", "domestication"))
        plots[[me+1]]<- ggplot(dfc, aes(x=dpa, y=ME, fill = domestication)) +
        geom_bar(stat="identity",position=position_dodge(), color="black", size=0.3) +
        geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.3,position=position_dodge(0.9)) +
        ggtitle(paste(which.module," ",module.color,", P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
        theme_bw() +
        theme(plot.title=element_text( size=11),legend.position = "none")
    }
    for(page in 1:ceiling(Nmodules/9))
    {
        if(Nmodules>(9*page))
        {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
        else
        {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    }
    
    plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
    
    dev.off()# end the big plot
}


# focus on AD1 netwrok module 5

# Q: What other genes were also found in this module
table(AD1net12$colors==5)
# FALSE  TRUE
# 34072  1508
#1508 genes in this module

# examine connectivity with this module
MM=signedKME(multiExpr[[1]]$data, AD1net12$MEs)
order(MM[genes,"kME5"], decreasing=TRUE)[genes==profilin[1]]
# ranked at 122


# Q: what GO terms are enriched
library(topGO)
load('~/Dropbox/Scripts/D5annotation.Rdata')
# "annotation221" "annot"
load('~/Dropbox/Scripts/cottonGOenrich.RData')
# all included in individual network
universe<-colnames(multiExpr[[1]]$data)
# genes of interest
genes<-universe[AD1net12$colors==5]
geneList <- factor(as.integer(universe %in% genes))
names(geneList) <- universe
pdf("s4.topGO_AD1_ME5.pdf")
for(on in c("MF","BP","CC"))
{
    print(on)
    # Make topGO object
    GOdata <- new("topGOdata", ontology = on, allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    # fisher test
    result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    results.table <- GenTable(GOdata, result, topNodes = length(result@score))
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
    results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
    # label ontology type
    results.table$ontology<-on
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
    keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
    if(exists("enrich")) enrich<- rbind(enrich, keep)
    if(!exists("enrich")) enrich<- keep
    
    # draw figure for GO terms pval<=0.05 before FDR correction
    if(is.na(sigNo<-length(keep$ontology))){next}
    showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
    mtext(on, line=-1)
}
dev.off()
enrich
# [1] GO.ID       Term        Annotated   Significant Expected    result1     qval.bh     ontology
# <0 rows> (or 0-length row.names)
# no significant enrichment


# Compare network with DE
data = read.table("AD1.fiber.counts",header=TRUE, sep="\t");
# Take a quick look at what is in the data set:
dim(data);  #37223   25
names(data);
names(data)<-gsub("AD1_|.sort.bam","",names(data))
names(data)
#  Make each row corresponds to a gene and column to a sample or auxiliary information.
datExprT = as.data.frame(t(data[, -1]));
names(datExprT) = data[,1];
ADdata<-datExprT
ADdpa<-as.numeric(gsub(".*_|dpa","",rownames(ADdata)))

# normalization using DESeq2 rlog
count<-as.data.frame( t(datExprT) )
coldata<-data.frame( sample = rownames(ADdata), accession = gsub("_.*","", rownames(ADdata)), dpa = gsub(".*_","", rownames(ADdata)) )
coldata<-merge(coldata,info, all.x=TRUE,by="accession")
library(DESeq2)
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ domestication + dpa)
dds<-DESeq(dds)
resultsNames(dds)
# "Intercept"                 "domesticationDomesticated" "domesticationWild"         "dpa10dpa"
# "dpa15dpa"                  "dpa20dpa"                  "dpa5dpa"

# how many DEs between wild and domesticated
summary( (res=results(dds, contrast=c("domestication", "Domesticated", "Wild")) ), alpha=0.05 )
# out of 36439 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 3811, 10%
# LFC < 0 (down)   : 4521, 12%
# outliers [1]     : 800, 2.2%
# low counts [2]   : 5285, 15%
# (mean count < 1.5)
rownames(res)<-rownames(count)
write.table(res, "s4.DE_domestication.txt",row.names=TRUE,sep="\t")
DE_genes<-rownames(res)[res$padj<0.05 & res$log2FoldChange>0]
length(DE_genes<-DE_genes[!is.na(DE_genes)] )

# Are DE genes enriched in any modules?
probes<-colnames(multiExpr[[1]]$data)
asDEs<-probes %in% DE_genes
moduleDEs<-as.data.frame(table(AD1net12$colors[asDEs]))
names(moduleDEs)<-c("moduleLabels","DEs")
moduleAll <-as.data.frame(table(AD1net12$colors))
names(moduleAll)<-c("moduleLabels","All")
moduleDEs<-merge(moduleDEs, moduleAll, by="moduleLabels" ,all.y=TRUE)
moduleDEs$DEs[is.na(moduleDEs$DEs)]<-0
## calculate enrichment
tt<-colSums(moduleDEs[,2:3])
#  DEs   All
#  3788 35580  3788 DE genes were included in network
moduleDEs$fisherP<-round( apply(moduleDEs[,2:3],1,
function(x)fisher.test(matrix(as.numeric(c(x,tt-x)), nrow = 2, dimnames = list( c( "DEs","all"),c("inModule", "out"))) ,alternative="greater" )$p.value) ,3)
# CCs enriched in below modules
moduleDEs[moduleDEs$fisherP<0.05,]
# moduleLabels  DEs  All fisherP
# 3             2  906 5963       0
# 6             5 1062 1508       0
# 7             6  570 1460       0
# 9             8  467  857       0
# 10            9  142  854       0
# 12           12  247  636       0
# 15           21   49  233       0
# 19           34   35   95       0



#putative prf1 interactors
prf1<-"Gorai.009G028500"
prf1_interactor<-c("Gorai.005G187300","Gorai.012G001900","Gorai.006G236200","Gorai.005G001300","Gorai.004G221600","Gorai.001G254600","Gorai.006G044300","Gorai.004G153600","Gorai.008G122800","Gorai.011G016600","Gorai.006G167800","Gorai.009G418900","Gorai.009G297000","Gorai.012G017200","Gorai.007G355300","Gorai.011G213600","Gorai.006G017200","Gorai.012G054300")
#where are these interactors in network
AD1net12$colors[match(prf1_interactor,colnames(multiExpr[[1]]$data))]
#  3  5  2 16  0 47 NA  1  5 30  3  3  3  3  9  1  3 19  1
# only Gorai.004G153600 also in ME5
# are they also differentially expressed

#how many considered in DE analsys
table(rownames(count) %in% prf1_interactor)
# FALSE  TRUE
# 37205    18
# how many are they also differentially expressed as prf1
table(DE_genes %in% prf1_interactor)
# FALSE  TRUE
#  3810     1
# "Gorai.004G153600" also differentially expressed



# JJ's curated list of cotton genes related to Cytoskeleton and Cell wall functions in fibers, as CCs
list<- read.xls("genelist_cytoskeleton&cellwall.xlsx")
dim(list)
# locate CCs in modules
probes<-colnames(multiExpr[[1]]$data)
asCCs<-probes %in% list$Gene.ID
moduleCCs<-as.data.frame(table(AD1net12$colors[asCCs]))
names(moduleCCs)<-c("moduleLabels","CCs")
moduleAll <-as.data.frame(table(AD1net12$colors))
names(moduleAll)<-c("moduleLabels","All")
moduleCCs<-merge(moduleCCs, moduleAll, by="moduleLabels" ,all.y=TRUE)
moduleCCs$CCs[is.na(moduleCCs$CCs)]<-0
## calculate enrichment
tt<-colSums(moduleCCs[,2:3])
#  CCs   All
# 506 35580   506 CC genes were included in network
moduleCCs$fisherP<-round( apply(moduleCCs[,2:3],1,
function(x)fisher.test(matrix(as.numeric(c(x,tt-x)), nrow = 2, dimnames = list( c( "CCs","all"),c("inModule", "out"))) ,alternative="greater" )$p.value) ,3)
# CCs enriched in below modules
moduleCCs[moduleCCs$fisherP<0.05,]
# only module 19,11,9 are enriched with JJ's CC genes
# for ME5, 16 CC genes were found.

#how many CCs considered in DE analsys
table(rownames(count) %in% list$Gene.ID)
# FALSE  TRUE
# 36699   524
# how many are they also differentially expressed
table(DE_genes %in% list$Gene.ID)
# FALSE  TRUE
# 3769    42
mm<-matrix(as.numeric(c(42,3811,(524-42),(37223-3811))), nrow = 2, dimnames = list( c( "CCs","all"),c("inDE", "out")))
fisher.test(mm ,alternative="greater" )$p.value # 0.962885, not enriched


## Transcriptions factors
TF<-read.table("~/Dropbox/Literature/Network Analysis/GraP/Gene families/Family_TF.txt", header=TRUE,sep="\t")
TFs<-TF$Gene_ID
length(TFs<-unique(TFs))
# 3275
# how many TF in networks
table(probes %in% TFs)
# FALSE  TRUE
# 32454  3126
table(probes[AD1net12$colors==5] %in% TFs)
# FALSE  TRUE
# 1437    71
mm<-matrix(as.numeric(c(71,1508,(3126-71),(35580-1508))), nrow = 2, dimnames = list( c( "TFs","all"),c("inDE", "out")))
fisher.test(mm ,alternative="greater" )$p.value # 1, not enriched
#how many TFs considered in DE analsys
table(rownames(count) %in% TFs)
# FALSE  TRUE
# 33964  3259
# how many are they also differentially expressed
table(DE_genes %in% TFs)
# FALSE  TRUE
# 3644   167
mm<-matrix(as.numeric(c(167,3811,(3259-167),(37223-3811))), nrow = 2, dimnames = list( c( "TFs","all"),c("inDE", "out")))
fisher.test(mm ,alternative="greater" )$p.value # 1, not enriched

###################################################
# more pairwise DE analysis
coldata$sample0<-paste(coldata$domestication,coldata$dpa, sep=".")
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample0)
# pairwise deseq workflow
batch<- rbind(
c("Wild.5dpa", "Wild.10dpa" ),
c("Wild.10dpa", "Wild.15dpa" ),
c("Wild.15dpa", "Wild.20dpa" ),

c("Domesticated.5dpa", "Domesticated.10dpa" ),
c("Domesticated.10dpa", "Domesticated.15dpa" ),
c("Domesticated.15dpa", "Domesticated.20dpa" ),

c("Wild.5dpa",  "Domesticated.5dpa" ),
c("Wild.10dpa", "Domesticated.10dpa" ),
c("Wild.15dpa", "Domesticated.15dpa" ),
c("Wild.20dpa", "Domesticated.20dpa" )
 )

pairwiseDE<-function(dds, contrast,savePath)
{
    # DE analysis
    print(contrast)
    ddsPW <-dds[,dds$sample0 %in% contrast]
    ddsPW$sample0<-droplevels(ddsPW$sample0)
    res <- results(DESeq(ddsPW))
    print( summary(res,alpha=.05) ) # print results
     write.table(res, file=paste(savePath,"DEs/",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
}
apply(batch,1,function(x) pairwiseDE(dds,x,savePath = ""))
# [1] "Wild.5dpa"  "Wild.10dpa"
# out of 34797 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 409, 1.2%
# LFC < 0 (down)   : 1772, 5.1%
# outliers [1]     : 1086, 3.1%
# low counts [2]   : 6828, 20%
# (mean count < 3.9)
#
# [1] "Wild.10dpa" "Wild.15dpa"
# out of 35108 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 1, 0.0028%
# LFC < 0 (down)   : 1, 0.0028%
# outliers [1]     : 1156, 3.3%
# low counts [2]   : 0, 0%
# (mean count < 0.1)
#
# [1] "Wild.15dpa" "Wild.20dpa"
# out of 34980 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 1509, 4.3%
# LFC < 0 (down)   : 1448, 4.1%
# outliers [1]     : 670, 1.9%
# low counts [2]   : 8621, 25%
# (mean count < 7.4)
#
# [1] "Domesticated.5dpa"  "Domesticated.10dpa"
# out of 34700 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 656, 1.9%
# LFC < 0 (down)   : 1124, 3.2%
# outliers [1]     : 275, 0.79%
# low counts [2]   : 8546, 25%
# (mean count < 8.1)
#
# [1] "Domesticated.10dpa" "Domesticated.15dpa"
# out of 34775 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 1182, 3.4%
# LFC < 0 (down)   : 1026, 3%
# outliers [1]     : 405, 1.2%
# low counts [2]   : 8527, 25%
# (mean count < 9.5)
#
# [1] "Domesticated.15dpa" "Domesticated.20dpa"
# out of 34656 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 3316, 9.6%
# LFC < 0 (down)   : 2657, 7.7%
# outliers [1]     : 457, 1.3%
# low counts [2]   : 6839, 20%
# (mean count < 3.7)
#
# [1] "Wild.5dpa"         "Domesticated.5dpa"
# out of 34459 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 711, 2.1%
# LFC < 0 (down)   : 425, 1.2%
# outliers [1]     : 280, 0.81%
# low counts [2]   : 8509, 25%
# (mean count < 6.5)
#
# [1] "Wild.10dpa"         "Domesticated.10dpa"
# out of 35001 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 427, 1.2%
# LFC < 0 (down)   : 374, 1.1%
# outliers [1]     : 1104, 3.2%
# low counts [2]   : 10165, 29%
# (mean count < 15.5)
#
# [1] "Wild.15dpa"         "Domesticated.15dpa"
# out of 35094 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 392, 1.1%
# LFC < 0 (down)   : 356, 1%
# outliers [1]     : 739, 2.1%
# low counts [2]   : 8580, 24%
# (mean count < 10.3)
#
# [1] "Wild.20dpa"         "Domesticated.20dpa"
# out of 34357 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 444, 1.3%
# LFC < 0 (down)   : 380, 1.1%
# outliers [1]     : 619, 1.8%
# low counts [2]   : 11849, 34%
# (mean count < 12)

profilin<- c( "Gorai.009G028500", "Gorai.009G028400", "Gorai.004G282400", "Gorai.013G229500", "Gorai.010G078400", "Gorai.001G025300", "Gorai.003G061200", "Gorai.010G078500")
DEfiles<-list.files("DEs/")
out<-data.frame()
for(file in DEfiles)
{
    x<-read.table(paste("DEs/", file, sep=""),header=TRUE, sep="\t")
    y<-x[profilin,]
    y$comparison<-file
    out<-rbind(out,y)
    # print(file)
    # print(y)
}

####################################################
######################################################
# Homoeolog ratio analysis
x<-read.table("all.AD1.fiber.counts",sep="\t",header=TRUE)
rownames(x)<-x$gene
x<-x[,-1]
y <- log2(x+1)
z<-y[,seq(1,47, by=2)] -y[,seq(2,48, by=2)]  # log(At/Dt)
gsub("..bam","",names(x)[seq(1,47, by=2)]) == gsub("..bam","",names(x)[seq(2,48, by=2)])
boxplot(z)  # log2 ratio center at zero for all samples, looking good
apply(z,2,quantile)
apply(z,2,function(x)quantile(x,c(0.05, 0.95)) )
names(z)<-gsub(".sort.A.bam|AD1_","",names(z))

info<-read.table("accesionInfo.txt",sep="\t",header=TRUE)
sample<-data.frame(name=names(z),accession=gsub("_.*","",names(z) ), dpa=gsub(".*_","",names(z) ) )
sample<-merge(sample,info,by="accession")

### how to determine homoeolog expression bias: At!=Dt
# 1. student t test on log2(ratio)!=0, with bonferroni correction
n=24
stat<-data.frame(mean=rowMeans(z), sd=apply(z,1,sd) )
stat$t<-stat$mean/(stat$sd/sqrt(n))
stat$p <- 2*pt(-abs(stat$t),df=n-1)
stat$p.bh<-p.adjust(stat$p, method="BH")
stat$p.bonferroni<-p.adjust(stat$p, method="bonferroni")
# valcano plot
plot(stat$mean,log2(stat$p ))
abline(h=log2(0.05))  # pvalue
abline(v=2.8)
abline(v=-2.8)

# 2. DESeq comparing At vs Dt read counts
library(DESeq2)
# get library factors from total counts
data = read.table("AD1.fiber.counts",header=TRUE, sep="\t");
names(data);
names(data)<-gsub("AD1_|.sort.bam","",names(data))
names(data)
count<-data[,-1]
rownames(count)<-data$gene
coldata <- data.frame(sample = names(count),accession = gsub("_.*","", names(count)))
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
dds<-estimateSizeFactors(dds)
sizeFactors(dds)
totalSizeFactors<-data.frame(sampleS=names(count),sf=sizeFactors(dds) )

# Take a quick look at what is in the data set:
dim(x);  #37223   25
names(x)<-gsub("AD1_|.sort|.bam","",names(x))
count<-x
coldata<-data.frame( sample = names(x), accession = gsub("_.*","", names(x)), dpa = gsub(".*_|dpa[.].","", names(x)) )
coldata<-merge(coldata,info, all.x=TRUE,by="accession")
coldata$homoeolog<-gsub(".*[.]","",coldata$sample)
coldata$sampleS<-gsub("[.].*","",coldata$sample)
coldata$sampleH<-paste(coldata$domestication,coldata$dpa,sep="." )
coldata<-merge(coldata, totalSizeFactors, by="sampleS")

# homoeolog foe wild and domesticated each dpa tests
for( i in c("Domesticated","Wild") )
{
    select<-which(coldata$domestication==i)
    dds <- DESeqDataSetFromMatrix( countData = count[,select], colData = coldata[select,], design = ~ homoeolog + dpa)
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    sizeFactors(dds)<-coldata$sf[select]
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds, c("homoeolog", "A","D"))
    print(i)
    print(res)
    summary(res, alpha=0.05)
    assign(i,res)
}
# "Domesticated.10" "Domesticated.15" "Domesticated.20" "Domesticated.5"
# "Wild.10"  "Wild.15"         "Wild.20"         "Wild.5"
# "Domesticated"27-27% "Wild" 23-23%
# Homoeolog bias is always balanced in fiber as shown here. Very interesting, seed transcriptome shows unbalance towards A bias


# pairwise at each dpa tests
for( i in unique(coldata$sampleH ) )
{
    select<-which(coldata$sampleH==i)
    dds <- DESeqDataSetFromMatrix( countData = count[,select], colData = coldata[select,], design = ~ homoeolog)
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    sizeFactors(dds)<-coldata$sf[select]
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)
    res <- results(dds)
    print(i)
    print(res)
    summary(res, alpha=0.05)
    assign(i,res)
}
# "Domesticated.10" "Domesticated.15" "Domesticated.20" "Domesticated.5"
# "Wild.10"  "Wild.15"         "Wild.20"         "Wild.5"

homoeo <- rbind()


load('R-04-buildNetwork.RData')
load('R-02-dataInput.RData')
probes<-colnames(multiExpr[[1]]$data)
net<-AD1net12
module<- net$colors
sampleH<-unique(coldata$sampleH )

homoeo<-data.frame(nodeName=rownames(Domesticated.10) )
for( j in sampleH )
{
    res<-get(j)
    tempt<-res[,c("log2FoldChange","padj")]
    tempt$log2FoldChange[tempt$padj>0.05 | is.na(tempt$padj)] <-0
    names(tempt)<-paste(j, names(tempt),sep="_" )
    homoeo<-cbind(homoeo, tempt)
    # Dt <- rownames(res)[res$log2FoldChange>log2(1.5) & res$padj<0.05 & !is.na(res$padj)]
    # At <- rownames(res)[res$log2FoldChange<log2(2/3) & res$padj<0.05 & !is.na(res$padj)]
    # print(j)
    # print("A bias")
    # print( table(module[At] ) )
    # print("D bias")
    # print( table(module[Dt] ) )
}
write.table(homoeo,"homoeologExpressionBias.txt",sep="\t",row.names=FALSE, quote=FALSE)


###################################
####### Output for cytoscape ######
load('R-04-buildNetwork.RData')
load('R-02-dataInput.RData')
probes<-colnames(multiExpr[[1]]$data)

probes = colnames(multiExpr[[2]]$data)
# locate genes in Module 5 only
asME5<-which(AD1net12$colors==5)
length(asME5)# 1508

# Export FA network to cytoscape
# connect to speedy, or whatever used to build network and store the TOM file (ISOLONE, biocrunch, etc.)
# Cytoscape [2] allows the user to input an edge file and a node file, allowing the user to specify for example the link weights and the node colors.
# Get topological overlap, "TOM"
load("AD1_power12_TOM-block.1.RData")
# Select the corresponding Topological Overlap
subTOM = as.matrix(TOM)[asME5, asME5];
# str(subTOM)
subProbes = probes[asME5];
dimnames(subTOM) = list(subProbes, subProbes)

# gather node infor
#homoeolog expression bias
homoeo<-read.table("homoeologExpressionBias.txt",sep="\t", header=TRUE)
ME5<-homoeo[homoeo$nodeName %in% subProbes,]
# transcription factors
TF<-read.table("Family_TF.txt", header=TRUE,sep="\t")
TF<-TF[,-1]
dim(TF<-unique(TF)) #3282
names(TF)[1]<-"nodeName"
# JJ's list
CC<-read.table("genelist_cytoskeleton&cellwall.txt",sep="\t",header=TRUE)
dim( CC<-unique(CC[,-2])  )
names(CC)<-c("nodeName","Family")
CC$Terms<-"Cell wall & cytoskeleton"
# because no overlap between CC and TF
intersect(TF$nodeName,CC$nodeName)
func<-rbind(TF,CC)
dim(ME5<-merge(ME5,func, by="nodeName", all.x=TRUE )  )
table(ME5$nodeName ==subProbes)


# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(subTOM,
    edgeFile = "AD1-edges-ME5.txt",
    nodeFile = "AD1-nodes-ME5.txt",
weighted = TRUE, threshold = quantile(subTOM, c(0.9)),  # change quantile to get 0.75
    nodeNames = subProbes,
    nodeAttr = ME5[,-1]   )

top10_neighbours<-c("Gorai.001G128800", "Gorai.002G121400", "Gorai.003G112600", "Gorai.003G116500", "Gorai.004G293100", "Gorai.005G073200", "Gorai.009G028500", "Gorai.011G000100", "Gorai.011G230600", "Gorai.011G277300", "Gorai.012G118400")


moduleLabels<-AD1net12$colors
names(moduleLabels)<-probes
homoeo$module<-moduleLabels[homoeo$nodeName]
for(i in names(homoeo)[c(2,4,6,8,10,12,14,16)] )
{   bias<-homoeo[,c("nodeName","module",i)]
    cutoff<-log2(1)
    all<-as.data.frame(table(bias$module))
    At<-as.data.frame(table(bias$module[bias[,i]>  cutoff] ))
    Dt<-as.data.frame(table(bias$module[bias[,i]< -cutoff] ))
    t<-merge(merge(all,At,by="Var1",all.x=TRUE,all.y=TRUE),Dt,by="Var1",all.x=TRUE,all.y=TRUE)
    names(t)<-c("module","all","At","Dt")
    t[is.na(t)]<-0
    t$towards<-ifelse(t$At==t$Dt ,"-",ifelse(t$At>t$Dt,"At","Dt"))
    t$pval<-apply(t[,3:4],1,function(x) {ee<-(x[1]+x[1])/2; tail<-ifelse(x[1]>x[2],"greater","less"); fisher.test(matrix(c(x[1],x[2],ee,ee), nrow = 2, dimnames = list( c( "Observed","Expected"),c("At", "Dt"))), alternative=tail )$p.value })
    table(t$padj<-p.adjust(t$pval,method="bonferroni")<0.05 )
    print(i)
    print(t)
    print( xtabs(padj~towards, data=t) )
}
# basically all unbalanced towards At bias



########### Export network for JJ's ##########################
# JJ's list
CC<-read.table("genelist_cytoskeleton&cellwall.txt",sep="\t",header=TRUE)
dim( CC<-unique(CC[,-2])  )  #585
aggregate(CC$Function, list(CC$Gene.ID), function(x)paste(unique(x),sep="/")) ->CC531
names(CC531)<-c("nodeName","Family")
CC531$Terms<-"Cell wall & cytoskeleton"
asCC<-which(probes %in% CC531$nodeName)
length(asCC)# 506


# get domestication
# 2. DESeq comparing At vs Dt read counts
library(DESeq2)
# get library factors from total counts
data = read.table("AD1.fiber.counts",header=TRUE, sep="\t");
names(data);
names(data)<-gsub("AD1_|.sort.bam","",names(data))
names(data)
count<-data[,-1]
rownames(count)<-data$gene
coldata <- data.frame(sample = names(count),accession = gsub("_.*","", names(count)))
info<-read.table("accesionInfo.txt",header=TRUE,sep="\t")
coldata<-merge(coldata,info,merge="accesion")
coldata$dpa<-gsub(".*_|dpa","",coldata$sample)
# deseq
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ domestication + dpa)
dds<-DESeq(dds)
res <- results(dds,c("domestication","Domesticated","Wild") )
summary(res, alpha=0.05)
AD1<-res
### about 23%
#######################
data = read.table("Adip.fiber.counts",header=TRUE, sep="\t");
names(data);
names(data)<-gsub("A1_|.sort.bam","",names(data))
names(data)
count<-data[,-1]
rownames(count)<-data$gene
coldata <- data.frame(sample = names(count),accession = gsub("_.*","", names(count)))
info<-read.table("accesionInfo.txt",header=TRUE,sep="\t")
coldata<-merge(coldata,info,merge="accesion")
coldata$dpa<-gsub(".*_|dpa.*","",coldata$sample)
# deseq
dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ domestication + dpa)
dds<-DESeq(dds)
res <- results(dds,c("domestication","Domesticated","Wild") )
summary(res, alpha=0.05)
A1<-res
### ??? only 2.5%



# Export CC network to cytoscape
# connect to speedy, or whatever used to build network and store the TOM file (ISOLONE, biocrunch, etc.)
# Cytoscape [2] allows the user to input an edge file and a node file, allowing the user to specify for example the link weights and the node colors.
# Get topological overlap, "TOM"
load("AD1_power12_TOM-block.1.RData")
# Select the corresponding Topological Overlap
subTOM = as.matrix(TOM)[asCC, asCC];
# str(subTOM)
subProbes = probes[asCC];
dimnames(subTOM) = list(subProbes, subProbes)

AD1$log2FoldChange[AD1$padj>0.05] =0
AD1$nodeName <-rownames(AD1)
node<-merge(CC531, as.data.frame(AD1[,c("log2FoldChange","nodeName")]),by="nodeName", all.x=TRUE)
node0<- node[na.omit(match(subProbes,node$nodeName)),]
table(node0$nodeName==subProbes)
node0<-as.data.frame(apply(node0,2,as.character))
node0$module<-labels2colors(AD1net12$colors[asCC])

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(subTOM,
edgeFile = "AD1-edges-CC.txt",
nodeFile = "AD1-nodes-CC.txt",
weighted = TRUE, threshold = quantile(subTOM, c(0.9)),  # change quantile to get 0.75
nodeNames = subProbes,
nodeAttr = node0[,-1]   )

