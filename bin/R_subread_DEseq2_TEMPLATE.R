#########################################################################################
#####Run Rsubread for counts
#install.packages("BiocManager")
#BiocManager::install("Rsubread")
#BiocManager::install("DESeq2")


library(Rsubread);


files=c(bamlist)

gtf=("temp_annotation")

###Sense Unique
senseunique=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext = gtf,
GTF.attrType = "BLANK",
GTF.featureType = "exon",
useMetaFeatures = METAFEATURE,
allowMultiOverlap = TRUE,
isPairedEnd=PAIRED,
nthreads = cpus,
strandSpecific = 1)

write.table(x=data.frame(
senseunique$annotation[,c("GeneID","Length")],
senseunique$counts,stringsAsFactors=FALSE),
file="outputprefix/Rsubread_sense.txt",
quote=FALSE,sep="\t",row.names=FALSE)


#########################################################################################
#########################################################################################
#######DESeq2

library(DESeq2)
#Read in the count matrix
cts <- as.matrix(read.csv("outputprefix/Rsubread_sense.txt",sep="\t",row.names=1))
coldata <- read.table("COL_DATA",row.names=1)
names(coldata) <- c("condition")


cts <- subset(cts, select=-c(Length)) #Get rid of length column


# print(cts)
# print(coldata)

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

#dds <- dds[ rowSums(counts(dds)) > 1,]   #Take everything over 1 because we added one in the normalization python steps
dds$condition <- factor(dds$condition, levels=c("CONTROL","TREATMENT"))
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, file="outputprefix/DESeq2_sense.csv")
#############################################
