#I imported the data using file import because it wasn't recognizing it as in the same folder for some reason. 
#Thus the code is not included. Just make sure to import the data you want to use as file and that row 1 is considered names.

DrugDesign = data.frame(row.names = colnames( data),
                       condition = c("untreated", "untreated", "untreated", "untreated", "untreated", "untreated", "untreated",
                                     "treated", "treated", "treated", "treated", "treated", 
                                     "untreated", "untreated", "untreated", "untreated", "untreated", "untreated", "untreated", 
                                     "treated", "treated", "treated", "treated", "treated"),
                       location = c("blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood",
                                    "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor")
                       )
DrugDesign

install.packages("BiocManager")
BiocManager::install("DESeq")
library("DESeq")

tumorSamples = DrugDesign$location =="tumor"
countTable = data[ , tumorSamples ]
condition = DrugDesign$condition[tumorSamples]

cds = newCountDataSet( countTable, condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

plotDispEsts(cds)

str(fitInfo(cds))

res = nbinomTest(cds, "untreated", "treated")
head(res)
plotMA(res)

plcdsFull = newCountDataSet(data, DrugDesign)
cdsFull = estimateSizeFactors(cdsFull)
cdsFull = estimateDispersions(cdsFull)

hist(res$pval, breaks = 100, col="skyblue", border="slateblue", main="")
resSig = res[res$padj<0.1, ]
head(resSig[order(resSig$pval),])

plotDispEsts(cdsFull)

fit1 = fitNbinomGLMs( cdsFull, count ~ location + condition)
fit0 = fitNbinomGLMs( cdsFull, count ~ location)

str(fit1)
str(fit0)

pvalsGlm = nbinomGLMTest(fit1, fit0)
padjGLM = p.adjust(pvalsGlm, method = "BH")

tab1 = table("paired-end only" = res$padj < .1, "all samples" = padjGlm < .1)
addmargin(tab1)