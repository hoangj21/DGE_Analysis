#I imported the data using file import because it wasn't recognizing it as in the same folder for some reason. 
#Thus the code is not included. Just make sure to import the data you want to use as file and that row 1 is considered names.

install.packages("BiocManager")
BiocManager::install("DESeq")
library("DESeq")

#Setting up metadata
DrugDesign = data.frame(row.names = colnames( data),
                       condition = c("untreated", "untreated", "untreated", "untreated", "untreated", "untreated", "untreated",
                                     "treated", "treated", "treated", "treated", "treated", 
                                     "untreated", "untreated", "untreated", "untreated", "untreated", "untreated", "untreated", 
                                     "treated", "treated", "treated", "treated", "treated"),
                       location = c("blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood", "blood",
                                    "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor", "tumor")
                       )
DrugDesign

#This is for a single parameter normalization currently
tumorSamples = DrugDesign$location =="tumor"
countTable = data[ , tumorSamples ]
condition = DrugDesign$condition[tumorSamples]

#This is only comparing untreated to treated in tumor cells. Comparing all four is later, also that doesn't work yet.
cds = newCountDataSet( countTable, condition)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

#plotted normalized data
plotDispEsts(cds)
#not quite sure how to interpret this but here it is.
str(fitInfo(cds))

#this is actually testing to find differential expression between treated and untreated tumor.
res = nbinomTest(cds, "untreated", "treated")
head(res)
#visualizations of the results
plotMA(res)
hist(res$pval, breaks = 100, col="skyblue", border="slateblue", main="")

#these are the statistically significant results in order of their lowest p-values. 
resSig = res[res$padj<0.1, ]
head(resSig[order(resSig$pval),])


#normalization with two factors
plcdsFull = newCountDataSet(data, DrugDesign)
cdsFull = estimateSizeFactors(cdsFull)
cdsFull = estimateDispersions(cdsFull)

plotDispEsts(cdsFull)

#!! This does not work yet, not sure why. Neither converges.
fit1 = fitNbinomGLMs( cdsFull, count ~ location + condition)
fit0 = fitNbinomGLMs( cdsFull, count ~ location)

str(fit1)
str(fit0)

pvalsGlm = nbinomGLMTest(fit1, fit0)
padjGLM = p.adjust(pvalsGlm, method = "BH")

tab1 = table("paired-end only" = res$padj < .1, "all samples" = padjGlm < .1)
addmargin(tab1)
