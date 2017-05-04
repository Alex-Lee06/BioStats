library(limma)
library(affy)
mydata = read.csv("GSE59741_series_matrix.txt")
pd = read.AnnotatedDataFrame("sample.txt",header=TRUE,sep=",",row.names=1)
mydata = ReadAffy(filenames=pd$filename,phenoData=pd,verbose = TRUE)
sampleNames(mydata)=row.names(pData(pd))
eset = rma(mydata)
expression_data = exprs(eset)

TS = paste(pd$treatment)
TS = factor(TS, levels = c("E", "A", "B"))

exp_design = model.matrix(~ 0 +TS)
colnames(exp_design)=levels(TS)

fit=lmFit(eset,exp_design)

f = "ftp://ftp.arabidopsis.org/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt"
annots = read.delim(f, na.strings ="", fill=TRUE, header = T, sep="\t")
annots = annots[,c(1,5,6)]

#(For now to see the gene names double click on annots in environment)
