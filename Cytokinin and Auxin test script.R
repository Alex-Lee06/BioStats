library(limma)
library(affy)

#uses the created text file to orginize data
pd = read.AnnotatedDataFrame("sample.txt",header=TRUE,sep=",",row.names=1)
mydata = ReadAffy(filenames=pd$filename,phenoData=pd,verbose = TRUE)
sampleNames(mydata)=row.names(pData(pd))
eset = rma(mydata)
expression_data = exprs(eset)

TS = factor(pd$treatment, levels = c("E", "A", "B"))
#TS = factor(c("E","E", "E", "A", "A", "A", "B", "B", "B"))

exp_design = model.matrix(~ 0 +TS)
colnames(exp_design)=levels(TS)

#linear model of the expression data.
fit=lmFit(eset,exp_design)

#compares the treatments to each other
contMatrix= makeContrasts(AvE = A - E, BvA = B - A, BvE = B - E, levels = exp_design)

#creates a matirx for the compared treatments
fit2 = contrasts.fit(fit, contMatrix)

diff_exp = eBayes(fit2)

#creates a venn diagram of the top genes with specified p.value
result = decideTests(diff_exp, p.value = 0.01)
vennDiagram(result)

#uses the affymetrix information from affy themselves
f = "ftp://ftp.arabidopsis.org/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt"
annots = read.delim(f, na.strings ="", fill=TRUE, header = T, sep="\t")
annots = annots[,c(1,5,6)]

#shows top data with the p.value of specified value. 
N = dim(eset)[1]
top.A.AvE = toptable(diff_exp, coef=1, number = N, p.value = 0.001)
top.B.BvA = toptable(diff_exp, coef=2, number = N, p.value = 0.001)
top.E.BvE = toptable(diff_exp, coef = 3, number = N, p.value = 0.05)

#shows the gene names of compared treatments.
row.names(annots) = annots$array_element_name
top.B.BvA.merged = merge(annots, top.B.BvA, by = "row.names")
top.A.AvE.merged = merge(annots, top.A.AvE, by = "row.names")
top.E.BvE.merged = merge(annots, top.E.BvE, by = "row.names")

#shows all genes that are compared and have specified p.value from above.
genes.Diff = unique(top.B.BvA.merged$locus)
genes.E.BvE = unique(top.E.BvE.merged$locus)
genes.A.AvE = unique(top.A.AvE.merged$locus)

#shows the log fold changes of each compared treatment
logFC = top.E.BvE.merged[,'logFC']
logFC_AvE = top.A.AvE.merged[,'logFC']
logFC_BvE = top.B.BvA.merged[,'logFC']

#creates a histogram of the logFC for top.E.BvE.merged
hist(logFC)

#scatter plot for logFC of each treatment. Not sure if you need this though
plot(logFC)
plot(logFC_AvE)
plot(logFC_BvE)
