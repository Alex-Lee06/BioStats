library(limma)
library(affy)
mydata = read.csv("GSE59741_series_matrix.txt")
pd = read.AnnotatedDataFrame("sample.txt",header=TRUE,sep=",",row.names=1)
