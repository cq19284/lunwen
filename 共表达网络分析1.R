library('WGCNA')
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
# 读取CSV文件\
fpkm <- read.csv("D:\\rcode\\rpm429.csv")
dim(fpkm)
rownames(fpkm) <- fpkm[,1]
fpkm<- fpkm[,-1]
WGCNA_matrix = t(fpkm[order(apply(fpkm,1,mad), decreasing = T)[1:42389],]) #按mad进行排序，取前30%的基因进行后续分析
datExpr0 = as.data.frame(WGCNA_matrix) 

#全部基因
datExpr_all = as.data.frame(t(fpkm))
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK


sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(120,20)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


# 画一条阈值线，根据实际情况人工设定，这里是15
abline(h = 230000, col = "red")
# 确定阈值线下的集群
clust = cutreeStatic(sampleTree, cutHeight = 230000, minSize = 10)
table(clust)
# clust 1包含想要留下的样本
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

datTraits <- read.csv("D:\\python_projects\\create_video_project\\pythonProject2\\new.csv",header=TRUE)
dim(datTraits)
rownames(datTraits) <- datTraits[,1]
rownames(datTraits)
datTraits <- datTraits[,-1]
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = "WGCNA0.3-dataInput.RData")