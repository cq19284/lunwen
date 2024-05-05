library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#载入第一步中的表达量和表型值
lnames = load(file = "WGCNA0.3-dataInput.RData")
lnames

# 设置网络构建参数选择范围
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 
powers

# 计算无尺度分布拓扑矩阵
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
#一页多图，一行2列
par(mfrow = c(1,2))
#字号
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power（左图）
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
abline(h=0.90,col="red")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red") 

sft$powerEstimate
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 200, 
  deepSplit = 1,  reassignThreshold = 0, 
  mergeCutHeight = 0.3,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F, 
  verbose = 3
)

table(net$colors)
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-auto.RData")