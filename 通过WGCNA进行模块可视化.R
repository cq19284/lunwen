library('WGCNA')
enableWGCNAThreads()
setwd("D:\\rcode")
# 允许R语言以最大线程运行
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 载入第一步中的表达量和表型值
lnames = load(file = "WGCNA0.3-dataInput.RData")
lnames
# 载入第二步的网络数据
lnames = load(file = "networkConstruction-auto.RData")
lnames
lnames = load(file = "TOM_object.RData")
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# 模块检测时的计算，重新算一次
dissTOM = 1-TOM
# 转换dissTOM，方便在热图中显示

nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
# Open a graphical window
pdf("Network heatmap plot_selected genes.pdf", width = 9, height = 9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot_selected genes")
dev.off()