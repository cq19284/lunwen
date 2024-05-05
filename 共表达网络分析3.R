
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
datTraits <- read.csv("new.csv",header=TRUE)

datTraits <- datTraits[,-1]
datTraits
sr = as.data.frame(datTraits)
# 明确基因和样本数
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# 用颜色标签重新计算MEs
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# 指定datTrait中感兴趣的一个性状，这里选择sr_g


names(sr) = "sr"



moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 查看每个模块的颜色及包含的基因数目
table(moduleColors)


#pdf("Module-trait associations.pdf",width = 8, height=10)
# 通过p值显示相关性
#textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
#dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3))
# 通过热图显示相关性
#labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(sr), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
# 想改变热图配色可以更改 colors = greenWhiteRed(50)
#dev.off()

# 指定datTrait中感兴趣的一个性状，这里选择sr_g


#  各基因模块的名字（颜色）
modNames = substring(names(MEs), 3)

# 计算MM的P值
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership 
), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# 计算性状和基因表达量之间的相关性（GS）
geneTraitSignificance = as.data.frame(cor(datExpr, sr, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                            nSamples))
names(geneTraitSignificance) = paste("GS.", names(sr), sep="")
names(GSPvalue) = paste("p.GS.", names(sr), sep="")

module = "violet"
column = match(module, modNames)
moduleGenes = moduleColors==module
pdf("Module membership vs gene significance.pdf",width = 7, height=7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, 
                                              column]),abs(geneTraitSignificance[moduleGenes, 1]), xlab = 
                       paste("Module Membership in", module, "module"), ylab = "Gene 
  significance for body sr", main = paste("Module membership 
  vs gene significance"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 
                       1.2, col = 'black')
dev.off()