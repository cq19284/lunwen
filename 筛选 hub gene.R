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
## 这部分在专题（3）中计算过一次，因为没有保存，所以这里又重新计算
nSamples = nrow(datExpr)
# 指定datTrait中感兴趣的一个性状，这里选择sr_g

names(sr) = "sr"
# 各基因模块的名字（颜色）
modNames = substring(names(MEs), 3)

module = "violet"
column = match(module, modNames)
moduleGenes = moduleColors==module
lightyellow_module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes])
names(lightyellow_module)="genename"
lightyellow_KME<-as.data.frame(datKME[moduleGenes,column]) 
names(lightyellow_KME)="KME"
rownames(lightyellow_KME)=lightyellow_module$genename
FilterGenes = abs(lightyellow_KME$KME) 0.8
table(FilterGenes)

lightyellow_hub<-subset(lightyellow_KME, abs(lightyellow_KME$KME)>0.8)
write.csv(lightyellow_hub, "hubgene_KME_lightyellow.csv")
