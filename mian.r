########################################expression data and network data
net <- read.csv(file=".../background.network.txt", head=F, sep="\t") #background network data
pathway <- load(file=".../reactom/reactom.pathway.total.names.Rdata") #pathway data from reactom
GTEx_290 <- read.table(file=".../GTEx_breast_TPM.txt", head=T, sep="\t") #expression matrix of normal samples
TCGA_1093 <- read.table(file=".../TCGA_breast_TPM.txt", head=T, sep="\t") #expression matrix of cancer samples


########################################conert expression matrix to delta rank matrix
genes <- intersect(row.names(GTEx_290), row.names(TCGA_1093))
TCGA_1093 <- TCGA_1093[which(row.names(TCGA_1093) %in% genes ), ]
GTEx_290 <- GTEx_290[which(row.names(GTEx_290) %in% genes), ]
breast_300 <- cbind(GTEx_290, apply(GTEx_290, 1, mean))

rank_300 <- rank.matrix(breast_300)
rank_TCGA <- rank.matrix(TCGA_1093)
n = 290
x = cbind(rank_300, rank_TCGA)
dim(x)

n_normal = 290
n_cancer = 1093
deltarank.result <- delta.rank(net, x, n_normal, n_cancer)
save(deltarank.result, file=".../deltarank.result.Rdata")

####################################### caculate the edge-perturbation matrix
EPm_normal <- EPm (deltarank.result[[2]], deltarank.result[[3]])
EPm_cancer <- EPm (deltarank.result[[2]], deltarank.result[[4]])  

save(EPm_normal, file="EPm_normal.Rdata")
save(EPm_cancer, file="EPm_cancer.Rdata")

########################################################### select features (Sd+ kw test)
data <- t(cbind(EPm_cancer, EPm_normal))
group <- as.factor(c(rep("cancer", 1093), rep("normal", 290)))

kw.test <- function(x){
  data<- as.data.frame(cbind(x, group))
  colnames(data) <- c("value", "group")
  p<- kruskal.test(value ~ group, data =data)$p.value
  return(p)
}
P.value <- apply(data, 2, kw.test)
f1 <- order(P.value)[1:30000]

sd = apply(EPm_cancer, 1, sd)   
f2 =rev(order(sd))[1:30000]

fea.loc <- intersect(f1, f2) #location of feature
feature <- EPm_cancer[fea.loc,] #feature matrix used for clustering

####################################################consensus cluster
library(ConsensusClusterPlus)
results <- ConsensusClusterPlus(
           feature, maxK = 6, reps=10, pItem=0.8, pFeature=0.8, clusterAlg="km",title="untitled_consensus_cluster",
           innerLinkage="average", finalLinkage="average", distance="pearson", ml=NULL,
           tmyPal=NULL,seed=NULL,plot=NULL,writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,verbose=F)





































#####################################################
setwd("1617net")
load(file="deltarank.TCGA.reactom.Rdata")
load(file="stablenet.edge.reactom.Rdata")
load(file="stablenet.data.reactom.Rdata")
load(file="deltarank.290.reactom.Rdata")

SDR.cancer <- SDRm(as.numeric(stablenet.data), deltarank.TCGA.reactom)
SDR.normal <- SDRm(as.numeric(stablenet.data), deltarank.290.reactom)

data <- t(cbind(SDR.cancer, SDR.normal))
group <- as.factor(c(rep("cancer", 1093), rep("normal", 290)))
 
kw.test <- function(x){
        data<- as.data.frame(cbind(x, group))
         colnames(data) <- c("value", "group")
         p<- kruskal.test(value ~ group, data =data)$p.value
         return(p)
     }
P.value <- apply(data, 2, kw.test)
f1 <- order(P.value)[1:30000]
sd = apply(SDR.cancer, 1, sd)   
f2 =rev(order(sd))[1:30000]
fea.loc <- intersect(f1, f2) #1911 features

SDR.cancer.1911 <- SDR.cancer[fea.loc, ]
SDR.cancer.1911.nor <- apply(SDR.cancer.1911, 2, scale)
pheatmap(SDR.cancer.1911.nor, show_colnames = F)



