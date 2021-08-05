#########################################################################
# k-means gene clustering
# Travis S. Johnson PhD
#
# This script performs k-means clustering on the gene expression results

#########################################################################
# Loading requried packages and functions
library(pheatmap)
library(ggplot2)
library(DEGAS)
#library(cowplot)

findElbow <- function(Nk,Ek){
  names(Nk) = as.character(Nk)
  names(Ek) = as.character(Nk)
  NkErrors = rep(NA,length(Nk))
  names(NkErrors) = as.character(Nk)
  for (k in Nk){
    Nk1 = Nk[as.character(min(Nk):k)]
    Ek1 = Ek[as.character(min(Nk):k)]
    Nk2 = Nk[as.character((k+1):max(Nk))]
    Ek2 = Ek[as.character((k+1):max(Nk))]
    if(!any(is.na(Nk1)) & length(Nk1)>2){
      linfit1 = lm(Ek1~Nk1)
      linfit1res = sum(abs(linfit1$residuals))
    }else{
      linfit1res = 0
    }
    if(!any(is.na(Nk2)) & length(Nk2)>2){
      linfit2 = lm(Ek2~Nk2)
      linfit2res = sum(abs(linfit2$residuals))
    }else{
      linfit2res = 0
    }
    NkErrors[as.character(k)] = sum(linfit1res,linfit2res)
  }
  return(as.numeric(names(NkErrors)[NkErrors==min(NkErrors)]))
}

#########################################################################
# Loading data

X = read.table("~/Desktop/IBRI_consulting/kalwat/GSE104714_EtOH_avg_expression.txt",row.names=1,header=TRUE)

#########################################################################
# identifying optimal clusters

nclusts = 2:20
WCSSP = rep(NA,length(nclusts))
names(WCSSP) = as.character(nclusts)
kmeans.res = list()

for (n in nclusts){
  tmp = pheatmap(log2(X+1),kmeans_k=n)
  kmeans.res[[as.character(n)]] = tmp
  WCSSP[as.character(n)] = tmp$kmeans$tot.withinss/sum(tmp$kmeans$tot.withinss,tmp$kmeans$betweenss)
}

optimal.k = findElbow(nclusts,WCSSP)
pdf("~/Desktop/IBRI_consulting/kalwat/EtOH_kmeans_output/optimal_clusters.pdf")
plot(nclusts,WCSSP,type="l")
abline(v=optimal.k,col="red")
dev.off()

for(i in 1:optimal.k){
tmp = as.data.frame(t(apply(log2(X[kmeans.res[[as.character(optimal.k)]]$kmeans$cluster==i,]),1,centerFunc)))
df = data.frame(gene=rep(row.names(tmp),dim(tmp)[2]),time=rep(c(1,3,5,7,9,11),each=dim(tmp[1])),expr=unlist(tmp))
pdf(paste0("~/Desktop/IBRI_consulting/kalwat/EtOH_kmeans_output/genes_cluster",i,".pdf"))
plot(ggplot(df,aes(x=time,y=expr,group=gene)) + geom_line(aes(colour=gene)) + theme(legend.position = "none"))
dev.off()
}


