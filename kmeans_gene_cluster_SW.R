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
  NkDiffs = rep(NA,length(Nk))
  names(NkErrors) = names(NkDiffs) = as.character(Nk)
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
    NkErrors[as.character(k)] = sum(linfit1res,linfit2res)/2
    NkDiffs[as.character(k)] = abs(mean(linfit1res)-mean(linfit2res))
  }
  NkErrorsDiffs = NkErrors+NkDiffs
  names(NkErrorsDiffs) = names(NkErrors)
  return(as.numeric(names(NkErrorsDiffs)[NkErrorsDiffs==min(NkErrorsDiffs)]))
}

#########################################################################
# Loading data
X = read.csv("~/Desktop/IBRI_consulting/kalwat/edgeR_TJ/exported(QLFT)_SW_TJ.csv")
row.names(X) = X$X
X = X[,c("logFC.time1h.treatSW","logFC.time2h.treatSW","logFC.time6h.treatSW","logFC.time24h.treatSW")]
colnames(X) = c("1","2","6","24")


#########################################################################
# identifying optimal clusters

nclusts = 2:20
WCSSP = rep(NA,length(nclusts))
names(WCSSP) = as.character(nclusts)
kmeans.res = list()

for (n in nclusts){
  tmp = pheatmap(X,kmeans_k=n)
  kmeans.res[[as.character(n)]] = tmp
  WCSSP[as.character(n)] = tmp$kmeans$tot.withinss/sum(tmp$kmeans$tot.withinss,tmp$kmeans$betweenss)
}

test.wcssp = c(15, 14, 13, 12, 11, 10, 9, 8, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3, 2.5, 2)
message(paste0("The correct elbow should be 10. The findElbow function produced: ",findElbow(nclusts,test.wcssp),"."))
optimal.k = findElbow(nclusts,WCSSP)
pdf("~/Desktop/IBRI_consulting/kalwat/SW_kmeans_output/optimal_clusters.pdf")
plot(nclusts,WCSSP,type="l")
abline(v=optimal.k,col="red")
dev.off()

pdf("~/Desktop/IBRI_consulting/kalwat/SW_kmeans_output/optimal_cluster_heatmap.pdf")
kmeans.res[[as.character(optimal.k)]]
dev.off()

for(i in 1:optimal.k){
tmp = as.data.frame(t(apply(X[kmeans.res[[as.character(optimal.k)]]$kmeans$cluster==i,],1,centerFunc)))
df = data.frame(gene=rep(row.names(tmp),dim(tmp)[2]),time=rep(as.numeric(colnames(X)),each=dim(tmp)[1]),expr=unlist(tmp))
pdf(paste0("~/Desktop/IBRI_consulting/kalwat/SW_kmeans_output/genes_cluster",i,".pdf"))
plot(ggplot(df,aes(x=time,y=expr,group=gene)) + geom_line(aes(colour=gene)) + theme(legend.position = "none"))
dev.off()
}

write.table(X,file="~/Desktop/IBRI_consulting/kalwat/SW_kmeans_output/final_data_matrix.txt",sep='\t',quote=FALSE)

# Move to DP_GP_cluster directory
# DP_GP_cluster.py -i /Users/johnstrs/Desktop/IBRI_consulting/kalwat/SW_kmeans_output/final_data_matrix_TJ.txt -o /Users/johnstrs/Desktop/IBRI_consulting/kalwat/SW_DPGPtruetimeoff_output/SWDPGP --plot -p png --fast
