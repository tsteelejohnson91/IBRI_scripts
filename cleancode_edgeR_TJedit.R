####################################################################
# This file performs DGE analysis with edgeR
# Gitanjali Roy and Travis S. Johnson

library(edgeR)
library(limma)
library(magrittr)
cnt_type = "DMSO"
trt_type = "SW"
file_path = "~/Desktop/IBRI_consulting/kalwat/edgeR_TJ/"

## Load data I had to change the naming convention here (this part will need to be edited for future experiments)
# This part is a little tricky because the baselines are duplicates for DMSO and treatment.
# This shouldn't matter because DMSO == treatment at baseline according to our conversations.
x <- read.csv(paste0(file_path,"sample_matrix_rawcount_",trt_type,".csv"), header=T, stringsAsFactors=FALSE,row.names=1)
newcols = unlist(lapply(colnames(x),function(x) if(substr(x,1,8)=="Baseline"){return(c(sub("Baseline",cnt_type,x),sub("Baseline",trt_type,x)))}else{return(x)}))
x = x[,unlist(lapply(colnames(x),function(x) if(substr(x,1,8)=="Baseline"){return(c(x,x))}else{return(x)}))]
colnames(x) = newcols
x = as.data.frame(x)

## Setting grouping variables
time<-factor(sub("[_].*","",sub(".*[.]","",colnames(x))),levels=c("0h","1h","2h","6h","24h"))
treat<-factor(sub("[.].*","",colnames(x)),levels=c(cnt_type,trt_type))
levels(time)
levels(treat)

## Creating DGEList object
y<-DGEList(counts=x, genes = row.names(x), group=paste0(as.character(treat),"_",as.character(time)))
y <- calcNormFactors(y)
y
plotMDS(y, method="bcv", col=ifelse(sub("[_].*","",y$samples$group)==trt_type,"red","blue"))


## Filter out lowly expressed genes:
keep <- filterByExpr(y)
table(keep)
y <- y[keep,,keep.lib.sizes=FALSE]

## TMM Normalization
y <- calcNormFactors(y)
y$samples
plotMDS(y, labels = time, method="bcv", col=ifelse(sub("[_].*","",y$samples$group)==trt_type,"red","blue"))
# multi-dimensional scaling (MDS) plots visualizes the differences between the expression profiles of different samples in two dimensions.

## I changed the model matrix so that we can get the time+treatment dependent effects
design <- model.matrix(~time+treat+time*treat, data = y$samples) # the 0+ in the model formula is an instruction not to include an intercept column and
#instead to include a column for each group.Removing the intercept (with ~0) is only used when you intend to form contrasts from the levels of a single factor and should not be used in any other circumstances.
design

## Estimate dispersion:
z <- estimateDisp(y,design)
sqrt(z$common.dispersion)
plotBCV(z)
plotMDS(z, labels = time, method="bcv", col=ifelse(sub("[_].*","",y$samples$group)==trt_type,"red","blue"))

#To perform quasi-likelihood F-tests:
fit <- glmQLFit(z,design)
plotQLDisp(fit)

## Time course trend analysis
# I think this is what we want. Using the new design matrix we are choosing genes that are both significantly
# altered over time and based on treatment. I.e. the interaction between Tg treatment and time is significant.
test_coefficients = colnames(design)[7:10] # For future analyses this may need to be changed
test_coefficients
qlf <- glmQLFTest(fit,coef=test_coefficients)
tab1 <- as.data.frame(topTags(qlf, n=Inf))
# requiring a FC > 2 and BH-FDR < 0.05 (we can also use FC > 1.5 for more genes)
keep = decideTests(qlf,lfc=log2(2),p.value=0.05)
summary(keep)
tab2 = tab1[row.names(keep)[keep==1],]
write.csv(tab2, paste0(file_path,"exported(QLFT)_",trt_type,"_TJ.csv"))
tab3 = tab2[,2:5]
colnames(tab3) = c("1","2","6","24")
write.table(tab3,file=paste0(file_path,"final_data_matrix_",trt_type,"_TJ.txt"),sep='\t',row.names=TRUE,quote=FALSE)
