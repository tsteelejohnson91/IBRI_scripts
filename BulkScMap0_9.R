#######################################################################################
# Bulk to scRNA-seq mapping (v0.9)
# Alan Rupp editted by Travis Johnson
# Main edits optimized script to run in roughly 1 hour

# This is an adaption of a script written by Alan Rupp to identify cell types enriched
# for disease associated genes. The code has been adapted to accept standardized input.

#######################################################################################
# Here you must supply the paths to the expression files and label files for the gene expression
# data. 

# 1) The bulk expression file should consist of a count matrix where genes are rows and samples are columns. The row names should be genes and the column names should be samples with no additional indexing. Sample names need to be unique. Duplicated row names will be removed based mean expression.

# 2) The bulk label file should consist of a label matrix where the samples are rows and the sample data are columns. The row names should contain the samples as the bulk expression matrix. They must be unique. The column names contain the sample characteristics. These will be used to calculate the differentially expressed genes and stratify the data.

# 3) The single cell expression file can either be in the same format as the bulk expression file ,i.e.
# genes are rows and cells are columns, or in matrix market format, i.e. barcodes.tsv.gz, genes.tsv.gz,
# and matrix.mtx.gz .

sc_expr_path = "~/jonathan_flak/cns_diabetes/data/sc_data/GSE113576"
sc_meta_path = "~/jonathan_flak/cns_diabetes/data/sc_data/aau53224_Moffitt_Table-S1_formatted.csv"
sc_meta_comp_column = "subtype"
sc_meta_subset_column = "type"
sc_meta_subset_condition = "%notin% c('Unstable','Ambiguous')"
bulk_expr_path = "~/jonathan_flak/cns_diabetes/data/bulk_data/ILMN_1019_Flak_mRNAseq16_Aprl2021.csv"
bulk_meta_path = "~/jonathan_flak/cns_diabetes/data/bulk_data/ILMN_1019_Flak_metadata_Aprl2021.csv"
bulk_meta_comp_column = "diet"
bulk_meta_subset_column = ""
bulk_meta_subset_condition = ""
out_name = "neuralSubtype_dietEnriched_regionControlled2"
bulkLog2FC_cutoff = log2(1.5)
bulkPvalue_cutoff = 0.05
scLog2FC_cutoff = log2(1.5)
scPvalue_cutoff = 0.05

#######################################################################################
library(knitr)
library(kableExtra)
library(tidyverse)
library(ggrepel)
library(readxl)
library(Seurat)
library(DESeq2)

`%notin%` <- Negate(`%in%`)

get_extension <- function(file){ 
    ex <- strsplit(basename(file), split="\\.")[[1]]
    return(ex[-1])
}

remove_duplicate_gene_names <- function(data){
  dups <- data[duplicated(data[,1]),1] %>% unique()
  drop <- c()
  if(length(dups)>0){
    for (dup in dups) {
      ids <- which(data[,1] == dup)
      expr <- rowSums(data[ids,2:dim(data)[2]])
      expr <- sort(expr, decreasing = TRUE)
      drop <- c(drop, names(expr)[2:length(expr)])
    }
    data <- data[-as.numeric(drop), ]
    data <- as.data.frame(data)
  }
  row.names(data) <- data[,1]
  data <- data[,-1]
  return(data)
}

load_bulk_expression <- function(file){
  ext <- get_extension(file)
  if(ext == "xlsx"){
    data <- read_xlsx(file)
  }else if (ext == "txt"){
    data <- read.table(file,sep = "\t")
  }else if (ext == "csv"){
    data <- read.csv(file)
  }else{
    message("ERROR: No bulk expression file found")
    break
  }
  data <- remove_duplicate_gene_names(data)
  return(data)
}

load_metadata <- function(file){
  ext <- get_extension(file)
  if(ext == "xlsx"){
    labels <- read_xlsx(file)
  }else if (ext == "txt"){
    labels <- read.table(file,sep = "\t")
  }else if (ext == "csv"){
    labels <- read.csv(file)
  }else{
    message("ERROR: No label file found")
    break
  }
  row.names(labels) <- labels[,1]
  labels = labels[,-1]
  return(labels)
}

load_sc_expression <- function(file){
  if(dir.exists(file)){
    f_list = list.files(file)
    if(length(grep("barcodes.tsv",f_list))>0 && 
       length(grep("genes.tsv",f_list))>0 && 
       length(grep("matrix.mtx",f_list))>0){
      data <- Read10X(file,gene.column=1)
      feats <- read.table(paste0(file,"/genes.tsv"),sep="\t")
      row.names(data) <- NULL
      data <- remove_duplicate_gene_names(cbind(feats[,2],as.data.frame(data)))
    }else{
      message("ERROR: No SC expression file found")
      break
    }
  }else if(file.exists(file)){
    ext = get_extension(file)
    if(ext == "xlsx"){
      data <- read_xlsx(file)
      data <- remove_duplicate_gene_names(data)
    }else if (ext == "txt"){
      data <- read.table(file,sep = "\t")
      data <- remove_duplicate_gene_names(data)
    }else if (ext == "csv"){
      data <- read.csv(file)
      data <- remove_duplicate_gene_names(data)
    }else{
      message("ERROR: No SC expression file found")
      break
    }
  }else{
    message("ERROR: No SC expression file found")
    break
  }
  return(data)
}

list_to_string <- function(str_list,delim){
  out = ""
  if(length(str_list)==1){
    out = str_list[1]
  }else{
    for (str in str_list){
      if (str != str_list[length(str_list)]){
        out = paste0(out,str,delim)
      }else{
        out = paste0(out,str)
      }
    }
  }
  return(out)
}

# Can only handle binary comparisons
get_enrichment <- function(expr,meta,subset_col,subset_cond,comp_col){
  meta = meta[colnames(expr),]
  if(subset_col==""){
    com_str <- paste0("DESeqDataSetFromMatrix(countData = expr,colData = meta,design = ~ ",
                      list_to_string(colnames(meta),"+"),")")
    dds <- eval(parse(text=com_str))
    dds <- DESeq(dds)
    res <- results(dds, name=resultsNames(dds)[grep(comp_col,resultsNames(dds))])
    res <- lfcShrink(dds,coef=resultsNames(dds)[grep(comp_col,resultsNames(dds))],type="apeglm")
  }else if(subset_cond!=""){
    com_str <- paste0("DESeqDataSetFromMatrix(countData = expr2,colData = meta2,design = ~ ",
                      comp_col,")")
    expr2 = eval(parse(text=paste0("expr[,meta[,'",subset_col,"']",subset_cond,"]")))
    meta2 = eval(parse(text=paste0("meta[meta[,'",subset_col,"']",subset_cond,",]")))
    dds <- eval(parse(text=com_str))
    dds <- DESeq(dds)
    res <- results(dds, name=resultsNames(dds)[grep(comp_col,resultsNames(dds))])
    res <- lfcShrink(dds,coef=resultsNames(dds)[grep(comp_col,resultsNames(dds))],type="apeglm")
  }else{
    message("ERROR: subset condition arguments missing")
  }
  res = as.data.frame(res)
  res = res[order(res$pvalue),]
  return(res)
}

hyper_test <- function(ident,markersSC,markersBLK,all_genes,deg_cutoffs) {
  genes <- eval(parse(text=paste0("filter(markersSC, cluster == ident & ",deg_cutoffs,") %>% pull(gene)")))
  p = phyper(
    q = length(intersect(genes,markersBLK)),
    m = length(genes),
    n = length(common_genes),
    k = length(markersBLK),
    lower.tail = FALSE
  )
  return(p)
}

#######################################################################################
# read in bulk data
bulk_metadata <- load_metadata(bulk_meta_path)
bulk_expression <- load_bulk_expression(bulk_expr_path)

# Load the sc data
sc_metadata <- load_metadata(sc_meta_path)
sc_metadata <- eval(parse(text=paste0("sc_metadata[sc_metadata[,sc_meta_subset_column]",sc_meta_subset_condition,",]")))
sc_expression <- load_sc_expression(sc_expr_path)
keep_cells <- intersect(row.names(sc_metadata),colnames(sc_expression))
sc_metadata <- sc_metadata[keep_cells,]
sc_expression <- sc_expression[,keep_cells]
sc_expression <- CreateSeuratObject(sc_expression,min.cells=0,min.features=0)

common_genes <- intersect(rownames(sc_expression), rownames(bulk_expression))

# Differentially expressed genes in bulk expression
bulk_deg_table <- get_enrichment(bulk_expression[common_genes,],bulk_metadata,bulk_meta_subset_column,
                               bulk_meta_subset_condition,bulk_meta_comp_column)

# enriched genes
bulk_degs <- row.names(subset(bulk_deg_table,pvalue < bulkPvalue_cutoff & log2FoldChange > bulkLog2FC_cutoff))

# Single cell cluster enrichment
sc_expression <- NormalizeData(sc_expression)
sc_expression$meta_label <- sc_metadata[,sc_meta_comp_column]
sc_expression@active.ident <- as.factor(sc_expression$meta_label)
names(sc_expression@active.ident) <- names(sc_expression$meta_label)

scdeg_data = paste0(sc_expr_path,sc_meta_path,sc_meta_comp_column,
                    sc_meta_subset_column,sc_meta_subset_condition,
                    as.character(scPvalue_cutoff),
                    as.character(scLog2FC_cutoff))
scdeg_data = paste0(gsub(" ","",gsub("[[:punct:]]","",scdeg_data)),".rds")

if(!file.exists(scdeg_data)){
  st = Sys.time()
  sc_degs <- FindAllMarkers(
    sc_expression,
    features = common_genes,
    verbose = TRUE,
    only.pos = TRUE,
    logfc.threshold=scLog2FC_cutoff,
    min.diff.pct=.2
  )
  Sys.time()-st
saveRDS(sc_degs,file=scdeg_data)
}else{
  sc_degs = readRDS(scdeg_data)
}

# Hypergeometric test
hyper_results <- data.frame(
  "cluster" = levels(sc_expression@active.ident),
  "p" = map_dbl(levels(sc_expression@active.ident), function(x) hyper_test(x,sc_degs,bulk_degs,common_genes,paste0("p_val<",scPvalue_cutoff)))
)

pdf(paste0(out_name,".pdf"),height = 7, width = 14)
ggplot(hyper_results, aes(x = cluster, y = -log10(p), color = cluster)) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_point() + theme(legend.position = "none",
                       axis.text.x = element_text(angle = 90,
                                                  hjust = 1, vjust = 0.5,
                                                  color = "black", size = 8)) + 
  geom_text_repel(data = filter(hyper_results, p < 0.05),
                  aes(label = cluster), fontface = "bold", size = 3,
                  ylim = c(-log10(0.05), NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, NA)) +
  xlab(NULL) +
  ylab(expression(italic("P")*" value (-log"[10]*")"))
dev.off()


# Output the significance table
write.csv(hyper_results,file=paste0(out_name,".csv"))
