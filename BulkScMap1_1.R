#######################################################################################
# Bulk to scRNA-seq mapping (v1.1)
# Alan Rupp editted by Travis Johnson
# Main edits optimized script to run in roughly 1 hour for no null group
# and roughly 15 min using a null group (immature oligo) of 1692 cells

# This is an adaption of a script written by Alan Rupp to identify cell types enriched
# for disease associated genes. The code has been adapted to accept standardized input.

#######################################################################################
# Here you must supply the paths to the expression files and label files for the gene expression
# data. 

# 1) The bulk expression file should consist of a count matrix where genes are rows and 
#    samples are columns. The row names should be genes and the column names should be 
#    samples with no additional indexing. Sample names need to be unique. Duplicated 
#    row names will be removed based mean expression.

# 2) The bulk label file should consist of a label matrix where the samples are rows and 
#    the sample data are columns. The row names should contain the samples as the bulk 
#    expression matrix. They must be unique. The column names contain the sample characteristics.
#    These will be used to calculate the differentially expressed genes and stratify the data.

# 3) The single cell expression file can either be in the same format as the bulk expression file ,i.e.
# genes are rows and cells are columns, or in matrix market format, i.e. barcodes.tsv.gz, genes.tsv.gz,
# and matrix.mtx.gz .

sc_expr_path = "~/jonathan_flak/cns_diabetes/data/sc_data/GSE113576"
sc_meta_path = "~/jonathan_flak/cns_diabetes/data/sc_data/aau53224_Moffitt_Table-S1_TJclust.csv"
sc_meta_comp_column = "newtype"
sc_meta_subset_column = "type"
sc_meta_subset_condition = "%notin% c('Unstable','Ambiguous')"
sc_meta_comparison_group = "c('Immature oligodendrocyte')"
bulk_expr_path = "~/jonathan_flak/cns_diabetes/data/bulk_data/ILMN_1019_Flak_mRNAseq16_Aprl2021.csv"
bulk_meta_path = "~/jonathan_flak/cns_diabetes/data/bulk_data/ILMN_1019_Flak_metadata_Aprl2021.csv"
bulk_meta_comp_column = "diet"
bulk_meta_subset_column = "brain_region"
bulk_meta_subset_condition = "=='POA'"
out_name = "neuronSubtypes_oligoNull_dietEnriched_regionPOA_newV2"
bulkLog2FC_cutoff = log2(1.25)
bulkPvalue_cutoff = 0.05
scLog2FC_cutoff = log2(1.5)
scPvalue_cutoff = 0.05
scMinDiffPerc = 0.2
#######################################################################################
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

FindAllMarkersToBaseline <- function(seurat.obj,features,verbose,only.pos,
                                     logfc.threshold,min.diff.pct,baseline.groups){
  tmpDEG = list()
  test.groups <- unique(seurat.obj$meta_label)[unique(seurat.obj$meta_label) %notin% baseline.groups]
  for (g in test.groups){
    tmpCells <- names(seurat.obj$meta_label)[seurat.obj$meta_label %in% c(g,baseline.groups)]
    tmpSeurat <- subset(seurat.obj,cells=tmpCells)
    tmpDEG[[g]] <- FindMarkers(tmpSeurat,
                               ident.1=g,
                               features=features,
                               verbose=verbose,
                               only.pos=only.pos,
                               logfc.threshold=logfc.threshold,
                               min.diff.pct=min.diff.pct)
    tmpDEG[[g]]$cluster <- g
    tmpDEG[[g]]$gene <- row.names(tmpDEG[[g]])
  }
  return(do.call(rbind,tmpDEG))
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

get_log2FC <- function(cluster) {
  result <- pbapply::pbsapply(bulk_degs, function(x) {
    log2(
      mean(sc_expression[["RNA"]]@data[x, sc_expression@active.ident == cluster])/
        mean(sc_expression[["RNA"]]@data[x, sc_expression@active.ident %in% nullGrp])
    )
  })
  data.frame("cluster" = cluster, "gene" = bulk_degs, "avg_log2FC" = result)
}

quantile_normalize <- function(df) {
  # code adapted from Dave Tang
  ranks <- apply(df, 2, rank, ties.method = "min")
  means <- apply(df, 2, sort) %>% apply(., 1, mean)
  get_values <- function(rank_values) {
    values <- means[rank_values]
    # if ties, average each value by its position
    if (any(duplicated(rank_values))) {
      dups <- unique(rank_values[duplicated(rank_values)])
      new_means <- sapply(dups, function(x) {
        missing_ranks <- seq(x, x+sum(rank_values == x)-1)
        mean(means[missing_ranks])
      })
      for (i in 1:length(dups)) {
        values[rank_values == dups[i]] <- new_means[i]
      }
    }
    return(values)
  }
  new_values <- apply(ranks, 2, get_values)
  rownames(new_values) <- rownames(df); colnames(new_values) <- colnames(df)
  return(new_values)
}

#######################################################################################
# Loading data
# Load the bulk data
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

################################################################################
# Differentially expressed genes in bulk expression
bulk_deg_table <- get_enrichment(bulk_expression[common_genes,],bulk_metadata,bulk_meta_subset_column,
                               bulk_meta_subset_condition,bulk_meta_comp_column)

# enriched genes
bulk_degs <- row.names(subset(bulk_deg_table,pvalue < bulkPvalue_cutoff & log2FoldChange > bulkLog2FC_cutoff))

################################################################################
# Differentially expressed genes in sc data
sc_expression <- NormalizeData(sc_expression)
sc_expression$meta_label <- sc_metadata[,sc_meta_comp_column]
sc_expression@active.ident <- as.factor(sc_expression$meta_label)
names(sc_expression@active.ident) <- names(sc_expression$meta_label)

scdeg_data = paste0(sc_expr_path,sc_meta_path,sc_meta_comp_column,
                    sc_meta_subset_column,sc_meta_subset_condition,
                    sc_meta_comparison_group,
                    as.character(scPvalue_cutoff),
                    as.character(scLog2FC_cutoff),
                    as.character(scMinDiffPerc))
scdeg_data = paste0(gsub(" ","",gsub("[[:punct:]]","",scdeg_data)),".rds")

if(!file.exists(scdeg_data)){
  st = Sys.time()
  if(sc_meta_comparison_group==""){
    sc_degs <- FindAllMarkers(
      sc_expression,
      features = common_genes,
      verbose = TRUE,
      only.pos = TRUE,
      logfc.threshold=scLog2FC_cutoff,
      min.diff.pct=scMinDiffPerc
    )
  }else if(!any(eval(parse(text=sc_meta_comparison_group)) %notin% unique(sc_expression$meta_label))){
    sc_degs <- FindAllMarkersToBaseline(
      sc_expression,
      features = common_genes,
      verbose = TRUE,
      only.pos = TRUE,
      logfc.threshold=scLog2FC_cutoff,
      min.diff.pct=scMinDiffPerc,
      baseline.groups=eval(parse(text=sc_meta_comparison_group))
    )
  }else{
    message("ERROR: incorrect sc_meta_comparison_group parameter specified")
  }
  Sys.time()-st
saveRDS(sc_degs,file=scdeg_data)
}else{
  sc_degs = readRDS(scdeg_data)
}

################################################################################
# Hypergeometric test
hyper_results <- data.frame(
  "cluster" = unique(sc_degs$cluster),
  "p" = map_dbl(unique(sc_degs$cluster), 
  				function(x) hyper_test(x,sc_degs,bulk_degs,common_genes,
  				paste0("p_val<",scPvalue_cutoff)))
)

################################################################################
# GSEA test

nullGrp = eval(parse(text=sc_meta_comparison_group))
log_enrichment <- map_dfr(unique(sc_degs$cluster), get_log2FC)

# put data.frame in wide format with rownames
wide_enrichment <- pivot_wider(
  log_enrichment,
  names_from = "cluster",
  values_from = "avg_log2FC"
) %>%
  as.data.frame() %>%
  column_to_rownames("gene")

# filling in missing values
min_value <- unlist(wide_enrichment) %>% .[is.finite(.)] %>% min()
wide_enrichment <- apply(wide_enrichment, c(1, 2), function(x) {
  ifelse(is.finite(x), x, min_value)
})

# find enrichment score by calculating the cumulative sum of the fold-enrichment
es <- apply(wide_enrichment, 2, cumsum)
es_max <- apply(es, 2, max)

# to get a P value for association, I use random permutations of the markers
null_markers <- replicate(
  1000,
  apply(wide_enrichment, 1, function(x) sample(x, 1))
)

# find the number of times that a random marker set exceeds max for each cluster
null_es <- apply(null_markers, 2, cumsum)
null_es_max <- apply(null_es, 2, max)

# calculate P value from null datasets
p_vals <- sapply(es_max, function(x) sum(null_es_max > x)/length(null_es_max))
gsea_results <- data.frame(cluster=names(p_vals),p=as.numeric(p_vals))

################################################################################
# Visualizing results
pdf(paste0(out_name,"_hyperEnriched.pdf"),height = 7, width = 14)
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

pdf(paste0(out_name,"_gseaEnriched.pdf"),height = 7, width = 14)
ggplot(gsea_results, aes(x = cluster, y = -log10(p), color = cluster)) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_point() + theme(legend.position = "none",
                       axis.text.x = element_text(angle = 90,
                                                  hjust = 1, vjust = 0.5,
                                                  color = "black", size = 8)) + 
  geom_text_repel(data = filter(gsea_results, p < 0.05),
                  aes(label = cluster), fontface = "bold", size = 3,
                  ylim = c(-log10(0.05), NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, NA)) +
  xlab(NULL) +
  ylab(expression(italic("P")*" value (-log"[10]*")"))
dev.off()


################################################################################
# Output the significance tables
write.csv(hyper_results,file=paste0(out_name,"_hyperEnriched.csv"))
write.csv(hyper_results,file=paste0(out_name,"_gseaEnriched.csv"))
