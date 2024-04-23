#Download data from NCBI using GEO accession GSE183947 

#Extract gz file using gunzip in terminal->Data has Genes in Y-axis and Samples(30 normal,30 cancerous tissue) in X-axis

#Script to run survival analysis using TCGA Data -> survival analysis is branch of statistics for analyzing the expected duration of time until one event occurs, such as death in biological organisms and failure in systems.


library(TCGAbiolinks) #for integrative analysis with TCGA data
library(survminer) #provides functions for facilitating survival visualization and interpretation of results.
library(survival) #provides functions and tools for survival analysis
library(SummarizedExperiment) #store and manipulate high-dimensional biological data
library(tidyverse) #contains packages for data cleaning, manipulation, visualization, and analysis
library(DESeq2) #to calculate the difference in gene expression between two or more groups of samples

#summary of data

getProjectSummary("TCGA-THCA")

#getting clinical data for TCGA-BRCA(breast cancer gene) cohort using TCGAbiolinks
clinical_thca<-GDCquery_clinic("TCGA-THCA")
any(colnames(clinical_thca) %in% c("vital_status","days_to_last_follow_up","days_to_death"))
which(colnames(clinical_thca) %in% c("vital_status","days_to_last_follow_up","days_to_death"))
clinical_thca[,c(9,38,44)]

#no of people dead and alive
table(clinical_thca$vital_status)

#for survival analysis, we need time of event, status information and event information
#change certain values the way they are encoded
#1- (status information)
clinical_thca$deceased <- ifelse(clinical_thca$vital_status=="Alive",FALSE, TRUE)
table(clinical_thca$deceased)

#2- time of event 
clinical_thca$overall_survival <- ifelse(clinical_thca$vital_status=="Alive",clinical_thca$days_to_last_follow_up,clinical_thca$days_to_death)

#3- event information
  #survival probabilities of 2 groups of patients high and low expression of TP53 (which needs gene expression of data)

#download gene expression data

query_THCA <- GDCquery(project ="TCGA-THCA", data.category = "Transcriptome Profiling",experimental.strategy = "RNA-Seq", workflow.type = "STAR - Counts",data.type = "Gene Expression Quantification", sample.type = "Primary Tumor", access = "open")

output_query_THCA <- getResults(query_THCA)

#get 30 datas
tumor<- output_query_THCA$cases[1:30]
tumor 

query_THCA <- GDCquery(project ="TCGA-THCA", data.category = "Transcriptome Profiling",experimental.strategy = "RNA-Seq", workflow.type = "STAR - Counts",data.type = "Gene Expression Quantification", sample.type = c("Primary Tumor","Solid Tissue Normal"), access = "open", barcode=tumor)

#download data
GDCdownload(query_THCA)

#prepare data
tcga_thca_data<-GDCprepare(query_THCA, summarizedExperiment = TRUE) #extract data from files
thca_matrix <-assay(tcga_thca_data,"unstranded") #count of data

#extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_thca_data))
coldata <- as.data.frame(colData(tcga_thca_data))

#variant stabilizing transformation (Gaussian shape) counts used in survival analysis  
dds <- DESeqDataSetFromMatrix(countData = thca_matrix, colData = coldata, design = ~1)

#removing genes with the sum total of 10 reads from all the samples(low expressing genes)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#vst
vsd <- vst(dds, blind = FALSE)
thca_matrix_vst<- assay(vsd)  #variant stabilized transform counts
thca_matrix_vst[1:10,1:10]

#get data for TP53 gene , add gene symbols(not so informative atm) to gene ids , change shape of the data
thca_tp53 <-thca_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = "counts", -gene_id) %>%#change shape 
  left_join(.,gene_metadata,by="gene_id") %>%
  filter(gene_name == "TP53")


#divide into low and high expression group 
#gene median value
median_value <- median(thca_tp53$counts)

#condition for dividing groups
thca_tp53$strata <- ifelse(thca_tp53$counts >= median_value, "HIGH", "LOW")

#now as we have 3 information for survival analysis
#let's merge all dataframes but case_id are different so change case_id
thca_tp53$case_id <- gsub('-01.*','',thca_tp53$case_id)
thca_tp53 <- merge(thca_tp53, clinical_thca,
      by.x = "case_id", by.y = "submitter_id", all = FALSE, 
      sort = TRUE, suffixes = c(".x",".y"), no.dups = TRUE,
      incomparables = NULL)

#fitting survival curve

fit <- survfit(Surv(overall_survival,deceased) ~ strata, data = thca_tp53)
fit
ggsurvplot(fit, data = thca_tp53, pval=T, risk.table = T)

fit2<- survdiff(Surv(overall_survival.y,deceased.y) ~ strata, data = thca_tp53)
fit2

