library(TCGA2STAT)
library(dplyr)
library(knitr)
data_dir="/Volumes/Share/TCGApublicdata/" #sever
results_dir = data_dir

#General settings
data.type = "Mutation"
type = "somatic"
type = "all"
clinical = TRUE

## 8 cancer data needed to be analyzed
# 1.Breast cancer
cancer = "BRCA"
# 2.Cervical and endocervical cancers
cancer ="CESC"
# 3.Ovarian cancer
cancer = "OV"
# 4.Colorectal adenocarcinoma
cancer ="COADREAD"
# 4.Colon adenocarcinoma
cancer = "COAD"
# 5.Glioblastoma multiforme
cancer="GBM"
# 5.Glioma
cancer="GBMLGG"
# 6.Lung adenocarcinoma
cancer = "LUAD"
# 6. Lung squamous cell carcinoma
cancer = "LUSC"
# 7.Ovarian serous cystadenocarcinoma
cancer ="OV"
# 8.Pancreatic adenocarcinoma
cancer = "PAAD"
# Prostate adenocarcinoma
cancer ="PRAD"
# Stomach adenocarcinoma 
cancer = "STAD"
# Stomach and Esophageal carcinoma 
cancer = "STES"

# function to load data from remote repository
load_data<-function(disease = cancer, data.type=data.type, type = type, data_dir=data_dir,force_reload=TRUE){
  FILE=paste0(data_dir,"/mtx_",disease,"_",data.type,"_",type,"_",".rda")
  if (all(file.exists(FILE),!(force_reload))){
    load(file=FILE)
  }else{
    mtx<-getTCGA(disease = disease, data.type = data.type, type = type, clinical = TRUE)
    save(file = FILE,list = c("mtx"))
  }
  return(mtx)
}

#funtion to overview downloaded data
summarize_data <- function(mtx = mtx) {
  print(paste0("Dimensions of clinical matrix, patients X parameters: ", paste(dim(mtx$clinical), collapse = " ")))
  print("Head of the clinical matrix")
  print(mtx$clinical[1:5, 1:7])
  print("List of clinical values, and frequency of each variable: ")
  clin_vars <- apply(mtx$clinical, 2, function(x) length(table(x[ !(is.na(x) & x != "" )]))) %>% as.data.frame()
  # Filter clinical variables to have at least 2, but no more than 10 categories,
  # And they are not dates
  clin_vars <- clin_vars[ as.numeric(clin_vars$.) > 1 & as.numeric(clin_vars$.) < 10 & !grepl("years|days|date|vital", rownames(clin_vars), perl = TRUE) , , drop = FALSE]
  print(kable(clin_vars))
  return(rownames(clin_vars))
}

# A function to create expression matrix
make_expression_matrix <- function(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir) {
  # transposed expression matrix (genes start at column 4) 
  mtx.expression <- mtx$merged.dat[, 4:ncol(mtx$merged.dat) ] %>% t 
  # Set column names as patient IDs
  colnames(mtx.expression) <- mtx$merged.dat$bcr 
  # Set row names as probe IDs
  rownames(mtx.expression) <- colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ] 
  # Save gzipped matrix
  fileName.gz <- gzfile(paste0(results_dir, "/mtx_", disease, "_", data.type, "_", type, "_1expression.txt.gz"), "w")
  write.table(mtx.expression, fileName.gz, sep = ";", quote = FALSE)
  close(fileName.gz)
}

# A function to create probe ID - gene symbol mapping
make_mapping_matrix <- function(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir) {
  mtx.mapping <- data.frame(probeID = colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ],
                            geneID = colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ])
  # Save gzipped matrix
  fileName.gz <- gzfile(paste0(results_dir, "/mtx_", disease, "_", data.type, "_", type, "_2mapping.txt.gz"), "w")
  write.table(mtx.mapping, fileName.gz, sep = ";", quote = FALSE, col.names = FALSE, row.names = FALSE)
  close(fileName.gz)
}

# A function to create sample annotation matrix
make_annotation_matrix <- function(mtx = mtx, disease = cancer, data.type = data.type, type = type, clinical_annotations = clinical_annotations, results_dir = results_dir) {
  mtx.sample <- mtx$merged.dat[, c("bcr", "OS", "status")] # First 3 columns are c("sample_id", "surv_time", "dead_1_alive_0")
  colnames(mtx.sample) <- c("sample_id", "surv_time", "dead_1_alive_0")
  # Append selected clinical annotations
  mtx.sample <- left_join(mtx.sample, data.frame(sample_id = rownames(mtx$clinical), mtx$clinical[, colnames(mtx$clinical) %in% clinical_annotations], stringsAsFactors = FALSE), by = "sample_id")
  mtx.sample[ is.na(mtx.sample) ] <- "N/A"
  mtx.sample <- mtx.sample[, apply(mtx.sample, 2, function(x) length(x[ x != "N/A"])) > 40] # Keep clinical annitations with at least 40 patients
  # Save gzipped matrix
  fileName.gz <- gzfile(paste0(results_dir, "/mtx_", disease, "_", data.type, "_", type, "_3sample.txt.gz"), "w")
  write.table(mtx.sample, fileName.gz, sep = ";", quote = FALSE, row.names = FALSE)
  close(fileName.gz)
}


# mtx <- load_data(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)

# clinical_annotations <- summarize_data(mtx = mtx)

# make_annotation_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, clinical_annotations = clinical_annotations, results_dir = results_dir)

# Get data for all cancers
get_data <- function(cancers, data.type = data.type, type = type, data_dir = data_dir, force_reload) {
  for (cancer in cancers) {
    print(paste0("Processing cancer ", cancer))
    mtx <- load_data(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload)
    clinical_annotations <- summarize_data(mtx = mtx)
    make_annotation_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, clinical_annotations = clinical_annotations, results_dir = results_dir)
    rm(list = c("mtx", "clinical_annotations"))
  }
}

# All cancers with Mutation data
data.type = "Mutation"; type = "" 
cancer_TCGA = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS")
# All cancers with miRNAseq data. Uncomment to get miRNAseq data
# data.type = "miRNASeq"; type = "rpmmm"
# cancer_TCGA = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS")

# sink("TCGA2SynTarget.txt", split = FALSE)
get_data(cancers = cancer_TCGA, data.type = data.type, type = type, data_dir = data_dir, force_reload = TRUE)
# sink(type = "message")
# sink()

# Cleanup intermediate files
unlink(paste0(data_dir, "/*.txt.gz"))
