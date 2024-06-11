#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("TCGAbiolinks")

dir <- "C:/Users/suade/Documents/Universidad/UOC/2/TFM/Código"
setwd(dir)
library(TCGAbiolinks)
library(tidyverse)

# Función para crear una consulta y guardar los datos
geneDataTCGA <- function(proj, path, sample_type) {
  # crear la consulta 
  query <- GDCquery(project = proj,
                            data.category = 'Transcriptome Profiling',
                            data.type = 'Gene Expression Quantification',
                            experimental.strategy = 'RNA-Seq',
                            workflow.type = 'STAR - Counts',
                            access = 'open',
                            sample.type = sample_type)
  # Descargar los datos
  GDCdownload(query, directory = path)
  # Leer los datos
  rna <- GDCprepare(query, directory = path)
  # Guardar los datos en un archivo RDS
  write_rds(rna, file=paste0(path, "/tpm_", proj, ".RDS"))
}

# Función para extraer los datos
extractAllGenes <- function(data, metric) {
  # extraer toda la expresión génica
  tpmDat <- as.data.frame(assays(data)[[metric]])
  # seleccionar tpm
  tpmGene <- as.data.frame(t(tpmDat)) %>%
  rownames_to_column("barcode")
  return(tpmGene)
}

filter_and_normalize <- function(data) {
  # Filtrar genes con baja expresión
  data_filtered <- data %>%
    select_if(~ sum(. > 1) > 0.1 * nrow(data))
  
  # Normalización logarítmica
  data_normalized <- data_filtered %>%
    mutate(across(-barcode, ~ log2(. + 1)))
  
  return(data_normalized)
}


# Definir los directorios de descarga
download_dir_TP <- file.path(getwd(), "BRCA_data/TP/GDCdata/")
download_dir_TM <- file.path(getwd(), "BRCA_data/TM/GDCdata/")

# Crear los directorios si no existen
if (!dir.exists(download_dir_TP)) {
  dir.create(download_dir_TP, recursive = TRUE)
}
if (!dir.exists(download_dir_TM)) {
  dir.create(download_dir_TM, recursive = TRUE)
}

# Lista de los proyectos disponibles
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-BRCA')
getProjectSummary('CMI-MBC')
getProjectSummary('CMI-ASC')

# Descargar datos
geneDataTCGA('TCGA-BRCA', download_dir_TP, c("Primary Tumor"))
geneDataTCGA('TCGA-BRCA', download_dir_TM, c("Metastatic"))
geneDataTCGA('CMI-MBC', download_dir_TM, c("Metastatic"))
geneDataTCGA('CMI-ASC', download_dir_TM, c("Metastatic"))


TM_genes_dat <- data.frame()
TP_genes_dat <- data.frame()

tpm <- read_rds(paste0(download_dir_TP, "/tpm_TCGA-BRCA", ".RDS"))
TP_genes_temp <- extractAllGenes(tpm, "tpm_unstrand")
write_csv(TP_genes_temp, file="BRCA_data/BRCA-TP_genes.csv")


tpm <- read_rds(paste0(download_dir_TM, "/tpm_TCGA-BRCA", ".RDS"))
TM_genes_temp <- extractAllGenes(tpm, "tpm_unstrand")
write_csv(TM_genes_temp, file="BRCA_data/BRCA-TM_genes_TCGA-BRCA.csv")
TM_genes_dat <- bind_rows(TM_genes_dat, TM_genes_temp)

tpm <- read_rds(paste0(download_dir_TM, "/tpm_CMI-MBC", ".RDS"))
TM_genes_temp <- extractAllGenes(tpm, "tpm_unstrand")
write_csv(TM_genes_temp, file="BRCA_data/BRCA-TM_genes_CMI-MBC.csv")
TM_genes_dat <- bind_rows(TM_genes_dat, TM_genes_temp)

tpm <- read_rds(paste0(download_dir_TM, "/tpm_CMI-ASC", ".RDS"))
TM_genes_temp <- extractAllGenes(tpm, "tpm_unstrand")
write_csv(TM_genes_temp, file="BRCA_data/BRCA-TM_genes_CMI-ASC.csv")
TM_genes_dat <- bind_rows(TM_genes_dat, TM_genes_temp)

write_csv(TM_genes_dat, "BRCA_data/BRCA-TM_genes.csv")






