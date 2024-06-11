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
download_dir_TP <- file.path(getwd(), "PRAD_data/TP/GDCdata/")
download_dir_TM <- file.path(getwd(), "PRAD_data/TM/GDCdata/")

# Crear los directorios si no existen
if (!dir.exists(download_dir_TP)) {
  dir.create(download_dir_TP, recursive = TRUE)
}
if (!dir.exists(download_dir_TM)) {
  dir.create(download_dir_TM, recursive = TRUE)
}

# Lista de los proyectos disponibles
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-PRAD')
getProjectSummary('CMI-MPC')
getProjectSummary('WCDT-MCRPC')

# Descargar datos
geneDataTCGA('TCGA-PRAD', download_dir_TP, c("Primary Tumor"))
geneDataTCGA('TCGA-PRAD', download_dir_TM, c("Metastatic"))
geneDataTCGA('CMI-MPC', download_dir_TP, c("Primary Tumor"))
geneDataTCGA('WCDT-MCRPC', download_dir_TM, c("Metastatic"))


TM_genes_dat <- data.frame()
TP_genes_dat <- data.frame()

tpm <- read_rds(paste0(download_dir_TP, "/tpm_TCGA-PRAD", ".RDS"))
TP_genes_temp <- extractAllGenes(tpm, "tpm_unstrand")
write_csv(TP_genes_temp, file="PRAD_data/PRAD-TP_genes_TCGA-PRAD.csv")
TP_genes_dat <- bind_rows(TP_genes_dat, TP_genes_temp)

tpm <- read_rds(paste0(download_dir_TP, "/tpm_CMI-MPC", ".RDS"))
TP_genes_temp <- extractAllGenes(tpm, "tpm_unstrand")
write_csv(TP_genes_temp, file="PRAD_data/PRAD-TP_genes-CMI-MPC.csv")
TP_genes_dat <- bind_rows(TP_genes_dat, TP_genes_temp)

write_csv(TP_genes_dat, "PRAD_data/PRAD-TP_genes.csv")


tpm <- read_rds(paste0(download_dir_TM, "/tpm_TCGA-PRAD", ".RDS"))
TM_genes_temp <- extractAllGenes(tpm, "tpm_unstrand")
write_csv(TM_genes_temp, file="PRAD_data/PRAD-TM_genes_TCGA-PRAD.csv")
TM_genes_dat <- bind_rows(TM_genes_dat, TM_genes_temp)

tpm <- read_rds(paste0(download_dir_TM, "/tpm_WCDT-MCRPC", ".RDS"))
TM_genes_temp <- extractAllGenes(tpm, "tpm_unstrand")
write_csv(TM_genes_temp, file="PRAD_data/PRAD-TM_genes_WCDT-MCRPC.csv")
TM_genes_dat <- bind_rows(TM_genes_dat, TM_genes_temp)

write_csv(TM_genes_dat, "PRAD_data/PRAD-TM_genes.csv")






