
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("GEOquery", "edgeR", "limma", "impute", "preprocessCore")
cran_pkgs <- c("WGCNA", "matrixStats")
extra_bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Attempting to install Bioconductor package:", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Attempting to install CRAN package:", pkg))
    install.packages(pkg)
  }
}

for (pkg in extra_bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Attempting to install additional Bioconductor package:", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(GEOquery)
library(edgeR)
library(limma)
library(WGCNA)
library(matrixStats)

allowWGCNAThreads()

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)

## ---------------------------
## 1) Descargar raw counts desde GEO
## ---------------------------

# Specify the GEO accession ID for your dataset
geo_accession_id <- "GSE68719"

# Download supplementary files for the specified GEO accession ID
# This function will download files into a directory named after the geo_accession_id
getGEOSuppFiles(geo_accession_id, makeDirectory = TRUE)

# Construct the path to the correct count file within the downloaded directory
# Based on the `list.files` output in the previous execution,
# 'GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz' seems to contain the counts.
count_file <- file.path(geo_accession_id, paste0(geo_accession_id, "_mlpd_PCG_DESeq2_norm_counts.txt.gz"))

# Read the raw counts into a data frame
counts_df <- tryCatch({
  read.delim(
    gzfile(count_file),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  message("Could not read the expected count file. Please check the directory contents.")
  message("Error details: ", e$message)
  # List files in the directory to help user find the correct file
  message("Files in the directory:")
  print(list.files(file.path(geo_accession_id)))
  return(NULL)
})


# Quick inspection if data was loaded successfully
if (!is.null(counts_df)) {
  print(dim(counts_df))
  print(colnames(counts_df)[1:min(10, ncol(counts_df))])
  print(head(counts_df[, 1:min(6, ncol(counts_df))]))
} else {
  message("Data frame 'counts_df' is NULL. Please manually locate the raw counts file and update the 'count_file' variable.")
}


# Extraer los símbolos de los genes (segunda columna) para usarlos como IDs de fila
gene_symbols <- counts_df[[2]]

# Asegurarse de que los símbolos de los genes sean únicos
unique_gene_symbols <- make.unique(gene_symbols)

# Crear la matriz de conteos, eliminando las dos primeras columnas (EnsemblID y symbol)
count_matrix <- as.matrix(counts_df[, -c(1, 2)])
rownames(count_matrix) <- unique_gene_symbols # Usar símbolos de genes únicos aquí
mode(count_matrix) <- "numeric"

# Verificaciones
stopifnot(is.matrix(count_matrix))
stopifnot(ncol(count_matrix) == 73) # Se espera 73 columnas después de remover 2 de las 75 originales

print(dim(count_matrix))
print(count_matrix[1:5, 1:min(5, ncol(count_matrix))])


## ---------------------------
## 4) Metadata de muestras (adaptado para Parkinson/Control)
## ---------------------------

# Get sample names from the count_matrix columns
samples <- colnames(count_matrix)

# Determine disease status based on sample name prefix
disease_status <- ifelse(startsWith(samples, "P"), "Parkinson",
                         ifelse(startsWith(samples, "C"), "Control", NA))

# Create the metadata data frame
meta <- data.frame(
  sample = samples,
  disease_status = disease_status,
  row.names = samples,
  stringsAsFactors = FALSE
)

# Convert disease_status to a factor with desired levels
meta$disease_status <- factor(meta$disease_status, levels = c("Control", "Parkinson"))

# Display the metadata
print(meta)

# Display a summary table of the disease status
print(table(meta$disease_status, useNA = "ifany"))


set.seed(123) # for reproducibility

# Identify samples for each group
parkinson_samples <- rownames(meta[meta$disease_status == "Parkinson", ])
control_samples <- rownames(meta[meta$disease_status == "Control", ])

# Determine the number of samples to keep for each group
target_n <- 29

# Subsample Parkinson samples if there are more than target_n
if (length(parkinson_samples) > target_n) {
  parkinson_samples_balanced <- sample(parkinson_samples, target_n)
} else {
  parkinson_samples_balanced <- parkinson_samples
}

# Subsample Control samples if there are more than target_n
if (length(control_samples) > target_n) {
  control_samples_balanced <- sample(control_samples, target_n)
} else {
  control_samples_balanced <- control_samples
}

# Combine the balanced sample names
balanced_sample_names <- c(parkinson_samples_balanced, control_samples_balanced)

# Filter count_matrix to include only the balanced samples
count_matrix <- count_matrix[, balanced_sample_names]

# Filter metadata to include only the balanced samples
meta <- meta[balanced_sample_names, ]

# Verify the new dimensions and sample counts
print(paste("New dimensions of count_matrix:", paste(dim(count_matrix), collapse = " x ")))
print("New disease status distribution:")
print(table(meta$disease_status, useNA = "ifany"))

## ---------------------------
## 5) Crear objeto edgeR
## ---------------------------
y <- DGEList(counts = count_matrix, samples = meta)

# Chequeos
class(y$counts)
dim(y$counts)
head(y$samples)


design <- model.matrix(~ disease_status, data = meta)
design
colnames(design)

## ---------------------------
## 7) Filtrado de genes poco expresados
## ---------------------------
keep <- filterByExpr(y, design = design, min.total.count = 15)
y <- y[keep, , keep.lib.sizes = FALSE]

dim(y)

## ---------------------------
## 8) Normalización TMM
## ---------------------------
y <- calcNormFactors(y, method = "TMM")
y$samples


## ---------------------------
## 9) logCPM para WGCNA
## ---------------------------
logCPM <- cpm(y, log = TRUE, prior.count = 2)

dim(logCPM)
logCPM[1:5, 1:5]


## ---------------------------
## 10) Filtrado por variabilidad
## ---------------------------
gene_var <- rowVars(logCPM)

n_top <- min(5000, nrow(logCPM))
top_idx <- order(gene_var, decreasing = TRUE)[seq_len(n_top)]

expr_filt <- logCPM[top_idx, ]

## ---------------------------
## 11) Matriz para WGCNA
##     filas = muestras, columnas = genes
## ---------------------------
datExpr <- t(expr_filt)

# Verificación
dim(datExpr)
head(rownames(datExpr))
head(colnames(datExpr))

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

## ---------------------------
## 12) Traits para heatmap módulo-rasgo
## ---------------------------
traitData <- data.frame(
  parkinson = ifelse(meta[rownames(datExpr), "disease_status"] == "Parkinson", 1, 0),
  control = ifelse(meta[rownames(datExpr), "disease_status"]  == "Control", 1, 0)
)
rownames(traitData) <- rownames(datExpr)

#group_dummies <- model.matrix(~ 0 + disease_status, data = meta[rownames(datExpr), ])
traitData_full <- cbind(traitData)

traitData_full


## ---------------------------
## 13) Clustering de muestras
## ---------------------------
sampleTree <- hclust(dist(datExpr), method = "average")

pdf("01_sample_clustering.pdf", width = 8, height = 5)
par(cex = 0.8, mar = c(5, 4, 2, 2))
plot(sampleTree, main = "Clustering de muestras", xlab = "", sub = "")
dev.off()


## ---------------------------
## 14) Elegir soft-threshold
## ---------------------------
powers <- c(1:10, seq(12, 20, by = 2))

sft <- pickSoftThreshold(
  datExpr,
  dataIsExpr = TRUE,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(use = "p", maxPOutliers = 0.1),
  powerVector = powers,
  verbose = 5
)

pdf("02_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Potencia",
     ylab = "Scale-free topology fit (signed R^2)",
     type = "n",
     main = "Selección de potencia")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     cex = 0.9)
abline(h = 0.80, col = "red", lty = 2)

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Potencia",
     ylab = "Conectividad media",
     type = "n",
     main = "Conectividad media")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = 0.9)

dev.off()

candidate_powers <- sft$fitIndices[, 1][
  (-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]) >= 0.80
]
softPower <- if (length(candidate_powers) > 0) candidate_powers[1] else 6
message("Potencia elegida: ", softPower)

## ---------------------------
## 15) Construcción de red y módulos
## ---------------------------
net <- blockwiseModules(
  datExpr,
  power = softPower,
  corType = "bicor",
  maxPOutliers = 0.1,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- net$colors
table(moduleColors)

pdf("03_gene_dendrogram_modules.pdf", width = 12, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Dendrograma de genes y módulos"
)
dev.off()

## ---------------------------
## 16) Asociación módulo-rasgo
## ---------------------------
MEs <- orderMEs(net$MEs)

moduleTraitCor <- cor(MEs, traitData_full, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

textMatrix <- paste0(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 1),
  ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

pdf("04_module_trait_heatmap.pdf", width = 10, height = 7)
par(mar = c(8, 10, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(moduleTraitCor),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Asociación módulo-rasgo"
)
dev.off()


## ---------------------------
## 17) GS y MM usando Parkinson como ejemplo
## ---------------------------
GS.Parkinson <- as.data.frame(cor(datExpr, traitData$parkinson, use = "p"))
names(GS.Parkinson) <- "GS.Parkinson"

GSPvalue.Parkinson <- as.data.frame(corPvalueStudent(as.matrix(GS.Parkinson), nSamples = nrow(datExpr)))
names(GSPvalue.Parkinson) <- "p.GS.Parkinson"

MM <- as.data.frame(cor(datExpr, MEs, use = "p"))
MM.pvalue <- as.data.frame(corPvalueStudent(as.matrix(MM), nSamples = nrow(datExpr)))

bestME <- rownames(moduleTraitCor)[which.max(abs(moduleTraitCor[, "parkinson"]))]
moduleOfInterest <- sub("^ME", "", bestME)

message("Módulo más asociado a Parkinson: ", moduleOfInterest)

moduleGenes <- moduleColors == moduleOfInterest

hubTable <- data.frame(
  Gene = colnames(datExpr)[moduleGenes],
  Module = moduleColors[moduleGenes],
  GS_Parkinson = GS.Parkinson[moduleGenes, "GS.Parkinson"],
  GS_Parkinson_p = GSPvalue.Parkinson[moduleGenes, "p.GS.Parkinson"],
  MM = MM[moduleGenes, bestME],
  MM_p = MM.pvalue[moduleGenes, bestME]
)

hubTable <- hubTable[order(abs(hubTable$MM), abs(hubTable$GS_Parkinson), decreasing = TRUE), ]

write.csv(
  hubTable,
  file = paste0("05_hub_table_", moduleOfInterest, "_Parkinson.csv"),
  row.names = FALSE
)

pdf(paste0("06_GS_vs_MM_", moduleOfInterest, "_Parkinson.pdf"), width = 6, height = 6)
verboseScatterplot(
  abs(hubTable$MM),
  abs(hubTable$GS_Parkinson),
  xlab = paste("Module Membership en", moduleOfInterest),
  ylab = "Gene Significance para Parkinson",
  main = paste("GS vs MM - módulo", moduleOfInterest),
  cex = 1.0,
  col = moduleOfInterest
)
dev.off()
## -----------------------------------
## 18) Red de coexpresión solo Parkinson
## -----------------------------------

# 1. Filtrar datExpr y meta para solo muestras de Parkinson
parkinson_sample_names <- rownames(meta[meta$disease_status == "Parkinson", ])
datExpr_parkinson <- datExpr[parkinson_sample_names, ]
meta_parkinson <- meta[parkinson_sample_names, ]

message("Dimensiones de datExpr_parkinson: ", paste(dim(datExpr_parkinson), collapse = " x "))

# 2. Elegir soft-threshold para las muestras de Parkinson
powers_parkinson <- c(1:10, seq(12, 20, by = 2))
sft_parkinson <- pickSoftThreshold(
  datExpr_parkinson,
  dataIsExpr = TRUE,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(use = "p", maxPOutliers = 0.1),
  powerVector = powers_parkinson,
  verbose = 5
)

pdf("07_soft_threshold_parkinson.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sft_parkinson$fitIndices[, 1],
     -sign(sft_parkinson$fitIndices[, 3]) * sft_parkinson$fitIndices[, 2],
     xlab = "Potencia",
     ylab = "Scale-free topology fit (signed R^2)",
     type = "n",
     main = "Selección de potencia (Parkinson)")
text(sft_parkinson$fitIndices[, 1],
     -sign(sft_parkinson$fitIndices[, 3]) * sft_parkinson$fitIndices[, 2],
     labels = powers_parkinson,
     cex = 0.9)
abline(h = 0.80, col = "red", lty = 2)

plot(sft_parkinson$fitIndices[, 1],
     sft_parkinson$fitIndices[, 5],
     xlab = "Potencia",
     ylab = "Conectividad media",
     type = "n",
     main = "Conectividad media (Parkinson)")
text(sft_parkinson$fitIndices[, 1],
     sft_parkinson$fitIndices[, 5],
     labels = powers_parkinson,
     cex = 0.9)

dev.off()

candidate_powers_parkinson <- sft_parkinson$fitIndices[, 1][
  (-sign(sft_parkinson$fitIndices[, 3]) * sft_parkinson$fitIndices[, 2]) >= 0.80
]
softPower_parkinson <- if (length(candidate_powers_parkinson) > 0) candidate_powers_parkinson[1] else 6
message("Potencia elegida para Parkinson: ", softPower_parkinson)

# 3. Construcción de red y módulos para Parkinson
net_parkinson <- blockwiseModules(
  datExpr_parkinson,
  power = softPower_parkinson,
  corType = "bicor",
  maxPOutliers = 0.1,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors_parkinson <- net_parkinson$colors
table(moduleColors_parkinson)

pdf("08_gene_dendrogram_modules_parkinson.pdf", width = 12, height = 6)
plotDendroAndColors(
  net_parkinson$dendrograms[[1]],
  moduleColors_parkinson[net_parkinson$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Dendrograma de genes y módulos (Parkinson)"
)
dev.off()

# Opcional: Calcular MEs para la red de Parkinson
MEs_parkinson <- orderMEs(net_parkinson$MEs)

message("Red de coexpresión para Parkinson construida. Número de módulos: ", length(table(moduleColors_parkinson)))

## -----------------------------------
## 19) Red de coexpresión solo Control
## -----------------------------------

# 1. Filtrar datExpr y meta para solo muestras de Control
control_sample_names <- rownames(meta[meta$disease_status == "Control", ])
datExpr_control <- datExpr[control_sample_names, ]
meta_control <- meta[control_sample_names, ]

message("Dimensiones de datExpr_control: ", paste(dim(datExpr_control), collapse = " x "))

# 2. Elegir soft-threshold para las muestras de Control
powers_control <- c(1:10, seq(12, 20, by = 2))
sft_control <- pickSoftThreshold(
  datExpr_control,
  dataIsExpr = TRUE,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(use = "p", maxPOutliers = 0.1),
  powerVector = powers_control,
  verbose = 5
)

pdf("09_soft_threshold_control.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sft_control$fitIndices[, 1],
     -sign(sft_control$fitIndices[, 3]) * sft_control$fitIndices[, 2],
     xlab = "Potencia",
     ylab = "Scale-free topology fit (signed R^2)",
     type = "n",
     main = "Selección de potencia (Control)")
text(sft_control$fitIndices[, 1],
     -sign(sft_control$fitIndices[, 3]) * sft_control$fitIndices[, 2],
     labels = powers_control,
     cex = 0.9)
abline(h = 0.80, col = "red", lty = 2)

plot(sft_control$fitIndices[, 1],
     sft_control$fitIndices[, 5],
     xlab = "Potencia",
     ylab = "Conectividad media",
     type = "n",
     main = "Conectividad media (Control)")
text(sft_control$fitIndices[, 1],
     sft_control$fitIndices[, 5],
     labels = powers_control,
     cex = 0.9)

dev.off()

candidate_powers_control <- sft_control$fitIndices[, 1][
  (-sign(sft_control$fitIndices[, 3]) * sft_control$fitIndices[, 2]) >= 0.80
]
softPower_control <- if (length(candidate_powers_control) > 0) candidate_powers_control[1] else 6
message("Potencia elegida para Control: ", softPower_control)

# 3. Construcción de red y módulos para Control
net_control <- blockwiseModules(
  datExpr_control,
  power = softPower_control,
  corType = "bicor",
  maxPOutliers = 0.1,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors_control <- net_control$colors
table(moduleColors_control)

pdf("10_gene_dendrogram_modules_control.pdf", width = 12, height = 6)
plotDendroAndColors(
  net_control$dendrograms[[1]],
  moduleColors_control[net_control$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Dendrograma de genes y módulos (Control)"
)
dev.off()

# Opcional: Calcular MEs para la red de Control
MEs_control <- orderMEs(net_control$MEs)

message("Red de coexpresión para Control construida. Número de módulos: ", length(table(moduleColors_control)))


## -----------------------------------------------------------## 18) Exportar módulos a Cytoscape para redes de Parkinson y Control## -----------------------------------------------------------# Función para exportar módulos a Cytoscape para una red específica
exportModulesToCytoscape <- function(datExpr_net, net_obj, moduleColors_net, softPower_net, network_label_prefix) {
  message(paste("Calculando TOM para la red", network_label_prefix, "..."))
  # Calcular TOM para toda la red si no está ya guardada en net_obj
  # (blockwiseModules con saveTOMs=FALSE no guarda el TOM completo)
  TOM_net <- TOMsimilarityFromExpr(
    datExpr_net,
    corType = "bicor",
    maxPOutliers = 0.1,
    networkType = "signed",
    TOMType = "signed",
    power = softPower_net,
    verbose = 0 # Silenciar salida detallada para este cálculo si no es necesario
  )
  dimnames(TOM_net) <- list(colnames(datExpr_net), colnames(datExpr_net))

  # Obtener colores de módulos únicos (excluyendo 'grey' para genes no asignados)
  unique_modules <- unique(moduleColors_net)
  unique_modules <- unique_modules[unique_modules != "grey"]

  message(paste("Exportando", length(unique_modules), "módulos para la red", network_label_prefix, "..."))

  # Crear directorio para los archivos de Cytoscape de esta red
  output_dir <- paste0(network_label_prefix, "_Cytoscape_Exports")
  dir.create(output_dir, showWarnings = FALSE) # showWarnings = FALSE para no advertir si ya existe

  for (module_color in unique_modules) {
    message(paste("  Exportando módulo:", module_color, "de la red", network_label_prefix))
    # Seleccionar genes en el módulo actual
    moduleGenes <- (moduleColors_net == module_color)

    if (sum(moduleGenes) == 0) {
      message(paste("    Módulo", module_color, "está vacío, saltando."))
      next
    }

    # Extraer la submatriz TOM para este módulo
    TOM_mod <- TOM_net[moduleGenes, moduleGenes]

    # Crear atributos de nodo para Cytoscape
    # Incluimos solo el color del módulo ya que GS/MM requerirían recalcularse para sub-redes con un trait específico.
    nodeAttr <- data.frame(
      module = rep(module_color, sum(moduleGenes))
    )
    # Asegurar que los nombres de los nodos sean únicos antes de asignarlos como nombres de fila
    rownames(nodeAttr) <- make.unique(colnames(TOM_mod))

    # Definir nombres de archivos con la ruta de la carpeta
    edge_file <- file.path(output_dir, paste0("cytoscape_edges_", network_label_prefix, "_", module_color, "_.txt"))
    node_file <- file.path(output_dir, paste0("cytoscape_nodes_", network_label_prefix, "_", module_color, "_.txt"))

    # Exportar a Cytoscape
    exportNetworkToCytoscape(
      TOM_mod,
      weighted = TRUE,
      threshold = 0.10, # Usando el umbral original
      nodeNames = colnames(TOM_mod),
      altNodeNames = colnames(TOM_mod),
      nodeAttr = nodeAttr,
      edgeFile = edge_file,
      nodeFile = node_file
    )
  }
}

# Exportar módulos para la red de Parkinson
exportModulesToCytoscape(datExpr_parkinson, net_parkinson, moduleColors_parkinson, softPower_parkinson, "Parkinson")

# Exportar módulos para la red de Control
exportModulesToCytoscape(datExpr_control, net_control, moduleColors_control, softPower_control, "Control")

## ----------------------------------------------------## 19) Guardar resultados para redes de Parkinson y Control## ----------------------------------------------------
# Crear directorios para los resultados de WGCNA
dir.create("Parkinson_WGCNA_Results", showWarnings = FALSE)
dir.create("Control_WGCNA_Results", showWarnings = FALSE)

# Guardar resultados para la red de Parkinson
message("Guardando resultados de WGCNA para la red de Parkinson...")
saveRDS(list(
  datExpr = datExpr_parkinson,
  meta = meta_parkinson,
  sft = sft_parkinson,
  net = net_parkinson,
  moduleColors = moduleColors_parkinson,
  MEs = MEs_parkinson), file = file.path("Parkinson_WGCNA_Results", "08_wgcna_Parkinson_results_GSE290452.rds"))

# Guardar resultados para la red de Control
message("Guardando resultados de WGCNA para la red de Control...")
saveRDS(list(
  datExpr = datExpr_control,
  meta = meta_control,
  sft = sft_control,
  net = net_control,
  moduleColors = moduleColors_control,
  MEs = MEs_control), file = file.path("Control_WGCNA_Results", "08_wgcna_Control_results_GSE290452.rds"))


## ---------------------------
## 18A) Enriquecimiento funcional de módulos (Red Global)
## ---------------------------

# Crear carpeta para resultados
dir.create("09_module_enrichment", showWarnings = FALSE)

# Genes usados en WGCNA (de la red global)
all_genes <- colnames(datExpr)

# Limpiar posibles versiones de Ensembl, por ejemplo ENSMUSG00000000001.1 -> ENSMUSG00000000001
all_genes_clean <- sub("\\..*$", "", all_genes)

# Inferir tipo de identificador
if (all(grepl("^ENSMUSG", all_genes_clean))) {
  keytype_used <- "ENSEMBL"
} else if (all(grepl("^[0-9]+$", all_genes_clean))) {
  keytype_used <- "ENTREZID"
} else {
  keytype_used <- "SYMBOL"
}

message("Tipo de identificador inferido: ", keytype_used)

# Tabla de mapeo gene ID -> Entrez ID
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(all_genes_clean),
  keytype = keytype_used,
  columns = c("ENTREZID", "SYMBOL")
)

gene_map <- gene_map[!is.na(gene_map$ENTREZID), ]
gene_map <- gene_map[!duplicated(gene_map[, keytype_used]), ]

# Universo: todos los genes usados para construir la red
# Esto es importante para que el enriquecimiento use como fondo los genes analizados,
# no todos los genes del genoma.
background_entrez <- unique(gene_map$ENTREZID)

# Crear una tabla gen-módulo (usando moduleColors de la red global)
gene_module_df <- data.frame(
  Gene_original = all_genes,
  Gene_clean = all_genes_clean,
  Module = moduleColors,
  stringsAsFactors = FALSE
)

# Agregar Entrez ID y símbolo
gene_module_df <- merge(
  gene_module_df,
  gene_map,
  by.x = "Gene_clean",
  by.y = keytype_used,
  all.x = TRUE
)

# Guardar tabla de genes por módulo
write.csv(
  gene_module_df,
  file = "09_module_enrichment/09_gene_module_annotation_global_network.csv",
  row.names = FALSE
)

# Lista de genes Entrez por módulo
module_gene_list <- split(
  gene_module_df$ENTREZID,
  gene_module_df$Module
)

# Eliminar NA y módulos muy pequeños
module_gene_list <- lapply(module_gene_list, function(x) unique(na.omit(x)))
module_gene_list <- module_gene_list[sapply(module_gene_list, length) >= 10]

# Opcional: excluir módulo grey, porque en WGCNA suele representar genes no asignados
module_gene_list <- module_gene_list[names(module_gene_list) != "grey"]

# Enriquecimiento GO Biological Process por módulo
ego_modules <- compareCluster(
  geneCluster = module_gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  universe = background_entrez,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10,
  readable = TRUE
)

# Convertir a tabla
ego_modules_df <- as.data.frame(ego_modules)

write.csv(
  ego_modules_df,
  file = "09_module_enrichment/10_GO_BP_enrichment_all_modules_global_network.csv",
  row.names = FALSE
)

# Dotplot comparativo de módulos
pdf("09_module_enrichment/11_GO_BP_dotplot_all_modules_global_network.pdf",
    width = 12, height = 8)

print(
  dotplot(ego_modules, showCategory = 10) +
    ggplot2::ggtitle("Enriquecimiento GO Biological Process por módulo (Red Global)")
)

dev.off()

## -----------------------------------------------------------------
## 19) Enriquecimiento del módulo más asociado a Parkinson (Red Global)
## -----------------------------------------------------------------

# Identificar el módulo más asociado al rasgo 'parkinson'
bestME <- rownames(moduleTraitCor)[which.max(abs(moduleTraitCor[, "parkinson"]))]
moduleOfInterest <- sub("^ME", "", bestME)

message("Módulo de interés para enriquecimiento (asociado a Parkinson): ", moduleOfInterest)

genes_interest <- gene_module_df$ENTREZID[
  gene_module_df$Module == moduleOfInterest
]

genes_interest <- unique(na.omit(genes_interest))

ego_interest <- enrichGO(
  gene = genes_interest,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  universe = background_entrez,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10,
  readable = TRUE
)

ego_interest_df <- as.data.frame(ego_interest)

write.csv(
  ego_interest_df,
  file = paste0(
    "09_module_enrichment/12_GO_BP_enrichment_module_",
    moduleOfInterest,
    "_Parkinson_global_network.csv"
  ),
  row.names = FALSE
)

pdf(
  paste0(
    "09_module_enrichment/13_GO_BP_dotplot_module_",
    moduleOfInterest,
    "_Parkinson_global_network.pdf"
  ),
  width = 8,
  height = 6
)

print(
  dotplot(ego_interest, showCategory = 15) +
    ggplot2::ggtitle(
      paste("GO BP enrichment - módulo", moduleOfInterest, "(asociado a Parkinson - Red Global)")
    )
)

dev.off()

## ----------------------------------------------------
## 18) Guardar resultados de la Red Global
## ----------------------------------------------------

# Crear directorio para los resultados de la red global
dir.create("Global_WGCNA_Results", showWarnings = FALSE)


# Función para exportar módulos a Cytoscape para una red específica
exportModulesToCytoscape <- function(datExpr_net, net_obj, moduleColors_net, softPower_net, network_label_prefix) {
  message(paste("Calculando TOM para la red", network_label_prefix, "..."))
  # Calcular TOM para toda la red si no está ya guardada en net_obj
  # (blockwiseModules con saveTOMs=FALSE no guarda el TOM completo)
  TOM_net <- TOMsimilarityFromExpr(
    datExpr_net,
    corType = "bicor",
    maxPOutliers = 0.1,
    networkType = "signed",
    TOMType = "signed",
    power = softPower_net,
    verbose = 0 # Silenciar salida detallada para este cálculo si no es necesario
  )
  dimnames(TOM_net) <- list(colnames(datExpr_net), colnames(datExpr_net))

  # Obtener colores de módulos únicos (excluyendo 'grey' para genes no asignados)
  unique_modules <- unique(moduleColors_net)
  unique_modules <- unique_modules[unique_modules != "grey"]

  message(paste("Exportando", length(unique_modules), "módulos para la red", network_label_prefix, "..."))

  # Crear directorio para los archivos de Cytoscape de esta red
  output_dir <- paste0(network_label_prefix, "_Cytoscape_Exports")
  dir.create(output_dir, showWarnings = FALSE) # showWarnings = FALSE para no advertir si ya existe

  for (module_color in unique_modules) {
    message(paste("  Exportando módulo:", module_color, "de la red", network_label_prefix))
    # Seleccionar genes en el módulo actual
    moduleGenes <- (moduleColors_net == module_color)

    if (sum(moduleGenes) == 0) {
      message(paste("    Módulo", module_color, "está vacío, saltando."))
      next
    }

    # Extraer la submatriz TOM para este módulo
    TOM_mod <- TOM_net[moduleGenes, moduleGenes]

    # Crear atributos de nodo para Cytoscape
    # Incluimos solo el color del módulo ya que GS/MM requerirían recalcularse para sub-redes con un trait específico.
    nodeAttr <- data.frame(
      module = rep(module_color, sum(moduleGenes))
    )
    # Asegurar que los nombres de los nodos sean únicos antes de asignarlos como nombres de fila
    rownames(nodeAttr) <- make.unique(colnames(TOM_mod))

    # Definir nombres de archivos con la ruta de la carpeta
    edge_file <- file.path(output_dir, paste0("cytoscape_edges_", network_label_prefix, "_", module_color, "_.txt"))
    node_file <- file.path(output_dir, paste0("cytoscape_nodes_", network_label_prefix, "_", module_color, "_.txt"))

    # Exportar a Cytoscape
    exportNetworkToCytoscape(
      TOM_mod,
      weighted = TRUE,
      threshold = 0.10, # Usando el umbral original
      nodeNames = colnames(TOM_mod),
      altNodeNames = colnames(TOM_mod),
      nodeAttr = nodeAttr,
      edgeFile = edge_file,
      nodeFile = node_file
    )
  }
}

# Exportar módulos para la red global
exportModulesToCytoscape(datExpr, net, moduleColors, softPower, "Global")

message("Guardando resultados de WGCNA para la red global...")
saveRDS(list(
  y = y,
  meta = meta,
  logCPM = logCPM,
  datExpr = datExpr,
  traitData = traitData_full,
  sft = sft,
  net = net,
  moduleColors = moduleColors,
  MEs = MEs,
  moduleTraitCor = moduleTraitCor,
  moduleTraitPvalue = moduleTraitPvalue,
  hubTable = hubTable,
  moduleOfInterest = moduleOfInterest
), file = file.path("Global_WGCNA_Results", "08_wgcna_Global_results_GSE68719.rds"))


## ---------------------------
## 20) Resumen General del Análisis WGCNA
## ---------------------------

cat("\n### Resumen del Análisis de Co-expresión de Genes (WGCNA) ###\n")
cat("-----------------------------------------------------------\n")

cat("\n#### Red de Co-expresión Global (Parkinson y Control Combinados) ####\n")
cat("-----------------------------------------------------------\n")
cat("Número de muestras en la red global: ", nrow(datExpr), "\n")
cat("Número de genes en la red global: ", ncol(datExpr), "\n")
cat("Potencia de soft-thresholding elegida para la red global: ", softPower, "\n")
cat("Módulos detectados en la red global:\n")
print(table(moduleColors))
cat("Módulo más asociado a Parkinson (red global): ", moduleOfInterest, "\n")
cat("Top 10 genes hub del módulo ", moduleOfInterest, " (red global):\n")
print(head(hubTable, 10))

cat("\n#### Red de Co-expresión Específica de Parkinson ####\n")
cat("-----------------------------------------------------\n")
cat("Número de muestras en la red Parkinson: ", nrow(datExpr_parkinson), "\n")
cat("Número de genes en la red Parkinson: ", ncol(datExpr_parkinson), "\n")
cat("Potencia de soft-thresholding elegida para Parkinson: ", softPower_parkinson, "\n")
cat("Módulos detectados en la red Parkinson:\n")
print(table(moduleColors_parkinson))

cat("\n#### Red de Co-expresión Específica de Control ####\n")
cat("---------------------------------------------------\n")
cat("Número de muestras en la red Control: ", nrow(datExpr_control), "\n")
cat("Número de genes en la red Control: ", ncol(datExpr_control), "\n")
cat("Potencia de soft-thresholding elegida para Control: ", softPower_control, "\n")
cat("Módulos detectados en la red Control:\n")
print(table(moduleColors_control))

cat("\n#### Análisis de Enriquecimiento Funcional ####\n")
cat("-------------------------------------------\n")
cat("Se realizó el enriquecimiento de términos GO Biological Process para los módulos de la red global.\n")
cat("Los resultados detallados se encuentran en la carpeta '09_module_enrichment'.\n")
cat("El módulo '", moduleOfInterest, "' (asociado a Parkinson en la red global) también fue analizado individualmente para enriquecimiento.\n")

cat("\n-----------------------------------------------------------\n")