#Set work directory
setwd("/Volumes/Macintosh HD/Users/anasofia/Master_UNIR/Primer_trimestre/Secuenciacion_Omicas/Actividad2")

#Instalar tximport
#BiocManager::install("tximport")

#Librerias necesarias
library(DESeq2)
library(tximport)
library(tidyverse)

#Para hacer la matriz de conteos
samples <- c("AbrahamSimpson", "HomerSimpson", "BartSimpson", "LisaSimpson", "MaggieSimpson")

#Crear el path para o ficheiro quant.sf
files <- file.path(
  paste0(samples, "_quant_creado"),
  "quant.sf"
)

#Dar nombre a los ficheros
names(files) <- samples

#Carga el archivo transcrito → gen
tx2gene <- read.delim(
  "Transcrito_a_Gen.tsv",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(tx2gene) <- c("transcript_id", "gene_id")
str(tx2gene)

#Importar los ficheros quant.sf
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene
)

#Guardar matriz de conteos
write.table(
  txi$counts,
  file = "matriz_conteos_genes.tsv",
  sep = "\t",
)

#Normalización de los datos
#Cargar el CSV en R
conteos <- read.csv2("Matriz_de_conteos.csv", sep = ";", dec = ".")

#Crear la tabla experimental
diseño <- read.csv2("Design.csv", sep = ",")

#Filtrar solo para las muestras que queremos
diseño_grupo8 <- diseño %>% 
  filter(Condition %in% c("Normopeso", "Sobrepeso/Obeso1"))

#Cambiar el el nombre de una de las observciones dentro de una variabel
diseño_grupo8 <- diseño_grupo8 %>%
  mutate(Condition = recode(Condition,
                            "Sobrepeso/Obeso1" = "Obeso1"))

#Colocar los nombres de las filas en la tabla de disño
rownames(diseño_grupo8) <- diseño_grupo8$Sample

#Eliminar la primera columna que lleva los nombre de las filas
diseño_grupo8 <- diseño_grupo8[, -1]

#Colocar los nombres de las filas en la table de los conteos
rownames(conteos) <- conteos$Gene.ID

#Eliminar la primera columna que lleva los nombre de las filas
conteos <- conteos[, -1]

#Redondear la matriz de conteos que si no no funciona DEseq2
conteos_final <- round(conteos)

#Comprobar que el nombre de las muestras es el mismo - tiene que dar lo mismo
colnames(conteos)
rownames(diseño_grupo8)

#Comprobar tipos
str(conteos)
str(diseño_grupo8)

#Asegurar factor y nivel de referencia
diseño_grupo8$Condition <- factor(
  diseño_grupo8$Condition,
  levels = c("Normopeso", "Obeso1")
)

#Crear el objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = conteos_final,
  colData = diseño_grupo8,
  design = ~ Condition
)

#Filtrado básico (RECOMENDADO) - Quita genes sin expresión real:
dds <- dds[rowSums(counts(dds)) > 10, ]

#Ejecutar DESeq2 (AQUÍ ocurre la normalización)
dds <- DESeq(dds)

# Obtener resultados Obeso1 vs Normopeso
resultados <- results(
  dds,
  contrast = c("Condition", "Obeso1", "Normopeso")
)

#Guardar resultados (muy recomendable)
write.csv(
  as.data.frame(resultados),
  file = "DESeq2_resultados_Obeso1_vs_Normopeso.csv"
)

#Conteos normalizados (opcional pero útil)
conteos_norm <- counts(dds, normalized = TRUE)

write.csv(
  conteos_norm,
  file = "conteos_normalizados_DESeq2.csv"
)

#Normalización para visualización (IMPORTANTE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vsd_mat <- assay(vsd)

#Guardar la tabla de normalización para visualización
write.csv(vsd_mat, "vsd_normalizado.csv")
