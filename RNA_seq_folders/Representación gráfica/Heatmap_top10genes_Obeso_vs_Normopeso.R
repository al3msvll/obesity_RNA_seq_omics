#Estalecer directorio de trabajo
setwd("C:/Users/albar/Downloads/UNIR/Secuenciación y Ómicas de Próxima Generación/mubio03_act2")

#Descarga de las librerías necesarias
if (!require("BiocManager", quietly = TRUE)) # Verificamos si el paquete está instalado sin cargarlo,
  install.packages("BiocManager") # Si no está, se instala.
BiocManager::install("DESeq2", force = TRUE) # Instalamos el programa que vamos a utilizar para realizar el análisis diferencial de genes
library(DESeq2)
library(ggplot2)
library(pheatmap)

# 1. CARGAR Y PREPARAR DATOS
# Cargamos el archivo, usamos 'row.names = 1' para decirle a R que la primera columna es la que contiene los nombres de genes.
df <- read.csv("DESeq2_resultados_Obeso1_vs_Normopeso.csv", row.names = 1, header = TRUE, check.names = FALSE)
vsd <- read.csv("vsd_normalizado.csv", row.names = 1, header = TRUE, check.names = FALSE)

# Ordenamos por padj (p-valor ajustado), nuestro p-valor estándar es 0.05 y seleccionamos los 10 primeros genes
# Eliminamos cualquier valor NA en padj antes de ordenar
df_results_orden <- df[!is.na(df$padj), ]
df_results_orden <- df_results_orden[order(df_results_orden$padj, decreasing = FALSE), ]
top_10_gene_names <- rownames(df_results_orden)[1:10]

# Filtramos la matriz VSD para solo los 10 genes seleccionados
mat_top10 <- as.matrix(vsd[top_10_gene_names, ])


#Comenzamos definiendo la condición de los 5 individuos
sample_groups <- data.frame(
  Grupo = factor(c("Obeso", "Obeso", "Normopeso", "Normopeso", "Normopeso"))
)
rownames(sample_groups) <- colnames(data_matrix)

# 2. HEATMAP DIFERENCIAL
#Realizamos un heatmap de la expresión génica de los 10 genes con mayor diferencia entre los dos grupos
mat_top10_escalada <- t(scale(t(mat_top10)))

pheatmap(
  mat_top10_escalada,
  scale = "row", 
  annotation_col = sample_groups,
  color = colorRampPalette(c("blue", "white", "red"))(100), 
  fontsize = 8,
  main = "Heatmap de los 10 Genes Más Significativos",
  treeheight_row = 50,
  treeheight_col = 30,
  filename = "Heatmap_top10genes_Obeso_vs_Normopeso.pdf" 
)
