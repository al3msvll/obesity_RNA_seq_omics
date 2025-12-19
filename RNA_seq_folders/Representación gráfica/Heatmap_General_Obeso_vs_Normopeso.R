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
df <- read.csv("DESeq2_resultados_Obeso1_vs_Normopeso.csv", row.names = 1)
vsd <- read.csv("vsd_normalizado.csv", row.names = 1, header = TRUE, check.names = FALSE)

data_matrix <- as.matrix(vsd)

# 2. HEATMAP GENERAL
#Realizamos un heatmap de la expresión génica general
#Comenzamos definiendo la condición de los 5 individuos
sample_groups <- data.frame(
  Grupo = factor(c("Obeso", "Obeso", "Normopeso", "Normopeso", "Normopeso"))
)
rownames(sample_groups) <- colnames(data_matrix)
mat_escalada <- t(scale(t(data_matrix)))

pheatmap(
  mat_escalada,
  scale = "row",
  annotation_col = sample_groups,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize = 8,
  main = "Heatmap General de Expresión",
  treeheight_row = 50, 
  treeheight_col = 30, 
  filename = "Heatmap_General_Obeso_vs_Normopeso.pdf" 
)
