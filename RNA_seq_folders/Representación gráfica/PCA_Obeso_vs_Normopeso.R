# Instalar los paquetes si no los tienes (descomentar si es necesario)
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("ggrepel")

# Cargar librerías
library(ggplot2)
library(dplyr)
library(ggrepel)

# 1. CARGAR Y PREPARAR DATOS

# Leer el archivo, usando la primera columna como nombres de fila (genes)
vsd_data <- read.csv("vsd_normalizado.csv", row.names = 1)

# Transponer la matriz: PCA se realiza sobre las muestras (filas) y genes (columnas)
vsd_t <- t(vsd_data)

# -------------------------------------------------------------------
# 2. DEFINIR LA AGRUPACIÓN (¡AJUSTA ESTO A TUS METADATOS!)

sample_group <- data.frame(
  Sample = rownames(vsd_t),
  Group = factor(c("Normopeso", "Normopeso", "Obeso", "Obeso", "Obeso"))
)

# -------------------------------------------------------------------
# 3. REALIZAR EL PCA

pca_result <- prcomp(vsd_t, center = TRUE, scale. = FALSE)

# 4. PREPARAR DATOS PARA EL GRÁFICO

# A) Datos de Muestras (Scores)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- rownames(pca_scores)
pca_plot_data <- pca_scores %>%
  left_join(sample_group, by = "Sample")

# B) Datos de Cargas de Genes (Loadings/Rotations)
gene_loadings <- as.data.frame(pca_result$rotation)
gene_loadings$GeneName <- rownames(gene_loadings)

# 5. IDENTIFICAR LOS GENES MÁS CONTRIBUYENTES
# Para saber qué genes etiquetar, seleccionamos los que tienen mayor impacto en PC1 y PC2.
# Usamos la suma de las contribuciones al cuadrado a PC1 y PC2.
gene_loadings$Contribution <- gene_loadings$PC1^2 + gene_loadings$PC2^2

# Seleccionar los 10 genes con mayor contribución total
top_genes_for_plot <- gene_loadings %>%
  arrange(desc(Contribution)) %>%
  head(10)

# Calcular la varianza explicada por PC1 y PC2
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# 6. GENERAR EL PCA PLOT (Muestras + Etiquetas de Genes)

pca_plot_with_genes <- ggplot() +
  
  # 1. Puntos de Muestras (Muestras como puntos grandes, coloreadas por grupo)
  geom_point(data = pca_plot_data, aes(x = PC1, y = PC2, color = Group), size = 4, alpha = 0.8) +
  # 2. Etiquetas de Muestras
  geom_text_repel(data = pca_plot_data, aes(x = PC1, y = PC2, label = Sample), color = "black", size = 3) +
  
  # 3. Etiquetas de Genes (Los 10 más contribuyentes)
  # Usamos un punto 'x' o '+' para distinguir las genes de las muestras
  geom_point(data = top_genes_for_plot, aes(x = PC1*5, y = PC2*5), shape = 3, size = 3, color = "darkgreen") +
  # Nota: Multiplicamos por 5 las coordenadas del gen para escalarlas y que aparezcan en el mismo rango visual que las muestras.
  geom_text_repel(data = top_genes_for_plot, 
                  aes(x = PC1*5, y = PC2*5, label = GeneName), 
                  color = "darkgreen", 
                  size = 3, 
                  max.overlaps = 50) +
  
  # 4. Líneas de los ejes (opcional)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # 5. Etiquetas y Título
  xlab(paste0("PC1: ", percentVar[1], "% varianza explicada")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza explicada")) +
  ggtitle("PCA de Muestras y Contribución de Genes (vsd_normalizado)") +
  
  # Temas y formato
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Mostrar el gráfico
print(pca_plot_with_genes)

ggsave("PCA2_vsd_normalizado.png", plot = pca_plot, width = 8, height = 7)


