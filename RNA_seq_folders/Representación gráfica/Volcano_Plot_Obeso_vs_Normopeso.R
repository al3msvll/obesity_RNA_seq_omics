
# Instalar los paquetes si no los tienes (descomenta si es necesario)
# install.packages("ggplot2")
# install.packages("dplyr")
install.packages("ggrepel") # Necesario para etiquetar sin superposición

# Cargar las librerías
library(ggplot2)
library(dplyr)
library(ggrepel)

# 1. CARGAR Y PREPARAR DATOS
# Cargamos el archivo. Notar que usamos 'row.names = 1' para decirle a R que la primera
# columna es la que contiene los IDs (nombres de genes).
df <- read.csv("DESeq2_resultados_Obeso1_vs_Normopeso.csv", row.names = 1)

# El nombre de la columna se llama 'log2FoldChange' y 'padj' en el archivo, pero
# por si acaso, renombramos las columnas si no tienen el nombre correcto.
# (Este paso es opcional si ya tienen el nombre correcto)
# colnames(df)[colnames(df) == "log2FoldChange"] <- "Log2FC"
# colnames(df)[colnames(df) == "padj"] <- "padj"

# AHORA: Extraemos los nombres de fila (que son los genes) a una columna llamada GeneName
df$GeneName <- rownames(df)

# --- Definir Umbrales ---
LFC_THRESHOLD <- 1.0     # |Log2FC| > 1
P_ADJ_THRESHOLD <- 0.05  # padj < 0.05

# 2. PREPARAR COLUMNAS PARA PLOT
df <- df %>%
  # 1. Calcular -log10(padj) para el eje Y
  mutate(neg_log10_padj = -log10(padj)) %>%
  
  # 2. Categorizar los genes por regulación
  mutate(Regulation = case_when(
    log2FoldChange >= LFC_THRESHOLD & padj < P_ADJ_THRESHOLD ~ "Upregulated (Obeso > Normopeso)",
    log2FoldChange <= -LFC_THRESHOLD & padj < P_ADJ_THRESHOLD ~ "Downregulated (Obeso < Normopeso)",
    TRUE ~ "Not Significant"
  ))

# 3. Identificar los 5 genes más importantes para etiquetar
top_genes <- df %>%
  filter(Regulation != "Not Significant") %>%
  arrange(padj) %>%
  head(5)

# 4. GENERAR EL GRÁFICO
volcano_plot <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  # Puntos de dispersión (coloreados por regulación)
  geom_point(aes(color = Regulation), alpha = 0.8, size = 2) +
  
  # Líneas de umbral
  geom_hline(yintercept = -log10(P_ADJ_THRESHOLD), linetype = "dashed", color = "#828282") +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", color = "#828282") +
  
  # Etiquetas de los 5 genes principales (usando ggrepel para que no se superpongan)
  geom_text_repel(data = top_genes,
                  aes(label = GeneName),
                  size = 3,
                  max.overlaps = 20) +
  
  # Etiquetas y Título
  labs(
    title = expression(paste("Volcano Plot: Obeso vs. Normopeso (", Log[2], "FC > 1, padj < 0.05)")),
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"P-valor Ajustado (padj)"),
    color = "Regulación"
  ) +
  
  # Colores personalizados
  scale_color_manual(values = c(
    "Upregulated (Obeso > Normopeso)" = "red",
    "Downregulated (Obeso < Normopeso)" = "blue",
    "Not Significant" = "black"
  )) +
  
  theme_minimal() +
  guides(color = guide_legend(nrow = 2)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 9)
        )

# Mostrar el gráfico
print(volcano_plot)

# Guardar el gráfico como un archivo PNG
ggsave("Volcano_Plot_Obeso_vs_Normopeso.png", plot = volcano_plot, width = 10, height = 8)

