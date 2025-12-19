# obesity_RNA_seq_omics
# Análisis de Expresión Diferencial: Obesidad vs. Normopeso (RNA-seq)
Este repositorio contiene el flujo de trabajo en R para el análisis de datos transcriptómicos de muestras de tejido de pacientes con obesidad frente a pacientes con normopeso. El objetivo es identificar genes diferencialmente expresados (DEGs) y visualizar la variabilidad biológica entre grupos. El proyecto analiza perfiles metabólicos de personajes de Los Simpson para comprender la base genética de la obesidad.

## 🎯 Objetivos del Proyecto
- Realizar el control de calidad y alineamiento de lecturas de archivos FASTQ.
- Cuantificar niveles de expresión y normalizar los datos para comparaciones justas.
- Identificar e interpretar genes diferencialmente expresados según el fenotipo (Obeso vs. Normopeso).
- Sintetizar los hallazgos en un póster científico para comunicación académica.

## 🛠️ Metodología de Análisis
- **Mapeo y Cuantificación:** Asignación de lecturas a genes para obtener la matriz de conteos básica.   
- **Normalización:** Ajuste de los datos para reducir variabilidad técnica y diferencias en el número de lecturas.   
- **Contraste:** Identificación de genes con mayor/menor expresión entre los grupos asignados.   
- **Interpretación:** Uso de bases de datos como GeneCards o PubMed para relacionar los genes con el fenotipo de obesidad.

## 📊 Visualizaciones Principales
El análisis genera tres tipos de gráficos clave para la interpretación de resultados:
- **PCA (Principal Component Analysis):** Visualización de la agrupación de las muestras y la contribución de los genes principales a la variancia.
- **Volcano Plot:** Identificación de genes significativos basados en su log2 Fold Change y valor p ajustado.
- **Heatmap:** Mapa de calor de los 10 genes más significativos para observar patrones de expresión entre condiciones.

## 📂 Estructura del Repositorio
```
obesity_RNA_seq_omics/
├── Control de Calidad/                 
├── Matriz de conteos/            
│   ├── Abraham_quant_creado/       
│   └── Bart_quant_creado/
│   └── Homer_quant_creado/       
│   └── Lisa_quant_creado/
│   └── Maggie_quant_creado/          
├── Normalización y Análisis diferencias/           
└── Representación gráfica         
```
## Autores del trabajo
**Autores:** Ainhoa Artetxe, Alba Xiaohe Elias, Alejandra Martin, Alejandro Pascual, Alicia Muñoz, Ana Sofia Santos Tedim

