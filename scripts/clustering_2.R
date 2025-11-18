rm(list=ls())

library(dplyr)
library(factoextra)
library(stats)
library(ggdendro)
library(cluster)
library(ggplot2) 
library(gridExtra)
library(pheatmap)

setwd("/home/samuel/Documentos/Documents/clases/algoritmos e IA/tema_4/tema__4_export/dataset")

# Lectura de datos
df <- read.csv("Dataset expresión genes.csv")

# Extraemos los genes que empiecen por AQ
df_genes <- df %>% dplyr::select(starts_with("AQ_"))
str(df_genes)

# Mapa de calor
pheatmap(df_genes, scale="row", clustering_method = "complete")

# Hacemos analisis de como están los datos
is.na(colSums(df_genes)) # ver si hay missing
df_genes_scale <- scale(df_genes)  # Normalización z-score



#### Clustering no jerarquico: K-means

# funcion kmeans
#   x: matriz numerica o dataframe
#   centers: numero de clusteres que se desean formar
#   iter.max: numero de iteraciones para el algoritmo
#   nstart: veces que se ejecuta el algoritmo con diferentes centroides (10 - 25)
#   algorithm: varias opciones (Hartingan-wong por defecto)


# k=5
kmeans.result <- kmeans(df_genes_scale, centers = 5, iter.max = 100, nstart = 25)
# Visualizacion
fviz_cluster(kmeans.result, df_genes_scale, xlab = '', ylab = '') +
  ggtitle("Cluster plot, centers = 5", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))

# k=3
kmeans.result <- kmeans(df_genes_scale, centers = 3, iter.max = 100, nstart = 25)
# Visualizacion
fviz_cluster(kmeans.result, df_genes_scale, xlab = '', ylab = '') +
  ggtitle("Cluster plot, centers = 3", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))

# n optimo de clusters
fviz_nbclust(df_genes_scale, kmeans, method = "wss") +
  ggtitle("Optimal number of clusters", subtitle = "") +
  theme_classic()





#### Clustering jerarquico aglomerativo

# Calcular la matriz de distancia
dist_matrix <- dist(df_genes_scale)

# Se ejecuta el algoritmo de clusterización jerárquica aglomerativa

# agrupa los clusters usando la distancia entre los puntos más CERCANOS
hclust_model_single <- hclust(dist_matrix, method = "single") 
# agrupa los clusters usando la distancia entre los puntos más ALEJADOS
hclust_model_complete <- hclust(dist_matrix, method = "complete") 
# agrupa los clusters usando el PROMEDIO de todas las distancias entre los puntos de ambos clusters
hclust_model_average <- hclust(dist_matrix, method = "average") 
# agrupa los clusters tratando de que sean lo más COMPACTOS posible minimizando la dispersión interna
hclust_model_ward <- hclust(dist_matrix, method = "ward.D") 



# Definimos los colores que vamos a utilizar en los dendogramas
colors <- rainbow(5)

# single: conecta puntos cercanos y puede encontrar clusters largos y delgados, 
#         PERO a veces conecta muchos puntos formando cadenas, lo que puede ser
#         poco útil
clust_single <- fviz_dend(hclust_model_single, 
                          cex = 0.5,
                          k = 5,
                          palette = colors,
                          main = "Single",
                          xlab = "Índice de Observaciones",
                          ylab = "Distancia") + theme_classic()
clust_single

# complete: crea clusters compactos y bien definidos, PERO es sensible a puntos 
#           extremos (outliers), que pueden distorsionar los clusters
clust_complete <- fviz_dend(hclust_model_complete, 
                          cex = 0.5,
                          k = 5,
                          palette = colors,
                          main = "Complete",
                          xlab = "Índice de Observaciones",
                          ylab = "Distancia") + 
  theme_classic()

clust_complete
# average:  encuentra un equilibrio entre single y complete, PERO a veces no es 
#           tan bueno para clusters con tamaños muy diferentes
clust_average <- fviz_dend(hclust_model_average, 
                          cex = 0.5,
                          k = 5,
                          palette = colors,
                          main = "Average",
                          xlab = "Índice de Observaciones",
                          ylab = "Distancia") + 
  theme_classic()

# Ward: crea clusters redondeados y homogéneos similares a los que genera 
#       k-means, PERO no funciona tan bien si los clusters tienen formas raras 
#       o tamaños muy diferentes
clust_ward <- fviz_dend(hclust_model_ward, 
                          cex = 0.5,
                          k = 5,
                          palette = colors,
                          main = "Ward",
                          xlab = "Índice de Observaciones",
                          ylab = "Distancia") + 
  theme_classic()

grid.arrange(clust_single, clust_complete, clust_average, clust_ward, nrow = 2)


# single: datos con patrones lineales -> "Creo que los datos tienen estructuras locales fuertes, y estoy más interesado en cómo se conectan los puntos cercanos."
# complete: datos donde esperas clusters compactos y bien definidos -> "Los datos están agrupados en regiones claramente separadas y no quiero que un punto extremo distorsione el análisis"
# average: cuando no tienes una forma clara en mente para los clusters -> "Espero una mezcla de clusters compactos y algo más dispersos, pero quiero un balance entre lo local y lo global."
# ward.D: datos en los que esperas clusters compactos y homogéneos -> "Los datos deberían agruparse en clusters compactos con baja variabilidad interna."


# Podemos asignar al cluster al que pertenece cada elemento
df$cluster_single <- as.factor(cutree(hclust_model_single, k = 5))
df$cluster_complete <- as.factor(cutree(hclust_model_complete, k = 5))
df$cluster_average <- as.factor(cutree(hclust_model_average, k = 5))
df$cluster_ward <- as.factor(cutree(hclust_model_ward, k = 5))





#### Clustering jerarquico divisivo

# Implementación del clustering divisivo

# ideal para datos donde las distancias más pequeñas
diana_euclidean <- diana(df_genes_scale, metric = "euclidean", stand = FALSE) 
# ideal para datos con diferentes escalas o datos categóricos
diana_manhattan <- diana(df_genes_scale, metric = "manhattan", stand = FALSE) 

colors <- rainbow(5)
clust_diana_euclidean <- fviz_dend(diana_euclidean, 
                                   cex = 0.5, 
                                   k = 5,
                                   palette = colors, 
                                   main = 'Euclidean',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()


colors <- rainbow(5)
clust_diana_manhattan <- fviz_dend(diana_manhattan, 
                                   cex = 0.5, 
                                   k = 5,
                                   palette = colors, 
                                   main = 'Manhattan',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()


grid.arrange(clust_diana_euclidean, clust_diana_manhattan, nrow = 2)




# Ejemplo de las diapositivas 
########################################

rm(list=ls())
# Cargar librerías
library(ggplot2)
library(cluster)

# Ejemplo de coordenadas (representación de plantas, clientes, tiendas, etc.)
data <- data.frame(
  x = c(2, 4, 6, 10),  # Coordenadas x (pueden ser precios, calificaciones, etc.)
  y = c(3, 4, 7, 10)   # Coordenadas y
)

# Calcular la matriz de distancias (usando Euclidiana)
dist_matrix <- dist(data)
dist_matrix

# Realizar el clustering con diferentes métodos de enlace
hclust_single <- hclust(dist_matrix, method = "single")
hclust_complete <- hclust(dist_matrix, method = "complete")
hclust_average <- hclust(dist_matrix, method = "average")
hclust_ward <- hclust(dist_matrix, method = "ward.D2")

hclust_single$order
hclust_complete$order
hclust_average$order
hclust_ward$order

# Graficar resultados
par(mfrow = c(2, 2))  # Dividir la ventana en 2x2 para graficar múltiples
plot(hclust_single, main = "Single Linkage")
plot(hclust_complete, main = "Complete Linkage")
plot(hclust_average, main = "Average Linkage")
plot(hclust_ward, main = "Ward's Linkage")







