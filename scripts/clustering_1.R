rm(list=ls())

library(ggplot2)
library(cluster)

iris_data <- iris[, 1:4] #se carga dataset que viene por defecto en R
 
# CLUSTERING NO JERARQUICO
set.seed(1234)
km <- kmeans(iris_data, centers = 3) #se asigno 3 cluster o k =3, ya que conocemos que el dataset tiene 3  grupos

# Resultados
table(km$cluster)  # distribución de clústeres
table(km$cluster, iris$Species)  # comparación con las clases reales

# Representacion grafica
pca <- prcomp(iris_data)
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Cluster = as.factor(km$cluster))

ggplot(df, aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Clustering no jerárquico con k-means (PCA)")


# Otra forma de hacerlo
library(cluster)
pam_result <- pam(iris_data, k = 3)
# Resultados
table(pam_result$clustering, iris$Species)
# Visualización
clusplot(pam_result, main = "Clustering con PAM", color = TRUE, shade = TRUE)


## CLUSTERING JERARQUICO

# Calcular la matriz de distancias
dist_matrix <- dist(iris_data) #se puede coocar como parametro que distancia o formula de D se va a usar
dist_matrix

# Realizar el clustering con diferentes metodos de agrupamiento
hclust_single <- hclust(dist_matrix, method = "single")
hclust_complete <- hclust(dist_matrix, method = "complete")
hclust_average <- hclust(dist_matrix, method = "average")
hclust_ward <- hclust(dist_matrix, method = "ward.D2")

# Variable de indices para el orden del dendograma
hclust_single$order
hclust_complete$order
hclust_average$order
hclust_ward$order

# Graficar resultados
plot(hclust_single, main = "Single Linkage")
plot(hclust_complete, main = "Complete Linkage")
plot(hclust_average, main = "Average Linkage")
plot(hclust_ward, main = "Ward's Linkage")

