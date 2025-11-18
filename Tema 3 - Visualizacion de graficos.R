# Material adicional para alumnos MU Bioinformática
# Tema 3. Visualización de gráficos


# Resetamos entorno de R ----
rm (list = ls())



# Librerías
library(ggplot2)
library(dplyr)
library(mice)
library(ggstatsplot)


# Importar CSV ----

## Establecemos directorio donde cogeremos el archivo ----
setwd("/Users/samuelsalazardiaz/Desktop/1er Cuatri/Estadística y R para Ciencias de la Salud/Scripts/Tema 3")

# Detectar el archivo por patrón (evita errores por espacios/acentos)
f <- list.files(pattern = "^Base de datos repaso tema contrastes de hipotesis.*\\.csv$", ignore.case = TRUE)[1]
df <- if (!is.na(f)) read.csv(f, fileEncoding = "UTF-8-BOM", check.names = FALSE) else stop("No se encontró el archivo .csv")

## depuracion de datos ----
str(df)
df$sexo <- as.factor(df$sexo)
df$trat <- as.factor(df$trat)
colSums(is.na(df))


## imputar metodo mice ----
imputed_data <- mice(df, m = 5, method = "pmm", seed = 123)
df <- complete(imputed_data)
colSums(is.na(df))
any(is.na(df))

## imputamos por mediana ----
df <- df %>% 
  mutate(across(where(is.numeric), 
                ~ifelse(is.na(.), 
                        median(., na.rm = TRUE), 
                        .)))

colSums(is.na(df))


## imputamos por media ----
df <- df %>% 
  mutate(across(where(is.numeric), 
                ~ifelse(is.na(.), 
                        mean(., na.rm = TRUE), 
                        .)))

colSums(is.na(df))



# Ejercicio 1: Histogramas ----
## Cuál es el rango de edad más frecuente en la muestra? ----
## ¿La distribución es simétrica o sesgada? ----


# Ejercicio 2: Diagramas de Barras ----
## ¿Hay más hombres o mujeres en la muestra? ----
## ¿Qué podrías inferir sobre la composición de la población estudiada? ---- 


# Ejercicio 3: Boxplots ----
## ¿Hay diferencias notables entre el IMC de hombres y mujeres? ----
## ¿Existen valores atípicos? ----


# Ejercicio 4: Diagramas de Dispersión ----
## ¿Existe una relación positiva entre el peso y el IMC? ----
## ¿Hay valores atípicos en la relación? ----
## ¿Existen diferencias por el sexo? ----


# Ejercicio 5: Heatmap de Correlación ----
## ¿Qué variables tienen una correlación fuerte entre estas: delta_peso_trat, delta_musculo_trat, delta_imc_trat, delta_grasa_trat, delta_agua_trat, delta_grasa_visceral_trat? ----
## ¿Existen correlaciones inesperadas? ----
## ¿Existen correlaciones diferentes entre tratamiento en mujeres? ----


















# Ejercicio 1: Histogramas ----
## Cuál es el rango de edad más frecuente en la muestra? ----
## ¿La distribución es simétrica o asimétrica? ----
ggplot(df, aes(x = edad)) +
  geom_histogram(aes(y = ..density..),fill = "blue", color = "black", alpha = 0.7, bins = 50) +
  geom_density(color = "red", size = 2) +
  scale_x_continuous(breaks = seq(20, 80, by = 2)) +
  labs(title = "Distribución de la edad", x = "Edad", y = "Frecuencia") +
  theme_minimal()

ggplot(df, aes(x = edad)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7, bins = 50) +
  labs(title = "Distribución de la edad", x = "Edad", y = "Frecuencia") +
  theme_minimal()


# Ejercicio 2: Diagramas de Barras ----
## ¿Hay más hombres o mujeres en la muestra? ----
## ¿Qué podrías inferir sobre la composición de la población estudiada? ---- 
ggplot(df, aes(x = sexo, fill = sexo)) +
  geom_bar(color = "black", alpha = 0.7) +
  labs(title = "Distribución por sexo", x = "Sexo", y = "Frecuencia") +
  theme_minimal()


# Ejercicio 3: Boxplots ----
## ¿Hay diferencias notables entre el IMC de hombres y mujeres? ----
## ¿Existen valores atípicos? ----
ggplot(df, aes(x = sexo, y = imc0, fill = sexo)) +
  geom_boxplot() +
  labs(title = "Distribución del IMC por sexo", x = "Sexo", y = "IMC") +
  theme_minimal()


# Ejercicio 4: Diagramas de Dispersión ----
## ¿Existe una relación positiva entre el peso y el IMC? ----
## ¿Hay valores atípicos en la relación? ----
ggplot(df, aes(x = peso0, y = imc0)) +
  geom_point(color = "darkgreen", alpha = 0.5) +
  labs(title = "Relación entre Peso e IMC", x = "Peso", y = "IMC") +
  theme_minimal()

ggplot(df, aes(x = peso0, y = imc0)) +
  geom_point(color = "darkgreen", alpha = 0.5) +
  geom_smooth(method = "lm", color= "red", se = TRUE) +
  labs(title = "Relación entre Peso e IMC", x = "Peso", y = "IMC") +
  theme_minimal()

## ¿Existen diferencias por el sexo? ----
levels(df$sexo)
ggplot(df, aes(x = peso0, y = imc0, color = sexo)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("H" = "lightgreen", "M" = "lightpink"))+
  geom_smooth(aes(x = peso0, y = imc0, fill = sexo), method = "lm", color= "red", se = TRUE) +
  labs(title = "Relación entre Peso e IMC", x = "Peso", y = "IMC") +
  theme_minimal()



# Ejercicio 5: Heatmap de Correlación ----
## ¿Qué variables tienen una correlación fuerte entre estas: delta_peso_trat, delta_musculo_trat, delta_imc_trat, delta_grasa_trat, delta_agua_trat, delta_grasa_visceral_trat? ----
## ¿Existen correlaciones inesperadas? ----

# Seleccionar solo las variables numéricas
colnames(df)
numeric_df <- df %>% select_if(is.numeric)
select_df <- df %>% select(delta_peso_trat, delta_musculo_trat, delta_imc_trat, delta_bmr_trat, delta_grasa_trat, delta_agua_trat, delta_grasa_visceral_trat)

# Calcular la matriz de correlación
cor_matrix <- cor(select_df)

# Crear el heatmap
pheatmap(cor_matrix)
pheatmap(scale(select_df))
table(df$imc0_cat)


## ¿Existen correlaciones diferentes entre tratamiento en mujeres? ----
select_trat1_df <- df %>% 
  filter(trat == "control" & sexo == "M" & imc0_cat =="obeso I basal") %>%
  select(delta_peso_trat, delta_musculo_trat, delta_imc_trat, delta_grasa_trat, delta_agua_trat, delta_grasa_visceral_trat)

pheatmap(select_trat1_df, 
         main = "Heatmap - Grupo Control")

select_trat2_df <- df %>% 
  filter(trat == "farmaco" & sexo == "M" & imc0_cat =="obeso I basal") %>%
  select(delta_peso_trat, delta_musculo_trat, delta_imc_trat, delta_grasa_trat, delta_agua_trat, delta_grasa_visceral_trat)

# Heatmap para el grupo fármaco
pheatmap(select_trat2_df, 
         main = "Heatmap - Grupo Fármaco")

