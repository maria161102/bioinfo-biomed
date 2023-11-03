
# TRABAJO R

# Mª Isabel Gálvez Castaño. 132557

##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")


# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/datos-trabajoR", followlocation = TRUE)
data = read.table(file = "datos-trabajoR.txt", head = T)


##############################################################################
# 1. Carga los datos y exáminalos en R. Emplea las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos?
# Hay 2 variables y 5 tratamientos
head(data) #vemos las primeras líneas de la tabla
summary(data) #resumen de los datos
dim(data) # esto mide las dimensiones
str(data) #permite visualizar la estructura interna de los datos


##############################################################################
# 2. Haz un boxplot para nuestros datos. Uno para cada variable. Colorea a Variable 1 y a Variable 2 de forma diferente 
boxplot(Variable1~Tratamiento, data = data, main="Boxplot Variable 1", col="hotpink", ylab = "Variable 1")
boxplot(Variable2~Tratamiento, data = data, main="Boxplot Variable 2", col="turquoise", ylab = "Variable 2")


##############################################################################
# 3. Haz un gráfico de dispersión con las dos variables. Cada tratamiento debe de ir de un color distinto.
plot(Variable2~Variable1, data = data, col=data$Tratamiento, main="Gráfico de dispersión - Variables 1 y 2", xlab="Variable 1", ylab = "Variable 2")


##############################################################################
# 4. Ponle leyenda al gráfico del apartado anterior. En el margen inferior derecho.
legend(x="bottomright", legend = c("Tto 1","Tto 2", "Tto 3", "Tto 4", "Tto 5"), fill = c("black", "red", "green", "lightblue", "blue"), title = "Tratamiento")


##############################################################################
# 5. Haz un histograma para cada variable. Recuerda mantener los colores.
variable1=data[,2]
variable2=data[,3]
tratamiento=data[,1]
hist(x=variable1, col="hotpink", main = "Histograma de la Variable 1", xlab = "Variable 1")
hist(x=variable2, col = "turquoise", main = "Histograma de la Variable 2", xlab = "Variable 2")


##############################################################################
# 6. Haz un factor en la columna tratamiento y guárdalo en una variable.
factor_tratamiento=factor(tratamiento)
levels(factor_tratamiento)


##############################################################################
# 7.Calcula la media y la desviación estándar para cada tratamiento
mean_variable1 <- tapply(data$Variable1, data$Tratamiento, mean)
sd_variable1 <- tapply(data$Variable1, data$Tratamiento, sd)
mean_variable1
sd_variable1

mean_variable2 <- tapply(data$Variable2, data$Tratamiento, mean)
sd_variable2 <- tapply(data$Variable2, data$Tratamiento, sd)
mean_variable2
sd_variable2


##############################################################################
# 8. Averigua cuántos elementos tiene cada tratamiento.
# Cada tratamiento tiene 10 elementos.
table(factor_tratamiento)


##############################################################################
# 9. Extrae los datos para el tratamiento 1 y el tratamiento 4 y guárdalos cada uno en una variable diferente.
datos_tratamiento1 = subset(data, tratamiento==1)
datos_tratamiento4 = subset(data, tratamiento==4)


##############################################################################
# 10. Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la Variable 1 son iguales.
# ¿Puedes comprobarlo? Para ello, necesitarás comprobar primero si los datos se distribuyen de forma normal.
# En función del resultado de la prueba de normalidad, ¿qué test usarías? ** En general, asumimos que las muestras son independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo a tus resultados.

shapiro.test(datos_tratamiento1$Variable1)
shapiro.test(datos_tratamiento4$Variable1)
# El resultado del shapiro test es que en ambos casos los datos siguen una distribucción normal.

# En función del resultado utilizo el método t-Student.
t.test(datos_tratamiento1$Variable1, datos_tratamiento4$Variable1)
# Según el resultado del t.test las medias de los tratamiento 1 y 4 no son iguales.


# Comprobamos si sus varianzas son iguales con Fisher’s F-test
var.test(datos_tratamiento1$Variable1, datos_tratamiento4$Variable1)
# Obtenemos un valor de p-value<0,05, se rechaza la hipóteis nula, lo que indica que sus varianzas no son iguales.

