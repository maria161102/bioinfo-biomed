#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 22 OCTUBRE 23:59
## Se requiere la entrega de un Jupyter Notebook con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) # esto mide las dimensiones
head(data) #vemos las primeras linesas de la tabla
tail(data) # vemos las últimas líneas de la tabla

# Hacemos un primer histograma para explorar los datos
hist(data) 
hist(data,col = "gray", main ="GSE5583 - Histogram") # con esto hemos cambiado el titulo del histograma y también el color

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
# Con la transformación conseguimos manipular los datos y que cambien en base al logaritmo. Sirve para que se puedan visualizar de una mejor manera en el histograma
data_log=log2(data)
hist(data_log)


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot?
#Estos parámetros sirven para cambiar los colores del boxplot. Utilizamos un vector para poder acceder a los elementos que queremos.En este caso los 3 primeros son WT y los 3 siguientes son KO.
#Le ponemos título con el parámetro main. 
#las=2 nos sirve para poner los títulos de los ejes en vertical
#Un boxplot es una gráfica que nos muestra una serie de datos en base a sus cuartiles. En él vemos los cuartiles, la mediana y valores anómalos

boxplot(data_log, col=c("turquoise", "turquoise", "turquoise", "pink", "pink", "pink"),
        main="GSE5583 - boxplot", las=2)

boxplot(data_log)

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
# Sí, es correcta
hc =hclust(as.dist(1-cor(data_log)))
plot(hc, main="GSE5583 - Hierarchical Clustering")


#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
# Hemos generado las columnas separadas de ambos grupos WT y KO. 
# Con head hemos generado una matrix de los datos llamados wt.

wt<-data[,1:3]
ko<-data[,4:6]
class(wt)
head(wt)

# Calcula las medias de las muestras para cada condición. Usa apply
# 1 indica hacer la media de las filas, 2 indica hacer la media de las columnas

wt.mean =apply(wt, 1, mean)
ko.mean=apply(ko, 1, mean)
head(wt.mean)
head(ko.mean)


# ¿Cuál es la media más alta?
# La media más alta es la de KO

max(wt.mean)
max(ko.mean)
 
# Ahora hacemos un scatter plot (gráfico de dispersión)

plot(ko.mean ~ wt.mean)
plot(ko.mean ~ wt.mean, 
     col="hot pink")
# Cambiamos el nombre de los ejes y le ponemos título
plot(ko.mean ~ wt.mean, 
  col="hotpink", 
  xlab ="WT", ylab = "KO", 
  main = "GSE5583 - Scatter Plot")

# Añadir una línea diagonal con abline
abline(0, 1, col="green")

#También se puede añadir una líena horizontal
abline(h=2, col="blue")

#O una vertical
abline(v=5, col= "purple")

# ¿Eres capaz de añadirle un grid?
grid() 

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.mean)
hist(diff.mean, col = "lightyellow")

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# Si hago un t-test con mis datos tranformados no son datos fiables, ya que los estoy manipulando. Los logaritmos sirven mejor en gráficas para poder visualizarlo mejor
# ¿Cuántas valores tiene cada muestra?
# Cada gen tiene dos condiciones, con 6 muestras (3 para cada condición).
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data)) { #Para cada gen 
  x = wt[i,] #gen wt número i
  y = ko[i,] #gen ko número i
  
  #Hacemos el test
  t = t.test(x, y)
  #Añadimos el p-value a la lista
  pvalue[i] = t$p.value
  #Añadimos las estadísticas a la lista
  tstat[i] = t$statistic
}
head(pvalue)
length(pvalue)
# Ahora comprobamos que hemos hecho TODOS los cálculos


# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?
# Que manipulamos los datos y podemos visualizarlos mejor en la gráfica, agrandandolo en escala 
hist(pvalue)
hist(-log10(pvalue), col = "purple")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano", col ="black")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v=diff.mean_cutoff, col = "blue", lwd = 3)
#abline(v=-diff.mean_cutoff, col ="red", lwd=3
abline(h=-log10(pvalue_cutoff), col = "green", lwd=3)

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main =" GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]), col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
plot(diff.mean, -log10(pvalue), main=" GSE5583 - Volcano #3")
points(diff.mean[filter_combined & diff.mean < 0],
       -log10 (pvalue[filter_combined & diff.mean <0]), col = "red")
points(diff.mean[filter_combined & diff.mean >0],
       -log10(pvalue[filter_combined & diff.mean >0]), col = "blue") 

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# Rowv y Colv --> Parámetros que controlan la agrupación y el orden de las filas y las columnas en el heatmap. 
# hclust --> Parámetro para especificar objetos de clustering
# as.dist --> Parámetro que se utiliza para convertir una matriz de datos en una matriz de distancias.
# as.dendogram --> Parámetro para personalizar cómo se agrupan y ordenan las filas o columnas.
# 1-cor --> Parámetro para representar la distancia o la disimilitud entre variables en lugar de su similitud.
# cexCol --> Parámetro para ajustar el tamaño de las etiquetas de las columnas.
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv = rowv, Colv = colv, cexCol = 0.7, labRow = FALSE)

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Hacemos nuestro heatmap


# Lo guardamos en un archivo PDF
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,col=redgreen(75), scale ="row", labRow=FALSE)

# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table(filtered, "GSE5583_DE.txt", sep = "\t", 
            quote = FALSE)
