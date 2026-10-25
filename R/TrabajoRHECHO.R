##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

#Primero voy a importar el archivo con los nombres de las columnas y lo guardamos en una variable que la voy a llamar datos
datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "\t")
 
#Ahora, vamos a utilizar las funciones head(), tail() y dim()
#head() nos va a mostrar las primeras seis filas
head(datos)

#tail() nos va a mostrar las últimas seis filas 
tail(datos)

#dim() nos muestra el número de filas y de columnas 
dim(datos)


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

#Primero creo las variables 
control =  datos[c("C2.TPM","C1.TPM","C3.TPM")]
condicion = datos[c( "X1W5.TPM","X1W6.TPM","X1W5.COV")]

#Ahora vamos a calcular la media de cada gen en cada variable utilizando apply()
media_control = apply(control, 1, mean)
media_condicion = apply(condicion, 1, mean)

#Ahora vamos a calcular la desviación estándar de cada gen en cada variable
de_control = apply(control, 1, sd)
de_condicion = apply(condicion, 1, sd)

#Añado los resultados como nuevas columnas
datos$media_control = media_control
datos$de_control = de_control
datos$media_condicion = media_condicion
datos$de_condicion = de_condicion

#Vamos a comprobar con head() y dim() que las columnas se hayan añadido
head(datos)
dim(datos)


### 3.- En una sola imagen, muestra un histograma para la media del control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

#Para poner ambos histogramas en una sola imagen, voy a utilizar la función par() de R,
#que nos permite cambiar los ajustes de la ventana donde se dibujan los gráficos. Además, dentro voy a utilizar
#mfrow para indicar en cuántos huecos quiero dividir la ventana. En este caso en una fila y dos columnas, por ello le añado c(1,2)

par(mfrow = c(1,2))

#Creo cada histograma
hist(datos$media_control, main = "Media de control", col = "lightblue", xlab = "Media")
hist(datos$media_condicion, main = "Media de condición", col = "lightgreen", xlab = "Media")

#Al ver los histogramas, se observa que casi todos los genes caen en una única barra pegada al cero.
#Esto pasa porque los datos tienen un sesgo positivo, es decir, casi todos los genes tienen valores bajos
#pero unos pocos tienen valores enormes y estiran todo el eje.
#Para ver mejor la forma de los datos aplico una transformación logarítmica en base 2, que comprime los valores
#altos y permite apreciar mejor la distribución de los valores bajos. Como el logaritmo de 0 no existe (da -Inf) y hay genes con media 0, sumo una constante
#a cada valor, en este caso 1, antes de calcular el logaritmo. Así todos los valores son mayores que cero
#y no pierdo ningún gen. Además, los genes con media 0 se quedan en 0 después de la transformación y
#no salen valores negativos. Guardo los resultados en variables nuevas para no modificar la tabla original.
#Esta transformación la uso únicamente para las gráficas. Los tests estadísticos los haré con los datos sin transformar.

control_log = log2(datos$media_control + 1)
condicion_log = log2(datos$media_condicion + 1)

#Vuelvo a poner par(mfrow = c(1,2)) para que se vean los dos histogramas en una imagen 
par(mfrow = c(1,2))

hist(control_log, main = "Media de control (log2)", col = "lightblue", xlab = "log2(Media + 1)")
hist(condicion_log, main = "Media de condición (log2)", col = "lightgreen", xlab = "log2(Media + 1)")

#Ahora se pueden ver los histogramas más bonitos y se aprecia mejor la forma de los datos

#Para que en el siguiente apartado no salga en solo medio hueco de la ventana, vuelvo a poner los ajustes iniciales, es decir, una fila y una columna
par(mfrow = c(1,1))


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

#Vamos a hacer cada boxplot con la leyenda arriba a la izquierda
boxplot(datos$media_control, datos$media_condicion, col = c("lightblue", "lightgreen"), names = c("Control", "Condición"), main = "Medias de control y condición", ylab = "Media")

legend("topleft", legend = c("Control", "Condición"), fill = c("lightblue", "lightgreen"))

#Al igual que pasaba con los histogramas, las cajas quedan aplastadas contra el cero y no se ve nada,
#por el sesgo positivo de los datos. Por eso hago ahora el boxplot con los datos transformados con log2
#que ya calculé en el anterior ejercicio (control_log y condicion_log)

#He añadido ylim para ajustar el eje vertical y dejar espacio para la leyenda, porque sin él la leyenda tapaba algunos puntos del control
boxplot(control_log, condicion_log, col = c("lightblue", "lightgreen"), names = c("Control", "Condición"), main = "Medias de control y condición (log2)", ylab = "log2(Media + 1)", ylim = c(0, 24))

legend("topleft", legend = c("Control", "Condición"), fill = c("lightblue", "lightgreen"))


 
### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### Realiza un test de comparativa de medias para cada gen entre condición y control. 
### Guarda el p-value de la comparativa y añádelo como una columna extra a tu tabla original.

#Primero voy a comprobar si las medias de expresión de los genes siguen una distribución normal
#utilizando el test de Shapiro-Wilk. La hipótesis nula es que los datos siguen una distribución normal
#Si el p-value es menor que 0.05, rechazo la hipótesis nula de normalidad.

shapiro.test(datos$media_control)
shapiro.test(datos$media_condicion)

#En ambos casos obtengo un p-value menor que 0.05, por lo que rechazo la hipótesis nula
#y concluyo que las medias de expresión de los genes no siguen una distribución normal.
#Sin embargo, este resultado se refiere a las medias de los 619 genes en conjunto y no
#permite saber si las réplicas de cada gen siguen una distribución normal.
#Además, como solo tengo tres réplicas por grupo, resulta difícil comprobar la normalidad
#de cada gen por separado.

#Para comparar las medias de control y condición de cada gen voy a utilizar el t- test,
#que permite comparar ambos grupos aunque tengan varianzas diferentes.
#Como solo tengo tres réplicas por grupo, tendré que interpretar los resultados con precaución.

#Creo un vector para guardar el p-value que obtenga para cada gen
pvalue = rep(NA_real_, nrow(datos))

#Ahora recorro todos los genes para comparar el control y la condición de cada uno de ellos
for(i in 1:nrow(datos)) {

  #Selecciono las tres réplicas del control del gen correspondiente
  x = unlist(control[i,])

  #Selecciono las tres réplicas de la condición del mismo gen
  y = unlist(condicion[i,])

  #Si los dos grupos tienen una desviación estándar igual a cero,
  #dejo el p-value como NA porque no se puede realizar correctamente el test
  if(sd(x) == 0 & sd(y) == 0) {
    pvalue[i] = NA_real_
  } else {

    #Realizo el t-test para comparar ambos grupos.
    #Utilizo tryCatch para que el bucle pueda continuar si algún test da error.
    test = tryCatch(t.test(x, y), error = function(e) NULL)

    #Si el test se ha realizado correctamente, guardo su p-value.
    if(!is.null(test)) {
      pvalue[i] = test$p.value
    }
  }
}

#Compruebo que he obtenido un resultado para cada uno de los 619 genes
length(pvalue)

#Compruebo si hay algún gen para el que no se haya podido calcular el p-value
sum(is.na(pvalue))

#Añado los p-values a la tabla original como una nueva columna
datos$pvalue = pvalue

#Compruebo que la columna se ha añadido correctamente
head(datos)
dim(datos)
#Compruebo que se ha generado un resultado para cada gen
length(pvalue)

#Compruebo cuántos genes no tienen un p-value calculado
sum(is.na(pvalue))

#Añado los p-values como una nueva columna a la tabla original
datos$pvalue = pvalue

#Compruebo que la nueva columna se ha añadido correctamente
head(datos)
dim(datos)


### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

#Primero establezco el límite de significación en 0.05
limite = 0.05

#Creo un filtro para seleccionar los genes cuyo p-value sea menor que 0.05.
#También compruebo que tengan un p-value calculado.
filtro = !is.na(datos$pvalue) & datos$pvalue < limite

#Utilizo el filtro para crear una tabla con los genes que cumplen este criterio
genes_filtrados = datos[filtro, ]

#Muestro el listado de genes diferencialmente expresados
#Incluyo el nombre del gen, las medias de ambos grupos y el p-value.
genes_filtrados[, c("Gene.Name", "media_control", "media_condicion", "pvalue")]

#Calculo cuántos genes cumplen el criterio de significación
n_genes = nrow(genes_filtrados)
n_genes

#Para saber cuáles están sobreexpresados y cuáles reprimidos,
#comparo la media de expresión de la condición con la del control.
#Si la media de condición es mayor, el gen está sobreexpresado.
#Si la media de condición es menor, el gen está reprimido.

#Cuento los genes sobreexpresados
sobreexpresados = sum(genes_filtrados$media_condicion > genes_filtrados$media_control)
sobreexpresados

#Cuento los genes reprimidos
reprimidos = sum(genes_filtrados$media_condicion < genes_filtrados$media_control)
reprimidos

#En total, hay 362 genes con un p-value menor que 0.05.
#De estos, 50 tienen una mayor expresión en la condición y 312 tienen
#una menor expresión que en el control.


### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

#Selecciono las seis muestras de los genes diferencialmente expresados
datos_heatmap = genes_filtrados[, c("C2.TPM", "C1.TPM", "C3.TPM", "X1W5.TPM", "X1W6.TPM", "X1W5.COV")]

#Convierto los datos seleccionados en una matriz numérica para poder realizar el heatmap
datos_heatmap = as.matrix(datos_heatmap)

#Transformo los datos utilizando log2(x + 1)
#Sumo 1 para poder transformar los valores iguales a cero.
datos_heatmap = log2(datos_heatmap + 1)

#Realizo un agrupamiento jerárquico de los genes y de las muestras
rowv = as.dendrogram(hclust(as.dist(1 - cor(t(datos_heatmap)))))
colv = as.dendrogram(hclust(as.dist(1 - cor(datos_heatmap))))

#Represento el heatmap con los genes y las muestras agrupados
heatmap(datos_heatmap, Rowv = rowv, Colv = colv,
	cexCol = 0.7, labRow = FALSE,
	col = colorRampPalette(c("blue", "green"))(50),
	scale = "row",
	main = "Heatmap de genes diferencialmente expresados")