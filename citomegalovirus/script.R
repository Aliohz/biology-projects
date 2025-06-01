
## Simulación

## simulación usando Poisson y Gamma
## Codigo para 10000 replicas de la simulación

## largo de adn
lengthdna <- 229354
## número de palindromes
numberpalin <- 296
## iteraciones
iterate <- 10000

#matris para guardar las simulaciones

sim <- matrix(nrow=iterate, ncol=numberpalin)
for(i in 1:iterate){
  sim[i,] <- sample(1:lengthdna, numberpalin, replace=F)
}

# Codigo para sumar el número de conteos
breakslength <- 230001
breaksize <- 500
mybreaks <- seq(1,breakslength,breaksize)

# Matriz generada para guardar la simulación
sim.hist <-matrix(nrow=iterate, ncol=length(mybreaks)-1)
# Bucle para obtener los conteos al azar generados 

for(i in 1:iterate){
  sim.hist[i,] <- hist(sim[i,], breaks=mybreaks, right=F, plot=F)$counts
}
# Compilar los resultados de conteos en densidades de las observasiones
# Revisar la ayuda para el comando (table)  and probar  table(sample(2:11, 1000, replace=T))
sim.dens<-table(sim.hist)/(nrow(sim.hist)*ncol(sim.hist))

# codigo para resumir las distancias entre los datos de las octetas 
distances<-apply(sim, 1, function(x){diff(sort(x), lag=7)})

# Análisis de los datos simulados y analisis de las probabilidades
## Conteos
## Grafica
plot(as.numeric(names(sim.dens)), sim.dens, ylab="Frequencia", xlab="Número de palindromes")
mtext(side=3, line=2.7, "Simulados (circulos)", cex=1.3)
mtext(side=3, line=1.2, "Distribución teorica (linea)", cex=1.3)
lines(as.numeric(names(sim.dens)), dpois(0:(length(sim.dens)-1),
                                         lambda=(numberpalin*breaksize)/lengthdna ), col="red")

# Cargar los datos de palindromos en CMV, sustituyendo uno de los datos de la simulación
real_data <- read.csv("cmv.csv")
cmv <- real_data$location

#### Analizar la distribución de palindromes por secciones de adn
## Aquí s importante analizar la acumulación de palindromes (lag) en el espacio del adn
## Entonces para esto usaremos la distribución gamma.
## la distribución gamma requiere dos momementos, llamados shape (alfa) y rate (beta).
## También analizaremos la acumulación de palindromes en secciones de diferentes largos

par(mfcol = c(3, 2))

## Cada grafica tiene en el eje horizontal la localización en la secuencia de ADN en acumulaciones de palindormes de 5 a 10 palindromes
## el eje vertical tiene la distribución de probabilidad de observar el espaciamiento en estos grupos de clusters bajo una distribución gamma.
# La probabilidades estan en escala de logaritmos naturales.  

plot(cmv[1:(296-5)], pgamma(diff(cmv, lag=5), shape=5, rate=296/229354, log=T), ylab = "", xlab = "", main = "Grupos de 5", pch = 19) ## lag corresponde al número de palindromes donde se tendra la distancia del primer palindrome al ultimo palindrome del grupo. 
abline (h = c ( log(0.01), log(0.05), log(0.1), log(0.5)), lty = 2) ##  Agregar las lineas para señalar las p = 0.01, p = 0.05, p = 0.1 and p = 0.5 lty = 2 hace lineas punteads. 
text(x=c(230000,230000,230000),y=c(log(0.01), log(0.05), log(0.1),  log(0.5)),labels=c("0.01","0.05","0.1", "0.5"),cex=.6) ## Agregar etiquetas a las lineas

plot(cmv[1:(296-6)], pgamma(diff(cmv, lag=6), shape=6, rate=296/229354, log=T), ylab = "Dist. Gamma (log p)", xlab = "", main = "Grupos de 6", pch = 19)
abline (h = c ( log(0.01), log(0.05), log(0.1), log(0.5)), lty = 2)    
text(x=c(230000,230000,230000),y=c(log(0.01), log(0.05), log(0.1), log(0.5)),labels=c("0.01","0.05","0.1", "0.5"),cex=.6) 

plot(cmv[1:(296-7)], pgamma(diff(cmv, lag=7), shape=7, rate=296/229354, log=T), ylab = "", xlab = "Localización en la secuencia de ADN (bp)", main = "Grupos de 7", pch = 19)
abline (h = c ( log(0.01), log(0.05), log(0.1), log(0.5)), lty = 2)   ## Adding lines which is pointing out where p = 0.01, p = 0.05, p = 0.1 and p = 0.5 lty = 2 makes a dashed line. 
text(x=c(230000,230000,230000),y=c(log(0.01), log(0.05), log(0.1), log(0.5)),labels=c("0.01","0.05","0.1", "0.5"),cex=.6) ## adding labels for cutting lines

plot(cmv[1:(296-8)], pgamma(diff(cmv, lag=8), shape=8, rate=296/229354, log=T), ylab = "", xlab = "", main = "Grupos de 8", pch = 19)
abline (h = c ( log(0.01), log(0.05), log(0.1), log(0.5)), lty = 2)   ## Adding lines which is pointing out where p = 0.01, p = 0.05, p = 0.1 and p = 0.5 lty = 2 makes a dashed line. 
text(x=c(230000,230000,230000),y=c(log(0.01), log(0.05), log(0.1), log(0.5)),labels=c("0.01","0.05","0.1", "0.5"),cex=.8) ## adding labels for cutting lines

plot(cmv[1:(296-9)], pgamma(diff(cmv, lag=9), shape=9, rate=296/229354, log=T), ylab = "", xlab = "", main = "Grupos de 9", pch = 19)
abline (h = c ( log(0.01), log(0.05), log(0.1), log(0.5)), lty = 2)   ## Adding lines which is pointing out where p = 0.01, p = 0.05, p = 0.1 and p = 0.5 lty = 2 makes a dashed line. 
text(x=c(230000,230000,230000),y=c(log(0.01), log(0.05), log(0.1), log(0.5)),labels=c("0.01","0.05","0.1", "0.5"),cex=.6) ## adding labels for cutting lines

plot(cmv[1:(296-10)], pgamma(diff(cmv, lag=10), shape=10, rate=296/229354, log=T), ylab = "", xlab = "Localización en la secuencia de ADN (bp)", main = "Groups of 10", pch = 19)
abline (h = c ( log(0.01), log(0.05), log(0.1), log(0.5)), lty = 2)   ## Adding lines which is pointing out where p = 0.01, p = 0.05, p = 0.1 and p = 0.5 lty = 2 makes a dashed line. 
text(x=c(230000,230000,230000),y=c(log(0.01), log(0.05), log(0.1), log(0.5)),labels=c("0.01","0.05","0.1", "0.5"),cex=.6) ## adding labels for cutting lines

##### Unimos esto con una distribución de frecuencias
par(mfrow=c(1,2))
hist750<-hist(cmv,breaks=seq(0, 230000, 750), ylim = c(0,8)
              ,col="grey", main= "a)", xlab="Localización en la secuencia de ADN (bp)", xlim=c(0,250000),ylab="Frecuencia de palindromes", cex = 0.8)
abline(h=qpois(c(0.99),lambda=((296/229354)*750)), lty=2) ## Adding dashed line to point out 99 % quantile. lty = 2 makes a dashed pattern in the line. 
abline(h=qpois(c(0.95),lambda=((296/229354)*750)), lty=2)  ## Adding dashed line to point out 95 % quantile.

expect.750<- (296/229354)*750 ## Calcular lambda esperado para intervalos de 500 pb
plot(dpois(0:8, lambda = 750*(296/229354)), 0:8, type ='l', main = "b)", ylab= "", xlab= "Distribución de Poisson ", cex = 0.8)
abline(h= expect.750, lty=2) ## Agregar linea discontinua para señalar lambda bajo una distribución de Poissson
text(0.20, 1.15,"Lambda, Promedio, y Esperanza", cex= 0.6) ## etiquetar linea.
abline(h=qpois(c(0.99), (296/229354)*750), lty=2) ## Agregar linea para señalar el cuantil 99%. lty = 2 hace que la linea sea discontinua
text(0.20,4.25,"Percentil 99%", cex = 0.6)
abline(h=qpois(c(0.95), (296/229354)*750), lty=2) ## Agregar una linea discontinua para señalar el cuantil 95% . lty = 2 hace que la linea sea discontinua
text(0.20,3.25,"Percentil 95%", cex = 0.6) 

par(mfrow=c(1,1))

## Función para obtener la localizacion de inicio y final de los agrupameintos de cluster bajo el tamaño de grupo (lag) debinido, con lambda y valores de p bajo una distribución Gamma
## 4 datos son requeridos en esta funcion. file = juego de datos con la localizacion de los palindromes,  lag = tamaño del grupo a calcular el distanciamiento entre palindromes, lambda = media esperada de palindromes para el intervalo de ADN, p = valores de probabilidad de interes, 
## regularmente nos interesamos en los datos que estan en la cola de la probabilidad de distribución.
# la salida es un objeto tipo dataframe.

location.gamma <- function( file, lag, lambda, p){
  location_gamma <-data.frame(index = 1:(296-lag),loc=file[1:(296-lag)],dist=diff(cmv, lag= lag), prob= pgamma(diff(file, lag= lag), shape=lag, rate = lambda, log=T))
  return (data.frame (Indice = location_gamma$index [location_gamma$prob < log(p)],
                      loc.O=location_gamma$loc [location_gamma$prob < log(p)],loc.lag= location_gamma$loc[location_gamma$prob < log(p)]+ lag,
                      p = signif (exp(location_gamma$prob[location_gamma$prob < log(p)]),3)))
}

## Sliding window analysis
##USAGE:
# cmv.counts<-slide.window(hcmv, 230000, blocksize=3000, incr=400)

slide.window<-function(seqdata, total.length, blocksize=500, incr=100){
  block.right<-seq(blocksize,total.length, by=incr)
  block.left<-block.right-blocksize+1
  block.mid<-(block.right+block.left)/2
  nblocks<-length(block.right)
  count<-numeric(nblocks)
  rate<-numeric(nblocks)
  
  for(i in 1:nblocks) {
    focal<-seqdata[seqdata>block.left[i] & seqdata < block.right[i]]
    count[i]<-length(focal)
    rate[i]<-count[i] / blocksize
  }
  data.frame(count, rate, block.left, block.mid, block.right)
}

Table3<- slide.window(cmv, 230000, blocksize=1500, incr=500) ## Saving the data.frame to get locations for the bigger peaks.
sorted_Table3 <- Table3[order(Table3$count, decreasing = TRUE), ]
sorted_Table3

# plot(count ~ block.mid, data=cmv.counts, type="l")
plot(count ~ block.mid, data=slide.window(cmv, 230000, blocksize=1000, incr=500), type="l", xlab="Localización en la secuencia de ADn (bp)", ylab = "Número de Palindromes")
abline (v = 91500, lty = 10, col= "blue") ## adding a line to mark the starting location in the first recomeded section to search for replication site. 
abline (v = 94000, lty = 10, col = "blue") ## adding a line to mark the ending location in the first recomeded section to search for replication site.
abline (v = 194000, lty = 10, col = "blue") ## adding a line to mark the starting location in the second recomeded section to search for replication site.
abline (v = 196500, lty = 10, col = "blue") ## adding a line to mark the starting location in the second recomeded section to search for replication site.
