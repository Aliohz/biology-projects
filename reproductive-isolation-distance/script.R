library("hierfstat") # Para calcular la matriz de distancias genéticas, Fst
library("adegenet") # Para importar las secuencias genéticas del archivo FASTA alineado
library("vegan") # Para hacer el test de mantel
library("geosphere") # Para calcular la matriz de distancias geográficas
library("MASS") # Para agregar la densidad de puntos en el gráfico de dispersión

# Se carga el archivo fasta con las 166 secuencias alineadas y se extraen los SNPs (Single nucleotide polymorphisms)
seq.snp = fasta2DNAbin('./mtDNA_aligned.fas')
obj = DNAbin2genind(seq.snp) # Se almacena en un objeto genind propio del paquete adegenet.

# Se carga la información geográfica y el nombre de los sitios muestreados. Se asume que cada coordenada corresponde a una subpoblación. 
meta <- read.table('./meta.csv', sep=',', header = T)

# Se agrega al objeto genind, pero no supe cómo trabajar las coordenadas dentro de este objeto para obtener las distancias geográficas.
obj$pop = as.factor(meta[match(rownames(obj$tab), meta$Genbank_accession),]$Site)

# Se convierte el objeto genind al objeto hierfstat para calcular el Fst.
obj.hf = genind2hierfstat(obj)

# Se calcula la matriz de Fst por el métdodo “Dch”: Cavalli-Sforza and Edwards Chord distance. This distance is used as default since Takezaki & Nei (1996) found that it was the best to retrieve the relation among samples.
# Como se trata de información genética mitocondrial, el atributo diploid se cambia a Falso, para trabajar con datos haploides.
gen_dist <- genet.dist(obj.hf, method='Dch', diploid=F)

# Calculo de la matriz de distancias geográficas.
# Se obtienen todos los sitios únicos (subpoblaciones) y se ordenan alfabéticamente, al igual que aparecen en la matriz de distancias genéticas.
sites <- unique(meta[, c("Site", "Lat", "Long")])
sites <- sites[order(sites$Site), ] # Este paso es CRÍTICO, dado que la correlación de Mantel revisa las posiciones de la matriz, asumiendo que tienen el mismo orden.

# Se obtiene la matriz de distancias en m.
geo_dist <- distm(sites[, c("Long", "Lat")], fun = distHaversine)
geo_dist <- geo_dist / 1000 # Se convierte a km.

# Cambia los nombres de filas y columnas de la matriz de distancias.
rownames(geo_dist) <- sites$Site
colnames(geo_dist) <- sites$Site

# Se convierte a una matriz que contiene los datos únicamente en la parte triangular de abajo, que es la que se usa en la prueba de Mantel.
geo_dist <- as.dist(geo_dist)

# Test de Mantel
mantel_test <- mantel(gen_dist, geo_dist, method = "pearson", permutations = 9999)
mantel_test

# Histograma de los resultados de Mantel
# Extrae las permunaciones para graficar la distribución
null_distribution <- mantel_test$perm
#null_distribution

# Gráfico con el valor calculado, las permutaciones y el IC de 95%
hist(null_distribution, 
     breaks = 30, 
     col = "#afc4b4", 
     main = "Distribución nula del test de Mantel",
     xlab = "Estadístico de Mantel (r)",
     ylab = "Frequencia",
     xlim = range(c(null_distribution, mantel_test$statistic)))

abline(v = mantel_test$statistic, col = "#9c2c33", lwd = 2) # r calculado

quantiles <- quantile(null_distribution, probs = c(0.025, 0.975)) # IC 95%
abline(v = quantiles, col = "#362d96", lwd = 2, lty = 2)

legend("topright", 
       legend = c(paste0("r = ", round(mantel_test$statistic, 3)), paste0("p = ", mantel_test$signif), "I.C. 95% del modelo nulo"), 
       col = c("#9c2c33", NA, "#362d96"),  
       lwd = c(2, NA, 2),
       lty = c(1, NA, 2),
       seg.len = 1,
       bty = "n")

# Gráfico de dispersión de la distancia geográfica y distancia genética

# Calcula la densidad de puntos dentro del gráfico, permite saber donde está concentrado la mayoría de nuestros datos.
dens <- kde2d(as.vector(geo_dist), as.vector(gen_dist), n=300) 

# Asigna diferentes colores para distinguir las densidades, de blanco a rojo que equivale de menor a mayor densidad.
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red")) 
plot(geo_dist, gen_dist, pch=20,cex=.5,
     xlab = "Distancia geográfica (km)",
     ylab = "Distancia genética (Fst)")
image(dens, col=transp(myPal(300),.7), add=TRUE)
title("Aislamiento por distancia")

# Regresión lineal entre las distancias geográficas y genéticas.
modelo <- lm(as.vector(gen_dist) ~ as.vector(geo_dist))
abline(modelo)
modelo