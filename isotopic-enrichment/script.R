# Cargar datos desde el archivo csv.
datos <- read.csv("./datos.csv")


#  ¿Qué organismo tiene el mayor valor isotopico de δ15N?
max_dN <- max(datos$d15N)

cat("El organismo con el mayor valor isotopico de δ15N es", datos$Taxon[datos$d15N == max_dN], "del gremio", datos$Guild[datos$d15N == max_dN], ".")

  
#  ¿Hay correlacion en la variación de δ15N y δ13C?

plot(datos$d13C, datos$d15N)

shapiro.test(datos$d13C)
shapiro.test(datos$d15N)

cor(datos$d13C, datos$d15N, method = "spearman")
cor.test(datos$d13C, datos$d15N, method = "spearman") #  suppressWarnings() can be used # No se rechaza Ho, no hay correlación.

cor.test(datos$d13C, datos$d15N)

#  ¿Cuál es el valor promedio y variación de cada isotopo por cada taxa y por cada gremio alimenticio?

taxa_sd <- taxa_medias <- gremio_medias <- gremio_sd <- c()

for (isotopo in c('d13C', 'd15N')) {
  
  for (taxa_no in unique(datos$taxa_no)) {
    media <- mean(datos[[isotopo]][datos$taxa_no == taxa_no])
    desv_est <- sd(datos[[isotopo]][datos$taxa_no == taxa_no])
    
    taxa_medias <- append(taxa_medias, media)
    taxa_sd <- append(taxa_sd, desv_est)
    
    #cat("El promedio de taxa", taxa_no, "para el isotopo", isotopo, "es:", media, "\n")
  }
  #cat("\n")
  
  for (guild in unique(datos$Guild)) {
    print(guild)
    media <- mean(datos[[isotopo]][datos$Guild == guild])
    desv_est <- sd(datos[[isotopo]][datos$Guild == guild])
    
    gremio_medias <- append(gremio_medias, media)
    gremio_sd <- append(gremio_sd, desv_est)
    
    #cat("El promedio del gremio", guild, "para el isotopo", isotopo, "es:", media, "\n")
  }
  #cat("\n")
  
}


taxa_df <- data.frame(
  taxa_no <- c(unique(datos$taxa_no)),
  Media_d13C <- taxa_medias[1:10],
  SD_d13C <- taxa_sd[1:10],
  Media_d15N <- taxa_medias[11:20],
  SD_d15N <- taxa_sd[11:20]
)

nombres<- c("taxa_no","Media_d13C","SD_d13C","Media_d15N","SD_d15N")
colnames(taxa_df)<- nombres

taxa_df

gremio_df <- data.frame(
  gremio <- unique(datos$Guild),
  Media_d13C <- gremio_medias[1:4],
  SD_d13C <- gremio_sd[1:4],
  Media_d15N <- gremio_medias[5:8],
  SD_d15N <- gremio_sd[5:8]
)

nombres<- c("Gremio","Media_d13C","SD_d13C","Media_d15N","SD_d15N")
colnames(gremio_df)<- nombres

gremio_df

#  ¿Qué gremios alimenticios tienen las las mayores diferencias entre los valores isotopicos de δ15N?

min_d15N <- min(gremio_df$Media_d15N)
max_d15N <- max(gremio_df$Media_d15N)

cat("Los gremios con mayor diferencia en sus valores isotópicos promedio de δ15N son:", gremio_df$Gremio[gremio_df$Media_d15N == min_d15N], 
    "(", min_d15N,  ") y", gremio_df$Gremio[gremio_df$Media_d15N == max_d15N], "(", max_d15N, ").")

#  ¿Cuanto mayor es el valor isotopico δ13C, en promedio, de los gremios tróficos con menor y mayor valor isotopico de δ15N?

proporcion_13C <- gremio_df$Media_d13C[gremio_df$Gremio == "Detritivore-A"] / gremio_df$Media_d13C[gremio_df$Gremio == "Omnivore-A"]
proporcion_13C

proporcion_15N <- gremio_df$Media_d15N[gremio_df$Gremio == "Omnivore-A"] / gremio_df$Media_d15N[gremio_df$Gremio == "Detritivore-A"]
proporcion_15N