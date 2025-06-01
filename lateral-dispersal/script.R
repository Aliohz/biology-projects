## cargar paquetes
library (rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 1)
library (tidybayes)

# Data
datos <- read.csv("datos.csv", skip=4)

## Convertir FFG a enteros
## Sc= 1, CF= 2, Cg = 3, Pr = 4, Sh = 5
datos$F_no <- ifelse(datos$FFG == "Sc", 1,
                     ifelse(datos$FFG == "CF", 2,
                            ifelse(datos$FFG == "Cg", 3,
                                   ifelse(datos$FFG == "Pr", 4,
                                          ifelse(datos$FFG == "Sh", 5, NA)))))

################################################################################
## Primer gráfico
frame();
plot.window(xlim=c(0, 40), ylim=c(0,110));

axis(side=1, at=seq(0, 40, by = 10), cex.axis=1);
axis(side=2, at=seq(0, 110, by= 10), cex.axis = 1); #note how you use the at argument to control location of tickmarks.

mtext("Abundancia", side = 2,line=2.5, cex=1);
mtext("Distancia (m)",side=1,cex=1, line = 2.5);

# Sc
points (datos$Distancia_m[datos$F_no == 1], datos$Abundancia[datos$F_no == 1], col= "#D81B60", pch= 16)

# CF
points (datos$Distancia_m[datos$F_no == 2], datos$Abundancia[datos$F_no == 2], col= "#1E88E5", pch= 16)

# Cg
points (datos$Distancia_m[datos$F_no == 3], datos$Abundancia[datos$F_no == 3], col= "#FFC107", pch= 16)

# Pr
points (datos$Distancia_m[datos$F_no == 4], datos$Abundancia[datos$F_no == 4], col= "#986cca", pch= 16)

# Sh
points (datos$Distancia_m[datos$F_no == 5], datos$Abundancia[datos$F_no == 5], col= "#54ad4d", pch= 16)

# Leyenda
legend(25,100, col = c("#D81B60","#1E88E5","#FFC107","#986cca","#54ad4d"), title ="Grupo funcional alimenticio",
       legend = c("Raspador (Sc)",  "Colector-filtrador (CF)", "Colector-recogedor (Cg)","Depredador (Pr)","Fragmentador (Sh)"), bty ="n", pch = c(16,16,16,16,16), cex= 1)

## Codigo en el lenguage de stan para generar el modelo bayesiano jerárquico
sink("poisson2.stan")
cat("
    data {
    
    int <lower=1> N; //numero de puntos
    int <lower=0> A[N];
    vector [N] D; // Distancia en donde fue muestreado (variable independiente)
    int  <lower=1> FF[N]; //FFG
    int  <lower=1> F_no; // No
    }
    
    parameters {
    vector [F_no] a;
    real <lower=0> sigma_a;
    real mu_a;
    
    vector [F_no] b;
    real mu_b;    
    real sigma_b;
    }
    
    model {
    //priors del intercepto. 
    mu_a ~ normal(0,10);
    sigma_a ~ normal(0,10);
    a ~ normal(mu_a, sigma_a);
    
    //priors de la pendiente 
    mu_b ~ normal(0,10);
    sigma_b ~ normal(0,10);
    b ~ normal(mu_b, sigma_b);

    //likelihood
    for(i in 1:N){
    A[i] ~ poisson(a[FF[i]]*exp(-b[FF[i]]*D[i]));
    }
    }
    
    "
    ,fill=TRUE)
sink()

## hacer una lista
abund_data <- list(A=datos$Abundancia, D=datos$Distancia_m, N=length(datos$Distancia_m), FF=datos$F_no, F_no=5)

## correr modelo

abund_fit2<-stan(file='poisson2.stan', data = abund_data, 
                 iter = 10000, chains = 4, warmup = 7500)

# revisar resultados
# revisar estimaciones

print(abund_fit2, digits_summary = 6)
## revisar la convergencia de las cadenas
traceplot(abund_fit2, pars = c("a"))
traceplot(abund_fit2, pars = c("b"))

################################################################################
## Tabla de los interceptos (a) y pendientes (b) con sus intervalos de credibilidad bayesiano
modelo <- extract(abund_fit2)
modelo

Sc <-data.frame (c(mean_a = mean(modelo$a[,1]), CI = quantile(modelo$a[,1], probs =c(0.025,0.975)), mean_b = mean(modelo$b[,1]), CI = quantile(modelo$b[,1], probs =c(0.025,0.975))))
row.names(Sc)<- c("a_estimado", "a_CI_2.5%", "a_CI_97.5%", 
                       "b_estimado", "b_CI_2.5%", "b_CI_97.5%" )
colnames(Sc)<- c("Sc")

CF <-data.frame (c(mean_a = mean(modelo$a[,2]), CI = quantile(modelo$a[,2], probs =c(0.025,0.975)), mean_b = mean(modelo$b[,2]), CI = quantile(modelo$b[,2], probs =c(0.025,0.975))))
row.names(CF)<- c("a_estimado", "a_CI_2.5%", "a_CI_97.5%", 
                       "b_estimado", "b_CI_2.5%", "b_CI_97.5%" )
colnames(CF)<- c("CF")

Cg <-data.frame (c(mean_a = mean(modelo$a[,3]), CI = quantile(modelo$a[,3], probs =c(0.025,0.975)), mean_b = mean(modelo$b[,3]), CI = quantile(modelo$b[,3], probs =c(0.025,0.975))))
row.names(Cg)<- c("a_estimado", "a_CI_2.5%", "a_CI_97.5%", 
                       "b_estimado", "b_CI_2.5%", "b_CI_97.5%" )
colnames(Cg)<- c("Cg")

Pr <-data.frame (c(mean_a = mean(modelo$a[,4]), CI = quantile(modelo$a[,4], probs =c(0.025,0.975)), mean_b = mean(modelo$b[,4]), CI = quantile(modelo$b[,4], probs =c(0.025,0.975))))
row.names(Pr)<- c("a_estimado", "a_CI_2.5%", "a_CI_97.5%", 
                       "b_estimado", "b_CI_2.5%", "b_CI_97.5%" )
colnames(Pr)<- c("Pr")

Sh <-data.frame (c(mean_a = mean(modelo$a[,5]), CI = quantile(modelo$a[,5], probs =c(0.025,0.975)), mean_b = mean(modelo$b[,5]), CI = quantile(modelo$b[,5], probs =c(0.025,0.975))))
row.names(Sh)<- c("a_estimado", "a_CI_2.5%", "a_CI_97.5%", 
                       "b_estimado", "b_CI_2.5%", "b_CI_97.5%" )
colnames(Sh)<- c("Sh")

estimados <- as.data.frame(cbind(Sc,CF,Cg,Pr,Sh))
estimados <- as.data.frame(t(estimados))

# Preparación para el siguiente bloque
estimados$F_no <- c(1,2,3,4,5)

print(estimados)

################################################################################
## Gráfica de los modelos dispersión de los distintos grupos funcionales
frame();
plot.window(xlim=c(0, 40), ylim=c(0,110));

axis(side=1, at=seq(0, 40, by = 10), cex.axis=1);
axis(side=2, at=seq(0, 110, by= 10), cex.axis = 1); #note how you use the at argument to control location of tickmarks.

mtext("Abundancia", side = 2,line=2.5, cex=1);
mtext("Distancia (m)",side=1,cex=1, line = 2.5);

# Sc
points (datos$Distancia_m[datos$F_no == 1], datos$Abundancia[datos$F_no == 1], col= "#D81B60", pch= 16)

Sc_fitted <- as.data.frame(sort(datos$Distancia_m[datos$F_no==1]))
colnames(Sc_fitted) <- c("Distancia_m")
Sc_fitted$abundancia <-  estimados$a_estimado[estimados$F_no==1] * exp(-(estimados$b_estimado[estimados$F_no==1]) * Sc_fitted$Distancia_m)
  
points (Sc_fitted$Distancia_m, Sc_fitted$abundancia,
        col= "#D81B60", pch= 16, type= "l", lwd = 2)

# CF
points (datos$Distancia_m[datos$F_no == 2], datos$Abundancia[datos$F_no == 2], col= "#1E88E5", pch= 16)

Sc_fitted <- as.data.frame(sort(datos$Distancia_m[datos$F_no==2]))
colnames(Sc_fitted) <- c("Distancia_m")
Sc_fitted$abundancia <-  estimados$a_estimado[estimados$F_no==2] * exp(-(estimados$b_estimado[estimados$F_no==2]) * Sc_fitted$Distancia_m)

points (Sc_fitted$Distancia_m, Sc_fitted$abundancia,
        col= "#1E88E5", pch= 16, type= "l", lwd = 2)

# Cg
points (datos$Distancia_m[datos$F_no == 3], datos$Abundancia[datos$F_no == 3], col= "#FFC107", pch= 16)

Sc_fitted <- as.data.frame(sort(datos$Distancia_m[datos$F_no==3]))
colnames(Sc_fitted) <- c("Distancia_m")
Sc_fitted$abundancia <-  estimados$a_estimado[estimados$F_no==3] * exp(-(estimados$b_estimado[estimados$F_no==3]) * Sc_fitted$Distancia_m)

points (Sc_fitted$Distancia_m, Sc_fitted$abundancia,
        col= "#FFC107", pch= 16, type= "l", lwd = 2)

# Pr
points (datos$Distancia_m[datos$F_no == 4], datos$Abundancia[datos$F_no == 4], col= "#986cca", pch= 16)

Sc_fitted <- as.data.frame(sort(datos$Distancia_m[datos$F_no==4]))
colnames(Sc_fitted) <- c("Distancia_m")
Sc_fitted$abundancia <-  estimados$a_estimado[estimados$F_no==4] * exp(-(estimados$b_estimado[estimados$F_no==4]) * Sc_fit4ted$Distancia_m)

points (Sc_fitted$Distancia_m, Sc_fitted$abundancia,
        col= "#986cca", pch= 16, type= "l", lwd = 2)

# Sh
points (datos$Distancia_m[datos$F_no == 5], datos$Abundancia[datos$F_no == 5], col= "#54ad4d", pch= 16)

Sc_fitted <- as.data.frame(sort(datos$Distancia_m[datos$F_no==5]))
colnames(Sc_fitted) <- c("Distancia_m")
Sc_fitted$abundancia <-  estimados$a_estimado[estimados$F_no==5] * exp(-(estimados$b_estimado[estimados$F_no==5]) * Sc_fitted$Distancia_m)

points (Sc_fitted$Distancia_m, Sc_fitted$abundancia,
        col= "#54ad4d", pch= 16, type= "l", lwd = 2)

# Leyenda
legend(25,100, col = c("#D81B60","#1E88E5","#FFC107","#986cca","#54ad4d"), title ="Grupo funcional alimenticio",
       legend = c("Raspador (Sc)",  "Colector-filtrador (CF)", "Colector-recogedor (Cg)","Depredador (Pr)","Fragmentador (Sh)"), bty ="n", pch = c(16,16,16,16,16), cex= 1)

################################################################################
## Matrix de pendientes para comparar las tasas de dispersión

dispersion <- matrix(nrow = nrow(estimados), ncol = nrow(estimados))

# Nombres de filas y columnas
rownames(dispersion) <- rownames(estimados)
colnames(dispersion) <- rownames(estimados)

# Llenar la matriz con los ratios b_i/b_j
for(i in 1:nrow(estimados)) {
  for(j in 1:nrow(estimados)) {
    dispersion[i, j] <- estimados$b_estimado[i] / estimados$b_estimado[j]
  }
}

# Convertir a dataframe
dispersion <- as.data.frame(dispersion)
dispersion

################################################################################
## Distancia media de dispersion es 1/b
distancia.lateral.media <- 1 / estimados$b_estimado
names(distancia.lateral.media) <- rownames(estimados)
print(distancia.lateral.media)
