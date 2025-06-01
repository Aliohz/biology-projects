
    data {
    
    int <lower=1> N; // numero de puntos
    int <lower=0> A[N]; // Abundancia (variable respuesta)
    vector [N] D; // Distancia en donde fue muestreado (variable independiente)
    int  <lower=1> FF[N]; //FFG, grupos funcionales
    int  <lower=1> F_no; // Número de grupos funcionales
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
    //priors del intercepto, no son informativos
    mu_a ~ normal(45,50);
    sigma_a ~ normal(0,10);
    a ~ normal(mu_a, sigma_a);
    
    //priors de la pendiente, no son informativos
    mu_b ~ normal(0,0.1);
    sigma_b ~ normal(0,0.1);
    b ~ normal(mu_b, sigma_b);

    //likelihood
    for(i in 1:N){
    A[i] ~ poisson(a[FF[i]]*exp(-b[FF[i]]*D[i])); // el modelo es y = a * exp(-b * x)
    }
    }
    
    
