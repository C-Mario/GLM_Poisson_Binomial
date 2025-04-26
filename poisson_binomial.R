setwd("C:/PERSONAL_MAPU/UNIVERSIDADES/NACIONAL/MODELOS_LINEALES_GENERALIZADOS")
rm(list = ls())
check_and_install <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    library(package_name, character.only = TRUE)
    message(paste("Paquete", package_name, "instalado y cargado."))
  } else {
    message(paste("Paquete", package_name, "ya está instalado y cargado."))
  }
}
# Cargar librerías necesarias
check_and_install("tibble")
check_and_install("dplyr")
check_and_install("readr")

## Población Poisson ####
players_202223 <- read_csv("players_202223.csv")
View(players_202223)
baloncesto_data <- players_202223 %>% select(AGE,MIN,GP,FGM,FGA,FG3M,FG3A,FTM,
                                               FTA,PTS)
View(baloncesto_data)

# Matriz X
x1 <- as.matrix(baloncesto_data %>% select(-last_col()))
X <- model.matrix(~x1)
Y <- baloncesto_data$PTS

set.seed(1040)


# Creación de la función
fisher_scoring_poisson <- function(y, X, beta_init, tol = 1e-4, max_iter = 100) {
  beta <- beta_init
  for (iter in 1:max_iter) {
    eta <- X %*% beta
    mu <- exp(eta)
    
    W <- diag(as.vector(mu))
    
    z <- eta + (y - mu) / mu
    beta_new <- solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% z)
    
    if (max(abs(beta_new - beta)) < tol) {
      message("Convergió en iteración ", iter)
      return(as.vector(beta_new))
    }
    
    beta <- beta_new
  }
  warning("No convergió")
  return(as.vector(beta))
}

# Ejecutamos
fisher_scoring_poisson(y = Y,X = X,beta_init = rep(0.01, ncol(X)),
                       max_iter = 100)
beta_est <- fisher_scoring_poisson(y = y, X = X, beta_init = c(0, 0),
                                   max_iter = 500)
beta_est


## Población binomial
# Parámetros
set.seed(1040)
n <- 100
beta <- c(1, 0.5)

# Simulación de datos
xbin <- runif(n, min = 0, max = 10)
X2 <- model.matrix(~ xbin)  # Incluye intercepto automáticamente
eta <- X2 %*% beta

# Dataset simulado
pi <- exp(eta)/(1+exp(eta))
yb <- c()
for (i in 1:n) {
  yb[i] <- rbinom(1,size = 1,prob = pi[i])
}

data_bin <- tibble(y = yb, x = x1)

# Creación de la función
fisher_scoring_binomial <- function(y, X, beta_init, tol = 1e-4, max_iter = 100) {
  beta <- beta_init
  for (iter in 1:max_iter) {
    eta <- X %*% beta
    eta <- pmin(eta, 700)
    eta <- pmax(eta, -700)
    
    pi <- exp(eta) / (1 + exp(eta))
    
    W <- diag(as.vector(pi * (1 - pi)))
    
    z <- eta + (y - pi) / (pi * (1 - pi))
    
    beta_new <- tryCatch({
      solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% z)
    }, error = function(e) {
      stop("Fallo al invertir la matriz. ¿X tiene multicolinealidad?")
    })
    
    if (max(abs(beta_new - beta), na.rm = TRUE) < tol) {
      message("Convergió en iteración ", iter)
      return(as.vector(beta_new))
    }
    
    beta <- beta_new
  }
  warning("No convergió")
  return(as.vector(beta))
}

beta_binomial <- fisher_scoring_binomial(y = yb, X = X2, 
                                         beta_init = rep(0.01, ncol(X2)))
