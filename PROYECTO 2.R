#####################################################
#####################################################
# Proyecto 2, Estadistica - Alejandro Borda 202020727
#####################################################
#####################################################


#---------------------------------------------
# PARTE I: NEWTON RAPHSON
#---------------------------------------------

#DEFINIMOS l(theta), sus derivadas parciales y la funcion f a resolver

log_likelihood = function(par1, par2, data) {
  shape = par1
  scale = par2
  n = length(data)
  sum(dgamma(data, shape, scale, log = TRUE))
}

f = function(alpha_k, beta_k, data){
  n = length(data)
  f1 = sum(log(data) - digamma(alpha_k) - log(beta_k))
  f2 = sum(data/beta_k^2 - alpha_k/beta_k)
  return (c(f1, f2))
}

# IMPLEMENTAMOS NEWTON-RAPHSON

Df_inverse = function(alpha, beta, data) {
  
  n = length(data)
  
  f1_alpha = function(alpha, beta, n) {
    -n * trigamma(alpha)
  }
  f1_beta = function(alpha, beta, n) {
    -n / beta
  }
  f2_alpha = function(alpha, beta, n) {
    -n / beta
  }
  f2_beta = function(alpha, beta, n) {
    -2*n*mean(data)/beta^3 + alpha*n/beta^2
  }
  
  # Creamos Df
  matrix_Df = matrix(c(f1_alpha(alpha, beta, n),
                       f1_beta(alpha, beta, n),
                       f2_alpha(alpha, beta, n),
                       f2_beta(alpha, beta, n)), nrow = 2, byrow = TRUE)
  
  # Inversa de Df
  matrix_inverse = solve(matrix_Df)
  return(matrix_inverse)
}


NR = function(alpha_k, beta_k, data) {
  if (!all(is.numeric(data))) {
    stop("Data vector contains non-numeric values.")
  }
  
  tol = 1e-6  # Threshold de convergencia (recomendado)
  max_iter = 500 
  iter = 0  
  converged = FALSE 
  
  while (!converged && iter < max_iter) {
    
    n = length(data)

    Dfinv = Df_inverse(alpha_k, beta_k, data)
    coord_n = c(alpha_k, beta_k) - Dfinv %*% f(alpha_k, beta_k, data)
    
    # Revisamos convergencia
    if (!any(is.na(coord_n))) {
      if (max(abs(coord_n - c(alpha_k, beta_k))) < tol) {
        converged = TRUE
      } else {
        alpha_k = coord_n[1]
        beta_k = coord_n[2]
        iter = iter + 1
      }
    } else {
      stop("Convergence failed: NaN values encountered.")
    }
    print(paste("Iteration:", iter))
  }
  
  return(coord_n)
}


######################################################
#--------PARTE II GRAFICAS DE LOS ESTIMADORES---------
######################################################

alpha_v = 2.5
beta_v = 4.6
n_list = list(100, 200, 500, 1000, 2000)
m = 500
moments = list(list(), list(), list(), list(), list()) #Cada lista interior corresponde
EMVs = list(list(), list(), list(), list(), list()) # a un n fijo

for (j in 1:5) { # iteramos sobre el numero de n's 
  
  n = n_list[[j]] 
  EMV_values = list()
  moments_values = list()
  
  for (i in 1:m) { # 500 iteraciones de un n fijo
    gammaData = rgamma(n, shape=alpha_v, scale=beta_v)
    est_alpha = mean(gammaData)^2 / var(gammaData) #alpha y beta por met. de momentos
    est_beta = var(gammaData) / mean(gammaData)
    est=list(est_alpha, est_beta)
    moments_values = c(moments_values, list(est))
    
    # Por NR
    EMV = NR(est_alpha, est_beta, gammaData)
    EMV_values = c(EMV_values, list(EMV))
  }
  
  EMVs[[j]] = EMV_values
  moments[[j]] = moments_values
}

###GRAFICAMOS para EMV###

par(mfrow=c(3, 2))
for (j in 1:5) {
  
  n = n_list[j]
  EMV_values = EMVs[[j]]
  x_values = sapply(EMV_values, "[[", 1)
  y_values = sapply(EMV_values, "[[", 2)
  
  plot(x_values, y_values, type="p", col="black", xlab="Alpha", ylab="Beta", main=paste("Sample Size:", n))
  points(alpha_v, beta_v, col = "blue", pch = 16)
}

#Boxplot
alpha_values = list()
beta_values = list()

for (j in 1:5) {
  
  EMV_values = EMVs[[j]]
  alpha_j = sapply(EMV_values, "[[", 1)
  beta_j = sapply(EMV_values, "[[", 2)
  
  alpha_values[[j]] = alpha_j
  beta_values[[j]] = beta_j
}

par(mfrow=c(2, 1))
boxplot(alpha_values, main="Boxplot of Alpha", xlab="Sample Size", ylab="Alpha")
boxplot(beta_values, main="Boxplot of Beta", xlab="Sample Size", ylab="Beta")



###GRAFICAMOS POR MET DE MOMENTS###

par(mfrow=c(3, 2))

for (j in 1:5) { 
  
  n = n_list[j]
  moment_values = moments[[j]]
  x_values = sapply(moment_values, "[[", 1)
  print(x_values)
  y_values = sapply(moment_values, "[[", 2)
  
  plot(x_values, y_values, type="p", col="orange", xlab="Alpha", ylab="Beta", main=paste("Sample Size:", n))
  points(alpha_v, beta_v, col = "blue", pch = 16)
}

#Boxplot
alpha_values = list()
beta_values = list()

for (j in 1:5) {  
  moment_values = moments[[j]]
  
  alpha_j = sapply(moment_values, "[[", 1)
  beta_j = sapply(moment_values, "[[", 2)
  
  alpha_values[[j]] = alpha_j
  beta_values[[j]] = beta_j
}

par(mfrow=c(2, 1))

boxplot(alpha_values, main="Boxplot of Alpha", xlab="Sample Size", ylab="Alpha")
boxplot(beta_values, main="Boxplot of Beta", xlab="Sample Size", ylab="Beta")



###################################################################
#-------------- PARTE III: Conjuntos de Confianza---------------
###################################################################

#Informacion de Fisher
If = function(alpha, beta) {
  matrix_If = matrix(c(trigamma(alpha),
                       1/beta,
                       1/beta,
                       alpha/beta^2), nrow = 2, byrow = TRUE)
  return (matrix_If)
  
}

theta_0 = list(alpha_v, beta_v) #Vector de valores reales

# Definimos la funcion que calcula el porcentaje de datos que caen en el set C

C = function(theta_0, data, n) {
  rho = qchisq(0.95, 2)
  counter = 0 
  
  for (i in 1:length(data)) {
    
    d = c(data[[i]][1], data[[i]][2]) #Formateamos los datos para calcular el rango
    thetaVec = c(theta_0[[1]], theta_0[[2]])
    F_info = If(data[[i]][1], data[[i]][2])
    mag = n*t(d - thetaVec) %*% F_info %*% (d - thetaVec) #Valor para revisar <= rho

    if (mag <= rho) { #Condicion de estar en C
      counter = counter + 1 
    }
  }
  
  return (counter / length(data))
}

print("Probabilidades para distintos n: ")
cat("Para n=100 la probabilidad es: ", C(theta_0, EMVs[[1]], 100), "\n")
cat("Para n=200 la probabilidad es: ", C(theta_0, EMVs[[2]], 200), "\n")
cat("Para n=500 la probabilidad es: ", C(theta_0, EMVs[[3]], 500), "\n")
cat("Para n=1000 la probabilidad es: ", C(theta_0, EMVs[[4]], 1000), "\n")
cat("Para n=2000 la probabilidad es: ", C(theta_0, EMVs[[5]], 2000), "\n")


#AREA DEL ELISPE: 

AE = function(theta_0, EMVs) {
  
  areas = list(list(), list(), list(), list(), list())
  rho = qchisq(0.95, 1)
  
  for (n in n_list){
    ind = which(n_list == n)
    data = EMVs[[ind]]
    A_values = list()
  
    for (i in 1:length(data)) {
      d = c(data[[i]][1], data[[i]][2]) #Formateamos los datos para calcular el rango
      thetaVec = c(theta_0[[1]], theta_0[[2]])
      F_info = If(theta_0[[1]], theta_0[[2]])
    
      U = eigen(F_info)$vectors #Obtenemos U para tener la variable Y
      lambda = eigen(F_info)$values
      D = diag(lambda)
      Y = t(U) %*% (d - thetaVec)
      cov_Y = cov(Y)
      eigen_Y = eigen(cov_Y)
      eigenvalues = eigen_Y$values
      eigenvectors = eigen_Y$vectors
        
      a = sqrt(max(eigenvalues)) #Semiejes de la elipse
      b = sqrt(min(eigenvalues))
      Area = pi * a * b #Area de la elipse
      A_values = c(A_values, list(Area)) #Agregamos a la lista
      
    }
    areas[[ind]] = A_values
  }
  return(areas)
  
}

#BOXPLOTS

areas = AE(theta_0, EMVs)
area_values = list()

for (j in 1:5) {  
  moment_values = areas[[j]]
  alpha_j = sapply(moment_values, "[[", 1)
  area_values[[j]] = alpha_j
}
par(mfrow=c(1, 1))
boxplot(area_values, main="Boxplot of Areas", xlab="Sample Size", ylab="area")


###################################################################
#------------ PARTE IV: Convergencia de la covarianza--------------
###################################################################


If_inv = solve(If(alpha_v, beta_v)) #Matriz inversa de If(theta_0)

covariances = list()
for (n in n_list){
  ind = which(n_list == n)
  EMV_val = unlist(EMVs[[ind]])
  x_list = list()
  y_list = list()
  
  for (k in 1:500){
    x = c(EMV_val[(k - 1) * 2 + 1]) #Extraemos alpha y beta correspondiente
    y = c(EMV_val[(k - 1) * 2 + 2])
    
    x_list[[k]] = x
    y_list[[k]] = y
    
    
  }
  x_matrix = matrix(c(unlist(x_list), unlist(y_list)), ncol = 2, byrow = FALSE)
  covariances[[ind]] = n*cov(x_matrix) #n*Cov a comparar con If(theta_0) inversa
  
}

print(covariances)
print(If_inv)

























