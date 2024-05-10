########################################################
#
#Proyecto 1 Estadistica
# Alejandro Borda y Sebastian Suescun
#
########################################################




########################################################
#PARTE 1, RAYLEIGH F^-1
########################################################

# Constantes y densidades

n=10000
thetaList=list(2.7,7.5)

FR=function(x,theta){1 - exp(-x^2 / (2 * theta^2))}
U=runif(n,min=0,max=1) #Muestra uniforme

#Calculamos inversa explicitamente
find_root = function(u,theta) {
  uniroot(function(x) FR(x,theta) - u, c(0, 100))$root
}

#Calculamos la muestra por el metodo F^-1 para ambos valores theta
FR_inverse1=sapply(U, find_root, theta = thetaList[[1]])
FR_inverse2=sapply(U, find_root, theta=thetaList[[2]])


#Definiendo la distribucion Rayleigh REAL

CDF_rayleigh = function(x, theta){1 - exp(-x^2 / (2 * theta^2))}

q_rayleigh = function(p, theta){q = sqrt(-2 * theta^2 * log(1 - p))
return(q)}

theo_quantiles1=q_rayleigh(ppoints(n), thetaList[[1]])
theo_quantiles2=q_rayleigh(ppoints(n), thetaList[[2]])

#HISTOGRAMAS

hist(FR_inverse1, freq = FALSE, main = "Hist. de densidad y densidad real Rayleigh(2.7)", xlab = "Rayleigh theta=2.7", ylab = "Densidad")
curve((x / thetaList[[1]]^2) * exp(-x^2 / (2 * thetaList[[1]]^2)), col = "black", lwd = 2, add = TRUE, n = 1000)

hist(FR_inverse2, freq = FALSE, main = "Hist. de densidad y densidad real Rayleigh(7.5)", xlab = "Rayleigh theta=7.5", ylab = "Densidad")
curve((x / thetaList[[2]]^2) * exp(-x^2 / (2 * thetaList[[2]]^2)), col = "black", lwd = 2, add = TRUE, n = 1000)

#QQ-PLOTS
qqplot(theo_quantiles1,FR_inverse1,
       main = "Q-Q Plot - Rayleigh(2.7)",
       xlab = "Cuantiles teoricos",
       ylab = "Cuantiles metodo F^-1")


qqplot(theo_quantiles2,FR_inverse2,
       main = "Q-Q Plot - Rayleigh(7.5)",
       xlab = "Cuantiles teoricos",
       ylab = "Cuantiles metodo F^-1s")


#########################################################
#PARTE 2 
#Convergencia de Histogramas normas L1, L2, Linf
#########################################################

#-------------------------------------------------------
#CAUCHY (0,1)
#-------------------------------------------------------


n_list =list(1000, 2000, 5000, 10000, 20000)
#Listas para guardar y promediar las diferentes normas en cada iteracion
l1_list =list()
l2_list= list()
linf_list = list()
l1_ave=list()
l2_ave=list()
linf_ave=list()

for (n in n_list) { #Iteramos sobre los tamaños de muestra
  l1_list=list()
  l2_list=list()
  linf_list=list()
  for (i in 1:20) { #Realizamos las 20 muestras i.i.d para promediar

    #Generar muestra para n fijo
    cauchyData=rcauchy(n,0,1)
    cauchyData = cauchyData[cauchyData >= -12.7062 & cauchyData <= 12.7062]#Filtramos los datos Cauchy
  
    breaks=sqrt(n)
    hist_data = histogram = hist(cauchyData, 
                                  seq(min(cauchyData), 
                                      max(cauchyData), 
                                      length.out = breaks), 
                                  prob = "TRUE",
                                  plot = "TRUE")
    lines(density(cauchyData))
  
    interval_starts = hist_data$breaks[-length(hist_data$breaks)]
    interval_ends = hist_data$breaks[-1]
  
    density_values=hist_data$density #vector de densidades de cada intervalo
  
  
    #Obtener comeinzo + 4 datos intermedios por intervalo de histograma y densidad
  
    delta=0.2
    #delta corresponde al peso de cada separacion del intervalo, como tomamos 5
    # elementos por intervalo cada uno tiene peso 1/5
  
    #VECTOR DEL HISTOGRAMA:
    #Dado que el histograma es constante en un intervalo para generar el vector
    #ordenado simplemente repetimos el valor de densidad 5 veces por intervalo
    hist_density_values = numeric()
    for (i in seq_along(density_values)) {
      hist_density_values = c(hist_density_values, rep(density_values[i], each = 5))
    }
    #VECTOR DE DENSIDAD REAL:
    #Calculamos explicitamente la densidad real para cada uno de los 5 puntos equidistantes del intervalo 
    #del histograma y formamos un vector ordenado para calcular la distancia de ambos vectores a continuacion
    density_vector = numeric()
    for (i in seq_along(interval_starts)) {
      start = interval_starts[i]
      end = interval_ends[i]
      d = end - start  # Longitud del intervalo
      q1 = start + 0.2 * d
      q2 = start + 0.4 * d
      q3 = start + 0.6 * d
      q4 = start + 0.8 * d
      tfactor=pcauchy(12.7062)-pcauchy(-12.7062) #factor de la funcion de densidad truncada
      density_q0 = dcauchy(start,0,1)/tfactor
      density_q1 = dcauchy(q1,0,1)/tfactor
      density_q2 = dcauchy(q2,0,1)/tfactor
      density_q3 = dcauchy(q3,0,1)/tfactor
      density_q4 = dcauchy(q4,0,1)/tfactor
      density_vector=c(density_vector, list(density_q0, density_q1, density_q2, density_q3, density_q4))
      
    }
  
    #Calculamos las diferentes normas
    
    #L1
    histogram=hist_data
    density_vector = unlist(density_vector)
    L1_vec=abs((density_vector - hist_density_values))*(histogram$breaks[2] - histogram$breaks[1])/5
    L1_norm=sum(L1_vec)
    #L INF
    linf_norm = max(L1_vec)
    #L2
    L2_vec=((density_vector - hist_density_values)^2)*(histogram$breaks[2] - histogram$breaks[1])/5
    L2_norm=sum(L2_vec)
    L2_norm=sqrt(L2_norm)
    
    #Formateamos los datos obtenidos
    l1_list=c(l1_list, L1_norm)
    l2_list=c(l2_list, L2_norm)
    linf_list=c(linf_list, linf_norm)
  }
  #Guardar y formatear normas obtenidas por iteracion de n
  l1_list=unlist(l1_list)
  l2_list=unlist(l2_list)
  linf_list=unlist(linf_list)
  ave1=sum(l1_list)/20
  ave2=sum(l2_list)/20
  aveInf=sum(linf_list)/20
  l1_ave=c(l1_ave, ave1)
  l2_ave=c(l2_ave, ave2)
  linf_ave=c(linf_ave, aveInf)

}
#Resultados CAUCHY(0,1)
print("Normas promediadas l1, l2 y linf para respectivos n y 20 muestras iid por n : ")
print("L1: ")
print(unlist(l1_ave))
print("L2: ")
print(unlist(l2_ave))
print("Linf: ")
print(unlist(linf_ave))


#log-log Plot
#Promedios calculados para las diferentes normas
n_list=unlist(n_list)
l1_ave=unlist(l1_ave)
l2_ave=unlist(l2_ave)
linf_ave=unlist(linf_ave)

#Transformadas a logaritmos para L1, regresion lineal y graficar
log_n_list = log(n_list)
log_l1_ave = log(l1_ave)

lm_model_l1 = lm(log_l1_ave ~ log_n_list)

plot(log_n_list, log_l1_ave, 
     xlab = "log(n)", ylab = "log(L1 promedio)", 
     main = "Log-Log Plot para Cauchy(0,1) en L1")
abline(lm_model_l1, col = "red")

eq_l1 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l1)[1], digits = 3)) + .(format(coef(lm_model_l1)[2], digits = 3)) * italic(x)))
text(x = 7, y = -2.7, labels = eq_l1, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

#Mismo proceso para L2
log_l2_ave = log(l2_ave)
lm_model_l2 = lm(log_l2_ave ~ log_n_list)

plot(log_n_list, log_l2_ave, 
     xlab = "log(n)", ylab = "log(L2 promedio)", 
     main = "Log-Log Plot para cauchy(0,1) en L2")
abline(lm_model_l2, col = "blue")

eq_l2 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l2)[1], digits = 3)) + .(format(coef(lm_model_l2)[2], digits = 3)) * italic(x)))
text(x = 7, y = -3.6, labels = eq_l2, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

# Mismo proceso para Linf
log_linf_ave = log(linf_ave)

lm_model_linf = lm(log_linf_ave ~ log_n_list)

plot(log_n_list, log_linf_ave, 
     xlab = "log(n)", ylab = "log(L_inf promedio)", 
     main = "Log-Log Plot para cauchy(0,1) en Linf")

abline(lm_model_linf, col = "green")

eq_l3 = as.expression(bquote(italic(y) == .(format(coef(lm_model_linf)[1], digits = 3)) + .(format(coef(lm_model_linf)[2], digits = 3)) * italic(x)))
text(x = 7, y = -6.5, labels = eq_l3, pos = 4, cex = 0.8) #Agregar ecuacion al grafico


#---------------------------------------------------------
#Gamma(1,3)
#---------------------------------------------------------

n_list =list(1000, 2000, 5000, 10000, 20000)
#Listas para guardar y promediar las diferentes normas en cada iteracion
l1_list =list()
l2_list= list()
linf_list = list()
l1_ave=list()
l2_ave=list()
linf_ave=list()

for (n in n_list) {
  l1_list=list()
  l2_list=list()
  linf_list=list()
  for (i in 1:20) {
    
    #Generar muestra para n fijo
    gammaData=rgamma(n,1,3)
    
    hist_data = hist(gammaData, freq=FALSE, breaks=1/4*sqrt(n), plot = TRUE)
    lines(density(gammaData))
    
    interval_starts = hist_data$breaks[-length(hist_data$breaks)]
    interval_ends = hist_data$breaks[-1]
    
    density_values=hist_data$density #vector de densidades de cada intervalo
    
    
    #Obtener comeinzo + 4 datos intermedios por intervalo de histograma y densidad
    
    delta=0.2
    #delta corresponde al peso de cada separacion del intervalo, como tomamos 5
    # elementos por intervalo cada uno tiene peso 1/5
    
    #VECTOR DEL HISTOGRAMA:
    #Dado que el histograma es constante en un intervalo para generar el vector
    #ordenado simplemente repetimos el valor de densidad 5 veces por intervalo
    hist_density_values = numeric()
    for (i in seq_along(density_values)) {
      hist_density_values = c(hist_density_values, rep(density_values[i], each = 5))
    }
    
    
    #VECTOR DE DENSIDAD REAL:
    
    #Calculamos explicitamente la densidad real para cada uno de los 5 puntos equidistantes del intervalo 
    #del histograma y formamos un vector ordenado para calcular la distancia de ambos vectores a continuacion
    
    density_vector = numeric()
    
    for (i in seq_along(interval_starts)) {
      start = interval_starts[i]
      end = interval_ends[i]
      d = end - start  # Interval length
      q0 = start
      q1 = start + 0.2 * d
      q2 = start + 0.4 * d
      q3 = start + 0.6 * d
      q4 = start + 0.8 * d
      tfactor=pcauchy(12.7062)-pcauchy(-12.7062) #factor de la funcion de densidad truncada
      density_q0 = dgamma(q0,1,3)
      density_q1 = dgamma(q1,1,3)
      density_q2 = dgamma(q2,1,3)
      density_q3 = dgamma(q3,1,3)
      density_q4 = dgamma(q4,1,3)
      density_vector=c(density_vector, list(density_q0, density_q1, density_q2, density_q3, density_q4))
      
    }
    
    #Calculamos las diferentes normas
    
    #L1
    histogram=hist_data
    density_vector = unlist(density_vector)
    L1_vec=abs((density_vector - hist_density_values))*(histogram$breaks[2] - histogram$breaks[1])/5
    L1_norm=sum(L1_vec)
    #L INF
    linf_norm = max(L1_vec)
    #L2
    L2_vec=((density_vector - hist_density_values)^2)*(histogram$breaks[2] - histogram$breaks[1])/5
    L2_norm=sum(L2_vec)
    L2_norm=sqrt(L2_norm)
    
    #Formateamos los datos obtenidos
    l1_list=c(l1_list, L1_norm)
    l2_list=c(l2_list, L2_norm)
    linf_list=c(linf_list, linf_norm)
  }
  #Guardar y formatear normas obtenidas por iteracion de n
  l1_list=unlist(l1_list)
  l2_list=unlist(l2_list)
  linf_list=unlist(linf_list)
  ave1=sum(l1_list)/20
  ave2=sum(l2_list)/20
  aveInf=sum(linf_list)/20
  l1_ave=c(l1_ave, ave1)
  l2_ave=c(l2_ave, ave2)
  linf_ave=c(linf_ave, aveInf)
  
}
#Resultados GAMMA(1,3)
print("Normas promediadas l1, l2 y linf para respectivos n y 20 muestras iid por n : ")
print("L1: ")
print(unlist(l1_ave))
print("L2: ")
print(unlist(l2_ave))
print("Linf: ")
print(unlist(linf_ave))


#log-log Plot
#Promedios calculados para las diferentes normas
n_list=unlist(n_list)
l1_ave=unlist(l1_ave)
l2_ave=unlist(l2_ave)
linf_ave=unlist(linf_ave)

#Transformadas a logaritmos para L1, regresion lineal y graficar
log_n_list = log(n_list)
log_l1_ave = log(l1_ave)

lm_model_l1 = lm(log_l1_ave ~ log_n_list)

plot(log_n_list, log_l1_ave, 
     xlab = "log(n)", ylab = "log(L1 promedio)", 
     main = "Log-Log Plot para gamma(1,3) en L1")
abline(lm_model_l1, col = "red")

eq_l1 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l1)[1], digits = 3)) + .(format(coef(lm_model_l1)[2], digits = 3)) * italic(x)))
text(x = 7, y = -2.7, labels = eq_l1, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

#Mismo proceso para L2
log_l2_ave = log(l2_ave)
lm_model_l2 = lm(log_l2_ave ~ log_n_list)

plot(log_n_list, log_l2_ave, 
     xlab = "log(n)", ylab = "log(L2 promedio)", 
     main = "Log-Log Plot para gamma(1,3) en L2")
abline(lm_model_l2, col = "blue")

eq_l2 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l2)[1], digits = 3)) + .(format(coef(lm_model_l2)[2], digits = 3)) * italic(x)))
text(x = 7, y = -3.6, labels = eq_l2, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

# Mismo proceso para Linf
log_linf_ave = log(linf_ave)

lm_model_linf = lm(log_linf_ave ~ log_n_list)

plot(log_n_list, log_linf_ave, 
     xlab = "log(n)", ylab = "log(L_inf promedio)", 
     main = "Log-Log Plot para gamma(1,3) en Linf")

abline(lm_model_linf, col = "green")

eq_l3 = as.expression(bquote(italic(y) == .(format(coef(lm_model_linf)[1], digits = 3)) + .(format(coef(lm_model_linf)[2], digits = 3)) * italic(x)))
text(x = 7, y = -6.5, labels = eq_l3, pos = 4, cex = 0.8) #Agregar ecuacion al grafico


--------------------------------------------------------------
#Gamma(4,3)
#-------------------------------------------------------------


n_list =list(1000, 2000, 5000, 10000, 20000)
#Listas para guardar y promediar las diferentes normas en cada iteracion
l1_list =list()
l2_list= list()
linf_list = list()
l1_ave=list()
l2_ave=list()
linf_ave=list()

for (n in n_list) {
  l1_list=list()
  l2_list=list()
  linf_list=list()
  for (i in 1:20) {
    
    #Generar muestra para n fijo
    gammaData=rgamma(n,4,3)
    
    hist_data = hist(gammaData, freq=FALSE, breaks=0.25*sqrt(n), plot = TRUE)
    lines(density(gammaData))
    
    interval_starts = hist_data$breaks[-length(hist_data$breaks)]
    interval_ends = hist_data$breaks[-1]
    
    density_values=hist_data$density #vector de densidades de cada intervalo
    
    
    #Obtener comeinzo + 4 datos intermedios por intervalo de histograma y densidad
    
    delta=0.2
    #delta corresponde al peso de cada separacion del intervalo, como tomamos 5
    # elementos por intervalo cada uno tiene peso 1/5
    
    #VECTOR DEL HISTOGRAMA:
    #Dado que el histograma es constante en un intervalo para generar el vector
    #ordenado simplemente repetimos el valor de densidad 5 veces por intervalo
    hist_density_values = numeric()
    for (i in seq_along(density_values)) {
      hist_density_values = c(hist_density_values, rep(density_values[i], each = 5))
    }
    
    
    #VECTOR DE DENSIDAD REAL:
    
    #Calculamos explicitamente la densidad real para cada uno de los 5 puntos equidistantes del intervalo 
    #del histograma y formamos un vector ordenado para calcular la distancia de ambos vectores a continuacion
    
    density_vector = numeric()
    
    for (i in seq_along(interval_starts)) {
      start = interval_starts[i]
      end = interval_ends[i]
      d = end - start  # Interval length
      q0 = start
      q1 = start + 0.2 * d
      q2 = start + 0.4 * d
      q3 = start + 0.6 * d
      q4 = start + 0.8 * d
      tfactor=pcauchy(12.7062)-pcauchy(-12.7062) #factor de la funcion de densidad truncada
      density_q0 = dgamma(q0,4,3)
      density_q1 = dgamma(q1,4,3)
      density_q2 = dgamma(q2,4,3)
      density_q3 = dgamma(q3,4,3)
      density_q4 = dgamma(q4,4,3)
      density_vector=c(density_vector, list(density_q0, density_q1, density_q2, density_q3, density_q4))
      
    }
    
    #Calculamos las diferentes normas
    
    #L1
    histogram=hist_data
    density_vector = unlist(density_vector)
    L1_vec=abs((density_vector - hist_density_values))*(histogram$breaks[2] - histogram$breaks[1])/5
    L1_norm=sum(L1_vec)
    #L INF
    linf_norm = max(L1_vec)
    #L2
    L2_vec=((density_vector - hist_density_values)^2)*(histogram$breaks[2] - histogram$breaks[1])/5
    L2_norm=sum(L2_vec)
    L2_norm=sqrt(L2_norm)
    
    #Formateamos los datos obtenidos
    l1_list=c(l1_list, L1_norm)
    l2_list=c(l2_list, L2_norm)
    linf_list=c(linf_list, linf_norm)
  }
  #Guardar y formatear normas obtenidas por iteracion de n
  l1_list=unlist(l1_list)
  l2_list=unlist(l2_list)
  linf_list=unlist(linf_list)
  ave1=sum(l1_list)/20
  ave2=sum(l2_list)/20
  aveInf=sum(linf_list)/20
  l1_ave=c(l1_ave, ave1)
  l2_ave=c(l2_ave, ave2)
  linf_ave=c(linf_ave, aveInf)
  
}
#Resultados GAMMA(4,3)

print("Normas promediadas l1, l2 y linf para respectivos n y 20 muestras iid por n : ")
print("L1: ")
print(unlist(l1_ave))
print("L2: ")
print(unlist(l2_ave))
print("Linf: ")
print(unlist(linf_ave))

#log-log Plot
#Promedios calculados para las diferentes normas
n_list=unlist(n_list)
l1_ave=unlist(l1_ave)
l2_ave=unlist(l2_ave)
linf_ave=unlist(linf_ave)

#Transformadas a logaritmos para L1, regresion lineal y graficar
log_n_list = log(n_list)
log_l1_ave = log(l1_ave)

lm_model_l1 = lm(log_l1_ave ~ log_n_list)

plot(log_n_list, log_l1_ave, 
     xlab = "log(n)", ylab = "log(L1 promedio)", 
     main = "Log-Log Plot para gamma(4,3) en L1")
abline(lm_model_l1, col = "red")

eq_l1 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l1)[1], digits = 3)) + .(format(coef(lm_model_l1)[2], digits = 3)) * italic(x)))
text(x = 7, y = -2.7, labels = eq_l1, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

#Mismo proceso para L2
log_l2_ave = log(l2_ave)
lm_model_l2 = lm(log_l2_ave ~ log_n_list)

plot(log_n_list, log_l2_ave, 
     xlab = "log(n)", ylab = "log(L2 promedio)", 
     main = "Log-Log Plot para gamma(4,3) en L2")
abline(lm_model_l2, col = "blue")

eq_l2 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l2)[1], digits = 3)) + .(format(coef(lm_model_l2)[2], digits = 3)) * italic(x)))
text(x = 7, y = -2.8, labels = eq_l2, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

# Mismo proceso para Linf
log_linf_ave = log(linf_ave)

lm_model_linf = lm(log_linf_ave ~ log_n_list)

plot(log_n_list, log_linf_ave, 
     xlab = "log(n)", ylab = "log(L_inf promedio)", 
     main = "Log-Log Plot para gamma(4,3) en Linf")

abline(lm_model_linf, col = "green")

eq_l3 = as.expression(bquote(italic(y) == .(format(coef(lm_model_linf)[1], digits = 3)) + .(format(coef(lm_model_linf)[2], digits = 3)) * italic(x)))
text(x = 7, y = -8, labels = eq_l3, pos = 4, cex = 0.8) #Agregar ecuacion al grafico



#-----------------------------------------------------------------
#Unif(0,1)
#-----------------------------------------------------------------

n_list =list(1000, 2000, 5000, 10000, 20000)
#Listas para guardar y promediar las diferentes normas en cada iteracion
l1_list =list()
l2_list= list()
linf_list = list()
l1_ave=list()
l2_ave=list()
linf_ave=list()

for (n in n_list) {
  l1_list=list()
  l2_list=list()
  linf_list=list()
  for (i in 1:20) {
    
    #Generar muestra para n fijo
    unifData=runif(n,0,1)
    
    hist_data = hist(unifData, freq=FALSE, breaks=1/4*sqrt(n), plot = TRUE)
    lines(density(unifData))
    
    interval_starts = hist_data$breaks[-length(hist_data$breaks)]
    interval_ends = hist_data$breaks[-1]
    
    density_values=hist_data$density #vector de densidades de cada intervalo
    
    
    #Obtener comeinzo + 4 datos intermedios por intervalo de histograma y densidad
    
    delta=0.2
    #delta corresponde al peso de cada separacion del intervalo, como tomamos 5
    # elementos por intervalo cada uno tiene peso 1/5
    
    #VECTOR DEL HISTOGRAMA:
    #Dado que el histograma es constante en un intervalo para generar el vector
    #ordenado simplemente repetimos el valor de densidad 5 veces por intervalo
    hist_density_values = numeric()
    for (i in seq_along(density_values)) {
      hist_density_values = c(hist_density_values, rep(density_values[i], each = 5))
    }
    
    
    #VECTOR DE DENSIDAD REAL:
    
    #Calculamos explicitamente la densidad real para cada uno de los 5 puntos equidistantes del intervalo 
    #del histograma y formamos un vector ordenado para calcular la distancia de ambos vectores a continuacion
    
    density_vector = numeric()
    
    for (i in seq_along(interval_starts)) {
      start = interval_starts[i]
      end = interval_ends[i]
      d = end - start  # Interval length
      q0 = start
      q1 = start + 0.2 * d
      q2 = start + 0.4 * d
      q3 = start + 0.6 * d
      q4 = start + 0.8 * d
      tfactor=pcauchy(12.7062)-pcauchy(-12.7062) #factor de la funcion de densidad truncada
      density_q0 = dunif(q0,0,1)
      density_q1 = dunif(q1,0,1)
      density_q2 = dunif(q2,0,1)
      density_q3 = dunif(q3,0,1)
      density_q4 = dunif(q4,0,1)
      density_vector=c(density_vector, list(density_q0, density_q1, density_q2, density_q3, density_q4))
      
    }
    
    #Calculamos las diferentes normas
    
    #L1
    histogram=hist_data
    density_vector = unlist(density_vector)
    L1_vec=abs((density_vector - hist_density_values))*(histogram$breaks[2] - histogram$breaks[1])/5
    L1_norm=sum(L1_vec)
    #L INF
    linf_norm = max(L1_vec)
    #L2
    L2_vec=((density_vector - hist_density_values)^2)*(histogram$breaks[2] - histogram$breaks[1])/5
    L2_norm=sum(L2_vec)
    L2_norm=sqrt(L2_norm)
    
    #Formateamos los datos obtenidos
    l1_list=c(l1_list, L1_norm)
    l2_list=c(l2_list, L2_norm)
    linf_list=c(linf_list, linf_norm)
  }
  #Guardar y formatear normas obtenidas por iteracion de n
  l1_list=unlist(l1_list)
  l2_list=unlist(l2_list)
  linf_list=unlist(linf_list)
  ave1=sum(l1_list)/20
  ave2=sum(l2_list)/20
  aveInf=sum(linf_list)/20
  l1_ave=c(l1_ave, ave1)
  l2_ave=c(l2_ave, ave2)
  linf_ave=c(linf_ave, aveInf)
  
}

#Resultados UNIF(0,1)
print("Normas promediadas l1, l2 y linf para respectivos n y 20 muestras iid por n : ")
print("L1: ")
print(unlist(l1_ave))
print("L2: ")
print(unlist(l2_ave))
print("Linf: ")
print(unlist(linf_ave))

#log-log Plot
#Promedios calculados para las diferentes normas
n_list=unlist(n_list)
l1_ave=unlist(l1_ave)
l2_ave=unlist(l2_ave)
linf_ave=unlist(linf_ave)

#Transformadas a logaritmos para L1, regresion lineal y graficar
log_n_list = log(n_list)
log_l1_ave = log(l1_ave)

lm_model_l1 = lm(log_l1_ave ~ log_n_list)

plot(log_n_list, log_l1_ave, 
     xlab = "log(n)", ylab = "log(L1 promedio)", 
     main = "Log-Log Plot para unif(0,1) en L1")
abline(lm_model_l1, col = "red")

eq_l1 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l1)[1], digits = 3)) + .(format(coef(lm_model_l1)[2], digits = 3)) * italic(x)))
text(x = 7, y = -2.7, labels = eq_l1, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

#Mismo proceso para L2
log_l2_ave = log(l2_ave)
lm_model_l2 = lm(log_l2_ave ~ log_n_list)

plot(log_n_list, log_l2_ave, 
     xlab = "log(n)", ylab = "log(L2 promedio)", 
     main = "Log-Log Plot para unif(0,1) en L2")
abline(lm_model_l2, col = "blue")

eq_l2 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l2)[1], digits = 3)) + .(format(coef(lm_model_l2)[2], digits = 3)) * italic(x)))
text(x = 7, y = -3.6, labels = eq_l2, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

# Mismo proceso para Linf
log_linf_ave = log(linf_ave)

lm_model_linf = lm(log_linf_ave ~ log_n_list)

plot(log_n_list, log_linf_ave, 
     xlab = "log(n)", ylab = "log(L_inf promedio)", 
     main = "Log-Log Plot para unif(0,1) en Linf")

abline(lm_model_linf, col = "green")

eq_l3 = as.expression(bquote(italic(y) == .(format(coef(lm_model_linf)[1], digits = 3)) + .(format(coef(lm_model_linf)[2], digits = 3)) * italic(x)))
text(x = 7, y = -6.5, labels = eq_l3, pos = 4, cex = 0.8) #Agregar ecuacion al grafico




#--------------------------------------------------------
#Rayleigh(2)
#--------------------------------------------------------
theta=2

rayleigh_density = function(x, theta) {(x / theta^2) * exp(-x^2 / (2 * theta^2))
}

FR=function(x,theta){1 - exp(-x^2 / (2 * theta^2))}


n_list =list(1000, 2000, 5000, 10000, 20000)
#Listas para guardar y promediar las diferentes normas en cada iteracion
l1_list =list()
l2_list= list()
linf_list = list()
l1_ave=list()
l2_ave=list()
linf_ave=list()

for (n in n_list) {
  l1_list=list()
  l2_list=list()
  linf_list=list()
  for (i in 1:20) {
    
    #Generar muestra para n fijo
    U=runif(n,min=0,max=1)
    
    #Calculamos inversa explicitamente
    find_root = function(u,theta) {
      uniroot(function(x) FR(x,theta) - u, c(0, 100))$root
    }
    
    rayleighData=sapply(U, find_root, theta = theta)
    
    hist_data = hist(rayleighData, freq=FALSE, breaks=1/4*sqrt(n), plot = FALSE)
    #lines(density(rayleighData))
    
    interval_starts = hist_data$breaks[-length(hist_data$breaks)]
    interval_ends = hist_data$breaks[-1]
    
    density_values=hist_data$density #vector de densidades de cada intervalo
    
    
    #Obtener comeinzo + 4 datos intermedios por intervalo de histograma y densidad
    
    delta=0.2
    #delta corresponde al peso de cada separacion del intervalo, como tomamos 5
    # elementos por intervalo cada uno tiene peso 1/5
    
    #VECTOR DEL HISTOGRAMA:
    #Dado que el histograma es constante en un intervalo para generar el vector
    #ordenado simplemente repetimos el valor de densidad 5 veces por intervalo
    hist_density_values = numeric()
    for (i in seq_along(density_values)) {
      hist_density_values = c(hist_density_values, rep(density_values[i], each = 5))
    }
    
    
    #VECTOR DE DENSIDAD REAL:
    
    #Calculamos explicitamente la densidad real para cada uno de los 5 puntos equidistantes del intervalo 
    #del histograma y formamos un vector ordenado para calcular la distancia de ambos vectores a continuacion
    
    density_vector = numeric()
    
    for (i in seq_along(interval_starts)) {
      start = interval_starts[i]
      end = interval_ends[i]
      d = end - start  # Interval length
      q0 = start
      q1 = start + 0.2 * d
      q2 = start + 0.4 * d
      q3 = start + 0.6 * d
      q4 = start + 0.8 * d
      tfactor=pcauchy(12.7062)-pcauchy(-12.7062) #factor de la funcion de densidad truncada
      density_q0 = rayleigh_density(q0,2)
      density_q1 = rayleigh_density(q1,2)
      density_q2 = rayleigh_density(q2,2)
      density_q3 = rayleigh_density(q3,2)
      density_q4 = rayleigh_density(q4,2)
      density_vector=c(density_vector, list(density_q0, density_q1, density_q2, density_q3, density_q4))
      
    }
    
    #Calculamos las diferentes normas
    
    #L1
    histogram=hist_data
    density_vector = unlist(density_vector)
    L1_vec=abs((density_vector - hist_density_values))*(histogram$breaks[2] - histogram$breaks[1])/5
    L1_norm=sum(L1_vec)
    #L INF
    linf_norm = max(L1_vec)
    #L2
    L2_vec=((density_vector - hist_density_values)^2)*(histogram$breaks[2] - histogram$breaks[1])/5
    L2_norm=sum(L2_vec)
    L2_norm=sqrt(L2_norm)
    
    #Formateamos los datos obtenidos
    l1_list=c(l1_list, L1_norm)
    l2_list=c(l2_list, L2_norm)
    linf_list=c(linf_list, linf_norm)
  }
  #Guardar y formatear normas obtenidas por iteracion de n
  l1_list=unlist(l1_list)
  l2_list=unlist(l2_list)
  linf_list=unlist(linf_list)
  ave1=sum(l1_list)/20
  ave2=sum(l2_list)/20
  aveInf=sum(linf_list)/20
  l1_ave=c(l1_ave, ave1)
  l2_ave=c(l2_ave, ave2)
  linf_ave=c(linf_ave, aveInf)
  
}

#Resultados RAYLEIGH(0,1)
print("Normas promediadas l1, l2 y linf para respectivos n y 20 muestras iid por n : ")
print("L1: ")
print(unlist(l1_ave))
print("L2: ")
print(unlist(l2_ave))
print("Linf: ")
print(unlist(linf_ave))

#log-log Plot
#Promedios calculados para las diferentes normas
n_list=unlist(n_list)
l1_ave=unlist(l1_ave)
l2_ave=unlist(l2_ave)
linf_ave=unlist(linf_ave)

#Transformadas a logaritmos para L1, regresion lineal y graficar
log_n_list = log(n_list)
log_l1_ave = log(l1_ave)

lm_model_l1 = lm(log_l1_ave ~ log_n_list)

plot(log_n_list, log_l1_ave, 
     xlab = "log(n)", ylab = "log(L1 promedio)", 
     main = "Log-Log Plot para Rayleigh(2) en L1")
abline(lm_model_l1, col = "red")

eq_l1 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l1)[1], digits = 3)) + .(format(coef(lm_model_l1)[2], digits = 3)) * italic(x)))
text(x = 7, y = -2.7, labels = eq_l1, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

#Mismo proceso para L2
log_l2_ave = log(l2_ave)
lm_model_l2 = lm(log_l2_ave ~ log_n_list)

plot(log_n_list, log_l2_ave, 
     xlab = "log(n)", ylab = "log(L2 promedio)", 
     main = "Log-Log Plot para Rayleigh(2) en L2")
abline(lm_model_l2, col = "blue")

eq_l2 = as.expression(bquote(italic(y) == .(format(coef(lm_model_l2)[1], digits = 3)) + .(format(coef(lm_model_l2)[2], digits = 3)) * italic(x)))
text(x = 7, y = -3.6, labels = eq_l2, pos = 4, cex = 0.8) #Agregar ecuacion al grafico

# Mismo proceso para Linf
log_linf_ave = log(linf_ave)

lm_model_linf = lm(log_linf_ave ~ log_n_list)

plot(log_n_list, log_linf_ave, 
     xlab = "log(n)", ylab = "log(L_inf promedio)", 
     main = "Log-Log Plot para Rayleigh(2) en Linf")

abline(lm_model_linf, col = "green")

eq_l3 = as.expression(bquote(italic(y) == .(format(coef(lm_model_linf)[1], digits = 3)) + .(format(coef(lm_model_linf)[2], digits = 3)) * italic(x)))
text(x = 7, y = -6.5, labels = eq_l3, pos = 4, cex = 0.8) #Agregar ecuacion al grafico














