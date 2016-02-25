library (rgdal)
library (alphahull)
library (sp)
library (maptools)
library (raster)
library (adehabitatHR)
library (rgeos)
library (tlocoh)
library(np)
library(fields)
library(stringr)
library (gridExtra)

## con estas funciones defines las formas, la x y la y de cada forma, y el área total=12

circulo<- circle_to_shape (2, 2, 12)
circulos_iguales<- circles_to_shape (c(0,4,4), c(2,0,4),c(4,4,4))
circulos_dif<- circles_to_shape (c(0,4,4), c(2,0,4),c(8,3,1))

rosco<- rosco_to_shape (2,2,12, 5)
roscos_iguales<- roscos_to_shape (c(0,4,4), c(2,0,4),c(4, 4, 4), c(2, 2, 2))
roscos_dif<- roscos_to_shape (c(0,4,4), c(2,0,4),c(8,3,1), c(2, 1, 0.5))

star<- estrella_to_shape (2,2, area=12, puntas= 5, relacion_radios=3)
stars_iguales<- estrellas_to_shape (c(0,4,4), c(2,0,4), area=c(4, 4, 4), puntas= 5, relacion_radios=3)
stars_dif<- estrellas_to_shape (c(0,4,4), c(2,0,4), area=c(8, 3, 1), puntas= 5, relacion_radios=3)

#' nsample = vector de 3 elementos, e.g. c(50, 100, 250). Si no pones nada esos serán los tamaños de muestra por defecto
#' n sample minima = 50
#' errores = numero entero. porcentaje de errores que queremos en nuestra muestra, 
#' % puntos al azar lanzados fuera del polígono de la especie. Por defecto errores=0
#' perc = vector de 5 elementos. porcetnajes de puntos que usa el convex hull. por defecto c(100, 95, 90, 85, 80) 
#' k= valores para k en el locohk. 
#' r= valores para r en el locohr
#' a= valores para a en el locoha
#' thresh= valores para cortar el mapa continuo que da el kernel. están definidos como percentiles de la prob predicha.
#' 1= el mínimo valor de predicción. sugerencia c(1, 5, 10, 15, 20)

experimentos<- function (obj= list(circulo, circulos_iguales, circulos_dif, rosco, roscos_iguales, roscos_dif, 
                                   star, stars_iguales, stars_dif), 
                         nsample = c(50, 100, 250), errores=0, 
                         perc=c(100, 95, 90, 85, 80), 
                         k=c(10, 15, 20, 25, 30), 
                         r=c(2, 2.2, 2.4, 2.8, 3),
                         a= c(5, 5.5, 6, 6.5, 7), 
                         thresh= c(1, 5, 10, 15, 20)){
  
  variables<- c("Shape", "N_frag", "Size_frag", 
                "N_sample", "Sampling_method", "Algorithm", "parameter", 
                "correct_negatives", "Correct_positives", "False_positive", 
                "False_negative", "Predicted_area", "Real_area", "Sensitivity", "Specificity", "TSS", "N_frag_pred")
  
  obj<- obj
  b<- c("circulo", "circulos_iguales", "circulos_dif", "rosco", "roscos_iguales", 
        "roscos_dif", "star", "stars_iguales", "stars_dif")
  
  modelos<- 5
  formas<- length (b)
  parametros_modelo<- 5
  sample_size<- 3
  sample_method<- 3
  
  total<- sample_size*sample_method*modelos*parametros_modelo
  filas<- total * formas
  n<- total- 1
  i<- 0
  
  res<- matrix (0,  total*formas, length (variables))
  colnames (res)<- variables
  
  for (j in 1:length (obj)){
    shape<- obj[[j]]
    nombre<- b[j]
    i<- i+1
    ## forma
    res [i:(i+n), 1]<- nombre
    ## numero de patches y tamaño relativo
    patches<- str_locate_all(pattern ='_',nombre)
    if (nrow (patches[[1]])==0){
      res [i:(i+n), 2]<-1 
      res [i:(i+n), 3]<-NA
    }
    if (nrow (patches[[1]])==1){
      res [i:(i+n), 2]<-3 
      patches<- str_locate_all(pattern ='iguales',nombre) 
      if (nrow (patches[[1]])==1){
        res [i:(i+n), 3]<- "equal" 
      } else {
        res [i:(i+n), 3]<- "heterog"
      }
    }
    
    ##muestreos
    if (errores==0){
    samples50<- muestreos (shape,n=nsample [1])
    samples100<- muestreos (shape, n= nsample [2])
    samples250<- muestreos (shape, n= nsample [3])
    } else {
    samples50<- muestreos_error (shape,n=nsample [1], p_error= errores)
    samples100<- muestreos_error (shape,n=nsample [2], p_error= errores)
    samples250<- muestreos_error (shape,n=nsample [3], p_error= errores)  
    }
    
    res[i:(i+n), 4]<- c(rep(nsample [1], total/3), rep(nsample [2], total/3), rep(nsample [3], total/3))
    res[i:(i+n), 5]<- rep (c(rep ("random", total/9), rep( "unif", total/9), 
                             rep("clustered", total/9)), 3)
    res[i:(i+n), 6]<- rep( c(rep("mcp", parametros_modelo), rep("locohk", parametros_modelo), 
                          rep( "locohr", parametros_modelo), rep( "locoha", parametros_modelo), 
                          rep( "kern", parametros_modelo)), 9) 
    
    ## algoritmos para predecir el área
    t<- i
    res[t:(t+n), 7]<- rep(c(paste ("perc=", perc [1], sep=""),
                         paste ("perc=", perc [2], sep=""),
                         paste ("perc=", perc [3], sep=""),
                         paste ("perc=", perc [4], sep=""),
                         paste ("perc=", perc [5], sep=""),
                         paste ("k=", k [1], sep=""),
                         paste ("k=", k [2], sep=""),
                         paste ("k=", k [3], sep=""),
                         paste ("k=", k [4], sep=""),
                         paste ("k=", k [5], sep=""),
                         paste ("r=", r [1], sep=""),
                         paste ("r=", r [2], sep=""),
                         paste ("r=", r [3], sep=""),
                         paste ("r=", r [4], sep=""),
                         paste ("r=", r [5], sep=""),
                         paste ("a=", a [1], sep=""),
                         paste ("a=", a [2], sep=""),
                         paste ("a=", a [3], sep=""),
                         paste ("a=", a [4], sep=""),
                         paste ("a=", a [5], sep=""),
                         paste ("thresh=", thresh [1], sep=""),
                         paste ("thresh=", thresh [2], sep=""),
                         paste ("thresh=", thresh [3], sep=""),
                         paste ("thresh=", thresh [4], sep=""),
                         paste ("thresh=", thresh [5], sep="")), 9)
                       
    puntos<- samples50$rand 
    pred_rand_50<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    res[t:(t+24), 17]<- n_frag_pred2 (pred_rand_50)
    
    puntos<- samples50$reg 
    pred_reg_50<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    t<- t+ (modelos* parametros_modelo)
    res[t:(t+24), 17]<- n_frag_pred2(pred_reg_50)
    
    puntos<- samples50$clus
    pred_clus_50<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    t<- t+ (modelos* parametros_modelo)
    res[t:(t+24), 17]<- n_frag_pred2(pred_clus_50)
    
    puntos<- samples100$rand 
    pred_rand_100<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    t<- t+ (modelos* parametros_modelo)
    res[t:(t+24), 17]<- n_frag_pred2(pred_rand_100)
    
    puntos<- samples100$reg 
    pred_reg_100<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    t<- t+ (modelos* parametros_modelo)
    res[t:(t+24), 17]<- n_frag_pred2(pred_reg_100)
    
    puntos<- samples100$clus
    pred_clus_100<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    t<- t+ (modelos* parametros_modelo)
    res[t:(t+24), 17]<- n_frag_pred2(pred_clus_100)
    
    puntos<- samples250$rand 
    pred_rand_250<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    t<- t+ (modelos* parametros_modelo)
    res[t:(t+24), 17]<- n_frag_pred2(pred_rand_250)
    
    puntos<- samples250$reg 
    pred_reg_250<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    t<- t+ (modelos* parametros_modelo)
    res[t:(t+24), 17]<- n_frag_pred2(pred_reg_250)
    
    puntos<- samples250$clus
    pred_clus_250<- shapefiles_pred2 (puntos, perc, k, r, a, thresh)
    t<- t+ (modelos* parametros_modelo)
    res[t:(t+24), 17]<- n_frag_pred2(pred_clus_250)
    
    c_real<- calc_area (shape, 100, 100)
    
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_rand_50 [[l]])){
        re<- NA
      }else{
      pred<- calc_area (pred_rand_50 [[l]], 100, 100)
      re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
 
    i<- i+ modelos* parametros_modelo
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_reg_50 [[l]])){
        re<- NA
      }else{
        pred<- calc_area (pred_reg_50 [[l]], 100, 100)
        re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
    
    i<- i+ modelos* parametros_modelo
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_clus_50 [[l]])){
        re<- NA
      }else{
        pred<- calc_area (pred_clus_50 [[l]], 100, 100)
        re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
    
    i<- i+ modelos* parametros_modelo
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_rand_100 [[l]])){
        re<- NA
      }else{
        pred<- calc_area (pred_rand_100 [[l]], 100, 100)
        re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
    
    i<- i+ modelos* parametros_modelo
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_reg_100 [[l]])){
        re<- NA
      }else{
        pred<- calc_area (pred_reg_100 [[l]], 100, 100)
        re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
    
    i<- i+ modelos* parametros_modelo
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_clus_100 [[l]])){
        re<- NA
      }else{
        pred<- calc_area (pred_clus_100 [[l]], 100, 100)
        re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
    
    i<- i+ modelos* parametros_modelo
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_rand_250 [[l]])){
        re<- NA
      }else{
        pred<- calc_area (pred_rand_250 [[l]], 100, 100)
        re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
    
    i<- i+ modelos* parametros_modelo
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_reg_250 [[l]])){
        re<- NA
      }else{
        pred<- calc_area (pred_reg_250 [[l]], 100, 100)
        re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
    
    i<- i+ modelos* parametros_modelo
    for (l in 1:(modelos* parametros_modelo)){
      if (is.na (pred_clus_250 [[l]])){
        re<- NA
      }else{
        pred<- calc_area (pred_clus_250 [[l]], 100, 100)
        re<- accuracy (c_real, pred)
      }
      res[i+(l-1), 8:11]<- re
    }
    
    res[, 13]<- as.numeric (res[, 9]) + as.numeric (res[, 11])
    res[, 12]<- as.numeric (res[, 9]) + as.numeric (res[, 10])
    res[, 14]<- round (as.numeric (res[, 9]) / as.numeric (res[, 13]), 2)
    res[, 15]<- round (as.numeric (res[, 8]) / (as.numeric (res[, 8])+  as.numeric (res[, 10])), 2)
    res[, 16]<- round (as.numeric (res [, 14]) + as.numeric (res [, 15]) -1, 2)
    
    i<- i+ (modelos* parametros_modelo -1)
  }
  
  write.table (res, paste ("resultados_", errores, ".csv", sep=""), sep=",", row.names = F)
  res
}


