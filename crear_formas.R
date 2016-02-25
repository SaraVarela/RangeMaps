library (maptools)
library (raster)
library (rgeos)
library (gridExtra)

## createCircle function modified from Roger Bivand 
createCircle <- function(x,y,r,start=0,
                         end=2*pi,
                         nsteps=100000,...){
  rs <- seq(start,end,len=nsteps)
  xc <- x+r*cos(rs)
  yc <- y+r*sin(rs)
  my.pol<-cbind(xc,yc)
  my.pol <- rbind(my.pol, my.pol[1,])
  my.pol
}

## pasar el circulo a un shapefile en proyección latlong
circle_to_shape<- function (x,y,area){
r<- sqrt (area/pi)
circle<-createCircle(x,y,r)
circle_pol <- Polygons(list(Polygon(circle)), ID="1")
circle_shape<- SpatialPolygons(list(circle_pol))
crs (circle_shape)<- CRS("+proj=longlat +datum=WGS84") 
circle_shape
}

## en este caso x, y, area son vectores
circles_to_shape<- function (x,y,area){
  r<- sqrt (area[1]/pi)
  circle<-createCircle(x[1],y[1],r)
  circle_pol1 <- Polygons(list(Polygon(circle)), ID="1")
  r<- sqrt (area[2]/pi)
  circle<-createCircle(x[2],y[2],r)
  circle_pol2 <- Polygons(list(Polygon(circle)), ID="2")
  r<- sqrt (area[3]/pi)
  circle<-createCircle(x[3],y[3],r)
  circle_pol3 <- Polygons(list(Polygon(circle)), ID="3")
  
  circle_shape<- SpatialPolygons(list(circle_pol1, circle_pol2, circle_pol3))
  crs (circle_shape)<- CRS("+proj=longlat +datum=WGS84") 
  circle_shape
}

## para un rosco
rosco_to_shape<- function (x,y,area_final, area2){
  area1<- area_final *2 + area2
  r<- sqrt (area1/pi)
  r2<- sqrt (area2/pi)
  circle<-createCircle(x,y,r)
  circle2<-createCircle(x,y,r2)
  
  circle_pol <- Polygons(list(Polygon(circle),  
                              Polygon(circle2, hole=TRUE)), ID="1")
  
  circle_shape<- SpatialPolygons(list(circle_pol))
  crs (circle_shape)<- CRS("+proj=longlat +datum=WGS84") 
 
  corte_x<- round ((extent (circle_shape)[1]+ extent (circle_shape)[2])/2, 2)
  b<- as(extent(corte_x, extent (circle_shape)[2],
         extent (circle_shape)[3], extent (circle_shape)[4]), 'SpatialPolygons')
  crs(b) <- crs(circle_shape)
  rosco<- crop(circle_shape, b)
  rosco
}

## para 3 roscos, 4 vectores de 3 elementos
roscos_to_shape<- function (x,y,area_final, area2){
  area1<- area_final[1]  *2 + area2[1]
  r<- sqrt (area1/pi)
  r2<- sqrt (area2[1]/pi)
  circle<-createCircle(x[1],y[1],r)
  circle2<-createCircle(x[1],y[1],r2)
  circle_pol1 <- Polygons(list(Polygon(circle),  
                              Polygon(circle2, hole=TRUE)), ID="1")
  
  circle_shape1<- SpatialPolygons(list(circle_pol1))
  crs (circle_shape1)<- CRS("+proj=longlat +datum=WGS84") 
  
  corte_x<- round ((extent (circle_shape1)[1]+ extent (circle_shape1)[2])/2, 2)
  b<- as(extent(corte_x, extent (circle_shape1)[2],
                extent (circle_shape1)[3], extent (circle_shape1)[4]), 'SpatialPolygons')
  crs(b) <- crs(circle_shape1)
  rosco1<- crop(circle_shape1, b)
 
  area1<- area_final[2]  *2 + area2[2]
  r<- sqrt (area1/pi)
  r2<- sqrt (area2[2]/pi)
  circle<-createCircle(x[2],y[2],r)
  circle2<-createCircle(x[2],y[2],r2)
  circle_pol2 <- Polygons(list(Polygon(circle),  
                               Polygon(circle2, hole=TRUE)), ID="2")
  circle_shape2<- SpatialPolygons(list(circle_pol2))
  crs (circle_shape2)<- CRS("+proj=longlat +datum=WGS84") 
  
  corte_x<- round ((extent (circle_shape2)[1]+ extent (circle_shape2)[2])/2, 2)
  b<- as(extent(corte_x, extent (circle_shape2)[2],
                extent (circle_shape2)[3], extent (circle_shape2)[4]), 'SpatialPolygons')
  crs(b) <- crs(circle_shape2)
  rosco2<- crop(circle_shape2, b)
  
  area1<- area_final[3]  *2 + area2[3]
  r<- sqrt (area1/pi)
  r2<- sqrt (area2[3]/pi)
  circle<-createCircle(x[3],y[3],r)
  circle2<-createCircle(x[3],y[3],r2)
  circle_pol3 <- Polygons(list(Polygon(circle),  
                               Polygon(circle2, hole=TRUE)), ID="3")
  circle_shape3<- SpatialPolygons(list(circle_pol3))
  crs (circle_shape3)<- CRS("+proj=longlat +datum=WGS84") 
  
  corte_x<- round ((extent (circle_shape3)[1]+ extent (circle_shape3)[2])/2, 2)
  b<- as(extent(corte_x, extent (circle_shape3)[2],
                extent (circle_shape3)[3], extent (circle_shape3)[4]), 'SpatialPolygons')
  crs(b) <- crs(circle_shape3)
  rosco3<- crop(circle_shape3, b)
  
  roscos<- gUnion(rosco1, rosco2, byid=TRUE) 
  roscos<- gUnion(roscos, rosco3, byid=TRUE) 
  roscos
}

## para crear la elispe
createEllipse <- function(x,y, a, b,
                         nsteps=10000,...){
  xc <- seq(-a,a,len=nsteps)
  yc <- c(sqrt(1- (xc^2/a^2)) * b, -sqrt(1- (xc^2/a^2)) * b)
  my.pol<-cbind(c(xc, rev (xc)), rev (yc))
  my.pol <- rbind(my.pol, my.pol[1,])
  my.pol<- cbind (my.pol [,1]+ x, my.pol [,2]+ y)
}


## ratio definirá la relación entre los radios, la forma de la elipse
ellipse_to_shape<- function (x,y,ratio, area){
  a<- sqrt (area/ratio*pi)
  b<- a/ratio
  elipse<- createEllipse(x,y, a, b)
  elipse_pol <- Polygons(list(Polygon(elipse)), ID="1")
  elipse_shape<- SpatialPolygons(list(elipse_pol))
  crs (elipse_shape)<- CRS("+proj=longlat +datum=WGS84") 
  elipse_shape
}

## en este caso x, y, area son vectores
ellipses_to_shape<- function (x,y, ratio, area){
  a<- sqrt (area[1]/ratio*pi)
  b<- a/ratio
  elipse<- createEllipse(x[1],y[1], a, b)
  elipse_pol1 <- Polygons(list(Polygon(elipse)), ID="1")
  
  a<- sqrt (area[2]/ratio*pi)
  b<- a/ratio
  elipse<- createEllipse(x[2],y[2], a, b)
  elipse_pol2 <- Polygons(list(Polygon(elipse)), ID="2")
  
  a<- sqrt (area[3]/ratio*pi)
  b<- a/ratio
  elipse<- createEllipse(x[3],y[3], a, b)
  elipse_pol3 <- Polygons(list(Polygon(elipse)), ID="3")
  
  elipse_shape<- SpatialPolygons(list(elipse_pol1, 
                                      elipse_pol2, elipse_pol3))
  crs (elipse_shape)<- CRS("+proj=longlat +datum=WGS84") 
  elipse_shape
}


## estrella
estrella<- function (x, y, area, relacion_radios, puntas){
  
  r_interior<- sqrt (area /(relacion_radios * puntas * sin (pi/puntas)))
  r_exterior<- relacion_radios * r_interior
  
  rs <- seq(0, 2*pi, len= (2* puntas+1))
  rs<- rs [-length (rs)]
  rs2<- rs [c(FALSE, TRUE)]
  rs3<- rs [c(TRUE, FALSE)]
  
  xc <- x+r_exterior*cos(rs3)
  yc <- y+r_exterior*sin(rs3)
  
  xxc <- x+r_interior*cos(rs2)
  yyc <- y+r_interior*sin(rs2)
  
  xc<- c(rbind (xc, xxc))
  yc<- c(rbind (yc, yyc))
  my.pol<-cbind(xc,yc)
  my.pol <- rbind(my.pol, my.pol[1,])
  my.pol
}


## pasar la estrella a un shapefile en proyección latlong
estrella_to_shape<- function (x,y, area, relacion_radios, puntas){
  circle<-estrella (x,y,area, relacion_radios, puntas)
  circle_pol <- Polygons(list(Polygon(circle)), ID="1")
  circle_shape<- SpatialPolygons(list(circle_pol))
  crs (circle_shape)<- CRS("+proj=longlat +datum=WGS84") 
  circle_shape
}

## en este caso x, y, area son vectores
estrellas_to_shape<- function (x,y,area, relacion_radios, puntas){
  circle<-estrella(x[1],y[1],area [1], relacion_radios, puntas)
  circle_pol1 <- Polygons(list(Polygon(circle)), ID="1")
  circle<-estrella(x[2],y[2],area [2], relacion_radios, puntas)
  circle_pol2 <- Polygons(list(Polygon(circle)), ID="2")
  circle<-estrella(x[3],y[3],area [3], relacion_radios, puntas)
  circle_pol3 <- Polygons(list(Polygon(circle)), ID="3")
  circle_shape<- SpatialPolygons(list(circle_pol1, circle_pol2, circle_pol3))
  crs (circle_shape)<- CRS("+proj=longlat +datum=WGS84") 
  circle_shape
}



## muestreos
muestreos<- function (shape, n){
  n_rand<- spsample (shape, n= n, type="random")
  n_reg<- spsample (shape, n= n, type="regular")
  n_clus<- spsample (shape, n= n, type="clustered", nclusters= n/2)
  list (rand=n_rand, reg=n_reg, clus=n_clus)
}


## funcion para usar npudens (que calcula la densidad de puntos)
## hacer una interpolación espacial con un thin plate spline
## y cortar con el threshold de minima presencia

kernel<- function (coord){
  #para mejorar el tiempo de computación
  options(np.tree=TRUE)
  f <- npudens(~x+y,ckertype="epanechnikov",data=coord)
  #### Thin plate spline model
  tps <- Tps(f$eval, f$dens)
  r <- raster(ncol=100, nrow=100)
  e<- extent (min (f$eval[,1]), max (f$eval[,1]), min (f$eval[,2]), max (f$eval[,2]))
  extent (r)<- e
  # use model to predict values at all locations
  p <- interpolate(r, tps)
  ## threshold con 0 omission, no nos dejamos ningún punto fuera
  thresh<- min (f$dens)
  map<- reclassify (p, c(-Inf, thres, 0, thres, +Inf, 1))
  pol <- rasterToPolygons(map, fun=function(x){x==1}, dissolve=TRUE)
  pol
}

## para poner diferentes thresholds al kernel

kernel2<- function (coord, thresh){
  #para mejorar el tiempo de computación
  options(np.tree=TRUE)
  f <- npudens(~x+y,ckertype="epanechnikov",data=coord)
  #### Thin plate spline model
  tps <- Tps(f$eval, f$dens)
  r <- raster(ncol=100, nrow=100)
  e<- extent (min (f$eval[,1]), max (f$eval[,1]), min (f$eval[,2]), max (f$eval[,2]))
  extent (r)<- e
  # use model to predict values at all locations
  p <- interpolate(r, tps)
  ## threshold con 0 omission, no nos dejamos ningún punto fuera
  densi<- sort (f$dens)
  t<- round ((length (densi) * thresh)/100, 0)
  t[1]<- 1
  thres<- densi [t]
  res_k<- list ()
  for (i in 1:5){
  map<- reclassify (p, c(-Inf, thres [i], 0, thres[i], +Inf, 1))
  pol <- rasterToPolygons(map, fun=function(x){x==1}, dissolve=TRUE)
  res_k [[i]]<- pol
  }
  res_k
}

## para que los modelos converjan en un solo parámetro.

shapefiles_pred<- function (puntos){
  
  mcp<- mcp(puntos, percent=100)
  
  k<- 10
  while (class (try (LoCoH.k(puntos, k=k, 
                             unin = "m", 
                             unout = "m2",
                             duplicates="remove"), silent=T))== "try-error"){ 
    k<- k+1 
  }
  locohk<- LoCoH.k(puntos, k=k, 
                   unin = "m", 
                   unout = "m2",
                   duplicates="remove")
  
  
  r<- 2
  while (class (try (LoCoH.r(puntos, r=r, 
                             unin = "m", 
                             unout = "m2",
                             duplicates="remove"), silent=T))== "try-error"){ 
    r<- r+0.1 
  }
  locohr<-LoCoH.r(puntos, r=r, 
                  unin = "m", 
                  unout = "m2",
                  duplicates="remove")
  
  a<- 5
  while (class (try (LoCoH.a(puntos, a=a, 
                             unin = "m", 
                             unout = "m2",
                             duplicates="remove"), silent=T))== "try-error"){ 
    a<- a+0.5 
  }
  locoha<-LoCoH.a(puntos, a=a, 
                  unin = "m", 
                  unout = "m2",
                  duplicates="remove")
  
  coord<- as.data.frame(puntos)
  names (coord)<- c("x", "y")
  kern<- kernel (coord)
 
  list (k=k, r=r, a=a, mcp=mcp,locohk=locohk,locohr=locohr, locoha=locoha, kern=kern)
}


## funcion para hacer los shapefiles con las predicciones. para 5 parámetros por modelo

shapefiles_pred2<- function (puntos, perc, k, r, a, thresh){

  resul<- list ()
  for (i in 1:5){
  resul [[i]]<- mcp(puntos, percent=perc[i])
  }
  
  
  for (i in 1:5){
    if (class (try (LoCoH.k(puntos, k=k[i], 
            unin = "m", 
            unout = "m2",
            duplicates="remove"), silent=T))== "try-error"){ 
      
            resul[[i+5]]<- NA 
            
      } else {
        
      resul[[i+5]]<- LoCoH.k(puntos, k=k[i], 
                         unin = "m", 
                         unout = "m2",
                         duplicates="remove")
    }
  }
    
  
  for (i in 1:5){
    if (class (try (LoCoH.r(puntos, r=r[i], 
                            unin = "m", 
                            unout = "m2",
                            duplicates="remove"), silent=T))== "try-error"){ 
      
      resul[[i+10]]<- NA 
      
    } else {
      
      resul[[i+10]]<- LoCoH.r(puntos, r=r[i], 
                           unin = "m", 
                           unout = "m2",
                           duplicates="remove")
    }
  }
  
  
  for (i in 1:5){
    if (class (try (LoCoH.a(puntos, a=a[i], 
                            unin = "m", 
                            unout = "m2",
                            duplicates="remove"), silent=T))== "try-error"){ 
      
      resul[[i+15]]<- NA 
      
    } else {
      
      resul[[i+15]]<- LoCoH.a(puntos, a=a[i], 
                           unin = "m", 
                           unout = "m2",
                           duplicates="remove")
    }
  }
  
  
  coord<- as.data.frame(puntos)
  names (coord)<- c("x", "y")
  re<- kernel2 (coord, thresh)
  resul <- c( resul, re)
  length (resul)
  resul
}

# área
calc_area<- function (shape, col, row){ 
  r <- raster(ncol=col, nrow=row)
  extent(r) <- extent(-2, 6, -2, 6)
  rp <- rasterize(shape, r)
  real<- reclassify (rp, c(-Inf, 0, 0, 1, +Inf, 1))
  real
}

## calcular áreas predichas vs. reales. accuracy
accuracy<- function (real, pred){
  pred [which (is.na (pred@data@values))]<- 0 
  real [which (is.na (real@data@values))]<- 0 
  res<- ((real + 1)* real) + pred
  a<- length (res@data@values [res@data@values==0])
  b<- length (res@data@values [res@data@values==3])
  d<- length (res@data@values [res@data@values==1])
  e<- length (res@data@values [res@data@values==2])
  res<- c(a, b, d, e)
  names (res)<- c("0sbien", "1sbien", "0smal", "1smal")
  res
}


## muestreos con errores, 
muestreos_error<- function (shape, n, p_error){
  m<- n * (1- p_error/100)
  er<- n * (p_error/100)
  
  cuadrado<- Polygon(cbind (x=c(-2,-2, 7, 7, -2), y=c(-3, 7, 7, -3, -3)))
  c_pol <- Polygons(list(Polygon(cuadrado)), ID="1")
  c_shape<- SpatialPolygons(list(c_pol))
  crs (c_shape)<- CRS("+proj=longlat +datum=WGS84") 
  shape_errores<- c_shape - shape
  
  n_rand1<- spsample (shape, n= m, type="random")
  n_rand2<- spsample (shape_errores, n= er, type="random")
  n_rand<- n_rand1 + n_rand2
  n_reg1<- spsample (shape, n= m, type="regular")
  n_reg<- n_reg1 + n_rand2
  n_clus1<- spsample (shape, n= m, type="clustered", nclusters= m/2)
  n_clus<- n_clus1 + n_rand2
  
list (rand=n_rand, reg=n_reg, clus=n_clus)
}

## contar el numero de fragmentos que se generan
## la función local convex hull tiene un bug y hay problemas con holes en las geometrías. 
## para evitar el error, si dan problemas les asigno NA.

n_frag_pred<- function (pred_shapes){
shapes<- pred_shapes [-c(1:3)]
res_n<- NULL
  for (i in 1:length (shapes)){
    if (class (try (gUnaryUnion (shapes[[i]]), silent=T))!= "try-error"){ 
      frag<- gUnaryUnion (shapes[[i]])
      frag2<- unlist (frag@polygons) [[1]]
      n<- length (frag2@Polygons)
      res_n<- c(res_n, n)  
    } else {
      n<- NA
      res_n<- c(res_n, n)  
    }
  }
res_n
}

## contar el numero de fragmentos que se generan cuando contamos con 5 parámetros por modelo
## la función local convex hull tiene un bug y hay problemas con holes en las geometrías. 
## para evitar el error, si dan problemas les asigno NA.

n_frag_pred2<- function (pred_shapes){
  shapes<- pred_shapes
  res_n<- NULL
  for (i in 1:length (shapes)){
    if (class (try (gUnaryUnion (shapes[[i]]), silent=T))!= "try-error"){ 
      frag<- gUnaryUnion (shapes[[i]])
      frag2<- unlist (frag@polygons) [[1]]
      n<- length (frag2@Polygons)
      res_n<- c(res_n, n)  
    } else {
      n<- NA
      res_n<- c(res_n, n)  
    }
  }
  res_n
}
