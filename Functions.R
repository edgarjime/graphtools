
if (lgLocal) {
Ruta <- "d:/AdministracionProyectos/Investigacion/LasagnaPlot/"
} else {
Ruta <- "/home/edgarcimat/Investigacion/Visualizacion/"
}

RutaD <-  paste(Ruta,"Data/",sep="")
RutaG <- paste(Ruta,"Borrador/",sep="")

library(reshape2)
library(ggplot2)
library(Rcpp)
library(VIM)
library(longCatEDA)
library(VIM)
library(RColorBrewer)
library(kohonen)
library(stringr)
library(data.table)
library(microbenchmark)

library(grid)
library(sp)

## Rcpp
sourceCpp("MatrixZero.cpp")

## Windows endings an begginnings localization
fnVentana <- function(nLargo,nVentana){
  Res <- nLargo%%nVentana
  Div <- nLargo%/%nVentana
  
  ## Result Table
  dtRes <- matrix(NA, nrow=nVentana, ncol=2)
  dtAux <- rep(0, nVentana)
  ## Ponemos los elementos de ajuste al final
  if (Res!=0) {
    dtAux[(length(dtAux)-Res+1):length(dtAux)] <- 1
  }
    
  #print(dtAux)
  
  dtRes[nVentana,2] <- nLargo
  dtRes[nVentana,1] <- nLargo-(Div-1) - dtAux[nVentana]
  
  for (i in (nVentana-1):1) {
    dtRes[i,2] <- dtRes[i+1,1]-1
    dtRes[i,1] <- dtRes[i,2]-(Div-1)-dtAux[i]
  }

  return(dtRes)
}



## Functions for analisis of missing data 
fnNASumUlt <- function(dtX,nLast) {
  return(sum(is.na(dtX[(length(dtX)-nLast+1):length(dtX)])))
}

fnNASumFirst <- function(dtX,nFirst) {
  return(sum(is.na(dtX[1:nFirst])))
}


fnNASum <- function(dtX) {
  return(sum(is.na(dtX)))
}

fnFullSum <- function(dtX) {
  return(sum(!is.na(dtX)))
}


fnLargoNA <- function (a){
  rl <- rle(is.na(c(a)))
  sal <- max(rl$lengths[rl$values])
  salF <- ifelse(is.infinite(sal),0,sal)
}

## Densities in specific points 
fnSumUltD <- function(dtX,nLast) {
  return(sum(dtX[(length(dtX)-nLast):length(dtX)], na.rm=TRUE))
}


## Extreme values labels 
fnAtipicoEtiq <- function(x,Bajo,Alto) {
    res <- ifelse(x<Bajo,-Inf,x)
    res <- ifelse(res>Alto,Inf,res)
    res
}


fnCeroSum <- function(dtX) {
        return(sum((dtX==0)))
}

## Conversion of result from values of square are function 
## to x,y coordinates
## dtVector has 5 elements, area, i, d[j], d2[j], d1[j]
fn_cr_zero <- function (dtVector) {
  dtSalida <- matrix(NA,ncol=2,nrow=4)
  i <- dtVector[2]
  d <- dtVector[3]
  d2 <- dtVector[4]
  d1 <- dtVector[5]
  
  ## Coordinates in 0
  xa <- xc <- d1 + 1
  xb <- xd <- d2 - 1
  yd <- yc <- i
  ya <- yb <- i-(i-d-1)
  
  ## Coordinates in 1
  dtSalida[1,] <- c(ya,xa)
  dtSalida[2,] <- c(yb,xb)
  dtSalida[3,] <- c(yd,xd)
  dtSalida[4,] <- c(yc,xc)
  dtSalida <- dtSalida + 1
  return(dtSalida)
}

## find square matrices of NA inside a matrix
## m1c is the matrix of data and minArea is a constant
fnPolCuadrados <- function(m1c,minArea) {
  nnren <- nrow(m1c)
  nncol <- ncol(m1c)
  m1c  <- as.matrix(m1c)
  
  ## Shadow matrix for squea algorithm (NA=0 y !NA=1)
  idx <- which(!is.na(m1c))
  m1c[idx] <- 1
  idx <- which(is.na(m1c))
  m1c[idx] <- 0
  
  
  ## Square algorithm
  conArea <- 0.005
  lsCuadradoTotal <- list()
  lsCuadradoArea <- list()
  for (i in 1:(nnren*nncol)) {
    ## i <- 1
    ## Find greatest submatrix (fn_zero_matrix is in MatrixZero.cpp)
    dtCuadrado <- fn_zero_matrix(m1c)
    
    if (dtCuadrado[1]<conArea*nncol*nnren)
      break(i)
    
    dtCoordenadas <- fn_cr_zero(dtCuadrado)
    colnames(dtCoordenadas) <- c("x", "y")
    
    idx  <- dtCoordenadas[1,"x"]:dtCoordenadas[3,"x"]
    idy  <- dtCoordenadas[1,"y"]:dtCoordenadas[2,"y"]
    
    m1c[idx,idy]  <- 1
    
    ## Save submatrix
    dtCoordenadas <- data.frame(dtCoordenadas,id=i)
    lsCuadradoTotal <- c(lsCuadradoTotal,list(dtCoordenadas))
    lsCuadradoArea <- c(lsCuadradoArea,dtCuadrado[1])
  }
  
  dtCuadradosTotal <- rbindlist(lsCuadradoTotal)
  dtCuadradosTotalArea <- c(unlist(lsCuadradoArea))
  
  return(list(Coor=dtCuadradosTotal,Areas=dtCuadradosTotalArea))
}

## Find area in squares areas
fnAreaPol <- function(m1, dtCuadradosTotal, dtCuadradosTotalArea) {
  ## Shadow matrix (NA=1, other data=0)
  idx <- which(!is.na(m1))
  m1[idx] <- 0
  idx <- which(is.na(m1))
  m1[idx] <- 1
  
  ## Realizar el procedimiento anterior para los cuadrados
  dtCooNACuad <- which(m1==1, arr.ind=TRUE)
  dtCooNoNACuad  <- which(m1==0, arr.ind=TRUE)
  iCuad <- length(table(dtCuadradosTotal$id))
  
  sumTCuad <- rep(0, iCuad)
  for (i in 1:iCuad) {
    sump <- sum(point.in.polygon(dtCooNACuad[,2], dtCooNACuad[,1], unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"y"]),  unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"x"]))>0)
    sumTCuad[i] <- sump
  }
  
  sumTNCuad <- rep(0, iCuad)
  for (i in 1:iCuad) {
    sump <- sum(point.in.polygon(dtCooNoNACuad[,2], dtCooNoNACuad[,1], unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"y"]),  unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"x"]))>0)
    sumTNCuad[i] <- sump
  }
  
  
  dtPor0Cuad <- sum(sumTCuad)/nrow(dtCooNACuad)
  dtPor0NoCuad  <- sum(sumTNCuad)/nrow(dtCooNoNACuad)
  dtPorCuadSep <- dtCuadradosTotalArea/nrow(dtCooNACuad)
  return(list(Por0=dtPor0Cuad, Por0N=dtPor0NoCuad, PorCuad=dtPorCuadSep ))
  
}


## Graphical functions
## dtDatosGra matrix of data with a colun of groups
## lgGroups: mark groups (horizontal) in the table
## dtVentana: mark vertical lines 
## lsResCuadF: localization of "squares areas" inside the graph
fnGraficaVacio <- function(dtDatosGra,IDSample, Tipo, nnCLuster, lgOrd, nnRepClus,AtipRetail, lgGrupos=FALSE,dtVentana=NULL, lsResCuadF=NULL) {
  
  lgCuadrados <- ifelse(is.null(lsResCuadF),FALSE,TRUE)
  
 
  ## los porcentajes solo funcionan para n=1000
  nnVis <- 1000
  dtEtiquetasPosY <- c(1,seq(from=200, to=nnVis, by=200))
  dtEtiquetasY <- paste(c(1,seq(from=20, to=100, by=20)),"%")
  
  nnSemanas <- ncol(dtDatosGra) - 3
  dtEtiquetasX <- c(1,seq(from=4, to=nnSemanas, by=4))
  
  dtDatosKGLim <- dtDatosGra
  nCLusters <- nnCLuster
  dtDatosKGLim[,"ClienteOrden"] <- 1:nrow(dtDatosKGLim)
  
  dtGrupos <- table(dtDatosKGLim[,"Grupo"])
  dtGrupos <- c(1, dtGrupos)
  dtGrupos <- cumsum(dtGrupos)
  
  ##dtDatosKGExp  <- melt(dtDatosKG[,c(-1,-3)][1:10,], id.vars = "ClienteOrden", measure.vars = colnames(dtDatosKG[,4:ncol(dtDatosKG)]))
  ##dtDatosKGLim  <- melt(dtDatosKGLim[,c(-2,-3)], id.vars = c("ClienteOrden"), measure.vars = colnames(dtDatosKGLim[,4:ncol(dtDatosKGLim)]))
  ##colnames(dtDatosKGLim) <- c("Cliente","Dia","Total")
  
  dtDatosKGLim  <- reshape2::melt(dtDatosKGLim[,c(-2,-3)], id.vars = c("ClienteOrden"), measure.vars = colnames(dtDatosKGLim[,4:ncol(dtDatosKGLim)]))
  colnames(dtDatosKGLim) <- c("Cliente","Dia","Total")
  
  ## Atypical data labels
  dtDatosKGLim[,"Total"] <- fnAtipicoEtiq(dtDatosKGLim[,"Total"],AtipRetail[1],AtipRetail[2])
  
  idxIB <- which(dtDatosKGLim[,"Total"]==-Inf)
  idxIA <- which(dtDatosKGLim[,"Total"]==Inf)
  dtDatosKGLim[idxIB,"Total"] <- NA
  dtDatosKGLim[idxIA,"Total"] <- NA
  
  nnCortes <- 5
  dtDatosKGLim[,"Total"] <- cut(dtDatosKGLim[,"Total"], breaks=nnCortes, dig.lab = 2)
  
  ## Correct label
  dtDatosKGLim$Total <- factor(dtDatosKGLim$Total, levels = c("Atip_Low", levels(dtDatosKGLim$Total) ,"Atip_High"))
  
  dtDatosKGLim[idxIB,"Total"] <- "Atip_Low"
  dtDatosKGLim[idxIA,"Total"] <- "Atip_High"
  
  dtColores <- c("red", brewer.pal(n = nnCortes, name = "YlGn"),"purple")
  
  p <- ggplot(dtDatosKGLim, aes(Dia, Cliente, z= Total)) + geom_tile(aes(fill = Total)) +
    scale_fill_manual(na.value = '#87ceeb', values = dtColores)
    #scale_fill_continuous(na.value = '#FFFFFF',low = "#717171", high = "#000000")
  
  ## Labels for graphing
  dtTitulo <- paste("Sales($) per client and week","-----" ,"Fraction of Missing Data in square areas-",round(lsPorCuad$Por0,3),sep="" )
  p <- p + labs(x="Week",y="Client",title=dtTitulo) + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  ## Percentage label needs more work (1000 only)
  p <- p + labs(fill = "Sales($)") + scale_x_discrete(breaks=dtEtiquetasX) + scale_y_continuous(labels = function(x) paste0(x/10, "%"), sec.axis = dup_axis())
  
  if (lgGrupos) {
    for (iLinea in 1:length(dtGrupos)){
      p <- p + geom_hline(yintercept=dtGrupos[iLinea], linetype="dashed", color = "black", size=0.25)
    }
  }
  
  if (!is.null(dtVentana)) {
      for (iVentana in 1:(nrow(dtVentana)-1)) {
        p <- p + geom_vline(xintercept=dtVentana[iVentana,2]+0.5, linetype="dashed", color = "black", size=0.25)
      }
  }
  
  
  if (lgCuadrados) {
    ##  i <- 1
    dtCuadradosTotal <- lsResCuadF$Coor
    iCuad <- length(table(dtCuadradosTotal$id))
    for (i in 1:iCuad) {
      dtPol <- data.frame(Dia=unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"y"])+c(-0.5,0.5,0.5,-0.5),Cliente=unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"x"]),Total=0)
      Base <- (dtPol[2,"Dia"]-dtPol[1,"Dia"]) 
      Altura <- (dtPol[3,"Cliente"]-dtPol[1,"Cliente"])
      
      conTam <- 1000
      conAltura <- 31
      if (Altura > conAltura) {
        PosBase <- dtPol[1,"Dia"] + Base/2
        PosAltura <- dtPol[1,"Cliente"]  + Altura/2
        dtEtiquetaCuad <- data.frame(Dia=PosBase, Cliente=PosAltura,Total=0)
        p <- p + geom_polygon(data=dtPol, linetype="longdash", color = "orange3", size=0.8, fill=NA)
        p <- p + geom_label(data=dtEtiquetaCuad,label=paste(Base, paste(round((Altura/conTam)*100,2),"%",sep=""), sep="-"), size=2.2)
      }
    }
  }
  
  p <- p + theme(legend.position="bottom",legend.key.size = unit(0.5,"line")) + theme(legend.background = element_rect(fill="lightblue",
                                                                                    size=0.3, linetype="solid", 
                                                                                    colour ="darkblue"))
  #p <- p + theme(legend.text = element_text(size=8))
  
  p <- p + theme(plot.title = element_text(size=7),
                 legend.title = element_text(size = 5), 
                 legend.text = element_text(size = 7))
  
  Nombre <- paste(RutaG, "Grafica8_SalesProductsRetailGroup_", IDSample, "_", Tipo,"_Clus_",  nnCLuster, "_Order_",lgOrd,"_Cuad_", lgCuadrados,"_Rep_" ,nnRepClus,".jpeg", sep="")
  jpeg(Nombre, width = 4800, height = 3000, units = "px", res = 800)
  print(p)
  dev.off()
  
}

