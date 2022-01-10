## Establecer rutas y Funciones

lgLocal =TRUE
source("Functions.R")

dtDatosCliM = read.csv("DataVignette.csv")
AtipOnline <- quantile(unlist(dtDatosCliM[,-1]), c(0.01,0.99),na.rm=TRUE)
nnSemanas <- ncol(dtDatosCliM) - 1
dtEtiquetasX <- seq(from=1, to=nnSemanas, by=4)

## Analysis of NA
## Total and percentiles
sum(is.na(dtDatosCliM[,-1]))
dtVacioRen  = apply(dtDatosCliM[,-1],1,fnNASum)
quantile(dtVacioRen)

## SOM Classfication

lgSOM = FALSE
dtCluster = NA
if (lgSOM) {
    dtDatosKGSom <- dtDatosCliM
    dtDatosKGSom[,-1] <- ifelse(is.na(dtDatosKGSom[,-1]),1,0)

    ## SOM Comparison
    nnGrid <- 5*sqrt(ncol(dtDatosKGSom[,-1])*nrow(dtDatosKGSom[,-1]))
    nnMapSize <- sqrt(nnGrid)

                                        # general topology
    som_grid <- somgrid(xdim = 30, ydim=30, topo="hexagonal")

                                        # train the SOM, number of iterations,
                                        # the learning rates, and the neighbourhood
    som_model <- som(as.matrix(dtDatosKGSom[,-1]),
                     grid=som_grid,
                     rlen=100,
                     alpha=c(0.05,0.01),
                     keep.data = TRUE )

    plot(som_model, type="changes")
    plot(som_model, type="count")
    plot(som_model, type="dist.neighbours")
    plot(som_model, type="codes")

    nCLusters <- 16
    som_cluster <- cutree(hclust(dist(getCodes(som_model))), nCLusters)

    ## plot these results
    pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')
    plot(som_model, type="mapping", bgcol = pretty_palette[som_cluster], main = "Clusters")
    add.cluster.boundaries(som_model, som_cluster)

    ## Hacer el match entre el cluster
    idx <- match(som_model$unit.classif,1:length(som_cluster))
    dtCluster <- som_cluster[idx]
}

## 3

## Conversion to graph using kmeans and sampling
dtDatosKG <- data.frame(Cliente=dtDatosCliM[,"Cliente"], ClienteOrden=1:nrow(dtDatosCliM),
                        Grupo=dtCluster,dtDatosCliM[,-1])
colnames(dtDatosKG) <- c("Cliente","ClienteOrden","Grupo",1:(ncol(dtDatosCliM)-1))

lgKmeans = TRUE
if (lgKmeans) {
    dtLimites <- fnVentana(ncol(dtDatosKG[,-c(1:3)]),4)
    conShift <- 3
    ## Kmeans Class
    dtSec4 <- apply(dtDatosKG[,(dtLimites[4,1]+conShift):(dtLimites[4,2]+conShift)],1,fnNASum)
    dtSec3 <- apply(dtDatosKG[,(dtLimites[3,1]+conShift):(dtLimites[3,2]+conShift)],1,fnNASum)
    dtSec2 <- apply(dtDatosKG[,(dtLimites[2,1]+conShift):(dtLimites[2,2]+conShift)],1,fnNASum)
    dtSec1 <- apply(dtDatosKG[,(dtLimites[1,1]+conShift):(dtLimites[1,2]+conShift)],1,fnNASum)

    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    indexT <- apply(cbind(dtSec1,dtSec2,dtSec3,dtSec4),2,range01)
    nCluster <- 16
    groupsKT <- kmeans(indexT,centers = nCluster,nstart = 50, iter.max = 50)
    dtDatosKG[,"Grupo"] <- groupsKT$cluster
    dtW <- c(1,4,16,64)

    lgIndex <- TRUE
    if (lgIndex) {

        groups <- dtDatosKG[,"Grupo"]
        lgOrdenGrupo <- FALSE
        lgEscalera <- TRUE
        if (lgOrdenGrupo) {
            ## Indizado por grupos
            dtSumNAT <- apply(dtDatosKG[,-c(1:3)],1,fnFullSum)
            indexgroups <- tapply(dtSumNAT,groups,sum,na.rm=TRUE)#/(table(groups)*nnTime)
            dtOrdenNuevo <- order(indexgroups, decreasing=c("TRUE","TRUE","TRUE","TRUE"))
        }

        ## Indizado por escalera
        if (lgEscalera) {
            index4 <- tapply(dtSec4, groups,sum,na.rm=TRUE)/(table(groups)*(dtLimites[4,2] - dtLimites[4,1] +1))
            index3 <- tapply(dtSec3, groups,sum,na.rm=TRUE)/(table(groups)*(dtLimites[3,2] - dtLimites[3,1] +1))
            index2 <- tapply(dtSec2, groups,sum,na.rm=TRUE)/(table(groups)*(dtLimites[2,2] - dtLimites[2,1] +1))
            index1 <- tapply(dtSec1, groups,sum,na.rm=TRUE)/(table(groups)*(dtLimites[1,2] - dtLimites[1,1] +1))
            indexgroups <- index4
            indexa <- cbind(index1,index2,index3,index4)
            dtOrdenNuevo <- order(apply(dtW*t(indexa),2,sum),method="radix", decreasing=c("FALSE"))
        }

        dtTabla <- cbind(Nuevo=1:nCluster,Viejo=dtOrdenNuevo ,Etiqueta=as.numeric(names(indexgroups[dtOrdenNuevo])))
        groupsP <- dtTabla[match(groups,dtTabla[,"Etiqueta"]),"Nuevo"]
        groups <- groupsP
        dtDatosKG[,"Grupo"] <- groups
    }
}

## Manual selection between Alg1 and Alg2 and the Alg1
idx <- order(dtSec4, dtSec3, dtSec2, dtSec1,
             apply(dtDatosKG[,-c(1:3)],1,fnNASum),
             apply(dtDatosKG[,-c(1:3)],1,sum,na.rm=TRUE),
             method="radix",
             decreasing=c("FALSE","FALSE","FALSE","FALSE","FALSE","TRUE"))

idx <- order(dtDatosKG[,"Grupo"], dtSec4, dtSec3, dtSec2, dtSec1,
             apply(dtDatosKG[,-c(1:3)],1,fnNASum),
             apply(dtDatosKG[,-c(1:3)],1,sum,na.rm=TRUE),
             method="radix",
             decreasing=c("FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","TRUE"))



## Veryfing elemts for plotting
dtDatosKGLim <- dtDatosKG[idx,]
dtDatosKGLim[,"ClienteOrden"] <- 1:nrow(dtDatosKGLim)


## Making a ramdom sample
lgRandom <- FALSE
if(lgRandom) {
  nnVis <- 1000
  idxEleccion <- sort(sample(1:nrow(dtDatosKG), nnVis))
} else {
  idxEleccion <- 1:nrow(dtDatosKG)
}


## Data limits
dtDatosKGLim <- dtDatosKGLim[idxEleccion,]
dtDatosKGLim[,"ClienteOrden"] <- 1:nrow(dtDatosKGLim)

dtGrupos <- table(dtDatosKGLim[,"Grupo"])
dtGrupos <- c(1, dtGrupos)
dtGrupos <- cumsum(dtGrupos)


## Square finding algoritm modification
m1 <- dtDatosKGLim[,4:ncol(dtDatosKGLim)]
conArea <- 0.001
lsResCuad <- fnPolCuadrados(m1,conArea)
lsPorCuad <- fnAreaPol(as.matrix(m1), lsResCuad$Coor, lsResCuad$Areas)
# plot(lsPorCuad$PorCuad)


## 4
dtDatosKGLim  <- melt(dtDatosKGLim[,c(-1,-3)], id.vars = c("ClienteOrden"), measure.vars = colnames(dtDatosKGLim[,4:ncol(dtDatosKGLim)]))
colnames(dtDatosKGLim) <- c("ClienteOrden","Dia","Total")

## Atypical labeling
dtDatosKGLim[,"Total"] <- fnAtipicoEtiq(dtDatosKGLim[,"Total"],AtipOnline[1],AtipOnline[2])

idxIB <- which(dtDatosKGLim[,"Total"]==-Inf)
idxIA <- which(dtDatosKGLim[,"Total"]==Inf)
dtDatosKGLim[idxIB,"Total"] <- NA
dtDatosKGLim[idxIA,"Total"] <- NA

nnCortes <- 5
dtDatosKGLim[,"Total"] <- cut(dtDatosKGLim[,"Total"], breaks=nnCortes)

## Right level for label
dtDatosKGLim$Total <- factor(dtDatosKGLim$Total, levels = c("Atip_Low", levels(dtDatosKGLim$Total) ,"Atip_High"))

dtDatosKGLim[idxIB,"Total"] <- "Atip_Low"
dtDatosKGLim[idxIA,"Total"] <- "Atip_High"
dtColores <- c("red", brewer.pal(n = nnCortes, name = "YlGn"),"purple")

## Original grapch
p <- ggplot(dtDatosKGLim, aes(Dia, ClienteOrden, z= Total)) + geom_tile(aes(fill = Total)) +
     scale_fill_manual(na.value = '#87ceeb', values = dtColores)
##  scale_fill_continuous(na.value = '#87ceeb',low = "darkgray", high = "red")

## Graph labels
dtTitulo <- paste("Sales($) per client and week","-----" ,"Fraction of Missing Data in square areas-",round(lsPorCuad$Por0,3),sep="" )
p <- p + labs(x="Week",y="Client",title=dtTitulo) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p <- p + labs(fill = "Sales($)") + scale_x_discrete(breaks=dtEtiquetasX)
p <- p + theme(plot.title = element_text(size=10))
p <- p + theme(legend.position="bottom") + theme(legend.background = element_rect(fill="lightblue",
                                                                                  size=0.5, linetype="solid",
                                                                                  colour ="darkblue"))

p <- p + theme(legend.text = element_text(size=7), legend.title = element_text(size=7), legend.key.size = unit(.25, "cm"))


dtVentana <- dtLimites
if (!is.null(dtVentana)) {
  for (iVentana in 1:(nrow(dtVentana)-1)) {
    p <- p + geom_vline(xintercept=dtVentana[iVentana,2]+0.5, linetype="dashed", color = "black", size=0.25)
  }
}


lgGrupos <- TRUE
if (lgGrupos) {
  for (iLinea in 1:length(dtGrupos)){
    p <- p + geom_hline(yintercept=dtGrupos[iLinea], linetype="dotted", color = "black", size=0.3)
  }
}

lgCuadrados <- TRUE
if (lgCuadrados) {
  ##  i <- 1
  dtCuadradosTotal <- lsResCuad$Coor
  iCuad <- length(table(dtCuadradosTotal$id))
  for (i in 1:iCuad) {
    dtPol <- data.frame(Dia=unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"y"])+c(-0.5,0.5,0.5,-0.5),ClienteOrden=unlist(dtCuadradosTotal[dtCuadradosTotal$id==i,"x"]),Total=0)
    Base <- (dtPol[2,"Dia"]-dtPol[1,"Dia"])
    Altura <- (dtPol[3,"ClienteOrden"]-dtPol[1,"ClienteOrden"])

    conAltura <- 31
    if (Altura > conAltura) {
      PosBase <- dtPol[1,"Dia"] + Base/2
      PosAltura <- dtPol[1,"ClienteOrden"]  + Altura/2
      dtEtiquetaCuad <- data.frame(Dia=PosBase, ClienteOrden=PosAltura,Total=0)
      p <- p + geom_polygon(data=dtPol, linetype="longdash", color = "orange3", size=0.2, fill=NA)
      p <- p + geom_label(data=dtEtiquetaCuad,label=paste(Base, Altura, sep="-"),size=2.2)
    }
  }
}


Nombre <- paste(RutaG, "FigureVignette_X", ".jpeg", sep="")
jpeg(Nombre, width = 4800, height = 3000, units = "px", res = 800)
print(p)
dev.off()



