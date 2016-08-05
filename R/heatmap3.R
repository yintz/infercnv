##' Please note this code is from the library GMD
##' All credit for this code goes to GMD's author's.
##' I do not recommend using this version of the code, which
##' has been poorly modified for our use but recommend using
##' the official version from the package GMD
##' https://cran.r-project.org/web/packages/GMD/index.html
##' A copy of gtools::invalid
##' 
##' see \code{invalid} in package:gtools for details
##' @title Test if a value is missing, empty, or contains only NA or NULL values
##' @param x value to be tested
.invalid <- 
  function(x) 
{
  if (missing(x) || is.null(x) || length(x) == 0) 
    return(TRUE)
  if (is.list(x)) 
    return(all(sapply(x, .invalid)))
  else if (is.vector(x)) 
    return(all(is.na(x)))
  else return(FALSE)
}

##' Please note this code is from the library GMD
##' All credit for this code goes to GMD's author's.
##' I do not recommend using this version of the code, which
##' has been poorly modified for our use but recommend using
##' the official version from the package GMD
##' https://cran.r-project.org/web/packages/GMD/index.html
##' Call a function with arguments
##'
##' Call a function with arguments
##' @title Call a function with arguments
##' @param FUN function or function name
##' @param ... unnameed function arguments
##' @param MoreArgs named (or unnameed) function arguments
.call.FUN <-
  function(FUN,...,MoreArgs)
{
  FUN <- match.fun(FUN)
  tmp.MoreArgs <- list(...)
  if (!.invalid(MoreArgs)){
    if (length(MoreArgs)>=1){
      for (i in 1:length(MoreArgs)) tmp.MoreArgs[[names(MoreArgs)[i]]] <- MoreArgs[[i]]
    }   
  }
  ret <- do.call(FUN, tmp.MoreArgs)
  if ("call" %in% names(ret)){
    ret$call <- match.call()
  }
  if ("call" %in% names(attributes(ret))){
    attr(ret,"call") <- match.call()
  }
  return(ret)
}

##' Please note this code is from the library GMD
##' All credit for this code goes to GMD's author's.
##' I do not recommend using this version of the code, which
##' has been poorly modified for our use but recommend using
##' the official version from the package GMD
##' https://cran.r-project.org/web/packages/GMD/index.html
##' Scale values to make them follow Standard Normal Distribution
##'
##' Scale values to make them follow Standard Normal Distribution
##' @title Scale values to make them follow Standard Normal Distribution
##' @param x numeric
##' @param scale character, indicating the type to scale.
##' @param na.rm logical
##' @return an object with the same dimention of `x'.
.scale.data <-
  function(x,scale,na.rm=TRUE)
{
  if(scale=="row"){
    x <- sweep(x,1,rowMeans(x,na.rm=na.rm),FUN="-")
    sx <- apply(x,1,sd,na.rm=na.rm)
    x <- sweep(x,1,sx,FUN="/")
  } else if(scale=="column"){
    x <- sweep(x,2,colMeans(x,na.rm=na.rm),FUN="-")
    sx <- apply(x,2,sd,na.rm=na.rm)
    x <- sweep(x,2,sx,,FUN="/")
  }
  x
}

##' Please note this code is from the library GMD
##' All credit for this code goes to GMD's author's.
##' I do not recommend using this version of the code, which
##' has been poorly modified for our use but recommend using
##' the official version from the package GMD
##' https://cran.r-project.org/web/packages/GMD/index.html
##' Scale values to a new range: c(low, high)
##' @title Scale values to a new range.
##' @param x numeric
##' @param low numeric, lower bound of target values
##' @param high numeric, higher bound of target values
##' @return an object with the same dimention of `x'.
.scale.x <-
  function(x,low=0,high=1,na.rm=TRUE)
{
  if(identical(max(x,na.rm=na.rm),min(x,na.rm=na.rm))) NA
  a <- 1/(max(x)-min(x))
  b <- -min(x)/(max(x)-min(x))
  a*x+b
}

##' Please note this code is from the library GMD
##' All credit for this code goes to GMD's author's.
##' I do not recommend using this version of the code, which
##' has been poorly modified for our use but recommend using
##' the official version from the package GMD
##' https://cran.r-project.org/web/packages/GMD/index.html
##' Plot text
##'
##' Plot text
##' @title Plot text
##' @param x character, text to plot
##' @param cex
##' @param forecolor color of foreground
##' @param bg color of background
##' @param bordercolor color of border
##' @param axes as in \code{graphics:::plot}
##' @param ... additional arguments for \code{graphics:::text}
.plot.text <- function(x,xlim=c(0,1),ylim=c(0,1),cex=1,forecolor=par("fg"),bg=par("bg"),bordercolor=NA,axes=FALSE,...){
  if (.invalid(x)){
    x <- NULL
  }
  if (is.null(x)){
    x <- ""
  } else if (is.na(x)){
    x <- 'NA'
  }

  plot(xlim,ylim,type="n",ylab="",xlab="",xaxt="n",yaxt="n",axes=axes)
  rect(xleft=0, ybottom=0, xright=1, ytop=1, col=bg, border=bordercolor)
  text(0.5,0.5,x,cex=cex,...)
}

##' Please note this code is from the library GMD
##' All credit for this code goes to GMD's author's.
##' I do not recommend using this version of the code, which
##' has been poorly modified for our use but recommend using
##' the official version from the package GMD
##' https://cran.r-project.org/web/packages/GMD/index.html
##' This was originally heatmap.3.
heatmap.cnv <-
  function(x,
           
           ## whether a dissimilarity matrix
           diss=inherits(x,"dist"),

           ## dendrogram control
           Rowv=TRUE,
           Colv=TRUE,
           dendrogram=c("both","row","column","none"),

           ## dist object
           dist.row,
           dist.col,
           dist.FUN=gdist,
           dist.FUN.MoreArgs=list(method="euclidean"),
           
           ## hclust object
           hclust.row,
           hclust.col,
           hclust.FUN=hclust,
           hclust.FUN.MoreArgs=list(method="ward"),
           
           ## data scaling
           scale=c("none","row","column"),
           na.rm=TRUE,

           ## clustering control
           cluster.by.row=FALSE,
           cluster.by.col=FALSE,
           kr=NA,      
           kc=NA,
           row.clusters=NA,
           col.clusters=NA,

           ## image plot
           revR=FALSE,
           revC=FALSE,
           add.expr,
           
           ## mapping data to colors
           breaks,
           ## centering colors to a value
           x.center,
           ## colors
           color.FUN=gplots::bluered,
           ##
           ## block sepration
           sepList=list(NULL,NULL),
           sep.color=c("gray45","gray45"),
           sep.lty=1,
           sep.lwd=2,
           

           ## cell labeling
           cellnote,
           cex.note=1.0,
           notecol="cyan",
           na.color=par("bg"),

           ## level trace
           trace=c("none","column","row","both"),
           tracecol="cyan",
           hline,
           vline,
           linecol=tracecol,

           ## Row/Column Labeling
           labRow=TRUE, ## shown by default
           labCol=TRUE, ## shown by default
           srtRow=NULL,
           srtCol=NULL,
           sideRow=4,
           sideCol=1,
           ##
           margin.for.labRow,
           margin.for.labCol,
           ColIndividualColors,
           RowIndividualColors,
           cexRow,
           cexCol,
           labRow.by.group=FALSE,
           labCol.by.group=FALSE,
           

           ## plot color key + density info
           key=TRUE,
           key.title="Color Key",
           key.xlab="Value",
           key.ylab="Count",
           
           keysize=1.5,
           mapsize=9,
           mapratio=4/3,
           sidesize=3,
           cex.key.main=0.75,
           cex.key.xlab=0.75,
           cex.key.ylab=0.75,
           density.info=c("histogram","density","none"),
           denscol=tracecol,
           densadj=0.25,

           ## plot titles/labels
           main="Heatmap",
           sub="",
           xlab="",
           ylab="",
           cex.main=2,
           cex.sub=1.5,
           font.main=2,
           font.sub=3,
           adj.main=0.5,
           mgp.main=c(1.5,0.5,0),
           mar.main=3,
           mar.sub=3,
           ## plot ##
           if.plot=TRUE,
           
           ## plot of partition (left/top of heatmap)
           plot.row.partition=FALSE,
           plot.col.partition=FALSE,
           cex.partition=1.25,
           color.partition.box="gray45",
           color.partition.border="#FFFFFF",

           ## plot of summary (right/bottom of heatmap)
           plot.row.individuals=FALSE,
           plot.col.individuals=FALSE,
           plot.row.clusters=FALSE,
           plot.col.clusters=FALSE,
           plot.row.clustering=FALSE,
           plot.col.clustering=FALSE,

           ##
           plot.row.individuals.list=FALSE,
           plot.col.individuals.list=FALSE,
           plot.row.clusters.list=FALSE,
           plot.col.clusters.list=FALSE,
           plot.row.clustering.list=FALSE,
           plot.col.clustering.list=FALSE,
           
           ## for plot of clusters - row
           row.data=FALSE,
           ## for plot of clusters - col
           col.data=FALSE,

           ##
           if.plot.info=FALSE,
           text.box,
           cex.text=1.0,

           ## Force in layout info
           force_lmat=NULL,
           force_lwid=NULL,
           force_lhei=NULL,
           force_add=FALSE,

           ## extras
           ...
           )
{
  ## check input - take1 ##
  if (is.data.frame(x)){
    x <- as.matrix(x)
  }
  x.ori <- x
  
  if(!inherits(x,"dist") & !is.matrix(x)){
    stop("`x' should either be a matrix, a data.frame or a `dist' object.")
  }

  if (! sideRow %in% c(2,4)){
    stop('sideRow must be either 2 or 4.')
  }
  
  if (! sideCol %in% c(1,3)){
    stop('sideCol must be either 1 or 3.')
  }
  
  ## store input
  Rowv.ori <- Rowv
  Colv.ori <- Colv
  
  ## check
  dendrogram <- match.arg(dendrogram)
  if ( (dendrogram %in% c("both","row")) & !inherits(Rowv,"dendrogram") ){
    warning("Discrepancy: row dendrogram is asked;  Rowv is set to `TRUE'.")
    Rowv <- TRUE
  }

  if ( (dendrogram %in% c("both","col")) & !inherits(Colv,"dendrogram") ){
    warning("Discrepancy: col dendrogram is asked;  Colv is set to `TRUE'.")
    Colv <- TRUE
  }

  if (identical(Rowv, FALSE) | missing(Rowv)){
    if(!identical(cluster.by.row,FALSE)){
      warning("Discrepancy: No row dendrogram is asked; cluster.by.row is set to `FALSE'.")
      cluster.by.row <- FALSE
    }
  } else {
    if(!identical(cluster.by.row,TRUE)){
      warning("Discrepancy: row dendrogram is asked; cluster.by.row is set to `TRUE'.")
      cluster.by.row <- TRUE
    }
  }

  if (identical(Colv, FALSE) | .invalid(Colv)){
    if(!identical(cluster.by.col,FALSE)){
      warning("Discrepancy: No col dendrogram is asked; cluster.by.col is set to `FALSE'.")
      cluster.by.col <- FALSE
    }
  } else {
    if(!identical(cluster.by.col,TRUE)){
      warning("Discrepancy: col dendrogram is asked; cluster.by.col is set to `TRUE'.")
      cluster.by.col <- TRUE
    }
  }
  
  if (!.invalid(kr)){
    if (is.numeric(kr)){
      if(!plot.row.partition){
        warning("Discrepancy: kr is set, therefore plot.row.partition is set to `TRUE'.")
        plot.row.partition <- TRUE
      }
    }
  }
  
  if (!.invalid(kc)){
    if (is.numeric(kc)){
      if(!plot.col.partition){
        warning("Discrepancy: kc is set, therefore plot.col.partition is set to `TRUE'.")
        plot.col.partition <- TRUE
      }
    }
  }
  
  ## generate dist.obj - row/col ##
  if (inherits(x,"dist")){
    dist.row <- dist.col <- x ## dist.obj
    x <- as.matrix(x)
    mat.row <- mat.col <- x ## mat.obj
    symm <- TRUE
  } else if (is.matrix(x)){
    symm <- isSymmetric(x)
    if (diss){
      if (!symm){
        stop("Dissimilary matrix should be symmetric. Please set `diss' to FALSE if `x' is not dissimilary matrix.")
      } else {
        flush.console()
      }
      mat.row <- mat.col <- x
      dist.row <- dist.col <- as.dist(x)
    } else{
      if (cluster.by.row) {
        if (.invalid(dist.row)){
          dist.row <- .call.FUN(dist.FUN,x,MoreArgs=dist.FUN.MoreArgs)
        }
        mat.row <- as.matrix(dist.row)
      } else {
        dist.row <- NULL
        mat.row <- NULL
      }
      if (cluster.by.col) {
        if (.invalid(dist.col)){
          dist.col <- .call.FUN(dist.FUN,t(x),MoreArgs=dist.FUN.MoreArgs)
        }
        mat.col <- as.matrix(dist.col)
      } else {
        dist.col <- NULL
        mat.col <- NULL
      }
    }
  }

  ## check input - take2: di ##
  di <- dim(x)
  
  if(length(di)!=2 || !is.numeric(x)){
    stop("`x' should only contain `numeric' values and can be converted to a 2-D matrix.")
  }

  ## parse param ##
  scale <- if(symm && .invalid(scale)) "none" else match.arg(scale) ## no scale on symmetric matrix
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  dist.FUN <- match.fun(dist.FUN)
  hclust.FUN <- match.fun(hclust.FUN)
  color.FUN <- match.fun(color.FUN)
  
  ## NG if both breaks and scale are specified ##
  if(!.invalid(breaks) & (scale!="none")){
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.",
            "Please consider using only one or the other.")
  }

  ## nr and nc ##
  nr <- di[1]
  nc <- di[2]

  ## check input - take3: nr,nc ##
  if(nr <=1 || nc <=1)
    stop("`x' must have at least 2 rows and 2 columns")

  ## font size of row/col labels ##
  cexRow0 <- 0.2+1/log10(nr)
  cexCol0 <- 0.2+1/log10(nc)

  if (.invalid(cexRow)) {
    cexRow <- cexRow0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for cexRow.')
    cexRow <- cexRow0*cexRow
  }
  if (.invalid(cexCol)) {
    cexCol <- cexCol0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for cexCol.')
    cexCol <- cexCol0*cexCol
  }

  ## cellnote ##
  ## ##if(.invalid(cellnote)) cellnote <- matrix("",ncol=ncol(x),nrow=nrow(x))

  ## ------------------------------------------------------------------------
  ## parse dendrogram ##
  ## ------------------------------------------------------------------------  
  
  if (missing(Rowv)) Rowv <- FALSE
  
  if (.invalid(Colv)) Colv <- if(symm) Rowv else FALSE
  if (Colv=="Rowv") {
    if ((!isTRUE(Rowv) | !symm) ){
      Colv <- FALSE
      warning("`Colv' is specified to use \"Rowv\", but either `Rowv' is invalid or `x' is not symmetric; Colv is suppressed.")
    } else{
      Colv <- Rowv
    }
  }
  
  ## ------------------------------------------------------------------------
  ## generate hclust.obj - row/col
  ## ------------------------------------------------------------------------
  flush.console()

  if ( (!inherits(Rowv,"dendrogram") & !identical(Rowv,FALSE)) | (cluster.by.row & .invalid(row.clusters))){
    if (.invalid(hclust.row)){
      hclust.row <- .call.FUN(hclust.FUN,dist.row,MoreArgs=hclust.FUN.MoreArgs)
    } else {
      if (length(hclust.row$order) != nr){
        stop("`hclust.row' should have equal size as the rows.")
      }
    }
  } else{
    hclust.row <- NULL
  }
  
  if(symm){
    hclust.col <- hclust.row
  }

  if ( !inherits(Colv,"dendrogram") & !identical(Colv,FALSE) | (cluster.by.col & .invalid(col.clusters))){
    if (.invalid(hclust.col)){
      hclust.col <- .call.FUN(hclust.FUN,dist.col,MoreArgs=hclust.FUN.MoreArgs)
    } else {
      if (length(hclust.col$order) != nc){
        stop("`hclust.col' should have equal size as the cols.")
      }
    }
  } else {
    hclust.col <- NULL
  }

  ## ------------------------------------------------------------------------
  ## generate hclust.obj - row/col
  ## ------------------------------------------------------------------------
  ddr <- ddc <- NULL
  ## get the dendrograms and reordering row/column indices - row ##
  if(inherits(Rowv,"dendrogram")){
    if (attr(Rowv,"members") != nr){
      stop("`Rowv' should contain equal size of members as the rows.")
    }
    ddr <- Rowv ## use Rowv 'as-is',when it is dendrogram
    rowInd <- order.dendrogram(ddr)
  } else if (is.integer(Rowv)){ ## compute dendrogram and do reordering based on given vector
    ddr <- as.dendrogram(hclust.row)
    ddr <-  reorder(ddr,Rowv) 
    rowInd <- order.dendrogram(ddr)
    if(nr != length(rowInd)){
      stop("`rowInd' is of wrong length.")
    }
  } else if (isTRUE(Rowv)){ ## if TRUE,compute dendrogram and do reordering based on rowMeans

    Rowv <- rowMeans(x,na.rm=TRUE)
    ddr <- as.dendrogram(hclust.row)
    ddr <- reorder(ddr,Rowv)
    rowInd <- order.dendrogram(ddr)
    if(nr !=length(rowInd)){
      stop("`rowInd' is of wrong length.")
    }
  } else{
    rowInd <- nr:1 ## from bottom.
  }
  
  ## get the dendrograms and reordering row/column indices - col ##
  if(inherits(Colv,"dendrogram")){
    if (attr(Colv,"members") != nc){
      stop("`Colv' should contain equal size of members as the cols.")
    }
    ddc <- Colv ## use Colv 'as-is',when it is dendrogram
    colInd <- order.dendrogram(ddc)
  } else if(identical(Colv,"Rowv")) {
    if(exists("ddr")){
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    } else{
      colInd <- rowInd
    }
  } else if(is.integer(Colv)){## compute dendrogram and do reordering based on given vector
    ddc <- as.dendrogram(hclust.col)
    ddc <- reorder(ddc,Colv)
    colInd <- order.dendrogram(ddc)
    if(nc != length(colInd))
      stop("`colInd' is of wrong length.")
  } else if (isTRUE(Colv)){## if TRUE,compute dendrogram and do reordering based on rowMeans
    Colv <- colMeans(x,na.rm=TRUE)
    ddc <- as.dendrogram(hclust.col)
    ddc <- reorder(ddc,Colv)
    colInd <- order.dendrogram(ddc)
    if(nc !=length(colInd))
      stop("`colInd' is of wrong length.")
  } else{
    colInd <- 1:nc ## from left
  }

  ## ------------------------------------------------------------------------
  ## check consistency
  ## ------------------------------------------------------------------------
  
  ## Xmisc::logme(dendrogram)
  ## Xmisc::logme(Colv)
  ## Xmisc::logme(Rowv)
  
  ## dendrogram - check consistency: Rowv ##
  if ( is.null(ddr) & (dendrogram %in% c("both","row"))){
    warning("Discrepancy: Rowv is invalid or FALSE, while dendrogram is `",
            dendrogram,"'. Omitting row dendogram.")
    if (is.logical(Colv) & (Colv.ori) & dendrogram=="both")
      dendrogram <- "column"
    else
      dendrogram <- "none"
  }
  
  ## dendrogram - check consistency: Colv ##
  if ( is.null(ddc) & (dendrogram %in% c("both","column"))){
    warning("Discrepancy: Colv is invalid or FALSE, while dendrogram is `",
            dendrogram,"'. Omitting column dendogram.")
    if (is.logical(Rowv) & (identical(Rowv.ori,TRUE) | is.numeric(Rowv.ori) | inherits(Rowv.ori,"dendrogram")) & dendrogram=="both")
      dendrogram <- "row"
    else
      dendrogram <- "none"
  }
  
  ## check consistency
  if (is.null(ddr)){
    if(isTRUE(cluster.by.row) | isTRUE(plot.row.partition) | isTRUE(plot.row.clusters) | isTRUE(plot.row.clustering) ){
      warning("Using invalid `Rowv' while allowing",
              "`cluster.by.row' or `plot.row.partition' or `plot.row.clusters' or `plot.row.clustering'",
              "can produce unpredictable results; Forced to be disabled.")
    }
  }
  
  if (is.null(ddc)){
    if(isTRUE(cluster.by.col) | isTRUE(plot.col.partition) | isTRUE(plot.col.clusters) | isTRUE(plot.col.clustering) ){
      warning("Using invalid `Colv' while allowing",
              "`cluster.by.col' or `plot.col.partition' or `plot.col.clusters' or `plot.col.clustering'",
              "can produce unpredictable results; Forced to be disabled.")
    }
  }

  if (is.null(ddr)) cluster.by.row <- plot.row.partition <- plot.row.clusters <- plot.row.clustering <- FALSE
  if (is.null(ddc)) cluster.by.col <- plot.col.partition <- plot.col.clusters <- plot.col.clustering <- FALSE

  ## ------------------------------------------------------------------------
  ## Reordering
  ## ------------------------------------------------------------------------
  flush.console()
  ## reorder x and cellnote ##
  x <- x[rowInd,colInd]
  
  if (!.invalid(cellnote)) cellnote <- cellnote[rowInd,colInd]
  
  ## reorder labels - row ##
  if(identical(labRow,TRUE)){ ## Note: x is already reorderred 
    labRow <- if (is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
  } else if(identical(labRow,FALSE) | .invalid(labRow)){
    labRow <- rep("",nrow(x))
  } else if(is.character(labRow)){
    labRow <- labRow[rowInd]
  }

  ## reorder cellnote/labels - col ##
  if (identical(labCol,TRUE)){
    labCol <- if(is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
  } else if(identical(labCol,FALSE) | .invalid(labCol)){
    labCol <- rep("",ncol(x))
  } else if(is.character(labCol)){
    labCol <- labCol[colInd]
  }
  
  ## ------------------------------------------------------------------------
  ## scale
  ## center to 0 and scale to 1 in row or col but not both! ##
  ## ------------------------------------------------------------------------
  flush.console()
  x <- .scale.data(x,scale,na.rm)

  ## ------------------------------------------------------------------------
  ## labels for observations/clusters/
  ## ------------------------------------------------------------------------
  ## margin for labels
  
  margin.for.labRow0 <- max(nchar(labRow))*0.75+0.2
  margin.for.labCol0 <- max(nchar(labCol))*0.75+0.2
  
  if (.invalid(margin.for.labRow)){
    margin.for.labRow <- margin.for.labRow0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labRow.')
    margin.for.labRow <- margin.for.labRow0*margin.for.labRow
  }
  
  if (.invalid(margin.for.labCol)){
    margin.for.labCol <- margin.for.labCol0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labCol.')
    margin.for.labCol <- margin.for.labCol0*margin.for.labCol    
  }
  
  ## group unique labels - row ## ##??check
  if (!.invalid(labRow.by.group) & !identical(labRow.by.group,FALSE)){
    group.value <- unique(labRow)
    group.index <- sapply(group.value,function(x,y) min(which(y==x)),y=labRow)
    labRow <- rep("",length(labRow))
    labRow[group.index] <- group.value
  }
  
  ## group unique labels - col ## ##??check
  if (!.invalid(labCol.by.group) & !identical(labCol.by.group,FALSE)){
    group.value <- unique(labCol)
    group.index <- sapply(group.value,function(x,y) min(which(y==x)),y=labCol)
    labCol <- rep("",length(labCol))
    labCol[group.index] <- group.value
  }
  
  ## ------------------------------------------------------------------------
  ## color breaks
  ## ------------------------------------------------------------------------
  flush.console()
  
  ## set breaks for binning x into colors ##
  if(.invalid(breaks)){
    breaks <- 16
  }
  
  ## get x.range according to the value of x.center ##
  if (!.invalid(x.center)){ ## enhanced
    if (is.numeric(x.center)){
      x.range.old <- range(x,na.rm=TRUE)
      dist.to.x.center <- max(abs(x.range.old-x.center))
      x.range <- c(x.center-dist.to.x.center,x.center+dist.to.x.center)
    } else {
      stop("`x.center' should be numeric.")
    } 
  } else{
    x.range <- range(x,na.rm=TRUE)
  }


  ## set breaks for centering colors to the value of x.center ##
  if(length(breaks)==1){
    breaks <-
      seq(min(min(x,na.rm=TRUE),x.range[1]),
          max(max(x,na.rm=TRUE),x.range[2]),
          length.out=breaks)
  }

  ## count of breaks and colors ##
  nbr <- length(breaks)
  ncolor <- length(breaks)-1

  ## generate colors ##
  colors <- color.FUN(ncolor)
  
  ## set up breaks and force values outside the range into the endmost bins ##
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[] <- ifelse(x<min.breaks,min.breaks,x)
  x[] <- ifelse(x>max.breaks,max.breaks,x)

  ## ------------------------------------------------------------------------
  ## check if it is sufficient to draw side plots ##
  ## ------------------------------------------------------------------------
  if (cluster.by.row){
    if (!.invalid(row.clusters)) {## suppress kr
      if(!is.numeric(row.clusters) | length(row.clusters)!=nr | !(.is.grouped(row.clusters))){
        warning("`row.clusters' is not a grouped numeric vector of length nrow(x); cluster.by.row is set to FALSE.")
        cluster.by.row <- FALSE
      } else{
        row.clusters <- row.clusters[rowInd]
        kr <- length(unique(row.clusters))
      }
    } else {
      if (.invalid(kr)) kr <- 2
      if (is.numeric(kr) & length(kr)==1){
        row.clusters <- cutree(hclust.row,k=kr)
        row.clusters <- row.clusters[rowInd]
      } else {
        warning("`kr' should be numeric of length one; cluster.by.row is set to FALSE.")
        cluster.by.row <- FALSE
      }
    }
  }

  if (cluster.by.col){
    if (!.invalid(col.clusters)) {## suppress kc
      if(!is.numeric(col.clusters) | length(col.clusters)!=nc | !(.is.grouped(col.clusters))){
        warning("`col.clusters' is not a grouped numeric vector of length ncol(x); cluster.by.col is set to FALSE.")
        cluster.by.col <- FALSE
      } else{
        col.clusters <- col.clusters[colInd]
        kc <- length(unique(col.clusters))
        if(revC){ ## x columns reversed
          col.clusters <- rev(col.clusters)
        }

      }
    } else {
      if (.invalid(kc)) kc <- 2
      if (is.numeric(kc) & length(kc)==1){
        col.clusters <- cutree(hclust.col,k=kc)
        col.clusters <- col.clusters[colInd]
        if(revC){ ## x columns reversed
          col.clusters <- rev(col.clusters)
        }

      } else {
        warning("`kc' should be numeric of length one; cluster.by.col is set to FALSE.")
        cluster.by.col <- FALSE
      }
    }
  }
  
  ## ------------------------------------------------------------------------
  ## Plotting
  ## ------------------------------------------------------------------------
  if (if.plot){

    ir <- length(plot.row.individuals.list)
    ic <- length(plot.col.individuals.list)
    cr <- length(plot.row.clustering.list)
    cc <- length(plot.col.clustering.list)

    flush.console()
    if(mapratio<=1){
      sr <- 12
      sc <- sr*mapratio
    } else {
      sc <- 12
      sr <- sc/mapratio
    }
    
    ## calculate the plot layout ##
    ## 1) for heatmap
    lmat <- matrix(1,nrow=sr,ncol=sc) 
    lwid <- c(rep(mapsize/sc,sc))
    lhei <- c(rep(mapsize/mapratio/sr,sr))

    ## 2) row.clusters
    if (plot.row.partition | plot.row.clusters){ 
      lmat <- cbind(max(lmat,na.rm=TRUE)+1,lmat) 
      lwid <- c(0.3,lwid) 
    } else {
      lmat <- cbind(NA,lmat)
      lwid <- c(0.02,lwid) 

    }
    
    ## 3) col.clusters
    if (plot.col.partition | plot.col.clusters){ 
      lmat <- rbind(c(NA,rep(max(lmat,na.rm=TRUE)+1,sc)),lmat) 
      lhei <- c(0.3/mapratio,lhei) 
    } else {
      lmat <- rbind(NA,lmat)
      lhei <- c(0.02/mapratio,lhei)
    }

    if(!.invalid(RowIndividualColors)) { ## 4) add middle column to layout for vertical sidebar ##??check
      if(!is.character(RowIndividualColors) || length(RowIndividualColors) !=nr)
        stop("'RowIndividualColors' must be a character vector of length nrow(x)")
      lmat <- cbind(c(rep(NA,nrow(lmat)-sr),rep(max(lmat,na.rm=TRUE)+1,sr)),lmat)
      lwid <- c(0.2,lwid) 
    } else {
      lmat <- cbind(NA,lmat)
      lwid <- c(0.02,lwid) 
    }
    
    if(!.invalid(ColIndividualColors)) { ## 5) add middle row to layout for horizontal sidebar ##??check
      if(!is.character(ColIndividualColors) || length(ColIndividualColors) !=nc){
        stop("'ColIndividualColors' must be a character vector of length ncol(x)")
      }
      lmat <- rbind(c(rep(NA,ncol(lmat)-sc),rep(max(lmat,na.rm=TRUE)+1,sc)),lmat) 
      lhei <- c(0.2/mapratio,lhei) 
    } else {
      lmat <- rbind(NA,lmat)
      lhei <- c(0.02/mapratio,lhei) 
    }

    ## 6) for row-dend
    lmat <- cbind(c(rep(NA,nrow(lmat)-sr),
                    rep(max(lmat,na.rm=TRUE)+1,sr)),
                  lmat
                  ) 
    lwid <- c(keysize,lwid)
    
    ## 7) for col-dend, 8) for kay
    lmat <- rbind(c(
                    max(lmat,na.rm=TRUE)+2,
                    rep(NA,ncol(lmat)-sc-1),
                    rep(max(lmat,na.rm=TRUE)+1,sc)
                    ),
                  lmat
                  )
    lhei <- c(keysize/mapratio,lhei)

    ## text.box##
    ## numbered 999 ##
    ## 9) for RowPlot (from bottom)
    if(.invalid(text.box)){
      text.box <- "made by\nFunction: heatmap.3\nPackage: GMD\nin R"
    }
    if(plot.row.individuals) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat,na.rm=TRUE)),nrow(lmat)-sr),# text
                      rep((ir:1)+max(lmat,na.rm=TRUE)+(1),each=floor(sr/ir)),rep(NA,sr%%ir)
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }
    
    ## 10) for ColPlot from right
    if(plot.col.individuals) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat,na.rm=TRUE)),ncol(lmat)-sc-1),# text
                      rep((1:ic)+max(lmat,na.rm=TRUE)+(1),each=floor(sc/ic)),rep(NA,sc%%ic),
                      ##NA # change to numeric if text.box
                      999
                      )
                    )
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }
    
    ## 11) for RowPlot (from bottom)
    if(plot.row.clusters) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),nrow(lmat)-sr-1), # text
                      rep((kr:1)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sr/kr)),rep(NA,sr%%kr),
                      ##NA
                      999
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }
    
    ## 12) for ColPlot from right 
    if(plot.col.clusters) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),ncol(lmat)-sc-2),# text
                      rep((1:kc)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sc/kc)),rep(NA,sc%%kc),
                      ##NA,NA # change to numeric if text.box
                      999,999
                      )
                    ) 
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }

    ## 13) for RowPlot (from bottom)
    if(plot.row.clustering) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),nrow(lmat)-sr-2), # text
                      rep(c((cr:1)+max(lmat[lmat!=999],na.rm=TRUE)+(1)),each=floor(sr/cr)),rep(NA,sr%%cr),
                      ##NA,NA
                      999,999
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }

    
    ## 14) for ColPlot from right
    if(plot.col.clustering) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),ncol(lmat)-sc-3),# text
                      rep((1:cc)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sc/cc)),rep(NA,sc%%cc),
                      ##NA,NA,NA # change to numeric if text.box
                      999,999,999
                      )
                    ) 
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }

    lmat[is.na(lmat)] <- 0
    if (any(lmat==999)) flag.text <- TRUE else flag.text <- FALSE
    lmat[lmat==999] <- max(lmat[lmat!=999])+1
    
    ## Graphics `output' ##
    ## layout
    if(!is.null(force_lmat)){
        lmat <- force_lmat
    }
    if(!is.null(force_lwid)){
        lwid <- force_lwid
    }
    if(!is.null(force_lhei)){
        lhei <- force_lhei
    }
    if(!force_add){
        layout(lmat,widths=lwid,heights=lhei,respect=FALSE)
    }

    ## reverse columns
    if(revC){ ## x columns reversed
      iy <- nr:1
      ddc <- rev(ddc)
      x <- x[iy,]
      if (!.invalid(cellnote)) cellnote <- cellnote[iy,]
    } else {
      iy <- 1:nr
    }

    ## reverse rows
    if(revR){ ## x columns reversed
      ix <- nc:1
      ddr <- rev(ddr)
      x <- x[,ix]
      if (!.invalid(cellnote)) cellnote <- cellnote[,ix]
    } else {
      ix <- 1:nc
    }
    
    ## 1) draw the main carpet/heatmap
    margins <- c(margin.for.labCol,0,0,margin.for.labRow)
    mgp <- c(3,1,0)
    par(mar=margins,mgp=mgp);outer=FALSE

    x.save <- x
    if(!symm || scale !="none"){ ##??
      x <- t(x)
      if (!.invalid(cellnote)) cellnote <- t(cellnote)
    }

    image(1:nc,1:nr,
          x,
          xlim=0.5+c(0,nc),ylim=0.5+c(0,nr),
          axes=FALSE,xlab="",ylab="",col=colors,breaks=breaks,
          ...)
    
    ## plot/color NAs
    if(!.invalid(na.color) & any(is.na(x))){
      mmat <- ifelse(is.na(x),1,NA)
      image(1:nc,1:nr,mmat,axes=FALSE,xlab="",ylab="",
            col=na.color,add=TRUE)
    }

    ##
    ## labCol (?)
    if ((dendrogram %in% c("both","col")) & sideCol==3) {
      warning("Discrepancy: col dendrogram is asked; srtCol is set to 1.")
      sideCol <- 1
    }
    if (!length(srtCol)) {
      axis(sideCol,1:nc,labels=labCol,las=2,line=-0.5,tick=0,cex.axis=cexCol,outer=outer)
    } else {
      if (sideCol==1){
        if (sideCol==1) .sideCol <- par("usr")[3]-0.5*srtCol/90 else .sideCol <- par("usr")[4]+0.5*srtCol/90
        text(1:nc,.sideCol,labels=labCol,srt=srtCol,pos=1,xpd=TRUE,cex=cexCol)
      }
    }
    
    if(!.invalid(xlab)) mtext(xlab,side=1,line=margins[1]-1.25)

    ## labRow (?)
    if ((dendrogram %in% c("both","row")) & sideRow==2) {
      warning("Discrepancy: row dendrogram is asked; sideRow is set to 4.")
      sideRow <- 4
    }
    if (!length(srtRow)) {
      axis(sideRow,iy,labels=labRow,las=2,line=-0.5,tick=0,cex.axis=cexRow,outer=outer)
    } else {
      if (sideRow==4){
        if (sideRow==4) .sideRow <- par("usr")[2]+0.5*srtRow/90 else .sideRow <- par("usr")[1]-0.5*srtRow/90
        text(.sideRow,iy,labels=labRow,srt=srtRow,pos=1,xpd=TRUE,cex=cexRow)
      }
    }
    
    if(!.invalid(ylab)) mtext(ylab,side=4,line=margins[4]-1.25)

    if (!.invalid(add.expr))
      eval(substitute(add.expr))
    
    ## Enhanced: add 'sep.color' colored spaces to visually separate sections
    if (plot.row.partition | plot.row.clusters){ ##??
      plot.row.partitionList <- get.sep(clusters=row.clusters,type="row")
    } else {
      plot.row.partitionList <- NULL
    }
    if (plot.col.partition | plot.col.clusters){ ##??
      plot.col.partitionList <- get.sep(clusters=col.clusters,type="column")
    } else {
      plot.col.partitionList <- NULL
    }

    row.sepList <- sepList[[1]]
    if (!.invalid(row.sepList)){
      for (i in 1:length(row.sepList)){
        i.sep <- row.sepList[[i]]
        rect(
             xleft=i.sep[1]+0.5,
             ybottom=i.sep[2]+0.5,
             xright=i.sep[3]+0.5,
             ytop=i.sep[4]+0.5,
             lty=sep.lty,
             lwd=sep.lwd,
             col=FALSE,
             border=sep.color[1]
             )
      }
    }

    col.sepList <- sepList[[2]]
    if (!.invalid(col.sepList)){
      for (i in 1:length(col.sepList)){
        i.sep <- col.sepList[[i]]
        rect(
             xleft=i.sep[1]+0.5,
             ybottom=i.sep[2]+0.5,
             xright=i.sep[3]+0.5,
             ytop=i.sep[4]+0.5,
             lty=sep.lty,
             lwd=sep.lwd,
             col=FALSE,
             border=sep.color[2]
             )
      }
    }
        
    ## show traces
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    
    x.scaled  <- .scale.x(t(x),min.scale,max.scale)
  
    if(.invalid(hline)) hline=median(breaks)
    if(.invalid(vline)) vline=median(breaks)
    
    if(trace %in% c("both","column")){
      for( i in colInd ){
        if(!.invalid(vline)){
          vline.vals <- .scale.x(vline,min.scale,max.scale)
          abline(v=i-0.5+vline.vals,col=linecol,lty=2)
        }
        xv <- rep(i,nrow(x.scaled))+x.scaled[,i]-0.5
        xv <- c(xv[1],xv)
        yv <- 1:length(xv)-0.5
        lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
      }
    }

    if(trace %in% c("both","row")){
      for( i in rowInd ){
        if(!.invalid(hline)){
          hline.vals <- .scale.x(hline,min.scale,max.scale)
          abline(h=i+hline,col=linecol,lty=2)
        }
        yv <- rep(i,ncol(x.scaled))+x.scaled[i,]-0.5
        yv <- rev(c(yv[1],yv))
        xv <- length(yv):1-0.5
        lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
      }
    }

    ## cellnote
    if(!.invalid(cellnote)){
      text(x=c(row(cellnote)),
           y=c(col(cellnote)),
           labels=c(cellnote),
           col=notecol,
           cex=cex.note)
    }
    
    ## 2) plot.row.partition
    if(plot.row.partition |plot.row.clusters) { ##row.clusters
      par(mar=c(margins[1],0.5,0,0.1))

      row.clusters.unique <- unique(row.clusters)
      row.clusters.unique <- row.clusters.unique[!is.na(row.clusters.unique)]
      
      image(rbind(1:nr),
            xlim=0.5+c(0,1),ylim=0.5+c(0,nr),
            col=par("bg"),
            axes=FALSE,add=force_add)

      if (!.invalid(plot.row.partitionList)){
        for (i in 1:length(plot.row.partitionList)){
          i.sep <- plot.row.partitionList[[i]]
          rect(
               xleft=0+0.5,
               ybottom=i.sep[2]+0.5,
               xright=1+0.5,
               ytop=i.sep[4]+0.5,
               lty=sep.lty,
               lwd=sep.lwd,
               col=color.partition.box,
               border=color.partition.border
               )
          g <- row.clusters.unique[i]
          s <- g
          text(x=1,y=(i.sep[2]+0.5+i.sep[4]+0.5)/2,labels=s,col=color.partition.border,
               cex=cex.partition,
               srt=90
               )
        }
      }
      
    } 
    
    ## 3) plot.col.partition
    if(plot.col.partition | plot.col.clusters) {
      par(mar=c(0.1,0,0,margins[4]))
      col.clusters.unique <- unique(col.clusters)
      col.clusters.unique <- col.clusters.unique[!is.na(col.clusters.unique)]
      
      image(cbind(1:nc),
            xlim=0.5+c(0,nc),ylim=0.5+c(0,1),
            col=par("bg"),
            axes=FALSE,add=force_add)
      
      if (!.invalid(plot.col.partitionList)){
        for (i in 1:length(plot.col.partitionList)){
          i.sep <- plot.col.partitionList[[i]]
          rect(
               xleft=i.sep[1]+0.5,
               ybottom=0+0.5,
               xright=i.sep[3]+0.5,
               ytop=1+0.5,
               lty=sep.lty,
               lwd=sep.lwd,
               col=color.partition.box,
               border=color.partition.border
               )
          g <- col.clusters.unique[i]
          s <- g
          text(x=(i.sep[1]+0.5+i.sep[3]+0.5)/2,y=1,labels=s,col=color.partition.border,
               cex=cex.partition,
               srt=0
               )
        }
      }

    }
    
    ## 4) draw the side color bars - for row
    if(!.invalid(RowIndividualColors)) {    
      par(mar=c(margins[1],0,0,0.5))
      image(rbind(1:nr),col=RowIndividualColors[rowInd],axes=FALSE,add=force_add)
    } 
    
    ## 5) draw the side color bars - for col
    if(!.invalid(ColIndividualColors)) {
      par(mar=c(0.5,0,0,margins[4]))
      image(cbind(1:nc),col=ColIndividualColors[colInd],axes=FALSE,add=force_add)
    }

    ## 6) row-dend
    par(mar=c(margins[1],0,0,0))
    if(dendrogram %in% c("both","row")){
      plot(ddr,horiz=TRUE,axes=FALSE,yaxs="i",leaflab="none")
    }else{
      .plot.text(ylim=range(iy))
      if (sideRow==2){
        .sideRow <- par("usr")[2]-0.5*srtCol/90
        text(.sideRow,iy,labels=labRow,srt=srtRow,pos=1,xpd=TRUE,cex=cexRow)        
      }
    }

    ## 7) col-dend and title
    mar3 <- (if(!is.null(main)) mar.main else 0) +
      (if(!is.null(sub)) mar.sub else 0)
    par(mar=c(0,0,mar3,margins[4]))
    
    if(dendrogram %in% c("both","column")){
      plot(ddc,axes=FALSE,xaxs="i",leaflab="none")
    } else{
      .plot.text(xlim=range(1:nc))
      if (sideCol==3){
        .sideCol <- par("usr")[3]+0.5*srtCol/90
        text(1:nc,.sideCol,labels=labCol,srt=srtCol,pos=1,xpd=TRUE,cex=cexCol)
      }
    }

    ## title
    if (is.null(sub)) main.line <- 1 else main.line <- 3
    if(!is.null(main)) title(main,cex.main=cex.main,adj=adj.main,mgp=mgp.main,font.main=font.main,line=main.line)
    if(!is.null(sub)) title(sub,cex.main=cex.sub,adj=adj.main,mgp=mgp.main,font.main=font.sub,line=0)
    ##if(!is.null(main)) title(main,cex.main=1.5*op[["cex.main"]])

    ## 8) plot the color-key
    if(key){
      cex.key <- 0.75
      op.ori <- par()

      par(mar=c(2,1.5,0.75,1)*keysize,cex=cex.key,mgp=c(0.75,0,0),tcl=-0.05)
      z <- seq(x.range[1],x.range[2],length=length(colors))
      
      image(z=matrix(z,ncol=1),
            col=colors,
            breaks=breaks,
            xaxt="n",
            yaxt="n",
            xlab=key.xlab,
            ylab="",
            main="",add=force_add
            )
      par(usr=c(0,1,0,1))
      lv <- pretty(breaks)
      xv <- .scale.x(as.numeric(lv),x.range[1],x.range[2])
      axis(1,at=xv,labels=lv,cex.axis=cex.key*1)
      
      if(density.info=="density"){
        ## Experimental : also plot density of data
        dens <- density(x,adjust=densadj,na.rm=TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- .scale.x(dens$x,x.range[1],x.range[2])
        lines(dens$x,dens$y / max(dens$y) * 0.95,col=denscol,lwd=1)
        axis(2,at=pretty(dens$y)/max(dens$y) * 0.95,pretty(dens$y),cex.axis=cex.key*1)
        ##title("Color Key and Density",cex.lab=cex.key*0.25)
        title(key.title,cex.main=cex.key,font.main=1)
        mtext(side=2,"Density",line=0.75,cex=cex.key)
      } else if(density.info=="histogram"){
        h <- hist(x,plot=FALSE,breaks=breaks)
        hx <- .scale.x(breaks,x.range[1],x.range[2])
        hy <- c(h$counts,h$counts[length(h$counts)])
        lines(hx,hy/max(hy)*0.95,lwd=1,type="s",col=denscol)
        axis(2,at=pretty(hy)/max(hy)*0.95,pretty(hy),cex.axis=cex.key*1)
        ##title("Color Key and Histogram",cex.main=cex.key*0.25)
        title(key.title,cex.main=cex.key,font.main=1)
        mtext(side=2,key.ylab,line=0.75,cex=cex.key)
      } else{
        title(key.title,cex.main=cex.key,font.main=1)
      }
    } else{
      if(!force_add){
      .plot.text()
      }
    }

    ## 9)
    if(plot.row.individuals) {
      .plot.text("Row\nIndividuals",cex=cex.text,bg="white")
      for (i in 1:ir) {
        ##.plot.text()
        tmp <- plot.row.individuals.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 10)
    if(plot.col.individuals) {
      .plot.text("Column\nIndividuals",cex=cex.text,bg="white",srt=90)
      for (i in 1:ic) {
        ##.plot.text()
        tmp <- plot.col.individuals.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }
    
    ## 11) for RowPlot from bottom
    if (plot.row.clusters){
      .plot.text("Row\nClusters",cex=cex.text,bg="white")
      
      tmp <- plot.row.clusters.list[[1]]
      row.data <- row.data[rowInd]
      for (i in unique(row.clusters)){
        i.x <- row.data[row.clusters==i]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
        i.main <- sprintf("Row group %s (n=%s)",i,length(i.x))
        title(i.main,cex.main=1,font.main=1)
      }
    }
    
    ## 12) for ColPlot from left
    if (plot.col.clusters){
      .plot.text("Col\nClusters",cex=cex.text,bg="white",srt=90)
      
      tmp <- plot.col.clusters.list[[1]]
      col.data <- if(revC) col.data[rev(colInd)] else col.data[colInd]
      for (i in unique(col.clusters)){
        i.x <- col.data[col.clusters==i]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
        i.main <- sprintf("Col group %s (n=%s)",i,length(i.x))
        title(i.main,cex.main=1,font.main=1)
      }
    }

    ## 13)
    if(plot.row.clustering) {
      .plot.text("Row\nClustering",cex=cex.text,bg="white")
      
      for (i in 1:cr) {
        ##.plot.text()
        tmp <- plot.row.clustering.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 14)
    if(plot.col.clustering) {
      .plot.text("Column\nClustering",cex=cex.text,bg="white",srt=90)
      
      for (i in 1:cc) {
        ##.plot.text()
        tmp <- plot.col.clustering.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }
    
    ## 15) text
    if (!.invalid(text.box) & if.plot.info){
      .plot.text(text.box,cex=cex.text,bg="gray75")
    } else{
      if (flag.text){
        .plot.text()
      }
    }
 
  }
  
  ret <-
    list(x.ori=x.ori,
         x=x.save,
         rowInd=rowInd,colInd=colInd,
         row.clusters=row.clusters,col.clusters=col.clusters,
         dist.row=dist.row,dist.col=dist.col,
         hclust.row=hclust.row,hclust.col=hclust.col,
         kr=kr,kc=kc
         )
  class(ret) <- c("hclustering",class(ret))
  invisible(ret)
}

##' Please note this code is from the library GMD
##' All credit for this code goes to GMD's author's.
##' I do not recommend using this version of the code, which
##' has been poorly modified for our use but recommend using
##' the official version from the package GMD
##' https://cran.r-project.org/web/packages/GMD/index.html
##' Get row or column lines of separation for \code{heatmap.3} according to clusters
##' @title Get row or column lines of separation for heatmap.3
##' @param clusters a numerical vector, indicating the cluster labels of observations.
##' @param type string, one of the following: \code{c("row","column","both")}
get.sep <-
  function(clusters,type=c("row","column","both"))
{
  ##   if(!is.numeric(clusters) | !(.is.grouped(clusters))){
  ## stop("`clusters' should be a grouped numeric vector.")
  ##   }
  tmp.whichna <- which(is.na(clusters))
  tmp.which <- which(!duplicated(clusters))
  
  tmp.sep <- data.frame(start=tmp.which,end=c(tmp.which[-1],length(clusters)+1)-1)
  tmp.sep2 <- tmp.sep[tmp.sep$start<=tmp.sep$end,]

  ## lines of separation 
  sepList <- list()
  if (type=="row"){
    xmax <- length(clusters)
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(0,tmp.sep2[i.s,1]-1,xmax,tmp.sep2[i.s,2])
    }
  } else if (type=="column"){
    ymax <- length(clusters)
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(tmp.sep2[i.s,1]-1,0,tmp.sep2[i.s,2],ymax)
    }
  } else if (type=="both"){
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(tmp.sep2[i.s,1]-1,tmp.sep2[i.s,1]-1,tmp.sep2[i.s,2],tmp.sep2[i.s,2])
    }
  }
  sepList
}


