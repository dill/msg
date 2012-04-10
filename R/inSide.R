# MVB's inSide
#"in.poly" <-function( xy, poly) {
#inSide.mvb<-function( poly, x,y) {
#  xy<-cbind(x,y)
#  poly<-matrix(c(poly$x[-(length(poly$x))],poly$y[-(length(poly$x))]),length(poly$x)-1,2)
#  if( ncol( poly)==2)
#    poly <- poly.for.testing( poly)
#  n <- nrow( xy); np <- nrow( poly); nnp <- rep( n,np)
#  #check1 <- xor( xy[,1]>=rep( poly[,1], nnp), xy[,1]>rep( poly[,3], nnp))
#  #check2 <- rep( poly[,2], nnp) + rep( poly[,4], nnp) * xy[,1] > xy[,2]
#  #a<-as.vector( rowSums( matrix( check1 & check2, n, np)) %% 2 > 0)
#  #check1 <- xor( xy[,1]>rep( poly[,1], nnp), xy[,1]>rep( poly[,3], nnp))
#  #check2 <- rep( poly[,2], nnp) + rep( poly[,4], nnp) * xy[,1] >= xy[,2]
#  #a<-as.vector( rowSums( matrix( check1 & check2, n, np)) %% 2 > 0)
#  #check1 <- xor( xy[,1]>=rep( poly[,1], nnp), xy[,1]>=rep( poly[,3], nnp))
#  #check2 <- rep( poly[,2], nnp) + rep( poly[,4], nnp) * xy[,1] > xy[,2]
#  #b<-as.vector( rowSums( matrix( check1 & check2, n, np)) %% 2 > 0)
#  #a|b
#  check1 <- xor( xy[,1]>=rep( poly[,1], nnp), xy[,1]>rep( poly[,3], nnp))
#  check2 <- rep( poly[,2], nnp) + rep( poly[,4], nnp) * xy[,1] >= xy[,2]
#  a<-as.vector( rowSums( matrix( check1 & check2, n, np)) %% 2 > 0)
#  check1 <- xor( xy[,1]>rep( poly[,1], nnp), xy[,1]>rep( poly[,3], nnp))
#  check2 <- rep( poly[,2], nnp) + rep( poly[,4], nnp) * xy[,1] > xy[,2]
#  b<-as.vector( rowSums( matrix( check1 & check2, n, np)) %% 2 > 0)
#  a|b
#}
#in.poly1<-function( poly, x,y){
#  if( ncol( poly)==2)
#    poly <- poly.for.testing( poly)
#  n <- nrow( xy); np <- nrow( poly); nnp <- rep( n,np)
#  check1 <- xor( xy[,1]>=rep( poly[,1], nnp), xy[,1]>=rep( poly[,3], nnp))
#  check2 <- rep( poly[,2], nnp) + rep( poly[,4], nnp) * xy[,1] > xy[,2]
#  as.vector( matrix( check1 & check2, n, np) %*% rep.int( 1, np) %% 2 > 0)
#}
#
#in.poly<-function( poly, x,y){
#   ind1<-in.poly1(poly,x,y)
#   ind2<-in.poly1(list(x=-bnd$x,y=-bnd$y),-x,-y)
#
#   ind1|ind2
#}
#
#"poly.for.testing" <-function( xy) {
#  if( is.list( xy))
#    xy <- as.data.frame(xy)
#  if( is.data.frame( xy))
#    xy <- as.matrix( xy)
#  if( !(is.matrix( xy) || is.data.frame( xy)) || ncol( xy)!=2 || !is.numeric( xy[2,1]) ||
#      !is.numeric( xy[1,2]))
#stop( "xy must by nX2 numeric matrix or data.frame")
#
## Columns of poly are start.x, a, end.x, b
#
## Duplicate points at end of polygons. No problem if done already,
## as will be dropped in last line.
#  na.rows <- c( 0, index( is.na( xy[,1])), nrow( xy)+1)
#  xy <- rbind( xy, xy[1,])
#  xy[ na.rows[-1],] <- xy[ na.rows[ -length( na.rows)]+1,]
#  poly <- matrix( c( xy, xy[ c( 2:nrow( xy), 1), ]), ncol=4)
#
#  poly[,4] <- (poly[,4]-poly[,2])/(poly[,3]-poly[,1])
#  poly[,2] <- poly[,2] - poly[,1]*poly[,4]
#  poly <- poly[ -na.rows[-1],]
#  poly[ poly[,1] != poly[,3], ]
#}
#
#"index" <-
#function (lvector)
#seq_along( lvector)[lvector]

#' Are points inside a polygon?
#'
#' This is Simon Wood's inSide from soap
#'
#' @param bnd a polygon boundary (a list with elements x and y)
#' @param x x coordinates of the points to test
#' @param y y coordinates of the points to test
#' @return logical vector
#' @export
# http://www.maths.bath.ac.uk/~sw283/simon/software.html
inSide <- function(bnd,x,y)
## tests whether each point x[i],y[i] is inside the boundary defined
## by bnd$x, bnd$y, or by multiple boundary loops in bnd[[1]]$x,
## bnd[[1]]$y, bnd[[2]]$x, ... etc.
{ if (is.null(bnd$x)) { ## translate into form that C routine needs
    bn <- list(x=bnd[[1]]$x,y=bnd[[1]]$y)
    if (length(bnd)>1) for (i in 2:length(bnd)) {
      bn$x <- c(bn$x,NA,bnd[[i]]$x);bn$y <- c(bn$y,NA,bnd[[i]]$y)
    }
    bnd <- bn
  }
  ## replace NA segment separators with a numeric code 
  lowLim <- min(c(bnd$x,bnd$y),na.rm=TRUE)-1
  ind <- is.na(bnd$x)|is.na(bnd$y)
  bnd$x[ind] <- bnd$y[ind] <- lowLim - 1
  n <- length(bnd$x)
  if (n != length(bnd$y)) stop("x and y must be same length")
  um <-.C("in_out1",bx=as.double(bnd$x),by=as.double(bnd$y),break.code=as.double(lowLim),
          x=as.double(x),y=as.double(y),inside=as.integer(y*0),nb=as.integer(n),
          n=as.integer(length(x)),PACKAGE="msg")
  as.logical(um$inside)
}

