# function to create a grid for soap
make_soap_grid<-function(bnd,n.grid,mat=FALSE,log=FALSE,delta=FALSE){

   # set the grid size, if the input is a 2-vec then it is m and n
   if(length(n.grid)==2){
      m<-n.grid[1]
      n<-n.grid[2]
   }else{
      m<-n<-n.grid
   }

   # min and max values of the boundary (but not on the boundary)
   xmin<-min(bnd$x,na.rm=TRUE)
   ymin<-min(bnd$y,na.rm=TRUE)
   xmax<-max(bnd$x,na.rm=TRUE)
   ymax<-max(bnd$y,na.rm=TRUE)

   # create the grid
   xm <- seq(xmin,xmax,length=m)
   yn<-seq(ymin,ymax,length=n)
   xx <- rep(xm,n)
   yy<-rep(yn,rep(m,n))

   onoff<-inSide(bnd,xx,yy)
   xx<-xx[onoff]
   yy<-yy[onoff]

   ret<-list(x=xx,y=yy)

   # if we want an image plot, return a matrix
   if(mat){
     mat<-matrix(NA,m,n)
     mat[onoff]<-0
     ret$mat<-mat
   }
   
   # return the logical for the grid
   if(log){
      ret$log<-onoff
   }

   # return the deltas
   if(delta){
      ret$deltax<-abs(xm[1]-xm[2])
      ret$deltay<-abs(yn[1]-yn[2])
   }

   return(ret)

}

