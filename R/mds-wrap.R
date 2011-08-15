# wrapper function for C code
woodpath<-function(xpoints,ypoints,bnd,start=NA,faster=0,debug=0){

   # put everything in the right format
   xbnd<-bnd$x
   ybnd<-bnd$y
   nbnd<-length(xbnd)
   len<-length(xpoints)

   if(is.na(start)){
      # when we do full MDS then we just compute the upper triangle
      pl<-(len*(len-1)/2)
      insert<-FALSE # used later
      start<-0
   }else{
      # in the insertion case, start gives the length of the old points
      # here we are calculating a length(old.points)*length(new.points) matrix
      pl<-len-start # length of new points
      pl<-pl*start # size of matrix
      insert<-TRUE # used later
   }

   # set up the reference grid in R, since it's much easier
   ref.grid<-create_refgrid(bnd)

   # load the library
   #dyn.load("wood.so")
   library.dynam("mdspack",package=c("mdspack"))

   ## code for running everything at once...
   wood_ret<-.C("wood_path",len=as.integer(len),start=as.integer(start), 
                x=as.double(xpoints),y=as.double(ypoints),
                nbnd=as.integer(nbnd),
                xbnd=as.double(xbnd),ybnd=as.double(ybnd),
                xref=as.double(ref.grid$x),yref=as.double(ref.grid$y),
                ngrid=as.integer(ref.grid$nrefx),
                refdelx=as.double(ref.grid$deltax),
                refdely=as.double(ref.grid$deltay),
                refio=as.integer(ref.grid$log),
                nref=as.integer(length(ref.grid$x)),
                xstart=as.double(ref.grid$x[1]),
                ystart=as.double(ref.grid$y[1]),
                pathlen=as.double(rep(0,pl)),
                faster=as.integer(faster),
                debug=as.integer(debug))

   # full MDS
   if(!insert){
      # get passed back an array which is the upper diagonal
      # create a matrix
      D<-matrix(0,len,len)
      # R fills columns first, so fill the lower triangle first
      # then take the transpose for the same effect   
      D[lower.tri(D)]<-wood_ret$pathlen
      D<-t(D)

      D<-D+t(D)

   # insertion
   }else{
      D<-matrix(wood_ret$pathlen,ncol=start,nrow=len-start)
      D<-t(D)
   }

   # unload the library
   if(R.version$os=="darwin9.8.0"){
      extraslash<-"/"
   }else{
      extraslash<-""
   }

   library.dynam.unload("mdspack",paste(.libPaths()[1],extraslash,"mdspack",sep=""))

   return(D)
}
