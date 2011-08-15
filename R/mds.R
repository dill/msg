# create the distance matrix
create_distance_matrix<-function(xpoints,ypoints,bnd,start=NA,faster=0,debug=0){
   # requires the library soap   
   # args:
   #  xpoints,ypoints      data points
   #  bnd                  boundary list(x=c(),y=c())
   #  res                  resolution to test if line is inside

   # make sure that x and y are the same length
   if(length(xpoints)!=length(ypoints)){
      cat("ERROR: Vector lengths not the same!\n")
      return(FALSE)
   }
   if(any(is.na(xpoints))|any(is.na(ypoints))){
      cat("ERROR: Some elements are NA!\n")
      return(FALSE)
   }

   # create a matrix to hold the distances

   D<-woodpath(xpoints,ypoints,bnd,start,faster,debug)

   return(D)
}
