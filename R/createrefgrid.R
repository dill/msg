create_refgrid<-function(bnd,dens=25){
   # create reference grid for the path finding

   # starting value
   res<-10


   if(length(dens)==2){

      grid<-make_soap_grid(bnd,dens,log=TRUE,delta=TRUE)

   }else{

      # using a square grid: make_soap_grid   
      grid<-make_soap_grid(bnd,res,log=TRUE,delta=TRUE)
      
      while(length(grid$x)<dens){
         res<-res+1
         grid<-make_soap_grid(bnd,res)
      }

   }

   grid$nrefx<-res
   grid$nrefy<-res

   return(grid)
}
