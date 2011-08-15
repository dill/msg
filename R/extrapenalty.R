# make some extra penalties for msg...
# want to add \int (f')^2 for each dimension 
extra.penalty<-function(object,data){

   # object has the Duchon object already, just need to 
   # add new elements to object$S

   eps<-(1e-15)^(1/4)

   k<-dim(object$S[[1]])[1]

   N<-100
   N<-20 #### ARBITRARY!!!!!!!

   if(!is.null(object$xt$int.res)){
      N<-object$xt$int.res
   }

   ### first need to create the mesh we want to integrate over
   # mesh function
   mesh <- function(x,d,w=1/length(x)+x*0) { 
      n <- length(x) 
      W <- X <- matrix(0,n^d,d) 
      for (i in 1:d) {
         X[,i] <- x;W[,i] <- w
         x<- rep(x,rep(n,length(x))) 
         w <- rep(w,rep(n,length(w)))
      } 
      w <- exp(rowSums(log(W))) ## column product of W gives weights 
      list(X=X,w=w) ## each row of X gives co-ordinates of a node
   }

   # extract the grid point locations in MDS space
   grid.points<-object$msg$mds.obj$points

   # find the extent of that
   grid.max<-max(grid.points)
   grid.min<-min(grid.points)

   # make a mesh to integrate over in the mds space
   ip <- mesh(grid.min+(1:N-.5)/N*(grid.max-grid.min),
              ncol(grid.points),
              rep(2/N,N))

   # this is going to be waaaaaaay to big, so trim it.


   # first just if the extent in each dimension is 
   # greater than the extent of the grid in that dimension

   # in each dimension
   for(d in 1:ncol(grid.points)){
      # find the max and min
      d.max<-max(grid.points[,d])
      d.min<-min(grid.points[,d])

      # if any points in those dimensions fall 
      # outside of the max and min
      ind<-ip$X[,d]>d.max | ip$X[,d]<d.min
      # remove the points and their corresponding weights
      ip$X<-ip$X[!ind,]
      ip$w<-ip$w[!ind]
   }

   # next step is to do something clever!

   # find the distances from even grid point to every mesh point,
   # put them in D
   D<-as.matrix(dist(rbind(ip$X,grid.points),diag=T,upper=T))
   D<-D[1:nrow(ip$X),]
   D<-D[,(nrow(ip$X)+1):ncol(D)]

   # want the minimum distance each mesh point has to any grid point
   rmin<-apply(D,1,min)

   ipd<-2*diff(grid.min+(1:N-.5)/N*(grid.max-grid.min))[1]
   too.far<-which(rmin>ipd)
   ind<-unique(too.far)

   # now we have a nicer mesh to integrate over!
   ip$X<-ip$X[-ind,]
   ip$w<-ip$w[-ind]

   # plot the integration grid
   #plot(ip$X,pch=19,cex=0.3,asp=1)
   #points(grid.points,pch=19,cex=0.3,col="red")

   # root the weights, since we square them in a bit
   ip$w<-sqrt(ip$w)


   # for prediction give make the mesh a data frame and give it
   # the right column names
   pmat<-as.data.frame(ip$X)
   names(pmat)<-paste("mds-",1:dim(pmat)[2],sep="")

   offset<-length(object$S)

   # for each dimension we're dealing with
   for(d in 1:ncol(grid.points)){
   
      # finite first differences

      # perturb forward
      pmatf<-pmat
      pmatf[,d]<-pmatf[,d]+2*eps
      ffd<-Predict.matrix(object,pmatf)
                                 
      # perturb backward
      pmatb<-pmat
      pmatb[,d]<-pmatb[,d]-2*eps
      bfd<-Predict.matrix(object,pmatb)

      # find the finite difference
      fd<-ip$w*(ffd-bfd)/(2*eps)

      # integrate!
      fd<-t(fd)%*%fd
      # enforce symmetry (from smooth.construct.tp...)
      fd <- (fd + t(fd))/2

      object$S[[offset+d]]<-fd

      object$rank<-c(object$rank,sum(eigen(fd)$values>1e-14))

      # clean up
      rm(pmatf); rm(pmatb); rm(fd); gc()

   }

   object
}
