smooth.construct.msg.smooth.spec<-function(object,data,knots){

   # this is the msg smooth.spec file
   # this does most of the work

   # for now TODO:
   #  just get spatial smoothing working
   #     do the clever stuff with previous objects as in gam.mds.R
   # THEN:
   #  if no bnd supplied do the general thing
   #  also do the s(.) thing
   #  select grid resolution sensibly

   # extract the boundary
   bnd<-object$xt$bnd

   # for the WAD case 
   if(!is.null(bnd)){
      if(length(names(data))!=2){
         stop("msg can only be used with 2D smooths!\n")
      }
      if(any(names(data)!=c("x","y"))){
         stop("Names of data are not \"x\" and \"y\"\n")
      }
   }

   if(is.null(object$xt$mds.dim)){
      stop("No MDS projection dimension supplied!\n")
   }

   # extract the MDS dimension
   mds.dim<-object$xt$mds.dim

   grid.res<-object$xt$mds.grid.res
   if(is.null(grid.res)){
      # pick a grid size...

      # just pick something for now
      #grid.res<-c(40,40)
      grid.res<-120
   }

   # if there was an old object in the extra stuff, use it
   #old.obj<-object$xt$old.obj
   old.obj<-NULL

   if(!is.null(old.obj)){
      # object to store all the results for later
      new.obj<-old.obj

      # pull out the grid D matrix
      D.grid<-old.obj$D
      my.grid<-old.obj$grid

      # also the pred and sample D matrices if
      # they are there
      if(!is.null(old.obj$D.samp)){
         D.samp<-old.obj$D.samp
      }else{
         D.samp<-NULL
      }
      if(!is.null(old.obj$D.pred)){
         D.pred<-old.obj$D.pred
      }else{
         D.pred<-NULL
      }

      if(!is.null(old.obj$m)){
         m<-old.obj$m
      }
      if(!is.null(old.obj$bs)){
         bs<-old.obj$bs
      }
      if(!is.null(old.obj$k)){
         k<-old.obj$k
      }

      if(!is.null(old.obj$mds.dim)){
         mds.dim<-old.obj$mds.dim
      }

   }else{
      # object to store all the results for later
      new.obj<-list()

      # in the WAD case, need to create a grid but not in the 
      # GDS case...
      if(!is.null(bnd)){
         # create the grid
         grid.obj<-create_refgrid(bnd,grid.res)
         D.grid<-create_distance_matrix(grid.obj$x,grid.obj$y,bnd,faster=0)
         grid.obj<-list(D=D.grid,grid=list(x=grid.obj$x,y=grid.obj$y))

         D.grid<-grid.obj$D
         my.grid<-grid.obj$grid
         # store!
         new.obj$D<-D.grid
         new.obj$grid<-my.grid

         object$msg<-new.obj
      }
   }

   # for the WAD case, insert the points into the grid
   if(!is.null(bnd)){
      # now we have the grid object, insert the data into that
      # and store it as the data
      grid.mds<-cmdscale(D.grid,eig=TRUE,k=mds.dim,x.ret=TRUE)
      object$msg$mds.obj<-grid.mds
      mds.data<-as.data.frame(insert.mds.generic(grid.mds,data,my.grid,bnd=bnd))

      object$msg$metric<-"WAD"

   # in the GDS case
   }else{

      new.obj<-list()

      # choose the distance metric
      if(is.null(object$xt$metric)){
         object$xt$metric<-"euclidean"
      }

      # the data in matrix form
      #data.names<-names(data)
      mdata<-as.matrix(as.data.frame(data))

      ## separate the predictors from response
      #ind<-rep(FALSE,ncol(mdata))
      #ind[match(object$term,data.names)]<-TRUE

      ## save the response
      #response.var<-mdata[,!ind]
      #mdata<-mdata[,ind]

      new.obj$metric<-object$xt$metric
      dist.metric<-object$xt$metric

      if(dist.metric=="mahalanobis"){
         # find the inverse covariance matrix
         cov.mat<-try(solve(cov(mdata)))

         # if something bad happened use the pseudo-inverse
         if(class(cov.mat)=="try-error"){
            cov.mat<-ginv(cov(mdata))
         }

         D<-apply(mdata,1,mahalanobis,x=mdata,cov=cov.mat,inverted=TRUE)
      }else{
      # Euclidean or anything else allowed by dist()
         D<-as.matrix(dist(mdata,method=dist.metric,diag=T,upper=T))
      }
      new.obj$D<-D

      mds.obj<-cmdscale(D,eig=TRUE,k=mds.dim,x.ret=TRUE)
      new.obj$mds.obj<-mds.obj
      mds.data<-as.data.frame(mds.obj$points)

      object$msg<-new.obj
   }

   # make some variable names up
   mds.names<-paste("mds-",1:dim(mds.data)[2],sep="")
   # remove any already in the data
   names(mds.data)<-mds.names

   # make sure there are the right stuff is in the object before passing
   # to Duchon, but save beforehand!
   save.dim<-object$dim
   save.term<-object$term
   save.data<-data

   object$term<-mds.names
   object$dim<-mds.dim
   data<-mds.data

   object$msg$term<-mds.names
   object$msg$dim<-mds.dim
   object$msg$data<-mds.data

   # if knots were supplied, they're going to be ignored, warn about that!
   if(length(knots)!=0){
      warning("Knots were supplied but will be ignored!\n")
      knots<-list()
   }

   # set the penalty order
   object$p.order<-c(2,mds.dim/2-1)

   # make the duchon splines object as usual
   object<-smooth.construct.ds.smooth.spec(object,data,knots)

   if(!is.null(object$xt$extra.penalty)){
      object<-extra.penalty(object)
   }

   # recover the stuff we want in the object
   object$term<-save.term
   object$dim<-save.dim
   data<-save.data

   class(object)<-"msg.smooth"
   object
}

Predict.matrix.msg.smooth<-function(object,data){

   save.dim<-object$dim
   save.term<-object$term

   object$term<-object$msg$term
   object$dim<-object$msg$dim

   #### MAGIC HAPPENS HERE!!!!

   if(!is.null(object$xt$bnd)){
      bnd<-object$xt$bnd
   }else{
      bnd<-NULL
   }

   mds.obj<-object$msg$mds.obj
   my.grid<-object$msg$grid
   mds.data<-as.data.frame(insert.mds(data,my.grid,mds.obj,bnd))


   # make some variable names up
   mds.names<-paste("mds-",1:dim(mds.data)[2],sep="")
   # remove any already in the data
   names(mds.data)<-mds.names

   Predict.matrix.duchon.spline(object,mds.data)
}
