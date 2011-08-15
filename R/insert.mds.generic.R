# generic MDS insertion routine
# will call the WAD version if bnd is supplied
insert.mds.generic<-function(mds.obj,new.points,old.points,dist.metric="euclidean",bnd=NULL){

   # if bnd supplied return what insert.mds() returns
   if(!is.null(bnd)){
     return(insert.mds(new.points,old.points,mds.obj,bnd))
   }

   # otherwise do the non-WAD thing...

   big.points<-rbind(old.points,new.points)
   ind<-1:nrow(old.points)
   ind2<-(nrow(old.points)+1):nrow(big.points)

   lambda.inverse<-diag(1/mds.obj$eig[1:dim(mds.obj$points)[2]])

   # mahalanobis distances
   if(dist.metric=="mahalanobis"){

      cov.mat<-solve(cov(big.points))
      new.dist<-apply(big.points,1,mahalanobis,x=big.points,
                      cov=cov.mat,inverted=TRUE)[ind,]
   }else{
   # Euclidean
      new.dist<-as.matrix(dist(big.points,method=dist.metric,diag=T,upper=T))[ind,]
   }

   new.dist<-new.dist[,ind2]
   S<- -1/2*mds.obj$x
   d<- -(new.dist^2-diag(S))
   mds.points<-t(1/2*(lambda.inverse %*% t(mds.obj$points) %*% d))

   return(mds.points)
}
