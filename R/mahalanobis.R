# replacement of the mahalanobis() function in R
# automatically does P-M pseudo-inverse when p>n

#mahalanobis<-function (x, center, cov=NULL, inverted = FALSE, xc.svd=NULL, ...){
mahalanobis<-function (x, center, cov=NULL, inverted = FALSE, cov.inv=NULL, ...){

   x<-if (is.vector(x)) 
         matrix(x, ncol = length(x))
      else as.matrix(x)

   center<-as.matrix(center)

   if(nrow(x)>=ncol(x) & !is.null(cov)){
      # use the function in the stats package
      res<-stats::mahalanobis(x,center,cov,inverted,...)
   }else{
      # number of data
      n<-nrow(x)

      # column centred t(x) 
      xc<-t(sweep(x,2,colMeans(x)))
      #xc<-sweep(x,2,colMeans(x))
      
      y <- sweep(x, 2, center)# = (x - center)
      y<-t(y)

      # find the singular value decomposition (if we need to)
      #if(is.null(xc.svd)){
      #   xc.svd<-svd(xc)
      #}
      # find cov.inv if we need to....
      if(is.null(cov.inv)){
         # find the singular value decomposition
         xc.svd<-svd(xc)
         # only use non-"zero" elements of d, as in MASS::ginv()
         nz<-xc.svd$d > (sqrt(.Machine$double.eps) * xc.svd$d[1])

         # calculate the inverse of the covariance matrix
         cov.inv<-xc.svd$u[,nz]%*%(((1/xc.svd$d[nz])^2)*t(xc.svd$u[,nz]))*(n-1)
         #cov.inv<-xc.svd$v[,nz]%*%(((1/xc.svd$d[nz])^2)*t(xc.svd$v[,nz]))*(n-1)
      }

      # calculate the Mahalanobis distance
      res<-rowSums((t(y)%*%cov.inv)*t(y))
      #res<-rowSums((y%*%cov.inv)*y)

      # quick way of doing the above
      #Z<-t(y)%*%xc.svd$u[,nz]%*%diag(sqrt(n-1)/xc.svd$d[nz])
      #Z<-t(y)%*%xc.svd$u[,nz]%*%diag(1/xc.svd$d[nz])
      #res<-rowMeans(Z%*%t(Z)*(n-1))

      # attach the svd to the result, to save time later
      #attr(res,"xc.svd")<-xc.svd
      attr(res,"cov.inv")<-cov.inv
   }
   return(res)
}
