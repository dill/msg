// Simon's algorithm for finding the path
// Copyright 2009-2010 David Lawrence Miller

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "wood.h"

// this doesn't actually do anything, eps gets set below.
double eps=1e-6;

void wood_path(int *len, int *start, double *x, double *y, int *nbnd, double *xbnd, double *ybnd, double *xref, double *yref, int *ngrid, double *refdelx, double *refdely, int *refio, int *nref, double *xstart, double *ystart, double *pathlen, int *faster, int *debug){
   // args:
   //   len         the length of x and y *
   //   start       point to start at *
   //   x,y         lists of x and y points
   //   nbnd        length of bnd
   //   bnd         the boundary that got in the way
   //   faster      0 = no speedup, 1= use cached paths
   //  return:
   //   length of the path

   // * these can be manipulated to control which we evaluate
   //  eg. in the insertion case
   double **bnd, p1[2], p2[2];
   int i,j,k,err=0, ilim, jstart;
   int app[2];

   // storage
   node** savedpaths;
   node* thispath=NULL;

   bnd=(double**)malloc(sizeof(double*)*(*nbnd));
   bnd[0]=(double*)malloc(sizeof(double)*(*nbnd)*2);
   // put bnd in the right format and allocate memory
   for(i=0; i<*nbnd; i++){
      bnd[i]=bnd[0]+i*2;
      bnd[i][0]=xbnd[i];
      bnd[i][1]=ybnd[i];
   }

   // first of all, set the epsilon to use...
   set_epsilon(*nbnd,xbnd,ybnd);

   // insertion counter
   k=0;

   if(*faster){
      create_refpaths(xref,yref,*nref,*refdelx, *refdely, *xstart, 
                      *ystart,refio,*ngrid,*nbnd,&savedpaths,bnd);
   }

   // first calculate all of the Euclidean paths
   get_euc_path(x, y, *nbnd, bnd, *len, pathlen, *start);

   // switch between insertion and full MDS 
   // insertion format is c(old,new)
   //  * so *start gives the index of the limit of the old points for insertion
   //    if we're not doing insertion then this is just the length of the 
   //    point vector
   
   if(*start != 0){
      // insertion loop variables 
      ilim=*start;
      jstart=*start;
   }else{
      /// full MDS loop variables
      ilim=*len;
      jstart=0;
   }

   // indexing here is rather sticky...
   // i,j   index the points
   // k     indexes the path lengths, only find pathlen[k] when we don't
   //        have a Euclidean path
   // l     counts the size of the saved path array
   // m     indexes the saved paths

   // #### Main for loops 
   for(i=0; i<ilim; i++){
      if(*start==0){ jstart=i+1;} // make sure that j is set right for full mds
      for(j=jstart; j<(*len); j++){
         // if no euclidean path was found, calculate the path
         if(pathlen[k] == (-1)){
            p1[0]=x[i]; p1[1]=y[i];
            p2[0]=x[j]; p2[1]=y[j];

            //printf("cat(\"p1=list(x=%.16f,y=%.16f);p2=list(x=%.16f,y=%.16f)\\n\")\n",
            //       p1[0],p1[1],p2[0],p2[1]);

            if(*debug==1){
               printf("# DEBUG: p1=list(x=%f,y=%f); p2=list(x=%f,y=%f);\n",p1[0],p1[1],p2[0],p2[1]);
            }

//            printf("cat(\"i=%d,j=%d\\n\")\n",i+1,j+1);


            if(*faster){
               // can we do an append?
               // do the append check for p1   
               append_check(p1, p2, *xstart, *ystart, *refdelx, *refdely, *ngrid, refio, *nbnd, bnd, app);

            }else{
               app[0]=0;
            }

            err=0;

            // if an append will work...
            if(app[0]!=0){
               err=+append_path(&savedpaths[app[1]],&thispath,p2,p1,app[0],*nbnd,bnd);
            }else{
               err=+make_bnd_path(p1,p2,*nbnd,bnd,&thispath,0);
            }

// DEBUG
//printf("cat(\"### first ###\\n\")\n");
//PrintPath(&thispath);

            if(err==1){
               FreeList(&thispath);
               err=make_bnd_path(p1,p2,*nbnd,bnd,&thispath,0);
            }

            // take the start path and optimize it...
            err=+iter_path(&thispath,*nbnd,bnd);

            if(err==1){
               FreeList(&thispath);
               err=make_bnd_path(p1,p2,*nbnd,bnd,&thispath,0);
               err=iter_path(&thispath,*nbnd,bnd);
            }

            // make sure that we can measure the length of the path
            if(thispath ==NULL){
               // use hypot() and warn!
               pathlen[k]=hypot(p2[0]-p1[0],p2[1]-p1[1]);
               printf("# ERROR: make_bnd_path FAILED. Using hypot()!\n");
               printf("# DEBUG: p1=list(x=%f,y=%f); p2=list(x=%f,y=%f);\n",p1[0],p1[1],p2[0],p2[1]);
            }else{
               // find the length of the path
               pathlen[k]=hull_length(&thispath);
            }

// DEBUG
//printf("cat(\"### final ###\\n\")\n");
//PrintPath(&thispath);
            FreeList(&thispath);

         }
         // increment pathlen counter
         k++;
      }    
   }

   // if we saved paths, free them!   
   if(*faster){
      // free all the saved paths
      for(i=0;i<((*ngrid)*(*ngrid)*((*ngrid)*(*ngrid)/2-1));i++){
         if(savedpaths[i]!=NULL){
            FreeList(&savedpaths[i]);
         }
      }
      free(savedpaths);
   }

   free(bnd[0]);
   free(bnd);
}

void get_euc_path(double x[], double y[], int nbnd, double **bnd, int npathlen,
                  double *pathlen, int full){
   /*
      get the within-area Euclidean distance, if possible
      args:
         x,y      points to check
         nbnd     length of boundary
         bnd      boundary
         pathlen  variable to store the lengths of the path in
         full     do full or insertion MDS?
   */
   int i,j,k, *retint,ilim,jstart;
   double p1[2], p2[2];

   k=0;

   // allocate memory for retint
   retint=(int*)malloc(sizeof(int)*(nbnd-1));
   for(i=0; i<(nbnd-1); i++){
      retint[i]=retint[0]+i;
   }

   // get loop variables
   if(full != 0){
      // insertion loop variables 
      ilim=full;
      jstart=full;
   }else{
      /// full MDS loop variables
      ilim=npathlen;
      jstart=0;
   }

   for(i=0; i<ilim; i++){
      if(full==0){
         jstart=i+1;
      }
      for(j=jstart; j<npathlen; j++){
         p1[0]=x[i]; p1[1]=y[i]; // set p1
         p2[0]=x[j]; p2[1]=y[j]; // set p2
 
         // check to see if we have to do the path finding or
         // just the hypotenuse 
         sp_do_intersect(p1, p2, nbnd, bnd, retint);
   
         //                           vvv just hypot when we touch only 1 vertex 
         if(iarrsum((nbnd-1),retint) > 0){
            pathlen[k] = (-1);
         }else{
            pathlen[k]=hypot(p2[0]-p1[0],p2[1]-p1[1]);
         }
         k++; // increment pathlen counter
      }
   }
   free(retint);
}

int iter_path(node** mypath,int nbnd, double **bnd){
   /*
      args:
         path        the shortest within-domain path
         nbnd        length of the boundary
         bnd         the boundary
   */

   int conv, conv_stop;
   node* prevpath=NULL;

   // convergence stop
   conv=0;
   conv_stop=10;

   // keep going until we don't remove any more points.
   do{
      // free up prevpath before we copy onto it
      if(conv>0){
         FreeList(&prevpath);
      }

      // save previous path
      //prev.path<-my.path
      CopyList(*mypath,&prevpath);

      // add new vertices
      alter_step(mypath,nbnd,bnd);
      // DEBUG
      //printf("cat(\"### alter_step ###\\n\")\n");
      //PrintPath(mypath);

      // delete step, remove anything that doesn't need to be there
      delete_step(mypath,nbnd,bnd);
      // DEBUG
      //printf("cat(\"### delete_step ###\\n\")\n");
      //PrintPath(mypath);

      // increment convergence stopper 
      conv++;

   } while( !(has_converged(prevpath,*mypath)) & (conv<conv_stop) );

   if(!(has_converged(prevpath,*mypath))){
      printf("# WARNING: path find finished without convergence!\n");
      printf("# conv = %d\n",conv);
      printf("# convergence = %d\n",has_converged(prevpath,*mypath) );
      FreeList(&prevpath);
      return 1;
   }

   FreeList(&prevpath);
   return 0;
}


// create the initial path:
// p1, p1 1st intersection, some of bnd, p2 1st intersection, p2
int make_bnd_path(double p1[], double p2[], int nbnd, double **bnd, node** path, int delfirst)
{
   /* Args:
   * p1, p2 points
   * nbnd length of boundary
   * bnd boundary
   * path head node of the linked list
   * delfirst delete before finding the lengths? 0=yes,1=no
   */
 
   double ip1[2],ip2[2], curr_insert[2];
   int intind[2],i,start;
   int err=0, sortind=0,indi;
   node* bnd1 = NULL;
   node* bnd2 = NULL;
 
   // find the first intersection between p1, p2 and the boundary side that
   // each point intersects
   err=first_ips(p1, p2, nbnd, bnd, ip1, ip2, intind);
 
   // if there are no errors
   if(err==0){
 
      // first sort the intersection indices
      indi=intind[0];
      itwosort(intind);
      // set a var here to know if we switched, needed later
      if(indi!=intind[0]){
         sortind=1;
      }

      // since we ordered intind first, we don't need to worry too much
      // we put the bit that is contiguous in first eg:
      // 0 1 2 3} [4 5 6] {7 8 9 <<<- do the [ ] bit first!
 
      // push everything in
      // vvvvvvvvvv <- since we want it to be inclusive
      for(i=(intind[0]+1);i<=(intind[1]);i++){
         curr_insert[0]=bnd[i][0];
         curr_insert[1]=bnd[i][1];
         Push(&bnd1,curr_insert);
      }
 
      // create the second boundary segment
      // want intind[1]+1:intind[0]
      if(intind[1]!=nbnd){
         start=intind[1];
      }else{
         start=1; // handle the case where start is actually the end
      }
 
      // insert until we hit the end
      for(i=(start+1);i<nbnd;i++){
         curr_insert[0]=bnd[i][0];
         curr_insert[1]=bnd[i][1];
         Push(&bnd2,curr_insert);
      }
 
      // insert from the start back to intend[0]
      if(intind[0]!=0){
         for(i=0;i<(intind[0]+1);i++){
            curr_insert[0]=bnd[i][0];
            curr_insert[1]=bnd[i][1];
            Push(&bnd2,curr_insert);
         }
      }


      if(!sortind){
         curr_insert[0]=ip1[0]; curr_insert[1]=ip1[1];
         AppendNode(&bnd1,curr_insert);
         Push(&bnd2,curr_insert);
   

         if(!((fabs(ip1[0]-p1[0]) <eps) & (fabs(ip1[1]-p1[1]) <eps))){
            curr_insert[0]=p1[0]; curr_insert[1]=p1[1];
            AppendNode(&bnd1,curr_insert);
            Push(&bnd2,curr_insert);
         }
 
         curr_insert[0]=ip2[0]; curr_insert[1]=ip2[1];
         Push(&bnd1,curr_insert);
         AppendNode(&bnd2,curr_insert);
   
         if(!((fabs(ip2[0]-p2[0]) <eps) & (fabs(ip2[1]-p2[1]) <eps)) ){
            curr_insert[0]=p2[0]; curr_insert[1]=p2[1];
            Push(&bnd1,curr_insert);
            AppendNode(&bnd2,curr_insert);
         }
      }else{
         curr_insert[0]=ip1[0]; curr_insert[1]=ip1[1];
         Push(&bnd1,curr_insert);
         AppendNode(&bnd2,curr_insert);
   
         if(!((fabs(ip1[0]-p1[0]) <eps) & (fabs(ip1[1]-p1[1]) <eps) )){
            curr_insert[0]=p1[0]; curr_insert[1]=p1[1];
            Push(&bnd1,curr_insert);
            AppendNode(&bnd2,curr_insert);
         }
  
         curr_insert[0]=ip2[0]; curr_insert[1]=ip2[1];
         AppendNode(&bnd1,curr_insert);
         Push(&bnd2,curr_insert);
   
         if(!((fabs(ip2[0]-p2[0]) <eps) & (fabs(ip2[1]-p2[1]) <eps) )){
            curr_insert[0]=p2[0]; curr_insert[1]=p2[1];
            AppendNode(&bnd1,curr_insert);
            Push(&bnd2,curr_insert);
         }
      }
//printf("cat(\"# pre-final makes ###\\n\")\n");
//PrintPath(&bnd1);
//PrintPath(&bnd2);
//printf("cat(\"# pre-final makes ###\\n\")\n");
 
// DEBUG
//printf("cat(\" path1 ###\\n\")\n");
//PrintPath(&bnd1);
//printf("# path1 ###\n");
//printf("ip1<-list(x=%f,y=%f)\n",ip1[0],ip1[1]);
//printf("ip2<-list(x=%f,y=%f)\n",ip2[0],ip2[1]);
//printf("p1<- list(x=%f,y=%f)\n",p1[0],p1[1]);
//printf("p2<- list(x=%f,y=%f)\n",p2[0],p2[1]);
//printf("points(ip1)\n");
//printf("points(ip2)\n");
//printf("points(p1,pch=19)\n");
//printf("points(p2,pch=19)\n");
//printf("scan()\n");
//printf("cat(\" path2 ###\\n\")\n");
//PrintPath(&bnd2);
//printf("ip1<-list(x=%f,y=%f)\n",ip1[0],ip1[1]);
//printf("ip2<-list(x=%f,y=%f)\n",ip2[0],ip2[1]);
//printf("p1<- list(x=%f,y=%f)\n",p1[0],p1[1]);
//printf("p2<- list(x=%f,y=%f)\n",p2[0],p2[1]);
//printf("points(ip1)\n");
//printf("points(ip2)\n");
//printf("points(p1,pch=19)\n");
//printf("points(p2,pch=19)\n");
//printf("scan()\n");
//printf("# path2 ###\n");
      // delete before testing length?
      if(delfirst==0){
            if(Length(&bnd1)>2){
               delete_step(&bnd1, nbnd, bnd);
            }
            if(Length(&bnd2)>2){
               delete_step(&bnd2, nbnd, bnd);
            }
      }

//printf("cat(\"# final makes ###\\n\")\n");
//PrintPath(&bnd1);
//PrintPath(&bnd2);
//printf("cat(\"# final makes ###\\n\")\n");
 
      // pick the shorter path to return
      if(hull_length(&bnd1)<hull_length(&bnd2)){
         CopyList(bnd1,path);
         FreeList(&bnd1);
         FreeList(&bnd2);
      }else{
         CopyList(bnd2,path);
         FreeList(&bnd1);
         FreeList(&bnd2);
      }
 
      return 0;
   }else{
      *path=NULL;
//      printf("# ERROR: make_bnd_path FAILED. Error returned from first_ips\n");
// printf("# DEBUG: p1=list(x=%f,y=%f); p2=list(x=%f,y=%f);\n",p1[0],p1[1],p2[0],p2[1]);
      return 1;
   }// end of error if()
 
}

// append one path to the end of another
int append_path(node** oldpath, node** newpath, double p1[2], double p2[2], int end,
                int nbnd, double **bnd){
   /*
    *
    *
   */

   if(Length(oldpath)<5){
      return 1;
   }

   // blank what is currently in newpath
   FreeList(newpath);

   // replace with what's in oldpath
   CopyList(*oldpath,newpath);

   // if the point->endpoint path is Euclidean in the domain then
   // just add that point
   if(end==1){
      AppendNode(newpath,p1);
      Push(newpath,p2);
   }else{
      Push(newpath,p1);
      AppendNode(newpath,p2);
   }

   delete_step(newpath,nbnd,bnd);

   // done
   return 0;
}


// iterate over the points in the path:
// delete as many points as possible, making sure that the
// path is still outside
void delete_step(node** path, int nbnd, double **bnd){
   // Args:
   //  path     the current path
   //  nbnd     # elements of boundary
   //  bnd      the boundary
   // Return:
   //           revised path with dropped vertices
   
   int i, *intbnd, in[1];
   double mytrip[3][2], p1[2], p2[2], xmp[1], ymp[1], *bx, *by, xmin, ymin, mina[2], break_code;
	int inout_n=1;
   // convergence stop
   int conv=0;
   int conv_stop=50;
   
   // some linked lists
   node* current=NULL;   // iterator
   node* prevpath=NULL;
   // needed for deletion
   node* start_ptr=NULL;
   node* end_ptr=NULL;

   // setup intbnd
   intbnd=(int*)malloc(sizeof(int)*(nbnd-1));
   for(i=0; i<(nbnd-1); i++){
      intbnd[i]=intbnd[0]+i;
   }

   ////// needed later on for the in_out() call
   bx=(double*)malloc(sizeof(double)*nbnd);
   by=(double*)malloc(sizeof(double)*nbnd);
   for(i=0; i<nbnd; i++){
      bx[i]=bx[0]+i;
      by[i]=by[0]+i;
      bx[i]=bnd[i][0];
      by[i]=bnd[i][1];
   }
   
   // to handle holes, multiple boundaries
   // ignore this at the moment
   xmin=minarr(nbnd,bx);
   ymin=minarr(nbnd,by);
   mina[0] = xmin; mina[1] = ymin;
   break_code=minarr(2,mina)-1;

   // keep going until we don't remove any more points.
   do{
      if(conv>0){
         FreeList(&prevpath);
      }

      // use prevpath to keep a copy of the previous path for comparison
      CopyList(*path,&prevpath);

      // start point for triplet selection
      current = *path;   // iterator

      // loop over the full path
      while(current!=NULL){

         // equivalent of some ANDs in the above, but doesn't cause a memory
         // problem that occurs when  the while() evaluates all of the 
         // conditions.
         if(current->next==NULL){
            break;
         }else if(current->next->next==NULL){
            break;
         }

         // given current is at point i here...

         // create the current triplet to inspect
         mytrip[0][0]=current->data[0];
         mytrip[0][1]=current->data[1];
        	current=current->next;
         mytrip[1][0]=current->data[0];
         mytrip[1][1]=current->data[1];
        	current=current->next;
         mytrip[2][0]=current->data[0];
         mytrip[2][1]=current->data[1];

         // pointer is now at i+2

         // if we are going forward and back again, just remove the points
         if( (fabs(mytrip[0][0]-mytrip[2][0])<eps) & 
            (fabs(mytrip[0][1]-mytrip[2][1])<eps)){
            // current is sitting at the 3rd entry
            // create a pointer to that
            end_ptr=current->next; // pointer to i+3

            // break out of the loop if that would crash
            if(current->next==NULL){
               break;
            }

            // go back twice
            current=current->prev; // pointer now at i+1
            current=current->prev; // pointer at i

            // change where next points to
            current->next=end_ptr; // point i next to i+3

            start_ptr=current; // pointer to i

            // go forward again (remember next has changed)
            current=current->next; // back to i+3

            // free memory of i+2 and i+1
            free(current->prev->prev);
            free(current->prev);

            // change previous
            current->prev=start_ptr; //set i+3 prev to i
         }else{
            //// This is all setup for the next if() vvvvv

            // start and end points of the trip in p1 and p2
            p1[0]=mytrip[0][0];  p1[1]=mytrip[0][1];
            p2[0]=mytrip[2][0];  p2[1]=mytrip[2][1];

            // see if the line between p1 and p2 intersects the boundary
            sp_do_intersect(p1,p2, nbnd,bnd,intbnd);

            // midpoints
            xmp[0]=(mytrip[2][0]+mytrip[0][0])/2;
            ymp[0]=(mytrip[2][1]+mytrip[0][1])/2;

            in[0]=0;
            in_out(bx,by,&break_code,xmp,ymp,in, &nbnd,&inout_n);
            // if deleting point i makes the resulting line cross the
            // the boundary then keep it in the path, else get rid of it 

            // first part asks if there are any intersections, if there are
            // none then that's okay. Second part asks if midpoints are inside.
            if((iarrsum((nbnd-1),intbnd)==0) & in[0]){
               // remove point i+1 by setting next pointer from i to 
               // i+2, and prev from i+2 to i

               // current is sitting at the 3rd entry
               // create a pointer to that
               end_ptr=current; // pointer to i+2

               // go back twice
               current=current->prev; // pointer now at i+1
               current=current->prev; // pointer at i

               // change where next points to
               current->next=end_ptr; // point i next to i+2
               start_ptr=current; // pointer to i

               // go forward again (remember next has changed)
               current=current->next; // current at i+2
               
               // free (i+1)'s memory
               free(current->prev);

               // change previous
               current->prev=start_ptr; //set i+2 prev to i

            } // end if on del middle 

         }// end of main if 
        	   current=current->prev; // go back to i, need this to catch all triplets

      } // end iteration over path
      conv++; // increment run counter

   } while( !(has_converged(prevpath,*path)) & (conv<conv_stop)); // end of do loop

   // free some memory
   free(bx); free(by);
   free(intbnd);
   FreeList(&prevpath);
}           

// alter the path
void alter_step(node** path, int nbnd, double **bnd)
{
   // Args:
   //  path     the current path
   //  nbnd     # elements of boundary
   //  bnd      the boundary
   // Return:
   //           revised path with added/ammended vertices

   double ep1[2], ep2[2], mid[2], triplen, tp1[2], tp2[2];
   int err;
   node* prevpath=NULL;
   node* newpath=NULL;
   node* end_ptr=NULL;
   node* end1=NULL;
   node* end2=NULL;
   node* current=NULL;

   // convergence stop
   int conv=0;
   int conv_stop=50;

   int oneshot;
   oneshot=Length(path);

   // iterate over the points in the path:
   // alter the path, until on two(?) consecutive runs there are
   // no changes to the path
   do{
      // iterator
      current = *path;

      if(conv>0){
         FreeList(&prevpath);
      }
      // use prevpath to keep a copy of the previous path for comparison            
      CopyList(*path,&prevpath);

      // start point for triplet selection
      while(current!=NULL){
         // equivalent of some ANDs in the above, but doesn't cause a memory
         // problem when the while() evaluates all of the conditions.
         if(current->next==NULL){
            break;
         }else if(current->next->next==NULL){
            break;
         }
         // for each point i, look at the line i to i+2

         // create the current triplet to inspect
         // make a copy of it and work out its length
         ep1[0]=current->data[0];
         ep1[1]=current->data[1];
        	current=current->next;
         mid[0]=current->data[0];
         mid[1]=current->data[1];
        	current=current->next;
         ep2[0]=current->data[0];
         ep2[1]=current->data[1];

         // current at i+2
         
         // calculate old trip length...
         triplen=hypot(mid[0]-ep1[0],mid[1]-ep1[1]);
         triplen=triplen+hypot(ep2[0]-mid[0],ep2[1]-mid[1]);

         // does it go inside-outside-inside?
         if(facing(ep1, ep2, nbnd, bnd)){
            // create a new path
            err=make_bnd_path(ep1,ep2,nbnd,bnd,&newpath,0);

            // provided there were no errors in the path making...
            if(err==0){
               // only insert the path if it's better!
               if((hull_length(&newpath)<=triplen)){

                  ////////// getting the path into the right format...
                  // only reverse the order if we need to...
                  // that is when the line joining current element of the full path
                  // and the first element of the new path is not inside the domain. 
                  tp1[0]=current->data[0];
                  tp1[1]=current->data[1];
                  tp2[0]=newpath->data[0];
                  tp2[1]=newpath->data[1];

                  if((fabs(tp2[0]-tp1[0]) <eps) & (fabs(tp2[1]-tp1[1]) <eps)){
                     ReverseList(&newpath);
                  }else{
                     end1=current;
                     end2=newpath;

                     while(end1->next !=NULL){
                        end1=end1->next;
                     }
                     while(end2->next !=NULL){
                        end2=end2->next;
                     }

                     tp1[0]=end1->data[0];
                     tp1[1]=end1->data[1];
                     tp2[0]=end2->data[0];
                     tp2[1]=end2->data[1];
                     
                     if(!((fabs(tp2[0]-tp1[0]) <eps) & (fabs(tp2[1]-tp1[1]) <eps))){
                        ReverseList(&newpath);
                     }
                  }

                  // remove the first and last entries in newpath, since otherwise
                  // we duplicated ep1 and ep2
                  if(Length(&newpath)>=3){
                     DelTopBot(&newpath);
                  }
                  /////////// done getting the path in the right format
//printf("cat(\"### path to insert\\n\")\n");
//PrintPath(&newpath);

                  // modify path
                  // from i, stitch in newpath up to i+2
                  // current is at i+2 now
                  
                  // create a pointer to that
                  end_ptr=current; // pointer to i+2
                  
                  // go back twice
                  current=current->prev; // pointer now at i+1
                  current=current->prev; // pointer at i
                  
                  // change where next points to
                  newpath->prev=current; // point head of newpath back 

                  free(current->next); // free the memory at i+1

                  current->next=newpath; // point i next to the head of newpath
      
                  // fast forward to the end of newpath
                  while(current->next!=NULL){
                     current=current->next;
                  }
                  // current now at end of newpath
   
                  // point the end of newpath to i+2
                  end_ptr->prev=current; // set previous of i+2 to be the end of newpath
                  current->next=end_ptr;
   
               }else{
                  // free the path we made
                  FreeList(&newpath);
                  current=current->prev;
               }// end insert if 
            }
            newpath=NULL; // make sure that newpath doesn't point anywhere
         }else{
            current=current->prev;
         } // end facing
      } // end of iteration over the path
      conv++;

      // if the initial path was just of length three, don't iterate
      // since we can't do any better...
      if(oneshot==3){
         break;
      }
   } while( !(has_converged(prevpath,*path)) & (conv<conv_stop) );
   //end of main do
   FreeList(&prevpath);
}

// check convergence
int has_converged(node* path1, node* path2)
{
   // args
   //    path1, path2      the two paths to compare
   // return
   //    1 if the paths are the same, 0 otherwise
   node* current1 = path1;
   node* current2 = path2;

   // first check length then their contents
   if(Length(&path1)==Length(&path2)){
      // check if the new and old paths are the same
      while(current1!=NULL){
         if(( (current1->data[0]) != (current2->data[0]) ) & 
            ( (current1->data[1]) != (current2->data[1]) )) return 0;
         current1=current1->next;
         current2=current2->next;
      }
   }else return 0;

   return 1;
}
