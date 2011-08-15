/* Code for soap film smoothing. Copyright Simon Wood 2006-2008.

   R CMD SHLIB soap.c creates appropriate soap.so from this, can
   then be loaded by dyn.load("soap.so") and called with .C() */

/* SNIPSNIPSNIP */
//void in_out(double *, double *, double *,double *,double *,int *, int *, int * );

#include <stdlib.h>
#include "utils.h"

void in_out(double *bx, double *by, double *break_code, double *x,double *y,int *in, int *nb, int *n)
{
   in_out1(bx,by,break_code,x,y,in,nb,n);
//   in_out2(bx,by,break_code,x,y,in,nb,n);


}

void in_out1(double *bx, double *by, double *break_code, double *x,double *y,int *in, int *nb, int *n)
/* finds out whether points in arrays x,y are inside boundary or outside, by counting boundary 
   crossings. The boundaries nodes are defined by bx, by.  bx[i] and by[i] less than or equal to 
   break_code signals a break in the boundary (e.g. between island and external boundary.) Each 
   section of boundary is assumed to be a closed loop. nb is dimenion of bx and by; n is dimension 
   of x and y. `in' will contain a 1 for an interior point and a 0 otherwise, on exit. 
   Both bx[i] and by[i] or neither must be less than the break_code.
*/ 
{ double xx,yy,dum,x0,x1,y0,y1;
  int i,j,count,start,swap;
  for (i=0;i<*n;i++) { /* loop through all test points */
      xx=x[i];yy=y[i]; /* the current test point */
      start=0; /* start of current boundary section */
      for (count=0,j=0;j<*nb-1;j++) { /* loop through entire boundary */
	x0 = bx[j]; /* start node */
        if (x0 <= *break_code) start=j+1; /* next segment start */
        else { /* not a new section start */
	  x1 = bx[j+1]; /* end node */
          if (x1 <= *break_code) x1=bx[start]; /* must join up segment end */
          if (x0!=x1) { /* x0==x1 => segment immaterial to decision */
	    if (x1<x0) { dum=x0;x0=x1;x1=dum;swap=1;} else swap=0; /* ordered */
            if (x0<xx&&x1>=xx) { /* might have a crossing */
	      y0 = by[j]; /* start node y co-ord */
              y1 = by[j+1]; /* end node y co-ord */
              if (y1 <= *break_code) y1=by[start]; /* must join up */
              if (y0<=yy&&y1<=yy) count++; /* definite crossing */
              else { /* more detail needed to determine crossing */
		if (!(y0>yy&&y1>yy)) { /* could still be one */
                  if (swap) {dum=y0;y0=y1;y1=dum;}
                  dum = (xx-x0)*(y1-y0)/(x1-x0)+y0; /* at what y does vertical cross segment */
                  if (yy>=dum) count++; /* it's a crossing */
		} /* end - could still be one */
              } /* end - more detail */
            } /* end - might be a crossing */
          } /* end - does seg matter */
        } /* end - not skipped because break */
      } /* end boundary loop */
     if (count%2) in[i]=1;else in[i]=0; /* record result */
  } /* end x,y test loop */
}

void in_out2(double *bx, double *by, double *break_code, double *x,double *y,int *in, int *nb, int *n)
{ double xx,yy,dum,x0,x1,y0,y1;
  int i,j,count,start,swap;
  for (i=0;i<*n;i++) { /* loop through all test points */
      xx=-x[i];yy=-y[i]; /* the current test point */
      start=0; /* start of current boundary section */
      for (count=0,j=0;j<*nb-1;j++) { /* loop through entire boundary */
	x0 = -bx[j]; /* start node */
        if (x0 <= *break_code) start=j+1; /* next segment start */
        else { /* not a new section start */
	  x1 = -bx[j+1]; /* end node */
          if (x1 <= *break_code) x1=bx[start]; /* must join up segment end */
          if (x0!=x1) { /* x0==x1 => segment immaterial to decision */
	    if (x1<x0) { dum=x0;x0=x1;x1=dum;swap=1;} else swap=0; /* ordered */
            if (x0<xx&&x1>=xx) { /* might have a crossing */
	      y0 = -by[j]; /* start node y co-ord */
              y1 = -by[j+1]; /* end node y co-ord */
              if (y1 <= *break_code) y1=by[start]; /* must join up */
              if (y0<=yy&&y1<=yy) count++; /* definite crossing */
              else { /* more detail needed to determine crossing */
		if (!(y0>yy&&y1>yy)) { /* could still be one */
                  if (swap) {dum=y0;y0=y1;y1=dum;}
                  dum = (xx-x0)*(y1-y0)/(x1-x0)+y0; /* at what y does vertical cross segment */
                  if (yy>=dum) count++; /* it's a crossing */
		} /* end - could still be one */
              } /* end - more detail */
            } /* end - might be a crossing */
          } /* end - does seg matter */
        } /* end - not skipped because break */
      } /* end boundary loop */
     if (count%2) in[i]=1|in[i];else in[i]=0|in[i]; /* record result */
  } /* end x,y test loop */
}

/* SNIPSNIPSNIP */
