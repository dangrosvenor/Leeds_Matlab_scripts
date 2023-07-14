/* F = lininterp1f_multidim(X,Y,XI,Ydefault)
   Performs 1D interpolation on an arrays of vectors.
   X or Y are arrays of vectors, e.g. of size [N M]
   At the moment it operates assuming that either X is an array (size [N M]) and Y is a vector (size [N]; i.e. the same Y-values for all of the X-vectors), OR vice versa (size(X)=[1 N]; size(Y)=[N M]; the same X-values for a series of Y-vectors)
   
  IMPORTANT - the vectors for the arrays must be ordered so that the inerpolation data is stored in the first dimension. I.e. each interpolation vector is of length N. The other dimension (of length M above) can actually be further split up into more dimensions if desired. So, could have an array of size [N P Q R]
   The X vectors must be increasing monotonically, but do not have to be linearly spaced
   XI = the x values at which the interpolated values are required
   Ydefault = the value to return when the requested values are out of range
   -----  Daniel Grosvenor (1D routine written by Umberto Picchini, see below)
*/

/*  
F = lininterp1f(X,Y,XI,Ydefault) returns the value of the 1-D function Y at the
    points XI using linear interpolation. Length(F)=length(XI).
    The vector X specifies the coordinates of the underlying interval.
    Ydefault is returned for values of XI outside the coordinates in X.
    For lininterp1f to work properly:
    X        must be a monotonically increasing array;
    Y        must be an array with length(Y)=length(X);
	XI       must be an array.
    Ydefault must be a scalar value or an empty matrix [].

	Warning:  not much in the way of error checking, since this slows
              things down, so pay attention to the argument passed to the function!!!

    function lininterp1f(), V. 1.1 Umberto Picchini (umberto.picchini@uniroma1.it), November 2005


    Installation:
	--------------

    Simply copy the trapzf.dll file in a directory recognizable by MATLAB. You won't
	need lininterp1f.c to perform computations, but it is useful if you want to
	read the user instructions or to customize your code.
	If you modify lininterp1f.c then you must run the following command

    mex lininterp1f.c

    in order to obtain the new dll file.


Ex1: 
>> x = [1:1:1000];
>> y =log(sqrt(x+1.001)-1.001);
>> xv =[5:.001:100];
>> yinterp =lininterp1f(x,y,xv,[]);

Ex2: 
>> x=[1:1:1000];
>> y=log(sqrt(x+1.001)-1.001);
>> xv=[5:.001:100];
>> tic; y1=interp1(x,y,xv,'linear'); eval1=toc;  % the MATLAB standard linear interpolator
>> eval1

eval1 =

    0.14

>> tic; y2=lininterp1f(x,y,xv,[]); eval2=toc;
>> eval2

eval2 =

    0.05


Ex3: (run the following script)
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
tic;
rand('seed',0)
for(i=1:50)
    x = [1:1:1000] + rand(1,1000);
    y = x.^(1./x);
    xv = [3:.001:100];
    yv = interp1(x,y,xv,'linear');   % the MATLAB standard linear interpolator
end
eval1 = toc;

tic;
rand('seed',0)
for(i=1:50)
    x = [1:1:1000] + rand(1,1000);
    y = x.^(1./x);
    xv = [3:.001:100];
    yv = lininterp1f(x,y,xv,[]);
end
eval2 = toc;
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


then you obtain:

eval1 =

     9.1640


eval2 =

    3.4750
*/


#include "mex.h"

/* Input arguments */
#define X_data      prhs[0]
#define Y_data      prhs[1]
#define X_interp    prhs[2]
#define Y_default   prhs[3]

/* Output arguments */
#define Y_out       plhs[0]





/* 1-D interpolation routine */
static double lininterp1f(double *yinterp, double *xv, double *yv, double *x, double *ydefault, int m, int minterp)
{
  double temp;    
  int i, j; 
	int nrowsinterp, nrowsdata;
	nrowsinterp = minterp;
	nrowsdata = m;

 /* loop over all of the values we want to interpolate for */
	for (i=0; i<nrowsinterp; i++)
	{
	  if((x[i] < xv[0]) || (x[i] > xv[nrowsdata-1])) {
				yinterp[i] = *ydefault;
	  }
	  else {

/* loop through all of the x-values */
			  for(j=1; j<nrowsdata; j++) {
			    if(x[i]<=xv[j]) {
 temp = (x[i]-xv[j-1]) / (xv[j]-xv[j-1]) * (yv[j]-yv[j-1]) + yv[j-1];
 yinterp[i] = temp;
                       break;
			    }
			  }
	  }
	}

    return *yinterp;
} 


/* -------------------------------------------------------------------------------
   Wrapper routine to run the 1D interpolation for each part of the N-D array 
       by Daniel Grosvenor
----------------------------------------------------------------------------------*/

static double run_lininterp1f(double *yinterp2, double *xv, double *yv, double *x, double *ydefault, int *dims_xv, int ndims_xv, int *dims_yv, int ndims_yv, int m, int minterp)
{
  /* ydefault is the default value returned if x is outside the range of xv */
    int i, j, iloc; 
    int nloops;
    double *x1D,*y1D;
    double *y_out;

    y_out = (double *)malloc(m * sizeof(double) );

*/ 
mexPrintf("\nEntered wrapper routine, ndims_yv=%d, dims_yv[0]=%d, dims_yv[1]=%d\n",ndims_yv,dims_yv[0],dims_yv[1]); */


 if (ndims_xv>ndims_yv) {  /*case where xv is an array - assume yv is a vector */
   nloops=1;
     for (i=1;i<ndims_xv;i++) {
		       nloops = nloops * dims_xv[i];
     }

     iloc=0;
     for (i=1;i<=nloops;i++) {
     /* Do the actual computations in the 1D subroutine. m is the number of elements in the vector for the x-data and minterp is the number of elements for the x values that we want to inerpolate for */

       x1D = &xv[iloc];
       lininterp1f(&y_out,x1D,yv,x,ydefault,m,minterp);

          for (j=1;j<=m;j++) { 
             yinterp2[iloc+j-1]=y_out[j-1];
          }

       iloc=iloc+m; /* increment the pointer location */
     } /* for */

   } */ if */


 else {  /* yv is an array, not a vector. Assume xv is a vector  */
     nloops=1;
     for (i=1;i<ndims_yv;i++) {
		       nloops = nloops * dims_yv[i];
     }
     iloc=0;
     for (i=1;i<=nloops;i++) {
     /* Do the actual computations in the 1D subroutine. m is the number of elements in the vector for the x-data and minterp is the number of elements for the x values that we want to inerpolate for */

            y1D = &yv[iloc];   	 
            lininterp1f(&y_out,xv,y1D,x,ydefault,m,minterp);

            /* mexPrintf("\nDone one interpolation, yinterp2[0]=%f, yinterp2[50]=%f\n",yinterp2[0],yinterp2[50]); */


            for (j=1;j<=m;j++) { 
               yinterp2[iloc+j-1]=y_out[j-1];
            }
            iloc=iloc+m;
     } /* end of for *
 }  /* end of if-else */




     /* mexPrintf("\nEntered wrapper routine, ndims=%d, dims[0]=%d, dims[1]=%d\n",ndims_xv,dims_xv[0],dims_xv[1]);

     mexErrMsgTxt("Halting here\n"); */

 return *yinterp2;

}  /* end of subroutine */



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int m, minterp,ndims_xv,ndims_yv,ndims_x,i,nloops,iloc,iloc2;
double *yinterp3;
double *xv, *yv, *x, *ydefault, *y_out;
int *dims_out, *ndims_out, *dims_xv, *dims_yv, *dims_x;



/* Find the dimension of the data */
/*   m = mxGetNumberOfElements(X_data);  */
/* Find the dimension of the x values for which interpolation is desired */
   minterp = mxGetNumberOfElements(X_interp);
/* Create an mxArray of real values for the output */


   dims_xv = mxGetDimensions(prhs[0]);
   ndims_xv = mxGetNumberOfDimensions(prhs[0]); 

   ndims_yv = mxGetNumberOfDimensions(prhs[1]); 
   dims_yv = mxGetDimensions(prhs[1]);  

   if (ndims_xv>ndims_yv) {
     dims_out = dims_xv;
     ndims_out = ndims_xv;
     m = dims_xv[0];
   }
   else {
     dims_out = dims_yv;
     ndims_out = ndims_yv;
     m = dims_yv[0];
   }

 dims_out[0] = minterp;

/* Get the data passed in */
   xv       = mxGetPr(X_data);
   yv       = mxGetPr(Y_data);
   x        = mxGetPr(X_interp);
   ydefault = mxGetPr(Y_default);

 




     Y_out = mxCreateNumericArray(ndims_out, dims_out,
        mxDOUBLE_CLASS, mxREAL); 

     yinterp3  = mxGetPr(Y_out);

     iloc=0;
     iloc2=0;

    if (ndims_xv>ndims_yv) {  /*case where xv is an array - assume yv is a vector */
       nloops=1;
       for (i=1;i<ndims_xv;i++) {
		       nloops = nloops * dims_xv[i];
       }


      for (i=1;i<=nloops;i++) {
       /* Do the actual computations in the 1D subroutine. m is the number of elements in the vector for the x-data and minterp is the number of elements for the x values that we want to inerpolate for */

   
       /*       y_out = &yinterp[iloc];   */
       lininterp1f(&yinterp3[iloc2],&xv[iloc],yv,x,ydefault,m,minterp);      



       iloc=iloc+m; /* increment the pointer location */
       iloc2=iloc2+minterp;
      }
    } /* end if */


   else {  /* yv is an array, not a vector. Assume xv is a vector  */
     nloops=1;
     for (i=1;i<ndims_yv;i++) {
		       nloops = nloops * dims_yv[i];
     }
     for (i=1;i<=nloops;i++) {
     /* Do the actual computations in the 1D subroutine. m is the number of elements in the vector for the x-data and minterp is the number of elements for the x values that we want to inerpolate for */

       lininterp1f(&yinterp3[iloc2],xv,&yv[iloc],x,ydefault,m,minterp);

       /*

       iloc=iloc+m;
       iloc2=iloc2+minterp;
     }
 }  /* end of if-else */


 /* Pass the full N-D arrays to the wrapper routine  
   run_lininterp1f(yinterp3,xv,yv,x,ydefault,dims_xv,ndims_xv,dims_yv,ndims_yv,m,minterp); */



   return;
}
