#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>

#define  Min(a,b)   ( (a) < (b) ? (a) : (b) )
#define  Max(a,b)   ( (a) > (b) ? (a) : (b) )
#define  _(i,j,M,N)    ( ((i)<0 || (i)>=(M) || (j)<0 || (j)>=(N)) ? _rngerr(i,j,M,N,__LINE__) : (i)+(j)*(M) )
int _rngerr( int i, int j,
          int M, int N,
          int line
        )
{
  char s[120];
  sprintf(s,"Index out of range (i=%i, j=%i, M=%i, N=%i) at line %i.",i,j,M,N,line);
  mexErrMsgTxt(s);
  return 0;
}

/**
  * Bilinear interpolation in 2D matrix
  **/
double               // returns X(i,j)
BilinearInterpolation(
  double  *X,        // matrix
  int      M,
  int      N,        // size if matrix
  double   ir,
  double   jr,        // indices to matrix, ranging from 0 to M-1, resp. N-1
  double   xInvalid  // value to return if i, j are too near rims of X that interpolation cannot be done
) {
  double x;
  int ii, ji;

  /* ir, jr set to float part of i(n), j(n),
     ii, ji set to floor of i(n), j(n) */
  ir -= ii = int(floor(ir));
  jr -= ji = int(floor(jr));

  /* perform bilinear interpolation */
  if ( ir==0 )
    if ( jr==0 )
      x = (ii>=0 && ji>=0 && ii<M   && ji<N  ) ?
           X[_(ii,ji,M,N)] :
           xInvalid;
    else
      x = (ii>=0 && ji>=0 && ii<M   && ji<N-1) ?
          (1-jr)*X[_(ii,ji  ,M,N)] +
             jr *X[_(ii,ji+1,M,N)] :
          xInvalid;
  else
    if ( jr==0 )
      x = (ii>=0 && ji>=0 && ii<M-1 && ji<N  ) ?
          (1-ir)*X[_(ii  ,ji,M,N)] +
             ir *X[_(ii+1,ji,M,N)] :
           xInvalid;
    else
      x = (ii>=0 && ji>=0 && ii<M-1 && ji<N-1) ?
          (1-jr)*((1-ir)*X[_(ii  ,ji  ,M,N)] +
                     ir *X[_(ii+1,ji  ,M,N)]) +
             jr *((1-ir)*X[_(ii  ,ji+1,M,N)] +
                     ir *X[_(ii+1,ji+1,M,N)]) :
          xInvalid;
  return x;
}


/**
  * One-dimensional subpixel detection from 3 points with coordinates
  * (-1,f1), (0,f2), (1,f3). It is assumed that f0 is local extrem among fn, fp.
  * Extrem is searched with subinteger precision by fitting parabola.
  **/
double   // returns function value in the extrem
SubpixDetect1D(
  double  f1,
  double  f2,
  double  f3,   // 3 successive points
  double *x,    // abscissa of local extrem
  double *dfxx  // 2nd derivative in the extrem
) {
  double q = f1-f3;
  *dfxx = f1-2*f2+f3;
  *x = 0.5*q/(*dfxx);
//  if ( fabs(*x)>0.5 ) printf("%f : %f %f %f\n",*x,f1,f2,f3);
  return q*q/8/(*dfxx)+f2;
}


void mexFunction( int dstn, mxArray **dstp, int srcn, const mxArray **srcp )
{
  mxArray *dst_L;
  #define src_I  srcp[0]
  #define src_S  srcp[1]

  double *I, *L, *tmpL, ci, cj, dfii, dfij, dfjj, Iij, Ikl, f;
  int M, N, SM, SN, offsL, i, j, k, l, notmax, nrowsL, maxnext;

  if ( srcn != 2 ) mexErrMsgTxt("Bad number of input arguments.");
  if ( mxGetM(src_S)*mxGetN(src_S) != 2 ) mexErrMsgTxt("Bad size of 2nd argument.");
  M = mxGetM(src_I);
  N = mxGetN(src_I);
  I = mxGetPr(src_I);
  SM = int(mxGetPr(src_S)[0]);
  SN = int(mxGetPr(src_S)[1]);

  nrowsL = 5; // number of rows of L
  maxnext = M*N/ ( ( SM &&  SN) ? 4 :
                   ( SM && !SN) ? 2 :
                   (!SM &&  SN) ? 2 :
                                  1 );  // maximum number of extrems
  tmpL = (double*)mxCalloc(maxnext*nrowsL,sizeof(double));
  if ( !tmpL ) mexErrMsgTxt("Out of memory.");

  // loop for all pixels of the input image except  rims
  offsL = 0;
  for ( i = SM; i < M-SM; i++ )
  for ( j = SN; j < N-SN; j++ ) {
    Iij = I[_(i,j,M,N)];

    // find if I(i,j) is a local maximum
//    notmax = notmin = 0;
    for ( k = -SM; k <= SM; k++ )
    for ( l = -SN; l <= SN; l++ )
      if ( k || l )
        if ( Iij <= I[_(i+k,j+l,M,N)] ) goto NO_EXTREM;
//        Ikl = I[_(i+k,j+l,M,N)];
//        notmax += (Iij <= Ikl);
//        notmin += (Iij >= Ikl);
//        if ( notmin && notmax ) goto NO_EXTREM;

    if ( SM )
      f = SubpixDetect1D(
        I[_(i-1,j  ,M,N)],
        Iij,
        I[_(i+1,j  ,M,N)],
        &ci,
        &dfii
      );
    else {
      ci = 0;
      dfii = mxGetNaN();
    }
    if ( SN )
      f = SubpixDetect1D(
        I[_(i  ,j-1,M,N)],
        Iij,
        I[_(i  ,j+1,M,N)],
        &cj,
        &dfjj
      );
    else {
      cj = 0;
      dfjj = mxGetNaN();
    }

    tmpL[offsL  ] = i + ci + 1;
    tmpL[offsL+1] = j + cj + 1; // +1 are present because MATLAB starts indexing arrays from 1
    tmpL[offsL+2] = (SM && SN) ? BilinearInterpolation(I,M,N,i+ci,j+cj,mxGetNaN()) : f;
    tmpL[offsL+3] = dfii;
    tmpL[offsL+4] = dfjj;
    offsL += 5;
    if ( offsL >= maxnext*nrowsL ) {
      printf("%i %i\n",offsL,maxnext*nrowsL);
      mexErrMsgTxt("Internal error, report to pajdla@cmp.felk.cvut.cz.");
    }

    NO_EXTREM: ;
  }

  // allocating the output list L and removing there the maxima from tmpL
  dst_L = mxCreateDoubleMatrix(nrowsL,offsL/nrowsL,mxREAL);
  if ( dst_L == 0 ) mexErrMsgTxt("Out of memory.");
  L = mxGetPr(dst_L);
  for ( i = 0; i < offsL; i++ ) L[i] = tmpL[i];
  mxFree(tmpL);

  dstp[0] = dst_L;
}
