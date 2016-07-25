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

void mexFunction( int dstn, mxArray **dstp, int srcn, const mxArray **srcp )
{
  mxArray *dst_L;
  #define src_I  srcp[0]
  #define src_S  srcp[1]

  double *I, *L, *tmpL, ci, cj, dfii, dfij, dfjj, Iij, Ikl, f, J;
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
                                  1 );  // maximum number of extremas
  tmpL = (double*)mxCalloc(maxnext*nrowsL,sizeof(double));
  if ( !tmpL ) mexErrMsgTxt("Out of memory.");

  // loop for all pixels of the input image except  rims
  offsL = 0;
  for ( i = SM; i < M-SM; i++ )
  for ( j = SN; j < N-SN; j++ ) 
    {
     Iij = I[_(i,j,M,N)];
     if (Iij==0) goto NO_EXTREM;
     for ( k = -SM; k <= SM; k++ )
      for ( l = -SN; l <= SN; l++ )
       if ( k || l )
         { 
          J = I[_(i+k,j+l,M,N)];
          if (Iij < J) goto NO_EXTREM;
          // if (J && (Iij > J)) I[_(i+k,j+l,M,N)] = 0; 
         }
     tmpL[offsL  ] = i + 1;
     tmpL[offsL+1] = j + 1; // +1 are present because MATLAB starts indexing arrays from 1
     tmpL[offsL+2] = Iij;
     tmpL[offsL+3] = 0;
     tmpL[offsL+4] = 0;
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
