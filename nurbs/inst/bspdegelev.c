/*
   Copyright (C) 2000 Mark Spink

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include "mexmat.h"

#define min(x,y) (x < y ? x : y) 
#define max(x,y) (x > y ? x : y) 

/* Compute logarithm of the gamma function */
/* Algorithm from 'Numerical Recipes in C, 2nd Edition' pg214. */
double gammaln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6] = {76.18009172947146,-86.50532032291677,
                          24.01409824083091,-1.231739572450155,
                          0.12086650973866179e-2, -0.5395239384953e-5};
  int j;
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j=0; j<=5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

/* computes ln(n!) */
/* Numerical Recipes in C */
/* Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215. */
double factln(int n)
{
  static int ntop = 0;
  static double a[101];
  
  if (n <= 1) return 0.0;
  while (n > ntop)
  {
    ++ntop;
    a[ntop] = gammaln(ntop+1.0);
  }
  return a[n];
}

/* Computes the binomial coefficient. */

/*     ( n )      n!                  */
/*     (   ) = --------               */
/*     ( k )   k!(n-k)!               */

/* Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215. */
double bincoeff(int n, int k)
{
  return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

/* Degree elevate a B-Spline t times. */
/* i.e. from d to d+t */

/* INPUT: */

/*   n,p,U,Pw,t  */

/* OUTPUT:  */

/*   nh,Uh,Qw  */

/* Modified version of Algorithm A5.9 from 'The NURBS BOOK' pg206. */
int bspdegelev(int d, double *c, int mc, int nc, double *k, int nk, 
               int t, int *nh, double *ic, double *ik)
{
  int row,col;

  int ierr = 0;
  int i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, cind, oldr, mul;
  int n, lbz, rbz, save, tr, kj, first, kind, last, bet, ii;
  double inv, ua, ub, numer, den, alf, gam;
  double **bezalfs, **bpts, **ebpts, **Nextbpts, *alfs; 

  double **ctrl  = vec2mat(c, mc, nc);
  double **ictrl = vec2mat(ic, mc, nc*(t+1));

  n = nc - 1;
  
  bezalfs = matrix(d+1,d+t+1);
  bpts = matrix(mc,d+1);
  ebpts = matrix(mc,d+t+1);
  Nextbpts = matrix(mc,d+1);
  alfs = (double *) mxMalloc(d*sizeof(double));

  m = n + d + 1;
  ph = d + t;
  ph2 = ph / 2;

  /* compute bezier degree elevation coefficients   */
  bezalfs[0][0] = bezalfs[ph][d] = 1.0;

  for (i = 1; i <= ph2; i++)
  {
    inv = 1.0 / bincoeff(ph,i);
    mpi = min(d,i);
    
    for (j = max(0,i-t); j <= mpi; j++)
      bezalfs[i][j] = inv * bincoeff(d,j) * bincoeff(t,i-j);
  }    
  
  for (i = ph2+1; i <= ph-1; i++)
  {
    mpi = min(d, i);
    for (j = max(0,i-t); j <= mpi; j++)
      bezalfs[i][j] = bezalfs[ph-i][d-j];
  } 
      
  mh = ph;
  kind = ph+1;
  r = -1;
  a = d;
  b = d+1;
  cind = 1;
  ua = k[0];

  for (ii = 0; ii < mc; ii++)
    ictrl[0][ii] = ctrl[0][ii];
  
  for (i = 0; i <= ph; i++)
    ik[i] = ua;
    
  /* initialise first bezier seg */
  for (i = 0; i <= d; i++)
    for (ii = 0; ii < mc; ii++)
      bpts[i][ii] = ctrl[i][ii];  

  /* big loop thru knot vector */
  while (b < m)
  {
    i = b;
    while (b < m && k[b] == k[b+1])
      b++;

    mul = b - i + 1;
    mh += mul + t;
    ub = k[b];
    oldr = r;
    r = d - mul;
    
    /* insert knot u(b) r times */
    if (oldr > 0)
      lbz = (oldr+2) / 2;
    else
      lbz = 1;

    if (r > 0)
      rbz = ph - (r+1)/2;
    else
      rbz = ph;  

    if (r > 0)
    {
      /* insert knot to get bezier segment */
      numer = ub - ua;
      for (q = d; q > mul; q--)
        alfs[q-mul-1] = numer / (k[a+q]-ua);
      for (j = 1; j <= r; j++)  
      {
        save = r - j;
        s = mul + j;            

        for (q = d; q >= s; q--)
          for (ii = 0; ii < mc; ii++)
            bpts[q][ii] = alfs[q-s]*bpts[q][ii]+(1.0-alfs[q-s])*bpts[q-1][ii];

        for (ii = 0; ii < mc; ii++)
          Nextbpts[save][ii] = bpts[d][ii];
      }  
    }
    /* end of insert knot */

    /* degree elevate bezier */
    for (i = lbz; i <= ph; i++)
    {
      for (ii = 0; ii < mc; ii++)
        ebpts[i][ii] = 0.0;
      mpi = min(d, i);
      for (j = max(0,i-t); j <= mpi; j++)
        for (ii = 0; ii < mc; ii++)
          ebpts[i][ii] = ebpts[i][ii] + bezalfs[i][j]*bpts[j][ii];
    }
    /* end of degree elevating bezier */

    if (oldr > 1)
    {
      /* must remove knot u=k[a] oldr times */
      first = kind - 2;
      last = kind;
      den = ub - ua;
      bet = (ub-ik[kind-1]) / den;
      
      /* knot removal loop */
      for (tr = 1; tr < oldr; tr++)
      {        
        i = first;
        j = last;
        kj = j - kind + 1;
        while (j - i > tr)
        {
          /* loop and compute the new control points */
	     /* for one removal step    */
          if (i < cind)
          {
            alf = (ub-ik[i])/(ua-ik[i]);
            for (ii = 0; ii < mc; ii++)
              ictrl[i][ii] = alf * ictrl[i][ii] + (1.0-alf) * ictrl[i-1][ii];
          }        
          if (j >= lbz)
          {
            if (j-tr <= kind-ph+oldr)
            {  
              gam = (ub-ik[j-tr]) / den;
              for (ii = 0; ii < mc; ii++)
                ebpts[kj][ii] = gam*ebpts[kj][ii] + (1.0-gam)*ebpts[kj+1][ii];
            }
            else
            {
              for (ii = 0; ii < mc; ii++)
                ebpts[kj][ii] = bet*ebpts[kj][ii] + (1.0-bet)*ebpts[kj+1][ii];
            }
          }
          i++;
          j--;
          kj--;
        }      
        
        first--;
        last++;
      }                    
    }
    /* end of removing knot n=k[a] */

    /* load the knot ua  */
    if (a != d)
      for (i = 0; i < ph-oldr; i++)
      {
        ik[kind] = ua;
        kind++;
      }

    /* load ctrl pts into ic  */
    for (j = lbz; j <= rbz; j++)
    {
      for (ii = 0; ii < mc; ii++)
        ictrl[cind][ii] = ebpts[j][ii];
      cind++;
    }
    
    if (b < m)
    {
      /* setup for next pass thru loop  */
      for (j = 0; j < r; j++)
        for (ii = 0; ii < mc; ii++)
          bpts[j][ii] = Nextbpts[j][ii];
      for (j = r; j <= d; j++)
        for (ii = 0; ii < mc; ii++)
          bpts[j][ii] = ctrl[b-d+j][ii];
      a = b;
      b++;
      ua = ub;
    }
    else
      /* end knot  */
      for (i = 0; i <= ph; i++)
        ik[kind+i] = ub;
  }                  
  /* end while loop  */

  *nh = mh - ph - 1;

  freevec2mat(ctrl);
  freevec2mat(ictrl);
  freematrix(bezalfs);
  freematrix(bpts);
  freematrix(ebpts);
  freematrix(Nextbpts);
  mxFree(alfs);

  return(ierr);
}

/* Matlab gateway function  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mc, nc;
  int mk, nk;
  int mu, nu;

  int d, t;
  double *c, *k;

  int nic, nik;
  double *ic, *ik;

  int nh, i, n;
  double *inc, *ink;

  if (nrhs != 4)
    mexErrMsgTxt("Four inputs required\n");

  if (nlhs > 2)
    mexErrMsgTxt("Maximum of two outputs required\n");

  /* spline degree */
  d = (int) mxGetScalar(prhs[0]);

  /* control points */
  c = mxGetPr(prhs[1]);
  mc = mxGetM(prhs[1]);
  nc = mxGetN(prhs[1]);

  /* knot sequence */
  k = mxGetPr(prhs[2]);
  mk = mxGetM(prhs[2]);
  nk = mxGetN(prhs[2]);

  /* raise degree t times */
  t = (int) mxGetScalar(prhs[3]);

  /* allocate work space t times larger than original number */
  /* of control points and knots */
  inc = (double *) mxMalloc((t+1)*mc*nc*sizeof(double));
  ink = (double *) mxMalloc((t+1)*nk*sizeof(double));

  bspdegelev(d, c, mc, nc, k, nk, t, &nh, inc, ink);

  /* new coefficients and knot sequence  */
  nic = nh + 1;
  nik = nic + d + t + 1;
  plhs[0] = mxCreateDoubleMatrix(mc, nic, mxREAL);
  ic = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(mk, nik, mxREAL);
  ik = mxGetPr(plhs[1]);

  /* copy new control points */
  n = mc * nic;
  for (i = 0; i < n; i++)
    ic[i] = inc[i];   
  
  /* copy new knots */
  for (i = 0; i < nik; i++)
    ik[i] = ink[i];      

  /* free working space */
  mxFree(ink);
  mxFree(inc);
}
