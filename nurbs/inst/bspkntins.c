/* Insert Knot into a B-Spline 

 INPUT: 

   d - spline degree             integer               
   c - control points            double  matrix(mc,nc)      
   k - knot sequence             double  vector(nk)    
   u - new knots                 double  vector(nu)    

 OUTPUT: 

   ic - new control points double  matrix(mc,nc+nu)    
   ik - new knot sequence  double  vector(nk+nu)       

 Modified version of Algorithm A5.4 from 'The NURBS BOOK' pg164. 

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

int bspkntins(int d, double *c, int mc, int nc, double *k, int nk, 
              double *u, int nu, double *ic, double *ik)
{
  int ierr = 0;
  int a, b, r, l, i, j, m, n, s, q, ind;
  double alfa;

  double **ctrl  = vec2mat(c, mc, nc);
  double **ictrl = vec2mat(ic, mc, nc+nu);

  n = nc - 1;
  r = nu - 1;

  m = n + d + 1;
  a = findspan(n, d, u[0], k);
  b = findspan(n, d, u[r], k);
  ++b;

  for (q = 0; q < mc; q++)
  {
    for (j = 0; j <= a-d; j++) ictrl[j][q] = ctrl[j][q];
    for (j = b-1; j <= n; j++) ictrl[j+r+1][q] = ctrl[j][q];
  }
  for (j = 0; j <= a; j++)   ik[j] = k[j];
  for (j = b+d; j <= m; j++) ik[j+r+1] = k[j];

  i = b + d - 1;
  s = b + d + r;
  for (j = r; j >= 0; j--)
  {
    while (u[j] <= k[i] && i > a)
    {
      for (q = 0; q < mc; q++)
        ictrl[s-d-1][q] = ctrl[i-d-1][q];
      ik[s] = k[i];
      --s;
      --i;
    }
    for (q = 0; q < mc; q++)
      ictrl[s-d-1][q] = ictrl[s-d][q];
    for (l = 1; l <= d; l++)
    {
      ind = s - d + l;
      alfa = ik[s+l] - u[j];
      if (fabs(alfa) == 0.0)
        for (q = 0; q < mc; q++)
          ictrl[ind-1][q] = ictrl[ind][q];
      else
      {
        alfa /= (ik[s+l] - k[i-d+l]);
        for (q = 0; q < mc; q++)
          ictrl[ind-1][q] = alfa*ictrl[ind-1][q]+(1.0-alfa)*ictrl[ind][q];
      }
    }

    ik[s] = u[j];
    --s;
  }

  freevec2mat(ctrl);
  freevec2mat(ictrl);

  return ierr;
}

/* Matlab gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mc, nc;
  int mk, nk;
  int mu, nu;

  int d;
  double *c, *k, *u;

  int nic, nik;
  double *ic, *ik;

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

  /* knots to be inserted */
  u = mxGetPr(prhs[3]);
  mu = mxGetM(prhs[3]);
  nu = mxGetN(prhs[3]);

  /* new coefficients and knot sequence */
  nic = nc + nu;
  nik = nk + nu;
  plhs[0] = mxCreateDoubleMatrix(mc, nic, mxREAL);
  ic = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(mk, nik, mxREAL);
  ik = mxGetPr(plhs[1]);

  bspkntins(d, c, mc, nc, k, nk, u, nu, ic, ik);
}
