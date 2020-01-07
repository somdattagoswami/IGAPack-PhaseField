/* Evaluation of univariate B-Spline.

 MATLAB SYNTAX: 

        p = bspeval(d,c,k,u) 
 
 INPUT: 

   d - degree of the B-Spline  integer 
   c - control points          double  matrix(mc,nc) 
   k - knot sequence           double  vector(nk)    
   u - parametric points       double  vector(nu)    

 OUTPUT: 

   p - evaluated points        double matrix(mc,nu)  

 Modified version of Algorithm A3.1 from 'The NURBS BOOK' pg82. 

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

int bspeval(int d, double *c, int mc, int nc, double *k, int nk, double *u,
            int nu, double *p)
{
  int ierr = 0;
  int i, s, tmp1, row, col;
  double tmp2;

  /* Construct the control points */
  double **ctrl = vec2mat(c,mc,nc);

  /* Contruct the evaluated points */
  double **pnt = vec2mat(p,mc,nu);
 
  /* space for the basis functions */
  double *N = (double*) mxMalloc((d+1)*sizeof(double));

  /* for each parametric point i */
  for (col = 0; col < nu; col++)
  {
    /* find the span of u[col] */
    s = findspan(nc-1, d, u[col], k);
    basisfun(s, u[col], d, k, N);
    
    tmp1 = s - d;
    for (row = 0; row < mc; row++)
    {
      tmp2 = 0.0;   
      for (i = 0; i <= d; i++)
	tmp2 += N[i] * ctrl[tmp1+i][row];
      pnt[col][row] = tmp2;
    }
  }

  mxFree(N);
  freevec2mat(pnt);
  freevec2mat(ctrl);

  return ierr;
} 

/* Matlab gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mc, nc;
  int mk, nk;
  int mu, nu;

  int d;
  double *c, *p, *k, *u;

  if (nrhs != 4)
    mexErrMsgTxt("Four inputs required\n");

  if (nlhs > 1)
    mexErrMsgTxt("One output required\n");

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

  /* parametric points */
  u = mxGetPr(prhs[3]);
  mu = mxGetM(prhs[3]);
  nu = mxGetN(prhs[3]);

  /* evaluated points */
  plhs[0] = mxCreateDoubleMatrix(mc, nu, mxREAL);
  p = mxGetPr(plhs[0]);

  bspeval(d, c, mc, nc, k, nk, u, nu, p);
}



