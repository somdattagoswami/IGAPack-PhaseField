/* Basis Function. 

 INPUT: 

   i - knot span  ( from FindSpan() ) 
   u - parametric point 
   p - spline degree 
   U - knot sequence 
    
 OUTPUT: 

   N - Basis functions vector[p+1] 

 Algorithm A2.2 from 'The NURBS BOOK' pg70. 

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

void basisfun(int i, double u, int p, double *U, double *N)
{
  int j,r;
  double saved, temp;

  /* work space */
  double *left  = (double*) mxMalloc((p+1)*sizeof(double));
  double *right = (double*) mxMalloc((p+1)*sizeof(double));
  
  N[0] = 1.0;
  for (j = 1; j <= p; j++)
  {
    left[j]  = u - U[i+1-j];
    right[j] = U[i+j] - u;
    saved = 0.0;
    
    for (r = 0; r < j; r++)
    {
      temp = N[r] / (right[r+1] + left[j-r]);
      N[r] = saved + right[r+1] * temp;
      saved = left[j-r] * temp;
    } 

    N[j] = saved;
  }
  
  mxFree(left);
  mxFree(right);
}
