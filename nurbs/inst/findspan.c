/* Find the knot span of the parametric point u. 

 INPUT: 

   n - number of control points - 1 
   p - spline degree                
   u - parametric point             
   U - knot sequence                

 RETURN: 

   s - knot span 

 Algorithm A2.1 from 'The NURBS BOOK' pg68 

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

int findspan(int n, int p, double u, double *U)
{
  int low, high, mid;
   
  /* special case */
  if (u == U[n+1]) return(n);
    
  /* do binary search */
  low = p;
  high = n + 1;
  mid = (low + high) / 2;
  while (u < U[mid] || u >= U[mid+1])
  {
    if (u < U[mid])
      high = mid;
    else
      low = mid;
    mid = (low + high) / 2;
  }  

  return(mid);
}

