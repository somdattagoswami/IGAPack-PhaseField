/* Useful functions for accessing and manipulating matlab data structures.
   ======================================================================= 

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

/* convert c vector to c matrix. */
double **vec2mat(double *vec, int nrows, int ncols) 
{
  int col;
  double **mat;

  mat = (double**) mxMalloc (ncols*sizeof(double*));
  mat[0] = vec;
  for (col = 1; col < ncols; col++)
    mat[col] = mat[col-1] + nrows;  
  return mat;
}

/* create a new c matrix */
double **matrix(int nrows, int ncols) 
{
  int col;
  double **mat;

  mat = (double**) mxMalloc (ncols*sizeof(double*));
  mat[0] = (double*) mxMalloc (nrows*ncols*sizeof(double));
  for (col = 1; col < ncols; col++)
    mat[col] = mat[col-1] + nrows;  
  return mat;
}

void freevec2mat(double **mat)
{
  mxFree(mat);
}

void freematrix(double **mat)
{
  mxFree(mat[0]);
  mxFree(mat);
}


