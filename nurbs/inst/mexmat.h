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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __MEXMAT_H
#define __MEXMAT_H

#include "math.h"
#include "mex.h"

double **vec2mat(double *vec, int nrows, int ncols);
double **matrix(int nrows, int ncols);
void freevec2mat(double **mat);
void freematrix(double **mat);

#endif /* __MEXMAT_H */

