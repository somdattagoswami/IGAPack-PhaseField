% Compile the C functions
  mex bspeval.c basisfun.c findspan.c mexmat.c
  mex bspderiv.c mexmat.c
  mex bspkntins.c findspan.c mexmat.c
  mex bspdegelev.c mexmat.c

% Add the NURBS toolbox to the path
  addpath (pwd);
  savepath()

