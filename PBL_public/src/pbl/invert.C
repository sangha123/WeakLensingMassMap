#include "invert.h"
#include "tnt.h"
#include "jama_lu.h"
#include <iostream>

using namespace TNT;

Array2D<double> invert(Array2D<double> A){
  
  if(A.dim1()!=A.dim2()) cerr<<"Cannot invert, Matrix is not square."<<endl;

  JAMA::LU<double> A_lu(A);
  TNT::Array2D<double> id(A.dim1(),A.dim2(),double(0));
  for (int i = 0; i < A.dim1(); i++) id[i][i] = 1;
  return A_lu.solve(id);
}
