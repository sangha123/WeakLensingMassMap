#include <vector.h>
#include <iostream>
#include <string>
#include "grid_field.h"

using namespace std;

grid_field::grid_field(vector<double> x_new,vector<double> y_new,vector<double> field_new,int nx_new,int ny_new){
  nx=nx_new;
  ny=ny_new;
  npts=x_new.size();
 
  for(int i=0;i<nx*ny;i++){
    psi_grid.push_back(0);
    n_grid.push_back(0);
  }
  
  for(int i=0;i<npts;i++){
    x.push_back(x_new[i]);
    y.push_back(y_new[i]);
    field.push_back(field_new[i]);
  }
 
}

vector<double> grid_field::bin_pbl(){
  double wtot;
  double xmin=min(x);
  double ymin=min(y);
  double xmax=max(x);
  double ymax=max(y);
  
  double dx=(xmax-xmin)/nx;
  double dy=dx;

  double h=sqrt((xmax-xmin)*(ymax-ymin)/npts);
  //cerr<<xmin<<"  "<<ymin<<"  "<<xmax<<"  "<<ymax<<"  "<<h<<endl;
  
  for(int i=0;i<nx;i++){
    double xc=xmin+(i+0.5)*dx;
    for(int j=0;j<ny;j++){
      double yc=ymin+(j+0.5)*dy;
      wtot=0;
	  //The weighting
      for(int k=0;k<npts;k++){
	double W=exp(-(pow((x[k]-xc),2)+pow((y[k]-yc),2))/(2*h*h));
	wtot+=W;
	psi_grid[i*ny+j]+=W*field[k];
      }
      psi_grid[i*ny+j]=psi_grid[i*ny+j]/wtot;
      
    } 
  }
  
  int i0,j0;
  for(int k=0;k<npts;k++){
    i0=floor(x[k] *nx);
    j0=floor(y[k] *ny);
    
    n_grid[i0*ny+j0]+=1;
  }
  return psi_grid;
}



double grid_field::min(vector<double> x){
  double min_val=0;
  for(int i=0;i<npts;i++){
    if(x[i] < min_val) min_val=x[i];

  }
  return min_val;
}

double grid_field::max(vector<double> x){
  double max_val=0;
  for(int i=0;i<npts;i++){
    if(x[i] > max_val) max_val=x[i];

  }
return max_val;
}
