#include "binnedgrid.h"
#include <vector.h>
#include <cmath>

using namespace std;

binnedgrid::binnedgrid() {
  nx=1;
  ny=1;
}

binnedgrid::binnedgrid(vector<double> x, vector<double> y, 
		       vector<double> data_arr, int nx_new,int ny_new){
  // Assume the domain is 0->1 in both directions
  nx=nx_new;
  ny=ny_new;
  
  for (int i=0; i < nx*ny; i++) {
    data.push_back(0);
    w.push_back(0);
    count.push_back(0);
  }
  
  int nelem=x.size();
  int i,j;
  double wx,wy;
  double dx=max(x)-min(x);
  double dy=max(y)-min(y);
  for (int idx=0; idx < nelem; idx++) {
    i=floor(x[idx]/dx*nx);
    j=floor(y[idx]/dx*nx);
    wx=abs(x[idx]*nx-i);
    wy=abs(y[idx]*nx-j);
    if ((i > 0)&&(i < nx)&&(j > 0)&&(j < ny)){
      data[i*ny+j]+=data_arr[idx]*(1-wx)*(1-wy);
      data[(i-1)*ny+j]+=data_arr[idx]*wx*(1-wy);
      data[i*ny+j-1]+=data_arr[idx]*(1-wx)*wy;
      data[(i-1)*ny+j-1]+=data_arr[idx]*wx*wy;
      w[i*ny+j]+=(1-wx)*(1-wy);
      w[(i-1)*ny+j]+=wx*(1-wy);
      w[i*ny+j-1]+=(1-wx)*wy;
      w[(i-1)*ny+j-1]+=wx*wy;
    }
    count[i*ny+j]++; 
  }
  for (int i=0; i < nx*ny; i++) {
    if (w[i] > 0) {
      data[i]=data[i]/w[i];
    }
  }
}

double binnedgrid::min(vector<double> x){
  double min_val=0;
  int npts=x.size();
  for(int i=0;i<npts;i++){
    if(x[i] < min_val) min_val=x[i];

  }
  return min_val;
}

double binnedgrid::max(vector<double> x){
  double max_val=0;
  int npts=x.size();
  for(int i=0;i<npts;i++){
    if(x[i] > max_val) max_val=x[i];

  }
return max_val;
}
