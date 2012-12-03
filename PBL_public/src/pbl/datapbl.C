#include <iostream>
#include "datapbl.h"
#include <cmath>
#include <vector.h>
#include "grid_field.h"
#define pi 3.14

datapbl::datapbl(){
  npts=1;
  e1.resize(1);
  e2.resize(1);
  w.resize(1);
}
datapbl::datapbl(int npts_new) {
  npts=npts_new;
  e1.resize(npts);
  e2.resize(npts);
  w.resize(npts);
}

datapbl::datapbl(vector<double> x_input, vector<double> y_input, vector<double> e1_input,vector<double> e2_input, vector<double> w_input) {
  npts=x_input.size();
  for (int i=0; i < npts; i++) {
    x.push_back(x_input[i]);
    y.push_back(y_input[i]);
    e1.push_back(e1_input[i]);
    e2.push_back(e2_input[i]);
    w.push_back(w_input[i]);
  }
  //compute rho...fill in w.  
}

double datapbl::get_data(string name, int i) {
  // 0 returns w
  // 1 returns e1
  // 2 returns e2
  if ((i >= 0)&&(i < npts)) {
  if (name.find("W") != string::npos) {
      return(w[i]);
    }
    else if (name.find("E1") != string::npos) {
      return(e1[i]);
    } 
    else if (name.find("E2") != string::npos) {
      return(e2[i]);
    }
  }
  return 0;

}

vector<double> datapbl::get_data(string name){
  if (name.find("W") != string::npos) {
    return w;
  }else if (name.find("E1") != string::npos){
    return e1;
  }else if (name.find("E2") != string::npos){
    return e2;
  } 
}

double datapbl::w_deriv(double q){
/* q=r/h
this an sph equation.
Review on SPH by J.J Monghan
eqn 2.6*/
  if (q < 1.) return (15./(14.*pi))*(pow((2.-q),3)-4.*pow((1.-q),3));
  if (q < 2.) return (15./(14.*pi))*pow((2.-q),3);
}

vector<double> datapbl::compute_density() {
  vector<double> r2,rho;
  double rho_ave=0.;
  for (int i=0; i < npts; i++) {
    r2.push_back(0.);
  }
  for (int i=0; i < npts; i++) {
    for (int j=0; j < npts; j++) {
      r2[j]=pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
    }
    sort(r2.begin(),r2.end());
    //rho.push_back(1/sqrt(r2[4]*r2[3]+0.00005)); // 3rd nearest neighbor...for simulations
    rho.push_back(1/r2[4]);
    rho_ave+=rho[i]/npts;
  }
  for (int i=0; i < npts; i++) {
    rho[i]=rho[i]/rho_ave;
  }
  return rho;


}


/*
vector<double> datapbl::compute_density() {
  // First, estimate the "typical" length scale.
  double rho0=1./npts;
  //double h=1./sqrt(rho0);
  double h=0.17;
  vector<double> rho_arr;
  double rho_ave=0;
  for (int i=0; i < npts; i++) {
    double rho=0.;
    for (int j=0; j < npts; j++) {
      double q=sqrt(pow((x[i]-x[j]),2)+pow((y[i]-y[j]),2))/h;
      rho+=2*pi*q*w_deriv(q)/h;
    }
    rho_arr.push_back(rho);
    rho_ave+=rho_arr[i]/npts;
  }
  for (int i=0; i < npts; i++) {
    rho_arr[i]=rho_arr[i]/rho_ave;
  }
  return rho_arr;
}

*/
