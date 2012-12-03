#include <vector.h>
#include <iostream>
#include <string>
#include "tnt.h"
#include "jama_lu.h"
#include "invert.h"

class datapbl;

using namespace std;

class modelpbl{
 public:
  modelpbl();
  modelpbl(vector<double> x, vector<double> y, vector<double> h,int order=3);
  double get_data(string name,int i);
  void set_eta(double val);
  int get_npts(){return npts;};
  double get_eta(){return eta;};
  vector<vector<double> > get_matrix(string name);
  vector<double> get_data(string name);
  void set_psi(int i, double val);
  void set_psi(vector<double> valarray);
  void create_matrix();
  void update_kappa();
  void update_gamma();
  void update_ellipticity();
  void update_all();
  void set_eta_step(double val);
  double get_eta_step(){return eta_step;};
  double chi2();
  void set_data(datapbl d);//This functions sets e1_data,e2_data to the measured ellipticity extracting them from datapbl
  vector<double> matmul(vector<vector<double> > A,vector<double> psi);
  Array1D<double> psi_newton_pbl();
 private:
  int npts;
  double eta;
  double eta_step;
  double w_weak;
  vector<double> x;
  vector<double> y;
  vector<double> h;
  int order;
  vector<double> e1;
  vector<double> e2;
  vector<double> w_data;
  vector<double> e1_data;
  vector<double> e2_data;
  vector<double> psi;
  vector<double> kappa;
  vector<double> gamma1;
  vector<double> gamma2;
  /*  The D object is as follows:
      D_0 = d/dx
      D_1 = d/dy
      D_2 = d^2/dx^2
      D_3 = d^2/dxdy
      D_4 = d^2/dy^2
      and so on...
  */
  vector<vector<vector<double> > > D;
  vector<vector<double> > kmat;
  vector<vector<double> > g1mat;
  vector<vector<double> > g2mat;
};





