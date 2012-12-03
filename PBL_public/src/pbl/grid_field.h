#include <vector.h>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

class grid_field{
 public:
  grid_field(vector<double> x_new, vector<double> y_new, vector<double> field_new,int nx_new,int ny_new);
  vector<double> bin_pbl();
  vector<int> get_n_grid(){return n_grid;};
  double min(vector<double> x);
  double max(vector<double> x);
 private:
  int nx;
  int ny;
  int npts;
  vector<double> psi_grid;
  vector<int> n_grid;
  vector<double> x;
  vector<double> y;
  vector<double> field;
  //Later on write a template that does these...for vectors,arrays,lists etc.
  
};
