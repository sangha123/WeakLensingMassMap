#ifndef DATAPBL_H
#define DATAPBL_H

#include <vector.h>
#include <string>

using namespace std;

class datapbl{
 public:
  datapbl();
  datapbl(int npts_new);
  datapbl(vector<double> x_input, vector<double> y_input, vector<double> e1_input,
	  vector<double> e2_input, vector<double> w_input);
  double get_data(string name, int i);
  vector<double> get_data(string name);
  vector<double> compute_density();
 private:
  int npts;
  vector<double> x;
  vector<double> y;
  vector<double> e1;
  vector<double> e2;
  vector<double> w;
  double w_deriv(double q);
};

#endif
