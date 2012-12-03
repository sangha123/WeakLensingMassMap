#include <iostream>
#include <vector.h>

class binnedgrid {
 public:
  binnedgrid();
  binnedgrid(vector<double> x, vector<double> y, vector<double> data_arr, 
	     int nx,int ny);
  vector<double> get_data(){return data;};
  vector<double> get_w(){return w;};
  vector<int> get_count(){return count;};
 private:
  int nx;
  int ny;
  vector<double> w;
  vector<double> data;
  vector<int> count;
  double max(vector<double> );
  double min(vector<double> );
};

