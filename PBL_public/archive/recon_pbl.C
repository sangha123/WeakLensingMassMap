#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "datapbl.h"
#include "modelpbl.h"
#include "invert.h"
#include "grid_field.h"
#include "binnedgrid.h"
#include "time.h"

int main(int argc, char* argv[]) {
  // This code will take a configuration file of the form:
  // 1. datafile (containing column delimited x y e1 e2 sigma_e)
  // 2. int constsigma (0 to use the parameters in the datafile, 
  //                    1 to use the parameter in line 3)
  // 3. double sigma_e (if constsigma=0, then this number is thrown away)
  // 4. Eta
  // 4a. niter
  // 5. output filename
  // 6. idx associated with h (e.g. h=h1*rho0^(-0.5)/(rho/rho0)^idx), idx=0.5 means constant noise
  //                                             idx=0.=constant length
  // 7. h1, the length scale associated with a density of "1".  
  // 8. The number of modeled isothermal spheres assumed for the initial conditions
  // 9-11, 12-14, ..., xc, yc, theta_e for each of the isothermal spheres
if (argc == 1) {
    cerr << "calling procedure is: kappapbl [configfile] \n";
    return(0);
  } 
 
 fstream datafile, configfile;
 configfile.open(argv[1],ios::in);
 if (! configfile.is_open()) {
   cerr << "Boo!  File does not exist.\n";
   return(0);
 } 

 char dataname[200],outpblfile[200];
 
 vector<double> x,y,w,e1,e2,psi0,ktrue;
 configfile >> dataname;
 datafile.open(dataname,ios::in);
 int constsigma,n_iso,niter;
 double sige,h1,idx,eta;
 
 configfile >> constsigma;
 configfile >> sige;
 configfile >> eta;
 configfile >> niter;
 configfile >> outpblfile;
 configfile >> idx;
 configfile >> h1;

 int n=0;
 
 string line;
 double d1,d2,d3,d4,d5,d6,d7,d8;
 while  (! (datafile.eof() || datafile.peek()== EOF)) {
   datafile >> d1 >> d2 >> d3 >> d4 >> d5;  // >> d5 >> d6 >> d7 >> d8;
   x.push_back(d1);
   y.push_back(d2);
   e1.push_back(d3);
   e2.push_back(d4);
   ktrue.push_back(d5);
   w.push_back(100);
   //w.push_back(1./(d5*d5));
   //psi0.push_back(d7);
   if (constsigma == 1) {
     w[n]=1./(sige*sige);
   }
   getline(datafile,line);
    n++;
  }
 datafile.close();
 cerr << "n=" << n << "\n";

 // Now, do we want to start the initial conditions with an isothermal sphere (or two)?
 configfile >> n_iso;
 cerr << "initializing with " << n_iso << " isothermal spheres\n";
 vector<double>psi_new;
 for (int i=0; i < n; i++) {
   psi_new.push_back(0.);
 }
 double xc,yc,te;
 for (int i_iso=0; i_iso < n_iso; i_iso++) {
   configfile >> xc;
   configfile >> yc;
   configfile >> te;
   for (int i=0; i < n; i++) {
     double theta=sqrt(pow(x[i]-xc,2)+pow(y[i]-yc,2));
     psi_new[i]+=te*theta;
   }
 }
 configfile.close();


 // Initialize the data
 datapbl mydata(x,y,e1,e2,w);
 vector<double> rho=mydata.compute_density();
 /*for(int i=0;i<n;i++){
   cout<<x[i]<<"  "<<y[i]<<"  "<<rho[i]<<endl;
   }*/
 double rho0=n;
 int nx0=20;
 vector<double> rho_grid;
 if (idx==1){
   grid_field rhobin(x,y,rho,nx0,nx0);
   rho_grid=rhobin.bin_pbl();
   for(int i=0; i < n; i++) {
     int i0=floor(x[i]*nx0);
     int j0=floor(y[i]*nx0);
     double wx=abs(x[i]*nx0-i0);
     double wy=abs(y[i]*nx0-j0);
     rho[i]=rho_grid[i0*nx0+j0]*(1-wx)*(1-wy)+rho_grid[i0*nx0+j0-1]*(1-wx)*wy
       +rho_grid[(i0-1)*nx0+j0]*wx*(1-wy)+rho_grid[(i0-1)*nx0+j0-1]*wx*wy;
   }
 }
 vector<double> h;
 for (int i=0; i < n; i++) {
   h.push_back(h1/sqrt(rho0)/pow(rho[i],idx));
 }
 modelpbl mymodel(x,y,h,3);
 mymodel.set_data(mydata);
 mymodel.set_eta_step(30.);
 mymodel.set_eta(eta);
 mymodel.set_psi(psi_new);
 mymodel.update_all();

 //***************************************************************
 // 
 // At this point, both the model and the data are configured
 // and initialized.
 // Below, we minimize chi^2
 //
 //***************************************************************


 double dof=n-1.;
 cerr << "chi2= " << mymodel.chi2()/dof << "\n";
 
 int npts=mymodel.get_npts();
 Array1D<double> dpsi(n);
 for(int j=0;j<n;j++){
     dpsi[j]=0;
 }
 double chi0;
 time_t seconds1,seconds2;
 double factor;
 factor=0.45;
 double dchi=1000;
 double eps=0.01;
 double chi=mymodel.chi2()/dof;
 while (dchi > eps){
   seconds1=time(NULL);
   
   dpsi=mymodel.psi_newton_pbl();
   for(int j=0;j<n;j++){
     psi_new[j]+=factor*dpsi[j];
   }
   seconds2=time(NULL);
   cerr<<"time for each iteration= "<<(seconds2-seconds1)<< " sec.\n";
   mymodel.set_psi(psi_new);
   cerr << "chi2= " << mymodel.chi2()/dof << "\n";
   chi0=mymodel.chi2()/dof;
   dchi=(chi-chi0);
   chi=chi0;
     
 }

 mymodel.set_eta(0.);
 cerr<<"final chi^2=  "<<mymodel.chi2()/dof<<endl;
 
 //****************************************************************
 // We have completed our minimization.  Now, time to export
 //****************************************************************


 // First, grid up the data...
 vector<double> gamma1;
 vector<double> gamma2;
 vector<double> kappa;
 
 gamma1=mymodel.get_data("GAMMA1");
 gamma2=mymodel.get_data("GAMMA2");
 kappa=mymodel.get_data("KAPPA");

 int nx=40;
 int ny=40;

 vector<double> g1_bin;
 vector<double> g2_bin;
 vector<double> k_bin;
 vector<double> ktrue_bin;
 vector<int> n_grid;
 binnedgrid kapbin(x,y,kappa,nx,ny);
 k_bin=kapbin.get_data();
 binnedgrid gam1bin(x,y,gamma1,nx,ny);
 g1_bin=gam1bin.get_data();
 binnedgrid gam2bin(x,y,gamma2,nx,ny);
 g2_bin=gam2bin.get_data();
 vector<double> weight=gam2bin.get_w();
 n_grid=gam2bin.get_count();

 binnedgrid ktruebin(x,y,ktrue,nx,ny);
 ktrue_bin=ktruebin.get_data();

 ofstream outbindata, outptdata;

 

 outbindata.open("gridded.dat");
 outptdata.open(outpblfile);

 

 // Print out the kappa field, g1, and g2 fields

 for(int i=0;i<nx;i++){
   for(int j=0;j<ny;j++){
     outbindata <<i<<"  "<<j<<"  "<<k_bin[i*ny+j]<<" " << ktrue_bin[i*ny+j]<<"  "<<g1_bin[i*ny+j]<<"  "<<g2_bin[i*ny+j]
		<<"  "<<n_grid[i*ny+j]<< endl;
   }
 }
 outbindata.close();
 
 // Finally, compute the best fit statistics

 // 1 Compute the nearest gridcell and compute delta kappa^2
  double stat1=0;
  double stat2=0;
  double stat3=0;
  double Sxx=0,Sxy=0,Sy=0,Sx=0;
  for (int i=0; i < n; i++) {
    double x=kappa[i];
    double y=ktrue[i];
    Sxx+=x*x;
    Sxy+=x*y;
    Sx+=x;
    Sy+=y;
  }
  
  double lambda=(Sxy-Sy-Sx+n)/(Sxx-2*Sx+n);
  for (int i=0; i < n; i++) {
    double x=kappa[i];
    double y=ktrue[i];
    stat1+=pow(y-(lambda*x+(1-lambda)),2)/n;
  }

  vector<double> kappa_pbl=mymodel.get_data("KAPPA");

  binnedgrid kbin(x,y,kappa_pbl,nx,ny);
  k_bin=kbin.get_data();

  Sxx=0;Sxy=0;Sy=0;Sx=0;
  for (int i=0; i < nx*nx; i++) {
    double x=k_bin[i];
    double y=ktrue_bin[i];
    Sxx+=x*x;
    Sxy+=x*y;
    Sx+=x;
    Sy+=y;
  }
  
  double lambda_grid=(Sxy-Sy-Sx+nx*nx)/(Sxx-2*Sx+nx*nx);
  int c=0;
  for (int i=0; i < nx*nx; i++) {
    double x=k_bin[i];
    double y=ktrue_bin[i];
    if (n_grid[i] > 0) {
      stat2+=pow(y-(lambda_grid*x+(1-lambda_grid)),2);
      c=c+1;
    }
    stat3+=pow(y-(lambda_grid*x+(1-lambda_grid)),2)*n_grid[i]/n;
  }

  stat2=stat2/c;

  cerr << "Awesomeness statistic1(particle) = " << stat1 << "\n";
  cerr << "Awesomeness statistic2 = " << stat2 << "\n";
 cerr << "Awesomeness statistic3(weighted bin) = " << stat3 << "\n";

  // Print out the individual data points
 for(int i=0;i<npts;i++){
   outptdata<<x[i]<<"  "<<y[i]<<"  "<<mymodel.get_data("PSI",i)<<"  "
	    <<kappa[i]*lambda+(1-lambda)<<"  "<<y[i] <<"\n";
 }
 outptdata.close();

 return(0);
   
}


 

 
 
