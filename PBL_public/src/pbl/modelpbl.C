#include "modelpbl.h"
#include "datapbl.h"
#include <cmath>
#include "tnt.h"
#include "jama_lu.h"
#include "invert.h"
#include "binnedgrid.h"
#define pi 3.14
using namespace std;

modelpbl::modelpbl()
{
  npts=1;
  eta=0.0001;
  eta_step=0.;
  psi.resize(1);
  kappa.resize(1);
  gamma1.resize(1);
  gamma2.resize(1);
  e1.resize(1);
  e2.resize(1);
  e1_data.resize(1);
  e2_data.resize(1);
  
}

modelpbl::modelpbl(vector<double> x_new, vector<double> y_new, 
		   vector<double> h_new, int order_new)
{
  eta=1.;
  eta_step=0.;
  npts=x_new.size();
  x=x_new;
  y=y_new;
  h=h_new;
  order=order_new;
  
  for (int i=0; i < npts; i++) {
    psi.push_back(0);
    kappa.push_back(0);
    gamma1.push_back(0);
    gamma2.push_back(0);
    e1.push_back(0);
    e2.push_back(0);
    e1_data.push_back(0);
    e2_data.push_back(0);
  }
  create_matrix();
}   


void modelpbl::set_eta(double val){
  eta=val;
  
}


double modelpbl::get_data(string name,int i) {
  // If we like, we can say that element -1, and nx yields a value... 0.
  // all names are in uppercase
 
  if (name.find("PSI") != string::npos) {
    if (i>0 && i<npts) {
      return(psi[i]);
    } else {
      return 0;
    }
  }
  if (name.find("KAPPA") != string::npos) {
    if (i>0 && i<npts) {
      return(kappa[i]);
    } else {
      return 0;
    }
  }
  if (name.find("GAMMA1") != string::npos) {
    if (i>0 && i<npts) {
      return(gamma1[i]);
    } else {
      return 0;
    }
    
  }
  if (name.find("GAMMA2") != string::npos) {
    if (i>0 && i<npts) {
      return(gamma2[i]);
    } else {
      return 0;
    }
  }

  if (name.find("E1") != string::npos) {
    if (i>0 && i<npts) {
      return(e1[i]);
    } else {
      return 0;
    }
  }

  if (name.find("E2") != string::npos) {
    if (i>0 && i<npts) {
      return(e2[i]);
    } else {
      return 0;
    }
  }
  cerr << "Warning! modelpbl.get_data called without correct datatype\n";
  return 0;
}

vector<double> modelpbl::get_data(string name)  {
  if (name == "PSI")
    return psi;
  else if (name == "KAPPA")
    return kappa;
  else if (name == "GAMMA1") 
    return gamma1;
   else if (name == "GAMMA2" ) 
    return gamma2;
   else if (name == "E1") 
    return e1; 
   else if (name == "E2") 
    return e2; 
  
   else {
    cerr << "Type does not exist\n";
    return vector<double>(0);
  }
}

void modelpbl::set_psi(int i,double val)
{
  //cerr<<"inside set_psi"<<psi.size()<<endl;
  if ((i>=0)&&(i < npts)) {
    psi[i]=val;
  } else {
    cerr << i <<  " outside the domain\n";
  }
}

void modelpbl::set_psi(vector<double> valarray) {
  if (valarray.size() == psi.size()) {
    psi=valarray;
    
  }
  update_all();
}

void modelpbl::create_matrix() {
  int pars=9;
  // For now, the order may only be 2 or 3.
  if (order==2) {
    pars=5;
  } else {
    order=3;
  }

  //  cerr << "pars= " << pars << "\n";
  // cerr << "npts= " << npts << "\n";

  // First, create the blank A array;
  
  vector<double> dum;
  for (int i=0; i < npts; i++) {
    dum.push_back(0.);
  }
  vector<vector<double> > dum1;
  for (int j=0; j < npts; j++) {
    dum1.push_back(dum);
  }
  kmat=dum1;
  g1mat=dum1;
  g2mat=dum1;
  // dum1 has dimension npts*npts
  for (int j=0; j < pars; j++) {
    D.push_back(dum1);
  }
  // D now has the appropriate dimensionality (pars,npts,npts);


  vector<vector <double> > X_arr, D_out;
  vector<double> W;
  W=dum;
  Array2D<double> A_arr (pars,pars),A_inv(pars,pars);
  // Allocate memory for X_arr
  for (int i=0; i < pars; i++) {
    X_arr.push_back(dum);
    D_out.push_back(dum);
  }
  
  
  for (int i=0; i< npts; i++) {
    
    for (int alpha=0; alpha < pars; alpha++) {
      for (int beta=0; beta < pars; beta++) {
	A_arr[alpha][beta]=0.;
      }
    }

    for (int j=0; j < npts; j++) {
      double dx=x[j]-x[i];
      double dy=y[j]-y[i];
      double drsq=(dx*dx+dy*dy)/(h[j]*h[j]);
      double dr=sqrt(drsq);
      W[j]=exp(-(pow(drsq,1)/2.));
      //if (dr < 1.) W[j]=(15./(14.*pi))*(pow((2.-dr),3)-4.*pow((1.-dr),3));
      //if (dr < 2.) W[j]=(15./(14.*pi))*pow((2.-dr),3);
      X_arr[0][j]=dx;
      X_arr[1][j]=dy;
      X_arr[2][j]=0.5*dx*dx;
      X_arr[3][j] = dx*dy;
      X_arr[4][j] = 0.5*dy*dy;
      if (order > 2){      
	X_arr[5][j] = 1./6. *pow(dx,3);
	X_arr[6][j] = 1./2. *dx*dx*dy;
	X_arr[7][j] = 1./2. *dx*dy*dy;
	X_arr[8][j] = 1./6. *pow(dy,3);
      }
      for (int alpha=0; alpha < pars; alpha++) {
	for (int beta=0; beta < pars; beta++) {
	  A_arr[alpha][beta]+=(X_arr[alpha][j]*X_arr[beta][j]*W[j]);
	}
      }
    }
    A_inv=invert(A_arr);
    for (int alpha=0;alpha < pars; alpha++) {
      for (int j=0; j < npts; j++) {
	D_out[alpha][j]=0;
	for (int beta=0; beta < pars; beta++) {
	  // Note: Make sure that A_inv is being called in the correct order
	  D_out[alpha][j]+=A_inv[alpha][beta]*X_arr[beta][j];
	}
      }
    }


    for (int alpha=0; alpha < pars; alpha++) {
      double tot=0.;
      for (int j=0; j < npts; j++) {
	D[alpha][i][j]=W[j]*D_out[alpha][j];
	tot+=D[alpha][i][j];
      }      
      D[alpha][i][i]-=tot;
    }
    // Now, create the specific matrices for kappa and gamma
    for (int j=0; j < npts; j++) {
      kmat[i][j]=(D[2][i][j]+D[4][i][j])/2.;
      g1mat[i][j]=(D[2][i][j]-D[4][i][j])/2.;
      g2mat[i][j]=(D[3][i][j]);
    }

  }
}

void modelpbl::update_all() {
  update_kappa();
  update_gamma();
  update_ellipticity();
}

void modelpbl::update_kappa(){
  for (int i=0; i < npts; i++) {
    kappa[i]=0.;
    for (int j=0; j < npts; j++) {
      kappa[i]+=psi[j]*kmat[i][j];
    }
  }
}

void modelpbl::update_gamma(){
  for (int i=0; i < npts; i++) {
    gamma1[i]=0.;
    gamma2[i]=0.;
    for (int j=0; j < npts; j++) {
      gamma1[i]+=psi[j]*g1mat[i][j];
      gamma2[i]+=psi[j]*g2mat[i][j];
    }
  }
} 
  



void modelpbl::update_ellipticity(){
  
  double g1,g2,g,u,H,gsq,gammasq;
  double epsilon=0.00001;
  for (int i=0; i < npts; i++) { 
    gsq=(pow(gamma1[i]/(1-kappa[i]),2)+pow(gamma2[i]/(1-kappa[i]),2))+epsilon;
    u=eta_step*((gsq*gsq-1.)/gsq);
    H=1./(1.+exp(-2*u));
 
    e1[i]=(1-H)*(gamma1[i]/(1-kappa[i]))+H*(gamma1[i]*(1-kappa[i]))/(gamma1[i]*gamma1[i]+gamma2[i]*gamma2[i]+epsilon);
    e2[i]=(1-H)*(gamma2[i]/(1-kappa[i]))+H*(gamma2[i]*(1-kappa[i]))/(gamma1[i]*gamma1[i]+gamma2[i]*gamma2[i]+epsilon);
    }
}

void modelpbl::set_eta_step(double value){
  eta_step=value;
}

// Compute chi^2 for a particular model
double modelpbl::chi2() {
  double chi2_out=0.;
  for (int i=0; i < npts; i++) {
    chi2_out+=w_data[i]*pow((e1[i]-e1_data[i]),2);
    chi2_out+=w_data[i]*pow((e2[i]-e2_data[i]),2);
    chi2_out+=eta*pow(kappa[i],2);
  } 
  return(chi2_out);
}



void modelpbl::set_data(datapbl d){

  e1_data=d.get_data("E1");
  e2_data=d.get_data("E2");
  w_data=d.get_data("W"); 
}

vector<vector<double> > modelpbl::get_matrix(string name)  {
  if (name.find("KAPPAMATRIX") != string::npos) {
    return kmat;
  } else if (name.find("GAMMA1MATRIX") != string::npos) {
    return g1mat;
  } else if (name.find("GAMMA2MATRIX") != string::npos) {
    return g2mat;
  }

}

Array1D<double> modelpbl::psi_newton_pbl()
{
  Array1D<double> b(npts),dpsi(npts);
  Array2D<double> M(npts,npts);
  Array2D<double> A(3*npts,npts);//number of constraints=3*npts;  //e1,e2,kappa
  vector<double> dyarr(3*npts,0);
  
  double epsilon=0.000001;
  double de1dgamma1,de1dgamma2,de1dkappa,de2dgamma1,de2dgamma2,de2dkappa,g1,g2,u,H;
 
  double dudkappa,dudgamma1,dudgamma2,dHdkappa,dHdgamma1,dHdgamma2,gsq,gammasq;
  //cerr<<"Before any loop"<<time(NULL)<<endl;
  for(int i=0;i<npts;i++){
    b[i]=0.;
    dpsi[i]=0;
    for(int j=0;j<npts;j++){
      M[i][j]=0.;
      A[i][j]=0;
      A[i+npts][j]=0;
      A[i+2*npts][j]=0.;
    }
  }
  //cerr<<"After initialization"<<time(NULL)<<endl;
  for(int i=0;i<npts;i++){
    
    gsq=(pow((gamma1[i]/(1-kappa[i])),2)+pow((gamma2[i]/(1-kappa[i])),2))+epsilon;
    gammasq=gamma1[i]*gamma1[i]+gamma2[i]*gamma2[i]+epsilon;
    //The step function
    u=eta_step*((gsq*gsq-1.)/gsq);
    H=1./(1.+exp(-2*u));
    //cerr<<"eta_step"<<eta_step<<endl;
    //Derivatives of the step function
    dudkappa=(2*eta_step/(1-kappa[i])) *(gsq*gsq+1)/(gsq);
    dudgamma1=(2*eta_step*gamma1[i]/pow((1-kappa[i]),2))*(gsq*gsq+1)/(gsq*gsq);
    dudgamma2=(2*eta_step*gamma2[i]/pow((1-kappa[i]),2)) *(gsq*gsq+1)/(gsq*gsq);
    
    dHdkappa=(2./(2+exp(2*u)+exp(-2*u))) *dudkappa;
    dHdgamma1=(2./(2+exp(2*u)+exp(-2*u))) *dudgamma1;
    dHdgamma2=(2./(2+exp(2*u)+exp(-2*u))) *dudgamma2;

    dyarr[i]=e1_data[i]-e1[i];
    dyarr[i+npts]=e2_data[i]-e2[i];
    dyarr[i+2*npts]=-kappa[i];

    //Derivatives of smoothed ellipticities

    de1dkappa=gamma1[i]*(1-H)/pow((1-kappa[i]),2)
              -H*gamma1[i]/(gammasq)
              -gamma1[i]*dHdkappa/(1-kappa[i])
              +gamma1[i]*(1-kappa[i])*dHdkappa/(gammasq);

    de1dgamma1=(1-H)/(1-kappa[i])
                -2*pow(gamma1[i],2)*(1-kappa[i])*H/pow(gammasq,2)
                +(1-kappa[i])*H/(gammasq)
                -gamma1[i]*dHdgamma1/(1-kappa[i])+gamma1[i]*(1-kappa[i])*dHdgamma1/(gammasq);

    de1dgamma2=-2*gamma1[i]*gamma2[i]*(1-kappa[i])*H/pow(gammasq,2)
               -gamma1[i]*dHdgamma2/(1-kappa[i])
               +gamma1[i]*(1-kappa[i])*dHdgamma2/gammasq;

    de2dkappa=gamma2[i]*(1-H)/pow((1-kappa[i]),2)
              -gamma2[i]*H/gammasq
              -gamma2[i]*dHdkappa/(1-kappa[i])
              +gamma2[i]*(1-kappa[i])*dHdkappa/gammasq;

    de2dgamma1=-2*gamma1[i]*gamma2[i]*(1-kappa[i])*H/pow(gammasq,2)
               -gamma2[i]*dHdgamma1/(1-kappa[i])
               +gamma2[i]*(1-kappa[i])*dHdgamma1/gammasq;

    de2dgamma2=(1-H)/(1-kappa[i])
               -2*pow(gamma2[i],2)*(1-kappa[i])*H/pow(gammasq,2)
               +(1-kappa[i])*H/gammasq-gamma2[i]*dHdgamma2/(1-kappa[i])
               +gamma2[i]*(1-kappa[i])*dHdgamma2/gammasq;

    
    for(int j=0;j<npts;j++){
      A[i][j]=kmat[i][j]*de1dkappa+g1mat[i][j]*de1dgamma1+g2mat[i][j]*de1dgamma2;
      A[i+npts][j]=kmat[i][j]*de2dkappa+g1mat[i][j]*de2dgamma1+g2mat[i][j]*de2dgamma2;
      A[i+2*npts][j]=kmat[i][j];
      b[j]+=w_data[i]*A[i][j]*dyarr[i]+w_data[i]*A[i+npts][j]*dyarr[i+npts]+
	eta*A[i+2*npts][j]*dyarr[i+2*npts];
      for(int k=0;k<=j;k++){
	M[j][k]+=w_data[i]*A[i][j]*A[i][k]+w_data[i]*A[i+npts][j]*A[i+npts][k]+
	  eta*A[i+2*npts][j]*A[i+2*npts][k];
	M[k][j]=M[j][k];
	if (j == k)  M[j][k]+=1.e-5;
      }
    }
  }
  //cerr<<"Before Matrix inversion"<<time(NULL)<<endl;
 JAMA::LU<double> A_lu(M);
 dpsi=A_lu.solve(b);
 //cerr<<"eta = "<<eta<<endl;
 // cerr<<"After Matrix inversion"<<time(NULL)<<endl;
 return dpsi;
}

