#include<exception>
#include<vector>
#include "defs.hpp"


void setBCs(int nx,int ny,int nvar,int nb,real**** Q) {
	for (int k=0; k<nvar; k++) {
                for (int i=0;i<nx;i++) {Q[i][0][k][nb]=Q[i][1][k][nb];
                                            Q[i][ny-1][k][nb]=Q[i][ny-2][k][nb];}
                for (int j=0;j<ny;j++) {Q[0][j][k][nb]=Q[1][j][k][nb];
                                            Q[nx-1][j][k][nb]=Q[nx-2][j][k][nb];}
        }
}

// Take the max 4 finite difference of neighbouring points
real max4diff(real**** W,int i,int j,int k,int nb) {
        real dif1=abs(W[i+1][j][k][nb]-W[i][j][k][nb]);
        real dif2=MAX(dif1,abs(W[i-1][j][k][nb]-W[i][j][k][nb]));
        dif1=MAX(dif2,abs(W[i][j+1][k][nb]-W[i][j][k][nb]));
        dif2=MAX(dif1,abs(W[i][j-1][k][nb]-W[i][j][k][nb]));
        return dif2;
}

// Sum up values in array using i- and j- range
real sum_irjr(real**** U, int nb, int k, 
		int i1, int i2, int j1, int j2) {
	real val=0.;
	for (int i=i1; i<=i2; i++) {
		for (int j=j1; j<=j2; j++) {
			val+=U[i][j][k][nb];
		}
	}
	if (i2-i1>1 or j2-j1>1) {
		cout<<" Note that sum_irjr contains summation of more than four cells!"<<endl;
		cout<<" i2-i1="<<i2-i1<<", j2-j1="<<j2-j1<<endl;
	}
	return val;
}

// Dot product of two vector
real dot_prod(vector<real> v1,int v1s,int v1e,vector<real> v2,int v2s,int v2e) {
	if (v2e-v2s!=v1e-v1s) {
		cout<<"Error as non-equal vectorlengths are used!"<<endl;
		throw exception();
	}
	real sum=0.;
	for (int i=0; i<=v1e-v1s; i++) {
		sum+=v1[v1s+i]*v2[v2s+i];
	}
	return sum;
}

// Produce sign of any value
real sign(real val) {
	if (val>0) {return 1.;}
	else if (val<0) {return -1.;}
	else {return 0.;}
}

// Compute wavespeeds (include magnetic fields)
void wavespeeds(vector<real> w,real gamma,int norm,real &a,real &ca,real &can,real &cf) {
	// norm is the normal direction 
	ca=sqrt(MAG(w[5],w[6],w[7])/w[0]);
        a=sqrt(gamma*w[4]/w[0]);
        can=sqrt(SQR(w[norm])/w[0]);
        real intm0=SQR(ca)+SQR(a);
        real intm1=SQR(ca)-SQR(a);
        cf=sqrt(0.5*(intm0+sqrt(intm1*intm1+4.*SQR(a)*(SQR(ca)-SQR(can)))));
}

// Convert primitive variables to conservative variables
void w2u(vector<real> w,vector<real> &u,real gamma) {
	u[0]=w[0];
	u[1]=w[0]*w[1];
	u[2]=w[0]*w[2];
	u[3]=w[0]*w[3];
	u[4]=w[4]/(gamma-1)+0.5*w[0]*MAG(w[1],w[2],w[3]);
	if (MAG_field) {
		u[4]+=0.5*MAG(w[5],w[6],w[7]);
		u[5]=w[5]; u[6]=w[6]; u[7]=w[7];
	}	
}

// Calculate Roe Averages (include magnetic field too)
void RoeAvg(vector<real> wL,vector<real> wR,real &RT,real &r,vector<real> &v,vector<real> &B,real &Xfac,real &Yfac) {
	for (int i=0;i<3;i++) {
		v[i]=(sqrt(wL[0])*wL[i+1]+sqrt(wR[0])*wR[i+1])/
			(sqrt(wL[0]+sqrt(wR[0])));
		B[i]=(sqrt(wL[0])*wL[i+5]+sqrt(wR[0])*wR[i+5])/
                        (sqrt(wL[0]+sqrt(wR[0])));
	}
	RT=sqrt(wR[0]/wL[0]);
	r=sqrt(wR[0]*wL[0]);
	Xfac=0.5*(pow(wL[6]-wR[6],2)+pow(wL[7]-wR[7],2))/
		(sqrt(wL[0])+sqrt(wR[0]));
	Yfac=0.5*(wL[0]+wR[0])/r;
}

// Calculate Roe eigenvalues (include magnetic field too)
void RoeEigen(real gamma,real d,vector<real> v,real h,vector<real> B,real Xfac,real Yfac,vector<real> &lambda) {
	real vsq=MAG(v[0],v[1],v[2]);
	real btsq=pow(B[1],2)+pow(B[2],2);
	real bt_starsq=(gamma-1.-(gamma-2)*Yfac)*btsq;
	real bt=sqrt(btsq);
	real bt_star=sqrt(bt_starsq);

	real vaxsq=pow(B[0],2)/d;
	real vax=sqrt(vaxsq);

	real hp=h-(vaxsq+btsq/d);
	real twid_asq=((gamma-1.)*(hp-0.5*vsq)-(gamma-2.)*Xfac);
	twid_asq=MAX(twid_asq,eps);
	real ct2=bt_starsq/d;
	real tsum=vaxsq+ct2+twid_asq;
	real tdif=vaxsq+ct2-twid_asq;
	real cf2_cs2=sqrt(tdif*tdif+4.*twid_asq*ct2);
	real cfsq=0.5*(tsum+cf2_cs2);
	real cfast=sqrt(cfsq);
	real cssq=twid_asq*vaxsq/cfsq;
	real cslow=sqrt(cssq);

	// Eigenvalues
	lambda[0]=v[0]-cfast;
	lambda[1]=v[0]-vax;
	lambda[2]=v[0]-cslow;
	lambda[3]=v[0];
	lambda[4]=v[0]+cslow;
	lambda[5]=v[0]+vax;
	lambda[6]=v[0]+cfast;
}
