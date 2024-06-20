#include <vector>
#include <cmath>
#include <exception>
#include<iomanip>
#include "defs.hpp"

void celledges(meshblock &dom,int nb,int step);
void riemannS(string fluxMth,vector<real> wL,vector<real> wR,real gamma,char direc,vector<real> &flux);
void savearray(meshblock &dom, real*** array, string arrname);
void WENO2D_3rd(meshblock &dom,string direc,int nb);
void WENO2D_5th(meshblock &dom,string direc,int nb);

void MUSCL2D(meshblock &dom, int nb, int step) {
	// Convert U to W
	if (step==1) {
		dom.U2W("MUSCL",nb);
	} else if (step==2) {
		dom.Us2W("MUSCL",nb);
	}

	if (dom.limiter=="WENO3") {
		WENO2D_3rd(dom,"x",nb);
		WENO2D_3rd(dom,"y",nb);
	} else if (dom.limiter=="WENO5") { 
		// Note that this directional splitting of WENO is 2nd order accurate, thus more robust method is required.
		WENO2D_5th(dom,"x",nb);
		WENO2D_5th(dom,"y",nb);
	} else {
		// Apply slope limiter
		celledges(dom,nb,step);
	}

	if (CT_mtd) {
		for (int i=0; i<dom.nx-1; i++) { 
                	for (int j=0; j<dom.ny-1; j++) {
                        	// Constrained transport method: Assumed magnetic field indices from 5 to 7
                                dom.wxL[i][j][5]= step==1 ? dom.Bi[i+1][j][0][nb] : dom.Bis[i+1][j][0][nb];
                                dom.wxR[i][j][5]= step==1 ? dom.Bi[i][j][0][nb] : dom.Bis[i][j][0][nb];
				dom.wyL[i][j][6]= step==1 ? dom.Bi[i][j+1][1][nb] : dom.Bis[i][j+1][1][nb];
                                dom.wyR[i][j][6]= step==1 ? dom.Bi[i][j][1][nb] : dom.Bis[i][j][1][nb];
			}
		}
	}
                                                                    

	for (int i=0;i<dom.nx;i++) {
		for (int j=0;j<dom.ny;j++) {
			for (int k=0;k<dom.nvar;k++) {
				dom.ff[i][j][k][nb]=0.; dom.gg[i][j][k][nb]=0.;
			}
		}
	}

	vector<real> wL(8),wR(8),flux(8);
	for (int i=dom.nxminb; i<=dom.nxmax; i++) {
		for (int j=dom.nyminb; j<=dom.nymax; j++) {
			int ip=i+1;
			int jp=j+1;
			// (a) In x-direction-----------------------------------------------------------------------
			for (int k=0; k<dom.nvar; k++) {
				wL[k]=dom.wxL[i][j][k];
				wR[k]=dom.wxR[ip][j][k];
			}
			if (wL[4]<0) {
				cout<<"Negative pressure at i="<<i<<", j="<<j<<", nb="<<nb<<endl;
				cout<<"pL="<<wL[4]<<", pR="<<wR[4]<<endl;
				throw exception();
			}
			// Compute flux at i+1/2
			riemannS(dom.fluxMth,wL,wR,dom.gamma,'x',flux);
			// Contribution to the residual of cell (i,j)	
			for (int k=0; k<dom.nvar; k++) {
				dom.ff[ip][j][k][nb]=flux[k];
			if (isnan(flux[k])) {
				real x=((dom.nxmin+dom.icoord[nb][0]-nghosts)+0.5)*dom.dx[dom.lp[nb][0]] - dom.dx[0];
				real y=((dom.nymin+dom.icoord[nb][1]-nghosts)+0.5)*dom.dy[dom.lp[nb][0]] - dom.dy[0];
				cout<<"At X: nb="<<nb<<", i="<<i<<", j="<<j<<", flux["<<k<<"]="<<flux[k]
					<<", (x,y)=("<<x<<","<<y<<")"<<endl;
				cout<<"inner_bounds=("<<dom.innerbounds[nb][0]<<","<<dom.innerbounds[nb][1]<<","<<
					dom.innerbounds[nb][2]<<")"<<endl;
				cout<<"W=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<dom.W[i][j][kk][nb]<<" ";} cout<<")"<<endl;
				cout<<"dwdx=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<dom.dwdx[i][j][kk]<<" ";} cout<<")"<<endl;
				cout<<"wL=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<wL[kk]<<" ";} cout<<")"<<endl;
				cout<<"wR=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<wR[kk]<<" ";} cout<<")"<<endl;
				cout<<"flux=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<flux[kk]<<" ";} cout<<")"<<endl;	
				throw exception(); }
			}
			// Ensure divergence free magnetic field
			if (CT_mtd and wL[5]-wR[5]!=0) {
				cout<<std::setprecision(Nprec)<<"diffBx not equal to zero! BxL="<<wL[5]<<", BxR="<<wR[5]
					<<", at i="<<i<<", j="<<j<<", nb="<<nb<<endl;
				throw exception();
			}

			// (b) In y-direction----------------------------------------------------------------------
			for (int k=0; k<dom.nvar; k++) {
				wL[k]=dom.wyL[i][j][k];
				wR[k]=dom.wyR[i][jp][k];
			}
			// Compute flux at j+1/2
			riemannS(dom.fluxMth,wL,wR,dom.gamma,'y',flux);
			// Contribution to the residual of cell (i,j)
			for (int k=0; k<dom.nvar; k++) {
				dom.gg[i][jp][k][nb]=flux[k];
			if (isnan(flux[k])) {
				real x=((dom.nxmin+dom.icoord[nb][0]-nghosts)+0.5)*dom.dx[dom.lp[nb][0]] - dom.dx[0];
				real y=((dom.nymin+dom.icoord[nb][1]-nghosts)+0.5)*dom.dy[dom.lp[nb][0]] - dom.dy[0];
				cout<<"At Y: nb="<<nb<<", i="<<i<<", j="<<j<<", flux["<<k<<"]="<<flux[k]
					<<", (x,y)=("<<x<<","<<y<<")"<<endl;
                                cout<<"inner_bounds=("<<dom.innerbounds[nb][0]<<","<<dom.innerbounds[nb][1]<<","<<
                                        dom.innerbounds[nb][2]<<")"<<endl;
                                cout<<"W=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<dom.W[i][j][kk][nb]<<" ";} cout<<")"<<endl;
                                cout<<"dwdx=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<dom.dwdy[i][j][kk]<<" ";} cout<<")"<<endl;
                                cout<<"wL=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<wL[kk]<<" ";} cout<<")"<<endl;
                                cout<<"wR=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<wR[kk]<<" ";} cout<<")"<<endl;
                                cout<<"flux=("; for (int kk=0;kk<dom.nvar;kk++) {cout<<flux[kk]<<" ";} cout<<")"<<endl;
				throw exception(); }
			}
                        if (CT_mtd and wL[6]-wR[6]!=0) {
                                cout<<std::setprecision(Nprec)<<"diffBy not equal to zero! ByL="<<wL[6]<<", ByR="<<wR[6]
                                        <<", at i="<<i<<", j="<<j<<", nb="<<nb<<endl;
                                throw exception();
                        }
		}
	}
}
