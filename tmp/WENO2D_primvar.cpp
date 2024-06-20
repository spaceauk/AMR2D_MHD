// This WENO scheme works well with the RS. Obtained from https://doi.org/10.1016/j.compfluid.2015.04.026
#include "defs.hpp"

void savearray(meshblock &dom,real*** array, string arrname);

// 3rd order WENO reconstruction
void WENO2D_3rd(meshblock &dom,string direc,int nb){
	int ii,jj;
#ifdef DOUBLE_PREC
	real epsilon=1.E-15;
#else
	real epsilon=1.E-6;
#endif
	for (int i=dom.nxminb; i<=dom.nxmaxb; i++) {
                for (int j=dom.nyminb; j<=dom.nymaxb; j++) {
                        for (int k=0; k<dom.nvar; k++) {
				// Shift Forward
                                if (direc=="x"){
                                        ii=i+1; jj=j;
                                } else if (direc=="y"){
                                        ii=i; jj=j+1;
                                }
				real wL=dom.W[i][j][k][nb];
				real wR=dom.W[i][j][k][nb];

                                real wLp=dom.W[ii][jj][k][nb];
                                real wRp=dom.W[ii][jj][k][nb];
	
				// Shift Backwards
                                if (direc=="x"){
                                        ii=i-1; jj=j;
                                } else if (direc=="y"){
                                        ii=i; jj=j-1;
                                }
                                real wLm=dom.W[ii][jj][k][nb];
                                real wRm=dom.W[ii][jj][k][nb];

				// Obtain right flux---------------------
				// Polynomials (from eqn (3.14): High-order WENO FVM on Cartesian Grids with AMR)
				real p0n = 0.5*wRm+0.5*wR;
				real p1n = 1.5*wR-0.5*wRp;
				// Smoothness indicators (Beta factors)
				real B0n = pow((wR-wRm),2);
				real B1n = pow((wRp-wR),2);
				// Constants
				real d0n=2./3.; 
				real d1n=1./3.; 
			        // Alpha weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real alpha0n = d0n/pow((epsilon+B0n),2);
				real alpha1n = d1n/pow((epsilon+B1n),2);
				real alphasumn = alpha0n + alpha1n;
				// ENO stencils weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real w0n = alpha0n/alphasumn;
				real w1n = alpha1n/alphasumn;
				// Numerical Flux at cell boundary
				if (direc=="x") {
					dom.wxR[i][j][k] = w0n*p0n + w1n*p1n;
				} else if (direc=="y") {
					dom.wyR[i][j][k] = w0n*p0n + w1n*p1n;
				}

				// Obtain left flux----------------------
				// Polynomials (from eqn (3.14): High-order WENO FVM on Cartesian Grids with AMR)
				real p0p = 1.5*wL-0.5*wLm;
				real p1p = 0.5*wLp+0.5*wL;
				// Smoothness indicators (Beta factors)
				real B0p = pow((wL-wLm),2);
				real B1p = pow((wLp-wL),2);
				// Constants
				real d0p=1./3.; 
				real d1p=2./3.; 
				// Alpha weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real alpha0p = d0p/pow((epsilon+B0p),2);
				real alpha1p = d1p/pow((epsilon+B1p),2);
				real alphasump = alpha0p + alpha1p;
				// ENO stencils weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real w0p = alpha0p/alphasump;
				real w1p = alpha1p/alphasump;
				// Numerical Flux at cell boundary
				if (direc=="x") {
					dom.wxL[i][j][k] = w0p*p0p + w1p*p1p;
				} else if (direc=="y") {
					dom.wyL[i][j][k] = w0p*p0p + w1p*p1p;
				}
                        }
                }
        }
}


// 5th order WENO reconstruction
void WENO2D_5th(meshblock &dom,string direc,int nb){
	int ii,jj,iii,jjj;
#ifdef DOUBLE_PREC
        real epsilon=1.E-15;
#else
        real epsilon=1.E-6;
#endif	
	for (int i=dom.nxminb; i<=dom.nxmaxb; i++) {
                for (int j=dom.nyminb; j<=dom.nymaxb; j++) {
                        for (int k=0; k<dom.nvar; k++) {
				// Shift Forward
                                if (direc=="x"){
                                        ii=i+1; jj=j;
			                iii=i+2; jjj=j;
                                } else if (direc=="y"){
                                        ii=i; jj=j+1;
                                        iii=i; jjj=j+2;
                                }
				real wL=dom.W[i][j][k][nb];
				real wR=dom.W[i][j][k][nb];

                                real wLp=dom.W[ii][jj][k][nb];
				real wLpp=dom.W[iii][jjj][k][nb];
                                real wRp=dom.W[ii][jj][k][nb];
				real wRpp=dom.W[iii][jjj][k][nb];
	
				// Shift Backwards
                                if (direc=="x"){
                                        ii=i-1; jj=j;
			                iii=i-2; jjj=j;
                                } else if (direc=="y"){
                                        ii=i; jj=j-1;
                                        iii=i; jjj=j-2;
                                }
                                real wLm=dom.W[ii][jj][k][nb];
				real wLmm=dom.W[iii][jjj][k][nb];
                                real wRm=dom.W[ii][jj][k][nb];
				real wRmm=dom.W[iii][jjj][k][nb];

				// Obtain right flux---------------------
				//Polynomials (from eqn (3.14): High-order WENO FVM on Cartesian Grids with AMR)
				real p0n = (-1*wRmm+5*wRm+2*wR)/6;
				real p1n = (2*wRm+5*wR-wRp)/6;
				real p2n = (11*wR-7*wRp+2*wRpp)/6;
				// Smoothness indicators (Beta factors)
				real B0n = (13./12.)*pow((wRmm-2*wRm+wR),2) + (1./4.)*pow((wRmm-4*wRm+3*wR),2);
				real B1n = (13./12.)*pow((wRmm-2*wR+wRp),2) + (1./4.)*pow((wRm-wRp),2);
				real B2n = (13./12.)*pow((wR-2*wRp+wRpp),2) + (1./4.)*pow((3*wR-4*wRp+wRpp),2);
				// Constants
				real d0n=1./10.; 
				real d1n=6./10.; 
				real d2n=3./10.; 
			        // Alpha weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real alpha0n = d0n/pow((epsilon+B0n),2);
				real alpha1n = d1n/pow((epsilon+B1n),2);
				real alpha2n = d2n/pow((epsilon+B2n),2);
				real alphasumn = alpha0n + alpha1n + alpha2n;
				// ENO stencils weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real w0n = alpha0n/alphasumn;
				real w1n = alpha1n/alphasumn;
				real w2n = alpha2n/alphasumn;
				// Numerical Flux at cell boundary
				if (direc=="x") {
					dom.wxR[i][j][k] = w0n*p0n + w1n*p1n +w2n*p2n;
				} else if (direc=="y") {
					dom.wyR[i][j][k] = w0n*p0n + w1n*p1n +w2n*p2n;
				}

				// Obtain left flux----------------------
				// Polynomials (from eqn (3.14): High-order WENO FVM on Cartesian Grids with AMR)
				real p0p = (2*wLmm-7*wLm+11*wL)/6;
				real p1p = (-1*wLm+5*wL+2*wLp)/6;
				real p2p = (2*wL+5*wLp-wLpp)/6;
				// Smoothness indicators (Beta factors)
				real B0p = (13./12.)*pow((wLmm-2*wLm+wL),2) + (1./4.)*pow((wLmm-4*wLm+3*wL),2);
				real B1p = (13./12.)*pow((wLmm-2*wL+wLp),2) + (1./4.)*pow((wLm-wLp),2);
				real B2p = (13./12.)*pow((wL-2*wLp+wLpp),2) + (1./4.)*pow((3*wL-4*wLp+wLpp),2);
				// Constants
				real d0p=3./10.; 
				real d1p=6./10.; 
				real d2p=1./10.; 
				// Alpha weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real alpha0p = d0p/pow((epsilon+B0p),2);
				real alpha1p = d1p/pow((epsilon+B1p),2);
				real alpha2p = d2p/pow((epsilon+B2p),2);
				real alphasump = alpha0p + alpha1p + alpha2p;
				// ENO stencils weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real w0p = alpha0p/alphasump;
				real w1p = alpha1p/alphasump;
				real w2p = alpha2p/alphasump;
				// Numerical Flux at cell boundary
				if (direc=="x") {
					dom.wxL[i][j][k] = w0p*p0p + w1p*p1p +w2p*p2p;
				} else if (direc=="y") {
					dom.wyL[i][j][k] = w0p*p0p + w1p*p1p +w2p*p2p;
				}
                        }
                }
        }
}



