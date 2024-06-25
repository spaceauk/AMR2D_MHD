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

#ifdef OPENMP
	dom.omp_time["celledge_tmp"]=omp_get_wtime();
#endif	
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
#ifdef OPENMP
        dom.omp_time["celledge_tmp"]-=omp_get_wtime();
	dom.omp_time["celledge"]-=dom.omp_time["celledge_tmp"];
#endif

	if (CT_mtd) {
#ifdef OPENMP
		dom.omp_time["CT2D_tmp"]=omp_get_wtime();
#endif		
		#pragma omp parallel for collapse(2) default(none) shared(dom,nb,step) 
		for (int i=0; i<dom.nx-1; i++) { 
                	for (int j=0; j<dom.ny-1; j++) {
                        	// Constrained transport method: Assumed magnetic field indices from 5 to 7
                                dom.wxL[i][j][5]= step==1 ? dom.Bi[i+1][j][0][nb] : dom.Bis[i+1][j][0][nb];
                                dom.wxR[i][j][5]= step==1 ? dom.Bi[i][j][0][nb] : dom.Bis[i][j][0][nb];
				dom.wyL[i][j][6]= step==1 ? dom.Bi[i][j+1][1][nb] : dom.Bis[i][j+1][1][nb];
                                dom.wyR[i][j][6]= step==1 ? dom.Bi[i][j][1][nb] : dom.Bis[i][j][1][nb];
			}
		}
#ifdef OPENMP
		dom.omp_time["CT2D_tmp"]-=omp_get_wtime();
		dom.omp_time["CT2D"]-=dom.omp_time["CT2D_tmp"];
#endif		
	}

#ifdef OPENMP
	dom.omp_time["RS_tmp"]=omp_get_wtime();
#endif	
	#pragma omp parallel for collapse(2) default(none) shared(dom,nb) 
	for (int i=dom.nxminb; i<=dom.nxmax; i++) {
		for (int j=dom.nyminb; j<=dom.nymax; j++) {
			vector<real> wL(8),wR(8),flux(8);
			int ip=i+1;
			int jp=j+1;
			// (a) In x-direction-----------------------------------------------------------------------
			for (int k=0; k<dom.nvar; k++) {
				wL[k]=dom.wxL[i][j][k];
				wR[k]=dom.wxR[ip][j][k];
			}
			// Compute flux at i+1/2
			riemannS(dom.fluxMth,wL,wR,dom.gamma,'x',flux);
			// Contribution to the residual of cell (i,j)	
			for (int k=0; k<dom.nvar; k++) {
				dom.ff[ip][j][k][nb]=flux[k];
			if (isnan(flux[k])) {
				real x=((dom.nxmin+dom.icoord[nb][0]-nghosts)+0.5)*dom.dx[dom.lp[nb][0]] - dom.dx[0];
				real y=((dom.nymin+dom.icoord[nb][1]-nghosts)+0.5)*dom.dy[dom.lp[nb][0]] - dom.dy[0];
				printf("At X: nb=%d, i=%d, j=%d, flux[%d]=%f, (x,y)=(%f,%f), pL=%f & pR=%f \n",
					nb,i,j,k,flux[k],x,y,wL[4],wR[4]);
				throw exception(); }
			}
			// Ensure divergence free magnetic field
			if (CT_mtd and wL[5]-wR[5]!=0) {
				printf("diffBx not equal to zero! BxL=%f, BxR=%f at i=%d, j=%d, nb=%d \n",
					wL[5],wR[5],i,j,nb);
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
				printf("At Y: nb=%d, i=%d, j=%d, flux[%d]=%f, (x,y)=(%f,%f) \n",
					nb,i,j,k,flux[k],x,y);
				throw exception(); }
			}
                        if (CT_mtd and wL[6]-wR[6]!=0) {
				printf("diffBy not equal to zero! ByL=%f, ByR=%f, at i=%d, j=%d, nb=%d \n",
					wL[6],wR[6],i,j,nb);
                                throw exception();
                        }
		}
	}
#ifdef OPENMP
	dom.omp_time["RS_tmp"]-=omp_get_wtime();
	dom.omp_time["RS"]-=dom.omp_time["RS_tmp"];
#endif
}
