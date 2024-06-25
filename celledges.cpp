#include "defs.hpp"
void savearray(meshblock &dom,real*** array,string arrname);
real slopelimiter(string limiter,real r,real eta);

void celledges(meshblock &dom,int nb,int step) {
	// Obtain cell differences
	int im,jm,ii,jj;
	#pragma omp parallel for collapse(3) default(none) shared(dom,nb) private(im,jm,ii,jj)
	for (int i=0; i<dom.nx; i++) {
		for (int j=0; j<dom.ny; j++) {
			for (int k=0; k<dom.nvar; k++) {
				im= i>0 ? i-1 : i;
				jm= j>0 ? j-1 : j;
				ii=im+1; jj=jm+1;
				dom.dwdx[i][j][k]=dom.W[ii][j][k][nb]-dom.W[im][j][k][nb];
				dom.dwdy[i][j][k]=dom.W[i][jj][k][nb]-dom.W[i][jm][k][nb];
				if (isnan(dom.dwdx[i][j][k]) or isnan(dom.dwdy[i][j][k])) {
					printf("NaN values encountered where dx=%f, dy=%f at i=%d, j=%d, k=%d, nb=%d \n",
					      dom.dwdx[i][j][k],dom.dwdy[i][j][k],i,j,k,nb);
					throw exception();					
				}				
			}
		}
	}

	// Obtain cell left and right edges
	real r, phi;
	real eta, smoothr=1; // Smoothness indicator for Compact 3rd order scheme (M.Cada)
	#pragma omp parallel for collapse(3) default(none) shared(dom,nb) private(r,phi,eta,smoothr)
	for (int i=0; i<dom.nx-1; i++) { // For x
		for (int j=0; j<dom.ny; j++) {
			for (int k=0; k<dom.nvar; k++) {
				if (CT_mtd and k==5) {continue;}
				if (dom.limiter=="Cada3rd") {
					eta=(pow(dom.dwdx[i+1][j][k],2)+pow(dom.dwdx[i][j][k],2))
						/pow(smoothr*dom.dx[dom.lp[nb][0]],2);
				} else {eta=0.;}
				r=dom.dwdx[i+1][j][k]/(dom.dwdx[i][j][k]+eps);
				phi=slopelimiter(dom.limiter,r,eta);
				dom.wxL[i][j][k]=dom.W[i][j][k][nb]+0.5*phi*dom.dwdx[i][j][k];
				r=dom.dwdx[i][j][k]/(dom.dwdx[i+1][j][k]+eps);
				phi=slopelimiter(dom.limiter,r,eta);
				dom.wxR[i][j][k]=dom.W[i][j][k][nb]-0.5*phi*dom.dwdx[i+1][j][k];
				if (isnan(dom.wxL[i][j][k])) {
					printf("celledge -> nan @ i+1=%d, j=%d, k=%d, nb=%d, dwdx=%f, W=%f \n",
						i+1,j,k,nb,dom.dwdx[i+1][j][k],dom.W[i+1][j][k][nb]);
					throw exception();}
			}
		}
	}

	#pragma omp parallel for collapse(3) default(none) shared(dom,nb) private(r,phi,eta,smoothr)
	for (int i=0; i<dom.nx; i++) { // For y
		for (int j=0; j<dom.ny-1; j++) {
			for (int k=0; k<dom.nvar; k++) {
				if (CT_mtd and k==6) {continue;}
				if (dom.limiter=="Cada3rd") {
                                        eta=(pow(dom.dwdy[i][j+1][k],2)+pow(dom.dwdy[i][j][k],2))
                                                /pow(smoothr*dom.dy[dom.lp[nb][0]],2);
                                } else {eta=0.;}
				r=dom.dwdy[i][j+1][k]/(dom.dwdy[i][j][k]+eps);
                                phi=slopelimiter(dom.limiter,r,eta);
				dom.wyL[i][j][k]=dom.W[i][j][k][nb]+0.5*phi*dom.dwdy[i][j][k];
				r=dom.dwdy[i][j][k]/(dom.dwdy[i][j+1][k]+eps);
				phi=slopelimiter(dom.limiter,r,eta);
				dom.wyR[i][j][k]=dom.W[i][j][k][nb]-0.5*phi*dom.dwdy[i][j+1][k];
			}
		}
	}

}	
