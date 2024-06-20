#include "defs.hpp"
void savearray(meshblock &dom,real*** array,string arrname);
real slopelimiter(string limiter,real r,real eta);

void celledges(meshblock &dom,int nb,int step) {
	// Obtain cell differences
	for (int i=0; i<dom.nx; i++) {
		for (int j=0; j<dom.ny; j++) {
			for (int k=0; k<dom.nvar; k++) {
				int im= i>0 ? i-1 : i;
				int jm= j>0 ? j-1 : j;
				dom.dwdx[i][j][k]=dom.W[i][j][k][nb]-dom.W[im][j][k][nb];
				dom.dwdy[i][j][k]=dom.W[i][j][k][nb]-dom.W[i][jm][k][nb];
				if (isnan(dom.dwdx[i][j][k]) or isnan(dom.dwdy[i][j][k])) {
					cout<<"Nan values encountered - dx="<<dom.dwdx[i][j][k]<<", dy="<<dom.dwdy[i][j][k]<<
						" at i="<<i<<", j="<<j<<", k="<<k<<", nb="<<nb<<endl;
					cout<<"Wim="<<dom.W[im][j][k][nb]<<", W="<<dom.W[i][j][k][nb]<<", Wjm="
						<<dom.W[i][jm][k][nb]<<", Us="<<dom.Us[i][j][k][nb]<<", U="<<dom.U[i][j][k][nb]<<endl;
					throw exception();					
				}				
			}
		}
	}

	// Obtain cell left and right edges
	real r, phi;
	real eta, smoothr=1; // Smoothness indicator for Compact 3rd order scheme (M.Cada)
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
					cout<<"celledge -> nan! @ i+1="<<i+1<<", j="<<j<<", k="<<k
					<<", nb="<<nb<<", dwdx="<<dom.dwdx[i+1][j][k]<<", W="<<dom.W[i+1][j][k][nb]<<endl;
					throw exception();}
			}
		}
	}

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
