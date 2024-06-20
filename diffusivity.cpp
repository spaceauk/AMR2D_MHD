// Applies diffusivity to all variables thus, not just momentum diffusivity (or viscosity) for this case.
// up(n)=up(n)+\eta\nabla^2(u(n))

#include "defs.hpp"

void diffusivity(real**** U,real eta,int nb,int nxmin,int nxmax,
		int nymin,int nymax,int nvar) {

	for (int i=nxmin;i<=nxmax;i++) {
		for (int j=nymin;j<=nymax;j++) {
			for (int k=0;k<nvar;k++) {
				U[i][j][k][nb]=U[i][j][k][nb]+eta*(U[i+1][j][k][nb]+U[i-1][j][k][nb]+
					U[i][j+1][k][nb]+U[i][j-1][k][nb]-4.*U[i][j][k][nb]);
			}
		}
	}
}
