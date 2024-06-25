// Refine to next level for face-centered variables
#include "defs.hpp"
real slopelimiter(string limiter,real r,real eta);

real refineSL_fc(real**** Bi,int nx,int ny,int i,int j,int k,int nb,int signx,int signy) {
	real r,phi;
	real result=Bi[i][j][k][nb];
	// Slope limiter for refining
	// Note that the sign is used to determine which side of the refined cell center lies on the original 
	// coarse cell center. E.g.,
	//          i          i+1
	//    |     .     |     .     |
	//    |  .  |  .  |  .  |  .  |
	//     i-1/2  i+1/2
	//       0     1      --> sign allocated here to tell which side
	real sx=0., sy=0.;
	r=(Bi[i][j][k][nb]-Bi[i-1][j][k][nb])/(Bi[i+1][j][k][nb]-Bi[i][j][k][nb] +eps);
	if (k==0) {
		sx= signx==0 ? 0. : 0.5;
	} else if (k==1) {
		sx= signx==0 ? -0.25 : 0.25;
	}
	phi=slopelimiter(rlimiter,r,0)*(Bi[i+1][j][k][nb]-Bi[i][j][k][nb]);
	result = result+sx*phi;

	r=(Bi[i][j][k][nb]-Bi[i][j-1][k][nb])/(Bi[i][j+1][k][nb]-Bi[i][j][k][nb] +eps);
	if (k==0) {
		sy= signy==0 ? -0.25 : 0.25;
	} else if (k==1) {
		sy= signy==0 ? 0. : 0.5;
	}
        phi=slopelimiter(rlimiter,r,0)*(Bi[i][j+1][k][nb]-Bi[i][j][k][nb]);
        result = result+sy*phi;
	return result;
}

void refinefc(meshblock &dom, real**** Bi, int nvar, int nb, 
	    int son1, int son2, int son3, int son4) {
	int ii,jj;
	int signx,signy;
	// Update U
	#pragma omp parallel for collapse(3) default(none) \
           shared(dom,Bi,nvar,nb,son1,son2,son3,son4) private(ii,jj,signx,signy)
	for (int i=0;i<dom.nx;i++) {
	  for (int j=0;j<dom.ny;j++) {
	    for (int k=0;k<nvar;k++) {
	      signx=(i+nghosts)%2; signy=(j+nghosts)%2;
	      ii=(i+nghosts)/2; jj=(j+nghosts)/2;
	      Bi[i][j][k][son1]=refineSL_fc(Bi,dom.nx,dom.ny,ii,jj,k,nb,signx,signy);
	      ii=(i-nghosts+2)/2+dom.nx2;	      
	      Bi[i][j][k][son2]=refineSL_fc(Bi,dom.nx,dom.ny,ii,jj,k,nb,signx,signy);
	      ii=(i+nghosts)/2; jj=(j-nghosts+2)/2+dom.ny2;
	      Bi[i][j][k][son3]=refineSL_fc(Bi,dom.nx,dom.ny,ii,jj,k,nb,signx,signy);
	      ii=(i-nghosts+2)/2+dom.nx2;
	      Bi[i][j][k][son4]=refineSL_fc(Bi,dom.nx,dom.ny,ii,jj,k,nb,signx,signy);
	    }		
	  }
	}
}
