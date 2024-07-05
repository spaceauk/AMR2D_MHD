#include "defs.hpp"
real slopelimiter(string limiter,real r,real eta);

real refineSL_Bi(real**** Bi,int nx,int ny,int i,int j,int k,int nb,int signx,int signy) {
	real r,DyBxm,DyBxp;
	real DxBym,DxByp;
	real result;
	string CT_rlimit=rlimiter;

	real Bx_jm=Bi[i][j-1][0][nb], Bx=Bi[i][j][0][nb], Bx_jp=Bi[i][j+1][0][nb];
	real Bx_ip=Bi[i+1][j][0][nb], Bx_ipjm=Bi[i+1][j-1][0][nb], Bx_ipjp=Bi[i+1][j+1][0][nb];
	real By_im=Bi[i-1][j][1][nb], By=Bi[i][j][1][nb], By_ip=Bi[i+1][j][1][nb];
        real By_imjp=Bi[i-1][j+1][1][nb], By_jp=Bi[i][j+1][1][nb], By_ipjp=Bi[i+1][j+1][1][nb]; 	

	// Slope limiter for refining
	r=(Bx-Bx_jm)/(Bx_jp-Bx +eps);
	DyBxm=slopelimiter(CT_rlimit,r,0)*(Bx_jp-Bx);
	r=(Bx_ip-Bx_ipjm)/(Bx_ipjp-Bx_ip +eps);
        DyBxp=slopelimiter(CT_rlimit,r,0)*(Bx_ipjp-Bx_ip);

	r=(By-By_im)/(By_ip-By +eps);
        DxBym=slopelimiter(CT_rlimit,r,0)*(By_ip-By);
	r=(By_jp-By_imjp)/(By_ipjp-By_jp +eps);
        DxByp=slopelimiter(CT_rlimit,r,0)*(By_ipjp-By_jp);

	real ax=(Bx_ip-Bx);
	real by=-1.*ax;
	real ay=0.5*(DyBxp+DyBxm);
	real bx=0.5*(DxByp+DxBym);
	real axy=(DyBxp-DyBxm);
	real byy=axy/(-2.);
	real bxy=(DxByp-DxBym);
	real axx=bxy/(-2.);
	real a0=0.5*(Bx_ip+Bx)-0.25*axx;
	real b0=0.5*(By_jp+By)-0.25*byy;
	real ayy=0., bxx=0.;

	if (k==0) {
		real sx= signx==0 ? 0 : 0.5;
        	real sy= signy==0 ? -0.25 : 0.25;
		result=a0+ax*sx+ay*sy+axx*sx*sx+axy*sx*sy+ayy*sy*sy;
	} else if (k==1) {
		real sx= signx==0 ? -0.25 : 0.25;
        	real sy= signy==0 ? 0 : 0.5;
		result=b0+bx*sx+by*sy+bxx*sx*sx+bxy*sx*sy+byy*sy*sy;
	}
	if (isnan(result)) {cout<<"Error at refine Bi... i="<<i<<", j="<<j<<", k="<<k<<endl; throw exception();} 
	
	return result;
}

void refineBi_CT(meshblock &dom, real**** Bi, int nvar, int nb, 
	    int son1, int son2, int son3, int son4) {
	int ii,jj;
	int signx,signy;
	// Update Bi after refinement
	#pragma omp parallel for collapse(3) default(none) \
           shared(dom,Bi,nvar,nb,son1,son2,son3,son4) private(ii,jj,signx,signy)
	for (int i=0;i<dom.nx;i++) {
	  for (int j=0;j<dom.ny;j++) {
	    for (int k=0;k<nvar;k++) {
	      signx=(i+nghosts)%2; signy=(j+nghosts)%2;
	      ii=(i+nghosts)/2; jj=(j+nghosts)/2;
	      Bi[i][j][k][son1]=refineSL_Bi(Bi,dom.nx,dom.ny,ii,jj,k,nb,signx,signy);
	      ii=(i-nghosts+2)/2+dom.nx2;	      
	      Bi[i][j][k][son2]=refineSL_Bi(Bi,dom.nx,dom.ny,ii,jj,k,nb,signx,signy);
	      ii=(i+nghosts)/2; jj=(j-nghosts+2)/2+dom.ny2;
	      Bi[i][j][k][son3]=refineSL_Bi(Bi,dom.nx,dom.ny,ii,jj,k,nb,signx,signy);
	      ii=(i-nghosts+2)/2+dom.nx2;
	      Bi[i][j][k][son4]=refineSL_Bi(Bi,dom.nx,dom.ny,ii,jj,k,nb,signx,signy);
	    }		
	  }
	}
}
