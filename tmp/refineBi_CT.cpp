// Refine to next level
// Note how the dad block is shared with son blocks: E.g., nx-1=35 (w/ ghosts) , nxmax=33, nxmin=2
//  i     (i+2)/2     i/2+nx2        -> This ensures that the internal cells of dad block are
//  0        1          17      (g)     all matched to the internal cells of son block.
//  1	     1          17      (b)
//  2        2          18      (i)  
//  3        2          18
//           :
//  32       17         33
//  33       17         33      (i)
//  34       18         34      (b)
//  35       18         34      (g)

#include "defs.hpp"
real slopelimiter(string limiter,real r,real eta);

real refineSL_Bi(real**** Bi,int nx,int ny,int i,int j,int k,int nb,int signx,int signy) {
	real r,DyBxm,DyBxp;
	real DxBym,DxByp;
	real result;
	string CT_rlimit=rlimiter;

	// Slope limiter for refining
	r=(Bi[i][j+1][0][nb]-Bi[i][j][0][nb])/(Bi[i][j][0][nb]-Bi[i][j-1][0][nb] +eps);
	DyBxm=slopelimiter(CT_rlimit,r,0)*(Bi[i][j][0][nb]-Bi[i][j-1][0][nb]);
	r=(Bi[i+1][j+1][0][nb]-Bi[i+1][j][0][nb])/(Bi[i+1][j][0][nb]-Bi[i+1][j-1][0][nb] +eps);
        DyBxp=slopelimiter(CT_rlimit,r,0)*(Bi[i+1][j][0][nb]-Bi[i+1][j-1][0][nb]);

	r=(Bi[i+1][j][1][nb]-Bi[i][j][1][nb])/(Bi[i][j][1][nb]-Bi[i-1][j][1][nb] +eps);
        DxBym=slopelimiter(CT_rlimit,r,0)*(Bi[i][j][1][nb]-Bi[i-1][j][1][nb]);
	r=(Bi[i+1][j+1][1][nb]-Bi[i][j+1][1][nb])/(Bi[i][j+1][1][nb]-Bi[i-1][j+1][1][nb] +eps);
        DxByp=slopelimiter(CT_rlimit,r,0)*(Bi[i][j+1][1][nb]-Bi[i-1][j+1][1][nb]);

	real ax=0.5*(Bi[i+1][j][0][nb]-Bi[i][j][0][nb]);
	real by=-1.*ax;
	real ay=0.25*(DyBxp+DyBxm);
	real bx=0.25*(DxByp+DxBym);
	real axy=0.25*(DyBxp-DyBxm);
	real byy=axy/(-2.);
	real bxy=0.25*(DxByp-DxBym);
	real axx=bxy/(-2.);
	real a0=0.5*(Bi[i+1][j][0][nb]+Bi[i][j][0][nb])-0.25*axx;
	real b0=0.5*(Bi[i][j+1][1][nb]+Bi[i][j][1][nb])-0.25*byy;
	real ayy=0., bxx=0.;

	int sx= signx==0 ? -1 : 1;
	int sy= signy==0 ? -1 : 1;
	if (k==0) {
		result=a0+ax*sx+ay*sy+axx+axy*sx*sy+ayy;
	} else if (k==1) {
		result=b0+bx*sx+by*sy+bxx+bxy*sx*sy+byy;
	}
	//if (abs(result-Bi[i][j][k][nb])>0.1) {
	if (i==0 or j==0 or i==nx-1 or j==ny-1) {
		cout<<endl<<"Bi="<<Bi[i][j][k][nb]<<", result="<<result<<"; "<<endl;
	        cout<<"Bi_ip="<<Bi[i+1][j][k][nb]<<", Bi_jp="<<Bi[i][j+1][k][nb]<<endl;
		cout<<" At i="<<i<<", j="<<j<<", k="<<k<<endl;
		throw exception();}
	if (isnan(result)) {cout<<"Error at refine Bi... i="<<i<<", j="<<j<<", k="<<k<<endl; throw exception();}
	return result;
}

void refineBi_CT(meshblock &dom, real**** Bi, int nvar, int nb, 
	    int son1, int son2, int son3, int son4) {
	int ii,jj;
	int signx,signy;
	// Update Bi after refinement
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
