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

real refineSL(real**** Q,int i,int j,int k,int nb,int signx,int signy) {
	real r,phi;
	real result=Q[i][j][k][nb];
	// Slope limiter for refining
	// Note that the sign is used to determine which side of the refined cell center lies on the original 
	// coarse cell center. E.g.,
	//          i          i+1
	//    |     .     |     .     |
	//    |  .  |  .  |  .  |  .  |
	//     i-1/2  i+1/2
	//       0     1      --> sign allocated here to tell which side
	r=(Q[i][j][k][nb]-Q[i-1][j][k][nb])/(Q[i+1][j][k][nb]-Q[i][j][k][nb] +eps);
	phi=0.25*slopelimiter(rlimiter,r,0)*(Q[i+1][j][k][nb]-Q[i][j][k][nb]);
        result = signx==0 ? result-phi : result+phi;

	r=(Q[i][j][k][nb]-Q[i][j-1][k][nb])/(Q[i][j+1][k][nb]-Q[i][j][k][nb] +eps);
        phi=0.25*slopelimiter(rlimiter,r,0)*(Q[i][j+1][k][nb]-Q[i][j][k][nb]);
        result = signy==0 ? result-phi : result+phi;
	return result;
}

void refine(meshblock &dom, real**** Q, int nvar, int nb, 
	    int son1, int son2, int son3, int son4) {
	int ii,jj;
	int signx,signy;
	// Update U
	#pragma omp parallel for collapse(3) default(none) \
	   shared(dom,Q,nvar,nb,son1,son2,son3,son4) private(ii,jj,signx,signy)
	for (int i=0;i<dom.nx;i++) {
	  for (int j=0;j<dom.ny;j++) {
	    for (int k=0;k<nvar;k++) {
	      signx=(i+nghosts)%2; signy=(j+nghosts)%2;
	      ii=(i+nghosts)/2; jj=(j+nghosts)/2;
	      Q[i][j][k][son1]=refineSL(Q,ii,jj,k,nb,signx,signy);
	      ii=(i-nghosts+2)/2+dom.nx2;	      
	      Q[i][j][k][son2]=refineSL(Q,ii,jj,k,nb,signx,signy);
	      ii=(i+nghosts)/2; jj=(j-nghosts+2)/2+dom.ny2;
	      Q[i][j][k][son3]=refineSL(Q,ii,jj,k,nb,signx,signy);
	      ii=(i-nghosts+2)/2+dom.nx2;
	      Q[i][j][k][son4]=refineSL(Q,ii,jj,k,nb,signx,signy);
	    }		
	  }
	}
	// Obtain block original icoords but refined to current resolution
	real nx0=dom.icoord[nb][0]*2;
	real ny0=dom.icoord[nb][1]*2;
	//  Note that icoord is written such that it is the x/y position
	//  of the block in terms of the number of cells which is at the 
	//  assumed resolution of the leaf block. It is written like this 
	//  cause of how the x & y are computed over here.
	
	// Update icoords
	cout<<" ["<<son1<<","<<son2<<","<<son3<<","<<son4<<"] ";
	dom.U2W("refine",son1);
	dom.icoord[son1][0]=nx0;
	dom.icoord[son1][1]=ny0;

	dom.U2W("refine",son2);
        dom.icoord[son2][0]=nx0+dom.nx-2*nghosts;
        dom.icoord[son2][1]=ny0;

	dom.U2W("refine",son3);
        dom.icoord[son3][0]=nx0;
        dom.icoord[son3][1]=ny0+dom.ny-2*nghosts;

	dom.U2W("refine",son4);
        dom.icoord[son4][0]=nx0+dom.nx-2*nghosts;
        dom.icoord[son4][1]=ny0+dom.ny-2*nghosts;
}
