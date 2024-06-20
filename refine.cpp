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

void refine(meshblock &dom, real**** Q, int nvar, int nb, 
	    int son1, int son2, int son3, int son4) {
	// Update U
	for (int i=0;i<dom.nx;i++) {
	  for (int j=0;j<dom.ny;j++) {
	    for (int k=0;k<nvar;k++) {
	      Q[i][j][k][son1]=Q[(i+nghosts)/2][(j+nghosts)/2][k][nb];
	      Q[i][j][k][son2]=Q[(i-nghosts+2)/2+dom.nx2][(j+nghosts)/2][k][nb];
	      Q[i][j][k][son3]=Q[(i+nghosts)/2][(j-nghosts+2)/2+dom.ny2][k][nb];
	      Q[i][j][k][son4]=Q[(i-nghosts+2)/2+dom.nx2][(j-nghosts+2)/2+dom.ny2][k][nb];
	    }		
	  }
	}
	// Obtain block original icoords but refined to current resolution
	real nx0=dom.icoord[nb][0]*2;
	real ny0=dom.icoord[nb][1]*2;
	//  Note that icoord is written such that it is the x/y position
	//  of the block in terms of the number of cells which is at the 
	//  assumed resolution of the leaf block. It is written like this 
	//  cause of how the x & y are computed over here. I can definitely
	//  correct this to make it less confusing.
	
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
