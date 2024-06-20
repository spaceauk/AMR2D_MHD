// Refinement criteria based on grad(P), grad(rho)
#include "defs.hpp"

real max4diff(real**** W,int i,int j,int k,int nb);

void criteria(meshblock &dom, int nb) {
	real grad=0;
	real rThresh=0.05; 
	real cThresh=1.0;
	
	for (int i=dom.nxmin;i<=dom.nxmax;i++) {
	  for (int j=dom.nymin;j<=dom.nymax;j++) {
	    real dm=MAX(eps,abs(dom.W[i][j][4][nb]));
	    real intm=max4diff(dom.W,i,j,4,nb);
	    grad=MAX(grad,intm/dm);
	    
	    dm=MAX(eps,abs(dom.W[i][j][0][nb]));
	    intm=max4diff(dom.W,i,j,0,nb);
	    grad=MAX(grad,intm/dm);

	    if (MAG_field) {
		    dm=MAX(eps,abs(dom.W[i][j][5][nb]));
            	    intm=max4diff(dom.W,i,j,5,nb);
            	    grad=MAX(grad,intm/dm);

		    dm=MAX(eps,abs(dom.W[i][j][6][nb]));
                    intm=max4diff(dom.W,i,j,6,nb);
                    grad=MAX(grad,intm/dm);
	    }

	    if (grad>=rThresh) {
	      // Only mark if nl<maxlevs-1 (-1 is due to cpp start from 0)
	      if (dom.lp[nb][0]<maxlevs-1) {
		cout<<"Refining at nb="<<nb<<" as grad="<<grad<<" at i="<<i<<", j="<<j<<"; # of blocks="<<dom.lastActive<<endl;
		dom.iref[nb]=true;
		dom.icoarse[nb]=false;
	      }
	      return;
	    }
	  }
	}
	if (grad<=cThresh) {dom.icoarse[nb]=true;}
}

