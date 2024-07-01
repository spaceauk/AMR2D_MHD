#include<exception>
#include "defs.hpp"


void extboundary(meshblock &dom,real**** Q,int nvar) {

	// Internal boundaries
	#pragma omp parallel for default(none) shared(dom,Q,nvar)
	for (int nb=0; nb<dom.nbounds; nb++){
	  int nbs=dom.innerbounds[nb][0];	 
	  int n1,n2;
	  int ii, jj;
	  int n3, n4;
	  switch(dom.innerbounds[nb][2]) {
		  case -1: // Left Box boundary
			  ii=dom.nxmin+(nghosts-1);
			  for (int i=0;i<nghosts;i++) {
 			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<nvar;k++) {
			        Q[i][j][k][nbs]=Q[ii][j][k][nbs];
			      } 
			    }
			    ii--;
			  }
			  break;
		  case -2: // Right Box boundary
			  ii=dom.nxmax;
			  for (int i=dom.nxp1;i<dom.nx;i++) {
 			    for (int j=0;j<dom.ny;j++) {
                              for (int k=0;k<nvar;k++) {
                                Q[i][j][k][nbs]=Q[ii][j][k][nbs];
                              }
                            }
			    ii--;
			  }
			  break;
		  case -3: // Bottom Box boundary
			  jj=dom.nymin+(nghosts-1);
			  for (int j=0;j<nghosts;j++) {
 			    for (int i=0;i<dom.nx;i++) {
		              for (int k=0;k<nvar;k++) {
		                Q[i][j][k][nbs]=Q[i][jj][k][nbs];
			      }
			    }
			    jj--;
			  }
			  break;
		  case -4: // Top Box boundary
			  jj=dom.nymax;
			  for (int j=dom.nyp1;j<dom.ny;j++) {
                            for (int i=0;i<dom.nx;i++) {				    			      
                              for (int k=0;k<nvar;k++) {
                                Q[i][j][k][nbs]=Q[i][jj][k][nbs];
                              }
                            }
			    jj--;
			  }
			  break;
	  }
	}
}
