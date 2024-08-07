#include<exception>
#include "defs.hpp"

real sum_irjr(real**** U, int nb, int k, int i1, int i2, int j1, int j2);
void extboundary(meshblock &dom,real**** Q,int nvar);

void boundary(meshblock &dom,real**** Q,int nvar) {
	int nbs;
	int n1,n2,n3,n4; 
	int ii,jj;

	// Internal boundaries
	for (int nb=0; nb<dom.nbounds; nb++){
	  nbs=dom.innerbounds[nb][0];	 
	  switch(dom.innerbounds[nb][2]) {
		//----------------------------------------------------LEFT-----------------------------------------------------
		  case 1: // Left @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int j=0;j<dom.ny;j++) {
		            for (int k=0;k<nvar;k++) {
			      // Own's left
			      ii=dom.nxmax-(nghosts-1);
			      for (int i=0;i<nghosts;i++) {
			        Q[i][j][k][nbs]=Q[ii][j][k][n1];
				ii++;
			      }
			      // Neighbor's right
			      ii=dom.nxmin;
			      for (int i=dom.nxp1;i<dom.nx;i++) {
			        Q[i][j][k][n1]=Q[ii][j][k][nbs];
			        ii++;
			      }
			    }
			  }
			  break;
		  case 2: // Left @ higher resolution 1
			  n1=dom.innerbounds[nb][1];
			  n2=n1+2;
			  // Own's left
			  for (int j=dom.nymin;j<=dom.ny2;j++) {
			    int ny2j=dom.ny2+j-1-(nghosts-2);
			    for (int i=0;i<nghosts;i++) {
			      int i1=2*(i-nghosts/2)-nghosts/2+dom.nxmax;
			      int i2=i1+1;
			      int j1=2*j-dom.nymin; 
			      int j2=j1+1;
 			      for (int k=0;k<nvar;k++) {
			        Q[i][j][k][nbs]=0.25*sum_irjr(Q,n1,k,i1,i2,j1,j2);
			        Q[i][ny2j][k][nbs]=0.25*sum_irjr(Q,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  n3=dom.lp[dom.lp[n1][1]][6]; // Bottom corner 
			  // 1. Note has to use the coarser version of bottom blocks as the upper block is less refined.
			  // 2. Do not need to update bottom corner if the neighbouring corner is coarser since the blocks already have all the data from coarser blocks.
			  if (n3!=-1) { // Ensure doesn't lie on grid boundary			    
			   if (dom.leafs[n3]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=0;j<nghosts;j++) {
                                    ii=dom.nxmax-nghosts+(i+1); jj=dom.nymax-nghosts+(j+1);
                                    Q[i][j][k][nbs]=Q[ii][jj][k][n3];
                                  }}}
			   } else {
			      n3=dom.lp[n3][2]+3;
			      if (dom.leafs[n3]) { // Means more refined corner
			        for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxmax-nghosts*2+2*i+1; jj=dom.nymax-nghosts*2+2*j+1;
                                      Q[i][j][k][nbs]=0.25*sum_irjr(Q,n3,k,ii,ii+1,jj,jj+1);
                                    }}}
			      } 
			    }
			  }
			  n4=dom.lp[dom.lp[n2][1]][7]; // Top corner
			  if (n4!=-1) { // Ensure does not lie at grid boundary
                            if (dom.leafs[n4]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=0;j<nghosts;j++) {
                                    ii=dom.nxmax-nghosts+(i+1); jj=dom.nymin+nghosts-(j+1);
                                    Q[i][dom.ny-1-j][k][nbs]=Q[ii][jj][k][n4];
                                  }}}
                            } else  {
                              n4=dom.lp[n4][2]+1;
			      if (dom.leafs[n4]) { // Means more refined corner
                                for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxmax-nghosts*2+2*i+1; jj=dom.nymin+nghosts*2-2*(j+1);
                                      Q[i][dom.ny-1-j][k][nbs]=0.25*sum_irjr(Q,n4,k,ii,ii+1,jj,jj+1);
                                    }}}
                              } 
			    }
                          }			  
			  // Neighbor's right
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts;i++) {
			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<nvar;k++) {
			        Q[i+dom.nxmax][j][k][n1]=Q[(i+2*nghosts-1)/2][(j+nghosts)/2][k][nbs];
			        Q[i+dom.nxmax][j][k][n2]=Q[(i+2*nghosts-1)/2][(j-nghosts+2)/2+dom.ny2][k][nbs];
		  	      }
			    }
			  }
			  break;
	  	//----------------------------------------------------RIGHT-----------------------------------------------------
		  case 4: // Right @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int j=0;j<dom.ny;j++) {
		            for (int k=0;k<nvar;k++) {
			      // Own's right
			      ii=dom.nxmin;
			      for (int i=dom.nxp1;i<dom.nx;i++) {
			        Q[i][j][k][nbs]=Q[ii][j][k][n1];
				ii++;
			      }
			      // Neighbor's left
			      ii=dom.nxmax-(nghosts-1);
			      for (int i=0;i<nghosts;i++) {
			        Q[i][j][k][n1]=Q[ii][j][k][nbs];
				ii++;
			      }
			    }
			  }
			  break;
		  case 5: // Right @ higher resolution 1
			  n1=dom.innerbounds[nb][1];
			  n2=n1+2;
			  // Own's right
			  for (int i=dom.nxp1;i<dom.nx;i++) {
			    for (int j=dom.nymin;j<=dom.ny2;j++) {
			      int i1=2*(i-dom.nxmax+(nghosts/2-1));
			      int i2=i1+1;
		              int j1=2*j-dom.nymin;
			      int j2=j1+1;
			      int ny2j=dom.ny2+j-1-(nghosts-2);
 			      for (int k=0;k<nvar;k++) {
			        Q[i][j][k][nbs]=0.25*sum_irjr(Q,n1,k,i1,i2,j1,j2);
			        Q[i][ny2j][k][nbs]=0.25*sum_irjr(Q,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  n3=dom.lp[dom.lp[n1][1]][6]; // Bottom corner
			  if (n3!=-1) { // Ensure not at grid boundary
			    if (dom.leafs[n3]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                  for (int j=0;j<nghosts;j++) {
                                    ii=dom.nxminb+(i-dom.nxmax); jj=dom.nymax-nghosts+(j+1);
                                    Q[i][j][k][nbs]=Q[ii][jj][k][n3];
                                  }}}
                            } else {
			      n3=dom.lp[n3][2]+2;
			      if (dom.leafs[n3]) { // Means more refined corner
			        for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxminb+2*(i-dom.nxmax)-1; jj=dom.nymax-nghosts*2+2*j+1;
                                      Q[i][j][k][nbs]=0.25*sum_irjr(Q,n3,k,ii,ii+1,jj,jj+1);
                                    }}}
			      } 
			    }			    
			  }			  
			  n4=dom.lp[dom.lp[n2][1]][7]; // Top corner
			  if (n4!=-1) { // ensure not at grid boundary
                            if (dom.leafs[n4]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                  for (int j=0;j<nghosts;j++) {
                                    ii=dom.nxminb+(i-dom.nxmax); jj=dom.nymin+nghosts-(j+1);
                                    Q[i][dom.ny-1-j][k][nbs]=Q[ii][jj][k][n4];
                                  }}}
                            } else {
                              n4=dom.lp[n4][2];
			      if (dom.leafs[n4]) { // Means more refined corner
                                for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxminb+2*(i-dom.nxmax)-1; jj=dom.nymin+nghosts*2-2*(j+1);
                                      Q[i][dom.ny-1-j][k][nbs]=0.25*sum_irjr(Q,n4,k,ii,ii+1,jj,jj+1);
                                    }}}
                              } 
			    }
                          }
			  // Neighbor's left
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts;i++) {
			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<nvar;k++) {
			        Q[i-1][j][k][n1]=Q[(i+1)/2+dom.nxmax-nghosts/2][(j+nghosts)/2][k][nbs];
			        Q[i-1][j][k][n2]=Q[(i+1)/2+dom.nxmax-nghosts/2][(j-nghosts+2)/2+dom.ny2][k][nbs];
		  	      }
			    }
			  }
			  break;
		//----------------------------------------------------BOTTOM-----------------------------------------------------
		  case 7: // Bottom @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int i=0;i<dom.nx;i++) {
		            for (int k=0;k<nvar;k++) {
			      // Own's bottom
			      jj=dom.nymax-(nghosts-1);
			      for (int j=0;j<nghosts;j++) {
			        Q[i][j][k][nbs]=Q[i][jj][k][n1];
				jj++;
			      }
			      // Neighbor's top
			      jj=dom.nymin;
			      for (int j=dom.nyp1;j<dom.ny;j++) {
			        Q[i][j][k][n1]=Q[i][jj][k][nbs];
				jj++;
			      }
			    }
			  }
			  break;
		  case 8: // Bottom @ higher resolution 1
			  n1=dom.innerbounds[nb][1];
			  n2=n1+1;
			  // Own's bottom
			  for (int i=dom.nxmin;i<=dom.nx2;i++) {
			    int nx2i=dom.nx2+i-1-(nghosts-2);
			    for (int j=0;j<nghosts;j++) {
			      int i1=2*i-dom.nxmin;
                              int i2=i1+1;
			      int j1=2*(j-nghosts/2)-nghosts/2+dom.nymax;
			      int j2=j1+1;
 			      for (int k=0;k<nvar;k++) {
			        Q[i][j][k][nbs]=0.25*sum_irjr(Q,n1,k,i1,i2,j1,j2);
			        Q[nx2i][j][k][nbs]=0.25*sum_irjr(Q,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  n3=dom.lp[dom.lp[n1][1]][4]; // left corner
			  if (n3!=-1) { // ensure does not lie on grid boundary
			    if (dom.leafs[n3]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=0;j<nghosts;j++) {
                                    ii=dom.nxmax-nghosts+(i+1); jj=dom.nymax-nghosts+(j+1);
                                    Q[i][j][k][nbs]=Q[ii][jj][k][n3];
                                  }}}
                            } else {
			      n3=dom.lp[n3][2]+3;
			      if (dom.leafs[n3]) { // Means more refined corner
			        for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxmax-nghosts*2+2*i+1; jj=dom.nymax-nghosts*2+2*j+1;
                                      Q[i][j][k][nbs]=0.25*sum_irjr(Q,n3,k,ii,ii+1,jj,jj+1);
                                    }}}
			      } 
			    }
			  }
			  n4=dom.lp[dom.lp[n2][1]][5]; // right corner
			  if (n4!=-1) { // ensure not at grid boundary
                            if (dom.leafs[n4]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=0;j<nghosts;j++) {
                                    ii=dom.nxmin+nghosts-(i+1); jj=dom.nymax-nghosts+(j+1);
                                    Q[dom.nx-1-i][j][k][nbs]=Q[ii][jj][k][n4];
                                  }}}
                            } else {
                              n4=dom.lp[n4][2]+2;
			      if (dom.leafs[n4]) { // Means more refined corner
                                for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
				      ii=dom.nxmin+nghosts*2-2*(i+1); jj=dom.nymax-nghosts*2+2*j+1;
                                      Q[dom.nx-1-i][j][k][nbs]=0.25*sum_irjr(Q,n4,k,ii,ii+1,jj,jj+1);
                                    }}}
                              } 
			    }
                          }
			  // Neighbor's top
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts;j++) {
			    for (int i=0;i<dom.nx;i++) {
			      for (int k=0;k<nvar;k++) {
			        Q[i][j+dom.nymax][k][n1]=Q[(i+nghosts)/2][(j+2*nghosts-1)/2][k][nbs];
			        Q[i][j+dom.nymax][k][n2]=Q[(i-nghosts+2)/2+dom.nx2][(j+2*nghosts-1)/2][k][nbs];
		  	      }
			    }
			  }
			  break;
		//----------------------------------------------------TOP-----------------------------------------------------
		  case 10: // TOP @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int i=0;i<dom.nx;i++) {
		            for (int k=0;k<nvar;k++) {
			      // Own's top
			      jj=dom.nymin;
			      for (int j=dom.nyp1;j<dom.ny;j++) {
			        Q[i][j][k][nbs]=Q[i][jj][k][n1];
				jj++;
			      }
			      // Neighbor's bottom
			      jj=dom.nymax-(nghosts-1);
			      for (int j=0;j<nghosts;j++) {
			        Q[i][j][k][n1]=Q[i][jj][k][nbs];
				jj++;
			      }
			    }
			  }
			  break;
		  case 11: // TOP @ higher resolution 1
			  n1=dom.innerbounds[nb][1];
			  n2=n1+1;
			  // Own's top
			  for (int i=dom.nxmin;i<=dom.nx2;i++) {
			    for (int j=dom.nyp1;j<dom.ny;j++) {
		              int i1=2*i-dom.nxmin;
			      int i2=i1+1;
			      int j1=2*(j-dom.nymax+(nghosts/2-1));
			      int j2=j1+1;
			      int nx2i=dom.nx2+i-1-(nghosts-2);
 			      for (int k=0;k<nvar;k++) {
			        Q[i][j][k][nbs]=0.25*sum_irjr(Q,n1,k,i1,i2,j1,j2);
			        Q[nx2i][j][k][nbs]=0.25*sum_irjr(Q,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  n3=dom.lp[dom.lp[n1][1]][4]; // left corner
			  if (n3!=-1) { // ensure not at grid boundary
			    if (dom.leafs[n3]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=dom.nymaxb;j<dom.ny;j++) {
                                    ii=dom.nxmax-nghosts+(i+1); jj=dom.nyminb+(j-dom.nymax);
                                    Q[i][j][k][nbs]=Q[ii][jj][k][n3];
                                  }}}
                            } else {
			      n3=dom.lp[n3][2]+1;
			      if (dom.leafs[n3]) { // Means more refined corner
			      for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=dom.nymaxb;j<dom.ny;j++) {
				    ii=dom.nxmax-nghosts*2+2*i+1; jj=dom.nyminb+2*(j-dom.nymax)-1;
                                    Q[i][j][k][nbs]=0.25*sum_irjr(Q,n3,k,ii,ii+1,jj,jj+1);
                                  }}}
			      } 
			    }			      
			  }
			  n4=dom.lp[dom.lp[n2][1]][5]; // right corner
			  if (n4!=-1) { // ensure not at grid boundary
                            if (dom.leafs[n4]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=dom.nymaxb;j<dom.ny;j++) {
                                    ii=dom.nxmin+nghosts-(i+1); jj=dom.nyminb+(j-dom.nymax);
                                    Q[dom.nx-1-i][j][k][nbs]=Q[ii][jj][k][n4];
                                  }}}
			    } else {
                              n4=dom.lp[n4][2];
			      if (dom.leafs[n4]) { // Means more refined corner
                                for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=dom.nymaxb;j<dom.ny;j++) {
                                      ii=dom.nxminb+2*(i-dom.nxmax)-1; jj=dom.nyminb+2*(j-dom.nymax)-1;
                                      Q[i][j][k][nbs]=0.25*sum_irjr(Q,n4,k,ii,ii+1,jj,jj+1);
                                    }}}
                              } 
			    }
                          }
			  // Neighbor's bottom
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts;j++) {
			    for (int i=0;i<dom.nx;i++) {
			      for (int k=0;k<nvar;k++) {
			        Q[i][j-1][k][n1]=Q[(i+nghosts)/2][(j+1)/2+dom.nymax-nghosts/2][k][nbs];
			        Q[i][j-1][k][n2]=Q[(i-nghosts+2)/2+dom.nx2][(j+1)/2+dom.nymax-nghosts/2][k][nbs];
		  	      }
			    }
			  }
			  break;
		  default:
			  if (dom.innerbounds[nb][2]<-4 and dom.innerbounds[nb][2]>-1) {
			  	printf("Boundary: Error as invalid direction obtained! direc=%d \n",dom.innerbounds[nb][2]);
			  	throw exception();
			  }
	  }
	}

	// External boundaries
	extboundary(dom,Q,nvar);
}
