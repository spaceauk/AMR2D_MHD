#include<exception>
#include "defs.hpp"

real sum_irjr(real**** U, int nb, int k, int i1, int i2, int j1, int j2);

void boundary(meshblock &dom,real**** Q) {

	// Internal boundaries
	for (int nb=0; nb<dom.nbounds; nb++){
	  int nbs=dom.innerbounds[nb][0];	 
	  int n1,n2;
	  int ii, jj;
	  // YH: note that the BC are separated and thus, you can modify each of the boundaries (e.g., outflow).
	  switch(dom.innerbounds[nb][2]) {
		  case -1: // Left Box boundary
			  ii=dom.nxmin+(nghosts-1);
			  for (int i=0;i<nghosts;i++) {
 			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<dom.nvar;k++) {
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
                              for (int k=0;k<dom.nvar;k++) {
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
		              for (int k=0;k<dom.nvar;k++) {
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
                              for (int k=0;k<dom.nvar;k++) {
                                Q[i][j][k][nbs]=Q[i][jj][k][nbs];
                              }
                            }
			    jj--;
			  }
			  break;
		//----------------------------------------------------LEFT-----------------------------------------------------
		  case 1: // Left @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int j=0;j<dom.ny;j++) {
		            for (int k=0;k<dom.nvar;k++) {
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
			  for (int j=nghosts/2;j<=dom.ny2+nghosts/2;j++) {
			    int ny2j=dom.ny2+j-1-(nghosts-2);
			    for (int i=0;i<nghosts;i++) {
			      int i1=2*(i-nghosts/2)-nghosts/2+dom.nxmax;
			      int i2=i1+1;
			      int j1=2*j-dom.nymin; // YH: The use of -2 is only valid when NG=2. Must ensure that j1=j=nymin for general NG.
					    // But I cannot just use -nymin as my j cannot start from 1 since it will result in -ve indices. Rather j must start from nghosts/2 instead. NG should also be even as odd number will leave behind one cell instead of two for sum.
			      int j2=j1+1;
 			      for (int k=0;k<dom.nvar;k++) {
			        if (j<=dom.ny2) Q[i][j][k][nbs]=0.25*sum_irjr(Q,n1,k,i1,i2,j1,j2);
			        if (j>=dom.nymin) Q[i][ny2j][k][nbs]=0.25*sum_irjr(Q,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  for (int k=0;k<dom.nvar;k++) {
			    for (int i=0;i<nghosts;i++) {
			      // Account for non-overlapping ghost cells due to neighbouring refined block
			      for (int j=0;j<nghosts/2;j++) {				
			        Q[i][j][k][nbs]=Q[i][nghosts/2][k][nbs];
			        Q[i][dom.ny-1-j][k][nbs]=Q[i][dom.ny-1-nghosts/2][k][nbs];
			      }
			    }
			  }
			  // Neighbor's right
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts;i++) {
			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<dom.nvar;k++) {
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
		            for (int k=0;k<dom.nvar;k++) {
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
			    for (int j=nghosts;j<=dom.ny2+nghosts/2;j++) {
			      int i1=2*(i-dom.nxmax+(nghosts/2-1));
			      int i2=i1+1;
		              int j1=2*j-dom.nymin;
			      int j2=j1+1;
			      int ny2j=dom.ny2+j-1-(nghosts-2);
 			      for (int k=0;k<dom.nvar;k++) {
			        if (j<=dom.ny2) Q[i][j][k][nbs]=0.25*sum_irjr(Q,n1,k,i1,i2,j1,j2);
			        if (j>=dom.nymin) Q[i][ny2j][k][nbs]=0.25*sum_irjr(Q,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  for (int k=0;k<dom.nvar;k++) {
			    for (int i=dom.nxp1;i<dom.nx;i++) {
			      for (int j=0;j<nghosts/2;j++) {
			        Q[i][j][k][nbs]=Q[i][nghosts/2][k][nbs];
			        Q[i][dom.ny-1-j][k][nbs]=Q[i][dom.ny-1-nghosts/2][k][nbs];
			      }
			    }
			  }
			  // Neighbor's left
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts;i++) {
			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<dom.nvar;k++) {
				// YH: Rmb to correct nxm1 to nghosts/2 for general b.c.s
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
		            for (int k=0;k<dom.nvar;k++) {
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
			  for (int i=nghosts/2;i<=dom.nx2+nghosts/2;i++) {
			    int nx2i=dom.nx2+i-1-(nghosts-2);
			    for (int j=0;j<nghosts;j++) {
			      int i1=2*i-dom.nxmin;
                              int i2=i1+1;
			      int j1=2*(j-nghosts/2)-nghosts/2+dom.nymax;
			      int j2=j1+1;
 			      for (int k=0;k<dom.nvar;k++) {
			        if (i<=dom.nx2) Q[i][j][k][nbs]=0.25*sum_irjr(Q,n1,k,i1,i2,j1,j2);
			        if (i>=dom.nxmin) Q[nx2i][j][k][nbs]=0.25*sum_irjr(Q,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  for (int k=0;k<dom.nvar;k++) {
			    for (int j=0;j<nghosts;j++) {
			      for (int i=0;i<nghosts/2;i++) {
			        Q[i][j][k][nbs]=Q[nghosts/2][j][k][nbs];
 			        Q[dom.nx-1-i][j][k][nbs]=Q[dom.nx-1-nghosts/2][j][k][nbs];
			      }
			    }
			  }
			  // Neighbor's top
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts;j++) {
			    for (int i=0;i<dom.nx;i++) {
			      for (int k=0;k<dom.nvar;k++) {
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
		            for (int k=0;k<dom.nvar;k++) {
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
			  for (int i=nghosts/2;i<=dom.nx2+nghosts/2;i++) {
			    for (int j=dom.nyp1;j<dom.ny;j++) {
		              int i1=2*i-dom.nxmin;
			      int i2=i1+1;
			      int j1=2*(j-dom.nymax+(nghosts/2-1));
			      int j2=j1+1;
			      int nx2i=dom.nx2+i-1-(nghosts-2);
 			      for (int k=0;k<dom.nvar;k++) {
			        if (i<=dom.nx2) Q[i][j][k][nbs]=0.25*sum_irjr(Q,n1,k,i1,i2,j1,j2);
			        if (i>=dom.nxmin) Q[nx2i][j][k][nbs]=0.25*sum_irjr(Q,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  for (int k=0;k<dom.nvar;k++) {
			    for (int j=dom.nyp1;j<dom.ny;j++) {
			      for (int i=0;i<nghosts/2;i++) {
			        Q[i][j][k][nbs]=Q[nghosts/2][j][k][nbs];
			        Q[dom.nx-1-i][j][k][nbs]=Q[dom.nx-1-nghosts/2][j][k][nbs];
			      }
			    }
			  }
			  // Neighbor's bottom
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts;j++) {
			    for (int i=0;i<dom.nx;i++) {
			      for (int k=0;k<dom.nvar;k++) {
			        Q[i][j-1][k][n1]=Q[(i+nghosts)/2][(j+1)/2+dom.nymax-nghosts/2][k][nbs];
			        Q[i][j-1][k][n2]=Q[(i-nghosts+2)/2+dom.nx2][(j+1)/2+dom.nymax-nghosts/2][k][nbs];
		  	      }
			    }
			  }
			  break;
		  default:
			  cout<<"Boundary: Error as invalid direction obtained! direc="<<dom.innerbounds[nb][2]<<endl;
			  throw exception();
	  }
	}
}
