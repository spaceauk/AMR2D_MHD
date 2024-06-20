#include<exception>
#include "defs.hpp"

real sum_irjr(real**** U, int nb, int k, int i1, int i2, int j1, int j2);

void boundary2(meshblock &dom) {

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
			        dom.Us[i][j][k][nbs]=dom.Us[ii][j][k][nbs];
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
                                dom.Us[i][j][k][nbs]=dom.Us[ii][j][k][nbs];
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
		                dom.Us[i][j][k][nbs]=dom.Us[i][jj][k][nbs];
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
                                dom.Us[i][j][k][nbs]=dom.Us[i][jj][k][nbs];
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
			        dom.Us[i][j][k][nbs]=dom.Us[ii][j][k][n1];
				ii++;
			      }
			      // Neighbor's right
			      ii=dom.nxmin;
			      for (int i=dom.nxp1;i<dom.nx;i++) {
			        dom.Us[i][j][k][n1]=dom.Us[ii][j][k][nbs];
			        ii++;
			      }
			    }
			  }
			  break;
		  case 2: // Left @ higher resolution 1
			  n1=dom.innerbounds[nb][1];
			  n2=n1+2;
			  // Own's left
			  for (int j=1;j<=dom.ny2+(nghosts-1);j++) {
			    int ny2j=dom.ny2+j-1;
			    for (int i=0;i<nghosts;i++) {
			      int i1=2*(i-1)-1+dom.nxmax;
			      int i2=i1+1;
			      int j1=2*j-dom.nymin;
			      int j2=j1+1;
 			      for (int k=0;k<dom.nvar;k++) {
			        if (j<=dom.ny2) dom.Us[i][j][k][nbs]=0.25*sum_irjr(dom.Us,n1,k,i1,i2,j1,j2);
			        if (j>=dom.nymin) dom.Us[i][ny2j][k][nbs]=0.25*sum_irjr(dom.Us,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  for (int k=0;k<dom.nvar;k++) {
			    for (int i=0;i<nghosts;i++) {
			      dom.Us[i][0][k][nbs]=dom.Us[i][1][k][nbs];
			      dom.Us[i][dom.ny-1][k][nbs]=dom.Us[i][dom.ny-2][k][nbs];
			    }
			  }
			  // Neighbor's right
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts;i++) {
			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<dom.nvar;k++) {
			        dom.Us[i+dom.nxmax][j][k][n1]=dom.Us[(i+3)/2][(j+2)/2][k][nbs];
			        dom.Us[i+dom.nxmax][j][k][n2]=dom.Us[(i+3)/2][(j)/2+dom.ny2][k][nbs];
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
			        dom.Us[i][j][k][nbs]=dom.Us[ii][j][k][n1];
				ii++;
			      }
			      // Neighbor's left
			      ii=dom.nxmax-(nghosts-1);
			      for (int i=0;i<nghosts;i++) {
			        dom.Us[i][j][k][n1]=dom.Us[ii][j][k][nbs];
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
			    for (int j=1;j<=dom.ny2+(nghosts-1);j++) {
			      int i1=2*(i-dom.nxmax);
			      int i2=i1+1;
		              int j1=2*j-dom.nymin;
			      int j2=j1+1;
			      int ny2j=dom.ny2+j-1;
 			      for (int k=0;k<dom.nvar;k++) {
			        if (j<=dom.ny2) dom.Us[i][j][k][nbs]=0.25*sum_irjr(dom.Us,n1,k,i1,i2,j1,j2);
			        if (j>=dom.nymin) dom.Us[i][ny2j][k][nbs]=0.25*sum_irjr(dom.Us,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  for (int k=0;k<dom.nvar;k++) {
			    for (int i=dom.nxp1;i<dom.nx;i++) {
			      dom.Us[i][0][k][nbs]=dom.Us[i][1][k][nbs];
			      dom.Us[i][dom.ny-1][k][nbs]=dom.Us[i][dom.ny-2][k][nbs];
			    }
			  }
			  // Neighbor's left
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts;i++) {
			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<dom.nvar;k++) {
			        dom.Us[i-1][j][k][n1]=dom.Us[(i+1)/2+dom.nxm1][(j+2)/2][k][nbs];
			        dom.Us[i-1][j][k][n2]=dom.Us[(i+1)/2+dom.nxm1][(j)/2+dom.ny2][k][nbs];
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
			        dom.Us[i][j][k][nbs]=dom.Us[i][jj][k][n1];
				jj++;
			      }
			      // Neighbor's top
			      jj=dom.nymin;
			      for (int j=dom.nyp1;j<dom.ny;j++) {
			        dom.Us[i][j][k][n1]=dom.Us[i][jj][k][nbs];
				jj++;
			      }
			    }
			  }
			  break;
		  case 8: // Bottom @ higher resolution 1
			  n1=dom.innerbounds[nb][1];
			  n2=n1+1;
			  // Own's bottom
			  for (int i=1;i<=dom.nx2+(nghosts-1);i++) {
			    int nx2i=dom.nx2+i-1;
			    for (int j=0;j<nghosts;j++) {
			      int i1=2*i-dom.nxmin;
                              int i2=i1+1;
			      int j1=2*(j-1)-1+dom.nymax;
			      int j2=j1+1;
 			      for (int k=0;k<dom.nvar;k++) {
			        if (i<=dom.nx2) dom.Us[i][j][k][nbs]=0.25*sum_irjr(dom.Us,n1,k,i1,i2,j1,j2);
			        if (i>=dom.nxmin) dom.Us[nx2i][j][k][nbs]=0.25*sum_irjr(dom.Us,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  for (int k=0;k<dom.nvar;k++) {
			    for (int j=0;j<nghosts;j++) {
			      dom.Us[0][j][k][nbs]=dom.Us[1][j][k][nbs];
 			      dom.Us[dom.nx-1][j][k][nbs]=dom.Us[dom.nx-2][j][k][nbs];
			    }
			  }
			  // Neighbor's top
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts;j++) {
			    for (int i=0;i<dom.nx;i++) {
			      for (int k=0;k<dom.nvar;k++) {
			        dom.Us[i][j+dom.nymax][k][n1]=dom.Us[(i+2)/2][(j+3)/2][k][nbs];
			        dom.Us[i][j+dom.nymax][k][n2]=dom.Us[(i)/2+dom.nx2][(j+3)/2][k][nbs];
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
			        dom.Us[i][j][k][nbs]=dom.Us[i][jj][k][n1];
				jj++;
			      }
			      // Neighbor's bottom
			      jj=dom.nymax-(nghosts-1);
			      for (int j=0;j<nghosts;j++) {
			        dom.Us[i][j][k][n1]=dom.Us[i][jj][k][nbs];
				jj++;
			      }
			    }
			  }
			  break;
		  case 11: // TOP @ higher resolution 1
			  n1=dom.innerbounds[nb][1];
			  n2=n1+1;
			  // Own's top
			  for (int i=1;i<=dom.nx2+(nghosts-1);i++) {
			    for (int j=dom.nyp1;j<dom.ny;j++) {
		              int i1=2*i-dom.nxmin;
			      int i2=i1+1;
			      int j1=2*(j-dom.nymax);
			      int j2=j1+1;
			      int nx2i=dom.nx2+i-1;
 			      for (int k=0;k<dom.nvar;k++) {
			        if (i<=dom.nx2) dom.Us[i][j][k][nbs]=0.25*sum_irjr(dom.Us,n1,k,i1,i2,j1,j2);
			        if (i>=dom.nxmin) dom.Us[nx2i][j][k][nbs]=0.25*sum_irjr(dom.Us,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  for (int k=0;k<dom.nvar;k++) {
			    for (int j=dom.nyp1;j<dom.ny;j++) {
			      dom.Us[0][j][k][nbs]=dom.Us[1][j][k][nbs];
			      dom.Us[dom.nx-1][j][k][nbs]=dom.Us[dom.nx-2][j][k][nbs];
			    }
			  }
			  // Neighbor's bottom
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts;j++) {
			    for (int i=0;i<dom.nx;i++) {
			      for (int k=0;k<dom.nvar;k++) {
			        dom.Us[i][j-1][k][n1]=dom.Us[(i+2)/2][(j+1)/2+dom.nym1][k][nbs];
			        dom.Us[i][j-1][k][n2]=dom.Us[(i)/2+dom.nx2][(j+1)/2+dom.nym1][k][nbs];
		  	      }
			    }
			  }
			  break;
		  default:
			  cout<<"Boundary2: Error as invalid direction obtained! direc="<<dom.innerbounds[nb][2]<<endl;
			  throw exception();
	  }
	}
}
