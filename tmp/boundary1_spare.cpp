#include<exception>
#include "defs.hpp"

real sum_irjr(real**** U, int nb, int k, int i1, int i2, int j1, int j2);

void boundary1(meshblock &dom) {

	// Internal boundaries
	for (int nb=0; nb<dom.nbounds; nb++){
	  int nbs=dom.innerbounds[nb][0];	 
	  int n1,n2;
	  int ii,jj;
	  // YH: note that the BC are separated and thus, you can modify each of the boundaries (e.g., outflow).
	  switch(dom.innerbounds[nb][2]) {
		  case -1: // Left Box boundary			  
			  for (int j=0;j<dom.ny;j++) {
			    for (int k=0;k<dom.nvar;k++) {
			      for (int i=0;i<nghosts;i++) {
				ii=dom.nxmin+nghosts-i-1;
			        dom.U[i][j][k][nbs]=dom.U[ii][j][k][nbs];
			      }
			    }
			  }
			  break;
		  case -2: // Right Box boundary
			  for (int j=0;j<dom.ny;j++) {
                            for (int k=0;k<dom.nvar;k++) {
			      for (int i=0;i<nghosts;i++) {
			        ii=dom.nxmax-i;
                                dom.U[dom.nxmaxb+i][j][k][nbs]=dom.U[ii][j][k][nbs];
			      }
                            }
                          }
			  break;
		  case -3: // Bottom Box boundary
			  for (int i=0;i<dom.nx;i++) {
		            for (int k=0;k<dom.nvar;k++) {
			      for (int j=0;j<nghosts;j++) {
				jj=dom.nymin+nghosts-j-1;
		                dom.U[i][j][k][nbs]=dom.U[i][jj][k][nbs];
			      }
			    }
			  }
			  break;
		  case -4: // Top Box boundary
                          for (int i=0;i<dom.nx;i++) {
                            for (int k=0;k<dom.nvar;k++) {
			      for (int j=0;j<nghosts;j++) {
				jj=dom.nymax-j;
                                dom.U[i][dom.nymaxb+j][k][nbs]=dom.U[i][jj][k][nbs];
			      }
                            }
                          }
			  break;
		//----------------------------------------------------LEFT-----------------------------------------------------
		  case 1: // Left @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int j=0;j<dom.ny;j++) {
		            for (int k=0;k<dom.nvar;k++) {
			      for (int i=0;i<nghosts;i++) {
			        // Own's left
				ii=dom.nxmax-(nghosts-1)+i;
			        dom.U[i][j][k][nbs]=dom.U[ii][j][k][n1];
			        // Neighbor's right
			        dom.U[dom.nxp1+i][j][k][n1]=dom.U[dom.nxmin+i][j][k][nbs];
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
                            for (int i=1;i<=nghosts-1;i++) {
                              int i1=2*(i-1)-1+dom.nxmax;
                              int i2=i1+1;
                              int j1=2*j-dom.nymin;
                              int j2=j1+1;
                              for (int k=0;k<dom.nvar;k++) {
                                if (j<=dom.ny2) dom.U[i][j][k][nbs]=0.25*sum_irjr(dom.U,n1,k,i1,i2,j1,j2);
                                if (j>=dom.nymin) dom.U[i][ny2j][k][nbs]=0.25*sum_irjr(dom.U,n2,k,i1,i2,j1,j2);
                              }
                            }
                          }
			  // Neighbor's right
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts-1;i++) {
                            for (int j=0;j<dom.ny;j++) {
                              for (int k=0;k<dom.nvar;k++) {
                                dom.U[i+dom.nxmax][j][k][n1]=dom.U[(i+3)/2][(j+2)/2][k][nbs];
                                dom.U[i+dom.nxmax][j][k][n2]=dom.U[(i+3)/2][(j)/2+dom.ny2][k][nbs];
                              }
                            }
                          }
			  break;
	  	//----------------------------------------------------RIGHT-----------------------------------------------------
		  case 4: // Right @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int j=0;j<dom.ny;j++) {
		            for (int k=0;k<dom.nvar;k++) {
			      for (int i=0;i<nghosts;i++) {
			        // Own's right
			        dom.U[dom.nxp1+i][j][k][nbs]=dom.U[dom.nxmin+i][j][k][n1];
			        // Neighbor's left
				ii=dom.nxmax-(nghosts-1)+i;
			        dom.U[i][j][k][n1]=dom.U[ii][j][k][nbs];
			      }
			    }
			  }
			  break;
		  case 5: // Right @ higher resolution 1
			  n1=dom.innerbounds[nb][1];      // nb | n2 = ID1+2
			  n2=n1+2;                        // nb | n1 = ID1
			  // Own's right
			  // issue is when set i<dom.nx 
			  for (int i=dom.nxp1;i<dom.nx-1;i++) {
                            for (int j=1;j<=dom.ny2+(nghosts-1);j++) {
                              int i1=2*(i-dom.nxmax);
                              int i2=i1+1;
                              int j1=2*j-dom.nymin;
                              int j2=j1+1;
                              int ny2j=dom.ny2+j-1;
                              for (int k=0;k<dom.nvar;k++) {
                                if (j<=dom.ny2) dom.U[i][j][k][nbs]=0.25*sum_irjr(dom.U,n1,k,i1,i2,j1,j2);
                                if (j>=dom.nymin) dom.U[i][ny2j][k][nbs]=0.25*sum_irjr(dom.U,n2,k,i1,i2,j1,j2);
                              }
                            }
                          }
			  // Neighbor's left
			  // to finer level (0th order interpolation)
			  for (int i=1+1;i<=nghosts;i++) {
                            for (int j=0;j<dom.ny;j++) {
                              for (int k=0;k<dom.nvar;k++) {
                                dom.U[i-1][j][k][n1]=dom.U[(i+1)/2+dom.nxm1][(j+2)/2][k][nbs];
                                dom.U[i-1][j][k][n2]=dom.U[(i+1)/2+dom.nxm1][(j)/2+dom.ny2][k][nbs];
                              }
                            }
                          }
			  break;
		//----------------------------------------------------BOTTOM-----------------------------------------------------
		  case 7: // Bottom @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int i=0;i<dom.nx;i++) {
		            for (int k=0;k<dom.nvar;k++) {
			      for (int j=0;j<nghosts;j++) {
			        // Own's bottom
				jj=dom.nymax-(nghosts-1)+j;
			        dom.U[i][j][k][nbs]=dom.U[i][jj][k][n1];
			        // Neighbor's top
			        dom.U[i][dom.nyp1+j][k][n1]=dom.U[i][dom.nymin+j][k][nbs];
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
                            for (int j=0+1;j<nghosts;j++) {
                              int i1=2*i-dom.nxmin;
                              int i2=i1+1;
                              int j1=2*(j-1)-1+dom.nymax;
                              int j2=j1+1;
                              for (int k=0;k<dom.nvar;k++) {
                                if (i<=dom.nx2) dom.U[i][j][k][nbs]=0.25*sum_irjr(dom.U,n1,k,i1,i2,j1,j2);
                                if (i>=dom.nxmin) dom.U[nx2i][j][k][nbs]=0.25*sum_irjr(dom.U,n2,k,i1,i2,j1,j2);
                              }
                            }
                          }
			  // Neighbor's top
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts-1;j++) {
                            for (int i=0;i<dom.nx;i++) {
                              for (int k=0;k<dom.nvar;k++) {
                                dom.U[i][j+dom.nymax][k][n1]=dom.U[(i+2)/2][(j+3)/2][k][nbs];
                                dom.U[i][j+dom.nymax][k][n2]=dom.U[(i)/2+dom.nx2][(j+3)/2][k][nbs];
                              }
                            }
                          }
			  break;
		//----------------------------------------------------TOP-----------------------------------------------------
		  case 10: // TOP @ same resolution
			  n1=dom.innerbounds[nb][1];	 
			  for (int i=0;i<dom.nx;i++) {
		            for (int k=0;k<dom.nvar;k++) {
			      for (int j=0;j<nghosts;j++) {
			        // Own's top
			        dom.U[i][dom.nyp1+j][k][nbs]=dom.U[i][dom.nymin+j][k][n1];
			        // Neighbor's bottom
				jj=dom.nymax-(nghosts-1)+j;
			        dom.U[i][j][k][n1]=dom.U[i][jj][k][nbs];
			      }
			    }
			  }
			  break;
		  case 11: // TOP @ higher resolution 1
			  n1=dom.innerbounds[nb][1];
			  n2=n1+1;
			  // Own's top
			  for (int i=1;i<=dom.nx2+(nghosts-1);i++) {
                            for (int j=dom.nyp1;j<dom.ny-1;j++) {
                              int i1=2*i-dom.nxmin;
                              int i2=i1+1;
                              int j1=2*(j-dom.nymax);
                              int j2=j1+1;
                              int nx2i=dom.nx2+i-1;
                              for (int k=0;k<dom.nvar;k++) {
                                if (i<=dom.nx2) dom.U[i][j][k][nbs]=0.25*sum_irjr(dom.U,n1,k,i1,i2,j1,j2);
                                if (i>=dom.nxmin) dom.U[nx2i][j][k][nbs]=0.25*sum_irjr(dom.U,n2,k,i1,i2,j1,j2);
                              }
                            }
                          }
			  // Neighbor's bottom
			  // to finer level (0th order interpolation)
			  for (int j=1+1;j<=nghosts;j++) {
                            for (int i=0;i<dom.nx;i++) {
                              for (int k=0;k<dom.nvar;k++) {
                                dom.U[i][j-1][k][n1]=dom.U[(i+2)/2][(j+1)/2+dom.nym1][k][nbs];
                                dom.U[i][j-1][k][n2]=dom.U[(i)/2+dom.nx2][(j+1)/2+dom.nym1][k][nbs];
                              }
                            }
                          }
			  break;
		  default:
			  cout<<"Boundary 1: Error as invalid direction obtained! direc="<<dom.innerbounds[nb][2]<<endl;
			  throw exception();
	  }
	}
}
