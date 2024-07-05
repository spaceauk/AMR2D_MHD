#include<exception>
#include "defs.hpp"

real sum_irjr(real**** U, int nb, int k, int i1, int i2, int j1, int j2);
void extboundary(meshblock &dom,real**** Bi,int nvar);

// fc for face centered
void boundaryfc(meshblock &dom,real**** Bi,int nvar) {
	int nbs;
	int n1,n2,n3,n4; 
	int ii,jj;
	int i1, i2, j1, j2;

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
			        Bi[i][j][k][nbs]=Bi[ii][j][k][n1];
				ii++;
			      }
			      // Neighbor's right
			      ii=dom.nxmin;
			      for (int i=dom.nxp1;i<dom.nx;i++) {
			        Bi[i][j][k][n1]=Bi[ii][j][k][nbs];
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
			      i1=2*(i-nghosts/2)-nghosts/2+dom.nxmax;
			      j1=2*j-dom.nymin; 
 			      for (int k=0;k<nvar;k++) {
				i2= k==0 ? i1 : i1+1;
				j2= k==0 ? j1+1 : j1;
			        Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n1,k,i1,i2,j1,j2);
			        Bi[i][ny2j][k][nbs]=0.5*sum_irjr(Bi,n2,k,i1,i2,j1,j2);
			      }
			    }
			  }
			  n3=dom.lp[dom.lp[n1][1]][6]; // Bottom corner
			  if (n3!=-1) { // Ensure doesn't lie on grid boundary			    
			   if (dom.leafs[n3]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=0;j<nghosts;j++) {
                                    ii=dom.nxmax-nghosts+(i+1); jj=dom.nymax-nghosts+(j+1);
                                    Bi[i][j][k][nbs]=Bi[ii][jj][k][n3];
                                  }}}
			   } else {
			      n3=dom.lp[n3][2]+3;
			      if (dom.leafs[n3]) { // Means more refined corner
			        for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      i1=dom.nxmax-nghosts*2+2*i+1; j1=dom.nymax-nghosts*2+2*j+1;
				      i2= k==0 ? i1 : i1+1;
				      j2= k==0 ? j1+1 : j1;
                                      Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n3,k,i1,i2,j1,j2);
                                    }}}
			      } else { // coarser corner
			        n3=dom.lp[dom.lp[dom.lp[n1][1]][1]][6];
                                for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxmax-nghosts/2+(i+2)/2; jj=dom.nymax-nghosts/2+(j+2)/2;
                                      Bi[i][j][k][nbs]=Bi[ii][jj][k][n3];
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
                                    Bi[i][dom.ny-1-j][k][nbs]=Bi[ii][jj][k][n4];
                                  }}}
                            } else  {
                              n4=dom.lp[n4][2]+1;
			      if (dom.leafs[n4]) { // Means more refined corner
                                for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      i1=dom.nxmax-nghosts*2+2*i+1; j1=dom.nymin+nghosts*2-2*(j+1);
				      i2= k==0 ? i1 : i1+1;
                                      j2= k==0 ? j1+1 : j1;
                                      Bi[i][dom.ny-1-j][k][nbs]=0.5*sum_irjr(Bi,n4,k,i1,i2,j1,j2);
                                    }}}
                              } else  { // coarser corner
			        n4=dom.lp[dom.lp[dom.lp[n2][1]][1]][7];
                                for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxmax-nghosts/2+(i+2)/2; jj=dom.nymin+nghosts/2-(j+2)/2;
                                      Bi[i][dom.ny-1-j][k][nbs]=Bi[ii][jj][k][n4];
                                    }}}
			      }
			    }
                          }			  
			  // Neighbor's right
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts;i++) {
			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<nvar;k++) {
			        Bi[i+dom.nxmax][j][k][n1]=Bi[(i+2*nghosts-1)/2][(j+nghosts)/2][k][nbs];
			        Bi[i+dom.nxmax][j][k][n2]=Bi[(i+2*nghosts-1)/2][(j-nghosts+2)/2+dom.ny2][k][nbs];
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
			        Bi[i][j][k][nbs]=Bi[ii][j][k][n1];
				ii++;
			      }
			      // Neighbor's left
			      ii=dom.nxmax-(nghosts-1);
			      for (int i=0;i<nghosts;i++) {
			        Bi[i][j][k][n1]=Bi[ii][j][k][nbs];
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
			      i1=2*(i-dom.nxmax+(nghosts/2-1));
		              j1=2*j-dom.nymin;
			      int ny2j=dom.ny2+j-1-(nghosts-2);
 			      for (int k=0;k<nvar;k++) {
				i2= k==0 ? i1 : i1+1;
                                j2= k==0 ? j1+1 : j1;
			        Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n1,k,i1,i2,j1,j2);
			        Bi[i][ny2j][k][nbs]=0.5*sum_irjr(Bi,n2,k,i1,i2,j1,j2);
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
                                    Bi[i][j][k][nbs]=Bi[ii][jj][k][n3];
                                  }}}
                            } else {
			      n3=dom.lp[n3][2]+2;
			      if (dom.leafs[n3]) { // Means more refined corner
			        for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      i1=dom.nxminb+2*(i-dom.nxmax)-1; j1=dom.nymax-nghosts*2+2*j+1;
				      i2= k==0 ? i1 : i1+1;
                                      j2= k==0 ? j1+1 : j1;
                                      Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n3,k,i1,i2,j1,j2);
                                    }}}
			      } else { // coarser corner
			        n3=dom.lp[dom.lp[dom.lp[n1][1]][1]][6];
                                for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxminb+(i-dom.nxmax+1)/2; jj=dom.nymax-nghosts/2+(j+2)/2;
                                      Bi[i][j][k][nbs]=Bi[ii][jj][k][n3];
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
                                    Bi[i][dom.ny-1-j][k][nbs]=Bi[ii][jj][k][n4];
                                  }}}
                            } else {
                              n4=dom.lp[n4][2];
			      if (dom.leafs[n4]) { // Means more refined corner
                                for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      i1=dom.nxminb+2*(i-dom.nxmax)-1; j1=dom.nymin+nghosts*2-2*(j+1);
				      i2= k==0 ? i1 : i1+1;
                                      j2= k==0 ? j1+1 : j1;
                                      Bi[i][dom.ny-1-j][k][nbs]=0.5*sum_irjr(Bi,n4,k,i1,i2,j1,j2);
                                    }}}
                              } else { // coarser corner
			        n4=dom.lp[dom.lp[dom.lp[n2][1]][1]][7];
                                for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxminb+(i-dom.nxmax+1)/2; jj=dom.nymin+nghosts/2-(j+2)/2;
                                      Bi[i][dom.ny-1-j][k][nbs]=Bi[ii][jj][k][n4];
                                    }}}
			      }
			    }
                          }
			  // Neighbor's left
			  // to finer level (0th order interpolation)
			  for (int i=1;i<=nghosts;i++) {
			    for (int j=0;j<dom.ny;j++) {
			      for (int k=0;k<nvar;k++) {
			        Bi[i-1][j][k][n1]=Bi[(i+1)/2+dom.nxmax-nghosts/2][(j+nghosts)/2][k][nbs];
			        Bi[i-1][j][k][n2]=Bi[(i+1)/2+dom.nxmax-nghosts/2][(j-nghosts+2)/2+dom.ny2][k][nbs];
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
			        Bi[i][j][k][nbs]=Bi[i][jj][k][n1];
				jj++;
			      }
			      // Neighbor's top
			      jj=dom.nymin;
			      for (int j=dom.nyp1;j<dom.ny;j++) {
			        Bi[i][j][k][n1]=Bi[i][jj][k][nbs];
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
			      i1=2*i-dom.nxmin;
			      j1=2*(j-nghosts/2)-nghosts/2+dom.nymax;
 			      for (int k=0;k<nvar;k++) {
				i2= k==0 ? i1 : i1+1;
                                j2= k==0 ? j1+1 : j1;
			        Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n1,k,i1,i2,j1,j2);
			        Bi[nx2i][j][k][nbs]=0.5*sum_irjr(Bi,n2,k,i1,i2,j1,j2);
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
                                    Bi[i][j][k][nbs]=Bi[ii][jj][k][n3];
                                  }}}
                            } else {
			      n3=dom.lp[n3][2]+3;
			      if (dom.leafs[n3]) { // Means more refined corner
			        for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      i1=dom.nxmax-nghosts*2+2*i+1; j1=dom.nymax-nghosts*2+2*j+1;
				      i2= k==0 ? i1 : i1+1;
                                      j2= k==0 ? j1+1 : j1;
                                      Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n3,k,i1,i2,j1,j2);
                                    }}}
			      } else { // coarser corner
			        n3=dom.lp[dom.lp[dom.lp[n1][1]][1]][4];
                                for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxmax-nghosts/2+(i+2)/2; jj=dom.nymax-nghosts/2+(j+2)/2;
                                      Bi[i][j][k][nbs]=Bi[ii][jj][k][n3];
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
                                    Bi[dom.nx-1-i][j][k][nbs]=Bi[ii][jj][k][n4];
                                  }}}
                            } else {
                              n4=dom.lp[n4][2]+2;
			      if (dom.leafs[n4]) { // Means more refined corner
                                for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
				      i1=dom.nxmin+nghosts*2-2*(i+1); j1=dom.nymax-nghosts*2+2*j+1;
				      i2= k==0 ? i1 : i1+1;
                                      j2= k==0 ? j1+1 : j1;
                                      Bi[dom.nx-1-i][j][k][nbs]=0.5*sum_irjr(Bi,n4,k,i1,i2,j1,j2);
                                    }}}
                              } else { // coarser corner
			        n4=dom.lp[dom.lp[dom.lp[n2][1]][1]][5];
                                for (int k=0;k<nvar;k++) {
                                  for (int i=0;i<nghosts;i++) {
                                    for (int j=0;j<nghosts;j++) {
                                      ii=dom.nxmin+nghosts/2-(i+2)/2; jj=dom.nymax-nghosts/2+(j+2)/2;
                                      Bi[dom.nx-1-i][j][k][nbs]=Bi[ii][jj][k][n4];
                                    }}} 
			      }
			    }
                          }
			  // Neighbor's top
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts;j++) {
			    for (int i=0;i<dom.nx;i++) {
			      for (int k=0;k<nvar;k++) {
			        Bi[i][j+dom.nymax][k][n1]=Bi[(i+nghosts)/2][(j+2*nghosts-1)/2][k][nbs];
			        Bi[i][j+dom.nymax][k][n2]=Bi[(i-nghosts+2)/2+dom.nx2][(j+2*nghosts-1)/2][k][nbs];
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
			        Bi[i][j][k][nbs]=Bi[i][jj][k][n1];
				jj++;
			      }
			      // Neighbor's bottom
			      jj=dom.nymax-(nghosts-1);
			      for (int j=0;j<nghosts;j++) {
			        Bi[i][j][k][n1]=Bi[i][jj][k][nbs];
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
		              i1=2*i-dom.nxmin;
			      j1=2*(j-dom.nymax+(nghosts/2-1));
			      int nx2i=dom.nx2+i-1-(nghosts-2);
 			      for (int k=0;k<nvar;k++) {
				i2= k==0 ? i1 : i1+1;
                                j2= k==0 ? j1+1 : j1;
			        Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n1,k,i1,i2,j1,j2);
			        Bi[nx2i][j][k][nbs]=0.5*sum_irjr(Bi,n2,k,i1,i2,j1,j2);
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
                                    Bi[i][j][k][nbs]=Bi[ii][jj][k][n3];
                                  }}}
                            } else {
			      n3=dom.lp[n3][2]+1;
			      if (dom.leafs[n3]) { // Means more refined corner
			      for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=dom.nymaxb;j<dom.ny;j++) {
				    i1=dom.nxmax-nghosts*2+2*i+1; j1=dom.nyminb+2*(j-dom.nymax)-1;
				    i2= k==0 ? i1 : i1+1;
                                    j2= k==0 ? j1+1 : j1;
                                    Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n3,k,i1,i2,j1,j2);
                                  }}}
			    } else { // coarser corner
			      n3=dom.lp[dom.lp[dom.lp[n1][1]][1]][4];
                              for (int k=0;k<nvar;k++) {
                                for (int i=0;i<nghosts;i++) {
                                  for (int j=dom.nymaxb;j<dom.ny;j++) {
                                    ii=dom.nxmax-nghosts/2+(i+2)/2; jj=dom.nyminb+(j-dom.nymax+1)/2;
                                    Bi[i][j][k][nbs]=Bi[ii][jj][k][n3];
                                  }}}
			      }
			    }			      
			  }
			  n4=dom.lp[dom.lp[n2][1]][5]; // right corner
			  if (n4!=-1) { // ensure not at grid boundary
                            if (dom.leafs[n4]) { // same res corner
                              for (int k=0;k<nvar;k++) {
                                for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                  for (int j=dom.nymaxb;j<dom.ny;j++) {
                                    ii=dom.nxminb+(i-dom.nxmax); jj=dom.nyminb+(j-dom.nymax);
                                    Bi[i][j][k][nbs]=Bi[ii][jj][k][n4];
                                  }}}
			    } else {
                              n4=dom.lp[n4][2];
			      if (dom.leafs[n4]) { // Means more refined corner
                                for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=dom.nymaxb;j<dom.ny;j++) {
                                      i1=dom.nxminb+2*(i-dom.nxmax)-1; j1=dom.nyminb+2*(j-dom.nymax)-1;
				      i2= k==0 ? i1 : i1+1;
                                      j2= k==0 ? j1+1 : j1;
                                      Bi[i][j][k][nbs]=0.5*sum_irjr(Bi,n4,k,i1,i2,j1,j2);
                                    }}}
                              } else { // coarser corner
			        n4=dom.lp[dom.lp[dom.lp[n2][1]][1]][5];
                                for (int k=0;k<nvar;k++) {
                                  for (int i=dom.nxmaxb;i<dom.nx;i++) {
                                    for (int j=dom.nymaxb;j<dom.ny;j++) {
                                      ii=dom.nxminb+(i-dom.nxmax+1)/2; jj=dom.nyminb+(j-dom.nymax+1)/2;
                                      Bi[i][j][k][nbs]=Bi[ii][jj][k][n4];
                                    }}}
			      }
			    }
                          }
			  // Neighbor's bottom
			  // to finer level (0th order interpolation)
			  for (int j=1;j<=nghosts;j++) {
			    for (int i=0;i<dom.nx;i++) {
			      for (int k=0;k<nvar;k++) {
			        Bi[i][j-1][k][n1]=Bi[(i+nghosts)/2][(j+1)/2+dom.nymax-nghosts/2][k][nbs];
			        Bi[i][j-1][k][n2]=Bi[(i-nghosts+2)/2+dom.nx2][(j+1)/2+dom.nymax-nghosts/2][k][nbs];
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
	extboundary(dom,Bi,nvar);
}
