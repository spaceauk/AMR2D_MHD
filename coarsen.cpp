// Coarsen to parent level: son1 to son4 --> dad
#include<exception>
#include "defs.hpp"

real sum_irjr(real**** U, int nb, int k, int i1, int i2, int j1, int j2);

void coarsen(meshblock &dom, real**** Q, int nvar, int dad, int son1, int son2,
	     int son3, int son4) {
	int i1, i2, j1, j2;
	int nx2i, ny2j;
	// The block itself
	for (int i=dom.nxmin;i<=dom.nx2;i++) {
	  for (int j=dom.nymin;j<=dom.ny2;j++) {
	    i1=2*i-dom.nxmin;
	    i2=i1+1;
	    j1=2*j-dom.nymin;
	    j2=j1+1;
	    nx2i=dom.nx2+i-1-(nghosts-2);
	    ny2j=dom.ny2+j-1-(nghosts-2);
	    for (int k=0;k<nvar;k++) {
	      Q[i][j][k][dad]=0.25*sum_irjr(Q,son1,k,i1,i2,j1,j2);
	      Q[nx2i][j][k][dad]=0.25*sum_irjr(Q,son2,k,i1,i2,j1,j2);
	      Q[i][ny2j][k][dad]=0.25*sum_irjr(Q,son3,k,i1,i2,j1,j2);
	      Q[nx2i][ny2j][k][dad]=0.25*sum_irjr(Q,son4,k,i1,i2,j1,j2);
	    }
	  }
	}
	cout<<" "<<dad<<" ";

	// The ghost cells of the block
	if (nghosts>=2) {
	  // Left
	  for (int k=0;k<nvar;k++) {
	    for (int i=0;i<nghosts/2;i++) {
	      for (int j=0;j<nghosts/2;j++) {
		i2=dom.nxminb-2*i;
		i1=i2-1;
		j2=dom.nyminb-2*j;
		j1=j2-1;
	        Q[dom.nxminb-i][dom.nyminb-j][k][dad]=0.25*sum_irjr(Q,son1,k,i1,i2,j1,j2);
		j1=dom.nymaxb+2*j;
		j2=j1+1;
	        Q[dom.nxminb-i][dom.nymaxb+j][k][dad]=0.25*sum_irjr(Q,son3,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  for (int i=0;i<=dom.nxminb;i++) {
	    for (int j=dom.nymin;j<=dom.ny2;j++) {
	      j1=2*j-dom.nymin;
	      j2=j1+1;
	      ny2j=j+dom.ny2-1-(nghosts-2);
	      i2=dom.nxminb-2*i; i1=i2-1;
	      if (i1<0) {i1=0; i2=1;}
	      for (int k=0;k<nvar;k++) {
	        Q[dom.nxminb-i][j][k][dad]=0.25*sum_irjr(Q,son1,k,i1,i2,j1,j2);
		Q[dom.nxminb-i][ny2j][k][dad]=0.25*sum_irjr(Q,son3,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  // Right
	  for (int k=0;k<nvar;k++) {
	    for (int i=0;i<nghosts/2;i++) {
	      for (int j=0;j<nghosts/2;j++) {
		i1=dom.nxmaxb+2*i;
		i2=i1+1;
		j2=dom.nyminb-2*j;
		j1=j2-1;
	        Q[dom.nxmaxb+i][dom.nyminb-j][k][dad]=0.25*sum_irjr(Q,son2,k,i1,i2,j1,j2);
		j1=dom.nymaxb+2*j;
		j2=j1+1;
	        Q[dom.nxmaxb+i][dom.nymaxb+j][k][dad]=0.25*sum_irjr(Q,son4,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  for (int i=0;i<nghosts;i++) {
	    for (int j=dom.nymin;j<=dom.ny2;j++) {
	      j1=2*j-dom.nymin;
	      j2=j1+1;
	      ny2j=j+dom.ny2-1-(nghosts-2);
	      i1=dom.nxmaxb+2*i;
	      i2=i1+1;
	      if (i2>dom.nx-1) {i1=dom.nx-2; i2=i1+1;}
	      for (int k=0;k<nvar;k++) {
	        Q[dom.nxmaxb+i][j][k][dad]=0.25*sum_irjr(Q,son2,k,i1,i2,j1,j2);
		Q[dom.nxmaxb+i][ny2j][k][dad]=0.25*sum_irjr(Q,son4,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  // Bottom
	  for (int k=0;k<nvar;k++) {
	    for (int j=0;j<nghosts/2;j++) {
	      for (int i=0;i<nghosts/2;i++) {
		j2=dom.nyminb-2*j;
		j1=j2-1;
		i2=dom.nxminb-2*i;
		i1=i2-1;
	        Q[dom.nxminb-i][dom.nyminb-j][k][dad]=0.25*sum_irjr(Q,son1,k,i1,i2,j1,j2);
		i1=dom.nxmaxb+2*i;
		i2=i1+1;
	        Q[dom.nxmaxb+i][dom.nyminb-j][k][dad]=0.25*sum_irjr(Q,son2,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  for (int i=dom.nxmin;i<=dom.nx2;i++) {
	    for (int j=0;j<=dom.nyminb;j++) {
	      i1=2*i-dom.nxmin;
	      i2=i1+1;
	      nx2i=i+dom.nx2-1-(nghosts-2);
	      j2=dom.nyminb-2*j;
	      j1=j2-1;
	      if (j1<0) {j1=0;j2=1;}
	      for (int k=0;k<nvar;k++) {
	        Q[i][dom.nyminb-j][k][dad]=0.25*sum_irjr(Q,son1,k,i1,i2,j1,j2);
		Q[nx2i][dom.nyminb-j][k][dad]=0.25*sum_irjr(Q,son2,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  // Top
	  for (int k=0;k<nvar;k++) {
	    for (int j=0;j<nghosts/2;j++) {
	      for (int i=0;i<nghosts/2;i++) {
		j1=dom.nymaxb+2*j;
		j2=j1+1;
		i2=dom.nxminb-2*i;
		i1=i2-1;
	        Q[dom.nxminb-i][dom.nymaxb+j][k][dad]=0.25*sum_irjr(Q,son3,k,i1,i2,j1,j2);
		i1=dom.nxmaxb+2*i;
		i2=i1+1;
	        Q[dom.nxmaxb+i][dom.nymaxb+j][k][dad]=0.25*sum_irjr(Q,son4,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  for (int i=dom.nxmin;i<=dom.nx2;i++) {
	    for (int j=0;j<nghosts;j++) {
	      i1=2*i-dom.nxmin;
	      i2=i1+1;
	      nx2i=i+dom.nx2-1-(nghosts-2);
	      j1=dom.nymaxb+2*j;
	      j2=j1+1;
	      if (j2>dom.ny-1) {j1=dom.ny-2; j2=j1+1;}
	      for (int k=0;k<nvar;k++) {
	        Q[i][dom.nymaxb+j][k][dad]=0.25*sum_irjr(Q,son3,k,i1,i2,j1,j2);
		Q[nx2i][dom.nymaxb+j][k][dad]=0.25*sum_irjr(Q,son4,k,i1,i2,j1,j2);
	      }
	    }
	  }

	  // Convert to primitive variables
	  dom.U2W("coarsen",dad);

	} else {
	  cout<<"Error at coarsening as wrong ghost cells chosen!"<<endl;
	  throw exception();
	}
}
