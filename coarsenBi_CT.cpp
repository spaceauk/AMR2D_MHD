// Coarsen to parent level: son1 to son4 --> dad
#include<exception>
#include "defs.hpp"

real sum_irjr(real**** U, int nb, int k, int i1, int i2, int j1, int j2);

void coarsenBi_CT(meshblock &dom, real**** Bi, int nvar, int dad, int son1, int son2,
	     int son3, int son4) {
	int i1, i2, j1, j2;
	int nx2i, ny2j;
	// The block itself
	#pragma omp parallel for collapse(2) default(none) \
	   shared(dom,Bi,nvar,dad,son1,son2,son3,son4) private(i1,i2,j1,j2,nx2i,ny2j)
	for (int i=dom.nxmin;i<=dom.nx2;i++) {
	  for (int j=dom.nymin;j<=dom.ny2;j++) {
	    i1=2*i-dom.nxmin;
	    j1=2*j-dom.nymin;
	    nx2i=dom.nx2+i-1-(nghosts-2);
	    ny2j=dom.ny2+j-1-(nghosts-2);
	    for (int k=0;k<nvar;k++) {
	      i2= k==0 ? i1 : i1+1;
	      j2= k==0 ? j1+1 : j1;
	      Bi[i][j][k][dad]=0.5*sum_irjr(Bi,son1,k,i1,i2,j1,j2);
	      Bi[nx2i][j][k][dad]=0.5*sum_irjr(Bi,son2,k,i1,i2,j1,j2);
	      Bi[i][ny2j][k][dad]=0.5*sum_irjr(Bi,son3,k,i1,i2,j1,j2);
	      Bi[nx2i][ny2j][k][dad]=0.5*sum_irjr(Bi,son4,k,i1,i2,j1,j2);
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
		i1=dom.nxminb-2*i-1;
		i2= k==0 ? i1 : i1+1;
		j1=dom.nyminb-2*j-1;
                j2= k==0 ? j1+1 : j1;
	        Bi[dom.nxminb-i][dom.nyminb-j][k][dad]=0.5*sum_irjr(Bi,son1,k,i1,i2,j1,j2);
		j1=dom.nymaxb+2*j;
		j2= k==0 ? j1+1 : j1;;
	        Bi[dom.nxminb-i][dom.nymaxb+j][k][dad]=0.5*sum_irjr(Bi,son3,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  for (int i=0;i<=dom.nxminb;i++) {
	    for (int j=dom.nymin;j<=dom.ny2;j++) {
	      j1=2*j-dom.nymin;
	      ny2j=j+dom.ny2-1-(nghosts-2);
	      i1=dom.nxminb-2*i-1;
	      if (i1<0) {i1=0;}
	      for (int k=0;k<nvar;k++) {
		i2= k==0 ? i1 : i1+1;
                j2= k==0 ? j1+1 : j1;
	        Bi[dom.nxminb-i][j][k][dad]=0.5*sum_irjr(Bi,son1,k,i1,i2,j1,j2);
		Bi[dom.nxminb-i][ny2j][k][dad]=0.5*sum_irjr(Bi,son3,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  // Right
	  for (int k=0;k<nvar;k++) {
	    for (int i=0;i<nghosts/2;i++) {
	      for (int j=0;j<nghosts/2;j++) {
		i1=dom.nxmaxb+2*i;
		i2= k==0 ? i1 : i1+1;
		j1=dom.nyminb-2*j-1;
                j2= k==0 ? j1+1 : j1;
	        Bi[dom.nxmaxb+i][dom.nyminb-j][k][dad]=0.5*sum_irjr(Bi,son2,k,i1,i2,j1,j2);
		j1=dom.nymaxb+2*j;
		j2= k==0 ? j1+1 : j1;
	        Bi[dom.nxmaxb+i][dom.nymaxb+j][k][dad]=0.5*sum_irjr(Bi,son4,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  for (int i=0;i<nghosts;i++) {
	    for (int j=dom.nymin;j<=dom.ny2;j++) {
	      j1=2*j-dom.nymin;
	      ny2j=j+dom.ny2-1-(nghosts-2);
	      i1=dom.nxmaxb+2*i;
	      if (i1+1>dom.nx-1) {i1=dom.nx-2;}
	      for (int k=0;k<nvar;k++) {
		i2= k==0 ? i1 : i1+1;
                j2= k==0 ? j1+1 : j1;
	        Bi[dom.nxmaxb+i][j][k][dad]=0.5*sum_irjr(Bi,son2,k,i1,i2,j1,j2);
		Bi[dom.nxmaxb+i][ny2j][k][dad]=0.5*sum_irjr(Bi,son4,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  // Bottom
	  for (int k=0;k<nvar;k++) {
	    for (int j=0;j<nghosts/2;j++) {
	      for (int i=0;i<nghosts/2;i++) {
		j1=dom.nyminb-2*j-1;
                j2= k==0 ? j1+1 : j1;
		i1=dom.nxminb-2*i-1;
		i2= k==0 ? i1 : i1+1;
	        Bi[dom.nxminb-i][dom.nyminb-j][k][dad]=0.5*sum_irjr(Bi,son1,k,i1,i2,j1,j2);
		i1=dom.nxmaxb+2*i;
		i2= k==0 ? i1 : i1+1;
	        Bi[dom.nxmaxb+i][dom.nyminb-j][k][dad]=0.5*sum_irjr(Bi,son2,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  for (int i=dom.nxmin;i<=dom.nx2;i++) {
	    for (int j=0;j<=dom.nyminb;j++) {
	      i1=2*i-dom.nxmin;
	      nx2i=i+dom.nx2-1-(nghosts-2);
	      j1=dom.nyminb-2*j-1;
	      if (j1<0) {j1=0;}
	      for (int k=0;k<nvar;k++) {
		i2= k==0 ? i1 : i1+1;
                j2= k==0 ? j1+1 : j1;
	        Bi[i][dom.nyminb-j][k][dad]=0.5*sum_irjr(Bi,son1,k,i1,i2,j1,j2);
		Bi[nx2i][dom.nyminb-j][k][dad]=0.5*sum_irjr(Bi,son2,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  // Top
	  for (int k=0;k<nvar;k++) {
	    for (int j=0;j<nghosts/2;j++) {
	      for (int i=0;i<nghosts/2;i++) {
		j1=dom.nymaxb+2*j;
                j2= k==0 ? j1+1 : j1;
		i1=dom.nxminb-2*i-1;
		i2= k==0 ? i1 : i1+1;
	        Bi[dom.nxminb-i][dom.nymaxb+j][k][dad]=0.5*sum_irjr(Bi,son3,k,i1,i2,j1,j2);
		i1=dom.nxmaxb+2*i;
		i2= k==0 ? i1 : i1+1;
	        Bi[dom.nxmaxb+i][dom.nymaxb+j][k][dad]=0.5*sum_irjr(Bi,son4,k,i1,i2,j1,j2);
	      }
	    }
	  }
	  for (int i=dom.nxmin;i<=dom.nx2;i++) {
	    for (int j=0;j<nghosts;j++) {
	      i1=2*i-dom.nxmin;
	      nx2i=i+dom.nx2-1-(nghosts-2);
	      j1=dom.nymaxb+2*j;
	      if (j1+1>dom.ny-1) {j1=dom.ny-2;}
	      for (int k=0;k<nvar;k++) {
		i2= k==0 ? i1 : i1+1;
                j2= k==0 ? j1+1 : j1;
	        Bi[i][dom.nymaxb+j][k][dad]=0.5*sum_irjr(Bi,son3,k,i1,i2,j1,j2);
		Bi[nx2i][dom.nymaxb+j][k][dad]=0.5*sum_irjr(Bi,son4,k,i1,i2,j1,j2);
	      }
	    }
	  }

	} else {
	  cout<<"Error at coarsening as wrong ghost cells chosen!"<<endl;
	  throw exception();
	}
}
