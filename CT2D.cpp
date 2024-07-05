#include<exception>
#include<iomanip>
#include "defs.hpp"
#include "timeIntegral.hpp"

void setBCs(int nx,int ny,int nvar,int nb,real**** Q);
void boundary(meshblock &dom,real**** Q,int nvar);
void boundaryfc(meshblock &dom,real**** Bi,int nvar);
void boundaryec(meshblock &dom,real**** EMFz,int nvar);

void CT2D(meshblock &dom,real**** q,real**** Bi,int step,real**** Binew) {
	real res;

	// Main loop
	if (CTtype==1) {
	for (int nb=0;nb<dom.lastActive;nb++) {
          if (dom.leafs[nb]) {
	    #pragma omp parallel for collapse(2) default(none) shared(dom,nb)
	    for (int i=1; i<dom.nx; i++) {
	      for (int j=1; j<dom.ny; j++) {
                dom.EMF[i][j][0][nb]=0.25*((dom.gg[i-1][j][5][nb]+dom.gg[i][j][5][nb])
				-(dom.ff[i][j-1][6][nb]+dom.ff[i][j][6][nb]));
	      }
	    }	
	  }
	}

	} else if (CTtype==2) { // Upwind type 
	real p[dom.nx][dom.ny];
	real C[dom.nx][dom.ny];

	 for (int nb=0;nb<dom.lastActive;nb++) {
          if (dom.leafs[nb]) {

	    #pragma omp parallel for collapse(2) default(none) shared(dom,p,C,q,nb)	  
	    for (int i=0; i<dom.nx; i++) {
	      for (int j=0; j<dom.ny; j++) {
	        p[i][j]=(dom.gamma-1)*(q[i][j][4][nb]-0.5*(MAG(q[i][j][1][nb],q[i][j][2][nb],q[i][j][3][nb])/q[i][j][0][nb]+MAG(q[i][j][5][nb],q[i][j][6][nb],q[i][j][7][nb])));
	        C[i][j]=sqrt((dom.gamma*p[i][j]+MAG(q[i][j][5][nb],q[i][j][6][nb],q[i][j][7][nb]))/q[i][j][0][nb]);
	      }
	    }

	    #pragma omp parallel for collapse(2) default(none) shared(dom,nb,p,C,q)
            for (int i=1; i<dom.nx; i++) {
              for (int j=1; j<dom.ny; j++) {
		real beta=0.5;
		real delta=0.1;

		real dPx=0.5*abs(p[i][j]+p[i][j-1]-p[i-1][j]-p[i-1][j-1]);
		real dPy=0.5*abs(p[i][j]+p[i-1][j]-p[i][j-1]-p[i-1][j-1]);
		real divV=0.5*(q[i][j][1][nb]/q[i][j][0][nb]+q[i][j-1][1][nb]/q[i][j-1][0][nb]-q[i-1][j][1][nb]/q[i-1][j][0][nb]-q[i-1][j-1][1][nb]/q[i-1][j-1][0][nb]) + 
		     0.5*(q[i][j][2][nb]/q[i][j][0][nb]+q[i-1][j][2][nb]/q[i-1][j][2][nb]-q[i][j-1][2][nb]/q[i][j-1][2][nb]-q[i-1][j-1][2][nb]/q[i-1][j-1][0][nb]);

		// Switch 1: pick out situations in vicinity of strong magnetosonic shock
		real tempP=MIN(p[i-1][j-1],p[i][j-1]);
		real tempP1=MIN(p[i-1][j],p[i][j]);
		tempP=MIN(tempP,tempP1);		
		bool sw1= (dPx+dPy)>beta*tempP ? true:false;

		// Switch 2: pick out strongly compressive motions
		real tempC=MIN(C[i-1][j-1],C[i][j-1]);
		real tempC1=MIN(C[i-1][j],C[i][j]);
		tempC=MIN(tempC,tempC1);
		bool sw2= (-delta*tempC)>divV ? true:false;

		real psi=0.5;
		if (sw1 and sw2) psi=dPx/(dPx+dPy);
		if (psi<0 or psi>1) {printf("Incorrect psi produced... \n"); throw exception();}
		dom.EMF[i][j][0][nb]=0.5*(1.-psi)*(dom.gg[i-1][j][5][nb]+dom.gg[i][j][5][nb])
				-0.5*psi*(dom.ff[i][j-1][6][nb]+dom.ff[i][j][6][nb]);
	      }
	    }	
	  }
	}
	}
	boundaryec(dom,dom.EMF,1);

	for (int nb=0;nb<dom.lastActive;nb++) {
	  if (dom.leafs[nb]) {
	    #pragma omp parallel for collapse(2) default(none) shared(dom,step,Binew,Bi,nb) private(res)
	    for (int i=0; i<dom.nx; i++) {
	      for (int j=0; j<dom.ny; j++) {
	        int ip = i<dom.nx-1 ? i+1 : i;
	        int jp = j<dom.ny-1 ? j+1 : j;
	        res=(1./dom.dy[dom.lp[nb][0]])*(dom.EMF[i][jp][0][nb]-dom.EMF[i][j][0][nb]);
	        Binew[i][j][0][nb]=RK2ndv1(dom.dt,Binew[i][j][0][nb],Bi[i][j][0][nb],res,step);	      
	        res=-(1./dom.dx[dom.lp[nb][0]])*(dom.EMF[ip][j][0][nb]-dom.EMF[i][j][0][nb]);
	        Binew[i][j][1][nb]=RK2ndv1(dom.dt,Binew[i][j][1][nb],Bi[i][j][1][nb],res,step);
	      }	  
	    }
	  }
	}	
	boundaryfc(dom,Binew,2);

	real divB;
	for (int nb=0;nb<dom.lastActive;nb++) {
          if (dom.leafs[nb]) {
	    #pragma omp parallel for collapse(2) default(none) shared(dom,q,Binew,nb) private(divB)
	    for (int i=0; i<dom.nx; i++) {
	      for (int j=0; j<dom.ny; j++) {
	        int ip = i<dom.nx-1 ? i+1 : i;
	        int ii = ip-1;
	        int jp = j<dom.ny-1 ? j+1 : j;
	        int jj = jp-1;
	        // Check if div(B) is at machine precision
	        divB=abs((Binew[ip][j][0][nb]-Binew[ii][j][0][nb])/dom.dx[dom.lp[nb][0]] 
			    +(Binew[i][jp][1][nb]-Binew[i][jj][1][nb])/dom.dy[dom.lp[nb][0]]);
	        if (false and divB>eps and (i>=dom.nxmin and i<=dom.nxmax) and (j>=dom.nymin and j<=dom.nymax)) {
		  printf("CT error as div(B)=%f at nb=%d, i=%d, j=%d \n",divB,nb,i,j);
	          throw exception();
	        }
	        // Update magnetic field with face centered values	    
	        real Bfcx=0.5*(Binew[ip][j][0][nb]+Binew[ii][j][0][nb]);
	        real Bfcy=0.5*(Binew[i][jp][1][nb]+Binew[i][jj][1][nb]);
	        // Apply energy correction for improve positivity
	        q[i][j][4][nb]=q[i][j][4][nb]+0.5*((Bfcx*Bfcx+Bfcy*Bfcy)-(pow(q[i][j][5][nb],2)+pow(q[i][j][6][nb],2)));
	        // Update cell-centered magnetic field with face-centered magnetic field
	        q[i][j][5][nb]=Bfcx;
                q[i][j][6][nb]=Bfcy;
	      }
	    }
	  }
	}
}
