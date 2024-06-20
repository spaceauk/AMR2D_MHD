#include<exception>
#include<iomanip>
#include "defs.hpp"
#include "timeIntegral.hpp"

void setBCs(int nx,int ny,int nvar,int nb,real**** Q);

void CT2D(meshblock &dom,real**** q,real**** Bi,int nb,int step,real**** Binew) {
	double EMFz[dom.nx][dom.ny];
	real res;

	// Main loop
	for (int i=1; i<dom.nx; i++) {
	  for (int j=1; j<dom.ny; j++) {
            EMFz[i][j]=0.25*((dom.gg[i-1][j][5][nb]+dom.gg[i][j][5][nb])-(dom.ff[i][j-1][6][nb]+dom.ff[i][j][6][nb]));
	  }
	}	
        for (int i=0;i<dom.nx;i++) {EMFz[i][0]=EMFz[i][1];} // Need this otherwise will have wrong values at count==2
        for (int j=0;j<dom.ny;j++) {EMFz[0][j]=EMFz[1][j];}

	for (int i=0; i<dom.nx; i++) {
	  for (int j=0; j<dom.ny; j++) {
	    int ip = i<dom.nx-1 ? i+1 : i;
	    int jp = j<dom.ny-1 ? j+1 : j;
	    // If div free, will remain div free as (assumed dx=dy here) 
            // -(EMFz[i+1][j+1]-EMFz[i+1][j])+(EMFz[i][j+1]-EMFz[i][j]) +
            //  (EMFz[i+1][j+1]-EMFz[i][j+1])-(EMFz[i+1][j]-EMFz[i][j]) = 0
	    res=(1./dom.dy[dom.lp[nb][0]])*(EMFz[i][jp]-EMFz[i][j]);
	    Binew[i][j][0][nb]=RK2ndv1(dom.dt,Binew[i][j][0][nb],Bi[i][j][0][nb],res,step);	      
	    res=-(1./dom.dx[dom.lp[nb][0]])*(EMFz[ip][j]-EMFz[i][j]);
	    Binew[i][j][1][nb]=RK2ndv1(dom.dt,Binew[i][j][1][nb],Bi[i][j][1][nb],res,step);
	  }
	}
	//setBCs(dom.nx,dom.ny,2,nb,Binew);
	real divB;
	for (int i=0; i<dom.nx; i++) {
	  for (int j=0; j<dom.ny; j++) {
	    int ip = i<dom.nx-1 ? i+1 : i;
	    int ii = ip-1;
	    int jp = j<dom.ny-1 ? j+1 : j;
	    int jj = jp-1;
	    // Check if div(B) is at machine precision
	    divB=abs((Binew[ip][j][0][nb]-Binew[ii][j][0][nb])/dom.dx[dom.lp[nb][0]] 
			    +(Binew[i][jp][1][nb]-Binew[i][jj][1][nb])/dom.dy[dom.lp[nb][0]]);
	    if (divB>1.E-8 and (i>=dom.nxmin and i<=dom.nxmax) and (j>=dom.nymin and j<=dom.nymax)) {
	      cout<<"CT error as div(B)="<<divB<<" at nb="<<nb<<", i="<<i<<", j="<<j<<endl;
	      cout<<std::setprecision(Nprec)<<"  where dBx="<<Binew[ip][j][0][nb]-Binew[ii][j][0][nb]<<", Bx1="<<Binew[ii][j][0][nb]<<" & Bx2="
		  <<Binew[ip][j][0][nb]<<", dx="<<dom.dx[dom.lp[nb][0]]<<endl;
	      cout<<std::setprecision(Nprec)<<"        dBy="<<Binew[i][jp][1][nb]-Binew[i][jj][1][nb]<<", By1="<<Binew[i][jj][1][nb]<<" & By2="
		  <<Binew[i][jp][1][nb]<<", dy="<<dom.dy[dom.lp[nb][0]]<<endl;
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
