#include <fstream>
#include <string>
#include "defs.hpp"
#include <exception>

void meshblock::setSize(int nx, int ny, int nvar, real lenx, real leny) {
	this->nx=nx;  
	this->ny=ny; 
	this->nvar=nvar;
	if (MAG_field and nvar<8) {
		cout<<"Insufficient number of variables used due to magnetic field!"<<endl;
		throw exception();
	} else if (!MAG_field and nvar>5) {
		cout<<"Reduce number of variables to save space."<<endl;
		nvar=5;
	}

	// Use new to allocate memory to array for conservative variable
	U = new real***[nx]; // [r,ru,rv,rw,E,Bx,By,Bz]
	Us = new real***[nx];
	W = new real***[nx]; // [r,u,v,w,p,Bx,By,Bz]
	dwdx = new real**[nx];			   
	dwdy = new real**[nx];			    
	wxL = new real**[nx];			    
	wxR = new real**[nx];			    
	wyL = new real**[nx];			    
	wyR = new real**[nx];			    
	ff = new real***[nx];
	gg = new real***[nx];	
	if  (CT_mtd) {
		Bi = new double***[nx]; 
		Bis = new double***[nx];
		EMF = new double***[nx];
	}

	for (int i=0; i<nx; i++) {
		U[i] = new real**[ny];
		Us[i] = new real**[ny];
		W[i] = new real**[ny];
		dwdx[i] = new real*[ny];
		dwdy[i] = new real*[ny];
		wxL[i] = new real*[ny];
		wxR[i] = new real*[ny];
		wyL[i] = new real*[ny];
		wyR[i] = new real*[ny];
		ff[i] = new real**[ny];
		gg[i] = new real**[ny];
		if (CT_mtd) {
			Bi[i] = new double**[ny]; 
			Bis[i] = new double**[ny];
			EMF[i] = new double**[ny];
		}

		for (int j =0; j<ny; j++) {
			U[i][j] = new real*[nvar];
			Us[i][j] = new real*[nvar];
			W[i][j] = new real*[nvar];
			dwdx[i][j] = new real[nvar];
			dwdy[i][j] = new real[nvar];
			wxL[i][j] = new real[nvar];
			wxR[i][j] = new real[nvar];
			wyL[i][j] = new real[nvar];
			wyR[i][j] = new real[nvar];
			ff[i][j] = new real*[nvar];
			gg[i][j] = new real*[nvar];
			if (CT_mtd) {
				Bi[i][j] = new double*[2];  // For 2D only 2 magnetic field direction with only 1 electric field
				Bis[i][j] = new double*[2]; // perpendicular to them.
				EMF[i][j] = new double*[1];
			}
			for (int k=0;k<nvar;k++) {
				U[i][j][k] = new real[maxblocks];
				Us[i][j][k] = new real[maxblocks];
                                W[i][j][k] = new real[maxblocks];
				ff[i][j][k] = new real[maxblocks];
				gg[i][j][k] = new real[maxblocks];
				if (CT_mtd and k<2) {
					Bi[i][j][k] = new double[maxblocks]; 
					Bis[i][j][k] = new double[maxblocks];
					if (k<1) EMF[i][j][k] = new double[maxblocks];
				}
			}
		}
	}

	// Calculate quantities (note boundary cells)
	dx = new real[maxlevs];
        dy = new real[maxlevs];
	// Assume that root blocks are distributed evenly along x & y
	dx[0] = lenx/(sqrt(nbroots)*(nx-2*nghosts-1));
	dy[0] = leny/(sqrt(nbroots)*(ny-2*nghosts-1));
	cout<<"dx[0]="<<dx[0]<<", dy[0]="<<dy[0]<<endl;
	for (int l=1; l<maxlevs; l++) {
                dx[l]=dx[l-1]/2.;
                dy[l]=dy[l-1]/2.;		
        }

	// Initializing AMR-related stuffs
        lp = new int*[maxblocks];
        icoord = new int*[maxblocks];
        leafs = new bool[maxblocks];
        iref = new bool[maxblocks];
        icoarse = new bool[maxblocks];

        for (int nb=0; nb<maxblocks; nb++) {
                lp[nb] = new int[8];
                icoord[nb] = new int[2];
        }

}

void meshblock::setType(string IC, string limiter, string fluxMth) {
	this->IC=IC;
	this->limiter=limiter;
	this->fluxMth=fluxMth;
}

// Conversion between U & W
void meshblock::U2W(string loc,int nb) {
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			W[i][j][0][nb]=U[i][j][0][nb];
			W[i][j][1][nb]=U[i][j][1][nb]/U[i][j][0][nb];
			W[i][j][2][nb]=U[i][j][2][nb]/U[i][j][0][nb];
			W[i][j][3][nb]=U[i][j][3][nb]/U[i][j][0][nb];
			W[i][j][4][nb]=(gamma-1)*(U[i][j][4][nb]-0.5*W[i][j][0][nb]*
				   MAG(W[i][j][1][nb],W[i][j][2][nb],W[i][j][3][nb]));
			if (MAG_field) {
				W[i][j][5][nb]=U[i][j][5][nb];
				W[i][j][6][nb]=U[i][j][6][nb];
				W[i][j][7][nb]=U[i][j][7][nb];
				W[i][j][4][nb]=W[i][j][4][nb]-(gamma-1)*(0.5*MAG(U[i][j][5][nb],U[i][j][6][nb],U[i][j][7][nb]));
			}
			// Check for negative pressure
			if (W[i][j][4][nb]<=0 and i>nxminb and i<nxmaxb and j>nyminb and j<nymaxb) {
				cout<<"Error as zero/negative pressure! P="<<W[i][j][4][nb]<<" at i="<<i<<" & j="<<j<<". \n";
				cout<<" - @ inner boundaries = "<<innerbounds[nb][2]<<endl;
				cout<<" - Bx="<<W[i][j][5][nb]<<", By="<<W[i][j][6][nb]<<", Bz="<<W[i][j][7][nb]<<endl;
				cout<<"U2W in function: "<<loc<<" at loop count="<<count<<endl;
				throw exception();
			}
		}
	}
}

void meshblock::Us2W(string loc,int nb) {
        for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                        W[i][j][0][nb]=Us[i][j][0][nb];
                        W[i][j][1][nb]=Us[i][j][1][nb]/Us[i][j][0][nb];
                        W[i][j][2][nb]=Us[i][j][2][nb]/Us[i][j][0][nb];
                        W[i][j][3][nb]=Us[i][j][3][nb]/Us[i][j][0][nb];
                        W[i][j][4][nb]=(gamma-1)*(Us[i][j][4][nb]-0.5*W[i][j][0][nb]*
                                   MAG(W[i][j][1][nb],W[i][j][2][nb],W[i][j][3][nb]));
                        if (MAG_field) {
                                W[i][j][5][nb]=Us[i][j][5][nb];
                                W[i][j][6][nb]=Us[i][j][6][nb];
                                W[i][j][7][nb]=Us[i][j][7][nb];
                                W[i][j][4][nb]=W[i][j][4][nb]-(gamma-1)*(0.5*MAG(Us[i][j][5][nb],Us[i][j][6][nb],Us[i][j][7][nb]));
                        }
                        // Check for negative pressure
                        if (W[i][j][4][nb]<=0 and i>nxminb and i<nxmaxb and j>nyminb and j<nymaxb) {
                                cout<<"Error as zero/negative pressure! P="<<W[i][j][4]<<" at i="<<i<<" & j="<<j<<". \n";
                                cout<<"Us2W in function: "<<loc<<" at loop count="<<count<<endl;
				cout<<" - Center: Bx="<<W[i][j][5][nb]<<", By="<<W[i][j][6][nb]<<", Bz="<<W[i][j][7][nb]<<endl;
				if (CT_mtd) cout<<" - Edge  : Bx="<<Bi[i][j][0][nb]<<", By="<<Bi[i][j][1][nb]<<endl;
                                cout<<"Inner domain:"<<nxmin<<", "<<nxmax<<", "<<nymin<<", "<<nymax<<endl;
                                throw exception();
                        }
			if (isnan(W[i][j][1][nb]) or isnan(W[i][j][2][nb]) or isnan(W[i][j][3][nb])) {
				cout<<"NaN values encountered in Us2W! nb="<<nb<<", i="<<i<<", j="<<j<<
					", Us="<<Us[i][j][0][nb]<<endl;
				throw exception();
				
			}
                }
        }
}

void meshblock::W2U(int nb) {
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			U[i][j][0][nb]=W[i][j][0][nb];
			U[i][j][1][nb]=W[i][j][1][nb]*W[i][j][0][nb];
			U[i][j][2][nb]=W[i][j][2][nb]*W[i][j][0][nb];
			U[i][j][3][nb]=W[i][j][3][nb]*W[i][j][0][nb];
			U[i][j][4][nb]=(W[i][j][4][nb]/(gamma-1))+0.5*W[i][j][0][nb]*
				   MAG(W[i][j][1][nb],W[i][j][2][nb],W[i][j][3][nb]);
			if (MAG_field) {
				U[i][j][5][nb]=W[i][j][5][nb];
				U[i][j][6][nb]=W[i][j][6][nb];
				U[i][j][7][nb]=W[i][j][7][nb];
				U[i][j][4][nb]=U[i][j][4][nb]+0.5*MAG(W[i][j][5][nb],W[i][j][6][nb],W[i][j][7][nb]);
			}
		}
	}
}


// Calculate wave speed
void meshblock::wavespeed(int i, int j,int nb) {
	speed=sqrt(MAG(W[i][j][1][nb],W[i][j][2][nb],W[i][j][3][nb]));	
	c=sqrt(gamma*W[i][j][4][nb]/W[i][j][0][nb]);
	if (MAG_field) {
		ca=sqrt((MAG(W[i][j][5][nb],W[i][j][6][nb],W[i][j][7][nb]))/W[i][j][0][nb]);
		cax=sqrt((SQR(W[i][j][5][nb]))/W[i][j][0][nb]);
		cay=sqrt((SQR(W[i][j][6][nb]))/W[i][j][0][nb]);
		real intm=SQR(ca)+SQR(c);		
		//cfx=sqrt(0.5*(SQR(c)+SQR(ca))+0.5*sqrt(SQR(intm)
		//		-4*(SQR(c)*SQR(cax))));		
		//cfy=sqrt(0.5*(SQR(c)+SQR(ca))+0.5*sqrt(SQR(intm)
		//		-4*(SQR(c)*SQR(cay))));	
		// Rewrite as +ve definite form
		real intm1=SQR(ca)-SQR(c);
		cfx=sqrt(0.5*(intm+sqrt(intm1*intm1+4.*SQR(ca)*
			(SQR(ca)-SQR(cax)))));
		cfy=sqrt(0.5*(intm+sqrt(intm1*intm1+4.*SQR(ca)*
                        (SQR(ca)-SQR(cay)))));
	} else {ca=c; cax=c; cay=c; cfx=c; cfy=c;}
}


// Delete memory from arrays

