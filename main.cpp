#include<fstream>
#include<iostream>
#include<sstream>   // Contains a set of types & functions for string
#include<cstdlib>
#include<math.h>
#include<exception> // Call error exception

#include "defs.hpp"
#include "timeIntegral.hpp"

void MUSCL2D(meshblock &dom, int nb, int step);
void savearray(meshblock &dom, real*** array, string arrname);
void save4py(meshblock &dom, int id, real t);
// Multiple static grids related stuffs
void locate_bounds(meshblock &dom);
void boundary(meshblock &dom,real**** Q,int nvar);
// AMR related stuffs
void admesh(meshblock &dom);
// Constrained transport for solenoidal magnetic field
void CT2D(meshblock &dom,real**** q,real**** Bi,int step,real**** Binew);
// Apply diffusivity to all variables for enhanced stability
void diffusivity(real**** U,real eta,int nb,int nxmin,int nxmax,int nymin,int nymax,int nvar);

int main() {
	meshblock dom; // Domain using arrays
	string IC,limiter,fluxMth;
	real CFL;
	int nx, ny, nvar;

	ifstream inputFile("input.txt");
	if (inputFile.is_open()) {
		// Input Domain-related properties
		inputFile >> nx >> ny >> nvar;
		real lenx, leny;
		inputFile >> lenx >> leny;
		inputFile >> CFL;
		// Must account for 'ghost' cells thus, +2
		dom.setSize(nx+2*nghosts,ny+2*nghosts,nvar,lenx,leny);
		cout<<"Domain size is "<<dom.nx-2<<" x "<<dom.ny-2<<" with nvar="<<nvar<<endl;
		cout<<"CFL="<<CFL<<", base dx="<<dom.dx[0]<<" & dy="<<dom.dy[0]<<endl;
		// Choose types of solver & problem
		if (!MAG_field) {
			cout<<"-----Hydrodynamics solver-----"<<endl;
		} else if (MAG_field) {
			cout<<"-----MHD solver-----"<<endl;
		}
		inputFile >> IC >> limiter >> fluxMth;
		dom.setType(IC,limiter,fluxMth);
		cout<<"Problems chosen: "<<dom.IC<<endl;
		cout<<"Limiter selected: "<<dom.limiter<<endl;
		cout<<"Types of Riemann Solver chosen: "<<dom.fluxMth<<endl;
	} else {
		cout<<"Unable to open the file..."<<endl;
		cout<<"***Terminating program***"<<endl;
		throw exception();
	}

	// Use existing data or not
	int choice;
	cout<<"For plotting, ";
	cout<<"(1) Use existing data (2) Run simulation and obtain new data"<<endl;
	cin>>choice;
	if (choice==1) {
		ostringstream ssx, ssy;
	        ssx<<nx; string n_x=ssx.str();
	        ssy<<ny; string n_y=ssy.str();
	        string plot_var;
		cout<<"0) Density; 1) u; 2) v; 3) w; 4) P; ";
	        if (MAG_field) {
        	        cout<<"5) Bx; 6) By; 7) Bz;"<<endl;
	        } else {
        	        cout<<endl;
	        }
	        cout<<"Choose variables for plotting: ";cin>>plot_var;
	        string command="python3 call_plot.py "+n_x+" "+ n_y+" "+plot_var+" "+dom.IC;
	        system(command.c_str());
		exit(0);	
	} else if (choice!=1 and choice!=2) {
		cout<<"Error as wrong selection!"<<endl;
		throw exception();
	}

	
	// Set up initial conditions
	dom.basegrid();
	for (int nb=0; nb<dom.lastActive; nb++) {
		if (dom.leafs[nb]) dom.W2U(nb);
	}

	// Initial adaption of mesh
	cout<<"--> leafs be4 admesh: ";
	for (int nb=0;nb<dom.lastActive;nb++) {cout<<dom.leafs[nb]<<" ";} cout<<endl;
	cout<<"--> Location of blocks be4 admesh: ";
	for (int nb=0;nb<dom.lastActive;nb++) { if (dom.leafs[nb]) cout<<"("<<dom.lp[nb][0]<<":"
		<<dom.icoord[nb][0]<<","<<dom.icoord[nb][1]<<") ";}
        cout<<endl;
	if (maxlevs>1) {admesh(dom);}
	cout<<"--> leafs aft admesh: ";
	for (int nb=0;nb<dom.lastActive;nb++) {cout<<dom.leafs[nb]<<" ";} cout<<endl;
	cout<<"--> Location of blocks aft admesh: ";
	for (int nb=0;nb<dom.lastActive;nb++) { if (dom.leafs[nb]) cout<<"("<<dom.lp[nb][0]<<":"
		<<dom.icoord[nb][0]<<","<<dom.icoord[nb][1]<<") ";} 
	cout<<endl;

	// Main loop--------------------------------------
	cout<<"Entering main loop..."<<endl;
	int count=0, simu_count=0; dom.count=count;
	real lambda1=0, lambda2=0, lambda;
	real t=0.0, dt=0.0;
        real res=0.;
	while (t<=dom.tEnd and count<=2000) {
		locate_bounds(dom);
	        boundary(dom,dom.U,dom.nvar);	

		if (t==0) {
		cout<<"0) # of blocks="<<dom.lastActive<<" where maxblocks="<<maxblocks<<", nbleafs="<<dom.nbleafs<<endl;
		cout<<"1) Calc. nbound="<<dom.nbounds<<"; ";
		for (int nb=0;nb<dom.nbounds;nb++) {
			cout<<"("<<dom.innerbounds[nb][0]<<","<<dom.innerbounds[nb][1]<<","<<dom.innerbounds[nb][2]<<") ";
		} cout<<" end;"<<endl;} 		

		// RK2 1st step: update Us
		for (int nb=0; nb<dom.lastActive; nb++) {
   		  if (dom.leafs[nb]) {
		    MUSCL2D(dom,nb,1);		
		    for (int i = 0; i < dom.nx; i++) {
		      for (int j = 0; j < dom.ny; j++) {
		        for (int k = 0; k < dom.nvar; k++) {
			  if ((i>=dom.nxminb and i<=dom.nxmax) and (j>=dom.nyminb and j<=dom.nymax)) {
				  res=(dom.ff[i+1][j][k][nb]-dom.ff[i][j][k][nb])/dom.dx[dom.lp[nb][0]]
					  +(dom.gg[i][j+1][k][nb]-dom.gg[i][j][k][nb])/dom.dy[dom.lp[nb][0]];
			  } else {res=0.;}
            		  dom.Us[i][j][k][nb] = RK2ndv1(dom.dt,dom.Us[i][j][k][nb],dom.U[i][j][k][nb],res,1);
		        }
		      }
		    }
		  }
		}  
		// Apply constrained transport
		if (CT_mtd) {
			boundary(dom,dom.Bi,2);
			boundary(dom,dom.ff,dom.nvar);
			boundary(dom,dom.gg,dom.nvar);
                        CT2D(dom,dom.Us,dom.Bi,1,dom.Bis);
		}
		// Boundary conditions on Us
                boundary(dom,dom.Us,dom.nvar);

		// RK2 2nd step: update U
		for (int nb=0; nb<dom.lastActive; nb++) {
		  if (dom.leafs[nb]) {
                    MUSCL2D(dom,nb,2);
                    for (int i = 0; i < dom.nx; i++) {
                      for (int j = 0; j < dom.ny; j++) {
                        for (int k = 0; k < dom.nvar; k++) {
			  if ((i>=dom.nxminb and i<=dom.nxmax) and (j>=dom.nyminb and j<=dom.nymax)) {
                                  res=(dom.ff[i+1][j][k][nb]-dom.ff[i][j][k][nb])/dom.dx[dom.lp[nb][0]]
					  +(dom.gg[i][j+1][k][nb]-dom.gg[i][j][k][nb])/dom.dy[dom.lp[nb][0]];
                          } else {res=0.;}
                          dom.U[i][j][k][nb] = RK2ndv1(dom.dt,dom.Us[i][j][k][nb],dom.U[i][j][k][nb],res,2);;
                        }
                      }
                    }
		    // Apply viscosity for stability (only at very refined regions)
		    if (diffusion) {if (dom.lp[nb][0]>=maxlevs-2) {
			    diffusivity(dom.U,dom.eta,nb,dom.nxmin,dom.nxmax,dom.nymin,dom.nymax,dom.nvar);}
		    }
		  }
		} 
		if (CT_mtd) {
			boundary(dom,dom.ff,dom.nvar);
                        boundary(dom,dom.gg,dom.nvar);
                        CT2D(dom,dom.U,dom.Bis,2,dom.Bi);
		}
		// Boundary conditions on U
                boundary(dom,dom.U,dom.nvar);

		// Calculate time for advancement		
		lambda=1.E6;
		for (int nb=0; nb<dom.lastActive; nb++) {
		  if (dom.leafs[nb]) {
		    dom.U2W("main",nb);
		    lambda1=0.; lambda2=0.; 
		    for (int i=dom.nxmin; i<=dom.nxmax; i++) {
		      for (int j=dom.nymin; j<=dom.nymax; j++) {
		        dom.wavespeed(i,j,nb);
			lambda1=max(lambda1,dom.speed+abs(dom.cfx));
			lambda2=max(lambda2,dom.speed+abs(dom.cfy));				
		      }
		    }		 		
		    lambda1=max(lambda1,lambda2);
		    lambda=min(lambda,min(dom.dx[dom.lp[nb][0]],dom.dy[dom.lp[nb][0]])/lambda1);
		  }
		}
		t=t+dt; dom.dt=dt;
		cout<<count<<") At t="<<t<<", dt="<<dt<<" and # of blocks are "<<dom.lastActive<<endl;
		count=count+1; dom.count=count;
		dt=CFL*lambda;
		if (dt<0 or isnan(dt)) {
			cout<<"Error as invalid dt produced!!!"<<endl;
			throw exception();
		}

		if (count%20==0) {
			simu_count+=1;
			save4py(dom,simu_count,t);
		}

		// Update mesh
		if (maxlevs>1) {admesh(dom);}
	}

	// Output data
	for (int nb=0;nb<dom.lastActive;nb++) {
		dom.W2U(nb);
	}
	save4py(dom,simu_count+1,t);	

	// For post-processing of data
	ostringstream ssx, ssy;
	ssx<<nx*sqrt(nbroots); string n_x=ssx.str();
	ssy<<ny*sqrt(nbroots); string n_y=ssy.str();
	string plot_var;
	cout<<"0) Density; 1) u; 2) v; 3) w; 4) P; ";
	if (MAG_field) {
		cout<<"5) Bx; 6) By; 7) Bz;"<<endl;
	} else {
		cout<<endl;
	}
	cout<<"Choose variables for plotting: ";cin>>plot_var;
	string command="python3 call_plot.py "+n_x+" "+ n_y+" "+plot_var+" "+dom.IC;
	system(command.c_str());

	return 0;
}
