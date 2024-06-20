#include<exception> // Call error exception
#include<vector>		    

#include "defs.hpp"

void meshblock::IC2Dtype(int nb) {
	real tEnd, gamma;
	vector<real> rinit(4),pinit(4),Einit(4);
	vector<real> uinit(4),vinit(4),winit(4);
       	vector<real> Bxinit(4),Byinit(4),Bzinit(4);

	// Initialize the vectors to zero first
	rinit={0}; pinit={0}; Einit={0};
	uinit={0}; vinit={0}; winit={0};
	Bxinit={0}; Byinit={0}; Bzinit={0};

	// Set conditions for specific I.C.s selected
	if (IC=="SST") {
		if (nb==0) cout<<"2D Sod Shock tube config along x."<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={0.1,1.0,1.0,0.1};
		rinit={0.125,1.0,1.0,0.125};

	} else if (IC=="SST1") {
		if (nb==0) cout<<"2D Sod Shock tube config along x (reverse direction)."<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={1.0,0.1,0.1,1.0};
		rinit={1.0,0.125,0.125,1.0};

	} else if (IC=="SSTy") {
                if (nb==0) cout<<"2D Sod Shock tube config along y."<<endl;
                tEnd=0.2;
                gamma=2.0;
                pinit={0.1,0.1,1.0,1.0};
                rinit={0.125,0.125,1.0,1.0};

	} else if (IC=="BWx") {
		if (nb==0) cout<<"2D Brio & Wu shocktube config along x."<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={0.1,1.0,1.0,0.1};
		rinit={0.125,1.0,1.0,0.125};
		Bxinit={0.75,0.75,0.75,0.75};
		Byinit={-1.0,1.0,1.0,-1.0};
	
	} else if (IC=="BWx1") {
		if (nb==0) cout<<"2D Brio & Wu shocktube config along x. (opposite direction)"<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={1.0,0.1,0.1,1.0};
		rinit={1.0,0.125,0.125,1.0};
		Bxinit={0.75,0.75,0.75,0.75};
		Byinit={1.0,-1.0,-1.0,1.0};

	} else if (IC=="BWy") {
		if (nb==0) cout<<"2D Brio & Wu shocktube config along y."<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={0.1,0.1,1.0,1.0};
		rinit={0.125,0.125,1.0,1.0};
		Byinit={0.75,0.75,0.75,0.75};
		Bxinit={-1.0,-1.0,1.0,1.0};

	} else if (IC=="BWy1") {
		if (nb==0) cout<<"2D Brio & Wu shocktube config along y. (opposite direction)"<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={1.0,1.0,0.1,0.1};
		rinit={1.0,1.0,0.125,0.125};
		Byinit={0.75,0.75,0.75,0.75};
		Bxinit={1.0,1.0,-1.0,-1.0};

	} else if (IC=="rotor") {
		if (nb==0) cout<<"MHD Rotor problem in 2D"<<endl;
		gamma=1.4;
		tEnd=0.295;
		pinit[0]=0.5;
		Bxinit[0]=2.5/(sqrt(4*pi));
	} else {
		cout<<"No valid initial condition chosen!"<<endl;
		cout<<"Look at IC2Dtype file to choose the correct I.C.s!"<<endl;
		cout<<"***Program is stopped!***"<<endl;
		throw exception();
	}

	// Assign values to arrays (Must account for boundary cells)
	real x, y, radius;
	if (IC=="rotor") {
		for (int i=0; i<nx; i++) {
			x=((i+icoord[nb][0]-nghosts)+0.5)*dx[lp[nb][0]] - dx[0];
			for (int j=0; j<ny; j++) {
				y=((j+icoord[nb][1]-nghosts)+0.5)*dy[lp[nb][0]] - dy[0];
				radius=sqrt(pow(x-0.5,2)+pow(y-0.5,2));
				real fr=(23.-200.*radius)/3.;
				for (int k=0;k<nvar;k++) {W[i][j][k][nb]=0.;}
				W[i][j][4][nb]=pinit[0];
				W[i][j][5][nb]=Bxinit[0];
				if (radius<=0.1) {
					W[i][j][0][nb]=10.;
					W[i][j][1][nb]=-(y-0.5)/0.1;
					W[i][j][2][nb]=(x-0.5)/0.1;
				} else if (radius>0.1 and radius<0.115) {
					W[i][j][0][nb]=(1+9.*fr);
					W[i][j][1][nb]=(-(y-0.5)*fr/0.1);
					W[i][j][2][nb]=(x-0.5)*fr/0.1;
				} else {
					W[i][j][0][nb]=1.;					
				}
			}
		}

	} else { 
		for (int i=0; i<nx; i++) {
			x=((i+icoord[nb][0]-nghosts)+0.5)*dx[lp[nb][0]] - dx[0];
			for (int j=0; j<ny; j++) {
				y=((j+icoord[nb][1]-nghosts)+0.5)*dy[lp[nb][0]] - dy[0];

				if (x>=0.5 and y>=0.5) {
					W[i][j][0][nb]=rinit[0];
					W[i][j][1][nb]=uinit[0];
					W[i][j][2][nb]=vinit[0];
					W[i][j][3][nb]=winit[0];
					W[i][j][4][nb]=pinit[0];
					if (MAG_field) {
						W[i][j][5][nb]=Bxinit[0];
						W[i][j][6][nb]=Byinit[0];
						W[i][j][7][nb]=Bzinit[0];				 	 
				       	}
				} else if (x<0.5 and y>=0.5) {
					W[i][j][0][nb]=rinit[1];
					W[i][j][1][nb]=uinit[1];
					W[i][j][2][nb]=vinit[1];
					W[i][j][3][nb]=winit[1];
					W[i][j][4][nb]=pinit[1];
					if (MAG_field) {
						W[i][j][5][nb]=Bxinit[1];
						W[i][j][6][nb]=Byinit[1];
						W[i][j][7][nb]=Bzinit[1];
					}
				} else if (x<0.5 and y<0.5) {
					W[i][j][0][nb]=rinit[2];
					W[i][j][1][nb]=uinit[2];
					W[i][j][2][nb]=vinit[2];
					W[i][j][3][nb]=winit[2];
					W[i][j][4][nb]=pinit[2];
					if (MAG_field) {
						W[i][j][5][nb]=Bxinit[2];
						W[i][j][6][nb]=Byinit[2];
						W[i][j][7][nb]=Bzinit[2];
					}
				} else if (x>=0.5 and y<0.5) {
					W[i][j][0][nb]=rinit[3];
					W[i][j][1][nb]=uinit[3];
					W[i][j][2][nb]=vinit[3];
					W[i][j][3][nb]=winit[3];
					W[i][j][4][nb]=pinit[3];
					if (MAG_field) {
						W[i][j][5][nb]=Bxinit[3];
						W[i][j][6][nb]=Byinit[3];
						W[i][j][7][nb]=Bzinit[3];	
					}
				}	
			}
		}
		
	}

	// Assign values from initial condition to class parameters
	setParam(gamma,tEnd);

	// Initialize magnetic field for constrained transport method
	if (CT_mtd) {
		for (int i=0; i<nx; i++) {
			for (int j=0; j<ny; j++) {
				// 20 -----------
				//    | 10 | 11 |
				// 10 -----------
				//    | 00 | 01 | 
				// 00 -----------  Thus, 01e=0.5*(00c+01c)
				//   00   01   02
				if (i<nx-1) {Bi[i+1][j][0][nb]=0.5*(W[i][j][5][nb]+W[i+1][j][5][nb]);}
				else {Bi[0][j][0][nb]=Bi[1][j][0][nb];}
				if (j<ny-1) {Bi[i][j+1][1][nb]=0.5*(W[i][j][6][nb]+W[i][j+1][6][nb]);}
				else {Bi[i][0][1][nb]=Bi[i][1][1][nb];}
			}
		}
	}
}
