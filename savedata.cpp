#include<fstream>  // input & output
#include<iomanip>  // setw()
#include<sstream>		    

#include"defs.hpp"

void save4py(meshblock &dom, int id, real t) {

	stringstream ss;
	ss<<setw(3)<<setfill('0')<<id;
	string idd=ss.str();
	dom.fname="./data/"+dom.IC+idd+"t_"+to_string(t)+".dat";

	ofstream Wdata;
	int width=16;
	Wdata.open(dom.fname);	
	for (int nb=0;nb<dom.lastActive;nb++) {
	if (dom.leafs[nb]) {
		
	Wdata<<"#"<<setw(width-1)<<"x"<<setw(width)<<"y"<<setw(width)<<"xdelta"<<setw(width)<<"ydelta"<<
    	       setw(width)<<"density"<<setw(width)<<"U"<<setw(width)<<"v"<<setw(width)<<"w"<<setw(width)<<"pressure";
        if (MAG_field) Wdata<<setw(width)<<"Bx"<<setw(width)<<"By"<<setw(width)<<"Bz";
        Wdata<<endl;

	for (int j=dom.nymin; j<=dom.nymax; j++) {
		// Note that the x & y here are curated for my gnuplot...
		real y=((j+dom.icoord[nb][1]-nghosts)+0.5)*dom.dy[dom.lp[nb][0]] - dom.dy[0];
		for (int i=dom.nxmin; i<=dom.nxmax; i++) {
			real x=((i+dom.icoord[nb][0]-nghosts)+0.5)*dom.dx[dom.lp[nb][0]] - dom.dx[0];
			// x,y,r,u,v,w,p 
			Wdata<<setw(width)<<x<<setw(width)<<y<<
				setw(width)<<dom.dx[dom.lp[nb][0]]/2.<<
				setw(width)<<dom.dy[dom.lp[nb][0]]/2.<<
				setw(width)<<dom.W[i][j][0][nb]<<
				setw(width)<<dom.W[i][j][1][nb]<<setw(width)<<dom.W[i][j][2][nb]<<setw(width)<<dom.W[i][j][3][nb]<<
				setw(width)<<dom.W[i][j][4][nb];
			if (MAG_field) { // Magnetic field (Bx,By,Bz)
				Wdata<<setw(width)<<dom.W[i][j][5][nb]<<setw(width)<<dom.W[i][j][6][nb]<<setw(width)<<dom.W[i][j][7][nb];
			}
			Wdata<<endl;
		}
	}
	}
	}
	Wdata.close();
}

void savearray(meshblock &dom,real*** array, string arrname) {

	string ARRname="savefile/"+arrname;	
	real intm;
	ofstream arr;
	arr.open(ARRname);
	for (int k=0; k<dom.nvar; k++) {
		arr<<"At ivar="<<k<<endl;
		for (int j=0; j<dom.ny; j++) {
			for (int i=0; i<dom.nx; i++) {
				intm=array[i][j][k];
				arr<<intm<<" ";
			}
			arr<<endl;
		}
		arr<<"\n \n";
	}
	arr.close();

}
