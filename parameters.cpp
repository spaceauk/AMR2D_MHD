#include<exception>
#include "defs.hpp"

void meshblock::setParam(real gamma, real tEnd) {
        this->gamma=gamma;
        this->tEnd=tEnd;

	// Note that nx and ny are the total cells including ghost cells 
	// Before including ghost cells, the user input to the nx/ny must be even due to the nature of the code written here
	nxmin=nghosts; nxmax=nx-1-nghosts;	
	nymin=nghosts; nymax=ny-1-nghosts;
	if (nxmax%2!=1 or nymax%2!=1) {
		cout<<"Domain size must be even!"<<endl;
		if (nghosts%2!=0) {
	        	cout<<"With AMR, # of ghost cells must be even #."<<endl; 
			cout<<"For e.g., nghosts should be 2 and 4 for 2nd order and 4th order spatial orders."<<endl;	
		}
		//throw exception();
	}
	nxminb=nxmin-1; nxmaxb=nxmax+1;
	nyminb=nymin-1; nymaxb=nymax+1;
	// Assuming only two ghost cells,
	//  -------------------------------
        //  | 0 | 1 | 2 | ... |n-3|n-2|n-1|
        //  -------------------------------
        //    G   G                 G   G
        //        b                 b
        //            i   ...   i
        // So I will call boundary cells (b) as one cell outside the inner cells (i) and ghost cells (G) will contain the b cells too.
        // Take note of the < or <= in my loop based on following def
        nxp1=nxmaxb, nyp1=nymaxb;    // p1 = nmax+1
        nxm1=nxmax-1, nym1=nymax-1;  // m1 = nmax-1
        nx2=(nx-1)/2, ny2=(ny-1)/2;
        nx2m1=nx2-1, ny2m1=ny2-1;
        nx2p1=nx2+1, ny2p1=ny2+1;
}
