// Include C++ headers
#include <cmath>
#include <cstdint>

#include <iostream>
using namespace std;

// Declarations or constants
#ifdef DOUBLE_PREC
using real = double;
#define eps 1.0E-15
#define Nprec 15
#else
using real = float;
#define eps 1.0E-6
#define Nprec 6
#endif
#define pi 3.14159265359
#define R 287.0 // Gas constant for air

// Global variables
const int nghosts=4;   // # of ghost cells
const int maxlevs=5;
const int nbroots=4;    // # of root blocks
const int maxblocks= nbroots>pow(nbroots,maxlevs-1) ? nbroots : pow(4,maxlevs)*2;
// (a) Additional features
const bool MAG_field=true;
const bool CT_mtd=true;                // Apply constrained transport method to ensure divergence free magnetic field
const bool diffusion=false;            // Apply diffusion on all variables for stability

// Forward declarations needed for class files
class meshblock {
	public: 
		// Domain properties
		int nx, ny, nvar;
		int nxmin, nxmax, nxminb, nxmaxb;
		int nymin, nymax, nyminb, nymaxb;
		int nxp1, nxm1, nx2, nx2p1, nx2m1;
		int nyp1, nym1, ny2, ny2p1, ny2m1;
		// Fluid Properties
		real gamma;
		real tEnd;
		real eta=0.001; // Diffusivity
		// AMR parameters
                int lastActive=nbroots;
		int nbounds; // Number of boundaries 
                int** lp;
                int** icoord;
                bool* iref;
                bool* icoarse;
                bool* leafs;
                real* dx;
                real* dy;
                int nbleafs, oldnbleafs;
                int** innerbounds;

		// Related arrays
		real**** U;      // Conservative variables 
		real**** Us;
	       	real**** W;      // Primitive variables

		real*** dwdx;
		real*** dwdy;
		real*** wxL;
		real*** wxR;
		real*** wyL;
		real*** wyR;
		real**** ff;       // Flux along x
		real**** gg;       // Flux along y
		double**** Bi;	  // Cell edge magnetic field for constrained transport
		double**** Bis;
		double**** EMF;
		// Types of solver & problem
		string IC;
		string limiter;
		string fluxMth;
	
		// Calulated quantities after input
		real dt;   // Constant timing used throughout regardless of different refinement level
		real speed, c;
		real ca, cax, cay, cfx, cfy, beta; // When there is M-field
		int count;
		// For post-processing
		string fname;	
		
	void setSize(int no_x, int no_y, int no_v, real lenx, real leny);
	void setParam (real gammaa, real tEnd0);
	void setType (string IC0, string limiter0, string fluxMth0);
	void U2W(string loc,int nb);
	void Us2W(string loc,int nb);
	void W2U(int nb);
	void wavespeed(int i, int j, int nb);

	// AMR-related stuffs
        void basegrid();
        void IC2Dtype(int nb);
};


// All functions
#define MIN(a, b) ((a)<(b) ? (a):(b))
#define MAX(a, b) ((a)>(b) ? (a):(b))
#define MAG(x,y,z) (x*x+y*y+z*z)
#define SQR(x) (x*x)
