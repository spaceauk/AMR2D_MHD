// Setup base grid for AMR. Start with 4 blocks.
#include "defs.hpp"

void meshblock::basegrid() {
	cout<<"Setting up basegrid..."<<endl;

	// ----------------1st level------------------
	// Sibling numbering follows this convention:
	//                      Y
	//    -------------     ^
	//    | ID2 | ID3 |     |
	//    -------------     |
	//    | ID0 | ID1 |     +-------> X
	//    -------------
	// Initialized root blocks here
	// Note that my ID all also must minus one since the ID used here are use as index for lp to form quadtree structure
	// - Note that son is not the sibling as sibling always start with 0 while son is based along the 0 to lastActive range.
	lp[0][0]=0;    // Level of refinement
	lp[0][1]=0;    // Father
	lp[0][2]=0;    // Son (only 1st=ID0)
	lp[0][3]=0;    // Sibling ID (not refined start from ID0)
	lp[0][4]=-1;   // Neighbor left
	lp[0][5]=1;    // Neighbor right
	lp[0][6]=-1;   // Neighbor down
	lp[0][7]=2;    // Neighbor up (note that -1=grid boundary)
	// Coordinates of each block
	icoord[0][0]=0;
	icoord[0][1]=0;
	leafs[0]=true;
	iref[0]=true;

	lp[1][0]=0;    // Level of refinement
        lp[1][1]=0;    // Father
        lp[1][2]=0;    // Son (only 1st=ID1)
        lp[1][3]=1;    // Sibling ID
        lp[1][4]=0;    // Neighbor left
        lp[1][5]=-1;   // Neighbor right
        lp[1][6]=-1;   // Neighbor down
        lp[1][7]=3;    // Neighbor up (note that -1=grid boundary)
        // Coordinates of each block
        icoord[1][0]=nx-2*nghosts;
        icoord[1][1]=0;
        leafs[1]=true;
        iref[1]=true;

	lp[2][0]=0;    // Level of refinement
        lp[2][1]=0;    // Father
        lp[2][2]=0;    // Son (only 1st=ID1)
        lp[2][3]=0;    // Sibling ID
        lp[2][4]=-1;   // Neighbor left
        lp[2][5]=3;    // Neighbor right
        lp[2][6]=0;    // Neighbor down
        lp[2][7]=-1;   // Neighbor up (note that -1=grid boundary)
        // Coordinates of each block
        icoord[2][0]=0;
        icoord[2][1]=ny-2*nghosts;
        leafs[2]=true;
        iref[2]=true;

	lp[3][0]=0;    // Level of refinement
        lp[3][1]=0;    // Father
        lp[3][2]=0;    // Son (only 1st=ID1)
        lp[3][3]=0;    // Sibling ID
        lp[3][4]=2;    // Neighbor left
        lp[3][5]=-1;   // Neighbor right
        lp[3][6]=1;    // Neighbor down
        lp[3][7]=-1;   // Neighbor up (note that -1=grid boundary)
        // Coordinates of each block
        icoord[3][0]=nx-2*nghosts;
        icoord[3][1]=ny-2*nghosts;
        leafs[3]=true;
        iref[3]=true;

	// The resulting number of blocks and leafs
	nbleafs=4;
	oldnbleafs=nbleafs;
	innerbounds = new int*[nbleafs*4]; 
        for (int i=0;i<nbleafs*4;i++) {
                innerbounds[i] = new int[3];
        }

	cout<<"Setting up initial conditions..."<<endl;
	for (int nb=0; nb<lastActive; nb++) {
		IC2Dtype(nb);
	        W2U(nb);
	}

	// Refine basegrid ...
}
