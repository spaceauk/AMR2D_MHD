// Update pointers after coarsening

#include "defs.hpp"

void updatelpdown(meshblock &dom, int dad, 
		int son1, int son2, int son3, int son4) {
	int nb1;

	// Update the father
	dom.leafs[dad]=true;

	// Clear the sons
	dom.leafs[son1]=false;
	dom.leafs[son2]=false;
	dom.leafs[son3]=false;
	dom.leafs[son4]=false;
	dom.nbleafs-=3; 

	// Update neighbors
	// Left
	nb1=dom.lp[dad][4];
	if (nb1!=-1 and !dom.leafs[nb1]) { 
		// if neighbours more refined (assumed that adjacent block refinement diff is at most 1)
		// E.g.,  ----------------------
		//        | n2 | n3 |          |
		//        -----------   dad    |
		//        | n0 | n1 |          |
		//        ----------------------   where ni is neighbouring more refined cells
		nb1=dom.lp[dom.lp[nb1][2]][5]; // n1
		dom.lp[nb1][5]=dad;           
		dom.lp[dom.lp[nb1][7]][5]=dad; // n3
	}
	// Right
        nb1=dom.lp[dad][5];
        if (nb1!=-1 and !dom.leafs[nb1]) {
                nb1=dom.lp[nb1][2];
                dom.lp[nb1][4]=dad;
                dom.lp[dom.lp[nb1][7]][4]=dad;
        }
	// Bottom
        nb1=dom.lp[dad][6];
        if (nb1!=-1 and !dom.leafs[nb1]) {
                nb1=dom.lp[dom.lp[nb1][2]][7];
                dom.lp[nb1][7]=dad;
                dom.lp[dom.lp[nb1][5]][7]=dad;
        }
	// Top
        nb1=dom.lp[dad][7];
        if (nb1!=-1 and !dom.leafs[nb1]) {
                nb1=dom.lp[nb1][2];
                dom.lp[nb1][6]=dad;
                dom.lp[dom.lp[nb1][5]][6]=dad;
        }
}
