// Update pointers after refining
#include "defs.hpp"

void updatelpup(meshblock &dom, int nb,
		int son1, int son2, int son3, int son4) {
	int nsr, nlev2;
	// Level of refinement of the sons
	nlev2=dom.lp[nb][0]+1;

	// Son 1-----------------------------------------------------
	dom.lp[son1][0]=nlev2;  // Refinment level
	dom.lp[son1][1]=nb;     // Father ID
	dom.lp[son1][3]=0;      // Sibling ID
	dom.lp[son1][5]=son2;   // Right
	dom.lp[son1][7]=son3;	// Up
	// Left ****
	if (dom.lp[nb][4]==-1 or dom.leafs[dom.lp[nb][4]]) { 
		// If left neighbor is at grid boundary or less refined
		dom.lp[son1][4]=dom.lp[nb][4];
	} else {
		nsr=dom.lp[dom.lp[nb][4]][2]+1; // | n0 | n1 | s0 | s1 | : thus need to +1 as n1 is beside son1(s0) instead of n0.
		dom.lp[son1][4]=nsr;
		dom.lp[nsr][5]=son1;
	}
	// Down ****
	if (dom.lp[nb][6]==-1 or dom.leafs[dom.lp[nb][6]]) {
		dom.lp[son1][6]=dom.lp[nb][6];
	} else {
		nsr=dom.lp[dom.lp[nb][6]][2]+2;
		dom.lp[son1][6]=nsr;
		dom.lp[nsr][7]=son1;
	}
	// Son 2------------------------------------------------------
	dom.lp[son2][0]=nlev2;
        dom.lp[son2][1]=nb;
        dom.lp[son2][3]=1;
        dom.lp[son2][4]=son1;  // Left
        dom.lp[son2][7]=son4;   // Up
        // Right ****
        if (dom.lp[nb][5]==-1 or dom.leafs[dom.lp[nb][5]]) {
                dom.lp[son2][5]=dom.lp[nb][5];
        } else {
                nsr=dom.lp[dom.lp[nb][5]][2]; 
                dom.lp[son2][5]=nsr;
                dom.lp[nsr][4]=son2;
        }
        // Down ****
        if (dom.lp[nb][6]==-1 or dom.leafs[dom.lp[nb][6]]) {
                dom.lp[son2][6]=dom.lp[nb][6];
        } else {
                nsr=dom.lp[dom.lp[nb][6]][2]+3;
                dom.lp[son2][6]=nsr;
                dom.lp[nsr][7]=son2;
        }
	// Son 3------------------------------------------------------
        dom.lp[son3][0]=nlev2;
        dom.lp[son3][1]=nb;
        dom.lp[son3][3]=2;
        dom.lp[son3][5]=son4;  // Right
        dom.lp[son3][6]=son1;   // Down
        // Left ****
        if (dom.lp[nb][4]==-1 or dom.leafs[dom.lp[nb][4]]) {
                dom.lp[son3][4]=dom.lp[nb][4];
        } else {
                nsr=dom.lp[dom.lp[nb][4]][2]+3;
                dom.lp[son3][4]=nsr;
                dom.lp[nsr][5]=son3;
        }
        // Up ****
        if (dom.lp[nb][7]==-1 or dom.leafs[dom.lp[nb][7]]) {
                dom.lp[son3][7]=dom.lp[nb][7];
        } else {
                nsr=dom.lp[dom.lp[nb][7]][2];
                dom.lp[son3][7]=nsr;
                dom.lp[nsr][6]=son3;
        }
	// Son 4------------------------------------------------------
        dom.lp[son4][0]=nlev2;
        dom.lp[son4][1]=nb;
        dom.lp[son4][3]=3;
        dom.lp[son4][4]=son3;  // Left
        dom.lp[son4][6]=son2;   // Down
        // Right ****
        if (dom.lp[nb][5]==-1 or dom.leafs[dom.lp[nb][5]]) {
                dom.lp[son4][5]=dom.lp[nb][5];
        } else {
                nsr=dom.lp[dom.lp[nb][5]][2]+2;
                dom.lp[son4][5]=nsr;
                dom.lp[nsr][4]=son4;
        }
        // Up ****
        if (dom.lp[nb][7]==-1 or dom.leafs[dom.lp[nb][7]]) {
                dom.lp[son4][7]=dom.lp[nb][7];
        } else {
                nsr=dom.lp[dom.lp[nb][7]][2]+1;
                dom.lp[son4][7]=nsr;
                dom.lp[nsr][6]=son4;
        }
	// Note that the neighbors marked **** are given on a coarser level
	// if not available at same level
	
	// Father no longer leaf
	dom.leafs[nb]=false;
	// All sons are leafs
	dom.leafs[son1]=true;
	dom.leafs[son2]=true;
	dom.leafs[son3]=true;
	dom.leafs[son4]=true;
	dom.nbleafs+=4; 
	
	// De-mark for refining/coarsening
	dom.iref[nb]=false;
	dom.icoarse[nb]=false;
	dom.icoarse[son1]=false;
	dom.icoarse[son2]=false;
	dom.icoarse[son3]=false;
	dom.icoarse[son4]=false;
}
