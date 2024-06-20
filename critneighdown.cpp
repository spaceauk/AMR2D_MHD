// Refinement criteria to allow maximum difference of 1 level betw
// neighbouring blocks when coarsening (i.e. only with icoarse)

#include "defs.hpp"

void critneighdown(int** lp, bool* leafs, bool* icoarse, int lastActive) {
	int n1;

	for (int nb=0; nb<lastActive; nb++) {
	  if (icoarse[nb]) {
	    // Left
	    n1=lp[nb][4];
	    if (n1!=-1) {
	      if (!leafs[n1]) icoarse[nb]=false;
	    }
	    // Right
	    n1=lp[nb][5];
	    if (n1!=-1) {
              if (!leafs[n1]) icoarse[nb]=false;
            }
	    // Bottom
	    n1=lp[nb][6];
            if (n1!=-1) {
              if (!leafs[n1]) icoarse[nb]=false;
            }
	    // Top
	    n1=lp[nb][7];
            if (n1!=-1) {
              if (!leafs[n1]) icoarse[nb]=false;
            }
	  }
	}
}
