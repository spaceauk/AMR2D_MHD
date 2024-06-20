// Mark all the leaf blocks for refinement or coarsening
#include "defs.hpp"

void criteria(meshblock &dom, int nb);

void markref(meshblock &dom) {
	// Initialize 
	for (int nb=0; nb<dom.lastActive; nb++) {
		dom.iref[nb]=false;
		dom.icoarse[nb]=false;
	}
	// Mark refine or coarsen based on criteria
	for (int nb=0; nb<dom.lastActive; nb++) {
		if (dom.leafs[nb] and dom.lp[nb][0]>=0) {
			criteria(dom,nb);
		}
	}
}
