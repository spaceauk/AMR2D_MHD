//==============================================================
// Updates the mesh by refining and coarsening all leaf blocks
//==============================================================
#include<exception>
#include "defs.hpp"

void markref(meshblock &dom);
void critneighup(int** lp, bool* iref, int lastActive);
void critneighdown(int** lp, bool* leafs, bool* icoarse, int lastActive);
void refine(meshblock &dom,real**** Q,int nvar,int dad,int son1,int son2,int son3,int son4);
void updatelpup(meshblock &dom,int dad,int son1,int son2,int son3,int son4);
void coarsen(meshblock &dom,real**** Q,int nvar,int dad,int son1,int son2,int son3,int son4);
void updatelpdown(meshblock &dom,int dad,int son1,int son2,int son3,int son4);

void admesh(meshblock &dom) {
	int brother[3];
	int dad, son1, son2, son3, son4;
	bool print=true;
	// Flag for refinement/coarsen
	markref(dom);
	// Mark for refining on proximity of higher resolution grid
	critneighup(dom.lp,dom.iref,dom.lastActive); 

	// Refine marked blocks by order of levels to ensure neighbours are updated correctly
	for (int nl=0; nl<maxlevs-1; nl++) {
	  for (int nb=0; nb<dom.lastActive; nb++) {
	    if (dom.lp[nb][0]==nl and dom.iref[nb]) {
	      if (dom.lp[nb][2]==0) {
	        // 1st time a son is produced
		son1=dom.lastActive; 
		son2=son1+1;
		son3=son1+2;
		son4=son1+3;
		dom.lp[nb][2]=son1;
		// Update total number of available blocks
		dom.lastActive+=4;
		if (dom.lastActive-1>=maxblocks) {
			cout<<"nblocks exceeded max limit="<<maxblocks<<" where lastActive="<<dom.lastActive<<endl;
			throw exception();
		} 
	      } else {
	        son1=dom.lp[nb][2];
		son2=dom.lp[son1][5]; // right  rmb ID:   3|4
		son3=dom.lp[son1][7]; // up               1|2   
		son4=dom.lp[son3][5]; // right
	      }
	      // Apply refinement
	      if (print) {cout<<"Refining to give son block: "; print=false;}
	      refine(dom,dom.U,dom.nvar,nb,son1,son2,son3,son4);
	      updatelpup(dom,nb,son1,son2,son3,son4);
	      if (CT_mtd) refine(dom,dom.Bi,2,nb,son1,son2,son3,son4);
	    }
	  }
	} 
	if (!print) cout<<endl;

	// Verify that the blocks to be coarsen will have max of 1 level difference only
	critneighdown(dom.lp,dom.leafs,dom.icoarse,dom.lastActive);

	print=true;
	// Coarsen in order of levels	
	for (int nl=maxlevs-1;nl>=1;nl--) {
	  for (int nb=0;nb<dom.lastActive;nb++) {
	    if (dom.icoarse[nb] and dom.lp[nb][3]==0 and 
			    dom.lp[nb][0]==nl) {
	      dad=dom.lp[nb][1];
	      son1=nb;
	      son2=dom.lp[nb][5];
	      son3=dom.lp[nb][7];
	      son4=dom.lp[son2][7]; 
	      brother[0]=dom.icoarse[son2];
	      brother[1]=dom.icoarse[son3];
	      brother[2]=dom.icoarse[son4];
	      // Only when entire group of 4 siblings are marked for coarsening is it allowed
	      if (brother[0] and brother[1] and brother[2]) {
		if (print) {cout<<"Coarsening to give dad block: "; print=false;}
	        coarsen(dom,dom.U,dom.nvar,dad,son1,son2,son3,son4);
		updatelpdown(dom,dad,son1,son2,son3,son4);
		if (CT_mtd) coarsen(dom,dom.Bi,2,dad,son1,son2,son3,son4);
	      }
	    }
	  }
	} 
	if (!print) cout<<endl;

}
