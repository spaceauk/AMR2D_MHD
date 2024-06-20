// Refinement criteria to allow maximum difference of 1 level betw
// neighbouring blocks when refining (i.e. only with irefup)

#include "defs.hpp"

void critneighup(int** lp, bool* iref, int lastActive) {
	int n1;

	for (int nl=maxlevs-1; nl>=0; nl--) {
	  for (int nb=0; nb<lastActive; nb++) {
	    if (iref[nb] and lp[nb][0]==nl) {
	      // Left
	      n1=lp[nb][4];
	      if (n1!=-1) {
		if (lp[n1][0]<lp[nb][0]) iref[n1]=true;	
	      }
	      // Right
	      n1=lp[nb][5];
	      if (n1!=-1) {
                if (lp[n1][0]<lp[nb][0]) iref[n1]=true;
              }
	      // Bottom
	      n1=lp[nb][6];
              if (n1!=-1) {
                if (lp[n1][0]<lp[nb][0]) iref[n1]=true;
              }
	      // Top
	      n1=lp[nb][7];
              if (n1!=-1) {
                if (lp[n1][0]<lp[nb][0]) iref[n1]=true;
              }
	    }
	  }
	}
}
