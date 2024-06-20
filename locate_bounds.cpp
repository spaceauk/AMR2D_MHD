//================================================================
// Locate neighbors to pass B.C.s.
// It makes a list for the inner boundaries and a list with the 
// outer boundaries according to:
// - bounds = 0: Left, 1: Right, 2: Bottom, 3: Top
// - innerbounds = (ninternal,3)  0: source ID, 1: dest. ID, 
//                                2: direc.
// - extbounds =   (nextertot,5)  0: source ID, 1: sorce rank,
//   			          2: dest. ID,  3: dest. rank,
//  			          4: direc.
// Note that if dest. ID = -1 here, it means that the boundary
// lies on the grid boundary and thus no destination block.
//
// Note that for the direc here will not start from 0!
// direc. = -1 : left box boundary
//        = -2 : right box boundary
//        = -3 : bottom box boundary
//        = -4 : top box boundary
//        =  1 : left boundary at same res
//        =  2 : left boundary at higher res 1
//        =  3 : left boundary at higher res 2 (only MPI)
//        =  4 : right boundary at same res
//        =  5 : right boundary at higher res 1
//        =  6 : right boundary at higher res 2 (only MPI)
//        =  7 : bottom boundary at same res
//        =  8 : bottom boundary at higher res 1
//        =  9 : bottom boundary at higher res 2 (only MPI)
//        = 10 : top boundary at same res
//        = 11 : top boundary at higher res 1
//        = 12 : top boundary at higher res 2 (only MPI)
//===============================================================
#include "defs.hpp"
#include<iostream>

void locate_bounds(meshblock &dom) {
	int n1,n2;
	bool bounds[dom.lastActive][4]; // Marked boundary that have been evaluated to prevent re-evaluated in subsequent iteration
	// Reset boundary flag, lists and their size
	for (int nb=0;nb<dom.lastActive;nb++) {
		for (int i=0;i<4;i++) {bounds[nb][i]=false;}
	}

	// Deallocate old innerbounds for differentt size
	for (int i=0;i<dom.oldnbleafs*4;i++) {
		delete [] dom.innerbounds[i]; 
	}
	delete [] dom.innerbounds;		
	dom.innerbounds = new int*[dom.nbleafs*4]; 
	for (int i=0;i<dom.nbleafs*4;i++) {
		dom.innerbounds[i] = new int[3];
	}
	dom.oldnbleafs=dom.nbleafs;
	int ninternal=-1; 

	// Locate internal boundaries
	for (int nb=0;nb<dom.lastActive;nb++) {
  	  if (dom.leafs[nb]) {		
	    //-----------------
	    //  Left boundary - thus only look at lp[nb][4] here
	    //-----------------
	    if (!bounds[nb][0]) {
	      // If box boundary
	      if (dom.lp[nb][4]==-1) {
	        bounds[nb][0]=true;
		ninternal+=1;
		dom.innerbounds[ninternal][0]=nb;
		dom.innerbounds[ninternal][1]=-1; 
		dom.innerbounds[ninternal][2]=-1;					
		// Only if neighbor is at least at same level
	      } else if (dom.lp[dom.lp[nb][4]][0]>=dom.lp[nb][0]) {
	        n1=dom.lp[nb][4]; // ID of left neighbor   					
		if (dom.leafs[n1]) { // If neighbor is at same resolution
		  ninternal+=1;
		  dom.innerbounds[ninternal][0]=nb;
		  dom.innerbounds[ninternal][1]=n1;
		  dom.innerbounds[ninternal][2]=1;
		  //  Left boundary of nb is shared with the right boundary of n1 cell. Thus marked both to avoid revisiting this
		  //  same boundary in the subsequent iteration.
		  bounds[nb][0]=true;
		  bounds[n1][1]=true;
		} else { // If neighbor is at higher resolution
		  bounds[nb][0]=true;
		  n1=dom.lp[n1][2]+1; // Need to +1 as the left boundary of nb is shared by ID2 of n1.
		  //  The vertical boundaries have ID difference of 2.       -------------      
		  //  For e.g., left boundary contains ID1 & ID3 while       | ID3 | ID4 |
		  //  right boundary contains ID2 & ID4.                     -------------
		  //  Thus, n2=n1+2 as n1 & n2 are the ID for the            | ID1 | ID2 |
		  //  vertical boundaries.                                   -------------
		  n2=n1+2; 
		  ninternal+=1;
		  dom.innerbounds[ninternal][0]=nb;
		  dom.innerbounds[ninternal][1]=n1;
		  dom.innerbounds[ninternal][2]=2;
		  bounds[n1][1]=true;
		  bounds[n2][1]=true;
		}
	      }
	    }

	    //-----------------
            //  Right boundary - thus only look at lp[nb][5] here
            //-----------------
            if (!bounds[nb][1]) {
	      // If box boundary
	      if (dom.lp[nb][5]==-1) {
	        bounds[nb][1]=true;
		ninternal+=1;
		dom.innerbounds[ninternal][0]=nb;
		dom.innerbounds[ninternal][1]=-1;
		dom.innerbounds[ninternal][2]=-2;
	      // Only if neighbor is at least at same level
	      } else if (dom.lp[dom.lp[nb][5]][0]>=dom.lp[nb][0]) {
		n1=dom.lp[nb][5];
		if (dom.leafs[n1]) { // If neighbor is at same res
		  ninternal+=1;
		  dom.innerbounds[ninternal][0]=nb;
		  dom.innerbounds[ninternal][1]=n1;
		  dom.innerbounds[ninternal][2]=4;
		  bounds[nb][1]=true;
		  bounds[n1][0]=true;
		} else { // If neighbor is at higher res
		  bounds[nb][1]=true;
		  n1=dom.lp[n1][2]; // Do not need to +1 as the right boundary of nb is shared by ID1 of n1.
		  n2=n1+2;
		  ninternal+=1;
		  dom.innerbounds[ninternal][0]=nb;
		  dom.innerbounds[ninternal][1]=n1;
		  dom.innerbounds[ninternal][2]=5;
		  bounds[n1][0]=true;
		  bounds[n2][0]=true;
		}
	      } 
	    }

	    //-----------------
            //  Bottom boundary - thus only look at lp[nb][6] here
            //-----------------
	    if (!bounds[nb][2]) {
 	      // If box boundary
	      if (dom.lp[nb][6]==-1) {
	        bounds[nb][2]=true;
		ninternal+=1;
		dom.innerbounds[ninternal][0]=nb;
		dom.innerbounds[ninternal][1]=-1;
		dom.innerbounds[ninternal][2]=-3;
	      // Only if neighbor is at least at same level
	      } else if (dom.lp[dom.lp[nb][6]][0]>=dom.lp[nb][0]) {
	        n1=dom.lp[nb][6];
	        if (dom.leafs[n1]) { // If neighbor is at same res
		  ninternal+=1;
		  dom.innerbounds[ninternal][0]=nb;
		  dom.innerbounds[ninternal][1]=n1;
		  dom.innerbounds[ninternal][2]=7;
		  bounds[nb][2]=true;
		  bounds[n1][3]=true;
		} else { // If neighbor is at higher res
	 	  bounds[nb][2]=true;
		  n1=dom.lp[n1][2]+2; // Bottom boundary of ID1 of nb is shared by ID3 of n1 below nb.
		  // For horizontal boundary, the ID differs by 1. E.g., ID1-ID2 & ID3-ID4
		  n2=n1+1;
		  ninternal+=1;
		  dom.innerbounds[ninternal][0]=nb;
		  dom.innerbounds[ninternal][1]=n1;
		  dom.innerbounds[ninternal][2]=8;
		  bounds[n1][3]=true;
		  bounds[n2][3]=true;
		}
	      }
	    }
	    
	    //-----------------
            //  Top boundary - thus only look at lp[nb][6] here
            //-----------------
	    if (!bounds[nb][3]) {
	      // If box boundary
	      if (dom.lp[nb][7]==-1) {
	        bounds[nb][3]=true;
		ninternal+=1;
		dom.innerbounds[ninternal][0]=nb;
		dom.innerbounds[ninternal][1]=-1;
		dom.innerbounds[ninternal][2]=-4;
	      // Only if neighbor is at least at same level	
	      } else if (dom.lp[dom.lp[nb][7]][0]>=dom.lp[nb][0]) {
	        n1=dom.lp[nb][7];
		if (dom.leafs[n1]) { // If neighbor is at same res
		  ninternal+=1;
		  dom.innerbounds[ninternal][0]=nb;
		  dom.innerbounds[ninternal][1]=n1;
		  dom.innerbounds[ninternal][2]=10;
		  bounds[nb][3]=true;
		  bounds[n1][2]=true;
		} else { // If neighbor is at higher res
                  bounds[nb][3]=true;
		  n1=dom.lp[n1][2]; // Top boundary of nb is shared by ID1 of n1 above nb.
		  n2=n1+1;
		  ninternal+=1;
		  dom.innerbounds[ninternal][0]=nb;
		  dom.innerbounds[ninternal][1]=n1;
		  dom.innerbounds[ninternal][2]=11;
		  bounds[n1][2]=true;
		  bounds[n2][2]=true;
		}
	      }	    
	    }
	  }
	}
	dom.nbounds=ninternal+1; 
}
