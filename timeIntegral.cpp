#include "defs.hpp"

real RK2ndv1(real dt,real Qs,real Q,real res,int step) {
	if (step==1) {
		return Q-dt*res;
	} else if (step==2) {
		return 0.5*(Q+Qs-dt*res);
	}
}
