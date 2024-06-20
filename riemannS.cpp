#include<vector>
#include<exception>

#include "defs.hpp"
#include "riemannS.hpp"

void riemannS(string fluxMth,vector<real> wL,vector<real> wR,real gamma,char direc,vector<real> &flux) {
	vector<real> wL_int(8),wR_int(8);
	for (int i=0; i<size(wL); i++) {
		wL_int[i]=wL[i]; wR_int[i]=wR[i];
	}
	if (direc=='x') {
	} else if (direc=='y'){
		wL_int[1]=wL[2]; // Swapped velocity direction
		wL_int[2]=wL[1];
		wR_int[1]=wR[2]; 
		wR_int[2]=wR[1];
		if (MAG_field) {
			wL_int[5]=wL[6]; // Swapped magnetic direction
			wL_int[6]=wL[5];
			wR_int[5]=wR[6]; 
			wR_int[6]=wR[5];
		}
	} else {
		cout<<"Wrong direction("<<direc<<") selected!!! \n";
		throw exception();
	}

	vector<real> flux_int(size(flux));
	if (fluxMth=="RUSA") {
		RUSA(wL_int,wR_int,gamma,flux_int);
	} else if (fluxMth=="HLLE") {
                HLLE(wL_int,wR_int,gamma,flux_int);
	} else if (fluxMth=="HLLC") {
		HLLC(wL_int,wR_int,gamma,flux_int);
	} else if (fluxMth=="HLLD") {
		HLLD(wL_int,wR_int,gamma,flux_int);
	} else if (fluxMth=="ROE") {
		if (MAG_field) {cout<<"ROE not ready with magnetic field!"<<endl; throw exception();}
		ROE(wL_int,wR_int,gamma,flux_int);
	} else {
		cout<<"Error as no appropriate Riemann solver chosen! \n";
		throw exception();
	}

	for (int i=0; i<size(wL); i++){flux[i]=flux_int[i];}
	if (direc=='x') {
	} else if (direc=='y'){
		flux[1]=flux_int[2]; // Swapped velocity direction
		flux[2]=flux_int[1];
		if (MAG_field) {
			flux[5]=flux_int[6]; // Swapped magnetic firection
			flux[6]=flux_int[5];
		}
	}
}

void RUSA(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux){
	// Left states
	real vnL=wL[1];
	real BnL=wL[5];
	real caL,aL,canL,cfL;
	wavespeeds(wL,gamma,5,aL,caL,canL,cfL);
	real eL=(wL[4]/(wL[0]*(gamma-1)))+0.5*MAG(wL[1],wL[2],wL[3])
			+0.5*MAG(wL[5],wL[6],wL[7])/wL[0];
	eL=eL*wL[0];
	// Specific enthalpy
	real HL=(eL+wL[4]+0.5*MAG(wL[5],wL[6],wL[7]))/wL[0];
	vector<real> qL(size(wL)), FL(size(wL));
	qL=wL;
	qL[1]=wL[0]*wL[1]; qL[2]=wL[0]*wL[2]; qL[3]=wL[0]*wL[3];
	qL[4]=eL;

	// Right states
        real vnR=wR[1];
        real BnR=wR[5];
	real caR,aR,canR,cfR;
        wavespeeds(wR,gamma,5,aR,caR,canR,cfR);
        real eR=(wR[4]/(wR[0]*(gamma-1)))+0.5*MAG(wR[1],wR[2],wR[3])
                        +0.5*MAG(wR[5],wR[6],wR[7])/wR[0];
        eR=eR*wR[0];
        // Specific enthalpy
        real HR=(eR+wR[4]+0.5*MAG(wR[5],wR[6],wR[7]))/wR[0];
        vector<real> qR(size(wR)), FR(size(wR));
        qR=wR;
        qR[1]=wR[0]*wR[1]; qR[2]=wR[0]*wR[2]; qR[3]=wR[0]*wR[3];
        qR[4]=eR;

	// Left fluxes
	real ptL=wL[4]+0.5*MAG(wL[5],wL[6],wL[7]);
	FL[0]=wL[0]*vnL;
	FL[1]=wL[0]*vnL*wL[1]+ptL-BnL*wL[5];
	FL[2]=wL[0]*vnL*wL[2]-BnL*wL[6];
	FL[3]=wL[0]*vnL*wL[3]-BnL*wL[7];
	FL[4]=wL[0]*vnL*HL-BnL*(wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]);
	FL[5]=vnL*wL[5]-BnL*wL[1];
	FL[6]=vnL*wL[6]-BnL*wL[2];
	FL[7]=vnL*wL[7]-BnL*wL[3];
	// Right fluxes
        real ptR=wR[4]+0.5*MAG(wR[5],wR[6],wR[7]);
        FR[0]=wR[0]*vnR;
        FR[1]=wR[0]*vnR*wR[1]+ptR-BnR*wR[5];
        FR[2]=wR[0]*vnR*wR[2]-BnR*wR[6];
        FR[3]=wR[0]*vnR*wR[3]-BnR*wR[7];
        FR[4]=wR[0]*vnR*HR-BnR*(wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7]);
        FR[5]=vnR*wR[5]-BnR*wR[1];
        FR[6]=vnR*wR[6]-BnR*wR[2];
        FR[7]=vnR*wR[7]-BnR*wR[3];

	// Compute RUSA flux
	real intm2=cfL+abs(wL[1]), intm3=cfR+abs(wR[1]);
	real intm=max(intm2,intm3);
	for (int k=0; k<size(wL); k++) {
		flux[k]=0.5*(FL[k]+FR[k])-0.5*intm*(qR[k]-qL[k]);
	}
}

void HLLE(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux){
	// Left states
	real vnL=wL[1];
	real BnL=wL[5];
	real caL,aL,canL,cfL;
	wavespeeds(wL,gamma,5,aL,caL,canL,cfL);
	real eL=(wL[4]/(wL[0]*(gamma-1)))+0.5*MAG(wL[1],wL[2],wL[3])
			+0.5*MAG(wL[5],wL[6],wL[7])/wL[0];
	eL=eL*wL[0];
	// Specific enthalpy
	real HL=(eL+wL[4]+0.5*MAG(wL[5],wL[6],wL[7]))/wL[0];
	vector<real> qL(size(wL)), FL(size(wL));
	qL=wL;
	qL[1]=wL[0]*wL[1]; qL[2]=wL[0]*wL[2]; qL[3]=wL[0]*wL[3];
	qL[4]=eL;

	// Right states
        real vnR=wR[1];
        real BnR=wR[5];
	real caR,aR,canR,cfR;
        wavespeeds(wR,gamma,5,aR,caR,canR,cfR);
        real eR=(wR[4]/(wR[0]*(gamma-1)))+0.5*MAG(wR[1],wR[2],wR[3])
                        +0.5*MAG(wR[5],wR[6],wR[7])/wR[0];
        eR=eR*wR[0];
        // Specific enthalpy
        real HR=(eR+wR[4]+0.5*MAG(wR[5],wR[6],wR[7]))/wR[0];
        vector<real> qR(size(wR)), FR(size(wR));
        qR=wR;
        qR[1]=wR[0]*wR[1]; qR[2]=wR[0]*wR[2]; qR[3]=wR[0]*wR[3];
        qR[4]=eR;

	// Left fluxes
	real ptL=wL[4]+0.5*MAG(wL[5],wL[6],wL[7]);
	FL[0]=wL[0]*vnL;
	FL[1]=wL[0]*vnL*wL[1]+ptL-BnL*wL[5];
	FL[2]=wL[0]*vnL*wL[2]-BnL*wL[6];
	FL[3]=wL[0]*vnL*wL[3]-BnL*wL[7];
	FL[4]=wL[0]*vnL*HL-BnL*(wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]);
	FL[5]=vnL*wL[5]-BnL*wL[1];
	FL[6]=vnL*wL[6]-BnL*wL[2];
	FL[7]=vnL*wL[7]-BnL*wL[3];
	// Right fluxes
        real ptR=wR[4]+0.5*MAG(wR[5],wR[6],wR[7]);
        FR[0]=wR[0]*vnR;
        FR[1]=wR[0]*vnR*wR[1]+ptR-BnR*wR[5];
        FR[2]=wR[0]*vnR*wR[2]-BnR*wR[6];
        FR[3]=wR[0]*vnR*wR[3]-BnR*wR[7];
        FR[4]=wR[0]*vnR*HR-BnR*(wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7]);
        FR[5]=vnR*wR[5]-BnR*wR[1];
        FR[6]=vnR*wR[6]-BnR*wR[2];
        FR[7]=vnR*wR[7]-BnR*wR[3];

	// Compute the Roe Averages
	vector<real> v_roe(3), B_roe(3);
	real RT,r,Xfac,Yfac;
	RoeAvg(wL,wR,RT,r,v_roe,B_roe,Xfac,Yfac);
	real H=(HL+RT*HR)/(1+RT);
	// Calculate Roe Eigenvalues
	vector<real> lambda(7);
	RoeEigen(gamma,r,v_roe,H,B_roe,Xfac,Yfac,lambda);

	// Wave speed estimates
	real SL=min(vnL-cfL,lambda[0]);
	SL=MIN(SL,0.);
	real SR=max(vnR+cfR,lambda[6]);
	SR=MAX(SR,0.);

	// Compute HLLE flux
	for (int k=0; k<size(wL); k++) {
		flux[k]=(SR*FL[k]-SL*FR[k]+SL*SR*(qR[k]-qL[k]))/(SR-SL);
	}
}

// HLLC flux from T.J. Linde
void HLLC(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux){
	// Left states
	real vnL=wL[1];
	real BnL=wL[5];
	real caL,aL,canL,cfL;
	wavespeeds(wL,gamma,5,aL,caL,canL,cfL);
	real eL=(wL[4]/(wL[0]*(gamma-1)))+0.5*MAG(wL[1],wL[2],wL[3])
			+0.5*MAG(wL[5],wL[6],wL[7])/wL[0];
	eL=eL*wL[0];
	// Specific enthalpy
	real HL=(eL+wL[4]+0.5*MAG(wL[5],wL[6],wL[7]))/wL[0];
	vector<real> qL(size(wL)), FL(size(wL));
	qL=wL;
	qL[1]=wL[0]*wL[1]; qL[2]=wL[0]*wL[2]; qL[3]=wL[0]*wL[3];
	qL[4]=eL;

	// Right states
        real vnR=wR[1];
        real BnR=wR[5];
	real caR,aR,canR,cfR;
        wavespeeds(wR,gamma,5,aR,caR,canR,cfR);
        real eR=(wR[4]/(wR[0]*(gamma-1)))+0.5*MAG(wR[1],wR[2],wR[3])
                        +0.5*MAG(wR[5],wR[6],wR[7])/wR[0];
        eR=eR*wR[0];
        // Specific enthalpy
        real HR=(eR+wR[4]+0.5*MAG(wR[5],wR[6],wR[7]))/wR[0];
        vector<real> qR(size(wR)), FR(size(wR));
        qR=wR;
        qR[1]=wR[0]*wR[1]; qR[2]=wR[0]*wR[2]; qR[3]=wR[0]*wR[3];
        qR[4]=eR;

	// Left fluxes
	real ptL=wL[4]+0.5*MAG(wL[5],wL[6],wL[7]);
	FL[0]=wL[0]*vnL;
	FL[1]=wL[0]*vnL*wL[1]+ptL-BnL*wL[5];
	FL[2]=wL[0]*vnL*wL[2]-BnL*wL[6];
	FL[3]=wL[0]*vnL*wL[3]-BnL*wL[7];
	FL[4]=wL[0]*vnL*HL-BnL*(wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]);
	FL[5]=vnL*wL[5]-BnL*wL[1];
	FL[6]=vnL*wL[6]-BnL*wL[2];
	FL[7]=vnL*wL[7]-BnL*wL[3];
	// Right fluxes
        real ptR=wR[4]+0.5*MAG(wR[5],wR[6],wR[7]);
        FR[0]=wR[0]*vnR;
        FR[1]=wR[0]*vnR*wR[1]+ptR-BnR*wR[5];
        FR[2]=wR[0]*vnR*wR[2]-BnR*wR[6];
        FR[3]=wR[0]*vnR*wR[3]-BnR*wR[7];
        FR[4]=wR[0]*vnR*HR-BnR*(wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7]);
        FR[5]=vnR*wR[5]-BnR*wR[1];
        FR[6]=vnR*wR[6]-BnR*wR[2];
        FR[7]=vnR*wR[7]-BnR*wR[3];

	// Compute the Roe Averages
	vector<real> v_roe(3), B_roe(3);
	real RT,r,Xfac,Yfac;
	RoeAvg(wL,wR,RT,r,v_roe,B_roe,Xfac,Yfac);
	real H=(HL+RT*HR)/(1+RT);
	// Calculate Roe Eigenvalues
	vector<real> lambda(7);
	RoeEigen(gamma,r,v_roe,H,B_roe,Xfac,Yfac,lambda);

	// Wave speed estimates
	real SL=min(vnL-cfL,lambda[0]);
	real SR=max(vnR+cfR,lambda[6]);
	real cf_roe=lambda[6]-v_roe[0];

	// Compute HLLC flux (Linde 1998 version 1)
	// -> uses Roe avg and heuristic linear function w.r.t. gap of jump
	real SM=v_roe[0];
	real sumdif1=0.;
	for (int k=0;k<size(qL);k++) {sumdif1+=abs(qR[k]-qL[k]);}
	real sumdif2=0.;
	for (int k=0;k<size(qL);k++) {
		sumdif2+=abs((FR[k]-FL[k])-v_roe[0]*(qR[k]-qL[k]));
	}
	real interm=1.-sumdif2/(sumdif1+eps);
	real alpha=MAX(0.,interm);
	if (alpha>1 or alpha<0 or isnan(alpha)) {
		cout<<"Error as alpha="<<alpha<<", SM="<<SM<<", SD1="<<sumdif1<<", SD2="<<sumdif2<<endl;
		throw exception();
	}
	if (SR-SL==0) {
		cout<<"Error as SR-SL=0 will give NaN value!"<<endl;
		throw exception();
	}
	vector<real> FLs(8),FRs(8);
	for (int k=0;k<8;k++) {
		FLs[k]=((SL*((1-alpha)*SR+alpha*v_roe[0]))/(SR-SL))*(qR[k]-qL[k])-(SL/(SR-SL))*FR[k]+((SR)/(SR-SL))*FL[k];
		FRs[k]=((SR*((1-alpha)*SL+alpha*v_roe[0]))/(SR-SL))*(qR[k]-qL[k])-(SL/(SR-SL))*FR[k]+((SR)/(SR-SL))*FL[k];
	}

	// Compute HLLC flux
	for (int k=0; k<size(wL); k++) {
		if (SL>0.) {
			flux[k]=FL[k];
		} else if (SL<=0. and SM>=0.) {
			flux[k]=FLs[k];
		} else if (SM<=0. and SR>=0.) {
			flux[k]=FRs[k];
		} else if (SR<0.) {
			flux[k]=FR[k];
		}
	}
}

void ROE(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux){
	real xnorm=1.0, ynorm=0.0;
        // Left states
        real vnL=wL[1]*xnorm+wL[2]*ynorm;
        real BnL=wL[5]*xnorm+wL[6]*ynorm;
	real caL,aL,canL,cfL;
        wavespeeds(wL,gamma,5,aL,caL,canL,cfL);
        real eL=(wL[4]/(wL[0]*(gamma-1)))+0.5*MAG(wL[1],wL[2],wL[3])
                        +0.5*MAG(wL[5],wL[6],wL[7])/wL[0];
        eL=eL*wL[0];
        // Specific enthalpy
        real HL=(eL+wL[4]+0.5*MAG(wL[5],wL[6],wL[7]))/wL[0];
        vector<real> qL(size(wL)), FL(size(wL));
        qL=wL;
        qL[1]=wL[0]*wL[1]; qL[2]=wL[0]*wL[2]; qL[3]=wL[0]*wL[3];
        qL[4]=eL;

        // Right states
        real vnR=wR[1]*xnorm+wR[2]*ynorm;
        real BnR=wR[5]*xnorm+wR[6]*ynorm;
	real caR,aR,canR,cfR;
        wavespeeds(wR,gamma,5,aR,caR,canR,cfR);
        real eR=(wR[4]/(wR[0]*(gamma-1)))+0.5*MAG(wR[1],wR[2],wR[3])
                        +0.5*MAG(wR[5],wR[6],wR[7])/wR[0];
        eR=eR*wR[0];
        // Specific enthalpy
        real HR=(eR+wR[4]+0.5*MAG(wR[5],wR[6],wR[7]))/wR[0];
        vector<real> qR(size(wR)), FR(size(wR));
        qR=wR;
        qR[1]=wR[0]*wR[1]; qR[2]=wR[0]*wR[2]; qR[3]=wR[0]*wR[3];
        qR[4]=eR;

        // Steps from Pg. 376 of Toro book (Not suited for MHD!)
        // 1) Compute Roe Averages
        real RT=sqrt(wR[0]/wL[0]);
        real r=sqrt(wR[0]*wL[0]);
        real u=(wL[1]+RT*wR[1])/(1+RT);
        real v=(wL[2]+RT*wR[2])/(1+RT);
        real w=(wL[3]+RT*wR[3])/(1+RT);
        /*real Bx=(sqrt(wR[0])*wL[5]+sqrt(wL[0])*wR[5]) 
           /(sqrt(wL[0])+sqrt(wR[0]));   
        real By=(sqrt(wR[0])*wL[6]+sqrt(wL[0])*wR[6]) 
           /(sqrt(wL[0])+sqrt(wR[0]));
        real Bz=(sqrt(wR[0])*wL[7]+sqrt(wL[0])*wR[7]) 
           /(sqrt(wL[0])+sqrt(wR[0]));*/
        real H=(HL+RT*HR)/(1+RT);
        real a=sqrt((gamma-1)*(H-0.5*MAG(u,v,w)));

	// 2) Compute averaged eigenvalues according to (11.58)
        real lambda[5];
        lambda[0]=u-a; lambda[4]=u+a;
        lambda[1]=u; lambda[2]=u; lambda[3]=u;

        // 3) Compute the averaged right eigenvectors according to (11.59)
        real K[5][8]={{1.,u-a,v,w,H-u*a,0.,0.,0.},
                {1.,u,v,w,0.5*MAG(u,v,w),0.,0.,0.},
                {0.,0.,1.,0.,v,0.,0.,0.},
                {0.,0.,0.,1.,w,0.,0.,0.},
                {1.,u+a,v,w,H+u*a,0.,0.,0.}};

        // 4) Compute the wave strengths according to (11.68)-(11.70) - doesn't account for y flux?
        real alpha[5];
        alpha[2]=(qR[2]-qL[2])-v*(qR[0]-qL[0]);
        alpha[3]=(qR[3]-qL[3])-w*(qR[0]-qL[0]);
        real du5bar=(qR[4]-qL[4])-(alpha[2])*v-(alpha[3])*w;
        alpha[1]=((gamma-1)/pow(a,2))*((qR[0]-qL[0])*(H-u*u)+u*(qR[1]-qL[1])-du5bar);
        alpha[0]=0.5*(1/a)*((qR[0]-qL[0])*(u+a)-(qR[1]-qL[1])-a*alpha[1]);
        alpha[4]=(qR[0]-qL[0])-(alpha[0]+alpha[1]);


        // 5) COmpute flux according to (11.27)-(11.29)
        // Left fluxes
        real ptL=wL[4]+0.5*MAG(wL[5],wL[6],wL[7]);
        FL[0]=wL[0]*vnL;
        FL[1]=wL[0]*vnL*wL[1]+xnorm*ptL-BnL*wL[5];
        FL[2]=wL[0]*vnL*wL[2]+ynorm*ptL-BnL*wL[6];
        FL[3]=wL[0]*vnL*wL[3]-BnL*wL[7];
        FL[4]=wL[0]*vnL*HL-BnL*(wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]);
        FL[5]=vnL*wL[5]-BnL*wL[1];
        FL[6]=vnL*wL[6]-BnL*wL[2];
        FL[7]=vnL*wL[7]-BnL*wL[3];
        // Right fluxes
        real ptR=wR[4]+0.5*MAG(wR[5],wR[6],wR[7]);
        FR[0]=wR[0]*vnR;
        FR[1]=wR[0]*vnR*wR[1]+xnorm*ptR-BnR*wR[5];
        FR[2]=wR[0]*vnR*wR[2]+ynorm*ptR-BnR*wR[6];
        FR[3]=wR[0]*vnR*wR[3]-BnR*wR[7];
        FR[4]=wR[0]*vnR*HR-BnR*(wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7]);
        FR[5]=vnR*wR[5]-BnR*wR[1];
        FR[6]=vnR*wR[6]-BnR*wR[2];
        FR[7]=vnR*wR[7]-BnR*wR[3];

        // Compute ROE flux
        for (int k=0; k<size(wL); k++) {
                real summ=0.;
                for (int j=0; j<5; j++) {summ+=alpha[j]*abs(lambda[j])*K[j][k];}
                flux[k]=0.5*(FL[k]+FR[k])-0.5*summ;
        }
}

void HLLD(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux){
	vector<real> qL(size(wL)), FL(size(wL));
	vector<real> wLs(size(wL)), wLss(size(wL));
	vector<real> qLs(size(wL)), qLss(size(wL));
	vector<real> qR(size(wL)), FR(size(wL));
	vector<real> wRs(size(wL)), wRss(size(wL));
	vector<real> qRs(size(wL)), qRss(size(wL));

	// Left states
	real vnL=wL[1];
	real BnL=wL[5];
	real caL,aL,canL,cfL;
        wavespeeds(wL,gamma,5,aL,caL,canL,cfL);
	real eL=(wL[4]/(wL[0]*(gamma-1)))+0.5*MAG(wL[1],wL[2],wL[3])
			+0.5*MAG(wL[5],wL[6],wL[7])/wL[0];
	eL=eL*wL[0];
	// Specific enthalpy
	real HL=(eL+wL[4]+0.5*MAG(wL[5],wL[6],wL[7]))/wL[0];
	qL=wL;
	qL[1]=wL[0]*wL[1]; qL[2]=wL[0]*wL[2]; qL[3]=wL[0]*wL[3];
	qL[4]=eL;

	// Right states
        real vnR=wR[1];
        real BnR=wR[5];
	real caR,aR,canR,cfR;
        wavespeeds(wR,gamma,5,aR,caR,canR,cfR);
        real eR=(wR[4]/(wR[0]*(gamma-1)))+0.5*MAG(wR[1],wR[2],wR[3])
                        +0.5*MAG(wR[5],wR[6],wR[7])/wR[0];
        eR=eR*wR[0];
        // Specific enthalpy
        real HR=(eR+wR[4]+0.5*MAG(wR[5],wR[6],wR[7]))/wR[0];
        qR=wR;
        qR[1]=wR[0]*wR[1]; qR[2]=wR[0]*wR[2]; qR[3]=wR[0]*wR[3];
        qR[4]=eR;

	// Left fluxes
	real ptL=wL[4]+0.5*MAG(wL[5],wL[6],wL[7]);
	FL[0]=wL[0]*vnL;
	FL[1]=wL[0]*vnL*wL[1]+ptL-BnL*wL[5];
	FL[2]=wL[0]*vnL*wL[2]-BnL*wL[6];
	FL[3]=wL[0]*vnL*wL[3]-BnL*wL[7];
	FL[4]=wL[0]*vnL*HL-BnL*(wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]);
	FL[5]=vnL*wL[5]-BnL*wL[1];
	FL[6]=vnL*wL[6]-BnL*wL[2];
	FL[7]=vnL*wL[7]-BnL*wL[3];
	// Right fluxes
        real ptR=wR[4]+0.5*MAG(wR[5],wR[6],wR[7]);
        FR[0]=wR[0]*vnR;
        FR[1]=wR[0]*vnR*wR[1]+ptR-BnR*wR[5];
        FR[2]=wR[0]*vnR*wR[2]-BnR*wR[6];
        FR[3]=wR[0]*vnR*wR[3]-BnR*wR[7];
        FR[4]=wR[0]*vnR*HR-BnR*(wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7]);
        FR[5]=vnR*wR[5]-BnR*wR[1];
        FR[6]=vnR*wR[6]-BnR*wR[2];
        FR[7]=vnR*wR[7]-BnR*wR[3];

	// Obtain wavespeeds 
	real SL=min(vnL,vnR)-max(cfL,cfR);
	real SR=max(vnL,vnR)+max(cfL,cfR);
 	real SM=((SR-vnR)*wR[0]*vnR-(SL-vnL)*wL[0]*vnL-ptR+ptL)/
		((SR-vnR)*wR[0]-(SL-vnL)*wL[0]);
	
	// Average total pressure in Riemann fan (Eq. 41)
	real pts=((SR-vnR)*wR[0]*ptL-(SL-vnL)*wL[0]*ptR
		+wL[0]*wR[0]*(SR-vnR)*(SL-vnL)*(vnR-vnL))
		/((SR-vnR)*wR[0]-(SL-vnL)*wL[0]);	

	// Star quantities
	// (a) Left *
	wLs[0]=wL[0]*(SL-vnL)/(SL-SM);
	wLs[1]=SM;
	wLs[5]=wL[5];
	if ((wL[0]*(SL-wL[1])*(SL-SM)-pow(wL[5],2))!=0) {
		wLs[2]=wL[2]-wL[5]*wL[6]*((SM-wL[1])/(wL[0]*(SL-wL[1])*(SL-SM)-pow(wL[5],2)));
		wLs[3]=wL[3]-wL[5]*wL[7]*((SM-wL[1])/(wL[0]*(SL-wL[1])*(SL-SM)-pow(wL[5],2)));
		wLs[6]=wL[6]*((wL[0]*pow(SL-wL[1],2)-pow(wL[5],2))/(wL[0]*(SL-wL[1])*(SL-SM)-pow(wL[5],2)));
		wLs[7]=wL[7]*((wL[0]*pow(SL-wL[1],2)-pow(wL[5],2))/(wL[0]*(SL-wL[1])*(SL-SM)-pow(wL[5],2)));
	} else {
		wLs[2]=wL[2]; wLs[3]=wL[3];
		wLs[6]=0.; wLs[7]=0.;
	}
	real eLs=(((SL-wL[1])*eL-ptL*wL[1]+pts*SM+wL[5]*
		(dot_prod(wL,1,3,wL,5,7)-
		 dot_prod(wLs,1,3,wLs,5,7)))/(SL-SM));
	// (b) Right *
	wRs[0]=wR[0]*(SR-vnR)/(SR-SM);
        wRs[1]=SM;
        wRs[5]=wR[5];
        if ((wR[0]*(SR-wR[1])*(SR-SM)-pow(wR[5],2))!=0) {
                wRs[2]=wR[2]-wR[5]*wR[6]*((SM-wR[1])/(wR[0]*(SR-wR[1])*(SR-SM)-pow(wR[5],2)));
                wRs[3]=wR[3]-wR[5]*wR[7]*((SM-wR[1])/(wR[0]*(SR-wR[1])*(SR-SM)-pow(wR[5],2)));
                wRs[6]=wR[6]*((wR[0]*pow(SR-wR[1],2)-pow(wR[5],2))/(wR[0]*(SR-wR[1])*(SR-SM)-pow(wR[5],2)));
                wRs[7]=wR[7]*((wR[0]*pow(SR-wR[1],2)-pow(wR[5],2))/(wR[0]*(SR-wR[1])*(SR-SM)-pow(wR[5],2)));
        } else {
                wRs[2]=wR[2]; wRs[3]=wR[3];
                wRs[6]=0; wRs[7]=0;
        }
        real eRs=(((SR-wR[1])*eR-ptR*wR[1]+pts*SM+wR[5]*
                (dot_prod(wR,1,3,wR,5,7)-
                 dot_prod(wRs,1,3,wRs,5,7)))/(SR-SM));

	// Star-star quantities
	// (a) Left **
	wLss[0]=wLs[0];
	wLss[1]=SM;
	wLss[2]=(sqrt(wLs[0])*wLs[2]+sqrt(wRs[0])*wRs[2]+
		(wRs[6]-wLs[6])*sign(wL[5]))/(sqrt(wLs[0])+sqrt(wRs[0]));
	wLss[3]=(sqrt(wLs[0])*wLs[3]+sqrt(wRs[0])*wRs[3]+
		(wRs[7]-wLs[7])*sign(wL[5]))/(sqrt(wLs[0])+sqrt(wRs[0]));
	wLss[5]=wLs[5];
	wLss[6]=(sqrt(wLs[0])*wRs[6]+sqrt(wRs[0])*wLs[6]+
		sqrt(wLs[0]*wRs[0])*(wRs[2]-wLs[2])*sign(wL[5]))
		/(sqrt(wLs[0])+sqrt(wRs[0]));
	wLss[7]=(sqrt(wLs[0])*wRs[7]+sqrt(wRs[0])*wLs[7]+
		sqrt(wLs[0]*wRs[0])*(wRs[3]-wLs[3])*sign(wL[5]))
		/(sqrt(wLs[0])+sqrt(wRs[0]));
	real eLss=eLs-sqrt(wLs[0])*(dot_prod(wLs,1,3,wLs,5,7)-
		dot_prod(wLss,1,3,wLss,5,7))*sign(wL[5]);
	// (b) Right **
	for (int k=0; k<size(wL); k++) {wRss[k]=wLss[k];}
	real eRss=eRs+sqrt(wRs[0])*(dot_prod(wRs,1,3,wRs,5,7)-
		dot_prod(wRss,1,3,wRss,5,7))*sign(wR[5]);

	// Covert from primitive to conservative
	w2u(wLs,qLs,gamma);
	w2u(wLss,qLss,gamma); 
	qLs[4]=eLs; qLss[4]=eLss;
	w2u(wRs,qRs,gamma);
	w2u(wRss,qRss,gamma); 
	qRs[4]=eRs; qRss[4]=eRss;

	// Propagation of Alfven waves in intermediate states (Eq. 51)
	real SLs=SM-abs(wL[5])/sqrt(wLs[0]);
	real SRs=SM+abs(wR[5])/sqrt(wRs[0]);

	// Compute HLLD flux
	for (int k=0; k<size(wL); k++) {
		if (SL>0) {
			flux[k]=FL[k];
		} else if (SL<=0. and SLs>=0.) {
			flux[k]=FL[k]+SL*qLs[k]-SL*qL[k];
		} else if (SLs<=0. and SM>=0.) {
			flux[k]=FL[k]+SLs*qLss[k]-(SLs-SL)*qLs[k]-SL*qL[k];
		} else if (SM<=0. and SRs>=0.) {
			flux[k]=FR[k]+SRs*qRss[k]-(SRs-SR)*qRs[k]-SR*qR[k];
		} else if (SRs<=0. and SR>=0.) {
			flux[k]=FR[k]+SR*qRs[k]-SR*qR[k];
		} else if (SR<0.) {
			flux[k]=FR[k];
		}
	}
}
