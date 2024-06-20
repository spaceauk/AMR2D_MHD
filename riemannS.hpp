real sign(real val);
real dot_prod(vector<real> v1,int v1s,int v1e,vector<real> v2,int v2s,int v2e);
void wavespeeds(vector<real> w,real gamma,int norm,real &a,real &ca,real &can,real &cf);
void RoeAvg(vector<real> wL,vector<real> wR,real &RT,real &r,vector<real> &v,vector<real> &B,real &Xfac,real &Yfac);
void RoeEigen(real gamma,real d,vector<real> v,real h,vector<real> B,real Xfac,real Yfac,vector<real> &lambda);
void w2u(vector<real> w,vector<real> &u,real gamma);
void RUSA(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux);
void HLLE(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux);
void HLLC(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux);
void ROE(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux);
void HLLD(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux);
