/*
 *  blackoil.h -- header file for Black Oil Simulation
 *
 *  Created by peghoty 2011/09/09
 *
 *  Xiangtan University
 *  peghoty@163.com
 *
 */ 

#ifndef FSLS_BLACKOIL_HEADER
#define FSLS_BLACKOIL_HEADER

#define BlackOil_CSR  1
#define BlackOil_FULL 2

/*!
 * \struct fsls_OilGasParam
 * \brief Parameters for 2D 2 phase (Oil and Gas) simulator.
 */
typedef struct
{
   // parameters for mesh
   double            XX;      // [0,XX] is the scope in the X-direction of the solving domain
   double            YY;      // [0,YY] is the scope in the Y-direction of the solving domain
   double            ZZ;      // [0,ZZ] is the scope in the Z-direction of the solving domain
   int               Nx;      // number of grids in X-direction
   int               Ny;      // number of grids in Y-direction
   int               Nz;      // number of grids in Z-direction
   double            Dx;      // grid size in X-direction
   double            Dy;      // grid size in Y-direction
   double            Dz;      // grid size in Z-direction
   double            V;       // Volumn for each grid
   
   // parameters for porosity and permeability
   double            Phi;     // porosity
   double            Kx;      // permeability in X direction
   double            Ky;      // permeability in Y direction 
   
   // parameters for the Productio Well
   int               Wx;      // X Coordinate of Well Location
   int               Wy;      // Y Coordinate of Well Location
   double            WI;      // Well Index
   double            OR;      // Well Oil Production Rate Constraint (STB/day)
   double            BHP;     // Well BHP Constraint (Psia)
   double            Pw;      // The current well Pressure 
   double            Sw;      // The current well Saturation 

   // parameters for time control
   double            omega;
   double            yitap;
   double            yitas;
   double            dt0;     // The First time step
   double            dt;      // The Current time step
   double            Tmax;    // The maximum time step
   double            Tend;    // The Total simulation time 
   
   // parameters for Newton Iteration
   double            tolN;    // tolerance for Residual Norm
   double            tolP;    // tolerance for Pressure
   double            tolS;    // tolerance for Saturation
   
   // parameters for Number System Conversion
   double            alpha;   // British System -> Metric System 
   
   int               MaxStep; // The maximum simulation step allowed
   int               SucStep; // The number of successful   time steps
   int               UncStep; // The number of unconvergent time steps
   int               NCFlag;  // Mark whether the Newton iteration is convergenced. 1: yes; 0: no
   
   // variables saved for fast computing
   double            Tow;
   double            Tgw;
   
   int               postpro; // whether do the post-processing. 1: yes; 0: no
   
   int               NumUnconvType1; // number of unconvergence type 1 (Rnorm increses)
   int               NumUnconvType2; // number of unconvergence type 1 (reaches the maximal number of iteration)
   
} fsls_OilGasParam;

#define fsls_OilGasParamXX(oilgasparam)      ((oilgasparam) -> XX)
#define fsls_OilGasParamYY(oilgasparam)      ((oilgasparam) -> YY)
#define fsls_OilGasParamZZ(oilgasparam)      ((oilgasparam) -> ZZ)
#define fsls_OilGasParamNx(oilgasparam)      ((oilgasparam) -> Nx)
#define fsls_OilGasParamNy(oilgasparam)      ((oilgasparam) -> Ny)
#define fsls_OilGasParamNz(oilgasparam)      ((oilgasparam) -> Nz)
#define fsls_OilGasParamDx(oilgasparam)      ((oilgasparam) -> Dx)
#define fsls_OilGasParamDy(oilgasparam)      ((oilgasparam) -> Dy)
#define fsls_OilGasParamDz(oilgasparam)      ((oilgasparam) -> Dz)
#define fsls_OilGasParamV(oilgasparam)       ((oilgasparam) -> V)
#define fsls_OilGasParamPhi(oilgasparam)     ((oilgasparam) -> Phi)
#define fsls_OilGasParamKx(oilgasparam)      ((oilgasparam) -> Kx)
#define fsls_OilGasParamKy(oilgasparam)      ((oilgasparam) -> Ky)

#define fsls_OilGasParamWx(oilgasparam)      ((oilgasparam) -> Wx)
#define fsls_OilGasParamWy(oilgasparam)      ((oilgasparam) -> Wy)
#define fsls_OilGasParamWI(oilgasparam)      ((oilgasparam) -> WI)
#define fsls_OilGasParamOR(oilgasparam)      ((oilgasparam) -> OR)
#define fsls_OilGasParamBHP(oilgasparam)     ((oilgasparam) -> BHP)
#define fsls_OilGasParamPw(oilgasparam)      ((oilgasparam) -> Pw)
#define fsls_OilGasParamSw(oilgasparam)      ((oilgasparam) -> Sw)

#define fsls_OilGasParamOmega(oilgasparam)   ((oilgasparam) -> omega)
#define fsls_OilGasParamYitap(oilgasparam)   ((oilgasparam) -> yitap)
#define fsls_OilGasParamYitas(oilgasparam)   ((oilgasparam) -> yitas)
#define fsls_OilGasParamDt0(oilgasparam)     ((oilgasparam) -> dt0)
#define fsls_OilGasParamDt(oilgasparam)      ((oilgasparam) -> dt)
#define fsls_OilGasParamTmax(oilgasparam)    ((oilgasparam) -> Tmax)
#define fsls_OilGasParamTend(oilgasparam)    ((oilgasparam) -> Tend)

#define fsls_OilGasParamTolN(oilgasparam)    ((oilgasparam) -> tolN)
#define fsls_OilGasParamTolP(oilgasparam)    ((oilgasparam) -> tolP)
#define fsls_OilGasParamTolS(oilgasparam)    ((oilgasparam) -> tolS)

#define fsls_OilGasParamAlpha(oilgasparam)   ((oilgasparam) -> alpha)

#define fsls_OilGasParamMaxStep(oilgasparam) ((oilgasparam) -> MaxStep)
#define fsls_OilGasParamSucStep(oilgasparam) ((oilgasparam) -> SucStep)
#define fsls_OilGasParamUncStep(oilgasparam) ((oilgasparam) -> UncStep)
#define fsls_OilGasParamNCFlag(oilgasparam)  ((oilgasparam) -> NCFlag)

#define fsls_OilGasParamTow(oilgasparam)     ((oilgasparam) -> Tow)
#define fsls_OilGasParamTgw(oilgasparam)     ((oilgasparam) -> Tgw)
#define fsls_OilGasParamPostpro(oilgasparam) ((oilgasparam) -> postpro)

#define fsls_OilGasParamNumUnconvType1(oilgasparam) ((oilgasparam) -> NumUnconvType1)
#define fsls_OilGasParamNumUnconvType2(oilgasparam) ((oilgasparam) -> NumUnconvType2)


/*!
 * \struct fsls_OilGasData
 * \brief Data for 2D 2 phase (Oil and Gas) simulator.
 */
typedef struct
{
   fsls_OilGasParam  *oilgasparam;
   
   /* Oil-Gas properties which are depentdent of Pressure */
   
   // Mu relevant
   double *Muo;
   double *dMuo;
   double *Mug;
   double *dMug;
   
   // Rs relevant
   double *Rs;
   double *dRs;
   double *Rsx1;
   double *Rsx2;
   double *Rsy1;
   double *Rsy2;
   double *Rsnow;   // current Rs   
   
   // B relevant
   double *Bo;
   double *dBo;
   double *Bonow;   // current Bo
   double *Bg;
   double *dBg;
   double *Bgnow;   // current Bg   

   // relative permeabilities
   double *Kro;
   double *dKro;
   double *Krg;
   double *dKrg;
   
   /* Transmissibility and their derivatives */
   
   // x-direction, Oil
   double *Yox1;
   double *Yox2;
   double *dpYox1;
   double *dpYox2;
   double *dsYox1;
   double *dsYox2;
   double *ddpYox1;
   double *ddpYox2;
   double *ddsYox1;
   double *ddsYox2;
   
   // x-direction, Gas
   double *Ygx1;
   double *Ygx2;
   double *dpYgx1;
   double *dpYgx2;
   double *dsYgx1;
   double *dsYgx2;
   double *ddpYgx1;
   double *ddpYgx2;
   double *ddsYgx1;
   double *ddsYgx2;
   
   // y-direction, Oil
   double *Yoy1;
   double *Yoy2;
   double *dpYoy1;
   double *dpYoy2;
   double *dsYoy1;
   double *dsYoy2;
   double *ddpYoy1;
   double *ddpYoy2;
   double *ddsYoy1;
   double *ddsYoy2;  
   
   // y-direction, Gas
   double *Ygy1;
   double *Ygy2;
   double *dpYgy1;
   double *dpYgy2;
   double *dsYgy1;
   double *dsYgy2;
   double *ddpYgy1;
   double *ddpYgy2;
   double *ddsYgy1;
   double *ddsYgy2;

   double *Pbegin;    // previous time step Pressure
   double *Sbegin;    // previous time step Saturation
   double *Pnow;      // current time step Pressure
   double *Snow;      // current time step Saturation

   /* solve relevant */
   double *R;         // size of 2*Nx*Ny
   double *delta;     // size of 2*Nx*Ny
   double *deltap;    // correction of Pressure
   double *deltas;    // correction of Saturation
   double *Jacobi;    // Jacobian Matrix (Full format)
   
   /* data for post-processing */
   int     *It_array;  // number of Newton iteration array
   double  *Ts_array;  // Successful time Step array
   double  *Tn_array;  // T_now array
   double  *Pw_array;  // Well Pressure for each successful time step
   double  *Sw_array;  // Well Saturation for each successful time step
   double  *OR_array;  // Well Oil Production Rate for each successful time step   
   double **P_array;   // The Pressure   sequence for all the successful iteration 
   double **S_array;   // The Saturation sequence for all the successful iteration

} fsls_OilGasData;

#define fsls_OilGasDataOilgasparam(oilgasdata) ((oilgasdata) -> oilgasparam)
// Oil-Gas properties which are depentdent of Pressure
#define fsls_OilGasDataMuo(oilgasdata)         ((oilgasdata) -> Muo)
#define fsls_OilGasDataDMuo(oilgasdata)        ((oilgasdata) -> dMuo)
#define fsls_OilGasDataMug(oilgasdata)         ((oilgasdata) -> Mug)
#define fsls_OilGasDataDMug(oilgasdata)        ((oilgasdata) -> dMug)

#define fsls_OilGasDataRs(oilgasdata)          ((oilgasdata) -> Rs)
#define fsls_OilGasDataDRs(oilgasdata)         ((oilgasdata) -> dRs)
#define fsls_OilGasDataRsx1(oilgasdata)        ((oilgasdata) -> Rsx1)
#define fsls_OilGasDataRsx2(oilgasdata)        ((oilgasdata) -> Rsx2)
#define fsls_OilGasDataRsy1(oilgasdata)        ((oilgasdata) -> Rsy1)
#define fsls_OilGasDataRsy2(oilgasdata)        ((oilgasdata) -> Rsy2)
#define fsls_OilGasDataRsnow(oilgasdata)       ((oilgasdata) -> Rsnow)

#define fsls_OilGasDataBo(oilgasdata)          ((oilgasdata) -> Bo)
#define fsls_OilGasDataDBo(oilgasdata)         ((oilgasdata) -> dBo)
#define fsls_OilGasDataBonow(oilgasdata)       ((oilgasdata) -> Bonow)
#define fsls_OilGasDataBg(oilgasdata)          ((oilgasdata) -> Bg)
#define fsls_OilGasDataDBg(oilgasdata)         ((oilgasdata) -> dBg)
#define fsls_OilGasDataBgnow(oilgasdata)       ((oilgasdata) -> Bgnow)

// Oil and Gas relative permeabilities
#define fsls_OilGasDataKro(oilgasdata)         ((oilgasdata) -> Kro)
#define fsls_OilGasDataDKro(oilgasdata)        ((oilgasdata) -> dKro)
#define fsls_OilGasDataKrg(oilgasdata)         ((oilgasdata) -> Krg)
#define fsls_OilGasDataDKrg(oilgasdata)        ((oilgasdata) -> dKrg)

// Transmissibility and their derivatives   
/* x-direction, Oil */
#define fsls_OilGasDataYox1(oilgasdata)        ((oilgasdata) -> Yox1)
#define fsls_OilGasDataYox2(oilgasdata)        ((oilgasdata) -> Yox2)
#define fsls_OilGasDataDpYox1(oilgasdata)      ((oilgasdata) -> dpYox1)
#define fsls_OilGasDataDpYox2(oilgasdata)      ((oilgasdata) -> dpYox2)
#define fsls_OilGasDataDsYox1(oilgasdata)      ((oilgasdata) -> dsYox1)
#define fsls_OilGasDataDsYox2(oilgasdata)      ((oilgasdata) -> dsYox2)
#define fsls_OilGasDataDdpYox1(oilgasdata)     ((oilgasdata) -> ddpYox1)
#define fsls_OilGasDataDdpYox2(oilgasdata)     ((oilgasdata) -> ddpYox2)
#define fsls_OilGasDataDdsYox1(oilgasdata)     ((oilgasdata) -> ddsYox1)
#define fsls_OilGasDataDdsYox2(oilgasdata)     ((oilgasdata) -> ddsYox2)
/* x-direction, Gas */
#define fsls_OilGasDataYgx1(oilgasdata)        ((oilgasdata) -> Ygx1)
#define fsls_OilGasDataYgx2(oilgasdata)        ((oilgasdata) -> Ygx2)
#define fsls_OilGasDataDpYgx1(oilgasdata)      ((oilgasdata) -> dpYgx1)
#define fsls_OilGasDataDpYgx2(oilgasdata)      ((oilgasdata) -> dpYgx2)
#define fsls_OilGasDataDsYgx1(oilgasdata)      ((oilgasdata) -> dsYgx1)
#define fsls_OilGasDataDsYgx2(oilgasdata)      ((oilgasdata) -> dsYgx2)
#define fsls_OilGasDataDdpYgx1(oilgasdata)     ((oilgasdata) -> ddpYgx1)
#define fsls_OilGasDataDdpYgx2(oilgasdata)     ((oilgasdata) -> ddpYgx2)
#define fsls_OilGasDataDdsYgx1(oilgasdata)     ((oilgasdata) -> ddsYgx1)
#define fsls_OilGasDataDdsYgx2(oilgasdata)     ((oilgasdata) -> ddsYgx2)
/* y-direction, Oil */
#define fsls_OilGasDataYoy1(oilgasdata)        ((oilgasdata) -> Yoy1)
#define fsls_OilGasDataYoy2(oilgasdata)        ((oilgasdata) -> Yoy2)
#define fsls_OilGasDataDpYoy1(oilgasdata)      ((oilgasdata) -> dpYoy1)
#define fsls_OilGasDataDpYoy2(oilgasdata)      ((oilgasdata) -> dpYoy2)
#define fsls_OilGasDataDsYoy1(oilgasdata)      ((oilgasdata) -> dsYoy1)
#define fsls_OilGasDataDsYoy2(oilgasdata)      ((oilgasdata) -> dsYoy2)
#define fsls_OilGasDataDdpYoy1(oilgasdata)     ((oilgasdata) -> ddpYoy1)
#define fsls_OilGasDataDdpYoy2(oilgasdata)     ((oilgasdata) -> ddpYoy2)
#define fsls_OilGasDataDdsYoy1(oilgasdata)     ((oilgasdata) -> ddsYoy1)
#define fsls_OilGasDataDdsYoy2(oilgasdata)     ((oilgasdata) -> ddsYoy2)
/* y-direction, Gas */
#define fsls_OilGasDataYgy1(oilgasdata)        ((oilgasdata) -> Ygy1)
#define fsls_OilGasDataYgy2(oilgasdata)        ((oilgasdata) -> Ygy2)
#define fsls_OilGasDataDpYgy1(oilgasdata)      ((oilgasdata) -> dpYgy1)
#define fsls_OilGasDataDpYgy2(oilgasdata)      ((oilgasdata) -> dpYgy2)
#define fsls_OilGasDataDsYgy1(oilgasdata)      ((oilgasdata) -> dsYgy1)
#define fsls_OilGasDataDsYgy2(oilgasdata)      ((oilgasdata) -> dsYgy2)
#define fsls_OilGasDataDdpYgy1(oilgasdata)     ((oilgasdata) -> ddpYgy1)
#define fsls_OilGasDataDdpYgy2(oilgasdata)     ((oilgasdata) -> ddpYgy2)
#define fsls_OilGasDataDdsYgy1(oilgasdata)     ((oilgasdata) -> ddsYgy1)
#define fsls_OilGasDataDdsYgy2(oilgasdata)     ((oilgasdata) -> ddsYgy2)

#define fsls_OilGasDataPbegin(oilgasdata)      ((oilgasdata) -> Pbegin)
#define fsls_OilGasDataSbegin(oilgasdata)      ((oilgasdata) -> Sbegin)
#define fsls_OilGasDataPnow(oilgasdata)        ((oilgasdata) -> Pnow)
#define fsls_OilGasDataSnow(oilgasdata)        ((oilgasdata) -> Snow)

#define fsls_OilGasDataR(oilgasdata)           ((oilgasdata) -> R)
#define fsls_OilGasDataDelta(oilgasdata)       ((oilgasdata) -> delta)
#define fsls_OilGasDataDeltap(oilgasdata)      ((oilgasdata) -> deltap)
#define fsls_OilGasDataDeltas(oilgasdata)      ((oilgasdata) -> deltas)
#define fsls_OilGasDataJacobi(oilgasdata)      ((oilgasdata) -> Jacobi)

#define fsls_OilGasDataItArray(oilgasdata)     ((oilgasdata) -> It_array)
#define fsls_OilGasDataTsArray(oilgasdata)     ((oilgasdata) -> Ts_array)
#define fsls_OilGasDataTnArray(oilgasdata)     ((oilgasdata) -> Tn_array)
#define fsls_OilGasDataPwArray(oilgasdata)     ((oilgasdata) -> Pw_array)
#define fsls_OilGasDataSwArray(oilgasdata)     ((oilgasdata) -> Sw_array)
#define fsls_OilGasDataORArray(oilgasdata)     ((oilgasdata) -> OR_array)
#define fsls_OilGasDataPArray(oilgasdata)      ((oilgasdata) -> P_array)
#define fsls_OilGasDataSArray(oilgasdata)      ((oilgasdata) -> S_array)

/* oilgas.c */
fsls_OilGasParam *fsls_OilGasParamCreate( );
int fsls_OilGasParamInitialize( fsls_OilGasParam  *oilgasparam, int Nx, int Ny, int Nz );
double fsls_WellIndex( double Kx, double Ky, double Dx, double Dy, double Dz );
int fsls_OilGasParamOutput( fsls_OilGasParam  *oilgasparam );
fsls_OilGasData *fsls_OilGasDataCreate( fsls_OilGasParam *oilgasparam );
int fsls_OilGasDataInitialize( fsls_OilGasData *oilgasdata );
int fsls_JacobiMatNZ( int Nx, int Ny );
int fsls_JacobiMatIA( int Nx, int Ny, int *ia );
int fsls_OilGasDataFinalize( fsls_OilGasData *oilgasdata );
int fsls_InitialPressureSaturation( fsls_OilGasData *oilgasdata );
int fsls_GetNextDt( fsls_OilGasData *oilgasdata, double Tnow );
void fsls_FROMxy2TOxy1( int Nx, int Ny, double *vx2, double *vy2, double *vx1, double *vy1 );
void fsls_LastPressureSaturation( fsls_OilGasData *oilgasdata, char *filename_P, char *filename_S );
void fsls_BackUpTsArray( fsls_OilGasData *oilgasdata, int SucStep );
void fsls_BackUpItArray( fsls_OilGasData *oilgasdata, int SucStep, int numiter );
void fsls_BackUpTnArray( fsls_OilGasData *oilgasdata, int SucStep, double Tnow );
void fsls_BackUpPwArray( fsls_OilGasData *oilgasdata, int SucStep );
void fsls_BackUpSwArray( fsls_OilGasData *oilgasdata, int SucStep );
void fsls_BackUpORArray( fsls_OilGasData *oilgasdata, int SucStep );
void fsls_PostProcessWellInfoStep( fsls_OilGasData *oilgasdata, int SucStep );
void fsls_PostProcessWellInfoTnow( fsls_OilGasData *oilgasdata, int SucStep );
void fsls_PostProcessPSDistribute( fsls_OilGasData *oilgasdata, int step );
void fsls_PostProcessPSLast( fsls_OilGasData *oilgasdata, int step );

/* newton.c */
int fsls_BackUpPSData( fsls_OilGasData *oilgasdata, int SucStep );
double fsls_DelatSMax( fsls_OilGasData *oilgasdata );
double fsls_DelatPMax( fsls_OilGasData *oilgasdata );
int fsls_UpdateWellPressure( fsls_OilGasData *oilgasdata );
int fsls_NewtonIteration( fsls_OilGasData *oilgasdata );

/* iniguess.c */
int fsls_SetNewtonInitialGuess( fsls_OilGasData *oilgasdata );
int fsls_SetNewtonInitialGuessTL( fsls_OilGasData *oilgasdataH, fsls_OilGasData *oilgasdata );
void fsls_BlackOilRestrict( int Nx_coarse, int Ny_coarse, int Nx_fine, double *data_fine, double *data_coarse );
void fsls_BlackOilProlong( int Nx_coarse, int Ny_coarse, int Nx_fine, double *data_coarse, double *data_fine );
void
fsls_BlackOilProlongWeight( int      Nx_coarse, 
                            int      Ny_coarse, 
                            int      Nx_fine, 
                            double  *data_coarse, 
                            double  *data_fine, 
                            double  *data_weight );
/* pdterms.c */
int fsls_ComputeBoBgRs( fsls_OilGasData *oilgasdata, int flag );
int fsls_UpdatePressureDependentTerms( fsls_OilGasData *oilgasdata );

/* sdterms.c */
int fsls_UpdateSaturationDependentTerms( fsls_OilGasData *oilgasdata );

/* trans.c */
int fsls_UpdateTransmissibility( fsls_OilGasData *oilgasdata );

/* residual.c */
int fsls_FormResidual( fsls_OilGasData *oilgasdata ); 

/* jacobi.c */
int fsls_FormJacobiMatrixFull( fsls_OilGasData *oilgasdata );

/* solver.c */
int fsls_JacobiSolveFull( fsls_OilGasData *oilgasdata );

#endif
 
