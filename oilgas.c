/*!
 *   oilgas.c -- basic operations
 *
 *   Created by peghoty  2011/09/09
 *   Xiangtan University
 *   peghoty@163.com
 *
 */

#include "basic.h"
#include "util.h"
#include "blackoil.h"

/*!
 * \fn fsls_OilGasParam *fsls_OilGasParamCreate
 * \brief Create a fsls_OilGasParam object.
 * \return pointer to a fsls_OilGasParam object
 * \author peghoty
 * \date 2011/09/09
 */ 
fsls_OilGasParam *
fsls_OilGasParamCreate( )
{
   fsls_OilGasParam  *oilgasparam = NULL;
   oilgasparam = fsls_CTAlloc(fsls_OilGasParam, 1); 
   return oilgasparam;
}

/*!
 * \fn int fsls_OilGasParamInitialize
 * \brief Initialize a fsls_OilGasParam object.
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_OilGasParamInitialize( fsls_OilGasParam  *oilgasparam, int Nx, int Ny, int Nz )
{
   // parameters for mesh
   double Dx, Dy, Dz, V;  
   double XX, YY, ZZ;
   
   // parameters for porosity and permeability
   double Phi;     // porosity
   double Kx, Ky;  // permeability in X- and Y- direction
   
   // parameters for the Productio Well
   int    Wx, Wy;  // X- and Y- Coordinate of Well Location
   double WI;      // Well Index
   double OR;      // Well Oil Production Rate Constraint (STB/day)
   double BHP;     // Well BHP Constraint (Psia)
   double Pw;      // The current well Pressure
   double Sw;      // The current well Saturation

   // parameters for time control
   double omega;
   double yitap;
   double yitas;
   double dt0;     // The first time step
   double dt;      // The Current time step
   double Tmax;    // The maximum time step
   double Tend;    // The Total simulation time 
    
   // parameters for Newton Iteration
   double tolN;    // tolerance for Residual Norm
   double tolP;    // tolerance for Pressure
   double tolS;    // tolerance for Saturation
   
   // parameters for Number System Conversion
   double alpha;   // British System -> Metric System 

   int    MaxStep; // The maximum simulation step allowed
   int    SucStep; // The number of successful time steps
   int    NCFlag;  // Mark whether the Newton iteration is convergenced
   
   int    postpro; // whether do the post-processing
   
   int    NumUnconvType1; // number of unconvergence type 1 (Rnorm increses)
   int    NumUnconvType2; // number of unconvergence type 1 (reaches the maximal number of iteration)   
  
   XX = 2200.0;
   YY = 1100.0;
   ZZ = 300.0;

   Dx = XX / Nx; 
   Dy = YY / Ny; 
   Dz = ZZ / Nz;
   V  = Dx*Dy*Dz / 5.615;  
   
   Phi = 0.25;     
   Kx  = 10.0;     
   Ky  = 10.0;    

   Wx  = (Nx + 1) / 2;     
   Wy  = (Ny + 1) / 2; 
   OR  = 5000;  
   BHP = 2000; 
   Pw  = 0.0;
   Sw  = 0.0;

   WI = fsls_WellIndex(Kx, Ky, Dx, Dy, Dz);

   omega  = 0.5;
   yitap  = 100.0;
   yitas  = 0.1;
   dt0    = 0.3;
   dt     = dt0;        
   Tmax   = 25.0;    
   Tend   = 1000.0; 
   
   tolN = 0.1;  
   tolP = 0.0001;  
   tolS = 0.005; 
   
   alpha = 0.001127;
   
   MaxStep = 1500;
   SucStep = 0;
   
   NCFlag  = 0;
   postpro = 0;
   
   NumUnconvType1 = 0;
   NumUnconvType2 = 0;   

   fsls_OilGasParamXX(oilgasparam)     = XX;
   fsls_OilGasParamYY(oilgasparam)     = YY;
   fsls_OilGasParamZZ(oilgasparam)     = ZZ;
   fsls_OilGasParamNx(oilgasparam)     = Nx;
   fsls_OilGasParamNy(oilgasparam)     = Ny;
   fsls_OilGasParamNz(oilgasparam)     = Nz;
   fsls_OilGasParamDx(oilgasparam)     = Dx;
   fsls_OilGasParamDy(oilgasparam)     = Dy;
   fsls_OilGasParamDz(oilgasparam)     = Dz;
   fsls_OilGasParamV(oilgasparam)      = V;
   fsls_OilGasParamPhi(oilgasparam)    = Phi;
   fsls_OilGasParamKx(oilgasparam)     = Kx;
   fsls_OilGasParamKy(oilgasparam)     = Ky;

   fsls_OilGasParamWx(oilgasparam)     = Wx;
   fsls_OilGasParamWy(oilgasparam)     = Wy;
   fsls_OilGasParamWI(oilgasparam)     = WI;
   fsls_OilGasParamOR(oilgasparam)     = OR;
   fsls_OilGasParamBHP(oilgasparam)    = BHP;
   fsls_OilGasParamPw(oilgasparam)     = Pw;
   fsls_OilGasParamSw(oilgasparam)     = Sw;

   fsls_OilGasParamOmega(oilgasparam)  = omega;
   fsls_OilGasParamYitap(oilgasparam)  = yitap;
   fsls_OilGasParamYitas(oilgasparam)  = yitas;   
   fsls_OilGasParamDt0(oilgasparam)    = dt0;
   fsls_OilGasParamDt(oilgasparam)     = dt;
   fsls_OilGasParamTmax(oilgasparam)   = Tmax;
   fsls_OilGasParamTend(oilgasparam)   = Tend;
   
   fsls_OilGasParamTolN(oilgasparam)   = tolN;
   fsls_OilGasParamTolP(oilgasparam)   = tolP;
   fsls_OilGasParamTolS(oilgasparam)   = tolS;
   
   fsls_OilGasParamAlpha(oilgasparam)  = alpha;
   
   fsls_OilGasParamMaxStep(oilgasparam) = MaxStep; 
   fsls_OilGasParamSucStep(oilgasparam) = SucStep;
   fsls_OilGasParamNCFlag(oilgasparam)  = NCFlag;
   fsls_OilGasParamPostpro(oilgasparam) = postpro;
   
   fsls_OilGasParamNumUnconvType1(oilgasparam) = NumUnconvType1;
   fsls_OilGasParamNumUnconvType2(oilgasparam) = NumUnconvType2; 
   
   return 0;  
}

/*!
 * \fn double fsls_WellIndex
 * \brief Compute the Well Index.
 * \author peghoty
 * \date 2011/09/09
 */ 
double 
fsls_WellIndex( double Kx, double Ky, double Dx, double Dy, double Dz )
{
   double kxy   = Kx / Ky;
   double kyx   = Ky / Kx;
   double sqkxy = sqrt(kxy);
   double sqkyx = sqrt(kyx);
   double beta1 = sqrt( sqkyx*Dx*Dx + sqkxy*Dy*Dy );
   double beta2 = sqrt(sqkxy) + sqrt(sqkxy);
   double ro    = 0.28 * beta1 / beta2;
   double rw    = 0.25;
   double s     = 0.0;
   double WI;
   
   beta1 = 2 * PI * sqrt(Kx*Ky) * Dz;
   beta2 = log(ro/rw) + s;
   WI    = beta1 / beta2; 
   
   return (WI);  
}

/*!
 * \fn int fsls_OilGasParamOutput
 * \brief Print the members of a fsls_OilGasParam object.
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_OilGasParamOutput( fsls_OilGasParam  *oilgasparam )
{
   // parameters for mesh
   int    Nx, Ny, Nz;     // number of grids in X-, Y-, and Z- direction
   double Dx, Dy, Dz, V;  // grid size in X-, Y-, and Z- direction; Volumn for each grid
   double XX, YY, ZZ;     // solving domain: [0,XX]x[0,YY]X[0,ZZ]
   
   // parameters for porosity and permeability
   double Phi;            // porosity
   double Kx, Ky;         // permeability in X- and Y- direction
   
   // parameters for the Productio Well
   int    Wx, Wy;         // X- and Y- Coordinate of Well Location
   double WI;             // Well Index
   double OR;             // Well Oil Production Rate Constraint (STB/day)
   double BHP;            // Well BHP Constraint (Psia)

   // parameters for time control
   double omega;
   double yitap;
   double yitas;
   double dt0;            // The first time step
   double Tmax;           // The maximum time step
   double Tend;           // The Total simulation time 
   
   // parameters for Newton Iteration
   double tolN;           // tolerance for Residual Norm
   double tolP;           // tolerance for Pressure
   double tolS;           // tolerance for Saturation
   
   // parameters for Number System Conversion
   double alpha;          // British System -> Metric System 
   
   int    MaxStep;        // The maximum simulation step allowed
   
   int    postpro;        // whether do the post-processing   
   
   XX  = fsls_OilGasParamXX(oilgasparam);
   YY  = fsls_OilGasParamYY(oilgasparam);
   ZZ  = fsls_OilGasParamZZ(oilgasparam);
   
   Nx  = fsls_OilGasParamNx(oilgasparam);
   Ny  = fsls_OilGasParamNy(oilgasparam);
   Nz  = fsls_OilGasParamNz(oilgasparam);
   Dx  = fsls_OilGasParamDx(oilgasparam);
   Dy  = fsls_OilGasParamDy(oilgasparam);
   Dz  = fsls_OilGasParamDz(oilgasparam);
   V   = fsls_OilGasParamV(oilgasparam);
   Phi = fsls_OilGasParamPhi(oilgasparam);
   Kx  = fsls_OilGasParamKx(oilgasparam);
   Ky  = fsls_OilGasParamKy(oilgasparam);

   Wx = fsls_OilGasParamWx(oilgasparam);
   Wy = fsls_OilGasParamWy(oilgasparam);
   WI = fsls_OilGasParamWI(oilgasparam);
   OR = fsls_OilGasParamOR(oilgasparam);
   BHP = fsls_OilGasParamBHP(oilgasparam);

   omega = fsls_OilGasParamOmega(oilgasparam);
   yitap = fsls_OilGasParamYitap(oilgasparam);
   yitas = fsls_OilGasParamYitas(oilgasparam);   
   dt0   = fsls_OilGasParamDt0(oilgasparam);
   Tmax  = fsls_OilGasParamTmax(oilgasparam);
   Tend  = fsls_OilGasParamTend(oilgasparam);
   tolN  = fsls_OilGasParamTolN(oilgasparam);
   tolP  = fsls_OilGasParamTolP(oilgasparam);
   tolS  = fsls_OilGasParamTolS(oilgasparam);
   
   alpha   = fsls_OilGasParamAlpha(oilgasparam);
   MaxStep = fsls_OilGasParamMaxStep(oilgasparam);
   postpro = fsls_OilGasParamPostpro(oilgasparam); 
   
   printf("\n\n");
   printf(" ======================================================================\n");
   printf("           Information for Some Members of OilGasParam Object          \n");
   printf(" ======================================================================\n");

   printf(" solving domain                  Omega        = (0,%.2lf)X(0,%.2lf)X(0,%.2lf)\n", XX, YY, ZZ);
   printf(" number of grids                 (Nx, Ny, Nz) = (%d %d %d)\n", Nx, Ny, Nz);
   printf(" grid size                       (Dx, Dy, Dz) = (%.2lf %.2lf %.2lf)\n", Dx, Dy, Dz);
   printf(" Volumn for each grid                         = %.15le\n", V);
   printf(" porosity                                     = %le\n", Phi);
   printf(" permeability                        (Kx, Ky) = (%.2lf, %.2lf)\n", Kx, Ky);
   printf(" Well Location                       (Wx, Wy) = (%d, %d)\n", Wx, Wy);
   printf(" Well Index                                   = %.15le\n", WI);
   printf(" Well Oil Production Rate (STB/day)        OR = %le\n", OR);
   printf(" Well BHP Constraint (Psia)               BHP = %le\n", BHP);
   printf(" omega                                        = %le\n", omega);
   printf(" yitap                                        = %le\n", yitap);
   printf(" yitas                                        = %le\n", yitas);
   printf(" British System -> Metric System              = %le\n", alpha);
   printf(" The maximum simulation step allowed          = %d\n",  MaxStep);
   printf(" The first time step                          = %le\n", dt0);
   printf(" The maximum time step                        = %le\n", Tmax);
   printf(" The Total simulation time                    = %le\n", Tend);
   printf(" tolerance for Residual Norm                  = %le\n", tolN);
   printf(" tolerance for Pressure                       = %le\n", tolP);
   printf(" tolerance for Saturation                     = %le\n", tolS); 
   printf(" Marker for post-processing                   = %d\n",  postpro);
   printf("\n"); 
      
   return 0;  
}

/*!
 * \fn fsls_OilGasData *fsls_OilGasDataCreate
 * \brief Create a fsls_OilGasData object.
 * \author peghoty
 * \date 2011/09/09
 */ 
fsls_OilGasData *
fsls_OilGasDataCreate( fsls_OilGasParam *oilgasparam )
{ 
   fsls_OilGasData *oilgasdata = NULL; 
   oilgasdata = fsls_CTAlloc(fsls_OilGasData, 1);
   fsls_OilGasDataOilgasparam(oilgasdata) = oilgasparam; 
   return oilgasdata;
}

/*!
 * \fn int fsls_OilGasDataInitialize
 * \brief Initialize a fsls_OilGasData object.
 *        Allocate memories for all working arrays.
 * \note Which will be more efficient, one by one, or once for all?
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_OilGasDataInitialize( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam);
   int MaxStep = fsls_OilGasParamMaxStep(oilgasparam);
   int postpro = fsls_OilGasParamPostpro(oilgasparam);
   
   int i;
   int N  = Nx*Ny;
   int N2 = 2*N;
   
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
   
   /* x-direction, Oil */
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
   
   /* x-direction, Gas */
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

   /* y-direction, Oil */
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
   
   /* y-direction, Gas */
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
   
   /* Jacobian Matrix */
   double *Jacobi = NULL;

   int     *It_array;  // number of Newton iteration array
   double  *Ts_array;  // Successful time Step array
   double  *Tn_array;  // T_now array
   double  *Pw_array;  // Well Pressure for each successful time step
   double  *Sw_array;  // Well Saturation for each successful time step 
   double  *OR_array;  // Well Oil Production Rate for each successful time step 
   double **P_array;   // The Pressure   sequence for all the successful iteration 
   double **S_array;   // The Saturation sequence for all the successful iteration
   
   
   //==================================================//
   //   Allocate Memories for all double type arrays   //
   //==================================================//  
   
   Muo     = fsls_CTAlloc(double, N);
   dMuo    = fsls_CTAlloc(double, N);
   Mug     = fsls_CTAlloc(double, N);
   dMug    = fsls_CTAlloc(double, N);
   Rs      = fsls_CTAlloc(double, N);
   dRs     = fsls_CTAlloc(double, N);
   Rsx1    = fsls_CTAlloc(double, N);
   Rsx2    = fsls_CTAlloc(double, N);
   Rsy1    = fsls_CTAlloc(double, N);
   Rsy2    = fsls_CTAlloc(double, N);
   Rsnow   = fsls_CTAlloc(double, N);

   Bo      = fsls_CTAlloc(double, N);
   dBo     = fsls_CTAlloc(double, N);
   Bonow   = fsls_CTAlloc(double, N);
   Bg      = fsls_CTAlloc(double, N);
   dBg     = fsls_CTAlloc(double, N);
   Bgnow   = fsls_CTAlloc(double, N);
   
   Kro     = fsls_CTAlloc(double, N);
   dKro    = fsls_CTAlloc(double, N);
   Krg     = fsls_CTAlloc(double, N);
   dKrg    = fsls_CTAlloc(double, N);

   Yox1    = fsls_CTAlloc(double, N);
   Yox2    = fsls_CTAlloc(double, N);
   dpYox1  = fsls_CTAlloc(double, N);
   dpYox2  = fsls_CTAlloc(double, N);
   dsYox1  = fsls_CTAlloc(double, N);
   dsYox2  = fsls_CTAlloc(double, N);
   ddpYox1 = fsls_CTAlloc(double, N);
   ddpYox2 = fsls_CTAlloc(double, N);
   ddsYox1 = fsls_CTAlloc(double, N);
   ddsYox2 = fsls_CTAlloc(double, N);

   Ygx1    = fsls_CTAlloc(double, N);
   Ygx2    = fsls_CTAlloc(double, N);
   dpYgx1  = fsls_CTAlloc(double, N);
   dpYgx2  = fsls_CTAlloc(double, N);
   dsYgx1  = fsls_CTAlloc(double, N);
   dsYgx2  = fsls_CTAlloc(double, N);
   ddpYgx1 = fsls_CTAlloc(double, N);
   ddpYgx2 = fsls_CTAlloc(double, N);
   ddsYgx1 = fsls_CTAlloc(double, N);
   ddsYgx2 = fsls_CTAlloc(double, N);
   
   Yoy1    = fsls_CTAlloc(double, N);
   Yoy2    = fsls_CTAlloc(double, N);
   dpYoy1  = fsls_CTAlloc(double, N);
   dpYoy2  = fsls_CTAlloc(double, N);
   dsYoy1  = fsls_CTAlloc(double, N);
   dsYoy2  = fsls_CTAlloc(double, N);
   ddpYoy1 = fsls_CTAlloc(double, N);
   ddpYoy2 = fsls_CTAlloc(double, N);
   ddsYoy1 = fsls_CTAlloc(double, N);
   ddsYoy2 = fsls_CTAlloc(double, N);
     
   Ygy1    = fsls_CTAlloc(double, N);
   Ygy2    = fsls_CTAlloc(double, N);
   dpYgy1  = fsls_CTAlloc(double, N);
   dpYgy2  = fsls_CTAlloc(double, N);
   dsYgy1  = fsls_CTAlloc(double, N);
   dsYgy2  = fsls_CTAlloc(double, N);
   ddpYgy1 = fsls_CTAlloc(double, N);
   ddpYgy2 = fsls_CTAlloc(double, N);
   ddsYgy1 = fsls_CTAlloc(double, N);
   ddsYgy2 = fsls_CTAlloc(double, N);
      
   Pbegin   = fsls_CTAlloc(double, N);
   Sbegin   = fsls_CTAlloc(double, N);
   Pnow     = fsls_CTAlloc(double, N);
   Snow     = fsls_CTAlloc(double, N);
   R        = fsls_CTAlloc(double, N2);
   delta    = fsls_CTAlloc(double, N2);
   deltap   = fsls_CTAlloc(double, N);
   deltas   = fsls_CTAlloc(double, N);
   
   /* coefficient matrix */
   Jacobi = fsls_CTAlloc(double, N2*N2);
   fsls_OilGasDataJacobi(oilgasdata) = Jacobi;
   
   It_array = fsls_CTAlloc(int,    MaxStep+1);
   Ts_array = fsls_CTAlloc(double, MaxStep+1);
   Tn_array = fsls_CTAlloc(double, MaxStep+1);
   Pw_array = fsls_CTAlloc(double, MaxStep+1);
   Sw_array = fsls_CTAlloc(double, MaxStep+1); 
   OR_array = fsls_CTAlloc(double, MaxStep+1);  
   
   if (postpro == 1)
   {
      P_array = fsls_CTAlloc(double *, MaxStep);
      S_array = fsls_CTAlloc(double *, MaxStep);
   
      for (i = 0; i < MaxStep; i ++)
      {
         P_array[i] = fsls_CTAlloc(double, N);
         S_array[i] = fsls_CTAlloc(double, N);  
      } 

      fsls_OilGasDataPArray(oilgasdata) = P_array;
      fsls_OilGasDataSArray(oilgasdata) = S_array;   
   }


   //============================================//
   //       Fill the arrays to oilgasdata        //
   //============================================// 
   
   fsls_OilGasDataMuo(oilgasdata)      = Muo;
   fsls_OilGasDataDMuo(oilgasdata)     = dMuo;
   fsls_OilGasDataMug(oilgasdata)      = Mug;
   fsls_OilGasDataDMug(oilgasdata)     = dMug;
   
   fsls_OilGasDataRs(oilgasdata)       = Rs;
   fsls_OilGasDataDRs(oilgasdata)      = dRs;
   fsls_OilGasDataRsx1(oilgasdata)     = Rsx1;
   fsls_OilGasDataRsx2(oilgasdata)     = Rsx2;
   fsls_OilGasDataRsy1(oilgasdata)     = Rsy1;
   fsls_OilGasDataRsy2(oilgasdata)     = Rsy2;
   fsls_OilGasDataRsnow(oilgasdata)    = Rsnow;

   fsls_OilGasDataBo(oilgasdata)       = Bo;
   fsls_OilGasDataDBo(oilgasdata)      = dBo;
   fsls_OilGasDataBonow(oilgasdata)    = Bonow;
   fsls_OilGasDataBg(oilgasdata)       = Bg;
   fsls_OilGasDataDBg(oilgasdata)      = dBg;
   fsls_OilGasDataBgnow(oilgasdata)    = Bgnow;  
   
   // Oil and Gas relative permeabilities
   fsls_OilGasDataKro(oilgasdata)      = Kro;
   fsls_OilGasDataDKro(oilgasdata)     = dKro;
   fsls_OilGasDataKrg(oilgasdata)      = Krg;
   fsls_OilGasDataDKrg(oilgasdata)     = dKrg;
   
   // Transmissibility and their derivatives   
   /* x-direction, Oil */
   fsls_OilGasDataYox1(oilgasdata)     = Yox1;
   fsls_OilGasDataYox2(oilgasdata)     = Yox2;
   fsls_OilGasDataDpYox1(oilgasdata)   = dpYox1;
   fsls_OilGasDataDpYox2(oilgasdata)   = dpYox2;
   fsls_OilGasDataDsYox1(oilgasdata)   = dsYox1;
   fsls_OilGasDataDsYox2(oilgasdata)   = dsYox2;
   fsls_OilGasDataDdpYox1(oilgasdata)  = ddpYox1;
   fsls_OilGasDataDdpYox2(oilgasdata)  = ddpYox2;
   fsls_OilGasDataDdsYox1(oilgasdata)  = ddsYox1;
   fsls_OilGasDataDdsYox2(oilgasdata)  = ddsYox2;
   /* x-direction, Gas */
   fsls_OilGasDataYgx1(oilgasdata)     = Ygx1;
   fsls_OilGasDataYgx2(oilgasdata)     = Ygx2;
   fsls_OilGasDataDpYgx1(oilgasdata)   = dpYgx1;
   fsls_OilGasDataDpYgx2(oilgasdata)   = dpYgx2;
   fsls_OilGasDataDsYgx1(oilgasdata)   = dsYgx1;
   fsls_OilGasDataDsYgx2(oilgasdata)   = dsYgx2;
   fsls_OilGasDataDdpYgx1(oilgasdata)  = ddpYgx1;
   fsls_OilGasDataDdpYgx2(oilgasdata)  = ddpYgx2;
   fsls_OilGasDataDdsYgx1(oilgasdata)  = ddsYgx1;
   fsls_OilGasDataDdsYgx2(oilgasdata)  = ddsYgx2;
   /* y-direction, Oil */
   fsls_OilGasDataYoy1(oilgasdata)     = Yoy1;
   fsls_OilGasDataYoy2(oilgasdata)     = Yoy2;
   fsls_OilGasDataDpYoy1(oilgasdata)   = dpYoy1;
   fsls_OilGasDataDpYoy2(oilgasdata)   = dpYoy2;
   fsls_OilGasDataDsYoy1(oilgasdata)   = dsYoy1;
   fsls_OilGasDataDsYoy2(oilgasdata)   = dsYoy2;
   fsls_OilGasDataDdpYoy1(oilgasdata)  = ddpYoy1;
   fsls_OilGasDataDdpYoy2(oilgasdata)  = ddpYoy2;
   fsls_OilGasDataDdsYoy1(oilgasdata)  = ddsYoy1;
   fsls_OilGasDataDdsYoy2(oilgasdata)  = ddsYoy2;
   /* y-direction, Gas */
   fsls_OilGasDataYgy1(oilgasdata)     = Ygy1;
   fsls_OilGasDataYgy2(oilgasdata)     = Ygy2;
   fsls_OilGasDataDpYgy1(oilgasdata)   = dpYgy1;
   fsls_OilGasDataDpYgy2(oilgasdata)   = dpYgy2;
   fsls_OilGasDataDsYgy1(oilgasdata)   = dsYgy1;
   fsls_OilGasDataDsYgy2(oilgasdata)   = dsYgy2;
   fsls_OilGasDataDdpYgy1(oilgasdata)  = ddpYgy1;
   fsls_OilGasDataDdpYgy2(oilgasdata)  = ddpYgy2;
   fsls_OilGasDataDdsYgy1(oilgasdata)  = ddsYgy1;
   fsls_OilGasDataDdsYgy2(oilgasdata)  = ddsYgy2;
   
   fsls_OilGasDataPbegin(oilgasdata)   = Pbegin;
   fsls_OilGasDataSbegin(oilgasdata)   = Sbegin;
   fsls_OilGasDataPnow(oilgasdata)     = Pnow;
   fsls_OilGasDataSnow(oilgasdata)     = Snow;  
   
   fsls_OilGasDataR(oilgasdata)        = R;
   fsls_OilGasDataDelta(oilgasdata)    = delta;
   fsls_OilGasDataDeltap(oilgasdata)   = deltap;
   fsls_OilGasDataDeltas(oilgasdata)   = deltas;

   fsls_OilGasDataItArray(oilgasdata)  = It_array;
   fsls_OilGasDataTsArray(oilgasdata)  = Ts_array;
   fsls_OilGasDataTnArray(oilgasdata)  = Tn_array;
   fsls_OilGasDataPwArray(oilgasdata)  = Pw_array;
   fsls_OilGasDataSwArray(oilgasdata)  = Sw_array;
   fsls_OilGasDataORArray(oilgasdata)  = OR_array;
        
   return 0;
}


/*!
 * \fn int fsls_OilGasDataFinalize
 * \brief Finalize a fsls_OilGasData object.
 * \author peghoty
 * \date 2011/09/09
 */
int
fsls_OilGasDataFinalize( fsls_OilGasData *oilgasdata )
{
   fsls_OilGasParam *oilgasparam = NULL;
   int MaxStep;
   int postpro;
   int i;
   
   if (oilgasdata)
   {
      oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
      MaxStep = fsls_OilGasParamMaxStep(oilgasparam);
      postpro = fsls_OilGasParamPostpro(oilgasparam);
      
      /* Muo */
      if ( fsls_OilGasDataMuo(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataMuo(oilgasdata) );
      }
      
      /* dMuo */
      if ( fsls_OilGasDataDMuo(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDMuo(oilgasdata) );
      }

      /* Mug */
      if ( fsls_OilGasDataMug(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataMug(oilgasdata) );
      }

      /* dMug */
      if ( fsls_OilGasDataDMug(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDMug(oilgasdata) );
      }      

      /* Rs */
      if ( fsls_OilGasDataRs(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataRs(oilgasdata) );
      }

      /* dRs */
      if ( fsls_OilGasDataDRs(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDRs(oilgasdata) );
      }

      /* Rsx1 */
      if ( fsls_OilGasDataRsx1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataRsx1(oilgasdata) );
      }  

      /* Rsx2 */ 
      if ( fsls_OilGasDataRsx2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataRsx2(oilgasdata) );
      }

      /* Rsy1 */
      if ( fsls_OilGasDataRsy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataRsy1(oilgasdata) );
      }

      /* Rsy2 */      
      if ( fsls_OilGasDataRsy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataRsy2(oilgasdata) );
      }  

      /* Rsnow */ 
      if ( fsls_OilGasDataRsnow(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataRsnow(oilgasdata) );
      }

      /* Bo */
      if ( fsls_OilGasDataBo(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataBo(oilgasdata) );
      }

      /* dBo */      
      if ( fsls_OilGasDataDBo(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDBo(oilgasdata) );
      }  

      /* Bonow */ 
      if ( fsls_OilGasDataBonow(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataBonow(oilgasdata) );
      }

      /* Bg */
      if ( fsls_OilGasDataBg(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataBg(oilgasdata) );
      }

      /* dBg */      
      if ( fsls_OilGasDataDBg(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDBg(oilgasdata) );
      }  

      /* Bgnow */ 
      if ( fsls_OilGasDataBgnow(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataBgnow(oilgasdata) );
      }

      /* Kro */
      if ( fsls_OilGasDataKro(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataKro(oilgasdata) );
      }

      /* dKro */      
      if ( fsls_OilGasDataDKro(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDKro(oilgasdata) );
      }      

      /* Krg */      
      if ( fsls_OilGasDataKrg(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataKrg(oilgasdata) );
      }  

      /* dKrg */
      if ( fsls_OilGasDataDKrg(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDKrg(oilgasdata) );
      }

      /* Yox1 */      
      if ( fsls_OilGasDataYox1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataYox1(oilgasdata) );
      }  

      /* Yox2 */
      if ( fsls_OilGasDataYox2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataYox2(oilgasdata) );
      }

      /* dpYox1 */      
      if ( fsls_OilGasDataDpYox1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDpYox1(oilgasdata) );
      }  

      /* dpYox2 */
      if ( fsls_OilGasDataDpYox2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDpYox2(oilgasdata) );
      }

      /* dsYox1 */      
      if ( fsls_OilGasDataDsYox1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDsYox1(oilgasdata) );
      }  

      /* dsYox2 */
      if ( fsls_OilGasDataDsYox2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDsYox2(oilgasdata) );
      }

      /* ddpYox1 */      
      if ( fsls_OilGasDataDdpYox1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdpYox1(oilgasdata) );
      }  

      /* ddpYox2 */
      if ( fsls_OilGasDataDdpYox2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdpYox2(oilgasdata) );
      }

      /* ddsYox1 */      
      if ( fsls_OilGasDataDdsYox1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdsYox1(oilgasdata) );
      }  

      /* ddsYox2 */
      if ( fsls_OilGasDataDdsYox2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdsYox2(oilgasdata) );
      }

      /* Ygx1 */      
      if ( fsls_OilGasDataYgx1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataYgx1(oilgasdata) );
      }  

      /* Ygx2 */
      if ( fsls_OilGasDataYgx2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataYgx2(oilgasdata) );
      }
      
      /* dpYgx1 */
      if ( fsls_OilGasDataDpYgx1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDpYgx1(oilgasdata) );
      }  

      /* dpYgx2 */
      if ( fsls_OilGasDataDpYgx2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDpYgx2(oilgasdata) );
      }

      /* dsYgx1 */      
      if ( fsls_OilGasDataDsYgx1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDsYgx1(oilgasdata) );
      }  

      /* dsYgx2 */
      if ( fsls_OilGasDataDsYgx2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDsYgx2(oilgasdata) );
      }

      /* ddpYgx1 */      
      if ( fsls_OilGasDataDdpYgx1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdpYgx1(oilgasdata) );
      }  

      /* ddpYgx2 */
      if ( fsls_OilGasDataDdpYgx2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdpYgx2(oilgasdata) );
      }

      /* ddsYgx1 */      
      if ( fsls_OilGasDataDdsYgx1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdsYgx1(oilgasdata) );
      }  

      /* ddsYgx2 */
      if ( fsls_OilGasDataDdsYgx2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdsYgx2(oilgasdata) );
      }

      /* Yoy1 */      
      if ( fsls_OilGasDataYoy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataYoy1(oilgasdata) );
      }  

      /* Yoy2 */
      if ( fsls_OilGasDataYoy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataYoy2(oilgasdata) );
      }

      /* dpYoy1 */      
      if ( fsls_OilGasDataDpYoy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDpYoy1(oilgasdata) );
      }  

      /* dpYoy2 */
      if ( fsls_OilGasDataDpYoy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDpYoy2(oilgasdata) );
      }

      /* dsYoy1 */      
      if ( fsls_OilGasDataDsYoy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDsYoy1(oilgasdata) );
      }  

      /* dsYoy2 */
      if ( fsls_OilGasDataDsYoy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDsYoy2(oilgasdata) );
      }
     
      /* ddpYoy1 */ 
      if ( fsls_OilGasDataDdpYoy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdpYoy1(oilgasdata) );
      }  

      /* ddpYoy2 */
      if ( fsls_OilGasDataDdpYoy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdpYoy2(oilgasdata) );
      }

      /* Yoy1 */      
      if ( fsls_OilGasDataDdsYoy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdsYoy1(oilgasdata) );
      }  

      /* Yoy2 */
      if ( fsls_OilGasDataDdsYoy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdsYoy2(oilgasdata) );
      }

      /* Ygy1 */      
      if ( fsls_OilGasDataYgy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataYgy1(oilgasdata) );
      }  

      /* Ygy2 */
      if ( fsls_OilGasDataYgy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataYgy2(oilgasdata) );
      }

      /* dpYgy1 */      
      if ( fsls_OilGasDataDpYgy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDpYgy1(oilgasdata) );
      }  

      /* dpYgy2 */
      if ( fsls_OilGasDataDpYgy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDpYgy2(oilgasdata) );
      }  

      /* dsYgy1 */
      if ( fsls_OilGasDataDsYgy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDsYgy1(oilgasdata) );
      }

      /* dsYgy2 */  
      if ( fsls_OilGasDataDsYgy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDsYgy2(oilgasdata) );
      }  

      /* ddpYgy1 */
      if ( fsls_OilGasDataDdpYgy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdpYgy1(oilgasdata) );
      }  

      /* ddpYgy2 */
      if ( fsls_OilGasDataDdpYgy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdpYgy2(oilgasdata) );
      }

      /* ddsYgy1 */      
      if ( fsls_OilGasDataDdsYgy1(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdsYgy1(oilgasdata) );
      }  

      /* ddsYgy2 */
      if ( fsls_OilGasDataDdsYgy2(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDdsYgy2(oilgasdata) );
      }  

      /* Pbegin */
      if ( fsls_OilGasDataPbegin(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataPbegin(oilgasdata) );
      }

      /* Sbegin */      
      if ( fsls_OilGasDataSbegin(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataSbegin(oilgasdata) );
      }  

      /* Pnow */
      if ( fsls_OilGasDataPnow(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataPnow(oilgasdata) );
      }  

      /* Snow */
      if ( fsls_OilGasDataSnow(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataSnow(oilgasdata) );
      }

      /* R */      
      if ( fsls_OilGasDataR(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataR(oilgasdata) );
      }  

      /* Delta */
      if ( fsls_OilGasDataDelta(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDelta(oilgasdata) );
      }  

      /* Deltap */
      if ( fsls_OilGasDataDeltap(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDeltap(oilgasdata) );
      }  

      /* Deltas */      
      if ( fsls_OilGasDataDeltas(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataDeltas(oilgasdata) );
      }  

      /* Jacobi */      
      if ( fsls_OilGasDataJacobi(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataJacobi(oilgasdata) );
      }  

      /* It_array */      
      if ( fsls_OilGasDataItArray(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataItArray(oilgasdata) );
      }

      /* Ts_array */
      if ( fsls_OilGasDataTsArray(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataTsArray(oilgasdata) );
      }

      /* Tn_array */      
      if ( fsls_OilGasDataTnArray(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataTnArray(oilgasdata) );
      }          

      /* Pw_array */
      if ( fsls_OilGasDataPwArray(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataPwArray(oilgasdata) );
      }  

      /* Sw_array */
      if ( fsls_OilGasDataSwArray(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataSwArray(oilgasdata) );
      }  
      
      /* OR_array */      
      if ( fsls_OilGasDataORArray(oilgasdata) )
      {
         fsls_TFree( fsls_OilGasDataORArray(oilgasdata) );
      }

      /* P_array and S_array */
      if (postpro)
      {
         for (i = 0; i < MaxStep; i ++)
         {
            /* P_array */
            if ( fsls_OilGasDataPArray(oilgasdata)[i] )
            {
               fsls_TFree( fsls_OilGasDataPArray(oilgasdata)[i] );
            } 
            /* S_array */
            if ( fsls_OilGasDataSArray(oilgasdata)[i] )
            {
               fsls_TFree( fsls_OilGasDataSArray(oilgasdata)[i] );
            }
         }
         fsls_TFree( fsls_OilGasDataPArray(oilgasdata) );
         fsls_TFree( fsls_OilGasDataSArray(oilgasdata) );
      } 
      fsls_TFree(oilgasparam);
      fsls_TFree(oilgasdata);
   }

   return 0;
}

/*!
 * \fn int fsls_InitialPressureSaturation
 * \brief Initialize the Pressure and Saturation.
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_InitialPressureSaturation( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   int    Nx = fsls_OilGasParamNx(oilgasparam);
   int    Ny = fsls_OilGasParamNy(oilgasparam);
   double Dz = fsls_OilGasParamDz(oilgasparam);
   int    postpro = fsls_OilGasParamPostpro(oilgasparam);
   
   double *Pbegin = fsls_OilGasDataPbegin(oilgasdata); 
   double *Sbegin = fsls_OilGasDataSbegin(oilgasdata); 

   double **P_array = NULL;
   double **S_array = NULL;   
   
   double Ptop    = 5000.0;  // Reservoir top pressure (Psia)
   double P0      = Ptop;
   double rho_oil = 46.244;
   double rho_gas = 0.0647;
   double beta    = 1.0 / 144.0;
   double gas_co  = 178.0934;
   double Bo, Gammao;
   double a = -3.922264e-8;
   double b = 3.0e-4;
   double c = 1.0;
   int    n = (int) (Dz / 2);
   int    Nxy = Nx*Ny;
   int    i; 
   
   for (i = 0; i < n; i ++)
   {
      Bo = a*P0*P0 + b*P0 + c;
      Gammao = beta*(rho_oil + gas_co*rho_gas) / Bo; // Bo is always positive
      P0 += Gammao; 
   }
   fsls_OilGasParamPw(oilgasparam) = P0;  // Set oilgasparam -> Pw = P0
   fsls_OilGasParamSw(oilgasparam) = 0.0; // Set oilgasparam -> Sw = 0.0
   
   for (i = 0; i < Nxy; i ++)
   {
      Pbegin[i] = P0;
      Sbegin[i] = 0.0;
   }
   
   if (postpro == 1)
   {
      P_array = fsls_OilGasDataPArray(oilgasdata);
      S_array = fsls_OilGasDataSArray(oilgasdata); 
      memcpy(P_array[0], Pbegin, Nxy*sizeof(double));
      memcpy(S_array[0], Sbegin, Nxy*sizeof(double));
   }

   return 0;
}

/*!
 * \fn int fsls_GetNextDt
 * \brief Get the time step for next Newton Iteration.
 * \param oilgasdata pointer to the fsls_OilGasData object
 * \param Tnow the current time
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_GetNextDt( fsls_OilGasData *oilgasdata, double Tnow ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   int    NCFlag = fsls_OilGasParamNCFlag(oilgasparam);
   int    Nx     = fsls_OilGasParamNx(oilgasparam);
   int    Ny     = fsls_OilGasParamNy(oilgasparam);

   double dt     = fsls_OilGasParamDt(oilgasparam);
   double omega  = fsls_OilGasParamOmega(oilgasparam);   
   double yitap  = fsls_OilGasParamYitap(oilgasparam);
   double yitas  = fsls_OilGasParamYitas(oilgasparam);     
   double Tmax   = fsls_OilGasParamTmax(oilgasparam);
   double Tend   = fsls_OilGasParamTend(oilgasparam);
   
   double *deltap = fsls_OilGasDataDeltap(oilgasdata);
   double *deltas = fsls_OilGasDataDeltas(oilgasdata);

   /* auxiliary variables */
   int     Nxy    = Nx*Ny;  
   double  psmax  = 0.0;
   double  temp   = 0.0;
   int     i;

   if (NCFlag == 0)
   {
      dt = 0.5*dt;  // Cut time step by half
   }
   else
   {
      // find max { deltap[i] / yitap, deltas[i] / yitas }
      for (i = 0; i < Nxy; i ++)
      {
         temp = deltap[i] / yitap;
         if (temp > psmax) psmax = temp;
      }
      for (i = 0; i < Nxy; i ++)
      {
         temp = deltas[i] / yitas;
         if (temp > psmax) psmax = temp;
      }
      
      // Auto timestep control
      dt = dt * (1.0 + omega) / (psmax + omega);      

      if (dt > Tmax) 
      {
         dt = Tmax;
      }
      
      if (Tend - Tnow < dt)
      {
         dt = Tend - Tnow;
      }     
   }
   
   fsls_OilGasParamDt(oilgasparam) = dt;
      
   return 0;
}

/*!
 * \fn void fsls_FROMxy2TOxy1
 * \brief Generate vx1 and vy1 by displace the data in vx2 and vy2.
 * \author peghoty
 * \date 2011/09/09
 */ 
void 
fsls_FROMxy2TOxy1( int Nx, int Ny, double *vx2, double *vy2, double *vx1, double *vy1 )
{
   int i,j,k,size,jNx;
   
   for (j = 0; j < Ny; j ++)
   {
      jNx = j*Nx;
      for (i = 0; i < Nx; i ++)
      {
         k = jNx + i;
         if (i == 0)
         {
            vx1[k] = 0.0;
         }
         else
         {
            vx1[k] = vx2[k-1];
         }
      }
   }   
   size = Nx*(Ny-1);
   memcpy(vy1+Nx, vy2, size*sizeof(double));
}

/*!
 * \fn void fsls_BackUpItArray
 * \brief Back up the SucStep-th number of Newton iteration.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_BackUpItArray( fsls_OilGasData *oilgasdata, int SucStep, int numiter )
{
   int *It_array  = fsls_OilGasDataItArray(oilgasdata);
   It_array[SucStep] = numiter;
}

/*!
 * \fn void fsls_BackUpTsArray
 * \brief Back up the SucStep-th Time step.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_BackUpTsArray( fsls_OilGasData *oilgasdata, int SucStep )
{
   double *Ts_array;  // Successful time Step array
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   Ts_array = fsls_OilGasDataTsArray(oilgasdata);
   Ts_array[SucStep] = fsls_OilGasParamDt(oilgasparam);
}

/*!
 * \fn void fsls_BackUpTsArray
 * \brief Back up the SucStep-th Time step.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_BackUpTnArray( fsls_OilGasData *oilgasdata, int SucStep, double Tnow ) 
{
   double *Tn_array  = fsls_OilGasDataTnArray(oilgasdata);
   Tn_array[SucStep] = Tnow;
}
        
/*!
 * \fn void fsls_BackUpPwArray
 * \brief Back up the SucStep-th Well Pressure.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_BackUpPwArray( fsls_OilGasData *oilgasdata, int SucStep )
{
   double *Pw_array;  // Successful time Step array
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   Pw_array = fsls_OilGasDataPwArray(oilgasdata);
   Pw_array[SucStep] = fsls_OilGasParamPw(oilgasparam);
}

/*!
 * \fn void fsls_BackUpSwArray
 * \brief Back up the SucStep-th Well Saturation.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_BackUpSwArray( fsls_OilGasData *oilgasdata, int SucStep )
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   int    Nx = fsls_OilGasParamNx(oilgasparam);
   int    Wx = fsls_OilGasParamWx(oilgasparam);
   int    Wy = fsls_OilGasParamWy(oilgasparam);
   double *Snow     = fsls_OilGasDataSnow(oilgasdata);
   double *Sw_array = fsls_OilGasDataSwArray(oilgasdata);  // Successful time Step array   
   double Sw;
   
   Sw = Snow[(Wy-1)*Nx + (Wx-1)];   
   
   Sw_array[SucStep] = Sw;
   fsls_OilGasParamSw(oilgasparam) = Sw;
}

/*!
 * \fn void fsls_BackUpORArray
 * \brief Back up the SucStep-th Well Oil Production Rate.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_BackUpORArray( fsls_OilGasData *oilgasdata, int SucStep )
{
   double *OR_array;  // Well Oil Production Rate for each successful time step
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   OR_array = fsls_OilGasDataORArray(oilgasdata);
   OR_array[SucStep] = fsls_OilGasParamOR(oilgasparam);
}

/*!
 * \fn void fsls_LastPressureSaturation
 * \brief Print the last pressure and saturation to given files.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_LastPressureSaturation( fsls_OilGasData *oilgasdata, char *filename_P, char *filename_S )
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam);
   
   double *Pnow = fsls_OilGasDataPnow(oilgasdata);
   double *Snow = fsls_OilGasDataSnow(oilgasdata);
   
   int Nxy = Nx*Ny;
   int i;
   
   FILE *fp_P = NULL;
   FILE *fp_S = NULL;
   
   fp_P = fopen(filename_P, "w");
   fp_S = fopen(filename_S, "w");
   fprintf(fp_P, "%d\n", Nxy);
   fprintf(fp_S, "%d\n", Nxy);
   for (i = 0; i < Nxy; i ++)
   {
      fprintf(fp_P, "%.15le\n", Pnow[i]);
      fprintf(fp_S, "%.15le\n", Snow[i]);
   }
   fclose(fp_P);
   fclose(fp_S);
}

/*!
 * \fn void fsls_PostProcessWellInfoStep
 * \brief Post-processing.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_PostProcessWellInfoStep( fsls_OilGasData *oilgasdata, int SucStep )
{
   double  *Ts_array = fsls_OilGasDataTsArray(oilgasdata);
   double  *Tn_array = fsls_OilGasDataTnArray(oilgasdata);
   double  *Pw_array = fsls_OilGasDataPwArray(oilgasdata);
   double  *Sw_array = fsls_OilGasDataSwArray(oilgasdata);
   double  *OR_array = fsls_OilGasDataORArray(oilgasdata);
   int     *It_array = fsls_OilGasDataItArray(oilgasdata);
   
   char *filename_Pw = NULL; 
   char *filename_Sw = NULL;
   char *filename_OR = NULL;
   char *filename_Ts = NULL;
   char *filename_Tn = NULL;
   char *filename_It = NULL;
   
   FILE *fp_Pw = NULL;
   FILE *fp_Sw = NULL;
   FILE *fp_OR = NULL;
   FILE *fp_Ts = NULL;
   FILE *fp_Tn = NULL;
   FILE *fp_It = NULL;
   
   int i;
   
   filename_Pw = "./output/Pw_step.dat";
   filename_Sw = "./output/Sw_step.dat";
   filename_OR = "./output/OR_step.dat";
   filename_Ts = "./output/Ts_step.dat";
   filename_Tn = "./output/Tn_step.dat";
   filename_It = "./output/It_step.dat";

   fp_Pw = fopen(filename_Pw, "w");
   fp_Sw = fopen(filename_Sw, "w");
   fp_OR = fopen(filename_OR, "w");
   fp_Ts = fopen(filename_Ts, "w");
   fp_Tn = fopen(filename_Tn, "w");
   fp_It = fopen(filename_It, "w");
 
   for (i = 0; i <= SucStep; i ++)
   {
      fprintf(fp_Pw, "%3d  %8.2lf\n", i, Pw_array[i]);
      fprintf(fp_Sw, "%3d  %5.2lf\n", i, Sw_array[i]);
      fprintf(fp_OR, "%3d  %8.2lf\n", i, OR_array[i]);
      fprintf(fp_Ts, "%3d  %5.2lf\n", i, Ts_array[i]);
      fprintf(fp_Tn, "%3d  %5.2lf\n", i, Tn_array[i]);
      fprintf(fp_It, "%3d  %2d\n",    i, It_array[i]);
   }
  
   fclose(fp_Pw);
   fclose(fp_Sw);
   fclose(fp_OR);
   fclose(fp_Ts);
   fclose(fp_Tn);
   fclose(fp_It);
}

/*!
 * \fn void fsls_PostProcessWellInfoTnow
 * \brief Post-processing.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_PostProcessWellInfoTnow( fsls_OilGasData *oilgasdata, int SucStep )
{
   double  *Ts_array = fsls_OilGasDataTsArray(oilgasdata);
   double  *Tn_array = fsls_OilGasDataTnArray(oilgasdata);
   double  *Pw_array = fsls_OilGasDataPwArray(oilgasdata);
   double  *Sw_array = fsls_OilGasDataSwArray(oilgasdata);
   double  *OR_array = fsls_OilGasDataORArray(oilgasdata);
   int     *It_array = fsls_OilGasDataItArray(oilgasdata);
   
   char *filename_Pw = NULL; 
   char *filename_Sw = NULL;
   char *filename_OR = NULL;
   char *filename_Ts = NULL;
   char *filename_It = NULL;
   
   FILE *fp_Pw = NULL;
   FILE *fp_Sw = NULL;
   FILE *fp_OR = NULL;
   FILE *fp_Ts = NULL;
   FILE *fp_It = NULL;
   
   int i;
   
   filename_Pw = "./output/Pw_tnow.dat";
   filename_Sw = "./output/Sw_tnow.dat";
   filename_OR = "./output/OR_tnow.dat";
   filename_Ts = "./output/Ts_tnow.dat";
   filename_It = "./output/It_tnow.dat";
   
   fp_Pw = fopen(filename_Pw, "w");
   fp_Sw = fopen(filename_Sw, "w");
   fp_OR = fopen(filename_OR, "w");
   fp_Ts = fopen(filename_Ts, "w");
   fp_It = fopen(filename_It, "w");
   
   for (i = 0; i <= SucStep; i ++)
   {
      fprintf(fp_Pw, "%8.2lf  %8.2lf\n", Tn_array[i], Pw_array[i]);
      fprintf(fp_Sw, "%8.2lf  %5.2lf\n", Tn_array[i], Sw_array[i]);
      fprintf(fp_OR, "%8.2lf  %8.2lf\n", Tn_array[i], OR_array[i]);
      fprintf(fp_Ts, "%8.2lf  %5.2lf\n", Tn_array[i], Ts_array[i]);
      fprintf(fp_It, "%8.2lf  %2d\n",    Tn_array[i], It_array[i]);
   }
   
   fclose(fp_Pw);
   fclose(fp_Sw);
   fclose(fp_OR);
   fclose(fp_Ts);
   fclose(fp_It);
}

/*!
 * \fn void fsls_PostProcessPSDistribute
 * \brief Post-processing.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_PostProcessPSDistribute( fsls_OilGasData *oilgasdata, int step )
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam);

   double **P_array  = fsls_OilGasDataPArray(oilgasdata);
   double **S_array  = fsls_OilGasDataSArray(oilgasdata);
   
   char *filename_P = NULL; 
   char *filename_S = NULL;
   
   char  NewPFile[120];
   char  NewSFile[120];    
   
   FILE *fp_P = NULL;
   FILE *fp_S = NULL;
   
   int i, j, k, jNx;
   
   filename_P = "./output/P";
   filename_S = "./output/S";
 
   sprintf(NewPFile, "%s_%d", filename_P, step);
   sprintf(NewSFile, "%s_%d", filename_S, step);   
   
   fp_P = fopen(NewPFile, "w");
   fp_S = fopen(NewSFile, "w");
   
   for (j = 0; j < Ny; j ++)
   {
      jNx = j*Nx;
      for (i = 0; i < Nx; i ++)
      {
         k = jNx + i;

         fprintf(fp_P, "%.15le\n", P_array[step][k]);
         fprintf(fp_S, "%.15le\n", S_array[step][k]);     
      }
   }
   
   fclose(fp_P);
   fclose(fp_S);
}

/*!
 * \fn void fsls_PostProcessPSLast
 * \brief Post-processing.
 * \author peghoty
 * \date 2011/09/09
 */ 
void
fsls_PostProcessPSLast( fsls_OilGasData *oilgasdata, int step )
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam);

   double *Pnow = fsls_OilGasDataPnow(oilgasdata);
   double *Snow = fsls_OilGasDataSnow(oilgasdata);
   
   char *filename_P = NULL; 
   char *filename_S = NULL;
   
   char  NewPFile[120];
   char  NewSFile[120];    
   
   FILE *fp_P = NULL;
   FILE *fp_S = NULL;
   
   int i, j, k, jNx;
   
   filename_P = "./output/P_last";
   filename_S = "./output/S_last";
 
   sprintf(NewPFile, "%s_%d", filename_P, step);
   sprintf(NewSFile, "%s_%d", filename_S, step);   
   
   fp_P = fopen(NewPFile, "w");
   fp_S = fopen(NewSFile, "w");
   
   for (j = 0; j < Ny; j ++)
   {
      jNx = j*Nx;
      for (i = 0; i < Nx; i ++)
      {
         k = jNx + i;

         fprintf(fp_P, "%.15le\n", Pnow[k]);
         fprintf(fp_S, "%.15le\n", Snow[k]);     
      }
   }
   
   fclose(fp_P);
   fclose(fp_S);
}
