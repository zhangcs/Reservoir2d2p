/*!
 *   newton.c -- Newton Iteration Process
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
 * \fn int fsls_BackUpPSData
 * \brief Back up 'Pnow', 'Snow' to 'Pbegin', 'Sbegin' if Newton Iteration converges.
 * \param oilgasdata pointer to the fsls_OilGasData object.
 * \param SucStep The number of successfull iterations for present.
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_BackUpPSData( fsls_OilGasData *oilgasdata, int SucStep ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   
   int Nx      = fsls_OilGasParamNx(oilgasparam);
   int Ny      = fsls_OilGasParamNy(oilgasparam);
   int postpro = fsls_OilGasParamPostpro(oilgasparam);
    
   double *Pbegin = fsls_OilGasDataPbegin(oilgasdata);
   double *Sbegin = fsls_OilGasDataSbegin(oilgasdata);
   double *Pnow   = fsls_OilGasDataPnow(oilgasdata);
   double *Snow   = fsls_OilGasDataSnow(oilgasdata);

   double **P_array;   // The Pressure   sequence for all the successful iteration 
   double **S_array;   // The Saturation sequence for all the successful iteration
   
   int Nxy = Nx*Ny;

   memcpy(Pbegin, Pnow, Nxy*sizeof(double));
   memcpy(Sbegin, Snow, Nxy*sizeof(double));
   
   if (postpro == 1)
   {
      P_array = fsls_OilGasDataPArray(oilgasdata);
      S_array = fsls_OilGasDataSArray(oilgasdata); 
      memcpy(P_array[SucStep], Pnow, Nxy*sizeof(double));
      memcpy(S_array[SucStep], Snow, Nxy*sizeof(double));
   }
   
   return 0;
}

/*!
 * \fn double fsls_DelatSMax
 * \brief Find the largest absolute component in 'deltas'.
 * \param oilgasdata pointer to the fsls_OilGasData object
 * \author peghoty
 * \date 2011/04/10
 */
double
fsls_DelatSMax( fsls_OilGasData *oilgasdata )
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam); 

   double *deltas = fsls_OilGasDataDeltas(oilgasdata); 
   
   int    i;
   int    Nxy = Nx*Ny;
   double deltas_max = 0.0;
   
   deltas_max = fabs(deltas[0]);
   for (i = 1; i < Nxy; i ++)
   {
      if (fabs(deltas[i]) > deltas_max)
      {
         deltas_max = fabs(deltas[i]);
      } 
   }

   return (deltas_max);
}

/*!
 * \fn double fsls_DelatPMax
 * \brief Find the largest absolute component in 'deltap/ave_P'.
 * \param oilgasdata pointer to the fsls_OilGasData object
 * \author peghoty
 * \date 2011/04/10
 */
double
fsls_DelatPMax( fsls_OilGasData *oilgasdata )
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam); 

   double *deltap = fsls_OilGasDataDeltap(oilgasdata); 
   double *Pnow   = fsls_OilGasDataPnow(oilgasdata);
   
   int i;
   int Nxy = Nx*Ny;
   double deltap_max = 0.0;
   double meanpnow   = 0.0;
   
   for (i = 0; i < Nxy; i ++)
   {
      meanpnow += Pnow[i];  
   }
   meanpnow /= Nxy;
      
   deltap_max = fabs(deltap[0]);
   for (i = 1; i < Nxy; i ++)
   {
      if (fabs(deltap[i]) > deltap_max)
      {
         deltap_max = fabs(deltap[i]);
      } 
   }

   return (deltap_max/meanpnow);
}

/*!
 * \fn int fsls_UpdateWellPressure
 * \brief Update the pressure of Well.
 * \param oilgasdata pointer to the fsls_OilGasData object
 * \author peghoty
 * \date 2011/04/11
 */
int 
fsls_UpdateWellPressure( fsls_OilGasData *oilgasdata )
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int    Nx     = fsls_OilGasParamNx(oilgasparam);
   int    Wx     = fsls_OilGasParamWx(oilgasparam);
   int    Wy     = fsls_OilGasParamWy(oilgasparam);
   double OR     = fsls_OilGasParamOR(oilgasparam);
   double BHP    = fsls_OilGasParamBHP(oilgasparam);
   double Pw     = fsls_OilGasParamPw(oilgasparam);
   double Tow    = fsls_OilGasParamTow(oilgasparam);
   
   /* arrays from 'oilgasdata' */
   double *Pnow  = fsls_OilGasDataPnow(oilgasdata); 

   /* auxiliary variables */
   int Wx1 = Wx - 1;
   int Wy1 = Wy - 1;

   if (Pw > BHP)
   {
      Pw = Pnow[Wy1*Nx + Wx1] - OR/Tow;
   }
   else
   {
      Pw = BHP;
   }
   fsls_OilGasParamPw(oilgasparam) = Pw;
  
   return 0;
}

/*!
 * \fn int fsls_NewtonIteration
 * \brief Newton Iteration process.
 * \param oilgasdata pointer to the fsls_OilGasData object
 * \author peghoty
 * \date 2011/04/04
 */ 
int 
fsls_NewtonIteration( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   int    Nx      = fsls_OilGasParamNx(oilgasparam);
   int    Ny      = fsls_OilGasParamNy(oilgasparam);   
   double tolN    = fsls_OilGasParamTolN(oilgasparam);
   double tolP    = fsls_OilGasParamTolP(oilgasparam);
   double tolS    = fsls_OilGasParamTolS(oilgasparam);
      
   /* auxiliary variables */
   int    innerloop  = 0;   
   int    Nxy        = Nx*Ny;
   int    Nxy2       = 2*Nxy;
   int    NumStep    = 50;
   double Rnorm      = 0.0;
   double Rnorm_old  = 0.0;
   double deltap_max = 0.0;
   double deltas_max = 0.0;

   /* ------ Zero out 'deltap' and 'deltas' ----- */
   fsls_ArrayInitialize(fsls_OilGasDataDeltap(oilgasdata), Nxy); 
   fsls_ArrayInitialize(fsls_OilGasDataDeltas(oilgasdata), Nxy); 


   //IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII// 
   //               M A I N    N E W T O N    L O O P             //
   //IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII// 
     
   for (innerloop = 1; innerloop < NumStep; innerloop ++)
   {

      //---------------------------------------------//
      //  Step 1: Update Pressure Dependent Terms    //
      //---------------------------------------------//
      
      fsls_UpdatePressureDependentTerms(oilgasdata);
      

      //---------------------------------------------//
      //  Step 2: Update Saturation Dependent Terms  //
      //---------------------------------------------//

      fsls_UpdateSaturationDependentTerms(oilgasdata);
   

      //---------------------------------------------//
      //  Step 3: Update Transmissibility Terms      //
      //---------------------------------------------//

      fsls_UpdateTransmissibility(oilgasdata);
   

      //---------------------------------------------//
      //  Step 4: Form the Residual                  //
      //  (right hand side of Jacobian equations)    //
      //---------------------------------------------//

      fsls_FormResidual(oilgasdata);
      Rnorm_old = Rnorm;
      Rnorm = fsls_Arrayl2Norm(fsls_OilGasDataR(oilgasdata), Nxy2);
      

      //===========================================================//
      //  Step 5: Determine if unconverged of too many iterations  //
      //===========================================================//

      if (innerloop > 6)
      {
         if ( (Rnorm > Rnorm_old) || (innerloop > 20) )
         {
            if (Rnorm > Rnorm_old)
            {
               fsls_OilGasParamNumUnconvType1(oilgasparam) ++; 
            }
            
            if (innerloop > 20)
            {
               fsls_OilGasParamNumUnconvType2(oilgasparam) ++;
            } 
            
            /* Mark unconvergence of Newton Iteration */
            fsls_OilGasParamNCFlag(oilgasparam) = 0; 
            
            break;
         }
      }
        
      //===========================================================//
      //  Step 6: Determine the convergence(Using both Criterion)  //
      //===========================================================//  

      if ( Rnorm < tolN )
      {
         deltas_max = fsls_DelatSMax(oilgasdata);
         if ( deltas_max < tolS )
         {
            deltap_max = fsls_DelatPMax(oilgasdata);
            if ( deltap_max < tolP )
            {
               /* Mark convergence of Newton Iteration */
               fsls_OilGasParamNCFlag(oilgasparam) = 1; 
               
               /* Update the Well Pressure */
               fsls_UpdateWellPressure(oilgasdata); 
               
               /* Compute 'Bo', 'Bg' and 'Rs' for next iteratiton */
               fsls_ComputeBoBgRs(oilgasdata, 1);

               break;
            }
         }   
      }


      //---------------------------------------------//
      //  Step 7: Form the Jacobian Matrix           //
      //---------------------------------------------//

      fsls_FormJacobiMatrixFull(oilgasdata); 


      //---------------------------------------------//
      //  Step 8: Solve the Jacobian Linear System   //
      //---------------------------------------------//

      fsls_JacobiSolveFull(oilgasdata);


   } // end for 'innerloop'

   return (innerloop-1);
}

