/*!
 *   solver.c -- solve the Jacobian linear system.
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
 * \fn int fsls_JacobiSolveFull
 * \brief solve the (full) Jacobian linear system.
 * \param oilgasdata pointer to the fsls_OilGasData object 
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_JacobiSolveFull( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam);

   /* arrays from 'oilgasdata' */ 
   double *Jacobi = fsls_OilGasDataJacobi(oilgasdata);
   double *delta  = fsls_OilGasDataDelta(oilgasdata);
   double *deltap = fsls_OilGasDataDeltap(oilgasdata);
   double *deltas = fsls_OilGasDataDeltas(oilgasdata);   
   double *R      = fsls_OilGasDataR(oilgasdata); 
   double *Pnow   = fsls_OilGasDataPnow(oilgasdata);
   double *Snow   = fsls_OilGasDataSnow(oilgasdata);       
   
   /* auxiliary variables */  
   int i, j;
   int Nxy  = Nx*Ny;
   int Nxy2 = 2*Nxy;   

   //----------------------------------------//
   //  copy the right hand side to 'delta'   //
   //  of the Jacobian linear system         //
   //----------------------------------------//
   
   for (i = 0; i < Nxy2; i ++)
   {
      delta[i] = R[i]; 
   }

   //----------------------------------------//
   //    Solve the Jacobian linear system    //
   //----------------------------------------//
   
   fsls_gselim(Jacobi, delta, Nxy2);
   
   
   //----------------------------------------//
   //  Get 'deltap' and 'deltas' seperately  //
   //----------------------------------------//  
   
   for (i = 0; i < Nxy; i ++)
   {
      j = 2*i;
      deltap[i] = delta[j];
      deltas[i] = delta[j+1];
   }
   
   
   //----------------------------------------//
   //          Update 'P' and 'S'            //
   //----------------------------------------//  
   
   for (i = 0; i < Nxy; i ++)
   {
      Pnow[i] += deltap[i];
      Snow[i] += deltas[i];
   }

   return 0;
}
