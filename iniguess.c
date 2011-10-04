/*!
 *   iniguess.c -- Initial guess for Newton Iteration
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
 * \fn int fsls_SetNewtonInitialGuess
 * \brief Set initial guess for Newton Iteration.
 * \param oilgasdata pointer to the fsls_OilGasData object.
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_SetNewtonInitialGuess( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);
   
   int Nx  = fsls_OilGasParamNx(oilgasparam);
   int Ny  = fsls_OilGasParamNy(oilgasparam);
   int Nxy = Nx*Ny; 
   
   double *Pbegin = fsls_OilGasDataPbegin(oilgasdata);
   double *Sbegin = fsls_OilGasDataSbegin(oilgasdata);
   double *Pnow   = fsls_OilGasDataPnow(oilgasdata);
   double *Snow   = fsls_OilGasDataSnow(oilgasdata);

   memcpy(Pnow, Pbegin, Nxy*sizeof(double));
   memcpy(Snow, Sbegin, Nxy*sizeof(double));

   return 0;
} 

