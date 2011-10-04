/*!
 *   sdterms.c -- Update the Saturation Dependent Terms
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
 * \fn int fsls_UpdateSaturationDependentTerms
 * \brief Update the Saturation Dependent Terms.
 * \param oilgasdata pointer to the fsls_OilGasData object 
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_UpdateSaturationDependentTerms( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int  Nx = fsls_OilGasParamNx(oilgasparam);
   int  Ny = fsls_OilGasParamNy(oilgasparam);

   /* arrays from 'oilgasdata' */ 
   double *Snow = fsls_OilGasDataSnow(oilgasdata);
   double *Kro  = fsls_OilGasDataKro(oilgasdata);
   double *dKro = fsls_OilGasDataDKro(oilgasdata);
   double *Krg  = fsls_OilGasDataKrg(oilgasdata);
   double *dKrg = fsls_OilGasDataDKrg(oilgasdata);
   
   /* auxiliary variables */
   int i;
   int Nxy = Nx*Ny;
   double S_threshold = 0.85;
   double theta;
   double beta  = 1.5;
   double gamma = 0.5;
   double coef  = - beta / S_threshold;
   
   //=====================================//
   //  Oil phase relative permeabilities  //
   //=====================================//
   
   for (i = 0; i < Nxy; i ++)
   {
      if (Snow[i] > S_threshold)
      {
         Kro[i]  = 0.0;
         dKro[i] = 0.0;
      }
      else
      {
         theta   = (S_threshold - Snow[i]) / S_threshold; // S_threshold = 0.85 is nonzero!
         Kro[i]  = pow(theta, beta);
         dKro[i] = coef*pow(theta, gamma); 
      }  
   }
     
   //=====================================//
   //  Gas phase relative permeabilities  //
   //=====================================// 
   
   for (i = 0; i < Nxy; i ++)
   {
      if (Snow[i] <= 0)
      {
         Snow[i] = 0.0;
         Krg[i]  = 0.0;
         dKrg[i] = 0.0;
      }
      else 
      {
         if (Snow[i] >= S_threshold)
         {
            Krg[i]  = pow(S_threshold, beta);
            dKrg[i] = 0.0;
         }
         else
         {
            Krg[i]  = pow(Snow[i], beta);
            dKrg[i] = beta*pow(Snow[i], gamma);
         }           
      }
   }

   return 0;
}

