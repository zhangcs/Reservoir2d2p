/*!
 *   pdterms.c -- Update the Pressure Dependent Terms
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
 * \fn int fsls_ComputeBoBgRs
 * \brief Compute 'Bo', 'Bg', and 'Rs'.
 * \param oilgasdata pointer to the fsls_OilGasData object 
 * \param flag 0: use it for the first time; 1: otherwise
 * \author peghoty
 * \date 2011/09/09
 */ 
int
fsls_ComputeBoBgRs( fsls_OilGasData *oilgasdata, int flag )
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam);

   /* arrays from 'oilgasdata' */ 
   double *Rs   = fsls_OilGasDataRs(oilgasdata);
   double *Bo   = fsls_OilGasDataBo(oilgasdata);
   double *Bg   = fsls_OilGasDataBg(oilgasdata);
   double *Pnow = NULL;

   /* auxiliary variables */   
   int    i;
   int    Nxy = Nx*Ny; 
   double pv, pv2;  
   double uu  = -3.922264e-8;
   double vv  = 3.0e-4;
   double ww  = 1.0;
   double xx  = 3.5927e+3;
   double yy  = -1.0226;
   double pp  = 4.76e-5;
   double qq  = 1.206511; 

   double P_threshold = 3824.321712; 
   
   if (flag == 0)
   {
      Pnow = fsls_OilGasDataPbegin(oilgasdata);
   }
   else
   {
      Pnow = fsls_OilGasDataPnow(oilgasdata);
   } 
      
   for (i = 0; i < Nxy; i ++)
   {
      pv  = Pnow[i];
      pv2 = pv*pv;
  
      /* Oil FVF */
      Bo[i] = uu*pv2 + vv*pv + ww;

      /* Gas FVF */
      Bg[i] = xx*pow(pv, yy);   

      /* Gas Solubility */
      if (Pnow[i] > P_threshold)
      {
         Rs[i] = 1.0;
      }
      else
      {
         Rs[i] = pp*pow(pv, qq);
      }     
   }
   
   return 0;
}

/*!
 * \fn int fsls_UpdatePressureDependentTerms
 * \brief Update the Pressure Dependent Terms.
 * \param oilgasdata pointer to the fsls_OilGasData object 
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_UpdatePressureDependentTerms( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */
   int Nx = fsls_OilGasParamNx(oilgasparam);
   int Ny = fsls_OilGasParamNy(oilgasparam);

   /* arrays from 'oilgasdata' */ 
   double *Pnow  = fsls_OilGasDataPnow(oilgasdata);   
   double *Muo   = fsls_OilGasDataMuo(oilgasdata);
   double *dMuo  = fsls_OilGasDataDMuo(oilgasdata);
   double *Mug   = fsls_OilGasDataMug(oilgasdata);
   double *dMug  = fsls_OilGasDataDMug(oilgasdata);
   double *Rsnow = fsls_OilGasDataRsnow(oilgasdata);
   double *dRs   = fsls_OilGasDataDRs(oilgasdata);
   double *Bonow = fsls_OilGasDataBonow(oilgasdata);
   double *dBo   = fsls_OilGasDataDBo(oilgasdata);
   double *Bgnow = fsls_OilGasDataBgnow(oilgasdata);
   double *dBg   = fsls_OilGasDataDBg(oilgasdata);
   double *Rsx1  = fsls_OilGasDataRsx1(oilgasdata);
   double *Rsx2  = fsls_OilGasDataRsx2(oilgasdata);
   double *Rsy1  = fsls_OilGasDataRsy1(oilgasdata);
   double *Rsy2  = fsls_OilGasDataRsy2(oilgasdata);

   /* auxiliary variables */   
   int i, j, k;
   int Nxy = Nx*Ny; 
   int Nx1 = Nx - 1;
   int Ny1 = Ny - 1;
   int jNx;
   
   double pv, pv2;  
   double Bo2, Bg2;
   
   double ee  = 6.0e-8;
   double ee2 = 2.0*ee;
   double ff  = -4.589e-4;
   double gg  = 1.1179;
   
   double rr  = 3.0e-10;
   double rr2 = 2.0*rr;
   double ss  = 1.0e-6;
   double tt  = 1.33e-2;
   
   double uu  = -3.922264e-8;
   double uu2 = 2.0*uu;
   double vv  = 3.0e-4;
   double ww  = 1.0;
   
   double xx   = 3592.7;
   double yy   = -1.0226;
   double xy   = xx*yy;
   double yym1 = yy - 1.0;
   
   double pp   = 0.0000476;
   double qq   = 1.206511; 
   double pq   = pp*qq;
   double qqm1 = qq - 1.0;
   double P_threshold = 3824.321712;  
      
   for (i = 0; i < Nxy; i ++)
   {
      pv  = Pnow[i];
      pv2 = pv*pv;
  
      //================================================//
      //      Oil Viscousity and its derivative         //
      //================================================//
      
      Muo[i]  = ee*pv2 + ff*pv + gg;      
      dMuo[i] = ee2*pv + ff;   

      //================================================//
      //      Gas Viscousity and its derivative         //
      //================================================//
      
      Mug[i]  = rr*pv2 + ss*pv + tt;
      dMug[i] = rr2*pv + ss;

      //================================================//
      //  Oil FVF and derivative ( dBo = d(1/Bo) / dp ) //
      //================================================//
      
      Bonow[i] = uu*pv2 + vv*pv + ww;
      Bo2      = Bonow[i]*Bonow[i];
      dBo[i]   = -(uu2*pv + vv) / Bo2; // Bo2 is always posotive

      //================================================//
      //  Gas FVF and derivative ( dBg = d(1/Bg) / dp ) //
      //================================================//

      Bgnow[i] = xx*pow(pv, yy);
      Bg2      = Bgnow[i]*Bgnow[i];
      dBg[i]   = - (xy*pow(pv, yym1)) / Bg2;      

      //================================================//
      //       Gas Solubility and its derivative        //
      //================================================//
      
      if (Pnow[i] > P_threshold)
      {
         Rsnow[i] = 1.0;
         dRs[i]   = 0.0;
      }
      else
      {
         Rsnow[i] = pp*pow(pv, qq);
         dRs[i]   = pq*pow(pv, qqm1);
      }     
   }

   //===============================================//
   //   Rsx2[k] = (Rs)_{i+0.5,j}, k = (j-1)*Nx + i  //
   //   Rsy2[k] = (Rs)_{i,j+0.5}, k = (j-1)*Nx + i  //
   //===============================================// 
   /* remark: here, we just take the average of the  
      values at the the middle of the adjoint cells. */
    
   for (j = 0; j < Ny; j ++)
   {
      jNx = j*Nx;
      for (i = 0; i < Nx1; i ++)
      {
         k = jNx + i;
         Rsx2[k] = 0.5*(Rsnow[k] + Rsnow[k+1]); // average value of left and right
      }
   }
   for (j = 0; j < Ny1; j ++)
   {
      jNx = j*Nx;
      for (i = 0; i < Nx; i ++)
      {
         k = jNx + i;
         Rsy2[k] = 0.5*(Rsnow[k] + Rsnow[k+Nx]); // average value of up and down
      }
   }  
   
   //===============================================//
   //   Rsx1[k] = (Rs)_{i-0.5,j}, k = (j-1)*Nx + i  //
   //   Rsy1[k] = (Rs)_{i,j-0.5}, k = (j-1)*Nx + i  //
   //===============================================// 
   
   fsls_FROMxy2TOxy1(Nx, Ny, Rsx2, Rsy2, Rsx1, Rsy1); 

   return 0;
}
