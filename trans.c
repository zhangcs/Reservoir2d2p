/*!
 *   trans.c -- Update the Transmissibility
 * 
 *       This subroutine is to compare the pressure between neighboring  
 *   grids, and then calculate the Timissibility using upwinding shemes. 
 *   The derivatives of Transmissibility are also calculated.
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
 * \fn int fsls_UpdateTransmissibility
 * \brief Update the Transmissibility.
 * \param oilgasdata pointer to the fsls_OilGasData object 
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_UpdateTransmissibility( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int    Nx    = fsls_OilGasParamNx(oilgasparam);
   int    Ny    = fsls_OilGasParamNy(oilgasparam);
   double Dx    = fsls_OilGasParamDx(oilgasparam);
   double Dy    = fsls_OilGasParamDy(oilgasparam);
   double Dz    = fsls_OilGasParamDz(oilgasparam);
   double Kx    = fsls_OilGasParamKx(oilgasparam);
   double Ky    = fsls_OilGasParamKy(oilgasparam);
   double alpha = fsls_OilGasParamAlpha(oilgasparam);

   /* arrays from 'oilgasdata' */ 
   double *Muo  = fsls_OilGasDataMuo(oilgasdata);
   double *dMuo = fsls_OilGasDataDMuo(oilgasdata);
   double *Mug  = fsls_OilGasDataMug(oilgasdata);
   double *dMug = fsls_OilGasDataDMug(oilgasdata);
   double *Bo   = fsls_OilGasDataBo(oilgasdata);   // use '->Bo' not '->Bonow'
   double *dBo  = fsls_OilGasDataDBo(oilgasdata);
   double *Bg   = fsls_OilGasDataBg(oilgasdata);   // use '->Bg' not '->Bgnow'
   double *dBg  = fsls_OilGasDataDBg(oilgasdata);
   double *Kro  = fsls_OilGasDataKro(oilgasdata);
   double *dKro = fsls_OilGasDataDKro(oilgasdata);
   double *Krg  = fsls_OilGasDataKrg(oilgasdata);
   double *dKrg = fsls_OilGasDataDKrg(oilgasdata);  
   double *Pnow = fsls_OilGasDataPnow(oilgasdata);
   
   // Transmissibility and their derivatives   
   
   /* x-direction, Oil */
   double *Yox1    = fsls_OilGasDataYox1(oilgasdata);         
   double *Yox2    = fsls_OilGasDataYox2(oilgasdata);        
   double *dpYox1  = fsls_OilGasDataDpYox1(oilgasdata);       
   double *dpYox2  = fsls_OilGasDataDpYox2(oilgasdata);       
   double *dsYox1  = fsls_OilGasDataDsYox1(oilgasdata);       
   double *dsYox2  = fsls_OilGasDataDsYox2(oilgasdata);       
   double *ddpYox1 = fsls_OilGasDataDdpYox1(oilgasdata);      
   double *ddpYox2 = fsls_OilGasDataDdpYox2(oilgasdata);      
   double *ddsYox1 = fsls_OilGasDataDdsYox1(oilgasdata);      
   double *ddsYox2 = fsls_OilGasDataDdsYox2(oilgasdata); 
        
   /* x-direction, Gas */
   double *Ygx1    = fsls_OilGasDataYgx1(oilgasdata);         
   double *Ygx2    = fsls_OilGasDataYgx2(oilgasdata);         
   double *dpYgx1  = fsls_OilGasDataDpYgx1(oilgasdata);       
   double *dpYgx2  = fsls_OilGasDataDpYgx2(oilgasdata);       
   double *dsYgx1  = fsls_OilGasDataDsYgx1(oilgasdata);       
   double *dsYgx2  = fsls_OilGasDataDsYgx2(oilgasdata);       
   double *ddpYgx1 = fsls_OilGasDataDdpYgx1(oilgasdata);      
   double *ddpYgx2 = fsls_OilGasDataDdpYgx2(oilgasdata);      
   double *ddsYgx1 = fsls_OilGasDataDdsYgx1(oilgasdata);      
   double *ddsYgx2 = fsls_OilGasDataDdsYgx2(oilgasdata); 
        
   /* y-direction, Oil */
   double *Yoy1    = fsls_OilGasDataYoy1(oilgasdata);         
   double *Yoy2    = fsls_OilGasDataYoy2(oilgasdata);         
   double *dpYoy1  = fsls_OilGasDataDpYoy1(oilgasdata);       
   double *dpYoy2  = fsls_OilGasDataDpYoy2(oilgasdata);      
   double *dsYoy1  = fsls_OilGasDataDsYoy1(oilgasdata);       
   double *dsYoy2  = fsls_OilGasDataDsYoy2(oilgasdata);       
   double *ddpYoy1 = fsls_OilGasDataDdpYoy1(oilgasdata);      
   double *ddpYoy2 = fsls_OilGasDataDdpYoy2(oilgasdata);      
   double *ddsYoy1 = fsls_OilGasDataDdsYoy1(oilgasdata);      
   double *ddsYoy2 = fsls_OilGasDataDdsYoy2(oilgasdata);  
       
   /* y-direction, Gas */
   double *Ygy1    = fsls_OilGasDataYgy1(oilgasdata);         
   double *Ygy2    = fsls_OilGasDataYgy2(oilgasdata);         
   double *dpYgy1  = fsls_OilGasDataDpYgy1(oilgasdata);       
   double *dpYgy2  = fsls_OilGasDataDpYgy2(oilgasdata);       
   double *dsYgy1  = fsls_OilGasDataDsYgy1(oilgasdata);       
   double *dsYgy2  = fsls_OilGasDataDsYgy2(oilgasdata);       
   double *ddpYgy1 = fsls_OilGasDataDdpYgy1(oilgasdata);      
   double *ddpYgy2 = fsls_OilGasDataDdpYgy2(oilgasdata);      
   double *ddsYgy1 = fsls_OilGasDataDdsYgy1(oilgasdata);      
   double *ddsYgy2 = fsls_OilGasDataDdsYgy2(oilgasdata);      
   
   /* auxiliary variables */  
   int i, j, k, k1, k2;
   int Nx1 = Nx - 1;
   int Ny1 = Ny - 1;
   int jNx;
   
   // Geometric part of 'Y'
   double coefx = alpha*Kx*Dy*Dz / Dx;
   double coefy = alpha*Ky*Dx*Dz / Dy;
   
   double MuBo, MuBg;
   double coreo, coreg;

   
   //IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII//
   //      Compare the pressure between every neighboring grids     //
   //IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII//
   
   
   //=====================================//
   //       X   -   D I R E C T I O N     //
   //=====================================//
   
   for (j = 0; j < Ny; j ++)
   {
   
      jNx = j*Nx;
      
      for (i = 0; i < Nx1; i ++)
      {
      
         k  = jNx + i;
         k1 = k + 1;
         
         //----------------------------------------------//
         //        Calculate X trans for oil             //
         //----------------------------------------------//
         
         if (Pnow[k] > Pnow[k1])
         {       
            MuBo       = 1.0 / (Muo[k]*Bo[k]);
            coreo      = Kro[k] * MuBo;
            
            Yox2[k]    = coefx * coreo;
            dpYox2[k]  = coefx / Muo[k] * (dBo[k]*Kro[k] - coreo*dMuo[k]);
            dsYox2[k]  = coefx * MuBo * dKro[k];       
            ddpYox2[k] = 0.0;          
            ddsYox2[k] = 0.0;            
         }
         else
         {
            MuBo       = 1.0 / (Muo[k1]*Bo[k1]);
            coreo      = Kro[k1] * MuBo;
            
            Yox2[k]    = coefx * coreo;
            dpYox2[k]  = 0.0;       
            dsYox2[k]  = 0.0;       
            ddpYox2[k] = coefx / Muo[k1] * (dBo[k1]*Kro[k1] - coreo*dMuo[k1]);          
            ddsYox2[k] = coefx * MuBo * dKro[k1];
         }
         

         //----------------------------------------------//
         //        Calculate X trans for gas             //
         //----------------------------------------------//
            
         if (Pnow[k] > Pnow[k1])
         {      
            MuBg       = 1.0 / (Mug[k]*Bg[k]);
            coreg      = Krg[k] * MuBg;
            
            Ygx2[k]    = coefx * coreg;
            dpYgx2[k]  = coefx / Mug[k] * (dBg[k]*Krg[k] - coreg*dMug[k]);       
            dsYgx2[k]  = coefx * MuBg * dKrg[k];      
            ddpYgx2[k] = 0.0;          
            ddsYgx2[k] = 0.0;                            
         }
         else
         {
            MuBg       = 1.0 / (Mug[k1]*Bg[k1]);
            coreg      = Krg[k1] * MuBg;
            
            Ygx2[k]    = coefx * coreg;
            dpYgx2[k]  = 0.0;       
            dsYgx2[k]  = 0.0;       
            ddpYgx2[k] = coefx / Mug[k1] * (dBg[k1]*Krg[k1] - coreg*dMug[k1]);          
            ddsYgx2[k] = coefx * MuBg * dKrg[k1];   
         }

      } // end for i
      
   } // end for j


   //=====================================//
   //       Y   -   D I R E C T I O N     //
   //=====================================//
   
   for (j = 0; j < Ny1; j ++)
   {
   
      jNx = j*Nx;
      
      for (i = 0; i < Nx; i ++)
      {
      
         k  = jNx + i;
         k2 = k + Nx;
         
         //----------------------------------------------//
         //        Calculate Y trans for oil             //
         //----------------------------------------------//
         
         if (Pnow[k] > Pnow[k2])
         {  
            MuBo       = 1.0 / (Muo[k]*Bo[k]);
            coreo      = Kro[k] * MuBo;
            
            Yoy2[k]    = coefy * coreo;
            dpYoy2[k]  = coefy / Muo[k] * (dBo[k]*Kro[k] - coreo*dMuo[k]);
            dsYoy2[k]  = coefy * MuBo * dKro[k];       
            ddpYoy2[k] = 0.0;          
            ddsYoy2[k] = 0.0;                  
         }
         else
         {   
            MuBo       = 1.0 / (Muo[k2]*Bo[k2]);
            coreo      = Kro[k2] * MuBo;
            
            Yoy2[k]    = coefy * coreo;
            dpYoy2[k]  = 0.0;       
            dsYoy2[k]  = 0.0;       
            ddpYoy2[k] = coefy / Muo[k2] * (dBo[k2]*Kro[k2] - coreo*dMuo[k2]);          
            ddsYoy2[k] = coefy * MuBo * dKro[k2];  
         }
         

         //----------------------------------------------//
         //        Calculate Y trans for gas             //
         //----------------------------------------------//
            
         if (Pnow[k] > Pnow[k2])
         {  
            MuBg       = 1.0 / (Mug[k]*Bg[k]);
            coreg      = Krg[k] * MuBg;
            
            Ygy2[k]    = coefy * coreg;
            dpYgy2[k]  = coefy / Mug[k] * (dBg[k]*Krg[k] - coreg*dMug[k]);      
            dsYgy2[k]  = coefy * MuBg * dKrg[k];      
            ddpYgy2[k] = 0.0;          
            ddsYgy2[k] = 0.0;                 
         }
         else
         {  
            MuBg       = 1.0 / (Mug[k2]*Bg[k2]);
            coreg      = Krg[k2] * MuBg;
            Ygy2[k]    = coefy * coreg;
            dpYgy2[k]  = 0.0;       
            dsYgy2[k]  = 0.0;       
            ddpYgy2[k] = coefy / Mug[k2] * (dBg[k2]*Krg[k2] - coreg*dMug[k2]);          
            ddsYgy2[k] = coefy * MuBg * dKrg[k2];     
         }

      } // end for i
      
   } // end for j


   //==============================================================//
   //              F R O M    *x2, *y2    T O   *x1, *y1           //
   //==============================================================//

   fsls_FROMxy2TOxy1(Nx, Ny, Yox2, Yoy2, Yox1, Yoy1);
   fsls_FROMxy2TOxy1(Nx, Ny, dpYox2, dpYoy2, dpYox1, dpYoy1);
   fsls_FROMxy2TOxy1(Nx, Ny, dsYox2, dsYoy2, dsYox1, dsYoy1);
   fsls_FROMxy2TOxy1(Nx, Ny, ddpYox2, ddpYoy2, ddpYox1, ddpYoy1);
   fsls_FROMxy2TOxy1(Nx, Ny, ddsYox2, ddsYoy2, ddsYox1, ddsYoy1);

   fsls_FROMxy2TOxy1(Nx, Ny, Ygx2, Ygy2, Ygx1, Ygy1);
   fsls_FROMxy2TOxy1(Nx, Ny, dpYgx2, dpYgy2, dpYgx1, dpYgy1);
   fsls_FROMxy2TOxy1(Nx, Ny, dsYgx2, dsYgy2, dsYgx1, dsYgy1);
   fsls_FROMxy2TOxy1(Nx, Ny, ddpYgx2, ddpYgy2, ddpYgx1, ddpYgy1);
   fsls_FROMxy2TOxy1(Nx, Ny, ddsYgx2, ddsYgy2, ddsYgx1, ddsYgy1);
 
   return 0;
}

