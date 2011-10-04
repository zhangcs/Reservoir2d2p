/*!
 *   residual.c -- Form the Residual, i.e, the right hand 
 *                 side of Jacobian equations.
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
 * \fn int fsls_FormResidual
 * \brief Form the Residual, i.e, the right hand side of Jacobian equations.
 * \param oilgasdata pointer to the fsls_OilGasData object 
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_FormResidual( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int    Nx    = fsls_OilGasParamNx(oilgasparam);
   int    Ny    = fsls_OilGasParamNy(oilgasparam);
   double alpha = fsls_OilGasParamAlpha(oilgasparam);
   double dt    = fsls_OilGasParamDt(oilgasparam);
   double V     = fsls_OilGasParamV(oilgasparam);
   double Phi   = fsls_OilGasParamPhi(oilgasparam);
   int    Wx    = fsls_OilGasParamWx(oilgasparam);
   int    Wy    = fsls_OilGasParamWy(oilgasparam);
   double WI    = fsls_OilGasParamWI(oilgasparam);
   double OR    = fsls_OilGasParamOR(oilgasparam);
   double BHP   = fsls_OilGasParamBHP(oilgasparam);
   double Pw    = fsls_OilGasParamPw(oilgasparam);

   /* arrays from 'oilgasdata' */ 
   double *Muo     = fsls_OilGasDataMuo(oilgasdata);
   double *Mug     = fsls_OilGasDataMug(oilgasdata);
   double *Bo      = fsls_OilGasDataBo(oilgasdata);
   double *Bg      = fsls_OilGasDataBg(oilgasdata);
   double *Kro     = fsls_OilGasDataKro(oilgasdata);
   double *Krg     = fsls_OilGasDataKrg(oilgasdata);  
   double *Pnow    = fsls_OilGasDataPnow(oilgasdata);
   double *Snow    = fsls_OilGasDataSnow(oilgasdata);
   double *R       = fsls_OilGasDataR(oilgasdata);
   double *Sbegin  = fsls_OilGasDataSbegin(oilgasdata);   
   double *Bonow   = fsls_OilGasDataBonow(oilgasdata);
   double *Bgnow   = fsls_OilGasDataBgnow(oilgasdata);
   double *Rsnow   = fsls_OilGasDataRsnow(oilgasdata);   
   double *Rs      = fsls_OilGasDataRs(oilgasdata);
   double *Rsx1    = fsls_OilGasDataRsx1(oilgasdata);
   double *Rsx2    = fsls_OilGasDataRsx2(oilgasdata);
   double *Rsy1    = fsls_OilGasDataRsy1(oilgasdata);
   double *Rsy2    = fsls_OilGasDataRsy2(oilgasdata);
   double *Yox1    = fsls_OilGasDataYox1(oilgasdata);         
   double *Yox2    = fsls_OilGasDataYox2(oilgasdata);           
   double *Ygx1    = fsls_OilGasDataYgx1(oilgasdata);         
   double *Ygx2    = fsls_OilGasDataYgx2(oilgasdata);           
   double *Yoy1    = fsls_OilGasDataYoy1(oilgasdata);         
   double *Yoy2    = fsls_OilGasDataYoy2(oilgasdata);             
   double *Ygy1    = fsls_OilGasDataYgy1(oilgasdata);         
   double *Ygy2    = fsls_OilGasDataYgy2(oilgasdata);              

   /* auxiliary variables */  
   int i, j, k;
   int jNx;
   int ROW = 2*Nx*Ny;
   int Nx1 = Nx - 1;
   int Ny1 = Ny - 1;
   int Wx1 = Wx - 1;
   int Wy1 = Wy - 1;
   
   double P1, P2, P3, P4, Pn;
   double coef = V*Phi / dt;
   double alphaWI = alpha*WI;
   double Tow, Tgw;
    
   for (j = 0; j < Ny; j ++)
   {
      jNx = j*Nx;
      
      for (i = 0; i < Nx; i ++)
      {
      
         k = jNx + i;
         
         //-----------------------------------//
         // Determine if it's on the boundary //
         //-----------------------------------//
         
         if (i == 0)
         {
            P1 = Pnow[k];
         }
         else
         {
            P1 = Pnow[k-1];
         }
         
         if (i == Nx1)
         {
            P2 = Pnow[k];
         }
         else
         {
            P2 = Pnow[k+1];
         }
         
         if (j == 0)
         {
            P3 = Pnow[k];
         }
         else
         {
            P3 = Pnow[k-Nx];
         }
         
         if (j == Ny1)
         {
            P4 = Pnow[k];
         }
         else
         {
            P4 = Pnow[k+Nx];
         }
         
         Pn = Pnow[k];
         
         //----------------------------------------------------//
         //                R1:   Oil Equation                  //
         //----------------------------------------------------//

         /* here, Bonow[k] and Bo[k] are always positive */
         R[2*k]  =   Yox1[k]*(P1 - Pn) + Yox2[k]*(P2 - Pn) 
                   + Yoy1[k]*(P3 - Pn) + Yoy2[k]*(P4 - Pn)
                   - coef*(   (1.0 - Snow[k]) / Bonow[k] 
                            - (1.0 - Sbegin[k]) / Bo[k]  );                                   
                      
         //----------------------------------------------------//
         //                  R2:  Gas Equation                 //
         //----------------------------------------------------//
           
         R[2*k+1] =   (Ygx1[k] + Rsx1[k]*Yox1[k])*(P1 - Pn)
                    + (Ygx2[k] + Rsx2[k]*Yox2[k])*(P2 - Pn)
                    + (Ygy1[k] + Rsy1[k]*Yoy1[k])*(P3 - Pn)
                    + (Ygy2[k] + Rsy2[k]*Yoy2[k])*(P4 - Pn)
                    - coef*(    Rsnow[k]*(1.0 - Snow[k]) / Bonow[k]
                              + Snow[k] / Bgnow[k]
                              - Rs[k]*(1.0 - Sbegin[k]) / Bo[k]
                              - Sbegin[k] / Bg[k]  );

  
         //----------------------------------------------------------//
         //                 W E L L   T E R M                        //
         //----------------------------------------------------------//
         
         /* this should be very careful, not (Wx,Wy) but (Wx-1,Wy-1) 
            since the index starts from 0 in C language.  peghoty */
            
         if (i == Wx1 && j == Wy1) 
         {        
             Tow = alphaWI*Kro[k] / (Muo[k]*Bo[k]); 
             Tgw = alphaWI*Krg[k] / (Mug[k]*Bg[k]);
             fsls_OilGasParamTow(oilgasparam) = Tow;
             fsls_OilGasParamTgw(oilgasparam) = Tgw;              

             if (Pw > BHP)
             {
                 R[2*k]   -= OR;
                 R[2*k+1] -= OR*( Rsnow[k] + Tgw / Tow );
             }
             else
             {
                 Pw = BHP;                             
                 fsls_OilGasParamPw(oilgasparam) = Pw;
                 OR = Tow*(Pn - BHP);
                 fsls_OilGasParamOR(oilgasparam) = OR;

                 R[2*k]   -= Tow*(Pn - BHP);
                 R[2*k+1] -= (Rs[k]*Tow + Tgw)*(Pn - BHP);
             }
         }
         
      } // end for i
      
   } // end for j


   //-----------------------------------//
   //  Make 'R' be the right hand side  //
   //  of the Jacobian linear system    //
   //-----------------------------------//
         
   for (i = 0; i < ROW; i ++)
   {
      R[i] = - R[i];
   } 
     
   return 0;
}
