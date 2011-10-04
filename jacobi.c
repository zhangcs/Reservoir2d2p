/*!
 *   jacob.c -- Form Jacobian Matrix for Newton Iteration.
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
 * \fn int fsls_FormJacobiMatrixFull
 * \brief Form (Full) Jacobian Matrix for Newton Iteration.
 * \param oilgasdata pointer to the fsls_OilGasData object 
 * \author peghoty
 * \date 2011/09/09
 */ 
int 
fsls_FormJacobiMatrixFull( fsls_OilGasData *oilgasdata ) 
{
   fsls_OilGasParam *oilgasparam = fsls_OilGasDataOilgasparam(oilgasdata);

   /* parameters from 'oilgasparam' */ 
   int    Nx     = fsls_OilGasParamNx(oilgasparam);
   int    Ny     = fsls_OilGasParamNy(oilgasparam);
   double alpha  = fsls_OilGasParamAlpha(oilgasparam);
   double dt     = fsls_OilGasParamDt(oilgasparam);
   double V      = fsls_OilGasParamV(oilgasparam);
   double Phi    = fsls_OilGasParamPhi(oilgasparam);
   int    Wx     = fsls_OilGasParamWx(oilgasparam);
   int    Wy     = fsls_OilGasParamWy(oilgasparam);
   double WI     = fsls_OilGasParamWI(oilgasparam);
   double OR     = fsls_OilGasParamOR(oilgasparam);
   double BHP    = fsls_OilGasParamBHP(oilgasparam);
   double Pw     = fsls_OilGasParamPw(oilgasparam);

   /* arrays from 'oilgasdata' */ 
   /* Rs */
   double *Rs    = fsls_OilGasDataRs(oilgasdata);
   double *Rsnow = fsls_OilGasDataRsnow(oilgasdata);
   double *dRs   = fsls_OilGasDataDRs(oilgasdata);
   double *Rsx1  = fsls_OilGasDataRsx1(oilgasdata);
   double *Rsx2  = fsls_OilGasDataRsx2(oilgasdata);
   double *Rsy1  = fsls_OilGasDataRsy1(oilgasdata);
   double *Rsy2  = fsls_OilGasDataRsy2(oilgasdata);

   double *Muo   = fsls_OilGasDataMuo(oilgasdata);
   double *dMuo  = fsls_OilGasDataDMuo(oilgasdata);
   double *Mug   = fsls_OilGasDataMug(oilgasdata);
   double *dMug  = fsls_OilGasDataDMug(oilgasdata);
   double *Bo    = fsls_OilGasDataBo(oilgasdata);
   double *Bonow = fsls_OilGasDataBonow(oilgasdata);
   double *dBo   = fsls_OilGasDataDBo(oilgasdata);
   double *Bg    = fsls_OilGasDataBg(oilgasdata);
   double *Bgnow = fsls_OilGasDataBgnow(oilgasdata);
   double *dBg   = fsls_OilGasDataDBg(oilgasdata);
   double *Kro   = fsls_OilGasDataKro(oilgasdata);
   double *dKro  = fsls_OilGasDataDKro(oilgasdata);
   double *Krg   = fsls_OilGasDataKrg(oilgasdata);
   double *dKrg  = fsls_OilGasDataDKrg(oilgasdata);  
   double *Pnow  = fsls_OilGasDataPnow(oilgasdata);
   double *Snow  = fsls_OilGasDataSnow(oilgasdata);
   
   /* Transmissibility and their derivatives */
   // x-direction, Oil
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
   // x-direction, Gas
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
   // y-direction, Oil
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
   // y-direction, Gas
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
   
   double Tow = fsls_OilGasParamTow(oilgasparam);
   double Tgw = fsls_OilGasParamTgw(oilgasparam);   
   
   double *Jacobi  = fsls_OilGasDataJacobi(oilgasdata);    
   
   /* auxiliary variables */  
   int    i, j, k;
   int    k2, k2p1;
   int    k2ROW, k2p1ROW;
   int    jNx;
   int    ROW = 2*Nx*Ny;   
   int    Nx1 = Nx - 1;
   int    Ny1 = Ny - 1; 
   int    Wx1 = Wx - 1;
   int    Wy1 = Wy - 1;
   int    Nxx = 2*Nx;     
   double a11, a12, a21, a22;
   double b11, b12, b21, b22;
   double P1, P2, P3, P4, Pn;
   double coef    = Phi*V / dt;
   double alphaWI = alpha*WI;
   double ratio;
   double cco, ccg;
   double p1_pn, p2_pn, p3_pn, p4_pn, pn_pw;

   //-------------------------------------------------------//
   //  Zero out the Jacobian matrix at the very beginning.  //
   //  This is quite necessary since the entries will be    //
   //  modified during th Gaussian-Elimination.  peghoty    //
   //-------------------------------------------------------//
   
   fsls_ArrayInitialize(Jacobi, ROW*ROW); 
   

   //-------------------------------------------------------//
   //      Compute entries for the Jacobian Matrix          //
   //-------------------------------------------------------//

   for (j = 0; j < Ny; j ++)
   {
      jNx = j*Nx;
      
      for (i = 0; i < Nx; i ++)
      {
      
         k       = jNx + i;
         k2      = 2*k;
         k2p1    = k2 + 1;
         k2ROW   = k2*ROW;
         k2p1ROW = k2p1*ROW;

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
         
         p1_pn = P1 - Pn; 
         p2_pn = P2 - Pn;
         p3_pn = P3 - Pn;
         p4_pn = P4 - Pn;
         
         //=====================================================//
         //             M A I N     D I A G O N A L             //
         //=====================================================//

         a11 =   ddpYox1[k] * p1_pn - Yox1[k]
               + ddpYox2[k] * p2_pn - Yox2[k]
               + ddpYoy1[k] * p3_pn - Yoy1[k]
               + ddpYoy2[k] * p4_pn - Yoy2[k]
               - coef * (1.0 - Snow[k]) * dBo[k];
         
         a12 =   ddsYox1[k] * p1_pn
               + ddsYox2[k] * p2_pn
               + ddsYoy1[k] * p3_pn
               + ddsYoy2[k] * p4_pn
               + coef / Bonow[k];
               
         a21 =   (ddpYgx1[k] + 0.5*dRs[k]*Yox1[k] + Rsx1[k]*ddpYox1[k]) * p1_pn - Ygx1[k] - Rsx1[k]*Yox1[k]
               + (ddpYgx2[k] + 0.5*dRs[k]*Yox2[k] + Rsx2[k]*ddpYox2[k]) * p2_pn - Ygx2[k] - Rsx2[k]*Yox2[k]
               + (ddpYgy1[k] + 0.5*dRs[k]*Yoy1[k] + Rsy1[k]*ddpYoy1[k]) * p3_pn - Ygy1[k] - Rsy1[k]*Yoy1[k] 
               + (ddpYgy2[k] + 0.5*dRs[k]*Yoy2[k] + Rsy2[k]*ddpYoy2[k]) * p4_pn - Ygy2[k] - Rsy2[k]*Yoy2[k]
               - coef * ( Snow[k]*dBg[k] + (1.0 - Snow[k])*(dRs[k]/Bonow[k] + Rsnow[k]*dBo[k]) );
      
         a22 =   (ddsYgx1[k] + Rsx1[k]*ddsYox1[k]) * p1_pn
               + (ddsYgx2[k] + Rsx2[k]*ddsYox2[k]) * p2_pn
               + (ddsYgy1[k] + Rsy1[k]*ddsYoy1[k]) * p3_pn
               + (ddsYgy2[k] + Rsy2[k]*ddsYoy2[k]) * p4_pn
               - coef * ( 1.0/Bgnow[k] - Rsnow[k]/Bonow[k] );

         Jacobi[k2ROW + k2]     = a11;
         Jacobi[k2ROW + k2p1]   = a12;
         Jacobi[k2p1ROW + k2]   = a21;
         Jacobi[k2p1ROW + k2p1] = a22;
         

         //=====================================================//
         //                W E L L   T E R M S                  //
         //=====================================================//
         
         /* this should be very careful, not (Wx,Wy) but (Wx-1,Wy-1) 
            since the index starts from 0 in C language.  peghoty */
            
         if (i == Wx1 && j == Wy1) 
         {           
 
            if (Pw > BHP)
            {       
               ratio = Tgw / Tow;
               cco   = dMuo[k]/Muo[k] - Bo[k]*dBo[k];
               ccg   = dMug[k]/Mug[k] - Bg[k]*dBg[k];
               
               b21   = OR * ( dRs[k] + ratio * (cco - ccg) );
               
               b22   = OR * (Muo[k]*Bo[k] / (Kro[k]*Mug[k]*Bg[k])) * (dKrg[k] - dKro[k]*Krg[k]/Kro[k]);

               Jacobi[k2p1ROW + k2]   -= b21;
               Jacobi[k2p1ROW + k2p1] -= b22;           
            }
            else
            {
               Pw = BHP;
               fsls_OilGasParamPw(oilgasparam) = Pw;
               
               pn_pw = Pn - Pw;

               b11  =   Tow * ( 1.0 - pn_pw * (dMuo[k] / Muo[k] - Bo[k]*dBo[k]) );

               b12  =   alphaWI * pn_pw * dKro[k] / (Muo[k]*Bo[k]);
     
               b21  =   Tgw * ( 1.0 - pn_pw * (dMug[k] / Mug[k] - Bg[k]*dBg[k]) )
                      + dRs[k] * Tow * pn_pw
                      + Rs[k] * b11;      
                    
               b22  =   alphaWI * pn_pw * dKrg[k] / (Mug[k]*Bg[k]) + Rs[k] * b12;       

               Jacobi[k2ROW + k2]     -= b11;
               Jacobi[k2ROW + k2p1]   -= b12;
               Jacobi[k2p1ROW + k2]   -= b21;
               Jacobi[k2p1ROW + k2p1] -= b22;
            }
            
         } // end if (i == Wx1 && j == Wy1) 
         

         //=====================================================//
         //            U P P E R    D I A G O N A L             //
         //=====================================================//

         if (j != Ny1)
         {      
            a11 = Yoy2[k] + dpYoy2[k] * p4_pn;
            
            a12 = dsYoy2[k]*p4_pn;
         
            a21 =   (dpYgy2[k] + 0.5*dRs[k]*Yoy2[k] + Rsy2[k]*dpYoy2[k]) * p4_pn
                  + Ygy2[k] + Rsy2[k]*Yoy2[k]; 
            
            a22 = (dsYgy2[k] + Rsy2[k]*dsYoy2[k]) * p4_pn;

            Jacobi[k2ROW + k2+Nxx]     = a11;
            Jacobi[k2ROW + k2p1+Nxx]   = a12;
            Jacobi[k2p1ROW + k2+Nxx]   = a21;
            Jacobi[k2p1ROW + k2p1+Nxx] = a22;      
         }


         //=====================================================//
         //            L O W E R    D I A G O N A L             //
         //=====================================================//
         
         if (j != 0)
         {      
            a11 = Yoy1[k] + dpYoy1[k] * p3_pn;
            
            a12 = dsYoy1[k] * p3_pn;
            
            a21 =   (dpYgy1[k] + 0.5*dRs[k]*Yoy1[k] 
                  + Rsy1[k]*dpYoy1[k]) * p3_pn
                  + Ygy1[k] + Rsy1[k]*Yoy1[k];
                  
            a22 =   (dsYgy1[k] + Rsy1[k]*dsYoy1[k]) * p3_pn;

            Jacobi[k2ROW + k2-Nxx]     = a11;
            Jacobi[k2ROW + k2p1-Nxx]   = a12;
            Jacobi[k2p1ROW + k2-Nxx]   = a21;
            Jacobi[k2p1ROW + k2p1-Nxx] = a22;
         }


         //=====================================================//
         //         U P P E R   O F F   D I A G O N A L         //
         //=====================================================//

         if (i != Nx1)
         {
            a11 = Yox2[k] + dpYox2[k] * p2_pn;
            a12 = dsYox2[k] * p2_pn;
            a21 =   Ygx2[k] + Rsx2[k]*Yox2[k] 
                  + ( dpYgx2[k] + 0.5*dRs[k]*Yox2[k] + Rsx2[k]*dpYox2[k] ) * p2_pn;
            a22 = ( dsYgx2[k] + Rsx2[k]*dsYox2[k] )*p2_pn;    

            Jacobi[k2ROW + k2+2]     = a11;
            Jacobi[k2ROW + k2p1+2]   = a12;
            Jacobi[k2p1ROW + k2+2]   = a21;
            Jacobi[k2p1ROW + k2p1+2] = a22;
         }


         //=====================================================//
         //        L O W E R   O F F   D I A G O N A L          //
         //=====================================================//

         if (i != 0)
         {
            a11 = Yox1[k] + dpYox1[k] * p1_pn;
            a12 = dsYox1[k] * p1_pn;
            a21 =   Ygx1[k] + Rsx1[k]*Yox1[k] 
                  + ( dpYgx1[k] + 0.5*dRs[k]*Yox1[k] + Rsx1[k]*dpYox1[k] ) * p1_pn;
            a22 = ( dsYgx1[k] + Rsx1[k]*dsYox1[k] ) * p1_pn;       

            Jacobi[k2ROW + k2-2]     = a11;
            Jacobi[k2ROW + k2p1-2]   = a12;
            Jacobi[k2p1ROW + k2-2]   = a21;
            Jacobi[k2p1ROW + k2p1-2] = a22;           
         }
   
      } // end for i   
      
   } // end for j

   return 0;
}
