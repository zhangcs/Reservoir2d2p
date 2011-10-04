/*!
 *       [== Black Oil Reservior Simulation Project ==]
 *   
 *   This program is for 2d 2phase(oil and gas) simulation.
 *
 *      The main algorithm is based on the matlab code of RS_4A 
 *   from Chensong and Xiaozhe.
 *
 *   Created by peghoty  2011/09/09
 *   Xiangtan University
 *   peghoty@163.com
 *
 */

#include "basic.h"
#include "util.h"
#include "blackoil.h"

int 
main( int argc, char *argv[] )
{
   struct timeval tStart,tEnd;

   int arg_index   = 0;
   int print_usage = 0;   
   
   fsls_OilGasParam  *oilgasparam = NULL;
   fsls_OilGasData   *oilgasdata  = NULL; 
 
   double Tend;     // The Total simulation time
   double Tnow;     // The current simulation time
   
   int    NCFlag;   // Mark whether the Newton iteration is convergenced
   int    SucStep;  // The number of successful   time steps
   int    UncStep;  // The number of unconvergent time steps
   int    MaxStep;  // The maximum simulation step allowed 
   int    outloop;  // variable for outer loop
   int    numiter;  // number of nonlinear iterations 
   
   int    Time_Test   = 1; // whether timing the whole simulation
   int    print_level = 1; // whether output the parameters 
   
   int    nx = 11;
   int    ny = 11;
   
   while (arg_index < argc)
   {
       if ( strcmp(argv[arg_index], "-nx") == 0 )
       {
           arg_index ++;
           nx = atoi(argv[arg_index++]);
       }
       else if ( strcmp(argv[arg_index], "-ny") == 0 )
       {
           arg_index ++;
           ny = atoi(argv[arg_index++]);
       }
       else if ( strcmp(argv[arg_index], "-tt") == 0 )
       {
           arg_index ++;
           Time_Test = atoi(argv[arg_index++]);
       }        
       else if ( strcmp(argv[arg_index], "-pl") == 0 )
       {
           arg_index ++;
           print_level = atoi(argv[arg_index++]);
       }    
       else if ( strcmp(argv[arg_index], "-help") == 0 )
       {
           print_usage = 1;
           break;
       }
       else
       {
           arg_index ++;
       }         
   } 
 
   if (print_usage)
   {
      printf("\n");
      printf("  Usage: %s [<options>]\n", argv[0]);
      printf("\n");
      printf("  -nx    <val> : number of grids in X-direction [default: 11]\n");
      printf("  -ny    <val> : number of grids in Y-direction [default: 11]\n");
      printf("  -tt    <val> : whether timing the whole simulation [default: 1]\n");
      printf("  -pl    <val> : whether print the soling information [default: 1]\n");
      printf("  -help        : using help message\n\n");      
      exit(1);
   }  
   
    
   if (Time_Test) GetTime(tStart); 
    
    
   //=================================================// 
   // Create and Initialize a fsls_OilGasParam object //
   //=================================================// 
   oilgasparam = fsls_OilGasParamCreate();
   fsls_OilGasParamInitialize(oilgasparam, nx, ny, 1);
   if (print_level)
   {
      fsls_OilGasParamOutput(oilgasparam);
   }
   
   
   //=================================================// 
   // Create and Initialize a fsls_OilGasData object  //
   //=================================================// 
   oilgasdata = fsls_OilGasDataCreate(oilgasparam);
   fsls_OilGasDataInitialize(oilgasdata);


   //========================================================// 
   // Set initial Pressure, Saturation, 'Bo', 'Bg', and 'Rs' //
   //========================================================// 
   fsls_InitialPressureSaturation(oilgasdata);
   fsls_ComputeBoBgRs(oilgasdata, 0); 
      

   //IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII// 
   //                M A I N    L O O P                //
   //IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII// 
   MaxStep = fsls_OilGasParamMaxStep(oilgasparam);
   Tend    = fsls_OilGasParamTend(oilgasparam);

   /* zero out the caculators */
   Tnow    = 0.0;
   SucStep = 0;
   UncStep = 0;
   numiter = 0;

   fsls_BackUpItArray(oilgasdata, SucStep, numiter);
   fsls_BackUpTnArray(oilgasdata, SucStep, Tnow);
   fsls_BackUpPwArray(oilgasdata, SucStep);
   fsls_BackUpSwArray(oilgasdata, SucStep);
   fsls_BackUpORArray(oilgasdata, SucStep);  
 
   if (print_level)
   {
      printf(" ==================================\n");
      printf("       NL      Tnow          dt    \n");
      printf(" ==================================\n");
   }

   for (outloop = 0; outloop < MaxStep; outloop ++)
   {

      //==========================================================// 
      //     Step 1:  Set initial guess for Newton Iteration      //
      //==========================================================// 
      fsls_SetNewtonInitialGuess(oilgasdata);


      //=================================================// 
      //       Step 2:  Newton Iteration Process         //
      //=================================================// 
      numiter = fsls_NewtonIteration(oilgasdata);
  
      NCFlag = fsls_OilGasParamNCFlag(oilgasparam);
      if (NCFlag == 1)
      {
         SucStep ++; 
         Tnow += fsls_OilGasParamDt(oilgasparam);
         fsls_BackUpPSData(oilgasdata, SucStep);
         fsls_BackUpItArray(oilgasdata, SucStep, numiter);
         fsls_BackUpTsArray(oilgasdata, SucStep);
         fsls_BackUpTnArray(oilgasdata, SucStep, Tnow);
         fsls_BackUpPwArray(oilgasdata, SucStep);
         fsls_BackUpSwArray(oilgasdata, SucStep);
         fsls_BackUpORArray(oilgasdata, SucStep);
      }
      else
      {
         UncStep ++;
      }

      if (print_level)
      {
         printf(" %3d  %2d    %.4le   %.4le\n", 
                outloop+1, numiter, Tnow, fsls_OilGasParamDt(oilgasparam));
      }

      //===========================================================//
      //  Step 3:  Update the time step for next Newton Iteration  //
      //===========================================================//
      fsls_GetNextDt(oilgasdata, Tnow);


      //=================================================//
      //    Step 4:   Judge if the total time arrives    //
      //=================================================//
      if (Tnow >= Tend)
      {
         printf("\n >>> Simulation is Over!\n\n");
         printf(" >>> SucStep        = %d\n", SucStep);
         printf(" >>> UncStep        = %d\n", UncStep);
         printf(" >>> NumUnconvType1 = %d\n", fsls_OilGasParamNumUnconvType1(oilgasparam));
         printf(" >>> NumUnconvType2 = %d\n", fsls_OilGasParamNumUnconvType2(oilgasparam));
         break;
      }     

   } // end for 'outloop'
   
   
   if (Time_Test) 
   {
      GetTime(tEnd);
      printf("\n >>> \033[31mTotal Simulation Time:\033[00m %.3lf seconds\n\n",mytime(tStart,tEnd)); 
   } 
    
   fsls_PostProcessWellInfoStep(oilgasdata, SucStep);      
   fsls_PostProcessWellInfoTnow(oilgasdata, SucStep);
   fsls_PostProcessPSLast(oilgasdata, SucStep);
   
   if ( fsls_OilGasParamPostpro(oilgasparam) )
   {
      /* Save any 'P' and 'S' vectors of your interest. */
      fsls_PostProcessPSDistribute(oilgasdata, 0);
      fsls_PostProcessPSDistribute(oilgasdata, 9);
      fsls_PostProcessPSDistribute(oilgasdata, 10);
      fsls_PostProcessPSDistribute(oilgasdata, 11);
      fsls_PostProcessPSDistribute(oilgasdata, 12);
      fsls_PostProcessPSDistribute(oilgasdata, SucStep);  
   }

   fsls_OilGasDataFinalize(oilgasdata);

   return (0);
}
