/*!
 * tool.c -- some basic functions.
 *
 * Created by peghoty 2011/09/09
 * Xiangtan University
 * peghoty@163.com
 *  
 */

#include "basic.h"
#include "util.h"

/*!
 * \fn fsls_Free
 * \brief free memory
 * \date 2011/09/09
 */
void
fsls_Free( char *ptr )
{
   if (ptr) free(ptr);
}

/*!
 * \fn fsls_CAlloc
 * \brief allocate memory
 * \date 2011/09/09
 */
char *
fsls_CAlloc( size_t count, size_t elt_size )
{
   char *ptr  = NULL;
   int   size = count*elt_size;

   if (size > 0)
   {
      ptr = calloc(count, elt_size);
      if (ptr == NULL)
      {
         fsls_OutOfMemory(size);
      }
   }
   else
   {
      ptr = NULL;
   }

   return ptr;
}

/*!
 * \fn fsls_OutOfMemory
 * \brief print information when out of memory
 * \date 2011/09/09
 */
int
fsls_OutOfMemory( size_t size )
{
   printf("\n \033[31mOut of memory\033[00m trying to allocate \033[31m%d\033[00m bytes!\n", (int) size);
   fflush(stdout);
   return 0;
}

/*!
 * \fn int fsls_ArrayInitialize
 * \brief Zero out a double array.
 * \param *x pointer to the double array.
 * \param n size of the 'x' array.
 * \author peghoty
 * \date 2011/09/09
 */
void 
fsls_ArrayInitialize( double *x, int n )
{
   memset(x, 0X0, n*sizeof(double));
}

/*!
 * \fn double fsls_Arrayl2Norm
 * \brief Compute the l^2 norm of a given vector
 * \param *x pointer to the vector
 * \return l2norm, the l^2 norm of x
 * \author peghoty
 * \date 2011/09/09 
 */
double   
fsls_Arrayl2Norm( double *x, int n )
{
   int i;
   double l2norm = 0.0;
   
   for (i = 0; i < n; i ++)
   {
      l2norm += x[i]*x[i];
   }
   l2norm = sqrt(l2norm); 
     
   return l2norm;
}

/*!
 * \fn fsls_gselim
 * \brief Gaussian Elimination
 * \param *A the pointer to the matrix to be solved
 * \param *x the pointer to the rhs vector at first and the solution at last
 * \param n the size of the matrix
 * \return 0 if success, 1 if fail
 * \date 2011/09/09
 */
int 
fsls_gselim( double *A, double *x, int n )
{
   int    err_flag = 0;
   int    j,k,m;
   double factor = 0.0;
   
   if (n == 1)   /* A is 1x1 */  
   {
      if (A[0] != 0.0)
      {
         x[0] = x[0] / A[0];
         return(err_flag);
      }
      else
      {
         err_flag = 1;
         return(err_flag);
      }
   }
   else   /* A is nxn */   
   {
   
      /* Forward elimination */
      for (k = 0; k < n-1; k ++)
      {
          if (A[k*n+k] != 0.0)
          {          
             for (j = k+1; j < n; j ++)
             {
                 if (A[j*n+k] != 0.0)
                 {
                    factor = A[j*n+k] / A[k*n+k];
                    for (m = k+1; m < n; m ++)
                    {
                        A[j*n+m] -= factor * A[k*n+m];
                    }
                    /* Elimination step for rhs */ 
                    x[j] -= factor * x[k];              
                 }
             }
          }
       }
       
       /* Back Substitution */
       for (k = n-1; k > 0; -- k)
       {
           x[k] /= A[k*n+k];
           for (j = 0; j < k; j ++)
           {
               if (A[j*n+k] != 0.0)
               {
                  x[j] -= x[k] * A[j*n+k];
               }
           }
       }
       
       x[0] /= A[0];
       
       return(err_flag);
    }
}
