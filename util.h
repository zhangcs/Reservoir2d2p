/*!
 * header "util.h"
 *
 * Created by peghoty 2011/09/09
 * Xiangtan University
 * peghoty@163.com
 *  
 */

#ifndef FSLS_UTIL_HEADER
#define FSLS_UTIL_HEADER 

/*----------------------------------------------------------------*
 *                      Macro Definition                          *
 *----------------------------------------------------------------*/ 

/* time testing */  
#define GetTime(a) gettimeofday(&a,NULL)
#define mytime(a,b) ((b.tv_sec-a.tv_sec) + (float)(b.tv_usec-a.tv_usec)/1000000.0) 

/* ratio of the circumference of a circle to its diameter */
#define PI 3.141592653589793232

/* max and min  */ 
#define fsls_max(a,b)  (((a)<(b)) ? (b) : (a))
#define fsls_min(a,b)  (((a)<(b)) ? (a) : (b))

/* allocate and free memory */
#define fsls_TFree(ptr) ( fsls_Free((char *)ptr), ptr = NULL )
#define fsls_CTAlloc(type, count) ( (type *)fsls_CAlloc((size_t)(count), (size_t)sizeof(type)) )

/*----------------------------------------------------------------*
 *                  Functions Declaration                         *
 *----------------------------------------------------------------*/ 
 
/* fsls.c */
void   fsls_Free( char *ptr );
char  *fsls_CAlloc( size_t count, size_t elt_size );
int    fsls_OutOfMemory( size_t size );
void   fsls_ArrayInitialize( double *x, int n );                     
double fsls_Arrayl2Norm( double *x, int n );
int    fsls_gselim(double *A, double *x, int n);

#endif
