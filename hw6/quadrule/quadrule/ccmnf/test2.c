/*
 * =====================================================================================
 *
 *       Filename:  test2.c
 *
 *    Description:  test ccmfindroots
 *
 *        Version:  1.0
 *        Created:  12/08/2011 09:28:45 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ma Shichao (Ph.D), masc@lsec.cc.ac.cn
 *        Company:  PHG_group lsec
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ccmNF.h"


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  func
 *  Description:  a test nonlinear equation
 * =====================================================================================
 */
static void
func ( INT N, INT M, const DOUBLE *X, DOUBLE *Func, DOUBLE *J)
{
    int i;

    assert(N == M);

    if (Func != NULL) {
	Func[0] = X[0] * X[0] * X[0] * X[0] + X[1] - 8.0;
	Func[1] = X[1] * X[1] * X[1] * X[1] + 2 * X[2] - 7.0;
	Func[2] = X[2] * X[2] * X[2] * X[2] + 3 * X[3] - 6.0;
	Func[3] = X[3] * X[3] * X[3] * X[3] + 4 * X[0] - 5.0;
    }

    if (J == NULL)
	return;

    bzero(J, N * M * sizeof(*J));

    for(i = 0; i < 4; i++)
    {
        J[i*N+i] = 4. * X[i] * X[i] * X[i];
    }
    
    J[1*N+0] = 1.0;
    J[2*N+1] = 2.0;
    J[3*N+2] = 3.0;
    J[0*N+3] = 4.0;
}		/* -----  end of function func ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PrintStr
 *  Description:  Print Found Roots In Structure
 * =====================================================================================
 */
    void
PrintStr ( Result *str )
{
    int i = 0;
    int j = 0;

    if(str->n == 0)
    {
        printf ( "No Root Found In This Permutation!\n" );
    }
    else
    {
        printf ( "Found Roots = %d   In %d Permutations\n\n", str->n, str->perm );
        for(i = 1; i <= str->n; i++)
        {
            printf ( "X[%d] = [", i );
            for(j = (i-1)*str->N; j < i*str->N; j++)
            {
                printf ( "%10.6f  ", str->X[j] );
            }
            printf ( "]\n" );

            printf ( "Func[%d] = [", i );
            for(j = (i-1)*str->N; j < i*str->N; j++)
            {
                printf ( "%e  ", str->Fval[j] );
            }
            printf ( "]\n\n" );
        }
    }
}		/* -----  end of function PrintStr  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  test ccmfindroots
 * =====================================================================================
 */
    int
main ( int argc, char *argv[] )
{
    int N = 4;
    int K = 100;
    int MaxIt = 50;
    DOUBLE tol = 1e-14;

    Result *res;

    /* Initialize Root Structure */
    res = (Result *) malloc((unsigned int) sizeof(Result));
    res->N = N;
    res->n = 0;
    res->x_size = 8;
    res->block = 1;
    res->perm = K;
    res->X = (DOUBLE *) malloc((unsigned int) (res->N * res->x_size * 
                res->block * sizeof(DOUBLE)));
    res->Fval = (DOUBLE *) malloc((unsigned int) (res->N * res->x_size * 
                res->block * sizeof(DOUBLE)));
    res->Jval = (DOUBLE *) malloc((unsigned int) (res->N * res->N * 
                res->x_size * res->block * sizeof(DOUBLE)));
    
    /* calculation find roots! */
    ccmFindRoots(func, N, N, K, tol, MaxIt, NULL, res, FALSE);
    
    /* Output found roots! */
    PrintStr(res);

    free(res->X);
    free(res->Fval);
    free(res->Jval);
    free(res);

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
