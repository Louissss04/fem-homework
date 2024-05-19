/* $Id: test.c,v 1.3 2011/12/09 01:40:46 zlb Exp $ */

/*
 * =====================================================================================
 *
 *       Filename:  test.c
 *
 *    Description:  test ccmnf
 *
 *        Version:  1.0
 *        Created:  12/06/2011 04:38:29 PM
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
#include <time.h>
#include <strings.h>
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
 *         Name:  InitValue
 *  Description:  Initial value generator
 * =====================================================================================
 */
static void
InitValue ( INT N, DOUBLE LBC, DOUBLE RBC, DOUBLE *X )
{
    int i;
    static int flag = 0;
    int seed;

    seed = (unsigned int) time(NULL) + flag * 10;
    flag++;

    srand(seed);
    for(i = 0; i < N; i++)
    {
         X[i] = LBC + (RBC - LBC) * (rand() / (double)(RAND_MAX));
    }
}		/* -----  end of function InitValue  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PrintVec
 *  Description:  Print a Vector on screen
 * =====================================================================================
 */
static void
PrintVec ( DOUBLE *VEC, INT N )
{
    int i;

    for(i = 0; i < N; i++)
    {
        printf("%10.6f  ", VEC[i]);
    }
}		/* -----  end of function PrintVec  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  main function of test ccmNF
 * =====================================================================================
 */
    int
main ( int argc, char *argv[] )
{
    DOUBLE LBC = -2.0;
    DOUBLE RBC = 2.0;
    DOUBLE tol = 1e-14;
    INT MaxIt = 50;
    INT N = 4;

    RetTYPE Rflag;

    DOUBLE *X;
    DOUBLE *FuncValue;
    DOUBLE *JacValue;

    X = calloc(N, sizeof(DOUBLE));
    FuncValue = calloc(N, sizeof(DOUBLE));
    JacValue = calloc(N * N, sizeof(DOUBLE));

    /* Initialization */
    InitValue(N, LBC, RBC, X);

    /* Use ccmNF for calculation */
    ccmNFlow(func, N, N, X, tol, MaxIt, FuncValue, JacValue, &Rflag);

    /* Output the answers */
    printf("\n-----------------------------------------------------------\n");
    printf("        Test for ccmNF.c \n");

    switch ( Rflag ) {
        case ROOT:
            printf("Find a root!\n");
            printf("X = [");
            PrintVec(X, N);
            printf("]\n");
            break;

        case SINGULARITY:
            printf("Meet a singularity!\n");
            break;

        case MAXIT:
            printf("Excced Max Iteration Times!\n");
            break;

        default:
            printf("Error!\n");
            break;
    }				/* -----  end switch  ----- */
    printf("\n-----------------------------------------------------------\n");

     /* Clear the space */
    free(X);
    free(FuncValue);
    free(JacValue);

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
