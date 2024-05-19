/* A Newton flow algorithm for solving a system nonlinear equations.
 * This is a reimplementation in C of Chen Chuanmiao's algorithm.
 *
 * $Id: ccmNF.c,v 1.4 2011/12/09 01:40:46 zlb Exp $ */

/*
 * =====================================================================================
 *
 *       Filename:  ccmNF.c
 *
 *    Description:  CCMNF
 *
 *        Version:  1.0
 *        Created:  12/06/2011 01:39:06 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ma Shichao (Ph.D), masc@lsec.cc.ac.cn
 *        Company:  PHG_group lsec
 *
 * =====================================================================================
 */

#include "ccmNF.h"


/*-----------------------------------------------------------------------------
 *  Norm functions for vector
 *-----------------------------------------------------------------------------*/


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Norm2
 *  Description:  
 * =====================================================================================
 */
static DOUBLE
Norm2 ( DOUBLE *X, INT N )
{
    INT i = 0;
    DOUBLE sum = 0.0;

    for(i = 0; i < N; i++)
    {
        sum += X[i] * X[i];
    }

    return sqrt(sum);
}		/* -----  end of function Norm2  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  NormInf
 *  Description:  
 * =====================================================================================
 */
static DOUBLE
NormInf ( DOUBLE *X, INT N )
{
    INT i = 0;
    DOUBLE max = fabs(X[0]);
    
    for(i = 1; i < N; i++)
    {
        if(fabs(X[i]) > max)
        {
            max = fabs(X[i]);
        }
    }

    return max;
}		/* -----  end of function NormInf  ----- */



/*-----------------------------------------------------------------------------
 *  Initial Value Generator
 *-----------------------------------------------------------------------------*/


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  InitGene
 *  Description:  Initial value generator
 * =====================================================================================
 */
static void
InitGene ( INT N, DOUBLE LBC, DOUBLE RBC, DOUBLE *X )
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
}		/* -----  end of function InitGene  ----- */




/*-----------------------------------------------------------------------------
 *  Solve Linear System
 *-----------------------------------------------------------------------------*/


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ls_factor
 *  Description:  LU factorization of the matrix trans(A)[M][Lda] with full pivoting
 *    Arguments:  M N are the number of rows and columns of trans(A), Lda >= N
 *
 *                'pivot_ptr' is the pointer to the integer array of size 2*N,
 *                either provided by the calling function or alocated in this
 *                function (if *pivot_ptr == NULL), which stores the positions of
 *                pivots: the row and column number for the i-th pivot are stored
 *                as:  (*pivot_ptr)[2*i] and (*pivot_ptr)[2*i+1]
 * =====================================================================================
 */
static INT
ls_factor ( INT m, INT n, INT trans, INT lda, DOUBLE *A, INT **pivot_ptr )
{
    int i, j, k, p, q, nfree = 0, *pivot;
    DOUBLE d, t;
#define a(i,j) A[trans ? (j) * lda + (i) : (i) * lda + (j)]

    assert(pivot_ptr != NULL);

    if (*pivot_ptr == NULL)
	*pivot_ptr = pivot = malloc(2 * n * sizeof(*pivot));
    else
	pivot = *pivot_ptr;

    for (j = 0; j < n && j < m; j++) {
	/* find pivot */
	d = fabs(a(j,j));
	p = q = j;
	for (i = j; i < m; i++) {
	    for (k = j; k < n; k++) {
		if (d < (t = fabs(a(i,k)))) {
		    d = t;
		    p = i;
		    q = k;
		}
	    }
	}
	if (d == 0.0)
	{
	    nfree = n - j;
	    for (i = j; i < n; i++) {
		pivot[2 * i + 1] = i;
		if (i < m)
		    a(i,i) = 0.;
	    }
//fprintf(stderr, "ls_factor(): singular system!\n");
	    break;
	}
	pivot[2 * j] = p;
	pivot[2 * j + 1] = q;
	if (p != j) {
	    /* swap rows j and k */
	    for (i = j; i < n; i++) {
		t = a(j,i);
		a(j,i) = a(p,i);
		a(p,i) = t;
	    }
	}
	if (q != j) {
	    /* swap columns q and j */
	    for (k = 0; k < m; k++) {
		t = a(k, j);
		a(k, j) = a(k, q);
		a(k, q) = t;
	    }
	}
	a(j,j) = d = 1.0 / a(j,j);
	for (i = j + 1; i < n; i++)
	    a(j,i) *= d;
	for (i = j + 1; i < m; i++) {
	    d = a(i,j);
	    for (k = j + 1; k < n; k++)
		a(i,k) -= d * a(j,k);
	}
    }

    for (; j < n; j++)
	pivot[2 * j + 1] = j;

#if 0
	for (i = 0; i < m; i++) {
	    for (k = 0; k < n; k++)
		fprintf(stderr, " %lg", (double)a(i,k));
	    fprintf(stderr, "\n");
	}
	/* check consistency of the equations */
	for (i = j; i < m; i++) {
	    if (fabs(a(i, n)) >= EPS)
		return -1;
	}
#endif

#undef a
    return nfree;
}		/* -----  end of function Ls_Factor  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ls_solve
 *  Description:  solves trans(A)*X = B, A has been factorized by ls_factor()
 *    Arguments:  m and n are the number of rows and columns of trans(A), lda >= n
 *                X and B may overlap if m >= n
 * =====================================================================================
 */
static void
ls_solve ( INT m, INT n, INT trans, INT lda, const DOUBLE *A, const int *pivot,
        DOUBLE *B, DOUBLE *X )
{
    int i, j, p;
    DOUBLE t;
#define a(i,j) A[trans ? (j) * lda + (i) : (i) * lda + (j)]
#define b(i) (B[i])
#define x(i) (X[i])

    for (j = 0; j < n && j < m; j++) {
	if (a(j,j) == 0.)
	    break;
	p = pivot[2 * j];
	if (p != j) {
	    /* swap rows j and k */
	    t = b(j);
	    b(j) = b(p);
	    b(p) = t;
	}
	b(j) *= a(j,j);
	for (i = j + 1; i < m; i++)
	    b(i) -= a(i,j) * b(j);
    }

    for (; j < n && j < m; j++)
	b(j) = 0.;

    for (j = n - 1; j >= 0; j--) {
	x(j) = (j < m ? b(j) : 0.);
	if (j >= m || a(j,j) == 0.)
	    continue;
	for (i = j + 1; i < n; i++)
	    x(j) -= a(j,i) * x(i);
    }

    /* reorder the solution */
    for (j = n - 1; j >= 0; j--) {
	assert(pivot[2 * j + 1] >= j);
	if ((i = pivot[2 * j + 1]) == j)
	    continue;
	t = x(j);
	x(j) = x(i);
	x(i) = t;
    }

#undef a
#undef b
#undef x

    return ;
}		/* -----  end of function ls_solver  ----- */




/*-----------------------------------------------------------------------------
 *  ccmNF functions
 *-----------------------------------------------------------------------------*/


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ccmNFlow
 *  Description:  ccm's Newton Flow Algorithm
 *    Arguments:  func          I   Function for Nonlinear Equation and Jacobian
 *                N             I   Dimension of Variables
 *                M             I   Number of Equations
 *                X             I/O Initial point and Final iteration point
 *                tol           I   tolerance of calculaton
 *                MaxIt         I   Maximum Iteartion Times
 *                FuncValue     O   Function Value at Final point
 *                JacValue      O   Jacobian Matrix at Final point
 *                Rflag         O   Return Status Flag
 * =====================================================================================
 */
    void
ccmNFlow ( PFn func, int N, int M, DOUBLE *X, DOUBLE tol, int MaxIt,
        DOUBLE *FuncValue, DOUBLE *JacValue, RetTYPE *Rflag )
{
#define Func(N, X, F)		func(N, N, X, F, NULL)
#define Jacobian(N, X, JAC)	func(N, N, X, NULL, JAC)
    int i;
    int count = 0;
    int *pivot_ptr;
    DOUBLE Delta = 0.1;
    DOUBLE alpha = 0.5;
    DOUBLE beta = 2.0;
    DOUBLE *Nvf, *X1;
    DOUBLE K0, K1, D1, Q1, H0, G1, NFX0, NFX1;

    pivot_ptr = (int *) malloc((unsigned int) (2 * N * sizeof(DOUBLE)));
    X1 = (DOUBLE *) malloc((unsigned int) (N * sizeof(DOUBLE)));
    Nvf = (DOUBLE *) malloc((unsigned int) (N * sizeof(DOUBLE)));

    func(N, N, X, FuncValue, JacValue);
    NFX0 = Norm2(FuncValue, N);
    
    for(i = 0; i < N; i++)
    {
        FuncValue[i] *= -1.0;
    }

    /* solve linear system to get newton direction */
    ls_factor(M, N, 1, N, JacValue, &pivot_ptr);
    ls_solve(M, N, 1, N, JacValue, pivot_ptr, FuncValue, Nvf);

    K0 = NormInf(Nvf, N);
    H0 = (K0 <= 2.0) ? 0.1 : (Delta / K0);

    while(1)
    {
        /* Test Step Length */
        while(1)
        {
            for(i = 0; i < N; i++)
            {
                X1[i] = X[i];
            }
            for(i = 0; i < N; i++)
            {
                X1[i] += H0 * Nvf[i];
            }
            Func(N, X1, FuncValue);
            NFX1 = Norm2(FuncValue, N);
            G1 = NFX1 / NFX0;
            
            if( (H0 - 1.0) <= 1e-8)
            {
                break;
            }
            else
            {
                if(G1 <= 1.0)
                {
                    break;
                }
                else
                {
                    H0 /= 2.0;
                }
            }
        }

        /* Update Data */
        for(i = 0; i < N; i++)
        {
            X[i] = X1[i];
            FuncValue[i] *= -1.0;
        }
        Jacobian(N, X, JacValue);

        ls_factor(M, N, 1, N, JacValue, &pivot_ptr);
        ls_solve(M, N, 1, N, JacValue, pivot_ptr, FuncValue, Nvf);
        K1 = NormInf(Nvf, N);
        NFX0 = NFX1;
        D1 = (G1 + H0 - 1.0) / H0 / H0;
        Q1 = K1 / K0;
        K0 = K1;
        count++;

        /* test root */
        if( K0 < tol)
        {
            func(N, N, X, FuncValue, JacValue);
            *Rflag = ROOT;
            break;
        }

        /* test singularity */
        if((H0 < Delta/10) && (Q1 > 2.0))
        {
            *Rflag = SINGULARITY;
            break;
        }

        /* Exceed Max Iteration Time */
        if(count > MaxIt)
        {
            func(N, N, X, FuncValue, JacValue);
            *Rflag = MAXIT;
            break;
        }

        /* Set Next Step Length */
        if(K0 > 10 * Delta)
        {
            H0 = Delta / K0;
        }
        else
        {
            if(fabs(D1) <= beta)
            {
                if(fabs(D1) > alpha)
                {
                    H0 = (2.0*H0 > 0.5/fabs(D1)) ? 0.5/fabs(D1) : 2.0*H0;
                }
                else
                {
                    H0 = (2.0*H0 > 1.0) ? 1.0 : 2.0*H0; 
                }
            }
        }
    }

    free(pivot_ptr);
    free(X1);
    free(Nvf);
}		/* -----  end of function ccmNFlow  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ccmFindRoots
 *  Description:  ccm's Newton Flow Algorithm with permutations
 *    Arguments:  func       I  Function for Nonlinear Equations and Jacobian
 *                N          I  Dimension of Variables
 *                M          I  Number of Equations
 *                K          I  Permutation Times
 *                InitValue  I  Initialization Function can be NULL
 *                res        O  Found roots information
 * =====================================================================================
 */
    void
ccmFindRoots ( PFn func, int N, int M, int K, DOUBLE tol, int MaxIt, 
        PIn InitValue, Result *res, BOOLTYPE FindAll )
{
    int i = 0;
    int j = 0;
    int k = 0;
    int pertime = 0;
    int newflag = 0;

    DOUBLE *X;
    DOUBLE *Y;
    DOUBLE *FuncValue;
    DOUBLE *JacValue;
    RetTYPE Rflag;

    X = (DOUBLE *) malloc((unsigned int) (N * sizeof(DOUBLE)));
    Y = (DOUBLE *) malloc((unsigned int) (N * sizeof(DOUBLE)));
    FuncValue = (DOUBLE *) malloc((unsigned int) (N * sizeof(DOUBLE)));
    JacValue = (DOUBLE *) malloc((unsigned int) (N * N * sizeof(DOUBLE)));

    /* can be parallel */
    for(pertime = 0; pertime < K; pertime++)
    {
        /* Generate an initial value */
        if(InitValue != NULL)
        {
            InitValue(N, X);
        }
        else
        {
            InitGene(N, -2.0, 2.0, X);
        }

        /* ccmNFlow calculation */
        ccmNFlow(func, N, M, X, tol, MaxIt, FuncValue, JacValue, &Rflag);

        /* test roots, record in the structure */
        if(Rflag == ROOT)
        {
            if(res->n != 0)
            {
                newflag = 1;
                /* compare to the found roots */
                for(k = 1; k <= res->n; k++)
                {
                    /* residual of current root and found ones */
                    for(i = 0; i < N; i++)
                    {
                        Y[i] = res->X[(k-1)*res->N + i] - X[i];
                    }

                    if(NormInf(Y, N) < 1e-6)
                    {
                        newflag = 0;
                        break;
                    }
                }
            }

            if( res->n == 0 || newflag == 1 )
            {
                res->n++;
                if( res->n > res->x_size * res->block )
                {
                    /* add one block */
                    res->block++;
                    res->X = (DOUBLE *) realloc( res->X, 
                            (unsigned int)(res->N * res->x_size * 
                                res->block * sizeof(DOUBLE)));
                    res->Fval = (DOUBLE *) realloc( res->Fval,
                            (unsigned int)(res->N * res->x_size *
                                res->block * sizeof(DOUBLE)));
                    res->Jval = (DOUBLE *) realloc( res->Jval,
                            (unsigned int)(res->N * res->N * 
                                res->x_size * res->block * sizeof(DOUBLE)));
                }

                /* copy  root  funcvalue  jacvalue */
                for(i = 0; i < N; i++)
                {
                    res->X[(res->n-1)*res->N + i] = X[i];
                    res->Fval[(res->n-1) * res->N + i] = FuncValue[i];
                    for(j = 0; j < N; j++)
                    {
                        res->Jval[(res->n-1)*res->N*res->N+i*N+j] = 
                            JacValue[i*N+j];
                    }
                }
            }
        } /* end if(Rflag == ROOT) */
    } /* end for pertime */

    free(X);
    free(Y);
    free(FuncValue);
    free(JacValue);
}		/* -----  end of function ccmFindRoots  ----- */
