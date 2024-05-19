/* A Newton flow algorithm for solving a system nonlinear equations.
 * This is a reimplementation in C of Chen Chuanmiao's algorithm.
 *
 * $Id: ccmNF.h,v 1.4 2011/12/09 01:40:46 zlb Exp $ */

/*
 * =====================================================================================
 *
 *       Filename:  ccmNF.h
 *
 *    Description:  CCMNF
 *
 *        Version:  1.0
 *        Created:  12/06/2011 01:39:46 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ma Shichao (Ph.D), masc@lsec.cc.ac.cn
 *        Company:  PHG_group lsec
 *
 * =====================================================================================
 */
#ifndef CCMNF_H_
#define CCMNF_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

typedef double DOUBLE;
typedef int INT;

typedef enum { ROOT, SINGULARITY, MAXIT} RetTYPE;

typedef enum { TRUE, FALSE} BOOLTYPE;

typedef void (*PFn)(INT N, INT M, const DOUBLE *X, DOUBLE *F, DOUBLE *JAC);

typedef void (*PIn)(INT N, DOUBLE *X);

typedef struct Result_
{
    int N;               // dimension of variables
    int n;               // number of root
    int x_size;          // root array size
    int block;           // block number
    int perm;            // permutation times
    DOUBLE *X;           //found root
    DOUBLE *Fval;        //root function value
    DOUBLE *Jval;        //root Jacobian matrix
}Result;

/* Function Declearation */
void ccmNFlow ( PFn func, int N, int M, DOUBLE *X, DOUBLE tol, int MaxIt,
        DOUBLE *FuncValue, DOUBLE *JacValue, RetTYPE *Rflag );

void ccmFindRoots ( PFn func, int N, int M, int K, DOUBLE tol, int MaxIt, 
        PIn InitValue, Result *res, BOOLTYPE FindAll );

#endif
