/* $Id: quadrule.h,v 1.15 2019/08/28 02:20:24 zlb Exp $ */

#if !defined(QUADRULE_H)
#define QUADRULE_H

/*----------------------- nonlinear solvers ------------------------*/
#define NEWTON		0	/* the built-in Newton solver (NE) */
#define MINPACK		1	/* MINPACK, LS or NE */
#define TENSOLVE	2	/* TENSOLVE, least square */
#define PETSC		3	/* PETSc, nonlinear equations */
#define CCMNF		4	/* Chen Chuanmiao's Newton flow algorithm */
#define GSLSIMAN	5	/* GSL's simulated annealing algorithm */
/*------------------------------------------------------------------*/

/*------------------------- MINPACK solvers ------------------------*/
#define LMDIF	0	/* Levenberg-Marquardt algorithm, FD Jacobian */
#define LMDER	1	/* Levenberg-Marquardt algorithm, analytic Jacobian */
#define HYBRD	2	/* Powell hybrid method, FD Jacobian */
#define HYBRJ	3	/* Powell hybrid method, analytic Jacobian */
static const char *minpack_solvers[] = {"lmdif1", "lmder1", "hybrd1", "hybrj1"};
#ifndef MINPACK_MAXFEV
# define MINPACK_MAXFEV	500	/* func evals <= MINPACK_MAXFEV*(nabsc+1) */
#endif
/*------------------------------------------------------------------*/

/*--------------------- floating point types -----------------------*/
#define FT_FLOAT	0
#define FT_DOUBLE	1
#define FT_LONG_DOUBLE	2
#define FT___FLOAT128	3
/*------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
extern long double strtold (const char *nptr, char **endptr);
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>

/* TRAP_FPE!=0 ==> trap floating point exceptions (Linux only?) */
#ifndef TRAP_FPE
# define TRAP_FPE	0	/* TRAP_FPE!=0 not working! */
#endif
#if TRAP_FPE
# define __USE_GNU
# include <fenv.h>		/* needs fenableexcept() */
# include <signal.h>
/*# include <fpu_control.h>*/
#endif  /* TRAP_FPE */

#if NONLINEAR_SOLVER == PETSC
# include <petscsnes.h>
#endif	/* NONLINEAR_SOLVER == PETSC */

#if FT_PHG == FT___FLOAT128
#include <quadmath.h>
#endif	/* FT_PHG == ... */

#if USE_MPI
# include <mpi.h>
#endif	/* USE_MPI */

#ifndef Unused
# define Unused(v)	(void)(v)
#endif

#if FT_PHG == FT_FLOAT
# error wll not work!
typedef float FLOAT;
# define Exp(a)		expf((FLOAT)(a))
# define Log(a)		logf((FLOAT)(a))
# define Sqrt(a)	sqrtf((FLOAT)(a))
# define Pow(a,b)	powf((FLOAT)(a), (FLOAT)(b))
# define Atanh(a)	atanhf((FLOAT)(a))
# define Tanh(a)	tanhf((FLOAT)(a))
# define Fabs(a)	fabsf((FLOAT)(a))
# define Asin(a)	asinf((FLOAT)(a))
# define PI		M_PI
# define STR2FLOAT	strtof
#elif FT_PHG == FT_DOUBLE
typedef double FLOAT;
# define Exp(a)		exp((FLOAT)(a))
# define Log(a)		log((FLOAT)(a))
# define Sqrt(a)	sqrt((FLOAT)(a))
# define Pow(a,b)	pow((FLOAT)(a), (FLOAT)(b))
# define Atanh(a)	atanh((FLOAT)(a))
# define Tanh(a)	tanh((FLOAT)(a))
# define Fabs(a)	fabs((FLOAT)(a))
# define Asin(a)	asin((FLOAT)(a))
# define PI		M_PI
# define STR2FLOAT	strtod
#elif FT_PHG == FT_LONG_DOUBLE
typedef long double FLOAT;
# define Exp(a)		expl((FLOAT)(a))
# define Log(a)		logl((FLOAT)(a))
# define Sqrt(a)	sqrtl((FLOAT)(a))
# define Pow(a,b)	powl((FLOAT)(a), (FLOAT)(b))
# define Atanh(a)	atanhl((FLOAT)(a))
# define Tanh(a)	tanhl((FLOAT)(a))
# define Fabs(a)	fabsl((FLOAT)(a))
# define Asin(a)	asinl((FLOAT)(a))
# define PI		3.1415926535897932384626433832795029L
# define STR2FLOAT	strtold
#elif FT_PHG == FT___FLOAT128
typedef __float128 FLOAT;
# define Exp(a)		expq((FLOAT)(a))
# define Log(a)		logq((FLOAT)(a))
# define Sqrt(a)	sqrtq((FLOAT)(a))
# define Pow(a,b)	powq((FLOAT)(a), (FLOAT)(b))
# define Atanh(a)	atanhq((FLOAT)(a))
# define Tanh(a)	tanhq((FLOAT)(a))
# define Fabs(a)	fabsq((FLOAT)(a))
# define Asin(a)	asinq((FLOAT)(a))
# define PI		3.1415926535897932384626433832795029Q
# define STR2FLOAT	strtoflt128
#endif	/* FT_PHG == ... */

#if FT_PHG == FT_FLOAT || FT_PHG == FT_DOUBLE
# define _F(n)  n               /* no suffix */
#elif FT_PHG == FT_LONG_DOUBLE
# define _F(n)  n ## L          /* 'L' suffix */
#elif FT_PHG == FT___FLOAT128
# if HAVE_Q_SUFFIX
#  define _F(n) n ## Q          /* 'Q' suffix */
# else
#  define _F(n) n ## L          /* fallback to 'L' suffix */
# endif
#else
# error unexpected!
#endif

#endif	/* !defined(QUADRULE_H) */
