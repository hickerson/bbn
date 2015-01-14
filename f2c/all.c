/* all.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    real t9, hv, phie, y[26];
} evolp1_;

#define evolp1_1 evolp1_

struct {
    real dt9, dhv, dphie, dydt[26];
} evolp2_;

#define evolp2_1 evolp2_

struct {
    real t90, hv0, phie0, y0[26];
} evolp3_;

#define evolp3_1 evolp3_

struct {
    real cy, ct, t9i, t9f, ytmin;
    integer inc;
} compr_;

#define compr_1 compr_

struct {
    real t, dt, dlt9dt;
} ttime_;

#define ttime_1 ttime_

struct {
    integer ltime, is, ip, it, mbad;
} flags_;

#define flags_1 flags_

struct {
    integer itime;
} tcheck_;

#define tcheck_1 tcheck_

struct {
    integer irun, isize, jsize;
} runopt_;

#define runopt_1 runopt_

struct {
    real f[88], r__[88];
} rates_;

#define rates_1 rates_

struct {
    real g, tau, xnu, c__[3], cosmo, xi[3], b;
} modpr_;

#define modpr_1 modpr_

struct {
    real dt1, eta1;
} varpr_;

#define varpr_1 varpr_

struct {
    real rhone0, rhob0, rhob, rnb;
} endens_;

#define endens_1 endens_

struct {
    real bl1, bl2, bl3, bl4, bl5, bm1, bm2, bm3, bm4, bm5, bn1, bn2, bn3, bn4,
	     bn5;
} besselcb_;

#define besselcb_1 besselcb_

struct {
    real t9mev, tnmev, tnu, cnorm;
    integer nu;
    real rhonu;
} nupar_;

#define nupar_1 nupar_

struct {
    real thm[14], hubcst;
} thermcb_;

#define thermcb_1 thermcb_

struct nucdat_1_ {
    real am[26], zm[26], dm[26];
};

#define nucdat_1 (*(struct nucdat_1_ *) &nucdat_)

struct {
    real xout[1040]	/* was [40][26] */, thmout[240]	/* was [40][6] */, 
	    t9out[40], tout[40], dtout[40], etaout[40], hubout[40];
} outdat_;

#define outdat_1 outdat_

struct {
    real bk0, bk1, bk2, bk3, bk4;
} kays_;

#define kays_1 kays_

struct {
    integer iform[88], ii[88], jj[88], kk[88], ll[88];
    real rev[88], q9[88];
} recpr_;

#define recpr_1 recpr_

union {
    struct {
	doublereal a[676]	/* was [26][26] */;
	real bkevin[26], yx[26];
    } _1;
    struct {
	doublereal a[676]	/* was [26][26] */;
	real bkevin[26], y[26];
    } _2;
    struct {
	doublereal a[676]	/* was [26][26] */;
	real b___[26], yx[26];
    } _3;
} lncoef_;

#define lncoef_1 (lncoef_._1)
#define lncoef_2 (lncoef_._2)
#define lncoef_3 (lncoef_._3)

struct recpr0_1_ {
    real reacpr[704]	/* was [88][8] */;
};

#define recpr0_1 (*(struct recpr0_1_ *) &recpr0_)

struct compr0_1_ {
    real cy0, ct0, t9i0, t9f0, ytmin0;
    integer inc0;
};

#define compr0_1 (*(struct compr0_1_ *) &compr0_)

struct modpr0_1_ {
    real c0[3], cosmo0, xi0[3];
};

#define modpr0_1 (*(struct modpr0_1_ *) &modpr0_)

struct varpr0_1_ {
    real dt0, eta0;
};

#define varpr0_1 (*(struct varpr0_1_ *) &varpr0_)

struct {
    integer itime;
} checkcb_;

#define checkcb_1 checkcb_

struct {
    integer nout;
    logical outfile;
} outopt_;

#define outopt_1 outopt_

/* Initialized data */

struct {
    real e_1[78];
    } nucdat_ = { 1.f, 1.f, 2.f, 3.f, 3.f, 4.f, 6.f, 7.f, 7.f, 8.f, 8.f, 9.f, 
	    10.f, 11.f, 11.f, 12.f, 12.f, 12.f, 13.f, 13.f, 14.f, 14.f, 14.f, 
	    15.f, 15.f, 16.f, 0.f, 1.f, 1.f, 1.f, 2.f, 2.f, 3.f, 3.f, 4.f, 
	    3.f, 5.f, 4.f, 5.f, 5.f, 6.f, 5.f, 6.f, 7.f, 6.f, 7.f, 6.f, 7.f, 
	    8.f, 7.f, 8.f, 8.f, .008665f, .007825f, .014102f, .01605f, 
	    .01603f, .002603f, .015125f, .016004f, .016929f, .022487f, 
	    .024609f, .012186f, .012939f, .009305f, .011432f, .014354f, 0.f, 
	    .018641f, .003354f, .005738f, .003242f, .003074f, .008597f, 
	    1.08e-4f, .00307f, -.005085f };

struct {
    real e_1[704];
    } recpr0_ = { 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 
	    12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 21.f, 22.f, 
	    23.f, 24.f, 25.f, 26.f, 27.f, 28.f, 29.f, 30.f, 31.f, 32.f, 33.f, 
	    34.f, 35.f, 36.f, 37.f, 38.f, 39.f, 40.f, 41.f, 42.f, 43.f, 44.f, 
	    45.f, 46.f, 47.f, 48.f, 49.f, 50.f, 51.f, 52.f, 53.f, 54.f, 55.f, 
	    56.f, 57.f, 58.f, 59.f, 60.f, 61.f, 62.f, 63.f, 64.f, 65.f, 66.f, 
	    67.f, 68.f, 69.f, 70.f, 71.f, 72.f, 73.f, 74.f, 75.f, 76.f, 77.f, 
	    78.f, 79.f, 80.f, 81.f, 82.f, 83.f, 84.f, 85.f, 86.f, 87.f, 88.f, 
	    1.f, 1.f, 4.f, 1.f, 1.f, 4.f, 1.f, 1.f, 1.f, 1.f, 1.f, 2.f, 2.f, 
	    2.f, 2.f, 3.f, 3.f, 3.f, 5.f, 2.f, 2.f, 2.f, 3.f, 5.f, 2.f, 2.f, 
	    2.f, 6.f, 6.f, 3.f, 3.f, 11.f, 9.f, 9.f, 2.f, 2.f, 2.f, 3.f, 3.f, 
	    2.f, 2.f, 2.f, 2.f, 2.f, 3.f, 3.f, 3.f, 3.f, 2.f, 2.f, 2.f, 3.f, 
	    3.f, 3.f, 3.f, 3.f, 3.f, 8.f, 7.f, 9.f, 9.f, 9.f, 10.f, 10.f, 2.f,
	     2.f, 2.f, 3.f, 3.f, 3.f, 3.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 3.f, 
	    2.f, 3.f, 3.f, 3.f, 3.f, 3.f, 3.f, 3.f, 3.f, 3.f, 1.f, 4.f, 10.f, 
	    16.f, 21.f, 11.f, 15.f, 18.f, 20.f, 23.f, 25.f, 2.f, 3.f, 5.f, 
	    7.f, 5.f, 9.f, 7.f, 9.f, 3.f, 4.f, 7.f, 7.f, 8.f, 6.f, 6.f, 6.f, 
	    3.f, 3.f, 4.f, 5.f, 5.f, 8.f, 9.f, 8.f, 13.f, 14.f, 15.f, 13.f, 
	    9.f, 12.f, 13.f, 14.f, 15.f, 16.f, 12.f, 13.f, 16.f, 7.f, 8.f, 
	    9.f, 11.f, 10.f, 12.f, 12.f, 13.f, 14.f, 6.f, 6.f, 10.f, 11.f, 
	    12.f, 14.f, 15.f, 17.f, 19.f, 22.f, 20.f, 22.f, 25.f, 25.f, 17.f, 
	    19.f, 21.f, 20.f, 22.f, 24.f, 24.f, 17.f, 13.f, 14.f, 15.f, 18.f, 
	    20.f, 13.f, 14.f, 16.f, 19.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 
	    0.f, 0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 2.f, 
	    2.f, 2.f, 2.f, 2.f, 3.f, 4.f, 5.f, 0.f, 0.f, 3.f, 3.f, 0.f, 3.f, 
	    3.f, 1.f, 1.f, 1.f, 1.f, 1.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 
	    2.f, 2.f, 6.f, 6.f, 6.f, 6.f, 6.f, 6.f, 3.f, 3.f, 3.f, 1.f, 0.f, 
	    2.f, 1.f, 2.f, 2.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 2.f, 
	    2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 6.f, 6.f, 6.f, 6.f, 6.f, 6.f, 6.f, 
	    6.f, 6.f, 6.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 
	    0.f, 0.f, 0.f, 0.f, 0.f, 2.f, 2.f, 4.f, 0.f, 0.f, 0.f, 0.f, 5.f, 
	    0.f, 0.f, 0.f, 0.f, 1.f, 2.f, 1.f, 2.f, 2.f, 1.f, 2.f, 0.f, 0.f, 
	    0.f, 2.f, 6.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f, 6.f, 6.f, 6.f, 0.f, 
	    0.f, 0.f, 2.f, 1.f, 1.f, 1.f, 2.f, 1.f, 0.f, 0.f, 1.f, 2.f, 3.f, 
	    0.f, 0.f, 0.f, 0.f, 0.f, 2.f, 2.f, 2.f, 6.f, 0.f, 0.f, 0.f, 0.f, 
	    0.f, 0.f, 6.f, 0.f, 2.f, 2.f, 2.f, 2.f, 2.f, 1.f, 1.f, 1.f, 1.f, 
	    2.f, 5.f, 6.f, 17.f, 22.f, 6.f, 14.f, 17.f, 19.f, 22.f, 24.f, 3.f,
	     4.f, 6.f, 8.f, 4.f, 8.f, 6.f, 6.f, 5.f, 6.f, 9.f, 6.f, 6.f, 7.f, 
	    8.f, 9.f, 5.f, 4.f, 6.f, 6.f, 6.f, 6.f, 6.f, 10.f, 14.f, 16.f, 
	    14.f, 8.f, 11.f, 13.f, 15.f, 17.f, 18.f, 17.f, 7.f, 9.f, 12.f, 
	    13.f, 14.f, 15.f, 15.f, 14.f, 17.f, 13.f, 14.f, 17.f, 12.f, 17.f, 
	    6.f, 6.f, 6.f, 6.f, 6.f, 19.f, 21.f, 24.f, 19.f, 21.f, 24.f, 17.f,
	     20.f, 22.f, 24.f, 23.f, 25.f, 26.f, 17.f, 26.f, 19.f, 21.f, 22.f,
	     25.f, 26.f, 20.f, 22.f, 24.f, 26.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
	     0.f, 0.f, 0.f, 0.f, 0.f, .471f, 1.63f, 2.61f, 1.19f, 1.002f, 
	    .998f, 1.07f, 4.7f, 1.63f, 2.61f, 1.19f, 1.07f, 4.69f, 1.53f, 
	    1.11f, 1.11f, 1.73f, 1.73f, 5.54f, 5.55f, 3.39f, 9.95f, 9.97f, 
	    1.31f, 3.04f, 2.34f, 1.002f, .758f, 1.3f, .973f, 3.03f, 7.01f, 
	    2.33f, 3.f, .618f, .754f, .292f, 1.58f, 4.02f, 4.02f, 3.08f, 
	    3.07f, 10.3f, 2.07f, 6.44f, 14.9f, .584f, 2.f, 3.58f, 3.58f, 
	    .807f, 3.5f, 3.49f, .886f, 3.58f, 2.71f, 1.002f, 3.003f, 1.002f, 
	    .709f, .884f, 1.19f, .9f, 3.57f, 2.7f, 3.62f, .706f, 5.13f, 9.36f,
	     11.f, 3.68f, 4.26f, 5.81f, 9.34f, 3.67f, 4.25f, 5.79f, 0.f, 0.f, 
	    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 25.82f, 72.62f, 
	    238.81f, 84.17f, 8.863f, 19.081f, 55.494f, 220.39f, 63.75f, 
	    229.932f, 65.054f, 46.631f, 201.291f, 17.118f, 28.64f, 18.423f, 
	    37.935f, 46.798f, 204.117f, 212.98f, 149.23f, 175.476f, 194.557f, 
	    23.59f, 132.95f, 39.1f, 32.08f, 32.382f, 1.595f, 76.427f, 100.84f,
	     185.173f, 6.975f, 146.08f, 24.674f, 13.301f, 79.89f, 51.753f, 
	    100.538f, 87.539f, 86.f, 76.96f, 66.16f, 50.63f, 107.13f, 159.36f,
	     18.26f, 84.42f, 177.73f, 218.82f, 7.555f, 100.753f, 132.83f, 
	    57.41f, 94.88f, 125.74f, 34.846f, 7.263f, 41.037f, 98.661f, 
	    22.553f, 87.621f, 118.452f, 53.706f, 84.678f, 140.734f, 57.623f, 
	    83.111f, 47.16f, 9.098f, 33.915f, 111.87f, 60.557f, 12.287f, 
	    1.835f, 88.47f, 25.711f };

struct {
    real e_1[5];
    integer e_2;
    } compr0_ = { .3f, .03f, 100.f, .01f, 1e-25f, 30 };

struct {
    real e_1[7];
    } modpr0_ = { 1.f, 880.1f, 3.f, 0.f, 0.f, 0.f, 0.f };

struct {
    real e_1[2];
    } varpr0_ = { 1e-4f, 3.162e-10f };


/* Table of constant values */

static real c_b4 = 1e-5f;
static doublereal c_b5 = 1.5;
static doublereal c_b12 = .33333333333333331;
static doublereal c_b13 = 1.3333333333333333;
static real c_b17 = 0.f;
static integer c__50 = 50;
static doublereal c_b24 = .5;
static integer c__0 = 0;
static integer c__1 = 1;
static real c_b60 = 1.f;
static doublereal c_b69 = .33333333;
static doublereal c_b72 = .66666666666666663;
static doublereal c_b73 = .3333333;
static doublereal c_b74 = .8333333;
static doublereal c_b83 = .2;
static doublereal c_b84 = 1.25;
static doublereal c_b89 = .6;
static doublereal c_b90 = 1.2;
static doublereal c_b91 = .375;
static doublereal c_b94 = .83333333;
static integer c__4 = 4;
static integer c__3 = 3;
static doublereal c_b364 = 10.;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     Changes (to run inder DEC unix f77): */
/*     ----------------------------------- */
/*     COMMON /bessel/ -> COMMON /besselcb/ */
/*     COMMON /therm/  -> COMMON /thermcb/ */
/*     COMMON /time/   -> COMMON /ttime/ */
/*     ir=1 -> ir=5 */
/*     iw=1 -> iw=6 */
/*     All `entry' routines removed */
/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int driver_(void)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Builtin functions */
    double r_lg10(real *);

    /* Local variables */
    static integer i__;
#define v ((real *)&evolp1_1)
#define v0 ((real *)&evolp3_1)
    static real dtl;
#define dvdt ((real *)&evolp2_1)
    static integer mvar, loop;
    static real dvdt0[29];
    extern /* Subroutine */ int check_(void), accum_(void);
    static real dtmin;
    extern /* Subroutine */ int start_(void), derivs_(integer *);

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] run */
/*     CALLS     - [subroutine] start, derivs, accum */
/* ----------REMARKS. */
/*     Runge-Kutta computational routine */
/* ----------PARAMETERS. */
/* Number of variables to be evolved. */
/* Number of nuclides in calculation. */
/* ----------COMMON AREAS. */
/* Lower limit on size of time step. */
/* Evolution paramete */
/* Evolution paramete */
/* Evolution paramete */
/* Computation parame */
/* Time variables. */
/* Flags,counters. */
/* Computation locat */
/* ==========================DECLARATION DIVISION================================= */
/* ----------EVOLUTION PARAMETERS. */
/* Run options. */
/* Temperature (in units of 10**9 K). */
/* Defined by hv = M(atomic)n(baryon)/t */
/* Chemical potential for electron. */
/* ----------EVOLUTION PARAMETERS (DERIVATIVES). */
/* Relative number abundances. */
/* ----------EVOLUTION PARAMETERS (ORIGINAL VALUES). */
/* Change in rel number abundances. */
/* ----------COMPUTATION PARAMETERS. */
/* Rel # abund at beginning of iteratio */
/* Time step limiting constant on abund */
/* Time step limiting constant on tempe */
/* Final temperature (in 10**9 K). */
/* Smallest abundances allowed. */
/* ----------TIME AND TIME STEP VARIABLES. */
/* Accumulation increment. */
/* Time. */
/* Time step. */
/* ----------COUNTERS AND FLAGS. */
/* (1/t9)*d(t9)/d(t). */
/* Counts which Runge-Kutta loop. */
/* Indicates termination status. */
/* # total time steps for particular ru */
/* ----------COMPUTATION LOCATION. */
/* # time steps after outputting a line */
/* ----------RUN OPTION. */
/* Time check. */
/* ----------TIME AND TIME STEP VARIABLES. */
/* Number of nuclides in computation. */
/* Mininum time step. */
/* ----------LABELS FOR VARIABLES TO BE TIME EVOLVED. */
/* Time step from limitation on abund c */
/* Total number of variables to be evol */
/* Variables to be time evolved. */
/* Time derivatives. */
/* Value of variables at original point */
/* ----------EQUIVALENCE STATEMENTS. */
/* Value of derivatives at original poi */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------INPUT INITIALIZATION INFORMATION, RELABEL---------------------------- */
    flags_1.ltime = 0;
/* Set termination indicator to zero. */
    start_();
/* Input initialization information. */
    mvar = runopt_1.isize + 3;
/* 20--------LOOP ONE------------------------------------------------------------- */
/* Total number of variables to be evol */
L200:
/* Begin Runge-Kutta looping. */
    loop = 1;
/* ..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED. */
/* Loop indicator. */
    derivs_(&loop);
    tcheck_1.itime = 4;
/* Time = 1st R-K loop. */
    check_();
/* ..........ACCUMULATE. */
/* Check interface subroutine. */
    if (evolp1_1.t9 <= compr_1.t9f || ttime_1.dt < (r__1 = 1e-16f / 
	    ttime_1.dlt9dt, dabs(r__1)) || flags_1.ip == compr_1.inc) {
	accum_();
    }
/* ..........POSSIBLY TERMINATE COMPUTATION. */
/* Low temp. */
/* Small dt. */
/* Enough iterations. */
    if (flags_1.ltime == 1) {
/* Return to run selection. */
	return 0;
    }
/* ..........RESET COUNTERS. */
    if (flags_1.ip == compr_1.inc) {
/* Reset iteration counters. */
	flags_1.ip = 0;
    }
    ++flags_1.ip;
    ++flags_1.is;
/* ..........ADJUST TIME STEP. */
    if (flags_1.is > 3) {
/* Adjust time step after 3 iterations. */
	dtmin = (r__1 = 1.f / ttime_1.dlt9dt, dabs(r__1)) * compr_1.ct;
/* Trial value for minimum time step (R */
	i__1 = runopt_1.isize;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Go through all abundance changes. */
	    if (evolp2_1.dydt[i__ - 1] != 0.f && evolp1_1.y[i__ - 1] > 
		    compr_1.ytmin) {
/* Computing 2nd power */
		r__2 = r_lg10(&evolp1_1.y[i__ - 1]) / r_lg10(&compr_1.ytmin);
		dtl = (r__1 = evolp1_1.y[i__ - 1] / evolp2_1.dydt[i__ - 1], 
			dabs(r__1)) * compr_1.cy * (r__2 * r__2 + 1.f);
/* (Ref 2). */
		if (dtl < dtmin) {
		    dtmin = dtl;
		}
/* Find smallest time st */
	    }
	}
	if (dtmin > ttime_1.dt * 1.5f) {
	    dtmin = ttime_1.dt * 1.5f;
	}
/* Limit change in time */
	ttime_1.dt = dtmin;
/* Set new time step. */
    }
    ttime_1.t += ttime_1.dt;
/* ..........STORE AND INCREMENT VALUES (Ref 3). */
/* Increment time. */
    i__1 = mvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v0[i__ - 1] = v[i__ - 1];
	dvdt0[i__ - 1] = dvdt[i__ - 1];
	v[i__ - 1] = v0[i__ - 1] + dvdt0[i__ - 1] * ttime_1.dt;
	if (i__ >= 4 && v[i__ - 1] < compr_1.ytmin) {
	    v[i__ - 1] = compr_1.ytmin;
	}
/* Set at minimum */
    }
/* 30--------LOOP TWO------------------------------------------------------------- */
    loop = 2;
/* ..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED. */
/* Step up loop counter. */
    derivs_(&loop);
    tcheck_1.itime = 7;
/* Time = 2nd R-K loop. */
    check_();
/* ..........INCREMENT VALUES. */
/* Check interface subroutine. */
    i__1 = mvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ - 1] = v0[i__ - 1] + (dvdt[i__ - 1] + dvdt0[i__ - 1]) * .5f * 
		ttime_1.dt;
	if (i__ >= 4 && v[i__ - 1] < compr_1.ytmin) {
	    v[i__ - 1] = compr_1.ytmin;
	}
/* Set at minimum */
    }
    goto L200;
/* ----------REFERENCES----------------------------------------------------------- */
/*     1)  Constraint on dt from the requirement that */
/*                (d(t9)/dt)*(dt/t9) < ct */
/*           Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 293, equation C6. */
/*     2)  Constraint on dt from */
/*                dtl < y/(dy/dt)*cy*(1+(log(y)/log(ytmin))**2) */
/*          Wagoner, R.V. 1969, page 293, equation C7 but with log term squared. */
/*     3)  Wagoner, R.V. 1969, page 292, equations C1, C2. */
} /* driver_ */

#undef dvdt
#undef v0
#undef v


/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int start_(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static real z__;
    extern doublereal ex_(real *);
    extern /* Subroutine */ int rate0_(void), rate1_(real *), bessel_(real *);

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] driver */
/*     CALLS     - [subroutine] rate1, bessel, rate0 */
/*               - [function] ex */
/* ----------REMARKS. */
/*     Sets initial conditions. */
/* ----------PARAMETERS. */
/* Number of nuclear reactions. */
/* Number of nuclides in calculation. */
/* Relation between time and temperatur */
/* ----------COMMON AREAS. */
/* Gravitational constant. */
/* Reaction rates. */
/* Evolution paramete */
/* Evolution paramete */
/* Evolution paramete */
/* Computation parame */
/* Model parameters. */
/* Variational parame */
/* Time variables. */
/* Energy densities. */
/* Eval of function */
/* Eval of function b */
/* Eval of function b */
/* Flags,counters. */
/* Neutrino parameter */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION RATES. */
/* Run options. */
/* Forward reaction rate coefficients. */
/* ----------EVOLUTION PARAMETERS. */
/* Reverse reaction rate coefficients. */
/* Temperature (in units of 10**9 K). */
/* Defined by hv = M(atomic)n(baryon)/t */
/* Chemical potential of electron. */
/* ----------EVOLUTION PARAMETERS (ORIGINAL VALUES). */
/* Relative number abundances. */
/* ----------COMPUTATION SETTINGS. */
/* Rel # abund at start of iteration. */
/* Initial temperature (in 10**9 K). */
/* Smallest abundances allowed. */
/* ----------EARLY UNIVERSE MODEL PARAMETERS. */
/* Accumulation increment. */
/* Gravitational constant. */
/* Neutron lifetime. */
/* Number of neutrino species. */
/* c(1) is variation of grav. constant. */
/* c(2) is neutron lifetime (sec). */
/* c(3) is number of neutrino species. */
/* Neutrino degeneracy parameters. */
/* ----------VARIATIONAL PARAMETERS. */
/* Fierz parameter. */
/* Initial time step. */
/* ----------TIME VARIABLES. */
/* Baryon-to-photon ratio. */
/* Time. */
/* ----------ENERGY DENSITIES. */
/* Time step. */
/* Initial electron neutrino mass densi */
/* ----------EVALUATION OF FUNCTIONS bl,bm,bn. */
/* Initial baryon mass density. */
/* ----------COUNTERS AND FLAGS. */
/* Evaluation of function bl(z). */
/* Indicates if output buffer printed. */
/* # total time steps for particular ru */
/* # time steps after outputting a line */
/* # times accumulated in output buffer */
/* ----------NEUTRINO PARAMETERS. */
/* Indicates if gaussian elimination fa */
/* Neutrino temperature. */
/* ----------RUN OPTION. */
/* Normalizing constant. */
/* ----------LOCAL VARIABLES. */
/* Number of nuclides in computation. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------INITIALIZE FLAGS AND COUNTERS---------------------------------------- */
/* Defined by z = m(electron)*c**2/k*t9 */
    flags_1.ltime = 0;
/* No output yet. */
    flags_1.is = 1;
/* First iteration coming up. */
    flags_1.ip = compr_1.inc;
/* Set to maximum allowed # of iteratio */
    flags_1.it = 0;
/* No accumulation yet. */
    flags_1.mbad = 0;
/* 20--------SETTINGS------------------------------------------------------------- */
/* ..........COMPUTATIONAL SETTINGS. */
/* No computational errors. */
    evolp1_1.t9 = compr_1.t9i;
/* Initial temperature. */
    nupar_1.tnu = evolp1_1.t9;
/* Initial neutrino temperature. */
/* Computing 2nd power */
    r__1 = evolp1_1.t9 * .09615f;
    ttime_1.t = 1 / (r__1 * r__1);
/* Initial time (Ref 1). */
    ttime_1.dt = varpr_1.dt1;
/* ..........MODEL SETTINGS. */
/* Initial time step. */
    modpr_1.g = modpr_1.c__[0] * 6.67e-8f;
/* Modify gravitational constant. */
    modpr_1.tau = modpr_1.c__[1];
/* Convert n half-life (min) to lifetim */
    modpr_1.tau /= .98f;
/* Coulomb correction (Ref 2). */
    modpr_1.xnu = modpr_1.c__[2];
/* 30--------COMPUTE INITIAL ABUNDANCES FOR NEUTRON AND PROTON-------------------- */
/* Number of neutrino species. */
    if (15.011f / evolp1_1.t9 + modpr_1.xi[0] > 58.f) {
/* Overabundance of antineut */
	evolp1_1.y[0] = 1e-25f;
/* Very little of neutrons. */
	evolp1_1.y[1] = 1.f;
/* Essentially all protons. */
    } else {
	if (15.011f / evolp1_1.t9 + modpr_1.xi[0] < -58.f) {
/* Overabundance of neutrino */
	    evolp1_1.y[0] = 1.f;
/* Essentially all neutrons. */
	    evolp1_1.y[1] = 1e-25f;
/* Very little of protons. */
	} else {
	    r__1 = 15.011f / evolp1_1.t9 + modpr_1.xi[0];
	    evolp1_1.y[0] = 1.f / (ex_(&r__1) + 1.f);
/* Initial n abundance (Ref */
	    r__1 = -15.011f / evolp1_1.t9 - modpr_1.xi[0];
	    evolp1_1.y[1] = 1.f / (ex_(&r__1) + 1.f);
/* Initial p abundance (Ref */
	}
    }
/* change: set cnorm regardless of the nu_e chem potential */
/*      IF (xi(1).ne.0.) THEN        !Electron neutrino degeneracy. */
    nupar_1.cnorm = 1.f;
    nupar_1.tnu = 1e-5f;
/* Low temperature. */
    rate1_(&c_b4);
/* Find normalization constant at low t */
    nupar_1.cnorm = 1 / modpr_1.tau / rates_1.f[0];
/*      END IF */
    evolp3_1.y0[0] = evolp1_1.y[0];
    evolp3_1.y0[1] = evolp1_1.y[1];
/* 40--------FIND RATIO OF BARYON DENSITY TO TEMPERATURE CUBED-------------------- */
    z__ = 5.93f / evolp1_1.t9;
/* Inverse of temperature. */
    bessel_(&z__);
    evolp1_1.hv = varpr_1.eta1 * 33683.f * 2.75f;
/* (Ref 4 but with final eta). */
/* Computing 3rd power */
    r__1 = z__;
    evolp1_1.phie = evolp1_1.hv * (evolp1_1.y[1] * 1.784e-5f) / (r__1 * (r__1 
	    * r__1) * .5f * (besselcb_1.bl1 - besselcb_1.bl2 * 2.f + 
	    besselcb_1.bl3 * 3.f - besselcb_1.bl4 * 4.f + besselcb_1.bl5 * 
	    5.f));
/* Chemical potential of electron (Ref */
/* Computing 3rd power */
    r__1 = evolp1_1.t9;
    endens_1.rhob0 = evolp1_1.hv * (r__1 * (r__1 * r__1));
/* Baryon density. */
    if (modpr_1.xi[0] == 0.f && modpr_1.xi[1] == 0.f && modpr_1.xi[2] == 0.f) 
	    {
/* Nonde */
/* Computing 4th power */
	r__1 = evolp1_1.t9, r__1 *= r__1;
	endens_1.rhone0 = r__1 * r__1 * 7.366f;
/* Electron neutrino density (Ref 6). */
    }
/* 50--------SET ABUNDANCES FOR REST OF NUCLIDES---------------------------------- */
    r__1 = 25.82f / evolp1_1.t9;
    d__1 = (doublereal) evolp1_1.t9;
    evolp1_1.y[2] = evolp1_1.y[0] * evolp1_1.y[1] * endens_1.rhob0 * ex_(&
	    r__1) / (pow_dd(&d__1, &c_b5) * 4.71e9f);
/* (Ref 7). */
    evolp3_1.y0[2] = evolp1_1.y[2];
    i__1 = runopt_1.isize;
    for (i__ = 4; i__ <= i__1; ++i__) {
	evolp1_1.y[i__ - 1] = compr_1.ytmin;
/* Set rest to minimum abundance. */
	evolp3_1.y0[i__ - 1] = evolp1_1.y[i__ - 1];
/* Init abundances at beginning of iter */
    }
    rate0_();
/* Compute weak decay rates. */
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     1) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148, */
/*          page 44, equation A15. */
/*     2) Coulomb correction obtained by dividing by correction factor Fp(t9) */
/*               Fp(t9) = 1 - 0.5(pi/(137<v>/c)) */
/*          Wagoner, R.V. 1973, Ap. J. 179, page 358. */
/*     3) For the nondegenerate case: */
/*          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148, */
/*          page 4, equation 3. */
/*        For the case with neutrino degeneracy: */
/*          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, */
/*          page 417, equation 9. */
/*     4) Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 250, equation 4. */
/*          3.3683e+4 = Mu(ng/t9**3) with Mu the atomic mass, ng the */
/*          photon density.  2.75 is for the 11/4 factor difference */
/*          between the initial and final values of eta. */
/*     5) Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A, */
/*          Kellogg Radiation Lab preprint OAP-714. */
/*          equation D.2. */
/*     6) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148, */
/*          page 43, equation A4. */
/*          7.366 is used instead of 14.73 as the latter is the sum total */
/*          for 2 neutrino species. */
/*     7) Initial deuterium abundance from nuclear statistical equilibrium */
/*          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148, */
/*          page 19, equation 17. */
} /* start_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int derivs_(integer *loop)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static real bar;
    extern /* Subroutine */ int sol_(integer *);
    static real sumy;
    extern /* Subroutine */ int rate1_(real *), rate2_(void), rate3_(void), 
	    rate4_(void), therm_(void);
    static real sumdy, sumzy, dphdt9, dlndt9, dphdln, dphdzy, summdy, sumzdy;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] driver */
/*     CALLS     - [subroutine] therm, rate1, rate4, rate3, rate2, sol */
/* ----------REMARKS. */
/*     Computes derivatives of */
/*       - Temperature */
/*       - hv */
/*       - Chemical potential */
/*       - abundances */
/* ----------PARAMETERS. */
/* Number of variables to be evolved. */
/* Number of nuclides in calculation. */
/* ----------COMMON AREAS. */
/* Evolution paramete */
/* Evolution paramete */
/* Evolution paramete */
/* Model parameters. */
/* Time variables. */
/* Dynamic variables. */
/* Energy densities. */
/* Nuclide data. */
/* Flags,counters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------EVOLUTION PARAMETERS. */
/* Run options. */
/* Temperature (in units of 10**9 K). */
/* Defined by hv = M(atomic)n(baryon)/t */
/* Chemical potential for electron. */
/* ----------EVOLUTION PARAMETERS (DERIVATIVES). */
/* Relative number abundances. */
/* Change in temperature. */
/* Change in hv. */
/* Change in chemical potential. */
/* ----------EVOLUTION PARAMETERS (ORIGINAL VALUES). */
/* Change in rel number abundances. */
/* ----------MODEL PARAMETERS. */
/* Rel # abund at beginning of iteratio */
/* Gravitational constant. */
/* ----------TIME VARIABLES. */
/* Cosmological constant. */
/* ----------DYNAMIC VARIABLES. */
/* (1/t9)*d(t9)/d(t). */
/* Thermodynamic variables. */
/* ----------ENERGY DENSITIES. */
/* Expansion rate. */
/* Initial baryon mass density. */
/* Baryon mass density. */
/* ----------NUCLIDE DATA. */
/* Baryon mass density (ratio to init v */
/* Charge of nuclide. */
/* ----------COUNTERS AND FLAGS. */
/* Mass excess of nuclide. */
/* ----------RUN OPTION. */
/* Indicates if gaussian elimination fa */
/* Run network size. */
/* ----------SUMS. */
/* Number of nuclides in computation. */
/* Sum of abundances. */
/* Sum of charge*abundances. */
/* Sum of abundance flows. */
/* Sum of (mass excess)*(abundance flow */
/* ----------DERIVATIVES. */
/* Sum of (charge)*(abundance flows). */
/* d(phi e)/d(t9). */
/* d(phi e)/d(h). */
/* d(phi e)/d(sumzy). */
/* (1/h)*d(h)/d(t9). */
/* ----------LOCAL VARIABLES. */
/* Baryon density and pressure terms. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------COMPUTE DERIVATIVES FOR ABUNDANCES----------------------------------- */
/* Counts which Runge-Kutta loop. */
    endens_1.rnb = evolp1_1.hv * evolp1_1.t9 * evolp1_1.t9 * evolp1_1.t9 / 
	    endens_1.rhob0;
/* ..........VARIOUS THERMODYNAMIC QUANTITIES. */
/* Baryon mass density (ratio to init v */
    therm_();
    thermcb_1.hubcst = sqrt(modpr_1.g * 8.3775813333333318f * thermcb_1.thm[9]
	     + modpr_1.cosmo / 3.f);
/* Expansion rate. */
    endens_1.rhob = thermcb_1.thm[8];
/* ..........COMPUTE REACTION RATE COEFFICIENTS. */
/* Baryon mass density. */
    rate1_(&evolp1_1.t9);
    switch (runopt_1.irun) {
	case 1:  goto L100;
	case 2:  goto L110;
	case 3:  goto L120;
    }
/* Run network selection. */
L100:
    rate4_();
/* Forward rate for all of reactions. */
L110:
    rate3_();
/* Forward rate for reactions with A < */
L120:
    rate2_();
/* ..........SOLVE COUPLED DIFFERENTIAL EQUATIONS. */
/* Forward rate for reactions with A < */
    sol_(loop);
    if (flags_1.mbad > 0) {
	return 0;
    }
/* 20--------COMPUTE DERIVATIVES FOR TEMPERATURE, hv, AND CHEMICAL POTENTIAL------ */
/* ..........INITIALIZE SUMS TO ZERO. */
/* Abort in case matrix not invertible. */
    sumy = 0.f;
    sumzy = 0.f;
    sumdy = 0.f;
    summdy = 0.f;
    sumzdy = 0.f;
/* ..........ACCUMULATE TO GET SUM. */
    i__1 = runopt_1.isize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sumy += evolp1_1.y[i__ - 1];
/* Sum of abundance. */
	sumzy += nucdat_1.zm[i__ - 1] * evolp1_1.y[i__ - 1];
/* Sum of charge*abundance. */
	sumdy += evolp2_1.dydt[i__ - 1];
/* Sum of abundance flow. */
	summdy += nucdat_1.dm[i__ - 1] * evolp2_1.dydt[i__ - 1];
/* Sum of (mass excess)*(abundanc */
	sumzdy += nucdat_1.zm[i__ - 1] * evolp2_1.dydt[i__ - 1];
/* Sum of (charge)*(abundance flo */
    }
/* ..........CHANGES IN TEMPERATURE, hv, AND CHEMICAL POTENTIAL. */
    dphdt9 = thermcb_1.thm[11] * (evolp1_1.hv * -1.07e-4f * sumzy / 
	    evolp1_1.t9 - thermcb_1.thm[10]);
    dphdln = -thermcb_1.thm[11] * 3.568e-5f * evolp1_1.hv * sumzy;
    dphdzy = thermcb_1.thm[11] * 3.568e-5f * evolp1_1.hv;
    bar = evolp1_1.t9 * 9.25e-5f * sumy + evolp1_1.t9 * 1.388e-4f * sumdy / (
	    thermcb_1.hubcst * 3.f) + summdy / (thermcb_1.hubcst * 3.f);
    dlndt9 = -(thermcb_1.thm[1] + thermcb_1.thm[4] + thermcb_1.thm[5] * 
	    dphdt9 + thermcb_1.thm[8] * 1.388e-4f * sumy) / (thermcb_1.thm[0] 
	    + thermcb_1.thm[2] + thermcb_1.thm[3] + thermcb_1.thm[6] + 
	    thermcb_1.thm[8] * bar + thermcb_1.thm[5] * (dphdln + dphdzy * 
	    sumzdy / (thermcb_1.hubcst * 3.f)));
/* (Ref 1) */
    evolp2_1.dt9 = thermcb_1.hubcst * 3.f / dlndt9;
    ttime_1.dlt9dt = evolp2_1.dt9 / evolp1_1.t9;
    evolp2_1.dhv = -evolp1_1.hv * (thermcb_1.hubcst * 3.f + ttime_1.dlt9dt * 
	    3.f);
/* (Ref 2) */
    evolp2_1.dphie = dphdt9 * evolp2_1.dt9 + dphdln * (thermcb_1.hubcst * 3.f)
	     + dphdzy * sumzdy;
/* (Ref 3) */
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     1)  Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A, */
/*          Kellogg Radiation Lab preprint OAP-714, */
/*          equation D.35. */
/*     2)  Kawano, L., 1992, preprint, equation D.19. */
/*     3)  Kawano, L., 1992, preprint, equation D.20. */
} /* derivs_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int accum_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] driver */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Output accumulator. */
/* ----------PARAMETERS. */
/* Number of variables to be evolved. */
/* Number of nuclides in calculation. */
/* ----------COMMON AREAS. */
/* Maximum # of lines to be printed. */
/* Evolution paramete */
/* Computation parame */
/* Time variables. */
/* Dynamic variables. */
/* Nuclide data. */
/* Flags,counters. */
/* Output data. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------EVOLUTION PARAMETERS. */
/* Run options. */
/* Temperature (in units of 10**9 K). */
/* Defined by hv = M(atomic)n(baryon)/t */
/* Chemical potential for electron. */
/* ----------COMPUTATION PARAMETERS. */
/* Relative number abundances. */
/* ----------TIME PARAMETERS. */
/* Accumulation increment. */
/* Time. */
/* ----------DYNAMIC VARIABLES. */
/* Time step. */
/* Thermodynamic variables. */
/* ----------NUCLIDE DATA. */
/* Expansion rate. */
/* ----------COUNTERS AND FLAGS. */
/* Atomic number of nuclide. */
/* Indicates if output buffer printed. */
/* # times accumulated in output buffer */
/* ----------OUTPUT ARRAYS. */
/* # time steps after outputting a line */
/* Nuclide mass fractions. */
/* Thermodynamic variables. */
/* Temperature (in units of 10**9 K). */
/* Time. */
/* Time step. */
/* Baryon-to-photon ratio. */
/* ----------RUN OPTION. */
/* Expansion rate. */
/* ===========================PROCEDURE DIVISION================================== */
/* Number of nuclides in computation. */
    ++flags_1.it;
/* 10--------SET UP OUTPUT VARIABLES---------------------------------------------- */
/* ..........DIVIDE NUMBER FRACTION BY THAT OF PROTON. */
/* Set up accumulation counter. */
    i__1 = runopt_1.isize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	outdat_1.xout[flags_1.it + i__ * 40 - 41] = evolp1_1.y[i__ - 1] / 
		evolp1_1.y[1];
    }
    outdat_1.xout[flags_1.it + 39] = evolp1_1.y[1] * nucdat_1.am[1];
/* Exception for proton. */
    outdat_1.xout[flags_1.it + 199] = evolp1_1.y[5] * nucdat_1.am[5];
/* ..........SUM UP ABUNDANCES OF HEAVY NUCLIDES. */
/* Exception for helium. */
    outdat_1.xout[flags_1.it + 359] = outdat_1.xout[flags_1.it + 359] + 
	    outdat_1.xout[flags_1.it + 399] + outdat_1.xout[flags_1.it + 439] 
	    + outdat_1.xout[flags_1.it + 479] + outdat_1.xout[flags_1.it + 
	    519] + outdat_1.xout[flags_1.it + 559] + outdat_1.xout[flags_1.it 
	    + 599] + outdat_1.xout[flags_1.it + 639] + outdat_1.xout[
	    flags_1.it + 679] + outdat_1.xout[flags_1.it + 719] + 
	    outdat_1.xout[flags_1.it + 759] + outdat_1.xout[flags_1.it + 799] 
	    + outdat_1.xout[flags_1.it + 839] + outdat_1.xout[flags_1.it + 
	    879] + outdat_1.xout[flags_1.it + 919] + outdat_1.xout[flags_1.it 
	    + 959] + outdat_1.xout[flags_1.it + 999];
/* ..........RELABEL TEMPERATURE, TIME, THERMODYNAMIC VARIABLES, ETC. */
/* Li8 to O16. */
    outdat_1.t9out[flags_1.it - 1] = evolp1_1.t9;
/* Temperature. */
    outdat_1.tout[flags_1.it - 1] = ttime_1.t;
/* Time. */
    outdat_1.thmout[flags_1.it - 1] = thermcb_1.thm[0];
/* rho photon. */
    outdat_1.thmout[flags_1.it + 39] = thermcb_1.thm[3];
/* rho electron. */
    outdat_1.thmout[flags_1.it + 79] = thermcb_1.thm[7];
/* rho neutrino. */
    outdat_1.thmout[flags_1.it + 119] = thermcb_1.thm[8];
/* rho baryon. */
    outdat_1.thmout[flags_1.it + 159] = evolp1_1.phie;
/* Chemical potential. */
    outdat_1.thmout[flags_1.it + 199] = thermcb_1.thm[9];
/* rho total. */
    outdat_1.dtout[flags_1.it - 1] = ttime_1.dt;
/* Time step. */
    outdat_1.etaout[flags_1.it - 1] = evolp1_1.hv / 33683.f;
/* Baryon to photon ratio. */
    outdat_1.hubout[flags_1.it - 1] = thermcb_1.hubcst;
/* 20--------INDICATE TERMINATION OF ACCUMULATION IF APPROPRIATE------------------ */
/* Expansion rate. */
    if (flags_1.it == 40 || flags_1.ip < compr_1.inc) {
	flags_1.ltime = 1;
    }
    return 0;
} /* accum_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int therm_(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), cosh(doublereal), sinh(
	    doublereal);

    /* Local variables */
    static real z__, z1, z2, z3, z4, z5;
    extern doublereal ex_(real *);
    static real cosh1, cosh2, cosh3, cosh4, cosh5, sinh1, sinh2, sinh3, sinh4,
	     sinh5;
    extern /* Subroutine */ int bessel_(real *), nudens_(void);

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] derivs */
/*     CALLS     - [subroutine] bessel, nudens */
/*               - [function] ex */
/* ----------REMARKS. */
/*     Computes various temperature dependent thermodynamic quantities. */
/* ----------PARAMETER. */
/* Number of nuclides in calculation. */
/* ----------COMMON AREAS. */
/* (mass(neutron)-mass(proton))/m(elect */
/* Evolution paramete */
/* Computation parame */
/* Model parameters. */
/* Dynamic variables. */
/* Energy densities. */
/* Eval of function b */
/* Eval of function b */
/* Eval of function b */
/* ==========================DECLARATION DIVISION================================= */
/* ----------EVOLUTION PARAMETERS. */
/* Integration parame */
/* Temperature (in units of 10**9 K). */
/* ----------COMPUTATION PARAMETERS. */
/* Chemical potential for electron. */
/* ----------EARLY UNIVERSE MODEL PARAMETERS. */
/* Initial temperature (in 10**9 K). */
/* Number of neutrino species. */
/* ----------DYNAMIC VARIABLES. */
/* Neutrino degeneracy parameters. */
/* ----------ENERGY DENSITIES. */
/* Thermodynamic variables. */
/* Initial electron neutrino mass densi */
/* Initial baryon mass density. */
/* ----------EVALUATION OF FUNCTIONS bl,bm,bn. */
/* Baryon mass density (ratio to init v */
/* Evaluation of function bl(z). */
/* Evaluation of function bm(z). */
/* ----------NEUTRINO PARAMETERS. */
/* Evaluation of function bn(z). */
/* Neutrino temperature. */
/* Neutrino energy density. */
/* ----------LOCAL VARIABLE. */
/* Type of neutrino. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------COMPUTE FACTORS------------------------------------------------------ */
/* Defined by z = m(electron)*c**2/k*t9 */
    z__ = 5.93f / evolp1_1.t9;
/* z = m(electron)c**2/k(t9). */
    d__1 = (doublereal) endens_1.rnb;
    nupar_1.tnu = pow_dd(&d__1, &c_b12) * compr_1.t9i;
/* ..........FACTORS OF z. */
/* Neutrino temperature. */
    z1 = z__;
    z2 = z__ * z__;
    z3 = z__ * z__ * z__;
    z4 = z__ * z__ * z__ * z__;
    z5 = z__ * z__ * z__ * z__ * z__;
/* ..........TRIGNOMETRIC FUNCTION VALUES. */
    if (evolp1_1.phie <= 17.f) {
/* No chance of overflow. */
	cosh1 = cosh(evolp1_1.phie);
	cosh2 = cosh(evolp1_1.phie * 2.f);
	cosh3 = cosh(evolp1_1.phie * 3.f);
	cosh4 = cosh(evolp1_1.phie * 4.f);
	cosh5 = cosh(evolp1_1.phie * 5.f);
	sinh1 = sinh(evolp1_1.phie);
	sinh2 = sinh(evolp1_1.phie * 2.f);
	sinh3 = sinh(evolp1_1.phie * 3.f);
	sinh4 = sinh(evolp1_1.phie * 4.f);
	sinh5 = sinh(evolp1_1.phie * 5.f);
    } else {
	cosh1 = 0.f;
	cosh2 = 0.f;
	cosh3 = 0.f;
	cosh4 = 0.f;
	cosh5 = 0.f;
	sinh1 = 0.f;
	sinh2 = 0.f;
	sinh3 = 0.f;
	sinh4 = 0.f;
	sinh5 = 0.f;
    }
    bessel_(&z__);
/* 20--------COMPUTE THERMODYNAMIC VARIABLES-------------------------------------- */
    thermcb_1.thm[0] = evolp1_1.t9 * 8.418f * evolp1_1.t9 * evolp1_1.t9 * 
	    evolp1_1.t9;
/* (Ref 1) */
    thermcb_1.thm[1] = thermcb_1.thm[0] * 4.f / evolp1_1.t9;
/* (Ref 2) */
    thermcb_1.thm[2] = thermcb_1.thm[0] / 3.f;
/* (Ref 3) */
    thermcb_1.thm[3] = (besselcb_1.bm1 * cosh1 - besselcb_1.bm2 * cosh2 + 
	    besselcb_1.bm3 * cosh3 - besselcb_1.bm4 * cosh4 + besselcb_1.bm5 *
	     cosh5) * 3206.f;
/* (Ref 4) */
    thermcb_1.thm[4] = z__ / evolp1_1.t9 * 3206.f * (besselcb_1.bn1 * cosh1 - 
	    besselcb_1.bn2 * 2.f * cosh2 + besselcb_1.bn3 * 3.f * cosh3 - 
	    besselcb_1.bn4 * 4.f * cosh4 + besselcb_1.bn5 * 5.f * cosh5);
/* (Ref 5) */
    thermcb_1.thm[5] = (besselcb_1.bm1 * sinh1 - besselcb_1.bm2 * 2.f * sinh2 
	    + besselcb_1.bm3 * 3.f * sinh3 - besselcb_1.bm4 * 4.f * sinh4 + 
	    besselcb_1.bm5 * 5.f * sinh5) * 3206.f;
/* (Ref 6) */
    thermcb_1.thm[6] = (besselcb_1.bl1 * cosh1 / z__ - besselcb_1.bl2 * cosh2 
	    / (z__ * 2.f) + besselcb_1.bl3 * cosh3 / (z__ * 3.f) - 
	    besselcb_1.bl4 * cosh4 / (z__ * 4.f) + besselcb_1.bl5 * cosh5 / (
	    z__ * 5.f)) * 3206.f;
/* (Ref 7) */
    if (modpr_1.xi[0] == 0.f && modpr_1.xi[1] == 0.f && modpr_1.xi[2] == 0.f) 
	    {
/* Nonde */
	d__1 = (doublereal) endens_1.rnb;
	thermcb_1.thm[7] = modpr_1.xnu * endens_1.rhone0 * pow_dd(&d__1, &
		c_b13);
/* (Ref 8) */
    } else {
/* Include effects of neutrino degenera */
	thermcb_1.thm[7] = 0.f;
	i__1 = modpr_1.xnu;
	for (nupar_1.nu = 1; nupar_1.nu <= i__1; ++nupar_1.nu) {
/* For every neutrino family. */
	    nudens_();
/* Compute neutrino energy density. */
	    thermcb_1.thm[7] += nupar_1.rhonu * 12.79264f;
/* Have 12.79264 from units ch */
	}
    }
    thermcb_1.thm[8] = endens_1.rhob0 * endens_1.rnb;
/* (Ref 9) */
    thermcb_1.thm[9] = thermcb_1.thm[0] + thermcb_1.thm[3] + thermcb_1.thm[7] 
	    + thermcb_1.thm[8];
/* (Ref 10 */
/* Computing 3rd power */
    r__1 = z__;
    thermcb_1.thm[10] = -(r__1 * (r__1 * r__1) / evolp1_1.t9) * (sinh1 * (
	    besselcb_1.bl1 * 3.f - z__ * besselcb_1.bm1) - sinh2 * (
	    besselcb_1.bl2 * 3.f - z__ * 2.f * besselcb_1.bm2) + sinh3 * (
	    besselcb_1.bl3 * 3.f - z__ * 3.f * besselcb_1.bm3) - sinh4 * (
	    besselcb_1.bl4 * 3.f - z__ * 4.f * besselcb_1.bm4) + sinh5 * (
	    besselcb_1.bl5 * 3.f - z__ * 5.f * besselcb_1.bm5));
/* (Ref 11 */
/* Computing 3rd power */
    r__1 = z__;
    thermcb_1.thm[11] = r__1 * (r__1 * r__1) * (cosh1 * besselcb_1.bl1 - 
	    cosh2 * 2.f * besselcb_1.bl2 + cosh3 * 3.f * besselcb_1.bl3 - 
	    cosh4 * 4.f * besselcb_1.bl4 + cosh5 * 5.f * besselcb_1.bl5);
/* (Ref 12 */
    if (thermcb_1.thm[11] != 0.f) {
	thermcb_1.thm[11] = 1.f / thermcb_1.thm[11];
    }
    thermcb_1.thm[12] = .565f / z1 + 1.f - 6.382f / z2 + 11.108f / z3 + 
	    36.492f / z4 + 27.512f / z5;
/* (Ref 13 */
    r__1 = z__ * -2.531f;
    thermcb_1.thm[13] = (5.252f / z1 - 16.229f / z2 + 18.059f / z3 + 34.181f /
	     z4 + 27.617f / z5) * ex_(&r__1);
/* (Ref 14 */
    return 0;
/* ----------REFERENCES AND NOTES------------------------------------------------- */
/*     1)  thm(1)  = rho photon */
/*         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148, */
/*          page 43, equation A2.) */
/*     2)  thm(2)  = d(rho photon)/d(t9) */
/*     3)  thm(3)  = (p photon)/c**2 */
/*         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, */
/*          page 43, equation A3.) */
/*     4)  thm(4)  = rho electron+positron */
/*         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9, */
/*          page 281, equation B44.) */
/*     5)  thm(5)  = d(rho electron+positron)/d(t9) */
/*     6)  thm(6)  = d(rho electron+positron)/d(phi e) */
/*     7)  thm(7)  = (p electron+positron)/c**2 */
/*         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9, */
/*          page 279, equation B27.) */
/*     8)  thm(8)  = rho neutrino */
/*                 = # neutrino species x rho electron neutrino (nondegenerate) */
/*                 = rho nu(e) + rho nu(m) + rho nu(t)          (degenerate) */
/*     9)  thm(9)  = rho baryon */
/*     10) thm(10) = rho total */
/*                 = rho photon + rho electron+positron + rho neutrino */
/*                              + rho baryon */
/*     11) thm(11) = d     /pi**2(hbar*c)**3(ne- - ne+)*z**3\ */
/*                   d(t9) \  2  (mc**2)**3                 / */
/*     12) thm(12) = d        /pi**2(hbar*c)**3(ne- - ne+)*z**3\ */
/*                   d(phi e) \  2  (mc**2)**3                 / */
/*     13) thm(13) = rate for n->p */
/*     14) thm(14) = rate for p->n */
} /* therm_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int bessel_(real *z__)
{
    /* Local variables */
    static integer i__;
    static real r__;
#define blz ((real *)&besselcb_1)
#define bmz ((real *)&besselcb_1 + 5)
#define bnz ((real *)&besselcb_1 + 10)
    extern /* Subroutine */ int knux_(real *);

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] start, therm */
/*     CALLS     - [subroutine] knux */
/* ----------REMARKS. */
/*     Evaluates functions bl(z), bm(z), and bn(z) using solutions to */
/*     modified Bessel functions. */
/* ----------COMMON AREAS. */
/* Eval function bl(z */
/* Eval function bm(z */
/* Eval function bn(z */
/* ==========================DECLARATION DIVISION================================= */
/* ----------EVALUATION OF FUNCTIONS bl,bm,bn. */
/* Coefficients K. */
/* Single variables equivalenced to arr */
/* Single variables equivalenced to arr */
/* ----------EVALUATIION OF MODIFIED BESSEL FUNCTIONS. */
/* Single variables equivalenced to arr */
/* ----------LOCAL VARIABLES. */
/* Values k0(r),k1(r),k2(r),k3(r),k4(r) */
/* Array containing values from functio */
/* Array containing values from functio */
/* Array containing values from functio */
/* Defined by z = m(electron)*c**2/k*t9 */
/* ----------EQUIVALENCE STATEMENTS. */
/* Multiples of z. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------LOCALLY DEFINED FUNCTIONS-------------------------------------------- */
/* Function bl. */
/* Function bm. */
/* 20--------CALCULATE FOR 1 THRU 5 Z--------------------------------------------- */
/* Function bn. */
    for (i__ = 1; i__ <= 5; ++i__) {
	r__ = i__ * *z__;
/* Multiples of z. */
	knux_(&r__);
/* Get k0(r),k1(r),k2(r),k3(r),k4(r),k( */
	blz[i__ - 1] = kays_1.bk2 / r__;
/* Put value from function bl into arra */
	bmz[i__ - 1] = .25f * (3.f * kays_1.bk3 + kays_1.bk1) / r__;
/* Put value from function bm into arra */
	bnz[i__ - 1] = .5f * (kays_1.bk4 + kays_1.bk2) / r__;
/* Put value from function bn into arra */
    }
    return 0;
} /* bessel_ */

#undef bnz
#undef bmz
#undef blz


/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int knux_(real *z__)
{
    /* Initialized data */

    static real ci0[7] = { 1.f,3.5156229f,3.0899424f,1.2067492f,.2659732f,
	    .0360768f,.0045813f };
    static real ci1[7] = { .5f,.87890594f,.51498869f,.15084934f,.02658733f,
	    .00301532f,3.2411e-4f };
    static real ck0[7] = { -.57721566f,.4227842f,.23069756f,.0348859f,
	    .00262698f,1.075e-4f,7.4e-6f };
    static real ck1[7] = { 1.f,.15443144f,-.67278579f,-.18156897f,-.01919402f,
	    -.00110404f,-4.686e-5f };
    static real c0[7] = { 1.25331414f,-.07832358f,.02189568f,-.01062446f,
	    .00587872f,-.0025154f,5.3208e-4f };
    static real c1[7] = { 1.25331414f,.23498619f,-.0365562f,.01504268f,
	    -.00780353f,.00325614f,-6.8245e-4f };

    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    double log(doublereal), pow_ri(real *, integer *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static real t, y;
    extern doublereal ex_(real *);
    static real bi0, bi1, coeff;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] bessel */
/*     CALLS     - [function] exp */
/* ----------REMARKS. */
/*     A subroutine for modified bessel functions of the second kind */
/*     k-nu(z). */
/* ----------COMMON AREAS. */
/* ===========================DECLARATION DIVISION================================ */
/* -----------MODIFIED BESSEL FUNCTION VALUES. */
/* Coefficients K. */
/* Values k0(z),k1(z) */
/* Values i0(z),i1(z). */
/* -----------EXPANSION COEFFICIENTS. */
/* Values k2(z),k3(z),k4(z). */
/* Expansion coefficients for i0 (z.le. */
/* Expansion coefficients for i1 (z.le. */
/* Expansion coefficients for k0 (z.le. */
/* Expansion coefficients for k1 (z.le. */
/* Expansion coefficients for k0 (z.gt. */
/* -----------VARIABLES TO BE EVALUATED. */
/* Expansion coefficients for k1 (z.gt. */
/* Input variable. */
/* Expansion variable = z/2. */
/* Expansion variable = z/3.75. */
/* ==============================DATA DIVISION==================================== */
/* ----------EXPANSION COEFFICIENTS. */
/* Logrithmic or exponential coefficien */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------COMPUTE K0 AND K1---------------------------------------------------- */
    if (*z__ <= 2.f) {
/* ..........COMPUTE FACTORS. */
/* (Ref. 1). */
	t = *z__ / 3.75f;
	y = *z__ / 2;
	coeff = log(y);
/* ..........VALUES FOR i0(z) and i1(z). */
	bi0 = ci0[0];
	bi1 = ci1[0];
	kays_1.bk0 = ck0[0];
	kays_1.bk1 = ck1[0];
	for (i__ = 2; i__ <= 7; ++i__) {
	    i__1 = i__ - 1 << 1;
	    bi0 += ci0[i__ - 1] * pow_ri(&t, &i__1);
	    i__1 = i__ - 1 << 1;
	    bi1 += ci1[i__ - 1] * pow_ri(&t, &i__1);
	    i__1 = i__ - 1 << 1;
	    kays_1.bk0 += ck0[i__ - 1] * pow_ri(&y, &i__1);
	    i__1 = i__ - 1 << 1;
	    kays_1.bk1 += ck1[i__ - 1] * pow_ri(&y, &i__1);
	}
/* ..........VALUES FOR k0(z) and k1(z). */
	kays_1.bk0 = -coeff * bi0 + kays_1.bk0;
	kays_1.bk1 = coeff * bi1 * *z__ + kays_1.bk1 / *z__;
    } else {
/* ..........COMPUTE FACTORS. */
/* (z.le.2.)               !(Ref. 2). */
	y = 2.f / *z__;
	r__1 = -(*z__);
	coeff = ex_(&r__1) / sqrt(*z__);
/* ..........VALUES FOR k0(z) and k1(z). */
	kays_1.bk0 = c0[0];
	kays_1.bk1 = c1[0];
	for (i__ = 2; i__ <= 7; ++i__) {
	    i__1 = i__ - 1;
	    kays_1.bk0 += c0[i__ - 1] * pow_ri(&y, &i__1);
	    i__1 = i__ - 1;
	    kays_1.bk1 += c1[i__ - 1] * pow_ri(&y, &i__1);
	}
	kays_1.bk0 = coeff * kays_1.bk0;
	kays_1.bk1 = coeff * kays_1.bk1;
    }
/* 20--------FIND K2, K3, AND K4 BY ITERATION (Ref. 3)---------------------------- */
/* (z.le.2.) */
    kays_1.bk2 = kays_1.bk1 / *z__ * 2.f + kays_1.bk0;
/* k2(z). */
    kays_1.bk3 = kays_1.bk2 / *z__ * 4.f + kays_1.bk1;
/* k3(z). */
    kays_1.bk4 = kays_1.bk3 / *z__ * 6.f + kays_1.bk2;
/* k4(z). */
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     Handbook of Mathematical Functions (Abramowitz and Stegun), */
/*       Dover Publications, Inc., New York */
/*       1) Polynomial approximations for z.le.2 */
/*         page 378, equations 9.8.1 and 9.8.3. */
/*         page 379, equations 9.8.5 and 9.8.7. */
/*       2) Polynomial approximations for z > 2 */
/*         page 379, equations 9.8.6 and 9.8.8. */
/*       3) Recursion relation from 1st line of 9.6.26, page 376. */
} /* knux_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int nudens_(void)
{
    /* System generated locals */
    real r__1, r__2, r__3;

    /* Local variables */
    //TODO errase extern /* Subroutine */ int func5_(), func6_();
    extern /* Subroutine */ doublereal func5_(), func6_();
    //extern doublereal xintd_(real *, real *, U_fp, integer *);
    extern doublereal xintd_(real *, real *, E_fp, integer *);
    static real uplim1, uplim2;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] therm */
/*     CALLS     - [function] xintd, eval */
/* ----------REMARKS. */
/*     Computes energy density contribution from neutrinos. */
/* ----------PARAMTER. */
/* ----------COMMON AREAS. */
/* Number of gaussian quads. */
/* Model parameters. */
/* ----------EXTERNAL FUNCTIONS. */
/* Integration parame */
/* Integral for neutrinos. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------MODEL PARAMETERS. */
/* Integral for antineutrinos. */
/* ----------NEUTRINO PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* Neutrino temperature (units of 10**9 */
/* Neutrino energy density. */
/* ----------LOCAL VARIABLES. */
/* Which neutrino type. */
/* Upper limit for neutrino energy inte */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------COMPUTE NEUTRINO ENERGY DENSITIES------------------------------------ */
/* Upper limit for antineu energy integ */
    if ((r__1 = modpr_1.xi[nupar_1.nu - 1], dabs(r__1)) <= .03f) {
/* ..........SMALL xi APPROXIMATION. */
/* Computing 4th power */
	r__1 = nupar_1.tnu, r__1 *= r__1;
/* Computing 2nd power */
	r__2 = modpr_1.xi[nupar_1.nu - 1];
/* Computing 4th power */
	r__3 = modpr_1.xi[nupar_1.nu - 1], r__3 *= r__3;
	nupar_1.rhonu = r__1 * r__1 * .65797251520666655f * (r__2 * r__2 * 
		.37995508052714932f + .875f + r__3 * r__3 * 
		.01924878176245234f);
    } else {
	if ((r__1 = modpr_1.xi[nupar_1.nu - 1], dabs(r__1)) >= 30.f) {
/* ..........LARGE xi APPROXIMATION. */
/* Computing 4th power */
	    r__1 = nupar_1.tnu, r__1 *= r__1;
/* Computing 4th power */
	    r__2 = modpr_1.xi[nupar_1.nu - 1], r__2 *= r__2;
/* Computing 2nd power */
	    r__3 = modpr_1.xi[nupar_1.nu - 1];
	    nupar_1.rhonu = r__1 * r__1 / 78.956701824799993f * (r__2 * r__2) 
		    * (19.740000000000002f / (r__3 * r__3) + 1);
	} else {
/* ..........DO INTEGRATION */
	    uplim1 = (modpr_1.xi[nupar_1.nu - 1] + 88.029f) * nupar_1.tnu;
	    uplim2 = (88.029f - modpr_1.xi[nupar_1.nu - 1]) * nupar_1.tnu;
	    if (uplim2 <= 0.f) {
		nupar_1.rhonu = xintd_(&c_b17, &uplim1, (E_fp)func5_, &c__50);
	    } else {
		nupar_1.rhonu = xintd_(&c_b17, &uplim1, (E_fp)func5_, &c__50) 
			+ xintd_(&c_b17, &uplim2, (E_fp)func6_, &c__50);
	    }
	}
/* (abs(xi(nu)).ge.30.) */
    }
/* (abs(xi(nu)).le.0.03) */
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     Forms of the integrals involved can be found in */
/*       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415. */
/*       Freese, K., Kolb, E.W., Turner, M.S., 1983, Phys. Rev. D, 27, 1689. */
} /* nudens_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------1ST PART OF INTEGRAL FOR n->p RATE----------------------------------- */
/* *************************************************************************** */
doublereal func1_(real *x)
{
    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal ex_(real *);
    static real part1, part2;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] rate1, nudens */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Contains integrands to be integrated. */
/* ----------COMMON AREAS. */
/* Model parameters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------MODEL PARAMETERS. */
/* Integration parame */
/* ----------NEUTRINO PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* Temperature (in units of MeV). */
/* Neutrino temperature (in units of Me */
/* Neutrino temperature (units of 10**9 */
/* ----------LOCAL VARIABLES. */
/* Normalizing constant. */
/* Value at which function is evaluated */
/* Exponential expression with photon t */
/* Exponential expression with neutrino */
    if (*x <= 0.f) {
	ret_val = 0.f;
    } else {
	r__1 = *x * -.511f / nupar_1.t9mev;
	part1 = 1.f / (ex_(&r__1) + 1.f);
	r__1 = (*x - 2.531f) * (.511f / nupar_1.tnmev) - modpr_1.xi[0];
	part2 = 1.f / (ex_(&r__1) + 1.f);
/*       func1 = cnorm*x*(x-2.531)**2*(x**2-1)**.5*part1*part2 */
/*      IF (x.gt.(2.531) THEN */
/* Computing 2nd power */
	r__1 = *x - 2.531f;
/* Computing 2nd power */
	r__2 = *x;
	d__1 = (doublereal) (r__2 * r__2 - 1);
	ret_val = nupar_1.cnorm * (*x + modpr_1.b) * (r__1 * r__1) * pow_dd(&
		d__1, &c_b24) * part1 * part2;
/*      ELSE */
/*        func1 = cnorm*(x-b)*(x-2.531)**2*(x**2-1)**.5*part1*part2 */
/*      END IF */
    }
    return ret_val;
} /* func1_ */

/* 20--------2ND PART OF INTEGRAL FOR n->p RATE----------------------------------- */
doublereal func2_(real *x)
{
    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal ex_(real *);
    static real part1, part2;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] rate1, nudens */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Contains integrands to be integrated. */
/* ----------COMMON AREAS. */
/* Model parameters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------MODEL PARAMETERS. */
/* Integration parame */
/* ----------NEUTRINO PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* Temperature (in units of MeV). */
/* Neutrino temperature (in units of Me */
/* Neutrino temperature (units of 10**9 */
/* ----------LOCAL VARIABLES. */
/* Normalizing constant. */
/* Value at which function is evaluated */
/* Exponential expression with photon t */
/* Exponential expression with neutrino */
    if (*x <= 1.f) {
	ret_val = 0.f;
    } else {
	r__1 = *x * .511f / nupar_1.t9mev;
	part1 = 1.f / (ex_(&r__1) + 1.f);
	r__1 = -(*x + 2.531f) * (.511f / nupar_1.tnmev) - modpr_1.xi[0];
	part2 = 1.f / (ex_(&r__1) + 1.f);
/*       func2 = cnorm*(x+b*)*(x+2.531)**2*(x**2-1)**.5*part1*part2 */
/*      IF (x.gt.(2.531) THEN */
/* Computing 2nd power */
	r__1 = *x + 2.531f;
/* Computing 2nd power */
	r__2 = *x;
	d__1 = (doublereal) (r__2 * r__2 - 1);
	ret_val = nupar_1.cnorm * (*x - modpr_1.b) * (r__1 * r__1) * pow_dd(&
		d__1, &c_b24) * part1 * part2;
/*      ELSE */
/*        func2 = cnorm*(x+b)*(x+2.531)**2*(x**2-1)**.5*part1*part2 */
/*      END IF */
    }
    return ret_val;
} /* func2_ */

/* 30--------1ST PART OF INTEGRAL FOR p->n RATE----------------------------------- */
doublereal func3_(real *x)
{
    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal ex_(real *);
    static real part1, part2;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] rate1, nudens */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Contains integrands to be integrated. */
/* ----------COMMON AREAS. */
/* Model parameters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------MODEL PARAMETERS. */
/* Integration parame */
/* ----------NEUTRINO PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* Temperature (in units of MeV). */
/* Neutrino temperature (in units of Me */
/* Neutrino temperature (units of 10**9 */
/* ----------LOCAL VARIABLES. */
/* Normalizing constant. */
/* Value at which function is evaluated */
/* Exponential expression with photon t */
/* Exponential expression with neutrino */
    if (*x <= 1.f) {
	ret_val = 0.f;
    } else {
	r__1 = *x * -.511f / nupar_1.t9mev;
	part1 = 1.f / (ex_(&r__1) + 1.f);
	r__1 = (*x + 2.531f) * (.511f / nupar_1.tnmev) + modpr_1.xi[0];
	part2 = 1.f / (ex_(&r__1) + 1.f);
/*       func3 = cnorm*x*(x+2.531)**2*(x**2-1)**.5*part1*part2 */
/* Computing 2nd power */
	r__1 = *x + 2.531f;
/* Computing 2nd power */
	r__2 = *x;
	d__1 = (doublereal) (r__2 * r__2 - 1);
	ret_val = nupar_1.cnorm * (*x - modpr_1.b) * (r__1 * r__1) * pow_dd(&
		d__1, &c_b24) * part1 * part2;
    }
    return ret_val;
} /* func3_ */

/* 40--------2ND PART OF INTEGRAL FOR p->n RATE----------------------------------- */
doublereal func4_(real *x)
{
    /* System generated locals */
    real ret_val, r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal ex_(real *);
    static real part1, part2;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] rate1, nudens */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Contains integrands to be integrated. */
/* ----------COMMON AREAS. */
/* Model parameters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------MODEL PARAMETERS. */
/* Integration parame */
/* ----------NEUTRINO PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* Temperature (in units of MeV). */
/* Neutrino temperature (in units of Me */
/* Neutrino temperature (units of 10**9 */
/* ----------LOCAL VARIABLES. */
/* Normalizing constant. */
/* Value at which function is evaluated */
/* Exponential expression with photon t */
/* Exponential expression with neutrino */
    if (*x <= 1.f) {
	ret_val = 0.f;
    } else {
	r__1 = *x * .511f / nupar_1.t9mev;
	part1 = 1.f / (ex_(&r__1) + 1.f);
	r__1 = -(*x - 2.531f) * (.511f / nupar_1.tnmev) + modpr_1.xi[0];
	part2 = 1.f / (ex_(&r__1) + 1.f);
/*       func4 = cnorm*x*(x-2.531)**2*(x**2-1)**.5*part1*part2 */
/* Computing 2nd power */
	r__1 = *x - 2.531f;
/* Computing 2nd power */
	r__2 = *x;
	d__1 = (doublereal) (r__2 * r__2 - 1);
	ret_val = nupar_1.cnorm * (*x + modpr_1.b * .511f) * (r__1 * r__1) * 
		pow_dd(&d__1, &c_b24) * part1 * part2;
    }
    return ret_val;
} /* func4_ */

/* 50--------INTEGRAL FOR ENERGY DENSITY OF NEUTRINO------------------------------ */
doublereal func5_(real *x)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double exp(doublereal);

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] rate1, nudens */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Contains integrands to be integrated. */
/* ----------COMMON AREAS. */
/* Model parameters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------MODEL PARAMETERS. */
/* Integration parame */
/* ----------NEUTRINO PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* Temperature (in units of MeV). */
/* Neutrino temperature (in units of Me */
/* Neutrino temperature (units of 10**9 */
/* ----------LOCAL VARIABLES. */
/* Normalizing constant. */
/* Value at which function is evaluated */
/* Exponential expression with photon t */
/* Exponential expression with neutrino */
/* Computing 3rd power */
    r__1 = *x;
    ret_val = r__1 * (r__1 * r__1) * .05066067740361991f / (exp(*x / 
	    nupar_1.tnu - modpr_1.xi[nupar_1.nu - 1]) + 1.f);
    return ret_val;
} /* func5_ */

/* 60--------INTEGRAL FOR ENERGY DENSITY OF ANTINEUTRINO-------------------------- */
doublereal func6_(real *x)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double exp(doublereal);

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] rate1, nudens */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Contains integrands to be integrated. */
/* ----------COMMON AREAS. */
/* Model parameters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------MODEL PARAMETERS. */
/* Integration parame */
/* ----------NEUTRINO PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* Temperature (in units of MeV). */
/* Neutrino temperature (in units of Me */
/* Neutrino temperature (units of 10**9 */
/* ----------LOCAL VARIABLES. */
/* Normalizing constant. */
/* Value at which function is evaluated */
/* Exponential expression with photon t */
/* Exponential expression with neutrino */
/* Computing 3rd power */
    r__1 = *x;
    ret_val = r__1 * (r__1 * r__1) * .05066067740361991f / (exp(*x / 
	    nupar_1.tnu + modpr_1.xi[nupar_1.nu - 1]) + 1.f);
    return ret_val;
} /* func6_ */

/* ----------REFERENCES----------------------------------------------------------- */
/*     Forms of the integrals involved can be found in */
/*       Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683. */
/*       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415. */
/* ========================IDENTIFICATION DIVISION================================ */
doublereal xintd_(real *xlow, real *xhi, E_fp func, integer *nq)
{
    /* Initialized data */

    static real u[6] = { -.93246951420315f,-.66120938646627f,-.2386191860832f,
	    .2386191860832f,.66120938646627f,.93246951420315f };
    static real w[6] = { .17132449237917f,.36076157304814f,.46791393457269f,
	    .46791393457269f,.36076157304814f,.17132449237917f };
    static integer np = 6;

    /* System generated locals */
    integer i__1, i__2;
    real ret_val;

    /* Local variables */
    static real f, x, sum, cent, dist;
    static integer nint, npnt;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] rate1, nudens */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Computes the integral of the function "func". */
/* ==========================DECLARATION DIVISION================================= */
/* ----------INPUT VARIABLES. */
/* Array of low limits. */
/* Array of high limits. */
/* ----------COMPUTATION VARIABLES. */
/* Number of six point gaussian quads. */
/* Size of quad interval. */
/* Center of quad interval. */
/* Variables of integration. */
/* ----------COUNTERS. */
/* Summation of terms. */
/* Interval number. */
/* Point number. */
/* ----------ABSCISSAS AND WEIGHT FACTORS. */
/* Total number of points in interval. */
/* Abscissas. */
/* ==============================DATA DIVISION==================================== */
/* ----------ABSCISSAS AND WEIGHT FACTORS. */
/* Weight factor. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------DO INTEGRATION------------------------------------------------------- */
/* 6 point Gaussian integration. */
    sum = 0.f;
    dist = (*xhi - *xlow) / (real) (*nq);
/* Size of quad interval. */
    i__1 = *nq;
    for (nint = 1; nint <= i__1; ++nint) {
	cent = *xlow + ((real) nint - .5f) * dist;
/* Center of interval. */
	i__2 = np;
	for (npnt = 1; npnt <= i__2; ++npnt) {
	    x = cent + dist * .5f * u[npnt - 1];
/* Integration point. */
	    f = (*func)(&x);
/* Evaluate function x(1). */
	    sum += f * w[npnt - 1];
/* Add up sum. */
	}
    }
/* 20--------GET INTEGRAL VALUE--------------------------------------------------- */
    ret_val = sum * dist * .5f;
/* Do integral. */
    return ret_val;
} /* xintd_ */

/* ========================IDENTIFICATION DIVISION================================ */
doublereal ex_(real *x)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double exp(doublereal);

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] start, rate2, rate3, rate4, sol */
/*               - [function] eval */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Exponential function with underflow precaution. */
/* ===========================PROCEDURE DIVISION================================== */
    if (*x > 88.029f) {
/* In danger of overflow. */
	ret_val = exp(88.029f);
    } else {
	if (*x < -88.722f) {
/* In danger of underflow. */
	    ret_val = 0.f;
	} else {
/* Value of x in allowable range. */
	    ret_val = exp(*x);
	}
    }
    return ret_val;
/* ----------NOTE----------------------------------------------------------------- */
/*     The overflow limit for the VAX/VMS system is exp(88.029). */
/*     The underflow limit for the VAX/VMS system is exp(-88.722). */
} /* ex_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int sol_(integer *loop)
{
    /* Initialized data */

    static real si[11] = { 1.f,1.f,1.f,1.f,1.f,2.f,3.f,2.f,1.f,1.f,2.f };
    static real sj[11] = { 0.f,1.f,1.f,0.f,1.f,0.f,0.f,1.f,1.f,1.f,0.f };
    static real sk[11] = { 0.f,0.f,1.f,0.f,0.f,1.f,0.f,0.f,1.f,0.f,2.f };
    static real sl[11] = { 1.f,1.f,1.f,2.f,2.f,1.f,1.f,1.f,2.f,3.f,1.f };

    /* Format strings */
    static char fmt_6000[] = "(\002 \002,\002** y(\002,i2,\002) fails to con"
	    "verge **\002)";
    static char fmt_6002[] = "(\002 \002,\002** \002,i2,\002 th diagonal ter"
	    "m equals zero **\002)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k, l, n, i1, j1;
    static real ci, cj, ck, cl, ri, rj, rk, rl;
    extern doublereal ex_(real *);
    static real yy[26], t932;
    static integer ind;
    static real t9m32, bdln;
    static integer isize1;
    extern /* Subroutine */ int eqslin_(integer *, integer *);
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___103 = { 0, 6, 0, fmt_6000, 0 };
    static cilist io___104 = { 0, 6, 0, fmt_6002, 0 };


/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] derivs */
/*     CALLS     - [subroutine] eqslin */
/*               - [function] ex */
/* ----------REMARKS. */
/*     Computes reverse strong and electromagnetic reaction rates. */
/*     Fills and solves matrix equation for dydt(i). */
/* ----------PARAMETERS. */
/* Input unit number. */
/* Output unit number. */
/* Number of nuclear reactions. */
/* -----------COMMON AREAS. */
/* Number of nuclides in calculation. */
/* Reaction parameter */
/* Reaction rates. */
/* Evolution paramete */
/* Evolution paramete */
/* Evolution paramete */
/* Time varying param */
/* Dynamic variables. */
/* Energy densities. */
/* Linear eqn coeffic */
/* Flags,counters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION PARAMETERS. */
/* Run option. */
/* Reaction code number (1-88). */
/* Incoming nuclide type (1-26). */
/* Incoming light nuclide type (1-6). */
/* Outgoing light nuclide type (1-6). */
/* Outgoing nuclide type (1-26). */
/* Reverse reaction coefficient. */
/* ----------REACTION RATES. */
/* Energy released in reaction. */
/* Forward reaction rate coefficients. */
/* ----------EVOLUTION PARAMETERS. */
/* Reverse reaction rate coefficients. */
/* Temperature (in units of 10**9 K). */
/* ----------EVOLUTION PARAMETERS (DERIVATIVES). */
/* Relative number abundances. */
/* ----------EVOLUTION PARAMETERS (ORIGINAL VALUES). */
/* Change in rel number abundances. */
/* ----------TIME VARIABLES. */
/* Rel # abund at start of iteration. */
/* ----------DYNAMIC VARIABLES. */
/* Time step. */
/* ----------ENERGY DENSITIES. */
/* Expansion rate. */
/* ----------COMPONENTS OF MATRIX EQUATION. */
/* Baryon mass density. */
/* Relates y(t-dt) to y(t). */
/* Contains y0 in inverse order. */
/* ----------COUNTERS AND FLAGS. */
/* yy in reverse order. */
/* Counts which Runge-Kutta loop. */
/* # time steps after outputting a line */
/* ----------RUN OPTIONS. */
/* Indicates if gaussian elimination fa */
/* Number of nuclides in computation. */
/* Equals isize + 1. */
/* ----------EVOLUTION EQUATION COEFFICIENTS. */
/* Number of reactions in computation. */
/* Equate to ii,jj,kk,ll. */
/* Equate to si,sj,sk,sl. */
/* ----------LOCAL VARIABLES. */
/* Coefficients of rate equation. */
/* Abundances at end of iteration. */
/* # of nuclide i,j,k,l */
/* (10**(-5))*volume expansion rate. */
/* Equate to iform. */
/* ==============================DATA DIVISION==================================== */
/* ----------NUMBER OF NUCLIDES IN REACTION TYPES 1-11. */
/* Element which does not converge. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------TEMPERATURE FACTORS AND INITIAL VALUES------------------------------- */
/* ..........TEMPERATURE FACTORS. */
    d__1 = (doublereal) evolp1_1.t9;
    t932 = pow_dd(&d__1, &c_b5);
/* t9**(3/2). */
    t9m32 = 1.f / t932;
/* ..........MATRIX SIZE. */
/* t9**(-3/2). */
    isize1 = runopt_1.isize + 1;
/* ..........INITIALIZE A-MATRIX. */
    i__1 = runopt_1.isize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = runopt_1.isize;
	for (j = 1; j <= i__2; ++j) {
	    lncoef_1.a[j + i__ * 26 - 27] = 0.;
/* Set a-matrix to zero. */
	}
    }
/* 20--------COMPUTE FACTORS FOR THE A-MATRIX------------------------------------- */
    i__1 = runopt_1.jsize;
    for (n = 1; n <= i__1; ++n) {
/* ..........EQUATE VARIABLES TO ARRAYS. */
	ind = recpr_1.iform[n - 1];
/* Type of reaction. */
	i__ = recpr_1.ii[n - 1];
/* ID # of incoming nuclide i. */
	j = recpr_1.jj[n - 1];
/* ID # of incoming nuclide j. */
	k = recpr_1.kk[n - 1];
/* ID # of outgoing nuclide k. */
	l = recpr_1.ll[n - 1];
/* ID # of outgoing nuclide l. */
	if (ind != 0 && i__ <= runopt_1.isize && l <= runopt_1.isize) {
/* Reactio */
	    ri = si[ind - 1];
/* # of incoming nuclide i. */
	    rj = sj[ind - 1];
/* # of incoming nuclide j. */
	    rk = sk[ind - 1];
/* # of outgoing nuclide k. */
	    rl = sl[ind - 1];
/* ..........COMPUTE DIFFERENT REACTION RATES. */
/* # of outgoing nuclide l. */
	    switch (ind) {
		case 1:  goto L201;
		case 2:  goto L202;
		case 3:  goto L203;
		case 4:  goto L204;
		case 5:  goto L205;
		case 6:  goto L206;
		case 7:  goto L207;
		case 8:  goto L208;
		case 9:  goto L209;
		case 10:  goto L210;
		case 11:  goto L211;
	    }
L201:
/* 1-0-0-1 configuration. */
	    ci = rates_1.f[n - 1];
/* (Ref 1). */
	    cj = 0.f;
	    ck = 0.f;
	    cl = rates_1.r__[n - 1];
	    goto L212;
L202:
/* 1-1-0-1 configuration. */
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * 1e10f * t932 * ex_(&
		    r__1) * rates_1.f[n - 1];
/* (Ref 2). */
	    rates_1.f[n - 1] = endens_1.rhob * rates_1.f[n - 1];
	    ci = evolp1_1.y[j - 1] * rates_1.f[n - 1] / 2.f;
	    cj = evolp1_1.y[i__ - 1] * rates_1.f[n - 1] / 2.f;
	    ck = 0.f;
	    cl = rates_1.r__[n - 1];
	    goto L212;
L203:
/* 1-1-1-1 configuration. */
	    rates_1.f[n - 1] = endens_1.rhob * rates_1.f[n - 1];
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * ex_(&r__1) * rates_1.f[
		    n - 1];
/* (Ref 3). */
	    ci = evolp1_1.y[j - 1] * rates_1.f[n - 1] / 2.f;
	    cj = evolp1_1.y[i__ - 1] * rates_1.f[n - 1] / 2.f;
	    ck = evolp1_1.y[l - 1] * rates_1.r__[n - 1] / 2.f;
	    cl = evolp1_1.y[k - 1] * rates_1.r__[n - 1] / 2.f;
	    goto L212;
L204:
/* 1-0-0-2 configuration. */
	    ci = rates_1.f[n - 1];
	    cj = 0.f;
	    ck = 0.f;
	    cl = evolp1_1.y[l - 1] * rates_1.r__[n - 1] / 2.f;
	    goto L212;
L205:
/* 1-1-0-2 configuration. */
	    rates_1.f[n - 1] = endens_1.rhob * rates_1.f[n - 1];
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * ex_(&r__1) * rates_1.f[
		    n - 1];
/* (Ref 3). */
	    ci = evolp1_1.y[j - 1] * rates_1.f[n - 1] / 2.f;
	    cj = evolp1_1.y[i__ - 1] * rates_1.f[n - 1] / 2.f;
	    ck = 0.f;
	    cl = evolp1_1.y[l - 1] * rates_1.r__[n - 1] / 2.f;
	    goto L212;
L206:
/* 2-0-1-1 configuration. */
	    rates_1.f[n - 1] = endens_1.rhob * rates_1.f[n - 1];
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * ex_(&r__1) * rates_1.f[
		    n - 1];
/* (Ref 3). */
	    ci = evolp1_1.y[i__ - 1] * rates_1.f[n - 1] / 2.f;
	    cj = 0.f;
	    ck = evolp1_1.y[l - 1] * rates_1.r__[n - 1] / 2.f;
	    cl = evolp1_1.y[k - 1] * rates_1.r__[n - 1] / 2.f;
	    goto L212;
L207:
/* 3-0-0-1 configuration. */
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * 1e20f * t932 * t932 * 
		    ex_(&r__1) * rates_1.f[n - 1];
/* (Ref 4). */
	    rates_1.f[n - 1] = endens_1.rhob * endens_1.rhob * rates_1.f[n - 
		    1];
	    ci = evolp1_1.y[i__ - 1] * evolp1_1.y[i__ - 1] * rates_1.f[n - 1] 
		    / 6.f;
	    cj = 0.f;
	    ck = 0.f;
	    cl = rates_1.r__[n - 1];
	    goto L212;
L208:
/* 2-1-0-1 configuration. */
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * 1e20f * t932 * t932 * 
		    ex_(&r__1) * rates_1.f[n - 1];
/* (Ref 4). */
	    rates_1.f[n - 1] = endens_1.rhob * endens_1.rhob * rates_1.f[n - 
		    1];
	    ci = evolp1_1.y[j - 1] * evolp1_1.y[i__ - 1] * rates_1.f[n - 1] / 
		    3.f;
	    cj = evolp1_1.y[i__ - 1] * evolp1_1.y[i__ - 1] * rates_1.f[n - 1] 
		    / 6.f;
	    ck = 0.f;
	    cl = rates_1.r__[n - 1];
	    goto L212;
L209:
/* 1-1-1-2 configuration. */
	    rates_1.f[n - 1] = endens_1.rhob * rates_1.f[n - 1];
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * 1e-10f * t9m32 * 
		    endens_1.rhob * ex_(&r__1) * rates_1.f[n - 1];
/* (Ref 5) */
	    ci = evolp1_1.y[j - 1] * rates_1.f[n - 1] / 2.f;
	    cj = evolp1_1.y[i__ - 1] * rates_1.f[n - 1] / 2.f;
	    ck = evolp1_1.y[l - 1] * evolp1_1.y[l - 1] * rates_1.r__[n - 1] / 
		    6.f;
	    cl = evolp1_1.y[k - 1] * evolp1_1.y[l - 1] * rates_1.r__[n - 1] / 
		    3.f;
	    goto L212;
L210:
/* 1-1-0-3 configuration. */
	    rates_1.f[n - 1] = endens_1.rhob * rates_1.f[n - 1];
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * 1e-10f * t9m32 * 
		    endens_1.rhob * ex_(&r__1) * rates_1.f[n - 1];
/* (Ref 5) */
	    ci = evolp1_1.y[j - 1] * rates_1.f[n - 1] / 2.f;
	    cj = evolp1_1.y[i__ - 1] * rates_1.f[n - 1] / 2.f;
	    ck = 0.f;
	    cl = evolp1_1.y[l - 1] * evolp1_1.y[l - 1] * rates_1.r__[n - 1] / 
		    6.f;
	    goto L212;
L211:
/* 2-0-2-1 configuration. */
	    rates_1.f[n - 1] = endens_1.rhob * rates_1.f[n - 1];
	    r__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_1.r__[n - 1] = recpr_1.rev[n - 1] * 1e-10f * t9m32 * 
		    endens_1.rhob * ex_(&r__1) * rates_1.f[n - 1];
/* (Ref 5) */
	    ci = evolp1_1.y[i__ - 1] * rates_1.f[n - 1] / 2.f;
	    cj = 0.f;
	    ck = evolp1_1.y[l - 1] * evolp1_1.y[k - 1] * rates_1.r__[n - 1] / 
		    3.f;
	    cl = evolp1_1.y[k - 1] * evolp1_1.y[k - 1] * rates_1.r__[n - 1] / 
		    6.f;
L212:
/* 30--------CONSTRUCT THE A-MATRIX----------------------------------------------- */
	    i__ = isize1 - i__;
/* Invert i index. */
	    j = isize1 - j;
/* Invert j index. */
	    k = isize1 - k;
/* Invert k index. */
	    l = isize1 - l;
/* ..........FILL I NUCLIDE COLUMN. */
/* Invert l index. */
	    if (j <= runopt_1.isize) {
		lncoef_1.a[j + i__ * 26 - 27] += rj * ci;
	    }
	    if (k <= runopt_1.isize) {
		lncoef_1.a[k + i__ * 26 - 27] -= rk * ci;
	    }
	    lncoef_1.a[i__ + i__ * 26 - 27] += ri * ci;
	    lncoef_1.a[l + i__ * 26 - 27] -= rl * ci;
/* ..........FILL J NUCLIDE COLUMN. */
	    if (j <= runopt_1.isize) {
		lncoef_1.a[j + j * 26 - 27] += rj * cj;
		if (k <= runopt_1.isize) {
		    lncoef_1.a[k + j * 26 - 27] -= rk * cj;
		}
		lncoef_1.a[i__ + j * 26 - 27] += ri * cj;
		lncoef_1.a[l + j * 26 - 27] -= rl * cj;
	    }
/* ..........FILL K NUCLIDE COLUMN. */
	    if (k <= runopt_1.isize) {
		if (j <= runopt_1.isize) {
		    lncoef_1.a[j + k * 26 - 27] -= rj * ck;
		}
		lncoef_1.a[k + k * 26 - 27] += rk * ck;
		lncoef_1.a[i__ + k * 26 - 27] -= ri * ck;
		lncoef_1.a[l + k * 26 - 27] += rl * ck;
	    }
/* ..........FILL L NUCLIDE COLUMN. */
	    if (j <= runopt_1.isize) {
		lncoef_1.a[j + l * 26 - 27] -= rj * cl;
	    }
	    if (k <= runopt_1.isize) {
		lncoef_1.a[k + l * 26 - 27] += rk * cl;
	    }
	    lncoef_1.a[i__ + l * 26 - 27] -= ri * cl;
	    lncoef_1.a[l + l * 26 - 27] += rl * cl;
	}
/* ((ind.ne.0).and.(i.le.isize).and.(l.le.isize)) */
    }
/* 40--------PUT A-MATRIX AND B-VECTOR IN FINAL FORM OF MATRIX EQUATION----------- */
/* n = 1,jsize */
    bdln = thermcb_1.hubcst * 3.f * 1e-5f;
/* (10**(-5))*(Expansion rate). */
    i__1 = runopt_1.isize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = isize1 - i__;
/* Invert the rows. */
	i__2 = runopt_1.isize;
	for (j = 1; j <= i__2; ++j) {
	    j1 = isize1 - j;
/* Invert the columns. */
	    if ((d__1 = lncoef_1.a[j + i__ * 26 - 27], abs(d__1)) < bdln * 
		    evolp3_1.y0[j1 - 1] / evolp3_1.y0[i1 - 1]) {
		lncoef_1.a[j + i__ * 26 - 27] = 0.;
/* Set 0 if tiny. */
	    } else {
		lncoef_1.a[j + i__ * 26 - 27] *= ttime_1.dt;
/* Bring dt over to other side. */
	    }
	}
	lncoef_1.a[i__ + i__ * 26 - 27] += 1.;
/* Add identity matrix to a-matrix. */
	lncoef_1.bkevin[i1 - 1] = evolp3_1.y0[i__ - 1];
/* Initial abundances. */
    }
/* 50--------SOLVE EQUATIONS TO GET DERIVATIVE------------------------------------ */
/* ..........SET MONITOR FLAG AND SOLVE BY GAUSSIAN ELIMINATION. */
    if (*loop == 1) {
	eqslin_(&flags_1.ip, &ierror);
    } else {
	eqslin_(&c__0, &ierror);
    }
/* ..........OBTAIN DERIVATIVE. */
    i__1 = runopt_1.isize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yy[i__ - 1] = lncoef_1.yx[isize1 - i__ - 1];
/* Abundance at t+dt. */
	evolp2_1.dydt[i__ - 1] = (yy[i__ - 1] - evolp3_1.y0[i__ - 1]) / 
		ttime_1.dt;
/* Take derivative. */
    }
/* 60--------POSSIBLE ERROR MESSAGES AND EXIT------------------------------------- */
    if (flags_1.mbad != 0) {
/* Problem in gaussian elimination. */
	if (flags_1.mbad == -1) {
	    s_wsfe(&io___103);
	    do_fio(&c__1, (char *)&ierror, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
/* Error message. */
	if (flags_1.mbad >= 1) {
	    s_wsfe(&io___104);
	    do_fio(&c__1, (char *)&flags_1.mbad, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
/* Error message. */
    }
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     1) The coefficients are given in general as: */
/*             ci = ri*(y(j)**rj)*(y(i)**(ri-1)*f(n)/ */
/*                  ((ri+rj)*fac(ri)*fac(rj)) */
/*             cj = rj*(y(i)**ri)*(y(j)**(rj-1)*f(n)/ */
/*                  ((ri+rj)*fac(ri)*fac(rj)) */
/*             ck = rk*(y(l)**rl)*(y(k)**(rk-1)*f(n)/ */
/*                  ((rk+rl)*fac(rk)*fac(rl)) */
/*             cl = rl*(y(k)**rk)*(y(l)**(rl-1)*f(n)/ */
/*                  ((rk+rl)*fac(rk)*fac(rl)) */
/*        in which fac(x) is the factorial of x. */
/*     2) Form of reverse rate given in */
/*        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247, */
/*          tables 1B, 4B, 7B. */
/*     3) Form of reverse rate given in */
/*        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247, */
/*          tables 2B, 3B, 5B, 6B, 8B, 9B, 10B. */
/*     4) Form of reverse rate given in */
/*        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247, */
/*          table 11B. */
/*     5) Form of reverse rate given in */
/*        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247, */
/*          tables 12B, 13B. */
} /* sol_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int eqslin_(integer *icnvm, integer *ierror)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static real r__;
    static doublereal x[26], a0[676]	/* was [26][26] */, cx, sum;
    static real xdy;
    static integer nord;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] sol */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Solves for new abundances using gaussian elimination */
/*     with back substitution, no pivoting. */
/* ----------PARAMETERS. */
/* Rank of matrix. */
/* Higher order in correction. */
/* ----------COMMON AREAS. */
/* Tolerance for convergence (.ge. 1.e- */
/* Computation parame */
/* Lin eqn coefficien */
/* Flags, counters. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------COMPUTATION PARAMETER. */
/* Run options. */
/* ----------MATRIX COEFFICIENTS FOR LINEAR EQUATION. */
/* Accumulation increment. */
/* Coefficient array. */
/* Right-hand vector w/o manipulation. */
/* ----------COUNTERS AND FLAGS. */
/* Solution vector. */
/* ----------RUN OPTION. */
/* Indicates type of error. */
/* ----------LOCAL MATRICES AND VECTORS. */
/* Number of nuclides in computation. */
/* Coefficient array w/o manipulation. */
/* ----------LOCAL COMPUTATION VARIABLES. */
/* Right-hand vector. */
/* Scaling factor in triangularization. */
/* Sum for backsubstitution. */
/* ----------LOCAL COUNTERS. */
/* Relative error. */
/* Order of correction. */
/* Convergence monitor. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------INITIALIZE VECTOR---------------------------------------------------- */
/* ..........SET COUNTERS TO ZERO. */
/* ith nuclide fails to converge. */
    nord = 0;
/* No corrections yet. */
    flags_1.mbad = 0;
/* ..........SET RIGHT-HAND AND SOLUTION VECTORS TO INITIAL VALUES. */
/* No errors yet. */
    i__1 = runopt_1.isize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__ - 1] = lncoef_2.bkevin[i__ - 1];
/* Right-hand vector. */
	lncoef_2.y[i__ - 1] = 0.f;
/* Solution vector. */
    }
/* ..........SAVE MATRIX. */
    if (*icnvm == compr_1.inc) {
/* Monitor convergence. */
	i__1 = runopt_1.isize;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = runopt_1.isize;
	    for (j = 1; j <= i__2; ++j) {
		a0[j + i__ * 26 - 27] = lncoef_2.a[j + i__ * 26 - 27];
/* Initial value of coefficient array. */
	    }
	}
    }
/* 20--------TRIANGULARIZE MATRIX AND SAVE OPERATOR------------------------------- */
/* ..........CHECK TO SEE THAT THERE ARE NO ZEROES AT PIVOT POINTS. */
    i__1 = runopt_1.isize - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (lncoef_2.a[i__ + i__ * 26 - 27] == 0.) {
/* Don't want to divide by zero. */
	    flags_1.mbad = i__;
/* Position of zero coefficient. */
	    return 0;
/* Terminate matrix evaluation. */
	}
/* ..........TRIANGULARIZE MATRIX. */
	i__2 = runopt_1.isize;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (lncoef_2.a[j + i__ * 26 - 27] != 0.) {
/* Progress diagonally down the column. */
		cx = lncoef_2.a[j + i__ * 26 - 27] / lncoef_2.a[i__ + i__ * 
			26 - 27];
/* Scaling factor down the column. */
		i__3 = runopt_1.isize;
		for (k = i__ + 1; k <= i__3; ++k) {
/* Progress diagonally along row. */
		    lncoef_2.a[j + k * 26 - 27] -= cx * lncoef_2.a[i__ + k * 
			    26 - 27];
/* Subtract scaled coeff along */
		}
		lncoef_2.a[j + i__ * 26 - 27] = cx;
/* ..........OPERATE ON RIGHT-HAND VECTOR. */
/* Scaled coefficient. */
		x[j - 1] -= cx * x[i__ - 1];
/* Subtract off scaled coefficient. */
	    }
	}
    }
/* 30--------DO BACK SUBSTITUTION------------------------------------------------- */
L300:
    x[runopt_1.isize - 1] /= lncoef_2.a[runopt_1.isize + runopt_1.isize * 26 
	    - 27];
/* Solution for ultimate positi */
    lncoef_2.y[runopt_1.isize - 1] += x[runopt_1.isize - 1];
    for (i__ = runopt_1.isize - 1; i__ >= 1; --i__) {
/* From i = penultimate to i = 1. */
	sum = 0.;
	i__1 = runopt_1.isize;
	for (j = i__ + 1; j <= i__1; ++j) {
	    sum += lncoef_2.a[i__ + j * 26 - 27] * x[j - 1];
/* Sum up all previous terms. */
	}
	x[i__ - 1] = (x[i__ - 1] - sum) / lncoef_2.a[i__ + i__ * 26 - 27];
	lncoef_2.y[i__ - 1] += x[i__ - 1];
/* Add difference to initial value. */
    }
/* 40--------TESTS AND EXITS------------------------------------------------------ */
    if (*icnvm == compr_1.inc) {
	i__1 = runopt_1.isize;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (lncoef_2.y[i__ - 1] != 0.f) {
		xdy = (d__1 = x[i__ - 1] / lncoef_2.y[i__ - 1], abs(d__1));
/* Relative error. */
		if (xdy > 2e-4f) {
		    if (nord < 1) {
/* Continue to higher orders. */
			++nord;
/* ..........FIND ERROR IN RIGHT-HAND VECTOR. */
			i__2 = runopt_1.isize;
			for (j = 1; j <= i__2; ++j) {
			    r__ = 0.f;
/* Initialize r. */
			    i__3 = runopt_1.isize;
			    for (k = 1; k <= i__3; ++k) {
				r__ += a0[j + k * 26 - 27] * lncoef_2.y[k - 1]
					;
/* Left side with approximate sol */
			    }
			    x[j - 1] = lncoef_2.bkevin[j - 1] - r__;
/* Subtract difference from right */
			}
/* ..........OPERATE ON RIGHT-HAND VECTOR. */
			i__2 = runopt_1.isize - 1;
			for (j = 1; j <= i__2; ++j) {
			    i__3 = runopt_1.isize;
			    for (k = j + 1; k <= i__3; ++k) {
				x[k - 1] -= lncoef_2.a[k + j * 26 - 27] * x[j 
					- 1];
/* Subtract off scaled coef */
			    }
			}
			goto L300;
/* Go for another iteratiion. */
		    } else {
/* ..........NOT ENOUGH CONVERGENCE. */
			flags_1.mbad = -1;
/* Signal error problem. */
			*ierror = i__;
/* ith nuclide for which x/y checked. */
			return 0;
		    }
/* (nord.lt.mord) */
		}
/* (xdy.gt.eps) */
	    }
/* (y(i).ne.0) */
	}
/* i = 1,isize */
    }
/* (icnvm.eq.inc) */
    return 0;
/* No more iterations & relative error */
} /* eqslin_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     Changes (to run inder DEC unix f77): */
/*     ----------------------------------- */
/*     COMMON /therm/ -> COMMON /thermcb/ */
/*     COMMON /rates/ f,r(nrec) -> COMMON /rates/ f,r */

/*     Default neutron lifetime 888.54 -> 885.7 (+- 0.8) - PDG 2002/3 */
/* Subroutine */ int rate0_(void)
{
/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] start */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Generates weak decay rates. */
/* ----------PARAMETER. */
/* ----------COMMON AREA. */
/* Number of nuclear reactions. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION RATE COEFFICIENTS. */
/* Reaction rates. */
/* Forward reaction rate coefficients. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------SET DECAY RATE COEFFICIENTS------------------------------------------ */
/* .......H3 -> e- + v + He3.........(Tilly-Weller-Hasan 1987) */
    rates_1.f[1] = 1.79e-9f;
/* .......Li8 -> e- + v + 2He4.......(Ajzenberg-Selove 1988) */
    rates_1.f[2] = .827f;
/* .......B12 -> e- + B + C12........(Ajzenberg-Selove 1990) */
    rates_1.f[3] = 34.3f;
/* .......C14 -> e- + v + N14........(Ajzenberg-Selove 1986) */
    rates_1.f[4] = 3.834e-12f;
/* .......B8 -> e+ + v + 2He4........(Ajzenberg-Selove 1988) */
    rates_1.f[5] = .9f;
/* .......C11 -> e+ + v + B11........(Ajzenberg-Selove 1990) */
    rates_1.f[6] = 5.668e-4f;
/* .......N12 -> e+ + v + C12........(Ajzenberg-Selove 1990) */
    rates_1.f[7] = 63.01f;
/* .......N13 -> e+ + v + C13........(Ajzenberg-Selove 1986) */
    rates_1.f[8] = .001159f;
/* .......O14 -> e+ + v + N14........(Ajzenberg-Selove 1986) */
    rates_1.f[9] = .0098171f;
/* .......O15 -> e+ + v + N15........(Ajzenberg-Selove 1986) */
    rates_1.f[10] = .0056704f;
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     Ajzenberg-Selove, F., 1990, Nucl. Phys. A506, 1. */
/*     Ajzenberg-Selove, F., 1988, Nucl. Phys. A490, 1. */
/*     Ajzenberg-Selove, F., 1986, Nucl. Phys. A449, 1. */
/*     Tilley, D.R., Weller, H.R., and Hasan, H.H., 1987, Nucl. Phys. A474, 1. */
} /* rate0_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int rate1_(real *tph)
{
    static real w[2], x[2], y[2], z__[2];
    extern doublereal func1_(real *), func2_(real *), func3_(real *), func4_(
	    real *);
    static real part1, part2, part3, part4;
    extern doublereal xintd_(real *, real *, E_fp, integer *);
    static real uplim1, uplim2, uplim3, uplim4;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] start, derivs */
/*     CALLS     - [function] xintd, eval */
/* ----------REMARKS. */
/*     Generates rate coefficients for weak n->p and p->n reactions. */
/* ----------PARAMETERS. */
/* Number of nuclear reactions. */
/* ----------COMMON AREAS. */
/* Number of gaussian quads. */
/* Reaction rates. */
/* Model parameters. */
/* Dynamic variable */
/* ----------EXTERNAL FUNCTIONS. */
/* Integration parame */
/* Part 1 of n->p rate. */
/* Part 2 of n->p rate. */
/* Part 1 of p->n rate. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION RATES. */
/* Part 2 of p->n rate. */
/* Forward reaction rate coefficients. */
/* ----------EARLY UNIVERSE MODEL PARAMETERS. */
/* Reverse reaction rate coefficients. */
/* Neutron lifetime. */
/* ----------DYNAMIC VARIABLES. */
/* Neutrino degeneracy parameters. */
/* xi(1) is e neutrino degeneracy param */
/* xi(2) is m neutrino degeneracy param */
/* xi(3) is t neutrino degeneracy param */
/* ----------NEUTRINO PARAMETERS. */
/* Thermodynamic variables (energy dens */
/* Temperature (in units of MeV). */
/* ----------LOCAL VARIABLES. */
/* Neutrino temperature (in units of Me */
/* Photon temperature. */
/* Upper limits for exponentials, forwa */
/* Upper limits for exponentials, rever */
/* Upper limits for integrals for forwa */
/* Upper limits for integrals for rever */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------COMPUTE WEAK REACTION RATES (NONDEGENERATE)-------------------------- */
/*      IF (xi(1).eq.0.) THEN */
/*        f(1)  = thm(13)/tau        !Forward rate for weak np reaction. */
/*        r(1)  = thm(14)/tau        !Reverse rate for weak np reaction. */
/*      ELSE */
/* 20--------COMPUTE WEAK REACTION RATES (DEGENERATE)----------------------------- */
/* Parts of integrals for forward rate. */
/* Parts of integrals for reverse rate. */
    nupar_1.t9mev = *tph * .086171f;
/* Convert photon temp to units of MeV. */
    nupar_1.tnmev = nupar_1.tnu * .086171f;
/* ..........COMPUTE OVERFLOW LIMITS FOR LIMITS OF INTEGRATION (Ref 1 & 2). */
/* Convert neutrino temp to units of Me */
    w[0] = -(nupar_1.t9mev / .511f) * -88.722f;
    w[1] = nupar_1.tnmev / .511f * (modpr_1.xi[0] + 88.029f) + 2.531f;
    x[0] = nupar_1.t9mev / .511f * 88.029f;
    x[1] = -(nupar_1.tnmev / .511f) * (modpr_1.xi[0] - 88.722f) - 2.531f;
    y[0] = -(nupar_1.t9mev / .511f) * -88.722f;
    y[1] = nupar_1.tnmev / .511f * (88.029f - modpr_1.xi[0]) - 2.531f;
    z__[0] = nupar_1.t9mev / .511f * 88.029f;
    z__[1] = -(nupar_1.tnmev / .511f) * (-88.722f - modpr_1.xi[0]) + 2.531f;
/* ..........COMPARE LIMITS AND TAKE LARGER OF THE TWO. */
    uplim1 = dabs(w[0]);
    uplim2 = dabs(x[0]);
    uplim3 = dabs(y[0]);
    uplim4 = dabs(z__[0]);
    if (uplim1 < dabs(w[1])) {
	uplim1 = w[1];
    }
    if (uplim2 < dabs(x[1])) {
	uplim2 = x[1];
    }
    if (uplim3 < dabs(y[1])) {
	uplim3 = y[1];
    }
    if (uplim4 < dabs(z__[1])) {
	uplim4 = z__[1];
    }
/* ..........EVALUATE THE INTEGRALS NUMERICALLY. */
    part1 = xintd_(&c_b60, &uplim1, (E_fp)func1_, &c__50);
    part2 = xintd_(&c_b60, &uplim2, (E_fp)func2_, &c__50);
    part3 = xintd_(&c_b60, &uplim3, (E_fp)func3_, &c__50);
    part4 = xintd_(&c_b60, &uplim4, (E_fp)func4_, &c__50);
    rates_1.f[0] = part1 + part2;
/* Add 2 integrals to get forward rate. */
    rates_1.r__[0] = part3 + part4;
/*      END IF !(xi(1).eq.0.) */
/* Add 2 integrals to get reverse rate. */
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     1) Forms of the integrals involved can be found in */
/*          Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683. */
/*          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415. */

/*     2) The overflow limit for the VAX/VMS system is exp(88.029). */
/*        The underflow limit for the VAX/VMS system is exp(-88.722). */
} /* rate1_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int rate2_(void)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    extern doublereal ex_(real *);
    static real t9a, t9b, t912, t913, t923, t932, t943, t953, t9c, t9d, t9e, 
	    t9f, t9m1, t9a32, t9b32, t9c13, t9d13, t9e13, t9f13, t9c56, t9d56,
	     t9e56, t9f56, t9m23, t9m32;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] derivs */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Generates rate coefficients for reactions involving nuclides */
/*     up to A = 9. */
/* ----------PARAMETER. */
/* Number of nuclear reactions. */
/* ----------COMMON AREAS. */
/* Number of nuclides in calculation. */
/* Reaction rates. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION RATES. */
/* Evolution parameters. */
/* Forward reaction rate coefficients. */
/* ----------EVOLUTION PARAMETER. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------TEMPERATURE FACTORS-------------------------------------------------- */
/* Temperature of photons (units of 10* */
    d__1 = (doublereal) evolp1_1.t9;
    t913 = pow_dd(&d__1, &c_b69);
/* t9**(1/3) */
    t923 = t913 * t913;
/* t9**(2/3) */
    t943 = t923 * t923;
/* t9**(4/3) */
    t953 = evolp1_1.t9 * t923;
/* t9**(5/3) */
    t912 = sqrt(evolp1_1.t9);
/* t9**(1/2) */
    t932 = evolp1_1.t9 * t912;
/* t9**(3/2) */
    t9m1 = 1 / evolp1_1.t9;
/* t9**(-1) */
    t9m23 = 1.f / t923;
/* t9**(-2/3) */
    t9m32 = 1.f / t932;
/* t9**(-3/2) */
    t9a = evolp1_1.t9 / (evolp1_1.t9 * 13.076f + 1.f);
/* For reaction 17. */
    d__1 = (doublereal) t9a;
    t9a32 = pow_dd(&d__1, &c_b5);
/* t9a**(3/2) */
    t9b = evolp1_1.t9 / (evolp1_1.t9 * 49.18f + 1.f);
/* For reaction 18. */
    d__1 = (doublereal) t9b;
    t9b32 = pow_dd(&d__1, &c_b5);
/* t9b**(3/2) */
    if (evolp1_1.t9 > 10.f) {
/* For reaction 22. */
	t9c = 1.f;
    } else {
	d__1 = (doublereal) (1.f - evolp1_1.t9 * .0969f);
	t9c = evolp1_1.t9 / (1.f - evolp1_1.t9 * .0969f + t953 * .0284f / 
		pow_dd(&d__1, &c_b72));
    }
    d__1 = (doublereal) t9c;
    t9c13 = pow_dd(&d__1, &c_b73);
/* t9c**(1/3) */
    d__1 = (doublereal) t9c;
    t9c56 = pow_dd(&d__1, &c_b74);
/* t9c**(5/6) */
    t9d = evolp1_1.t9 / (evolp1_1.t9 * .759f + 1.f);
/* For reaction 24. */
    d__1 = (doublereal) t9d;
    t9d13 = pow_dd(&d__1, &c_b73);
/* t9d**(1/3) */
    d__1 = (doublereal) t9d;
    t9d56 = pow_dd(&d__1, &c_b74);
/* t9d**(5/6) */
    t9e = evolp1_1.t9 / (evolp1_1.t9 * .1378f + 1.f);
/* For reaction 26. */
    d__1 = (doublereal) t9e;
    t9e13 = pow_dd(&d__1, &c_b73);
/* t9e**(1/3) */
    d__1 = (doublereal) t9e;
    t9e56 = pow_dd(&d__1, &c_b74);
/* t9e**(5/6) */
    t9f = evolp1_1.t9 / (evolp1_1.t9 * .1071f + 1.f);
/* For reaction 27. */
    d__1 = (doublereal) t9f;
    t9f13 = pow_dd(&d__1, &c_b73);
/* t9f**(1/3) */
    d__1 = (doublereal) t9f;
    t9f56 = pow_dd(&d__1, &c_b74);
/* 20--------NEUTRON, PHOTON REACTIONS-------------------------------------------- */
/* .......H(n,g)H2...................(Smith-Kawano-Malaney 1992) */
/* t9f**(5/6) */
    rates_1.f[11] = (1.f - t912 * .8504f + evolp1_1.t9 * .4895f - t932 * 
	    .09623f + evolp1_1.t9 * .008471f * evolp1_1.t9 - evolp1_1.t9 * 
	    2.8e-4f * t932) * 47420.f;
/* .......H2(n,g)H3..................(Wagoner 1969) */
    rates_1.f[12] = (evolp1_1.t9 * 18.9f + 1.f) * 66.2f;
/* .......He3(n,g)He4................(Wagoner 1969) */
    rates_1.f[13] = (evolp1_1.t9 * 905.f + 1.f) * 6.62f;
/* .......Li6(n,g)Li7................(Malaney-Fowler 1989) */
    rates_1.f[14] = 5100.f;
/* 30--------NEUTRON, PROTON REACTIONS-------------------------------------------- */
/* .......He3(n,p)H3.................(Smith-Kawano-Malaney 1992) */
    rates_1.f[15] = (1.f - t912 * .508f + evolp1_1.t9 * .228f) * 7.21e8f;
/* .......Be7(n,p)Li7................(Smith-Kawano-Malaney 1992) */
    r__1 = -.07486f / evolp1_1.t9;
    rates_1.f[16] = (1.f - t912 * .56f + evolp1_1.t9 * .179f - t932 * .0283f 
	    + evolp1_1.t9 * .002214f * evolp1_1.t9 - evolp1_1.t9 * 6.851e-5f *
	     t932) * 2.675e9f + t9a32 * 9.391e8f * t9m32 + t9m32 * 4.467e7f * 
	    ex_(&r__1);
/* 40--------NEUTRON, ALPHA REACTIONS--------------------------------------------- */
/* .......Li6(n,a)H3.................(Caughlan-Fowler 1988) */
    r__1 = -2.39f / evolp1_1.t9;
    rates_1.f[17] = t9m32 * 2.54e9f * ex_(&r__1) + (1.f - t9b32 * .261f / 
	    t932) * 1.68e8f;
/* .......Be7(n,a)He4................(Wagoner 1969) */
    rates_1.f[18] = (evolp1_1.t9 * 3760.f + 1.f) * 20500.f;
/* 50--------PROTON, PHOTON REACTIONS--------------------------------------------- */
/* .......H2(p,g)He3.................(Smith-Kawano-Malaney 1992) */
    r__1 = -3.72f / t913;
    rates_1.f[19] = t9m23 * 2650.f * ex_(&r__1) * (t913 * .112f + 1.f + t923 *
	     1.99f + evolp1_1.t9 * 1.56f + t943 * .162f + t953 * .324f);
/* .......H3(p,g)He4.................(Caughlan-Fowler 1988) */
    r__1 = -3.869f / t913;
    rates_1.f[20] = t9m23 * 2.2e4f * ex_(&r__1) * (t913 * .108f + 1.f + t923 *
	     1.68f + evolp1_1.t9 * 1.26f + t943 * .551f + t953 * 1.06f);
/* .......Li6(p,g)Be7................(Caughlan-Fowler 1988) */
    r__1 = -8.413f / t9c13;
    rates_1.f[21] = t9c56 * 6.69e5f * t9m32 * ex_(&r__1);
/* 60--------PROTON, ALPHA REACTIONS---------------------------------------------- */
/* .......Li6(p,a)He3................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 5.5f;
    r__1 = -8.413f / t913 - r__2 * r__2;
    r__3 = -17.763f / evolp1_1.t9;
    r__4 = -21.82f / evolp1_1.t9;
    rates_1.f[22] = t9m23 * 3.73e10f * ex_(&r__1) * (t913 * .05f + 1.f - t923 
	    * .061f - evolp1_1.t9 * .021f + t943 * .006f + t953 * .005f) + 
	    t9m32 * 1.33e10f * ex_(&r__3) + t9m1 * 1.29e9f * ex_(&r__4);
/* .......Li7(p,a)He4................(Smith-Kawano-Malaney 1992) */
    r__1 = -8.472f / t913;
    r__2 = -8.472f / t9d13;
    r__3 = -30.442f / evolp1_1.t9;
/* Computing 2nd power */
    r__5 = evolp1_1.t9 / 1.696f;
    r__4 = -8.472f / t913 - r__5 * r__5;
    r__6 = -4.478f / evolp1_1.t9;
    rates_1.f[23] = t9m23 * 1.096e9f * ex_(&r__1) - t9d56 * 4.83e8f * t9m32 * 
	    ex_(&r__2) + t9m32 * 1.06e10f * ex_(&r__3) + t9m23 * 1.56e5f * 
	    ex_(&r__4) * (t913 * .049f + 1.f - t923 * 2.498f + evolp1_1.t9 * 
	    .86f + t943 * 3.518f + t953 * 3.08f) + t9m32 * 1.55e6f * ex_(&
	    r__6);
/* 70--------ALPHA, PHOTON REACTIONS---------------------------------------------- */
/* .......H2(a,g)Li6.................(Caughlan-Fowler 1988) */
    r__1 = -7.423f / t913;
    r__2 = -8.228f / evolp1_1.t9;
    rates_1.f[24] = t9m23 * 30.1f * ex_(&r__1) * (t913 * .056f + 1.f - t923 * 
	    4.85f + evolp1_1.t9 * 8.85f - t943 * .585f - t953 * .584f) + 
	    t9m32 * 85.5f * ex_(&r__2);
/* .......H3(a,g)Li7.................(Smith-Kawano-Malaney 1992) */
    r__1 = -8.09f / t913;
    r__2 = -8.068f / t9e13;
    rates_1.f[25] = t9m23 * 303200.f * ex_(&r__1) * (t913 * .0516f + 1.f + 
	    t923 * .0229f + evolp1_1.t9 * .00828f - t943 * 3.28e-4f - t953 * 
	    3.01e-4f) + t9e56 * 510900.f * t9m32 * ex_(&r__2);
/* .......He3(a,g)Be7................(Smith-Kawano-Malaney 1992) */
    r__1 = -14.964f / t913;
    r__2 = -12.859f / t9f13;
    rates_1.f[26] = t9m23 * 4.817e6f * ex_(&r__1) * (t913 * .0325f + 1.f - 
	    t923 * .00104f - evolp1_1.t9 * 2.37e-4f - t943 * 8.11e-5f - t953 *
	     4.69e-5f) + t9f56 * 5.938e6f * t9m32 * ex_(&r__2);
/* 80--------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS------------------- */
/* .......H2(d,n)He3.................(Smith-Kawano-Malaney 1992) */
    r__1 = -4.259f / t913;
    rates_1.f[27] = t9m23 * 3.95e8f * ex_(&r__1) * (t913 * .098f + 1.f + t923 
	    * .765f + evolp1_1.t9 * .525f + t943 * .00961f + t953 * .0167f);
/* .......H2(d,p)H3..................(Smith-Kawano-Malaney 1992) */
    r__1 = -4.258f / t913;
    rates_1.f[28] = t9m23 * 4.17e8f * ex_(&r__1) * (t913 * .098f + 1.f + t923 
	    * .518f + evolp1_1.t9 * .355f - t943 * .01f - t953 * .018f);
/* .......H3(d,n)He4.................(Smith-Kawano-Malaney 1992) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .0754f;
    r__1 = -4.559f / t913 - r__2 * r__2;
    r__3 = -.4857f / evolp1_1.t9;
    rates_1.f[29] = t9m23 * 1.063e11f * ex_(&r__1) * (t913 * .092f + 1.f - 
	    t923 * .375f - evolp1_1.t9 * .242f + t943 * 33.82f + t953 * 
	    55.42f) + t9m23 * 8.047e8f * ex_(&r__3);
/* .......He3(d,p)He4................(Smith-Kawano-Malaney 1992) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .27f;
    r__1 = -7.144f / t913 - r__2 * r__2;
    r__3 = -1.762f / evolp1_1.t9;
    rates_1.f[30] = t9m23 * 5.021e10f * ex_(&r__1) * (t913 * .058f + 1.f + 
	    t923 * .603f + evolp1_1.t9 * .245f + t943 * 6.97f + t953 * 7.19f) 
	    + 5.212e8f / t912 * ex_(&r__3);
/* 90--------THREE PARTICLE REACTIONS--------------------------------------------- */
/* .......He3(He3,2p)He4.............(Caughlan-Fowler 1988) */
    r__1 = -12.276f / t913;
    rates_1.f[31] = t9m23 * 6.04e10f * ex_(&r__1) * (t913 * .034f + 1.f - 
	    t923 * .522f - evolp1_1.t9 * .124f + t943 * .353f + t953 * .213f);
/* .......Li7(d,na)He4...............(Caughlan-Fowler 1988) */
    r__1 = -10.259f / t913;
    rates_1.f[32] = t9m23 * 2.92e11f * ex_(&r__1);
/* .......Be7(d,pa)He4...............(Caughlan-Fowler 1988) */
    r__1 = -12.428f / t913;
    rates_1.f[33] = t9m23 * 1.07e12f * ex_(&r__1);
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     Smith, M., Kawano, L.H., and Malaney, R.A., 1992, submitted to Ap. J. */
/*     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5. */
/*     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data */
/*       Tables, 40, 283. */
/*     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247. */
} /* rate2_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int rate3_(void)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    extern doublereal ex_(real *);
    static real t9a, t912, t913, t923, t932, t943, t953, t915, t954, t9m1, 
	    t9a13, t9a56, t9m23, t9m32, t9m34, t9m15, t9m54;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] derivs */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Generates rate coefficients for reactions involving nuclides */
/*     up to A = 18. */
/* ----------PARAMETER. */
/* Number of nuclear reactions. */
/* ----------COMMON AREAS. */
/* Number of nuclides in calculation. */
/* Reaction rates. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION RATES. */
/* Evolution parameters. */
/* Forward reaction rate coefficients. */
/* ----------EVOLUTION PARAMETER. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------TEMPERATURE FACTORS-------------------------------------------------- */
/* Temperature of photons (units of 10* */
    d__1 = (doublereal) evolp1_1.t9;
    t913 = pow_dd(&d__1, &c_b69);
/* t9**(1/3) */
    t923 = t913 * t913;
/* t9**(2/3) */
    t943 = t923 * t923;
/* t9**(4/3) */
    t953 = evolp1_1.t9 * t923;
/* t9**(5/3) */
    t912 = sqrt(evolp1_1.t9);
/* t9**(1/2) */
    t932 = evolp1_1.t9 * t912;
/* t9**(3/2) */
    d__1 = (doublereal) evolp1_1.t9;
    t915 = pow_dd(&d__1, &c_b83);
/* t9**(1/5) */
    d__1 = (doublereal) evolp1_1.t9;
    t954 = pow_dd(&d__1, &c_b84);
/* t9**(5/4) */
    t9m1 = 1.f / evolp1_1.t9;
/* t9**(-1) */
    t9m23 = 1.f / t923;
/* t9**(-2/3) */
    t9m32 = 1.f / t932;
/* t9**(-3/2) */
    t9m34 = sqrt(t9m32);
/* t9**(-3/4) */
    t9m15 = 1.f / t915;
/* t9**(-1/5) */
    t9m54 = 1.f / t954;
/* t9**(-5/4) */
    t9a = evolp1_1.t9 / (evolp1_1.t9 / 15.1f + 1.f);
/* For reaction 53. */
    d__1 = (doublereal) t9a;
    t9a13 = pow_dd(&d__1, &c_b73);
/* t9a**(1/3) */
    d__1 = (doublereal) t9a;
    t9a56 = pow_dd(&d__1, &c_b74);
/* 20--------NEUTRON, PHOTON REACTIONS-------------------------------------------- */
/* .......Li7(n,g)Li8................(Wagoner 1969) */
/* t9a**(5/6) */
    r__1 = -2.62f / evolp1_1.t9;
    rates_1.f[34] = t9m32 * 9960.f * ex_(&r__1) + 4900.f;
/* .......B10(n,g)B11................(Wagoner 1969) */
    rates_1.f[35] = 66200.f;
/* .......B11(n,g)B12................(Malaney-Fowler 1989) */
    r__1 = -.223f / evolp1_1.t9;
    rates_1.f[36] = t9m32 * 2400.f * ex_(&r__1) + 729.f;
/* 30--------NEUTRON, PROTON REACTIONS-------------------------------------------- */
/* .......C11(n,p)B11................(Caughlan-Fowler 1988) */
    rates_1.f[37] = (1.f - t912 * .048f + evolp1_1.t9 * .01f) * 1.69e8f;
/* 40--------NEUTRON, ALPHA REACTIONS--------------------------------------------- */
/* .......B10(n,a)Li7................(Caughlan-Fowler 1988) */
    rates_1.f[38] = 5.07e8f;
/* 50--------PROTON, PHOTON REACTIONS--------------------------------------------- */
/* .......Be7(p,g)B8.................(Caughlan-Fowler 1988) */
    r__1 = -10.262f / t913;
    r__2 = -7.306f / evolp1_1.t9;
    rates_1.f[39] = t9m23 * 3.11e5f * ex_(&r__1) + t9m32 * 2530.f * ex_(&r__2)
	    ;
/* .......Be9(p,g)B10................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .846f;
    r__1 = -10.359f / t913 - r__2 * r__2;
    r__3 = -3.445f / evolp1_1.t9;
    r__4 = -10.62f / evolp1_1.t9;
    rates_1.f[40] = t9m23 * 1.33e7f * ex_(&r__1) * (t913 * .04f + 1.f + t923 *
	     1.52f + evolp1_1.t9 * .428f + t943 * 2.15f + t953 * 1.54f) + 
	    t9m32 * 96400.f * ex_(&r__3) + t9m32 * 2.72e6f * ex_(&r__4);
/* .......B10(p,g)C11................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 4.402f;
    r__1 = -12.062f / t913 - r__2 * r__2;
    r__3 = -12.041f / evolp1_1.t9;
    r__4 = -16.164f / evolp1_1.t9;
    rates_1.f[41] = t9m23 * 4.61e5f * ex_(&r__1) * (t913 * .035f + 1.f + t923 
	    * .426f + evolp1_1.t9 * .103f + t943 * .281f + t953 * .173f) + 
	    t9m32 * 1.93e5f * ex_(&r__3) + t9m32 * 11400.f * ex_(&r__4);
/* .......B11(p,g)C12................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .239f;
    r__1 = -12.095f / t913 - r__2 * r__2;
    r__3 = -1.733f / evolp1_1.t9;
    r__4 = -5.617f / evolp1_1.t9;
    rates_1.f[42] = t9m23 * 4.62e7f * ex_(&r__1) * (t913 * .035f + 1.f + t923 
	    * 3.f + evolp1_1.t9 * .723f + t943 * 9.91f + t953 * 6.07f) + 
	    t9m32 * 7890.f * ex_(&r__3) + t9m15 * 96800.f * ex_(&r__4);
/* .......C11(p,g)N12................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 1.627f;
    r__1 = -13.658f / t913 - r__2 * r__2;
    r__3 = -7.021f / evolp1_1.t9;
    rates_1.f[43] = t9m23 * 42400.f * ex_(&r__1) * (t913 * .031f + 1.f + t923 
	    * 3.11f + evolp1_1.t9 * .665f + t943 * 4.61f + t953 * 2.5f) + 
	    t9m32 * 8840.f * ex_(&r__3);
/* 60--------PROTON, NEUTRON REACTIONS-------------------------------------------- */
/* .......B12(p,n)C12................(Wagoner 1969) */
    r__1 = -12.12f / t913;
    rates_1.f[44] = t9m23 * 4.02e11f * ex_(&r__1);
/* 70--------PROTON, ALPHA REACTIONS---------------------------------------------- */
/* .......Be9(p,a)Li6................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .52f;
    r__1 = -10.359f / t913 - r__2 * r__2;
    r__3 = -3.046f / evolp1_1.t9;
    r__4 = -5.16f / evolp1_1.t9;
    rates_1.f[45] = t9m23 * 2.11e11f * ex_(&r__1) * (t913 * .04f + 1.f + t923 
	    * 1.09f + evolp1_1.t9 * .307f + t943 * 3.21f + t953 * 2.3f) + 
	    t9m1 * 4.51e8f * ex_(&r__3) + t9m34 * 6.7e8f * ex_(&r__4);
/* .......B10(p,a)Be7................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 4.402f;
    r__1 = -12.062f / t913 - r__2 * r__2;
    r__3 = -12.26f / evolp1_1.t9;
    rates_1.f[46] = t9m23 * 1.26e11f * ex_(&r__1) * (t913 * .035f + 1.f - 
	    t923 * .498f - evolp1_1.t9 * .121f + t943 * .3f + t953 * .184f) + 
	    t9m1 * 2.59e9f * ex_(&r__3);
/* .......B12(p,a)Be9................(Wagoner 1969) */
    r__1 = -12.12f / t913;
    rates_1.f[47] = t9m23 * 2.01e11f * ex_(&r__1);
/* 80--------ALPHA, PHOTON REACTIONS---------------------------------------------- */
/* .......Li6(a,g)B10................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 1.326f;
    r__1 = -18.79f / t913 - r__2 * r__2;
    r__3 = -3.484f / evolp1_1.t9;
    r__4 = -7.269f / evolp1_1.t9;
    rates_1.f[48] = t9m23 * 4.06e6f * ex_(&r__1) * (t913 * .022f + 1.f + t923 
	    * 1.54f + evolp1_1.t9 * .239f + t943 * 2.2f + t953 * .869f) + 
	    t9m32 * 1910.f * ex_(&r__3) + t9m1 * 10100.f * ex_(&r__4);
/* .......Li7(a,g)B11................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 4.195f;
    r__1 = -19.161f / t913 - r__2 * r__2;
    r__3 = -2.977f / evolp1_1.t9;
    r__4 = -6.227f / evolp1_1.t9;
    rates_1.f[49] = t9m23 * 3.55e7f * ex_(&r__1) * (t913 * .022f + 1.f + t923 
	    * .775f + evolp1_1.t9 * .118f + t943 * .884f + t953 * .342f) + 
	    t9m32 * 333.f * ex_(&r__3) + t9m1 * 4.1e4f * ex_(&r__4);
/* .......Be7(a,g)C11................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 4.769f;
    r__1 = -23.212f / t913 - r__2 * r__2;
    r__3 = -6.51f / evolp1_1.t9;
    r__4 = -10.039f / evolp1_1.t9;
    rates_1.f[50] = t9m23 * 8.45e7f * ex_(&r__1) * (t913 * .018f + 1.f + t923 
	    * .488f + evolp1_1.t9 * .061f + t943 * .296f + t953 * .095f) + 
	    t9m32 * 12500.f * ex_(&r__3) + t9m54 * 1.29e5f * ex_(&r__4);
/* 90--------ALPHA, PROTON REACTIONS---------------------------------------------- */
/* .......B8(a,p)C11.................(Wagoner 1969) */
    r__1 = -27.36f / t913;
    rates_1.f[51] = t9m23 * 1.08e15f * ex_(&r__1);
/* 100-------ALPHA, NEUTRON REACTIONS--------------------------------------------- */
/* .......Li8(a,n)B11................(Malaney-Fowler 1989) */
    r__1 = -19.461f / t9a13;
    rates_1.f[52] = t9a56 * 8.62e13f * t9m32 * ex_(&r__1);
/* .......Be9(a,n)C12................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .049f;
    r__1 = -23.87f / t913 - r__2 * r__2;
    r__3 = -1.184f / evolp1_1.t9;
    r__4 = -1.834f / evolp1_1.t9;
    r__5 = -4.179f / evolp1_1.t9;
    r__6 = -12.732f / evolp1_1.t9;
    rates_1.f[53] = t9m23 * 4.62e13f * ex_(&r__1) * (t913 * .017f + 1.f + 
	    t923 * 8.57f + evolp1_1.t9 * 1.05f + t943 * 74.51f + t953 * 
	    23.15f) + t9m32 * 7.34e-5f * ex_(&r__3) + t9m32 * .227f * ex_(&
	    r__4) + t9m32 * 1.26e5f * ex_(&r__5) + ex_(&r__6) * 2.4e8f;
/* 110-------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS------------------- */
/* .......Be9(d,n)B10................(original Wagoner code) */
    r__1 = 6.44f - 12.6f / t913;
    rates_1.f[54] = t9m23 * 7.16e8f * ex_(&r__1);
/* .......B10(d,p)B11................(original Wagoner code) */
    r__1 = 7.3f - 14.8f / t913;
    rates_1.f[55] = t9m23 * 9.53e8f * ex_(&r__1);
/* .......B11(d,n)C12................(original Wagoner code) */
    r__1 = 7.4f - 14.8f / t913;
    rates_1.f[56] = t9m23 * 1.41e9f * ex_(&r__1);
/* 120-------THREE PARTICLE REACTIONS--------------------------------------------- */
/* .......He4(an,g)Be9...............(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__1 = evolp1_1.t9;
    r__2 = -1.062f / evolp1_1.t9;
    rates_1.f[57] = 2.59e-6f / ((evolp1_1.t9 * .344f + 1.f) * (r__1 * r__1)) *
	     ex_(&r__2);
/* .......He4(2a,g)C12...............(Caughlan-Fowler 1988) */
    r__1 = -4.4027f / evolp1_1.t9;
    r__2 = -24.811f / evolp1_1.t9;
    rates_1.f[58] = t9m32 * 2.79e-8f * t9m32 * ex_(&r__1) + t9m32 * 1.35e-8f *
	     ex_(&r__2);
/* .......Li8(p,na)He4...............(original Wagoner code) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 2.53f;
    r__1 = -8.52f / t913 - r__2 * r__2;
    r__3 = -4.64f / evolp1_1.t9;
    rates_1.f[59] = t9m23 * 8.65e9f * ex_(&r__1) + t9m32 * 2.31e9f * ex_(&
	    r__3);
/* .......B8(n,pa)He4................(original Wagoner code) */
    rates_1.f[60] = 4.02e8f;
/* .......Be9(p,da)He4...............(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .52f;
    r__1 = -10.359f / t913 - r__2 * r__2;
    r__3 = -3.046f / evolp1_1.t9;
    r__4 = -5.8f / evolp1_1.t9;
    rates_1.f[61] = t9m23 * 2.11e11f * ex_(&r__1) * (t913 * .04f + 1.f + t923 
	    * 1.09f + evolp1_1.t9 * .307f + t943 * 3.21f + t953 * 2.3f) + 
	    t9m1 * 5.79e8f * ex_(&r__3) + t9m34 * 8.5e8f * ex_(&r__4);
/* .......B11(p,2a)He4...............(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 1.644f;
    r__1 = -12.095f / t913 - r__2 * r__2;
    r__3 = -1.734f / evolp1_1.t9;
    r__4 = -6.262f / evolp1_1.t9;
    r__5 = -14.154f / evolp1_1.t9;
    rates_1.f[62] = t9m23 * 2.2e12f * ex_(&r__1) * (t913 * .034f + 1.f + t923 
	    * .14f + evolp1_1.t9 * .034f + t943 * .19f + t953 * .116f) + 
	    t9m32 * 4.03e6f * ex_(&r__3) + t9m32 * 6.73e9f * ex_(&r__4) + 
	    t9m1 * 3.88e9f * ex_(&r__5);
/* .......C11(n,2a)He4...............(Wagoner 1969) */
    rates_1.f[63] = 1.58e8f;
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5. */
/*     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data */
/*       Tables, 40, 283. */
/*     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247. */
} /* rate3_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int rate4_(void)
{
    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10, r__11;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    extern doublereal ex_(real *);
    static real t9a, t9b, t912, t913, t923, t932, t943, t953, t935, t965, 
	    t938, t9a13, t9b13, t9a56, t9b56, t9m13, t9m23, t9m32, t9m65;

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] derivs */
/*     CALLS     - [function] ex */
/* ----------REMARKS. */
/*     Generates rate coefficients for rest of reactions. */
/* ----------PARAMETER. */
/* Number of nuclear reactions. */
/* ----------COMMON AREAS. */
/* Number of nuclides in calculation. */
/* Reaction rates. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION RATES. */
/* Evolution parameters. */
/* Forward reaction rate coefficients. */
/* ----------EVOLUTION PARAMETER. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------TEMPERATURE FACTORS-------------------------------------------------- */
/* Temperature of photons (units of 10* */
    d__1 = (doublereal) evolp1_1.t9;
    t913 = pow_dd(&d__1, &c_b69);
/* t9**(1/3) */
    t923 = t913 * t913;
/* t9**(2/3) */
    t943 = t923 * t923;
/* t9**(4/3) */
    t953 = evolp1_1.t9 * t923;
/* t9**(5/3) */
    t912 = sqrt(evolp1_1.t9);
/* t9**(1/2) */
    t932 = evolp1_1.t9 * t912;
/* t9**(3/2) */
    d__1 = (doublereal) evolp1_1.t9;
    t935 = pow_dd(&d__1, &c_b89);
/* t9**(3/5) */
    d__1 = (doublereal) evolp1_1.t9;
    t965 = pow_dd(&d__1, &c_b90);
/* t9**(6/5) */
    d__1 = (doublereal) evolp1_1.t9;
    t938 = pow_dd(&d__1, &c_b91);
/* t9**(3/8) */
    t9m13 = 1.f / t913;
/* t9**(1/3) */
    t9m23 = 1.f / t923;
/* t9**(-2/3) */
    t9m32 = 1.f / t932;
/* t9**(-3/2) */
    t9m65 = 1.f / t965;
/* t9**(-6/5) */
    d__1 = (doublereal) (evolp1_1.t9 * .0478f + 1.f);
    t9a = evolp1_1.t9 / (evolp1_1.t9 * .0478f + 1.f + t953 * .00756f / pow_dd(
	    &d__1, &c_b72));
/* For reaction 82. */
    d__1 = (doublereal) t9a;
    t9a13 = pow_dd(&d__1, &c_b69);
/* t9a**(1/3) */
    d__1 = (doublereal) t9a;
    t9a56 = pow_dd(&d__1, &c_b94);
/* t9a**(5/6) */
    d__1 = (doublereal) (evolp1_1.t9 * .0776f + 1.f);
    t9b = evolp1_1.t9 / (evolp1_1.t9 * .0776f + 1.f + t953 * .0264f / pow_dd(&
	    d__1, &c_b72));
/* For reaction 84. */
    d__1 = (doublereal) t9b;
    t9b13 = pow_dd(&d__1, &c_b69);
/* t9b**(1/3) */
    d__1 = (doublereal) t9b;
    t9b56 = pow_dd(&d__1, &c_b94);
/* 20--------NEUTRON, PHOTON REACTIONS-------------------------------------------- */
/* .......C12(n,g)C13................(Wagoner 1969) */
/* t9b**(5/6) */
    rates_1.f[64] = 450.f;
/* .......C13(n,g)C14................(Wagoner 1969) */
    r__1 = -1.67f / evolp1_1.t9;
    rates_1.f[65] = t9m32 * 2.38e5f * ex_(&r__1) + 119.f;
/* .......N14(n,g)N15................(Wagoner 1969) */
    rates_1.f[66] = 9940.f;
/* 30--------NEUTRON, PROTON REACTIONS-------------------------------------------- */
/* .......N13(n,p)C13................(Caughlan-Fowler 1988) */
    rates_1.f[67] = (1.f - t912 * .167f + evolp1_1.t9 * .037f) * 1.88e8f;
/* .......N14(n,p)C14................(Caughlan-Fowler 1988) */
    r__1 = -4.983f / evolp1_1.t9;
    rates_1.f[68] = (t912 * .361f + 1.f + evolp1_1.t9 * .502f) * 2.39e5f + 
	    1.112e8f / t912 * ex_(&r__1);
/* .......O15(n,p)N15................(Caughlan-Fowler 1988) */
    rates_1.f[69] = (t912 * .452f + 1.f - evolp1_1.t9 * .191f) * 3.5e8f;
/* 40--------NEUTRON, ALPHA REACTIONS--------------------------------------------- */
/* .......O15(n,a)C12................(Caughlan-Fowler 1988) */
    rates_1.f[70] = (t912 * .188f + 1.f + evolp1_1.t9 * .015f) * 3.5e7f;
/* 50--------PROTON, PHOTON REACTIONS--------------------------------------------- */
/* .......C12(p,g)N13................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 1.5f;
    r__1 = -13.69f / t913 - r__2 * r__2;
    r__3 = -4.925f / evolp1_1.t9;
    r__4 = -18.179f / evolp1_1.t9;
    rates_1.f[71] = t9m23 * 2.04e7f * ex_(&r__1) * (t913 * .03f + 1.f + t923 *
	     1.19f + evolp1_1.t9 * .254f + t943 * 2.06f + t953 * 1.12f) + 
	    t9m32 * 1.08e5f * ex_(&r__3) + t9m32 * 2.15e5f * ex_(&r__4);
/* .......C13(p,g)N14................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 2.f;
    r__1 = -13.717f / t913 - r__2 * r__2;
    r__3 = -5.701f / evolp1_1.t9;
    rates_1.f[72] = t9m23 * 8.01e7f * ex_(&r__1) * (t913 * .03f + 1.f + t923 *
	     .958f + evolp1_1.t9 * .204f + t943 * 1.39f + t953 * .753f) + 
	    t9m65 * 1.21e6f * ex_(&r__3);
/* .......C14(p,g)N15................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 5.721f;
    r__1 = -13.741f / t913 - r__2 * r__2;
    r__3 = -3.811f / evolp1_1.t9;
    r__4 = -4.739f / evolp1_1.t9;
    rates_1.f[73] = t9m23 * 6.8e6f * ex_(&r__1) * (t913 * .03f + 1.f + t923 * 
	    .503f + evolp1_1.t9 * .107f + t943 * .213f + t953 * .115f) + 
	    t9m32 * 5360.f * ex_(&r__3) + t9m13 * 98200.f * ex_(&r__4);
/* .......N13(p,g)O14................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 1.191f;
    r__1 = -15.202f / t913 - r__2 * r__2;
    r__3 = -6.348f / evolp1_1.t9;
    rates_1.f[74] = t9m23 * 4.04e7f * ex_(&r__1) * (t913 * .027f + 1.f - t923 
	    * .803f - evolp1_1.t9 * .154f + t943 * 5.f + t953 * 2.44f) + 
	    t9m32 * 2.43e5f * ex_(&r__3);
/* .......N14(p,g)O15................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 3.294f;
    r__1 = -15.228f / t913 - r__2 * r__2;
    r__3 = -3.011f / evolp1_1.t9;
    r__4 = -12.53f / evolp1_1.t9;
    rates_1.f[75] = t9m23 * 4.9e7f * ex_(&r__1) * (t913 * .027f + 1.f - t923 *
	     .778f - evolp1_1.t9 * .149f + t943 * .261f + t953 * .127f) + 
	    t9m32 * 2370.f * ex_(&r__3) + ex_(&r__4) * 21900.f;
/* .......N15(p,g)O16................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .45f;
    r__1 = -15.251f / t913 - r__2 * r__2;
    r__3 = -3.328f / evolp1_1.t9;
    r__4 = -4.665f / evolp1_1.t9;
    r__5 = -11.048f / evolp1_1.t9;
    rates_1.f[76] = t9m23 * 9.78e8f * ex_(&r__1) * (t913 * .027f + 1.f + t923 
	    * .219f + evolp1_1.t9 * .042f + t943 * 6.83f + t953 * 3.32f) + 
	    t9m32 * 11100.f * ex_(&r__3) + t9m32 * 14900.f * ex_(&r__4) + 
	    t9m32 * 3.8e6f * ex_(&r__5);
/* 60--------PROTON, ALPHA REACTIONS---------------------------------------------- */
/* .......N15(p,a)C12................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .522f;
    r__1 = -15.251f / t913 - r__2 * r__2;
    r__3 = -3.676f / evolp1_1.t9;
    r__4 = -8.926f / evolp1_1.t9;
    r__5 = -7.721f / evolp1_1.t9;
    r__6 = -11.418f / evolp1_1.t9;
    rates_1.f[77] = t9m23 * 1.08e12f * ex_(&r__1) * (t913 * .027f + 1.f + 
	    t923 * 2.62f + evolp1_1.t9 * .501f + t943 * 5.36f + t953 * 2.6f) 
	    + t9m32 * 1.19e8f * ex_(&r__3) + 5.41e8f / t912 * ex_(&r__4) + 
	    t9m32 * 4.72e7f * ex_(&r__5) + t9m32 * 2.2e8f * ex_(&r__6);
/* 70--------ALPHA, PHOTON REACTIONS---------------------------------------------- */
/* .......C12(a,g)O16................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__1 = evolp1_1.t9;
/* Computing 2nd power */
    r__3 = evolp1_1.t9 / 3.496f;
    r__2 = -32.12f / t913 - r__3 * r__3;
/* Computing 2nd power */
    r__4 = t9m23 * .0489f + 1.f;
/* Computing 2nd power */
    r__5 = evolp1_1.t9;
/* Computing 2nd power */
    r__6 = t9m23 * .2654f + 1.f;
    r__7 = -32.12f / t913;
    r__8 = -27.499f / evolp1_1.t9;
/* Computing 5th power */
    r__9 = evolp1_1.t9, r__10 = r__9, r__9 *= r__9;
    r__11 = -15.541f / evolp1_1.t9;
    rates_1.f[78] = 1.04e8f / (r__1 * r__1) * ex_(&r__2) / (r__4 * r__4) + 
	    1.76e8f / (r__5 * r__5) / (r__6 * r__6) * ex_(&r__7) + t9m32 * 
	    1250.f * ex_(&r__8) + r__10 * (r__9 * r__9) * .0143f * ex_(&r__11)
	    ;
/* 80--------ALPHA, PROTON REACTIONS---------------------------------------------- */
/* .......B10(a,p)C13................(Wagoner 1969) */
    r__1 = -27.99f / t913;
    rates_1.f[79] = t9m23 * 9.6e14f * ex_(&r__1);
/* .......B11(a,p)C14................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .347f;
    r__1 = -28.234f / t913 - r__2 * r__2;
    r__3 = -2.827f / evolp1_1.t9;
    r__4 = -5.178f / evolp1_1.t9;
    r__5 = -11.617f / evolp1_1.t9;
    rates_1.f[80] = t9m23 * 5.37e11f * ex_(&r__1) * (t913 * .015f + 1.f + 
	    t923 * 5.575f + evolp1_1.t9 * .576f + t943 * 15.888f + t953 * 
	    4.174f) + t9m32 * .00544f * ex_(&r__3) + t9m32 * 336.f * ex_(&
	    r__4) + 5.32e6f / t938 * ex_(&r__5);
/* .......C11(a,p)N14................(Caughlan-Fowler 1988) */
    r__1 = -31.883f / t9a13;
    rates_1.f[81] = t9a56 * 7.15e15f * t9m32 * ex_(&r__1);
/* .......N12(a,p)O15................(Caughlan-Fowler 1988) */
    r__1 = -35.6f / t913;
    rates_1.f[82] = t9m23 * 5.59e16f * ex_(&r__1);
/* .......N13(a,p)O16................(Caughlan-Fowler 1988) */
    r__1 = -35.829f / t9b13;
    rates_1.f[83] = t9b56 * 3.23e17f * t9m32 * ex_(&r__1);
/* 90--------ALPHA, NEUTRON REACTIONS--------------------------------------------- */
/* .......B10(a,n)N13................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 9.589f;
    r__1 = -27.989f / t913 - r__2 * r__2;
    rates_1.f[84] = t9m23 * 1.2e13f * ex_(&r__1);
/* .......B11(a,n)N14................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / .14f;
    r__1 = -28.234f / t913 - r__2 * r__2;
    r__3 = -2.827f / evolp1_1.t9;
    r__4 = -5.178f / evolp1_1.t9;
    r__5 = -8.596f / evolp1_1.t9;
    rates_1.f[85] = t9m23 * 6.97e12f * ex_(&r__1) * (t913 * .015f + 1.f + 
	    t923 * 8.115f + evolp1_1.t9 * .838f + t943 * 39.804f + t953 * 
	    10.456f) + t9m32 * 1.79f * ex_(&r__3) + t9m32 * 1710.f * ex_(&
	    r__4) + t935 * 4.49e6f * ex_(&r__5);
/* .......B12(a,n)N15................(Wagoner 1969) */
    r__1 = -28.45f / t913;
    rates_1.f[86] = t9m23 * 3.04e15f * ex_(&r__1);
/* .......C13(a,n)O16................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    r__2 = evolp1_1.t9 / 1.284f;
    r__1 = -32.329f / t913 - r__2 * r__2;
    r__3 = -9.373f / evolp1_1.t9;
    r__4 = -11.873f / evolp1_1.t9;
    r__5 = -20.409f / evolp1_1.t9;
    r__6 = -29.283f / evolp1_1.t9;
    rates_1.f[87] = t9m23 * 6.77e15f * ex_(&r__1) * (t913 * .013f + 1.f + 
	    t923 * 2.04f + evolp1_1.t9 * .184f) + t9m32 * 3.82e5f * ex_(&r__3)
	     + t9m32 * 1.41e6f * ex_(&r__4) + t9m32 * 2e9f * ex_(&r__5) + 
	    t9m32 * 2.92e9f * ex_(&r__6);
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data */
/*       Tables, 40, 283. */
/*     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247. */
} /* rate4_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* ----------PARAMETERS. */
/* Number of nuclear reactions. */
/* ----------COMMON AREAS. */
/* Number of nuclides in calculation. */
/* Reaction parameter */
/* Default comp param */
/* Default model para */
/* Default variationl */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION PARAMETERS VALUES. */
/* Nuclide data. */
/* ----------DEFAULT COMPUTATION PARAMETERS. */
/* Reaction parameters. */
/* Default time step limiting constant. */
/* Default time step limiting constant. */
/* Default initial temperature (in 10** */
/* Default final temperature (in 10**9 */
/* Default smallest abundances allowed. */
/* ----------DEFAULT MODEL PARAMETERS. */
/* Default accumulation increment. */
/* c0(1) is default variation of grav c */
/* c0(2) is default neutron half-life. */
/* c0(3) is default number of neutrinos */
/* Default cosmological constant. */
/* ----------DEFAULT VARIATIONAL PARAMETERS. */
/* Default neutrino degeneracy paramete */
/* Default initial time step. */
/* ----------NUCLIDE DATA. */
/* Default baryon-to-photon ratio. */
/* Atomic number of nuclide. */
/* Charge of nuclide. */
/* ==============================DATA DIVISION==================================== */
/*    Nuclide and corresponding number */
/*    -------------------------------- */
/*    1) N         7) Li6      13) B10      19) C13      25) O15 */
/*    2) P         8) Li7      14) B11      20) N13      26) O16 */
/*    3) H2        9) Be7      15) C11      21) C14 */
/*    4) H3       10) Li8      16) B12      22) N14 */
/*    5) He3      11) B8       17) C12      23) O14 */
/*    6) He4      12) Be9      18) N12      24) N15 */
/* -----------NUCLIDE DATA. */
/* Mass excess of nuclide. */
/* ----------REACTION RATE COEFFICIENTS (Ref 1). */
/*              reac# type n1 n2 n3 n4 rev-coeff q-value */
/*              ----  ---- -- -- -- -- --------- ------- */
/* N->P */
/* H3->He3 */
/* Li8->2He4 */
/* B12->C12 */
/* C14->N14 */
/* B8->2He4 */
/* C11->B11 */
/* N12->C12 */
/* N13->C13 */
/* O14->N14 */
/* O15->N15 */
/*              reac# type n1 n2 n3 n4 rev-coeff q-value */
/*              ----  ---- -- -- -- -- --------- ------- */
/* H(n,g)H2 */
/* H2(n,g)H3 */
/* He3(n,g)He4 */
/* Li6(n,g)Li7 */
/* He3(n,p)H3 */
/* Be7(n,p)Li7 */
/* Li6(n,a)H3 */
/* Be7(n,a)He4 */
/* H2(p,g)He3 */
/* H3(p,g)He4 */
/* Li6(p,g)Be7 */
/*              reac# type n1 n2 n3 n4 rev-coeff q-value */
/*              ----  ---- -- -- -- -- --------- ------- */
/* Li6(p,a)He3 */
/* Li7(p,a)He4 */
/* H2(a,p)Li6 */
/* H3(a,p)Li7 */
/* He3(a,p)Be7 */
/* H2(d,p)He3 */
/* H2(d,n)H3 */
/* H3(d,n)He4 */
/* He3(d,p)He4 */
/* He3(He3,2p)H */
/* Li7(d,na)He4 */
/*              reac# type n1 n2 n3 n4 rev-coeff q-value */
/*              ----  ---- -- -- -- -- --------- ------- */
/* Be7(d,pa)He4 */
/* Li7(n,g)Li8 */
/* B10(n,g)B11 */
/* B11(n,g)B12 */
/* C11(n,p)B11 */
/* B10(n,a)Li7 */
/* Be7(p,g)B8 */
/* Be9(p,g)B10 */
/* B10(p,g)C11 */
/* B11(p,g)C12 */
/* C11(p,g)N12 */
/*              reac# type n1 n2 n3 n4 rev-coeff q-value */
/*              ----  ---- -- -- -- -- --------- ------- */
/* B12(p,n)C12 */
/* Be9(p,a)Li6 */
/* B10(p,a)Be7 */
/* B12(p,a)Be9 */
/* Li6(a,g)B10 */
/* Li7(a,g)B11 */
/* Be7(a,g)C11 */
/* B8(a,p)C11 */
/* Li8(a,n)B11 */
/* Be9(a,n)C12 */
/* Be9(d,n)B10 */
/*              reac# type n1 n2 n3 n4 rev-coeff q-value */
/*              ----  ---- -- -- -- -- --------- ------- */
/* B10(d,p)B11 */
/* B11(d,n)C12 */
/* He4(an,g)Be9 */
/* He4(2a,g)C12 */
/* Li8(p,na)He4 */
/* B8(n,pa)He4 */
/* Be9(p,da)He4 */
/* B11(p,2a)Be4 */
/* C11(n,2a)He4 */
/* C12(n,g)C13 */
/* C13(n,g)C14 */
/*              reac# type n1 n2 n3 n4 rev-coeff q-value */
/*              ----  ---- -- -- -- -- --------- ------- */
/* N14(n,g)N15 */
/* N13(n,p)C13 */
/* N14(n,p)C14 */
/* O15(n,p)N15 */
/* O15(n,a)C12 */
/* C12(p,g)N13 */
/* C13(p,g)N14 */
/* C14(p,g)N15 */
/* N13(p,g)O14 */
/* N14(p,g)O15 */
/* N15(p,g)O16 */
/*              reac# type n1 n2 n3 n4 rev-coeff q-value */
/*              ----  ---- -- -- -- -- --------- ------- */
/* ----------DEFAULT COMPUTATION PARAMETERS. */
/* N15(p,a)C12 */
/* C12(a,g)O16 */
/* B10(a,p)C13 */
/* B11(a,p)C14 */
/* C11(a,p)N14 */
/* N12(a,p)O15 */
/* N13(a,p)O16 */
/* B10(a,n)N13 */
/* B11(a,n)N14 */
/* B12(a,n)N15 */
/* C13(a,n)O16 */
/* Default time step limiting constant. */
/* Default time step limiting constant. */
/* Default initial temperature. */
/* Default final temperature. */
/* Default smallest abundances allowed. */
/* -----------DEFAULT MODEL PARAMETERS. */
/* Default accumulation increment. */
/* Default variation of 3 parameters. */
/* Default cosmological constant. */
/* -----------DEFAULT VARIATIONAL PARAMETERS. */
/* Default neutrino degeneracy paramete */
/* Default initial time step. */
/* Default baryon-to-photon ratio. */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     Changes (to run inder unix f77): */
/*     ------------------------------- */
/*     Programme renamed from nucint.for to newint.f */
/*     COMMON /therm/  -> COMMON /thermcb/ */
/*     COMMON /bessel/ -> COMMON /besselcb/ */
/*     COMMON /check/  -> COMMON /checkcb/ */
/*     COMMON /time/   -> COMMON /ttime/ */

/*     Updated Correction to helium-4 abundance (using fitted rates and */
/*     *smallest* step-size) - Sarkar, Rep. Prog Phys. 59, 1493 (1996): */

/*     Y_p = Y_p - 0.0025 -> Y_p = Y_p - 0.0003 */


/*     output nucint.dat  -> newint.dat */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int check_(void)
{
    /* Format strings */
    static char fmt_200[] = "(7(e13.5,\002 \002))";

    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);

    /* Fortran I/O blocks */
    static cilist io___188 = { 0, 3, 0, fmt_200, 0 };


/* ----------REMARKS. */
/*     This is an interface subroutine, */
/*     a flexible module which allows user to manipulate physical quantities */
/*     of interest at certain key points during the computer run. */
/*     Included within this subroutine is a roster of all global variables */
/*     and their respective COMMON areas. */
/*     Current instructions accumulate abundances of deuterium, helium-3, */
/*     helium-4, and lithium-7 for eventual plotting, taking into account */
/*     the contribution of beryllium-7 to lithium-7 and tritium to helium-3. */
/* ----------PARAMETERS. */
/* Number of nuclear reactions. */
/* Number of variables to be evolved. */
/* Number of nuclides in calculation. */
/* ----------COMMON AREAS. */
/* Maximum # of lines to be printed. */
/* Reaction parameter */
/* Reaction parameter */
/* Reaction rates. */
/* Evolution paramete */
/* Evolution paramete */
/* Evolution paramete */
/* Default comp param */
/* Computation parame */
/* Default model para */
/* Model parameters. */
/* Default variationl */
/* Variational parame */
/* Time variables. */
/* Dynamic variables. */
/* Energy densities. */
/* Linear eqn coeffi */
/* Nuclide data. */
/* Eval function bl */
/* Eval function bm(z */
/* Eval function bn(z */
/* Coefficients K. */
/* Flags,counters. */
/* Computation locati */
/* Output data. */
/* Neutrino parameter */
/* Run options. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION PARAMETER VALUES. */
/* Output option. */
/* ----------REACTION PARAMETER NAMES. */
/* Reaction parameters. */
/* Reaction type code (1-11). */
/* Incoming nuclide type (1-26). */
/* Incoming light nuclide type (1-6). */
/* Outgoing light nuclide type (1-6). */
/* Outgoing nuclide type (1-26). */
/* Reverse reaction coefficient. */
/* ----------REACTION RATES. */
/* Energy released in reaction (in 10** */
/* Forward reaction rate coefficients. */
/* ----------EVOLUTION PARAMETERS. */
/* Reverse reaction rate coefficients. */
/* Temperature of photons (units of 10* */
/* Defined by hv = M(atomic)n(baryon)/t */
/* Chemical potential of electron. */
/* ----------EVOLUTION PARAMETERS (DERIVATIVES). */
/* Relative number abundances. */
/* Change in temperature. */
/* Change in hv. */
/* Change in chemical potential. */
/* ----------EVOLUTION PARAMETERS (ORIGINAL VALUES). */
/* Change in relative number abundances */
/* ----------DEFAULT COMPUTATION PARAMETERS. */
/* Rel # abundances at end of 1st R-K l */
/* Default cy. */
/* Default ct. */
/* Default t9i. */
/* Default t9f. */
/* Default ytmin. */
/* ----------COMPUTATION PARAMETERS. */
/* Default accumulation increment. */
/* Time step limiting constant on abund */
/* Time step limiting constant on tempe */
/* Initial temperature (in 10**9 K). */
/* Final temperature (in 10**9 K). */
/* Smallest abundances allowed. */
/* ----------DEFAULT MODEL PARAMETERS. */
/* Accumulation increment. */
/* Default c. */
/* Default cosmological constant. */
/* ----------EARLY UNIVERSE MODEL PARAMETERS. */
/* Default neutrino degeneracy paramete */
/* Gravitational constant. */
/* Neutron lifetime (sec). */
/* Number of neutrino species. */
/* c(1) is variation of gravitational c */
/* c(2) is neutron half-life (min). */
/* c(3) is number of neutrino species. */
/* Cosmological constant. */
/* Neutrino degeneracy parameters. */
/* xi(1) is e neutrino degeneracy param */
/* xi(2) is m neutrino degeneracy param */
/* xi(3) is t neutrino degeneracy param */
/* ----------DEFAULT VARIATIONAL PARAMETERS. */
/* Fierz parameter */
/* Default initial time step. */
/* ----------VARIATIONAL PARAMETERS. */
/* Default baryon-to-photon ratio. */
/* Initial time step. */
/* ----------TIME VARIABLES. */
/* Baryon-to-photon ratio. */
/* Time. */
/* Time step. */
/* ----------DYNAMIC VARIABLES. */
/* (1/t9)*d(t9)/d(t). */
/* Thermodynamic variables (energy dens */
/* ----------ENERGY DENSITIES. */
/* Expansion rate of the universe. */
/* Initial electron neutrino energy den */
/* Initial baryon energy density. */
/* Baryon energy density. */
/* ----------MATRIX COEFFICIENTS FOR LINEAR EQUATION. */
/* Baryon energy density (ratio to init */
/* Relates y(t+dt) to y(t). */
/* Contains y0 in inverse order. */
/* ----------NUCLIDE DATA. */
/* yy in reverse order. */
/* Atomic number of nuclide. */
/* Charge of nuclide. */
/* ----------EVALUATION OF FUNCTIONS bl,bm,bn. */
/* Mass excess of nuclide. */
/* Evaluation of function bl(z). */
/* Evaluation of function bm(z). */
/* ----------EVALUATION OF MODIFIED BESSEL FUNCTIONS. */
/* Evaluation of function bn(z). */
/* ----------FLAGS AND COUNTERS. */
/* Values k0(r),k1(r),k2(r),k3(r),k4(r) */
/* Indicates if output buffer printed. */
/* # total iterations for particular mo */
/* # iterations after outputing a line. */
/* # times accumulated in output buffer */
/* ----------COMPUTATION LOCATION. */
/* Indicates if gaussian elimination fa */
/* ----------OUTPUT ARRAYS. */
/* Time check. */
/* Nuclide mass fractions. */
/* Thermodynamic variables. */
/* Temperature (in units of 10**9 K). */
/* Time. */
/* Time step. */
/* Baryon to photon ratio. */
/* ----------NEUTRINO PARAMETERS. */
/* Expansion rate. */
/* Temperature (in units of MeV). */
/* Neutrino temperature (in units of Me */
/* Neutrino temperature. */
/* Normalizing constant. */
/* Neutrino energy density. */
/* ----------RUN OPTION. */
/* Type of neutrino. */
/* Run network size. */
/* Number of nuclides in computation. */
/* ----------OUTPUT FILE STATUS. */
/* Number of reactions in computation. */
/* Number of output requests. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------OPEN FILE------------------------------------------------------------ */
/* Indicates if output file used. */
    if (checkcb_1.itime == 1) {
/* Beginning of program. */
	o__1.oerr = 0;
	o__1.ounit = 3;
	o__1.ofnmlen = 10;
	o__1.ofnm = "newint.dat";
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }
/* 20--------WRITE INTO FILE------------------------------------------------------ */
    if (checkcb_1.itime == 8) {
/* Right after a run. */
	outdat_1.xout[flags_1.it + 279] += outdat_1.xout[flags_1.it + 319];
/* Add beryllium to lithium. */
	outdat_1.xout[flags_1.it + 159] += outdat_1.xout[flags_1.it + 119];
/* Add tritium to helium-3. */
	outdat_1.xout[flags_1.it + 199] += -3e-4f;
/* my correction for fitted rates+coarse steps */
	s_wsfe(&io___188);
	do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.etaout[flags_1.it - 1], (ftnlen)
		sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.xout[flags_1.it + 79], (ftnlen)sizeof(
		real));
	do_fio(&c__1, (char *)&outdat_1.xout[flags_1.it + 159], (ftnlen)
		sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.xout[flags_1.it + 199], (ftnlen)
		sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.xout[flags_1.it + 279], (ftnlen)
		sizeof(real));
	e_wsfe();
/* Output N_nu, tau_n, eta, H2, He3, He4, an */
    }
/* 30--------CLOSE FILE----------------------------------------------------------- */
    if (checkcb_1.itime == 10) {
/* End of program. */
	cl__1.cerr = 0;
	cl__1.cunit = 3;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     1) D.A. Dicus, E.W. Kolb, A.M. Gleeson, E.C.G. Sudarshan, V.L. Teplitz, */
/*        M.S. Turner, Phys. Rev. D26 (1982) 2694. (Rad corr, Coulomb Corr) */
/*     2) D. Seckel, Bartol preprint BA-93-16; G. Guyk and M.S. Turner, */
/*         FERMILAB preprint FERMILAB-Pub-93/181-A. (Nucleon mass) */
/*     4) S. Dodelson and M.S. Turner, Phys. Rev. D46 (1992) 3372; B. Fields, */
/*         S. Dodelson and M.S. Turner, Phys. Rev. D47 (1993) 4309. (Nu heating) */
/*     5) S. Sarkar, Rep. Prog Phys. 59 (1996) 1493 (review) */
} /* check_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*     Changes (to run inder unix f77): */
/*     ------------------------------- */
/* Programme renamed from nuc123.for to new123.f */
/* ir=1 -> ir=5 ... input unit number */
/* iw=1 -> iw=6 ... output unit number */
/* COMMON /check/ itime -> COMMON /checkcb/ itime */
/* COMMON /time/        -> COMMON /ttime/ */
/* output nuc123.dat    -> new123.dat */

/* ========================IDENTIFICATION DIVISION================================ */
/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_1000[] = "(6(/),2(\002 \002,4x,\002NN\002,6x,\002NN  U"
	    "U\002,6x,\002UU\002,4x,8(\002C\002),6x,\00211\002,8x,6(\0022\002"
	    "),6x,6(\0023\002),/),2(\002 \002,4x,\002NN\002,6x,\002NN  UU\002"
	    ",6x,\002UU  CC\002,12x,\0021111\002,6x,\00222\002,6x,\00222  3"
	    "3\002,6x,\00233\002,/),2(\002 \002,4x,\002NNNN    NN  UU\002,6x"
	    ",\002UU  CC\002,14x,\00211\002,14x,\00222\002,10x,\00233\002,/),"
	    "2(\002 \002,4x,\002NN  NN  NN  UU\002,6x,\002UU  CC\002,14x,\002"
	    "11\002,12x,\00222\002,10x,\00233\002,/),2(\002 \002,4x,\002NN   "
	    " NNNN  UU\002,6x,\002UU  CC\002,14x,\00211\002,10x,\00222\002,14"
	    "x,\00233\002,/),2(\002 \002,4x,\002NN\002,6x,\002NN  UU\002,6x"
	    ",\002UU  CC\002,14x,\00211\002,8x,\00222\002,8x,\00233\002,6x"
	    ",\00233\002,/),2(\002 \002,4x,\002NN\002,6x,\002NN  \002,10(\002U"
	    "\002),4x,8(\002C\002),4x,6(\0021\002),4x,10(\0022\002),4x,6(\0023"
	    "\002),/),/,\002 \002,26x,\002MODIFIED APRIL 1994\002,///,\002"
	    " \002,\002(Press <RETURN> to continue): \002,$)";
    static char fmt_3000[] = "(8(/),\002 \002,32x,\002MENU SELECTION\002,/"
	    ",\002 \002,32x,\002---- ---------\002,//,\002 \002,24x,\0021. HE"
	    "LP\002,/,\002 \002,24x,\0022. SET COMPUTATION PARAMETERS\002,/"
	    ",\002 \002,24x,\0023. SET MODEL PARAMETERS\002,/,\002 \002,24x"
	    ",\0024. RUN\002,/,\002 \002,24x,\0025. OUTPUT\002,/,\002 \002,24"
	    "x,\0026. EXIT\002,8(/),\002 \002,24x,\002Enter selection (1-6):"
	    " \002,$)";
    static char fmt_3001[] = "(i1)";

    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_wsfe(cilist *), e_wsfe(void), s_rsle(cilist *),
	     e_rsle(void), s_rsfe(cilist *), do_fio(integer *, char *, ftnlen)
	    , e_rsfe(void), f_clos(cllist *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int run_(void), help_(void);
    static integer inum;
    extern /* Subroutine */ int check_(void), setcom_(void), setmod_(void), 
	    output_(void);

    /* Fortran I/O blocks */
    static cilist io___189 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___191 = { 0, 5, 0, 0, 0 };
    static cilist io___192 = { 0, 6, 0, fmt_3000, 0 };
    static cilist io___193 = { 0, 5, 0, fmt_3001, 0 };


/* ----------LINKAGES. */
/*     CALLED BY - none */
/*     CALLS     - [subroutine] help, setcom, setmod, run, output */
/* ----------REMARKS. */
/*     Control program - */
/*       Offers user the main menu and channels through to various options. */
/*     Implementation - */
/*       To run this program, new123.f must be linked with nuccom.f */
/*       (containing the computation subroutines), nucrat.f (with the */
/*       reaction rates), and newint.f (with an interface subroutine). */
/*       This program has been written to be compatible with */
/*       ANSI FORTRAN-77 with the exception of the END DO statement */
/*       used to limit the number of statement labels. */
/*       This code was modified on the DEC/AXP system. */
/*     Notes - */
/*       The program utilizes Wagoner's code as the core of the computational */
/*       routines. */
/*     Documentation - */
/*       Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A, */
/*       Kellogg Radiation Lab preprint OAP-714. */
/*     Copy - */
/*       Version 4.1 (December 1991) */
/* ----------PARAMETERS. */
/* Input unit number (previous value = */
/* Output unit number (previous value = */
/* Number of nuclear reactions. */
/* ----------COMMON AREAS. */
/* Number of nuclides in calculation. */
/* Reaction parameter */
/* Reaction parameter */
/* Reaction rates. */
/* Default comp param */
/* Computation parame */
/* Default model para */
/* Model parameters. */
/* Default variationl */
/* Variational parame */
/* Computation locati */
/* Run options. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------REACTION PARAMETERS FROM BLOCK DATA. */
/* Output option. */
/* ----------REACTION PARAMETERS. */
/* Reaction parameters. */
/* Reaction type code (1-11). */
/* Incoming nuclide type (1-26). */
/* Incoming light nuclide type (1-6). */
/* Outgoing light nuclide type (1-6). */
/* Outgoing nuclide type (1-26). */
/* Reverse reaction coefficient. */
/* ----------REACTION RATES. */
/* Energy released in reaction. */
/* Forward reaction rate coefficients. */
/* ----------DEFAULT COMPUTATION PARAMETERS. */
/* Reverse reaction rate coefficients. */
/* Default cy. */
/* Default ct. */
/* Default t9i. */
/* Default t9f. */
/* Default ytmin. */
/* ----------COMPUTATIONAL PARAMETERS. */
/* Default accumulation increment. */
/* Time step limiting constant on abund */
/* Time step limiting constant on tempe */
/* Initial temperature (in 10**9 K). */
/* Final temperature (in 10**9 k). */
/* Smallest abundances allowed. */
/* ----------DEFAULT MODEL PARAMETERS. */
/* Accumulation increment. */
/* Default c. */
/* Default cosmological constant. */
/* ----------EARLY UNIVERSE MODEL PARAMETERS. */
/* Default neutrino degeneracy paramete */
/* c(1) is variation of gravitational c */
/* c(2) is neutron lifetime (sec). */
/* c(3) is number of neutrino species. */
/* Cosmological constant. */
/* ----------DEFAULT VARIATIONAL PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* Default initial time step. */
/* ----------VARIATIONAL PARAMETERS. */
/* Default baryon-to-photon ratio. */
/* Initial time step. */
/* ----------COMPUTATION LOCATION. */
/* Baryon-to-photon ratio. */
/* ----------RUN OPTION. */
/* Time check. */
/* Run network size. */
/* Number of nuclides in computation. */
/* ----------OUTPUT FILE STATUS. */
/* Number of reactions in computation. */
/* Number of output requests. */
/* ----------USER RESPONSE VARIABLES. */
/* Indicates if output file used. */
/* declare b */
/* Selection number. */
/* ===========================PROCEDURE DIVISION================================== */
/* set b Fierz term */
    modpr_1.b = -.2f;
/* 10--------OPEN FILES AND PRINT GREETING---------------------------------------- */
/* Fierz parameter. Extern. */
    o__1.oerr = 0;
    o__1.ounit = 2;
    o__1.ofnmlen = 10;
    o__1.ofnm = "new123.dat";
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/* Output file. */
    checkcb_1.itime = 1;
/* Time = beginning of program. */
    check_();
/* Check interface subroutine. */
    s_wsfe(&io___189);
    e_wsfe();
/* 20--------INPUT INITIALIZATION INFORMATION AND PAUSE--------------------------- */
    for (i__ = 1; i__ <= 88; ++i__) {
/* ..........READ IN REACTION PARAMETERS. */
	recpr_1.iform[i__ - 1] = (integer) recpr0_1.reacpr[i__ + 87];
/* Reaction type. */
	recpr_1.ii[i__ - 1] = (integer) recpr0_1.reacpr[i__ + 175];
/* Incoming nuclide type. */
	recpr_1.jj[i__ - 1] = (integer) recpr0_1.reacpr[i__ + 263];
/* Incoming nuclide type. */
	recpr_1.kk[i__ - 1] = (integer) recpr0_1.reacpr[i__ + 351];
/* Outgoing nuclide type. */
	recpr_1.ll[i__ - 1] = (integer) recpr0_1.reacpr[i__ + 439];
/* Outgoing nuclide type. */
	recpr_1.rev[i__ - 1] = recpr0_1.reacpr[i__ + 527];
/* Reverse reaction coefficient. */
	recpr_1.q9[i__ - 1] = recpr0_1.reacpr[i__ + 615];
/* ..........INITIALIZE REACTION RATES. */
/* Energy released. */
	rates_1.f[i__ - 1] = 0.f;
/* Forward rate coeff. */
	rates_1.r__[i__ - 1] = 0.f;
/* ..........SET RUN OPTIONS TO DEFAULT. */
/* Reverse rate coeff. */
    }
    runopt_1.irun = 1;
/* Do full run. */
    runopt_1.isize = 26;
/* Use all 26 nuclides. */
    runopt_1.jsize = 88;
/* ..........SET OUTPUT OPTION TO DEFAULT. */
/* Use all 88 reactions. */
    outopt_1.nout = 0;
/* No output requests. */
    outopt_1.outfile = FALSE_;
/* ..........SET VALUES TO DEFAULT. */
/* Output file not used. */
    compr_1.cy = compr0_1.cy0;
/* Time step limiting constant on abund */
    compr_1.ct = compr0_1.ct0;
/* Time step limiting constant on tempe */
    compr_1.t9i = compr0_1.t9i0;
/* Initial temperature. */
    compr_1.t9f = compr0_1.t9f0;
/* Final temperature. */
    compr_1.ytmin = compr0_1.ytmin0;
/* Smallest abundances allowed. */
    compr_1.inc = compr0_1.inc0;
/* Accumulation increment. */
    modpr_1.c__[0] = modpr0_1.c0[0];
/* Variation of gravitational constant. */
    modpr_1.c__[1] = modpr0_1.c0[1];
/* Neutron lifetime. */
    modpr_1.c__[2] = modpr0_1.c0[2];
/* Number of neutrino species. */
    modpr_1.cosmo = modpr0_1.cosmo0;
/* Cosmological constant. */
    modpr_1.xi[0] = modpr0_1.xi0[0];
/* Electron degeneracy parameter. */
    modpr_1.xi[1] = modpr0_1.xi0[1];
/* Muon degeneray parameter. */
    modpr_1.xi[2] = modpr0_1.xi0[2];
/* Tauon degeneracy parameter. */
    varpr_1.dt1 = varpr0_1.dt0;
/* Initial time step. */
    varpr_1.eta1 = varpr0_1.eta0;
/* ..........ACCEPT RETURN TO CONTINUE. */
/* Baryon-to-photon ratio. */
    s_rsle(&io___191);
    e_rsle();
/* 30--------PRINT MENU AND AWAIT RESPONSE---------------------------------------- */
/* ..........RETURN FROM LOOPING. */
/* Pause. */
L300:
/* ..........DISPLAY MENU. */
    s_wsfe(&io___192);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___193);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
/* 40--------BRANCH TO APPROPRIATE SECTION---------------------------------------- */
    switch (inum) {
	case 1:  goto L410;
	case 2:  goto L420;
	case 3:  goto L430;
	case 4:  goto L440;
	case 5:  goto L450;
	case 6:  goto L460;
    }
    goto L460;
/* Improper input or <RETURN>. */
L410:
/* Help section. */
    help_();
    goto L500;
L420:
/* Set computation parameters section. */
    setcom_();
    goto L500;
L430:
/* Set model parameters section. */
    setmod_();
    goto L500;
L440:
/* Run section. */
    checkcb_1.itime = 2;
/* Time = beginning of run section. */
    check_();
/* Check interface subroutine. */
    run_();
    checkcb_1.itime = 9;
/* Time = end of run section. */
    check_();
/* Check interface subroutine. */
    goto L500;
L450:
/* Output section. */
    output_();
    goto L500;
L460:
/* Exit section. */
    if (outopt_1.outfile) {
	cl__1.cerr = 0;
	cl__1.cunit = 2;
	cl__1.csta = "keep";
	f_clos(&cl__1);
/* Close output file. */
    } else {
	cl__1.cerr = 0;
	cl__1.cunit = 2;
	cl__1.csta = "delete";
	f_clos(&cl__1);
/* File not used - dispose. */
    }
    checkcb_1.itime = 10;
/* Time = end of program. */
    check_();
/* Check interface subroutine. */
    s_stop("", (ftnlen)0);
/* 50---------GO BACK TO MENU----------------------------------------------------- */
L500:
    goto L300;
} /* MAIN__ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int help_(void)
{
    /* Format strings */
    static char fmt_1000[] = "(8(/),\002 \002,32x,\002HELP SELECTION\002,/"
	    ",\002 \002,32x,\002---- ---------\002,//,\002 \002,24x,\0021. IN"
	    "TRODUCTION\002,/,\002 \002,24x,\0022. SETTING UP A RUN\002,/,"
	    "\002 \002,24x,\0023. RUNNING THE PROGRAM\002,/,\002 \002,24x,"
	    "\0024. OUTPUT OPTIONS\002,/,\002 \002,24x,\0025. GENERAL METHOD "
	    "OF COMPUTATION\002,/,\002 \002,24x,\0026. USING THE INTERFACE SU"
	    "BROUTINE\002,/,\002 \002,24x,\0027. EXIT\002,7(/),\002 \002,24x"
	    ",\002Enter selection (1-7): \002,$)";
    static char fmt_1001[] = "(i1)";
    static char fmt_2100[] = "(/,\002 \002,31x,\002INTRODUCTION\002,/,\002"
	    " \002,31x,\002------------\002,2(/),\002 \002,\002Welcome to the"
	    " wonderful world of primor\002,\002dial nucleosynthesis.  NUC123"
	    " is a      \002,/,\002 \002,\002FORTRAN program designed to prov"
	    "ide the \002,\002early universe researcher with the tools\002,/"
	    ",\002 \002,\002necessary for the investigation of primo\002,\002"
	    "rdial nucleosynthesis.  Its menu-driven \002,/,\002 \002,\002int"
	    "erface allows the user to first set c\002,\002omputation paramet"
	    "ers (such as the time \002,/,\002 \002,\002step) and model param"
	    "eters (such as the \002,\002neutron lifetime and number of neutr"
	    "i-  \002,/,\002 \002,\002nos) before doing single runs or multipl"
	    "\002,\002e runs (in which desired model parame-  \002,/,\002 "
	    "\002,\002ters are varied over a desired range.)  \002,\002After "
	    "the run, the user can utilize the \002,/,\002 \002,\002menu to e"
	    "ither produce an output file or\002,\002 to view the most recent"
	    " run on the     \002,/,\002 \002,\002screen.  The program comes "
	    "with an empty\002,\002 subroutine CHECK into which the user  "
	    " \002,/,\002 \002,\002may wish to put additional code to add "
	    "t\002,\002o the computation in an original manner.\002,10(/),"
	    "\002 \002,\002(Enter <RETURN> to go back to help menu): \002,$)";
    static char fmt_2200[] = "(/,\002 \002,29x,\002SETTING UP A RUN\002,/"
	    ",\002 \002,29x,\002------- -- - ---\002,2(/),\002 \002,\002I. Se"
	    "tting computation parameters.      \002,/,\002 \002,\002   The a"
	    "ccuracy of the computation and t\002,\002he relevant temperature"
	    " region can be   \002,/,\002 \002,\002   set by the following pa"
	    "rameters:     \002,/,\002 \002,\002    A. Time step limiting con"
	    "stant 1  (d\002,\002efault value of 0.03)                    "
	    "\002,/,\002 \002,\002    B. Time step limiting constant 2  (d"
	    "\002,\002efault value of 0.003)                   \002,/,\002"
	    " \002,\002    C. Initial time step  (default value\002,\002 of 1"
	    "0**-4)                             \002,/,\002 \002,\002    D. I"
	    "nitial temperature  (default val\002,\002ue of 10**2)           "
	    "                 \002,/,\002 \002,\002       This is the tempera"
	    "ture at the be\002,\002ginning of the run in units of 10**9 K "
	    " \002,/,\002 \002,\002    E. Final temperature  (default valu"
	    "e\002,\002 of 10**-2)                             \002,/,\002"
	    " \002,\002       This is the termination temperatu\002,\002re of"
	    " the run in units of 10**9 K       \002,/,\002 \002,\002    F. S"
	    "mallest abundances allowed  (def\002,\002ault value of 10**-25) "
	    "                 \002,/,\002 \002,\002       Elemental abundance"
	    "s are not allo\002,\002wed to drop below this value           "
	    " \002,/,\002 \002,\002    G. # of iterations for each accumul"
	    "a\002,\002tion  (default value of 300)             \002,/,\002"
	    " \002,\002       This is the number of iterations \002,\002befor"
	    "e values are put in an output array\002,6(/),\002 \002,\002(Ente"
	    "r 1 to continue, <RETURN> to end): \002,$)";
    static char fmt_2202[] = "(/,\002 \002,\002II. Setting model parameters."
	    "           \002,/,\002 \002,\002   Default values here give what"
	    " is know\002,\002n as the standard model with best guess \002,/"
	    ",\002 \002,\002   figure on the neutron lifetime of 889\002,\002"
	    ".1 seconds.  Nonstandard scenarios can\002,/,\002 \002,\002   be"
	    " investigated by varying the follow\002,\002ing parameters:     "
	    "                    \002,/,\002 \002,\002    A. The gravitationa"
	    "l constant       \002,/,\002 \002,\002       (The default value "
	    "of one here gi\002,\002ves the usual 6.6720e-8 dyne*cm**2/g**2"
	    ")\002,/,\002 \002,\002    B. Neutron life-time  (default valu"
	    "e\002,\002 of 889.1 seconds)                    \002,/,\002 \002,"
	    "\002    C. Number of neutrino species  (defa\002,\002ult value o"
	    "f 3 light neutrinos)         \002,/,\002 \002,\002    D. Final b"
	    "aryon-to-photon ratio  (se\002,\002t to log(eta) = -9.5)        "
	    "           \002,/,\002 \002,\002    E. Cosmological constant  (d"
	    "efault v\002,\002alue of 0)                              \002,/"
	    ",\002 \002,\002    F. Neutrino degeneracy parameters  (\002,\002"
	    "default values all 0)                   \002,/,\002 \002,\002   "
	    "    There are 3 separate parameters f\002,\002or the electron, m"
	    "uon, and tau neutrinos\002,11(/),\002 \002,\002(Enter <RETURN> t"
	    "o go back to help menu): \002,$)";
    static char fmt_2300[] = "(/,\002 \002,28x,\002RUNNING THE PROGRAM\002"
	    ",/,\002 \002,28x,\002------- --- -------\002,2(/),\002 \002,\002"
	    "I. Setting run speed.                   \002,/,\002 \002,\002   "
	    "The code can be run at 3 different se\002,\002ttings of speed.  "
	    "The running of the    \002,/,\002 \002,\002   code can be speede"
	    "d up by reducing th\002,\002e number of nuclides and reactions. "
	    " The\002,/,\002 \002,\002   complete computation takes into accou"
	    "\002,\002nt the following nuclides: n, p, d, t,  \002,/,\002 "
	    "\002,\002   He3, He4, Li6, Li7, Be7, Li8, B8, Be9\002,\002,B10, "
	    "B11, C11, B12, C12, N12, C13, N13,\002,/,\002 \002,\002   C14, N"
	    "14, O14, N15, O15, and O16.    \002,/,\002 \002,\002   The given"
	    " CPU percentages and abundan\002,\002ce variations are with rega"
	    "rd to a      \002,/,\002 \002,\002   single run with all default"
	    " parameter\002,\002 values.                                \002,"
	    "/,\002 \002,\002    A. 26 nuclides, 88 reactions (defaul\002,"
	    "\002t)                                      \002,/,\002 \002,"
	    "\002       nuclides from n to O16           \002,/,\002 \002,"
	    "\002    B. 18 nuclides, 60 reactions        \002,/,\002 \002,"
	    "\002       nuclides from n to N12           \002,/,\002 \002,"
	    "\002       (63% CPU time, variation = .1%)  \002,/,\002 \002,"
	    "\002    C.  9 nuclides, 25 reactions        \002,/,\002 \002,"
	    "\002       nuclides from n to Be7           \002,/,\002 \002,"
	    "\002       (20% CPU time, variation = .5%)  \002,4(/),\002 \002"
	    ",\002(Enter 1 to continue, <RETURN> to end): \002,$)";
    static char fmt_2302[] = "(/,\002 \002,\002II. Do single run.           "
	    "           \002,/,\002 \002,\002    A. Interactive.             "
	    "        \002,/,\002 \002,\002       In an interactive session, t"
	    "he us\002,\002er can readily input the computational  \002,/,"
	    "\002 \002,\002       and model parameters and begin th\002,\002e"
	    " computation process.  The run itself  \002,/,\002 \002,\002    "
	    "   is commenced when option 2, \"GO\",\002,\002 in the \"RUN\" s"
	    "ection is requested.     \002,//,\002 \002,\002    B. Batch.    "
	    "                       \002,/,\002 \002,\002       To run the pr"
	    "ogram in a batch mod\002,\002e, it must be altered slightly so t"
	    "hat  \002,/,\002 \002,\002       the I/O takes place with files "
	    "in\002,\002stead of a terminal.  This is done by   \002,/,\002"
	    " \002,\002       setting different values for the \002,\002input"
	    " and output unit number parameters \002,/,\002 \002,\002      "
	    " \"ir\" and \"iw\" and assigning them \002,\002to different file"
	    "s in NUC123.  In the   \002,/,\002 \002,\002       file assigned"
	    " the \"ir\" unit numbe\002,\002r, one must place the responses t"
	    "o the  \002,/,\002 \002,\002       queries of the program.      "
	    "    \002,10(/),\002 \002,\002(Enter 1 to continue, <RETURN> to e"
	    "nd): \002,$)";
    static char fmt_2304[] = "(/,\002 \002,\002III. Do multiple runs.       "
	    "          \002,/,\002 \002,\002   A wide range of early universe"
	    " model\002,\002s can be covered by doing many runs    \002,/,"
	    "\002 \002,\002   while one or more parameters are var\002,\002ie"
	    "d over a range of interest.  The     \002,/,\002 \002,\002   par"
	    "ameters that can be varied are th\002,\002e following:          "
	    "                 \002,/,\002 \002,\002    A. Eta                "
	    "             \002,\002       - Logrithmic variation          "
	    "\002,/,\002 \002,\002    B. Gravitational constant          \002,"
	    "\002       - Linear variation              \002,/,\002 \002,\002"
	    "    C. Neutron lifetime                \002,\002       - Linear "
	    "variation              \002,/,\002 \002,\002    D. Number of neu"
	    "trino species      \002,\002       - Linear variation           "
	    "   \002,/,\002 \002,\002    E. Cosmological constant          "
	    " \002,\002       - Linear variation              \002,/,\002 "
	    "\002,\002    F. Neutrino degeneracy parameters  \002,\002       "
	    "- Linear variation              \002,/,\002 \002,\002        1. "
	    "Electron neutrino           \002,/,\002 \002,\002        2. Muon"
	    " neutrino               \002,/,\002 \002,\002        3. Tauon ne"
	    "utrino              \002,/,\002 \002,\002   At most 3 parameters"
	    " can be varied. \002,\002 The first parameter inputted will be "
	    " \002,/,\002 \002,\002   will be varied in the outermost loop"
	    "\002,\002 and the third parameter inputted will \002,/,\002 \002,"
	    "\002   be varied in the innermost loop.    \002,7(/),\002 \002"
	    ",\002(Enter <RETURN> to go back to help menu): \002,$)";
    static char fmt_2400[] = "(/,\002 \002,30x,\002OUTPUT OPTIONS\002,/"
	    ",\002 \002,30x,\002------ -------\002,2(/),\002 \002,\002I.  Req"
	    "uest output file.                \002,/,\002 \002,\002   After a"
	    " run, the user can request the\002,\002 program to put the resul"
	    "ting numbers   \002,/,\002 \002,\002   into an output file.  Thi"
	    "s can be don\002,\002e as many times as desired and all the  "
	    "\002,/,\002 \002,\002   information will be put in one new fi"
	    "\002,\002le under the name of \"NUC123.DAT.\"  If  \002,/,\002"
	    " \002,\002   there is no request during the entire\002,\002 runn"
	    "ing of the program, this file is   \002,/,\002 \002,\002   not c"
	    "reated.  If an output file is re\002,\002quested after a multipl"
	    "e run, only the  \002,/,\002 \002,\002   information from the ve"
	    "ry last run wi\002,\002ll be given.  The output file will give"
	    " \002,/,\002 \002,\002   the computational and model paramete"
	    "r\002,\002s for each run and will contain the     \002,/,\002"
	    " \002,\002   following information:               \002,/,\002"
	    " \002,\002    A. Temperatures in decreasing order \002,/,\002"
	    " \002,\002    B. Abundances for n, p, d, t, He3, H\002,\002e4, L"
	    "i6, Li7, Be7, and Li8 & up         \002,/,\002 \002,\002       ("
	    "p and He4 are in mass fraction, \002,\002the rest in ratios to t"
	    "he p abundance)  \002,/,\002 \002,\002    C. Time, time interval"
	    ", chemical pot\002,\002ential of the electron                 "
	    " \002,/,\002 \002,\002    D. Energy densities for photons, el"
	    "e\002,\002ctrons, electron neutrinos, and baryons \002,/,\002"
	    " \002,\002    E. Baryon-to-photon ratio, expansion\002,\002 rate"
	    " of the universe                   \002,5(/),\002 \002,\002(Ente"
	    "r 1 to continue, <RETURN> to end): \002,$)";
    static char fmt_2402[] = "(/,\002 \002,\002II.  Request output on screen"
	    ".         \002,/,\002 \002,\002   Instead of waiting to print ou"
	    "t an o\002,\002utput file, the user can immediately   \002,/,"
	    "\002 \002,\002   access the results of the latest run\002,\002 b"
	    "y requesting the output on the       \002,/,\002 \002,\002   scr"
	    "een.  There are four screens on e\002,\002ach of which are displ"
	    "ayed the         \002,/,\002 \002,\002   computational and model"
	    " parameters a\002,\002nd the temperature:                    "
	    "\002,/,\002 \002,\002    A. Abundances for d, t, He3, He4, a\002,"
	    "\002nd Li7                                 \002,/,\002 \002,\002"
	    "       (He4 in mass fraction, rest as a\002,\002 ratio with the "
	    "p abundance)           \002,/,\002 \002,\002    B. Abundances fo"
	    "r n, p, Li6, Be7, a\002,\002nd Li8 & up                         "
	    "   \002,/,\002 \002,\002       (p in mass fraction, rest as a "
	    "r\002,\002atio with the p abundance)             \002,/,\002 "
	    "\002,\002    C. Energy densities for photons, el\002,\002ectrons"
	    ", electron neutrinos, & baryons \002,/,\002 \002,\002    D. Time"
	    ", time interval, chemical po\002,\002tential of the electron,   "
	    "            \002,/,\002 \002,\002       baryon-to-photon ratio, "
	    "and expa\002,\002nsion rate of the universe             \002,11("
	    "/),\002 \002,\002(Enter <RETURN> to go back to help menu): \002,"
	    "$)";
    static char fmt_2500[] = "(/,\002 \002,22x,\002GENERAL METHOD OF COMPUTA"
	    "TION\002,/,\002 \002,22x,\002------- ------ -- -----------\002,2"
	    "(/),\002 \002,\002I. Time evolution algorithm.            \002,/,"
	    "\002 \002,\002   The program utilizes a 2-point Runge-\002,\002K"
	    "utta scheme (located in subroutine     \002,/,\002 \002,\002   D"
	    "RIVER) to time-evolve the temperatur\002,\002e, the quantity hv "
	    "(the ratio of the    \002,/,\002 \002,\002   baryon density to T"
	    "**3), the chemical\002,\002 potential of the electron, and the  "
	    "   \002,/,\002 \002,\002   nuclide abundances.  In the 2-point "
	    "R\002,\002unge-Kutta routine, a variable v at time\002,/,\002"
	    " \002,\002   t0 (= v0) is evolved to a time t1 by \002,\002addin"
	    "g to v0 the average of the         \002,/,\002 \002,\002   deriv"
	    "atives evaluated at t0 and at t1\002,\002 multiplied by dt:     "
	    "                 \002,/,\002 \002,\002       v1 = v0 + 0.5(dvdt("
	    "t0)+dvdt(t1)) \002,/,\002 \002,\002   where dvdt(t1) is gotten b"
	    "y first fin\002,\002ding v1' = v0 + dvdt(t0).  The         \002,"
	    "/,\002 \002,\002   derivatives of the nuclide abundances\002,"
	    "\002 are first computed and these are used  \002,/,\002 \002,"
	    "\002   to find the derivatives of t9, hv, an\002,\002d phie (thi"
	    "s is done in subroutine      \002,/,\002 \002,\002   DERIVS).  T"
	    "o compute the time derivat\002,\002ives of the nuclide abundance"
	    "s, a matrix\002,/,\002 \002,\002   equation is set up (in subrou"
	    "tine SOL\002,\002) and is solved (in subroutine EQSLIN)  \002,/"
	    ",\002 \002,\002   by gaussian elimination utilizing imp\002,\002"
	    "licit differentiation.                  \002,6(/),\002 \002,\002"
	    "(Enter 1 to continue, <RETURN> to end): \002,$)";
    static char fmt_2502[] = "(/\002 \002,\002II. Hierarchy of Subroutines. "
	    "  \002,/,\002 \002,\002    NUC123                       \002,"
	    "\002     Main program (main menu)    \002,/,\002 \002,\002      "
	    "  HELP                     \002,\002     Help option            "
	    "     \002,/,\002 \002,\002        SETCOM                   \002"
	    ",\002     Set computational parameters\002,/,\002 \002,\002     "
	    "   SETMOD                   \002,\002     Set model parameters  "
	    "      \002,/,\002 \002,\002        RUN                      \002,"
	    "\002     Run computation code        \002,/,\002 \002,\002      "
	    "      DRIVER               \002,\002     Main routine (Runge-Kut"
	    "ta loop)    \002,/,\002 \002,\002                START          "
	    "  \002,\002     Initialization routine      \002,/,\002 \002,"
	    "\002                    RATE0        \002,\002     Computes weak"
	    " decay rates   \002,/,\002 \002,\002                DERIVS      "
	    "     \002,\002     Computes time derivatives   \002,/,\002 \002"
	    ",\002                    THERM        \002,\002     Computes ene"
	    "rgy densities   \002,/,\002 \002,\002                        BES"
	    "SEL   \002,\002     Gives functions of Kn       \002,/,\002 \002,"
	    "\002                            KNUX \002,\002     Computes modi"
	    "fied Bessel fcn Kn    \002,/,\002 \002,\002                     "
	    "   NUDENS   \002,\002     Computes neutrino energy density   "
	    "\002,/,\002 \002,\002                    RATE1-4      \002,\002 "
	    "    Computes rates for reactions\002,/,\002 \002,\002           "
	    "         SOL          \002,\002     Builds A matrix for eqn dy/d"
	    "t = Ay \002,/,\002 \002,\002                        EQSLIN   "
	    "\002,\002     Solves dy/dt=Ay by gaussian elim   \002,/,\002 "
	    "\002,\002                ACCUM            \002,\002     Output a"
	    "ccumulator          \002,/,\002 \002,\002        OUTPUT         "
	    "          \002,\002     Allows user to output result\002,4(/)"
	    ",\002 \002,\002(Enter <RETURN> to go back to help menu): \002,$)";
    static char fmt_2600[] = "(/,\002 \002,22x,\002USING THE INTERFACE SUBRO"
	    "UTINE\002,/,\002 \002,22x,\002----- --- --------- ----------\002"
	    ",2(/),\002 \002,\002I. Purpose.                             \002"
	    ",/,\002 \002,\002   The interface subroutine CHECK is des\002"
	    ",\002igned to be an outlet of the program    \002,/,\002 \002"
	    ",\002   into which alterations can be easily \002,\002plugged.  "
	    "Programs are normally modified\002,/,\002 \002,\002   by searchi"
	    "ng through the program, ide\002,\002ntifying the appropriate are"
	    "as for      \002,/,\002 \002,\002   alterations, and interspersi"
	    "ng new co\002,\002mmands while deleting some old ones.    \002,/,"
	    "\002 \002,\002   This process can get tricky unless on\002,\002e"
	    " actively documents the alterations:   \002,/,\002 \002,\002   o"
	    "ne might lose track of all of the mo\002,\002difications and del"
	    "etions.  Thus, it is \002,/,\002 \002,\002   worthwhile to put m"
	    "ost if not all of \002,\002the necessary changes into one       "
	    "   \002,/,\002 \002,\002   subroutine which is to be called fro"
	    "m\002,\002 strategic locations in the main        \002,/,\002"
	    " \002,\002   program.  Furthermore, by putting cha\002,\002nges "
	    "into one small subroutine, one need\002,/,\002 \002,\002   only "
	    "to compile the subroutine CHECK \002,\002each time instead of th"
	    "e entire nucleo- \002,/,\002 \002,\002   synthesis code.        "
	    "              \002,8(/),\002 \002,\002(Enter 1 to continue, <RET"
	    "URN> to end): \002,$)";
    static char fmt_2602[] = "(/,\002 \002,\002II. Description.             "
	    "           \002,/,\002 \002,\002   Subroutine CHECK is an empty "
	    "subrouti\002,\002ne with a large COMMON area, giving the \002,/"
	    ",\002 \002,\002   user ready access to all of the impor\002,\002"
	    "tant variables in the computations.  The\002,/,\002 \002,\002   "
	    "routine is called from various locati\002,\002ons in the main pr"
	    "ogram and the location\002,/,\002 \002,\002   spot in the progra"
	    "m is labeled by the flag \"itime\".  The set call locations  "
	    "\002,/,\002 \002,\002   are given below:                     "
	    "\002,/,\002 \002,\002    A. itime = 1 (NUC123, very beginning"
	    "\002,\002 of program run)                        \002,/,\002 "
	    "\002,\002       (appropriate for opening files, i\002,\002nitial"
	    "izing variables)                  \002,/,\002 \002,\002    B. it"
	    "ime = 2 (NUC123, right before g\002,\002oing into the RUN sectio"
	    "n)              \002,/,\002 \002,\002    C. itime = 3 (RUN, righ"
	    "t before goin\002,\002g into DRIVER to do the computations)  "
	    " \002,/,\002 \002,\002    D. itime = 4 (DRIVER, in 1st R-K lo"
	    "o\002,\002p after computing derivatives in DERIVS)\002,/,\002"
	    " \002,\002    E. itime = 7 (DRIVER, in 2nd R-K loo\002,\002p aft"
	    "er computing derivatives in DERIVS)\002,/,\002 \002,\002    F. i"
	    "time = 8 (RUN, right after comin\002,\002g back from DRIVER)    "
	    "                 \002,/,\002 \002,\002    G. itime = 9 (NUC123, "
	    "right after co\002,\002ming back from the RUN section)        "
	    " \002,/,\002 \002,\002    H. itime =10 (NUC123, very end of p"
	    "r\002,\002ogram run)                              \002,/,\002"
	    " \002,\002       (appropriate for closing files)  \002,/,\002"
	    " \002,\002   The difference between the (2,9) pair\002,\002ing a"
	    "nd the (3,8) pairing is that for a \002,/,\002 \002,\002   multi"
	    "ple run, the (3,8) pairing would\002,\002 be called before and a"
	    "fter every run   \002,/,\002 \002,\002   but the (2,9) pairing w"
	    "ould be called\002,\002 before and after the entire sequence. "
	    " \002,4(/),\002 \002,\002(Enter 1 to continue, <RETURN> to end): "
	    "\002,$)";
    static char fmt_2604[] = "(/,\002 \002,\002III. Implementation.         "
	    "          \002,/,\002 \002,\002   The additional program stateme"
	    "nts ar\002,\002e needed in the subroutine CHECK.  If a\002,/,"
	    "\002 \002,\002   particular command is to be executed\002,\002 w"
	    "hen the computer is at a certain     \002,/,\002 \002,\002   loc"
	    "ation in the program -- say label\002,\002ed by itime = 8 -- the"
	    "n in CHECK, one  \002,/,\002 \002,\002   must place the command "
	    "under the sta\002,\002tement, IF (itime.eq.8)....  The user  "
	    "\002,/,\002 \002,\002   is at leisure to place his own locat\002,"
	    "\002ion indicators (5,6) and CALL CHECK    \002,/,\002 \002,\002"
	    "   statements anywhere in the program a\002,\002s long as there "
	    "is a COMMON /checkcb/    \002,/,\002 \002,\002   statement in th"
	    "e particular subrouti\002,\002ne to carry the value of itime alo"
	    "ng.  \002,15(/),\002 \002,\002(Enter <RETURN> to go back to help"
	    " menu): \002,$)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsle(cilist *), e_rsle(void);

    /* Local variables */
    static integer inum;

    /* Fortran I/O blocks */
    static cilist io___195 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___196 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___198 = { 0, 6, 0, fmt_2100, 0 };
    static cilist io___199 = { 0, 5, 0, 0, 0 };
    static cilist io___200 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___201 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___202 = { 0, 6, 0, fmt_2202, 0 };
    static cilist io___203 = { 0, 5, 0, 0, 0 };
    static cilist io___204 = { 0, 6, 0, fmt_2300, 0 };
    static cilist io___205 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___206 = { 0, 6, 0, fmt_2302, 0 };
    static cilist io___207 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___208 = { 0, 6, 0, fmt_2304, 0 };
    static cilist io___209 = { 0, 5, 0, 0, 0 };
    static cilist io___210 = { 0, 6, 0, fmt_2400, 0 };
    static cilist io___211 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___212 = { 0, 6, 0, fmt_2402, 0 };
    static cilist io___213 = { 0, 5, 0, 0, 0 };
    static cilist io___214 = { 0, 6, 0, fmt_2500, 0 };
    static cilist io___215 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___216 = { 0, 6, 0, fmt_2502, 0 };
    static cilist io___217 = { 0, 5, 0, 0, 0 };
    static cilist io___218 = { 0, 6, 0, fmt_2600, 0 };
    static cilist io___219 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___220 = { 0, 6, 0, fmt_2602, 0 };
    static cilist io___221 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___222 = { 0, 6, 0, fmt_2604, 0 };
    static cilist io___223 = { 0, 5, 0, 0, 0 };


/* ----------LINKAGES. */
/*     CALLED BY - [program] nuc123 */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Displays description and workings of the program. */
/* ----------PARAMETERS. */
/* Input unit number (previous value = */
/* ==========================DECLARATION DIVISION================================= */
/* ----------USER RESPONSE VARIABLES. */
/* Output unit number (previous value = */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------PRINT HELP SELECTION------------------------------------------------- */
/* ..........RETURN FROM LOOPING. */
/* Selection number. */
L100:
/* ..........DISPLAY MENU. */
    s_wsfe(&io___195);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___196);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
/* 20--------BRANCH TO APPROPRIATE SECTION---------------------------------------- */
    switch (inum) {
	case 1:  goto L210;
	case 2:  goto L220;
	case 3:  goto L230;
	case 4:  goto L240;
	case 5:  goto L250;
	case 6:  goto L260;
	case 7:  goto L270;
    }
    goto L270;
/* 21--------INTRODUCTION SECTION------------------------------------------------- */
/* Improper input or <RETURN>. */
L210:
/* Setting up a run section. */
    s_wsfe(&io___198);
    e_wsfe();
    s_rsle(&io___199);
    e_rsle();
    goto L300;
/* 22--------SET UP RUN SECTION--------------------------------------------------- */
L220:
/* Setting up a run section. */
    s_wsfe(&io___200);
    e_wsfe();
    s_rsfe(&io___201);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___202);
	e_wsfe();
	s_rsle(&io___203);
	e_rsle();
	goto L300;
    } else {
	goto L300;
    }
/* 23--------RUN PROGRAM SECTION-------------------------------------------------- */
/* (inum.eq.1) */
L230:
/* Running the program section. */
    s_wsfe(&io___204);
    e_wsfe();
    s_rsfe(&io___205);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___206);
	e_wsfe();
	s_rsfe(&io___207);
	do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
	e_rsfe();
	if (inum == 1) {
	    s_wsfe(&io___208);
	    e_wsfe();
	    s_rsle(&io___209);
	    e_rsle();
	    goto L300;
	} else {
	    goto L300;
	}
/* (inum.eq.1) */
    } else {
	goto L300;
    }
/* 24--------OUTPUT OPTIONS SECTION----------------------------------------------- */
/* (inum.eq.1) */
L240:
/* Output options section. */
    s_wsfe(&io___210);
    e_wsfe();
    s_rsfe(&io___211);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___212);
	e_wsfe();
	s_rsle(&io___213);
	e_rsle();
	goto L300;
    } else {
	goto L300;
    }
/* 25--------METHOD OF COMPUTATION SECTION---------------------------------------- */
/* (inum.eq.1) */
L250:
/* General method of computation sectio */
    s_wsfe(&io___214);
    e_wsfe();
    s_rsfe(&io___215);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___216);
	e_wsfe();
	s_rsle(&io___217);
	e_rsle();
	goto L300;
    } else {
	goto L300;
    }
/* 26--------USING INTERFACE SUBROUTINE SECTION. */
/* (inum.eq.1) */
L260:
/* Using the interface subroutine secti */
    s_wsfe(&io___218);
    e_wsfe();
    s_rsfe(&io___219);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___220);
	e_wsfe();
	s_rsfe(&io___221);
	do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
	e_rsfe();
	if (inum == 1) {
	    s_wsfe(&io___222);
	    e_wsfe();
	    s_rsle(&io___223);
	    e_rsle();
	    goto L300;
	} else {
	    goto L300;
	}
/* (inum.eq.1) */
    } else {
	goto L300;
    }
/* 27--------EXIT SECTION--------------------------------------------------------- */
/* (inum.eq.1) */
L270:
/* Exit section. */
    return 0;
/* 30--------GO BACK TO MAIN MENU------------------------------------------------- */
L300:
    goto L100;
} /* help_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int setcom_(void)
{
    /* Format strings */
    static char fmt_1000[] = "(8(/),\002 \002,21x,\002SET COMPUTATION PARAME"
	    "TERS SELECTION\002,/,\002 \002,21x,\002--- ----------- ---------"
	    "- ---------\002,//,\002 \002,10x,\002 1. CHANGE TIME-STEP LIMITI"
	    "NG CONSTANT 1  FROM \002,f5.3,/,\002 \002,10x,\002 2. CHANGE TIM"
	    "E-STEP LIMITING CONSTANT 2  FROM \002,f5.3,/,\002 \002,10x,\002 "
	    "3. CHANGE INITIAL TIME-STEP              FROM \002,1pe8.2,\002 S"
	    "ECONDS\002,/,\002 \002,10x,\002 4. CHANGE INITIAL TEMPERATURE   "
	    "         FROM \002,1pe8.2,\002 (10**9 K)\002,/,\002 \002,10x,"
	    "\002 5. CHANGE FINAL TEMPERATURE              FROM \002,1pe8.2"
	    ",\002 (10**9 K)\002,/,\002 \002,10x,\002 6. CHANGE SMALLEST ABUN"
	    "DANCES ALLOWED    FROM \002,1pe8.2,/,\002 \002,10x,\002 7. CHANG"
	    "E ACCUMULATION INCREMENT         FROM \002,1pe8.2,\002 ITERATIONS"
	    "\002,/,\002 \002,10x,\002 8. RESET ALL TO DEFAULT VALUES\002,/"
	    ",\002 \002,10x,\002 9. EXIT\002,5(/),\002 \002,10x,\002Enter sel"
	    "ection (1-9): \002,$)";
    static char fmt_1001[] = "(i1)";
    static char fmt_2100[] = "(\002 \002,\002Enter value for time step limit"
	    "ing constant 1: \002,$)";
    static char fmt_2200[] = "(\002 \002,\002Enter value for time step limit"
	    "ing constant 2: \002,$)";
    static char fmt_2300[] = "(\002 \002,\002Enter value for initial time st"
	    "ep: \002,$)";
    static char fmt_2400[] = "(\002 \002,\002Enter value for initial tempera"
	    "ture: \002,$)";
    static char fmt_2500[] = "(\002 \002,\002Enter value for final temperatu"
	    "re: \002,$)";
    static char fmt_2600[] = "(\002 \002,\002Enter value for smallest abunda"
	    "nces allowed: \002,$)";
    static char fmt_2700[] = "(\002 \002,\002Enter value for accumulation in"
	    "crement: \002,$)";
    static char fmt_2800[] = "(\002 \002,\002All values reset to default - P"
	    "ress <RETURN> to continue: \002,$)";

    /* System generated locals */
    real r__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsfe(cilist *), e_rsfe(void), s_rsle(cilist *), do_lio(integer 
	    *, integer *, char *, ftnlen), e_rsle(void);

    /* Local variables */
    static integer inum;

    /* Fortran I/O blocks */
    static cilist io___224 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___225 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___227 = { 0, 6, 0, fmt_2100, 0 };
    static cilist io___228 = { 0, 5, 0, 0, 0 };
    static cilist io___229 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___230 = { 0, 5, 0, 0, 0 };
    static cilist io___231 = { 0, 6, 0, fmt_2300, 0 };
    static cilist io___232 = { 0, 5, 0, 0, 0 };
    static cilist io___233 = { 0, 6, 0, fmt_2400, 0 };
    static cilist io___234 = { 0, 5, 0, 0, 0 };
    static cilist io___235 = { 0, 6, 0, fmt_2500, 0 };
    static cilist io___236 = { 0, 5, 0, 0, 0 };
    static cilist io___237 = { 0, 6, 0, fmt_2600, 0 };
    static cilist io___238 = { 0, 5, 0, 0, 0 };
    static cilist io___239 = { 0, 6, 0, fmt_2700, 0 };
    static cilist io___240 = { 0, 5, 0, 0, 0 };
    static cilist io___241 = { 0, 6, 0, fmt_2800, 0 };
    static cilist io___242 = { 0, 5, 0, 0, 0 };


/* ----------LINKAGES. */
/*     CALLED BY - [program] nuc123 */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Allows resetting of computation parameters. */
/* ----------PARAMETERS. */
/* Input unit number (previous value = */
/* ----------COMMON AREAS. */
/* Output unit number (previous value = */
/* Default comp param */
/* Computation parame */
/* Default variationl */
/* ==========================DECLARATION DIVISION================================= */
/* ----------DEFAULT COMPUTATION PARAMETERS. */
/* Variational parame */
/* Default cy. */
/* Default ct. */
/* Default t9i. */
/* Default t9f. */
/* Default ytmin. */
/* ----------COMPUTATION PARAMETERS. */
/* Default accumulation increment. */
/* Time step limiting constant on abund */
/* Time step limiting constant on tempe */
/* Initial temperature (in 10**9 K). */
/* Final temperature (in 10**9 K). */
/* Smallest abundances allowed. */
/* ----------DEFAULT VARIATIONAL  PARAMETERS. */
/* Accumulation increment. */
/* ----------VARIATIONAL  PARAMETERS. */
/* Default initial dt. */
/* ----------LOCAL VARIABLES. */
/* Initial time step. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------PRINT RESET SELECTION AND AWAIT RESPONSE----------------------------- */
/* ..........RETURN FROM LOOPING. */
/* Selection number. */
L100:
/* ..........DISPLAY RESET SELECTIONS. */
    s_wsfe(&io___224);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&varpr_1.dt1, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    r__1 = (real) compr_1.inc;
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___225);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
/* 20--------BRANCH TO APPROPRIATE SECTION---------------------------------------- */
    switch (inum) {
	case 1:  goto L210;
	case 2:  goto L220;
	case 3:  goto L230;
	case 4:  goto L240;
	case 5:  goto L250;
	case 6:  goto L260;
	case 7:  goto L270;
	case 8:  goto L280;
	case 9:  goto L300;
    }
    goto L300;
/* Improper input or <RETURN>. */
L210:
/* Change time step limiting const 1 se */
    s_wsfe(&io___227);
    e_wsfe();
    s_rsle(&io___228);
    do_lio(&c__4, &c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    e_rsle();
/* L2101: */
    goto L400;
L220:
/* Change time step limiting const 2 se */
    s_wsfe(&io___229);
    e_wsfe();
    s_rsle(&io___230);
    do_lio(&c__4, &c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L230:
/* Change initial time step section. */
    s_wsfe(&io___231);
    e_wsfe();
    s_rsle(&io___232);
    do_lio(&c__4, &c__1, (char *)&varpr_1.dt1, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L240:
/* Change initial temperature section. */
    s_wsfe(&io___233);
    e_wsfe();
    s_rsle(&io___234);
    do_lio(&c__4, &c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L250:
/* Change final temperature section. */
    s_wsfe(&io___235);
    e_wsfe();
    s_rsle(&io___236);
    do_lio(&c__4, &c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L260:
/* Change smallest abundances allowed s */
    s_wsfe(&io___237);
    e_wsfe();
    s_rsle(&io___238);
    do_lio(&c__4, &c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L270:
/* Change accumulation increment sectio */
    s_wsfe(&io___239);
    e_wsfe();
    s_rsle(&io___240);
    do_lio(&c__3, &c__1, (char *)&compr_1.inc, (ftnlen)sizeof(integer));
    e_rsle();
    goto L400;
L280:
/* Reset all to default values section. */
    compr_1.cy = compr0_1.cy0;
/* Time step limiting constant on abund */
    compr_1.ct = compr0_1.ct0;
/* Time step limiting constant on tempe */
    varpr_1.dt1 = varpr0_1.dt0;
/* Time step. */
    compr_1.t9i = compr0_1.t9i0;
/* Initial temperature. */
    compr_1.t9f = compr0_1.t9f0;
/* Final temperature. */
    compr_1.ytmin = compr0_1.ytmin0;
/* Smallest abundances allowed. */
    compr_1.inc = compr0_1.inc0;
/* Accumulation increment. */
    s_wsfe(&io___241);
    e_wsfe();
    s_rsle(&io___242);
    e_rsle();
    goto L400;
L300:
/* Exit section. */
    return 0;
/* 40--------GO BACK TO MENU------------------------------------------------------ */
L400:
    goto L100;
} /* setcom_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int setmod_(void)
{
    /* Format strings */
    static char fmt_1000[] = "(8(/),\002 \002,24x,\002SET MODEL PARAMETERS S"
	    "ELECTION\002,/,\002 \002,24x,\002--- ----- ---------- --------"
	    "-\002,//,\002 \002,10x,\002 1. CHANGE GRAVITATIONAL CONSTANT    "
	    "     FROM \002,1pe10.3,/,\002 \002,10x,\002 2. CHANGE NEUTRON LI"
	    "FETIME               FROM \002,1pe10.3,\002 SECONDS\002,/,\002"
	    " \002,10x,\002 3. CHANGE NUMBER OF NEUTRINO SPECIES     FROM "
	    "\002,1pe10.3,/,\002 \002,10x,\002 4. CHANGE FINAL BARYON-TO-PHOT"
	    "ON RATIO   FROM \002,1pe10.3,/,\002 \002,10x,\002 5. CHANGE COSM"
	    "OLOGICAL CONSTANT          FROM \002,1pe10.3,/,\002 \002,10x,"
	    "\002 6. CHANGE XI-ELECTRON                    FROM \002,1pe10.3,"
	    "/,\002 \002,10x,\002 7. CHANGE XI-MUON                        FR"
	    "OM \002,1pe10.3,/,\002 \002,10x,\002 8. CHANGE XI-TAUON         "
	    "              FROM \002,1pe10.3,/,\002 \002,10x,\002 9. RESET AL"
	    "L TO DEFAULT VALUES\002,/,\002 \002,10x,\00210. EXIT\002,4(/)"
	    ",\002 \002,10x,\002 Enter selection (1-10): \002,$)";
    static char fmt_1001[] = "(i2)";
    static char fmt_2100[] = "(\002 \002,\002Enter value for variation of gr"
	    "avitational \002,\002constant: \002,$)";
    static char fmt_2200[] = "(\002 \002,\002Enter value for neutron lifetim"
	    "e (sec): \002,$)";
    static char fmt_2300[] = "(\002 \002,\002Enter value for number of neutr"
	    "ino species: \002,$)";
    static char fmt_2400[] = "(\002 \002,\002Enter value for baryon-to-photo"
	    "n ratio: \002,$)";
    static char fmt_2500[] = "(\002 \002,\002Enter value for cosmological co"
	    "nstant: \002,$)";
    static char fmt_2600[] = "(\002 \002,\002Enter value for xi electron:"
	    " \002,$)";
    static char fmt_2700[] = "(\002 \002,\002Enter value for xi muon: \002,$)"
	    ;
    static char fmt_2800[] = "(\002 \002,\002Enter value for xi tauon: \002,"
	    "$)";
    static char fmt_2802[] = "(\002 \002,\002Number of neutrinos set to 3"
	    "\002)";
    static char fmt_2804[] = "(\002 \002,\002Press <RETURN> to continue: "
	    "\002,$)";
    static char fmt_2900[] = "(\002 \002,\002All values reset to default - P"
	    "ress <RETURN> to continue: \002,$)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsfe(cilist *), e_rsfe(void), s_rsle(cilist *), do_lio(integer 
	    *, integer *, char *, ftnlen), e_rsle(void);

    /* Local variables */
    static integer inum;

    /* Fortran I/O blocks */
    static cilist io___243 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___244 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___246 = { 0, 6, 0, fmt_2100, 0 };
    static cilist io___247 = { 0, 5, 0, 0, 0 };
    static cilist io___248 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___249 = { 0, 5, 0, 0, 0 };
    static cilist io___250 = { 0, 6, 0, fmt_2300, 0 };
    static cilist io___251 = { 0, 5, 0, 0, 0 };
    static cilist io___252 = { 0, 6, 0, fmt_2400, 0 };
    static cilist io___253 = { 0, 5, 0, 0, 0 };
    static cilist io___254 = { 0, 6, 0, fmt_2500, 0 };
    static cilist io___255 = { 0, 5, 0, 0, 0 };
    static cilist io___256 = { 0, 6, 0, fmt_2600, 0 };
    static cilist io___257 = { 0, 5, 0, 0, 0 };
    static cilist io___258 = { 0, 6, 0, fmt_2700, 0 };
    static cilist io___259 = { 0, 5, 0, 0, 0 };
    static cilist io___260 = { 0, 6, 0, fmt_2800, 0 };
    static cilist io___261 = { 0, 5, 0, 0, 0 };
    static cilist io___262 = { 0, 6, 0, fmt_2802, 0 };
    static cilist io___263 = { 0, 6, 0, fmt_2804, 0 };
    static cilist io___264 = { 0, 5, 0, 0, 0 };
    static cilist io___265 = { 0, 6, 0, fmt_2900, 0 };
    static cilist io___266 = { 0, 5, 0, 0, 0 };


/* ----------LINKAGES. */
/*     CALLED BY - [program] nuc123 */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Allows resetting of model parameters. */
/* ----------PARAMETERS. */
/* Input unit number (previous value = */
/* ----------COMMON AREAS. */
/* Output unit number (previous value = */
/* Default model para */
/* Model parameters. */
/* Default variationl */
/* ==========================DECLARATION DIVISION================================= */
/* ----------DEFAULT MODEL PARAMETERS. */
/* Variational parame */
/* Default c. */
/* Default cosmological constant. */
/* ----------EARLY UNIVERSE MODEL PARAMETERS. */
/* Default neutrino degeneracy paramete */
/* c(1) is variation of gravitational c */
/* c(2) is neutron lifetime (sec). */
/* c(3) is number of neutrino species. */
/* Cosmological constant. */
/* ----------DEFAULT VARIATIONAL PARAMETERS. */
/* Neutrino degeneracy parameters. */
/* ----------VARIATIONAL PARAMETERS. */
/* Default eta. */
/* ----------USER RESPONSE VARIABLES. */
/* Intial baryon-to-photon ratio. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------PRINT RESET SELECTION AND AWAIT RESPONSE----------------------------- */
/* ..........RETURN FROM LOOPING. */
/* Selection number. */
L100:
/* ..........DISPLAY RESET SELECTIONS. */
    s_wsfe(&io___243);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&varpr_1.eta1, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___244);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
/* 20--------BRANCH TO APPROPRIATE SECTION---------------------------------------- */
    switch (inum) {
	case 1:  goto L210;
	case 2:  goto L220;
	case 3:  goto L230;
	case 4:  goto L240;
	case 5:  goto L250;
	case 6:  goto L260;
	case 7:  goto L270;
	case 8:  goto L280;
	case 9:  goto L290;
	case 10:  goto L300;
    }
    goto L300;
/* Improper input or <RETURN>. */
L210:
/* Change gravitational constant sectio */
    s_wsfe(&io___246);
    e_wsfe();
    s_rsle(&io___247);
    do_lio(&c__4, &c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L220:
/* Change neutron lifetime section. */
    s_wsfe(&io___248);
    e_wsfe();
    s_rsle(&io___249);
    do_lio(&c__4, &c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L230:
/* Change number of neutrino species se */
    s_wsfe(&io___250);
    e_wsfe();
    s_rsle(&io___251);
    do_lio(&c__4, &c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L240:
/* Change baryon-to-photon ratio sectio */
    s_wsfe(&io___252);
    e_wsfe();
    s_rsle(&io___253);
    do_lio(&c__4, &c__1, (char *)&varpr_1.eta1, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L250:
/* Change cosmological constant section */
    s_wsfe(&io___254);
    e_wsfe();
    s_rsle(&io___255);
    do_lio(&c__4, &c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L260:
/* Change neutrino degeneracy section. */
    s_wsfe(&io___256);
    e_wsfe();
    s_rsle(&io___257);
    do_lio(&c__4, &c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L270:
/* Change neutrino degeneracy section. */
    s_wsfe(&io___258);
    e_wsfe();
    s_rsle(&io___259);
    do_lio(&c__4, &c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L280:
/* Change neutrino degeneracy section. */
    s_wsfe(&io___260);
    e_wsfe();
    s_rsle(&io___261);
    do_lio(&c__4, &c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_rsle();
    if (modpr_1.xi[2] != 0.f && modpr_1.c__[2] < 3.f) {
	modpr_1.c__[2] = 3.f;
	s_wsfe(&io___262);
	e_wsfe();
	s_wsfe(&io___263);
	e_wsfe();
	s_rsle(&io___264);
	e_rsle();
    }
    goto L400;
L290:
/* Reset all to default values section. */
    modpr_1.c__[0] = modpr0_1.c0[0];
    modpr_1.c__[1] = modpr0_1.c0[1];
    modpr_1.c__[2] = modpr0_1.c0[2];
    modpr_1.cosmo = modpr0_1.cosmo0;
    modpr_1.xi[0] = modpr0_1.xi0[0];
    modpr_1.xi[1] = modpr0_1.xi0[1];
    modpr_1.xi[2] = modpr0_1.xi0[2];
    varpr_1.eta1 = varpr0_1.eta0;
    s_wsfe(&io___265);
    e_wsfe();
    s_rsle(&io___266);
    e_rsle();
    goto L400;
L300:
/* Exit section. */
    return 0;
/* 40--------GO BACK TO MENU------------------------------------------------------ */
L400:
    goto L100;
} /* setmod_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int run_(void)
{
    /* Initialized data */

    static char vtype[22*8] = "baryon/photon ratio   " "gravitational consta"
	    "nt" "neutron lifetime      " "# of neutrino species " "cosmologi"
	    "cal constant " "xi-electron           " "xi-muon               " 
	    "xi-tauon              ";

    /* Format strings */
    static char fmt_1000[] = "(8(/),\002 \002,32x,\002RUN SELECTION\002,/"
	    ",\002 \002,32x,\002--- ---------\002,//,\002 \002,27x,\002 1. SE"
	    "T RUN NETWORK\002,/,\002 \002,27x,\002 2. GO\002,/,\002 \002,27x,"
	    "\002 3. DO MULTIPLE RUNS\002,/,\002 \002,27x,\002 4. EXIT\002,10"
	    "(/),\002 \002,27x,\002 Enter selection (1-4): \002,$)";
    static char fmt_1001[] = "(i1)";
    static char fmt_2100[] = "(\002 \002,\002Enter network size (1-26 nuclid"
	    "es (default); \002,\0022-18; 3-9): \002,$)";
    static char fmt_2104[] = "(\002 \002,\002Run network set to \002,i1,\002"
	    " - Press <RETURN> \002,\002to continue: \002,$)";
    static char fmt_2200[] = "(\002 \002,\002Begin computation run....\002)";
    static char fmt_2202[] = "(\002 \002,\002Computation completed - Press <"
	    "RETURN> to continue: \002,$)";
    static char fmt_2300[] = "(\002 \002,\002Enter the number of loopings to"
	    " be done (1 \002,\002(default); 2; 3): \002,$)";
    static char fmt_2302[] = "(8(/),\002 \002,30x,\002QUANTITY TO VARY\002"
	    ",/,\002 \002,30x,\002-------- -- ----\002,//,\002 \002,25x,\002 "
	    "1.  ETA     (LOGRITHMIC VARIATION)\002,/,\002 \002,25x,\002 2.  "
	    "G           (LINEAR VARIATION)\002,/,\002 \002,25x,\002 3.  TAU "
	    "        (LINEAR VARIATION)\002,/,\002 \002,25x,\002 4.  # NEUTRI"
	    "NOS (LINEAR VARIATION)\002,/,\002 \002,25x,\002 5.  LAMBDA      "
	    "(LINEAR VARIATION)\002,/,\002 \002,25x,\002 6.  XI-ELECTRON (LIN"
	    "EAR VARIATION)\002,/,\002 \002,25x,\002 7.  XI-MUON     (LINEAR "
	    "VARIATION)\002,/,\002 \002,25x,\002 8.  XI-TAUON    (LINEAR VARI"
	    "ATION)\002,/,\002 \002,25x,\002 9.  NO SELECTION\002,5(/),\002"
	    " \002,25x,\002 Enter selection (1-9): \002,$)";
    static char fmt_2304[] = "(\002 \002,\002No selection made - Reduce numb"
	    "er of \002,\002loopings by one\002,/,\002 \002,\002Press <RETURN"
	    "> to continue: \002,$)";
    static char fmt_2306[] = "(\002 \002,\002Enter minimum value: \002,$)";
    static char fmt_2308[] = "(\002 \002,\002Enter maximum value: \002,$)";
    static char fmt_2310[] = "(\002 \002,\002Enter increment: \002,$)";
    static char fmt_2312[] = "(\002 \002,\002Zero increment not allowed: tro"
	    "uble with \002,\002dividing by zero\002)";
    static char fmt_2314[] = "(\002 \002,\002Run from \002,1pe12.5,\002 to"
	    " \002,1pe12.5,\002 in increments of \002,1pe12.5)";
    static char fmt_2316[] = "(\002 \002,\002Confirm these values (Y or N):"
	    " \002,$)";
    static char fmt_2301[] = "(a1)";
    static char fmt_2318[] = "(\002 \002,\002Run \002,a22,/,\002    from "
	    "\002,1pe12.5,\002 to \002,1pe12.5,\002 in increments of \002,1pe"
	    "12.5)";
    static char fmt_2320[] = "(\002 \002,\002No selection made - \002,\002Pr"
	    "ess <RETURN> to continue: \002,$)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsle(cilist *), do_lio(integer 
	    *, integer *, char *, ftnlen), e_rsle(void);
    double r_lg10(real *);
    integer i_nint(real *);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, l, inum[3], jnum, knum, lnum[3];
    static real rnum1[3], rnum2[3], rnum3[3];
    extern /* Subroutine */ int check_(void);
    static integer inumb;
#define qvary ((real *)&modpr_1 + 3)
    static integer lnumb1, lnumb2, lnumb3;
    static real rnumb1, rnumb2, rnumb3;
    static char lchose[1];
    extern /* Subroutine */ int driver_(void);

    /* Fortran I/O blocks */
    static cilist io___269 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___270 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___272 = { 0, 6, 0, fmt_2100, 0 };
    static cilist io___273 = { 0, 5, 0, 0, 0 };
    static cilist io___274 = { 0, 6, 0, fmt_2104, 0 };
    static cilist io___275 = { 0, 5, 0, 0, 0 };
    static cilist io___276 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___277 = { 0, 6, 0, fmt_2202, 0 };
    static cilist io___278 = { 0, 5, 0, 0, 0 };
    static cilist io___279 = { 0, 6, 0, fmt_2300, 0 };
    static cilist io___280 = { 0, 5, 0, 0, 0 };
    static cilist io___288 = { 0, 6, 0, fmt_2302, 0 };
    static cilist io___289 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___290 = { 0, 6, 0, fmt_2304, 0 };
    static cilist io___291 = { 0, 5, 0, 0, 0 };
    static cilist io___292 = { 0, 6, 0, fmt_2306, 0 };
    static cilist io___293 = { 0, 5, 0, 0, 0 };
    static cilist io___294 = { 0, 6, 0, fmt_2308, 0 };
    static cilist io___295 = { 0, 5, 0, 0, 0 };
    static cilist io___296 = { 0, 6, 0, fmt_2310, 0 };
    static cilist io___297 = { 0, 5, 0, 0, 0 };
    static cilist io___298 = { 0, 6, 0, fmt_2312, 0 };
    static cilist io___299 = { 0, 6, 0, fmt_2314, 0 };
    static cilist io___300 = { 0, 6, 0, fmt_2316, 0 };
    static cilist io___301 = { 0, 5, 0, fmt_2301, 0 };
    static cilist io___304 = { 0, 6, 0, fmt_2318, 0 };
    static cilist io___306 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___313 = { 0, 6, 0, fmt_2202, 0 };
    static cilist io___314 = { 0, 6, 0, fmt_2320, 0 };
    static cilist io___315 = { 0, 5, 0, 0, 0 };


/* ----------LINKAGES. */
/*     CALLED BY - [program] nuc123 */
/*     CALLS     - [subroutine] driver */
/* ----------REMARKS. */
/*     Activates computation routine. */
/* ----------PARAMETERS. */
/* Input unit number (previous value = */
/* Output unit number (previous value = */
/* Number of nuclear reactions. */
/* Total # of nuclear reactions for iru */
/* Total # of nuclear reactions for iru */
/* Number of nuclides in calculation. */
/* Total # of nuclides for irun = 2. */
/* ----------COMMON AREAS. */
/* Total # of nuclides for irun = 3. */
/* Model parameters. */
/* Variational parame */
/* Computation locati */
/* ==========================DECLARATION DIVISION================================= */
/* ----------MODEL PARAMETERS. */
/* Run options. */
/* Baryon-to-photon ratio. */
/* c(1) is variation of gravitational c */
/* c(2) is neutron lifetime (sec). */
/* c(3) is number of neutrino species. */
/* Cosmological constant. */
/* ----------RUN OPTION. */
/* Neutrino degeneracy parameters. */
/* Run network size. */
/* Number of nuclides in computation. */
/* ----------USER INTERACTION VARIABLES. */
/* Number of reactions in computation. */
/* Run parameter for outer loop. */
/* Run parameter for middle loop. */
/* Run parameter for inner loop. */
/* Run parameter starting value. */
/* Run parameter end value. */
/* Run parameter increment. */
/* Selection number. */
/* Selection number. */
/* Number of loopings to be done. */
/* Number of loopings rejected. */
/* Run parameter for outer loop. */
/* Run parameter for middle loop. */
/* Run parameter for inner loop. */
/* Run parameter end value. */
/* ----------FLAG AND LABELS. */
/* User response (alphanumeric). */
/* Computation location. */
/* ----------EQUIVALENCE VARIABLE. */
/* Label for quantities being varied. */
/* ----------EQUIVALENCE STATEMENTS. */
/* Array set equal to c, cosmo, and xi. */
/* ==============================DATA DIVISION==================================== */
/* ----------LABELS FOR QUANTITIES BEING VARIED. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------PRINT RUN SELECTION AND AWAIT RESPONSE------------------------------- */
/* ..........RETURN FROM LOOPING. */
L100:
/* ..........DISPLAY RUN SELECTIONS. */
    s_wsfe(&io___269);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___270);
    do_fio(&c__1, (char *)&inumb, (ftnlen)sizeof(integer));
    e_rsfe();
/* 20--------BRANCH TO APPROPRIATE SECTION---------------------------------------- */
    switch (inumb) {
	case 1:  goto L210;
	case 2:  goto L220;
	case 3:  goto L230;
	case 4:  goto L240;
    }
    goto L240;
/* 21--------SET RUN NETWORK SECTION---------------------------------------------- */
/* Improper input or <RETURN>. */
L210:
    s_wsfe(&io___272);
    e_wsfe();
    s_rsle(&io___273);
    do_lio(&c__3, &c__1, (char *)&inumb, (ftnlen)sizeof(integer));
    e_rsle();
/* Read in selection number. */
    if (inumb != 1 && inumb != 2 && inumb != 3) {
	inumb = 1;
    }
    if (inumb != runopt_1.irun) {
/* Run network changed from previously. */
	runopt_1.irun = inumb;
/* Run network size selection. */
    }
    if (runopt_1.irun == 1) {
/* Maximal network size. */
	runopt_1.isize = 26;
	runopt_1.jsize = 88;
    } else {
	if (runopt_1.irun == 2) {
/* Abridged network size. */
	    runopt_1.isize = 18;
	    runopt_1.jsize = 64;
	} else {
	    if (runopt_1.irun == 3) {
/* Minimal network size. */
		runopt_1.isize = 9;
		runopt_1.jsize = 34;
	    }
	}
    }
/* (irun.eq.1) */
    s_wsfe(&io___274);
    do_fio(&c__1, (char *)&runopt_1.irun, (ftnlen)sizeof(integer));
    e_wsfe();
    s_rsle(&io___275);
    e_rsle();
    goto L300;
/* 22--------GO SECTION----------------------------------------------------------- */
L220:
    s_wsfe(&io___276);
    e_wsfe();
    checkcb_1.itime = 3;
    check_();
/* Call interface subr before computati */
    driver_();
/* Do nucleosynthesis computation. */
    checkcb_1.itime = 8;
    check_();
/* Call interface subr after computatio */
    s_wsfe(&io___277);
    e_wsfe();
    s_rsle(&io___278);
    e_rsle();
    goto L300;
/* 23--------DO MULTIPLE RUNS SECTION--------------------------------------------- */
/* ..........GET NUMBER OF LOOPINGS. */
L230:
    s_wsfe(&io___279);
    e_wsfe();
    s_rsle(&io___280);
    do_lio(&c__3, &c__1, (char *)&jnum, (ftnlen)sizeof(integer));
    e_rsle();
/* Read in number of loopings to be don */
    if (jnum != 1 && jnum != 2 && jnum != 3) {
	jnum = 1;
/* Default number of loopings. */
    }
    knum = 0.f;
/* No loopings rejected for now. */
    for (i__ = 1; i__ <= 3; ++i__) {
	if (i__ > jnum) {
	    rnum1[i__ - 1] = 0.f;
/* Initialize initial parameter. */
	    rnum2[i__ - 1] = 0.f;
/* Initialize terminal parameter. */
	    rnum3[i__ - 1] = 1.f;
/* Initialize incremental parameter. */
	    inum[i__ - 1] = 0;
/* Initialize selection number. */
	} else {
/* ..........OBTAIN QUANTITY TO VARY. */
	    s_wsfe(&io___288);
	    e_wsfe();
	    s_rsfe(&io___289);
	    do_fio(&c__1, (char *)&inum[i__ - 1], (ftnlen)sizeof(integer));
	    e_rsfe();
	    if (inum[i__ - 1] < 1 || inum[i__ - 1] > 8) {
/* No selection ma */
		s_wsfe(&io___290);
		e_wsfe();
		s_rsle(&io___291);
		e_rsle();
		++knum;
/* Step up number of loopings rejected. */
		rnum1[i__ - 1] = 0.f;
/* Initialize initial parameter. */
		rnum2[i__ - 1] = 0.f;
/* Initialize terminal parameter. */
		rnum3[i__ - 1] = 1.f;
/* Initialize incremental parameter. */
		inum[i__ - 1] = 0;
/* Initialize selection number. */
	    } else {
/* ..........INPUT RUN SPECIFICATIONS. */
/* ((inum(i).ge.1).and.(inum(i).le.8)) */
L231:
		s_wsfe(&io___292);
		e_wsfe();
		s_rsle(&io___293);
		do_lio(&c__4, &c__1, (char *)&rnum1[i__ - 1], (ftnlen)sizeof(
			real));
		e_rsle();
/* Read in starting value. */
		s_wsfe(&io___294);
		e_wsfe();
		s_rsle(&io___295);
		do_lio(&c__4, &c__1, (char *)&rnum2[i__ - 1], (ftnlen)sizeof(
			real));
		e_rsle();
/* Read in terminating value. */
L232:
		s_wsfe(&io___296);
		e_wsfe();
		s_rsle(&io___297);
		do_lio(&c__4, &c__1, (char *)&rnum3[i__ - 1], (ftnlen)sizeof(
			real));
		e_rsle();
/* Read in incremental value. */
		if (rnum3[i__ - 1] == 0.f) {
/* Trouble with 0 division later on */
		    s_wsfe(&io___298);
		    e_wsfe();
		    goto L232;
		}
		s_wsfe(&io___299);
		do_fio(&c__1, (char *)&rnum1[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&rnum2[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&rnum3[i__ - 1], (ftnlen)sizeof(real));
		e_wsfe();
/* Display inp */
		s_wsfe(&io___300);
		e_wsfe();
		s_rsfe(&io___301);
		do_fio(&c__1, lchose, (ftnlen)1);
		e_rsfe();
/* Get confirm */
		if (*(unsigned char *)lchose != 'Y' && *(unsigned char *)
			lchose != 'y') {
		    goto L231;
		}
	    }
/* ((inum(i).lt.1).or.(inum(i).gt.8)) */
	}
/* (i.gt.jnum) */
    }
/* i = 1,3 */
    jnum -= knum;
/* Number of valid loopings. */
    if (jnum != 0) {
/* ..........WRITE OUT QUANTITY TO VARY, RUN SPECIFICATIONS. */
/* Run requested. */
	i__1 = jnum + knum;
	for (l = 1; l <= i__1; ++l) {
/* Check all loopings. */
	    if (inum[l - 1] != 0) {
/* Proper selection was made. */
		s_wsfe(&io___304);
		do_fio(&c__1, vtype + (inum[l - 1] - 1) * 22, (ftnlen)22);
		do_fio(&c__1, (char *)&rnum1[l - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&rnum2[l - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&rnum3[l - 1], (ftnlen)sizeof(real));
		e_wsfe();
/* Display run */
/* ..........GET LOGS OF eta VALUES FOR LOGRITHMIC INCREMENTATION. */
		if (inum[l - 1] == 1) {
/* Work with exponents for eta incre */
		    rnum1[l - 1] = r_lg10(&rnum1[l - 1]);
		    rnum2[l - 1] = r_lg10(&rnum2[l - 1]);
		}
	    }
	}
/* ..........COMPUTE NUMBER OF RUNS FOR EACH LOOPING. */
	for (l = 1; l <= 3; ++l) {
	    r__1 = (rnum2[l - 1] - rnum1[l - 1] + rnum3[l - 1]) / rnum3[l - 1]
		    ;
	    lnum[l - 1] = i_nint(&r__1);
	}
/* ..........DO MULTIPLE RUNS. */
	s_wsfe(&io___306);
	e_wsfe();
/* Inform user of beginning of computat */
	i__1 = lnum[0] - 1;
	for (lnumb1 = 0; lnumb1 <= i__1; ++lnumb1) {
/* Outer loop. */
	    rnumb1 = rnum1[0] + (real) lnumb1 * rnum3[0];
/* Value of param fo */
	    if (inum[0] >= 1 && inum[0] <= 8) {
		if (inum[0] == 1) {
		    d__1 = (doublereal) rnumb1;
		    varpr_1.eta1 = pow_dd(&c_b364, &d__1);
/* Vary baryon-to-photon ratio. */
		} else {
		    qvary[inum[0] - 2] = rnumb1;
/* Vary other quantities. */
		}
	    }
	    i__2 = lnum[1] - 1;
	    for (lnumb2 = 0; lnumb2 <= i__2; ++lnumb2) {
/* Middle loop. */
		rnumb2 = rnum1[1] + (real) lnumb2 * rnum3[1];
/* Value of param */
		if (inum[1] >= 1 && inum[1] <= 8) {
		    if (inum[1] == 1) {
			d__1 = (doublereal) rnumb2;
			varpr_1.eta1 = pow_dd(&c_b364, &d__1);
/* Vary baryon-to-photon ratio. */
		    } else {
			qvary[inum[1] - 2] = rnumb2;
/* Vary other quantities. */
		    }
		}
		i__3 = lnum[2] - 1;
		for (lnumb3 = 0; lnumb3 <= i__3; ++lnumb3) {
/* Inner loop. */
		    rnumb3 = rnum1[2] + (real) lnumb3 * rnum3[2];
/* Value of para */
		    if (inum[2] >= 1 && inum[2] <= 8) {
			if (inum[2] == 1) {
			    d__1 = (doublereal) rnumb3;
			    varpr_1.eta1 = pow_dd(&c_b364, &d__1);
/* Vary baryon-to-photon ratio. */
			} else {
			    qvary[inum[2] - 2] = rnumb3;
/* Vary other quantities. */
			}
		    }
		    checkcb_1.itime = 3;
		    check_();
/* Check interface subr before computat */
		    driver_();
/* Do nucleosynthesis computation. */
		    checkcb_1.itime = 8;
		    check_();
/* Check interface subroutine after com */
		}
/* lnumb3 = 0,lnum(3)-1 */
	    }
/* lnumb2 = 0,lnum(2)-1 */
	}
/* lnumb1 = 0,lnum(1)-1 */
	s_wsfe(&io___313);
	e_wsfe();
/* Inform user of completion of computa */
    } else {
/* (jnum.eq.0) */
	s_wsfe(&io___314);
	e_wsfe();
    }
/* (jnum.ne.0) */
    s_rsle(&io___315);
    e_rsle();
    goto L300;
/* 24--------EXIT SECTION--------------------------------------------------------- */
L240:
    return 0;
/* 30--------GO BACK TO MENU------------------------------------------------------ */
L300:
    goto L100;
} /* run_ */

#undef qvary


/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int output_(void)
{
    /* Format strings */
    static char fmt_1000[] = "(8(/),\002 \002,30x,\002OUTPUT SELECTION\002"
	    ",/,\002 \002,30x,\002------ ---------\002,//,\002 \002,25x,\002 "
	    "1. REQUEST OUTPUT FILE\002,/,\002 \002,25x,\002 2. REQUEST OUTPU"
	    "T ON SCREEN\002,/,\002 \002,25x,\002 3. EXIT\002,11(/),\002 \002"
	    ",25x,\002 Enter selection (1-3): \002,$)";
    static char fmt_1001[] = "(i1)";
    static char fmt_2000[] = "(54x,\002NUCLIDE ABUNDANCE YIELDS\002,/,54x"
	    ",\002------- --------- ------\002,//)";
    static char fmt_2002[] = "(\002 Computational parameters:\002,/,\002   c"
	    "y = \002,f5.3,\002/  ct = \002,f5.3,\002/  initial temp = \002,1"
	    "pe8.2,\002/  final temp = \002,1pe8.2,\002/  smallest abundances"
	    " allowed = \002,1pe8.2)";
    static char fmt_2004[] = "(\002 Model parameters:\002,/,\002   g = \002,"
	    "f5.2,\002/  tau = \002,f6.2,\002/  # nu = \002,f5.2,\002/  lambd"
	    "a = \002,1pe10.3,\002/  xi-e = \002,e10.3,\002/  xi-m = \002,e10"
	    ".3,\002/  xi-t = \002,e10.3,/)";
    static char fmt_2006[] = "(4x,\002Temp\002,8x,\002N/H\002,10x,\002P\002,"
	    "10x,\002D/H\002,9x,\002T/H\002,8x,\002He3/H\002,8x,\002He4\002,8"
	    "x,\002Li6/H\002,7x,\002Li7/H\002,7x,\002Be7/H\002,6x,\002Li8/H&up"
	    "\002,/,132(\002-\002))";
    static char fmt_2008[] = "(1pe10.3,1p10e12.3)";
    static char fmt_2010[] = "(\002 \002,/,4x,\002Temp\002,9x,\002T\002,10x"
	    ",\002rhog\002,8x,\002rhoe\002,7x,\002rhone\002,8x,\002rhob\002,8"
	    "x,\002phie\002,9x,\002dt\002,9x,\002eta\002,10x,\002H\002,/,132"
	    "(\002-\002))";
    static char fmt_2012[] = "(1pe10.3,9e12.3)";
    static char fmt_2014[] = "(///)";
    static char fmt_2016[] = "(\002 \002,\002Output file requested - Press <"
	    "RETURN> to continue: \002,$)";
    static char fmt_3000[] = "(8(/),\002 \002,26x,\002SCREEN OUTPUT SELECT"
	    "ION\002,/,\002 \002,26x,\002------ ------ ---------\002,//,\002"
	    " \002,25x,\002 1. DISPLAY D,T,HE3,HE4,LI7\002,/,\002 \002,25x"
	    ",\002 2. DISPLAY N,P,LI6,BE7,LI8&UP\002,/,\002 \002,25x,\002 3. "
	    "DISPLAY RHOG,RHOE,RHONE,RHOB\002,/,\002 \002,25x,\002 4. DISPLAY"
	    " T,DT,PHIE,ETA,H\002,/,\002 \002,25x,\002 5. EXIT\002,9(/),\002"
	    " \002,25x,\002 Enter selection (1-5): \002,$)";
    static char fmt_3100[] = "(\002 \002,\002Computational parameters:\002"
	    ",/,\002 \002,\002   cy = \002,f5.3,\002/ ct = \002,f5.3,\002/ in"
	    "itial temp = \002,1pe8.2,\002/ final temp = \002,1pe8.2,/,\002"
	    " \002,\002   smallest abundances allowed = \002,1pe8.2)";
    static char fmt_3102[] = "(\002 \002,\002Model parameters:\002,/,\002"
	    " \002,\002   g = \002,f5.2,\002/ tau = \002,f6.2,\002/ # nu ="
	    " \002,f5.2,\002/ lambda = \002,1pe10.3,/,\002 \002,\002   xi-e = "
	    "\002,e10.3,\002/ xi-m = \002,e10.3,\002/ xi-t = \002,e10.3,/)";
    static char fmt_3104[] = "(4x,\002Temp\002,8x,\002D/H\002,9x,\002T/H\002"
	    ",8x,\002He3/H\002,8x,\002He4\002,8x,\002Li7/H\002,/,\002 \002,80("
	    "\002-\002))";
    static char fmt_3106[] = "(1pe10.3,1p5e12.3)";
    static char fmt_3108[] = "(\002 \002,\002Press <RETURN> to continue: "
	    "\002,$)";
    static char fmt_3204[] = "(4x,\002Temp\002,8x,\002N/H\002,10x,\002P\002,"
	    "9x,\002Li6/H\002,7x,\002Be7/H\002,6x,\002Li8/H&up\002,/,\002 "
	    "\002,80(\002-\002))";
    static char fmt_3304[] = "(4x,\002Temp\002,8x,\002rhog\002,8x,\002rho"
	    "e\002,7x,\002rhone\002,8x,\002rhob\002,/,\002 \002,80(\002-\002))"
	    ;
    static char fmt_3306[] = "(1pe10.3,4e12.3)";
    static char fmt_3404[] = "(4x,\002Temp\002,8x,\002time\002,8x,\002phi"
	    "e\002,9x,\002dt\002,9x,\002eta\002,10x,\002H\002,/,\002 \002,80"
	    "(\002-\002))";
    static char fmt_3406[] = "(1pe10.3,5e12.3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsle(cilist *), e_rsle(void);

    /* Local variables */
    static integer i__, j, inum;

    /* Fortran I/O blocks */
    static cilist io___316 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___317 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___319 = { 0, 2, 0, fmt_2000, 0 };
    static cilist io___320 = { 0, 2, 0, fmt_2002, 0 };
    static cilist io___321 = { 0, 2, 0, fmt_2004, 0 };
    static cilist io___322 = { 0, 2, 0, fmt_2006, 0 };
    static cilist io___324 = { 0, 2, 0, fmt_2008, 0 };
    static cilist io___326 = { 0, 2, 0, fmt_2010, 0 };
    static cilist io___327 = { 0, 2, 0, fmt_2012, 0 };
    static cilist io___328 = { 0, 2, 0, fmt_2014, 0 };
    static cilist io___329 = { 0, 6, 0, fmt_2016, 0 };
    static cilist io___330 = { 0, 5, 0, 0, 0 };
    static cilist io___331 = { 0, 6, 0, fmt_3000, 0 };
    static cilist io___332 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___333 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___334 = { 0, 6, 0, fmt_3100, 0 };
    static cilist io___335 = { 0, 6, 0, fmt_3102, 0 };
    static cilist io___336 = { 0, 6, 0, fmt_3104, 0 };
    static cilist io___337 = { 0, 6, 0, fmt_3106, 0 };
    static cilist io___338 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___339 = { 0, 6, 0, fmt_3108, 0 };
    static cilist io___340 = { 0, 5, 0, 0, 0 };
    static cilist io___341 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___342 = { 0, 6, 0, fmt_3100, 0 };
    static cilist io___343 = { 0, 6, 0, fmt_3102, 0 };
    static cilist io___344 = { 0, 6, 0, fmt_3204, 0 };
    static cilist io___345 = { 0, 6, 0, fmt_3106, 0 };
    static cilist io___346 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___347 = { 0, 6, 0, fmt_3108, 0 };
    static cilist io___348 = { 0, 5, 0, 0, 0 };
    static cilist io___349 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___350 = { 0, 6, 0, fmt_3100, 0 };
    static cilist io___351 = { 0, 6, 0, fmt_3102, 0 };
    static cilist io___352 = { 0, 6, 0, fmt_3304, 0 };
    static cilist io___353 = { 0, 6, 0, fmt_3306, 0 };
    static cilist io___354 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___355 = { 0, 6, 0, fmt_3108, 0 };
    static cilist io___356 = { 0, 5, 0, 0, 0 };
    static cilist io___357 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___358 = { 0, 6, 0, fmt_3100, 0 };
    static cilist io___359 = { 0, 6, 0, fmt_3102, 0 };
    static cilist io___360 = { 0, 6, 0, fmt_3404, 0 };
    static cilist io___361 = { 0, 6, 0, fmt_3406, 0 };
    static cilist io___362 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___363 = { 0, 6, 0, fmt_3108, 0 };
    static cilist io___364 = { 0, 5, 0, 0, 0 };


/* ----------LINKAGES. */
/*     CALLED BY - [program] nuc123 */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Outputs computational results either into an output file or onto */
/*     the screen */
/* ----------PARAMETERS. */
/* Input unit number (previous value = */
/* Output unit number (previous value = */
/* Number of nuclides in calculation. */
/* ----------COMMON AREAS. */
/* Maximum # of line to be printed. */
/* Computation parame */
/* Model parameters. */
/* Flags, counters. */
/* Output data. */
/* ==========================DECLARATION DIVISION================================= */
/* ----------COMPUTATION SETTINGS. */
/* Output option. */
/* Time step limiting constant on abund */
/* Time step limiting constant on tempe */
/* Initial temperature (in 10**9 K). */
/* Final temperature (in 10**9 K). */
/* ----------EARLY UNIVERSE MODEL PARAMETERS. */
/* Smallest abundances allowed. */
/* c(1) is variation of gravitational c */
/* c(2) is neutron lifetime (sec). */
/* c(3) is number of neutrino species. */
/* Cosmological constant. */
/* ----------COUNTER. */
/* Neutrino degeneracy parameters. */
/* ----------OUTPUT ARRAYS. */
/* # times accumulated in output buffer */
/* Nuclide mass fractions. */
/* Thermodynamic variables. */
/* Temperature (in units of 10**9 K). */
/* Time. */
/* Time step. */
/* Baryon-to-photon ratio. */
/* ----------OUTPUT FILE STATUS. */
/* Expansion rate. */
/* Number of output requests. */
/* ----------USER INTERACTION VARIABLES. */
/* Indicates if output file used. */
/* ===========================PROCEDURE DIVISION================================== */
/* 10--------PRINT OUTPUT SELECTION AND AWAIT RESPONSE---------------------------- */
/* ..........RETURN FROM LOOPING. */
/* Selection number. */
L100:
/* ..........DISPLAY OUTPUT SELECTIONS. */
    s_wsfe(&io___316);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___317);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
/* ..........BRANCH TO APPROPRIATE SECTION. */
    switch (inum) {
	case 1:  goto L200;
	case 2:  goto L300;
	case 3:  goto L400;
    }
    goto L400;
/* 20--------REQUEST OUTPUT SECTION----------------------------------------------- */
/* Improper input or <RETURN>. */
L200:
/*      DO j = 1,it                  !Temperature in MeV. */
/*        t9out(j) = t9out(j)*.08617 */
/*      END DO */
/*      DO j = 1,it                  !Energy density as fraction of total. */
/*        thmout(j,1) = thmout(j,1)/thmout(j,6)  !Rhog. */
/*        thmout(j,2) = thmout(j,2)/thmout(j,6)  !Rhoe. */
/*        thmout(j,3) = thmout(j,3)/thmout(j,6)  !Rhone. */
/*        thmout(j,4) = thmout(j,4)/thmout(j,6)  !Rhob. */
/*      END DO */
/* ..........PRINT CAPTION. */
    ++outopt_1.nout;
/* Keep track of number of output reque */
    if (outopt_1.nout == 1) {
	s_wsfe(&io___319);
	e_wsfe();
    }
    s_wsfe(&io___320);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___321);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT HEADINGS, ABUNDANCES FOR NEUTRON TO LI8. */
    s_wsfe(&io___322);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___324);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&outdat_1.xout[j + i__ * 40 - 41], (ftnlen)
		    sizeof(real));
	}
	e_wsfe();
    }
/* ..........PRINT THERMODYNAMIC QUANTITIES. */
    s_wsfe(&io___326);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___327);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.tout[j - 1], (ftnlen)sizeof(real));
	for (i__ = 1; i__ <= 5; ++i__) {
	    do_fio(&c__1, (char *)&outdat_1.thmout[j + i__ * 40 - 41], (
		    ftnlen)sizeof(real));
	}
	do_fio(&c__1, (char *)&outdat_1.dtout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.etaout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.hubout[j - 1], (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___328);
    e_wsfe();
    outopt_1.outfile = TRUE_;
/* Output file requested. */
    s_wsfe(&io___329);
    e_wsfe();
    s_rsle(&io___330);
    e_rsle();
    goto L500;
/* 30--------REQUEST OUTPUT ON SCREEN SECTION------------------------------------- */
/* ..........RETURN FROM LOOPING. */
L300:
/*      DO j = 1,it                  !Temperature in MeV. */
/*        t9out(j) = t9out(j)*.08617 */
/*      END DO */
/*      DO j = 1,it                  !Energy density as fraction of total. */
/*        thmout(j,1) = thmout(j,1)/thmout(j,6)  !Rhog. */
/*        thmout(j,2) = thmout(j,2)/thmout(j,6)  !Rhoe. */
/*        thmout(j,3) = thmout(j,3)/thmout(j,6)  !Rhone. */
/*        thmout(j,4) = thmout(j,4)/thmout(j,6)  !Rhob. */
/*      END DO */
/* ..........DISPLAY SCREEN OUTPUT SELECTIONS. */
    s_wsfe(&io___331);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___332);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    switch (inum) {
	case 1:  goto L310;
	case 2:  goto L320;
	case 3:  goto L330;
	case 4:  goto L340;
	case 5:  goto L350;
    }
    goto L350;
/* Improper input or <RETURN>. */
L310:
/* ..........PRINT CAPTION. */
/* Display d,t,he3,he4,li7. */
    s_wsfe(&io___333);
    e_wsfe();
    s_wsfe(&io___334);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___335);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT HEADINGS, ABUNDANCES FOR D,T,HE3,HE4,LI7. */
    s_wsfe(&io___336);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___337);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	for (i__ = 3; i__ <= 6; ++i__) {
	    do_fio(&c__1, (char *)&outdat_1.xout[j + i__ * 40 - 41], (ftnlen)
		    sizeof(real));
	}
	do_fio(&c__1, (char *)&outdat_1.xout[j + 279], (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___338);
    e_wsfe();
    s_wsfe(&io___339);
    e_wsfe();
    s_rsle(&io___340);
    e_rsle();
    goto L360;
L320:
/* ..........PRINT CAPTION. */
/* Display n,p,li6,be7,li8&up. */
    s_wsfe(&io___341);
    e_wsfe();
    s_wsfe(&io___342);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___343);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT HEADINGS, ABUNDANCES FOR N,P,LI6,BE7,LI8&UP. */
    s_wsfe(&io___344);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___345);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	for (i__ = 1; i__ <= 2; ++i__) {
	    do_fio(&c__1, (char *)&outdat_1.xout[j + i__ * 40 - 41], (ftnlen)
		    sizeof(real));
	}
	do_fio(&c__1, (char *)&outdat_1.xout[j + 239], (ftnlen)sizeof(real));
	for (i__ = 9; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&outdat_1.xout[j + i__ * 40 - 41], (ftnlen)
		    sizeof(real));
	}
	e_wsfe();
    }
    s_wsfe(&io___346);
    e_wsfe();
    s_wsfe(&io___347);
    e_wsfe();
    s_rsle(&io___348);
    e_rsle();
    goto L360;
L330:
/* ..........PRINT CAPTION. */
/* Display rhog,rhoe,rhone,rhob. */
    s_wsfe(&io___349);
    e_wsfe();
    s_wsfe(&io___350);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___351);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT ENERGY DENSITIES. */
    s_wsfe(&io___352);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___353);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	for (i__ = 1; i__ <= 4; ++i__) {
	    do_fio(&c__1, (char *)&outdat_1.thmout[j + i__ * 40 - 41], (
		    ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    s_wsfe(&io___354);
    e_wsfe();
    s_wsfe(&io___355);
    e_wsfe();
    s_rsle(&io___356);
    e_rsle();
    goto L360;
L340:
/* ..........PRINT CAPTION. */
/* Display t,dt,phie,eta,hubcst. */
    s_wsfe(&io___357);
    e_wsfe();
    s_wsfe(&io___358);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___359);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT THERMODYNAMIC QUANTITIES. */
    s_wsfe(&io___360);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___361);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.tout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.thmout[j + 159], (ftnlen)sizeof(real))
		;
	do_fio(&c__1, (char *)&outdat_1.dtout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.etaout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.hubout[j - 1], (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___362);
    e_wsfe();
    s_wsfe(&io___363);
    e_wsfe();
    s_rsle(&io___364);
    e_rsle();
    goto L360;
L350:
/* Exit. */
    goto L500;
L360:
    goto L300;
/* 40--------EXIT SECTION--------------------------------------------------------- */
L400:
    return 0;
/* 50--------GO BACK TO MENU------------------------------------------------------ */
L500:
    goto L100;
} /* output_ */

/* Main program alias */ int new123_ () { MAIN__ (); return 0; }
