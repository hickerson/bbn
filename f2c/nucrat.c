/* nucrat.f -- translated by f2c (version 20100827).
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
    real f[88], r__[88];
} rates_;

#define rates_1 rates_

struct {
    real g, tau, xnu, c__[3], cosmo, xi[3];
} modpr_;

#define modpr_1 modpr_

struct {
    real thm[14], hubcst;
} thermcb_;

#define thermcb_1 thermcb_

struct {
    real t9mev, tnmev, tnu, cnorm, rhonu;
    integer nu;
} nupar_;

#define nupar_1 nupar_

struct {
    real t9, hv, phie, y[26];
} evolp1_;

#define evolp1_1 evolp1_

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

struct nucdat_1_ {
    real am[26], zm[26], dm[26];
};

#define nucdat_1 (*(struct nucdat_1_ *) &nucdat_)

/* Initialized data */

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
    } modpr0_ = { 1.f, 885.7f, 3.f, 0.f, 0.f, 0.f, 0.f };

struct {
    real e_1[2];
    } varpr0_ = { 1e-4f, 3.162e-10f };

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


/* Table of constant values */

static real c_b3 = 1.f;
static integer c__50 = 50;
static doublereal c_b12 = .33333333;
static doublereal c_b13 = 1.5;
static doublereal c_b15 = .66666666666666663;
static doublereal c_b16 = .3333333;
static doublereal c_b17 = .8333333;
static doublereal c_b26 = .2;
static doublereal c_b27 = 1.25;
static doublereal c_b32 = .6;
static doublereal c_b33 = 1.2;
static doublereal c_b34 = .375;
static doublereal c_b37 = .83333333;

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
    extern /* Subroutine */ int func1_(), func2_(), func3_(), func4_();
    static real part1, part2, part3, part4;
    extern doublereal xintd_(real *, real *, U_fp, integer *);
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
/* Dynamic variables. */
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
/* Parts of integrals for forward rate. */
/* Parts of integrals for reverse rate. */
    if (modpr_1.xi[0] == 0.f) {
	rates_1.f[0] = thermcb_1.thm[12] / modpr_1.tau;
/* Forward rate for weak np reaction. */
	rates_1.r__[0] = thermcb_1.thm[13] / modpr_1.tau;
/* Reverse rate for weak np reaction. */
    } else {
/* 20--------COMPUTE WEAK REACTION RATES (DEGENERATE)----------------------------- */
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
	z__[1] = -(nupar_1.tnmev / .511f) * (-88.722f - modpr_1.xi[0]) + 
		2.531f;
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
	part1 = xintd_(&c_b3, &uplim1, (U_fp)func1_, &c__50);
	part2 = xintd_(&c_b3, &uplim2, (U_fp)func2_, &c__50);
	part3 = xintd_(&c_b3, &uplim3, (U_fp)func3_, &c__50);
	part4 = xintd_(&c_b3, &uplim4, (U_fp)func4_, &c__50);
	rates_1.f[0] = part1 + part2;
/* Add 2 integrals to get forward rate. */
	rates_1.r__[0] = part3 + part4;
/* Add 2 integrals to get reverse rate. */
    }
/* (xi(1).eq.0.) */
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
    t913 = pow_dd(&d__1, &c_b12);
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
    t9a32 = pow_dd(&d__1, &c_b13);
/* t9a**(3/2) */
    t9b = evolp1_1.t9 / (evolp1_1.t9 * 49.18f + 1.f);
/* For reaction 18. */
    d__1 = (doublereal) t9b;
    t9b32 = pow_dd(&d__1, &c_b13);
/* t9b**(3/2) */
    if (evolp1_1.t9 > 10.f) {
/* For reaction 22. */
	t9c = 1.f;
    } else {
	d__1 = (doublereal) (1.f - evolp1_1.t9 * .0969f);
	t9c = evolp1_1.t9 / (1.f - evolp1_1.t9 * .0969f + t953 * .0284f / 
		pow_dd(&d__1, &c_b15));
    }
    d__1 = (doublereal) t9c;
    t9c13 = pow_dd(&d__1, &c_b16);
/* t9c**(1/3) */
    d__1 = (doublereal) t9c;
    t9c56 = pow_dd(&d__1, &c_b17);
/* t9c**(5/6) */
    t9d = evolp1_1.t9 / (evolp1_1.t9 * .759f + 1.f);
/* For reaction 24. */
    d__1 = (doublereal) t9d;
    t9d13 = pow_dd(&d__1, &c_b16);
/* t9d**(1/3) */
    d__1 = (doublereal) t9d;
    t9d56 = pow_dd(&d__1, &c_b17);
/* t9d**(5/6) */
    t9e = evolp1_1.t9 / (evolp1_1.t9 * .1378f + 1.f);
/* For reaction 26. */
    d__1 = (doublereal) t9e;
    t9e13 = pow_dd(&d__1, &c_b16);
/* t9e**(1/3) */
    d__1 = (doublereal) t9e;
    t9e56 = pow_dd(&d__1, &c_b17);
/* t9e**(5/6) */
    t9f = evolp1_1.t9 / (evolp1_1.t9 * .1071f + 1.f);
/* For reaction 27. */
    d__1 = (doublereal) t9f;
    t9f13 = pow_dd(&d__1, &c_b16);
/* t9f**(1/3) */
    d__1 = (doublereal) t9f;
    t9f56 = pow_dd(&d__1, &c_b17);
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
    t913 = pow_dd(&d__1, &c_b12);
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
    t915 = pow_dd(&d__1, &c_b26);
/* t9**(1/5) */
    d__1 = (doublereal) evolp1_1.t9;
    t954 = pow_dd(&d__1, &c_b27);
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
    t9a13 = pow_dd(&d__1, &c_b16);
/* t9a**(1/3) */
    d__1 = (doublereal) t9a;
    t9a56 = pow_dd(&d__1, &c_b17);
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
    t913 = pow_dd(&d__1, &c_b12);
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
    t935 = pow_dd(&d__1, &c_b32);
/* t9**(3/5) */
    d__1 = (doublereal) evolp1_1.t9;
    t965 = pow_dd(&d__1, &c_b33);
/* t9**(6/5) */
    d__1 = (doublereal) evolp1_1.t9;
    t938 = pow_dd(&d__1, &c_b34);
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
	    &d__1, &c_b15));
/* For reaction 82. */
    d__1 = (doublereal) t9a;
    t9a13 = pow_dd(&d__1, &c_b12);
/* t9a**(1/3) */
    d__1 = (doublereal) t9a;
    t9a56 = pow_dd(&d__1, &c_b37);
/* t9a**(5/6) */
    d__1 = (doublereal) (evolp1_1.t9 * .0776f + 1.f);
    t9b = evolp1_1.t9 / (evolp1_1.t9 * .0776f + 1.f + t953 * .0264f / pow_dd(&
	    d__1, &c_b15));
/* For reaction 84. */
    d__1 = (doublereal) t9b;
    t9b13 = pow_dd(&d__1, &c_b12);
/* t9b**(1/3) */
    d__1 = (doublereal) t9b;
    t9b56 = pow_dd(&d__1, &c_b37);
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

