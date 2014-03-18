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

extern struct {
    doublereal f[88], r__[88];
} rates_;

#define rates_1 rates_

extern struct {
    doublereal g, tau, xnu, c__[3], cosmo, xi[3];
} modpr_;

#define modpr_1 modpr_

extern struct {
    doublereal thm[14], hubcst;
} thermcb_;

#define thermcb_1 thermcb_

extern struct {
    doublereal t9mev, tnmev, tnu, cnorm, rhonu;
    shortint nu;
} nupar_;

#define nupar_1 nupar_

extern struct {
    doublereal t9, hv, phie, y[26];
} evolp1_;

#define evolp1_1 evolp1_

struct recpr0_1_ {
    doublereal reacpr[704]	/* was [88][8] */;
};

#define recpr0_1 (*(struct recpr0_1_ *) &recpr0_)

struct compr0_1_ {
    doublereal cy0, ct0, t9i0, t9f0, ytmin0;
    shortint inc0;
};

#define compr0_1 (*(struct compr0_1_ *) &compr0_)

struct modpr0_1_ {
    doublereal c0[3], cosmo0, xi0[3];
};

#define modpr0_1 (*(struct modpr0_1_ *) &modpr0_)

struct varpr0_1_ {
    doublereal dt0, eta0;
};

#define varpr0_1 (*(struct varpr0_1_ *) &varpr0_)

struct nucdat_1_ {
    doublereal am[26], zm[26], dm[26];
};

#define nucdat_1 (*(struct nucdat_1_ *) &nucdat_)

/* Initialized data */

struct {
    doublereal e_1[704];
    } recpr0_ = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.,
	     15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 
	    28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 
	    41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 
	    54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 
	    67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 
	    80., 81., 82., 83., 84., 85., 86., 87., 88., 1., 1., 4., 1., 1., 
	    4., 1., 1., 1., 1., 1., 2., 2., 2., 2., 3., 3., 3., 5., 2., 2., 
	    2., 3., 5., 2., 2., 2., 6., 6., 3., 3., 11., 9., 9., 2., 2., 2., 
	    3., 3., 2., 2., 2., 2., 2., 3., 3., 3., 3., 2., 2., 2., 3., 3., 
	    3., 3., 3., 3., 8., 7., 9., 9., 9., 10., 10., 2., 2., 2., 3., 3., 
	    3., 3., 2., 2., 2., 2., 2., 2., 3., 2., 3., 3., 3., 3., 3., 3., 
	    3., 3., 3., 1., 4., 10., 16., 21., 11., 15., 18., 20., 23., 25., 
	    2., 3., 5., 7., 5., 9., 7., 9., 3., 4., 7., 7., 8., 6., 6., 6., 
	    3., 3., 4., 5., 5., 8., 9., 8., 13., 14., 15., 13., 9., 12., 13., 
	    14., 15., 16., 12., 13., 16., 7., 8., 9., 11., 10., 12., 12., 13.,
	     14., 6., 6., 10., 11., 12., 14., 15., 17., 19., 22., 20., 22., 
	    25., 25., 17., 19., 21., 20., 22., 24., 24., 17., 13., 14., 15., 
	    18., 20., 13., 14., 16., 19., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 2., 2., 2., 2., 2., 3., 
	    4., 5., 0., 0., 3., 3., 0., 3., 3., 1., 1., 1., 1., 1., 2., 2., 
	    2., 2., 2., 2., 2., 2., 2., 6., 6., 6., 6., 6., 6., 3., 3., 3., 
	    1., 0., 2., 1., 2., 2., 1., 1., 1., 1., 1., 1., 1., 1., 2., 2., 
	    2., 2., 2., 2., 2., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 
	    4., 0., 0., 0., 0., 5., 0., 0., 0., 0., 1., 2., 1., 2., 2., 1., 
	    2., 0., 0., 0., 2., 6., 0., 0., 0., 0., 0., 1., 6., 6., 6., 0., 
	    0., 0., 2., 1., 1., 1., 2., 1., 0., 0., 1., 2., 3., 0., 0., 0., 
	    0., 0., 2., 2., 2., 6., 0., 0., 0., 0., 0., 0., 6., 0., 2., 2., 
	    2., 2., 2., 1., 1., 1., 1., 2., 5., 6., 17., 22., 6., 14., 17., 
	    19., 22., 24., 3., 4., 6., 8., 4., 8., 6., 6., 5., 6., 9., 6., 6.,
	     7., 8., 9., 5., 4., 6., 6., 6., 6., 6., 10., 14., 16., 14., 8., 
	    11., 13., 15., 17., 18., 17., 7., 9., 12., 13., 14., 15., 15., 
	    14., 17., 13., 14., 17., 12., 17., 6., 6., 6., 6., 6., 19., 21., 
	    24., 19., 21., 24., 17., 20., 22., 24., 23., 25., 26., 17., 26., 
	    19., 21., 22., 25., 26., 20., 22., 24., 26., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., .471, 1.63, 2.61, 1.19, 1.002, .998, 1.07,
	     4.7, 1.63, 2.61, 1.19, 1.07, 4.69, 1.53, 1.11, 1.11, 1.73, 1.73, 
	    5.54, 5.55, 3.39, 9.95, 9.97, 1.31, 3.04, 2.34, 1.002, .758, 1.3, 
	    .973, 3.03, 7.01, 2.33, 3., .618, .754, .292, 1.58, 4.02, 4.02, 
	    3.08, 3.07, 10.3, 2.07, 6.44, 14.9, .584, 2., 3.58, 3.58, .807, 
	    3.5, 3.49, .886, 3.58, 2.71, 1.002, 3.003, 1.002, .709, .884, 
	    1.19, .9, 3.57, 2.7, 3.62, .706, 5.13, 9.36, 11., 3.68, 4.26, 
	    5.81, 9.34, 3.67, 4.25, 5.79, 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 25.82, 72.62, 238.81, 84.17, 8.863, 19.081, 55.494, 
	    220.39, 63.75, 229.932, 65.054, 46.631, 201.291, 17.118, 28.64, 
	    18.423, 37.935, 46.798, 204.117, 212.98, 149.23, 175.476, 194.557,
	     23.59, 132.95, 39.1, 32.08, 32.382, 1.595, 76.427, 100.84, 
	    185.173, 6.975, 146.08, 24.674, 13.301, 79.89, 51.753, 100.538, 
	    87.539, 86., 76.96, 66.16, 50.63, 107.13, 159.36, 18.26, 84.42, 
	    177.73, 218.82, 7.555, 100.753, 132.83, 57.41, 94.88, 125.74, 
	    34.846, 7.263, 41.037, 98.661, 22.553, 87.621, 118.452, 53.706, 
	    84.678, 140.734, 57.623, 83.111, 47.16, 9.098, 33.915, 111.87, 
	    60.557, 12.287, 1.835, 88.47, 25.711 };

struct {
    doublereal e_1[5];
    shortint e_2;
    } compr0_ = { .3, .03, 100., .01, 1e-25, 30 };

struct {
    doublereal e_1[7];
    } modpr0_ = { 1., 885.7, 3., 0., 0., 0., 0. };

struct {
    doublereal e_1[2];
    } varpr0_ = { 1e-4, 3.162e-10 };

struct {
    doublereal e_1[78];
    } nucdat_ = { 1., 1., 2., 3., 3., 4., 6., 7., 7., 8., 8., 9., 10., 11., 
	    11., 12., 12., 12., 13., 13., 14., 14., 14., 15., 15., 16., 0., 
	    1., 1., 1., 2., 2., 3., 3., 4., 3., 5., 4., 5., 5., 6., 5., 6., 
	    7., 6., 7., 6., 7., 8., 7., 8., 8., .008665, .007825, .014102, 
	    .01605, .01603, .002603, .015125, .016004, .016929, .022487, 
	    .024609, .012186, .012939, .009305, .011432, .014354, 0., .018641,
	     .003354, .005738, .003242, .003074, .008597, 1.08e-4, .00307, 
	    -.005085 };


/* Table of constant values */

static doublereal c_b3 = 1.;
static shortint cs__50 = 50;
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
    rates_1.f[1] = 1.79e-9;
/* .......Li8 -> e- + v + 2He4.......(Ajzenberg-Selove 1988) */
    rates_1.f[2] = .827;
/* .......B12 -> e- + B + C12........(Ajzenberg-Selove 1990) */
    rates_1.f[3] = 34.3;
/* .......C14 -> e- + v + N14........(Ajzenberg-Selove 1986) */
    rates_1.f[4] = 3.834e-12;
/* .......B8 -> e+ + v + 2He4........(Ajzenberg-Selove 1988) */
    rates_1.f[5] = .9;
/* .......C11 -> e+ + v + B11........(Ajzenberg-Selove 1990) */
    rates_1.f[6] = 5.668e-4;
/* .......N12 -> e+ + v + C12........(Ajzenberg-Selove 1990) */
    rates_1.f[7] = 63.01;
/* .......N13 -> e+ + v + C13........(Ajzenberg-Selove 1986) */
    rates_1.f[8] = .001159;
/* .......O14 -> e+ + v + N14........(Ajzenberg-Selove 1986) */
    rates_1.f[9] = .0098171;
/* .......O15 -> e+ + v + N15........(Ajzenberg-Selove 1986) */
    rates_1.f[10] = .0056704;
    return 0;
/* ----------REFERENCES----------------------------------------------------------- */
/*     Ajzenberg-Selove, F., 1990, Nucl. Phys. A506, 1. */
/*     Ajzenberg-Selove, F., 1988, Nucl. Phys. A490, 1. */
/*     Ajzenberg-Selove, F., 1986, Nucl. Phys. A449, 1. */
/*     Tilley, D.R., Weller, H.R., and Hasan, H.H., 1987, Nucl. Phys. A474, 1. */
} /* rate0_ */

/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int rate1_(doublereal *tph)
{
    static doublereal w[2], x[2], y[2], z__[2];
    extern /* Subroutine */ int func1_(), func2_(), func3_(), func4_();
    static doublereal part1, part2, part3, part4;
    extern doublereal xintd_(doublereal *, doublereal *, U_fp, shortint *);
    static doublereal uplim1, uplim2, uplim3, uplim4;

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
    if (modpr_1.xi[0] == 0.) {
	rates_1.f[0] = thermcb_1.thm[12] / modpr_1.tau;
/* Forward rate for weak np reaction. */
	rates_1.r__[0] = thermcb_1.thm[13] / modpr_1.tau;
/* Reverse rate for weak np reaction. */
    } else {
/* 20--------COMPUTE WEAK REACTION RATES (DEGENERATE)----------------------------- */
	nupar_1.t9mev = *tph * .086171;
/* Convert photon temp to units of MeV. */
	nupar_1.tnmev = nupar_1.tnu * .086171;
/* ..........COMPUTE OVERFLOW LIMITS FOR LIMITS OF INTEGRATION (Ref 1 & 2). */
/* Convert neutrino temp to units of Me */
	w[0] = -(nupar_1.t9mev / .511) * -88.722;
	w[1] = nupar_1.tnmev / .511 * (modpr_1.xi[0] + 88.029) + 2.531;
	x[0] = nupar_1.t9mev / .511 * 88.029;
	x[1] = -(nupar_1.tnmev / .511) * (modpr_1.xi[0] - 88.722) - 2.531;
	y[0] = -(nupar_1.t9mev / .511) * -88.722;
	y[1] = nupar_1.tnmev / .511 * (88.029 - modpr_1.xi[0]) - 2.531;
	z__[0] = nupar_1.t9mev / .511 * 88.029;
	z__[1] = -(nupar_1.tnmev / .511) * (-88.722 - modpr_1.xi[0]) + 2.531;
/* ..........COMPARE LIMITS AND TAKE LARGER OF THE TWO. */
	uplim1 = abs(w[0]);
	uplim2 = abs(x[0]);
	uplim3 = abs(y[0]);
	uplim4 = abs(z__[0]);
	if (uplim1 < abs(w[1])) {
	    uplim1 = w[1];
	}
	if (uplim2 < abs(x[1])) {
	    uplim2 = x[1];
	}
	if (uplim3 < abs(y[1])) {
	    uplim3 = y[1];
	}
	if (uplim4 < abs(z__[1])) {
	    uplim4 = z__[1];
	}
/* ..........EVALUATE THE INTEGRALS NUMERICALLY. */
	part1 = xintd_(&c_b3, &uplim1, (U_fp)func1_, &cs__50);
	part2 = xintd_(&c_b3, &uplim2, (U_fp)func2_, &cs__50);
	part3 = xintd_(&c_b3, &uplim3, (U_fp)func3_, &cs__50);
	part4 = xintd_(&c_b3, &uplim4, (U_fp)func4_, &cs__50);
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
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    extern doublereal ex_(doublereal *);
    static doublereal t9a, t9b, t912, t913, t923, t932, t943, t953, t9c, t9d, 
	    t9e, t9f, t9m1, t9a32, t9b32, t9c13, t9d13, t9e13, t9f13, t9c56, 
	    t9d56, t9e56, t9f56, t9m23, t9m32;

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
    t913 = pow_dd(&evolp1_1.t9, &c_b12);
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
    t9m23 = 1. / t923;
/* t9**(-2/3) */
    t9m32 = 1. / t932;
/* t9**(-3/2) */
    t9a = evolp1_1.t9 / (evolp1_1.t9 * 13.076 + 1.);
/* For reaction 17. */
    t9a32 = pow_dd(&t9a, &c_b13);
/* t9a**(3/2) */
    t9b = evolp1_1.t9 / (evolp1_1.t9 * 49.18 + 1.);
/* For reaction 18. */
    t9b32 = pow_dd(&t9b, &c_b13);
/* t9b**(3/2) */
    if (evolp1_1.t9 > 10.) {
/* For reaction 22. */
	t9c = 1.;
    } else {
	d__1 = 1. - evolp1_1.t9 * .0969;
	t9c = evolp1_1.t9 / (1. - evolp1_1.t9 * .0969 + t953 * .0284 / pow_dd(
		&d__1, &c_b15));
    }
    t9c13 = pow_dd(&t9c, &c_b16);
/* t9c**(1/3) */
    t9c56 = pow_dd(&t9c, &c_b17);
/* t9c**(5/6) */
    t9d = evolp1_1.t9 / (evolp1_1.t9 * .759 + 1.);
/* For reaction 24. */
    t9d13 = pow_dd(&t9d, &c_b16);
/* t9d**(1/3) */
    t9d56 = pow_dd(&t9d, &c_b17);
/* t9d**(5/6) */
    t9e = evolp1_1.t9 / (evolp1_1.t9 * .1378 + 1.);
/* For reaction 26. */
    t9e13 = pow_dd(&t9e, &c_b16);
/* t9e**(1/3) */
    t9e56 = pow_dd(&t9e, &c_b17);
/* t9e**(5/6) */
    t9f = evolp1_1.t9 / (evolp1_1.t9 * .1071 + 1.);
/* For reaction 27. */
    t9f13 = pow_dd(&t9f, &c_b16);
/* t9f**(1/3) */
    t9f56 = pow_dd(&t9f, &c_b17);
/* 20--------NEUTRON, PHOTON REACTIONS-------------------------------------------- */
/* .......H(n,g)H2...................(Smith-Kawano-Malaney 1992) */
/* t9f**(5/6) */
    rates_1.f[11] = (1. - t912 * .8504 + evolp1_1.t9 * .4895 - t932 * .09623 
	    + evolp1_1.t9 * .008471 * evolp1_1.t9 - evolp1_1.t9 * 2.8e-4 * 
	    t932) * 47420.;
/* .......H2(n,g)H3..................(Wagoner 1969) */
    rates_1.f[12] = (evolp1_1.t9 * 18.9 + 1.) * 66.2;
/* .......He3(n,g)He4................(Wagoner 1969) */
    rates_1.f[13] = (evolp1_1.t9 * 905. + 1.) * 6.62;
/* .......Li6(n,g)Li7................(Malaney-Fowler 1989) */
    rates_1.f[14] = 5100.;
/* 30--------NEUTRON, PROTON REACTIONS-------------------------------------------- */
/* .......He3(n,p)H3.................(Smith-Kawano-Malaney 1992) */
    rates_1.f[15] = (1. - t912 * .508 + evolp1_1.t9 * .228) * 7.21e8;
/* .......Be7(n,p)Li7................(Smith-Kawano-Malaney 1992) */
    d__1 = -.07486 / evolp1_1.t9;
    rates_1.f[16] = (1. - t912 * .56 + evolp1_1.t9 * .179 - t932 * .0283 + 
	    evolp1_1.t9 * .002214 * evolp1_1.t9 - evolp1_1.t9 * 6.851e-5 * 
	    t932) * 2.675e9 + t9a32 * 9.391e8 * t9m32 + t9m32 * 4.467e7 * ex_(
	    &d__1);
/* 40--------NEUTRON, ALPHA REACTIONS--------------------------------------------- */
/* .......Li6(n,a)H3.................(Caughlan-Fowler 1988) */
    d__1 = -2.39 / evolp1_1.t9;
    rates_1.f[17] = t9m32 * 2.54e9 * ex_(&d__1) + (1. - t9b32 * .261 / t932) *
	     1.68e8;
/* .......Be7(n,a)He4................(Wagoner 1969) */
    rates_1.f[18] = (evolp1_1.t9 * 3760. + 1.) * 20500.;
/* 50--------PROTON, PHOTON REACTIONS--------------------------------------------- */
/* .......H2(p,g)He3.................(Smith-Kawano-Malaney 1992) */
    d__1 = -3.72 / t913;
    rates_1.f[19] = t9m23 * 2650. * ex_(&d__1) * (t913 * .112 + 1. + t923 * 
	    1.99 + evolp1_1.t9 * 1.56 + t943 * .162 + t953 * .324);
/* .......H3(p,g)He4.................(Caughlan-Fowler 1988) */
    d__1 = -3.869 / t913;
    rates_1.f[20] = t9m23 * 2.2e4 * ex_(&d__1) * (t913 * .108 + 1. + t923 * 
	    1.68 + evolp1_1.t9 * 1.26 + t943 * .551 + t953 * 1.06);
/* .......Li6(p,g)Be7................(Caughlan-Fowler 1988) */
    d__1 = -8.413 / t9c13;
    rates_1.f[21] = t9c56 * 6.69e5 * t9m32 * ex_(&d__1);
/* 60--------PROTON, ALPHA REACTIONS---------------------------------------------- */
/* .......Li6(p,a)He3................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 5.5;
    d__1 = -8.413 / t913 - d__2 * d__2;
    d__3 = -17.763 / evolp1_1.t9;
    d__4 = -21.82 / evolp1_1.t9;
    rates_1.f[22] = t9m23 * 3.73e10 * ex_(&d__1) * (t913 * .05 + 1. - t923 * 
	    .061 - evolp1_1.t9 * .021 + t943 * .006 + t953 * .005) + t9m32 * 
	    1.33e10 * ex_(&d__3) + t9m1 * 1.29e9 * ex_(&d__4);
/* .......Li7(p,a)He4................(Smith-Kawano-Malaney 1992) */
    d__1 = -8.472 / t913;
    d__2 = -8.472 / t9d13;
    d__3 = -30.442 / evolp1_1.t9;
/* Computing 2nd power */
    d__5 = evolp1_1.t9 / 1.696;
    d__4 = -8.472 / t913 - d__5 * d__5;
    d__6 = -4.478 / evolp1_1.t9;
    rates_1.f[23] = t9m23 * 1.096e9 * ex_(&d__1) - t9d56 * 4.83e8 * t9m32 * 
	    ex_(&d__2) + t9m32 * 1.06e10 * ex_(&d__3) + t9m23 * 1.56e5 * ex_(&
	    d__4) * (t913 * .049 + 1. - t923 * 2.498 + evolp1_1.t9 * .86 + 
	    t943 * 3.518 + t953 * 3.08) + t9m32 * 1.55e6 * ex_(&d__6);
/* 70--------ALPHA, PHOTON REACTIONS---------------------------------------------- */
/* .......H2(a,g)Li6.................(Caughlan-Fowler 1988) */
    d__1 = -7.423 / t913;
    d__2 = -8.228 / evolp1_1.t9;
    rates_1.f[24] = t9m23 * 30.1 * ex_(&d__1) * (t913 * .056 + 1. - t923 * 
	    4.85 + evolp1_1.t9 * 8.85 - t943 * .585 - t953 * .584) + t9m32 * 
	    85.5 * ex_(&d__2);
/* .......H3(a,g)Li7.................(Smith-Kawano-Malaney 1992) */
    d__1 = -8.09 / t913;
    d__2 = -8.068 / t9e13;
    rates_1.f[25] = t9m23 * 303200. * ex_(&d__1) * (t913 * .0516 + 1. + t923 *
	     .0229 + evolp1_1.t9 * .00828 - t943 * 3.28e-4 - t953 * 3.01e-4) 
	    + t9e56 * 510900. * t9m32 * ex_(&d__2);
/* .......He3(a,g)Be7................(Smith-Kawano-Malaney 1992) */
    d__1 = -14.964 / t913;
    d__2 = -12.859 / t9f13;
    rates_1.f[26] = t9m23 * 4.817e6 * ex_(&d__1) * (t913 * .0325 + 1. - t923 *
	     .00104 - evolp1_1.t9 * 2.37e-4 - t943 * 8.11e-5 - t953 * 4.69e-5)
	     + t9f56 * 5.938e6 * t9m32 * ex_(&d__2);
/* 80--------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS------------------- */
/* .......H2(d,n)He3.................(Smith-Kawano-Malaney 1992) */
    d__1 = -4.259 / t913;
    rates_1.f[27] = t9m23 * 3.95e8 * ex_(&d__1) * (t913 * .098 + 1. + t923 * 
	    .765 + evolp1_1.t9 * .525 + t943 * .00961 + t953 * .0167);
/* .......H2(d,p)H3..................(Smith-Kawano-Malaney 1992) */
    d__1 = -4.258 / t913;
    rates_1.f[28] = t9m23 * 4.17e8 * ex_(&d__1) * (t913 * .098 + 1. + t923 * 
	    .518 + evolp1_1.t9 * .355 - t943 * .01 - t953 * .018);
/* .......H3(d,n)He4.................(Smith-Kawano-Malaney 1992) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .0754;
    d__1 = -4.559 / t913 - d__2 * d__2;
    d__3 = -.4857 / evolp1_1.t9;
    rates_1.f[29] = t9m23 * 1.063e11 * ex_(&d__1) * (t913 * .092 + 1. - t923 *
	     .375 - evolp1_1.t9 * .242 + t943 * 33.82 + t953 * 55.42) + t9m23 
	    * 8.047e8 * ex_(&d__3);
/* .......He3(d,p)He4................(Smith-Kawano-Malaney 1992) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .27;
    d__1 = -7.144 / t913 - d__2 * d__2;
    d__3 = -1.762 / evolp1_1.t9;
    rates_1.f[30] = t9m23 * 5.021e10 * ex_(&d__1) * (t913 * .058 + 1. + t923 *
	     .603 + evolp1_1.t9 * .245 + t943 * 6.97 + t953 * 7.19) + 5.212e8 
	    / t912 * ex_(&d__3);
/* 90--------THREE PARTICLE REACTIONS--------------------------------------------- */
/* .......He3(He3,2p)He4.............(Caughlan-Fowler 1988) */
    d__1 = -12.276 / t913;
    rates_1.f[31] = t9m23 * 6.04e10 * ex_(&d__1) * (t913 * .034 + 1. - t923 * 
	    .522 - evolp1_1.t9 * .124 + t943 * .353 + t953 * .213);
/* .......Li7(d,na)He4...............(Caughlan-Fowler 1988) */
    d__1 = -10.259 / t913;
    rates_1.f[32] = t9m23 * 2.92e11 * ex_(&d__1);
/* .......Be7(d,pa)He4...............(Caughlan-Fowler 1988) */
    d__1 = -12.428 / t913;
    rates_1.f[33] = t9m23 * 1.07e12 * ex_(&d__1);
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
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    extern doublereal ex_(doublereal *);
    static doublereal t9a, t912, t913, t923, t932, t943, t953, t915, t954, 
	    t9m1, t9a13, t9a56, t9m23, t9m32, t9m34, t9m15, t9m54;

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
    t913 = pow_dd(&evolp1_1.t9, &c_b12);
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
    t915 = pow_dd(&evolp1_1.t9, &c_b26);
/* t9**(1/5) */
    t954 = pow_dd(&evolp1_1.t9, &c_b27);
/* t9**(5/4) */
    t9m1 = 1. / evolp1_1.t9;
/* t9**(-1) */
    t9m23 = 1. / t923;
/* t9**(-2/3) */
    t9m32 = 1. / t932;
/* t9**(-3/2) */
    t9m34 = sqrt(t9m32);
/* t9**(-3/4) */
    t9m15 = 1. / t915;
/* t9**(-1/5) */
    t9m54 = 1. / t954;
/* t9**(-5/4) */
    t9a = evolp1_1.t9 / (evolp1_1.t9 / 15.1 + 1.);
/* For reaction 53. */
    t9a13 = pow_dd(&t9a, &c_b16);
/* t9a**(1/3) */
    t9a56 = pow_dd(&t9a, &c_b17);
/* 20--------NEUTRON, PHOTON REACTIONS-------------------------------------------- */
/* .......Li7(n,g)Li8................(Wagoner 1969) */
/* t9a**(5/6) */
    d__1 = -2.62 / evolp1_1.t9;
    rates_1.f[34] = t9m32 * 9960. * ex_(&d__1) + 4900.;
/* .......B10(n,g)B11................(Wagoner 1969) */
    rates_1.f[35] = 66200.;
/* .......B11(n,g)B12................(Malaney-Fowler 1989) */
    d__1 = -.223 / evolp1_1.t9;
    rates_1.f[36] = t9m32 * 2400. * ex_(&d__1) + 729.;
/* 30--------NEUTRON, PROTON REACTIONS-------------------------------------------- */
/* .......C11(n,p)B11................(Caughlan-Fowler 1988) */
    rates_1.f[37] = (1. - t912 * .048 + evolp1_1.t9 * .01) * 1.69e8;
/* 40--------NEUTRON, ALPHA REACTIONS--------------------------------------------- */
/* .......B10(n,a)Li7................(Caughlan-Fowler 1988) */
    rates_1.f[38] = 5.07e8;
/* 50--------PROTON, PHOTON REACTIONS--------------------------------------------- */
/* .......Be7(p,g)B8.................(Caughlan-Fowler 1988) */
    d__1 = -10.262 / t913;
    d__2 = -7.306 / evolp1_1.t9;
    rates_1.f[39] = t9m23 * 3.11e5 * ex_(&d__1) + t9m32 * 2530. * ex_(&d__2);
/* .......Be9(p,g)B10................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .846;
    d__1 = -10.359 / t913 - d__2 * d__2;
    d__3 = -3.445 / evolp1_1.t9;
    d__4 = -10.62 / evolp1_1.t9;
    rates_1.f[40] = t9m23 * 1.33e7 * ex_(&d__1) * (t913 * .04 + 1. + t923 * 
	    1.52 + evolp1_1.t9 * .428 + t943 * 2.15 + t953 * 1.54) + t9m32 * 
	    96400. * ex_(&d__3) + t9m32 * 2.72e6 * ex_(&d__4);
/* .......B10(p,g)C11................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 4.402;
    d__1 = -12.062 / t913 - d__2 * d__2;
    d__3 = -12.041 / evolp1_1.t9;
    d__4 = -16.164 / evolp1_1.t9;
    rates_1.f[41] = t9m23 * 4.61e5 * ex_(&d__1) * (t913 * .035 + 1. + t923 * 
	    .426 + evolp1_1.t9 * .103 + t943 * .281 + t953 * .173) + t9m32 * 
	    1.93e5 * ex_(&d__3) + t9m32 * 11400. * ex_(&d__4);
/* .......B11(p,g)C12................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .239;
    d__1 = -12.095 / t913 - d__2 * d__2;
    d__3 = -1.733 / evolp1_1.t9;
    d__4 = -5.617 / evolp1_1.t9;
    rates_1.f[42] = t9m23 * 4.62e7 * ex_(&d__1) * (t913 * .035 + 1. + t923 * 
	    3. + evolp1_1.t9 * .723 + t943 * 9.91 + t953 * 6.07) + t9m32 * 
	    7890. * ex_(&d__3) + t9m15 * 96800. * ex_(&d__4);
/* .......C11(p,g)N12................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 1.627;
    d__1 = -13.658 / t913 - d__2 * d__2;
    d__3 = -7.021 / evolp1_1.t9;
    rates_1.f[43] = t9m23 * 42400. * ex_(&d__1) * (t913 * .031 + 1. + t923 * 
	    3.11 + evolp1_1.t9 * .665 + t943 * 4.61 + t953 * 2.5) + t9m32 * 
	    8840. * ex_(&d__3);
/* 60--------PROTON, NEUTRON REACTIONS-------------------------------------------- */
/* .......B12(p,n)C12................(Wagoner 1969) */
    d__1 = -12.12 / t913;
    rates_1.f[44] = t9m23 * 4.02e11 * ex_(&d__1);
/* 70--------PROTON, ALPHA REACTIONS---------------------------------------------- */
/* .......Be9(p,a)Li6................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .52;
    d__1 = -10.359 / t913 - d__2 * d__2;
    d__3 = -3.046 / evolp1_1.t9;
    d__4 = -5.16 / evolp1_1.t9;
    rates_1.f[45] = t9m23 * 2.11e11 * ex_(&d__1) * (t913 * .04 + 1. + t923 * 
	    1.09 + evolp1_1.t9 * .307 + t943 * 3.21 + t953 * 2.3) + t9m1 * 
	    4.51e8 * ex_(&d__3) + t9m34 * 6.7e8 * ex_(&d__4);
/* .......B10(p,a)Be7................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 4.402;
    d__1 = -12.062 / t913 - d__2 * d__2;
    d__3 = -12.26 / evolp1_1.t9;
    rates_1.f[46] = t9m23 * 1.26e11 * ex_(&d__1) * (t913 * .035 + 1. - t923 * 
	    .498 - evolp1_1.t9 * .121 + t943 * .3 + t953 * .184) + t9m1 * 
	    2.59e9 * ex_(&d__3);
/* .......B12(p,a)Be9................(Wagoner 1969) */
    d__1 = -12.12 / t913;
    rates_1.f[47] = t9m23 * 2.01e11 * ex_(&d__1);
/* 80--------ALPHA, PHOTON REACTIONS---------------------------------------------- */
/* .......Li6(a,g)B10................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 1.326;
    d__1 = -18.79 / t913 - d__2 * d__2;
    d__3 = -3.484 / evolp1_1.t9;
    d__4 = -7.269 / evolp1_1.t9;
    rates_1.f[48] = t9m23 * 4.06e6 * ex_(&d__1) * (t913 * .022 + 1. + t923 * 
	    1.54 + evolp1_1.t9 * .239 + t943 * 2.2 + t953 * .869) + t9m32 * 
	    1910. * ex_(&d__3) + t9m1 * 10100. * ex_(&d__4);
/* .......Li7(a,g)B11................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 4.195;
    d__1 = -19.161 / t913 - d__2 * d__2;
    d__3 = -2.977 / evolp1_1.t9;
    d__4 = -6.227 / evolp1_1.t9;
    rates_1.f[49] = t9m23 * 3.55e7 * ex_(&d__1) * (t913 * .022 + 1. + t923 * 
	    .775 + evolp1_1.t9 * .118 + t943 * .884 + t953 * .342) + t9m32 * 
	    333. * ex_(&d__3) + t9m1 * 4.1e4 * ex_(&d__4);
/* .......Be7(a,g)C11................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 4.769;
    d__1 = -23.212 / t913 - d__2 * d__2;
    d__3 = -6.51 / evolp1_1.t9;
    d__4 = -10.039 / evolp1_1.t9;
    rates_1.f[50] = t9m23 * 8.45e7 * ex_(&d__1) * (t913 * .018 + 1. + t923 * 
	    .488 + evolp1_1.t9 * .061 + t943 * .296 + t953 * .095) + t9m32 * 
	    12500. * ex_(&d__3) + t9m54 * 1.29e5 * ex_(&d__4);
/* 90--------ALPHA, PROTON REACTIONS---------------------------------------------- */
/* .......B8(a,p)C11.................(Wagoner 1969) */
    d__1 = -27.36 / t913;
    rates_1.f[51] = t9m23 * 1.08e15 * ex_(&d__1);
/* 100-------ALPHA, NEUTRON REACTIONS--------------------------------------------- */
/* .......Li8(a,n)B11................(Malaney-Fowler 1989) */
    d__1 = -19.461 / t9a13;
    rates_1.f[52] = t9a56 * 8.62e13 * t9m32 * ex_(&d__1);
/* .......Be9(a,n)C12................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .049;
    d__1 = -23.87 / t913 - d__2 * d__2;
    d__3 = -1.184 / evolp1_1.t9;
    d__4 = -1.834 / evolp1_1.t9;
    d__5 = -4.179 / evolp1_1.t9;
    d__6 = -12.732 / evolp1_1.t9;
    rates_1.f[53] = t9m23 * 4.62e13 * ex_(&d__1) * (t913 * .017 + 1. + t923 * 
	    8.57 + evolp1_1.t9 * 1.05 + t943 * 74.51 + t953 * 23.15) + t9m32 *
	     7.34e-5 * ex_(&d__3) + t9m32 * .227 * ex_(&d__4) + t9m32 * 
	    1.26e5 * ex_(&d__5) + ex_(&d__6) * 2.4e8;
/* 110-------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS------------------- */
/* .......Be9(d,n)B10................(original Wagoner code) */
    d__1 = 6.44 - 12.6 / t913;
    rates_1.f[54] = t9m23 * 7.16e8 * ex_(&d__1);
/* .......B10(d,p)B11................(original Wagoner code) */
    d__1 = 7.3 - 14.8 / t913;
    rates_1.f[55] = t9m23 * 9.53e8 * ex_(&d__1);
/* .......B11(d,n)C12................(original Wagoner code) */
    d__1 = 7.4 - 14.8 / t913;
    rates_1.f[56] = t9m23 * 1.41e9 * ex_(&d__1);
/* 120-------THREE PARTICLE REACTIONS--------------------------------------------- */
/* .......He4(an,g)Be9...............(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__1 = evolp1_1.t9;
    d__2 = -1.062 / evolp1_1.t9;
    rates_1.f[57] = 2.59e-6 / ((evolp1_1.t9 * .344 + 1.) * (d__1 * d__1)) * 
	    ex_(&d__2);
/* .......He4(2a,g)C12...............(Caughlan-Fowler 1988) */
    d__1 = -4.4027 / evolp1_1.t9;
    d__2 = -24.811 / evolp1_1.t9;
    rates_1.f[58] = t9m32 * 2.79e-8 * t9m32 * ex_(&d__1) + t9m32 * 1.35e-8 * 
	    ex_(&d__2);
/* .......Li8(p,na)He4...............(original Wagoner code) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 2.53;
    d__1 = -8.52 / t913 - d__2 * d__2;
    d__3 = -4.64 / evolp1_1.t9;
    rates_1.f[59] = t9m23 * 8.65e9 * ex_(&d__1) + t9m32 * 2.31e9 * ex_(&d__3);
/* .......B8(n,pa)He4................(original Wagoner code) */
    rates_1.f[60] = 4.02e8;
/* .......Be9(p,da)He4...............(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .52;
    d__1 = -10.359 / t913 - d__2 * d__2;
    d__3 = -3.046 / evolp1_1.t9;
    d__4 = -5.8 / evolp1_1.t9;
    rates_1.f[61] = t9m23 * 2.11e11 * ex_(&d__1) * (t913 * .04 + 1. + t923 * 
	    1.09 + evolp1_1.t9 * .307 + t943 * 3.21 + t953 * 2.3) + t9m1 * 
	    5.79e8 * ex_(&d__3) + t9m34 * 8.5e8 * ex_(&d__4);
/* .......B11(p,2a)He4...............(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 1.644;
    d__1 = -12.095 / t913 - d__2 * d__2;
    d__3 = -1.734 / evolp1_1.t9;
    d__4 = -6.262 / evolp1_1.t9;
    d__5 = -14.154 / evolp1_1.t9;
    rates_1.f[62] = t9m23 * 2.2e12 * ex_(&d__1) * (t913 * .034 + 1. + t923 * 
	    .14 + evolp1_1.t9 * .034 + t943 * .19 + t953 * .116) + t9m32 * 
	    4.03e6 * ex_(&d__3) + t9m32 * 6.73e9 * ex_(&d__4) + t9m1 * 3.88e9 
	    * ex_(&d__5);
/* .......C11(n,2a)He4...............(Wagoner 1969) */
    rates_1.f[63] = 1.58e8;
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
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    extern doublereal ex_(doublereal *);
    static doublereal t9a, t9b, t912, t913, t923, t932, t943, t953, t935, 
	    t965, t938, t9a13, t9b13, t9a56, t9b56, t9m13, t9m23, t9m32, 
	    t9m65;

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
    t913 = pow_dd(&evolp1_1.t9, &c_b12);
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
    t935 = pow_dd(&evolp1_1.t9, &c_b32);
/* t9**(3/5) */
    t965 = pow_dd(&evolp1_1.t9, &c_b33);
/* t9**(6/5) */
    t938 = pow_dd(&evolp1_1.t9, &c_b34);
/* t9**(3/8) */
    t9m13 = 1. / t913;
/* t9**(1/3) */
    t9m23 = 1. / t923;
/* t9**(-2/3) */
    t9m32 = 1. / t932;
/* t9**(-3/2) */
    t9m65 = 1. / t965;
/* t9**(-6/5) */
    d__1 = evolp1_1.t9 * .0478 + 1.;
    t9a = evolp1_1.t9 / (evolp1_1.t9 * .0478 + 1. + t953 * .00756 / pow_dd(&
	    d__1, &c_b15));
/* For reaction 82. */
    t9a13 = pow_dd(&t9a, &c_b12);
/* t9a**(1/3) */
    t9a56 = pow_dd(&t9a, &c_b37);
/* t9a**(5/6) */
    d__1 = evolp1_1.t9 * .0776 + 1.;
    t9b = evolp1_1.t9 / (evolp1_1.t9 * .0776 + 1. + t953 * .0264 / pow_dd(&
	    d__1, &c_b15));
/* For reaction 84. */
    t9b13 = pow_dd(&t9b, &c_b12);
/* t9b**(1/3) */
    t9b56 = pow_dd(&t9b, &c_b37);
/* 20--------NEUTRON, PHOTON REACTIONS-------------------------------------------- */
/* .......C12(n,g)C13................(Wagoner 1969) */
/* t9b**(5/6) */
    rates_1.f[64] = 450.;
/* .......C13(n,g)C14................(Wagoner 1969) */
    d__1 = -1.67 / evolp1_1.t9;
    rates_1.f[65] = t9m32 * 2.38e5 * ex_(&d__1) + 119.;
/* .......N14(n,g)N15................(Wagoner 1969) */
    rates_1.f[66] = 9940.;
/* 30--------NEUTRON, PROTON REACTIONS-------------------------------------------- */
/* .......N13(n,p)C13................(Caughlan-Fowler 1988) */
    rates_1.f[67] = (1. - t912 * .167 + evolp1_1.t9 * .037) * 1.88e8;
/* .......N14(n,p)C14................(Caughlan-Fowler 1988) */
    d__1 = -4.983 / evolp1_1.t9;
    rates_1.f[68] = (t912 * .361 + 1. + evolp1_1.t9 * .502) * 2.39e5 + 
	    1.112e8 / t912 * ex_(&d__1);
/* .......O15(n,p)N15................(Caughlan-Fowler 1988) */
    rates_1.f[69] = (t912 * .452 + 1. - evolp1_1.t9 * .191) * 3.5e8;
/* 40--------NEUTRON, ALPHA REACTIONS--------------------------------------------- */
/* .......O15(n,a)C12................(Caughlan-Fowler 1988) */
    rates_1.f[70] = (t912 * .188 + 1. + evolp1_1.t9 * .015) * 3.5e7;
/* 50--------PROTON, PHOTON REACTIONS--------------------------------------------- */
/* .......C12(p,g)N13................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 1.5;
    d__1 = -13.69 / t913 - d__2 * d__2;
    d__3 = -4.925 / evolp1_1.t9;
    d__4 = -18.179 / evolp1_1.t9;
    rates_1.f[71] = t9m23 * 2.04e7 * ex_(&d__1) * (t913 * .03 + 1. + t923 * 
	    1.19 + evolp1_1.t9 * .254 + t943 * 2.06 + t953 * 1.12) + t9m32 * 
	    1.08e5 * ex_(&d__3) + t9m32 * 2.15e5 * ex_(&d__4);
/* .......C13(p,g)N14................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 2.;
    d__1 = -13.717 / t913 - d__2 * d__2;
    d__3 = -5.701 / evolp1_1.t9;
    rates_1.f[72] = t9m23 * 8.01e7 * ex_(&d__1) * (t913 * .03 + 1. + t923 * 
	    .958 + evolp1_1.t9 * .204 + t943 * 1.39 + t953 * .753) + t9m65 * 
	    1.21e6 * ex_(&d__3);
/* .......C14(p,g)N15................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 5.721;
    d__1 = -13.741 / t913 - d__2 * d__2;
    d__3 = -3.811 / evolp1_1.t9;
    d__4 = -4.739 / evolp1_1.t9;
    rates_1.f[73] = t9m23 * 6.8e6 * ex_(&d__1) * (t913 * .03 + 1. + t923 * 
	    .503 + evolp1_1.t9 * .107 + t943 * .213 + t953 * .115) + t9m32 * 
	    5360. * ex_(&d__3) + t9m13 * 98200. * ex_(&d__4);
/* .......N13(p,g)O14................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 1.191;
    d__1 = -15.202 / t913 - d__2 * d__2;
    d__3 = -6.348 / evolp1_1.t9;
    rates_1.f[74] = t9m23 * 4.04e7 * ex_(&d__1) * (t913 * .027 + 1. - t923 * 
	    .803 - evolp1_1.t9 * .154 + t943 * 5. + t953 * 2.44) + t9m32 * 
	    2.43e5 * ex_(&d__3);
/* .......N14(p,g)O15................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 3.294;
    d__1 = -15.228 / t913 - d__2 * d__2;
    d__3 = -3.011 / evolp1_1.t9;
    d__4 = -12.53 / evolp1_1.t9;
    rates_1.f[75] = t9m23 * 4.9e7 * ex_(&d__1) * (t913 * .027 + 1. - t923 * 
	    .778 - evolp1_1.t9 * .149 + t943 * .261 + t953 * .127) + t9m32 * 
	    2370. * ex_(&d__3) + ex_(&d__4) * 21900.;
/* .......N15(p,g)O16................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .45;
    d__1 = -15.251 / t913 - d__2 * d__2;
    d__3 = -3.328 / evolp1_1.t9;
    d__4 = -4.665 / evolp1_1.t9;
    d__5 = -11.048 / evolp1_1.t9;
    rates_1.f[76] = t9m23 * 9.78e8 * ex_(&d__1) * (t913 * .027 + 1. + t923 * 
	    .219 + evolp1_1.t9 * .042 + t943 * 6.83 + t953 * 3.32) + t9m32 * 
	    11100. * ex_(&d__3) + t9m32 * 14900. * ex_(&d__4) + t9m32 * 3.8e6 
	    * ex_(&d__5);
/* 60--------PROTON, ALPHA REACTIONS---------------------------------------------- */
/* .......N15(p,a)C12................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .522;
    d__1 = -15.251 / t913 - d__2 * d__2;
    d__3 = -3.676 / evolp1_1.t9;
    d__4 = -8.926 / evolp1_1.t9;
    d__5 = -7.721 / evolp1_1.t9;
    d__6 = -11.418 / evolp1_1.t9;
    rates_1.f[77] = t9m23 * 1.08e12 * ex_(&d__1) * (t913 * .027 + 1. + t923 * 
	    2.62 + evolp1_1.t9 * .501 + t943 * 5.36 + t953 * 2.6) + t9m32 * 
	    1.19e8 * ex_(&d__3) + 5.41e8 / t912 * ex_(&d__4) + t9m32 * 4.72e7 
	    * ex_(&d__5) + t9m32 * 2.2e8 * ex_(&d__6);
/* 70--------ALPHA, PHOTON REACTIONS---------------------------------------------- */
/* .......C12(a,g)O16................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__1 = evolp1_1.t9;
/* Computing 2nd power */
    d__3 = evolp1_1.t9 / 3.496;
    d__2 = -32.12 / t913 - d__3 * d__3;
/* Computing 2nd power */
    d__4 = t9m23 * .0489 + 1.;
/* Computing 2nd power */
    d__5 = evolp1_1.t9;
/* Computing 2nd power */
    d__6 = t9m23 * .2654 + 1.;
    d__7 = -32.12 / t913;
    d__8 = -27.499 / evolp1_1.t9;
/* Computing 5th power */
    d__9 = evolp1_1.t9, d__10 = d__9, d__9 *= d__9;
    d__11 = -15.541 / evolp1_1.t9;
    rates_1.f[78] = 1.04e8 / (d__1 * d__1) * ex_(&d__2) / (d__4 * d__4) + 
	    1.76e8 / (d__5 * d__5) / (d__6 * d__6) * ex_(&d__7) + t9m32 * 
	    1250. * ex_(&d__8) + d__10 * (d__9 * d__9) * .0143 * ex_(&d__11);
/* 80--------ALPHA, PROTON REACTIONS---------------------------------------------- */
/* .......B10(a,p)C13................(Wagoner 1969) */
    d__1 = -27.99 / t913;
    rates_1.f[79] = t9m23 * 9.6e14 * ex_(&d__1);
/* .......B11(a,p)C14................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .347;
    d__1 = -28.234 / t913 - d__2 * d__2;
    d__3 = -2.827 / evolp1_1.t9;
    d__4 = -5.178 / evolp1_1.t9;
    d__5 = -11.617 / evolp1_1.t9;
    rates_1.f[80] = t9m23 * 5.37e11 * ex_(&d__1) * (t913 * .015 + 1. + t923 * 
	    5.575 + evolp1_1.t9 * .576 + t943 * 15.888 + t953 * 4.174) + 
	    t9m32 * .00544 * ex_(&d__3) + t9m32 * 336. * ex_(&d__4) + 5.32e6 /
	     t938 * ex_(&d__5);
/* .......C11(a,p)N14................(Caughlan-Fowler 1988) */
    d__1 = -31.883 / t9a13;
    rates_1.f[81] = t9a56 * 7.15e15 * t9m32 * ex_(&d__1);
/* .......N12(a,p)O15................(Caughlan-Fowler 1988) */
    d__1 = -35.6 / t913;
    rates_1.f[82] = t9m23 * 5.59e16 * ex_(&d__1);
/* .......N13(a,p)O16................(Caughlan-Fowler 1988) */
    d__1 = -35.829 / t9b13;
    rates_1.f[83] = t9b56 * 3.23e17 * t9m32 * ex_(&d__1);
/* 90--------ALPHA, NEUTRON REACTIONS--------------------------------------------- */
/* .......B10(a,n)N13................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 9.589;
    d__1 = -27.989 / t913 - d__2 * d__2;
    rates_1.f[84] = t9m23 * 1.2e13 * ex_(&d__1);
/* .......B11(a,n)N14................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / .14;
    d__1 = -28.234 / t913 - d__2 * d__2;
    d__3 = -2.827 / evolp1_1.t9;
    d__4 = -5.178 / evolp1_1.t9;
    d__5 = -8.596 / evolp1_1.t9;
    rates_1.f[85] = t9m23 * 6.97e12 * ex_(&d__1) * (t913 * .015 + 1. + t923 * 
	    8.115 + evolp1_1.t9 * .838 + t943 * 39.804 + t953 * 10.456) + 
	    t9m32 * 1.79 * ex_(&d__3) + t9m32 * 1710. * ex_(&d__4) + t935 * 
	    4.49e6 * ex_(&d__5);
/* .......B12(a,n)N15................(Wagoner 1969) */
    d__1 = -28.45 / t913;
    rates_1.f[86] = t9m23 * 3.04e15 * ex_(&d__1);
/* .......C13(a,n)O16................(Caughlan-Fowler 1988) */
/* Computing 2nd power */
    d__2 = evolp1_1.t9 / 1.284;
    d__1 = -32.329 / t913 - d__2 * d__2;
    d__3 = -9.373 / evolp1_1.t9;
    d__4 = -11.873 / evolp1_1.t9;
    d__5 = -20.409 / evolp1_1.t9;
    d__6 = -29.283 / evolp1_1.t9;
    rates_1.f[87] = t9m23 * 6.77e15 * ex_(&d__1) * (t913 * .013 + 1. + t923 * 
	    2.04 + evolp1_1.t9 * .184) + t9m32 * 3.82e5 * ex_(&d__3) + t9m32 *
	     1.41e6 * ex_(&d__4) + t9m32 * 2e9 * ex_(&d__5) + t9m32 * 2.92e9 *
	     ex_(&d__6);
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

