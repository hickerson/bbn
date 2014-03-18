/* nuccom.f -- translated by f2c (version 20100827).
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
    doublereal t9, hv, phie, y[26];
} evolp1_;

#define evolp1_1 evolp1_

extern struct {
    doublereal dt9, dhv, dphie, dydt[26];
} evolp2_;

#define evolp2_1 evolp2_

extern struct {
    doublereal t90, hv0, phie0, y0[26];
} evolp3_;

#define evolp3_1 evolp3_

extern struct {
    doublereal cy, ct, t9i, t9f, ytmin;
    shortint inc;
} compr_;

#define compr_1 compr_

extern struct {
    doublereal t, dt, dlt9dt;
} ttime_;

#define ttime_1 ttime_

extern struct {
    shortint ltime, is, ip, it, mbad;
} flags_;

#define flags_1 flags_

extern struct {
    shortint itime;
} tcheck_;

#define tcheck_1 tcheck_

extern struct {
    shortint irun, isize, jsize;
} runopt_;

#define runopt_1 runopt_

extern union {
    struct {
	doublereal f[88], r__;
    } _1;
    struct {
	doublereal f[88], r__[88];
    } _2;
} rates_;

#define rates_1 (rates_._1)
#define rates_2 (rates_._2)

extern struct {
    doublereal g, tau, xnu, c__[3], cosmo, xi[3];
} modpr_;

#define modpr_1 modpr_

extern struct {
    doublereal dt1, eta1;
} varpr_;

#define varpr_1 varpr_

extern struct {
    doublereal rhone0, rhob0, rhob, rnb;
} endens_;

#define endens_1 endens_

extern struct {
    doublereal bl1, bl2, bl3, bl4, bl5, bm1, bm2, bm3, bm4, bm5, bn1, bn2, 
	    bn3, bn4, bn5;
} besselcb_;

#define besselcb_1 besselcb_

extern struct {
    doublereal t9mev, tnmev, tnu, cnorm, rhonu;
    shortint nu;
} nupar_;

#define nupar_1 nupar_

extern struct {
    doublereal thm[14], hubcst;
} thermcb_;

#define thermcb_1 thermcb_

extern struct {
    doublereal am[26], zm[26], dm[26];
} nucdat_;

#define nucdat_1 nucdat_

extern struct {
    doublereal xout[1040]	/* was [40][26] */, thmout[240]	/* was [40][6]
	     */, t9out[40], tout[40], dtout[40], etaout[40], hubout[40];
} outdat_;

#define outdat_1 outdat_

extern struct {
    doublereal bk0, bk1, bk2, bk3, bk4;
} kays_;

#define kays_1 kays_

extern struct {
    shortint iform[88], ii[88], jj[88], kk[88], ll[88];
    doublereal rev[88], q9[88];
} recpr_;

#define recpr_1 recpr_

extern union {
    struct {
	doublereal a[676]	/* was [26][26] */, b[26], yx[26];
    } _1;
    struct {
	doublereal a[676]	/* was [26][26] */, b[26], y[26];
    } _2;
} lncoef_;

#define lncoef_1 (lncoef_._1)
#define lncoef_2 (lncoef_._2)

/* Table of constant values */

static doublereal c_b4 = 1e-5;
static doublereal c_b5 = 1.5;
static doublereal c_b12 = .33333333333333331;
static doublereal c_b13 = 1.3333333333333333;
static doublereal c_b17 = 0.;
static shortint cs__50 = 50;
static doublereal c_b24 = .5;
static shortint cs__0 = 0;
static integer c__1 = 1;

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
    shortint s__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_lg10(doublereal *);

    /* Local variables */
    static shortint i__;
#define v ((doublereal *)&evolp1_1)
#define v0 ((doublereal *)&evolp3_1)
    static doublereal dtl;
#define dvdt ((doublereal *)&evolp2_1)
    static shortint mvar, loop;
    static doublereal dvdt0[29];
    extern /* Subroutine */ int check_(void), accum_(void);
    static doublereal dtmin;
    extern /* Subroutine */ int start_(void), derivs_(shortint *);

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
    if (evolp1_1.t9 <= compr_1.t9f || ttime_1.dt < (d__1 = 1e-16 / 
	    ttime_1.dlt9dt, abs(d__1)) || flags_1.ip == compr_1.inc) {
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
	dtmin = (d__1 = 1. / ttime_1.dlt9dt, abs(d__1)) * compr_1.ct;
/* Trial value for minimum time step (R */
	s__1 = runopt_1.isize;
	for (i__ = 1; i__ <= s__1; ++i__) {
/* Go through all abundance changes. */
	    if (evolp2_1.dydt[i__ - 1] != 0. && evolp1_1.y[i__ - 1] > 
		    compr_1.ytmin) {
/* Computing 2nd power */
		d__2 = d_lg10(&evolp1_1.y[i__ - 1]) / d_lg10(&compr_1.ytmin);
		dtl = (d__1 = evolp1_1.y[i__ - 1] / evolp2_1.dydt[i__ - 1], 
			abs(d__1)) * compr_1.cy * (d__2 * d__2 + 1.);
/* (Ref 2). */
		if (dtl < dtmin) {
		    dtmin = dtl;
		}
/* Find smallest time st */
	    }
	}
	if (dtmin > ttime_1.dt * 1.5) {
	    dtmin = ttime_1.dt * 1.5;
	}
/* Limit change in time */
	ttime_1.dt = dtmin;
/* Set new time step. */
    }
    ttime_1.t += ttime_1.dt;
/* ..........STORE AND INCREMENT VALUES (Ref 3). */
/* Increment time. */
    s__1 = mvar;
    for (i__ = 1; i__ <= s__1; ++i__) {
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
    s__1 = mvar;
    for (i__ = 1; i__ <= s__1; ++i__) {
	v[i__ - 1] = v0[i__ - 1] + (dvdt[i__ - 1] + dvdt0[i__ - 1]) * .5 * 
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
    shortint s__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static shortint i__;
    static doublereal z__;
    extern doublereal ex_(doublereal *);
    extern /* Subroutine */ int rate0_(void), rate1_(doublereal *), bessel_(
	    doublereal *);

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
/* ----------EVOLUTION PARAMETERS. */
/* Forward reaction rate coefficients. */
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
/* ----------VARIATIONAL PARAMETERS. */
/* Neutrino degeneracy parameters. */
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
    d__1 = evolp1_1.t9 * .09615;
    ttime_1.t = 1 / (d__1 * d__1);
/* Initial time (Ref 1). */
    ttime_1.dt = varpr_1.dt1;
/* ..........MODEL SETTINGS. */
/* Initial time step. */
    modpr_1.g = modpr_1.c__[0] * 6.67e-8;
/* Modify gravitational constant. */
    modpr_1.tau = modpr_1.c__[1];
/* Convert n half-life (min) to lifetim */
    modpr_1.tau /= .98;
/* Coulomb correction (Ref 2). */
    modpr_1.xnu = modpr_1.c__[2];
/* 30--------COMPUTE INITIAL ABUNDANCES FOR NEUTRON AND PROTON-------------------- */
/* Number of neutrino species. */
    if (15.011 / evolp1_1.t9 + modpr_1.xi[0] > 58.) {
/* Overabundance of antineut */
	evolp1_1.y[0] = 1e-25;
/* Very little of neutrons. */
	evolp1_1.y[1] = 1.;
/* Essentially all protons. */
    } else {
	if (15.011 / evolp1_1.t9 + modpr_1.xi[0] < -58.) {
/* Overabundance of neutrino */
	    evolp1_1.y[0] = 1.;
/* Essentially all neutrons. */
	    evolp1_1.y[1] = 1e-25;
/* Very little of protons. */
	} else {
	    d__1 = 15.011 / evolp1_1.t9 + modpr_1.xi[0];
	    evolp1_1.y[0] = 1. / (ex_(&d__1) + 1.);
/* Initial n abundance (Ref */
	    d__1 = -15.011 / evolp1_1.t9 - modpr_1.xi[0];
	    evolp1_1.y[1] = 1. / (ex_(&d__1) + 1.);
/* Initial p abundance (Ref */
	}
    }
    if (modpr_1.xi[0] != 0.) {
/* Electron neutrino degeneracy. */
	nupar_1.cnorm = 1.;
	nupar_1.tnu = 1e-5;
/* Low temperature. */
	rate1_(&c_b4);
/* Find normalization constant at low t */
	nupar_1.cnorm = 1 / modpr_1.tau / rates_1.f[0];
    }
    evolp3_1.y0[0] = evolp1_1.y[0];
    evolp3_1.y0[1] = evolp1_1.y[1];
/* 40--------FIND RATIO OF BARYON DENSITY TO TEMPERATURE CUBED-------------------- */
    z__ = 5.93 / evolp1_1.t9;
/* Inverse of temperature. */
    bessel_(&z__);
    evolp1_1.hv = varpr_1.eta1 * 33683. * 2.75;
/* (Ref 4 but with final eta). */
/* Computing 3rd power */
    d__1 = z__;
    evolp1_1.phie = evolp1_1.hv * (evolp1_1.y[1] * 1.784e-5) / (d__1 * (d__1 *
	     d__1) * .5 * (besselcb_1.bl1 - besselcb_1.bl2 * 2. + 
	    besselcb_1.bl3 * 3. - besselcb_1.bl4 * 4. + besselcb_1.bl5 * 5.));
/* Chemical potential of electron (Ref */
/* Computing 3rd power */
    d__1 = evolp1_1.t9;
    endens_1.rhob0 = evolp1_1.hv * (d__1 * (d__1 * d__1));
/* Baryon density. */
    if (modpr_1.xi[0] == 0. && modpr_1.xi[1] == 0. && modpr_1.xi[2] == 0.) {
/* Nonde */
/* Computing 4th power */
	d__1 = evolp1_1.t9, d__1 *= d__1;
	endens_1.rhone0 = d__1 * d__1 * 7.366;
/* Electron neutrino density (Ref 6). */
    }
/* 50--------SET ABUNDANCES FOR REST OF NUCLIDES---------------------------------- */
    d__1 = 25.82 / evolp1_1.t9;
    evolp1_1.y[2] = evolp1_1.y[0] * evolp1_1.y[1] * endens_1.rhob0 * ex_(&
	    d__1) / (pow_dd(&evolp1_1.t9, &c_b5) * 4.71e9);
/* (Ref 7). */
    evolp3_1.y0[2] = evolp1_1.y[2];
    s__1 = runopt_1.isize;
    for (i__ = 4; i__ <= s__1; ++i__) {
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
/* Subroutine */ int derivs_(shortint *loop)
{
    /* System generated locals */
    shortint s__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static shortint i__;
    static doublereal bar;
    extern /* Subroutine */ int sol_(shortint *);
    static doublereal sumy;
    extern /* Subroutine */ int rate1_(doublereal *), rate2_(void), rate3_(
	    void), rate4_(void), therm_(void);
    static doublereal sumdy, sumzy, dphdt9, dlndt9, dphdln, dphdzy, summdy, 
	    sumzdy;

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
    thermcb_1.hubcst = sqrt(modpr_1.g * 8.3775813333333318 * thermcb_1.thm[9] 
	    + modpr_1.cosmo / 3.);
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
    sumy = 0.;
    sumzy = 0.;
    sumdy = 0.;
    summdy = 0.;
    sumzdy = 0.;
/* ..........ACCUMULATE TO GET SUM. */
    s__1 = runopt_1.isize;
    for (i__ = 1; i__ <= s__1; ++i__) {
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
    dphdt9 = thermcb_1.thm[11] * (evolp1_1.hv * -1.07e-4 * sumzy / 
	    evolp1_1.t9 - thermcb_1.thm[10]);
    dphdln = -thermcb_1.thm[11] * 3.568e-5 * evolp1_1.hv * sumzy;
    dphdzy = thermcb_1.thm[11] * 3.568e-5 * evolp1_1.hv;
    bar = evolp1_1.t9 * 9.25e-5 * sumy + evolp1_1.t9 * 1.388e-4 * sumdy / (
	    thermcb_1.hubcst * 3.) + summdy / (thermcb_1.hubcst * 3.);
    dlndt9 = -(thermcb_1.thm[1] + thermcb_1.thm[4] + thermcb_1.thm[5] * 
	    dphdt9 + thermcb_1.thm[8] * 1.388e-4 * sumy) / (thermcb_1.thm[0] 
	    + thermcb_1.thm[2] + thermcb_1.thm[3] + thermcb_1.thm[6] + 
	    thermcb_1.thm[8] * bar + thermcb_1.thm[5] * (dphdln + dphdzy * 
	    sumzdy / (thermcb_1.hubcst * 3.)));
/* (Ref 1) */
    evolp2_1.dt9 = thermcb_1.hubcst * 3. / dlndt9;
    ttime_1.dlt9dt = evolp2_1.dt9 / evolp1_1.t9;
    evolp2_1.dhv = -evolp1_1.hv * (thermcb_1.hubcst * 3. + ttime_1.dlt9dt * 
	    3.);
/* (Ref 2) */
    evolp2_1.dphie = dphdt9 * evolp2_1.dt9 + dphdln * (thermcb_1.hubcst * 3.) 
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
    shortint s__1;

    /* Local variables */
    static shortint i__;


#define xout_ref(a_1,a_2) outdat_1.xout[(a_2)*40 + a_1 - 41]
#define thmout_ref(a_1,a_2) outdat_1.thmout[(a_2)*40 + a_1 - 41]

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
    s__1 = runopt_1.isize;
    for (i__ = 1; i__ <= s__1; ++i__) {
	xout_ref(flags_1.it, i__) = evolp1_1.y[i__ - 1] / evolp1_1.y[1];
    }
    xout_ref(flags_1.it, 2) = evolp1_1.y[1] * nucdat_1.am[1];
/* Exception for proton. */
    xout_ref(flags_1.it, 6) = evolp1_1.y[5] * nucdat_1.am[5];
/* ..........SUM UP ABUNDANCES OF HEAVY NUCLIDES. */
/* Exception for helium. */
    xout_ref(flags_1.it, 10) = xout_ref(flags_1.it, 10) + xout_ref(flags_1.it,
	     11) + xout_ref(flags_1.it, 12) + xout_ref(flags_1.it, 13) + 
	    xout_ref(flags_1.it, 14) + xout_ref(flags_1.it, 15) + xout_ref(
	    flags_1.it, 16) + xout_ref(flags_1.it, 17) + xout_ref(flags_1.it, 
	    18) + xout_ref(flags_1.it, 19) + xout_ref(flags_1.it, 20) + 
	    xout_ref(flags_1.it, 21) + xout_ref(flags_1.it, 22) + xout_ref(
	    flags_1.it, 23) + xout_ref(flags_1.it, 24) + xout_ref(flags_1.it, 
	    25) + xout_ref(flags_1.it, 26);
/* ..........RELABEL TEMPERATURE, TIME, THERMODYNAMIC VARIABLES, ETC. */
/* Li8 to O16. */
    outdat_1.t9out[flags_1.it - 1] = evolp1_1.t9;
/* Temperature. */
    outdat_1.tout[flags_1.it - 1] = ttime_1.t;
/* Time. */
    thmout_ref(flags_1.it, 1) = thermcb_1.thm[0];
/* rho photon. */
    thmout_ref(flags_1.it, 2) = thermcb_1.thm[3];
/* rho electron. */
    thmout_ref(flags_1.it, 3) = thermcb_1.thm[7];
/* rho neutrino. */
    thmout_ref(flags_1.it, 4) = thermcb_1.thm[8];
/* rho baryon. */
    thmout_ref(flags_1.it, 5) = evolp1_1.phie;
/* Chemical potential. */
    thmout_ref(flags_1.it, 6) = thermcb_1.thm[9];
/* rho total. */
    outdat_1.dtout[flags_1.it - 1] = ttime_1.dt;
/* Time step. */
    outdat_1.etaout[flags_1.it - 1] = evolp1_1.hv / 33683.;
/* Baryon to photon ratio. */
    outdat_1.hubout[flags_1.it - 1] = thermcb_1.hubcst;
/* 20--------INDICATE TERMINATION OF ACCUMULATION IF APPROPRIATE------------------ */
/* Expansion rate. */
    if (flags_1.it == 40 || flags_1.ip < compr_1.inc) {
	flags_1.ltime = 1;
    }
    return 0;
} /* accum_ */

#undef thmout_ref
#undef xout_ref


/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int therm_(void)
{
    /* System generated locals */
    shortint s__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), cosh(doublereal), sinh(
	    doublereal);

    /* Local variables */
    static doublereal z__, z1, z2, z3, z4, z5;
    extern doublereal ex_(doublereal *);
    static doublereal cosh1, cosh2, cosh3, cosh4, cosh5, sinh1, sinh2, sinh3, 
	    sinh4, sinh5;
    extern /* Subroutine */ int bessel_(doublereal *), nudens_(void);

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
    z__ = 5.93 / evolp1_1.t9;
/* z = m(electron)c**2/k(t9). */
    nupar_1.tnu = pow_dd(&endens_1.rnb, &c_b12) * compr_1.t9i;
/* ..........FACTORS OF z. */
/* Neutrino temperature. */
    z1 = z__;
    z2 = z__ * z__;
    z3 = z__ * z__ * z__;
    z4 = z__ * z__ * z__ * z__;
    z5 = z__ * z__ * z__ * z__ * z__;
/* ..........TRIGNOMETRIC FUNCTION VALUES. */
    if (evolp1_1.phie <= 17.) {
/* No chance of overflow. */
	cosh1 = cosh(evolp1_1.phie);
	cosh2 = cosh(evolp1_1.phie * 2.);
	cosh3 = cosh(evolp1_1.phie * 3.);
	cosh4 = cosh(evolp1_1.phie * 4.);
	cosh5 = cosh(evolp1_1.phie * 5.);
	sinh1 = sinh(evolp1_1.phie);
	sinh2 = sinh(evolp1_1.phie * 2.);
	sinh3 = sinh(evolp1_1.phie * 3.);
	sinh4 = sinh(evolp1_1.phie * 4.);
	sinh5 = sinh(evolp1_1.phie * 5.);
    } else {
	cosh1 = 0.;
	cosh2 = 0.;
	cosh3 = 0.;
	cosh4 = 0.;
	cosh5 = 0.;
	sinh1 = 0.;
	sinh2 = 0.;
	sinh3 = 0.;
	sinh4 = 0.;
	sinh5 = 0.;
    }
    bessel_(&z__);
/* 20--------COMPUTE THERMODYNAMIC VARIABLES-------------------------------------- */
    thermcb_1.thm[0] = evolp1_1.t9 * 8.418 * evolp1_1.t9 * evolp1_1.t9 * 
	    evolp1_1.t9;
/* (Ref 1) */
    thermcb_1.thm[1] = thermcb_1.thm[0] * 4. / evolp1_1.t9;
/* (Ref 2) */
    thermcb_1.thm[2] = thermcb_1.thm[0] / 3.;
/* (Ref 3) */
    thermcb_1.thm[3] = (besselcb_1.bm1 * cosh1 - besselcb_1.bm2 * cosh2 + 
	    besselcb_1.bm3 * cosh3 - besselcb_1.bm4 * cosh4 + besselcb_1.bm5 *
	     cosh5) * 3206.;
/* (Ref 4) */
    thermcb_1.thm[4] = z__ / evolp1_1.t9 * 3206. * (besselcb_1.bn1 * cosh1 - 
	    besselcb_1.bn2 * 2. * cosh2 + besselcb_1.bn3 * 3. * cosh3 - 
	    besselcb_1.bn4 * 4. * cosh4 + besselcb_1.bn5 * 5. * cosh5);
/* (Ref 5) */
    thermcb_1.thm[5] = (besselcb_1.bm1 * sinh1 - besselcb_1.bm2 * 2. * sinh2 
	    + besselcb_1.bm3 * 3. * sinh3 - besselcb_1.bm4 * 4. * sinh4 + 
	    besselcb_1.bm5 * 5. * sinh5) * 3206.;
/* (Ref 6) */
    thermcb_1.thm[6] = (besselcb_1.bl1 * cosh1 / z__ - besselcb_1.bl2 * cosh2 
	    / (z__ * 2.) + besselcb_1.bl3 * cosh3 / (z__ * 3.) - 
	    besselcb_1.bl4 * cosh4 / (z__ * 4.) + besselcb_1.bl5 * cosh5 / (
	    z__ * 5.)) * 3206.;
/* (Ref 7) */
    if (modpr_1.xi[0] == 0. && modpr_1.xi[1] == 0. && modpr_1.xi[2] == 0.) {
/* Nonde */
	thermcb_1.thm[7] = modpr_1.xnu * endens_1.rhone0 * pow_dd(&
		endens_1.rnb, &c_b13);
/* (Ref 8) */
    } else {
/* Include effects of neutrino degenera */
	thermcb_1.thm[7] = 0.;
	s__1 = (shortint) modpr_1.xnu;
	for (nupar_1.nu = 1; nupar_1.nu <= s__1; ++nupar_1.nu) {
/* For every neutrino family. */
	    nudens_();
/* Compute neutrino energy density. */
	    thermcb_1.thm[7] += nupar_1.rhonu * 12.79264;
/* Have 12.79264 from units ch */
	}
    }
    thermcb_1.thm[8] = endens_1.rhob0 * endens_1.rnb;
/* (Ref 9) */
    thermcb_1.thm[9] = thermcb_1.thm[0] + thermcb_1.thm[3] + thermcb_1.thm[7] 
	    + thermcb_1.thm[8];
/* (Ref 10 */
/* Computing 3rd power */
    d__1 = z__;
    thermcb_1.thm[10] = -(d__1 * (d__1 * d__1) / evolp1_1.t9) * (sinh1 * (
	    besselcb_1.bl1 * 3. - z__ * besselcb_1.bm1) - sinh2 * (
	    besselcb_1.bl2 * 3. - z__ * 2. * besselcb_1.bm2) + sinh3 * (
	    besselcb_1.bl3 * 3. - z__ * 3. * besselcb_1.bm3) - sinh4 * (
	    besselcb_1.bl4 * 3. - z__ * 4. * besselcb_1.bm4) + sinh5 * (
	    besselcb_1.bl5 * 3. - z__ * 5. * besselcb_1.bm5));
/* (Ref 11 */
/* Computing 3rd power */
    d__1 = z__;
    thermcb_1.thm[11] = d__1 * (d__1 * d__1) * (cosh1 * besselcb_1.bl1 - 
	    cosh2 * 2. * besselcb_1.bl2 + cosh3 * 3. * besselcb_1.bl3 - cosh4 
	    * 4. * besselcb_1.bl4 + cosh5 * 5. * besselcb_1.bl5);
/* (Ref 12 */
    if (thermcb_1.thm[11] != 0.) {
	thermcb_1.thm[11] = 1. / thermcb_1.thm[11];
    }
    thermcb_1.thm[12] = .565 / z1 + 1. - 6.382 / z2 + 11.108 / z3 + 36.492 / 
	    z4 + 27.512 / z5;
/* (Ref 13 */
    d__1 = z__ * -2.531;
    thermcb_1.thm[13] = (5.252 / z1 - 16.229 / z2 + 18.059 / z3 + 34.181 / z4 
	    + 27.617 / z5) * ex_(&d__1);
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
/* Subroutine */ int bessel_(doublereal *z__)
{
    /* Local variables */
    static shortint i__;
    static doublereal r__;
#define blz ((doublereal *)&besselcb_1)
#define bmz ((doublereal *)&besselcb_1 + 5)
#define bnz ((doublereal *)&besselcb_1 + 10)
    extern /* Subroutine */ int knux_(doublereal *);

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
	bmz[i__ - 1] = .25 * (3. * kays_1.bk3 + kays_1.bk1) / r__;
/* Put value from function bm into arra */
	bnz[i__ - 1] = .5 * (kays_1.bk4 + kays_1.bk2) / r__;
/* Put value from function bn into arra */
    }
    return 0;
} /* bessel_ */

#undef bnz
#undef bmz
#undef blz


/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int knux_(doublereal *z__)
{
    /* Initialized data */

    static doublereal ci0[7] = { 1.,3.5156229,3.0899424,1.2067492,.2659732,
	    .0360768,.0045813 };
    static doublereal ci1[7] = { .5,.87890594,.51498869,.15084934,.02658733,
	    .00301532,3.2411e-4 };
    static doublereal ck0[7] = { -.57721566,.4227842,.23069756,.0348859,
	    .00262698,1.075e-4,7.4e-6 };
    static doublereal ck1[7] = { 1.,.15443144,-.67278579,-.18156897,
	    -.01919402,-.00110404,-4.686e-5 };
    static doublereal c0[7] = { 1.25331414,-.07832358,.02189568,-.01062446,
	    .00587872,-.0025154,5.3208e-4 };
    static doublereal c1[7] = { 1.25331414,.23498619,-.0365562,.01504268,
	    -.00780353,.00325614,-6.8245e-4 };

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), sqrt(doublereal);

    /* Local variables */
    static shortint i__;
    static doublereal t, y;
    extern doublereal ex_(doublereal *);
    static doublereal bi0, bi1, coeff;

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
    if (*z__ <= 2.) {
/* ..........COMPUTE FACTORS. */
/* (Ref. 1). */
	t = *z__ / 3.75;
	y = *z__ / 2;
	coeff = log(y);
/* ..........VALUES FOR i0(z) and i1(z). */
	bi0 = ci0[0];
	bi1 = ci1[0];
	kays_1.bk0 = ck0[0];
	kays_1.bk1 = ck1[0];
	for (i__ = 2; i__ <= 7; ++i__) {
	    i__1 = i__ - 1 << 1;
	    bi0 += ci0[i__ - 1] * pow_di(&t, &i__1);
	    i__1 = i__ - 1 << 1;
	    bi1 += ci1[i__ - 1] * pow_di(&t, &i__1);
	    i__1 = i__ - 1 << 1;
	    kays_1.bk0 += ck0[i__ - 1] * pow_di(&y, &i__1);
	    i__1 = i__ - 1 << 1;
	    kays_1.bk1 += ck1[i__ - 1] * pow_di(&y, &i__1);
	}
/* ..........VALUES FOR k0(z) and k1(z). */
	kays_1.bk0 = -coeff * bi0 + kays_1.bk0;
	kays_1.bk1 = coeff * bi1 * *z__ + kays_1.bk1 / *z__;
    } else {
/* ..........COMPUTE FACTORS. */
/* (z.le.2.)               !(Ref. 2). */
	y = 2. / *z__;
	d__1 = -(*z__);
	coeff = ex_(&d__1) / sqrt(*z__);
/* ..........VALUES FOR k0(z) and k1(z). */
	kays_1.bk0 = c0[0];
	kays_1.bk1 = c1[0];
	for (i__ = 2; i__ <= 7; ++i__) {
	    i__1 = i__ - 1;
	    kays_1.bk0 += c0[i__ - 1] * pow_di(&y, &i__1);
	    i__1 = i__ - 1;
	    kays_1.bk1 += c1[i__ - 1] * pow_di(&y, &i__1);
	}
	kays_1.bk0 = coeff * kays_1.bk0;
	kays_1.bk1 = coeff * kays_1.bk1;
    }
/* 20--------FIND K2, K3, AND K4 BY ITERATION (Ref. 3)---------------------------- */
/* (z.le.2.) */
    kays_1.bk2 = kays_1.bk1 / *z__ * 2. + kays_1.bk0;
/* k2(z). */
    kays_1.bk3 = kays_1.bk2 / *z__ * 4. + kays_1.bk1;
/* k3(z). */
    kays_1.bk4 = kays_1.bk3 / *z__ * 6. + kays_1.bk2;
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
    doublereal d__1, d__2, d__3;

    /* Local variables */
    extern /* Subroutine */ int func5_(), func6_();
    extern doublereal xintd_(doublereal *, doublereal *, U_fp, shortint *);
    static doublereal uplim1, uplim2;

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
    if ((d__1 = modpr_1.xi[nupar_1.nu - 1], abs(d__1)) <= .03) {
/* ..........SMALL xi APPROXIMATION. */
/* Computing 4th power */
	d__1 = nupar_1.tnu, d__1 *= d__1;
/* Computing 2nd power */
	d__2 = modpr_1.xi[nupar_1.nu - 1];
/* Computing 4th power */
	d__3 = modpr_1.xi[nupar_1.nu - 1], d__3 *= d__3;
	nupar_1.rhonu = d__1 * d__1 * .65797251520666655 * (d__2 * d__2 * 
		.37995508052714932 + .875 + d__3 * d__3 * .01924878176245234);
    } else {
	if ((d__1 = modpr_1.xi[nupar_1.nu - 1], abs(d__1)) >= 30.) {
/* ..........LARGE xi APPROXIMATION. */
/* Computing 4th power */
	    d__1 = nupar_1.tnu, d__1 *= d__1;
/* Computing 4th power */
	    d__2 = modpr_1.xi[nupar_1.nu - 1], d__2 *= d__2;
/* Computing 2nd power */
	    d__3 = modpr_1.xi[nupar_1.nu - 1];
	    nupar_1.rhonu = d__1 * d__1 / 78.956701824799993 * (d__2 * d__2) *
		     (19.740000000000002 / (d__3 * d__3) + 1);
	} else {
/* ..........DO INTEGRATION */
	    uplim1 = (modpr_1.xi[nupar_1.nu - 1] + 88.029) * nupar_1.tnu;
	    uplim2 = (88.029 - modpr_1.xi[nupar_1.nu - 1]) * nupar_1.tnu;
	    if (uplim2 <= 0.) {
		nupar_1.rhonu = xintd_(&c_b17, &uplim1, (U_fp)func5_, &cs__50)
			;
	    } else {
		nupar_1.rhonu = xintd_(&c_b17, &uplim1, (U_fp)func5_, &cs__50)
			 + xintd_(&c_b17, &uplim2, (U_fp)func6_, &cs__50);
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
doublereal func1_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal ex_(doublereal *);
    static doublereal part1, part2;

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
    if (*x <= 0.) {
	ret_val = 0.;
    } else {
	d__1 = *x * -.511 / nupar_1.t9mev;
	part1 = 1. / (ex_(&d__1) + 1.);
	d__1 = (*x - 2.531) * (.511 / nupar_1.tnmev) - modpr_1.xi[0];
	part2 = 1. / (ex_(&d__1) + 1.);
/* Computing 2nd power */
	d__1 = *x - 2.531;
/* Computing 2nd power */
	d__3 = *x;
	d__2 = d__3 * d__3 - 1;
	ret_val = nupar_1.cnorm * *x * (d__1 * d__1) * pow_dd(&d__2, &c_b24) *
		 part1 * part2;
    }
    return ret_val;
} /* func1_ */

/* 20--------2ND PART OF INTEGRAL FOR n->p RATE----------------------------------- */
doublereal func2_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal ex_(doublereal *);
    static doublereal part1, part2;

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
    if (*x <= 1.) {
	ret_val = 0.;
    } else {
	d__1 = *x * .511 / nupar_1.t9mev;
	part1 = 1. / (ex_(&d__1) + 1.);
	d__1 = -(*x + 2.531) * (.511 / nupar_1.tnmev) - modpr_1.xi[0];
	part2 = 1. / (ex_(&d__1) + 1.);
/* Computing 2nd power */
	d__1 = *x + 2.531;
/* Computing 2nd power */
	d__3 = *x;
	d__2 = d__3 * d__3 - 1;
	ret_val = nupar_1.cnorm * *x * (d__1 * d__1) * pow_dd(&d__2, &c_b24) *
		 part1 * part2;
    }
    return ret_val;
} /* func2_ */

/* 30--------1ST PART OF INTEGRAL FOR p->n RATE----------------------------------- */
doublereal func3_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal ex_(doublereal *);
    static doublereal part1, part2;

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
    if (*x <= 1.) {
	ret_val = 0.;
    } else {
	d__1 = *x * -.511 / nupar_1.t9mev;
	part1 = 1. / (ex_(&d__1) + 1.);
	d__1 = (*x + 2.531) * (.511 / nupar_1.tnmev) + modpr_1.xi[0];
	part2 = 1. / (ex_(&d__1) + 1.);
/* Computing 2nd power */
	d__1 = *x + 2.531;
/* Computing 2nd power */
	d__3 = *x;
	d__2 = d__3 * d__3 - 1;
	ret_val = nupar_1.cnorm * *x * (d__1 * d__1) * pow_dd(&d__2, &c_b24) *
		 part1 * part2;
    }
    return ret_val;
} /* func3_ */

/* 40--------2ND PART OF INTEGRAL FOR p->n RATE----------------------------------- */
doublereal func4_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal ex_(doublereal *);
    static doublereal part1, part2;

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
    if (*x <= 1.) {
	ret_val = 0.;
    } else {
	d__1 = *x * .511 / nupar_1.t9mev;
	part1 = 1. / (ex_(&d__1) + 1.);
	d__1 = -(*x - 2.531) * (.511 / nupar_1.tnmev) + modpr_1.xi[0];
	part2 = 1. / (ex_(&d__1) + 1.);
/* Computing 2nd power */
	d__1 = *x - 2.531;
/* Computing 2nd power */
	d__3 = *x;
	d__2 = d__3 * d__3 - 1;
	ret_val = nupar_1.cnorm * *x * (d__1 * d__1) * pow_dd(&d__2, &c_b24) *
		 part1 * part2;
    }
    return ret_val;
} /* func4_ */

/* 50--------INTEGRAL FOR ENERGY DENSITY OF NEUTRINO------------------------------ */
doublereal func5_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

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
    d__1 = *x;
    ret_val = d__1 * (d__1 * d__1) * .05066067740361991 / (exp(*x / 
	    nupar_1.tnu - modpr_1.xi[nupar_1.nu - 1]) + 1.);
    return ret_val;
} /* func5_ */

/* 60--------INTEGRAL FOR ENERGY DENSITY OF ANTINEUTRINO-------------------------- */
doublereal func6_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

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
    d__1 = *x;
    ret_val = d__1 * (d__1 * d__1) * .05066067740361991 / (exp(*x / 
	    nupar_1.tnu + modpr_1.xi[nupar_1.nu - 1]) + 1.);
    return ret_val;
} /* func6_ */

/* ----------REFERENCES----------------------------------------------------------- */
/*     Forms of the integrals involved can be found in */
/*       Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683. */
/*       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415. */
/* ========================IDENTIFICATION DIVISION================================ */
doublereal xintd_(doublereal *xlow, doublereal *xhi, D_fp func, shortint *nq)
{
    /* Initialized data */

    static doublereal u[6] = { -.93246951420315,-.66120938646627,
	    -.2386191860832,.2386191860832,.66120938646627,.93246951420315 };
    static doublereal w[6] = { .17132449237917,.36076157304814,
	    .46791393457269,.46791393457269,.36076157304814,.17132449237917 };
    static shortint np = 6;

    /* System generated locals */
    shortint s__1, s__2;
    doublereal ret_val;

    /* Local variables */
    static doublereal f, x, sum, cent, dist;
    static shortint nint, npnt;

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
    sum = 0.;
    dist = (*xhi - *xlow) / (doublereal) (*nq);
/* Size of quad interval. */
    s__1 = *nq;
    for (nint = 1; nint <= s__1; ++nint) {
	cent = *xlow + ((doublereal) nint - .5) * dist;
/* Center of interval. */
	s__2 = np;
	for (npnt = 1; npnt <= s__2; ++npnt) {
	    x = cent + dist * .5 * u[npnt - 1];
/* Integration point. */
	    f = (*func)(&x);
/* Evaluate function x(1). */
	    sum += f * w[npnt - 1];
/* Add up sum. */
	}
    }
/* 20--------GET INTEGRAL VALUE--------------------------------------------------- */
    ret_val = sum * dist * .5;
/* Do integral. */
    return ret_val;
} /* xintd_ */

/* ========================IDENTIFICATION DIVISION================================ */
doublereal ex_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double exp(doublereal);

/* ----------LINKAGES. */
/*     CALLED BY - [subroutine] start, rate2, rate3, rate4, sol */
/*               - [function] eval */
/*     CALLS     - none */
/* ----------REMARKS. */
/*     Exponential function with underflow precaution. */
/* ===========================PROCEDURE DIVISION================================== */
    if (*x > 88.029) {
/* In danger of overflow. */
	ret_val = exp(88.029);
    } else {
	if (*x < -88.722) {
/* In danger of underflow. */
	    ret_val = 0.;
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
/* Subroutine */ int sol_(shortint *loop)
{
    /* Initialized data */

    static doublereal si[11] = { 1.,1.,1.,1.,1.,2.,3.,2.,1.,1.,2. };
    static doublereal sj[11] = { 0.,1.,1.,0.,1.,0.,0.,1.,1.,1.,0. };
    static doublereal sk[11] = { 0.,0.,1.,0.,0.,1.,0.,0.,1.,0.,2. };
    static doublereal sl[11] = { 1.,1.,1.,2.,2.,1.,1.,1.,2.,3.,1. };

    /* Format strings */
    static char fmt_6000[] = "(\002 \002,\002** y(\002,i2,\002) fails to con"
	    "verge **\002)";
    static char fmt_6002[] = "(\002 \002,\002** \002,i2,\002 th diagonal ter"
	    "m equals zero **\002)";

    /* System generated locals */
    shortint s__1, s__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static shortint i__, j, k, l, n, i1, j1;
    static doublereal ci, cj, ck, cl, ri, rj, rk, rl;
    extern doublereal ex_(doublereal *);
    static doublereal yy[26], t932;
    static shortint ind;
    static doublereal t9m32, bdln;
    static shortint isize1;
    extern /* Subroutine */ int eqslin_(shortint *, shortint *);
    static shortint ierror;

    /* Fortran I/O blocks */
    static cilist io___103 = { 0, 6, 0, fmt_6000, 0 };
    static cilist io___104 = { 0, 6, 0, fmt_6002, 0 };



#define a_ref(a_1,a_2) lncoef_1.a[(a_2)*26 + a_1 - 27]

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
/* Dynamic variable */
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
    t932 = pow_dd(&evolp1_1.t9, &c_b5);
/* t9**(3/2). */
    t9m32 = 1. / t932;
/* ..........MATRIX SIZE. */
/* t9**(-3/2). */
    isize1 = runopt_1.isize + 1;
/* ..........INITIALIZE A-MATRIX. */
    s__1 = runopt_1.isize;
    for (i__ = 1; i__ <= s__1; ++i__) {
	s__2 = runopt_1.isize;
	for (j = 1; j <= s__2; ++j) {
	    a_ref(j, i__) = 0.;
/* Set a-matrix to zero. */
	}
    }
/* 20--------COMPUTE FACTORS FOR THE A-MATRIX------------------------------------- */
    s__1 = runopt_1.jsize;
    for (n = 1; n <= s__1; ++n) {
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
	    ci = rates_2.f[n - 1];
/* (Ref 1). */
	    cj = 0.;
	    ck = 0.;
	    cl = rates_2.r__[n - 1];
	    goto L212;
L202:
/* 1-1-0-1 configuration. */
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * 1e10 * t932 * ex_(&d__1)
		     * rates_2.f[n - 1];
/* (Ref 2). */
	    rates_2.f[n - 1] = endens_1.rhob * rates_2.f[n - 1];
	    ci = evolp1_1.y[j - 1] * rates_2.f[n - 1] / 2.;
	    cj = evolp1_1.y[i__ - 1] * rates_2.f[n - 1] / 2.;
	    ck = 0.;
	    cl = rates_2.r__[n - 1];
	    goto L212;
L203:
/* 1-1-1-1 configuration. */
	    rates_2.f[n - 1] = endens_1.rhob * rates_2.f[n - 1];
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * ex_(&d__1) * rates_2.f[
		    n - 1];
/* (Ref 3). */
	    ci = evolp1_1.y[j - 1] * rates_2.f[n - 1] / 2.;
	    cj = evolp1_1.y[i__ - 1] * rates_2.f[n - 1] / 2.;
	    ck = evolp1_1.y[l - 1] * rates_2.r__[n - 1] / 2.;
	    cl = evolp1_1.y[k - 1] * rates_2.r__[n - 1] / 2.;
	    goto L212;
L204:
/* 1-0-0-2 configuration. */
	    ci = rates_2.f[n - 1];
	    cj = 0.;
	    ck = 0.;
	    cl = evolp1_1.y[l - 1] * rates_2.r__[n - 1] / 2.;
	    goto L212;
L205:
/* 1-1-0-2 configuration. */
	    rates_2.f[n - 1] = endens_1.rhob * rates_2.f[n - 1];
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * ex_(&d__1) * rates_2.f[
		    n - 1];
/* (Ref 3). */
	    ci = evolp1_1.y[j - 1] * rates_2.f[n - 1] / 2.;
	    cj = evolp1_1.y[i__ - 1] * rates_2.f[n - 1] / 2.;
	    ck = 0.;
	    cl = evolp1_1.y[l - 1] * rates_2.r__[n - 1] / 2.;
	    goto L212;
L206:
/* 2-0-1-1 configuration. */
	    rates_2.f[n - 1] = endens_1.rhob * rates_2.f[n - 1];
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * ex_(&d__1) * rates_2.f[
		    n - 1];
/* (Ref 3). */
	    ci = evolp1_1.y[i__ - 1] * rates_2.f[n - 1] / 2.;
	    cj = 0.;
	    ck = evolp1_1.y[l - 1] * rates_2.r__[n - 1] / 2.;
	    cl = evolp1_1.y[k - 1] * rates_2.r__[n - 1] / 2.;
	    goto L212;
L207:
/* 3-0-0-1 configuration. */
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * 1e20 * t932 * t932 * 
		    ex_(&d__1) * rates_2.f[n - 1];
/* (Ref 4). */
	    rates_2.f[n - 1] = endens_1.rhob * endens_1.rhob * rates_2.f[n - 
		    1];
	    ci = evolp1_1.y[i__ - 1] * evolp1_1.y[i__ - 1] * rates_2.f[n - 1] 
		    / 6.;
	    cj = 0.;
	    ck = 0.;
	    cl = rates_2.r__[n - 1];
	    goto L212;
L208:
/* 2-1-0-1 configuration. */
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * 1e20 * t932 * t932 * 
		    ex_(&d__1) * rates_2.f[n - 1];
/* (Ref 4). */
	    rates_2.f[n - 1] = endens_1.rhob * endens_1.rhob * rates_2.f[n - 
		    1];
	    ci = evolp1_1.y[j - 1] * evolp1_1.y[i__ - 1] * rates_2.f[n - 1] / 
		    3.;
	    cj = evolp1_1.y[i__ - 1] * evolp1_1.y[i__ - 1] * rates_2.f[n - 1] 
		    / 6.;
	    ck = 0.;
	    cl = rates_2.r__[n - 1];
	    goto L212;
L209:
/* 1-1-1-2 configuration. */
	    rates_2.f[n - 1] = endens_1.rhob * rates_2.f[n - 1];
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * 1e-10 * t9m32 * 
		    endens_1.rhob * ex_(&d__1) * rates_2.f[n - 1];
/* (Ref 5) */
	    ci = evolp1_1.y[j - 1] * rates_2.f[n - 1] / 2.;
	    cj = evolp1_1.y[i__ - 1] * rates_2.f[n - 1] / 2.;
	    ck = evolp1_1.y[l - 1] * evolp1_1.y[l - 1] * rates_2.r__[n - 1] / 
		    6.;
	    cl = evolp1_1.y[k - 1] * evolp1_1.y[l - 1] * rates_2.r__[n - 1] / 
		    3.;
	    goto L212;
L210:
/* 1-1-0-3 configuration. */
	    rates_2.f[n - 1] = endens_1.rhob * rates_2.f[n - 1];
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * 1e-10 * t9m32 * 
		    endens_1.rhob * ex_(&d__1) * rates_2.f[n - 1];
/* (Ref 5) */
	    ci = evolp1_1.y[j - 1] * rates_2.f[n - 1] / 2.;
	    cj = evolp1_1.y[i__ - 1] * rates_2.f[n - 1] / 2.;
	    ck = 0.;
	    cl = evolp1_1.y[l - 1] * evolp1_1.y[l - 1] * rates_2.r__[n - 1] / 
		    6.;
	    goto L212;
L211:
/* 2-0-2-1 configuration. */
	    rates_2.f[n - 1] = endens_1.rhob * rates_2.f[n - 1];
	    d__1 = -recpr_1.q9[n - 1] / evolp1_1.t9;
	    rates_2.r__[n - 1] = recpr_1.rev[n - 1] * 1e-10 * t9m32 * 
		    endens_1.rhob * ex_(&d__1) * rates_2.f[n - 1];
/* (Ref 5) */
	    ci = evolp1_1.y[i__ - 1] * rates_2.f[n - 1] / 2.;
	    cj = 0.;
	    ck = evolp1_1.y[l - 1] * evolp1_1.y[k - 1] * rates_2.r__[n - 1] / 
		    3.;
	    cl = evolp1_1.y[k - 1] * evolp1_1.y[k - 1] * rates_2.r__[n - 1] / 
		    6.;
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
		a_ref(j, i__) = a_ref(j, i__) + rj * ci;
	    }
	    if (k <= runopt_1.isize) {
		a_ref(k, i__) = a_ref(k, i__) - rk * ci;
	    }
	    a_ref(i__, i__) = a_ref(i__, i__) + ri * ci;
	    a_ref(l, i__) = a_ref(l, i__) - rl * ci;
/* ..........FILL J NUCLIDE COLUMN. */
	    if (j <= runopt_1.isize) {
		a_ref(j, j) = a_ref(j, j) + rj * cj;
		if (k <= runopt_1.isize) {
		    a_ref(k, j) = a_ref(k, j) - rk * cj;
		}
		a_ref(i__, j) = a_ref(i__, j) + ri * cj;
		a_ref(l, j) = a_ref(l, j) - rl * cj;
	    }
/* ..........FILL K NUCLIDE COLUMN. */
	    if (k <= runopt_1.isize) {
		if (j <= runopt_1.isize) {
		    a_ref(j, k) = a_ref(j, k) - rj * ck;
		}
		a_ref(k, k) = a_ref(k, k) + rk * ck;
		a_ref(i__, k) = a_ref(i__, k) - ri * ck;
		a_ref(l, k) = a_ref(l, k) + rl * ck;
	    }
/* ..........FILL L NUCLIDE COLUMN. */
	    if (j <= runopt_1.isize) {
		a_ref(j, l) = a_ref(j, l) - rj * cl;
	    }
	    if (k <= runopt_1.isize) {
		a_ref(k, l) = a_ref(k, l) + rk * cl;
	    }
	    a_ref(i__, l) = a_ref(i__, l) - ri * cl;
	    a_ref(l, l) = a_ref(l, l) + rl * cl;
	}
/* ((ind.ne.0).and.(i.le.isize).and.(l.le.isize)) */
    }
/* 40--------PUT A-MATRIX AND B-VECTOR IN FINAL FORM OF MATRIX EQUATION----------- */
/* n = 1,jsize */
    bdln = thermcb_1.hubcst * 3. * 1e-5;
/* (10**(-5))*(Expansion rate). */
    s__1 = runopt_1.isize;
    for (i__ = 1; i__ <= s__1; ++i__) {
	i1 = isize1 - i__;
/* Invert the rows. */
	s__2 = runopt_1.isize;
	for (j = 1; j <= s__2; ++j) {
	    j1 = isize1 - j;
/* Invert the columns. */
	    if ((d__1 = a_ref(j, i__), abs(d__1)) < bdln * evolp3_1.y0[j1 - 1]
		     / evolp3_1.y0[i1 - 1]) {
		a_ref(j, i__) = 0.;
/* Set 0 if tiny. */
	    } else {
		a_ref(j, i__) = a_ref(j, i__) * ttime_1.dt;
/* Bring dt over to other side. */
	    }
	}
	a_ref(i__, i__) = a_ref(i__, i__) + 1.;
/* Add identity matrix to a-matrix. */
	lncoef_1.b[i1 - 1] = evolp3_1.y0[i__ - 1];
/* Initial abundances. */
    }
/* 50--------SOLVE EQUATIONS TO GET DERIVATIVE------------------------------------ */
/* ..........SET MONITOR FLAG AND SOLVE BY GAUSSIAN ELIMINATION. */
    if (*loop == 1) {
	eqslin_(&flags_1.ip, &ierror);
    } else {
	eqslin_(&cs__0, &ierror);
    }
/* ..........OBTAIN DERIVATIVE. */
    s__1 = runopt_1.isize;
    for (i__ = 1; i__ <= s__1; ++i__) {
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
	    do_fio(&c__1, (char *)&ierror, (ftnlen)sizeof(shortint));
	    e_wsfe();
	}
/* Error message. */
	if (flags_1.mbad >= 1) {
	    s_wsfe(&io___104);
	    do_fio(&c__1, (char *)&flags_1.mbad, (ftnlen)sizeof(shortint));
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

#undef a_ref


/* ========================IDENTIFICATION DIVISION================================ */
/* Subroutine */ int eqslin_(shortint *icnvm, shortint *ierror)
{
    /* System generated locals */
    shortint s__1, s__2, s__3;
    doublereal d__1;

    /* Local variables */
    static shortint i__, j, k;
    static doublereal r__, x[26], a0[676]	/* was [26][26] */, cx, sum, 
	    xdy;
    static shortint nord;


#define a_ref(a_1,a_2) lncoef_2.a[(a_2)*26 + a_1 - 27]
#define a0_ref(a_1,a_2) a0[(a_2)*26 + a_1 - 27]

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
    s__1 = runopt_1.isize;
    for (i__ = 1; i__ <= s__1; ++i__) {
	x[i__ - 1] = lncoef_2.b[i__ - 1];
/* Right-hand vector. */
	lncoef_2.y[i__ - 1] = 0.;
/* Solution vector. */
    }
/* ..........SAVE MATRIX. */
    if (*icnvm == compr_1.inc) {
/* Monitor convergence. */
	s__1 = runopt_1.isize;
	for (i__ = 1; i__ <= s__1; ++i__) {
	    s__2 = runopt_1.isize;
	    for (j = 1; j <= s__2; ++j) {
		a0_ref(j, i__) = a_ref(j, i__);
/* Initial value of coefficient array. */
	    }
	}
    }
/* 20--------TRIANGULARIZE MATRIX AND SAVE OPERATOR------------------------------- */
/* ..........CHECK TO SEE THAT THERE ARE NO ZEROES AT PIVOT POINTS. */
    s__1 = runopt_1.isize - 1;
    for (i__ = 1; i__ <= s__1; ++i__) {
	if (a_ref(i__, i__) == 0.) {
/* Don't want to divide by zero. */
	    flags_1.mbad = i__;
/* Position of zero coefficient. */
	    return 0;
/* Terminate matrix evaluation. */
	}
/* ..........TRIANGULARIZE MATRIX. */
	s__2 = runopt_1.isize;
	for (j = i__ + 1; j <= s__2; ++j) {
	    if (a_ref(j, i__) != 0.) {
/* Progress diagonally down the column. */
		cx = a_ref(j, i__) / a_ref(i__, i__);
/* Scaling factor down the column. */
		s__3 = runopt_1.isize;
		for (k = i__ + 1; k <= s__3; ++k) {
/* Progress diagonally along row. */
		    a_ref(j, k) = a_ref(j, k) - cx * a_ref(i__, k);
/* Subtract scaled coeff along */
		}
		a_ref(j, i__) = cx;
/* ..........OPERATE ON RIGHT-HAND VECTOR. */
/* Scaled coefficient. */
		x[j - 1] -= cx * x[i__ - 1];
/* Subtract off scaled coefficient. */
	    }
	}
    }
/* 30--------DO BACK SUBSTITUTION------------------------------------------------- */
L300:
    x[runopt_1.isize - 1] /= a_ref(runopt_1.isize, runopt_1.isize);
/* Solution for ultimate positi */
    lncoef_2.y[runopt_1.isize - 1] += x[runopt_1.isize - 1];
    for (i__ = runopt_1.isize - 1; i__ >= 1; --i__) {
/* From i = penultimate to i = 1. */
	sum = 0.;
	s__1 = runopt_1.isize;
	for (j = i__ + 1; j <= s__1; ++j) {
	    sum += a_ref(i__, j) * x[j - 1];
/* Sum up all previous terms. */
	}
	x[i__ - 1] = (x[i__ - 1] - sum) / a_ref(i__, i__);
	lncoef_2.y[i__ - 1] += x[i__ - 1];
/* Add difference to initial value. */
    }
/* 40--------TESTS AND EXITS------------------------------------------------------ */
    if (*icnvm == compr_1.inc) {
	s__1 = runopt_1.isize;
	for (i__ = 1; i__ <= s__1; ++i__) {
	    if (lncoef_2.y[i__ - 1] != 0.) {
		xdy = (d__1 = x[i__ - 1] / lncoef_2.y[i__ - 1], abs(d__1));
/* Relative error. */
		if (xdy > 2e-4) {
		    if (nord < 1) {
/* Continue to higher orders. */
			++nord;
/* ..........FIND ERROR IN RIGHT-HAND VECTOR. */
			s__2 = runopt_1.isize;
			for (j = 1; j <= s__2; ++j) {
			    r__ = 0.;
/* Initialize r. */
			    s__3 = runopt_1.isize;
			    for (k = 1; k <= s__3; ++k) {
				r__ += a0_ref(j, k) * lncoef_2.y[k - 1];
/* Left side with approximate sol */
			    }
			    x[j - 1] = lncoef_2.b[j - 1] - r__;
/* Subtract difference from right side. */
			}
/* ..........OPERATE ON RIGHT-HAND VECTOR. */
			s__2 = runopt_1.isize - 1;
			for (j = 1; j <= s__2; ++j) {
			    s__3 = runopt_1.isize;
			    for (k = j + 1; k <= s__3; ++k) {
				x[k - 1] -= a_ref(k, j) * x[j - 1];
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

#undef a0_ref
#undef a_ref


