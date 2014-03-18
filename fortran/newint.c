/* newint.f -- translated by f2c (version 20090411).
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
    real reacpr[704]	/* was [88][8] */;
} recpr0_;

#define recpr0_1 recpr0_

struct {
    integer iform[88], ii[88], jj[88], kk[88], ll[88];
    real rev[88], q9[88];
} recpr_;

#define recpr_1 recpr_

struct {
    real f[88], r__[88];
} rates_;

#define rates_1 rates_

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
    real cy0, ct0, t9i0, t9f0, ytmin0;
    integer inc0;
} compr0_;

#define compr0_1 compr0_

struct {
    real cy, ct, t9i, t9f, ytmin;
    integer inc;
} compr_;

#define compr_1 compr_

struct {
    real c0[3], cosmo0, xi0[3];
} modpr0_;

#define modpr0_1 modpr0_

struct {
    real g, tau, xnu, c__[3], cosmo, xi[3], b;
} modpr_;

#define modpr_1 modpr_

struct {
    real dt0, eta0;
} varpr0_;

#define varpr0_1 varpr0_

struct {
    real dt1, eta1;
} varpr_;

#define varpr_1 varpr_

struct {
    real t, dt, dlt9dt;
} ttime_;

#define ttime_1 ttime_

struct {
    real thm[14], hubcst;
} thermcb_;

#define thermcb_1 thermcb_

struct {
    real rhone0, rhob0, rhob, rnb;
} endens_;

#define endens_1 endens_

struct {
    doublereal a[676]	/* was [26][26] */;
    real b___[26], yx[26];
} lncoef_;

#define lncoef_1 lncoef_

struct {
    real am[26], zm[26], dm[26];
} nucdat_;

#define nucdat_1 nucdat_

struct {
    real bl1, bl2, bl3, bl4, bl5, bm1, bm2, bm3, bm4, bm5, bn1, bn2, bn3, bn4,
	     bn5;
} besselcb_;

#define besselcb_1 besselcb_

struct {
    real bk0, bk1, bk2, bk3, bk4;
} kays_;

#define kays_1 kays_

struct {
    integer ltime, is, ip, it, mbad;
} flags_;

#define flags_1 flags_

struct {
    integer itime;
} checkcb_;

#define checkcb_1 checkcb_

struct {
    real xout[1040]	/* was [40][26] */, thmout[240]	/* was [40][6] */, 
	    t9out[40], tout[40], dtout[40], etaout[40], hubout[40];
} outdat_;

#define outdat_1 outdat_

struct {
    real t9mev, tnmev, tnu, cnorm;
    integer nu;
    real rhonu;
} nupar_;

#define nupar_1 nupar_

struct {
    integer irun, isize, jsize;
} runopt_;

#define runopt_1 runopt_

struct {
    integer nout;
    logical outfile;
} outopt_;

#define outopt_1 outopt_

/* Table of constant values */

static integer c__1 = 1;

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
    static cilist io___1 = { 0, 3, 0, fmt_200, 0 };


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
	s_wsfe(&io___1);
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

