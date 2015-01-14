/* new123.f -- translated by f2c (version 20100827).
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
    integer itime;
} checkcb_;

#define checkcb_1 checkcb_

struct {
    integer irun, isize, jsize;
} runopt_;

#define runopt_1 runopt_

struct {
    integer nout;
    logical outfile;
} outopt_;

#define outopt_1 outopt_

struct {
    integer ltime, is, ip, it, mbad;
} flags_;

#define flags_1 flags_

struct {
    real xout[1040]	/* was [40][26] */, thmout[240]	/* was [40][6] */, 
	    t9out[40], tout[40], dtout[40], etaout[40], hubout[40];
} outdat_;

#define outdat_1 outdat_

/* Table of constant values */

static integer c__1 = 1;
static integer c__4 = 4;
static integer c__3 = 3;
static doublereal c_b255 = 10.;

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
    static cilist io___1 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___3 = { 0, 5, 0, 0, 0 };
    static cilist io___4 = { 0, 6, 0, fmt_3000, 0 };
    static cilist io___5 = { 0, 5, 0, fmt_3001, 0 };



#define reacpr_ref(a_1,a_2) recpr0_1.reacpr[(a_2)*88 + a_1 - 89]

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
    s_wsfe(&io___1);
    e_wsfe();
/* 20--------INPUT INITIALIZATION INFORMATION AND PAUSE--------------------------- */
    for (i__ = 1; i__ <= 88; ++i__) {
/* ..........READ IN REACTION PARAMETERS. */
	recpr_1.iform[i__ - 1] = (integer) reacpr_ref(i__, 2);
/* Reaction type. */
	recpr_1.ii[i__ - 1] = (integer) reacpr_ref(i__, 3);
/* Incoming nuclide type. */
	recpr_1.jj[i__ - 1] = (integer) reacpr_ref(i__, 4);
/* Incoming nuclide type. */
	recpr_1.kk[i__ - 1] = (integer) reacpr_ref(i__, 5);
/* Outgoing nuclide type. */
	recpr_1.ll[i__ - 1] = (integer) reacpr_ref(i__, 6);
/* Outgoing nuclide type. */
	recpr_1.rev[i__ - 1] = reacpr_ref(i__, 7);
/* Reverse reaction coefficient. */
	recpr_1.q9[i__ - 1] = reacpr_ref(i__, 8);
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
    s_rsle(&io___3);
    e_rsle();
/* 30--------PRINT MENU AND AWAIT RESPONSE---------------------------------------- */
/* ..........RETURN FROM LOOPING. */
/* Pause. */
L300:
/* ..........DISPLAY MENU. */
    s_wsfe(&io___4);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___5);
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

#undef reacpr_ref


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
    static cilist io___7 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___8 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___10 = { 0, 6, 0, fmt_2100, 0 };
    static cilist io___11 = { 0, 5, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___13 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_2202, 0 };
    static cilist io___15 = { 0, 5, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, fmt_2300, 0 };
    static cilist io___17 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_2302, 0 };
    static cilist io___19 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_2304, 0 };
    static cilist io___21 = { 0, 5, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_2400, 0 };
    static cilist io___23 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___24 = { 0, 6, 0, fmt_2402, 0 };
    static cilist io___25 = { 0, 5, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_2500, 0 };
    static cilist io___27 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_2502, 0 };
    static cilist io___29 = { 0, 5, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_2600, 0 };
    static cilist io___31 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_2602, 0 };
    static cilist io___33 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_2604, 0 };
    static cilist io___35 = { 0, 5, 0, 0, 0 };


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
    s_wsfe(&io___7);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___8);
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
    s_wsfe(&io___10);
    e_wsfe();
    s_rsle(&io___11);
    e_rsle();
    goto L300;
/* 22--------SET UP RUN SECTION--------------------------------------------------- */
L220:
/* Setting up a run section. */
    s_wsfe(&io___12);
    e_wsfe();
    s_rsfe(&io___13);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___14);
	e_wsfe();
	s_rsle(&io___15);
	e_rsle();
	goto L300;
    } else {
	goto L300;
    }
/* 23--------RUN PROGRAM SECTION-------------------------------------------------- */
/* (inum.eq.1) */
L230:
/* Running the program section. */
    s_wsfe(&io___16);
    e_wsfe();
    s_rsfe(&io___17);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___18);
	e_wsfe();
	s_rsfe(&io___19);
	do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
	e_rsfe();
	if (inum == 1) {
	    s_wsfe(&io___20);
	    e_wsfe();
	    s_rsle(&io___21);
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
    s_wsfe(&io___22);
    e_wsfe();
    s_rsfe(&io___23);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___24);
	e_wsfe();
	s_rsle(&io___25);
	e_rsle();
	goto L300;
    } else {
	goto L300;
    }
/* 25--------METHOD OF COMPUTATION SECTION---------------------------------------- */
/* (inum.eq.1) */
L250:
/* General method of computation sectio */
    s_wsfe(&io___26);
    e_wsfe();
    s_rsfe(&io___27);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___28);
	e_wsfe();
	s_rsle(&io___29);
	e_rsle();
	goto L300;
    } else {
	goto L300;
    }
/* 26--------USING INTERFACE SUBROUTINE SECTION. */
/* (inum.eq.1) */
L260:
/* Using the interface subroutine secti */
    s_wsfe(&io___30);
    e_wsfe();
    s_rsfe(&io___31);
    do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
    e_rsfe();
    if (inum == 1) {
	s_wsfe(&io___32);
	e_wsfe();
	s_rsfe(&io___33);
	do_fio(&c__1, (char *)&inum, (ftnlen)sizeof(integer));
	e_rsfe();
	if (inum == 1) {
	    s_wsfe(&io___34);
	    e_wsfe();
	    s_rsle(&io___35);
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
    static cilist io___36 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___37 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_2100, 0 };
    static cilist io___40 = { 0, 5, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___42 = { 0, 5, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_2300, 0 };
    static cilist io___44 = { 0, 5, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_2400, 0 };
    static cilist io___46 = { 0, 5, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_2500, 0 };
    static cilist io___48 = { 0, 5, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_2600, 0 };
    static cilist io___50 = { 0, 5, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_2700, 0 };
    static cilist io___52 = { 0, 5, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_2800, 0 };
    static cilist io___54 = { 0, 5, 0, 0, 0 };


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
    s_wsfe(&io___36);
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
    s_rsfe(&io___37);
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
    s_wsfe(&io___39);
    e_wsfe();
    s_rsle(&io___40);
    do_lio(&c__4, &c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    e_rsle();
/* L2101: */
    goto L400;
L220:
/* Change time step limiting const 2 se */
    s_wsfe(&io___41);
    e_wsfe();
    s_rsle(&io___42);
    do_lio(&c__4, &c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L230:
/* Change initial time step section. */
    s_wsfe(&io___43);
    e_wsfe();
    s_rsle(&io___44);
    do_lio(&c__4, &c__1, (char *)&varpr_1.dt1, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L240:
/* Change initial temperature section. */
    s_wsfe(&io___45);
    e_wsfe();
    s_rsle(&io___46);
    do_lio(&c__4, &c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L250:
/* Change final temperature section. */
    s_wsfe(&io___47);
    e_wsfe();
    s_rsle(&io___48);
    do_lio(&c__4, &c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L260:
/* Change smallest abundances allowed s */
    s_wsfe(&io___49);
    e_wsfe();
    s_rsle(&io___50);
    do_lio(&c__4, &c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L270:
/* Change accumulation increment sectio */
    s_wsfe(&io___51);
    e_wsfe();
    s_rsle(&io___52);
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
    s_wsfe(&io___53);
    e_wsfe();
    s_rsle(&io___54);
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
    static cilist io___55 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___56 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___58 = { 0, 6, 0, fmt_2100, 0 };
    static cilist io___59 = { 0, 5, 0, 0, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___61 = { 0, 5, 0, 0, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_2300, 0 };
    static cilist io___63 = { 0, 5, 0, 0, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_2400, 0 };
    static cilist io___65 = { 0, 5, 0, 0, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_2500, 0 };
    static cilist io___67 = { 0, 5, 0, 0, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_2600, 0 };
    static cilist io___69 = { 0, 5, 0, 0, 0 };
    static cilist io___70 = { 0, 6, 0, fmt_2700, 0 };
    static cilist io___71 = { 0, 5, 0, 0, 0 };
    static cilist io___72 = { 0, 6, 0, fmt_2800, 0 };
    static cilist io___73 = { 0, 5, 0, 0, 0 };
    static cilist io___74 = { 0, 6, 0, fmt_2802, 0 };
    static cilist io___75 = { 0, 6, 0, fmt_2804, 0 };
    static cilist io___76 = { 0, 5, 0, 0, 0 };
    static cilist io___77 = { 0, 6, 0, fmt_2900, 0 };
    static cilist io___78 = { 0, 5, 0, 0, 0 };


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
    s_wsfe(&io___55);
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
    s_rsfe(&io___56);
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
    s_wsfe(&io___58);
    e_wsfe();
    s_rsle(&io___59);
    do_lio(&c__4, &c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L220:
/* Change neutron lifetime section. */
    s_wsfe(&io___60);
    e_wsfe();
    s_rsle(&io___61);
    do_lio(&c__4, &c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L230:
/* Change number of neutrino species se */
    s_wsfe(&io___62);
    e_wsfe();
    s_rsle(&io___63);
    do_lio(&c__4, &c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L240:
/* Change baryon-to-photon ratio sectio */
    s_wsfe(&io___64);
    e_wsfe();
    s_rsle(&io___65);
    do_lio(&c__4, &c__1, (char *)&varpr_1.eta1, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L250:
/* Change cosmological constant section */
    s_wsfe(&io___66);
    e_wsfe();
    s_rsle(&io___67);
    do_lio(&c__4, &c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L260:
/* Change neutrino degeneracy section. */
    s_wsfe(&io___68);
    e_wsfe();
    s_rsle(&io___69);
    do_lio(&c__4, &c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L270:
/* Change neutrino degeneracy section. */
    s_wsfe(&io___70);
    e_wsfe();
    s_rsle(&io___71);
    do_lio(&c__4, &c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    e_rsle();
    goto L400;
L280:
/* Change neutrino degeneracy section. */
    s_wsfe(&io___72);
    e_wsfe();
    s_rsle(&io___73);
    do_lio(&c__4, &c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_rsle();
    if (modpr_1.xi[2] != 0.f && modpr_1.c__[2] < 3.f) {
	modpr_1.c__[2] = 3.f;
	s_wsfe(&io___74);
	e_wsfe();
	s_wsfe(&io___75);
	e_wsfe();
	s_rsle(&io___76);
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
    s_wsfe(&io___77);
    e_wsfe();
    s_rsle(&io___78);
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
    static cilist io___81 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___82 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___84 = { 0, 6, 0, fmt_2100, 0 };
    static cilist io___85 = { 0, 5, 0, 0, 0 };
    static cilist io___86 = { 0, 6, 0, fmt_2104, 0 };
    static cilist io___87 = { 0, 5, 0, 0, 0 };
    static cilist io___88 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___89 = { 0, 6, 0, fmt_2202, 0 };
    static cilist io___90 = { 0, 5, 0, 0, 0 };
    static cilist io___91 = { 0, 6, 0, fmt_2300, 0 };
    static cilist io___92 = { 0, 5, 0, 0, 0 };
    static cilist io___100 = { 0, 6, 0, fmt_2302, 0 };
    static cilist io___101 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___102 = { 0, 6, 0, fmt_2304, 0 };
    static cilist io___103 = { 0, 5, 0, 0, 0 };
    static cilist io___104 = { 0, 6, 0, fmt_2306, 0 };
    static cilist io___105 = { 0, 5, 0, 0, 0 };
    static cilist io___106 = { 0, 6, 0, fmt_2308, 0 };
    static cilist io___107 = { 0, 5, 0, 0, 0 };
    static cilist io___108 = { 0, 6, 0, fmt_2310, 0 };
    static cilist io___109 = { 0, 5, 0, 0, 0 };
    static cilist io___110 = { 0, 6, 0, fmt_2312, 0 };
    static cilist io___111 = { 0, 6, 0, fmt_2314, 0 };
    static cilist io___112 = { 0, 6, 0, fmt_2316, 0 };
    static cilist io___113 = { 0, 5, 0, fmt_2301, 0 };
    static cilist io___116 = { 0, 6, 0, fmt_2318, 0 };
    static cilist io___118 = { 0, 6, 0, fmt_2200, 0 };
    static cilist io___125 = { 0, 6, 0, fmt_2202, 0 };
    static cilist io___126 = { 0, 6, 0, fmt_2320, 0 };
    static cilist io___127 = { 0, 5, 0, 0, 0 };



#define vtype_ref(a_0,a_1) &vtype[(a_1)*22 + a_0 - 22]

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
    s_wsfe(&io___81);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___82);
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
    s_wsfe(&io___84);
    e_wsfe();
    s_rsle(&io___85);
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
    s_wsfe(&io___86);
    do_fio(&c__1, (char *)&runopt_1.irun, (ftnlen)sizeof(integer));
    e_wsfe();
    s_rsle(&io___87);
    e_rsle();
    goto L300;
/* 22--------GO SECTION----------------------------------------------------------- */
L220:
    s_wsfe(&io___88);
    e_wsfe();
    checkcb_1.itime = 3;
    check_();
/* Call interface subr before computati */
    driver_();
/* Do nucleosynthesis computation. */
    checkcb_1.itime = 8;
    check_();
/* Call interface subr after computatio */
    s_wsfe(&io___89);
    e_wsfe();
    s_rsle(&io___90);
    e_rsle();
    goto L300;
/* 23--------DO MULTIPLE RUNS SECTION--------------------------------------------- */
/* ..........GET NUMBER OF LOOPINGS. */
L230:
    s_wsfe(&io___91);
    e_wsfe();
    s_rsle(&io___92);
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
	    s_wsfe(&io___100);
	    e_wsfe();
	    s_rsfe(&io___101);
	    do_fio(&c__1, (char *)&inum[i__ - 1], (ftnlen)sizeof(integer));
	    e_rsfe();
	    if (inum[i__ - 1] < 1 || inum[i__ - 1] > 8) {
/* No selection ma */
		s_wsfe(&io___102);
		e_wsfe();
		s_rsle(&io___103);
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
		s_wsfe(&io___104);
		e_wsfe();
		s_rsle(&io___105);
		do_lio(&c__4, &c__1, (char *)&rnum1[i__ - 1], (ftnlen)sizeof(
			real));
		e_rsle();
/* Read in starting value. */
		s_wsfe(&io___106);
		e_wsfe();
		s_rsle(&io___107);
		do_lio(&c__4, &c__1, (char *)&rnum2[i__ - 1], (ftnlen)sizeof(
			real));
		e_rsle();
/* Read in terminating value. */
L232:
		s_wsfe(&io___108);
		e_wsfe();
		s_rsle(&io___109);
		do_lio(&c__4, &c__1, (char *)&rnum3[i__ - 1], (ftnlen)sizeof(
			real));
		e_rsle();
/* Read in incremental value. */
		if (rnum3[i__ - 1] == 0.f) {
/* Trouble with 0 division later on */
		    s_wsfe(&io___110);
		    e_wsfe();
		    goto L232;
		}
		s_wsfe(&io___111);
		do_fio(&c__1, (char *)&rnum1[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&rnum2[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&rnum3[i__ - 1], (ftnlen)sizeof(real));
		e_wsfe();
/* Display inp */
		s_wsfe(&io___112);
		e_wsfe();
		s_rsfe(&io___113);
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
		s_wsfe(&io___116);
		do_fio(&c__1, vtype_ref(0, inum[l - 1]), (ftnlen)22);
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
	s_wsfe(&io___118);
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
		    varpr_1.eta1 = pow_dd(&c_b255, &d__1);
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
			varpr_1.eta1 = pow_dd(&c_b255, &d__1);
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
			    varpr_1.eta1 = pow_dd(&c_b255, &d__1);
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
	s_wsfe(&io___125);
	e_wsfe();
/* Inform user of completion of computa */
    } else {
/* (jnum.eq.0) */
	s_wsfe(&io___126);
	e_wsfe();
    }
/* (jnum.ne.0) */
    s_rsle(&io___127);
    e_rsle();
    goto L300;
/* 24--------EXIT SECTION--------------------------------------------------------- */
L240:
    return 0;
/* 30--------GO BACK TO MENU------------------------------------------------------ */
L300:
    goto L100;
} /* run_ */

#undef vtype_ref
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
    static cilist io___128 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___129 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___131 = { 0, 2, 0, fmt_2000, 0 };
    static cilist io___132 = { 0, 2, 0, fmt_2002, 0 };
    static cilist io___133 = { 0, 2, 0, fmt_2004, 0 };
    static cilist io___134 = { 0, 2, 0, fmt_2006, 0 };
    static cilist io___136 = { 0, 2, 0, fmt_2008, 0 };
    static cilist io___138 = { 0, 2, 0, fmt_2010, 0 };
    static cilist io___139 = { 0, 2, 0, fmt_2012, 0 };
    static cilist io___140 = { 0, 2, 0, fmt_2014, 0 };
    static cilist io___141 = { 0, 6, 0, fmt_2016, 0 };
    static cilist io___142 = { 0, 5, 0, 0, 0 };
    static cilist io___143 = { 0, 6, 0, fmt_3000, 0 };
    static cilist io___144 = { 0, 5, 0, fmt_1001, 0 };
    static cilist io___145 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___146 = { 0, 6, 0, fmt_3100, 0 };
    static cilist io___147 = { 0, 6, 0, fmt_3102, 0 };
    static cilist io___148 = { 0, 6, 0, fmt_3104, 0 };
    static cilist io___149 = { 0, 6, 0, fmt_3106, 0 };
    static cilist io___150 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___151 = { 0, 6, 0, fmt_3108, 0 };
    static cilist io___152 = { 0, 5, 0, 0, 0 };
    static cilist io___153 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___154 = { 0, 6, 0, fmt_3100, 0 };
    static cilist io___155 = { 0, 6, 0, fmt_3102, 0 };
    static cilist io___156 = { 0, 6, 0, fmt_3204, 0 };
    static cilist io___157 = { 0, 6, 0, fmt_3106, 0 };
    static cilist io___158 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___159 = { 0, 6, 0, fmt_3108, 0 };
    static cilist io___160 = { 0, 5, 0, 0, 0 };
    static cilist io___161 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___162 = { 0, 6, 0, fmt_3100, 0 };
    static cilist io___163 = { 0, 6, 0, fmt_3102, 0 };
    static cilist io___164 = { 0, 6, 0, fmt_3304, 0 };
    static cilist io___165 = { 0, 6, 0, fmt_3306, 0 };
    static cilist io___166 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___167 = { 0, 6, 0, fmt_3108, 0 };
    static cilist io___168 = { 0, 5, 0, 0, 0 };
    static cilist io___169 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___170 = { 0, 6, 0, fmt_3100, 0 };
    static cilist io___171 = { 0, 6, 0, fmt_3102, 0 };
    static cilist io___172 = { 0, 6, 0, fmt_3404, 0 };
    static cilist io___173 = { 0, 6, 0, fmt_3406, 0 };
    static cilist io___174 = { 0, 6, 0, fmt_2014, 0 };
    static cilist io___175 = { 0, 6, 0, fmt_3108, 0 };
    static cilist io___176 = { 0, 5, 0, 0, 0 };



#define xout_ref(a_1,a_2) outdat_1.xout[(a_2)*40 + a_1 - 41]
#define thmout_ref(a_1,a_2) outdat_1.thmout[(a_2)*40 + a_1 - 41]

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
    s_wsfe(&io___128);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___129);
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
	s_wsfe(&io___131);
	e_wsfe();
    }
    s_wsfe(&io___132);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___133);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT HEADINGS, ABUNDANCES FOR NEUTRON TO LI8. */
    s_wsfe(&io___134);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___136);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&xout_ref(j, i__), (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
/* ..........PRINT THERMODYNAMIC QUANTITIES. */
    s_wsfe(&io___138);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___139);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.tout[j - 1], (ftnlen)sizeof(real));
	for (i__ = 1; i__ <= 5; ++i__) {
	    do_fio(&c__1, (char *)&thmout_ref(j, i__), (ftnlen)sizeof(real));
	}
	do_fio(&c__1, (char *)&outdat_1.dtout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.etaout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.hubout[j - 1], (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___140);
    e_wsfe();
    outopt_1.outfile = TRUE_;
/* Output file requested. */
    s_wsfe(&io___141);
    e_wsfe();
    s_rsle(&io___142);
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
    s_wsfe(&io___143);
    e_wsfe();
/* ..........READ IN SELECTION NUMBER. */
    s_rsfe(&io___144);
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
    s_wsfe(&io___145);
    e_wsfe();
    s_wsfe(&io___146);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___147);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT HEADINGS, ABUNDANCES FOR D,T,HE3,HE4,LI7. */
    s_wsfe(&io___148);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___149);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	for (i__ = 3; i__ <= 6; ++i__) {
	    do_fio(&c__1, (char *)&xout_ref(j, i__), (ftnlen)sizeof(real));
	}
	do_fio(&c__1, (char *)&xout_ref(j, 8), (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___150);
    e_wsfe();
    s_wsfe(&io___151);
    e_wsfe();
    s_rsle(&io___152);
    e_rsle();
    goto L360;
L320:
/* ..........PRINT CAPTION. */
/* Display n,p,li6,be7,li8&up. */
    s_wsfe(&io___153);
    e_wsfe();
    s_wsfe(&io___154);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___155);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT HEADINGS, ABUNDANCES FOR N,P,LI6,BE7,LI8&UP. */
    s_wsfe(&io___156);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___157);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	for (i__ = 1; i__ <= 2; ++i__) {
	    do_fio(&c__1, (char *)&xout_ref(j, i__), (ftnlen)sizeof(real));
	}
	do_fio(&c__1, (char *)&xout_ref(j, 7), (ftnlen)sizeof(real));
	for (i__ = 9; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&xout_ref(j, i__), (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    s_wsfe(&io___158);
    e_wsfe();
    s_wsfe(&io___159);
    e_wsfe();
    s_rsle(&io___160);
    e_rsle();
    goto L360;
L330:
/* ..........PRINT CAPTION. */
/* Display rhog,rhoe,rhone,rhob. */
    s_wsfe(&io___161);
    e_wsfe();
    s_wsfe(&io___162);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___163);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT ENERGY DENSITIES. */
    s_wsfe(&io___164);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___165);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	for (i__ = 1; i__ <= 4; ++i__) {
	    do_fio(&c__1, (char *)&thmout_ref(j, i__), (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    s_wsfe(&io___166);
    e_wsfe();
    s_wsfe(&io___167);
    e_wsfe();
    s_rsle(&io___168);
    e_rsle();
    goto L360;
L340:
/* ..........PRINT CAPTION. */
/* Display t,dt,phie,eta,hubcst. */
    s_wsfe(&io___169);
    e_wsfe();
    s_wsfe(&io___170);
    do_fio(&c__1, (char *)&compr_1.cy, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ct, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9i, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.t9f, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&compr_1.ytmin, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___171);
    do_fio(&c__1, (char *)&modpr_1.c__[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.c__[2], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.cosmo, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[0], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[1], (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&modpr_1.xi[2], (ftnlen)sizeof(real));
    e_wsfe();
/* ..........PRINT THERMODYNAMIC QUANTITIES. */
    s_wsfe(&io___172);
    e_wsfe();
    i__1 = flags_1.it;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___173);
	do_fio(&c__1, (char *)&outdat_1.t9out[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.tout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&thmout_ref(j, 5), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.dtout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.etaout[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&outdat_1.hubout[j - 1], (ftnlen)sizeof(real));
	e_wsfe();
    }
    s_wsfe(&io___174);
    e_wsfe();
    s_wsfe(&io___175);
    e_wsfe();
    s_rsle(&io___176);
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

#undef thmout_ref
#undef xout_ref


/* Main program alias */ int new123_ () { MAIN__ (); return 0; }
