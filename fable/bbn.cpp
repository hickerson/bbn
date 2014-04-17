#define FEM_TRANSLATION_UNIT_WITH_MAIN

//#include "bbn.hpp"
#include "BigBangModel.hpp"
#include <gsl/gsl_sf_bessel.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>

namespace bbn {

using namespace fem;

//
//========================IDENTIFICATION DIVISION================================
//
void
help(common& cmn)
{
  common_read read(cmn);
  common_write write(cmn);
  const int iw = 6;
  const int ir = 5;
  int inum = 0;
  static const char* format_1001 = "(i1)";
  //
  //----------LINKAGES.
  //     CALLED BY - [program] nuc123
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Displays description and workings of the program.
  //
  //----------PARAMETERS.
  //Input unit number (previous value =
  //Output unit number (previous value =
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------USER RESPONSE VARIABLES.
  //Selection number.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------PRINT HELP SELECTION-------------------------------------------------
  //
  //..........RETURN FROM LOOPING.
  statement_100:
  //..........DISPLAY MENU.
  write(iw,
    "(8(/),32x,"
	"'HELP SELECTION',/,32x,"
	"'---- ---------',/,/,24x,"
    "'1. INTRODUCTION',/,24x,"
	"'2. SETTING UP A RUN',/,24x,"
    "'3. RUNNING THE PROGRAM',/,24x,"
	"'4. OUTPUT OPTIONS',/,24x,"
    "'5. GENERAL METHOD OF COMPUTATION',/,24x,"
    "'6. USING THE INTERFACE SUBROUTINE',/,24x,"
	"'7. EXIT',7(/),24x,"
    "'Enter selection (1-7): ',$)");
  //..........READ IN SELECTION NUMBER.
  read(ir, format_1001), inum;
  //
  //20--------BRANCH TO APPROPRIATE SECTION----------------------------------------
  //
  switch (inum) {
    case 1: goto statement_210;
    case 2: goto statement_220;
    case 3: goto statement_230;
    case 4: goto statement_240;
    case 5: goto statement_250;
    case 6: goto statement_260;
    case 7: goto statement_270;
    default: break;
  }
  //Improper input or <RETURN>.
  goto statement_270;
  //
  //21--------INTRODUCTION SECTION-------------------------------------------------
  //
  //Setting up a run section.
  statement_210:
  write(iw,
    "(/,31x,"
	"'INTRODUCTION',/,31x,"
	"'------------',2(/),"
    "'Welcome to the wonderful world of primordial nucleosynthesis.  NUC123 is a      ',/,"
    "'FORTRAN program designed to provide the early universe researcher with the tools',/,"
    "'necessary for the investigation of primordial nucleosynthesis.  Its menu-driven ',/,"
    "'interface allows the user to first set computation parameters (such as the time ',/,"
    "'step) and model parameters (such as the neutron lifetime and number of neutri-  ',/,"
    "'nos) before doing single runs or multiple runs (in which desired model parame-  ',/,"
    "'ters are varied over a desired range.)  After the run, the user can utilize the ',/,"
    "'menu to either produce an output file or to view the most recent run on the     ',/,"
    "'screen.  The program comes with an empty subroutine CHECK into which the user   ',/,"
    "'may wish to put additional code to add to the computation in an original manner.',10(/),"
    "'(Enter <RETURN> to go back to help menu): ',$)");
  read(ir, star);
  goto statement_300;
  //
  //22--------SET UP RUN SECTION---------------------------------------------------
  //
  //Setting up a run section.
  statement_220:
  write(iw,
    "(/,29x,"
	"'SETTING UP A RUN',/,29x,"
	"'------- -- - ---',2(/),"
    "'I. Setting computation parameters.      										',/,"
    "'   The accuracy of the computation and the relevant temperature region can be	',/,"
    "'   set by the following parameters:     										',/,"
    "'    A. Time step limiting constant 1  (default value of 0.03)					',/,"
    "'    B. Time step limiting constant 2  (default value of 0.003)				',/,"
    "'    C. Initial time step  (default value'," "' of 10**-4)						',/,"
    "'    D. Initial temperature  (default value of 10**2)							',/,"
    "'       This is the temperature at the beginning of the run in units of 10**9 K',/,"
    "'    E. Final temperature  (default value of 10**-2)							',/,"
    "'       This is the termination temperature of the run in units of 10**9 K		',/,"
    "'    F. Smallest abundances allowed  (default value of 10**-25)				',/,"
    "'       Elemental abundances are not allowed to drop below this value			',/,"
    "'    G. # of iterations for each accumulation  (default value of 300)			',/,"
    "'       This is the number of iterations before values are put in an output array',6(/),"
    "'(Enter 1 to continue, <RETURN> to end): ',$)");
  read(ir, format_1001), inum;
  if (inum == 1) {
    write(iw,
      "(/,'II. Setting model parameters.												',/,"
      "'   Default values here give what is known as the standard model with best guess ',/,"
      "'   figure on the neutron lifetime of 889.1 seconds.  Nonstandard scenarios can  ',/,"
      "'   be investigated by varying the following parameters:							',/,"
      "'    A. The gravitational constant       										',/,"
      "'       (The default value of one here gives the usual 6.6720e-8 dyne*cm**2/g**2)',/,"
      "'    B. Neutron life-time  (default value of 889.1 seconds)						',/,"
      "'    C. Number of neutrino species  (default value of 3 light neutrinos)         ',/,"
      "'    D. Final baryon-to-photon ratio  (set to log(eta) = -9.5)                   ',/,"
      "'    E. Cosmological constant  (default value of 0)                              ',/,"
      "'    F. Neutrino degeneracy parameters  (default values all 0)                   ',/,"
      "'       There are 3 separate parameters for the electron, muon, and tau neutrinos',11(/),"
      "'(Enter <RETURN> to go back to help menu): ',$)");
    read(ir, star);
    goto statement_300;
  }
  else {
    goto statement_300;
    //(inum.eq.1)
  }
  //
  //23--------RUN PROGRAM SECTION--------------------------------------------------
  //
  //Running the program section.
  statement_230:
  write(iw,
    "(/,28x,"
	"'RUNNING THE PROGRAM',/,28x,"
	"'------- --- -------',2(/),"
    "' ','I. Setting run speed.                   ',/,"
    "'   The code can be run at 3 different settings of speed.  The running of the    	',/,"
    "'   code can be speeded up by reducing the number of nuclides and reactions.  The	',/,"
    "'   complete computation takes into account the following nuclides: n, p, d, t,  	',/,"
    "'   He3, He4, Li6, Li7, Be7, Li8, B8, Be9, B10, B11, C11, B12, C12, N12, C13, N13,	',/,"
    "'   C14, N14, O14, N15, O15, and O16.    											',/,"
    "'   The given CPU percentages and abundance variations are with regard to a      	',/,"
    "'   single run with all default parameter values.                                	',/,"
    "'    A. 26 nuclides, 88 reactions (default)                                      	',/,"
    "'       nuclides from n to O16           											',/,"
    "'    B. 18 nuclides, 60 reactions        											',/,"
    "'       nuclides from n to N12           											',/,"
    "'       (63% CPU time, variation = .1%)  											',/,"
    "'    C. 9 nuclides, 25 reactions        											',/,"
    "'       nuclides from n to Be7           											',/,"
    "'       (20% CPU time, variation = .5%) ',4(/),"
    "'(Enter 1 to continue, <RETURN> to end): ',$)");
  read(ir, format_1001), inum;
  if (inum == 1) {
    write(iw,
      "(/,"
	  "'II. Do single run.                      										',/,"
      "'    A. Interactive.                     										',/,"
      "'       In an interactive session, the user can readily input the computational  ',/,"
      "'       and model parameters and begin the computation process.  The run itself  ',/,"
      "'       is commenced when option 2, \"GO\", in the \"RUN\" section is requested. ',/,/,"
      "'    B. Batch.                           										',/,"
      "'       To run the program in a batch mode, it must be altered slightly so that  ',/,"
      "'       the I/O takes place with files instead of a terminal.  This is done by   ',/,"
      "'       setting different values for the input and output unit number parameters ',/,"
      "'       \"ir\" and \"iw\" and assigning them to different files in NUC123.  In the   ',/,"
      "'       file assigned the \"ir\" unit number, one must place the responses to the ',/,"
      "'       queries of the program.          			',10(/),"
      "'(Enter 1 to continue, <RETURN> to end): ',$)");
    read(ir, format_1001), inum;
    if (inum == 1) {
      write(iw,
        "(/,'III. Do multiple runs.          											',/,"
        "'   A wide range of early universe models can be covered by doing many runs  	',/,"
        "'   while one or more parameters are varied over a range of interest.  The    	',/,"
        "'   parameters that can be varied are the following:                          	',/,"
        "'    A. Eta                             										',"
        "'       - Logarithmic variation         										',/,"
        "'    B. Gravitational constant          										',"
        "'       - Linear variation              										',/,"
        "'    C. Neutron lifetime                										',"
        "'       - Linear variation              										',/,"
        "'    D. Number of neutrino species      										',"
        "'       - Linear variation              										',/,"
        "'    E. Cosmological constant           										',"
        "'       - Linear variation              										',/,"
        "'    F. Neutrino degeneracy parameters  										',"
        "'       - Linear variation              										',/,"
        "'        1. Electron neutrino           										',/,"
        "'        2. Muon neutrino               										',/,"
        "'        3. Tauon neutrino              										',/,"
        "'   At most 3 parameters can be varied.  The first parameter inputted will be  ',/,"
        "'   will be varied in the outermost loop and the third parameter inputted will ',/,"
        "'   be varied in the innermost loop.    										',7(/),"
        "'(Enter <RETURN> to go back to help menu): ',$)");
      read(ir, star);
      goto statement_300;
    }
    else {
      goto statement_300;
      //(inum.eq.1)
    }
  }
  else {
    goto statement_300;
    //(inum.eq.1)
  }
  //
  //24--------OUTPUT OPTIONS SECTION-----------------------------------------------
  //
  //Output options section.
  statement_240:
  write(iw,
    "(/,30x,"
	"'OUTPUT OPTIONS',/,30x,"
	"'------ -------',2(/),"
    "'I.  Request output file.                											',/,"
    "'   After a run, the user can request the program to put the resulting numbers   	',/,"
    "'   into an output file.  This can be done as many times as desired and all the  	',/,"
    "'   information will be put in one new file under the name of \"NUC123.DAT.\"  If  ',/,"
    "'   there is no request during the entire running of the program, this file is   	',/,"
    "'   not created.  If an output file is requested after a multiple run, only the  	',/,"
    "'   information from the very last run will be given.  The output file will give 	',/,"
    "'   the computational and model parameters for each run and will contain the     	',/,"
    "'   following information:               											',/,"
    "'    A. Temperatures in decreasing order 											',/,"
    "'    B. Abundances for n, p, d, t, He3, He4, Li6, Li7, Be7, and Li8 & up         	',/,"
    "'       (p and He4 are in mass fraction, the rest in ratios to the p abundance)  	',/,"
    "'    C. Time, time interval, chemical potential of the electron                  	',/,"
    "'    D. Energy densities for photons, electrons, electron neutrinos, and baryons 	',/,"
    "'    E. Baryon-to-photon ratio, expansion rate of the universe                   	',5(/),"
    "'(Enter 1 to continue, <RETURN> to end): ',$)");
  read(ir, format_1001), inum;
  if (inum == 1) {
    write(iw,
      "(/,'II.  Request output on screen.         										',/,"
      "'   Instead of waiting to print out an output file, the user can immediately 	',/,"
      "'   access the results of the latest run by requesting the output on the     	',/,"
      "'   screen.  There are four screens on each of which are displayed the       	',/,"
      "'   computational and model parameters and the temperature:                  	',/,"
      "'    A. Abundances for d, t, He3, He4, and Li7                               	',/,"
      "'       (He4 in mass fraction, rest as a ratio with the p abundance)         	',/,"
      "'    B. Abundances for n, p, Li6, Be7, and Li8 & up                          	',/,"
      "'       (p in mass fraction, rest as a ratio with the p abundance)           	',/,"
      "'    C. Energy densities for photons, electrons, electron neutrinos, & baryons 	',/,"
      "'    D. Time, time interval, chemical potential of the electron,             	',/,"
      "'       baryon-to-photon ratio, and expansion rate of the universe           	',11(/),"
      "'(Enter <RETURN> to go back to help menu): ',$)");
    read(ir, star);
    goto statement_300;
  }
  else {
    goto statement_300;
    //(inum.eq.1)
  }
  //
  //25--------METHOD OF COMPUTATION SECTION----------------------------------------
  //
  //General method of computation sectio
  statement_250:
  write(iw,
    "(/,22x,"
	"'GENERAL METHOD OF COMPUTATION',/,22x,"
    "'------- ------ -- -----------',2(/),"
    "'I. Time evolution algorithm.            ',/,"
    "'   The program utilizes a 2-point Runge-Kutta scheme (located in subroutine     ',/,"
    "'   DRIVER) to time-evolve the temperature, the quantity hv (the ratio of the    ',/,"
    "'   baryon density to T**3), the chemical potential of the electron, and the     ',/,"
    "'   nuclide abundances.  In the 2-point Runge-Kutta routine, a variable v at time',/,"
    "'   t0 (= v0) is evolved to a time t1 by adding to v0 the average of the         ',/,"
    "'   derivatives evaluated at t0 and at t1 multiplied by dt:                      ',/,"
    "'       v1 = v0 + 0.5(dvdt(t0)+dvdt(t1)) ',/,"
    "'   where dvdt(t1) is gotten by first finding v1'' = v0 + dvdt(t0).  The         ',/,"
    "'   derivatives of the nuclide abundances are first computed and these are used  ',/,"
    "'   to find the derivatives of t9, hv, and phie (this is done in subroutine      ',/,"
    "'   DERIVS).  To compute the time derivatives of the nuclide abundances, a matrix',/,"
    "'   equation is set up (in subroutine SOL) and is solved (in subroutine EQSLIN)  ',/,"
    "'   by gaussian elimination utilizing implicit differentiation.                  ',6(/),"
    "'(Enter 1 to continue, <RETURN> to end): ',$)");
  read(ir, format_1001), inum;
  if (inum == 1) {
    write(iw,
      "(/,"
	  "'II. Hierarchy of Subroutines.   ',/,"
      "'        NUC123    Main program (main menu)    ',/,"
      "'        HELP      Help option                 ',/,"
      "'        SETCOM    Set computational parameters',/,"
      "'        SETMOD    Set model parameters        ',/,"
      "'        RUN       Run computation code        ',/,"
      "'        DRIVER    Main routine (Runge-Kutta loop)    ',/,"
      "'        START     Initialization routine      ',/,"
      "'        RATE0     Computes weak decay rates   ',/,"
      "'        DERIVS    Computes time derivatives   ',/,"
      "'        THERM     Computes energy densities   ',/,"
      "'        BESSEL    Gives functions of Kn       ',/,"
      "'        KNUX      Computes modified Bessel fcn Kn    ',/,"
      "'        NUDENS    Computes neutrino energy density   ',/,"
      "'        RATE1-4   Computes rates for reactions',/,"
      "'        SOL       Builds A matrix for eqn dy/dt = Ay ',/,"
      "'        EQSLIN    Solves dy/dt=Ay by gaussian elim   ',/,"
      "'        ACCUM     Output accumulator          ',/,"
      "'        OUTPUT    Allows user to output result',4(/),"
      "'(Enter <RETURN> to go back to help menu): ',$)");
    read(ir, star);
    goto statement_300;
  }
  else {
    goto statement_300;
    //(inum.eq.1)
  }
  //
  //26--------USING INTERFACE SUBROUTINE SECTION.
  //
  //Using the interface subroutine secti
  statement_260:
  write(iw,
    "(/,22x,"
	"'USING THE INTERFACE SUBROUTINE',/,22x,"
    "'----- --- --------- ----------',2(/),"
    "'I. Purpose.                             ',/,"
    "'   The interface subroutine CHECK is designed to be an outlet of the program    ',/,"
    "'   into which alterations can be easily plugged.  Programs are normally modified',/,"
    "'   by searching through the program, identifying the appropriate areas for      ',/,"
    "'   alterations, and interspersing new commands while deleting some old ones.    ',/,"
    "'   This process can get tricky unless one actively documents the alterations:   ',/,"
    "'   one might lose track of all of the modifications and deletions.  Thus, it is ',/,"
    "'   worthwhile to put most if not all of the necessary changes into one          ',/,"
    "'   subroutine which is to be called from strategic locations in the main        ',/,"
    "'   program.  Furthermore, by putting changes into one small subroutine, one need',/,"
    "'   only to compile the subroutine CHECK each time instead of the entire nucleo- ',/,"
    "'   synthesis code.                      ',8(/),"
    "'(Enter 1 to continue, <RETURN> to end): ',$)");
  read(ir, format_1001), inum;
  if (inum == 1) {
    write(iw,
      "(/,"
	  "'II. Description.                        ',/,"
      "'   Subroutine CHECK is an empty subroutine with a large COMMON area, giving the ',/,"
      "'   user ready access to all of the important variables in the computations.  The',/,"
      "'   routine is called from various locations in the main program and the location',/,"
      "'   spot in the program is labeled by the flag \"itime\".  The set call locations  ',/,"
      "'   are given below:                     ',/,"
      "'    A. itime = 1 (NUC123, very beginning of program run)                        ',/,"
      "'       (appropriate for opening files, initializing variables)                  ',/,"
      "'    B. itime = 2 (NUC123, right before going into the RUN section)              ',/,"
      "'    C. itime = 3 (RUN, right before going into DRIVER to do the computations)   ',/,"
      "'    D. itime = 4 (DRIVER, in 1st R-K loop after computing derivatives in DERIVS)',/,"
      "'    E. itime = 7 (DRIVER, in 2nd R-K loop after computing derivatives in DERIVS)',/,"
      "'    F. itime = 8 (RUN, right after coming back from DRIVER)                     ',/,"
      "'    G. itime = 9 (NUC123, right after coming back from the RUN section)         ',/,"
      "'    H. itime =10 (NUC123, very end of program run)                              ',/,"
      "'       (appropriate for closing files)  ',/,"
      "'   The difference between the (2,9) pairing and the (3,8) pairing is that for a ',/,"
      "'   multiple run, the (3,8) pairing would be called before and after every run   ',/,"
      "'   but the (2,9) pairing would be called before and after the entire sequence.  ',4(/),"
      "'(Enter 1 to continue, <RETURN> to end): ',$)");
    read(ir, format_1001), inum;
    if (inum == 1) {
      write(iw,
        "(/,"
		"'III. Implementation.                   ',/,"
        "'   The additional program statements are needed in the subroutine CHECK.  If a',/,"
        "'   particular command is to be executed when the computer is at a certain     ',/,"
        "'   location in the program -- say labeled by itime = 8 -- then in CHECK, one  ',/,"
        "'   must place the command under the statement, IF (itime.eq.8)....  The user  ',/,"
        "'   is at leisure to place his own location indicators (5,6) and CALL CHECK    ',/,"
        "'   statements anywhere in the program as long as there is a COMMON /checkcb/    ',/,"
        "'   statement in the particular subroutine to carry the value of itime along.  ',15(/),"
        "'(Enter <RETURN> to go back to help menu): ',$)");
      read(ir, star);
      goto statement_300;
    }
    else {
      goto statement_300;
      //(inum.eq.1)
    }
  }
  else {
    goto statement_300;
    //(inum.eq.1)
  }
  //
  //27--------EXIT SECTION---------------------------------------------------------
  //
  //Exit section.
  statement_270:
  return;
  //
  //30--------GO BACK TO MAIN MENU-------------------------------------------------
  //
  statement_300:
  goto statement_100;
  //
}

//
//========================IDENTIFICATION DIVISION================================
//
void
common::setcom(common& cmn)
{
  common_read read(cmn);
  common_write write(cmn);
  /*
  float& cy = cmn.cy;
  float& ct = cmn.ct;
  float& t9i = cmn.t9i;
  float& t9f = cmn.t9f;
  float& ytmin = cmn.ytmin;
  int& inc = cmn.inc;
  float& dt1 = cmn.dt1;
  */
  //
  const int iw = 6;
  const int ir = 5;
  int inum = 0;
  //
  //----------LINKAGES.
  //     CALLED BY - [program] nuc123
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Allows resetting of computation parameters.
  //
  //----------PARAMETERS.
  //Input unit number (previous value =
  //Output unit number (previous value =
  //
  //----------COMMON AREAS.
  //Default comp param
  //Computation parame
  //Default variationl
  //Variational parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------DEFAULT COMPUTATION PARAMETERS.
  //Default cy.
  //Default ct.
  //Default t9i.
  //Default t9f.
  //Default ytmin.
  //Default accumulation increment.
  //
  //----------COMPUTATION PARAMETERS.
  //Time step limiting constant on abund
  //Time step limiting constant on tempe
  //Initial temperature (in 10**9 K).
  //Final temperature (in 10**9 K).
  //Smallest abundances allowed.
  //Accumulation increment.
  //
  //----------DEFAULT VARIATIONAL  PARAMETERS.
  //Default initial dt.
  //
  //----------VARIATIONAL  PARAMETERS.
  //Initial time step.
  //
  //----------LOCAL VARIABLES.
  //Selection number.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------PRINT RESET SELECTION AND AWAIT RESPONSE-----------------------------
  //
  //..........RETURN FROM LOOPING.
  statement_100:
  //..........DISPLAY RESET SELECTIONS.
  write(iw,
    "(8(/),21x,"
	"'SET COMPUTATION PARAMETERS SELECTION',/,21x,"
    "'--- ----------- ---------- ---------',/,/,10x,"
    "' 1. CHANGE TIME-STEP LIMITING CONSTANT 1  FROM ',f5.3,/,10x,"
    "' 2. CHANGE TIME-STEP LIMITING CONSTANT 2  FROM ',f5.3,/,10x,"
    "' 3. CHANGE INITIAL TIME-STEP              FROM ',1p,e8.2,' SECONDS',/,10x,"
	"' 4. CHANGE INITIAL TEMPERATURE            FROM ',1p,e8.2,' (10**9 K)',/,10x,"
	"' 5. CHANGE FINAL TEMPERATURE              FROM ',1p,e8.2,' (10**9 K)',/,10x,"
    "' 6. CHANGE SMALLEST ABUNDANCES ALLOWED    FROM ',1p,e8.2,/,10x,"
    "' 7. CHANGE ACCUMULATION INCREMENT         FROM ',1p,e8.2,' ITERATIONS',/,10x,"
	"' 8. RESET ALL TO DEFAULT VALUES',/,10x,"
	"' 9. EXIT',5(/),10x,"
	"'Enter selection (1-9): ',$)"),
    cy, ct, dt1, t9i, t9f, ytmin, fem::ffloat(inc);
  //..........READ IN SELECTION NUMBER.
  read(ir, "(i1)"), inum;
  //
  //20--------BRANCH TO APPROPRIATE SECTION----------------------------------------
  //
  switch (inum) {
    case 1: goto statement_210;
    case 2: goto statement_220;
    case 3: goto statement_230;
    case 4: goto statement_240;
    case 5: goto statement_250;
    case 6: goto statement_260;
    case 7: goto statement_270;
    case 8: goto statement_280;
    case 9: goto statement_300;
    default: break;
  }
  //Improper input or <RETURN>.
  goto statement_300;
  //Change time step limiting const 1 se
  statement_210:
  write(iw, "(' ','Enter value for time step limiting constant 1: ',$)");
  read(ir, star), cy;
  goto statement_400;
  //Change time step limiting const 2 se
  statement_220:
  write(iw, "(' ','Enter value for time step limiting constant 2: ',$)");
  read(ir, star), ct;
  goto statement_400;
  //Change initial time step section.
  statement_230:
  write(iw, "(' ','Enter value for initial time step: ',$)");
  read(ir, star), dt1;
  goto statement_400;
  //Change initial temperature section.
  statement_240:
  write(iw, "(' ','Enter value for initial temperature: ',$)");
  read(ir, star), t9i;
  goto statement_400;
  //Change final temperature section.
  statement_250:
  write(iw, "(' ','Enter value for final temperature: ',$)");
  read(ir, star), t9f;
  goto statement_400;
  //Change smallest abundances allowed s
  statement_260:
  write(iw, "(' ','Enter value for smallest abundances allowed: ',$)");
  read(ir, star), ytmin;
  goto statement_400;
  //Change accumulation increment sectio
  statement_270:
  write(iw, "(' ','Enter value for accumulation increment: ',$)");
  read(ir, star), inc;
  goto statement_400;
  //Reset all to default values section.
  statement_280:
  //Time step limiting constant on abundance
  cy = cmn.cy0;
  //Time step limiting constant on temperature
  ct = cmn.ct0;
  //Time step.
  dt1 = cmn.dt0;
  //Initial temperature.
  t9i = cmn.t9i0;
  //Final temperature.
  t9f = cmn.t9f0;
  //Smallest abundances allowed.
  ytmin = cmn.ytmin0;
  //Accumulation increment.
  inc = cmn.inc0;
  write(iw,
    "(' ','All values reset to default - Press <RETURN> ','to continue: ',$)");
  read(ir, star);
  goto statement_400;
  //Exit section.
  statement_300:
  return;
  //
  //40--------GO BACK TO MENU------------------------------------------------------
  //
  statement_400:
  goto statement_100;
  //
}

//
//========================IDENTIFICATION DIVISION================================
//
void
common::setmod(
  common& cmn)
{
  common_read read(cmn);
  common_write write(cmn);
/*
  arr_cref<float> c0(cmn.c0, dimension(3));
  arr_cref<float> xi0(cmn.xi0, dimension(3));
  arr_ref<float> c(cmn.c, dimension(3));
  float& cosmo = cmn.cosmo;
  arr_ref<float> xi(cmn.xi, dimension(3));
  float& eta1 = cmn.eta1;
  */
  //
  const int iw = 6;
  const int ir = 5;
  int inum = 0;
  //
  //----------LINKAGES.
  //     CALLED BY - [program] nuc123
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Allows resetting of model parameters.
  //
  //----------PARAMETERS.
  //Input unit number (previous value =
  //Output unit number (previous value =
  //
  //----------COMMON AREAS.
  //Default model para
  //Model parameters.
  //Default variationl
  //Variational parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------DEFAULT MODEL PARAMETERS.
  //Default c.
  //Default cosmological constant.
  //Default neutrino degeneracy parameter
  //
  //----------EARLY UNIVERSE MODEL PARAMETERS.
  //c(1) is variation of gravitational con
  //c(2) is neutron lifetime (sec).
  //c(3) is number of neutrino species.
  //Cosmological constant.
  //Neutrino degeneracy parameters.
  //
  //----------DEFAULT VARIATIONAL PARAMETERS.
  //Default eta.
  //
  //----------VARIATIONAL PARAMETERS.
  //Initial baryon-to-photon ratio.
  //
  //----------USER RESPONSE VARIABLES.
  //Selection number.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------PRINT RESET SELECTION AND AWAIT RESPONSE-----------------------------
  //
  //..........RETURN FROM LOOPING.
  statement_100:
  //..........DISPLAY RESET SELECTIONS.
  write(iw,
    "(8(/),24x,'SET MODEL PARAMETERS SELECTION',/,24x,"
    "'--- ----- ---------- ---------',/,/,10x,"
    "' 1. CHANGE GRAVITATIONAL CONSTANT         FROM ',1p,e10.3,/,10x,"
    "' 2. CHANGE NEUTRON LIFETIME               FROM ',1p,e10.3,' SECONDS',/,10x,"
	"' 3. CHANGE NUMBER OF NEUTRINO SPECIES     FROM ',1p,e10.3,/,10x,"
	"' 4. CHANGE FINAL BARYON-TO-PHOTON RATIO   FROM ',1p,e10.3,/,10x,"
    "' 5. CHANGE COSMOLOGICAL CONSTANT          FROM ',1p,e10.3,/,10x,"
    "' 6. CHANGE XI-ELECTRON                    FROM ',1p,e10.3,/,10x,"
    "' 7. CHANGE XI-MUON                        FROM ',1p,e10.3,/,10x,"
    "' 8. CHANGE XI-TAUON                       FROM ',1p,e10.3,/,10x,"
    "' 9. RESET ALL TO DEFAULT VALUES',/,10x,'10. EXIT',4(/),10x,"
    "' Enter selection (1-10): ',$)"),
    cmn.c(1), cmn.c(2), cmn.c(3), cmn.eta1, cmn.cosmo, cmn.xi(1), cmn.xi(2), cmn.xi(3);
  //..........READ IN SELECTION NUMBER.
  read(ir, "(i2)"), inum;
  //
  //20--------BRANCH TO APPROPRIATE SECTION----------------------------------------
  //
  switch (inum) {
    case 1: goto statement_210;
    case 2: goto statement_220;
    case 3: goto statement_230;
    case 4: goto statement_240;
    case 5: goto statement_250;
    case 6: goto statement_260;
    case 7: goto statement_270;
    case 8: goto statement_280;
    case 9: goto statement_290;
    case 10: goto statement_300;
    default: break;
  }
  //Improper input or <RETURN>.
  goto statement_300;
  //Change gravitational constant section.
  statement_210:
  write(iw,
    "(' ','Enter value for variation of gravitational ','constant: ',$)");
  read(ir, star), cmn.c(1);
  goto statement_400;
  //Change neutron lifetime section.
  statement_220:
  write(iw, "(' ','Enter value for neutron lifetime (sec): ',$)");
  read(ir, star), cmn.c(2);
  goto statement_400;
  //Change number of neutrino species section.
  statement_230:
  write(iw, "(' ','Enter value for number of neutrino species: ',$)");
  read(ir, star), cmn.c(3);
  goto statement_400;
  //Change baryon-to-photon ratio section.
  statement_240:
  write(iw, "(' ','Enter value for baryon-to-photon ratio: ',$)");
  read(ir, star), cmn.eta1;
  goto statement_400;
  //Change cosmological constant section.
  statement_250:
  write(iw, "(' ','Enter value for cosmological constant: ',$)");
  read(ir, star), cmn.cosmo;
  goto statement_400;
  //Change neutrino degeneracy section.
  statement_260:
  write(iw, "(' ','Enter value for xi electron: ',$)");
  read(ir, star), cmn.xi(1);
  goto statement_400;
  //Change neutrino degeneracy section.
  statement_270:
  write(iw, "(' ','Enter value for xi muon: ',$)");
  read(ir, star), cmn.xi(2);
  goto statement_400;
  //Change neutrino degeneracy section.
  statement_280:
  write(iw, "(' ','Enter value for xi tauon: ',$)");
  read(ir, star), cmn.xi(3);
  if ((cmn.xi(3) != 0.f) && (cmn.c(3) < 3.f)) {
    cmn.c(3) = 3.f;
    write(iw, "(' ','Number of neutrinos set to 3')");
    write(iw, "(' ','Press <RETURN> to continue: ',$)");
    read(ir, star);
  }
  goto statement_400;
  //Reset all to default values section.
  statement_290:
  cmn.c(1) = cmn.c0(1);
  cmn.c(2) = cmn.c0(2);
  cmn.c(3) = cmn.c0(3);
  cmn.cosmo = cmn.cosmo0;
  cmn.xi(1) = cmn.xi0(1);
  cmn.xi(2) = cmn.xi0(2);
  cmn.xi(3) = cmn.xi0(3);
  cmn.eta1 = cmn.eta0;
  write(iw,
    "(' ','All values reset to default - Press <RETURN> ','to continue: ',$)");
  read(ir, star);
  goto statement_400;
  //Exit section.
  statement_300:
  return;
  //
  //40--------GO BACK TO MENU------------------------------------------------------
  //
  statement_400:
  goto statement_100;
  //
}

/*
struct check_save
{
  fem::variant_bindings lncoef_bindings;
  fem::variant_bindings rates_bindings;
};
*/

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//     Changes (to run inder unix f77):
//     -------------------------------
//     Programme renamed from nucint.for to newint.f
//     COMMON /therm/  -> COMMON /thermcb/
//     COMMON /bessel/ -> COMMON /besselcb/
//     COMMON /check/  -> COMMON /checkcb/
//     COMMON /time/   -> COMMON /ttime/
//
//     Updated Correction to helium-4 abundance (using fitted rates and
//     *smallest* step-size) - Sarkar, Rep. Prog Phys. 59, 1493 (1996):
//
//     Y_p = Y_p - 0.0025 -> Y_p = Y_p - 0.0003
//
//     output nucint.dat  -> newint.dat
//
//========================IDENTIFICATION DIVISION================================
//
void
common::check(common& cmn)
{
  //FEM_CMN_SVE(check);
  common_write write(cmn);
  // COMMON modpr
  //arr_cref<float> c(cmn.c, dimension(3));
  // COMMON flags
  //int& it = cmn.it;
  // COMMON checkcb
  //int& itime = static_cast<common_checkcb&>(cmn).itime;
  // COMMON outdat
  const int itmax = 40;
  const int nnuc = 26;
  const int nrec = 88;
  //arr_ref<float, 2> xout(cmn.xout, dimension(itmax, nnuc));
  //arr_cref<float> etaout(cmn.etaout, dimension(itmax));
  //
  /*
  common_variant rates(cmn.common_rates, sve.rates_bindings);
  common_variant lncoef(cmn.common_lncoef, sve.lncoef_bindings);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r(dimension(nrec));
      rates.allocate(), f, r;
    }
    {
      mbr<double> a(dimension(nnuc, nnuc));
      mbr<float> b(dimension(nnuc));
      mbr<float> yx(dimension(nnuc));
      lncoef.allocate(), a, b, yx;
    }
  }
  */
  //* arr_cref<float> f( */ rates.bind<float>() /* , dimension(nrec)) */ ;
  //* arr_cref<float> r( */ rates.bind<float>() /* , dimension(nrec)) */ ;
  //* arr_cref<double, 2> a( */ lncoef.bind<double>() /* , dimension(nnuc, nnuc)) */ ;
  //* arr_cref<float> b( */ lncoef.bind<float>() /* , dimension(nnuc)) */ ;
  //* arr_cref<float> yx( */ lncoef.bind<float>() /* , dimension(nnuc)) */ ;
  //
  //----------REMARKS.
  //     This is an interface subroutine,
  //     a flexible module which allows user to manipulate physical quantities
  //     of interest at certain key points during the computer run.
  //     Included within this subroutine is a roster of all global variables
  //     and their respective COMMON areas.
  //     Current instructions accumulate abundances of deuterium, helium-3,
  //     helium-4, and lithium-7 for eventual plotting, taking into account
  //     the contribution of beryllium-7 to lithium-7 and tritium to helium-3.
  //
  //----------PARAMETERS.
  //Number of nuclear reactions.
  //Number of variables to be evolved.
  //Number of nuclides in calculation.
  //Maximum # of lines to be printed.
  //
  //----------COMMON AREAS.
  //Reaction parameter
  //Reaction rates.
  //Evolution parameter.
  //Default comp parameter.
  //Computation parameter.
  //Default model para
  //Model parameters.
  //Default variational parameter.
  //Variational parameter.
  //Time variables.
  //Dynamic variables.
  //Energy densities.
  //Linear eqn coefficient
  //Nuclide data.
  //Eval function bl
  //Eval function bm(z
  //Eval function bn(z
  //Coefficients K.
  //Flags,counters.
  //Computation location.
  //Output data.
  //Neutrino parameter.
  //Run options.
  //Output option.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION PARAMETER VALUES.
  //Reaction parameters.
  //
  //----------REACTION PARAMETER NAMES.
  //Reaction type code (1-11).
  //Incoming nuclide type (1-26).
  //Incoming light nuclide type (1-6).
  //Outgoing light nuclide type (1-6).
  //Outgoing nuclide type (1-26).
  //Reverse reaction coefficient.
  //Energy released in reaction (in 10**
  //
  //----------REACTION RATES.
  //Forward reaction rate coefficients.
  //Reverse reaction rate coefficients.
  //
  //----------EVOLUTION PARAMETERS.
  //Temperature of photons (units of 10*
  //Defined by hv = M(atomic)n(baryon)/t
  //Chemical potential of electron.
  //Relative number abundances.
  //
  //----------EVOLUTION PARAMETERS (DERIVATIVES).
  //Change in temperature.
  //Change in hv.
  //Change in chemical potential.
  //Change in relative number abundances
  //
  //----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
  //Rel # abundances at end of 1st R-K l
  //
  //----------DEFAULT COMPUTATION PARAMETERS.
  //Default cy.
  //Default ct.
  //Default t9i.
  //Default t9f.
  //Default ytmin.
  //Default accumulation increment.
  //
  //----------COMPUTATION PARAMETERS.
  //Time step limiting constant on abund
  //Time step limiting constant on tempe
  //Initial temperature (in 10**9 K).
  //Final temperature (in 10**9 K).
  //Smallest abundances allowed.
  //Accumulation increment.
  //
  //----------DEFAULT MODEL PARAMETERS.
  //Default constant.
  //Default cosmological constant.
  //Default neutrino degeneracy parameter.
  //
  //----------EARLY UNIVERSE MODEL PARAMETERS.
  //Gravitational constant.
  //Neutron lifetime (sec).
  //Number of neutrino species.
  //c(1) is variation of gravitational constant.
  //c(2) is neutron half-life (min).
  //c(3) is number of neutrino species.
  //Cosmological constant.
  //Neutrino degeneracy parameters.
  //xi(1) is e neutrino degeneracy parameter.
  //xi(2) is m neutrino degeneracy parameter.
  //xi(3) is t neutrino degeneracy parameter.
  //
  //----------DEFAULT VARIATIONAL PARAMETERS.
  //Default initial time step.
  //Default baryon-to-photon ratio.
  //
  //----------VARIATIONAL PARAMETERS.
  //Initial time step.
  //Baryon-to-photon ratio.
  //
  //----------TIME VARIABLES.
  //Time.
  //Time step.
  //(1/t9)*d(t9)/d(t).
  //
  //----------DYNAMIC VARIABLES.
  //Thermodynamic variables (energy density).
  //Expansion rate of the universe.
  //
  //----------ENERGY DENSITIES.
  //Initial electron neutrino energy density.
  //Initial baryon energy density.
  //Baryon energy density.
  //Baryon energy density (ratio to initi
  //
  //----------MATRIX COEFFICIENTS FOR LINEAR EQUATION.
  //Relates y(t+dt) to y(t).
  //Contains y0 in inverse order.
  //yy in reverse order.
  //
  //----------NUCLIDE DATA.
  //Atomic number of nuclide.
  //Charge of nuclide.
  //Mass excess of nuclide.
  //
  //----------EVALUATION OF FUNCTIONS bl,bm,bn.
  //Evaluation of function bl(z).
  //Evaluation of function bm(z).
  //Evaluation of function bn(z).
  //
  //----------EVALUATION OF MODIFIED BESSEL FUNCTIONS.
  //Values k0(r),k1(r),k2(r),k3(r),k4(r)
  //
  //----------FLAGS AND COUNTERS.
  //Indicates if output buffer printed.
  //# total iterations for particular mo
  //# iterations after outputing a line.
  //# times accumulated in output buffer
  //Indicates if gaussian elimination fa
  //
  //----------COMPUTATION LOCATION.
  //Time check.
  //
  //----------OUTPUT ARRAYS.
  //Nuclide mass fractions.
  //Thermodynamic variables.
  //Temperature (in units of 10**9 K).
  //Time.
  //Time step.
  //Baryon to photon ratio.
  //Expansion rate.
  //
  //----------NEUTRINO PARAMETERS.
  //Temperature (in units of MeV).
  //Neutrino temperature (in units of Me
  //Neutrino temperature.
  //Normalizing constant.
  //Neutrino energy density.
  //Type of neutrino.
  //
  //----------RUN OPTION.
  //Run network size.
  //Number of nuclides in computation.
  //Number of reactions in computation.
  //
  //----------OUTPUT FILE STATUS.
  //Number of output requests.
  //Indicates if output file used.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------OPEN FILE------------------------------------------------------------
  //
  //Beginning of program.
  if (itime == 1) {
    cmn.io.open(3, "newint.dat").status("new");
  }
  //
  //20--------WRITE INTO FILE------------------------------------------------------
  //
  //Right after a run.
  if (itime == 8) {
    //Add beryllium to lithium.
    xout(it, 8) += xout(it, 9);
    //Add tritium to helium-3.
    xout(it, 5) += xout(it, 4);
    xout(it, 6) = xout(it, 6) - 0.0003f;
    // my correction for fitted rates+coarse steps
    //Output N_nu, tau_n, eta, H2, He3, He4, an
    write(3, "(7(e13.5))"), c(3), c(2), etaout(it), xout(it, 3),
      xout(it, 5), xout(it, 6), xout(it, 8);
  }
  //
  //30--------CLOSE FILE-----------------------------------------------------------
  //
  //End of program.
  if (itime == 10) {
    cmn.io.close(3);
  }
  //
  //----------REFERENCES-----------------------------------------------------------
  //     1) D.A. Dicus, E.W. Kolb, A.M. Gleeson, E.C.G. Sudarshan, V.L. Teplitz,
  //        M.S. Turner, Phys. Rev. D26 (1982) 2694. (Rad corr, Coulomb Corr)
  //     2) D. Seckel, Bartol preprint BA-93-16; G. Guyk and M.S. Turner,
  //         FERMILAB preprint FERMILAB-Pub-93/181-A. (Nucleon mass)
  //     4) S. Dodelson and M.S. Turner, Phys. Rev. D46 (1992) 3372; B. Fields,
  //         S. Dodelson and M.S. Turner, Phys. Rev. D47 (1993) 4309. (Nu heating)
  //     5) S. Sarkar, Rep. Prog Phys. 59 (1996) 1493 (review)
  //
}

//
//========================IDENTIFICATION DIVISION================================
//
float
ex(
  float const& x)
{
  float return_value = 0;
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] start, rate2, rate3, rate4, sol
  //               - [function] eval
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Exponential function with underflow precaution.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //In danger of overflow.
  if (x > 88.029f) {
    return_value = exp(88.029f);
  }
  else {
    //In danger of underflow.
    if (x <  - 88.722f) {
      return_value = 0.f;
      //Value of x in allowable range.
    }
    else {
      return_value = exp(x);
    }
  }
  return return_value;
  //
  //----------NOTE-----------------------------------------------------------------
  //     The overflow limit for the VAX/VMS system is exp(88.029).
  //     The underflow limit for the VAX/VMS system is exp(-88.722).
  //
}

/*
struct knux_save
{
  arr<float> c0;
  arr<float> c1;
  arr<float> ci0;
  arr<float> ci1;
  arr<float> ck0;
  arr<float> ck1;

  knux_save() :
    c0(dimension(7), fem::fill0),
    c1(dimension(7), fem::fill0),
    ci0(dimension(7), fem::fill0),
    ci1(dimension(7), fem::fill0),
    ck0(dimension(7), fem::fill0),
    ck1(dimension(7), fem::fill0)
  {}
};
*/

//
//========================IDENTIFICATION DIVISION================================
//
#if 0
void
knux(
  common& cmn,
  float const& z)
{
  FEM_CMN_SVE(knux);
  // COMMON kays
  float& bk0 = cmn.bk0;
  float& bk1 = cmn.bk1;
  float& bk2 = cmn.bk2;
  float& bk3 = cmn.bk3;
  //
  // SAVE
  arr_ref<float> c0(sve.c0, dimension(7));
  arr_ref<float> c1(sve.c1, dimension(7));
  arr_ref<float> ci0(sve.ci0, dimension(7));
  arr_ref<float> ci1(sve.ci1, dimension(7));
  arr_ref<float> ck0(sve.ck0, dimension(7));
  arr_ref<float> ck1(sve.ck1, dimension(7));
  //
  if (is_called_first_time) {
    {
      static const float values[] = {
        1.f, 3.5156229f, 3.0899424f, 1.2067492f, 0.2659732f,
          0.0360768f, 0.0045813f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        ci0;
    }
    {
      static const float values[] = {
        0.5f, 0.87890594f, 0.51498869f, 0.15084934f, 0.02658733f,
          0.00301532f, 0.00032411f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        ci1;
    }
    {
      static const float values[] = {
        -0.57721566f, 0.42278420f, 0.23069756f, 0.03488590f,
          0.00262698f, 0.00010750f, 0.00000740f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        ck0;
    }
    {
      static const float values[] = {
        1.f, 0.15443144f, -0.67278579f, -0.18156897f, -0.01919402f,
          -0.00110404f, -0.00004686f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        ck1;
    }
    {
      static const float values[] = {
        1.25331414f, -0.07832358f, 0.02189568f, -0.01062446f,
          0.00587872f, -0.00251540f, 0.00053208f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        c0;
    }
    {
      static const float values[] = {
        1.25331414f, 0.23498619f, -0.03655620f, 0.01504268f,
          -0.00780353f, 0.00325614f, -0.00068245f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        c1;
    }
  }
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] bessel
  //     CALLS     - [function] exp
  //
  //----------REMARKS.
  //     A subroutine for modified bessel functions of the second kind
  //     k-nu(z).
  //
  //----------COMMON AREAS.
  //Coefficients K.
  //
  //===========================DECLARATION DIVISION================================
  //
  //-----------MODIFIED BESSEL FUNCTION VALUES.
  //Values k0(z),k1(z)
  //Values i0(z),i1(z).
  //Values k2(z),k3(z),k4(z).
  //
  //-----------EXPANSION COEFFICIENTS.
  //Expansion coefficients for i0 (z.le.2).
  //Expansion coefficients for i1 (z.le.2).
  //Expansion coefficients for k0 (z.le.2).
  //Expansion coefficients for k1 (z.le.2).
  //Expansion coefficients for k0 (z.gt.2).
  //Expansion coefficients for k1 (z.gt.2).
  //
  //-----------VARIABLES TO BE EVALUATED.
  //Input variable.
  //Expansion variable = z/2.
  //Expansion variable = z/3.75.
  //Logarithmic or exponential coefficient
  //
  //==============================DATA DIVISION====================================
  //
  //----------EXPANSION COEFFICIENTS.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------COMPUTE K0 AND K1----------------------------------------------------
  //
  //(Ref. 1).
  float t = 0;
  float y = 0;
  float coeff = 0;
  float bi0 = 0;
  float bi1 = 0;
  int i = 0;
  if (z <= 2.f) {
    //..........COMPUTE FACTORS.
    t = (z / 3.75f);
    y = (z / 2);
    coeff = fem::alog(y);
    //..........VALUES FOR i0(z) and i1(z).
    bi0 = ci0(1);
    bi1 = ci1(1);
    bk0 = ck0(1);
    bk1 = ck1(1);
    FEM_DO_SAFE(i, 2, 7) {
      int n = 2 * (i - 1);
      float p = pow(t,n);
      bi0 += ci0(i) * p;
      bi1 += ci1(i) * p;
      bk0 += ck0(i) * p;
      bk1 += ck1(i) * p;
    }
    //..........VALUES FOR k0(z) and k1(z).
    bk0 += -coeff * bi0;
    bk1 = coeff * bi1 * z + bk1 / z;
    //(z.le.2.)               !(Ref. 2).
  }
  else {
    //..........COMPUTE FACTORS.
    y = (2.0f / z);
    coeff = (ex(-z) / sqrt(z));
    //..........VALUES FOR k0(z) and k1(z).
    bk0 = c0(1);
    bk1 = c1(1);
    FEM_DO_SAFE(i, 2, 7) {
      bk0 += c0(i) * pow(y, (i - 1));
      bk1 += c1(i) * pow(y, (i - 1));
    }
    bk0 = coeff * bk0;
    bk1 = coeff * bk1;
    //(z.le.2.)
  }
  //
  //20--------FIND K2, K3, AND K4 BY ITERATION (Ref. 3)----------------------------
  //
  //k2(z).
  bk2 = 2.f * (bk1 / z) + bk0;
  //k3(z).
  bk3 = 4.f * (bk2 / z) + bk1;
  //k4(z).
  cmn.bk4 = 6.f * (bk3 / z) + bk2;
  //
  //----------REFERENCES-----------------------------------------------------------
  //     Handbook of Mathematical Functions (Abramowitz and Stegun),
  //       Dover Publications, Inc., New York
  //       1) Polynomial approximations for z.le.2
  //         page 378, equations 9.8.1 and 9.8.3.
  //         page 379, equations 9.8.5 and 9.8.7.
  //       2) Polynomial approximations for z > 2
  //         page 379, equations 9.8.6 and 9.8.8.
  //       3) Recursion relation from 1st line of 9.6.26, page 376.
  //
}
#endif

//
//========================IDENTIFICATION DIVISION================================
//
#if 0
void
bessel(
  common& cmn,
  float const& z)
{
  // COMMON kays
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] start, therm
  //     CALLS     - [subroutine] knux
  //
  //----------REMARKS.
  //     Evaluates functions bl(z), bm(z), and bn(z) using solutions to
  //     modified Bessel functions.
  //
  //----------COMMON AREAS.
  //Evaluate function bl(z).
  //Evaluate function bm(z).
  //Evaluate function bn(z).
  //Coefficients K.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------EVALUATION OF FUNCTIONS bl,bm,bn.
  //Single variables equivalenced to arr
  //
  //----------EVALUATION OF MODIFIED BESSEL FUNCTIONS.
  //Values k0(r),k1(r),k2(r),k3(r),k4(r).
  //
  //----------LOCAL VARIABLES.
  //Array containing values from function bl.
  //Defined by z = m(electron)*c**2/k*t9
  //Multiples of z.
  //
  //----------EQUIVALENCE STATEMENTS.
  //     EQUIVALENCE (blz(1),bl1),(blz(2),bl2),(blz(3),bl3),(blz(4),bl4),
  //    |            (blz(5),bl5)
  //     EQUIVALENCE (bmz(1),bm1),(bmz(2),bm2),(bmz(3),bm3),(bmz(4),bm4),
  //    |            (bmz(5),bm5)
  //     EQUIVALENCE (bnz(1),bn1),(bnz(2),bn2),(bnz(3),bn3),(bnz(4),bn4),
  //    |            (bnz(5),bn5)
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------LOCALLY DEFINED FUNCTIONS--------------------------------------------
  //
  //Function bl.
  //bl(z) = bk2 / z;
  //Function bm.
  //bm(z) = 0.25f * (3.f * cmn.bk3 + cmn.bk1) / z;
  //Function bn.
  //bn(z) = 0.5f * (cmn.bk4 + bk2) / z;
  //
  //20--------CALCULATE FOR 1 THRU 5 Z---------------------------------------------
  //
  int i = 0;
  arr_1d<5, float> blz(fem::fill0); //TODO change to c array
  arr_1d<5, float> bmz(fem::fill0); //TODO change to c array
  arr_1d<5, float> bnz(fem::fill0); //TODO change to c array
  FEM_DO_SAFE(i, 1, 5) {
    float r = i * z; 							//Multiples of z.
    knux(cmn, r); 								//Get k0(r),k1(r),k2(r),k3(r),k4(r),k(r)
    blz(i) = cmn.bk2/r; 					//Put value from function bl into array
    bmz(i) = 0.25*(3*cmn.bk3 + cmn.bk1)/r; 	//Put value from function bm into array
    bnz(i) = 0.5*(cmn.bk4 + cmn.bk2)/r; 	//Put value from function bn into array
  }

  // TODO remove this fix 
}
#endif

float getBesselL(float r)
{
  	float K2r = gsl_sf_bessel_Kn(2,r); 		/// Irregular modified cylindrical Bessel functions.
	return K2r/r;
}

float getBesselM(float r)
{
  	const float K3r = gsl_sf_bessel_Kn(3,r);/// Irregular modified cylindrical Bessel functions.
  	const float K1r = gsl_sf_bessel_Kn(1,r);/// Irregular modified cylindrical Bessel functions.
    return 0.25*(3*K3r + K1r)/r;			/// (Ref ?).
}

float getBesselN(float r)
{
  	const float K4r = gsl_sf_bessel_Kn(4,r);/// Irregular modified cylindrical Bessel functions.
  	const float K2r = gsl_sf_bessel_Kn(2,r);/// Irregular modified cylindrical Bessel functions.
    return 0.5*(K4r + K2r)/r; 				/// (Ref ?).
}


/*
struct rate0_save
{
  fem::variant_bindings rates_bindings;
};
*/

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//     Changes (to run inder DEC unix f77):
//     -----------------------------------
//     COMMON /therm/ -> COMMON /thermcb/
//     COMMON /rates/ f,r(nrec) -> COMMON /rates/ f,r
//
//     Default neutron lifetime 888.54 -> 885.7 (+- 0.8) - PDG 2002/3
//
void
rate0(
  common& cmn)
{
  //FEM_CMN_SVE(rate0);
  //common_variant rates(cmn.common_rates, sve.rates_bindings);
  const int nrec = 88;
  arr<float> f(dimension(nrec));
  /*
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r(dimension(nrec));
      rates.allocate(), f, r;
    }
  }
  */
  //arr_ref<float> f(rates.bind<float>(), dimension(nrec));
  //arr_cref<float> r( *//* rates.bind<float>() *//* , dimension(nrec)) *//* ;
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] start
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Generates weak decay rates.
  //
  //----------PARAMETER.
  //Number of nuclear reactions.
  //
  //----------COMMON AREA.
  //Reaction rates.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION RATE COEFFICIENTS.
  //Forward reaction rate coefficients.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------SET DECAY RATE COEFFICIENTS------------------------------------------
  //
  //.......H3 -> e- + v + He3.........(Tilly-Weller-Hasan 1987)
  f(2) = 1.79e-9f;
  //
  //.......Li8 -> e- + v + 2He4.......(Ajzenberg-Selove 1988)
  f(3) = 8.27e-1f;
  //
  //.......B12 -> e- + B + C12........(Ajzenberg-Selove 1990)
  f(4) = 3.43e+1f;
  //
  //.......C14 -> e- + v + N14........(Ajzenberg-Selove 1986)
  f(5) = 3.834e-12f;
  //
  //.......B8 -> e+ + v + 2He4........(Ajzenberg-Selove 1988)
  f(6) = 9.00e-1f;
  //
  //.......C11 -> e+ + v + B11........(Ajzenberg-Selove 1990)
  f(7) = 5.668e-4f;
  //
  //.......N12 -> e+ + v + C12........(Ajzenberg-Selove 1990)
  f(8) = 6.301e+1f;
  //
  //.......N13 -> e+ + v + C13........(Ajzenberg-Selove 1986)
  f(9) = 1.159e-3f;
  //
  //.......O14 -> e+ + v + N14........(Ajzenberg-Selove 1986)
  f(10) = 9.8171e-3f;
  //
  //.......O15 -> e+ + v + N15........(Ajzenberg-Selove 1986)
  f(11) = 5.6704e-3f;
  //
  //----------REFERENCES-----------------------------------------------------------
  //     Ajzenberg-Selove, F., 1990, Nucl. Phys. A506, 1.
  //     Ajzenberg-Selove, F., 1988, Nucl. Phys. A490, 1.
  //     Ajzenberg-Selove, F., 1986, Nucl. Phys. A449, 1.
  //     Tilley, D.R., Weller, H.R., and Hasan, H.H., 1987, Nucl. Phys. A474, 1.
  //
}

//typedef float (*func1_function_pointer)(common&, float const&);

//
//========================IDENTIFICATION DIVISION================================
//
//===========================PROCEDURE DIVISION==================================
//
//10--------1ST PART OF INTEGRAL FOR n->p RATE-----------------------------------
//***************************************************************************
//
float
common::func1(
  common& cmn,
  float const& x)
{
  float return_value = 0;
  // COMMON modpr
  //arr_cref<float> xi(cmn.xi, dimension(3));
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] rate1, nudens
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Contains integrands to be integrated.
  //
  //----------COMMON AREAS.
  //Model parameters.
  //Integration parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------MODEL PARAMETERS.
  //Neutrino degeneracy parameters.
  //
  //----------NEUTRINO PARAMETERS.
  //Temperature (in units of MeV).
  //Neutrino temperature (in units of Me
  //Neutrino temperature (units of 10**9
  //Normalizing constant.
  //
  //----------LOCAL VARIABLES.
  //Value at which function is evaluated
  //Exponential expression with photon t
  //Exponential expression with neutrino
  //
  float part1 = 0; //TODO move down inside loop
  float part2 = 0;
  if (x <= 0)
    return 0;
  else {
	// TODO don't recompute 
    part1 = 1 / (1 + ex(-.511f * x / cmn.t9mev));
    part2 = 1 / (1 + ex(+(x - 2.531f) * (.511f / cmn.tnmev) - cmn.xi(1)));
    return cmn.cnorm * x * pow2((x - 2.531f)) * pow(
      (pow2(x) - 1), 0.5) * part1 * part2;
  }
}

//typedef float (*func2_function_pointer)(common&, float const&);

//
//20--------2ND PART OF INTEGRAL FOR n->p RATE-----------------------------------
//
float
common::func2(
  common& cmn,
  float const& x)
{
  float return_value = 0;
  // COMMON modpr
  //arr_cref<float> xi(cmn.xi, dimension(3));
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] rate1, nudens
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Contains integrands to be integrated.
  //
  //----------COMMON AREAS.
  //Model parameters.
  //Integration parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------MODEL PARAMETERS.
  //Neutrino degeneracy parameters.
  //
  //----------NEUTRINO PARAMETERS.
  //Temperature (in units of MeV).
  //Neutrino temperature (in units of Me
  //Neutrino temperature (units of 10**9
  //Normalizing constant.
  //
  //----------LOCAL VARIABLES.
  //Value at which function is evaluated
  //Exponential expression with photon t
  //Exponential expression with neutrino
  //
  float part1 = 0; // TODO move down into loop
  float part2 = 0;
  if (x <= 1.f) {
    return_value = 0.f;
  }
  else {
    part1 = 1.f / (1.f + ex(+.511f * x / cmn.t9mev));
    part2 = 1.f / (1.f + ex(-(x + 2.531f) * (.511f / cmn.tnmev) - cmn.xi(1)));
    return_value = cmn.cnorm * x * pow2((x + 2.531f)) * pow(
      (pow2(x) - 1), .5f) * part1 * part2;
  }
  return return_value;
}

//typedef float (*func3_function_pointer)(common&, float const&);

//
//30--------1ST PART OF INTEGRAL FOR p->n RATE-----------------------------------
//
float
common::func3(
  common& cmn,
  float const& x)
{
  float return_value = 0; // TODO move inside?
  // COMMON modpr
  //arr_cref<float> xi(cmn.xi, dimension(3));
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] rate1, nudens
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Contains integrands to be integrated.
  //
  //----------COMMON AREAS.
  //Model parameters.
  //Integration parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------MODEL PARAMETERS.
  //Neutrino degeneracy parameters.
  //
  //----------NEUTRINO PARAMETERS.
  //Temperature (in units of MeV).
  //Neutrino temperature (in units of Me
  //Neutrino temperature (units of 10**9
  //Normalizing constant.
  //
  //----------LOCAL VARIABLES.
  //Value at which function is evaluated
  //Exponential expression with photon t
  //Exponential expression with neutrino
  //
  float part1 = 0;
  float part2 = 0;
  if (x <= 1.f) {
    return_value = 0.f;
  }
  else {
    part1 = 1.f / (1.f + ex(-.511f * x / cmn.t9mev));
    part2 = 1.f / (1.f + ex(+(x + 2.531f) * (.511f / cmn.tnmev) + cmn.xi(1)));
    return_value = cmn.cnorm * x * pow2((x + 2.531f)) * pow(
      (pow2(x) - 1), .5f) * part1 * part2;
  }
  return return_value;
}

//typedef float (*func4_function_pointer)(common&, float const&);

//
//40--------2ND PART OF INTEGRAL FOR p->n RATE-----------------------------------
//
float
common::func4(
  common& cmn,
  float const& x)
{
  float return_value = 0;
  // COMMON modpr
  //arr_cref<float> xi(cmn.xi, dimension(3));
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] rate1, nudens
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Contains integrands to be integrated.
  //
  //----------COMMON AREAS.
  //Model parameters.
  //Integration parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------MODEL PARAMETERS.
  //Neutrino degeneracy parameters.
  //
  //----------NEUTRINO PARAMETERS.
  //Temperature (in units of MeV).
  //Neutrino temperature (in units of Me
  //Neutrino temperature (units of 10**9
  //Normalizing constant.
  //
  //----------LOCAL VARIABLES.
  //Value at which function is evaluated
  //Exponential expression with photon t
  //Exponential expression with neutrino
  //
  float part1 = 0;
  float part2 = 0;
  if (x <= 1.f) {
    return_value = 0.f;
  	return 0;
  }
  else {
    part1 = 1.f / (1.f + ex(+.511f * x / cmn.t9mev));
    part2 = 1.f / (1.f + ex(-(x - 2.531f) * (.511f / cmn.tnmev) + cmn.xi(1)));
    return_value = cmn.cnorm * x * pow2((x - 2.531f)) * pow(
      (pow2(x) - 1), .5f) * part1 * part2;
  }
  return return_value;
  // TODO test and uncomment
  /*
  if (x <= 1.f) {
    return_value = 0.f;
  	return 0;
  }
  else {
	const float me = 0.511;
	const float K = x - 2.531;
    const float part1 = 1 + exp(+me * x / cmn.t9mev);
    const float part2 = 1 + exp(-K*(.511f / cmn.tnmev) + xi(1));
    return cmn.cnorm * x * K * K * sqrt(x*x - 1) / (part1 * part2);
  }
  */
}

//typedef float (*func5_function_pointer)(common&, float const&);

//
//50--------INTEGRAL FOR ENERGY DENSITY OF NEUTRINO------------------------------
//
float
common::func5(
  common& cmn,
  float const& x)
{
  float return_value = 0;
  // COMMON modpr
  //arr_cref<float> xi(cmn.xi, dimension(3));
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] rate1, nudens
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Contains integrands to be integrated.
  //
  //----------COMMON AREAS.
  //Model parameters.
  //Integration parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------MODEL PARAMETERS.
  //Neutrino degeneracy parameters.
  //
  //----------NEUTRINO PARAMETERS.
  //Temperature (in units of MeV).
  //Neutrino temperature (in units of Me
  //Neutrino temperature (units of 10**9
  //Normalizing constant.
  //
  //----------LOCAL VARIABLES.
  //Value at which function is evaluated
  //Exponential expression with photon t
  //Exponential expression with neutrino
  //
  return_value = 1.f / (2 * pow2(3.14159f)) * pow3(x) / (
    1.f + exp(x / cmn.tnu - cmn.xi(cmn.nu)));
  return return_value;
}

//typedef float (*func6_function_pointer)(common&, float const&);

//
//60--------INTEGRAL FOR ENERGY DENSITY OF ANTINEUTRINO--------------------------
//
float
common::func6(
  common& cmn,
  float const& x)
{
  float return_value = 0;
  // COMMON modpr
  //arr_cref<float> xi(cmn.xi, dimension(3));
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] rate1, nudens
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Contains integrands to be integrated.
  //
  //----------COMMON AREAS.
  //Model parameters.
  //Integration parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------MODEL PARAMETERS.
  //Neutrino degeneracy parameters.
  //
  //----------NEUTRINO PARAMETERS.
  //Temperature (in units of MeV).
  //Neutrino temperature (in units of Me
  //Neutrino temperature (units of 10**9
  //Normalizing constant.
  //
  //----------LOCAL VARIABLES.
  //Value at which function is evaluated
  //Exponential expression with photon t
  //Exponential expression with neutrino
  //
  return_value = 1.f / (2 * pow2(3.14159f)) * pow3(x) / (
    1.f + exp(x / cmn.tnu + cmn.xi(cmn.nu)));
  return return_value;
}

/*
struct xintd_save
{
  int np;
  arr<float> u;
  arr<float> w;

  xintd_save() :
    np(0),
    u(dimension(6), fem::fill0),
    w(dimension(6), fem::fill0)
  {}
};
*/

//
//----------REFERENCES-----------------------------------------------------------
//     Forms of the integrals involved can be found in
//       Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
//       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
//
//========================IDENTIFICATION DIVISION================================
//
float
common::xintd(
  common& cmn,
  float const& xlow,
  float const& xhi,
  //func1_function_pointer func,
  float(* func)(common &, const float&),
  int const& nq)
{
  //FEM_CMN_SVE(xintd);
  // SAVE
  //
  const int np = 6;
  static const float u[] = {
        -.93246951420315f, -.66120938646627f, -.23861918608320f,
          .23861918608320f, .66120938646627f, .93246951420315f
      };
  static const float w[] = {
	.17132449237917f, .36076157304814f, .46791393457269f,
	  .46791393457269f, .36076157304814f, .17132449237917f
  };
  /*
  //int& np = sve.np;
  //arr_ref<float> u(sve.u, dimension(6));
  //arr_ref<float> w(sve.w, dimension(6));
  if (is_called_first_time) {
    {
      static const float values[] = {
        -.93246951420315f, -.66120938646627f, -.23861918608320f,
          .23861918608320f, .66120938646627f, .93246951420315f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        u;
    }
    {
      static const float values[] = {
        .17132449237917f, .36076157304814f, .46791393457269f,
          .46791393457269f, .36076157304814f, .17132449237917f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        w;
    }
  }
  */
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] rate1, nudens
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Computes the integral of the function "func".
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------INPUT VARIABLES.
  //Array of low limits.
  //Array of high limits.
  //Number of six point gaussian quads.
  //
  //----------COMPUTATION VARIABLES.
  //Size of quad interval.
  //Center of quad interval.
  //Variables of integration.
  //Summation of terms.
  //
  //----------COUNTERS.
  //Interval number.
  //Point number.
  //Total number of points in interval.
  //
  //----------ABSCISSAS AND WEIGHT FACTORS.
  //Abscissas.
  //Weight factor.
  //
  //==============================DATA DIVISION====================================
  //
  //----------ABSCISSAS AND WEIGHT FACTORS.
  //6 point Gaussian integration.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------DO INTEGRATION-------------------------------------------------------
  //
  float sum = 0;
  int nint = 0;
  int npnt = 0;
  float dist = (xhi - xlow) / fem::ffloat(nq); 				/// Size of quad interval.
  FEM_DO_SAFE(nint, 1, nq) {
    float cent = xlow + (fem::ffloat(nint) - 0.5) * dist; 	/// Center of interval.
    FEM_DO_SAFE(npnt, 1, np) {
      //x = cent + 0.5f * dist * u(npnt); 					/// Integration point.
      float x = cent + 0.5 * dist * u[npnt-1]; 				/// Integration point.
      float y = func(cmn, x); 								/// Evaluate function x(1).
      //sum += f * w(npnt); 								/// Add up sum.
      sum += y * w[npnt-1]; 								/// Add up sum.
    }
  }
  //
  //20--------GET INTEGRAL VALUE---------------------------------------------------
  //
  //Do integral.
  return 0.5 * sum * dist;
}

/*
struct rate1_save
{
  fem::variant_bindings rates_bindings;
};
*/

//
//========================IDENTIFICATION DIVISION================================
//
void
common::rate1(
  common& cmn,
  float const& tph)
{
/*
  FEM_CMN_SVE(rate1);
  // COMMON modpr
  float& tau = cmn.tau;
  arr_cref<float> xi(cmn.xi, dimension(3));
  // COMMON thermcb
  arr_cref<float> thm(cmn.thm, dimension(14));
  // COMMON nupar
  float& t9mev = cmn.t9mev;
  float& tnmev = cmn.tnmev;
  //
  //common_variant rates(cmn.common_rates, sve.rates_bindings);
  const int nrec = 88;
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r(dimension(nrec));
      rates.allocate(), f, r;
    }
  }
  arr_ref<float> f(rates.bind<float>(), dimension(nrec));
  arr_ref<float> r(rates.bind<float>(), dimension(nrec));
  */
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] start, derivs
  //     CALLS     - [function] xintd, eval
  //
  //----------REMARKS.
  //     Generates rate coefficients for weak n->p and p->n reactions.
  //
  //----------PARAMETERS.
  //Number of nuclear reactions.
  //Number of gaussian quads.
  //
  //----------COMMON AREAS.
  //Reaction rates.
  //Model parameters.
  //Dynamic variables.
  //Integration parame
  //
  //----------EXTERNAL FUNCTIONS.
  //Part 1 of n->p rate.
  //Part 2 of n->p rate.
  //Part 1 of p->n rate.
  //Part 2 of p->n rate.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION RATES.
  //Forward reaction rate coefficients.
  //Reverse reaction rate coefficients.
  //
  //----------EARLY UNIVERSE MODEL PARAMETERS.
  //Neutron lifetime.
  //Neutrino degeneracy parameters.
  //xi(1) is e neutrino degeneracy param
  //xi(2) is m neutrino degeneracy param
  //xi(3) is t neutrino degeneracy param
  //
  //----------DYNAMIC VARIABLES.
  //Thermodynamic variables (energy dens
  //
  //----------NEUTRINO PARAMETERS.
  //Temperature (in units of MeV).
  //Neutrino temperature (in units of Me
  //
  //----------LOCAL VARIABLES.
  //Photon temperature.
  //Upper limits for exponentials, forwa
  //Upper limits for exponentials, rever
  //Upper limits for integrals for forwa
  //Upper limits for integrals for rever
  //Parts of integrals for forward rate.
  //Parts of integrals for reverse rate.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------COMPUTE WEAK REACTION RATES (NONDEGENERATE)--------------------------
  //
  arr_1d<2, float> w(fem::fill0);
  arr_1d<2, float> x(fem::fill0);
  arr_1d<2, float> y(fem::fill0);
  arr_1d<2, float> z(fem::fill0);
  float uplim1 = 0;
  float uplim2 = 0;
  float uplim3 = 0;
  float uplim4 = 0;
  const int iter = 50;
  float part1 = 0;
  float part2 = 0;
  float part3 = 0;
  float part4 = 0;
  if (xi(1) == 0.f) {
    //f(1) = thm(13) / tau; 	/// Forward rate for weak np reaction.
    f[0] = thm(13) / tau; 	/// Forward rate for weak np reaction.
    //r(1) = thm(14) / tau; 	/// Reverse rate for weak np reaction.
    r[0] = thm(14) / tau; 	/// Reverse rate for weak np reaction.
  }
  else {
    //
    //20--------COMPUTE WEAK REACTION RATES (DEGENERATE)-----------------------------
    //
    //Convert photon temp to units of MeV.
    t9mev = tph * .086171f;
    //Convert neutrino temp to units of Me
    tnmev = cmn.tnu * .086171f;
    //..........COMPUTE OVERFLOW LIMITS FOR LIMITS OF INTEGRATION (Ref 1 & 2).
    w(1) = (-(t9mev / .511f) * (-88.722f));
    w(2) = ((tnmev / .511f) * (88.029f + xi(1)) + 2.531f);
    x(1) = ((t9mev / .511f) * (88.029f));
    x(2) = (-(tnmev / .511f) * (-88.722f + xi(1)) - 2.531f);
    y(1) = (-(t9mev / .511f) * (-88.722f));
    y(2) = ((tnmev / .511f) * (88.029f - xi(1)) - 2.531f);
    z(1) = ((t9mev / .511f) * (88.029f));
    z(2) = (-(tnmev / .511f) * (-88.722f - xi(1)) + 2.531f);
    //..........COMPARE LIMITS AND TAKE LARGER OF THE TWO.
    uplim1 = fem::abs(w(1));
    uplim2 = fem::abs(x(1));
    uplim3 = fem::abs(y(1));
    uplim4 = fem::abs(z(1));
    if (uplim1 < fem::abs(w(2))) {
      uplim1 = w(2);
    }
    if (uplim2 < fem::abs(x(2))) {
      uplim2 = x(2);
    }
    if (uplim3 < fem::abs(y(2))) {
      uplim3 = y(2);
    }
    if (uplim4 < fem::abs(z(2))) {
      uplim4 = z(2);
    }
    //..........EVALUATE THE INTEGRALS NUMERICALLY.
    part1 = xintd(cmn, 1., uplim1, func1, iter);
    part2 = xintd(cmn, 1., uplim2, func2, iter);
    part3 = xintd(cmn, 1., uplim3, func3, iter);
    part4 = xintd(cmn, 1., uplim4, func4, iter);
    //f(1) = part1 + part2; 		/// Add 2 integrals to get forward rate.
    f[0] = part1 + part2; 		/// Add 2 integrals to get forward rate.
    //r(1) = part3 + part4; 		/// Add 2 integrals to get reverse rate.
    r[0] = part3 + part4; 		/// Add 2 integrals to get reverse rate.
    //(xi(1).eq.0.)
  }
  //
  //----------REFERENCES-----------------------------------------------------------
  //     1) Forms of the integrals involved can be found in
  //          Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
  //          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
  //
  //     2) The overflow limit for the VAX/VMS system is exp(88.029).
  //        The underflow limit for the VAX/VMS system is exp(-88.722).
  //
}

struct start_save
{
  fem::variant_bindings rates_bindings;
};

//
//========================IDENTIFICATION DIVISION================================
//
void
common::start(common& cmn)
{
/*
  FEM_CMN_SVE(start);
  // COMMON evolp1
  float& t9 = cmn.t9;
  float& hv = cmn.hv;
  const int nnuc = 26;
  arr_ref<float> y(static_cast<common_evolp1&>(cmn).y, dimension(nnuc));
  // COMMON evolp3
  arr_ref<float> y0(cmn.y0, dimension(nnuc));
  // COMMON modpr
  float& tau = cmn.tau;
  arr_cref<float> c(cmn.c, dimension(3));
  arr_cref<float> xi(cmn.xi, dimension(3));
  // COMMON endens
  float& rhob0 = cmn.rhob0;
  // COMMON nupar
  float& tnu = cmn.tnu;
  float& cnorm = cmn.cnorm;
  //
  common_variant rates(cmn.common_rates, sve.rates_bindings);
  const int nrec = 88;
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r;
      rates.allocate(), f, r;
    }
  }
  */
  //arr_cref<float> f(rates.bind<float>(), dimension(nrec));
  ///* float const& r */ rates.bind<float>();
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] driver
  //     CALLS     - [subroutine] rate1, bessel, rate0
  //               - [function] ex
  //
  //----------REMARKS.
  //     Sets initial conditions.
  //
  //----------PARAMETERS.
  //Number of nuclear reactions.
  //Number of nuclides in calculation.
  //Relation between time and temperatur
  //Gravitational constant.
  //
  //----------COMMON AREAS.
  //Reaction rates.
  //Evolution paramete
  //Computation parame
  //Model parameters.
  //Variational parame
  //Time variables.
  //Energy densities.
  //Eval of function
  //Eval of function b
  //Flags,counters.
  //Neutrino parameter
  //Run options.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION RATES.
  //Forward reaction rate coefficients.
  //
  //----------EVOLUTION PARAMETERS.
  //Temperature (in units of 10**9 K).
  //Defined by hv = M(atomic)n(baryon)/t
  //Chemical potential of electron.
  //Relative number abundances.
  //
  //----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
  //Rel # abund at start of iteration.
  //
  //----------COMPUTATION SETTINGS.
  //Initial temperature (in 10**9 K).
  //Smallest abundances allowed.
  //Accumulation increment.
  //
  //----------EARLY UNIVERSE MODEL PARAMETERS.
  //Gravitational constant.
  //Neutron lifetime.
  //Number of neutrino species.
  //c(1) is variation of grav. constant.
  //c(2) is neutron lifetime (sec).
  //c(3) is number of neutrino species.
  //Neutrino degeneracy parameters.
  //
  //----------VARIATIONAL PARAMETERS.
  //Initial time step.
  //Baryon-to-photon ratio.
  //
  //----------TIME VARIABLES.
  //Time.
  //Time step.
  //
  //----------ENERGY DENSITIES.
  //Initial electron neutrino mass densi
  //Initial baryon mass density.
  //
  //----------EVALUATION OF FUNCTIONS bl,bm,bn.
  //Evaluation of function bl(z).
  //
  //----------COUNTERS AND FLAGS.
  //Indicates if output buffer printed.
  //# total time steps for particular ru
  //# time steps after outputting a line
  //# times accumulated in output buffer
  //Indicates if gaussian elimination fa
  //
  //----------NEUTRINO PARAMETERS.
  //Neutrino temperature.
  //Normalizing constant.
  //
  //----------RUN OPTION.
  //Number of nuclides in computation.
  //
  //----------LOCAL VARIABLES.
  //Defined by z = m(electron)*c**2/k*t9
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------INITIALIZE FLAGS AND COUNTERS----------------------------------------
  //
  //No output yet.
  cmn.ltime = 0;
  //First iteration coming up.
  cmn.is = 1;
  //Set to maximum allowed # of iteratio
  cmn.ip = cmn.inc;
  //No accumulation yet.
  cmn.it = 0;
  //No computational errors.
  cmn.mbad = 0;
  //
  //20--------SETTINGS-------------------------------------------------------------
  //
  //..........COMPUTATIONAL SETTINGS.
  //Initial temperature.
  t9 = cmn.t9i;
  //Initial neutrino temperature.
  tnu = t9;
  //Initial time (Ref 1).
  const float const1 = 0.09615f;
  cmn.t = 1 / pow2((const1 * t9));
  //Initial time step.
  cmn.dt = cmn.dt1;
  //..........MODEL SETTINGS.
  //Modify gravitational constant.
  const float const2 = 6.6700e-8f;
  cmn.g = const2 * c(1);
  //Convert n half-life (min) to lifetim
  tau = c(2);
  //Coulomb correction (Ref 2).
  tau = tau / 0.98f;
  //Number of neutrino species.
  cmn.xnu = c(3);
  //
  //30--------COMPUTE INITIAL ABUNDANCES FOR NEUTRON AND PROTON--------------------
  //
  //Overabundance of antineut
  if ((15.011f / t9 + xi(1)) > 58.f) {
    y(1) = 1.e-25f; 			/// Very little of neutrons.
    y(2) = 1.f; 				/// Essentially all protons.
  }
  else {
    //Overabundance of neutrino
    if ((15.011f / t9 + xi(1)) <  - 58.f) {
      y(1) = 1.f; 				/// Essentially all neutrons.
      y(2) = 1.e-25f; 			/// Very little of protons.
    }
    else {
      //Initial n abundance (Ref
      y(1) = 1 / (ex(15.011f / t9 + xi(1)) + 1);
      //Initial p abundance (Ref
      y(2) = 1 / (ex(-15.011f / t9 - xi(1)) + 1);
    }
  }
  //Electron neutrino degeneracy.
  if (xi(1) != 0.f) {
    cnorm = 1.;
    tnu = .00001f; 				/// Low temperature.
    rate1(cmn, 0.00001f); 		/// Find normalization constant at low temperature.
    //cnorm = 1 / tau / f(1);
    cnorm = 1 / tau / f[0];
  }
  y0(1) = y(1);
  y0(2) = y(2);
  //
  //40--------FIND RATIO OF BARYON DENSITY TO TEMPERATURE CUBED--------------------
  //
  float z = 5.930f / t9; 				//Inverse of temperature.
  //bessel(cmn, z); // old call
  float bl1 = getBesselL(z);
  float bl2 = getBesselL(2*z);
  float bl3 = getBesselL(3*z);
  float bl4 = getBesselL(4*z);
  float bl5 = getBesselL(5*z);
  /*
  float bm1 = getBesselM(z);
  float bm2 = getBesselM(2*z);
  float bm3 = getBesselM(3*z);
  float bm4 = getBesselM(4*z);
  float bm5 = getBesselM(5*z);
  float bn1 = getBesselN(z);
  float bn2 = getBesselN(2*z);
  float bn3 = getBesselN(3*z);
  float bn4 = getBesselN(4*z);
  float bn5 = getBesselN(5*z);
  */
  hv = 3.3683e+4f * cmn.eta1 * 2.75; 		/// (Ref 4 but with final eta).
  cmn.phie = hv * (1.784e-5f * y(2)) / 
  	(0.5*z*z*z*(bl1 - 2*bl2 + 3*bl3 - 4*bl4 + 5*bl5));
											/// Chemical potential of electron (Ref 5).
  rhob0 = hv * pow3(t9); 					/// TODO Baryon density. 
  //Nonde
  if ((xi(1) == 0) && (xi(2) == 0) && (xi(3) == 0)) {
    cmn.rhone0 = 7.366f * pow4(t9); 		/// Electron neutrino density (Ref 6).
  }
  //
  //50--------SET ABUNDANCES FOR REST OF NUCLIDES----------------------------------
  //
  y(3) = y(1) * y(2) * rhob0 * ex(25.82f / t9) / (.471e+10f * pow(t9, 1.5f)); /// (Ref 7).
  y0(3) = y(3);
  int i = 0;
  FEM_DO_SAFE(i, 4, cmn.isize) {
    y(i) = cmn.ytmin; 						/// Set rest to minimum abundance.
    y0(i) = y(i); 							/// Initialize abundances at beginning of iter
  }
  //Compute weak decay rates.
  rate0(cmn);
  //
  //----------REFERENCES-----------------------------------------------------------
  //     1) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
  //          page 44, equation A15.
  //     2) Coulomb correction obtained by dividing by correction factor Fp(t9)
  //               Fp(t9) = 1 - 0.5(pi/(137<v>/c))
  //          Wagoner, R.V. 1973, Ap. J. 179, page 358.
  //     3) For the nondegenerate case:
  //          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
  //          page 4, equation 3.
  //        For the case with neutrino degeneracy:
  //          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49,
  //          page 417, equation 9.
  //     4) Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 250, equation 4.
  //          3.3683e+4 = Mu(ng/t9**3) with Mu the atomic mass, ng the
  //          photon density.  2.75 is for the 11/4 factor difference
  //          between the initial and final values of eta.
  //     5) Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
  //          Kellogg Radiation Lab preprint OAP-714.
  //          equation D.2.
  //     6) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
  //          page 43, equation A4.
  //          7.366 is used instead of 14.73 as the latter is the sum total
  //          for 2 neutrino species.
  //     7) Initial deuterium abundance from nuclear statistical equilibrium
  //          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
  //          page 19, equation 17.
  //
}

//
//========================IDENTIFICATION DIVISION================================
//
void
common::nudens(
  common& cmn)
{
	/*
  // COMMON modpr
  arr_cref<float> xi(cmn.xi, dimension(3));
  // COMMON nupar
  float& tnu = cmn.tnu;
  float& rhonu = cmn.rhonu;
  int& nu = cmn.nu;
  */
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] therm
  //     CALLS     - [function] xintd, eval
  //
  //----------REMARKS.
  //     Computes energy density contribution from neutrinos.
  //
  //----------PARAMTER.
  //Number of gaussian quads.
  //
  //----------COMMON AREAS.
  //Model parameters.
  //Integration parame
  //
  //----------EXTERNAL FUNCTIONS.
  //Integral for neutrinos.
  //Integral for antineutrinos.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------MODEL PARAMETERS.
  //Neutrino degeneracy parameters.
  //
  //----------NEUTRINO PARAMETERS.
  //Neutrino temperature (units of 10**9
  //Neutrino energy density.
  //Which neutrino type.
  //
  //----------LOCAL VARIABLES.
  //Upper limit for neutrino energy inte
  //Upper limit for antineu energy integ
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------COMPUTE NEUTRINO ENERGY DENSITIES------------------------------------
  //
  float uplim1 = 0;
  float uplim2 = 0;
  const int iter = 50;
  if (fem::abs(xi(nu)) <= 0.03f) {
    //..........SMALL xi APPROXIMATION.
    rhonu = 2.f * (pow2(3.14159f) / 30.f) * pow4((tnu)) * (
      7.f / 8.f + (15.f / (4 * pow2(3.14159f))) * pow2(xi(
      nu)) + (15.f / (8.f * pow4(3.14159f))) * pow4(xi(
      nu)));
  }
  else {
    if (fem::abs(xi(nu)) >= 30.f) {
      //..........LARGE xi APPROXIMATION.
      rhonu = (pow4((tnu))) / (8.f * pow2(3.14159f)) *
        pow4(xi(nu)) * (1 + 12.f * 1.645f / pow2(xi(nu)));
    }
    else {
      //..........DO INTEGRATION
      uplim1 = (88.029f + xi(nu)) * tnu;
      uplim2 = (88.029f - xi(nu)) * tnu;
      if (uplim2 <= 0.) {
        rhonu = xintd(cmn, 0, uplim1, func5, iter);
      }
      else {
        rhonu = xintd(cmn, 0, uplim1, func5, iter) 
			  + xintd(cmn, 0, uplim2, func6, iter);
      }
      //(abs(xi(nu)).ge.30.)
    }
    //(abs(xi(nu)).le.0.03)
  }
  //
  //----------REFERENCES-----------------------------------------------------------
  //     Forms of the integrals involved can be found in
  //       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
  //       Freese, K., Kolb, E.W., Turner, M.S., 1983, Phys. Rev. D, 27, 1689.
  //
}

//
//========================IDENTIFICATION DIVISION================================
//
void
common::therm(
  common& cmn)
{
  /*
  // COMMON evolp1
  float& t9 = cmn.t9;
  float& phie = cmn.phie;
  // COMMON modpr
  float& xnu = cmn.xnu;
  arr_cref<float> xi(cmn.xi, dimension(3));
  // COMMON thermcb
  arr_ref<float> thm(cmn.thm, dimension(14));
  // COMMON endens
  float& rnb = cmn.rnb;
  // COMMON besselcb
  float& bl1 = cmn.bl1;
  float& bl2 = cmn.bl2;
  float& bl3 = cmn.bl3;
  float& bl4 = cmn.bl4;
  float& bl5 = cmn.bl5;
  float& bm1 = cmn.bm1;
  float& bm2 = cmn.bm2;
  float& bm3 = cmn.bm3;
  float& bm4 = cmn.bm4;
  float& bm5 = cmn.bm5;
  */
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] derivs
  //     CALLS     - [subroutine] bessel, nudens
  //               - [function] ex
  //
  //----------REMARKS.
  //     Computes various temperature dependent thermodynamic quantities.
  //
  //----------PARAMETER.
  //Number of nuclides in calculation.
  //(mass(neutron)-mass(proton))/m(elect
  //
  //----------COMMON AREAS.
  //Evolution paramete
  //Computation parame
  //Model parameters.
  //Dynamic variables.
  //Energy densities.
  //Eval of function b
  //Integration parame
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------EVOLUTION PARAMETERS.
  //Temperature (in units of 10**9 K).
  //Chemical potential for electron.
  //
  //----------COMPUTATION PARAMETERS.
  //Initial temperature (in 10**9 K).
  //
  //----------EARLY UNIVERSE MODEL PARAMETERS.
  //Number of neutrino species.
  //Neutrino degeneracy parameters.
  //
  //----------DYNAMIC VARIABLES.
  //Thermodynamic variables.
  //
  //----------ENERGY DENSITIES.
  //Initial electron neutrino mass densi
  //Initial baryon mass density.
  //Baryon mass density (ratio to init v
  //
  //----------EVALUATION OF FUNCTIONS bl,bm,bn.
  //Evaluation of function bl(z).
  //Evaluation of function bm(z).
  //Evaluation of function bn(z).
  //
  //----------NEUTRINO PARAMETERS.
  //Neutrino temperature.
  //Neutrino energy density.
  //Type of neutrino.
  //
  //----------LOCAL VARIABLE.
  //Defined by z = m(electron)*c**2/k*t9
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------COMPUTE FACTORS------------------------------------------------------
  //
  //z = m(electron)c**2/k(t9).
  float z = 5.930f / t9;
  float bl1 = getBesselL(z);
  float bl2 = getBesselL(2*z);
  float bl3 = getBesselL(3*z);
  float bl4 = getBesselL(4*z);
  float bl5 = getBesselL(5*z);
  float bm1 = getBesselM(z);
  float bm2 = getBesselM(2*z);
  float bm3 = getBesselM(3*z);
  float bm4 = getBesselM(4*z);
  float bm5 = getBesselM(5*z);
  float bn1 = getBesselN(z);
  float bn2 = getBesselN(2*z);
  float bn3 = getBesselN(3*z);
  float bn4 = getBesselN(4*z);
  float bn5 = getBesselN(5*z);
  //Neutrino temperature.
  cmn.tnu = (pow((rnb), (1.f / 3.f))) * cmn.t9i;
  //..........FACTORS OF z.
  float z1 = z;
  float z2 = z*z;
  float z3 = z*z*z;
  float z4 = z*z*z*z;
  float z5 = z*z*z*z*z;
  //..........TRIGONOMETRIC FUNCTION VALUES.
  //No chance of overflow.
  float cosh1 = 0;
  float cosh2 = 0;
  float cosh3 = 0;
  float cosh4 = 0;
  float cosh5 = 0;
  float sinh1 = 0;
  float sinh2 = 0;
  float sinh3 = 0;
  float sinh4 = 0;
  float sinh5 = 0;
  if (phie <= 17.f) {
    cosh1 = cosh(phie);
    cosh2 = cosh(2*phie);
    cosh3 = cosh(3*phie);
    cosh4 = cosh(4*phie);
    cosh5 = cosh(5*phie);
    sinh1 = sinh(phie);
    sinh2 = sinh(2*phie);
    sinh3 = sinh(3*phie);
    sinh4 = sinh(4*phie);
    sinh5 = sinh(5*phie);
  }
  else {
    cosh1 = 0;
    cosh2 = 0;
    cosh3 = 0;
    cosh4 = 0;
    cosh5 = 0;
    sinh1 = 0;
    sinh2 = 0;
    sinh3 = 0;
    sinh4 = 0;
    sinh5 = 0;
  }
  //bessel(cmn, z);
  //
  //20--------COMPUTE THERMODYNAMIC VARIABLES--------------------------------------
  //
  //(Ref 1)
  thm(1) = 8.418f * t9 * t9 * t9 * t9;
  //(Ref 2)
  thm(2) = 4.f * thm(1) / t9;
  //(Ref 3)
  thm(3) = thm(1) / 3.f;
  //(Ref 4)
  thm(4) = 3206.f * (bm1 * cosh1 - bm2 * cosh2 + bm3 * cosh3 - bm4 *
    cosh4 + bm5 * cosh5);
  //(Ref 5)
  thm(5) = 3206.f * (z / t9) * (bn1 * cosh1 - 2.f * bn2 *
    cosh2 + 3.f * bn3 * cosh3 - 4.f * bn4 * cosh4 + 5.f *
    bn5 * cosh5);
  //(Ref 6)
  thm(6) = 3206.f * (bm1 * sinh1 - 2.f * bm2 * sinh2 + 3.f * bm3 *
    sinh3 - 4.f * bm4 * sinh4 + 5.f * bm5 * sinh5);
  //(Ref 7)
  thm(7) = 3206.f * (bl1 * cosh1 / z - bl2 * cosh2 / (2.f * z) +
    bl3 * cosh3 / (3.f * z) - bl4 * cosh4 / (4.f * z) + bl5 * cosh5 /
    (5.f * z));
  //Nonde
  if ((xi(1) == 0.f) && (xi(2) == 0.f) && (xi(3) == 0)) {
    //(Ref 8)
    thm(8) = xnu * cmn.rhone0 * (pow(rnb, (4.f / 3.f)));
    //Include effects of neutrino degenera
  }
  else {
    thm(8) = 0.f;
    //For every neutrino family.
    FEM_DO_SAFE(cmn.nu, 1, xnu) {
      //Compute neutrino energy density.
      nudens(cmn);
      //Have 12.79264 from units ch
      thm(8) += 12.79264f * cmn.rhonu;
    }
  }
  //(Ref 9)
  thm(9) = cmn.rhob0 * rnb;
  //(Ref 10
  thm(10) = thm(1) + thm(4) + thm(8) + thm(9);
  //(Ref 11
  thm(11) = -(pow3(z) / t9) * (sinh1 * (3.f * bl1 - z * bm1) - // TODO change to z3
    sinh2 * (3.f * bl2 - 2.f * z * bm2) + sinh3 * (3.f * bl3 - 3.f *
    z * bm3) - sinh4 * (3.f * bl4 - 4.f * z * bm4) + sinh5 * (3.f *
    bl5 - 5.f * z * bm5));
  //(Ref 12
  thm(12) = pow3(z) * (cosh1 * bl1 - 2.f * cosh2 * bl2 + 3.f * // TODO change to z3
    cosh3 * bl3 - 4.f * cosh4 * bl4 + 5.f * cosh5 * bl5);
  if (thm(12) != 0.f) {
    thm(12) = 1.f / thm(12);
  }
  //(Ref 13
  thm(13) = 1.000f + 0.565f / z1 - 6.382f / z2 + 11.108f / z3 +
    36.492f / z4 + 27.512f / z5;
  //(Ref 14
  const float q = 2.531f;
  thm(14) = (5.252f / z1 - 16.229f / z2 + 18.059f / z3 + 34.181f /
    z4 + 27.617f / z5) * ex(-q * z);
  //
  //----------REFERENCES AND NOTES-------------------------------------------------
  //     1)  thm(1)  = rho photon
  //         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
  //          page 43, equation A2.)
  //     2)  thm(2)  = d(rho photon)/d(t9)
  //     3)  thm(3)  = (p photon)/c**2
  //         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967,
  //          page 43, equation A3.)
  //     4)  thm(4)  = rho electron+positron
  //         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9,
  //          page 281, equation B44.)
  //     5)  thm(5)  = d(rho electron+positron)/d(t9)
  //     6)  thm(6)  = d(rho electron+positron)/d(phi e)
  //     7)  thm(7)  = (p electron+positron)/c**2
  //         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9,
  //          page 279, equation B27.)
  //     8)  thm(8)  = rho neutrino
  //                 = # neutrino species x rho electron neutrino (nondegenerate)
  //                 = rho nu(e) + rho nu(m) + rho nu(t)          (degenerate)
  //     9)  thm(9)  = rho baryon
  //     10) thm(10) = rho total
  //                 = rho photon + rho electron+positron + rho neutrino
  //                              + rho baryon
  //     11) thm(11) = d     /pi**2(hbar*c)**3(ne- - ne+)*z**3\
  //                   d(t9) \  2  (mc**2)**3                 /
  //     12) thm(12) = d        /pi**2(hbar*c)**3(ne- - ne+)*z**3\
  //                   d(phi e) \  2  (mc**2)**3                 /
  //     13) thm(13) = rate for n->p
  //     14) thm(14) = rate for p->n
  //
}

/*
struct eqslin_save
{
  fem::variant_bindings lncoef_bindings;
};
*/

//
//========================IDENTIFICATION DIVISION================================
//
void
common::eqslin(
  common& cmn,
  int const& icnvm,
  int& ierror)
{
  //FEM_CMN_SVE(eqslin);
  /*
  int& inc = cmn.inc;
  int& mbad = cmn.mbad;
  int& isize = cmn.isize;
  //
  const int nnuc = 26;
  common_variant lncoef(cmn.common_lncoef, sve.lncoef_bindings);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> a(dimension(nnuc, nnuc));
      mbr<float> b(dimension(nnuc));
      mbr<float> y(dimension(nnuc));
      lncoef.allocate(), a, b, y;
    }
  }
  arr_ref<double, 2> a(lncoef.bind<double>(), dimension(nnuc, nnuc));
  arr_cref<float> b(lncoef.bind<float>(), dimension(nnuc));
  arr_ref<float> y(lncoef.bind<float>(), dimension(nnuc));
  */
  int nord = 0;
  int i = 0;
  arr_1d<nnuc, float> x(fem::fill0);
  int j = 0;
  arr<float, 2> a0(dimension(nnuc, nnuc), fem::fill0);
  float cx = 0;
  int k = 0;
  float sum = 0;
  float xdy = 0;
  const float eps = 2.e-4f;
  const int mord = 1;
  float r = 0;
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] sol
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Solves for new abundances using gaussian elimination
  //     with back substitution, no pivoting.
  //
  //----------PARAMETERS.
  //Rank of matrix.
  //Higher order in correction.
  //Tolerance for convergence (.ge. 1.e-
  //
  //----------COMMON AREAS.
  //Computation parame
  //Lin eqn coefficien
  //Flags, counters.
  //Run options.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------COMPUTATION PARAMETER.
  //Accumulation increment.
  //
  //----------MATRIX COEFFICIENTS FOR LINEAR EQUATION.
  //Coefficient array.
  //Right-hand vector w/o manipulation.
  //Solution vector.
  //
  //----------COUNTERS AND FLAGS.
  //Indicates type of error.
  //
  //----------RUN OPTION.
  //Number of nuclides in computation.
  //
  //----------LOCAL MATRICES AND VECTORS.
  //Coefficient array w/o manipulation.
  //Right-hand vector.
  //
  //----------LOCAL COMPUTATION VARIABLES.
  //Scaling factor in triangularization.
  //Sum for backsubstitution.
  //Relative error.
  //
  //----------LOCAL COUNTERS.
  //Order of correction.
  //Convergence monitor.
  //ith nuclide fails to converge.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------INITIALIZE VECTOR----------------------------------------------------
  //
  //..........SET COUNTERS TO ZERO.
  //No corrections yet.
  nord = 0;
  //No errors yet.
  mbad = 0;
  //..........SET RIGHT-HAND AND SOLUTION VECTORS TO INITIAL VALUES.
  FEM_DO_SAFE(i, 1, isize) {
    x(i) = b(i); 		/// Right-hand vector.
    y(i) = 0.f; 		/// Solution vector.
  }
  //..........SAVE MATRIX.
  //Monitor convergence.
  if (icnvm == inc) {
    FEM_DO_SAFE(i, 1, isize) {
      FEM_DO_SAFE(j, 1, isize) {
        //Initial value of coefficient array.
        a0(j, i) = a(j, i);
      }
    }
  }
  //
  //20--------TRIANGULARIZE MATRIX AND SAVE OPERATOR-------------------------------
  //
  //..........CHECK TO SEE THAT THERE ARE NO ZEROES AT PIVOT POINTS.
  FEM_DO_SAFE(i, 1, isize - 1) {
    //Don't want to divide by zero.
    if (a(i, i) == 0) {
      //Position of zero coefficient.
      mbad = i;
      //Terminate matrix evaluation.
      return;
    }
    //..........TRIANGULARIZE MATRIX.
    FEM_DO_SAFE(j, i + 1, isize) {
      //Progress diagonally down the column.
      if (a(j, i) != 0) {
        //Scaling factor down the column.
        cx = a(j, i) / a(i, i);
        //Progress diagonally along row.
        FEM_DO_SAFE(k, i + 1, isize) {
          //Subtract scaled coeff along
          a(j, k) = a(j, k) - cx * a(i, k);
        }
        //Scaled coefficient.
        a(j, i) = cx;
        //..........OPERATE ON RIGHT-HAND VECTOR.
        //Subtract off scaled coefficient.
        x(j) = x(j) - cx * x(i);
      }
    }
  }
  //
  //30--------DO BACK SUBSTITUTION-------------------------------------------------
  //
  statement_300:
  //Solution for ultimate positi
  x(isize) = x(isize) / a(isize, isize);
  y(isize) += x(isize);
  //From i = penultimate to i = 1.
  FEM_DOSTEP(i, isize - 1, 1, -1) {
    sum = 0.e0;
    FEM_DO_SAFE(j, i + 1, isize) {
      //Sum up all previous terms.
      sum += a(i, j) * x(j);
    }
    x(i) = (x(i) - sum) / a(i, i);
    //Add difference to initial value.
    y(i) += x(i);
  }
  //
  //40--------TESTS AND EXITS------------------------------------------------------
  //
  if (icnvm == inc) {
    FEM_DO_SAFE(i, 1, isize) {
      if (y(i) != 0.f) {
        //Relative error.
        xdy = fem::dabs(x(i) / y(i));
        if (xdy > eps) {
          //Continue to higher orders.
          if (nord < mord) {
            nord++;
            //..........FIND ERROR IN RIGHT-HAND VECTOR.
            FEM_DO_SAFE(j, 1, isize) {
              //Initialize r.
              r = 0.e0;
              FEM_DO_SAFE(k, 1, isize) {
                //Left side with approximate sol
                r += a0(j, k) * y(k);
              }
              //Subtract difference from right side.
              x(j) = b(j) - r;
            }
            //..........OPERATE ON RIGHT-HAND VECTOR.
            FEM_DO_SAFE(j, 1, isize - 1) {
              FEM_DO_SAFE(k, j + 1, isize) {
                //Subtract off scaled coef
                x(k) = x(k) - a(k, j) * x(j);
              }
            }
            //Go for another iteratiion.
            goto statement_300;
          }
          else {
            //..........NOT ENOUGH CONVERGENCE.
            //Signal error problem.
            mbad = -1;
            //ith nuclide for which x/y checked.
            ierror = i;
            return;
            //(nord.lt.mord)
          }
          //(xdy.gt.eps)
        }
        //(y(i).ne.0)
      }
      //i = 1,isize
    }
    //(icnvm.eq.inc)
  }
  //No more iterations & relative error
  //
}

/*
struct sol_save
{
  fem::variant_bindings lncoef_bindings;
  fem::variant_bindings rates_bindings;
  arr<float> si;
  arr<float> sj;
  arr<float> sk;
  arr<float> sl;

  sol_save() :
    si(dimension(11), fem::fill0),
    sj(dimension(11), fem::fill0),
    sk(dimension(11), fem::fill0),
    sl(dimension(11), fem::fill0)
  {}
};
*/

//
//========================IDENTIFICATION DIVISION================================
//
void
common::sol(
  common& cmn,
  int const& loop)
{
  //FEM_CMN_SVE(sol);
  common_write write(cmn);
  /*
  const int nrec = 88;
  arr_cref<int> iform(cmn.iform, dimension(nrec));
  arr_cref<int> ii(cmn.ii, dimension(nrec));
  arr_cref<int> jj(cmn.jj, dimension(nrec));
  arr_cref<int> kk(cmn.kk, dimension(nrec));
  arr_cref<int> ll(cmn.ll, dimension(nrec));
  arr_cref<float> rev(cmn.rev, dimension(nrec));
  arr_cref<float> q9(cmn.q9, dimension(nrec));
  float& t9 = cmn.t9;
  const int nnuc = 26;
  arr_cref<float> y(static_cast<common_evolp1&>(cmn).y, dimension(nnuc));
  arr_ref<float> dydt(cmn.dydt, dimension(nnuc));
  arr_cref<float> y0(cmn.y0, dimension(nnuc));
  float& dt = cmn.dt;
  float& rhob = cmn.rhob;
  int& mbad = cmn.mbad;
  int& isize = cmn.isize;
  //
  common_variant rates(cmn.common_rates, sve.rates_bindings);
  common_variant lncoef(cmn.common_lncoef, sve.lncoef_bindings);
  arr_ref<float> si(sve.si, dimension(11));
  arr_ref<float> sj(sve.sj, dimension(11));
  arr_ref<float> sk(sve.sk, dimension(11));
  arr_ref<float> sl(sve.sl, dimension(11));
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r(dimension(nrec));
      rates.allocate(), f, r;
    }
    {
      mbr<double> a(dimension(nnuc, nnuc));
      mbr<float> b(dimension(nnuc));
      mbr<float> yx(dimension(nnuc));
      lncoef.allocate(), a, b, yx;
    }
  }
  arr_ref<float> f(rates.bind<float>(), dimension(nrec));
  arr_ref<float> r(rates.bind<float>(), dimension(nrec));
  arr_ref<double, 2> a(lncoef.bind<double>(), dimension(nnuc, nnuc));
  arr_ref<float> b(lncoef.bind<float>(), dimension(nnuc));
  arr_cref<float> yx(lncoef.bind<float>(), dimension(nnuc));
  if (is_called_first_time) {
    {
      static const float values[] = {
        1.f, 1.f, 1.f, 1.f, 1.f, 2.f, 3.f, 2.f, 1.f, 1.f, 2.f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        si;
    }
    {
      static const float values[] = {
        0.f, 1.f, 1.f, 0.f, 1.f, 0.f, 0.f, 1.f, 1.f, 1.f, 0.f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        sj;
    }
    {
      static const float values[] = {
        0.f, 0.f, 1.f, 0.f, 0.f, 1.f, 0.f, 0.f, 1.f, 0.f, 2.f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        sk;
    }
    {
      static const float values[] = {
        1.f, 1.f, 1.f, 2.f, 2.f, 1.f, 1.f, 1.f, 2.f, 3.f, 1.f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        sl;
    }
  }
  */
	const float si[] = {1, 1, 1, 1, 1, 2, 3, 2, 1, 1, 2};
	const float sj[] = {0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0};
	const float sk[] = {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2};
	const float sl[] = {1, 1, 1, 2, 2, 1, 1, 1, 2, 3, 1};
  float t932 = 0;
  float t9m32 = 0;
  int isize1 = 0;
  int i = 0;
  int j = 0;
  int n = 0;
  int ind = 0;
  int k = 0;
  int l = 0;
  float ri = 0;
  float rj = 0;
  float rk = 0;
  float rl = 0;
  float ci = 0;
  float cj = 0;
  float ck = 0;
  float cl = 0;
  float bdln = 0;
  int i1 = 0;
  int j1 = 0;
  int ierror = 0;
  arr_1d<nnuc, float> yy(fem::fill0);
  const int iw = 6;
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] derivs
  //     CALLS     - [subroutine] eqslin
  //               - [function] ex
  //
  //----------REMARKS.
  //     Computes reverse strong and electromagnetic reaction rates.
  //     Fills and solves matrix equation for dydt(i).
  //
  //----------PARAMETERS.
  //Input unit number.
  //Output unit number.
  //Number of nuclear reactions.
  //Number of nuclides in calculation.
  //
  //-----------COMMON AREAS.
  //Reaction parameter
  //Reaction rates.
  //Evolution paramete
  //Time varying param
  //Dynamic variable
  //Energy densities.
  //Linear eqn coeffic
  //Flags,counters.
  //Run option.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION PARAMETERS.
  //Reaction code number (1-88).
  //Incoming nuclide type (1-26).
  //Incoming light nuclide type (1-6).
  //Outgoing light nuclide type (1-6).
  //Outgoing nuclide type (1-26).
  //Reverse reaction coefficient.
  //Energy released in reaction.
  //
  //----------REACTION RATES.
  //Forward reaction rate coefficients.
  //Reverse reaction rate coefficients.
  //
  //----------EVOLUTION PARAMETERS.
  //Temperature (in units of 10**9 K).
  //Relative number abundances.
  //
  //----------EVOLUTION PARAMETERS (DERIVATIVES).
  //Change in rel number abundances.
  //
  //----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
  //Rel # abund at start of iteration.
  //
  //----------TIME VARIABLES.
  //Time step.
  //
  //----------DYNAMIC VARIABLES.
  //Expansion rate.
  //
  //----------ENERGY DENSITIES.
  //Baryon mass density.
  //
  //----------COMPONENTS OF MATRIX EQUATION.
  //Relates y(t-dt) to y(t).
  //Contains y0 in inverse order.
  //yy in reverse order.
  //
  //----------COUNTERS AND FLAGS.
  //Counts which Runge-Kutta loop.
  //# time steps after outputting a line
  //Indicates if gaussian elimination fa
  //
  //----------RUN OPTIONS.
  //Number of nuclides in computation.
  //Equals isize + 1.
  //Number of reactions in computation.
  //
  //----------EVOLUTION EQUATION COEFFICIENTS.
  //Equate to ii,jj,kk,ll.
  //Equate to si,sj,sk,sl.
  //Coefficients of rate equation.
  //
  //----------LOCAL VARIABLES.
  //Abundances at end of iteration.
  //# of nuclide i,j,k,l
  //(10**(-5))*volume expansion rate.
  //Equate to iform.
  //Element which does not converge.
  //
  //==============================DATA DIVISION====================================
  //
  //----------NUMBER OF NUCLIDES IN REACTION TYPES 1-11.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------TEMPERATURE FACTORS AND INITIAL VALUES-------------------------------
  //
  //..........TEMPERATURE FACTORS.
  //t9**(3/2).
  t932 = pow(t9, 1.5f);
  //t9**(-3/2).
  t9m32 = 1.f / t932;
  //..........MATRIX SIZE.
  isize1 = isize + 1;
  //..........INITIALIZE A-MATRIX.
  FEM_DO_SAFE(i, 1, isize) {
    FEM_DO_SAFE(j, 1, isize) {
      //Set a-matrix to zero.
      a(j, i) = 0.e0;
    }
  }
  //
  //20--------COMPUTE FACTORS FOR THE A-MATRIX-------------------------------------
  //
  FEM_DO_SAFE(n, 1, cmn.jsize) {
    //..........EQUATE VARIABLES TO ARRAYS.
    //Type of reaction.
    ind = iform(n);
    //ID # of incoming nuclide i.
    i = ii(n);
    //ID # of incoming nuclide j.
    j = jj(n);
    //ID # of outgoing nuclide k.
    k = kk(n);
    //ID # of outgoing nuclide l.
    l = ll(n);
    //Reactio
    if ((ind != 0) && (i <= isize) && (l <= isize)) {
      //# of incoming nuclide i.
      //ri = si(ind);
      ri = si[ind-1];
      //# of incoming nuclide j.
      //rj = sj(ind);
      rj = sj[ind-1];
      //# of outgoing nuclide k.
      //rk = sk(ind);
      rk = sk[ind-1];
      //# of outgoing nuclide l.
      //rl = sl(ind);
      rl = sl[ind-1];
      //..........COMPUTE DIFFERENT REACTION RATES.
      switch (ind) {
        case 1: goto statement_201;
        case 2: goto statement_202;
        case 3: goto statement_203;
        case 4: goto statement_204;
        case 5: goto statement_205;
        case 6: goto statement_206;
        case 7: goto statement_207;
        case 8: goto statement_208;
        case 9: goto statement_209;
        case 10: goto statement_210;
        case 11: goto statement_211;
        default: break;
      }
      statement_201: 			/// 1-0-0-1 configuration.
      //ci = f(n); 				/// (Ref 1).
      ci = f[n-1]; 				/// (Ref 1).
      cj = 0.f;
      ck = 0.f;
      //cl = r(n);
      cl = r[n-1];
      goto statement_212;
      statement_202: 			/// 1-1-0-1 configuration.
      //r(n) = rev(n) * 1.e+10f * t932 * ex(-q9(n) / t9) * f(n); 	/// (Ref 2).
      r[n-1] = rev(n) * 1.e+10f * t932 * ex(-q9(n) / t9) * f[n-1]; 	/// (Ref 2).
      //f(n) = rhob * f(n);
      f[n-1] = rhob * f[n-1];
      //ci = y(j) * f(n) / 2.f;
      ci = y(j) * [n-1]) / 2.f;
      //cj = y(i) * f(n) / 2.f;
      cj = y(i) * f[n-1] / 2.f;
      ck = 0;
      //cl = r(n);
      cl = r[n-1];
      goto statement_212;
      statement_203: /// 1-1-1-1 configuration.
      f[n-1] = rhob * f[n-1]; /// (Ref 3).
      r[n-1] = rev(n) * ex(-q9(n) / t9) * f[n-1];
      ci = y(j) * f[n-1] / 2.f;
      cj = y(i) * f[n-1] / 2.f;
      ck = y(l) * r[n-1] / 2.f;
      cl = y(k) * r[n-1] / 2.f;
      goto statement_212;
      //1-0-0-2 configuration.
      statement_204:
      ci = f[n-1];
      cj = 0.f;
      ck = 0.f;
      cl = y(l) * r[n-1] / 2.f;
      goto statement_212;
      //1-1-0-2 configuration.
      statement_205:
      f[n-1] = rhob * f[n-1];
      //(Ref 3).
      r[n-1] = rev(n) * ex(-q9(n) / t9) * f[n-1];
      ci = y(j) * f[n-1] / 2.f;
      cj = y(i) * f[n-1] / 2.f;
      ck = 0.f;
      cl = y(l) * r[n-1] / 2.f;
      goto statement_212;
      //2-0-1-1 configuration.
      statement_206:
      f[n-1] = rhob * f[n-1];
      //(Ref 3).
      r[n-1] = rev(n) * ex(-q9(n) / t9) * f[n-1];
      ci = y(i) * f[n-1] / 2.f;
      cj = 0.f;
      ck = y(l) * r[n-1] / 2.f;
      cl = y(k) * r[n-1] / 2.f;
      goto statement_212;
      //3-0-0-1 configuration.
      statement_207:
      //(Ref 4).
      r[n-1] = rev(n) * 1.e+20f * t932 * t932 * ex(-q9(n) / t9) * f[n-1];
      f[n-1] = rhob * rhob * f[n-1];
      ci = y(i) * y(i) * f[n-1] / 6.f;
      cj = 0.f;
      ck = 0.f;
      cl = r[n-1];
      goto statement_212;
      //2-1-0-1 configuration.
      statement_208:
      //(Ref 4).
      r[n-1] = rev(n) * 1.e+20f * t932 * t932 * ex(-q9(n) / t9) * f[n-1];
      f[n-1] = rhob * rhob * f[n-1];
      ci = y(j) * y(i) * f[n-1] / 3.f;
      cj = y(i) * y(i) * f[n-1] / 6.f;
      ck = 0.f;
      cl = r[n-1];
      goto statement_212;
      //1-1-1-2 configuration.
      statement_209:
      f[n-1] = rhob * f[n-1];
      //(Ref 5)
      r[n-1] = rev(n) * 1.e-10f * t9m32 * rhob * ex(-q9(n) / t9) * f[n-1];
      ci = y(j) * f[n-1] / 2.f;
      cj = y(i) * f[n-1] / 2.f;
      ck = y(l) * y(l) * r[n-1] / 6.f;
      cl = y(k) * y(l) * r[n-1] / 3.f;
      goto statement_212;
      //1-1-0-3 configuration.
      statement_210:
      f[n-1] = rhob * f[n-1];
      //(Ref 5)
      r[n-1] = rev(n) * 1.e-10f * t9m32 * rhob * ex(-q9(n) / t9) * f[n-1];
      ci = y(j) * f[n-1] / 2.f;
      cj = y(i) * f[n-1] / 2.f;
      ck = 0.f;
      cl = y(l) * y(l) * r[n-1] / 6.f;
      goto statement_212;
      //2-0-2-1 configuration.
      statement_211:
      f[n-1] = rhob * f[n-1];
      //(Ref 5)
      r[n-1] = rev(n) * 1.e-10f * t9m32 * rhob * ex(-q9(n) / t9) * f[n-1];
      ci = y(i) * f[n-1] / 2.f;
      cj = 0.f;
      ck = y(l) * y(k) * r[n-1] / 3.f;
      cl = y(k) * y(k) * r[n-1] / 6.f;
      statement_212:
      //
      //30--------CONSTRUCT THE A-MATRIX-----------------------------------------------
      //
      //Invert i index.
      i = isize1 - i;
      //Invert j index.
      j = isize1 - j;
      //Invert k index.
      k = isize1 - k;
      //Invert l index.
      l = isize1 - l;
      //..........FILL I NUCLIDE COLUMN.
      if (j <= isize) {
        a(j, i) += rj * ci;
      }
      if (k <= isize) {
        a(k, i) = a(k, i) - rk * ci;
      }
      a(i, i) += ri * ci;
      a(l, i) = a(l, i) - rl * ci;
      //..........FILL J NUCLIDE COLUMN.
      if (j <= isize) {
        a(j, j) += rj * cj;
        if (k <= isize) {
          a(k, j) = a(k, j) - rk * cj;
        }
        a(i, j) += ri * cj;
        a(l, j) = a(l, j) - rl * cj;
      }
      //..........FILL K NUCLIDE COLUMN.
      if (k <= isize) {
        if (j <= isize) {
          a(j, k) = a(j, k) - rj * ck;
        }
        a(k, k) += rk * ck;
        a(i, k) = a(i, k) - ri * ck;
        a(l, k) += rl * ck;
      }
      //..........FILL L NUCLIDE COLUMN.
      if (j <= isize) {
        a(j, l) = a(j, l) - rj * cl;
      }
      if (k <= isize) {
        a(k, l) += rk * cl;
      }
      a(i, l) = a(i, l) - ri * cl;
      a(l, l) += rl * cl;
      //((ind.ne.0).and.(i.le.isize).and.(l.le.isize))
    }
    //n = 1,jsize
  }
  //
  //40--------PUT A-MATRIX AND B-VECTOR IN FINAL FORM OF MATRIX EQUATION-----------
  //
  //(10**(-5))*(Expansion rate).
  bdln = 1.e-5f * (3.f * cmn.hubcst);
  FEM_DO_SAFE(i, 1, isize) {
    //Invert the rows.
    i1 = isize1 - i;
    FEM_DO_SAFE(j, 1, isize) {
      //Invert the columns.
      j1 = isize1 - j;
      if (fem::dabs(a(j, i)) < bdln * y0(j1) / y0(i1)) {
        //Set 0 if tiny.
        a(j, i) = 0.e0;
      }
      else {
        //Bring dt over to other side.
        a(j, i) = a(j, i) * dt;
      }
    }
    //Add identity matrix to a-matrix.
    a(i, i) += 1.e0;
    //Initial abundances.
    b(i1) = y0(i);
  }
  //
  //50--------SOLVE EQUATIONS TO GET DERIVATIVE------------------------------------
  //
  //..........SET MONITOR FLAG AND SOLVE BY GAUSSIAN ELIMINATION.
  if (loop == 1) {
    eqslin(cmn, cmn.ip, ierror);
  }
  else {
    eqslin(cmn, 0, ierror);
  }
  //..........OBTAIN DERIVATIVE.
  FEM_DO_SAFE(i, 1, isize) {
    yy(i) = yx(isize1 - i); 				/// Abundance at t+dt.
    dydt(i) = (yy(i) - y0(i)) / dt; 		/// Take derivative.
  }
  //
  //60--------POSSIBLE ERROR MESSAGES AND EXIT-------------------------------------
  //
  //Problem in gaussian elimination.
  if (mbad != 0) {
    //Error message.
    if (mbad ==  - 1) {
      write(iw, "(' ','** y(',i2,') fails to converge **')"), ierror;
    }
    //Error message.
    if (mbad >= 1) {
      write(iw, "(' ','** ',i2,' th diagonal term equals zero **')"), mbad;
    }
  }
  //
  //----------REFERENCES-----------------------------------------------------------
  //     1) The coefficients are given in general as:
  //             ci = ri*(y(j)**rj)*(y(i)**(ri-1)*f(n)/
  //                  ((ri+rj)*fac(ri)*fac(rj))
  //             cj = rj*(y(i)**ri)*(y(j)**(rj-1)*f(n)/
  //                  ((ri+rj)*fac(ri)*fac(rj))
  //             ck = rk*(y(l)**rl)*(y(k)**(rk-1)*f(n)/
  //                  ((rk+rl)*fac(rk)*fac(rl))
  //             cl = rl*(y(k)**rk)*(y(l)**(rl-1)*f(n)/
  //                  ((rk+rl)*fac(rk)*fac(rl))
  //        in which fac(x) is the factorial of x.
  //     2) Form of reverse rate given in
  //        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
  //          tables 1B, 4B, 7B.
  //     3) Form of reverse rate given in
  //        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
  //          tables 2B, 3B, 5B, 6B, 8B, 9B, 10B.
  //     4) Form of reverse rate given in
  //        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
  //          table 11B.
  //     5) Form of reverse rate given in
  //        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
  //          tables 12B, 13B.
  //
}

/*
struct rate2_save
{
  fem::variant_bindings rates_bindings;
};
*/

//
//========================IDENTIFICATION DIVISION================================
//
void
common::rate2(
  common& cmn)
{
  //FEM_CMN_SVE(rate2);
  // COMMON evolp1
  /*
  float& t9 = cmn.t9;
  //
  common_variant rates(cmn.common_rates, sve.rates_bindings);
  const int nrec = 88;
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r(dimension(nrec));
      rates.allocate(), f, r;
    }
  }
  arr_ref<float> f(rates.bind<float>(), dimension(nrec));
  */
  //* arr_cref<float> r( */ rates.bind<float>() /* , dimension(nrec)) */ ;
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] derivs
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Generates rate coefficients for reactions involving nuclides
  //     up to A = 9.
  //
  //----------PARAMETER.
  //Number of nuclear reactions.
  //Number of nuclides in calculation.
  //
  //----------COMMON AREAS.
  //Reaction rates.
  //Evolution parameters.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION RATES.
  //Forward reaction rate coefficients.
  //
  //----------EVOLUTION PARAMETER.
  //Temperature of photons (units of GK)
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------TEMPERATURE FACTORS--------------------------------------------------
  //
  //t9**(1/3)
  float t913 = pow(t9, (.33333333f));
  //t9**(2/3)
  float t923 = t913 * t913;
  //t9**(4/3)
  float t943 = t923 * t923;
  //t9**(5/3)
  float t953 = t9 * t923;
  //t9**(1/2)
  float t912 = sqrt(t9);
  //t9**(3/2)
  float t932 = t9 * t912;
  //t9**(-1)
  float t9m1 = 1 / t9;
  //t9**(-2/3)
  float t9m23 = 1.0f / t923;
  //t9**(-3/2)
  float t9m32 = 1.0f / t932;
  //For reaction 17.
  float t9a = t9 / (1.0f + 13.076f * t9);
  //t9a**(3/2)
  float t9a32 = pow(t9a, (1.5f));
  //For reaction 18.
  float t9b = t9 / (1.f + 49.18f * t9);
  //t9b**(3/2)
  float t9b32 = pow(t9b, (1.5f));
  //For reaction 22.
  float t9c = 0;
  if (t9 > 10.f) {
    t9c = 1.f;
  }
  else {
    t9c = t9 / (1.f - 9.69e-2f * t9 + 2.84e-2f * t953 / pow((
      1.f - 9.69e-2f * t9), (2.f / 3.f)));
  }
  //t9c**(1/3)
  float t9c13 = pow(t9c, (.3333333f));
  //t9c**(5/6)
  float t9c56 = pow(t9c, (.8333333f));
  //For reaction 24.
  float t9d = t9 / (1.f + 0.759f * t9);
  //t9d**(1/3)
  float t9d13 = pow(t9d, (.3333333f));
  //t9d**(5/6)
  float t9d56 = pow(t9d, (.8333333f));
  //For reaction 26.
  float t9e = t9 / (1.f + 0.1378f * t9);
  //t9e**(1/3)
  float t9e13 = pow(t9e, (.3333333f));
  //t9e**(5/6)
  float t9e56 = pow(t9e, (.8333333f));
  //For reaction 27.
  float t9f = t9 / (1.f + 0.1071f * t9);
  //t9f**(1/3)
  float t9f13 = pow(t9f, (.3333333f));
  //t9f**(5/6)
  float t9f56 = pow(t9f, (.8333333f));
  //
  //20--------NEUTRON, PHOTON REACTIONS--------------------------------------------
  //
  //.......H(n,g)H2...................(Smith-Kawano-Malaney 1992)
  f(12) = 4.742e+4f * (1.f - .8504f * t912 + .4895f * t9 - .09623f *
    t932 + 8.471e-3f * t9 * t9 - 2.80e-4f * t9 * t932);
  //
  //.......H2(n,g)H3..................(Wagoner 1969)
  f(13) = 6.62e+1f * (1.f + 18.9f * t9);
  //
  //.......He3(n,g)He4................(Wagoner 1969)
  f(14) = 6.62e+0f * (1.f + 905.f * t9);
  //
  //.......Li6(n,g)Li7................(Malaney-Fowler 1989)
  f(15) = 5.10e+3f;
  //
  //30--------NEUTRON, PROTON REACTIONS--------------------------------------------
  //
  //.......He3(n,p)H3.................(Smith-Kawano-Malaney 1992)
  f(16) = 7.21e+8f * (1.f - .508f * t912 + .228f * t9);
  //
  //.......Be7(n,p)Li7................(Smith-Kawano-Malaney 1992)
  f(17) = 2.675e+9f * (1.f - .560f * t912 + .179f * t9 - .0283f *
    t932 + 2.214e-3f * t9 * t9 - 6.851e-5f * t9 * t932) + 9.391e+8f *
    t9a32 * t9m32 + 4.467e+7f * t9m32 * ex(-0.07486f / t9);
  //
  //40--------NEUTRON, ALPHA REACTIONS---------------------------------------------
  //
  //.......Li6(n,a)H3.................(Caughlan-Fowler 1988)
  f(18) = 2.54e+9f * t9m32 * ex(-2.39f / t9) + 1.68e+8f * (1.f -
    .261f * t9b32 / t932);
  //
  //.......Be7(n,a)He4................(Wagoner 1969)
  f(19) = 2.05e+4f * (1.f + 3760.f * t9);
  //
  //50--------PROTON, PHOTON REACTIONS---------------------------------------------
  //
  //.......H2(p,g)He3.................(Smith-Kawano-Malaney 1992)
  f(20) = 2.65e+3f * t9m23 * ex(-3.720f / t913) * (1.f + .112f *
    t913 + 1.99f * t923 + 1.56f * t9 + .162f * t943 + .324f * t953);
  //
  //.......H3(p,g)He4.................(Caughlan-Fowler 1988)
  f(21) = 2.20e+4f * t9m23 * ex(-3.869f / t913) * (1.f + .108f *
    t913 + 1.68f * t923 + 1.26f * t9 + .551f * t943 + 1.06f * t953);
  //
  //.......Li6(p,g)Be7................(Caughlan-Fowler 1988)
  f(22) = 6.69e+5f * t9c56 * t9m32 * ex(-8.413f / t9c13);
  //
  //60--------PROTON, ALPHA REACTIONS----------------------------------------------
  //
  //.......Li6(p,a)He3................(Caughlan-Fowler 1988)
  f(23) = 3.73e+10f * t9m23 * ex(-8.413f / t913 - pow2((t9 /
    5.50f))) * (1.f + .050f * t913 - .061f * t923 - .021f * t9 +
    .006f * t943 + .005f * t953) + 1.33e+10f * t9m32 * ex(-17.763f /
    t9) + 1.29e+09f * t9m1 * ex(-21.820f / t9);
  //
  //.......Li7(p,a)He4................(Smith-Kawano-Malaney 1992)
  f(24) = 1.096e+9f * t9m23 * ex(-8.472f / t913) - 4.830e+8f *
    t9d56 * t9m32 * ex(-8.472f / t9d13) + 1.06e+10f * t9m32 * ex(
    -30.442f / t9) + 1.56e+5f * t9m23 * ex((-8.472f / t913) -
    pow2((t9 / 1.696f))) * (1.f + .049f * t913 - 2.498f * t923 +
    .860f * t9 + 3.518f * t943 + 3.08f * t953) + 1.55e+6f * t9m32 *
    ex(-4.478f / t9);
  //
  //70--------ALPHA, PHOTON REACTIONS----------------------------------------------
  //
  //.......H2(a,g)Li6.................(Caughlan-Fowler 1988)
  f(25) = 3.01e+01f * t9m23 * ex(-7.423f / t913) * (1.f + .056f *
    t913 - 4.85f * t923 + 8.85f * t9 - .585f * t943 - .584f * t953) +
    8.55e+1f * t9m32 * ex(-8.228f / t9);
  //
  //.......H3(a,g)Li7.................(Smith-Kawano-Malaney 1992)
  f(26) = 3.032e+5f * t9m23 * ex(-8.090f / t913) * (1.f + .0516f *
    t913 + .0229f * t923 + 8.28e-3f * t9 - 3.28e-4f * t943 -
    3.01e-4f * t953) + 5.109e+5f * t9e56 * t9m32 * ex(-8.068f /
    t9e13);
  //
  //.......He3(a,g)Be7................(Smith-Kawano-Malaney 1992)
  f(27) = 4.817e+6f * t9m23 * ex(-14.964f / t913) * (1.f + .0325f *
    t913 - 1.04e-3f * t923 - 2.37e-4f * t9 - 8.11e-5f * t943 -
    4.69e-5f * t953) + 5.938e+6f * t9f56 * t9m32 * ex(-12.859f /
    t9f13);
  //
  //80--------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------------
  //
  //.......H2(d,n)He3.................(Smith-Kawano-Malaney 1992)
  f(28) = 3.95e+8f * t9m23 * ex(-4.259f / t913) * (1.f + .098f *
    t913 + .765f * t923 + .525f * t9 + 9.61e-3f * t943 + .0167f *
    t953);
  //
  //.......H2(d,p)H3..................(Smith-Kawano-Malaney 1992)
  f(29) = 4.17e+8f * t9m23 * ex(-4.258f / t913) * (1.f + .098f *
    t913 + .518f * t923 + .355f * t9 - .010f * t943 - .018f * t953);
  //
  //.......H3(d,n)He4.................(Smith-Kawano-Malaney 1992)
  f(30) = 1.063e+11f * t9m23 * ex(-4.559f / t913 - pow2((t9 /
    .0754f))) * (1.f + .092f * t913 - .375f * t923 - .242f * t9 +
    33.82f * t943 + 55.42f * t953) + 8.047e+8f * t9m23 * ex(
    -0.4857f / t9);
  //
  //.......He3(d,p)He4................(Smith-Kawano-Malaney 1992)
  f(31) = 5.021e+10f * t9m23 * ex(-7.144f / t913 - pow2((t9 /
    .270f))) * (1.f + .058f * t913 + .603f * t923 + .245f * t9 +
    6.97f * t943 + 7.19f * t953) + 5.212e+8f / t912 * ex(-1.762f /
    t9);
  //
  //90--------THREE PARTICLE REACTIONS---------------------------------------------
  //
  //.......He3(He3,2p)He4.............(Caughlan-Fowler 1988)
  f(32) = 6.04e+10f * t9m23 * ex(-12.276f / t913) * (1.f + .034f *
    t913 - .522f * t923 - .124f * t9 + .353f * t943 + .213f * t953);
  //
  //.......Li7(d,na)He4...............(Caughlan-Fowler 1988)
  f(33) = 2.92e+11f * t9m23 * ex(-10.259f / t913);
  //
  //.......Be7(d,pa)He4...............(Caughlan-Fowler 1988)
  f(34) = 1.07e+12f * t9m23 * ex(-12.428f / t913);
  //
  //----------REFERENCES-----------------------------------------------------------
  //     Smith, M., Kawano, L.H., and Malaney, R.A., 1992, submitted to Ap. J.
  //     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5.
  //     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data
  //       Tables, 40, 283.
  //     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.
  //
}

/*
struct rate3_save
{
  fem::variant_bindings rates_bindings;
};
*/

//
//========================IDENTIFICATION DIVISION================================
//
void
common::rate3(
  common& cmn)
{
	/*
  FEM_CMN_SVE(rate3);
  // COMMON evolp1
  float& t9 = cmn.t9;
  //
  common_variant rates(cmn.common_rates, sve.rates_bindings);
  const int nrec = 88;
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r(dimension(nrec));
      rates.allocate(), f, r;
    }
  }
  arr_ref<float> f(rates.bind<float>(), dimension(nrec));
  */
  //* arr_cref<float> r( */ rates.bind<float>() /* , dimension(nrec)) */ ;
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] derivs
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Generates rate coefficients for reactions involving nuclides
  //     up to A = 18.
  //
  //----------PARAMETER.
  //Number of nuclear reactions.
  //Number of nuclides in calculation.
  //
  //----------COMMON AREAS.
  //Reaction rates.
  //Evolution parameters.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION RATES.
  //Forward reaction rate coefficients.
  //
  //----------EVOLUTION PARAMETER.
  //Temperature of photons (units of 10*
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------TEMPERATURE FACTORS--------------------------------------------------
  //
  //t9**(1/3)
  float t913 = pow(t9, (.33333333f));
  //t9**(2/3)
  float t923 = t913 * t913;
  //t9**(4/3)
  float t943 = t923 * t923;
  //t9**(5/3)
  float t953 = t9 * t923;
  //t9**(1/2)
  float t912 = sqrt(t9);
  //t9**(3/2)
  float t932 = t9 * t912;
  //t9**(1/5)
  float t915 = pow(t9, (.2f));
  //t9**(5/4)
  float t954 = pow(t9, (1.25f));
  //t9**(-1)
  float t9m1 = 1.0f / t9;
  //t9**(-2/3)
  float t9m23 = 1.0f / t923;
  //t9**(-3/2)
  float t9m32 = 1.0f / t932;
  //t9**(-3/4)
  float t9m34 = sqrt(t9m32);
  //t9**(-1/5)
  float t9m15 = 1.0f / t915;
  //t9**(-5/4)
  float t9m54 = 1.0f / t954;
  //For reaction 53.
  float t9a = t9 / (1.f + t9 / 15.1f);
  //t9a**(1/3)
  float t9a13 = pow(t9a, (.3333333f));
  //t9a**(5/6)
  float t9a56 = pow(t9a, (.8333333f));
  //
  //20--------NEUTRON, PHOTON REACTIONS--------------------------------------------
  //
  //.......Li7(n,g)Li8................(Wagoner 1969)
  f(35) = 4.90e+3f + 9.96e+3f * t9m32 * ex(-2.62f / t9);
  //
  //.......B10(n,g)B11................(Wagoner 1969)
  f(36) = 6.62e+4f;
  //
  //.......B11(n,g)B12................(Malaney-Fowler 1989)
  f(37) = 7.29e+2f + 2.40e+3f * t9m32 * ex(-0.223f / t9);
  //
  //30--------NEUTRON, PROTON REACTIONS--------------------------------------------
  //
  //.......C11(n,p)B11................(Caughlan-Fowler 1988)
  f(38) = 1.69e+8f * (1.f - .048f * t912 + .010f * t9);
  //
  //40--------NEUTRON, ALPHA REACTIONS---------------------------------------------
  //
  //.......B10(n,a)Li7................(Caughlan-Fowler 1988)
  f(39) = 5.07e+8f;
  //
  //50--------PROTON, PHOTON REACTIONS---------------------------------------------
  //
  //.......Be7(p,g)B8.................(Caughlan-Fowler 1988)
  f(40) = 3.11e+5f * t9m23 * ex(-10.262f / t913) + 2.53e+3f * t9m32 *
    ex(-7.306f / t9);
  //
  //.......Be9(p,g)B10................(Caughlan-Fowler 1988)
  f(41) = 1.33e+7f * t9m23 * ex(-10.359f / t913 - pow2((t9 /
    .846f))) * (1.f + .040f * t913 + 1.52f * t923 + .428f * t9 +
    2.15f * t943 + 1.54f * t953) + 9.64e+4f * t9m32 * ex(-3.445f /
    t9) + 2.72e+6f * t9m32 * ex(-10.620f / t9);
  //
  //.......B10(p,g)C11................(Caughlan-Fowler 1988)
  f(42) = 4.61e+5f * t9m23 * ex(-12.062f / t913 - pow2((t9 /
    4.402f))) * (1.f + .035f * t913 + .426f * t923 + .103f * t9 +
    .281f * t943 + .173f * t953) + 1.93e+5f * t9m32 * ex(-12.041f /
    t9) + 1.14e+4f * t9m32 * ex(-16.164f / t9);
  //
  //.......B11(p,g)C12................(Caughlan-Fowler 1988)
  f(43) = 4.62e+7f * t9m23 * ex(-12.095f / t913 - pow2((t9 /
    .239f))) * (1.f + .035f * t913 + 3.00f * t923 + .723f * t9 +
    9.91f * t943 + 6.07f * t953) + 7.89e+3f * t9m32 * ex(-1.733f /
    t9) + 9.68e+4f * t9m15 * ex(-5.617f / t9);
  //
  //.......C11(p,g)N12................(Caughlan-Fowler 1988)
  f(44) = 4.24e+4f * t9m23 * ex(-13.658f / t913 - pow2((t9 /
    1.627f))) * (1.f + .031f * t913 + 3.11f * t923 + .665f * t9 +
    4.61f * t943 + 2.50f * t953) + 8.84e+3f * t9m32 * ex(-7.021f /
    t9);
  //
  //60--------PROTON, NEUTRON REACTIONS--------------------------------------------
  //
  //.......B12(p,n)C12................(Wagoner 1969)
  f(45) = 4.02e+11f * t9m23 * ex(-12.12f / t913);
  //
  //70--------PROTON, ALPHA REACTIONS----------------------------------------------
  //
  //.......Be9(p,a)Li6................(Caughlan-Fowler 1988)
  f(46) = 2.11e+11f * t9m23 * ex(-10.359f / t913 - pow2((t9 /
    .520f))) * (1.f + .040f * t913 + 1.09f * t923 + .307f * t9 +
    3.21f * t943 + 2.30f * t953) + 4.51e+8f * t9m1 * ex(-3.046f /
    t9) + 6.70e+8f * t9m34 * ex(-5.160f / t9);
  //
  //.......B10(p,a)Be7................(Caughlan-Fowler 1988)
  f(47) = 1.26e+11f * t9m23 * ex(-12.062f / t913 - pow2((t9 /
    4.402f))) * (1.f + .035f * t913 - .498f * t923 - .121f * t9 +
    .300f * t943 + .184f * t953) + 2.59e+9f * t9m1 * ex(-12.260f /
    t9);
  //
  //.......B12(p,a)Be9................(Wagoner 1969)
  f(48) = 2.01e+11f * t9m23 * ex(-12.12f / t913);
  //
  //80--------ALPHA, PHOTON REACTIONS----------------------------------------------
  //
  //.......Li6(a,g)B10................(Caughlan-Fowler 1988)
  f(49) = 4.06e+6f * t9m23 * ex(-18.790f / t913 - pow2((t9 /
    1.326f))) * (1.f + .022f * t913 + 1.54f * t923 + .239f * t9 +
    2.20f * t943 + .869f * t953) + 1.91e+3f * t9m32 * ex(-3.484f /
    t9) + 1.01e+4f * t9m1 * ex(-7.269f / t9);
  //
  //.......Li7(a,g)B11................(Caughlan-Fowler 1988)
  f(50) = 3.55e+7f * t9m23 * ex(-19.161f / t913 - pow2((t9 /
    4.195f))) * (1.f + .022f * t913 + .775f * t923 + .118f * t9 +
    .884f * t943 + .342f * t953) + 3.33e+2f * t9m32 * ex(-2.977f /
    t9) + 4.10e+4f * t9m1 * ex(-6.227f / t9);
  //
  //.......Be7(a,g)C11................(Caughlan-Fowler 1988)
  f(51) = 8.45e+7f * t9m23 * ex(-23.212f / t913 - pow2((t9 /
    4.769f))) * (1.f + .018f * t913 + .488f * t923 + .061f * t9 +
    .296f * t943 + .095f * t953) + 1.25e+4f * t9m32 * ex(-6.510f /
    t9) + 1.29e+5f * t9m54 * ex(-10.039f / t9);
  //
  //90--------ALPHA, PROTON REACTIONS----------------------------------------------
  //
  //.......B8(a,p)C11.................(Wagoner 1969)
  f(52) = 1.08e+15f * t9m23 * ex(-27.36f / t913);
  //
  //100-------ALPHA, NEUTRON REACTIONS---------------------------------------------
  //
  //.......Li8(a,n)B11................(Malaney-Fowler 1989)
  f(53) = 8.62e+13f * t9a56 * t9m32 * ex(-19.461f / t9a13);
  //
  //.......Be9(a,n)C12................(Caughlan-Fowler 1988)
  f(54) = 4.62e+13f * t9m23 * ex(-23.870f / t913 - pow2((t9 /
    .049f))) * (1.f + .017f * t913 + 8.57f * t923 + 1.05f * t9 + 74.51f *
    t943 + 23.15f * t953) + 7.34e-5f * t9m32 * ex(-1.184f / t9) +
    2.27e-1f * t9m32 * ex(-1.834f / t9) + 1.26e+5f * t9m32 * ex(-4.179f /
    t9) + 2.40e+8f * ex(-12.732f / t9);
  //
  //110-------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------------
  //
  //.......Be9(d,n)B10................(original Wagoner code)
  f(55) = 7.16e+8f * t9m23 * ex(6.44f - 12.6f / t913);
  //
  //.......B10(d,p)B11................(original Wagoner code)
  f(56) = 9.53e+8f * t9m23 * ex(7.30f - 14.8f / t913);
  //
  //.......B11(d,n)C12................(original Wagoner code)
  f(57) = 1.41e+9f * t9m23 * ex(7.40f - 14.8f / t913);
  //
  //120-------THREE PARTICLE REACTIONS---------------------------------------------
  //
  //.......He4(an,g)Be9...............(Caughlan-Fowler 1988)
  f(58) = (2.59e-6f / ((1.f + .344f * t9) * pow2(t9))) * ex(-1.062f / t9);
  //
  //.......He4(2a,g)C12...............(Caughlan-Fowler 1988)
  f(59) = 2.79e-8f * t9m32 * t9m32 * ex(-4.4027f / t9) + 1.35e-8f *
    t9m32 * ex(-24.811f / t9);
  //
  //.......Li8(p,na)He4...............(original Wagoner code)
  f(60) = 8.65e+9f * t9m23 * ex(-8.52f / t913 - pow2((t9 /
    2.53f))) + 2.31e+9f * t9m32 * ex(-4.64f / t9);
  //
  //.......B8(n,pa)He4................(original Wagoner code)
  f(61) = 4.02e+8f;
  //
  //.......Be9(p,da)He4...............(Caughlan-Fowler 1988)
  f(62) = 2.11e+11f * t9m23 * ex(-10.359f / t913 - pow2((t9 /
    .520f))) * (1.f + .040f * t913 + 1.09f * t923 + .307f * t9 +
    3.21f * t943 + 2.30f * t953) + 5.79e+8f * t9m1 * ex(-3.046f /
    t9) + 8.50e+8f * t9m34 * ex(-5.800f / t9);
  //
  //.......B11(p,2a)He4...............(Caughlan-Fowler 1988)
  f(63) = 2.20e+12f * t9m23 * ex(-12.095f / t913 - pow2((t9 /
    1.644f))) * (1.f + .034f * t913 + .140f * t923 + .034f * t9 +
    .190f * t943 + .116f * t953) + 4.03e+6f * t9m32 * ex(-1.734f /
    t9) + 6.73e+9f * t9m32 * ex(-6.262f / t9) + 3.88e+9f * t9m1 * ex(
    -14.154f / t9);
  //
  //.......C11(n,2a)He4...............(Wagoner 1969)
  f(64) = 1.58e+8f;
  //
  //----------REFERENCES-----------------------------------------------------------
  //     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5.
  //     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data
  //       Tables, 40, 283.
  //     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.
  //
}

/*
struct rate4_save
{
  fem::variant_bindings rates_bindings;
};
*/

//
//========================IDENTIFICATION DIVISION================================
//
void
common::rate4(
  common& cmn)
{
	/*
  FEM_CMN_SVE(rate4);
  // COMMON evolp1
  float& t9 = cmn.t9;
  //
  common_variant rates(cmn.common_rates, sve.rates_bindings);
  const int nrec = 88;
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r(dimension(nrec));
      rates.allocate(), f, r;
    }
  }
  arr_ref<float> f(rates.bind<float>(), dimension(nrec));
  */
  //* arr_cref<float> r( */ rates.bind<float>() /* , dimension(nrec)) */ ;
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] derivs
  //     CALLS     - [function] ex
  //
  //----------REMARKS.
  //     Generates rate coefficients for rest of reactions.
  //
  //----------PARAMETER.
  //Number of nuclear reactions.
  //Number of nuclides in calculation.
  //
  //----------COMMON AREAS.
  //Reaction rates.
  //Evolution parameters.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION RATES.
  //Forward reaction rate coefficients.
  //
  //----------EVOLUTION PARAMETER.
  //Temperature of photons (units of 10*
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------TEMPERATURE FACTORS--------------------------------------------------
  //
  //t9**(1/3)
  float t913 = pow(t9, (.33333333f));
  //t9**(2/3)
  float t923 = t913 * t913;
  //t9**(4/3)
  float t943 = t923 * t923;
  //t9**(5/3)
  float t953 = t9 * t923;
  //t9**(1/2)
  float t912 = sqrt(t9);
  //t9**(3/2)
  float t932 = t9 * t912;
  //t9**(3/5)
  float t935 = pow(t9, (.6f));
  //t9**(6/5)
  float t965 = pow(t9, (1.2f));
  //t9**(3/8)
  float t938 = pow(t9, (.375f));
  //t9**(1/3)
  float t9m13 = 1.0f / t913;
  //t9**(-2/3)
  float t9m23 = 1.0f / t923;
  //t9**(-3/2)
  float t9m32 = 1.0f / t932;
  //t9**(-6/5)
  float t9m65 = 1.0f / t965;
  //For reaction 82.
  float t9a = t9 / (1.f + 4.78e-2f * t9 + 7.56e-3f * t953 / pow((
    1.f + 4.78e-2f * t9), (2.f / 3.f)));
  //t9a**(1/3)
  float t9a13 = pow(t9a, (.33333333f));
  //t9a**(5/6)
  float t9a56 = pow(t9a, (.83333333f));
  //For reaction 84.
  float t9b = t9 / (1.f + 7.76e-2f * t9 + 2.64e-2f * t953 / pow((
    1.f + 7.76e-2f * t9), (2.f / 3.f)));
  //t9b**(1/3)
  float t9b13 = pow(t9b, (.33333333f));
  //t9b**(5/6)
  float t9b56 = pow(t9b, (.83333333f));
  //
  //20--------NEUTRON, PHOTON REACTIONS--------------------------------------------
  //
  //.......C12(n,g)C13................(Wagoner 1969)
  f(65) = 4.50e+2f;
  //
  //.......C13(n,g)C14................(Wagoner 1969)
  f(66) = 1.19e+2f + 2.38e+5f * t9m32 * ex(-1.67f / t9);
  //
  //.......N14(n,g)N15................(Wagoner 1969)
  f(67) = 9.94e+3f;
  //
  //30--------NEUTRON, PROTON REACTIONS--------------------------------------------
  //
  //.......N13(n,p)C13................(Caughlan-Fowler 1988)
  f(68) = 1.88e+8f * (1.f - .167f * t912 + .037f * t9);
  //
  //.......N14(n,p)C14................(Caughlan-Fowler 1988)
  f(69) = 2.39e+5f * (1.f + .361f * t912 + .502f * t9) + 1.112e+8f /
    t912 * ex(-4.983f / t9);
  //
  //.......O15(n,p)N15................(Caughlan-Fowler 1988)
  f(70) = 3.50e+8f * (1.f + .452f * t912 - .191f * t9);
  //
  //40--------NEUTRON, ALPHA REACTIONS---------------------------------------------
  //
  //.......O15(n,a)C12................(Caughlan-Fowler 1988)
  f(71) = 3.50e+7f * (1.f + .188f * t912 + .015f * t9);
  //
  //50--------PROTON, PHOTON REACTIONS---------------------------------------------
  //
  //.......C12(p,g)N13................(Caughlan-Fowler 1988)
  f(72) = 2.04e+7f * t9m23 * ex(-13.690f / t913 - pow2((t9 /
    1.500f))) * (1.f + .030f * t913 + 1.19f * t923 + .254f * t9 +
    2.06f * t943 + 1.12f * t953) + 1.08e+5f * t9m32 * ex(-4.925f /
    t9) + 2.15e+5f * t9m32 * ex(-18.179f / t9);
  //
  //.......C13(p,g)N14................(Caughlan-Fowler 1988)
  f(73) = 8.01e+7f * t9m23 * ex(-13.717f / t913 - pow2((t9 /
    2.000f))) * (1.f + .030f * t913 + .958f * t923 + .204f * t9 +
    1.39f * t943 + .753f * t953) + 1.21e+6f * t9m65 * ex(-5.701f /
    t9);
  //
  //.......C14(p,g)N15................(Caughlan-Fowler 1988)
  f(74) = 6.80e+6f * t9m23 * ex(-13.741f / t913 - pow2((t9 /
    5.721f))) * (1.f + .030f * t913 + .503f * t923 + .107f * t9 +
    .213f * t943 + .115f * t953) + 5.36e+3f * t9m32 * ex(-3.811f /
    t9) + 9.82e+4f * t9m13 * ex(-4.739f / t9);
  //
  //.......N13(p,g)O14................(Caughlan-Fowler 1988)
  f(75) = 4.04e+7f * t9m23 * ex(-15.202f / t913 - pow2((t9 /
    1.191f))) * (1.f + .027f * t913 - .803f * t923 - .154f * t9 +
    5.00f * t943 + 2.44f * t953) + 2.43e+5f * t9m32 * ex(-6.348f /
    t9);
  //
  //.......N14(p,g)O15................(Caughlan-Fowler 1988)
  f(76) = 4.90e+7f * t9m23 * ex(-15.228f / t913 - pow2((t9 /
    3.294f))) * (1.f + .027f * t913 - .778f * t923 - .149f * t9 + .261f *
    t943 + .127f * t953) + 2.37e+3f * t9m32 * ex(-3.011f / t9) +
    2.19e+4f * ex(-12.530f / t9);
  //
  //.......N15(p,g)O16................(Caughlan-Fowler 1988)
  f(77) = 9.78e+8f * t9m23 * ex(-15.251f / t913 - pow2((t9 /
    .450f))) * (1.f + .027f * t913 + .219f * t923 + .042f * t9 +
    6.83f * t943 + 3.32f * t953) + 1.11e+4f * t9m32 * ex(-3.328f /
    t9) + 1.49e+4f * t9m32 * ex(-4.665f / t9) + 3.80e+6f * t9m32 * ex(
    -11.048f / t9);
  //
  //60--------PROTON, ALPHA REACTIONS----------------------------------------------
  //
  //.......N15(p,a)C12................(Caughlan-Fowler 1988)
  f(78) = 1.08e+12f * t9m23 * ex(-15.251f / t913 - pow2((t9 /
    .522f))) * (1.f + .027f * t913 + 2.62f * t923 + .501f * t9 +
    5.36f * t943 + 2.60f * t953) + 1.19e+8f * t9m32 * ex(-3.676f /
    t9) + 5.41e+8f / t912 * ex(-8.926f / t9) + 4.72e+7f * t9m32 * ex(
    -7.721f / t9) + 2.20e+8f * t9m32 * ex(-11.418f / t9);
  //
  //70--------ALPHA, PHOTON REACTIONS----------------------------------------------
  //
  //.......C12(a,g)O16................(Caughlan-Fowler 1988)
  f(79) = 1.04e+8f / pow2(t9) * ex(-32.120f / t913 - pow2((t9 /
    3.496f))) / pow2((1.f + .0489f * t9m23)) + 1.76e+8f / pow2((
    t9)) / pow2((1.f + .2654f * t9m23)) * ex(-32.120f / t913) +
    1.25e+3f * t9m32 * ex(-27.499f / t9) + 1.43e-2f * pow((t9), 5) *
    ex(-15.541f / t9);
  //
  //80--------ALPHA, PROTON REACTIONS----------------------------------------------
  //
  //.......B10(a,p)C13................(Wagoner 1969)
  f(80) = 9.60e+14f * t9m23 * ex(-27.99f / t913);
  //
  //.......B11(a,p)C14................(Caughlan-Fowler 1988)
  f(81) = 5.37e+11f * t9m23 * ex(-28.234f / t913 - pow2((t9 /
    0.347f))) * (1.f + .015f * t913 + 5.575f * t923 + .576f * t9 +
    15.888f * t943 + 4.174f * t953) + 5.44e-3f * t9m32 * ex(-2.827f /
    t9) + 3.36e+2f * t9m32 * ex(-5.178f / t9) + 5.32e+6f / t938 * ex(
    -11.617f / t9);
  //
  //.......C11(a,p)N14................(Caughlan-Fowler 1988)
  f(82) = 7.15e+15f * t9a56 * t9m32 * ex(-31.883f / t9a13);
  //
  //.......N12(a,p)O15................(Caughlan-Fowler 1988)
  f(83) = 5.59e+16f * t9m23 * ex(-35.60f / t913);
  //
  //.......N13(a,p)O16................(Caughlan-Fowler 1988)
  f(84) = 3.23e+17f * t9b56 * t9m32 * ex(-35.829f / t9b13);
  //
  //90--------ALPHA, NEUTRON REACTIONS---------------------------------------------
  //
  //.......B10(a,n)N13................(Caughlan-Fowler 1988)
  f(85) = 1.20e+13f * t9m23 * ex(-27.989f / t913 - pow2((t9 / 9.589f)));
  //
  //.......B11(a,n)N14................(Caughlan-Fowler 1988)
  f(86) = 6.97e+12f * t9m23 * ex(-28.234f / t913 - pow2((t9 /
    0.140f))) * (1.f + .015f * t913 + 8.115f * t923 + .838f * t9 +
    39.804f * t943 + 10.456f * t953) + 1.79e+0f * t9m32 * ex(
    -2.827f / t9) + 1.71e+3f * t9m32 * ex(-5.178f / t9) + 4.49e+6f *
    t935 * ex(-8.596f / t9);
  //
  //.......B12(a,n)N15................(Wagoner 1969)
  f(87) = 3.04e+15f * t9m23 * ex(-28.45f / t913);
  //
  //.......C13(a,n)O16................(Caughlan-Fowler 1988)
  f(88) = 6.77e+15f * t9m23 * ex(-32.329f / t913 - pow2((t9 /
    1.284f))) * (1.f + .013f * t913 + 2.04f * t923 + .184f * t9) +
    3.82e+5f * t9m32 * ex(-9.373f / t9) + 1.41e+6f * t9m32 * ex(
    -11.873f / t9) + 2.00e+9f * t9m32 * ex(-20.409f / t9) +
    2.92e+9f * t9m32 * ex(-29.283f / t9);
  //
  //----------REFERENCES-----------------------------------------------------------
  //     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data
  //       Tables, 40, 283.
  //     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.
  //
}

//
//========================IDENTIFICATION DIVISION================================
//
void
common::derivs(
  common& cmn,
  int const& loop)
{
	/*
  float& t9 = cmn.t9;
  float& hv = cmn.hv;
  const int nnuc = 26;
  arr_cref<float> y(static_cast<common_evolp1&>(cmn).y, dimension(nnuc));
  float& dt9 = cmn.dt9;
  arr_cref<float> dydt(cmn.dydt, dimension(nnuc));
  float& dlt9dt = cmn.dlt9dt;
  arr_cref<float> thm(cmn.thm, dimension(14));
  float& hubcst = cmn.hubcst;
  arr_cref<float> zm(cmn.zm, dimension(nnuc));
  arr_cref<float> dm(cmn.dm, dimension(nnuc));
  */
  //
  const float pi = boost::math::constants::pi<float>();
  float sumy = 0;
  float sumzy = 0;
  float sumdy = 0;
  float summdy = 0;
  float sumzdy = 0;
  int i = 0;
  float dphdt9 = 0;
  float dphdln = 0;
  float dphdzy = 0;
  float bar = 0;
  float dlndt9 = 0;
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] driver
  //     CALLS     - [subroutine] therm, rate1, rate4, rate3, rate2, sol
  //
  //----------REMARKS.
  //     Computes derivatives of
  //       - Temperature
  //       - hv
  //       - Chemical potential
  //       - abundances
  //
  //----------PARAMETERS.
  //Number of variables to be evolved.
  //Number of nuclides in calculation.
  //
  //----------COMMON AREAS.
  //Evolution paramete
  //Model parameters.
  //Time variables.
  //Dynamic variables.
  //Energy densities.
  //Nuclide data.
  //Flags,counters.
  //Run options.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------EVOLUTION PARAMETERS.
  //Temperature (in units of 10**9 K).
  //Defined by hv = M(atomic)n(baryon)/t
  //Chemical potential for electron.
  //Relative number abundances.
  //
  //----------EVOLUTION PARAMETERS (DERIVATIVES).
  //Change in temperature.
  //Change in hv.
  //Change in chemical potential.
  //Change in rel number abundances.
  //
  //----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
  //Rel # abund at beginning of iteratio
  //
  //----------MODEL PARAMETERS.
  //Gravitational constant.
  //Cosmological constant.
  //
  //----------TIME VARIABLES.
  //(1/t9)*d(t9)/d(t).
  //
  //----------DYNAMIC VARIABLES.
  //Thermodynamic variables.
  //Expansion rate.
  //
  //----------ENERGY DENSITIES.
  //Initial baryon mass density.
  //Baryon mass density.
  //Baryon mass density (ratio to init v
  //
  //----------NUCLIDE DATA.
  //Charge of nuclide.
  //Mass excess of nuclide.
  //
  //----------COUNTERS AND FLAGS.
  //Indicates if gaussian elimination fa
  //
  //----------RUN OPTION.
  //Run network size.
  //Number of nuclides in computation.
  //
  //----------SUMS.
  //Sum of abundances.
  //Sum of charge*abundances.
  //Sum of abundance flows.
  //Sum of (mass excess)*(abundance flow
  //Sum of (charge)*(abundance flows).
  //
  //----------DERIVATIVES.
  //d(phi e)/d(t9).
  //d(phi e)/d(h).
  //d(phi e)/d(sumzy).
  //(1/h)*d(h)/d(t9).
  //Baryon density and pressure terms.
  //
  //----------LOCAL VARIABLES.
  //Counts which Runge-Kutta loop.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------COMPUTE DERIVATIVES FOR ABUNDANCES-----------------------------------
  //
  //Baryon mass density (ratio to init v
  cmn.rnb = hv * t9 * t9 * t9 / cmn.rhob0;
  //..........VARIOUS THERMODYNAMIC QUANTITIES.
  therm(cmn);
  //Expansion rate.
  hubcst = sqrt((8.f / 3.f) * pi * cmn.g * (thm(10)) + (cmn.cosmo / 3.f));
  //Baryon mass density.
  cmn.rhob = thm(9);
  //..........COMPUTE REACTION RATE COEFFICIENTS.
  rate1(cmn, t9);
  //Run network selection.
  switch (cmn.irun) {
    case 1: goto statement_100;
    case 2: goto statement_110;
    case 3: goto statement_120;
    default: break;
  }
  statement_100:
  //Forward rate for all of reactions.
  rate4(cmn);
  statement_110:
  //Forward rate for reactions with A <
  rate3(cmn);
  statement_120:
  //Forward rate for reactions with A <
  rate2(cmn);
  //..........SOLVE COUPLED DIFFERENTIAL EQUATIONS.
  sol(cmn, loop);
  //Abort in case matrix not invertible.
  if (cmn.mbad > 0) {
    return;
  }
  //
  //20--------COMPUTE DERIVATIVES FOR TEMPERATURE, hv, AND CHEMICAL POTENTIAL------
  //
  //..........INITIALIZE SUMS TO ZERO.
  sumy = 0.f;
  sumzy = 0.f;
  sumdy = 0.f;
  summdy = 0.f;
  sumzdy = 0.f;
  //..........ACCUMULATE TO GET SUM.
  FEM_DO_SAFE(i, 1, cmn.isize) {
    sumy += y(i); 					/// Sum of abundance.
    sumzy += zm[i-1] * y(i); 			/// Sum of charge*abundance.
    sumdy += dydt(i); 				/// Sum of abundance flow.
    summdy += dm[i-1] * dydt(i); 		/// Sum of (mass excess)*(abundanc
    sumzdy += zm[i-1] * dydt(i); 		/// Sum of (charge)*(abundance flo
  }
  //..........CHANGES IN TEMPERATURE, hv, AND CHEMICAL POTENTIAL.
  dphdt9 = thm(12) * (-1.070e-4f * hv * sumzy / t9 - thm(11));
  dphdln = -thm(12) * 3.568e-5f * hv * sumzy;
  dphdzy = thm(12) * 3.568e-5f * hv;
  bar = 9.25e-5f * t9 * sumy + 1.388e-4f * t9 * sumdy / (3.f *
    hubcst) + summdy / (3.f * hubcst);
  //(Ref 1)
  dlndt9 = -(thm(2) + thm(5) + thm(6) * dphdt9 + thm(9) * 1.388e-4f *
    sumy) / (thm(1) + thm(3) + thm(4) + thm(7) + thm(9) * bar + thm(
    6) * (dphdln + dphdzy * sumzdy / (3.f * hubcst)));
  dt9 = (3.f * hubcst) / dlndt9;
  dlt9dt = dt9 / t9;
  //(Ref 2)
  cmn.dhv = -hv * 3*(hubcst + dlt9dt);
  //(Ref 3)
  cmn.dphie = dphdt9 * dt9 + dphdln * (3.f * hubcst) + dphdzy * sumzdy;
  //
  //----------REFERENCES-----------------------------------------------------------
  //     1)  Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
  //          Kellogg Radiation Lab preprint OAP-714,
  //          equation D.35.
  //     2)  Kawano, L., 1992, preprint, equation D.19.
  //     3)  Kawano, L., 1992, preprint, equation D.20.
  //
}

//
//========================IDENTIFICATION DIVISION================================
//
void
common::accum(
  common& cmn)
{
	/*
  // COMMON evolp1
  const int nnuc = 26;
  arr_cref<float> y(static_cast<common_evolp1&>(cmn).y, dimension(nnuc));
  // COMMON thermcb
  arr_cref<float> thm(cmn.thm, dimension(14));
  // COMMON nucdat
  arr_cref<float> am(cmn.am, dimension(nnuc));
  // COMMON flags
  int& it = cmn.it;
  // COMMON outdat
  const int itmax = 40;
  arr_ref<float, 2> xout(cmn.xout, dimension(itmax, nnuc));
  arr_ref<float, 2> thmout(cmn.thmout, dimension(itmax, 6));
  arr_ref<float> t9out(cmn.t9out, dimension(itmax));
  arr_ref<float> tout(cmn.tout, dimension(itmax));
  arr_ref<float> dtout(cmn.dtout, dimension(itmax));
  arr_ref<float> etaout(cmn.etaout, dimension(itmax));
  arr_ref<float> hubout(cmn.hubout, dimension(itmax));
  */
  //
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] driver
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Output accumulator.
  //
  //----------PARAMETERS.
  //Number of variables to be evolved.
  //Number of nuclides in calculation.
  //Maximum # of lines to be printed.
  //
  //----------COMMON AREAS.
  //Evolution paramete
  //Computation parame
  //Time variables.
  //Dynamic variables.
  //Nuclide data.
  //Flags,counters.
  //Output data.
  //Run options.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------EVOLUTION PARAMETERS.
  //Temperature (in units of 10**9 K).
  //Defined by hv = M(atomic)n(baryon)/t
  //Chemical potential for electron.
  //Relative number abundances.
  //
  //----------COMPUTATION PARAMETERS.
  //Accumulation increment.
  //
  //----------TIME PARAMETERS.
  //Time.
  //Time step.
  //
  //----------DYNAMIC VARIABLES.
  //Thermodynamic variables.
  //Expansion rate.
  //
  //----------NUCLIDE DATA.
  //Atomic number of nuclide.
  //
  //----------COUNTERS AND FLAGS.
  //Indicates if output buffer printed.
  //# times accumulated in output buffer
  //# time steps after outputting a line
  //
  //----------OUTPUT ARRAYS.
  //Nuclide mass fractions.
  //Thermodynamic variables.
  //Temperature (in units of 10**9 K).
  //Time.
  //Time step.
  //Baryon-to-photon ratio.
  //Expansion rate.
  //
  //----------RUN OPTION.
  //Number of nuclides in computation.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //Set up accumulation counter.
  it++;
  //
  //10--------SET UP OUTPUT VARIABLES----------------------------------------------
  //
  //..........DIVIDE NUMBER FRACTION BY THAT OF PROTON.
  int i = 0;
  FEM_DO_SAFE(i, 1, cmn.isize) {
    xout(it, i) = y(i) / y(2);
  }
  //Exception for proton.
  //xout(it, 2) = y(2) * am(2);
  xout(it, 2) = y(2) * am[1];
  //Exception for helium.
  //xout(it, 6) = y(6) * am(6);
  xout(it, 6) = y(6) * am[5];
  //..........SUM UP ABUNDANCES OF HEAVY NUCLIDES.
  //Li8 to O16.
  xout(it, 10) += xout(it, 11) + xout(it, 12) + xout(it, 13) + xout(it,
    14) + xout(it, 15) + xout(it, 16) + xout(it, 17) + xout(it, 18) + xout(it,
    19) + xout(it, 20) + xout(it, 21) + xout(it, 22) + xout(it, 23) + xout(it,
    24) + xout(it, 25) + xout(it, 26);
  //..........RELABEL TEMPERATURE, TIME, THERMODYNAMIC VARIABLES, ETC.
  //Temperature.
  t9out(it) = cmn.t9;
  //Time.
  tout(it) = cmn.t;
  //rho photon.
  thmout(it, 1) = thm(1);
  //rho electron.
  thmout(it, 2) = thm(4);
  //rho neutrino.
  thmout(it, 3) = thm(8);
  //rho baryon.
  thmout(it, 4) = thm(9);
  //Chemical potential.
  thmout(it, 5) = cmn.phie;
  //rho total.
  thmout(it, 6) = thm(10);
  //Time step.
  dtout(it) = cmn.dt;
  //Baryon to photon ratio.
  etaout(it) = cmn.hv / (3.3683e+4f);
  //Expansion rate.
  hubout(it) = cmn.hubcst;
  //
  //20--------INDICATE TERMINATION OF ACCUMULATION IF APPROPRIATE------------------
  //
  if ((it == itmax) || (cmn.ip < cmn.inc)) {
    cmn.ltime = 1;
  }
  //
}

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//     Changes (to run inder DEC unix f77):
//     -----------------------------------
//     COMMON /bessel/ -> COMMON /besselcb/
//     COMMON /therm/  -> COMMON /thermcb/
//     COMMON /time/   -> COMMON /ttime/
//     ir=1 -> ir=5
//     iw=1 -> iw=6
//     All `entry' routines removed
//
//========================IDENTIFICATION DIVISION================================
//
void
common::driver(
  common& cmn)
{
  const int nnuc = 26;
  /*
  arr_cref<float> y(static_cast<common_evolp1&>(cmn).y, dimension(nnuc));
  arr_cref<float> dydt(cmn.dydt, dimension(nnuc));
  float& ytmin = cmn.ytmin;
  int& inc = cmn.inc;
  float& t = cmn.t;
  float& dt = cmn.dt;
  float& dlt9dt = cmn.dlt9dt;
  int& ltime = cmn.ltime;
  int& is = cmn.is;
  int& ip = cmn.ip;
  int& itime = static_cast<common_tcheck&>(cmn).itime;
  int& isize = cmn.isize;
  */
  //
  int mvar = 0;
  int loop = 0;
  const float cl = 1.e-16f;
  float dtmin = 0;
  int i = 0;
  float dtl = 0;
  const int nvar = 29;
  arr_1d<nvar, float> v(fem::fill0);
  arr_1d<nvar, float> v0(fem::fill0);
  arr_1d<nvar, float> dvdt(fem::fill0);
  arr_1d<nvar, float> dvdt0(fem::fill0);
  //
  //----------LINKAGES.
  //     CALLED BY - [subroutine] run
  //     CALLS     - [subroutine] start, derivs, accum
  //
  //----------REMARKS.
  //     Runge-Kutta computational routine
  //
  //----------PARAMETERS.
  //Number of variables to be evolved.
  //Number of nuclides in calculation.
  //Lower limit on size of time step.
  //
  //----------COMMON AREAS.
  //Evolution paramete
  //Computation parame
  //Time variables.
  //Flags,counters.
  //Computation locat
  //Run options.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------EVOLUTION PARAMETERS.
  //Temperature (in units of 10**9 K).
  //Defined by hv = M(atomic)n(baryon)/t
  //Chemical potential for electron.
  //Relative number abundances.
  //
  //----------EVOLUTION PARAMETERS (DERIVATIVES).
  //Change in rel number abundances.
  //
  //----------EVOLUTION PARAMETERS (ORIGINAL VALUES).
  //Rel # abund at beginning of iteratio
  //
  //----------COMPUTATION PARAMETERS.
  //Time step limiting constant on abund
  //Time step limiting constant on tempe
  //Final temperature (in 10**9 K).
  //Smallest abundances allowed.
  //Accumulation increment.
  //
  //----------TIME AND TIME STEP VARIABLES.
  //Time.
  //Time step.
  //(1/t9)*d(t9)/d(t).
  //
  //----------COUNTERS AND FLAGS.
  //Counts which Runge-Kutta loop.
  //Indicates termination status.
  //# total time steps for particular ru
  //# time steps after outputting a line
  //
  //----------COMPUTATION LOCATION.
  //Time check.
  //
  //----------RUN OPTION.
  //Number of nuclides in computation.
  //
  //----------TIME AND TIME STEP VARIABLES.
  //Mininum time step.
  //Time step from limitation on abund c
  //
  //----------LABELS FOR VARIABLES TO BE TIME EVOLVED.
  //Total number of variables to be evol
  //Variables to be time evolved.
  //Time derivatives.
  //Value of variables at original point
  //Value of derivatives at original poi
  //
  //----------EQUIVALENCE STATEMENTS.
  //     EQUIVALENCE (v(4),y(1)),(dvdt(4),dydt(1)),(v0(4),y0(1))
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------INPUT INITIALIZATION INFORMATION, RELABEL----------------------------
  //
  //Set termination indicator to zero.
  ltime = 0;
  //Input initialization information.
  start(cmn);
  //Total number of variables to be evol
  mvar = isize + 3;
  //
  //20--------LOOP ONE-------------------------------------------------------------
  //
  //Begin Runge-Kutta looping.
  statement_200:
  //Loop indicator.
  loop = 1;
  //..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
  derivs(cmn, loop);
  //Time = 1st R-K loop.
  itime = 4;
  //Check interface subroutine.
  check(cmn);
  //..........ACCUMULATE.
  //Low temp.
  //Small dt.
  //Enough iterations.
  if ((cmn.t9 <= cmn.t9f) || (dt < fem::abs(cl / dlt9dt)) || (ip == inc)) {
    accum(cmn);
  }
  //..........POSSIBLY TERMINATE COMPUTATION.
  //Return to run selection.
  if (ltime == 1) {
    return;
  }
  //..........RESET COUNTERS.
  //Reset iteration counters.
  if (ip == inc) {
    ip = 0;
  }
  ip++;
  is++;
  //..........ADJUST TIME STEP.
  //Adjust time step after 3 iterations.
  if (is > 3) {
    //Trial value for minimum time step (R
    dtmin = fem::abs(1.f / dlt9dt) * cmn.ct;
    //Go through all abundance changes.
    FEM_DO_SAFE(i, 1, isize) {
      if ((dydt(i) != 0.f) && (y(i) > ytmin)) {
        //(Ref 2).
        dtl = fem::abs(y(i) / dydt(i)) * cmn.cy * (1.f + pow2((
          fem::alog10(y(i)) / fem::alog10(ytmin))));
        //Find smallest time st
        if (dtl < dtmin) {
          dtmin = dtl;
        }
      }
    }
    //Limit change in time
    if (dtmin > 1.5f * dt) {
      dtmin = 1.5f * dt;
    }
    //Set new time step.
    dt = dtmin;
  }
  //Increment time.
  t += dt;
  //..........STORE AND INCREMENT VALUES (Ref 3).
  FEM_DO_SAFE(i, 1, mvar) {
    v0(i) = v(i);
    dvdt0(i) = dvdt(i);
    v(i) = v0(i) + dvdt0(i) * dt;
    //Set at minimum
    if ((i >= 4) && (v(i) < ytmin)) {
      v(i) = ytmin;
    }
  }
  //
  //30--------LOOP TWO-------------------------------------------------------------
  //
  //Step up loop counter.
  loop = 2;
  //..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
  derivs(cmn, loop);
  //Time = 2nd R-K loop.
  itime = 7;
  //Check interface subroutine.
  check(cmn);
  //..........INCREMENT VALUES.
  FEM_DO_SAFE(i, 1, mvar) {
    v(i) = v0(i) + .5f * (dvdt(i) + dvdt0(i)) * dt;
    //Set at minimum
    if ((i >= 4) && (v(i) < ytmin)) {
      v(i) = ytmin;
    }
  }
  goto statement_200;
  //
  //----------REFERENCES-----------------------------------------------------------
  //     1)  Constraint on dt from the requirement that
  //                (d(t9)/dt)*(dt/t9) < ct
  //           Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 293, equation C6.
  //     2)  Constraint on dt from
  //                dtl < y/(dy/dt)*cy*(1+(log(y)/log(ytmin))**2)
  //          Wagoner, R.V. 1969, page 293, equation C7 but with log term squared.
  //     3)  Wagoner, R.V. 1969, page 292, equations C1, C2.
  //
}

/*
struct run_save
{
  arr<fem::str<22> > vtype;

  run_save() :
    vtype(dimension(8), fem::fill0)
  {}
};
*/

// Replaces the equivalence memory sharing used in the original Fortran. 
// TODO add to class
// TODO make private
void qvary(common& cmn, int index, float value)
{
	//----------EQUIVALENCE VARIABLE.
	//     REAL    qvary(7)             !Array set equal to c, cosmo, and xi.
	//
	//----------EQUIVALENCE STATEMENTS.
	//     EQUIVALENCE (qvary(1),c(1)), (qvary(4),cosmo), (qvary(5),xi(1))
	//
	if (index >= 1 && index <= 3)
		cmn.c(index) = value;
	else if (index == 4)
		cmn.cosmo = value;
	else if (index >= 5 && index <= 7)
		cmn.xi(index) = value;
	else
	{
		std::cerr << "index out of bounds." << std::endl;
		exit(1);
	}
}


//
//========================IDENTIFICATION DIVISION================================
//
void
common::run(common& cmn)
{
  //FEM_CMN_SVE(run);
  common_read read(cmn);
  common_write write(cmn);
  /*
  float& eta1 = cmn.eta1;
  //int& itime = static_cast<common_checkcb&>(cmn).itime;
  int& itime = cmn.itime;
  int& irun = cmn.irun;
  int& isize = cmn.isize;
  int& jsize = cmn.jsize;
  */
  //
  //str_arr_ref<1> vtype(sve.vtype, dimension(8));
  //if (is_called_first_time) {
    static const char* vtype[] = {
      	"baryon/photon ratio   ", 
		"gravitational constant",
        "neutron lifetime      ", 
		"# of neutrino species ",
        "cosmological constant ", 
		"xi-electron           ",
        "xi-muon               ", 
		"xi-tauon              "
	};
  //  fem::data_of_type_str(FEM_VALUES_AND_SIZE),
  //    vtype;
  //}
  const int iw = 6;
  const int ir = 5;
  int inumb = 0;
  const int nnuc = 26;
  const int nrec = 88;
  const int lnuc = 18;
  const int lrec = 64;
  const int knuc = 9;
  const int krec = 34;
  int jnum = 0;
  int knum = 0;
  int i = 0;
  arr_1d<3, float> rnum1(fem::fill0);
  arr_1d<3, float> rnum2(fem::fill0);
  arr_1d<3, float> rnum3(fem::fill0);
  arr_1d<3, int> inum(fem::fill0);
  fem::str<1> lchose;
  int l = 0;
  arr_1d<3, int> lnum(fem::fill0);
  int lnumb1 = 0;
  float rnumb1 = 0;
  int lnumb2 = 0;
  float rnumb2 = 0;
  int lnumb3 = 0;
  float rnumb3 = 0;
  static const char* format_1001 = "(i1)";
  static const char* format_2200 = "(' ','Begin computation run....')";
  static const char* format_2202 =
    "(' ','Computation completed - Press <RETURN> to ','continue: ',$)";
  //
  //----------LINKAGES.
  //     CALLED BY - [program] nuc123
  //     CALLS     - [subroutine] driver
  //
  //----------REMARKS.
  //     Activates computation routine.
  //
  //----------PARAMETERS.
  //Input unit number (previous value =
  //Output unit number (previous value =
  //Number of nuclear reactions.
  //Total # of nuclear reactions for iru
  //Number of nuclides in calculation.
  //Total # of nuclides for irun = 2.
  //Total # of nuclides for irun = 3.
  //
  //----------COMMON AREAS.
  //Model parameters.
  //Variational parame
  //Computation locati
  //Run options.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------MODEL PARAMETERS.
  //Baryon-to-photon ratio.
  //c(1) is variation of gravitational constant.
  //c(2) is neutron lifetime (sec).
  //c(3) is number of neutrino species.
  //Cosmological constant.
  //Neutrino degeneracy parameters.
  //
  //----------RUN OPTION.
  //Run network size.
  //Number of nuclides in computation.
  //Number of reactions in computation.
  //
  //----------USER INTERACTION VARIABLES.
  //Run parameter for outer loop.
  //Run parameter for middle loop.
  //Run parameter for inner loop.
  //Run parameter starting value.
  //Run parameter end value.
  //Run parameter increment.
  //Selection number.
  //Number of loops to be done.
  //Number of loops rejected.
  //Run parameter for outer loop.
  //Run parameter for middle loop.
  //Run parameter for inner loop.
  //Run parameter end value.
  //User response (alphanumeric).
  //
  //----------FLAG AND LABELS.
  //Computation location.
  //Label for quantities being varied.
  //
  //----------EQUIVALENCE VARIABLE.
  //     REAL    qvary(7)             !Array set equal to c, cosmo, and xi.
  //
  //----------EQUIVALENCE STATEMENTS.
  //     EQUIVALENCE (qvary(1),c(1)), (qvary(4),cosmo), (qvary(5),xi(1))
  //
  //==============================DATA DIVISION====================================
  //
  //----------LABELS FOR QUANTITIES BEING VARIED.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------PRINT RUN SELECTION AND AWAIT RESPONSE-------------------------------
  //
  //..........RETURN FROM LOOPING.
  statement_100:
  //..........DISPLAY RUN SELECTIONS.
  write(iw,
    "(8(/),32x,"
	"'RUN SELECTION',/,32x,"
	"'--- ---------',/,/,27x,"
    "' 1. SET RUN NETWORK',/,27x,"
	"' 2. GO',/,27x,"
    "' 3. DO MULTIPLE RUNS',/,27x,"
	"' 4. EXIT',10(/),27x,"
    "' Enter selection (1-4): ',$)");
  //..........READ IN SELECTION NUMBER.
  read(ir, format_1001), inumb;
  //
  //20--------BRANCH TO APPROPRIATE SECTION----------------------------------------
  //
  switch (inumb) {
    case 1: goto statement_210;
    case 2: goto statement_220;
    case 3: goto statement_230;
    case 4: goto statement_240;
    default: break;
  }
  //Improper input or <RETURN>.
  goto statement_240;
  //
  //21--------SET RUN NETWORK SECTION----------------------------------------------
  //
  statement_210:
  write(iw,
    "(' ','Enter network size (1-26 nuclides (default); ','2-18; 3-9): ',$)");
  //Read in selection number.
  read(ir, star), inumb;
  //
  if ((inumb != 1) && (inumb != 2) && (inumb != 3)) {
    inumb = 1;
  }
  //Run network changed from previously.
  if (inumb != irun) {
    //Run network size selection.
    irun = inumb;
  }
  //Maximal network size.
  if (irun == 1) {
    isize = nnuc;
    jsize = nrec;
  }
  else {
    //Abridged network size.
    if (irun == 2) {
      isize = lnuc;
      jsize = lrec;
    }
    else {
      //Minimal network size.
      if (irun == 3) {
        isize = knuc;
        jsize = krec;
      }
    }
    //(irun.eq.1)
  }
  write(iw,
    "(' ','Run network set to ',i1,' - Press <RETURN> ','to continue: ',$)"),
    irun;
  read(ir, star);
  goto statement_300;
  //
  //22--------GO SECTION-----------------------------------------------------------
  //
  statement_220:
  write(iw, format_2200);
  itime = 3;
  //Call interface subr before computati
  check(cmn);
  //Do nucleosynthesis computation.
  driver(cmn);
  itime = 8;
  //Call interface subr after computatio
  check(cmn);
  write(iw, format_2202);
  read(ir, star);
  goto statement_300;
  //
  //23--------DO MULTIPLE RUNS SECTION---------------------------------------------
  //
  //..........GET NUMBER OF LOOPINGS.
  statement_230:
  write(iw,
    "(' ','Enter the number of loops to be done (1 ','(default); 2; 3): ',"
    "$)");
  //Read in number of loops to be don
  read(ir, star), jnum;
  if ((jnum != 1) && (jnum != 2) && (jnum != 3)) {
    //Default number of loops.
    jnum = 1;
  }
  //No loops rejected for now.
  knum = 0.f;
  FEM_DO_SAFE(i, 1, 3) {
    if (i > jnum) {
      //Initialize initial parameter.
      rnum1(i) = 0.f;
      //Initialize terminal parameter.
      rnum2(i) = 0.f;
      //Initialize incremental parameter.
      rnum3(i) = 1.f;
      //Initialize selection number.
      inum(i) = 0;
    }
    else {
      //..........OBTAIN QUANTITY TO VARY.
      write(iw,
        "(8(/),30x,'QUANTITY TO VARY',/,30x,'-------- -- ----',/,/,"
        "' ',25x,' 1.  ETA     (LOGRITHMIC VARIATION)',/,25x,"
        "' 2.  G           (LINEAR VARIATION)',/,25x,"
        "' 3.  TAU         (LINEAR VARIATION)',/,25x,"
        "' 4.  # NEUTRINOS (LINEAR VARIATION)',/,25x,"
        "' 5.  LAMBDA      (LINEAR VARIATION)',/,25x,"
        "' 6.  XI-ELECTRON (LINEAR VARIATION)',/,25x,"
        "' 7.  XI-MUON     (LINEAR VARIATION)',/,25x,"
        "' 8.  XI-TAUON    (LINEAR VARIATION)',/,25x,"
		"' 9.  NO SELECTION',"
        "5(/),25x,' Enter selection (1-9): ',$)");
      read(ir, format_1001), inum(i);
      //No selection ma
      if ((inum(i) < 1) || (inum(i) > 8)) {
        write(iw,
          "(' ','No selection made - Reduce number of ','loops by one',/,"
          "' ','Press <RETURN> to continue: ',$)");
        read(ir, star);
        //Step up number of loops rejected.
        knum++;
        //Initialize initial parameter.
        rnum1(i) = 0.f;
        //Initialize terminal parameter.
        rnum2(i) = 0.f;
        //Initialize incremental parameter.
        rnum3(i) = 1.f;
        //Initialize selection number.
        inum(i) = 0;
        //((inum(i).ge.1).and.(inum(i).le.8))
      }
      else {
        //..........INPUT RUN SPECIFICATIONS.
        statement_231:
        write(iw, "(' ','Enter minimum value: ',$)");
        //Read in starting value.
        read(ir, star), rnum1(i);
        write(iw, "(' ','Enter maximum value: ',$)");
        //Read in terminating value.
        read(ir, star), rnum2(i);
        statement_232:
        write(iw, "(' ','Enter increment: ',$)");
        //Read in incremental value.
        read(ir, star), rnum3(i);
        //Trouble with 0 division later on
        if (rnum3(i) == 0.f) {
          write(iw,
            "(' ','Zero increment not allowed: trouble with dividing by zero')");
          goto statement_232;
        }
        //Display inp
        write(iw,
          "(' ','Run from ',1p,e12.5,' to ',1p,e12.5,' in increments of ',1p,"
          "e12.5)"),
          rnum1(i), rnum2(i), rnum3(i);
        write(iw, "(' ','Confirm these values (Y or N): ',$)");
        //Get confirm
        read(ir, "(a1)"), lchose;
        if ((lchose != "Y") && (lchose != "y")) {
          goto statement_231;
        }
        //((inum(i).lt.1).or.(inum(i).gt.8))
      }
      //(i.gt.jnum)
    }
    //i = 1,3
  }
  //Number of valid loops.
  jnum = jnum - knum;
  //Run requested.
  if (jnum != 0) {
    //..........WRITE OUT QUANTITY TO VARY, RUN SPECIFICATIONS.
    //Check all loops.
    FEM_DO_SAFE(l, 1, jnum + knum) {
      //Proper selection was made.
      if (inum(l) != 0) {
        //Display run
        write(iw,
          "(' ','Run ',a22,/,'    from ',1p,e12.5,' to ',1p,e12.5,"
          "' in increments of ',1p,e12.5)"),
          //vtype(inum(l)), rnum1(l), rnum2(l), rnum3(l);
          vtype[inum(l)-1], rnum1(l), rnum2(l), rnum3(l);
        //..........GET LOGS OF eta VALUES FOR LOGARITHMIC INCREMENTATION.
        //Work with exponents for eta incre
        if (inum(l) == 1) {
          rnum1(l) = log10(rnum1(l));
          rnum2(l) = log10(rnum2(l));
        }
      }
    }
    //..........COMPUTE NUMBER OF RUNS FOR EACH LOOPING.
    FEM_DO_SAFE(l, 1, 3) {
      lnum(l) = fem::nint((rnum2(l) - rnum1(l) + rnum3(l)) / rnum3(l));
    }
    //..........DO MULTIPLE RUNS.
    //Inform user of beginning of computation.
    write(iw, format_2200);
    //Outer loop.
    FEM_DO_SAFE(lnumb1, 0, lnum(1) - 1) {
      //Value of param fo
      rnumb1 = rnum1(1) + fem::ffloat(lnumb1) * rnum3(1);
      if ((inum(1) >= 1) && (inum(1) <= 8)) {
        if (inum(1) == 1) {
          //Vary baryon-to-photon ratio.
          eta1 = pow(10, rnumb1);
        }
        else {
          //Vary other quantities.
		  qvary(cmn, inum(1)-1, rnumb1);
        }
      }
      //Middle loop.
      FEM_DO_SAFE(lnumb2, 0, lnum(2) - 1) {
        //Value of param
        rnumb2 = rnum1(2) + fem::ffloat(lnumb2) * rnum3(2);
        if ((inum(2) >= 1) && (inum(2) <= 8)) {
          if (inum(2) == 1) {
            //Vary baryon-to-photon ratio.
            eta1 = pow(10, rnumb2);
          }
          else {
            //Vary other quantities.
		  	qvary(cmn, inum(2)-1, rnumb2);
          }
        }
        //Inner loop.
        FEM_DO_SAFE(lnumb3, 0, lnum(3) - 1) {
          //Value of parameter.
          rnumb3 = rnum1(3) + fem::ffloat(lnumb3) * rnum3(3);
          if ((inum(3) >= 1) && (inum(3) <= 8)) {
            if (inum(3) == 1) {
              //Vary baryon-to-photon ratio.
              eta1 = pow(10, rnumb3);
            }
            else {
              //Vary other quantities.
		  	  qvary(cmn, inum(3)-1, rnumb3);
            }
          }
          itime = 3;
          //Check interface subroutine before computation.
          check(cmn);
          //Do nucleosynthesis computation.
          driver(cmn);
          itime = 8;
          //Check interface subroutine after computation.
          check(cmn);
          //lnumb3 = 0,lnum(3)-1
        }
        //lnumb2 = 0,lnum(2)-1
      }
      //lnumb1 = 0,lnum(1)-1
    }
    //Inform user of completion of computation.
    write(iw, format_2202);
    //(jnum == 0)
  }
  else {
    write(iw, "(' ','No selection made - ','Press <RETURN> to continue: ',$)");
    //(jnum != 0)
  }
  read(ir, star);
  goto statement_300;
  //
  //24--------EXIT SECTION---------------------------------------------------------
  //
  statement_240:
  return;
  //
  //30--------GO BACK TO MENU------------------------------------------------------
  //
  statement_300:
  goto statement_100;
  //
}

//
//========================IDENTIFICATION DIVISION================================
//
void
output(common& cmn)
{
  common_read read(cmn);
  common_write write(cmn);
  float& cy = cmn.cy;
  float& ct = cmn.ct;
  float& t9i = cmn.t9i;
  float& t9f = cmn.t9f;
  float& ytmin = cmn.ytmin;
  arr_cref<float> c(cmn.c, dimension(3));
  float& cosmo = cmn.cosmo;
  arr_cref<float> xi(cmn.xi, dimension(3));
  int& it = cmn.it;
  const int itmax = 40;
  const int nnuc = 26;
  arr_cref<float, 2> xout(cmn.xout, dimension(itmax, nnuc));
  arr_cref<float, 2> thmout(cmn.thmout, dimension(itmax, 6));
  arr_cref<float> t9out(cmn.t9out, dimension(itmax));
  arr_cref<float> tout(cmn.tout, dimension(itmax));
  arr_cref<float> dtout(cmn.dtout, dimension(itmax));
  arr_cref<float> etaout(cmn.etaout, dimension(itmax));
  arr_cref<float> hubout(cmn.hubout, dimension(itmax));
  int& nout = cmn.nout;
  //
  const int iw = 6;
  const int ir = 5;
  int inum = 0;
  int j = 0;
  int i = 0;
  static const char* format_1001 = "(i1)";
  static const char* format_2014 = "(/,/,/)";
  static const char* format_3100 =
    "(' ','Computational parameters:',/,'   cy = ',f5.3,'/ ct = ',f5.3,"
    "'/ initial temp = ',1p,e8.2,'/ final temp = ',1p,e8.2,/,"
    "'   smallest abundances allowed = ',1p,e8.2)";
  static const char* format_3102 =
    "(' ','Model parameters:',/,'   g = ',f5.2,'/ tau = ',f6.2,"
    "'/ # nu = ',f5.2,'/ lambda = ',1p,e10.3,/,'   xi-e = ',e10.3,"
    "'/ xi-m = ',e10.3,'/ xi-t = ',e10.3,/)";
  static const char* format_3106 = "(1p,e10.3,1p,5e12.3)";
  static const char* format_3108 = "(' ','Press <RETURN> to continue: ',$)";
  //
  //----------LINKAGES.
  //     CALLED BY - [program] nuc123
  //     CALLS     - none
  //
  //----------REMARKS.
  //     Outputs computational results either into an output file or onto
  //     the screen
  //
  //----------PARAMETERS.
  //Input unit number (previous value =
  //Output unit number (previous value =
  //Number of nuclides in calculation.
  //Maximum # of line to be printed.
  //
  //----------COMMON AREAS.
  //Computation parame
  //Model parameters.
  //Flags, counters.
  //Output data.
  //Output option.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------COMPUTATION SETTINGS.
  //Time step limiting constant on abund
  //Time step limiting constant on tempe
  //Initial temperature (in 10**9 K).
  //Final temperature (in 10**9 K).
  //Smallest abundances allowed.
  //
  //----------EARLY UNIVERSE MODEL PARAMETERS.
  //c(1) is variation of gravitational c
  //c(2) is neutron lifetime (sec).
  //c(3) is number of neutrino species.
  //Cosmological constant.
  //Neutrino degeneracy parameters.
  //
  //----------COUNTER.
  //# times accumulated in output buffer
  //
  //----------OUTPUT ARRAYS.
  //Nuclide mass fractions.
  //Thermodynamic variables.
  //Temperature (in units of 10**9 K).
  //Time.
  //Time step.
  //Baryon-to-photon ratio.
  //Expansion rate.
  //
  //----------OUTPUT FILE STATUS.
  //Number of output requests.
  //Indicates if output file used.
  //
  //----------USER INTERACTION VARIABLES.
  //Selection number.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------PRINT OUTPUT SELECTION AND AWAIT RESPONSE----------------------------
  //
  //..........RETURN FROM LOOPING.
  statement_100:
  //..........DISPLAY OUTPUT SELECTIONS.
  write(iw,
    "(8(/),30x,'OUTPUT SELECTION',/,30x,'------ ---------',/,/,"
    "25x,' 1. REQUEST OUTPUT FILE',/,25x,' 2. REQUEST OUTPUT ON SCREEN',/,"
    "' ',25x,' 3. EXIT',11(/),25x,' Enter selection (1-3): ',$)");
  //..........READ IN SELECTION NUMBER.
  read(ir, format_1001), inum;
  //..........BRANCH TO APPROPRIATE SECTION.
  switch (inum) {
    case 1: goto statement_200;
    case 2: goto statement_300;
    case 3: goto statement_400;
    default: break;
  }
  //Improper input or <RETURN>.
  goto statement_400;
  //
  //20--------REQUEST OUTPUT SECTION-----------------------------------------------
  //
  statement_200:
  //      DO j = 1,it                  !Temperature in MeV.
  //        t9out(j) = t9out(j)*.08617
  //      END DO
  //      DO j = 1,it                  !Energy density as fraction of total.
  //        thmout(j,1) = thmout(j,1)/thmout(j,6)  !Rhog.
  //        thmout(j,2) = thmout(j,2)/thmout(j,6)  !Rhoe.
  //        thmout(j,3) = thmout(j,3)/thmout(j,6)  !Rhone.
  //        thmout(j,4) = thmout(j,4)/thmout(j,6)  !Rhob.
  //      END DO
  //..........PRINT CAPTION.
  //Keep track of number of output reque
  nout++;
  if (nout == 1) {
    write(2,
      "(54x,'NUCLIDE ABUNDANCE YIELDS',/,54x,'------- --------- ------',/,/)");
  }
  write(2,
    "(' Computational parameters:',/,'   cy = ',f5.3,'/  ct = ',f5.3,"
    "'/  initial temp = ',1p,e8.2,'/  final temp = ',1p,e8.2,"
    "'/  smallest abundances allowed = ',1p,e8.2)"),
    cy, ct, t9i, t9f, ytmin;
  write(2,
    "(' Model parameters:',/,'   g = ',f5.2,'/  tau = ',f6.2,'/  # nu = ',"
    "f5.2,'/  lambda = ',1p,e10.3,'/  xi-e = ',e10.3,'/  xi-m = ',e10.3,"
    "'/  xi-t = ',e10.3,/)"),
    c(1), c(2), c(3), cosmo, xi(1), xi(2), xi(3);
  //..........PRINT HEADINGS, ABUNDANCES FOR NEUTRON TO LI8.
  write(2,
    "(4x,'Temp',8x,'N/H',10x,'P',10x,'D/H',9x,'T/H',8x,'He3/H',8x,'He4',8x,"
    "'Li6/H',7x,'Li7/H',7x,'Be7/H',6x,'Li8/H&up',/,132('-'))");
  FEM_DO_SAFE(j, 1, it) {
    {
      write_loop wloop(cmn, 2, "(1p,e10.3,1p,10e12.3)");
      wloop, t9out(j);
      FEM_DO_SAFE(i, 1, 10) {
        wloop, xout(j, i);
      }
    }
  }
  //..........PRINT THERMODYNAMIC QUANTITIES.
  write(2,
    "(' ',/,4x,'Temp',9x,'T',10x,'rhog',8x,'rhoe',7x,'rhone',8x,'rhob',8x,"
    "'phie',9x,'dt',9x,'eta',10x,'H',/,132('-'))");
  FEM_DO_SAFE(j, 1, it) {
    {
      write_loop wloop(cmn, 2, "(1p,e10.3,9e12.3)");
      wloop, t9out(j), tout(j);
      FEM_DO_SAFE(i, 1, 5) {
        wloop, thmout(j, i);
      }
      wloop, dtout(j), etaout(j), hubout(j);
    }
  }
  write(2, format_2014);
  //Output file requested.
  cmn.outfile = true;
  write(iw,
    "(' ','Output file requested - Press <RETURN> to ','continue: ',$)");
  read(ir, star);
  goto statement_500;
  //
  //30--------REQUEST OUTPUT ON SCREEN SECTION-------------------------------------
  //
  //..........RETURN FROM LOOPING.
  statement_300:
  //      DO j = 1,it                  !Temperature in MeV.
  //        t9out(j) = t9out(j)*.08617
  //      END DO
  //      DO j = 1,it                  !Energy density as fraction of total.
  //        thmout(j,1) = thmout(j,1)/thmout(j,6)  !Rhog.
  //        thmout(j,2) = thmout(j,2)/thmout(j,6)  !Rhoe.
  //        thmout(j,3) = thmout(j,3)/thmout(j,6)  !Rhone.
  //        thmout(j,4) = thmout(j,4)/thmout(j,6)  !Rhob.
  //      END DO
  //..........DISPLAY SCREEN OUTPUT SELECTIONS.
  write(iw,
    "(8(/),26x,'SCREEN OUTPUT SELECTION',/,26x,"
    "'------ ------ ---------',/,/,25x,' 1. DISPLAY D,T,HE3,HE4,LI7',/,"
    "' ',25x,' 2. DISPLAY N,P,LI6,BE7,LI8&UP',/,25x,"
    "' 3. DISPLAY RHOG,RHOE,RHONE,RHOB',/,25x,"
    "' 4. DISPLAY T,DT,PHIE,ETA,H',/,25x,' 5. EXIT',9(/),25x,"
    "' Enter selection (1-5): ',$)");
  //..........READ IN SELECTION NUMBER.
  read(ir, format_1001), inum;
  switch (inum) {
    case 1: goto statement_310;
    case 2: goto statement_320;
    case 3: goto statement_330;
    case 4: goto statement_340;
    case 5: goto statement_350;
    default: break;
  }
  //Improper input or <RETURN>.
  goto statement_350;
  //Display d,t,he3,he4,li7.
  statement_310:
  //..........PRINT CAPTION.
  write(iw, format_2014);
  write(iw, format_3100), cy, ct, t9i, t9f, ytmin;
  write(iw, format_3102), c(1), c(2), c(3), cosmo, xi(1), xi(2), xi(3);
  //..........PRINT HEADINGS, ABUNDANCES FOR D,T,HE3,HE4,LI7.
  write(iw,
    "(4x,'Temp',8x,'D/H',9x,'T/H',8x,'He3/H',8x,'He4',8x,'Li7/H',/,"
    "80('-'))");
  FEM_DO_SAFE(j, 1, it) {
    {
      write_loop wloop(cmn, iw, format_3106);
      wloop, t9out(j);
      FEM_DO_SAFE(i, 3, 6) {
        wloop, xout(j, i);
      }
      wloop, xout(j, 8);
    }
  }
  write(iw, format_2014);
  write(iw, format_3108);
  read(ir, star);
  goto statement_360;
  //Display n,p,li6,be7,li8&up.
  statement_320:
  //..........PRINT CAPTION.
  write(iw, format_2014);
  write(iw, format_3100), cy, ct, t9i, t9f, ytmin;
  write(iw, format_3102), c(1), c(2), c(3), cosmo, xi(1), xi(2), xi(3);
  //..........PRINT HEADINGS, ABUNDANCES FOR N,P,LI6,BE7,LI8&UP.
  write(iw,
    "(4x,'Temp',8x,'N/H',10x,'P',9x,'Li6/H',7x,'Be7/H',6x,'Li8/H&up',/,"
    "80('-'))");
  FEM_DO_SAFE(j, 1, it) {
    {
      write_loop wloop(cmn, iw, format_3106);
      wloop, t9out(j);
      FEM_DO_SAFE(i, 1, 2) {
        wloop, xout(j, i);
      }
      wloop, xout(j, 7);
      FEM_DO_SAFE(i, 9, 10) {
        wloop, xout(j, i);
      }
    }
  }
  write(iw, format_2014);
  write(iw, format_3108);
  read(ir, star);
  goto statement_360;
  //Display rhog,rhoe,rhone,rhob.
  statement_330:
  //..........PRINT CAPTION.
  write(iw, format_2014);
  write(iw, format_3100), cy, ct, t9i, t9f, ytmin;
  write(iw, format_3102), c(1), c(2), c(3), cosmo, xi(1), xi(2), xi(3);
  //..........PRINT ENERGY DENSITIES.
  write(iw,
    "(4x,'Temp',8x,'rhog',8x,'rhoe',7x,'rhone',8x,'rhob',/,80('-'))");
  FEM_DO_SAFE(j, 1, it) {
    {
      write_loop wloop(cmn, iw, "(1p,e10.3,4e12.3)");
      wloop, t9out(j);
      FEM_DO_SAFE(i, 1, 4) {
        wloop, thmout(j, i);
      }
    }
  }
  write(iw, format_2014);
  write(iw, format_3108);
  read(ir, star);
  goto statement_360;
  //Display t,dt,phie,eta,hubcst.
  statement_340:
  //..........PRINT CAPTION.
  write(iw, format_2014);
  write(iw, format_3100), cy, ct, t9i, t9f, ytmin;
  write(iw, format_3102), c(1), c(2), c(3), cosmo, xi(1), xi(2), xi(3);
  //..........PRINT THERMODYNAMIC QUANTITIES.
  write(iw,
    "(4x,'Temp',8x,'time',8x,'phie',9x,'dt',9x,'eta',10x,'H',/,80('-'))");
  FEM_DO_SAFE(j, 1, it) {
    write(iw, "(1p,e10.3,5e12.3)"), t9out(j), tout(j), thmout(j, 5),
      dtout(j), etaout(j), hubout(j);
  }
  write(iw, format_2014);
  write(iw, format_3108);
  read(ir, star);
  goto statement_360;
  //Exit.
  statement_350:
  goto statement_500;
  statement_360:
  goto statement_300;
  //
  //40--------EXIT SECTION---------------------------------------------------------
  //
  statement_400:
  return;
  //
  //50--------GO BACK TO MENU------------------------------------------------------
  //
  statement_500:
  goto statement_100;
  //
}

/*
struct blockdata_unnamed_save
{
};
*/

//
//========================IDENTIFICATION DIVISION================================
//
common::common(
		int argc,
		char const* argv[])
		:
		fem::common(argc, argv),
		f(dimension(nrec)),
		r(dimension(nrec)),
		a(nrec,nrec),
		b(nrec,nrec),
		yx(nrec,nrec)
{
	/*
  FEM_CMN_SVE(blockdata_unnamed);
  // COMMON recpr0
  const int nrec = 88;
  arr_ref<float, 2> reacpr(cmn.reacpr, dimension(nrec, 8));
  // COMMON modpr0
  arr_ref<float> c0(cmn.c0, dimension(3));
  arr_ref<float> xi0(cmn.xi0, dimension(3));
  // COMMON nucdat
  const int nnuc = 26;
  arr_ref<float> am(cmn.am, dimension(nnuc));
  arr_ref<float> zm(cmn.zm, dimension(nnuc));
  arr_ref<float> dm(cmn.dm, dimension(nnuc));
  //
  */
  /*
  int i = 0;
  int j = 0;
  if (is_called_first_time) {
    {
      static const float values[] = {
        1.f, 1.f, 2.f, 3.f, 3.f, 4.f, 6.f, 7.f, 7.f, 8.f, 8.f, 9.f,
          10.f, 11.f, 11.f, 12.f, 12.f, 12.f, 13.f, 13.f, 14.f, 14.f,
          14.f, 15.f, 15.f, 16.f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        am;
    }
    {
      static const float values[] = {
        0.f, 1.f, 1.f, 1.f, 2.f, 2.f, 3.f, 3.f, 4.f, 3.f, 5.f, 4.f,
          5.f, 5.f, 6.f, 5.f, 6.f, 7.f, 6.f, 7.f, 6.f, 7.f, 8.f, 7.f,
          8.f, 8.f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        zm;
    }
    {
      static const float values[] = {
        .008665f, .007825f, .014102f, .016050f, .016030f, .002603f,
          .015125f, .016004f, .016929f, .022487f, .024609f, .012186f,
          .012939f, .009305f, .011432f, .014354f, .000000f, .018641f,
          .003354f, .005738f, .003242f, .003074f, .008597f, .000108f,
          .003070f, -.005085f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        dm;
    }
    {
      static const float values[] = {
        1.f, 1.f, 1.f, 0.f, 0.f, 2.f, 0.0f, 0.0f, 2.f, 1.f, 4.f, 0.f,
          0.f, 5.f, 0.0f, 0.0f, 3.f, 4.f, 10.f, 0.f, 0.f, 6.f, 0.0f,
          0.0f, 4.f, 1.f, 16.f, 0.f, 0.f, 17.f, 0.0f, 0.0f, 5.f, 1.f,
          21.f, 0.f, 0.f, 22.f, 0.0f, 0.0f, 6.f, 4.f, 11.f, 0.f, 0.f,
          6.f, 0.0f, 0.0f, 7.f, 1.f, 15.f, 0.f, 0.f, 14.f, 0.0f,
          0.0f, 8.f, 1.f, 18.f, 0.f, 0.f, 17.f, 0.0f, 0.0f, 9.f, 1.f,
          20.f, 0.f, 0.f, 19.f, 0.0f, 0.0f, 10.f, 1.f, 23.f, 0.f,
          0.f, 22.f, 0.0f, 0.0f, 11.f, 1.f, 25.f, 0.f, 0.f, 24.f,
          0.0f, 0.0f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 11) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        12.f, 2.f, 2.f, 1.f, 0.f, 3.f, 0.471f, 25.82f, 13.f, 2.f,
          3.f, 1.f, 0.f, 4.f, 1.63f, 72.62f, 14.f, 2.f, 5.f, 1.f,
          0.f, 6.f, 2.61f, 238.81f, 15.f, 2.f, 7.f, 1.f, 0.f, 8.f,
          1.19f, 84.17f, 16.f, 3.f, 5.f, 1.f, 2.f, 4.f, 1.002f,
          8.863f, 17.f, 3.f, 9.f, 1.f, 2.f, 8.f, 0.998f, 19.081f,
          18.f, 3.f, 7.f, 1.f, 4.f, 6.f, 1.070f, 55.494f, 19.f, 5.f,
          9.f, 1.f, 0.f, 6.f, 4.70f, 220.39f, 20.f, 2.f, 3.f, 2.f,
          0.f, 5.f, 1.63f, 63.750f, 21.f, 2.f, 4.f, 2.f, 0.f, 6.f,
          2.61f, 229.932f, 22.f, 2.f, 7.f, 2.f, 0.f, 9.f, 1.19f,
          65.054f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 12, 22) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        23.f, 3.f, 7.f, 2.f, 5.f, 6.f, 1.07f, 46.631f, 24.f, 5.f,
          8.f, 2.f, 0.f, 6.f, 4.69f, 201.291f, 25.f, 2.f, 6.f, 3.f,
          0.f, 7.f, 1.53f, 17.118f, 26.f, 2.f, 6.f, 4.f, 0.f, 8.f,
          1.11f, 28.640f, 27.f, 2.f, 6.f, 5.f, 0.f, 9.f, 1.11f,
          18.423f, 28.f, 6.f, 3.f, 0.f, 1.f, 5.f, 1.73f, 37.935f,
          29.f, 6.f, 3.f, 0.f, 2.f, 4.f, 1.73f, 46.798f, 30.f, 3.f,
          4.f, 3.f, 1.f, 6.f, 5.54f, 204.117f, 31.f, 3.f, 5.f, 3.f,
          2.f, 6.f, 5.55f, 212.980f, 32.f, 11.f, 5.f, 0.f, 2.f, 6.f,
          3.39f, 149.230f, 33.f, 9.f, 8.f, 3.f, 1.f, 6.f, 9.95f,
          175.476f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 23, 33) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        34.f, 9.f, 9.f, 3.f, 2.f, 6.f, 9.97f, 194.557f, 35.f, 2.f,
          8.f, 1.f, 0.f, 10.f, 1.31f, 23.59f, 36.f, 2.f, 13.f, 1.f,
          0.f, 14.f, 3.04f, 132.95f, 37.f, 2.f, 14.f, 1.f, 0.f, 16.f,
          2.34f, 39.10f, 38.f, 3.f, 15.f, 1.f, 2.f, 14.f, 1.002f,
          32.080f, 39.f, 3.f, 13.f, 1.f, 6.f, 8.f, 0.758f, 32.382f,
          40.f, 2.f, 9.f, 2.f, 0.f, 11.f, 1.30f, 1.595f, 41.f, 2.f,
          12.f, 2.f, 0.f, 13.f, 0.973f, 76.427f, 42.f, 2.f, 13.f,
          2.f, 0.f, 15.f, 3.03f, 100.840f, 43.f, 2.f, 14.f, 2.f, 0.f,
          17.f, 7.01f, 185.173f, 44.f, 2.f, 15.f, 2.f, 0.f, 18.f,
          2.33f, 6.975f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 34, 44) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        45.f, 3.f, 16.f, 2.f, 1.f, 17.f, 3.00f, 146.08f, 46.f, 3.f,
          12.f, 2.f, 6.f, 7.f, 0.618f, 24.674f, 47.f, 3.f, 13.f, 2.f,
          6.f, 9.f, 0.754f, 13.301f, 48.f, 3.f, 16.f, 2.f, 6.f, 12.f,
          0.292f, 79.89f, 49.f, 2.f, 7.f, 6.f, 0.f, 13.f, 1.58f,
          51.753f, 50.f, 2.f, 8.f, 6.f, 0.f, 14.f, 4.02f, 100.538f,
          51.f, 2.f, 9.f, 6.f, 0.f, 15.f, 4.02f, 87.539f, 52.f, 3.f,
          11.f, 6.f, 2.f, 15.f, 3.08f, 86.00f, 53.f, 3.f, 10.f, 6.f,
          1.f, 14.f, 3.07f, 76.96f, 54.f, 3.f, 12.f, 6.f, 1.f, 17.f,
          10.3f, 66.160f, 55.f, 3.f, 12.f, 3.f, 1.f, 13.f, 2.07f,
          50.63f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 45, 55) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        56.f, 3.f, 13.f, 3.f, 2.f, 14.f, 6.44f, 107.13f, 57.f, 3.f,
          14.f, 3.f, 1.f, 17.f, 14.9f, 159.36f, 58.f, 8.f, 6.f, 1.f,
          0.f, 12.f, 0.584f, 18.260f, 59.f, 7.f, 6.f, 0.f, 0.f, 17.f,
          2.00f, 84.420f, 60.f, 9.f, 10.f, 2.f, 1.f, 6.f, 3.58f,
          177.73f, 61.f, 9.f, 11.f, 1.f, 2.f, 6.f, 3.58f, 218.82f,
          62.f, 9.f, 12.f, 2.f, 3.f, 6.f, 0.807f, 7.555f, 63.f, 10.f,
          14.f, 2.f, 0.f, 6.f, 3.50f, 100.753f, 64.f, 10.f, 15.f,
          1.f, 0.f, 6.f, 3.49f, 132.83f, 65.f, 2.f, 17.f, 1.f, 0.f,
          19.f, 0.886f, 57.41f, 66.f, 2.f, 19.f, 1.f, 0.f, 21.f,
          3.58f, 94.88f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 56, 66) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        67.f, 2.f, 22.f, 1.f, 0.f, 24.f, 2.71f, 125.74f, 68.f, 3.f,
          20.f, 1.f, 2.f, 19.f, 1.002f, 34.846f, 69.f, 3.f, 22.f,
          1.f, 2.f, 21.f, 3.003f, 7.263f, 70.f, 3.f, 25.f, 1.f, 2.f,
          24.f, 1.002f, 41.037f, 71.f, 3.f, 25.f, 1.f, 6.f, 17.f,
          0.709f, 98.661f, 72.f, 2.f, 17.f, 2.f, 0.f, 20.f, 0.884f,
          22.553f, 73.f, 2.f, 19.f, 2.f, 0.f, 22.f, 1.19f, 87.621f,
          74.f, 2.f, 21.f, 2.f, 0.f, 24.f, 0.900f, 118.452f, 75.f,
          2.f, 20.f, 2.f, 0.f, 23.f, 3.57f, 53.706f, 76.f, 2.f, 22.f,
          2.f, 0.f, 25.f, 2.70f, 84.678f, 77.f, 2.f, 24.f, 2.f, 0.f,
          26.f, 3.62f, 140.734f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 67, 77) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        78.f, 3.f, 24.f, 2.f, 6.f, 17.f, 0.706f, 57.623f, 79.f, 2.f,
          17.f, 6.f, 0.f, 26.f, 5.13f, 83.111f, 80.f, 3.f, 13.f, 6.f,
          2.f, 19.f, 9.36f, 47.16f, 81.f, 3.f, 14.f, 6.f, 2.f, 21.f,
          11.0f, 9.098f, 82.f, 3.f, 15.f, 6.f, 2.f, 22.f, 3.68f,
          33.915f, 83.f, 3.f, 18.f, 6.f, 2.f, 25.f, 4.26f, 111.87f,
          84.f, 3.f, 20.f, 6.f, 2.f, 26.f, 5.81f, 60.557f, 85.f, 3.f,
          13.f, 6.f, 1.f, 20.f, 9.34f, 12.287f, 86.f, 3.f, 14.f, 6.f,
          1.f, 22.f, 3.67f, 1.835f, 87.f, 3.f, 16.f, 6.f, 1.f, 24.f,
          4.25f, 88.47f, 88.f, 3.f, 19.f, 6.f, 1.f, 26.f, 5.79f,
          25.711f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 78, 88) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
	*/

	/*
	am = {
		1.f, 1.f, 2.f, 3.f, 3.f, 4.f, 6.f, 7.f, 7.f, 8.f, 8.f, 9.f,
		10.f, 11.f, 11.f, 12.f, 12.f, 12.f, 13.f, 13.f, 14.f, 14.f,
		14.f, 15.f, 15.f, 16.f
	};
	zm = {
		0.f, 1.f, 1.f, 1.f, 2.f, 2.f, 3.f, 3.f, 4.f, 3.f, 5.f, 4.f,
		5.f, 5.f, 6.f, 5.f, 6.f, 7.f, 6.f, 7.f, 6.f, 7.f, 8.f, 7.f,
		8.f, 8.f
	};
	dm = {
		.008665f, .007825f, .014102f, .016050f, .016030f, .002603f,
		.015125f, .016004f, .016929f, .022487f, .024609f, .012186f,
		.012939f, .009305f, .011432f, .014354f, .000000f, .018641f,
		.003354f, .005738f, .003242f, .003074f, .008597f, .000108f,
		.003070f, -.005085f
	};
	reacpr = {
		{1,	1,	1,	0,	0,	2,	0,	0},
		{2,	1,	4,	0,	0,	5,	0,	0},
		{3,	4,	10,	0,	0,	6,	0,	0},
		{4,	1,	16,	0,	0,	17,	0,	0},
		{5,	1,	21,	0,	0,	22,	0,	0},
		{6,	4,	11,	0,	0,	6,	0,	0},
		{7,	1,	15,	0,	0,	14,	0,	0},
		{8,	1,	18,	0,	0,	17,	0,	0},
		{9,	1,	20,	0,	0,	19,	0,	0},
		{10,	1,	23,	0,	0,	22,	0,	0},
		{11,	1,	25,	0,	0,	24,	0,	0},
		{12,	2,	2,	1,	0,	3,	0.471,	25.82},
		{13,	2,	3,	1,	0,	4,	1.63,	72.62},
		{14,	2,	5,	1,	0,	6,	2.61,	238.81},
		{15,	2,	7,	1,	0,	8,	1.19,	84.17},
		{16,	3,	5,	1,	2,	4,	1.002,	8.863},
		{17,	3,	9,	1,	2,	8,	0.998,	19.081},
		{18,	3,	7,	1,	4,	6,	1.070,	55.494},
		{19,	5,	9,	1,	0,	6,	4.70,	220.39},
		{20,	2,	3,	2,	0,	5,	1.63,	63.750},
		{21,	2,	4,	2,	0,	6,	2.61,	229.932},
		{22,	2,	7,	2,	0,	9,	1.19,	65.054},
		{23,	3,	7,	2,	5,	6,	1.07,	46.631},
		{24,	5,	8,	2,	0,	6,	4.69,	201.291},
		{25,	2,	6,	3,	0,	7,	1.53,	17.118},
		{26,	2,	6,	4,	0,	8,	1.11,	28.640},
		{27,	2,	6,	5,	0,	9,	1.11,	18.423},
		{28,	6,	3,	0,	1,	5,	1.73,	37.935},
		{29,	6,	3,	0,	2,	4,	1.73,	46.798},
		{30,	3,	4,	3,	1,	6,	5.54,	204.117},
		{31,	3,	5,	3,	2,	6,	5.55,	212.980},
		{32,	11,	5,	0,	2,	6,	3.39,	149.230},
		{33,	9,	8,	3,	1,	6,	9.95,	175.476},
		{34,	9,	9,	3,	2,	6,	9.97,	194.557},
		{35,	2,	8,	1,	0,	10,	1.31,	23.59},
		{36,	2,	13,	1,	0,	14,	3.04,	132.95},
		{37,	2,	14,	1,	0,	16,	2.34,	39.10},
		{38,	3,	15,	1,	2,	14,	1.002,	32.080},
		{39,	3,	13,	1,	6,	8,	0.758,	32.382},
		{40,	2,	9,	2,	0,	11,	1.30,	1.595},
		{41,	2,	12,	2,	0,	13,	0.973,	76.427},
		{42,	2,	13,	2,	0,	15,	3.03,	100.840},
		{43,	2,	14,	2,	0,	17,	7.01,	185.173},
		{44,	2,	15,	2,	0,	18,	2.33,	6.975},
		{45,	3,	16,	2,	1,	17,	3.00,	146.08},
		{46,	3,	12,	2,	6,	7,	0.618,	24.674},
		{47,	3,	13,	2,	6,	9,	0.754,	13.301},
		{48,	3,	16,	2,	6,	12,	0.292,	79.89},
		{49,	2,	7,	6,	0,	13,	1.58,	51.753},
		{50,	2,	8,	6,	0,	14,	4.02,	100.538},
		{51,	2,	9,	6,	0,	15,	4.02,	87.539},
		{52,	3,	11,	6,	2,	15,	3.08,	86.00},
		{53,	3,	10,	6,	1,	14,	3.07,	76.96},
		{54,	3,	12,	6,	1,	17,	10.3,	66.160},
		{55,	3,	12,	3,	1,	13,	2.07,	50.63},
		{56,	3,	13,	3,	2,	14,	6.44,	107.13},
		{57,	3,	14,	3,	1,	17,	14.9,	159.36},
		{58,	8,	6,	1,	0,	12,	0.584,	18.260},
		{59,	7,	6,	0,	0,	17,	2.00,	84.420},
		{60,	9,	10,	2,	1,	6,	3.58,	177.73},
		{61,	9,	11,	1,	2,	6,	3.58,	218.82},
		{62,	9,	12,	2,	3,	6,	0.807,	7.555},
		{63,	10,	14,	2,	0,	6,	3.50,	100.753},
		{64,	10,	15,	1,	0,	6,	3.49,	132.83},
		{65,	2,	17,	1,	0,	19,	0.886,	57.41},
		{66,	2,	19,	1,	0,	21,	3.58,	94.88},
		{67,	2,	22,	1,	0,	24,	2.71,	125.74},
		{68,	3,	20,	1,	2,	19,	1.002,	34.846},
		{69,	3,	22,	1,	2,	21,	3.003,	7.263},
		{70,	3,	25,	1,	2,	24,	1.002,	41.037},
		{71,	3,	25,	1,	6,	17,	0.709,	98.661},
		{72,	2,	17,	2,	0,	20,	0.884,	22.553},
		{73,	2,	19,	2,	0,	22,	1.19,	87.621},
		{74,	2,	21,	2,	0,	24,	0.900,	118.452},
		{75,	2,	20,	2,	0,	23,	3.57,	53.706},
		{76,	2,	22,	2,	0,	25,	2.70,	84.678},
		{77,	2,	24,	2,	0,	26,	3.62,	140.734},
		{78,	3,	24,	2,	6,	17,	0.706,	57.623},
		{79,	2,	17,	6,	0,	26,	5.13,	83.111},
		{80,	3,	13,	6,	2,	19,	9.36,	47.16},
		{81,	3,	14,	6,	2,	21,	11.0,	9.098},
		{82,	3,	15,	6,	2,	22,	3.68,	33.915},
		{83,	3,	18,	6,	2,	25,	4.26,	111.87},
		{84,	3,	20,	6,	2,	26,	5.81,	60.557},
		{85,	3,	13,	6,	1,	20,	9.34,	12.287},
		{86,	3,	14,	6,	1,	22,	3.67,	1.835},
		{87,	3,	16,	6,	1,	24,	4.25,	88.47},
		{88,	3,	19,	6,	1,	26,	5.79,	25.711}
	};
	*/


		/*
	static const float reacpr[][] = {
		{ 	1, 	1, 	1, 	0, 	0, 	2, 	0, 	0 },
		{ 	2, 	1, 	4, 	0, 	0, 	5, 	0, 	0 },
		{ 	3, 	4, 	10,	0, 	0, 	6, 	0, 	0 },
		{ 	4, 	1, 	16,	0, 	0, 	17,	0, 	0 },
		{ 	5, 	1, 	21,	0, 	0, 	22, 0, 	0 },
		{ 	6, 	4, 	11, 0, 	0, 	6, 	0, 	0 },
		{ 	7, 	1, 	15, 0, 	0, 	14, 0, 	0 },
		{ 	8, 	1, 	18, 0, 	0, 	17, 0, 	0 },
		{ 	9, 	1, 	20, 0, 	0, 	19, 0, 	0 },
		{ 	10, 1, 	23, 0, 	0, 	22, 0, 	0 },
		{ 	11, 1, 	25, 0, 	0, 	24, 0, 	0 }
        { 	12, 2, 	2, 	1, 	0, 	3, 	0.471f, 25.82f },
		{ 	13,	2, 	3, 	1, 	0, 	4, 	1.63f, 	72.62f },
		{ 	14, 2, 	5, 	1, 	0, 	6, 	2.61f, 	238.81f },
		{ 	15, 2, 	7, 	1, 	0, 	8, 	1.19f, 	84.17f },
		{ 	16, 3, 	5, 	1, 	2, 	4, 	1.002f, 8.863f },
		{ 	17, 3, 	9, 	1, 	2, 	8, 	0.998f, 19.081f },
		{ 	18, 3, 	7, 	1, 	4, 	6, 	1.070f, 55.494f },
		{ 	19 	5, 	9, 	1, 	0, 	6, 	4.70f, 	220.39f },
		{ 	20, 2, 	3, 	2, 	0, 	5, 	1.63f, 	63.750f },
		{ 	21, 2, 	4, 	2, 	0, 	6, 	2.61f, 	229.932f},
		{	22, 2, 	7, 	2, 	0, 	9, 	1.19f, 65.054f };
      };
    {
      static const float values[] = {
        23, 3, 7, 2, 5, 6, 1.07f, 46.631f, 24, 5,
          8, 2, 0, 6, 4.69f, 201.291f, 25, 2, 6, 3,
          0, 7, 1.53f, 17.118f, 26, 2, 6, 4, 0, 8,
          1.11f, 28.640f, 27, 2, 6, 5, 0, 9, 1.11f,
          18.423f, 28, 6, 3, 0, 1, 5, 1.73f, 37.935f,
          29, 6, 3, 0, 2, 4, 1.73f, 46.798f, 30, 3,
          4, 3, 1, 6, 5.54f, 204.117f, 31, 3, 5, 3,
          2, 6, 5.55f, 212.980f, 32, 11, 5, 0, 2, 6,
          3.39f, 149.230f, 33, 9, 8, 3, 1, 6, 9.95f,
          175.476f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 23, 33) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        34, 9, 9, 3, 2, 6, 9.97f, 194.557f, 35, 2,
          8, 1, 0, 10, 1.31f, 23.59f, 36, 2, 13, 1,
          0, 14, 3.04f, 132.95f, 37, 2, 14, 1, 0, 16,
          2.34f, 39.10f, 38, 3, 15, 1, 2, 14, 1.002f,
          32.080f, 39, 3, 13, 1, 6, 8, 0.758f, 32.382f,
          40, 2, 9, 2, 0, 11, 1.30f, 1.595f, 41, 2,
          12, 2, 0, 13, 0.973f, 76.427f, 42, 2, 13,
          2, 0, 15, 3.03f, 100.840f, 43, 2, 14, 2, 0,
          17, 7.01f, 185.173f, 44, 2, 15, 2, 0, 18,
          2.33f, 6.975f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 34, 44) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        45, 3, 16, 2, 1, 17, 3.00f, 146.08f, 46, 3,
          12, 2, 6, 7, 0.618f, 24.674f, 47, 3, 13, 2,
          6, 9, 0.754f, 13.301f, 48, 3, 16, 2, 6, 12,
          0.292f, 79.89f, 49, 2, 7, 6, 0, 13, 1.58f,
          51.753f, 50, 2, 8, 6, 0, 14, 4.02f, 100.538f,
          51, 2, 9, 6, 0, 15, 4.02f, 87.539f, 52, 3,
          11, 6, 2, 15, 3.08f, 86.00f, 53, 3, 10, 6,
          1, 14, 3.07f, 76.96f, 54, 3, 12, 6, 1, 17,
          10.3f, 66.160f, 55, 3, 12, 3, 1, 13, 2.07f,
          50.63f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 45, 55) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        56, 3, 13, 3, 2, 14, 6.44f, 107.13f, 57, 3,
          14, 3, 1, 17, 14.9f, 159.36f, 58, 8, 6, 1,
          0, 12, 0.584f, 18.260f, 59, 7, 6, 0, 0, 17,
          2.00f, 84.420f, 60, 9, 10, 2, 1, 6, 3.58f,
          177.73f, 61, 9, 11, 1, 2, 6, 3.58f, 218.82f,
          62, 9, 12, 2, 3, 6, 0.807f, 7.555f, 63, 10,
          14, 2, 0, 6, 3.50f, 100.753f, 64, 10, 15,
          1, 0, 6, 3.49f, 132.83f, 65, 2, 17, 1, 0,
          19, 0.886f, 57.41f, 66, 2, 19, 1, 0, 21,
          3.58f, 94.88f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 56, 66) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        67, 2, 22, 1, 0, 24, 2.71f, 125.74f, 68, 3,
          20, 1, 2, 19, 1.002f, 34.846f, 69, 3, 22,
          1, 2, 21, 3.003f, 7.263f, 70, 3, 25, 1, 2,
          24, 1.002f, 41.037f, 71, 3, 25, 1, 6, 17,
          0.709f, 98.661f, 72, 2, 17, 2, 0, 20, 0.884f,
          22.553f, 73, 2, 19, 2, 0, 22, 1.19f, 87.621f,
          74, 2, 21, 2, 0, 24, 0.900f, 118.452f, 75,
          2, 20, 2, 0, 23, 3.57f, 53.706f, 76, 2, 22,
          2, 0, 25, 2.70f, 84.678f, 77, 2, 24, 2, 0,
          26, 3.62f, 140.734f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 67, 77) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
    {
      static const float values[] = {
        78, 3, 24, 2, 6, 17, 0.706f, 57.623f, 79, 2,
          17, 6, 0, 26, 5.13f, 83.111f, 80, 3, 13, 6,
          2, 19, 9.36f, 47.16f, 81, 3, 14, 6, 2, 21,
          11.0f, 9.098f, 82, 3, 15, 6, 2, 22, 3.68f,
          33.915f, 83, 3, 18, 6, 2, 25, 4.26f, 111.87f,
          84, 3, 20, 6, 2, 26, 5.81f, 60.557f, 85, 3,
          13, 6, 1, 20, 9.34f, 12.287f, 86, 3, 14, 6,
          1, 22, 3.67f, 1.835f, 87, 3, 16, 6, 1, 24,
          4.25f, 88.47f, 88, 3, 19, 6, 1, 26, 5.79f,
          25.711f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 78, 88) {
        FEM_DO_SAFE(j, 1, 8) {
          data, reacpr(i, j);
        }
      }
    }
	*/
    const float cy0 = .300f;
    const float ct0 = .030f;
    const float t9i0 = 1.00e+02f;
    const float t9f0 = 1.00e-02f;
    const float ytmin0 = 1.00e-25f;
    const float inc0 = 30;
    const float c0[] = { 1.00f, 885.7f, 3.0f };
	/*
    {
      static const float values[] = {
        1.00f, 885.7f, 3.0f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        cmn.c0;
    }
	*/
    const float cosmo0 = 0.00f;
	const float xi0[] = {0,0,0};
	/*
    {
      static const float values[] = {
        0.00f, 0.00f, 0.00f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        cmn.xi0;
    }
	*/
    dt0 = 1.00e-04f;
    eta0 = 3.162e-10f;
  // } if called for first time
  //
  //----------PARAMETERS.
  //Number of nuclear reactions.
  //Number of nuclides in calculation.
  //
  //----------COMMON AREAS.
  //Reaction parameter
  //Default comp param
  //Default model para
  //Default variationl
  //Nuclide data.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION PARAMETERS VALUES.
  //Reaction parameters.
  //
  //----------DEFAULT COMPUTATION PARAMETERS.
  //Default time step limiting constant.
  //Default initial temperature (in 10**
  //Default final temperature (in 10**9
  //Default smallest abundances allowed.
  //Default accumulation increment.
  //
  //----------DEFAULT MODEL PARAMETERS.
  //c0(1) is default variation of grav c
  //c0(2) is default neutron half-life.
  //c0(3) is default number of neutrinos
  //Default cosmological constant.
  //Default neutrino degeneracy paramete
  //
  //----------DEFAULT VARIATIONAL PARAMETERS.
  //Default initial time step.
  //Default baryon-to-photon ratio.
  //
  //----------NUCLIDE DATA.
  //Atomic number of nuclide.
  //Charge of nuclide.
  //Mass excess of nuclide.
  //
  //==============================DATA DIVISION====================================
  //
  //    Nuclide and corresponding number
  //    --------------------------------
  //    1) N         7) Li6      13) B10      19) C13      25) O15
  //    2) P         8) Li7      14) B11      20) N13      26) O16
  //    3) H2        9) Be7      15) C11      21) C14
  //    4) H3       10) Li8      16) B12      22) N14
  //    5) He3      11) B8       17) C12      23) O14
  //    6) He4      12) Be9      18) N12      24) N15
  //
  //-----------NUCLIDE DATA.
  //
  //----------REACTION RATE COEFFICIENTS (Ref 1).
  //              reac# type n1 n2 n3 n4 rev-coeff q-value
  //              ----  ---- -- -- -- -- --------- -------
  //N->P
  //H3->He3
  //Li8->2He4
  //B12->C12
  //C14->N14
  //B8->2He4
  //C11->B11
  //N12->C12
  //N13->C13
  //O14->N14
  //O15->N15
  //              reac# type n1 n2 n3 n4 rev-coeff q-value
  //              ----  ---- -- -- -- -- --------- -------
  //H(n,g)H2
  //H2(n,g)H3
  //He3(n,g)He4
  //Li6(n,g)Li7
  //He3(n,p)H3
  //Be7(n,p)Li7
  //Li6(n,a)H3
  //Be7(n,a)He4
  //H2(p,g)He3
  //H3(p,g)He4
  //Li6(p,g)Be7
  //              reac# type n1 n2 n3 n4 rev-coeff q-value
  //              ----  ---- -- -- -- -- --------- -------
  //Li6(p,a)He3
  //Li7(p,a)He4
  //H2(a,p)Li6
  //H3(a,p)Li7
  //He3(a,p)Be7
  //H2(d,p)He3
  //H2(d,n)H3
  //H3(d,n)He4
  //He3(d,p)He4
  //He3(He3,2p)H
  //Li7(d,na)He4
  //              reac# type n1 n2 n3 n4 rev-coeff q-value
  //              ----  ---- -- -- -- -- --------- -------
  //Be7(d,pa)He4
  //Li7(n,g)Li8
  //B10(n,g)B11
  //B11(n,g)B12
  //C11(n,p)B11
  //B10(n,a)Li7
  //Be7(p,g)B8
  //Be9(p,g)B10
  //B10(p,g)C11
  //B11(p,g)C12
  //C11(p,g)N12
  //              reac# type n1 n2 n3 n4 rev-coeff q-value
  //              ----  ---- -- -- -- -- --------- -------
  //B12(p,n)C12
  //Be9(p,a)Li6
  //B10(p,a)Be7
  //B12(p,a)Be9
  //Li6(a,g)B10
  //Li7(a,g)B11
  //Be7(a,g)C11
  //B8(a,p)C11
  //Li8(a,n)B11
  //Be9(a,n)C12
  //Be9(d,n)B10
  //              reac# type n1 n2 n3 n4 rev-coeff q-value
  //              ----  ---- -- -- -- -- --------- -------
  //B10(d,p)B11
  //B11(d,n)C12
  //He4(an,g)Be9
  //He4(2a,g)C12
  //Li8(p,na)He4
  //B8(n,pa)He4
  //Be9(p,da)He4
  //B11(p,2a)Be4
  //C11(n,2a)He4
  //C12(n,g)C13
  //C13(n,g)C14
  //              reac# type n1 n2 n3 n4 rev-coeff q-value
  //              ----  ---- -- -- -- -- --------- -------
  //N14(n,g)N15
  //N13(n,p)C13
  //N14(n,p)C14
  //O15(n,p)N15
  //O15(n,a)C12
  //C12(p,g)N13
  //C13(p,g)N14
  //C14(p,g)N15
  //N13(p,g)O14
  //N14(p,g)O15
  //N15(p,g)O16
  //              reac# type n1 n2 n3 n4 rev-coeff q-value
  //              ----  ---- -- -- -- -- --------- -------
  //N15(p,a)C12
  //C12(a,g)O16
  //B10(a,p)C13
  //B11(a,p)C14
  //C11(a,p)N14
  //N12(a,p)O15
  //N13(a,p)O16
  //B10(a,n)N13
  //B11(a,n)N14
  //B12(a,n)N15
  //C13(a,n)O16
  //
  //----------DEFAULT COMPUTATION PARAMETERS.
  //Default time step limiting constant.
  //Default initial temperature.
  //Default final temperature.
  //Default smallest abundances allowed.
  //Default accumulation increment.
  //-----------DEFAULT MODEL PARAMETERS.
  //Default variation of 3 parameters.
  //Default cosmological constant.
  //Default neutrino degeneracy paramete
  //
  //-----------DEFAULT VARIATIONAL PARAMETERS.
  //Default initial time step.
  //Default baryon-to-photon ratio.
  //
}

/*
struct program_new123_save
{
  fem::variant_bindings rates_bindings;
};
*/

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//     Changes (to run inder unix f77):
//     -------------------------------
// Programme renamed from nuc123.for to new123.f
// ir=1 -> ir=5 ... input unit number
// iw=1 -> iw=6 ... output unit number
// COMMON /check/ itime -> COMMON /checkcb/ itime
// COMMON /time/        -> COMMON /ttime/
// output nuc123.dat    -> new123.dat
//
//========================IDENTIFICATION DIVISION================================
//
void
program_new123(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  //blockdata_unnamed(cmn);
  //FEM_CMN_SVE(program_new123);
  common_read read(cmn);
  common_write write(cmn);
  /*
  const int nrec = 88;
  arr_cref<float, 2> reacpr(cmn.reacpr, dimension(nrec, 8));
  arr_ref<int> iform(cmn.iform, dimension(nrec));
  arr_ref<int> ii(cmn.ii, dimension(nrec));
  arr_ref<int> jj(cmn.jj, dimension(nrec));
  arr_ref<int> kk(cmn.kk, dimension(nrec));
  arr_ref<int> ll(cmn.ll, dimension(nrec));
  arr_ref<float> rev(cmn.rev, dimension(nrec));
  arr_ref<float> q9(cmn.q9, dimension(nrec));
  arr_cref<float> c0(cmn.c0, dimension(3));
  arr_cref<float> xi0(cmn.xi0, dimension(3));
  arr_ref<float> c(cmn.c, dimension(3));
  arr_ref<float> xi(cmn.xi, dimension(3));
  //int& itime = static_cast<common_checkcb&>(cmn).itime;
  int& itime = cmn.itime;
  bool& outfile = cmn.outfile;
  */
  //
  //common_variant rates(cmn.common_rates, sve.rates_bindings);
  /*
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> f(dimension(nrec));
      mbr<float> r(dimension(nrec));
      rates.allocate(), f, r;
    }
  }
  arr_ref<float> f(rates.bind<float>(), dimension(nrec));
  arr_ref<float> r(rates.bind<float>(), dimension(nrec));
  arr<float> f(dimension(nrec));
  arr<float> r(dimension(nrec));
  const int nnuc = 26;
  */ 
  const int iw = 6;
  int i = 0;
  const int ir = 5;
  int inum = 0;
  //
  //----------LINKAGES.
  //     CALLED BY - none
  //     CALLS     - [subroutine] help, setcom, setmod, run, output
  //
  //----------REMARKS.
  //     Control program -
  //       Offers user the main menu and channels through to various options.
  //     Implementation -
  //       To run this program, new123.f must be linked with nuccom.f
  //       (containing the computation subroutines), nucrat.f (with the
  //       reaction rates), and newint.f (with an interface subroutine).
  //       This program has been written to be compatible with
  //       ANSI FORTRAN-77 with the exception of the END DO statement
  //       used to limit the number of statement labels.
  //       This code was modified on the DEC/AXP system.
  //     Notes -
  //       The program utilizes Wagoner's code as the core of the computational
  //       routines.
  //     Documentation -
  //       Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
  //       Kellogg Radiation Lab preprint OAP-714.
  //     Copy -
  //       Version 4.1 (December 1991)
  //
  //----------PARAMETERS.
  //Input unit number (previous value =
  //Output unit number (previous value =
  //Number of nuclear reactions.
  //Number of nuclides in calculation.
  //
  //----------COMMON AREAS.
  //Reaction parameter
  //Reaction rates.
  //Default comp param
  //Computation parame
  //Default model para
  //Model parameters.
  //Default variationl
  //Variational parame
  //Computation locati
  //Run options.
  //Output option.
  //
  //==========================DECLARATION DIVISION=================================
  //
  //----------REACTION PARAMETERS FROM BLOCK DATA.
  //Reaction parameters.
  //
  //----------REACTION PARAMETERS.
  //Reaction type code (1-11).
  //Incoming nuclide type (1-26).
  //Incoming light nuclide type (1-6).
  //Outgoing light nuclide type (1-6).
  //Outgoing nuclide type (1-26).
  //Reverse reaction coefficient.
  //Energy released in reaction.
  //
  //----------REACTION RATES.
  //Forward reaction rate coefficients.
  //Reverse reaction rate coefficients.
  //
  //----------DEFAULT COMPUTATION PARAMETERS.
  //Default cy.
  //Default ct.
  //Default t9i.
  //Default t9f.
  //Default ytmin.
  //Default accumulation increment.
  //
  //----------COMPUTATIONAL PARAMETERS.
  //Time step limiting constant on abund
  //Time step limiting constant on tempe
  //Initial temperature (in 10**9 K).
  //Final temperature (in 10**9 k).
  //Smallest abundances allowed.
  //Accumulation increment.
  //
  //----------DEFAULT MODEL PARAMETERS.
  //Default c.
  //Default cosmological constant.
  //Default neutrino degeneracy paramete
  //
  //----------EARLY UNIVERSE MODEL PARAMETERS.
  //c(1) is variation of gravitational c
  //c(2) is neutron lifetime (sec).
  //c(3) is number of neutrino species.
  //Cosmological constant.
  //Neutrino degeneracy parameters.
  //
  //----------DEFAULT VARIATIONAL PARAMETERS.
  //Default initial time step.
  //Default baryon-to-photon ratio.
  //
  //----------VARIATIONAL PARAMETERS.
  //Initial time step.
  //Baryon-to-photon ratio.
  //
  //----------COMPUTATION LOCATION.
  //Time check.
  //
  //----------RUN OPTION.
  //Run network size.
  //Number of nuclides in computation.
  //Number of reactions in computation.
  //
  //----------OUTPUT FILE STATUS.
  //Number of output requests.
  //Indicates if output file used.
  //
  //----------USER RESPONSE VARIABLES.
  //Selection number.
  //
  //===========================PROCEDURE DIVISION==================================
  //
  //10--------OPEN FILES AND PRINT GREETING----------------------------------------
  //
  //Output file.
  cmn.io.open(2, "new123.dat")
    .status("new");
  //Time = beginning of program.
  cmn.itime = 1;
  //Check interface subroutine.
  cmn.check(cmn);
  write(iw,
    "(6(/),2(' ',4x,'NN',6x,'NN  UU',6x,'UU',4x,8('C'),6x,'11',8x,6('2'),6x,"
    "6('3'),/),2(' ',4x,'NN',6x,'NN  UU',6x,'UU  CC',12x,'1111',6x,'22',6x,"
    "'22  33',6x,'33',/),2(' ',4x,'NNNN    NN  UU',6x,'UU  CC',14x,'11',14x,"
    "'22',10x,'33',/),2(' ',4x,'NN  NN  NN  UU',6x,'UU  CC',14x,'11',12x,'22',"
    "10x,'33',/),2(' ',4x,'NN    NNNN  UU',6x,'UU  CC',14x,'11',10x,'22',14x,"
    "'33',/),2(' ',4x,'NN',6x,'NN  UU',6x,'UU  CC',14x,'11',8x,'22',8x,'33',"
    "6x,'33',/),2(' ',4x,'NN',6x,'NN  ',10('U'),4x,8('C'),4x,6('1'),4x,"
    "10('2'),4x,6('3'),/),/,26x,'MODIFIED APRIL 1994',/,/,/,"
    "'(Press <RETURN> to continue): ',$)");
  //
  //20--------INPUT INITIALIZATION INFORMATION AND PAUSE---------------------------
  //
  FEM_DO_SAFE(i, 1, cmn.nrec) {
    //..........READ IN REACTION PARAMETERS.
	/*
    cmn.iform(i) = cmn.reacpr(i, 2); 		/// Reaction type.
    cmn.ii(i) = cmn.reacpr(i, 3); 			/// Incoming nuclide type.
    cmn.jj(i) = cmn.reacpr(i, 4); 			/// Incoming nuclide type.
    cmn.kk(i) = cmn.reacpr(i, 5); 			/// Outgoing nuclide type.
    cmn.ll(i) = cmn.reacpr(i, 6); 			/// Outgoing nuclide type.
    cmn.rev(i) = cmn.reacpr(i, 7); 			/// Reverse reaction coefficient.
    cmn.q9(i) = cmn.reacpr(i, 8); 			/// Energy released.
	*/
    cmn.iform(i) = cmn.reacpr[i-1][2-1]; 		/// Reaction type.
    cmn.ii(i) = cmn.reacpr[i-1][3-1]; 			/// Incoming nuclide type.
    cmn.jj(i) = cmn.reacpr[i-1][4-1]; 			/// Incoming nuclide type.
    cmn.kk(i) = cmn.reacpr[i-1][5-1]; 			/// Outgoing nuclide type.
    cmn.ll(i) = cmn.reacpr[i-1][6-1]; 			/// Outgoing nuclide type.
    cmn.rev(i) = cmn.reacpr[i-1][7-1]; 			/// Reverse reaction coefficient.
    cmn.q9(i) = cmn.reacpr[i-1][8-1]; 			/// Energy released.
    //..........INITIALIZE REACTION RATES.
    cmn.f(i) = 0; 							/// Forward rate coeff.
    cmn.r(i) = 0; 							/// Reverse rate coeff.
    //..........SET RUN OPTIONS TO DEFAULT.
  }
  cmn.irun = 1; 							/// Do full run.
  cmn.isize = cmn.nnuc; 					/// Use all 26 nuclides.
  cmn.jsize = cmn.nrec; 					/// Use all 88 reactions.
  //..........SET OUTPUT OPTION TO DEFAULT.
  cmn.nout = 0; 							/// No output requests.
  cmn.outfile = false; 						/// Output file not used.
  //..........SET VALUES TO DEFAULT.
  cmn.cy = cmn.cy0; 						/// Time step limiting constant on abundance.
  cmn.ct = cmn.ct0; 						/// Time step limiting constant on temperature.
  cmn.t9i = cmn.t9i0; 						/// Initial temperature.
  cmn.t9f = cmn.t9f0; 						/// Final temperature.
  cmn.ytmin = cmn.ytmin0; 					/// Smallest abundances allowed.
  cmn.inc = cmn.inc0; 						/// Accumulation increment.
  cmn.c(1) = cmn.c0(1); 					/// Variation of gravitational constant.
  cmn.c(2) = cmn.c0(2); 					/// Neutron lifetime.
  cmn.c(3) = cmn.c0(3); 					/// Number of neutrino species.
  cmn.cosmo = cmn.cosmo0; 					/// Cosmological constant.
  cmn.xi(1) = cmn.xi0(1); 					/// Electron degeneracy parameter.
  cmn.xi(2) = cmn.xi0(2); 					/// Muon degeneracy parameter.
  cmn.xi(3) = cmn.xi0(3); 					/// Tau degeneracy parameter.
  cmn.dt1 = cmn.dt0; 						/// Initial time step.
  cmn.eta1 = cmn.eta0; 						/// Baryon-to-photon ratio.
  //..........ACCEPT RETURN TO CONTINUE.
  read(ir, star); 							/// Pause.
  //
  //30--------PRINT MENU AND AWAIT RESPONSE----------------------------------------
  //
  //..........RETURN FROM LOOPING.
  statement_300:
  //..........DISPLAY MENU.
  write(iw,
    "(8(/),32x,"
	"'MENU SELECTION',/,32x,"
	"'---- ---------',/,/,24x,"
    "'1. HELP',/,24x,"
	"'2. SET COMPUTATION PARAMETERS',/,24x,"
    "'3. SET MODEL PARAMETERS',/,24x,"
	"'4. RUN',/,24x,"
	"'5. OUTPUT',/,24x,"
	"'6. EXIT',8(/),24x,"
	"'Enter selection (1-6): ',$)");
  //..........READ IN SELECTION NUMBER.
  read(ir, "(i1)"), inum;
  //
  //40--------BRANCH TO APPROPRIATE SECTION----------------------------------------
  //
  switch (inum) {
    case 1: goto statement_410;
    case 2: goto statement_420;
    case 3: goto statement_430;
    case 4: goto statement_440;
    case 5: goto statement_450;
    case 6: goto statement_460;
    default: break;
  }
  //Improper input or <RETURN>.
  goto statement_460;
  //Help section.
  statement_410:
  help(cmn);
  goto statement_500;
  //Set computation parameters section.
  statement_420:
  cmn.setcom(cmn);
  goto statement_500;
  //Set model parameters section.
  statement_430:
  cmn.setmod(cmn);
  goto statement_500;
  //Run section.
  statement_440:
  //Time = beginning of run section.
  cmn.itime = 2;
  //Check interface subroutine.
  cmn.check(cmn);
  cmn.run(cmn);
  //Time = end of run section.
  cmn.itime = 9;
  //Check interface subroutine.
  cmn.check(cmn);
  goto statement_500;
  //Output section.
  statement_450:
  output(cmn);
  goto statement_500;
  //Exit section.
  statement_460:
  if (cmn.outfile) {
    //Close output file.
    cmn.io.close(2).status("keep");
  }
  else {
    //File not used - dispose.
    cmn.io.close(2).status("delete");
  }
  //Time = end of program.
  cmn.itime = 10;
  //Check interface subroutine.
  cmn.check(cmn);
  //FEM_STOP(0);
  //
  //50---------GO BACK TO MENU-----------------------------------------------------
  //
  statement_500:
  goto statement_300;
  //
}

} // namespace bbn

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    bbn::program_new123);
}

