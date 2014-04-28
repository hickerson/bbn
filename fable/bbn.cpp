#define FEM_TRANSLATION_UNIT_WITH_MAIN

#include "BigBangModel.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>

// namespace bbn {
using namespace bbn;
using namespace boost::math;

//using namespace fem;
//using namespace std;

//
//========================IDENTIFICATION DIVISION================================
//
void common::help()
{
	common_read read(*this);
	common_write write(*this);
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
	//..........RETURN FROM LOOP.
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
			"'   to find the derivatives of T9, hv, and phie (this is done in subroutine      ',/,"
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
void common::setcom()
{
	common_read read(*this);
	common_write write(*this);
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
	//Default T9i.
	//Default T9f.
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
	//..........RETURN FROM LOOP.
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
		cy, ct, dt1, T9i, T9f, ytmin, fem::ffloat(inc);
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
	read(ir, star), T9i;
	goto statement_400;
	//Change final temperature section.
statement_250:
	write(iw, "(' ','Enter value for final temperature: ',$)");
	read(ir, star), T9f;
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
	cy = cy0;
	//Time step limiting constant on temperature
	ct = ct0;
	//Time step.
	dt1 = dt0;
	//Initial temperature.
	T9i = T9i0;
	//Final temperature.
	T9f = T9f0;
	//Smallest abundances allowed.
	ytmin = ytmin0;
	//Accumulation increment.
	inc = inc0;
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
void common::setmod()
{
	common_read read(*this);
	common_write write(*this);
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
	//c[1] is variation of gravitational con
	//c[2] is neutron lifetime (sec).
	//c[3] is number of neutrino species.
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
	//..........RETURN FROM LOOP.
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
		c[1], c[2], c[3], eta1, cosmo, xi[1], xi[2], xi[3];
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
	read(ir, star), c[1];
	goto statement_400;
	//Change neutron lifetime section.
statement_220:
	write(iw, "(' ','Enter value for neutron lifetime (sec): ',$)");
	read(ir, star), c[2];
	goto statement_400;
	//Change number of neutrino species section.
statement_230:
	write(iw, "(' ','Enter value for number of neutrino species: ',$)");
	read(ir, star), c[3];
	goto statement_400;
	//Change baryon-to-photon ratio section.
statement_240:
	write(iw, "(' ','Enter value for baryon-to-photon ratio: ',$)");
	read(ir, star), eta1;
	goto statement_400;
	//Change cosmological constant section.
statement_250:
	write(iw, "(' ','Enter value for cosmological constant: ',$)");
	read(ir, star), cosmo;
	goto statement_400;
	//Change neutrino degeneracy section.
statement_260:
	write(iw, "(' ','Enter value for xi electron: ',$)");
	read(ir, star), xi[1];
	goto statement_400;
	//Change neutrino degeneracy section.
statement_270:
	write(iw, "(' ','Enter value for xi muon: ',$)");
	read(ir, star), xi[2];
	goto statement_400;
	//Change neutrino degeneracy section.
statement_280:
	write(iw, "(' ','Enter value for xi tauon: ',$)");
	read(ir, star), xi[3];
	if ((xi[3] != 0.f) && (c[3] < 3.f)) {
		c[3] = 3.f;
		write(iw, "(' ','Number of neutrinos set to 3')");
		write(iw, "(' ','Press <RETURN> to continue: ',$)");
		read(ir, star);
	}
	goto statement_400;
	//Reset all to default values section.
statement_290:
	c[1] = c0[1];
	c[2] = c0[2];
	c[3] = c0[3];
	cosmo = cosmo0;
	xi[1] = xi0[1];
	xi[2] = xi0[2];
	xi[3] = xi0[3];
	eta1 = eta0;
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
void common::check()
{
	common_write write(*this);
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
	//Default T9i.
	//Default T9f.
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
	//c[1] is variation of gravitational constant.
	//c[2] is neutron half-life (min).
	//c[3] is number of neutrino species.
	//Cosmological constant.
	//Neutrino degeneracy parameters.
	//xi[1] is e neutrino degeneracy parameter.
	//xi[2] is m neutrino degeneracy parameter.
	//xi[3] is t neutrino degeneracy parameter.
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
	//(1/T9)*d(T9)/d(t).
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
		io.open(3, "newint.dat").status("new");
	}
	//
	//20--------WRITE INTO FILE------------------------------------------------------
	//
	if (itime == 8) { 							/// Right after a run.
		xout(it,8) += xout(it,9); 				/// Add beryllium to lithium.
		xout(it,5) += xout(it,4); 				/// Add tritium to helium-3.
		xout(it,6) -= 0.0003f; 					/// my correction for fitted rates+coarse steps
		write(3, "(7(e13.5))"), c[3], c[2], etaout(it), xout(it, 3),
			xout(it,5), xout(it,6), xout(it,8);	/// Output N_nu, tau_n, eta, H2, He3, He4, an Li7.
	}
	//
	//30--------CLOSE FILE-----------------------------------------------------------
	//
	//End of program.
	if (itime == 10) {
		io.close(2);
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
double ex( double const& x)
{
	double return_value = 0;
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


double common::getBesselL(double r)
{
	double K2r = cyl_bessel_k(2,r); 			/// Irregular modified cylindrical Bessel functions.
	return K2r/r;
}

double common::getBesselM(double r)
{
	double K3r = cyl_bessel_k(3,r);			/// Irregular modified cylindrical Bessel functions.
	double K1r = cyl_bessel_k(1,r);			/// Irregular modified cylindrical Bessel functions.
	return (3*K3r + K1r)/4/r;				/// (Ref ?).
}

double common::getBesselN(double r)
{
	double K4r = cyl_bessel_k(4,r);			/// Irregular modified cylindrical Bessel functions.
	double K2r = cyl_bessel_k(2,r);			/// Irregular modified cylindrical Bessel functions.
	return (K4r + K2r)/2/r; 				/// (Ref ?).
}


//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//     Changes (to run inder DEC unix f77):
//     -----------------------------------
//     COMMON /therm/ -> COMMON /thermcb/
//     COMMON /rates/ f,r(nrec) -> COMMON /rates/ f,r
//
//     Default neutron lifetime 888.54 -> 885.7 (+- 0.8) - PDG 2002/3
//
void common::rate0()
{
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
	f[2] = 1.79e-9f;
	//
	//.......Li8 -> e- + v + 2He4.......(Ajzenberg-Selove 1988)
	f[3] = 8.27e-1f;
	//
	//.......B12 -> e- + B + C12........(Ajzenberg-Selove 1990)
	f[4] = 3.43e+1f;
	//
	//.......C14 -> e- + v + N14........(Ajzenberg-Selove 1986)
	f[5] = 3.834e-12f;
	//
	//.......B8 -> e+ + v + 2He4........(Ajzenberg-Selove 1988)
	f[6] = 9.00e-1f;
	//
	//.......C11 -> e+ + v + B11........(Ajzenberg-Selove 1990)
	f[7] = 5.668e-4f;
	//
	//.......N12 -> e+ + v + C12........(Ajzenberg-Selove 1990)
	f[8] = 6.301e+1f;
	//
	//.......N13 -> e+ + v + C13........(Ajzenberg-Selove 1986)
	f[9] = 1.159e-3f;
	//
	//.......O14 -> e+ + v + N14........(Ajzenberg-Selove 1986)
	f[10] = 9.8171e-3f;
	//
	//.......O15 -> e+ + v + N15........(Ajzenberg-Selove 1986)
	f[11] = 5.6704e-3f;
	//
	//----------REFERENCES-----------------------------------------------------------
	//     Ajzenberg-Selove, F., 1990, Nucl. Phys. A506, 1.
	//     Ajzenberg-Selove, F., 1988, Nucl. Phys. A490, 1.
	//     Ajzenberg-Selove, F., 1986, Nucl. Phys. A449, 1.
	//     Tilley, D.R., Weller, H.R., and Hasan, H.H., 1987, Nucl. Phys. A474, 1.
	//
}


//
//========================IDENTIFICATION DIVISION================================
//
//===========================PROCEDURE DIVISION==================================
//
//10--------1ST PART OF INTEGRAL FOR n->p RATE-----------------------------------
//***************************************************************************
//
	double
common::func1(
		double const& x)
{
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
	double part1 = 0; //TODO move down inside loop
	double part2 = 0;
	if (x <= 0)
		return 0;
	else {
		// TODO don't recompute 
		part1 = 1 / (1 + ex(-.511f * x / T9mev));
		part2 = 1 / (1 + ex(+(x - 2.531f) * (.511f / tnmev) - xi[1]));
		return cnorm * x * fem::pow2((x - 2.531f)) * pow(
				(fem::pow2(x) - 1), 0.5) * part1 * part2;
	}
}


//
//20--------2ND PART OF INTEGRAL FOR n->p RATE-----------------------------------
//
	double
common::func2(
		double const& x)
{
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
	if (x <= 1.f) {
		return 0;
	}
	else {		// TODO remove cmn
		double part1 = 1.f / (1.f + ex(+.511f * x / T9mev));
		double part2 = 1.f / (1.f + ex(-(x + 2.531f) * (.511f / tnmev) - xi[1]));
		return cnorm * x * fem::pow2((x + 2.531f)) 
				* pow((fem::pow2(x) - 1), .5f) * part1 * part2;
	}
}


//
//30--------1ST PART OF INTEGRAL FOR p->n RATE-----------------------------------
//
	double
common::func3(
		double const& x)
{
	double return_value = 0; // TODO move inside?
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
	double part1 = 0;
	double part2 = 0;
	if (x <= 1.f) {
		return_value = 0.f;
	}
	else {
		part1 = 1.f / (1.f + ex(-.511f * x / T9mev));
		part2 = 1.f / (1.f + ex(+(x + 2.531f) * (.511f / tnmev) + xi[1]));
		return_value = cnorm * x * fem::pow2((x + 2.531f)) * pow(
				(fem::pow2(x) - 1), .5f) * part1 * part2;	// TODO change to sqrt.
	}
	return return_value;
}


//
//40--------2ND PART OF INTEGRAL FOR p->n RATE-----------------------------------
//
	double
common::func4(
		double const& x)
{
	double return_value = 0;
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
	double part1 = 0;
	double part2 = 0;
	if (x <= 1.f) {
		return_value = 0.f;
		return 0;
	}
	else {
		part1 = 1.f / (1.f + ex(+.511f * x / T9mev));
		part2 = 1.f / (1.f + ex(-(x - 2.531f) * (.511f / tnmev) + xi[1]));
		return_value = cnorm * x * fem::pow2((x - 2.531f)) * pow(
				(fem::pow2(x) - 1), .5f) * part1 * part2;	// TODO change to sqrt.
	}
	return return_value;
	// TODO test and uncomment
	/*
	   if (x <= 1.f) {
	   return 0;
	   }
	   else {
	   const double me = 0.511;
	   const double K = x - 2.531;
	   const double part1 = 1 + exp(+me * x / cmn.T9mev);
	   const double part2 = 1 + exp(-K*(.511f / cmn.tnmev) + xi[1]);
	   return cmn.cnorm * x * K * K * sqrt(x*x - 1) / (part1 * part2);
	   }
	 */
}


//
//50--------INTEGRAL FOR ENERGY DENSITY OF NEUTRINO------------------------------
//
	double
common::func5(
		double const& x)
{
	double return_value = 0;
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
	return_value = 1.f / (2 * fem::pow2(3.14159f)) * x*x*x / (
			1.f + exp(x / tnu - xi[nu]));
	return return_value;
}


//
//60--------INTEGRAL FOR ENERGY DENSITY OF ANTINEUTRINO--------------------------
//
	double
common::func6(
		double const& x)
{
	double return_value = 0;
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
	return_value = 1.f / (2 * fem::pow2(3.14159f)) * x*x*x / (
			1.f + exp(x / tnu + xi[nu]));
	return return_value;
}


//
//----------REFERENCES-----------------------------------------------------------
//     Forms of the integrals involved can be found in
//       Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
//       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
//
//========================IDENTIFICATION DIVISION================================
//
double
common::xintd(
		double const& xlow,
		double const& xhi,
		//func1_function_pointer func,
		double(* func)(common &, const double&),
		int const& nq)
{
	const int np = 6;
	static const double u[] = {
		0,
		-.93246951420315f, -.66120938646627f, -.23861918608320f,
		.23861918608320f, .66120938646627f, .93246951420315f
	};
	static const double w[] = {
		0,
		.17132449237917f, .36076157304814f, .46791393457269f,
		.46791393457269f, .36076157304814f, .17132449237917f
	};
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
	double sum = 0;
	int nint = 0;
	int npnt = 0;
	double dist = (xhi - xlow) / fem::ffloat(nq); 				/// Size of quad interval.
	FEM_DO_SAFE(nint, 1, nq) {
		double cent = xlow + (fem::ffloat(nint) - 0.5) * dist; 	/// Center of interval.
		FEM_DO_SAFE(npnt, 1, np) {
			//x = cent + 0.5f * dist * u(npnt); 					/// Integration point.
			double x = cent + 0.5 * dist * u[npnt]; 				/// Integration point.
			double y = func( x); 								/// Evaluate function x(1).
			//sum += f * w(npnt); 								/// Add up sum.
			sum += y * w[npnt]; 								/// Add up sum.
		}
	}
	//
	//20--------GET INTEGRAL VALUE---------------------------------------------------
	//
	//Do integral.
	return 0.5 * sum * dist;
}


//
//========================IDENTIFICATION DIVISION================================
//
void common::rate1(
		double const& tph)
{
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
	//xi[1] is e neutrino degeneracy param
	//xi[2] is m neutrino degeneracy param
	//xi[3] is t neutrino degeneracy param
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
	//arr_1d<2, double> w(fem::fill0);
	//arr_1d<2, double> x(fem::fill0);
	//arr_1d<2, double> y(fem::fill0);
	//arr_1d<2, double> z(fem::fill0);
	double _w[2+1];
	double _x[2+1];
	double _y[2+1];
	double _z[2+1];
	double uplim1 = 0;
	double uplim2 = 0;
	double uplim3 = 0;
	double uplim4 = 0;
	const int iter = 50;
	double part1 = 0;
	double part2 = 0;
	double part3 = 0;
	double part4 = 0;
	if (xi[1] == 0) {
		//f(1) = thm(13) / tau; 	/// Forward rate for weak np reaction.
		f[1] = thm(13) / tau; 	/// Forward rate for weak np reaction.
		//r(1) = thm(14) / tau; 	/// Reverse rate for weak np reaction.
		r[1] = thm(14) / tau; 	/// Reverse rate for weak np reaction.
	}
	else {
		//
		//20--------COMPUTE WEAK REACTION RATES (DEGENERATE)-----------------------------
		//
		T9mev = tph * .086171f; //Convert photon temp to units of MeV.
		tnmev = tnu * .086171f; //Convert neutrino temp to units of Me
		//..........COMPUTE OVERFLOW LIMITS FOR LIMITS OF INTEGRATION (Ref 1 & 2).
		_w[1] = (-(T9mev / .511f) * (-88.722f));
		_w[2] = ((tnmev / .511f) * (88.029f + xi[1]) + 2.531f);
		_x[1] = ((T9mev / .511f) * (88.029f));
		_x[2] = (-(tnmev / .511f) * (-88.722f + xi[1]) - 2.531f);
		_y[1] = (-(T9mev / .511f) * (-88.722f));
		_y[2] = ((tnmev / .511f) * (88.029f - xi[1]) - 2.531f);
		_z[1] = ((T9mev / .511f) * (88.029f));
		_z[2] = (-(tnmev / .511f) * (-88.722f - xi[1]) + 2.531f);
		//..........COMPARE LIMITS AND TAKE LARGER OF THE TWO.
		uplim1 = abs(_w[1]);
		uplim2 = abs(_x[1]);
		uplim3 = abs(_y[1]);
		uplim4 = abs(_z[1]);
		if (uplim1 < abs(_w[2])) {
			uplim1 = _w[2];
		}
		if (uplim2 < abs(_x[2])) {
			uplim2 = _x[2];
		}
		if (uplim3 < abs(_y[2])) {
			uplim3 = _y[2];
		}
		if (uplim4 < abs(_z[2])) {
			uplim4 = _z[2];
		}
		//..........EVALUATE THE INTEGRALS NUMERICALLY.
		part1 = xintd( 1., uplim1, func1, iter);
		part2 = xintd( 1., uplim2, func2, iter);
		part3 = xintd( 1., uplim3, func3, iter);
		part4 = xintd( 1., uplim4, func4, iter);
		//f(1) = part1 + part2; 		/// Add 2 integrals to get forward rate.
		f[1] = part1 + part2; 		/// Add 2 integrals to get forward rate.
		//r(1) = part3 + part4; 		/// Add 2 integrals to get reverse rate.
		r[1] = part3 + part4; 		/// Add 2 integrals to get reverse rate.
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
void common::start()
{
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
	//c[1] is variation of gravitational constant.
	//c[2] is neutron lifetime (sec).
	//c[3] is number of neutrino species.
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
	//Defined by z = m(electron)*c**2/k*T9
	//
	//===========================PROCEDURE DIVISION==================================
	//
	//10--------INITIALIZE FLAGS AND COUNTERS----------------------------------------
	//
	ltime = 0; 									/// No output yet.
	is = 1; 									/// First iteration coming up.
	ip = inc; 								/// Set to maximum allowed # of iterations.
	it = 0; 									/// No accumulation yet.
	mbad = 0; 									/// No computational errors.
	//
	//20--------SETTINGS-------------------------------------------------------------
	//
	//..........COMPUTATIONAL SETTINGS.
	T9 = T9i; 									/// Initial temperature.
	tnu = T9; 									/// Initial neutrino temperature.
	const double const1 = 0.09615f; 				/// Initial time (Ref 1).
	t = 1 / fem::pow2((const1 * T9));
	dt = dt1; 								/// Initial time step.
	//..........MODEL SETTINGS.
	const double const2 = 6.6700e-8f; 			/// Modify gravitational constant.
	g = const2 * c[1];
	tau = c[2]; 								/// Convert n half-life (min) to lifetime (secs).
	tau = tau / 0.98f; 							/// Coulomb correction (Ref 2). 
												//  TODO <-- check this!
												/// This does not have enough digits 
												/// for today's lifetime measured values.
	xnu = c[3]; 								/// Number of neutrino species.
	//
	//30--------COMPUTE INITIAL ABUNDANCES FOR NEUTRON AND PROTON--------------------
	//
	if ((15.011f / T9 + xi[1]) > 58.f) { 		/// Overabundance of anti-neutrinos.
		y[1] = 1.e-25f; 						/// Very little of neutrons.
		y[2] = 1.f; 							/// Essentially all protons.
	}
	else {
		if ((15.011f / T9 + xi[1]) <  - 58.f) { 	/// Overabundance of neutrinos.
			y[1] = 1.f; 							/// Essentially all neutrons.
			y[2] = 1.e-25f; 						/// Very little of protons.
		}
		else {
			y[1] = 1/(ex(15.011/T9 + xi[1]) + 1); 	/// Initial n abundance (Ref 3).
			y[2] = 1/(ex(-15.011/T9 - xi[1]) + 1);	/// Initial p abundance (Ref 3).
		}
	}
	if (xi[1] != 0) { 								/// Electron neutrino degeneracy.
		cnorm = 1.;
		tnu = .00001f; 								/// Low temperature.
		rate1( 0.00001f); 						/// Find normalization constant at low temperature.
		//cnorm = 1 / tau / f(1);
		cnorm = 1 / tau / f[1];
	}
	y0[1] = y[1];
	y0[2] = y[2];
	//
	//40--------FIND RATIO OF BARYON DENSITY TO TEMPERATURE CUBED--------------------
	//
	double z = 5.930f / T9; 				//Inverse of temperature.
	//bessel( z); // old call
	double bl1 = getBesselL(z);
	double bl2 = getBesselL(2*z);
	double bl3 = getBesselL(3*z);
	double bl4 = getBesselL(4*z);
	double bl5 = getBesselL(5*z);
	/*
	   double bm1 = getBesselM(z);
	   double bm2 = getBesselM(2*z);
	   double bm3 = getBesselM(3*z);
	   double bm4 = getBesselM(4*z);
	   double bm5 = getBesselM(5*z);
	   double bn1 = getBesselN(z);
	   double bn2 = getBesselN(2*z);
	   double bn3 = getBesselN(3*z);
	   double bn4 = getBesselN(4*z);
	   double bn5 = getBesselN(5*z);
	 */
	hv = 3.3683e+4f * eta1 * 2.75; 		/// (Ref 4 but with final eta).
	phie = hv * (1.784e-5f * y[2]) / 
		(0.5*z*z*z*(bl1 - 2*bl2 + 3*bl3 - 4*bl4 + 5*bl5));
	/// Chemical potential of electron (Ref 5).
	rhob0 = hv * fem::pow3(T9); 					/// TODO Baryon density. 
	//Nonde
	if ((xi[1] == 0) && (xi[2] == 0) && (xi[3] == 0)) {
		rhone0 = 7.366f * fem::pow4(T9); 	/// Electron neutrino density (Ref 6).
	}
	//
	//50--------SET ABUNDANCES FOR REST OF NUCLIDES----------------------------------
	//
	y[3] = y[1] * y[2] * rhob0 * ex(25.82f / T9) / (.471e+10f * pow(T9, 1.5f)); /// (Ref 7).
	y0[3] = y[3];
	int i = 0;
	FEM_DO_SAFE(i, 4, isize) {
		y[i] = ytmin; 						/// Set rest to minimum abundance.
		y0[i] = y[i]; 						/// Initialize abundances at beginning of iter
	}
	//Compute weak decay rates.
	rate0();
	//
	//----------REFERENCES-----------------------------------------------------------
	//     1) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
	//          page 44, equation A15.
	//     2) Coulomb correction obtained by dividing by correction factor Fp(T9)
	//               Fp(T9) = 1 - 0.5(pi/(137<v>/c))
	//          Wagoner, R.V. 1973, Ap. J. 179, page 358.
	//     3) For the nondegenerate case:
	//          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
	//          page 4, equation 3.
	//        For the case with neutrino degeneracy:
	//          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49,
	//          page 417, equation 9.
	//     4) Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 250, equation 4.
	//          3.3683e+4 = Mu(ng/T9**3) with Mu the atomic mass, ng the
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
void common::nudens()
{
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
	double uplim1 = 0;
	double uplim2 = 0;
	const int iter = 50;
	if (abs(xi[nu]) <= 0.03f) {
		//..........SMALL xi APPROXIMATION.
		rhonu = 2.f * (fem::pow2(3.14159f) / 30.f) * fem::pow4((tnu))
			* (7.f / 8.f + (15.f / (4 * fem::pow2(3.14159f))) * fem::pow2(xi[nu]) 
					+ (15.f / (8.f * fem::pow4(3.14159f))) * fem::pow4(xi[nu]));
	}
	else {
		if (abs(xi[nu]) >= 30.f) {
			//..........LARGE xi APPROXIMATION.
			rhonu = (fem::pow4((tnu))) / (8.f * fem::pow2(3.14159f)) *
				fem::pow4(xi[nu]) * (1 + 12.f * 1.645f / fem::pow2(xi[nu]));
		}
		else {
			//..........DO INTEGRATION
			uplim1 = (88.029f + xi[nu]) * tnu;
			uplim2 = (88.029f - xi[nu]) * tnu;
			if (uplim2 <= 0.) {
				rhonu = xintd( 0, uplim1, func5, iter);
			}
			else {
				rhonu = xintd(0, uplim1, func5, iter) 
					+ xintd(0, uplim2, func6, iter);
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
void common::therm()
{
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
	//Defined by z = m(electron)*c**2/k*T9
	//
	//===========================PROCEDURE DIVISION==================================
	//
	//10--------COMPUTE FACTORS------------------------------------------------------
	//
	//z = m(electron)c**2/k(T9).
	double z = 5.930f / T9;
	double bl1 = getBesselL(z);
	double bl2 = getBesselL(2*z);
	double bl3 = getBesselL(3*z);
	double bl4 = getBesselL(4*z);
	double bl5 = getBesselL(5*z);
	double bm1 = getBesselM(z);
	double bm2 = getBesselM(2*z);
	double bm3 = getBesselM(3*z);
	double bm4 = getBesselM(4*z);
	double bm5 = getBesselM(5*z);
	double bn1 = getBesselN(z);
	double bn2 = getBesselN(2*z);
	double bn3 = getBesselN(3*z);
	double bn4 = getBesselN(4*z);
	double bn5 = getBesselN(5*z);
	//Neutrino temperature.
	tnu = (pow((rnb), (1.f / 3.f))) * T9i;
	//..........FACTORS OF z.
	double z1 = z;
	double z2 = z*z;
	double z3 = z*z*z;
	double z4 = z*z*z*z;
	double z5 = z*z*z*z*z;
	//..........TRIGONOMETRIC FUNCTION VALUES.
	//No chance of overflow.
	double cosh1 = 0;
	double cosh2 = 0;
	double cosh3 = 0;
	double cosh4 = 0;
	double cosh5 = 0;
	double sinh1 = 0;
	double sinh2 = 0;
	double sinh3 = 0;
	double sinh4 = 0;
	double sinh5 = 0;
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
	//
	//20--------COMPUTE THERMODYNAMIC VARIABLES--------------------------------------
	//
	//(Ref 1)
	thm(1) = 8.418f * T9 * T9 * T9 * T9;
	//(Ref 2)
	thm(2) = 4.f * thm(1) / T9;
	//(Ref 3)
	thm(3) = thm(1) / 3.f;
	//(Ref 4)
	thm(4) = 3206.f * (bm1 * cosh1 - bm2 * cosh2 + bm3 * cosh3 - bm4 *
			cosh4 + bm5 * cosh5);
	//(Ref 5)
	thm(5) = 3206.f * (z / T9) * (bn1 * cosh1 - 2.f * bn2 *
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
	if ((xi[1] == 0) && (xi[2] == 0) && (xi[3] == 0)) {
		//(Ref 8)
		thm(8) = xnu * rhone0 * (pow(rnb, (4.f / 3.f)));
		//Include effects of neutrino degenera
	}
	else {
		thm(8) = 0.f;
		//For every neutrino family.
		FEM_DO_SAFE(nu, 1, xnu) {
			//Compute neutrino energy density.
			nudens();
			//Have 12.79264 from units ch
			thm(8) += 12.79264f * rhonu;
		}
	}
	//(Ref 9)
	thm(9) = rhob0 * rnb;
	//(Ref 10
	thm(10) = thm(1) + thm(4) + thm(8) + thm(9);
	//(Ref 11
	thm(11) = -(fem::pow3(z) / T9) * (sinh1 * (3.f * bl1 - z * bm1) - // TODO change to z3
			sinh2 * (3.f * bl2 - 2.f * z * bm2) + sinh3 * (3.f * bl3 - 3.f *
				z * bm3) - sinh4 * (3.f * bl4 - 4.f * z * bm4) + sinh5 * (3.f *
					bl5 - 5.f * z * bm5));
	//(Ref 12
	thm(12) = fem::pow3(z) * (cosh1 * bl1 - 2.f * cosh2 * bl2 + 3.f * // TODO change to z3
			cosh3 * bl3 - 4.f * cosh4 * bl4 + 5.f * cosh5 * bl5);
	if (thm(12) != 0.f) {
		thm(12) = 1.f / thm(12);
	}
	//(Ref 13
	thm(13) = 1.000f + 0.565f / z1 - 6.382f / z2 + 11.108f / z3 +
		36.492f / z4 + 27.512f / z5;
	//(Ref 14
	const double q = 2.531f;
	thm(14) = (5.252f / z1 - 16.229f / z2 + 18.059f / z3 + 34.181f /
			z4 + 27.617f / z5) * ex(-q * z);
	//
	//----------REFERENCES AND NOTES-------------------------------------------------
	//     1)  thm(1)  = rho photon
	//         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
	//          page 43, equation A2.)
	//     2)  thm(2)  = d(rho photon)/d(T9)
	//     3)  thm(3)  = (p photon)/c**2
	//         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967,
	//          page 43, equation A3.)
	//     4)  thm(4)  = rho electron+positron
	//         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9,
	//          page 281, equation B44.)
	//     5)  thm(5)  = d(rho electron+positron)/d(T9)
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
	//     11) thm(11) = d     [pi**2(hbar*c)**3(ne- - ne+)*z**3]
	//                   d(T9) [  2  (mc**2)**3                 ]
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
void common::eqslin(
		int const& icnvm,
		int& ierror)
{
	int nord = 0;
	int i = 0;
	double x[nnuc+1];
	int j = 0;
	double a0[nnuc+1][nnuc+1];			// TODO fix to zero indexing.
	double cx = 0;
	int k = 0;
	double sum = 0;
	double xdy = 0;
	const double eps = 2.e-4f;
	const int mord = 1;
	double r = 0;
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
		x[i] = b[i]; 		/// Right-hand vector.
		y[i] = 0; 			/// Solution vector.
	}
	//..........SAVE MATRIX.
	//Monitor convergence.
	if (icnvm == inc) {
		FEM_DO_SAFE(i, 1, isize) {
			FEM_DO_SAFE(j, 1, isize) {
				a0[j][i] = a[j][i]; 			/// Initial value of coefficient array.
			}
		}
	}
	//
	//20--------TRIANGULARIZE MATRIX AND SAVE OPERATOR-------------------------------
	//
	//..........CHECK TO SEE THAT THERE ARE NO ZEROS AT PIVOT POINTS.
	FEM_DO_SAFE(i, 1, isize - 1) {
		if (a[i][i] == 0) { 					/// Don't want to divide by zero.
			mbad = i; 							/// Position of zero coefficient.
			std::cout << "Found divide by zero on 'a' diagonal.\n";
			return; 							/// Terminate matrix evaluation.
		}
	//..........TRIANGULARIZE MATRIX.
		FEM_DO_SAFE(j, i + 1, isize) {
			if (a[j][i] != 0) { 				/// Progress diagonally down the column.
				cx = a[j][i] / a[i][i]; 		/// Scaling factor down the column.
				FEM_DO_SAFE(k, i + 1, isize) { 	/// Progress diagonally along row.
					a[j][k] -= cx * a[i][k]; 	/// Subtract scaled coefficient along row.
				}
				a[j][i] = cx; 					/// Scaled coefficient.
	//..........OPERATE ON RIGHT-HAND VECTOR.
				x[j] -= cx * x[i]; 				/// Subtract off scaled coefficient.
			}
		}
	}
	//
	//30--------DO BACK SUBSTITUTION-------------------------------------------------
	//
statement_300:
	x[isize] = x[isize] / a[isize][isize]; 		/// Solution for ultimate position.
	y[isize] += x[isize];
	FEM_DOSTEP(i, isize - 1, 1, -1) { 			/// From i = penultimate to i = 1.
		sum = 0;
		FEM_DO_SAFE(j, i + 1, isize) {
			sum += a[i][j] * x[j]; 				/// Sum up all previous terms.
		}
		x[i] = (x[i] - sum) / a[i][i];
		y[i] += x[i]; 							/// Add difference to initial value.
	}
	//
	//40--------TESTS AND EXITS------------------------------------------------------
	//
	if (icnvm == inc) {
		FEM_DO_SAFE(i, 1, isize) {
			if (y[i] != 0.) {
				xdy = fem::dabs(x[i] / y[i]); 			/// Relative error.
				if (xdy > eps) {
					if (nord < mord) { 					/// Continue to higher orders.
						nord++;
						//..........FIND ERROR IN RIGHT-HAND VECTOR.
						FEM_DO_SAFE(j, 1, isize) {
							r = 0; 						/// Initialize r.
							FEM_DO_SAFE(k, 1, isize) {
								r += a0[j][k] * y[k]; 	/// Left side with approximate sol
							}
							x[j] = b[j] - r; 			/// Subtract difference from right side.
						}
						//..........OPERATE ON RIGHT-HAND VECTOR.
						FEM_DO_SAFE(j, 1, isize - 1) {
							FEM_DO_SAFE(k, j + 1, isize) {
								x[k] -= a[k][j] * x[j]; /// Subtract off scaled coefficient.
							}
						}
						goto statement_300; 			/// Go for another iteration.
					}
					else {
						//..........NOT ENOUGH CONVERGENCE.
						mbad = -1; 						/// Signal error problem.
						ierror = i; 					/// i'th nuclide for which x/y checked.
						return;
					} //(nord < mord)
				} //(xdy > eps)
			} //(y(i) != 0)
		} //i = 1...isize
	} //(icnvm == inc)
	//No more iterations & relative error
	//
}


//
//========================IDENTIFICATION DIVISION================================
//
void common::sol(
		int const& loop)
{
	common_write write(*this);
	const double si[] = {NOT_USED, 1, 1, 1, 1, 1, 2, 3, 2, 1, 1, 2};
	const double sj[] = {NOT_USED, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0};
	const double sk[] = {NOT_USED, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2};
	const double sl[] = {NOT_USED, 1, 1, 1, 2, 2, 1, 1, 1, 2, 3, 1};
	double T932 = 0;
	double T9m32 = 0;
	int isize1 = 0;
	int i = 0;
	int j = 0;
	int n = 0;
	int ind = 0;
	int k = 0;
	int l = 0;
	double ri = 0;
	double rj = 0;
	double rk = 0;
	double rl = 0;
	double ci = 0;
	double cj = 0;
	double ck = 0;
	double cl = 0;
	double bdln = 0;
	int ierror = 0;
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
	//T9**(3/2).
	T932 = pow(T9, 1.5);
	//T9**(-3/2).
	T9m32 = 1. / T932;
	//..........MATRIX SIZE.
	isize1 = isize + 1;
	//..........INITIALIZE A-MATRIX.
	FEM_DO_SAFE(i, 1, isize) {
		FEM_DO_SAFE(j, 1, isize) {
			a[j][i] = 0; 					/// Set a-matrix to zero.
		}
	}
	//
	//20--------COMPUTE FACTORS FOR THE A-MATRIX-------------------------------------
	//
	FEM_DO_SAFE(n, 1, jsize) {
		//..........EQUATE VARIABLES TO ARRAYS.
		ind = iform(n); 					/// Type of reaction.
		i = ii(n); 							/// ID # of incoming nuclide i.
		j = jj(n); 							/// ID # of incoming nuclide j.
		k = kk(n); 							/// ID # of outgoing nuclide k.
		l = ll(n); 							/// ID # of outgoing nuclide l.
		//Reactio
		if ((ind != 0) && (i <= isize) && (l <= isize)) {
			ri = si[ind]; 					/// # of incoming nuclide i.
			rj = sj[ind]; 					/// # of incoming nuclide j.
			rk = sk[ind]; 					/// # of outgoing nuclide k.
			rl = sl[ind]; 					/// # of outgoing nuclide l.
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
			//configuration1001:
statement_201: 			/// 1-0-0-1 configuration.
			//ci = f(n); 				/// (Ref 1).
			ci = f[n]; 				/// (Ref 1).
			cj = 0;
			ck = 0;
			//cl = r(n);
			cl = r[n];
			goto statement_212;
statement_202: 			/// 1-1-0-1 configuration.
			//r(n) = rev(n) * 1.e+10f * T932 * ex(-q9(n) / T9) * f(n); 	/// (Ref 2).
			r[n] = rev(n) * 1.e+10f * T932 * ex(-q9(n) / T9) * f[n]; 	/// (Ref 2).
			//f(n) = rhob * f(n);
			f[n] = rhob * f[n];
			//ci = y(j) * f(n) / 2.;
			ci = y[j] * f[n] / 2.;
			//cj = y(i) * f(n) / 2.;
			cj = y[i] * f[n] / 2.;
			ck = 0;
			//cl = r(n);
			cl = r[n];
			goto statement_212;
statement_203: 				/// 1-1-1-1 configuration.
			f[n] = rhob * f[n]; 			/// (Ref 3).
			r[n] = rev(n) * ex(-q9(n) / T9) * f[n];
			ci = y[j] * f[n] / 2.;
			cj = y[i] * f[n] / 2.;
			ck = y[l] * r[n] / 2.;
			cl = y[k] * r[n] / 2.;
			goto statement_212;
statement_204: 				/// 1-0-0-2 configuration.
			ci = f[n];
			cj = 0.;
			ck = 0.;
			cl = y[l] * r[n] / 2.;
			goto statement_212;
statement_205: 				/// 1-1-0-2 configuration.
			f[n] = rhob * f[n];
			r[n] = rev(n) * ex(-q9(n) / T9) * f[n]; 	/// (Ref 3).
			ci = y[j] * f[n] / 2.;
			cj = y[i] * f[n] / 2.;
			ck = 0.;
			cl = y[l] * r[n] / 2.;
			goto statement_212;
statement_206: 				/// 2-0-1-1 configuration.
			f[n] = rhob * f[n];
			r[n] = rev(n) * ex(-q9(n) / T9) * f[n]; 	/// (Ref 3).
			ci = y[i] * f[n] / 2.;
			cj = 0.;
			ck = y[l] * r[n] / 2.;
			cl = y[k] * r[n] / 2.;
			goto statement_212;
			//3-0-0-1 configuration.
statement_207:
			//(Ref 4).
			r[n] = rev(n) * 1.e+20f * T932 * T932 * ex(-q9(n) / T9) * f[n];
			f[n] = rhob * rhob * f[n];
			ci = y[i] * y[i] * f[n] / 6.;
			cj = 0.;
			ck = 0.;
			cl = r[n];
			goto statement_212;
			//2-1-0-1 configuration.
statement_208:
			//(Ref 4).
			r[n] = rev(n) * 1.e+20f * T932 * T932 * ex(-q9(n) / T9) * f[n];
			f[n] = rhob * rhob * f[n];
			ci = y[j] * y[i] * f[n] / 3.;
			cj = y[i] * y[i] * f[n] / 6.;
			ck = 0.;
			cl = r[n];
			goto statement_212;
			//1-1-1-2 configuration.
statement_209:
			f[n] = rhob * f[n];
			//(Ref 5)
			r[n] = rev(n) * 1.e-10f * T9m32 * rhob * ex(-q9(n) / T9) * f[n];
			ci = y[j] * f[n] / 2.;
			cj = y[i] * f[n] / 2.;
			ck = y[l] * y[l] * r[n] / 6.;
			cl = y[k] * y[l] * r[n] / 3.;
			goto statement_212;
			//1-1-0-3 configuration.
statement_210:
			f[n] = rhob * f[n];
			//(Ref 5)
			r[n] = rev(n) * 1.e-10f * T9m32 * rhob * ex(-q9(n) / T9) * f[n];
			ci = y[j] * f[n] / 2.;
			cj = y[i] * f[n] / 2.;
			ck = 0.;
			cl = y[l] * y[l] * r[n] / 6.;
			goto statement_212;
			//2-0-2-1 configuration.
statement_211:
			f[n] = rhob * f[n];
			//(Ref 5)
			r[n] = rev(n) * 1.e-10f * T9m32 * rhob * ex(-q9(n) / T9) * f[n];
			ci = y[i] * f[n] / 2.;
			cj = 0.;
			ck = y[l] * y[k] * r[n] / 3.;
			cl = y[k] * y[k] * r[n] / 6.;
statement_212:
			//
			//30--------CONSTRUCT THE A-MATRIX-----------------------------------------------
			//
			i = isize1 - i; 				/// Invert i index.
			j = isize1 - j; 				/// Invert j index.
			k = isize1 - k; 				/// Invert k index.
			l = isize1 - l; 				/// Invert l index.
			//..........FILL I NUCLIDE COLUMN.
			if (j <= isize) {
				a[j][i] += rj * ci;
			}
			if (k <= isize) {
				a[k][i] -= rk * ci;
			}
			a[i][i] += ri * ci;
			a[l][i] -= rl * ci;
			//..........FILL J NUCLIDE COLUMN.
			if (j <= isize) {
				a[j][j] += rj * cj;
				if (k <= isize) {
					a[k][j] -= rk * cj;
				}
				a[i][j] += ri * cj;
				a[l][j] -= rl * cj;
			}
			//..........FILL K NUCLIDE COLUMN.
			if (k <= isize) {
				if (j <= isize) {
					a[j][k] -= rj * ck;
				}
				a[k][k] += rk * ck;
				a[i][k] -= ri * ck;
				a[l][k] += rl * ck;
			}
			//..........FILL L NUCLIDE COLUMN.
			if (j <= isize) {
				a[j][l] -= rj * cl;
			}
			if (k <= isize) {
				a[k][l] += rk * cl;
			}
			a[i][l] -= ri * cl;
			a[l][l] += rl * cl;
		} //((ind.ne.0).and.(i.le.isize).and.(l.le.isize))
	} //n = 1,jsize
	//
	//40--------PUT A-MATRIX AND B-VECTOR IN FINAL FORM OF MATRIX EQUATION-----------
	//
	//(10**(-5))*(Expansion rate).
	bdln = 1.e-5f * (3. * hubcst);
	FEM_DO_SAFE(i, 1, isize) {
		int i1 = isize1 - i; 						/// Invert the rows.
		FEM_DO_SAFE(j, 1, isize) {
			int j1 = isize1 - j; 						/// Invert the columns.
			//std::cout << "i:"<<i<<" j:"<<j<<" i1:"<<i1<<" j1:"<<j1<<std::endl;
			if (y0[i1] == NOT_USED) exit(0);
			if (y0[j1] == NOT_USED) exit(0);
			if (a[j][i] == NOT_USED) exit(0);

			if (abs(a[j][i]) < bdln * y0[j1] / y0[i1])
				a[j][i] = 0; 							/// Set 0 if tiny.
			else
				a[j][i] *= dt; 							/// Bring dt over to other side.
		}
		a[i][i] += 1; 								/// Add identity matrix to a-matrix.
		b[i1] = y0[i]; 								/// Initial abundances.
	}
	//
	//50--------SOLVE EQUATIONS TO GET DERIVATIVE------------------------------------
	//
	//..........SET MONITOR FLAG AND SOLVE BY GAUSSIAN ELIMINATION.
	if (loop == 1) {
		eqslin(ip, ierror);
	}
	else {
		eqslin(0, ierror);
	}
	//..........OBTAIN DERIVATIVE.
	FEM_DO_SAFE(i, 1, isize) {
		//yy[i] = yx[isize1 - i]; 					/// Abundance at t+dt.
		//dydt[i] = (yy[i] - y0[i]) / dt; 			/// Take derivative.
		//dydt[i] = (y[isize1-i] - y0[i]) / dt; 		/// Take derivative.
		dydt[i] = (y[i] - y0[i]) / dt; 		/// Take derivative.
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
void common::rate2()
{
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
	//T9**(1/3)
	double T913 = pow(T9, (.33333333f));
	//T9**(2/3)
	double T923 = T913 * T913;
	//T9**(4/3)
	double T943 = T923 * T923;
	//T9**(5/3)
	double T953 = T9 * T923;
	//T9**(1/2)
	double T912 = sqrt(T9);
	//T9**(3/2)
	double T932 = T9 * T912;
	//T9**(-1)
	double T9m1 = 1 / T9;
	//T9**(-2/3)
	double T9m23 = 1.0f / T923;
	//T9**(-3/2)
	double T9m32 = 1.0f / T932;
	//For reaction 17.
	double T9a = T9 / (1.0f + 13.076f * T9);
	//T9a**(3/2)
	double T9a32 = pow(T9a, (1.5f));
	//For reaction 18.
	double T9b = T9 / (1. + 49.18f * T9);
	//T9b**(3/2)
	double T9b32 = pow(T9b, (1.5f));
	//For reaction 22.
	double T9c = 0;
	if (T9 > 10.) {
		T9c = 1.;
	}
	else {
		T9c = T9 / (1. - 9.69e-2f * T9 + 2.84e-2f * T953 / 
				pow(1.f - 9.69e-2f * T9, (2.f / 3.f)));
	}
	//T9c**(1/3)
	double T9c13 = pow(T9c, (.3333333f));
	//T9c**(5/6)
	double T9c56 = pow(T9c, (.8333333f));
	//For reaction 24.
	double T9d = T9 / (1. + 0.759f * T9);
	//T9d**(1/3)
	double T9d13 = pow(T9d, (.3333333f));
	//T9d**(5/6)
	double T9d56 = pow(T9d, (.8333333f));
	//For reaction 26.
	double T9e = T9 / (1. + 0.1378f * T9);
	//T9e**(1/3)
	double T9e13 = pow(T9e, (.3333333f));
	//T9e**(5/6)
	double T9e56 = pow(T9e, (.8333333f));
	//For reaction 27.
	double T9f = T9 / (1. + 0.1071f * T9);
	//T9f**(1/3)
	double T9f13 = pow(T9f, (.3333333f));
	//T9f**(5/6)
	double T9f56 = pow(T9f, (.8333333f));
	//
	//20--------NEUTRON, PHOTON REACTIONS--------------------------------------------
	//
	//.......H(n,g)H2...................(Smith-Kawano-Malaney 1992)
	f[12] = 4.742e+4f * (1. - .8504f * T912 + .4895f * T9 - .09623f *
			T932 + 8.471e-3f * T9 * T9 - 2.80e-4f * T9 * T932);
	//
	//.......H2(n,g)H3..................(Wagoner 1969)
	f[13] = 6.62e+1f * (1. + 18.9f * T9);
	//
	//.......He3(n,g)He4................(Wagoner 1969)
	f[14] = 6.62e+0f * (1. + 905. * T9);
	//
	//.......Li6(n,g)Li7................(Malaney-Fowler 1989)
	f[15] = 5.10e+3f;
	//
	//30--------NEUTRON, PROTON REACTIONS--------------------------------------------
	//
	//.......He3(n,p)H3.................(Smith-Kawano-Malaney 1992)
	f[16] = 7.21e+8f * (1. - .508f * T912 + .228f * T9);
	//
	//.......Be7(n,p)Li7................(Smith-Kawano-Malaney 1992)
	f[17] = 2.675e+9f * (1. - .560f * T912 + .179f * T9 - .0283f *
			T932 + 2.214e-3f * T9 * T9 - 6.851e-5f * T9 * T932) + 9.391e+8f *
		T9a32 * T9m32 + 4.467e+7f * T9m32 * ex(-0.07486f / T9);
	//
	//40--------NEUTRON, ALPHA REACTIONS---------------------------------------------
	//
	//.......Li6(n,a)H3.................(Caughlan-Fowler 1988)
	f[18] = 2.54e+9f * T9m32 * ex(-2.39f / T9) + 1.68e+8f * (1. -
			.261f * T9b32 / T932);
	//
	//.......Be7(n,a)He4................(Wagoner 1969)
	f[19] = 2.05e+4f * (1. + 3760. * T9);
	//
	//50--------PROTON, PHOTON REACTIONS---------------------------------------------
	//
	//.......H2(p,g)He3.................(Smith-Kawano-Malaney 1992)
	f[20] = 2.65e+3f * T9m23 * ex(-3.720f / T913) * (1. + .112f *
			T913 + 1.99f * T923 + 1.56f * T9 + .162f * T943 + .324f * T953);
	//
	//.......H3(p,g)He4.................(Caughlan-Fowler 1988)
	f[21] = 2.20e+4f * T9m23 * ex(-3.869f / T913) * (1. + .108f *
			T913 + 1.68f * T923 + 1.26f * T9 + .551f * T943 + 1.06f * T953);
	//
	//.......Li6(p,g)Be7................(Caughlan-Fowler 1988)
	f[22] = 6.69e+5f * T9c56 * T9m32 * ex(-8.413f / T9c13);
	//
	//60--------PROTON, ALPHA REACTIONS----------------------------------------------
	//
	//.......Li6(p,a)He3................(Caughlan-Fowler 1988)
	f[23] = 3.73e+10f * T9m23 * ex(-8.413f / T913 - fem::pow2((T9 /
					5.50f))) * (1. + .050f * T913 - .061f * T923 - .021f * T9 +
				.006f * T943 + .005f * T953) + 1.33e+10f * T9m32 * ex(-17.763f /
					T9) + 1.29e+09f * T9m1 * ex(-21.820f / T9);
	//
	//.......Li7(p,a)He4................(Smith-Kawano-Malaney 1992)
	f[24] = 1.096e+9f * T9m23 * ex(-8.472f / T913) - 4.830e+8f *
		T9d56 * T9m32 * ex(-8.472f / T9d13) + 1.06e+10f * T9m32 * ex(
				-30.442f / T9) + 1.56e+5f * T9m23 * ex((-8.472f / T913) -
					fem::pow2((T9 / 1.696f))) * (1. + .049f * T913 - 2.498f * T923 +
					.860f * T9 + 3.518f * T943 + 3.08f * T953) + 1.55e+6f * T9m32 *
					ex(-4.478f / T9);
	//
	//70--------ALPHA, PHOTON REACTIONS----------------------------------------------
	//
	//.......H2(a,g)Li6.................(Caughlan-Fowler 1988)
	f[25] = 3.01e+01f * T9m23 * ex(-7.423f / T913) * (1. + .056f *
			T913 - 4.85f * T923 + 8.85f * T9 - .585f * T943 - .584f * T953) +
		8.55e+1f * T9m32 * ex(-8.228f / T9);
	//
	//.......H3(a,g)Li7.................(Smith-Kawano-Malaney 1992)
	f[26] = 3.032e+5f * T9m23 * ex(-8.090f / T913) * (1. + .0516f *
			T913 + .0229f * T923 + 8.28e-3f * T9 - 3.28e-4f * T943 -
			3.01e-4f * T953) + 5.109e+5f * T9e56 * T9m32 * ex(-8.068f /
				T9e13);
	//
	//.......He3(a,g)Be7................(Smith-Kawano-Malaney 1992)
	f[27] = 4.817e+6f * T9m23 * ex(-14.964f / T913) * (1. + .0325f *
			T913 - 1.04e-3f * T923 - 2.37e-4f * T9 - 8.11e-5f * T943 -
			4.69e-5f * T953) + 5.938e+6f * T9f56 * T9m32 * ex(-12.859f /
				T9f13);
	//
	//80--------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------------
	//
	//.......H2(d,n)He3.................(Smith-Kawano-Malaney 1992)
	f[28] = 3.95e+8f * T9m23 * ex(-4.259f / T913) * (1. + .098f *
			T913 + .765f * T923 + .525f * T9 + 9.61e-3f * T943 + .0167f *
			T953);
	//
	//.......H2(d,p)H3..................(Smith-Kawano-Malaney 1992)
	f[29] = 4.17e+8f * T9m23 * ex(-4.258f / T913) * (1. + .098f *
			T913 + .518f * T923 + .355f * T9 - .010f * T943 - .018f * T953);
	//
	//.......H3(d,n)He4.................(Smith-Kawano-Malaney 1992)
	f[30] = 1.063e+11f * T9m23 * ex(-4.559f / T913 - fem::pow2((T9 /
					.0754f))) * (1. + .092f * T913 - .375f * T923 - .242f * T9 +
				33.82f * T943 + 55.42f * T953) + 8.047e+8f * T9m23 * ex(
					-0.4857f / T9);
	//
	//.......He3(d,p)He4................(Smith-Kawano-Malaney 1992)
	f[31] = 5.021e+10f * T9m23 * ex(-7.144f / T913 - fem::pow2((T9 /
					.270f))) * (1. + .058f * T913 + .603f * T923 + .245f * T9 +
				6.97f * T943 + 7.19f * T953) + 5.212e+8f / T912 * ex(-1.762f /
					T9);
	//
	//90--------THREE PARTICLE REACTIONS---------------------------------------------
	//
	//.......He3(He3,2p)He4.............(Caughlan-Fowler 1988)
	f[32] = 6.04e+10f * T9m23 * ex(-12.276f / T913) * (1. + .034f *
			T913 - .522f * T923 - .124f * T9 + .353f * T943 + .213f * T953);
	//
	//.......Li7(d,na)He4...............(Caughlan-Fowler 1988)
	f[33] = 2.92e+11f * T9m23 * ex(-10.259f / T913);
	//
	//.......Be7(d,pa)He4...............(Caughlan-Fowler 1988)
	f[34] = 1.07e+12f * T9m23 * ex(-12.428f / T913);
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
void common::rate3()
{
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
	//T9**(1/3)
	double T913 = pow(T9, (.33333333f));
	//T9**(2/3)
	double T923 = T913 * T913;
	//T9**(4/3)
	double T943 = T923 * T923;
	//T9**(5/3)
	double T953 = T9 * T923;
	//T9**(1/2)
	double T912 = sqrt(T9);
	//T9**(3/2)
	double T932 = T9 * T912;
	//T9**(1/5)
	double T915 = pow(T9, (.2f));
	//T9**(5/4)
	double T954 = pow(T9, (1.25f));
	//T9**(-1)
	double T9m1 = 1.0f / T9;
	//T9**(-2/3)
	double T9m23 = 1.0f / T923;
	//T9**(-3/2)
	double T9m32 = 1.0f / T932;
	//T9**(-3/4)
	double T9m34 = sqrt(T9m32);
	//T9**(-1/5)
	double T9m15 = 1.0f / T915;
	//T9**(-5/4)
	double T9m54 = 1.0f / T954;
	//For reaction 53.
	double T9a = T9 / (1. + T9 / 15.1f);
	//T9a**(1/3)
	double T9a13 = pow(T9a, (.3333333f));
	//T9a**(5/6)
	double T9a56 = pow(T9a, (.8333333f));
	//
	//20--------NEUTRON, PHOTON REACTIONS--------------------------------------------
	//
	//.......Li7(n,g)Li8................(Wagoner 1969)
	f[35] = 4.90e+3f + 9.96e+3f * T9m32 * ex(-2.62f / T9);
	//
	//.......B10(n,g)B11................(Wagoner 1969)
	f[36] = 6.62e+4f;
	//
	//.......B11(n,g)B12................(Malaney-Fowler 1989)
	f[37] = 7.29e+2f + 2.40e+3f * T9m32 * ex(-0.223f / T9);
	//
	//30--------NEUTRON, PROTON REACTIONS--------------------------------------------
	//
	//.......C11(n,p)B11................(Caughlan-Fowler 1988)
	f[38] = 1.69e+8f * (1. - .048f * T912 + .010f * T9);
	//
	//40--------NEUTRON, ALPHA REACTIONS---------------------------------------------
	//
	//.......B10(n,a)Li7................(Caughlan-Fowler 1988)
	f[39] = 5.07e+8f;
	//
	//50--------PROTON, PHOTON REACTIONS---------------------------------------------
	//
	//.......Be7(p,g)B8.................(Caughlan-Fowler 1988)
	f[40] = 3.11e+5f * T9m23 * ex(-10.262f / T913) + 2.53e+3f * T9m32 *
		ex(-7.306f / T9);
	//
	//.......Be9(p,g)B10................(Caughlan-Fowler 1988)
	f[41] = 1.33e+7f * T9m23 * ex(-10.359f / T913 - fem::pow2((T9 /
					.846f))) * (1. + .040f * T913 + 1.52f * T923 + .428f * T9 +
				2.15f * T943 + 1.54f * T953) + 9.64e+4f * T9m32 * ex(-3.445f /
					T9) + 2.72e+6f * T9m32 * ex(-10.620f / T9);
	//
	//.......B10(p,g)C11................(Caughlan-Fowler 1988)
	f[42] = 4.61e+5f * T9m23 * ex(-12.062f / T913 - fem::pow2((T9 /
					4.402f))) * (1. + .035f * T913 + .426f * T923 + .103f * T9 +
				.281f * T943 + .173f * T953) + 1.93e+5f * T9m32 * ex(-12.041f /
					T9) + 1.14e+4f * T9m32 * ex(-16.164f / T9);
	//
	//.......B11(p,g)C12................(Caughlan-Fowler 1988)
	f[43] = 4.62e+7f * T9m23 * ex(-12.095f / T913 - fem::pow2((T9 /
					.239f))) * (1. + .035f * T913 + 3.00f * T923 + .723f * T9 +
				9.91f * T943 + 6.07f * T953) + 7.89e+3f * T9m32 * ex(-1.733f /
					T9) + 9.68e+4f * T9m15 * ex(-5.617f / T9);
	//
	//.......C11(p,g)N12................(Caughlan-Fowler 1988)
	f[44] = 4.24e+4f * T9m23 * ex(-13.658f / T913 - fem::pow2((T9 /
					1.627f))) * (1. + .031f * T913 + 3.11f * T923 + .665f * T9 +
				4.61f * T943 + 2.50f * T953) + 8.84e+3f * T9m32 * ex(-7.021f /
					T9);
	//
	//60--------PROTON, NEUTRON REACTIONS--------------------------------------------
	//
	//.......B12(p,n)C12................(Wagoner 1969)
	f[45] = 4.02e+11f * T9m23 * ex(-12.12f / T913);
	//
	//70--------PROTON, ALPHA REACTIONS----------------------------------------------
	//
	//.......Be9(p,a)Li6................(Caughlan-Fowler 1988)
	f[46] = 2.11e+11f * T9m23 * ex(-10.359f / T913 - fem::pow2((T9 /
					.520f))) * (1. + .040f * T913 + 1.09f * T923 + .307f * T9 +
				3.21f * T943 + 2.30f * T953) + 4.51e+8f * T9m1 * ex(-3.046f /
					T9) + 6.70e+8f * T9m34 * ex(-5.160f / T9);
	//
	//.......B10(p,a)Be7................(Caughlan-Fowler 1988)
	f[47] = 1.26e+11f * T9m23 * ex(-12.062f / T913 - fem::pow2((T9 /
					4.402f))) * (1. + .035f * T913 - .498f * T923 - .121f * T9 +
				.300f * T943 + .184f * T953) + 2.59e+9f * T9m1 * ex(-12.260f /
					T9);
	//
	//.......B12(p,a)Be9................(Wagoner 1969)
	f[48] = 2.01e+11f * T9m23 * ex(-12.12f / T913);
	//
	//80--------ALPHA, PHOTON REACTIONS----------------------------------------------
	//
	//.......Li6(a,g)B10................(Caughlan-Fowler 1988)
	f[49] = 4.06e+6f * T9m23 * ex(-18.790f / T913 - fem::pow2((T9 /
					1.326f))) * (1. + .022f * T913 + 1.54f * T923 + .239f * T9 +
				2.20f * T943 + .869f * T953) + 1.91e+3f * T9m32 * ex(-3.484f /
					T9) + 1.01e+4f * T9m1 * ex(-7.269f / T9);
	//
	//.......Li7(a,g)B11................(Caughlan-Fowler 1988)
	f[50] = 3.55e+7f * T9m23 * ex(-19.161f / T913 - fem::pow2((T9 /
					4.195f))) * (1. + .022f * T913 + .775f * T923 + .118f * T9 +
				.884f * T943 + .342f * T953) + 3.33e+2f * T9m32 * ex(-2.977f /
					T9) + 4.10e+4f * T9m1 * ex(-6.227f / T9);
	//
	//.......Be7(a,g)C11................(Caughlan-Fowler 1988)
	f[51] = 8.45e+7f * T9m23 * ex(-23.212f / T913 - fem::pow2((T9 /
					4.769f))) * (1. + .018f * T913 + .488f * T923 + .061f * T9 +
				.296f * T943 + .095f * T953) + 1.25e+4f * T9m32 * ex(-6.510f /
					T9) + 1.29e+5f * T9m54 * ex(-10.039f / T9);
	//
	//90--------ALPHA, PROTON REACTIONS----------------------------------------------
	//
	//.......B8(a,p)C11.................(Wagoner 1969)
	f[52] = 1.08e+15f * T9m23 * ex(-27.36f / T913);
	//
	//100-------ALPHA, NEUTRON REACTIONS---------------------------------------------
	//
	//.......Li8(a,n)B11................(Malaney-Fowler 1989)
	f[53] = 8.62e+13f * T9a56 * T9m32 * ex(-19.461f / T9a13);
	//
	//.......Be9(a,n)C12................(Caughlan-Fowler 1988)
	f[54] = 4.62e+13f * T9m23 * ex(-23.870f / T913 - fem::pow2((T9 /
					.049f))) * (1. + .017f * T913 + 8.57f * T923 + 1.05f * T9 + 74.51f *
				T943 + 23.15f * T953) + 7.34e-5f * T9m32 * ex(-1.184f / T9) +
					2.27e-1f * T9m32 * ex(-1.834f / T9) + 1.26e+5f * T9m32 * ex(-4.179f /
							T9) + 2.40e+8f * ex(-12.732f / T9);
	//
	//110-------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------------
	//
	//.......Be9(d,n)B10................(original Wagoner code)
	f[55] = 7.16e+8f * T9m23 * ex(6.44f - 12.6f / T913);
	//
	//.......B10(d,p)B11................(original Wagoner code)
	f[56] = 9.53e+8f * T9m23 * ex(7.30f - 14.8f / T913);
	//
	//.......B11(d,n)C12................(original Wagoner code)
	f[57] = 1.41e+9f * T9m23 * ex(7.40f - 14.8f / T913);
	//
	//120-------THREE PARTICLE REACTIONS---------------------------------------------
	//
	//.......He4(an,g)Be9...............(Caughlan-Fowler 1988)
	f[58] = (2.59e-6f / ((1. + .344f * T9) * fem::pow2(T9))) * ex(-1.062f / T9);
	//
	//.......He4(2a,g)C12...............(Caughlan-Fowler 1988)
	f[59] = 2.79e-8f * T9m32 * T9m32 * ex(-4.4027f / T9) + 1.35e-8f *
		T9m32 * ex(-24.811f / T9);
	//
	//.......Li8(p,na)He4...............(original Wagoner code)
	f[60] = 8.65e+9f * T9m23 * ex(-8.52f / T913 - fem::pow2((T9 /
					2.53f))) + 2.31e+9f * T9m32 * ex(-4.64f / T9);
	//
	//.......B8(n,pa)He4................(original Wagoner code)
	f[61] = 4.02e+8f;
	//
	//.......Be9(p,da)He4...............(Caughlan-Fowler 1988)
	f[62] = 2.11e+11f * T9m23 * ex(-10.359f / T913 - fem::pow2((T9 /
					.520f))) * (1. + .040f * T913 + 1.09f * T923 + .307f * T9 +
				3.21f * T943 + 2.30f * T953) + 5.79e+8f * T9m1 * ex(-3.046f /
					T9) + 8.50e+8f * T9m34 * ex(-5.800f / T9);
	//
	//.......B11(p,2a)He4...............(Caughlan-Fowler 1988)
	f[63] = 2.20e+12f * T9m23 * ex(-12.095f / T913 - fem::pow2((T9 /
					1.644f))) * (1. + .034f * T913 + .140f * T923 + .034f * T9 +
				.190f * T943 + .116f * T953) + 4.03e+6f * T9m32 * ex(-1.734f /
					T9) + 6.73e+9f * T9m32 * ex(-6.262f / T9) + 3.88e+9f * T9m1 * ex(
						-14.154f / T9);
	//
	//.......C11(n,2a)He4...............(Wagoner 1969)
	f[64] = 1.58e+8f;
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
void common::rate4()
{
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
	//T9**(1/3)
	double T913 = pow(T9, (.33333333f));
	//T9**(2/3)
	double T923 = T913 * T913;
	//T9**(4/3)
	double T943 = T923 * T923;
	//T9**(5/3)
	double T953 = T9 * T923;
	//T9**(1/2)
	double T912 = sqrt(T9);
	//T9**(3/2)
	double T932 = T9 * T912;
	//T9**(3/5)
	double T935 = pow(T9, (.6f));
	//T9**(6/5)
	double T965 = pow(T9, (1.2f));
	//T9**(3/8)
	double T938 = pow(T9, (.375f));
	//T9**(1/3)
	double T9m13 = 1.0f / T913;
	//T9**(-2/3)
	double T9m23 = 1.0f / T923;
	//T9**(-3/2)
	double T9m32 = 1.0f / T932;
	//T9**(-6/5)
	double T9m65 = 1.0f / T965;
	//For reaction 82.
	double T9a = T9 / (1. + 4.78e-2f * T9 + 7.56e-3f * T953 / pow((
					1.f + 4.78e-2f * T9), (2.f / 3.f)));
	//T9a**(1/3)
	double T9a13 = pow(T9a, (.33333333f));
	//T9a**(5/6)
	double T9a56 = pow(T9a, (.83333333f));
	//For reaction 84.
	double T9b = T9 / (1. + 7.76e-2f * T9 + 2.64e-2f * T953 / pow((
					1.f + 7.76e-2f * T9), (2.f / 3.f)));
	//T9b**(1/3)
	double T9b13 = pow(T9b, (.33333333f));
	//T9b**(5/6)
	double T9b56 = pow(T9b, (.83333333f));
	//
	//20--------NEUTRON, PHOTON REACTIONS--------------------------------------------
	//
	//.......C12(n,g)C13................(Wagoner 1969)
	f[65] = 4.50e+2f;
	//
	//.......C13(n,g)C14................(Wagoner 1969)
	f[66] = 1.19e+2f + 2.38e+5f * T9m32 * ex(-1.67f / T9);
	//
	//.......N14(n,g)N15................(Wagoner 1969)
	f[67] = 9.94e+3f;
	//
	//30--------NEUTRON, PROTON REACTIONS--------------------------------------------
	//
	//.......N13(n,p)C13................(Caughlan-Fowler 1988)
	f[68] = 1.88e+8f * (1. - .167f * T912 + .037f * T9);
	//
	//.......N14(n,p)C14................(Caughlan-Fowler 1988)
	f[69] = 2.39e+5f * (1.f + .361f * T912 + .502f * T9) + 1.112e+8f /
		T912 * ex(-4.983f / T9);
	//
	//.......O15(n,p)N15................(Caughlan-Fowler 1988)
	f[70] = 3.50e+8f * (1.f + .452f * T912 - .191f * T9);
	//
	//40--------NEUTRON, ALPHA REACTIONS---------------------------------------------
	//
	//.......O15(n,a)C12................(Caughlan-Fowler 1988)
	f[71] = 3.50e+7f * (1.f + .188f * T912 + .015f * T9);
	//
	//50--------PROTON, PHOTON REACTIONS---------------------------------------------
	//
	//.......C12(p,g)N13................(Caughlan-Fowler 1988)
	f[72] = 2.04e+7f * T9m23 * ex(-13.690f / T913 - fem::pow2((T9 /
					1.500f))) * (1.f + .030f * T913 + 1.19f * T923 + .254f * T9 +
				2.06f * T943 + 1.12f * T953) + 1.08e+5f * T9m32 * ex(-4.925f /
					T9) + 2.15e+5f * T9m32 * ex(-18.179f / T9);
	//
	//.......C13(p,g)N14................(Caughlan-Fowler 1988)
	f[73] = 8.01e+7f * T9m23 * ex(-13.717f / T913 - fem::pow2((T9 /
					2.000f))) * (1.f + .030f * T913 + .958f * T923 + .204f * T9 +
				1.39f * T943 + .753f * T953) + 1.21e+6f * T9m65 * ex(-5.701f /
					T9);
	//
	//.......C14(p,g)N15................(Caughlan-Fowler 1988)
	f[74] = 6.80e+6f * T9m23 * ex(-13.741f / T913 - fem::pow2((T9 /
					5.721f))) * (1.f + .030f * T913 + .503f * T923 + .107f * T9 +
				.213f * T943 + .115f * T953) + 5.36e+3f * T9m32 * ex(-3.811f /
					T9) + 9.82e+4f * T9m13 * ex(-4.739f / T9);
	//
	//.......N13(p,g)O14................(Caughlan-Fowler 1988)
	f[75] = 4.04e+7f * T9m23 * ex(-15.202f / T913 - fem::pow2((T9 /
					1.191f))) * (1.f + .027f * T913 - .803f * T923 - .154f * T9 +
				5.00f * T943 + 2.44f * T953) + 2.43e+5f * T9m32 * ex(-6.348f /
					T9);
	//
	//.......N14(p,g)O15................(Caughlan-Fowler 1988)
	f[76] = 4.90e+7f * T9m23 * ex(-15.228f / T913 - fem::pow2((T9 /
					3.294f))) * (1.f + .027f * T913 - .778f * T923 - .149f * T9 + .261f *
				T943 + .127f * T953) + 2.37e+3f * T9m32 * ex(-3.011f / T9) +
					2.19e+4f * ex(-12.530f / T9);
	//
	//.......N15(p,g)O16................(Caughlan-Fowler 1988)
	f[77] = 9.78e+8f * T9m23 * ex(-15.251f / T913 - fem::pow2((T9 /
					.450f))) * (1.f + .027f * T913 + .219f * T923 + .042f * T9 +
				6.83f * T943 + 3.32f * T953) + 1.11e+4f * T9m32 * ex(-3.328f /
					T9) + 1.49e+4f * T9m32 * ex(-4.665f / T9) + 3.80e+6f * T9m32 * ex(
						-11.048f / T9);
	//
	//60--------PROTON, ALPHA REACTIONS----------------------------------------------
	//
	//.......N15(p,a)C12................(Caughlan-Fowler 1988)
	f[78] = 1.08e+12f * T9m23 * ex(-15.251f / T913 - fem::pow2((T9 /
					.522f))) * (1.f + .027f * T913 + 2.62f * T923 + .501f * T9 +
				5.36f * T943 + 2.60f * T953) + 1.19e+8f * T9m32 * ex(-3.676f /
					T9) + 5.41e+8f / T912 * ex(-8.926f / T9) + 4.72e+7f * T9m32 * ex(
						-7.721f / T9) + 2.20e+8f * T9m32 * ex(-11.418f / T9);
	//
	//70--------ALPHA, PHOTON REACTIONS----------------------------------------------
	//
	//.......C12(a,g)O16................(Caughlan-Fowler 1988)
	f[79] = 1.04e+8f / fem::pow2(T9) * ex(-32.120f / T913 - fem::pow2((T9 /
					3.496f))) / fem::pow2((1.f + .0489f * T9m23)) + 1.76e+8f / fem::pow2((
					T9)) / fem::pow2((1.f + .2654f * T9m23)) * ex(-32.120f / T913) +
					1.25e+3f * T9m32 * ex(-27.499f / T9) + 1.43e-2f * pow((T9), 5) *
					ex(-15.541f / T9);
	//
	//80--------ALPHA, PROTON REACTIONS----------------------------------------------
	//
	//.......B10(a,p)C13................(Wagoner 1969)
	f[80] = 9.60e+14f * T9m23 * ex(-27.99f / T913);
	//
	//.......B11(a,p)C14................(Caughlan-Fowler 1988)
	f[81] = 5.37e+11f * T9m23 * ex(-28.234f / T913 - fem::pow2((T9 /
					0.347f))) * (1.f + .015f * T913 + 5.575f * T923 + .576f * T9 +
				15.888f * T943 + 4.174f * T953) + 5.44e-3f * T9m32 * ex(-2.827f /
					T9) + 3.36e+2f * T9m32 * ex(-5.178f / T9) + 5.32e+6f / T938 * ex(
						-11.617f / T9);
	//
	//.......C11(a,p)N14................(Caughlan-Fowler 1988)
	f[82] = 7.15e+15f * T9a56 * T9m32 * ex(-31.883f / T9a13);
	//
	//.......N12(a,p)O15................(Caughlan-Fowler 1988)
	f[83] = 5.59e+16f * T9m23 * ex(-35.60f / T913);
	//
	//.......N13(a,p)O16................(Caughlan-Fowler 1988)
	f[84] = 3.23e+17f * T9b56 * T9m32 * ex(-35.829f / T9b13);
	//
	//90--------ALPHA, NEUTRON REACTIONS---------------------------------------------
	//
	//.......B10(a,n)N13................(Caughlan-Fowler 1988)
	f[85] = 1.20e+13f * T9m23 * ex(-27.989f / T913 - fem::pow2((T9 / 9.589f)));
	//
	//.......B11(a,n)N14................(Caughlan-Fowler 1988)
	f[86] = 6.97e+12f * T9m23 * ex(-28.234f / T913 - fem::pow2((T9 /
					0.140f))) * (1.f + .015f * T913 + 8.115f * T923 + .838f * T9 +
				39.804f * T943 + 10.456f * T953) + 1.79e+0f * T9m32 * ex(
					-2.827f / T9) + 1.71e+3f * T9m32 * ex(-5.178f / T9) + 4.49e+6f *
				T935 * ex(-8.596f / T9);
	//
	//.......B12(a,n)N15................(Wagoner 1969)
	f[87] = 3.04e+15f * T9m23 * ex(-28.45f / T913);
	//
	//.......C13(a,n)O16................(Caughlan-Fowler 1988)
	f[88] = 6.77e+15f * T9m23 * ex(-32.329f / T913 - fem::pow2((T9 /
					1.284f))) * (1.f + .013f * T913 + 2.04f * T923 + .184f * T9) +
		3.82e+5f * T9m32 * ex(-9.373f / T9) + 1.41e+6f * T9m32 * ex(
				-11.873f / T9) + 2.00e+9f * T9m32 * ex(-20.409f / T9) +
		2.92e+9f * T9m32 * ex(-29.283f / T9);
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
void common::derivs(
		int const& loop)
{
	const double pi = constants::pi<double>();
	int i = 0;
	double dphdT9 = 0;
	double dphdln = 0;
	double dphdzy = 0;
	double bar = 0;
	double dlndT9 = 0;
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
	//(1/T9)*d(T9)/d(t).
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
	//d(phi e)/d(T9).
	//d(phi e)/d(h).
	//d(phi e)/d(sumzy).
	//(1/h)*d(h)/d(T9).
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
	rnb = hv * T9 * T9 * T9 / rhob0;
	//..........VARIOUS THERMODYNAMIC QUANTITIES.
	therm();
	//Expansion rate.
	hubcst = sqrt((8.f / 3.f) * pi * g * (thm(10)) + (cosmo / 3.f));
	//Baryon mass density.
	rhob = thm(9);
	//..........COMPUTE REACTION RATE COEFFICIENTS.
	rate1(T9);
	//Run network selection.
	switch (irun) {
		case 1: goto statement_100;
		case 2: goto statement_110;
		case 3: goto statement_120;
		default: break;
	}
statement_100:
	//Forward rate for all of reactions.
	rate4();
statement_110:
	//Forward rate for reactions with A <
	rate3();
statement_120:
	//Forward rate for reactions with A <
	rate2();
	//..........SOLVE COUPLED DIFFERENTIAL EQUATIONS.
	sol(loop);
	//Abort in case matrix not invertible.
	if (mbad > 0) {
		return;
	}
	//
	//20--------COMPUTE DERIVATIVES FOR TEMPERATURE, hv, AND CHEMICAL POTENTIAL------
	//
	//..........INITIALIZE SUMS TO ZERO.
	double sumy = 0;
	double sumzy = 0;
	double sumdy = 0;
	double summdy = 0;
	double sumzdy = 0;
	//..........ACCUMULATE TO GET SUM.
	FEM_DO_SAFE(i, 1, isize) {
		sumy += y[i]; 					/// Sum of abundance.
		sumzy += zm[i] * y[i]; 			/// Sum of charge*abundance.
		sumdy += dydt[i]; 				/// Sum of abundance flow.
		summdy += dm[i] * dydt[i]; 		/// Sum of (mass excess)*(abundanc
		sumzdy += zm[i] * dydt[i]; 		/// Sum of (charge)*(abundance flo
	}
	//..........CHANGES IN TEMPERATURE, hv, AND CHEMICAL POTENTIAL.
	dphdT9 = thm(12) * (-1.070e-4f * hv * sumzy / T9 - thm(11));
	dphdln = -thm(12) * 3.568e-5f * hv * sumzy;
	dphdzy = thm(12) * 3.568e-5f * hv;
	bar = 9.25e-5f * T9 * sumy + 1.388e-4f * T9 * sumdy 
		/ (3.f * hubcst) + summdy / (3.f * hubcst);
	//(Ref 1)
	dlndT9 = -(thm(2) + thm(5) + thm(6) * dphdT9 + thm(9) * 1.388e-4f * sumy) 
		/ (thm(1) + thm(3) + thm(4) + thm(7) + thm(9) * bar + thm(6) 
				* (dphdln + dphdzy * sumzdy / (3.f * hubcst)));
	dT9 = (3.f * hubcst) / dlndT9;
	dlT9dt = dT9 / T9;
	//(Ref 2)
	dhv = -hv * 3*(hubcst + dlT9dt);
	//(Ref 3)
	dphie = dphdT9 * dT9 + dphdln * (3.f * hubcst) + dphdzy * sumzdy;
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
void common::accum()
{
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
	FEM_DO_SAFE(i, 1, isize) {
		xout(it, i) = y[i] / y[2];
	}
	//xout(it, 2) = y(2) * am(2); 		/// Exception for proton.
	xout(it, 2) = y[2] * am[2]; 			/// Exception for proton.
	//xout(it, 6) = y(6) * am(6); 		/// Exception for helium.
	xout(it, 6) = y[6] * am[6]; 			/// Exception for helium.
	//..........SUM UP ABUNDANCES OF HEAVY NUCLIDES.
	//Li8 to O16.
	xout(it, 10) += xout(it,11) + xout(it,12) + xout(it,13) + xout(it,14) 
		+ xout(it,15) + xout(it,16) + xout(it,17) + xout(it,18) 
		+ xout(it,19) + xout(it,20) + xout(it,21) + xout(it,22) 
		+ xout(it,23) + xout(it,24) + xout(it,25) + xout(it,26);
	//..........RELABEL TEMPERATURE, TIME, THERMODYNAMIC VARIABLES, ETC.
	//Temperature.
	T9out(it) = T9;
	//Time.
	tout(it) = t;
	//rho photon.
	thmout(it, 1) = thm(1);
	//rho electron.
	thmout(it, 2) = thm(4);
	//rho neutrino.
	thmout(it, 3) = thm(8);
	//rho baryon.
	thmout(it, 4) = thm(9);
	//Chemical potential.
	thmout(it, 5) = phie;
	//rho total.
	thmout(it, 6) = thm(10);
	//Time step.
	dtout(it) = dt;
	//Baryon to photon ratio.
	etaout(it) = hv / (3.3683e+4f);
	//Expansion rate.
	hubout(it) = hubcst;
	//
	//20--------INDICATE TERMINATION OF ACCUMULATION IF APPROPRIATE------------------
	//
	if ((it == itmax) || (ip < inc)) {
		ltime = 1;
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
void common::driver()
{
	int mvar = 0;
	int loop = 0;
	const double cl = 1.e-16f;
	double dtmin = 0;
	int i = 0;
	double dtl = 0;
	const int nvar = 29;
	double v[nvar+1];
	double v0[nvar+1];
	double dvdt[nvar+1];
	double dvdt0[nvar+1];
	FEM_DO_SAFE(i, 1, nvar) {
		v[i] = 0;
		v0[i] = 0;
		dvdt[i] = 0;
		dvdt0[i] = 0;
	}
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
	//(1/T9)*d(T9)/d(t).
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
	start();
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
	derivs(loop);
	//Time = 1st R-K loop.
	itime = 4;
	//Check interface subroutine.
	check();
	//..........ACCUMULATE.
	//Low temp.
	//Small dt.
	//Enough iterations.
	//if ((cmn.T9 <= cmn.T9f) || (dt < abs(cl / dlT9dt)) || (ip == inc)) {
	if ((T9 <= T9f) || (dt < abs(cl / dlT9dt)) || (ip == inc)) {
		accum();
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

		// test -->
		FEM_DO_SAFE(i, 1, isize)
			std::cout<<"dydt["<<i<<"]:"<<dydt[i]<<" y["<<i<<"]:"<<y[i]<<"\n";
		// <-- test

		//Trial value for minimum time step (R
		dtmin = abs(1. / dlT9dt) * ct;
		//Go through all abundance changes.
		FEM_DO_SAFE(i, 1, isize) {
			if ((dydt[i] != 0) && (y[i] > ytmin)) {
				dtl = abs(y[i] / dydt[i]) * cy 
					* (1.f + fem::pow2(( fem::alog10(y[i]) / fem::alog10(ytmin))));  /// (Ref 2).
				if (dtl < dtmin) { 												/// Find smallest time st
					dtmin = dtl;
				}
			}
		}
		//Limit change in time
		if (dtmin > 1.5f * dt) {
			dtmin = 1.5f * dt;
		}
		dt = dtmin; 					/// Set new time step.
	}
	t += dt; 							/// Increment time.
	//..........STORE AND INCREMENT VALUES (Ref 3).
	/*
	   FEM_DO_SAFE(i, 1, mvar) {
	   v0(i) = v(i);
	   dvdt0(i) = dvdt(i);
	   v(i) = v0(i) + dvdt0(i) * dt;
	//Set at minimum
	if ((i >= 4) && (v(i) < ytmin)) {
	v(i) = ytmin;
	}
	}
	 */
	FEM_DO_SAFE(i, 1, mvar) {
		v0[i] = v[i];
		dvdt0[i] = dvdt[i];
		v[i] = v0[i] + dvdt0[i] * dt;
		//Set at minimum
		if ((i >= 4) && (v[i] < ytmin)) {
			v[i] = ytmin;
		}
	}
	//
	//30--------LOOP TWO-------------------------------------------------------------
	//
	//Step up loop counter.
	loop = 2;
	//..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
	derivs(loop);
	//Time = 2nd R-K loop.
	itime = 7;
	//Check interface subroutine.
	check();
	//..........INCREMENT VALUES.
	FEM_DO_SAFE(i, 1, mvar) {
		v[i] = v0[i] + .5f * (dvdt[i] + dvdt0[i]) * dt;
		//Set at minimum
		if ((i >= 4) && (v[i] < ytmin)) {
			v[i] = ytmin;
		}
	}
	goto statement_200;
	//
	//----------REFERENCES-----------------------------------------------------------
	//     1)  Constraint on dt from the requirement that
	//                (d(T9)/dt)*(dt/T9) < ct
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
void common::qvary(int index, double value)
{
	//----------EQUIVALENCE VARIABLE.
	//     REAL    qvary(7)             !Array set equal to c, cosmo, and xi.
	//
	//----------EQUIVALENCE STATEMENTS.
	//     EQUIVALENCE (qvary(1),c(1)), (qvary(4),cosmo), (qvary(5),xi(1))
	//
	if (index >= 1 && index <= 3)
		c[index] = value;
	else if (index == 4)
		cosmo = value;
	else if (index >= 5 && index <= 7)
		xi[index] = value;
	else
	{
		std::cerr << "index out of bounds." << std::endl;
		exit(1);
	}
}


//
//========================IDENTIFICATION DIVISION================================
//
void common::run()
{
	common_read read(*this);
	common_write write(*this);
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
	arr_1d<3, double> rnum1(fem::fill0);
	arr_1d<3, double> rnum2(fem::fill0);
	arr_1d<3, double> rnum3(fem::fill0);
	arr_1d<3, int> inum(fem::fill0);
	fem::str<1> lchose;
	int l = 0;
	arr_1d<3, int> lnum(fem::fill0);
	int lnumb1 = 0;
	double rnumb1 = 0;
	int lnumb2 = 0;
	double rnumb2 = 0;
	int lnumb3 = 0;
	double rnumb3 = 0;
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
	//c[1] is variation of gravitational constant.
	//c[2] is neutron lifetime (sec).
	//c[3] is number of neutrino species.
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
	//..........RETURN FROM LOOP.
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
	check(); 					/// Call interface subroutine before computation.
	driver(); 					/// Do nucleosynthesis computation.
	itime = 8;
	check(); 					/// Call interface subroutine after computation.
	write(iw, format_2202);
	read(ir, star);
	goto statement_300;
	//
	//23--------DO MULTIPLE RUNS SECTION---------------------------------------------
	//
	//..........GET NUMBER OF LOOPS.
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
				} //((inum(i) < 1) or (inum(i) > 8))
			} //(i > jnum)
		} //i = 1,3
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
		//..........COMPUTE NUMBER OF RUNS FOR EACH LOOP.
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
					qvary(inum(1)-1, rnumb1);
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
						qvary(inum(2)-1, rnumb2);
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
							qvary(inum(3)-1, rnumb3);
						}
					}
					itime = 3;
					//Check interface subroutine before computation.
					check();
					//Do nucleosynthesis computation.
					driver();
					itime = 8;
					//Check interface subroutine after computation.
					check();
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
void common::output()
{
	common_read read(*this);
	common_write write(*this);
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
	//c[1] is variation of gravitational c
	//c[2] is neutron lifetime (sec).
	//c[3] is number of neutrino species.
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
	//..........RETURN FROM LOOP.
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
	//        T9out(j) = T9out(j)*.08617
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
		cy, ct, T9i, T9f, ytmin;
	write(2,
			"(' Model parameters:',/,'   g = ',f5.2,'/  tau = ',f6.2,'/  # nu = ',"
			"f5.2,'/  lambda = ',1p,e10.3,'/  xi-e = ',e10.3,'/  xi-m = ',e10.3,"
			"'/  xi-t = ',e10.3,/)"),
		c[1], c[2], c[3], cosmo, xi[1], xi[2], xi[3];
	//..........PRINT HEADINGS, ABUNDANCES FOR NEUTRON TO LI8.
	write(2,
			"(4x,'Temp',8x,'N/H',10x,'P',10x,'D/H',9x,'T/H',8x,'He3/H',8x,'He4',8x,"
			"'Li6/H',7x,'Li7/H',7x,'Be7/H',6x,'Li8/H&up',/,132('-'))");
	FEM_DO_SAFE(j, 1, it) {
		{
			write_loop wloop(*this, 2, "(1p,e10.3,1p,10e12.3)");
			wloop, T9out(j);
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
			write_loop wloop(*this, 2, "(1p,e10.3,9e12.3)");
			wloop, T9out(j), tout(j);
			FEM_DO_SAFE(i, 1, 5) {
				wloop, thmout(j, i);
			}
			wloop, dtout(j), etaout(j), hubout(j);
		}
	}
	write(2, format_2014);
	//Output file requested.
	outfile = true;
	write(iw,
			"(' ','Output file requested - Press <RETURN> to ','continue: ',$)");
	read(ir, star);
	goto statement_500;
	//
	//30--------REQUEST OUTPUT ON SCREEN SECTION-------------------------------------
	//
	//..........RETURN FROM LOOP.
statement_300:
	//      DO j = 1,it                  !Temperature in MeV.
	//        T9out(j) = T9out(j)*.08617
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
	write(iw, format_3100), cy, ct, T9i, T9f, ytmin;
	write(iw, format_3102), c[1], c[2], c[3], cosmo, xi[1], xi[2], xi[3];
	//..........PRINT HEADINGS, ABUNDANCES FOR D,T,HE3,HE4,LI7.
	write(iw,
			"(4x,'Temp',8x,'D/H',9x,'T/H',8x,'He3/H',8x,'He4',8x,'Li7/H',/,"
			"80('-'))");
	FEM_DO_SAFE(j, 1, it) {
		write_loop wloop(*this, iw, format_3106);
		wloop, T9out(j);
		FEM_DO_SAFE(i, 3, 6) {
			wloop, xout(j, i);
		}
		wloop, xout(j, 8);
	}
	write(iw, format_2014);
	write(iw, format_3108);
	read(ir, star);
	goto statement_360;
	//Display n,p,li6,be7,li8&up.
statement_320:
	//..........PRINT CAPTION.
	write(iw, format_2014);
	write(iw, format_3100), cy, ct, T9i, T9f, ytmin;
	write(iw, format_3102), c[1], c[2], c[3], cosmo, xi[1], xi[2], xi[3];
	//..........PRINT HEADINGS, ABUNDANCES FOR N,P,LI6,BE7,LI8&UP.
	write(iw,
			"(4x,'Temp',8x,'N/H',10x,'P',9x,'Li6/H',7x,'Be7/H',6x,'Li8/H&up',/,"
			"80('-'))");
	FEM_DO_SAFE(j, 1, it) {
		{
			write_loop wloop(iw, format_3106);
			wloop, T9out(j);
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
	write(iw, format_3100), cy, ct, T9i, T9f, ytmin;
	write(iw, format_3102), c[1], c[2], c[3], cosmo, xi[1], xi[2], xi[3];
	//..........PRINT ENERGY DENSITIES.
	write(iw,
			"(4x,'Temp',8x,'rhog',8x,'rhoe',7x,'rhone',8x,'rhob',/,80('-'))");
	FEM_DO_SAFE(j, 1, it) {
		{
			write_loop wloop(*this, iw, "(1p,e10.3,4e12.3)");
			wloop, T9out(j);
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
	write(iw, format_3100), cy, ct, T9i, T9f, ytmin;
	write(iw, format_3102), c[1], c[2], c[3], cosmo, xi[1], xi[2], xi[3];
	//..........PRINT THERMODYNAMIC QUANTITIES.
	write(iw,
			"(4x,'Temp',8x,'time',8x,'phie',9x,'dt',9x,'eta',10x,'H',/,80('-'))");
	FEM_DO_SAFE(j, 1, it) {
		write(iw, "(1p,e10.3,5e12.3)"), T9out(j), tout(j), thmout(j, 5),
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
common::common() :
	fem::common(),
	common_compr0(),
	common_compr(),
	common_varpr0(),
	common_varpr(),
	common_modpr0(),
	common_modpr(),
	//common_recpr0(),
	common_recpr(),
	common_evolp1(),
	common_evolp2(),
	common_evolp3(),
	common_ttime(),
	common_thermcb(),
	common_endens(),
	//common_nucdat(),
	//common_besselcb(),
	//common_kays(),
	common_flags(),
	//common_checkcb(),
	common_outdat(),
	common_nupar(),
	common_runopt(),
	common_outopt(),
	common_tcheck()
{
	cy0 = .300f;
	ct0 = .030f;
	T9i0 = 1.00e+02f;
	T9f0 = 1.00e-02f;
	ytmin0 = 1.00e-25f;			// TODO make smaller with double
	inc0 = 30;
	c0[1] = 1.00;
	c0[2] = 885.7;
	c0[3] = 3.0;
	cosmo0 = 0.00f;
	xi0[1] = 0;
	xi0[2] = 0;
	xi0[3] = 0;
	dt0 = 1.00e-04f;
	eta0 = 3.162e-10f;
	for (int i = 0; i < nnuc+1; i++)
	for (int j = 0; j < nnuc+1; j++)
		if (i == 0 or j == 0)
			a[i][j] = NOT_USED;
		else
			a[i][j] = 0;

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
	//c0[1] is default variation of grav c
	//c0[2] is default neutron half-life.
	//c0[3] is default number of neutrinos
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
void common::program_new123()
{
	common_read read(*this);
	common_write write(*this);
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
	//Default T9i.
	//Default T9f.
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
	//c[1] is variation of gravitational c
	//c[2] is neutron lifetime (sec).
	//c[3] is number of neutrino species.
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
	io.open(2, "new123.dat")
		.status("new");
	//Time = beginning of program.
	itime = 1;
	//Check interface subroutine.
	check();
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
	FEM_DO_SAFE(i, 1, nrec) {
		//..........READ IN REACTION PARAMETERS.
		/*
		   iform(i) = reacpr(i, 2); 		/// Reaction type.
		   ii(i) = reacpr(i, 3); 			/// Incoming nuclide type.
		   jj(i) = reacpr(i, 4); 			/// Incoming nuclide type.
		   kk(i) = reacpr(i, 5); 			/// Outgoing nuclide type.
		   ll(i) = reacpr(i, 6); 			/// Outgoing nuclide type.
		   rev(i) = reacpr(i, 7); 			/// Reverse reaction coefficient.
		   q9(i) = reacpr(i, 8); 			/// Energy released.
		 */
		iform(i) = reacpr[i-1][2-1]; 		/// Reaction type.
		ii(i) = reacpr[i-1][3-1]; 			/// Incoming nuclide type.
		jj(i) = reacpr[i-1][4-1]; 			/// Incoming nuclide type.
		kk(i) = reacpr[i-1][5-1]; 			/// Outgoing nuclide type.
		ll(i) = reacpr[i-1][6-1]; 			/// Outgoing nuclide type.
		rev(i) = reacpr[i-1][7-1]; 			/// Reverse reaction coefficient.
		q9(i) = reacpr[i-1][8-1]; 			/// Energy released.
		//..........INITIALIZE REACTION RATES.
		f[i] = 0; 							/// Forward rate coeff.
		r[i] = 0; 							/// Reverse rate coeff.
		//..........SET RUN OPTIONS TO DEFAULT.
	}
	irun = 1; 							/// Do full run.
	isize = nnuc; 					/// Use all 26 nuclides.
	jsize = nrec; 					/// Use all 88 reactions.
	//..........SET OUTPUT OPTION TO DEFAULT.
	nout = 0; 							/// No output requests.
	outfile = false; 						/// Output file not used.
	//..........SET VALUES TO DEFAULT.
	cy = cy0; 						/// Time step limiting constant on abundance.
	ct = ct0; 						/// Time step limiting constant on temperature.
	T9i = T9i0; 						/// Initial temperature.
	T9f = T9f0; 						/// Final temperature.
	ytmin = ytmin0; 					/// Smallest abundances allowed.
	inc = inc0; 						/// Accumulation increment.
	c[1] = c0[1]; 					/// Variation of gravitational constant.
	c[2] = c0[2]; 					/// Neutron lifetime.
	std::cout << "ntau:"<<c[2]<<std::endl; 
	c[3] = c0[3]; 					/// Number of neutrino species.
	cosmo = cosmo0; 					/// Cosmological constant.
	xi[1] = xi0[1]; 					/// Electron degeneracy parameter.
	xi[2] = xi0[2]; 					/// Muon degeneracy parameter.
	xi[3] = xi0[3]; 					/// Tau degeneracy parameter.
	dt1 = dt0; 						/// Initial time step.
	eta1 = eta0; 						/// Baryon-to-photon ratio.
	//..........ACCEPT RETURN TO CONTINUE.
	read(ir, star); 							/// Pause.
	//
	//30--------PRINT MENU AND AWAIT RESPONSE----------------------------------------
	//
	//..........RETURN FROM LOOP.
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
	help();
	goto statement_500;
	//Set computation parameters section.
statement_420:
	setcom();
	goto statement_500;
	//Set model parameters section.
statement_430:
	setmod();
	goto statement_500;
	//Run section.
statement_440:
	//Time = beginning of run section.
	itime = 2;
	//Check interface subroutine.
	check();
	run();
	//Time = end of run section.
	itime = 9;
	//Check interface subroutine.
	check();
	goto statement_500;
	//Output section.
statement_450:
	output();
	goto statement_500;
	//Exit section.
statement_460:
	if (outfile) {
		//Close output file.
		io.close(2).status("keep");
	}
	else {
		//File not used - dispose.
		io.close(2).status("delete");
	}
	//Time = end of program.
	itime = 10;
	//Check interface subroutine.
	check();
	//FEM_STOP(0);
	//
	//50---------GO BACK TO MENU-----------------------------------------------------
	//
statement_500:
	goto statement_300;
	//
}

//} // namespace bbn

	int
main(
		int argc,
		char const* argv[])
{
	common cmn(argc, argv);
	cmn.program_new123();
	//return fem::main_with_catch(
	//		argc, argv,
	//		program_new123);
}

