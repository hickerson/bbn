CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Changes (to run inder unix f77):
C     -------------------------------
C Programme renamed from nuc123.for to new123.f
C ir=1 -> ir=5 ... input unit number
C iw=1 -> iw=6 ... output unit number 
C COMMON /check/ itime -> COMMON /checkcb/ itime
C COMMON /time/        -> COMMON /ttime/
C output nuc123.dat    -> new123.dat
C
C========================IDENTIFICATION DIVISION================================

      PROGRAM new123 

C----------LINKAGES.
C     CALLED BY - none
C     CALLS     - [subroutine] help, setcom, setmod, run, output

C----------REMARKS.
C     Control program -
C       Offers user the main menu and channels through to various options.
C     Implementation -
C       To run this program, new123.f must be linked with nuccom.f 
C       (containing the computation subroutines), nucrat.f (with the
C       reaction rates), and newint.f (with an interface subroutine).
C       This program has been written to be compatible with
C       ANSI FORTRAN-77 with the exception of the END DO statement
C       used to limit the number of statement labels.
C       This code was modified on the DEC/AXP system.
C     Notes -
C       The program utilizes Wagoner's code as the core of the computational
C       routines.
C     Documentation -
C       Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
C       Kellogg Radiation Lab preprint OAP-714.
C     Copy -
C       Version 4.1 (December 1991)

C----------PARAMETERS.
      PARAMETER (ir=5)             !Input unit number (previous value = 1).
      PARAMETER (iw=6)             !Output unit number (previous value = 1).
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C----------COMMON AREAS.
      COMMON /recpr0/ reacpr                         !Reaction parameter values.
      COMMON /recpr/  iform,ii,jj,kk,ll,rev,q9       !Reaction parameter names.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /compr0/ cy0,ct0,t9i0,t9f0,ytmin0,inc0  !Default comp parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /modpr0/ c0,cosmo0,xi0                  !Default model parameters.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /varpr0/ dt0,eta0                       !Default variationl params.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.
      COMMON /checkcb/  itime                        !Computation location.
      COMMON /runopt/ irun,isize,jsize               !Run options.
      COMMON /outopt/ nout,outfile                   !Output option.


C==========================DECLARATION DIVISION=================================

C----------REACTION PARAMETERS FROM BLOCK DATA.
      REAL    reacpr(nrec,8)       !Reaction parameters.

C----------REACTION PARAMETERS.
      INTEGER iform(nrec)          !Reaction type code (1-11).
      INTEGER ii(nrec)             !Incoming nuclide type (1-26).
      INTEGER jj(nrec)             !Incoming light nuclide type (1-6).
      INTEGER kk(nrec)             !Outgoing light nuclide type (1-6).
      INTEGER ll(nrec)             !Outgoing nuclide type (1-26).
      REAL    rev(nrec)            !Reverse reaction coefficient.
      REAL    q9(nrec)             !Energy released in reaction.

C----------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      REAL    r(nrec)              !Reverse reaction rate coefficients.

C----------DEFAULT COMPUTATION PARAMETERS.
      REAL    cy0                  !Default cy.
      REAL    ct0                  !Default ct.
      REAL    t9i0                 !Default t9i.
      REAL    t9f0                 !Default t9f.
      REAL    ytmin0               !Default ytmin.
      INTEGER inc0                 !Default accumulation increment.

C----------COMPUTATIONAL PARAMETERS.
      REAL    cy                   !Time step limiting constant on abundances.
      REAL    ct                   !Time step limiting constant on temperature.
      REAL    t9i                  !Initial temperature (in 10**9 K).
      REAL    t9f                  !Final temperature (in 10**9 k).
      REAL    ytmin                !Smallest abundances allowed.
      INTEGER inc                  !Accumulation increment.

C----------DEFAULT MODEL PARAMETERS.
      REAL    c0(3)                !Default c.
      REAL    cosmo0               !Default cosmological constant.
      REAL    xi0(3)               !Default neutrino degeneracy parameters.

C----------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    c(3)                 !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C----------DEFAULT VARIATIONAL PARAMETERS.
      REAL    dt0                  !Default initial time step.
      REAL    eta0                 !Default baryon-to-photon ratio.

C----------VARIATIONAL PARAMETERS.
      REAL    dt1                  !Initial time step.
      REAL    eta1                 !Baryon-to-photon ratio.

C----------COMPUTATION LOCATION.
      INTEGER itime                !Time check.

C----------RUN OPTION.
      INTEGER irun                 !Run network size.
      INTEGER isize                !Number of nuclides in computation.
      INTEGER jsize                !Number of reactions in computation.

C----------OUTPUT FILE STATUS.
      INTEGER nout                 !Number of output requests.
      LOGICAL outfile              !Indicates if output file used.

C----------USER RESPONSE VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION==================================

C10--------OPEN FILES AND PRINT GREETING----------------------------------------

      OPEN (unit=2, file='new123.dat', status='new')  !Output file.
      itime = 1                    !Time = beginning of program.
      CALL check                   !Check interface subroutine.
      WRITE (iw,1000)
 1000 FORMAT (6(/),
     |        2(' ',4x,'NN',6x,'NN  UU',6x,'UU',4x,8('C'),6x,'11',8x,
     |        6('2'),6x,6('3'),/),
     |        2(' ',4x,'NN',6x,'NN  UU',6x,'UU  CC',12x,'1111',6x,
     |        '22',6x,'22  33',6x,'33',/),
     |        2(' ',4x,'NNNN    NN  UU',6x,'UU  CC',14x,'11',14x,
     |        '22',10x,'33',/),
     |        2(' ',4x,'NN  NN  NN  UU',6x,'UU  CC',14x,'11',12x,
     |        '22',10x,'33',/),
     |        2(' ',4x,'NN    NNNN  UU',6x,'UU  CC',14x,'11',10x,
     |        '22',14x,'33',/),
     |        2(' ',4x,'NN',6x,'NN  UU',6x,'UU  CC',14x,'11',8x,
     |        '22',8x,'33',6x,'33',/),
     |        2(' ',4x,'NN',6x,'NN  ',10('U'),4x,8('C'),4x,6('1'),4x,
     |        10('2'),4x,6('3'),/),/,
     |        ' ',26x,'MODIFIED APRIL 1994',///,
     |        ' ','(Press <RETURN> to continue): ',$)

C20--------INPUT INITIALIZATION INFORMATION AND PAUSE---------------------------

      DO i  = 1,nrec
C..........READ IN REACTION PARAMETERS.
        iform(i) = int(reacpr(i,2))!Reaction type.
        ii(i)    = int(reacpr(i,3))!Incoming nuclide type.
        jj(i)    = int(reacpr(i,4))!Incoming nuclide type.
        kk(i)    = int(reacpr(i,5))!Outgoing nuclide type.
        ll(i)    = int(reacpr(i,6))!Outgoing nuclide type.
        rev(i)   = reacpr(i,7)     !Reverse reaction coefficient.
        q9(i)    = reacpr(i,8)     !Energy released.
C..........INITIALIZE REACTION RATES.
        f(i)  = 0.                 !Forward rate coeff.
        r(i)  = 0.                 !Reverse rate coeff.
C..........SET RUN OPTIONS TO DEFAULT.
      END DO
      irun       = 1               !Do full run.
      isize      = nnuc            !Use all 26 nuclides.
      jsize      = nrec            !Use all 88 reactions.
C..........SET OUTPUT OPTION TO DEFAULT.
      nout    = 0                  !No output requests.
      outfile = .false.            !Output file not used.
C..........SET VALUES TO DEFAULT.
      cy    = cy0                  !Time step limiting constant on abundances.
      ct    = ct0                  !Time step limiting constant on temperature.
      t9i   = t9i0                 !Initial temperature.
      t9f   = t9f0                 !Final temperature.
      ytmin = ytmin0               !Smallest abundances allowed.
      inc   = inc0                 !Accumulation increment.
      c(1)  = c0(1)                !Variation of gravitational constant.
      c(2)  = c0(2)                !Neutron lifetime.
      c(3)  = c0(3)                !Number of neutrino species.
      cosmo = cosmo0               !Cosmological constant.
      xi(1) = xi0(1)               !Electron degeneracy parameter.
      xi(2) = xi0(2)               !Muon degeneray parameter.
      xi(3) = xi0(3)               !Tauon degeneracy parameter.
      dt1   = dt0                  !Initial time step.
      eta1  = eta0                 !Baryon-to-photon ratio.
C..........ACCEPT RETURN TO CONTINUE.
      READ (ir,*)                  !Pause.

C30--------PRINT MENU AND AWAIT RESPONSE----------------------------------------

C..........RETURN FROM LOOPING.
 300  CONTINUE
C..........DISPLAY MENU.
      WRITE (iw,3000)
 3000 FORMAT (8(/),
     |        ' ',32x,'MENU SELECTION',/,
     |        ' ',32x,'---- ---------',//,
     |        ' ',24x,'1. HELP',/,
     |        ' ',24x,'2. SET COMPUTATION PARAMETERS',/,
     |        ' ',24x,'3. SET MODEL PARAMETERS',/,
     |        ' ',24x,'4. RUN',/,
     |        ' ',24x,'5. OUTPUT',/,
     |        ' ',24x,'6. EXIT',8(/),
     |        ' ',24x,'Enter selection (1-6): ',$)
C..........READ IN SELECTION NUMBER.
      READ (ir,3001) inum
 3001 FORMAT(i1)

C40--------BRANCH TO APPROPRIATE SECTION----------------------------------------

      GO TO (410,420,430,440,450,460),inum
      GO TO 460                    !Improper input or <RETURN>.
 410  CONTINUE                     !Help section.
        CALL help
        GO TO 500
 420  CONTINUE                     !Set computation parameters section.
        CALL setcom
        GO TO 500
 430  CONTINUE                     !Set model parameters section.
        CALL setmod
        GO TO 500
 440  CONTINUE                     !Run section.
        itime = 2                  !Time = beginning of run section.
        CALL check                 !Check interface subroutine.
        CALL run
        itime = 9                  !Time = end of run section.
        CALL check                 !Check interface subroutine.
        GO TO 500
 450  CONTINUE                     !Output section.
        CALL output
        GO TO 500
 460  CONTINUE                     !Exit section.
        IF (outfile) THEN 
          CLOSE (unit=2,status='keep')   !Close output file.
        ELSE
          CLOSE (unit=2,status='delete') !File not used - dispose.
        END IF
        itime = 10                 !Time = end of program.
        CALL check                 !Check interface subroutine.
        STOP

C50---------GO BACK TO MENU-----------------------------------------------------

 500  CONTINUE
      GO TO 300

      END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE help

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - none

C----------REMARKS.
C     Displays description and workings of the program.

C----------PARAMETERS.
      PARAMETER (ir=5)             !Input unit number (previous value = 1).
      PARAMETER (iw=6)             !Output unit number (previous value = 1).


C==========================DECLARATION DIVISION=================================

C----------USER RESPONSE VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION==================================

C10--------PRINT HELP SELECTION-------------------------------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY MENU.
      WRITE (iw,1000)
 1000 FORMAT (8(/),
     |        ' ',32x,'HELP SELECTION',/,
     |        ' ',32x,'---- ---------',//,
     |        ' ',24x,'1. INTRODUCTION',/,
     |        ' ',24x,'2. SETTING UP A RUN',/,
     |        ' ',24x,'3. RUNNING THE PROGRAM',/,
     |        ' ',24x,'4. OUTPUT OPTIONS',/,
     |        ' ',24x,'5. GENERAL METHOD OF COMPUTATION',/,
     |        ' ',24x,'6. USING THE INTERFACE SUBROUTINE',/,
     |        ' ',24x,'7. EXIT',7(/),
     |        ' ',24x,'Enter selection (1-7): ',$)
C..........READ IN SELECTION NUMBER.
      READ (ir,1001) inum
 1001 FORMAT (i1)

C20--------BRANCH TO APPROPRIATE SECTION----------------------------------------

      GO TO (210,220,230,240,250,260,270),inum
      GO TO 270                    !Improper input or <RETURN>.

C21--------INTRODUCTION SECTION-------------------------------------------------

 210  CONTINUE                     !Setting up a run section.
        WRITE (iw,2100)
 2100   FORMAT (/,
     |          ' ',31x,'INTRODUCTION',/,
     |          ' ',31x,'------------',2(/),
     |          ' ','Welcome to the wonderful world of primor',
     |              'dial nucleosynthesis.  NUC123 is a      ',/,
     |          ' ','FORTRAN program designed to provide the ',
     |              'early universe researcher with the tools',/,
     |          ' ','necessary for the investigation of primo',
     |              'rdial nucleosynthesis.  Its menu-driven ',/,
     |          ' ','interface allows the user to first set c',
     |              'omputation parameters (such as the time ',/,
     |          ' ','step) and model parameters (such as the ',
     |              'neutron lifetime and number of neutri-  ',/,
     |          ' ','nos) before doing single runs or multipl',
     |              'e runs (in which desired model parame-  ',/,
     |          ' ','ters are varied over a desired range.)  ',
     |              'After the run, the user can utilize the ',/,
     |          ' ','menu to either produce an output file or',
     |              ' to view the most recent run on the     ',/,
     |          ' ','screen.  The program comes with an empty',
     |              ' subroutine CHECK into which the user   ',/,
     |          ' ','may wish to put additional code to add t',
     |              'o the computation in an original manner.',10(/),
     |          ' ','(Enter <RETURN> to go back to help menu): ',$)
        READ (ir,*) 
        GO TO 300

C22--------SET UP RUN SECTION---------------------------------------------------

 220  CONTINUE                     !Setting up a run section.
        WRITE (iw,2200)
 2200   FORMAT (/,
     |          ' ',29x,'SETTING UP A RUN',/,
     |          ' ',29x,'------- -- - ---',2(/),
     |          ' ','I. Setting computation parameters.      ',/,
     |          ' ','   The accuracy of the computation and t',
     |              'he relevant temperature region can be   ',/,
     |          ' ','   set by the following parameters:     ',/,
     |          ' ','    A. Time step limiting constant 1  (d',
     |              'efault value of 0.03)                    ',/,
     |          ' ','    B. Time step limiting constant 2  (d',
     |              'efault value of 0.003)                   ',/,
     |          ' ','    C. Initial time step  (default value',
     |              ' of 10**-4)                             ',/,
     |          ' ','    D. Initial temperature  (default val',
     |              'ue of 10**2)                            ',/,
     |          ' ','       This is the temperature at the be',
     |              'ginning of the run in units of 10**9 K  ',/,
     |          ' ','    E. Final temperature  (default value',
     |              ' of 10**-2)                             ',/,
     |          ' ','       This is the termination temperatu',
     |              're of the run in units of 10**9 K       ',/,
     |          ' ','    F. Smallest abundances allowed  (def',
     |              'ault value of 10**-25)                  ',/,
     |          ' ','       Elemental abundances are not allo',
     |              'wed to drop below this value            ',/,
     |          ' ','    G. # of iterations for each accumula',
     |              'tion  (default value of 300)             ',/,
     |          ' ','       This is the number of iterations ',
     |              'before values are put in an output array',6(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (ir,1001) inum
        IF (inum.eq.1) THEN
          WRITE (iw,2202)
 2202     FORMAT (/,
     |            ' ','II. Setting model parameters.           ',/,
     |            ' ','   Default values here give what is know',
     |                'n as the standard model with best guess ',/,
     |            ' ','   figure on the neutron lifetime of 889',
     |                '.1 seconds.  Nonstandard scenarios can',/,
     |            ' ','   be investigated by varying the follow',
     |                'ing parameters:                         ',/,
     |            ' ','    A. The gravitational constant       ',/,
     |            ' ','       (The default value of one here gi',
     |                'ves the usual 6.6720e-8 dyne*cm**2/g**2)',/,
     |            ' ','    B. Neutron life-time  (default value',
     |                ' of 889.1 seconds)                    ',/,
     |            ' ','    C. Number of neutrino species  (defa',
     |                'ult value of 3 light neutrinos)         ',/,
     |            ' ','    D. Final baryon-to-photon ratio  (se',
     |                't to log(eta) = -9.5)                   ',/,
     |            ' ','    E. Cosmological constant  (default v',
     |                'alue of 0)                              ',/,
     |            ' ','    F. Neutrino degeneracy parameters  (',
     |                'default values all 0)                   ',/,
     |            ' ','       There are 3 separate parameters f',
     |                'or the electron, muon, and tau neutrinos',11(/),
     |            ' ','(Enter <RETURN> to go back to help menu): ',$)
          READ (ir,*)
          GO TO 300
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C23--------RUN PROGRAM SECTION--------------------------------------------------

 230  CONTINUE                     !Running the program section.
        WRITE (iw,2300)
 2300   FORMAT (/,
     |          ' ',28x,'RUNNING THE PROGRAM',/,
     |          ' ',28x,'------- --- -------',2(/),
     |          ' ','I. Setting run speed.                   ',/,
     |          ' ','   The code can be run at 3 different se',
     |              'ttings of speed.  The running of the    ',/,
     |          ' ','   code can be speeded up by reducing th',
     |              'e number of nuclides and reactions.  The',/,
     |          ' ','   complete computation takes into accou',
     |              'nt the following nuclides: n, p, d, t,  ',/,
     |          ' ','   He3, He4, Li6, Li7, Be7, Li8, B8, Be9',    
     |              ',B10, B11, C11, B12, C12, N12, C13, N13,',/,
     |          ' ','   C14, N14, O14, N15, O15, and O16.    ',/,
     |          ' ','   The given CPU percentages and abundan',
     |              'ce variations are with regard to a      ',/,
     |          ' ','   single run with all default parameter',
     |              ' values.                                ',/,
     |          ' ','    A. 26 nuclides, 88 reactions (defaul',
     |              't)                                      ',/,
     |          ' ','       nuclides from n to O16           ',/,
     |          ' ','    B. 18 nuclides, 60 reactions        ',/,
     |          ' ','       nuclides from n to N12           ',/,
     |          ' ','       (63% CPU time, variation = .1%)  ',/,
     |          ' ','    C.  9 nuclides, 25 reactions        ',/,
     |          ' ','       nuclides from n to Be7           ',/,
     |          ' ','       (20% CPU time, variation = .5%)  ',4(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (ir,1001) inum
        IF (inum.eq.1) THEN
          WRITE (iw,2302)
 2302     FORMAT (/,
     |            ' ','II. Do single run.                      ',/,
     |            ' ','    A. Interactive.                     ',/,
     |            ' ','       In an interactive session, the us',
     |                'er can readily input the computational  ',/,
     |            ' ','       and model parameters and begin th',
     |                'e computation process.  The run itself  ',/,
     |            ' ','       is commenced when option 2, "GO",',
     |                ' in the "RUN" section is requested.     ',//,
     |            ' ','    B. Batch.                           ',/,
     |            ' ','       To run the program in a batch mod',
     |                'e, it must be altered slightly so that  ',/,
     |            ' ','       the I/O takes place with files in',
     |                'stead of a terminal.  This is done by   ',/,
     |            ' ','       setting different values for the ',
     |                'input and output unit number parameters ',/,
     |            ' ','       "ir" and "iw" and assigning them ',
     |                'to different files in NUC123.  In the   ',/,
     |            ' ','       file assigned the "ir" unit numbe',
     |                'r, one must place the responses to the  ',/,
     |            ' ','       queries of the program.          ',10(/),
     |            ' ','(Enter 1 to continue, <RETURN> to end): ',$)
          READ (ir,1001) inum
          IF (inum.eq.1) THEN
            WRITE (iw,2304)
 2304       FORMAT (/,
     |              ' ','III. Do multiple runs.                 ',/,
     |              ' ','   A wide range of early universe model',
     |                  's can be covered by doing many runs    ',/,
     |              ' ','   while one or more parameters are var',
     |                  'ied over a range of interest.  The     ',/,
     |              ' ','   parameters that can be varied are th',
     |                  'e following:                           ',/,
     |              ' ','    A. Eta                             ',
     |                  '       - Logrithmic variation          ',/,
     |              ' ','    B. Gravitational constant          ',
     |                  '       - Linear variation              ',/,
     |              ' ','    C. Neutron lifetime                ',
     |                  '       - Linear variation              ',/,
     |              ' ','    D. Number of neutrino species      ',
     |                  '       - Linear variation              ',/,
     |              ' ','    E. Cosmological constant           ',
     |                  '       - Linear variation              ',/,
     |              ' ','    F. Neutrino degeneracy parameters  ',
     |                  '       - Linear variation              ',/,
     |              ' ','        1. Electron neutrino           ',/,
     |              ' ','        2. Muon neutrino               ',/,
     |              ' ','        3. Tauon neutrino              ',/,
     |              ' ','   At most 3 parameters can be varied. ',
     |                  ' The first parameter inputted will be  ',/,
     |              ' ','   will be varied in the outermost loop',
     |                  ' and the third parameter inputted will ',/,
     |              ' ','   be varied in the innermost loop.    ',7(/),
     |              ' ','(Enter <RETURN> to go back to help menu): ',$)
            READ (ir,*)
            GO TO 300
          ELSE
            GO TO 300
          END IF !(inum.eq.1)
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C24--------OUTPUT OPTIONS SECTION-----------------------------------------------

 240  CONTINUE                     !Output options section.
        WRITE (iw,2400)
 2400   FORMAT (/,
     |          ' ',30x,'OUTPUT OPTIONS',/,
     |          ' ',30x,'------ -------',2(/),
     |          ' ','I.  Request output file.                ',/,
     |          ' ','   After a run, the user can request the',
     |              ' program to put the resulting numbers   ',/,
     |          ' ','   into an output file.  This can be don',
     |              'e as many times as desired and all the  ',/,
     |          ' ','   information will be put in one new fi',
     |              'le under the name of "NUC123.DAT."  If  ',/,
     |          ' ','   there is no request during the entire',
     |              ' running of the program, this file is   ',/,
     |          ' ','   not created.  If an output file is re',
     |              'quested after a multiple run, only the  ',/,
     |          ' ','   information from the very last run wi',
     |              'll be given.  The output file will give ',/,
     |          ' ','   the computational and model parameter',
     |              's for each run and will contain the     ',/,
     |          ' ','   following information:               ',/,
     |          ' ','    A. Temperatures in decreasing order ',/,
     |          ' ','    B. Abundances for n, p, d, t, He3, H',
     |              'e4, Li6, Li7, Be7, and Li8 & up         ',/,
     |          ' ','       (p and He4 are in mass fraction, ',
     |              'the rest in ratios to the p abundance)  ',/,
     |          ' ','    C. Time, time interval, chemical pot',
     |              'ential of the electron                  ',/,
     |          ' ','    D. Energy densities for photons, ele',
     |              'ctrons, electron neutrinos, and baryons ',/,
     |          ' ','    E. Baryon-to-photon ratio, expansion',
     |              ' rate of the universe                   ',5(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (ir,1001) inum
        IF (inum.eq.1) THEN
          WRITE (iw,2402)
 2402     FORMAT (/,
     |            ' ','II.  Request output on screen.         ',/,
     |            ' ','   Instead of waiting to print out an o',
     |                'utput file, the user can immediately   ',/,
     |            ' ','   access the results of the latest run',
     |                ' by requesting the output on the       ',/,
     |            ' ','   screen.  There are four screens on e',
     |                'ach of which are displayed the         ',/,
     |            ' ','   computational and model parameters a',
     |                'nd the temperature:                    ',/,
     |            ' ','    A. Abundances for d, t, He3, He4, a',
     |                'nd Li7                                 ',/,
     |            ' ','       (He4 in mass fraction, rest as a',
     |                ' ratio with the p abundance)           ',/,
     |            ' ','    B. Abundances for n, p, Li6, Be7, a',
     |                'nd Li8 & up                            ',/,
     |            ' ','       (p in mass fraction, rest as a r',
     |                'atio with the p abundance)             ',/,
     |            ' ','    C. Energy densities for photons, el',
     |                'ectrons, electron neutrinos, & baryons ',/,
     |            ' ','    D. Time, time interval, chemical po',
     |                'tential of the electron,               ',/,
     |            ' ','       baryon-to-photon ratio, and expa',
     |                'nsion rate of the universe             ',11(/),
     |            ' ','(Enter <RETURN> to go back to help menu): ',$)
          READ (ir,*)
          GO TO 300
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C25--------METHOD OF COMPUTATION SECTION----------------------------------------

 250  CONTINUE                     !General method of computation section.
        WRITE (iw,2500)
 2500   FORMAT (/,
     |          ' ',22x,'GENERAL METHOD OF COMPUTATION',/,
     |          ' ',22x,'------- ------ -- -----------',2(/),
     |          ' ','I. Time evolution algorithm.            ',/,
     |          ' ','   The program utilizes a 2-point Runge-',
     |              'Kutta scheme (located in subroutine     ',/,
     |          ' ','   DRIVER) to time-evolve the temperatur',
     |              'e, the quantity hv (the ratio of the    ',/,
     |          ' ','   baryon density to T**3), the chemical',
     |              ' potential of the electron, and the     ',/,
     |          ' ','   nuclide abundances.  In the 2-point R',
     |              'unge-Kutta routine, a variable v at time',/,
     |          ' ','   t0 (= v0) is evolved to a time t1 by ',
     |              'adding to v0 the average of the         ',/,
     |          ' ','   derivatives evaluated at t0 and at t1',
     |              ' multiplied by dt:                      ',/,
     |          ' ','       v1 = v0 + 0.5(dvdt(t0)+dvdt(t1)) ',/,
     |          ' ','   where dvdt(t1) is gotten by first fin',
     |              'ding v1'' = v0 + dvdt(t0).  The         ',/,
     |          ' ','   derivatives of the nuclide abundances',
     |              ' are first computed and these are used  ',/,
     |          ' ','   to find the derivatives of t9, hv, an',
     |              'd phie (this is done in subroutine      ',/,
     |          ' ','   DERIVS).  To compute the time derivat',
     |              'ives of the nuclide abundances, a matrix',/,
     |          ' ','   equation is set up (in subroutine SOL',
     |              ') and is solved (in subroutine EQSLIN)  ',/,
     |          ' ','   by gaussian elimination utilizing imp',
     |              'licit differentiation.                  ',6(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (ir,1001) inum
        IF (inum.eq.1) THEN
          WRITE (iw,2502)
 2502     FORMAT (/
     |            ' ','II. Hierarchy of Subroutines.   ',/,
     |            ' ','    NUC123                       ',
     |                '     Main program (main menu)    ',/,
     |            ' ','        HELP                     ',
     |                '     Help option                 ',/,
     |            ' ','        SETCOM                   ',
     |                '     Set computational parameters',/,
     |            ' ','        SETMOD                   ',
     |                '     Set model parameters        ',/,
     |            ' ','        RUN                      ',
     |                '     Run computation code        ',/,
     |            ' ','            DRIVER               ',
     |                '     Main routine (Runge-Kutta loop)    ',/,
     |            ' ','                START            ',
     |                '     Initialization routine      ',/,
     |            ' ','                    RATE0        ',
     |                '     Computes weak decay rates   ',/,
     |            ' ','                DERIVS           ',
     |                '     Computes time derivatives   ',/,
     |            ' ','                    THERM        ',
     |                '     Computes energy densities   ',/,
     |            ' ','                        BESSEL   ',
     |                '     Gives functions of Kn       ',/,
     |            ' ','                            KNUX ',
     |                '     Computes modified Bessel fcn Kn    ',/,
     |            ' ','                        NUDENS   ',
     |                '     Computes neutrino energy density   ',/,
     |            ' ','                    RATE1-4      ',
     |                '     Computes rates for reactions',/,
     |            ' ','                    SOL          ',
     |                '     Builds A matrix for eqn dy/dt = Ay ',/,
     |            ' ','                        EQSLIN   ',
     |                '     Solves dy/dt=Ay by gaussian elim   ',/,
     |            ' ','                ACCUM            ',
     |                '     Output accumulator          ',/,
     |            ' ','        OUTPUT                   ',
     |                '     Allows user to output result',4(/),
     |            ' ','(Enter <RETURN> to go back to help menu): ',$)
          READ (ir,*)
          GO TO 300
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C26--------USING INTERFACE SUBROUTINE SECTION.

 260  CONTINUE                     !Using the interface subroutine section.
        WRITE (iw,2600)
 2600   FORMAT (/,
     |          ' ',22x,'USING THE INTERFACE SUBROUTINE',/,
     |          ' ',22x,'----- --- --------- ----------',2(/),
     |          ' ','I. Purpose.                             ',/,
     |          ' ','   The interface subroutine CHECK is des',
     |              'igned to be an outlet of the program    ',/,
     |          ' ','   into which alterations can be easily ',
     |              'plugged.  Programs are normally modified',/,
     |          ' ','   by searching through the program, ide',
     |              'ntifying the appropriate areas for      ',/,
     |          ' ','   alterations, and interspersing new co',
     |              'mmands while deleting some old ones.    ',/,
     |          ' ','   This process can get tricky unless on',
     |              'e actively documents the alterations:   ',/,
     |          ' ','   one might lose track of all of the mo',
     |              'difications and deletions.  Thus, it is ',/,
     |          ' ','   worthwhile to put most if not all of ',
     |              'the necessary changes into one          ',/,
     |          ' ','   subroutine which is to be called from',
     |              ' strategic locations in the main        ',/,
     |          ' ','   program.  Furthermore, by putting cha',
     |              'nges into one small subroutine, one need',/,
     |          ' ','   only to compile the subroutine CHECK ',
     |              'each time instead of the entire nucleo- ',/,
     |          ' ','   synthesis code.                      ',8(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (ir,1001) inum
        IF (inum.eq.1) THEN
          WRITE (iw,2602)
 2602     FORMAT (/,
     |            ' ','II. Description.                        ',/,
     |            ' ','   Subroutine CHECK is an empty subrouti',
     |                'ne with a large COMMON area, giving the ',/,
     |            ' ','   user ready access to all of the impor',
     |                'tant variables in the computations.  The',/,
     |            ' ','   routine is called from various locati',
     |                'ons in the main program and the location',/,
     |            ' ','   spot in the program is labeled by the'
     |                ' flag "itime".  The set call locations  ',/,
     |            ' ','   are given below:                     ',/,
     |            ' ','    A. itime = 1 (NUC123, very beginning',
     |                ' of program run)                        ',/,
     |            ' ','       (appropriate for opening files, i',
     |                'nitializing variables)                  ',/,
     |            ' ','    B. itime = 2 (NUC123, right before g',
     |                'oing into the RUN section)              ',/,
     |            ' ','    C. itime = 3 (RUN, right before goin',
     |                'g into DRIVER to do the computations)   ',/,
     |            ' ','    D. itime = 4 (DRIVER, in 1st R-K loo',
     |                'p after computing derivatives in DERIVS)',/,
     |            ' ','    E. itime = 7 (DRIVER, in 2nd R-K loo',
     |                'p after computing derivatives in DERIVS)',/,
     |            ' ','    F. itime = 8 (RUN, right after comin',
     |                'g back from DRIVER)                     ',/,
     |            ' ','    G. itime = 9 (NUC123, right after co',
     |                'ming back from the RUN section)         ',/,
     |            ' ','    H. itime =10 (NUC123, very end of pr',
     |                'ogram run)                              ',/,
     |            ' ','       (appropriate for closing files)  ',/,
     |            ' ','   The difference between the (2,9) pair',
     |                'ing and the (3,8) pairing is that for a ',/,
     |            ' ','   multiple run, the (3,8) pairing would',
     |                ' be called before and after every run   ',/,
     |            ' ','   but the (2,9) pairing would be called',
     |                ' before and after the entire sequence.  ',4(/),
     |            ' ','(Enter 1 to continue, <RETURN> to end): ',$)
          READ (ir,1001) inum
          IF (inum.eq.1) THEN
            WRITE (iw,2604)
 2604       FORMAT (/,
     |              ' ','III. Implementation.                   ',/,
     |              ' ','   The additional program statements ar',
     |                  'e needed in the subroutine CHECK.  If a',/,
     |              ' ','   particular command is to be executed',
     |                  ' when the computer is at a certain     ',/,
     |              ' ','   location in the program -- say label',
     |                  'ed by itime = 8 -- then in CHECK, one  ',/,
     |              ' ','   must place the command under the sta',
     |                  'tement, IF (itime.eq.8)....  The user  ',/,
     |              ' ','   is at leisure to place his own locat',
     |                  'ion indicators (5,6) and CALL CHECK    ',/,
     |              ' ','   statements anywhere in the program a',
     |                  's long as there is a COMMON /checkcb/    ',/,
     |              ' ','   statement in the particular subrouti',
     |                  'ne to carry the value of itime along.  ',15(/),
     |              ' ','(Enter <RETURN> to go back to help menu): ',$)
            READ (ir,*)
            GO TO 300
          ELSE
            GO TO 300
          END IF !(inum.eq.1)
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C27--------EXIT SECTION---------------------------------------------------------

 270  CONTINUE                     !Exit section.
        RETURN

C30--------GO BACK TO MAIN MENU-------------------------------------------------

 300  CONTINUE
      GO TO 100

      END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE setcom

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - none

C----------REMARKS.
C     Allows resetting of computation parameters.

C----------PARAMETERS.
      PARAMETER (ir=5)             !Input unit number (previous value = 1).
      PARAMETER (iw=6)             !Output unit number (previous value = 1).

C----------COMMON AREAS.
      COMMON /compr0/ cy0,ct0,t9i0,t9f0,ytmin0,inc0  !Default comp parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /varpr0/ dt0,eta0                       !Default variationl params.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.


C==========================DECLARATION DIVISION=================================

C----------DEFAULT COMPUTATION PARAMETERS.
      REAL    cy0                  !Default cy.
      REAL    ct0                  !Default ct.
      REAL    t9i0                 !Default t9i.
      REAL    t9f0                 !Default t9f.
      REAL    ytmin0               !Default ytmin.
      INTEGER inc0                 !Default accumulation increment.

C----------COMPUTATION PARAMETERS.
      REAL    cy                   !Time step limiting constant on abundances.
      REAL    ct                   !Time step limiting constant on temperature.
      REAL    t9i                  !Initial temperature (in 10**9 K).
      REAL    t9f                  !Final temperature (in 10**9 K).
      REAL    ytmin                !Smallest abundances allowed.
      INTEGER inc                  !Accumulation increment.

C----------DEFAULT VARIATIONAL  PARAMETERS.
      REAL    dt0                  !Default initial dt.

C----------VARIATIONAL  PARAMETERS.
      REAL    dt1                  !Initial time step.

C----------LOCAL VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION==================================

C10--------PRINT RESET SELECTION AND AWAIT RESPONSE-----------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY RESET SELECTIONS.
      WRITE (iw,1000) cy,ct,dt1,t9i,t9f,ytmin,float(inc)
 1000 FORMAT (8(/),
     |        ' ',21x,'SET COMPUTATION PARAMETERS SELECTION',/,
     |        ' ',21x,'--- ----------- ---------- ---------',//,
     |        ' ',10x,' 1. CHANGE TIME-STEP LIMITING CONSTANT 1  FROM ',
     |            f5.3,/,
     |        ' ',10x,' 2. CHANGE TIME-STEP LIMITING CONSTANT 2  FROM ',
     |            f5.3,/,
     |        ' ',10x,' 3. CHANGE INITIAL TIME-STEP              FROM ',
     |            1pe8.2,' SECONDS',/,
     |        ' ',10x,' 4. CHANGE INITIAL TEMPERATURE            FROM ',
     |            1pe8.2,' (10**9 K)',/,
     |        ' ',10x,' 5. CHANGE FINAL TEMPERATURE              FROM ',
     |            1pe8.2,' (10**9 K)',/,
     |        ' ',10x,' 6. CHANGE SMALLEST ABUNDANCES ALLOWED    FROM ',
     |            1pe8.2,/,
     |        ' ',10x,' 7. CHANGE ACCUMULATION INCREMENT         FROM ',
     |            1pe8.2,' ITERATIONS',/,
     |        ' ',10x,' 8. RESET ALL TO DEFAULT VALUES',/,
     |        ' ',10x,' 9. EXIT',5(/),
     |        ' ',10x,'Enter selection (1-9): ',$)
C..........READ IN SELECTION NUMBER.
      READ (ir,1001) inum
1001  FORMAT (i1)

C20--------BRANCH TO APPROPRIATE SECTION----------------------------------------

      GO TO (210,220,230,240,250,260,270,280,300),inum
      GO TO 300                    !Improper input or <RETURN>.
 210  CONTINUE                     !Change time step limiting const 1 section.
        WRITE (iw,2100)
 2100   FORMAT (' ','Enter value for time step limiting constant 1: ',$)
        READ (ir,*) cy
 2101   FORMAT (f5.3)
        GO TO 400
 220  CONTINUE                     !Change time step limiting const 2 section.
        WRITE (iw,2200)
 2200   FORMAT (' ','Enter value for time step limiting constant 2: ',$)
        READ (ir,*) ct
        GO TO 400
 230  CONTINUE                     !Change initial time step section.
        WRITE (iw,2300)
 2300   FORMAT (' ','Enter value for initial time step: ',$)
        READ (ir,*) dt1
        GO TO 400
 240  CONTINUE                     !Change initial temperature section.
        WRITE (iw,2400)
 2400   FORMAT (' ','Enter value for initial temperature: ',$)
        READ (ir,*) t9i
        GO TO 400
 250  CONTINUE                     !Change final temperature section.
        WRITE (iw,2500)
 2500   FORMAT (' ','Enter value for final temperature: ',$)
        READ (ir,*) t9f
        GO TO 400
 260  CONTINUE                     !Change smallest abundances allowed section.
        WRITE (iw,2600)
 2600   FORMAT (' ','Enter value for smallest abundances allowed: ',$)
        READ (ir,*) ytmin
        GO TO 400
 270  CONTINUE                     !Change accumulation increment section.
        WRITE (iw,2700)
 2700   FORMAT (' ','Enter value for accumulation increment: ',$)
        READ (ir,*) inc
        GO TO 400
 280  CONTINUE                     !Reset all to default values section.
        cy    = cy0                !Time step limiting constant on abundances.
        ct    = ct0                !Time step limiting constant on temperature.
        dt1   = dt0                !Time step.
        t9i   = t9i0               !Initial temperature.
        t9f   = t9f0               !Final temperature.
        ytmin = ytmin0             !Smallest abundances allowed.
        inc   = inc0               !Accumulation increment.
        WRITE (iw,2800)
 2800   FORMAT (' ','All values reset to default - Press <RETURN> '
     |              'to continue: ',$)
        READ (ir,*)
        GO TO 400
 300  CONTINUE                     !Exit section.
        RETURN

C40--------GO BACK TO MENU------------------------------------------------------

 400  CONTINUE
      GO TO 100

      END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE setmod

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - none

C----------REMARKS.
C     Allows resetting of model parameters.

C----------PARAMETERS.
      PARAMETER (ir=5)             !Input unit number (previous value = 1).
      PARAMETER (iw=6)             !Output unit number (previous value = 1).

C----------COMMON AREAS.
      COMMON /modpr0/ c0,cosmo0,xi0                  !Default model parameters.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /varpr0/ dt0,eta0                       !Default variationl params.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.


C==========================DECLARATION DIVISION=================================

C----------DEFAULT MODEL PARAMETERS.
      REAL    c0(3)                !Default c.
      REAL    cosmo0               !Default cosmological constant.
      REAL    xi0(3)               !Default neutrino degeneracy parameters.

C----------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    c(3)                 !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C----------DEFAULT VARIATIONAL PARAMETERS.
      REAL    eta0                 !Default eta.

C----------VARIATIONAL PARAMETERS.
      REAL    eta1                 !Intial baryon-to-photon ratio.

C----------USER RESPONSE VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION==================================

C10--------PRINT RESET SELECTION AND AWAIT RESPONSE-----------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY RESET SELECTIONS.
      WRITE (iw,1000) c(1),c(2),c(3),eta1,cosmo,xi(1),xi(2),xi(3)
 1000 FORMAT (8(/),
     |        ' ',24x,'SET MODEL PARAMETERS SELECTION',/,
     |        ' ',24x,'--- ----- ---------- ---------',//,
     |        ' ',10x,' 1. CHANGE GRAVITATIONAL CONSTANT         FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 2. CHANGE NEUTRON LIFETIME               FROM ',
     |            1pe10.3,' SECONDS',/,
     |        ' ',10x,' 3. CHANGE NUMBER OF NEUTRINO SPECIES     FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 4. CHANGE FINAL BARYON-TO-PHOTON RATIO   FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 5. CHANGE COSMOLOGICAL CONSTANT          FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 6. CHANGE XI-ELECTRON                    FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 7. CHANGE XI-MUON                        FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 8. CHANGE XI-TAUON                       FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 9. RESET ALL TO DEFAULT VALUES',/,
     |        ' ',10x,'10. EXIT',4(/),
     |        ' ',10x,' Enter selection (1-10): ',$)
C..........READ IN SELECTION NUMBER.
      READ (ir,1001) inum
 1001 FORMAT (i2)

C20--------BRANCH TO APPROPRIATE SECTION----------------------------------------

      GO TO (210,220,230,240,250,260,270,280,290,300),inum
      GO TO 300                    !Improper input or <RETURN>.
 210  CONTINUE                     !Change gravitational constant section.
        WRITE (iw,2100)
 2100   FORMAT (' ','Enter value for variation of gravitational ',
     |              'constant: ',$)
        READ (ir,*) c(1)
        GO TO 400
 220  CONTINUE                     !Change neutron lifetime section.
        WRITE (iw,2200)
 2200   FORMAT (' ','Enter value for neutron lifetime (sec): ',$)
        READ (ir,*) c(2)
        GO TO 400
 230  CONTINUE                     !Change number of neutrino species section.
        WRITE (iw,2300)
 2300   FORMAT (' ','Enter value for number of neutrino species: ',$)
        READ (ir,*) c(3)
        GO TO 400
 240  CONTINUE                     !Change baryon-to-photon ratio section.
        WRITE (iw,2400)
 2400   FORMAT (' ','Enter value for baryon-to-photon ratio: ',$)
        READ (ir,*) eta1
        GO TO 400
 250  CONTINUE                     !Change cosmological constant section.
        WRITE (iw,2500)
 2500   FORMAT (' ','Enter value for cosmological constant: ',$)
        READ (ir,*) cosmo
        GO TO 400
 260  CONTINUE                     !Change neutrino degeneracy section.
        WRITE (iw,2600)
 2600   FORMAT (' ','Enter value for xi electron: ',$)
        READ (ir,*) xi(1)
        GO TO 400
 270  CONTINUE                     !Change neutrino degeneracy section.
        WRITE (iw,2700)
 2700   FORMAT (' ','Enter value for xi muon: ',$)
        READ (ir,*) xi(2)
        GO TO 400
 280  CONTINUE                     !Change neutrino degeneracy section.
        WRITE (iw,2800)
 2800   FORMAT (' ','Enter value for xi tauon: ',$)
        READ (ir,*) xi(3)
        IF ((xi(3).ne.0.).and.(c(3).lt.3.)) THEN
          c(3) = 3.
          WRITE (iw,2802)
 2802     FORMAT (' ','Number of neutrinos set to 3')
          WRITE (iw,2804)
 2804     FORMAT (' ','Press <RETURN> to continue: ',$)
          READ (ir,*)
        END IF
        GO TO 400
 290  CONTINUE                     !Reset all to default values section.
        c(1)   = c0(1)
        c(2)   = c0(2)
        c(3)   = c0(3)
        cosmo  = cosmo0
        xi(1)  = xi0(1)
        xi(2)  = xi0(2)
        xi(3)  = xi0(3)
        eta1   = eta0
        WRITE (iw,2900)
 2900   FORMAT (' ','All values reset to default - Press <RETURN> '
     |              'to continue: ',$)
        READ (ir,*)
        GO TO 400
 300  CONTINUE                     !Exit section.
        RETURN

C40--------GO BACK TO MENU------------------------------------------------------

 400  CONTINUE
      GO TO 100

      END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE run

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - [subroutine] driver

C----------REMARKS.
C     Activates computation routine.

C----------PARAMETERS.
      PARAMETER (ir=5)             !Input unit number (previous value = 1).
      PARAMETER (iw=6)             !Output unit number (previous value = 1).
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (lrec=64)          !Total # of nuclear reactions for irun = 2.
      PARAMETER (krec=34)          !Total # of nuclear reactions for irun = 3.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (lnuc=18)          !Total # of nuclides for irun = 2.
      PARAMETER (knuc=9)           !Total # of nuclides for irun = 3.

C----------COMMON AREAS.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.
      COMMON /checkcb/  itime                        !Computation location.
      COMMON /runopt/ irun,isize,jsize               !Run options.


C==========================DECLARATION DIVISION=================================

C----------MODEL PARAMETERS.
      REAL    eta1                 !Baryon-to-photon ratio.
      REAL    c(3)                 !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C----------RUN OPTION.
      INTEGER irun                 !Run network size.
      INTEGER isize                !Number of nuclides in computation.
      INTEGER jsize                !Number of reactions in computation.

C----------USER INTERACTION VARIABLES.
      REAL    rnumb1               !Run parameter for outer loop.
      REAL    rnumb2               !Run parameter for middle loop.
      REAL    rnumb3               !Run parameter for inner loop.
      REAL    rnum1(3)             !Run parameter starting value.
      REAL    rnum2(3)             !Run parameter end value.
      REAL    rnum3(3)             !Run parameter increment.
      INTEGER inumb                !Selection number.
      INTEGER inum(3)              !Selection number.
      INTEGER jnum                 !Number of loopings to be done.
      INTEGER knum                 !Number of loopings rejected.
      INTEGER lnumb1               !Run parameter for outer loop.
      INTEGER lnumb2               !Run parameter for middle loop.
      INTEGER lnumb3               !Run parameter for inner loop.
      INTEGER lnum(3)              !Run parameter end value.
      CHARACTER lchose             !User response (alphanumeric).

C----------FLAG AND LABELS.
      INTEGER itime                !Computation location.
      CHARACTER*22 vtype(8)        !Label for quantities being varied.

C----------EQUIVALENCE VARIABLE.
C     REAL    qvary(7)             !Array set equal to c, cosmo, and xi.

C----------EQUIVALENCE STATEMENTS.
C     EQUIVALENCE (qvary(1),c(1)), (qvary(4),cosmo), (qvary(5),xi(1))


C==============================DATA DIVISION====================================

C----------LABELS FOR QUANTITIES BEING VARIED.
      DATA vtype /'baryon/photon ratio   ',
     |            'gravitational constant',
     |            'neutron lifetime      ',
     |            '# of neutrino species ',
     |            'cosmological constant ',
     |            'xi-electron           ',
     |            'xi-muon               ',
     |            'xi-tauon              '/


C===========================PROCEDURE DIVISION==================================

C10--------PRINT RUN SELECTION AND AWAIT RESPONSE-------------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY RUN SELECTIONS.
      WRITE (iw,1000)
 1000 FORMAT (8(/),
     |        ' ',32x,'RUN SELECTION',/,
     |        ' ',32x,'--- ---------',//,
     |        ' ',27x,' 1. SET RUN NETWORK',/,
     |        ' ',27x,' 2. GO',/,
     |        ' ',27x,' 3. DO MULTIPLE RUNS',/,
     |        ' ',27x,' 4. EXIT',10(/),
     |        ' ',27x,' Enter selection (1-4): ',$)
C..........READ IN SELECTION NUMBER.
      READ (ir,1001) inumb
 1001 FORMAT (i1)

C20--------BRANCH TO APPROPRIATE SECTION----------------------------------------

      GO TO (210,220,230,240),inumb
      GO TO 240                    !Improper input or <RETURN>.

C21--------SET RUN NETWORK SECTION----------------------------------------------

 210  CONTINUE
        WRITE (iw,2100)
 2100   FORMAT (' ','Enter network size (1-26 nuclides (default); ',
     |              '2-18; 3-9): ',$)
        READ (ir,*) inumb          !Read in selection number.
        IF ((inumb.ne.1).and.(inumb.ne.2).and.(inumb.ne.3)) inumb = 1  !Default.
        IF (inumb.ne.irun) THEN    !Run network changed from previously.
          irun = inumb             !Run network size selection.
        END IF
        IF (irun.eq.1) THEN        !Maximal network size.
          isize = nnuc
          jsize = nrec
        ELSE
          IF (irun.eq.2) THEN      !Abridged network size.
            isize = lnuc
            jsize = lrec
          ELSE
            IF (irun.eq.3) THEN    !Minimal network size.
              isize = knuc
              jsize = krec
            END IF
          END IF
        END IF !(irun.eq.1) 
        WRITE (iw,2104) irun
 2104   FORMAT (' ','Run network set to ',i1,' - Press <RETURN> ',
     |              'to continue: ',$)
        READ (ir,*)
        GO TO 300

C22--------GO SECTION-----------------------------------------------------------

 220  CONTINUE
        WRITE (iw,2200)
 2200   FORMAT (' ','Begin computation run....')
        itime = 3
        CALL check                 !Call interface subr before computation.
        CALL driver                !Do nucleosynthesis computation.
        itime = 8
        CALL check                 !Call interface subr after computation.
        WRITE (iw,2202)
 2202   FORMAT (' ','Computation completed - Press <RETURN> to '
     |              'continue: ',$)
        READ (ir,*)
        GO TO 300

C23--------DO MULTIPLE RUNS SECTION---------------------------------------------

C..........GET NUMBER OF LOOPINGS.
 230  CONTINUE
        WRITE (iw,2300)
 2300   FORMAT (' ','Enter the number of loopings to be done (1 ',
     |              '(default); 2; 3): ',$)
        READ (ir,*) jnum           !Read in number of loopings to be done.
        IF ((jnum.ne.1).and.(jnum.ne.2).and.(jnum.ne.3)) THEN
          jnum = 1                 !Default number of loopings.
        END IF
        knum = 0.                  !No loopings rejected for now.
        DO i = 1,3
          IF (i.gt.jnum) THEN
            rnum1(i) = 0.          !Initialize initial parameter.
            rnum2(i) = 0.          !Initialize terminal parameter.
            rnum3(i) = 1.          !Initialize incremental parameter.
            inum(i)  = 0           !Initialize selection number.
          ELSE
C..........OBTAIN QUANTITY TO VARY.
            WRITE (iw,2302)
 2302       FORMAT (8(/),
     |              ' ',30x,'QUANTITY TO VARY',/,
     |              ' ',30x,'-------- -- ----',//,
     |              ' ',25x,' 1.  ETA     (LOGRITHMIC VARIATION)',/,
     |              ' ',25x,' 2.  G           (LINEAR VARIATION)',/,
     |              ' ',25x,' 3.  TAU         (LINEAR VARIATION)',/,
     |              ' ',25x,' 4.  # NEUTRINOS (LINEAR VARIATION)',/,
     |              ' ',25x,' 5.  LAMBDA      (LINEAR VARIATION)',/,
     |              ' ',25x,' 6.  XI-ELECTRON (LINEAR VARIATION)',/,
     |              ' ',25x,' 7.  XI-MUON     (LINEAR VARIATION)',/,
     |              ' ',25x,' 8.  XI-TAUON    (LINEAR VARIATION)',/,
     |              ' ',25x,' 9.  NO SELECTION',5(/),
     |              ' ',25x,' Enter selection (1-9): ',$)
            READ (ir,1001) inum(i)
            IF ((inum(i).lt.1).or.(inum(i).gt.8)) THEN  !No selection made.
              WRITE (iw,2304)
 2304         FORMAT (' ','No selection made - Reduce number of ',
     |                    'loopings by one',/,
     |                ' ','Press <RETURN> to continue: ',$)
              READ (ir,*)
              knum = knum + 1   !Step up number of loopings rejected.
              rnum1(i) = 0.     !Initialize initial parameter.
              rnum2(i) = 0.     !Initialize terminal parameter.
              rnum3(i) = 1.     !Initialize incremental parameter.
              inum(i)  = 0      !Initialize selection number.
            ELSE !((inum(i).ge.1).and.(inum(i).le.8)) 
C..........INPUT RUN SPECIFICATIONS.
 231          CONTINUE
              WRITE (iw,2306)
 2306         FORMAT (' ','Enter minimum value: ',$)
              READ (ir,*) rnum1(i)  !Read in starting value.
              WRITE (iw,2308)
 2308         FORMAT (' ','Enter maximum value: ',$)
              READ (ir,*) rnum2(i)  !Read in terminating value.
 232          CONTINUE
              WRITE (iw,2310)
 2310         FORMAT (' ','Enter increment: ',$)
              READ (ir,*) rnum3(i)  !Read in incremental value.
              IF (rnum3(i).eq.0.) THEN !Trouble with 0 division later on.
                WRITE (iw,2312)
 2312           FORMAT (' ','Zero increment not allowed: trouble with ',
     |                      'dividing by zero')
                GO TO 232
              END IF 
              WRITE (iw,2314) rnum1(i), rnum2(i), rnum3(i)  !Display input info.
 2314         FORMAT (' ','Run from ',1pe12.5,' to ',1pe12.5,
     |                    ' in increments of ',1pe12.5)
              WRITE (iw,2316)
 2316         FORMAT (' ','Confirm these values (Y or N): ',$)
              READ (ir,2301) lchose                         !Get confirmation.
 2301         FORMAT (a1)
              IF ((lchose.ne.'Y').and.(lchose.ne.'y')) GO TO 231
            END IF !((inum(i).lt.1).or.(inum(i).gt.8)) 
          END IF !(i.gt.jnum) 
        END DO !i = 1,3
        jnum = jnum-knum           !Number of valid loopings.
        IF (jnum.ne.0) THEN        !Run requested.
C..........WRITE OUT QUANTITY TO VARY, RUN SPECIFICATIONS.
          DO l = 1,jnum+knum       !Check all loopings.
            IF (inum(l).ne.0) THEN !Proper selection was made.
              WRITE (iw,2318) vtype(inum(l)),rnum1(l),      !Display run params.
     |                        rnum2(l), rnum3(l)
 2318         FORMAT (' ','Run ',a22,/,
     |                    '    from ',1pe12.5,' to ',1pe12.5,
     |                    ' in increments of ',1pe12.5)
C..........GET LOGS OF eta VALUES FOR LOGRITHMIC INCREMENTATION.
              IF (inum(l).eq.1) THEN  !Work with exponents for eta increments.
                rnum1(l) = log10(rnum1(l))
                rnum2(l) = log10(rnum2(l))
              END IF
            END IF
          END DO
C..........COMPUTE NUMBER OF RUNS FOR EACH LOOPING.
          DO l = 1,3
            lnum(l) = nint((rnum2(l)-rnum1(l)+rnum3(l))/rnum3(l))
          END DO
C..........DO MULTIPLE RUNS.
          WRITE (iw,2200)          !Inform user of beginning of computation.
          DO lnumb1 = 0,lnum(1)-1  !Outer loop.
            rnumb1 = rnum1(1)+float(lnumb1)*rnum3(1)  !Value of param for run.
            IF ((inum(1).ge.1).and.(inum(1).le.8)) THEN
              IF (inum(1).eq.1) THEN
                eta1 = 10**rnumb1    !Vary baryon-to-photon ratio.
              ELSE 
                qvary(inum(1)-1) = rnumb1  !Vary other quantities.
              END IF
            END IF
            DO lnumb2 = 0,lnum(2)-1  !Middle loop.
              rnumb2 = rnum1(2)+float(lnumb2)*rnum3(2)  !Value of param for run.
              IF ((inum(2).ge.1).and.(inum(2).le.8)) THEN
                IF (inum(2).eq.1) THEN
                  eta1 = 10**rnumb2  !Vary baryon-to-photon ratio.
                ELSE 
                  qvary(inum(2)-1) = rnumb2  !Vary other quantities.
                END IF
              END IF
              DO lnumb3 = 0,lnum(3)-1  !Inner loop.
                rnumb3 = rnum1(3)+float(lnumb3)*rnum3(3)  !Value of parameter.
                IF ((inum(3).ge.1).and.(inum(3).le.8)) THEN
                  IF (inum(3).eq.1) THEN
                    eta1 = 10**rnumb3  !Vary baryon-to-photon ratio.
                  ELSE 
                    qvary(inum(3)-1) = rnumb3  !Vary other quantities.
                  END IF
                END IF
                itime = 3
                CALL check         !Check interface subr before computation.
                CALL driver        !Do nucleosynthesis computation.
                itime = 8
                CALL check         !Check interface subroutine after computation
              END DO !lnumb3 = 0,lnum(3)-1  
            END DO !lnumb2 = 0,lnum(2)-1  
          END DO !lnumb1 = 0,lnum(1)-1  
          WRITE (iw,2202)          !Inform user of completion of computation.
        ELSE !(jnum.eq.0)
          WRITE (iw,2320)
 2320     FORMAT (' ','No selection made - ',
     |                'Press <RETURN> to continue: ',$)   
        END IF !(jnum.ne.0) 
        READ (ir,*)
        GO TO 300

C24--------EXIT SECTION---------------------------------------------------------

 240  CONTINUE
        RETURN

C30--------GO BACK TO MENU------------------------------------------------------

 300  CONTINUE
      GO TO 100

       END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE output

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - none

C----------REMARKS.
C     Outputs computational results either into an output file or onto 
C     the screen

C----------PARAMETERS.
      PARAMETER (ir=5)             !Input unit number (previous value = 1).
      PARAMETER (iw=6)             !Output unit number (previous value = 1).
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (itmax=40)         !Maximum # of line to be printed.

C----------COMMON AREAS.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags, counters.
      COMMON /outdat/ xout,thmout,t9out,tout,dtout,  !Output data.
     |                etaout,hubout
      COMMON /outopt/ nout,outfile                   !Output option.


C==========================DECLARATION DIVISION=================================

C----------COMPUTATION SETTINGS.
      REAL    cy                   !Time step limiting constant on abundances.
      REAL    ct                   !Time step limiting constant on temperature.
      REAL    t9i                  !Initial temperature (in 10**9 K).
      REAL    t9f                  !Final temperature (in 10**9 K).
      REAL    ytmin                !Smallest abundances allowed.

C----------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    c(3)                 !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C----------COUNTER.
      INTEGER it                   !# times accumulated in output buffer.

C----------OUTPUT ARRAYS.
      REAL    xout(itmax,nnuc)     !Nuclide mass fractions.
      REAL    thmout(itmax,6)      !Thermodynamic variables.
      REAL    t9out(itmax)         !Temperature (in units of 10**9 K).
      REAL    tout(itmax)          !Time.
      REAL    dtout(itmax)         !Time step.
      REAL    etaout(itmax)        !Baryon-to-photon ratio.
      REAL    hubout(itmax)        !Expansion rate.

C----------OUTPUT FILE STATUS.
      INTEGER nout                 !Number of output requests.
      LOGICAL outfile              !Indicates if output file used.

C----------USER INTERACTION VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION==================================

C10--------PRINT OUTPUT SELECTION AND AWAIT RESPONSE----------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY OUTPUT SELECTIONS.
      WRITE (iw,1000)
 1000 FORMAT (8(/),
     |        ' ',30x,'OUTPUT SELECTION',/,
     |        ' ',30x,'------ ---------',//,
     |        ' ',25x,' 1. REQUEST OUTPUT FILE',/,
     |        ' ',25x,' 2. REQUEST OUTPUT ON SCREEN',/,
     |        ' ',25x,' 3. EXIT',11(/),
     |        ' ',25x,' Enter selection (1-3): ',$)
C..........READ IN SELECTION NUMBER.
      READ (ir,1001) inum
 1001 FORMAT (i1)
C..........BRANCH TO APPROPRIATE SECTION.
      GO TO (200,300,400),inum
      GO TO 400                    !Improper input or <RETURN>.

C20--------REQUEST OUTPUT SECTION-----------------------------------------------

 200  CONTINUE
c      DO j = 1,it                  !Temperature in MeV.
c        t9out(j) = t9out(j)*.08617
c      END DO
c      DO j = 1,it                  !Energy density as fraction of total.
c        thmout(j,1) = thmout(j,1)/thmout(j,6)  !Rhog.
c        thmout(j,2) = thmout(j,2)/thmout(j,6)  !Rhoe.
c        thmout(j,3) = thmout(j,3)/thmout(j,6)  !Rhone.
c        thmout(j,4) = thmout(j,4)/thmout(j,6)  !Rhob.
c      END DO
C..........PRINT CAPTION.         
        nout = nout + 1            !Keep track of number of output requests.
        IF (nout.eq.1) THEN
          WRITE (2,2000)
 2000     FORMAT (54x,'NUCLIDE ABUNDANCE YIELDS',/,
     |            54x,'------- --------- ------',//)
        END IF
        WRITE (2,2002) cy,ct,t9i,t9f,ytmin
 2002   FORMAT (' Computational parameters:',/,
     |          '   cy = ',f5.3,'/  ct = ',f5.3,
     |          '/  initial temp = ',1pe8.2,
     |          '/  final temp = ',1pe8.2,
     |          '/  smallest abundances allowed = ',1pe8.2)
        WRITE (2,2004) c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
 2004   FORMAT (' Model parameters:',/,
     |          '   g = ',f5.2,'/  tau = ',f6.2,
     |          '/  # nu = ',f5.2,'/  lambda = ',1pe10.3,
     |          '/  xi-e = ',e10.3,'/  xi-m = ',e10.3,
     |          '/  xi-t = ',e10.3,/)
C..........PRINT HEADINGS, ABUNDANCES FOR NEUTRON TO LI8.         
        WRITE (2,2006)
 2006   FORMAT (4x,'Temp',8x,'N/H',10x,'P',10x,'D/H',9x,'T/H',8x,
     |          'He3/H',8x,'He4',8x,'Li6/H',7x,'Li7/H',7x,
     |          'Be7/H',6x,'Li8/H&up',/,132('-'))
        DO j = 1,it
          WRITE (2,2008) t9out(j),(xout(j,i),i=1,10)
 2008     FORMAT (1pe10.3,1p10e12.3)
        END DO
C..........PRINT THERMODYNAMIC QUANTITIES.         
        WRITE (2,2010)
 2010   FORMAT (' ',/,4x,'Temp',9x,'T',10x,'rhog',8x,'rhoe',7x,
     |              'rhone',8x,'rhob',8x,'phie',9x,'dt',9x,
     |              'eta',10x,'H',/,132('-')) 
        DO j = 1,it
          WRITE (2,2012) t9out(j),tout(j),(thmout(j,i),i=1,5),dtout(j),
     |                   etaout(j),hubout(j)
 2012     FORMAT (1pe10.3,9e12.3)
        END DO
        WRITE (2,2014)
 2014   FORMAT (///)
        outfile = .true.           !Output file requested.      
        WRITE (iw,2016)
 2016   FORMAT (' ','Output file requested - Press <RETURN> to '
     |              'continue: ',$)
        READ (ir,*)
        GO TO 500

C30--------REQUEST OUTPUT ON SCREEN SECTION-------------------------------------

C..........RETURN FROM LOOPING.
 300  CONTINUE
c      DO j = 1,it                  !Temperature in MeV.
c        t9out(j) = t9out(j)*.08617 
c      END DO
c      DO j = 1,it                  !Energy density as fraction of total.
c        thmout(j,1) = thmout(j,1)/thmout(j,6)  !Rhog.
c        thmout(j,2) = thmout(j,2)/thmout(j,6)  !Rhoe.
c        thmout(j,3) = thmout(j,3)/thmout(j,6)  !Rhone.
c        thmout(j,4) = thmout(j,4)/thmout(j,6)  !Rhob.
c      END DO
C..........DISPLAY SCREEN OUTPUT SELECTIONS.
        WRITE (iw,3000)
 3000   FORMAT (8(/),
     |          ' ',26x,'SCREEN OUTPUT SELECTION',/,
     |          ' ',26x,'------ ------ ---------',//,
     |          ' ',25x,' 1. DISPLAY D,T,HE3,HE4,LI7',/,
     |          ' ',25x,' 2. DISPLAY N,P,LI6,BE7,LI8&UP',/,
     |          ' ',25x,' 3. DISPLAY RHOG,RHOE,RHONE,RHOB',/,
     |          ' ',25x,' 4. DISPLAY T,DT,PHIE,ETA,H',/,
     |          ' ',25x,' 5. EXIT',9(/),
     |          ' ',25x,' Enter selection (1-5): ',$)
C..........READ IN SELECTION NUMBER.
        READ (ir,1001) inum
        GO TO (310,320,330,340,350),inum
        GO TO 350                  !Improper input or <RETURN>.
 310    CONTINUE                   !Display d,t,he3,he4,li7.
C..........PRINT CAPTION.
          WRITE (iw,2014)
          WRITE (iw,3100) cy,ct,t9i,t9f,ytmin
 3100     FORMAT (' ','Computational parameters:',/,
     |            ' ','   cy = ',f5.3,'/ ct = ',f5.3,
     |                '/ initial temp = ',1pe8.2,
     |                '/ final temp = ',1pe8.2,/,
     |            ' ','   smallest abundances allowed = ',1pe8.2)
          WRITE (iw,3102) c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
 3102     FORMAT (' ','Model parameters:',/,
     |            ' ','   g = ',f5.2,'/ tau = ',f6.2,
     |                '/ # nu = ',f5.2,'/ lambda = ',1pe10.3,/,
     |            ' ','   xi-e = ',e10.3,'/ xi-m = ',e10.3,
     |                '/ xi-t = ',e10.3,/)
C..........PRINT HEADINGS, ABUNDANCES FOR D,T,HE3,HE4,LI7.         
          WRITE (iw,3104)
 3104     FORMAT (4x,'Temp',8x,'D/H',9x,'T/H',8x,'He3/H',8x,
     |            'He4',8x,'Li7/H',/,' ',80('-'))
          DO j = 1,it
            WRITE (iw,3106) t9out(j),(xout(j,i),i=3,6),xout(j,8)
 3106       FORMAT (1pe10.3,1p5e12.3)
          END DO
          WRITE (iw,2014)
          WRITE (iw,3108)
 3108     FORMAT (' ','Press <RETURN> to continue: ',$)
          READ (ir,*)
          GO TO 360
 320    CONTINUE                   !Display n,p,li6,be7,li8&up.
C..........PRINT CAPTION.
          WRITE (iw,2014)
          WRITE (iw,3100) cy,ct,t9i,t9f,ytmin
          WRITE (iw,3102) c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
C..........PRINT HEADINGS, ABUNDANCES FOR N,P,LI6,BE7,LI8&UP.
          WRITE (iw,3204)
 3204     FORMAT (4x,'Temp',8x,'N/H',10x,'P',9x,
     |            'Li6/H',7x,'Be7/H',6x,'Li8/H&up',/,' ',80('-'))
          DO j = 1,it
            WRITE (iw,3106) t9out(j),(xout(j,i),i=1,2),xout(j,7),
     |                      (xout(j,i),i=9,10)
          END DO
          WRITE (iw,2014)
          WRITE (iw,3108)
          READ (ir,*)
          GO TO 360
 330    CONTINUE                   !Display rhog,rhoe,rhone,rhob.
C..........PRINT CAPTION.
          WRITE (iw,2014)
          WRITE (iw,3100) cy,ct,t9i,t9f,ytmin
          WRITE (iw,3102) c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
C..........PRINT ENERGY DENSITIES.
          WRITE (iw,3304)
 3304     FORMAT (4x,'Temp',8x,'rhog',8x,'rhoe',7x,'rhone',8x,'rhob',
     |            /,' ',80('-'))
          DO j = 1,it
            WRITE (iw,3306) t9out(j),(thmout(j,i),i=1,4)
 3306       FORMAT (1pe10.3,4e12.3)
          END DO
          WRITE (iw,2014)
          WRITE (iw,3108)
          READ (ir,*)
          GO TO 360
 340    CONTINUE                   !Display t,dt,phie,eta,hubcst.
C..........PRINT CAPTION.
          WRITE (iw,2014)
          WRITE (iw,3100) cy,ct,t9i,t9f,ytmin
          WRITE (iw,3102) c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
C..........PRINT THERMODYNAMIC QUANTITIES.         
          WRITE (iw,3404)
 3404     FORMAT (4x,'Temp',8x,'time',8x,'phie',9x,'dt',9x,'eta',10x,
     |            'H',/,' ',80('-'))
          DO j = 1,it
            WRITE (iw,3406) t9out(j),tout(j),thmout(j,5),dtout(j),
     |                      etaout(j),hubout(j)
 3406       FORMAT (1pe10.3,5e12.3)
          END DO
          WRITE (iw,2014)
          WRITE (iw,3108)
          READ (ir,*)
          GO TO 360
 350    CONTINUE                   !Exit.
          GO TO 500
 360    CONTINUE
        GO TO 300

C40--------EXIT SECTION---------------------------------------------------------

 400  CONTINUE
        RETURN

C50--------GO BACK TO MENU------------------------------------------------------

 500  CONTINUE
      GO TO 100

      END
