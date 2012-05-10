CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Changes (to run inder DEC unix f77): 
C     -----------------------------------
C     COMMON /therm/ -> COMMON /thermcb/
C     COMMON /rates/ f,r(nrec) -> COMMON /rates/ f,r
C
C     Default neutron lifetime 888.54 -> 885.7 (+- 0.8) - PDG 2002/3

      SUBROUTINE rate0

C----------LINKAGES.
C     CALLED BY - [subroutine] start
C     CALLS     - none

C----------REMARKS.
C     Generates weak decay rates.

C----------PARAMETER.
      PARAMETER (nrec=88)          !Number of nuclear reactions.

C----------COMMON AREA.
      COMMON /rates/  f,r    !Reaction rates.


C==========================DECLARATION DIVISION=================================

C----------REACTION RATE COEFFICIENTS.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      REAL    r(nrec)

C===========================PROCEDURE DIVISION==================================

C10--------SET DECAY RATE COEFFICIENTS------------------------------------------

C.......H3 -> e- + v + He3.........(Tilly-Weller-Hasan 1987)
      f(2)  = 1.79e-9

C.......Li8 -> e- + v + 2He4.......(Ajzenberg-Selove 1988)
      f(3)  = 8.27e-1

C.......B12 -> e- + B + C12........(Ajzenberg-Selove 1990)
      f(4)  = 3.43e+1

C.......C14 -> e- + v + N14........(Ajzenberg-Selove 1986)
      f(5)  = 3.834e-12

C.......B8 -> e+ + v + 2He4........(Ajzenberg-Selove 1988)
      f(6)  = 9.00e-1

C.......C11 -> e+ + v + B11........(Ajzenberg-Selove 1990)
      f(7)  = 5.668e-4

C.......N12 -> e+ + v + C12........(Ajzenberg-Selove 1990)
      f(8)  = 6.301e+1

C.......N13 -> e+ + v + C13........(Ajzenberg-Selove 1986)
      f(9)  = 1.159e-3

C.......O14 -> e+ + v + N14........(Ajzenberg-Selove 1986)
      f(10) = 9.8171e-3

C.......O15 -> e+ + v + N15........(Ajzenberg-Selove 1986)
      f(11) = 5.6704e-3

      RETURN

C----------REFERENCES-----------------------------------------------------------
C     Ajzenberg-Selove, F., 1990, Nucl. Phys. A506, 1.
C     Ajzenberg-Selove, F., 1988, Nucl. Phys. A490, 1.
C     Ajzenberg-Selove, F., 1986, Nucl. Phys. A449, 1.
C     Tilley, D.R., Weller, H.R., and Hasan, H.H., 1987, Nucl. Phys. A474, 1.

      END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE rate1(tph)

C----------LINKAGES.
C     CALLED BY - [subroutine] start, derivs
C     CALLS     - [function] xintd, eval

C----------REMARKS.
C     Generates rate coefficients for weak n->p and p->n reactions.

C----------PARAMETERS.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (iter=50)          !Number of gaussian quads.

C----------COMMON AREAS.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /modpr/  g,tau,xnu,c(3),cosmo,xi        !Model parameters.
      COMMON /thermcb/  thm,hubcst                     !Dynamic variables.
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Integration parameters.

C----------EXTERNAL FUNCTIONS.
      EXTERNAL func1               !Part 1 of n->p rate.
      EXTERNAL func2               !Part 2 of n->p rate.
      EXTERNAL func3               !Part 1 of p->n rate.
      EXTERNAL func4               !Part 2 of p->n rate.


C==========================DECLARATION DIVISION=================================

C----------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      REAL    r(nrec)              !Reverse reaction rate coefficients.

C----------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    tau                  !Neutron lifetime.
      REAL    xi(3)                !Neutrino degeneracy parameters.
     |                             !xi(1) is e neutrino degeneracy parameter.
     |                             !xi(2) is m neutrino degeneracy parameter.
     |                             !xi(3) is t neutrino degeneracy parameter.

C----------DYNAMIC VARIABLES.
      REAL    thm(14)              !Thermodynamic variables (energy densities).

C----------NEUTRINO PARAMETERS.
      REAL    t9mev                !Temperature (in units of MeV).
      REAL    tnmev                !Neutrino temperature (in units of MeV).

C----------LOCAL VARIABLES.
      REAL    tph                  !Photon temperature.
      REAL    w(2),x(2),           !Upper limits for exponentials, forward rate.
     |        y(2),z(2)            !Upper limits for exponentials, reverse rate.
      REAL    uplim1,uplim2,       !Upper limits for integrals for forward rate.
     |        uplim3,uplim4        !Upper limits for integrals for reverse rate.
      REAL    part1,part2,         !Parts of integrals for forward rate.
     |        part3,part4          !Parts of integrals for reverse rate.


C===========================PROCEDURE DIVISION==================================

C10--------COMPUTE WEAK REACTION RATES (NONDEGENERATE)--------------------------

C      IF (xi(1).eq.0.) THEN
C        f(1)  = thm(13)/tau        !Forward rate for weak np reaction.
C        r(1)  = thm(14)/tau        !Reverse rate for weak np reaction.
C      ELSE

C20--------COMPUTE WEAK REACTION RATES (DEGENERATE)-----------------------------

        t9mev = tph*.086171        !Convert photon temp to units of MeV.
        tnmev = tnu*.086171        !Convert neutrino temp to units of MeV.
C..........COMPUTE OVERFLOW LIMITS FOR LIMITS OF INTEGRATION (Ref 1 & 2).
        w(1) = (-(t9mev/.511)*(-88.722))
        w(2) = ((tnmev/.511)*(88.029+xi(1))+2.531)
        x(1) = ((t9mev/.511)*(88.029))
        x(2) = (-(tnmev/.511)*(-88.722+xi(1))-2.531)
        y(1) = (-(t9mev/.511)*(-88.722))
        y(2) = ((tnmev/.511)*(88.029-xi(1))-2.531)
        z(1) = ((t9mev/.511)*(88.029))
        z(2) = (-(tnmev/.511)*(-88.722-xi(1))+2.531)
C..........COMPARE LIMITS AND TAKE LARGER OF THE TWO.
        uplim1 = abs(w(1))
        uplim2 = abs(x(1))
        uplim3 = abs(y(1))
        uplim4 = abs(z(1))
        IF (uplim1.lt.abs(w(2))) uplim1 = w(2)
        IF (uplim2.lt.abs(x(2))) uplim2 = x(2)
        IF (uplim3.lt.abs(y(2))) uplim3 = y(2)
        IF (uplim4.lt.abs(z(2))) uplim4 = z(2)
C..........EVALUATE THE INTEGRALS NUMERICALLY.
        part1 = xintd(1.,uplim1,func1,iter)
        part2 = xintd(1.,uplim2,func2,iter)
        part3 = xintd(1.,uplim3,func3,iter)
        part4 = xintd(1.,uplim4,func4,iter)
        f(1) = part1 + part2       !Add 2 integrals to get forward rate.
        r(1) = part3 + part4       !Add 2 integrals to get reverse rate.
C      END IF !(xi(1).eq.0.)
      RETURN

C----------REFERENCES-----------------------------------------------------------
C     1) Forms of the integrals involved can be found in
C          Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
C          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
C
C     2) The overflow limit for the VAX/VMS system is exp(88.029).
C        The underflow limit for the VAX/VMS system is exp(-88.722).

      END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE rate2

C----------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [function] ex

C----------REMARKS.
C     Generates rate coefficients for reactions involving nuclides
C     up to A = 9.

C----------PARAMETER.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C----------COMMON AREAS.
      COMMON /rates/  f,r           !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y(nnuc)  !Evolution parameters.


C==========================DECLARATION DIVISION=================================

C----------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      real    r(nrec)

C----------EVOLUTION PARAMETER.
      REAL    t9                   !Temperature of photons (units of 10**9 K).


C===========================PROCEDURE DIVISION==================================

C10--------TEMPERATURE FACTORS--------------------------------------------------

      t913  = t9**(.33333333)      !t9**(1/3)
      t923  = t913*t913            !t9**(2/3)
      t943  = t923*t923            !t9**(4/3)
      t953  = t9*t923              !t9**(5/3)
      t912  = sqrt(t9)             !t9**(1/2)
      t932  = t9*t912              !t9**(3/2)
      t9m1  = 1/t9                 !t9**(-1)
      t9m23 = 1.0/t923             !t9**(-2/3)
      t9m32 = 1.0/t932             !t9**(-3/2)
      t9a   = t9/(1.0+13.076*t9)   !For reaction 17.
      t9a32 = t9a**(1.5)           !t9a**(3/2)
      t9b   = t9/(1.+49.18*t9)     !For reaction 18.
      t9b32 = t9b**(1.5)           !t9b**(3/2)
      IF (t9.gt.10.) THEN          !For reaction 22.
        t9c = 1.
      ELSE
        t9c = t9/(1.-9.69e-2*t9+2.84e-2*t953/(1.-9.69e-2*t9)**(2./3.))
      END IF
      t9c13 = t9c**(.3333333)      !t9c**(1/3)
      t9c56 = t9c**(.8333333)      !t9c**(5/6)
      t9d   = t9/(1.+0.759*t9)     !For reaction 24.
      t9d13 = t9d**(.3333333)      !t9d**(1/3)
      t9d56 = t9d**(.8333333)      !t9d**(5/6)
      t9e   = t9/(1.+0.1378*t9)    !For reaction 26.
      t9e13 = t9e**(.3333333)      !t9e**(1/3)
      t9e56 = t9e**(.8333333)      !t9e**(5/6)
      t9f   = t9/(1.+0.1071*t9)    !For reaction 27.
      t9f13 = t9f**(.3333333)      !t9f**(1/3)
      t9f56 = t9f**(.8333333)      !t9f**(5/6)


C20--------NEUTRON, PHOTON REACTIONS--------------------------------------------

C.......H(n,g)H2...................(Smith-Kawano-Malaney 1992)
      f(12)  = 4.742e+4*(1.-.8504*t912+.4895*t9-.09623*t932
     |                     +8.471e-3*t9*t9-2.80e-4*t9*t932)

C.......H2(n,g)H3..................(Wagoner 1969)
      f(13)  = 6.62e+1*(1.+18.9*t9)

C.......He3(n,g)He4................(Wagoner 1969)
      f(14)  = 6.62e+0*(1.+905.*t9)

C.......Li6(n,g)Li7................(Malaney-Fowler 1989)
      f(15)  = 5.10e+3

C30--------NEUTRON, PROTON REACTIONS--------------------------------------------

C.......He3(n,p)H3.................(Smith-Kawano-Malaney 1992)
      f(16)  = 7.21e+8*(1.-.508*t912+.228*t9)

C.......Be7(n,p)Li7................(Smith-Kawano-Malaney 1992)
      f(17)  = 2.675e+9*(1.-.560*t912+.179*t9-.0283*t932
     |        + 2.214e-3*t9*t9-6.851e-5*t9*t932)
     |        + 9.391e+8*t9a32*t9m32
     |        + 4.467e+7*t9m32*ex(-0.07486/t9)

C40--------NEUTRON, ALPHA REACTIONS---------------------------------------------

C.......Li6(n,a)H3.................(Caughlan-Fowler 1988)
      f(18)  = 2.54e+9*t9m32*ex(-2.39/t9)
     |         + 1.68e+8*(1.-.261*t9b32/t932)

C.......Be7(n,a)He4................(Wagoner 1969)
      f(19)  = 2.05e+4*(1.+3760.*t9)

C50--------PROTON, PHOTON REACTIONS---------------------------------------------

C.......H2(p,g)He3.................(Smith-Kawano-Malaney 1992)
      f(20)  = 2.65e+3*t9m23*ex(-3.720/t913)
     |         *(1.+.112*t913+1.99*t923+1.56*t9+.162*t943+.324*t953)

C.......H3(p,g)He4.................(Caughlan-Fowler 1988)
      f(21)  = 2.20e+4*t9m23*ex(-3.869/t913)
     |         *(1.+.108*t913+1.68*t923+1.26*t9+.551*t943+1.06*t953)

C.......Li6(p,g)Be7................(Caughlan-Fowler 1988)
      f(22)  = 6.69e+5*t9c56*t9m32*ex(-8.413/t9c13)

C60--------PROTON, ALPHA REACTIONS----------------------------------------------

C.......Li6(p,a)He3................(Caughlan-Fowler 1988)
      f(23)  = 3.73e+10*t9m23*ex(-8.413/t913-(t9/5.50)**2)
     |         *(1.+.050*t913-.061*t923-.021*t9+.006*t943+.005*t953)
     |         + 1.33e+10*t9m32*ex(-17.763/t9)
     |         + 1.29e+09*t9m1*ex(-21.820/t9)

C.......Li7(p,a)He4................(Smith-Kawano-Malaney 1992)
      f(24)  = 1.096e+9*t9m23*ex(-8.472/t913)
     |         - 4.830e+8*t9d56*t9m32*ex(-8.472/t9d13)
     |         + 1.06e+10*t9m32*ex(-30.442/t9)
     |         + 1.56e+5*t9m23*ex((-8.472/t913)-(t9/1.696)**2)
     |           *(1.+.049*t913-2.498*t923+.860*t9+3.518*t943+3.08*t953)
     |         + 1.55e+6*t9m32*ex(-4.478/t9)

C70--------ALPHA, PHOTON REACTIONS----------------------------------------------

C.......H2(a,g)Li6.................(Caughlan-Fowler 1988)
      f(25)  = 3.01e+01*t9m23*ex(-7.423/t913)
     |         *(1.+.056*t913-4.85*t923+8.85*t9-.585*t943-.584*t953)
     |         + 8.55e+1*t9m32*ex(-8.228/t9)

C.......H3(a,g)Li7.................(Smith-Kawano-Malaney 1992)
      f(26)  = 3.032e+5*t9m23*ex(-8.090/t913)
     |         *(1.+.0516*t913+.0229*t923+8.28e-3*t9
     |             -3.28e-4*t943-3.01e-4*t953)
     |         + 5.109e+5*t9e56*t9m32*ex(-8.068/t9e13)

C.......He3(a,g)Be7................(Smith-Kawano-Malaney 1992)
      f(27)  = 4.817e+6*t9m23*ex(-14.964/t913)
     |         *(1.+.0325*t913-1.04e-3*t923-2.37e-4*t9
     |             -8.11e-5*t943-4.69e-5*t953)
     |         + 5.938e+6*t9f56*t9m32*ex(-12.859/t9f13)

C80--------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------------

C.......H2(d,n)He3.................(Smith-Kawano-Malaney 1992)
      f(28)  = 3.95e+8*t9m23*ex(-4.259/t913)
     |         *(1.+.098*t913+.765*t923+.525*t9+9.61e-3*t943+.0167*t953)

C.......H2(d,p)H3..................(Smith-Kawano-Malaney 1992)
      f(29)  = 4.17e+8*t9m23*ex(-4.258/t913)
     |         *(1.+.098*t913+.518*t923+.355*t9-.010*t943-.018*t953)

C.......H3(d,n)He4.................(Smith-Kawano-Malaney 1992)
      f(30)  = 1.063e+11*t9m23*ex(-4.559/t913-(t9/.0754)**2)
     |         *(1.+.092*t913-.375*t923-.242*t9+33.82*t943+55.42*t953)
     |         + 8.047e+8*t9m23*ex(-0.4857/t9)

C.......He3(d,p)He4................(Smith-Kawano-Malaney 1992)
      f(31)  = 5.021e+10*t9m23*ex(-7.144/t913-(t9/.270)**2)
     |         *(1.+.058*t913+.603*t923+.245*t9+6.97*t943+7.19*t953)
     |         + 5.212e+8/t912*ex(-1.762/t9)

C90--------THREE PARTICLE REACTIONS---------------------------------------------

C.......He3(He3,2p)He4.............(Caughlan-Fowler 1988)
      f(32)  = 6.04e+10*t9m23*ex(-12.276/t913)
     |         *(1.+.034*t913-.522*t923-.124*t9+.353*t943+.213*t953)

C.......Li7(d,na)He4...............(Caughlan-Fowler 1988)
      f(33)  = 2.92e+11*t9m23*ex(-10.259/t913)

C.......Be7(d,pa)He4...............(Caughlan-Fowler 1988)
      f(34)  = 1.07e+12*t9m23*ex(-12.428/t913)

      RETURN

C----------REFERENCES-----------------------------------------------------------
C     Smith, M., Kawano, L.H., and Malaney, R.A., 1992, submitted to Ap. J.
C     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5.
C     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data 
C       Tables, 40, 283.
C     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.

      END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE rate3

C----------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [function] ex

C----------REMARKS.
C     Generates rate coefficients for reactions involving nuclides 
C     up to A = 18.

C----------PARAMETER.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C----------COMMON AREAS.
      COMMON /rates/  f,r           !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y(nnuc)  !Evolution parameters.


C==========================DECLARATION DIVISION=================================

C----------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      real    r(nrec)

C----------EVOLUTION PARAMETER.
      REAL    t9                   !Temperature of photons (units of 10**9 K).


C===========================PROCEDURE DIVISION==================================

C10--------TEMPERATURE FACTORS--------------------------------------------------

      t913  = t9**(.33333333)      !t9**(1/3)
      t923  = t913*t913            !t9**(2/3)
      t943  = t923*t923            !t9**(4/3)
      t953  = t9*t923              !t9**(5/3)
      t912  = sqrt(t9)             !t9**(1/2)
      t932  = t9*t912              !t9**(3/2)
      t915  = t9**(.2)             !t9**(1/5)
      t954  = t9**(1.25)           !t9**(5/4)
      t9m1  = 1.0/t9               !t9**(-1)
      t9m23 = 1.0/t923             !t9**(-2/3)
      t9m32 = 1.0/t932             !t9**(-3/2)
      t9m34 = sqrt(t9m32)          !t9**(-3/4)
      t9m15 = 1.0/t915             !t9**(-1/5)
      t9m54 = 1.0/t954             !t9**(-5/4)
      t9a   = t9/(1.+t9/15.1)      !For reaction 53.
      t9a13 = t9a**(.3333333)      !t9a**(1/3)
      t9a56 = t9a**(.8333333)      !t9a**(5/6)

C20--------NEUTRON, PHOTON REACTIONS--------------------------------------------

C.......Li7(n,g)Li8................(Wagoner 1969)
      f(35)  = 4.90e+3 + 9.96e+3*t9m32*ex(-2.62/t9)

C.......B10(n,g)B11................(Wagoner 1969)
      f(36)  = 6.62e+4

C.......B11(n,g)B12................(Malaney-Fowler 1989)
      f(37)  = 7.29e+2 + 2.40e+3*t9m32*ex(-0.223/t9)

C30--------NEUTRON, PROTON REACTIONS--------------------------------------------

C.......C11(n,p)B11................(Caughlan-Fowler 1988)
      f(38)  = 1.69e+8*(1.-.048*t912+.010*t9)

C40--------NEUTRON, ALPHA REACTIONS---------------------------------------------

C.......B10(n,a)Li7................(Caughlan-Fowler 1988)
      f(39)  = 5.07e+8

C50--------PROTON, PHOTON REACTIONS---------------------------------------------

C.......Be7(p,g)B8.................(Caughlan-Fowler 1988)
      f(40)  = 3.11e+5*t9m23*ex(-10.262/t913)
     |         + 2.53e+3*t9m32*ex(-7.306/t9)

C.......Be9(p,g)B10................(Caughlan-Fowler 1988)
      f(41)  = 1.33e+7*t9m23*ex(-10.359/t913-(t9/.846)**2)
     |         *(1.+.040*t913+1.52*t923+.428*t9+2.15*t943+1.54*t953)
     |         + 9.64e+4*t9m32*ex(-3.445/t9)
     |         + 2.72e+6*t9m32*ex(-10.620/t9)

C.......B10(p,g)C11................(Caughlan-Fowler 1988)
      f(42)  = 4.61e+5*t9m23*ex(-12.062/t913-(t9/4.402)**2)
     |         *(1.+.035*t913+.426*t923+.103*t9+.281*t943+.173*t953)
     |         + 1.93e+5*t9m32*ex(-12.041/t9)
     |         + 1.14e+4*t9m32*ex(-16.164/t9)

C.......B11(p,g)C12................(Caughlan-Fowler 1988)
      f(43)  = 4.62e+7*t9m23*ex(-12.095/t913-(t9/.239)**2)
     |         *(1.+.035*t913+3.00*t923+.723*t9+9.91*t943+6.07*t953)
     |         + 7.89e+3*t9m32*ex(-1.733/t9)
     |         + 9.68e+4*t9m15*ex(-5.617/t9)

C.......C11(p,g)N12................(Caughlan-Fowler 1988)
      f(44)  = 4.24e+4*t9m23*ex(-13.658/t913-(t9/1.627)**2)
     |         *(1.+.031*t913+3.11*t923+.665*t9+4.61*t943+2.50*t953)
     |         + 8.84e+3*t9m32*ex(-7.021/t9)

C60--------PROTON, NEUTRON REACTIONS--------------------------------------------

C.......B12(p,n)C12................(Wagoner 1969)
      f(45)  = 4.02e+11*t9m23*ex(-12.12/t913)

C70--------PROTON, ALPHA REACTIONS----------------------------------------------

C.......Be9(p,a)Li6................(Caughlan-Fowler 1988)
      f(46)  = 2.11e+11*t9m23*ex(-10.359/t913-(t9/.520)**2)
     |         *(1.+.040*t913+1.09*t923+.307*t9+3.21*t943+2.30*t953)
     |         + 4.51e+8*t9m1*ex(-3.046/t9)
     |         + 6.70e+8*t9m34*ex(-5.160/t9)

C.......B10(p,a)Be7................(Caughlan-Fowler 1988)
      f(47)  = 1.26e+11*t9m23*ex(-12.062/t913-(t9/4.402)**2)
     |         *(1.+.035*t913-.498*t923-.121*t9+.300*t943+.184*t953)
     |         + 2.59e+9*t9m1*ex(-12.260/t9)

C.......B12(p,a)Be9................(Wagoner 1969)
      f(48)  = 2.01e+11*t9m23*ex(-12.12/t913)

C80--------ALPHA, PHOTON REACTIONS----------------------------------------------

C.......Li6(a,g)B10................(Caughlan-Fowler 1988)
      f(49)  = 4.06e+6*t9m23*ex(-18.790/t913-(t9/1.326)**2)
     |         *(1.+.022*t913+1.54*t923+.239*t9+2.20*t943+.869*t953)
     |         + 1.91e+3*t9m32*ex(-3.484/t9)
     |         + 1.01e+4*t9m1*ex(-7.269/t9)

C.......Li7(a,g)B11................(Caughlan-Fowler 1988)
      f(50)  = 3.55e+7*t9m23*ex(-19.161/t913-(t9/4.195)**2)
     |         *(1.+.022*t913+.775*t923+.118*t9+.884*t943+.342*t953)
     |         + 3.33e+2*t9m32*ex(-2.977/t9)
     |         + 4.10e+4*t9m1*ex(-6.227/t9)

C.......Be7(a,g)C11................(Caughlan-Fowler 1988)
      f(51)  = 8.45e+7*t9m23*ex(-23.212/t913-(t9/4.769)**2)
     |         *(1.+.018*t913+.488*t923+.061*t9+.296*t943+.095*t953)
     |         + 1.25e+4*t9m32*ex(-6.510/t9)
     |         + 1.29e+5*t9m54*ex(-10.039/t9)

C90--------ALPHA, PROTON REACTIONS----------------------------------------------

C.......B8(a,p)C11.................(Wagoner 1969)
      f(52)  = 1.08e+15*t9m23*ex(-27.36/t913)

C100-------ALPHA, NEUTRON REACTIONS---------------------------------------------

C.......Li8(a,n)B11................(Malaney-Fowler 1989)
      f(53)  = 8.62e+13*t9a56*t9m32*ex(-19.461/t9a13)

C.......Be9(a,n)C12................(Caughlan-Fowler 1988)
      f(54)  = 4.62e+13*t9m23*ex(-23.870/t913-(t9/.049)**2)
     |         *(1.+.017*t913+8.57*t923+1.05*t9+74.51*t943+23.15*t953)
     |         + 7.34e-5*t9m32*ex(-1.184/t9)
     |         + 2.27e-1*t9m32*ex(-1.834/t9)
     |         + 1.26e+5*t9m32*ex(-4.179/t9)
     |         + 2.40e+8*ex(-12.732/t9)

C110-------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------------

C.......Be9(d,n)B10................(original Wagoner code)
      f(55)  = 7.16e+8*t9m23*ex(6.44-12.6/t913)

C.......B10(d,p)B11................(original Wagoner code)
      f(56)  = 9.53e+8*t9m23*ex(7.30-14.8/t913)

C.......B11(d,n)C12................(original Wagoner code)
      f(57)  = 1.41e+9*t9m23*ex(7.40-14.8/t913)

C120-------THREE PARTICLE REACTIONS---------------------------------------------

C.......He4(an,g)Be9...............(Caughlan-Fowler 1988)
      f(58)  = (2.59e-6/((1.+.344*t9)*t9**2))*ex(-1.062/t9)

C.......He4(2a,g)C12...............(Caughlan-Fowler 1988)
      f(59)  = 2.79e-8*t9m32*t9m32*ex(-4.4027/t9)
     |         + 1.35e-8*t9m32*ex(-24.811/t9)

C.......Li8(p,na)He4...............(original Wagoner code)
      f(60)  = 8.65e+9*t9m23*ex(-8.52/t913-(t9/2.53)**2)
     |         + 2.31e+9*t9m32*ex(-4.64/t9)

C.......B8(n,pa)He4................(original Wagoner code)
      f(61)  = 4.02e+8

C.......Be9(p,da)He4...............(Caughlan-Fowler 1988)
      f(62)  = 2.11e+11*t9m23*ex(-10.359/t913-(t9/.520)**2)
     |         *(1.+.040*t913+1.09*t923+.307*t9+3.21*t943+2.30*t953)
     |         + 5.79e+8*t9m1*ex(-3.046/t9)
     |         + 8.50e+8*t9m34*ex(-5.800/t9)

C.......B11(p,2a)He4...............(Caughlan-Fowler 1988)
      f(63)  = 2.20e+12*t9m23*ex(-12.095/t913-(t9/1.644)**2)
     |         *(1.+.034*t913+.140*t923+.034*t9+.190*t943+.116*t953)
     |         + 4.03e+6*t9m32*ex(-1.734/t9)
     |         + 6.73e+9*t9m32*ex(-6.262/t9)
     |         + 3.88e+9*t9m1*ex(-14.154/t9)

C.......C11(n,2a)He4...............(Wagoner 1969)
      f(64)  = 1.58e+8

      RETURN

C----------REFERENCES-----------------------------------------------------------
C     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5.
C     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data 
C       Tables, 40, 283.
C     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.

      END



C========================IDENTIFICATION DIVISION================================

      SUBROUTINE rate4

C----------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [function] ex

C----------REMARKS.
C     Generates rate coefficients for rest of reactions.

C----------PARAMETER.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C----------COMMON AREAS.
      COMMON /rates/  f,r           !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y(nnuc)  !Evolution parameters.


C==========================DECLARATION DIVISION=================================

C----------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      real    r(nrec)

C----------EVOLUTION PARAMETER.
      REAL    t9                   !Temperature of photons (units of 10**9 K).


C===========================PROCEDURE DIVISION==================================

C10--------TEMPERATURE FACTORS--------------------------------------------------

      t913  = t9**(.33333333)      !t9**(1/3)
      t923  = t913*t913            !t9**(2/3)
      t943  = t923*t923            !t9**(4/3)
      t953  = t9*t923              !t9**(5/3)
      t912  = sqrt(t9)             !t9**(1/2)
      t932  = t9*t912              !t9**(3/2)
      t935  = t9**(.6)             !t9**(3/5)
      t965  = t9**(1.2)            !t9**(6/5)
      t938  = t9**(.375)           !t9**(3/8)
      t9m13 = 1.0/t913             !t9**(1/3)
      t9m23 = 1.0/t923             !t9**(-2/3)
      t9m32 = 1.0/t932             !t9**(-3/2)
      t9m65 = 1.0/t965             !t9**(-6/5)
      t9a   = t9                   !For reaction 82.
     |          /(1.+4.78e-2*t9+7.56e-3*t953/(1.+4.78e-2*t9)**(2./3.))  
      t9a13 = t9a**(.33333333)     !t9a**(1/3)
      t9a56 = t9a**(.83333333)     !t9a**(5/6)
      t9b   = t9                   !For reaction 84.
     |          /(1.+7.76e-2*t9+2.64e-2*t953/(1.+7.76e-2*t9)**(2./3.))
      t9b13 = t9b**(.33333333)     !t9b**(1/3)
      t9b56 = t9b**(.83333333)     !t9b**(5/6)

C20--------NEUTRON, PHOTON REACTIONS--------------------------------------------

C.......C12(n,g)C13................(Wagoner 1969)
      f(65)  = 4.50e+2

C.......C13(n,g)C14................(Wagoner 1969)
      f(66)  = 1.19e+2 + 2.38e+5*t9m32*ex(-1.67/t9)

C.......N14(n,g)N15................(Wagoner 1969)
      f(67)  = 9.94e+3

C30--------NEUTRON, PROTON REACTIONS--------------------------------------------

C.......N13(n,p)C13................(Caughlan-Fowler 1988)
      f(68)  = 1.88e+8*(1.-.167*t912+.037*t9)

C.......N14(n,p)C14................(Caughlan-Fowler 1988)
      f(69)  = 2.39e+5*(1.+.361*t912+.502*t9)
     |         + 1.112e+8/t912*ex(-4.983/t9)

C.......O15(n,p)N15................(Caughlan-Fowler 1988)
      f(70)  = 3.50e+8*(1.+.452*t912-.191*t9)

C40--------NEUTRON, ALPHA REACTIONS---------------------------------------------

C.......O15(n,a)C12................(Caughlan-Fowler 1988)
      f(71)  = 3.50e+7*(1.+.188*t912+.015*t9)

C50--------PROTON, PHOTON REACTIONS---------------------------------------------

C.......C12(p,g)N13................(Caughlan-Fowler 1988)
      f(72)  = 2.04e+7*t9m23*ex(-13.690/t913-(t9/1.500)**2)
     |         *(1.+.030*t913+1.19*t923+.254*t9+2.06*t943+1.12*t953)
     |         + 1.08e+5*t9m32*ex(-4.925/t9)
     |         + 2.15e+5*t9m32*ex(-18.179/t9)

C.......C13(p,g)N14................(Caughlan-Fowler 1988)
      f(73)  = 8.01e+7*t9m23*ex(-13.717/t913-(t9/2.000)**2)
     |         *(1.+.030*t913+.958*t923+.204*t9+1.39*t943+.753*t953)
     |         + 1.21e+6*t9m65*ex(-5.701/t9)

C.......C14(p,g)N15................(Caughlan-Fowler 1988)
      f(74)  = 6.80e+6*t9m23*ex(-13.741/t913-(t9/5.721)**2)
     |         *(1.+.030*t913+.503*t923+.107*t9+.213*t943+.115*t953)
     |         + 5.36e+3*t9m32*ex(-3.811/t9)
     |         + 9.82e+4*t9m13*ex(-4.739/t9)

C.......N13(p,g)O14................(Caughlan-Fowler 1988)
      f(75)  = 4.04e+7*t9m23*ex(-15.202/t913-(t9/1.191)**2)
     |         *(1.+.027*t913-.803*t923-.154*t9+5.00*t943+2.44*t953)
     |         + 2.43e+5*t9m32*ex(-6.348/t9)

C.......N14(p,g)O15................(Caughlan-Fowler 1988)
      f(76)  = 4.90e+7*t9m23*ex(-15.228/t913-(t9/3.294)**2)
     |         *(1.+.027*t913-.778*t923-.149*t9+.261*t943+.127*t953)
     |         + 2.37e+3*t9m32*ex(-3.011/t9)
     |         + 2.19e+4*ex(-12.530/t9)

C.......N15(p,g)O16................(Caughlan-Fowler 1988)
      f(77)  = 9.78e+8*t9m23*ex(-15.251/t913-(t9/.450)**2)
     |         *(1.+.027*t913+.219*t923+.042*t9+6.83*t943+3.32*t953)
     |         + 1.11e+4*t9m32*ex(-3.328/t9)
     |         + 1.49e+4*t9m32*ex(-4.665/t9)
     |         + 3.80e+6*t9m32*ex(-11.048/t9)

C60--------PROTON, ALPHA REACTIONS----------------------------------------------

C.......N15(p,a)C12................(Caughlan-Fowler 1988)
      f(78)  = 1.08e+12*t9m23*ex(-15.251/t913-(t9/.522)**2)
     |         *(1.+.027*t913+2.62*t923+.501*t9+5.36*t943+2.60*t953)
     |         + 1.19e+8*t9m32*ex(-3.676/t9)
     |         + 5.41e+8/t912*ex(-8.926/t9)
     |         + 4.72e+7*t9m32*ex(-7.721/t9)
     |         + 2.20e+8*t9m32*ex(-11.418/t9)

C70--------ALPHA, PHOTON REACTIONS----------------------------------------------

C.......C12(a,g)O16................(Caughlan-Fowler 1988)
      f(79)  = 1.04e+8/t9**2*ex(-32.120/t913-(t9/3.496)**2)
     |         /(1.+.0489*t9m23)**2
     |         + 1.76e+8/(t9)**2/(1.+.2654*t9m23)**2*ex(-32.120/t913)
     |         + 1.25e+3*t9m32*ex(-27.499/t9)
     |         + 1.43e-2*(t9)**5*ex(-15.541/t9)

C80--------ALPHA, PROTON REACTIONS----------------------------------------------

C.......B10(a,p)C13................(Wagoner 1969)
      f(80)  = 9.60e+14*t9m23*ex(-27.99/t913)

C.......B11(a,p)C14................(Caughlan-Fowler 1988)
      f(81)  = 5.37e+11*t9m23*ex(-28.234/t913-(t9/0.347)**2)
     |         *(1.+.015*t913+5.575*t923+.576*t9+15.888*t943+4.174*t953)
     |         + 5.44e-3*t9m32*ex(-2.827/t9)
     |         + 3.36e+2*t9m32*ex(-5.178/t9)
     |         + 5.32e+6/t938*ex(-11.617/t9)

C.......C11(a,p)N14................(Caughlan-Fowler 1988)
      f(82)  = 7.15e+15*t9a56*t9m32*ex(-31.883/t9a13)

C.......N12(a,p)O15................(Caughlan-Fowler 1988)
      f(83)  = 5.59e+16*t9m23*ex(-35.60/t913)

C.......N13(a,p)O16................(Caughlan-Fowler 1988)
      f(84)  = 3.23e+17*t9b56*t9m32*ex(-35.829/t9b13)

C90--------ALPHA, NEUTRON REACTIONS---------------------------------------------

C.......B10(a,n)N13................(Caughlan-Fowler 1988)
      f(85)  = 1.20e+13*t9m23*ex(-27.989/t913-(t9/9.589)**2)

C.......B11(a,n)N14................(Caughlan-Fowler 1988)
      f(86)  = 6.97e+12*t9m23*ex(-28.234/t913-(t9/0.140)**2)
     |         *(1.+.015*t913+8.115*t923+.838*t9+39.804*t943
     |             +10.456*t953)
     |         + 1.79e+0*t9m32*ex(-2.827/t9)
     |         + 1.71e+3*t9m32*ex(-5.178/t9)
     |         + 4.49e+6*t935*ex(-8.596/t9)

C.......B12(a,n)N15................(Wagoner 1969)
      f(87)  = 3.04e+15*t9m23*ex(-28.45/t913)

C.......C13(a,n)O16................(Caughlan-Fowler 1988)
      f(88)  = 6.77e+15*t9m23*ex(-32.329/t913-(t9/1.284)**2)
     |         *(1.+.013*t913+2.04*t923+.184*t9)
     |         + 3.82e+5*t9m32*ex(-9.373/t9)
     |         + 1.41e+6*t9m32*ex(-11.873/t9)
     |         + 2.00e+9*t9m32*ex(-20.409/t9)
     |         + 2.92e+9*t9m32*ex(-29.283/t9)

      RETURN

C----------REFERENCES-----------------------------------------------------------
C     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data
C       Tables, 40, 283.
C     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.

      END



C========================IDENTIFICATION DIVISION================================

      BLOCK DATA   

C----------PARAMETERS.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C----------COMMON AREAS.
      COMMON /recpr0/ reacpr                         !Reaction parameter values.
      COMMON /compr0/ cy0,ct0,t9i0,t9f0,ytmin0,inc0  !Default comp parameters.
      COMMON /modpr0/ c0,cosmo0,xi0                  !Default model parameters.
      COMMON /varpr0/ dt0,eta0                       !Default variationl params.
      COMMON /nucdat/ am,zm,dm                       !Nuclide data.


C==========================DECLARATION DIVISION=================================

C----------REACTION PARAMETERS VALUES.
      REAL    reacpr(nrec,8)       !Reaction parameters.

C----------DEFAULT COMPUTATION PARAMETERS.
      REAL    cy0                  !Default time step limiting constant.
      REAL    ct0                  !Default time step limiting constant.
      REAL    t9i0                 !Default initial temperature (in 10**9 K).
      REAL    t9f0                 !Default final temperature (in 10**9 K).
      REAL    ytmin0               !Default smallest abundances allowed.
      INTEGER inc0                 !Default accumulation increment.

C----------DEFAULT MODEL PARAMETERS.
      REAL    c0(3)                !c0(1) is default variation of grav constant.
     |                             !c0(2) is default neutron half-life.
     |                             !c0(3) is default number of neutrinos.
      REAL    cosmo0               !Default cosmological constant.
      REAL    xi0(3)               !Default neutrino degeneracy parameters.

C----------DEFAULT VARIATIONAL PARAMETERS.
      REAL    dt0                  !Default initial time step.
      REAL    eta0                 !Default baryon-to-photon ratio.      

C----------NUCLIDE DATA.
      REAL    am(nnuc)             !Atomic number of nuclide.
      REAL    zm(nnuc)             !Charge of nuclide.
      REAL    dm(nnuc)             !Mass excess of nuclide.


C==============================DATA DIVISION====================================

C    Nuclide and corresponding number
C    --------------------------------
C    1) N         7) Li6      13) B10      19) C13      25) O15
C    2) P         8) Li7      14) B11      20) N13      26) O16
C    3) H2        9) Be7      15) C11      21) C14
C    4) H3       10) Li8      16) B12      22) N14
C    5) He3      11) B8       17) C12      23) O14
C    6) He4      12) Be9      18) N12      24) N15

C-----------NUCLIDE DATA.
      DATA am /1.,1.,2.,3.,3.,4.,6.,7.,7.,8.,8.,9.,10.,11.,11.,12.,
     |         12.,12.,13.,13.,14.,14.,14.,15.,15.,16./
      DATA zm /0.,1.,1.,1.,2.,2.,3.,3.,4.,3.,5.,4.,5.,5.,6.,5.,
     |         6.,7.,6.,7.,6.,7.,8.,7.,8.,8./
      DATA dm /.008665,.007825,.014102,.016050,.016030,.002603,.015125,
     |         .016004,.016929,.022487,.024609,.012186,.012939,.009305,
     |         .011432,.014354,.000000,.018641,.003354,.005738,.003242,
     |         .003074,.008597,.000108,.003070,-.005085/

C----------REACTION RATE COEFFICIENTS (Ref 1).
      DATA ((reacpr(i,j),j=1,8),i=1,11) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- --------- -------
     |             1.,1., 1.,0.,0., 2., 0.0  ,   0.0 ,     !N->P         
     |             2.,1., 4.,0.,0., 5., 0.0  ,   0.0 ,     !H3->He3
     |             3.,4.,10.,0.,0., 6., 0.0  ,   0.0 ,     !Li8->2He4
     |             4.,1.,16.,0.,0.,17., 0.0  ,   0.0 ,     !B12->C12
     |             5.,1.,21.,0.,0.,22., 0.0  ,   0.0 ,     !C14->N14
     |             6.,4.,11.,0.,0., 6., 0.0  ,   0.0 ,     !B8->2He4
     |             7.,1.,15.,0.,0.,14., 0.0  ,   0.0 ,     !C11->B11
     |             8.,1.,18.,0.,0.,17., 0.0  ,   0.0 ,     !N12->C12
     |             9.,1.,20.,0.,0.,19., 0.0  ,   0.0 ,     !N13->C13
     |            10.,1.,23.,0.,0.,22., 0.0  ,   0.0 ,     !O14->N14
     |            11.,1.,25.,0.,0.,24., 0.0  ,   0.0 /     !O15->N15
      DATA ((reacpr(i,j),j=1,8),i=12,22) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- --------- -------
     |            12.,2., 2.,1.,0., 3., 0.471,  25.82,     !H(n,g)H2
     |            13.,2., 3.,1.,0., 4., 1.63 ,  72.62,     !H2(n,g)H3
     |            14.,2., 5.,1.,0., 6., 2.61 , 238.81,     !He3(n,g)He4
     |            15.,2., 7.,1.,0., 8., 1.19 ,  84.17,     !Li6(n,g)Li7
     |            16.,3., 5.,1.,2., 4., 1.002,   8.863,    !He3(n,p)H3
     |            17.,3., 9.,1.,2., 8., 0.998,  19.081,    !Be7(n,p)Li7
     |            18.,3., 7.,1.,4., 6., 1.070,  55.494,    !Li6(n,a)H3
     |            19.,5., 9.,1.,0., 6., 4.70 , 220.39,     !Be7(n,a)He4
     |            20.,2., 3.,2.,0., 5., 1.63 ,  63.750,    !H2(p,g)He3
     |            21.,2., 4.,2.,0., 6., 2.61 , 229.932,    !H3(p,g)He4
     |            22.,2., 7.,2.,0., 9., 1.19 ,  65.054/    !Li6(p,g)Be7
      DATA ((reacpr(i,j),j=1,8),i=23,33) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- --------- -------
     |            23.,3., 7.,2.,5., 6., 1.07 ,  46.631,    !Li6(p,a)He3
     |            24.,5., 8.,2.,0., 6., 4.69 , 201.291,    !Li7(p,a)He4
     |            25.,2., 6.,3.,0., 7., 1.53 ,  17.118,    !H2(a,p)Li6
     |            26.,2., 6.,4.,0., 8., 1.11 ,  28.640,    !H3(a,p)Li7
     |            27.,2., 6.,5.,0., 9., 1.11 ,  18.423,    !He3(a,p)Be7
     |            28.,6., 3.,0.,1., 5., 1.73 ,  37.935,    !H2(d,p)He3
     |            29.,6., 3.,0.,2., 4., 1.73 ,  46.798,    !H2(d,n)H3
     |            30.,3., 4.,3.,1., 6., 5.54 , 204.117,    !H3(d,n)He4
     |            31.,3., 5.,3.,2., 6., 5.55 , 212.980,    !He3(d,p)He4
     |            32.,11.,5.,0.,2., 6., 3.39 , 149.230,    !He3(He3,2p)He4
     |            33.,9., 8.,3.,1., 6., 9.95 , 175.476/    !Li7(d,na)He4
      DATA ((reacpr(i,j),j=1,8),i=34,44) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- --------- -------
     |            34.,9., 9.,3.,2., 6., 9.97 , 194.557,    !Be7(d,pa)He4
     |            35.,2., 8.,1.,0.,10., 1.31 ,  23.59,     !Li7(n,g)Li8
     |            36.,2.,13.,1.,0.,14., 3.04 , 132.95,     !B10(n,g)B11
     |            37.,2.,14.,1.,0.,16., 2.34 ,  39.10,     !B11(n,g)B12
     |            38.,3.,15.,1.,2.,14., 1.002,  32.080,    !C11(n,p)B11
     |            39.,3.,13.,1.,6., 8., 0.758,  32.382,    !B10(n,a)Li7
     |            40.,2., 9.,2.,0.,11., 1.30 ,   1.595,    !Be7(p,g)B8
     |            41.,2.,12.,2.,0.,13., 0.973,  76.427,    !Be9(p,g)B10
     |            42.,2.,13.,2.,0.,15., 3.03 , 100.840,    !B10(p,g)C11
     |            43.,2.,14.,2.,0.,17., 7.01 , 185.173,    !B11(p,g)C12
     |            44.,2.,15.,2.,0.,18., 2.33 ,   6.975/    !C11(p,g)N12
      DATA ((reacpr(i,j),j=1,8),i=45,55) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- --------- -------
     |            45.,3.,16.,2.,1.,17., 3.00 , 146.08,     !B12(p,n)C12
     |            46.,3.,12.,2.,6., 7., 0.618,  24.674,    !Be9(p,a)Li6
     |            47.,3.,13.,2.,6., 9., 0.754,  13.301,    !B10(p,a)Be7
     |            48.,3.,16.,2.,6.,12., 0.292,  79.89,     !B12(p,a)Be9
     |            49.,2., 7.,6.,0.,13., 1.58 ,  51.753,    !Li6(a,g)B10
     |            50.,2., 8.,6.,0.,14., 4.02 , 100.538,    !Li7(a,g)B11
     |            51.,2., 9.,6.,0.,15., 4.02 ,  87.539,    !Be7(a,g)C11
     |            52.,3.,11.,6.,2.,15., 3.08 ,  86.00,     !B8(a,p)C11
     |            53.,3.,10.,6.,1.,14., 3.07 ,  76.96,     !Li8(a,n)B11   
     |            54.,3.,12.,6.,1.,17.,10.3  ,  66.160,    !Be9(a,n)C12
     |            55.,3.,12.,3.,1.,13., 2.07 ,  50.63/     !Be9(d,n)B10
      DATA ((reacpr(i,j),j=1,8),i=56,66) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- --------- -------
     |            56.,3.,13.,3.,2.,14., 6.44 , 107.13,     !B10(d,p)B11
     |            57.,3.,14.,3.,1.,17.,14.9  , 159.36,     !B11(d,n)C12
     |            58.,8., 6.,1.,0.,12., 0.584,  18.260,    !He4(an,g)Be9
     |            59.,7., 6.,0.,0.,17., 2.00 ,  84.420,    !He4(2a,g)C12
     |            60.,9.,10.,2.,1., 6., 3.58 , 177.73,     !Li8(p,na)He4
     |            61.,9.,11.,1.,2., 6., 3.58 , 218.82,     !B8(n,pa)He4
     |            62.,9.,12.,2.,3., 6., 0.807,   7.555,    !Be9(p,da)He4
     |            63.,10.,14.,2.,0.,6., 3.50 , 100.753,    !B11(p,2a)Be4
     |            64.,10.,15.,1.,0.,6., 3.49 , 132.83,     !C11(n,2a)He4
     |            65.,2.,17.,1.,0.,19., 0.886,  57.41,     !C12(n,g)C13
     |            66.,2.,19.,1.,0.,21., 3.58 ,  94.88/     !C13(n,g)C14
      DATA ((reacpr(i,j),j=1,8),i=67,77) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- --------- -------
     |            67.,2.,22.,1.,0.,24., 2.71 , 125.74,     !N14(n,g)N15
     |            68.,3.,20.,1.,2.,19., 1.002,  34.846,    !N13(n,p)C13
     |            69.,3.,22.,1.,2.,21., 3.003,   7.263,    !N14(n,p)C14
     |            70.,3.,25.,1.,2.,24., 1.002,  41.037,    !O15(n,p)N15
     |            71.,3.,25.,1.,6.,17., 0.709,  98.661,    !O15(n,a)C12
     |            72.,2.,17.,2.,0.,20., 0.884,  22.553,    !C12(p,g)N13
     |            73.,2.,19.,2.,0.,22., 1.19 ,  87.621,    !C13(p,g)N14
     |            74.,2.,21.,2.,0.,24., 0.900, 118.452,    !C14(p,g)N15
     |            75.,2.,20.,2.,0.,23., 3.57 ,  53.706,    !N13(p,g)O14
     |            76.,2.,22.,2.,0.,25., 2.70 ,  84.678,    !N14(p,g)O15
     |            77.,2.,24.,2.,0.,26., 3.62 , 140.734/    !N15(p,g)O16
      DATA ((reacpr(i,j),j=1,8),i=78,88) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- --------- -------
     |            78.,3.,24.,2.,6.,17., 0.706,  57.623,    !N15(p,a)C12
     |            79.,2.,17.,6.,0.,26., 5.13 ,  83.111,    !C12(a,g)O16
     |            80.,3.,13.,6.,2.,19., 9.36 ,  47.16,     !B10(a,p)C13  
     |            81.,3.,14.,6.,2.,21.,11.0  ,   9.098,    !B11(a,p)C14  
     |            82.,3.,15.,6.,2.,22., 3.68 ,  33.915,    !C11(a,p)N14  
     |            83.,3.,18.,6.,2.,25., 4.26 , 111.87,     !N12(a,p)O15  
     |            84.,3.,20.,6.,2.,26., 5.81 ,  60.557,    !N13(a,p)O16  
     |            85.,3.,13.,6.,1.,20., 9.34 ,  12.287,    !B10(a,n)N13  
     |            86.,3.,14.,6.,1.,22., 3.67 ,   1.835,    !B11(a,n)N14  
     |            87.,3.,16.,6.,1.,24., 4.25 ,  88.47,     !B12(a,n)N15  
     |            88.,3.,19.,6.,1.,26., 5.79 ,  25.711/    !C13(a,n)O16  

C----------DEFAULT COMPUTATION PARAMETERS.
      DATA cy0    /.300/           !Default time step limiting constant.
      DATA ct0    /.030/           !Default time step limiting constant.
      DATA t9i0   /1.00e+02/       !Default initial temperature.
      DATA t9f0   /1.00e-02/       !Default final temperature.
      DATA ytmin0 /1.00e-25/       !Default smallest abundances allowed.
      DATA inc0   /30/             !Default accumulation increment.
C-----------DEFAULT MODEL PARAMETERS.
      DATA c0     /1.00,885.7,3.0/ !Default variation of 3 parameters.
      DATA cosmo0 /0.00/           !Default cosmological constant.
      DATA xi0    /0.00,0.00,0.00/ !Default neutrino degeneracy parameter.

C-----------DEFAULT VARIATIONAL PARAMETERS.
      DATA dt0    /1.00e-04/       !Default initial time step.
      DATA eta0   /3.162e-10/      !Default baryon-to-photon ratio.

      END







