#ifndef __INCLUDE_H__
#define __INCLUDE_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <complex.h>
#include <string.h>
#include <strings.h>

/*--------------------------------------------------------------------*/

/*#define DEBUG*/

/*--------------------------------------------------------------------*/

#define pi      3.1415926535897932385
#define zeta3   1.2020569031595942855

#define hbar    6.58211889e-25    /// Planck constant in GeV.s
#define Gn      6.67428e-8        /// Gravitational constant in m^3.g^-1.s^-2
#define me      0.510998928       /// electron mass in MeV
#define kB      0.0861733238      /// Boltzmann constant in MeV/GK

/*--------------------------------------------------------------------*/

/* structure containing the cosmological model parameters */
// TODO make name not suck
typedef struct relicparam {
	int dummy;
	int entropy_model;
	double dd0,ndd,Tdend;
	double sd0,nsd,Tsend;
	double Sigmad0,nSigmad,TSigmaend;
	double nt0,nnt,Tnend;
	double mgravitino;
	double table_eff[276][3];
	double eta0;
	double nbnu;
	double life_neutron;
	double xinu1,xinu2,xinu3;
    double fierz;
} relicparam;

/*--------------------------------------------------------------------*/
/* Prototypes */

#define NNUCREAC    88
#define NNUC        26
#define NBETA       11

/*
    enum ReactionIndex
    Using this way, reactions can added or removed 
    without interfering with loops.

    The naming convention for the reaction

        X -> Y + Z, X + Y -> Z + W, etc...

    would be labeled
        
        X_YZ, XY_ZW, X_YZW, etc...
    
    where X, Y, Z, W, can be one g, v, e, n, p, d, t, h, a, 
    or, and atomic isotope, Xn, such as Li7, C12, etc...
    The special cases for light interchange isotopes are 
    g - gamma
    v - electron neutrino (or anti-)
    e - electron, positron, beta+/-
    n - neutron (Nu)
    p - proton (H1)
    d - deuteron (H2)
    t - triton (H3)
    h - helion (He3)
    a - alpha (He4)
*/
typedef enum ReactionIndex {
	n_p = 1,    /// n <-> p
	H3_evHe3,   /// H3 -> e- + v + He3
	Li8_evaHe3, /// Li8 -> e- + v + 2He4
	B12_evC12,  /// B12 -> e- + v + C12
	C14_evN14,  /// C14 -> e- + v + N14
	B8_evaHe3,  /// B8 -> e+ + v + 2He4
	C11_evB11,  /// C11 -> e+ + v + B11
	N12_evC12,  /// N12 -> e+ + v + C12
	N13_evC13,  /// N13 -> e+ + v + C13
	O14_evN14,  /// O14 -> e+ + v + N14
	O15_evN15,  /// O15 -> e+ + v + N15
	H1n_gH2,    /// H1 + n -> g + H2
	H2n_gH3,    /// H2 + n -> g + H3
	He3n_gHe3,  /// He3 + n -> g + He4
	Li6n_gLi7,  /// Li6 + n -> g + Li7
	He3n_pH3,   /// He3 + n -> p + H3
	Be7n_pLi7,  /// Be7 + n -> p + Li7
	Li6n_aH3,   /// Li6 + n -> a + H3
	Be7n_aHe3,  /// Be7 + n -> a + He4
	H2p_gHe3,   /// H2 + p -> g + He3
	H3p_gHe3,   /// H3 + p -> g + He4
	Li6p_gBe7,  /// Li6 + p -> g + Be7
	Li6p_aHe3,  /// Li6 + p -> a + He3
	Li7p_aHe3,  /// Li7 + p -> a + He4
	H2a_gLi6,   /// H2 + a -> p + Li6   changed to H2 + a -> g + Li6
	H3a_gLi7,   /// H3 + a -> p + Li7   changed to H3 + a -> g + Li7
	He3a_gBe7,  /// He3 + a -> p + Be7  changed to He3 + a -> g + Be7
	H2d_nHe3,   /// H2 + d -> p + He3   changed to H2 + d -> n + He3   
	H2d_pH3,    /// H2 + d -> n + H3    changed to H2 + d -> p + H3
	H3d_nHe3,   /// H3 + d -> n + He4
	He3d_pHe4,  /// He3 + d -> p + He4
	He3h_ppa,   /// He3 + He3 -> 2p + He4
	Li7d_naa,   /// Li7 + d -> n + a + He4
	Be7d_paa,   /// Be7 + d -> p + a + He4
	Li7n_gLi8,  /// Li7 + n -> g + Li8
	B10n_gB11,  /// B10 + n -> g + B11
	B11n_gB12,  /// B11 + n -> g + B12
	C11n_pB11,  /// C11 + n -> p + B11
	B10n_aLi7,  /// B10 + n -> a + Li7
	Be7p_gB8,   /// Be7 + p -> g + B8
	Be9p_gB10,  /// Be9 + p -> g + B10
	B10p_gC11,  /// B10 + p -> g + C11
	B11p_gC12,  /// B11 + p -> g + C12
	C11p_gN12,  /// C11 + p -> g + N12
	B12p_nC12,  /// B12 + p -> n + C12
	Be9p_aLi6,  /// Be9 + p -> a + Li6
	B10p_aBe7,  /// B10 + p -> a + Be7
	B12p_aBe9,  /// B12 + p -> a + Be9
	Li6a_gB10,  /// Li6 + a -> g + B10
	Li7a_gB11,  /// Li7 + a -> g + B11
	Be7a_gC11,  /// Be7 + a -> g + C11
	B8a_pC11,   /// B8 + a -> p + C11
	Li8a_nB11,  /// Li8 + a -> n + B11
	Be9a_nC12,  /// Be9 + a -> n + C12
	Be9d_nB10,  /// Be9 + d -> n + B10
	B10d_pB11,  /// B10 + d -> p + B11
	B11d_nC12,  /// B11 + d -> n + C12
	He4an_gBe9, /// He4 + a + n -> g + Be9
	He4aa_gC12, /// He4 + 2a -> g + C12
	Li8p_naHe4, /// Li8 + p -> n + a + He4
	B8n_paHe4,  /// B8 + n -> p + a + He4
	Be9p_daHe4, /// Be9 + p -> d + a + He4
	B11p_aaHe4, /// B11 + p -> 2a + Be4    changed to  B11 + p -> 2a + He4
	C11n_aaHe4, /// C11 + n -> 2a + He4
	C12n_gC13,  /// C12 + n -> g + C13
	C13n_gC14,  /// C13 + n -> g + C14
	N14n_gN15,  /// N14 + n -> g + N15
	N13n_pC13,  /// N13 + n -> p + C13
	N14n_pC14,  /// N14 + n -> p + C14
	O15n_pN15,  /// O15 + n -> p + N15
	O15n_aC12,  /// O15 + n -> a + C12
	C12p_gN13,  /// C12 + p -> g + N13
	C13p_gN14,  /// C13 + p -> g + N14
	C14p_gN15,  /// C14 + p -> g + N15
	N13p_gO14,  /// N13 + p -> g + O14
	N14p_gO15,  /// N14 + p -> g + O15
	N15p_gO16,  /// N15 + p -> g + O16
	N15p_aC12,  /// N15 + p -> a + C12
	C12a_gO16,  /// C12 + a -> g + O16
	B10a_pC13,  /// B10 + a -> p + C13
	B11a_pC14,  /// B11 + a -> p + C14
	C11a_pN14,  /// C11 + a -> p + N14
	N12a_pO15,  /// N12 + a -> p + O15
	N13a_pO16,  /// N13 + a -> p + O16
	B10a_nN13,  /// B10 + a -> n + N13
	B11a_nN14,  /// B11 + a -> n + N14
	B12a_nN15,  /// B12 + a -> n + N15
	C13a_nO16,  /// C13 + a -> n + O16
    REAC_MAX
} ReactionIndex;


typedef struct Reaction {
    ReactionIndex n;
    // TODO ReationType type;
    int type;
} Reaction;


/********************************************
 *  NuclideIndex
 *  enum can be added or removed 
 *  with interfering with loops.
 *
 *  Nuclide naming convention:
 *  g - gamma
 *  v - electron neutrino (or anti)
 *  e - electron, positron, or beta
 *  n - neutron 
 *  p - proton (H)
 *  D - deuteron (H2)
 *  T - triton (H3)
 *  h - helion (He3)
 *  a - alpha (He4)
 *  SA - S is the element symbol and 
 *
 *  Nuclides: 
 *  1:  n,    2: p, 
 *  3:  H2,   4: H3,   5: He3,  6: He4, 
 *  7:  Li6,  8: Li7,  9: Be7, 10: Li8, 
 *  11: B8,  12: Be9, 13: B10, 14: B11, 
 *  15: C11, 16: B12, 17: C12, 18: N12, 
 *  19: C13, 20: N13, 21: C14, 22: N14, 
 *  23: O14, 24: N15, 25: O15, 26: O16 
 *******************************************/
typedef enum NuclideIndex {
    //g=0, e=-1, v=0, 
    //n=1, p, 
    //d,   t,   h,   a, 
    Nu = 1,   H1,
    H2,  H3,  He3, He4,    
    Li6, Li7, Be7, Li8,
    B8,  Be9, B10, B11,
    C11, B12, C12, N12,
    C13, N13, C14, N14,
    O14, N15, O15, O16,
    NUC_MAX
} NuclideIndex;


typedef struct Nuclide {
    NuclideIndex i;     /// Isotopic index
    //const char *S;      /// Symbol name
    int A;              /// Atomic number
    int Z;              /// Proton number
    int N;              /// Neutron number 
    double dm;          /// Mass excess [MeV]
} Nuclide;


/* general.c */
double max(double x, double y);
double min(double x, double y);
double complex polylog(int n, int m, double x);
double Li2(double x);
double Li3(double x);
double Li4(double x);
double complex CLi2(double complex x);
double complex CLi3(double complex x);
double Cl2(double x);
double Cl3(double x);
double I0(double x);
double I1(double x);
double K0(double x);
double K1(double x);
double K2(double x);
double K3(double x);
double K4(double x);
double Lbessel(double x);
double Mbessel(double x);
double Nbessel(double x);
double K0exp(double x,double z);
double K1exp(double x,double z);
double K2exp(double x,double z);
double expcor(double x);
double kron(int x, int y);
double rand_gauss(void);
int test_integer(char name[]);
int test_file(char *name);


/* omega.c */
void Init_modeleff(int model_eff, struct relicparam* paramrelic);
double heff(double Temp, struct relicparam paramrelic);
double sgStar(double Temp, struct relicparam paramrelic);
double geff(double Temp, struct relicparam paramrelic);
void Init_cosmomodel(struct relicparam* paramrelic);
void Init_cosmomodel_param(double eta, double nbnu, double life_neutron, double xinu1, double xinu2, double xinu3, struct relicparam* paramrelic);
void Init_dark_density(double dd0, double ndd, double T_end, struct relicparam* paramrelic);
void Init_dark_entropy(double sd0, double nsd, double T_end, struct relicparam* paramrelic);
void Init_dark_entropySigmaD(double Sigmad0, double nSigmad, double T_end, struct relicparam* paramrelic);
void Init_nonthermal(double nt0, double nnt, double T_end, struct relicparam* paramrelic);
double dark_density(double T, struct relicparam paramrelic);
double dark_entropy(double T, struct relicparam paramrelic);
double dark_entropy_derivative(double T, struct relicparam paramrelic);
double dark_entropy_Sigmad(double T, struct relicparam paramrelic);
double nonthermal(double T, struct relicparam paramrelic);
double neutdens(double Tnu, struct relicparam paramrelic);
void Init_fierz(double eta, double nbnu, double life_neutron, double fierz, struct relicparam* paramrelic);


/* bbnrate.c */
void rate_weak(int err, double f[]);
void rate_pn(int err, struct relicparam paramrelic, double f[], double r[], double T9, double Tnu);
void rate_all(int err, double f[], double T9);


/* bbn.c */
// TODO void setup_reactions(Reaction[] reaction);
void setup_reactions(double reacparam[][8]);
//void setup_nuclides(Nuclide nuclide[]);
void setup_nuclides(int A[], int Z[], double dm[]);
int linearize(double T9, double reacparam[][8], double f[], double r[], int loop, int inc, int ip, double dt, double y0[], double y[], double dydt[], double H, double rhob);
int nucl(int err, struct relicparam paramrelic, double ratioH[]);
int nucl_failsafe(int err, struct relicparam paramrelic, double ratioH[]);
int nucl_witherrors(int err, struct relicparam paramrelic, double ratioH[], double sigma_ratioH[]);
int bbn_excluded(int err, struct relicparam paramrelic);

#endif  /// __INCLUDE_H__
