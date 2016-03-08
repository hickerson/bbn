#ifndef __REACTION_H__
#define __REACTION_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <map>

/*--------------------------------------------------------------------*/

/*#define DEBUG*/

/*--------------------------------------------------------------------*/
/* Prototypes */

/*
    enum ReactionIndex
    Using this way, reactions can added or removed 
    without interfering with loops.

    The naming convention for the reaction

        X -> x + Y, X + x -> y + Y, X + x + y -> z + 2Y etc...

    would be labeled
        
        X_xY, Xx_yY, X_xyY, Xx_yzY, Xxy_z2Y etc...
    
    where X, Y, Z, can be one of and atomic 
    isotope, such as H1, Li7, C12, etc...
    and x, y, z, are special cases for light 
    interchange isotopes can be one of
    g - gamma
    v - electron neutrino (nue)
    u - electron anti-neutrino (nuz)
    e - electron (beta-)
    z - positron (beta+)
    n - neutron (Nu1)
    p - proton (H1)
    d - deuteron (H2)
    t - triton (H3)
    h - helion (He3)
    a - alpha (He4)
*/
enum ReactionIndex {
	n_p = 961,    /// n <-> p
	H3_evHe3,   /// H3 -> e- + v + He3
	Li8_ev2He4, /// Li8 -> e- + v + 2He4
	B12_evC12,  /// B12 -> e- + v + C12
	C14_evN14,  /// C14 -> e- + v + N14
	B8_ev2He4,  /// B8 -> e+ + v + 2He4
	C11_evB11,  /// C11 -> e+ + v + B11
	N12_evC12,  /// N12 -> e+ + v + C12
	N13_evC13,  /// N13 -> e+ + v + C13
	O14_evN14,  /// O14 -> e+ + v + N14
	O15_evN15,  /// O15 -> e+ + v + N15
	H1n_gH2,    /// H1 + n -> g + H2
	H2n_gH3,    /// H2 + n -> g + H3
	He3n_gHe4,  /// He3 + n -> g + He4
	Li6n_gLi7,  /// Li6 + n -> g + Li7
	He3n_pH3,   /// He3 + n -> p + H3
	Be7n_pLi7,  /// Be7 + n -> p + Li7
	Li6n_tHe4,   /// Li6 + n -> t + He4
	Be7n_2He4,  /// Be7 + n -> 2He4
	H2p_gHe3,   /// H2 + p -> g + He3
	H3p_gHe4,   /// H3 + p -> g + He4
	Li6p_gBe7,  /// Li6 + p -> g + Be7
	Li6p_hHe4,  /// Li6 + p -> h + He4
	Li7p_2He4,  /// Li7 + p -> a + He4
	He4d_gLi6,  /// He4 + d -> g + Li6
	He4t_gLi7,  /// He4 + t -> g + Li7
	He4h_gBe7,  /// He4 + h -> g + Be7
	H2d_nHe3,   /// H2 + d -> n + He3   
	H2d_pH3,    /// H2 + d -> p + H3
	H3d_nHe4,   /// H3 + d -> n + He4
	He3d_pHe4,  /// He3 + d -> p + He4
	He3h_2pHe4, /// He3 + h -> 2p + He4
	Li7d_n2He4, /// Li7 + d -> n + 2He4
	Be7d_p2He4, /// Be7 + d -> p + 2He4
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
	Li8p_n2He4, /// Li8 + p -> n + 2He4
	B8n_p2He4,  /// B8 + n -> p + 2He4
	Be9p_d2He4, /// Be9 + p -> d + 2He4
	B11p_3He4,  /// B11 + p -> 3He4
	C11n_3He4,  /// C11 + n -> 3He4
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
    ReactionIndexOverflow
};

ReactionIndex operator++(ReactionIndex &,int);	// R++
//ReactionIndex operator--(ReactionIndex);
//ReactionIndex operator+(ReactionIndex a, int b);
/*
ReactionIndex &operator++(ReactionIndex index) {
	assert(index != ReactionIndexOverflow);
	index = static_cast<ReactionIndex>(index + 1);
	return index;
}
*/

#define REACMIN n_p
#define REACMAX C13a_nO16
#define REACBUF ReactionIndexOverflow

#define BETAMIN H3_evHe3
#define BETAMAX O15_evN15
#define NUMBETA (BETAMAX-BETAMIN+1)
#define BETABUF	(METAMAX+1) 


/********************************************
 *  NuclideIndex
 *  enum can be added or removed 
 *  with interfering with loops.
 *
 *  Nuclide naming convention:
 *  SyA - Sy is the element symbol of
 *        Nu, H, He, Li Be, B, C, N, O.
 *      - A is the atomic mass
 *  Nuclides: 
 *  0:  Nu0   1: Nu1,  
 *  2:  H1,   3: H2,   4: H3,   
 *  5:  He3,  6: He4, 
 *  7:  Li6,  8: Li7,  9: Be7, 10: Li8, 
 *  11: B8,  12: Be9, 13: B10, 14: B11, 
 *  15: C11, 16: B12, 17: C12, 18: N12, 
 *  19: C13, 20: N13, 21: C14, 22: N14, 
 *  23: O14, 24: N15, 25: O15, 26: O16 
 *******************************************/
enum NuclideIndex {
    //g=0, e=-1, v=0, 
    //n=1, p, 
    //d,   t,   h,   a, 
    NuclideIndexUnderflow=-2,
	None=-1,
    Nu0=30, Nu1, // Nu2, Nu4
    H1,  H2,  H3,  
    He3, He4, //He6, He8,
    Li6, Li7, Be7, Li8,
    B8,  Be9, B10, B11,
    C11, B12, C12, N12,
    C13, N13, C14, N14,
    O14, N15, O15, O16,
    NuclideIndexOverflow
};

#define NNUCREAC     	(ReactionIndexOverflow-n_p)
#define NNUC          	(NuclideIndexOverflow-Nu1)
//#define NNUC        	26 
//#define NNUCREAC   	88
//#define NBETA      	11
#define NUCBUF 			NuclideIndexOverflow
#define MinNuclideIndex Nu1
#define MaxNuclideIndex O16

//NuclideIndex & operator++(NuclideIndex &); 	/// prefix operator ++i
//NuclideIndex & operator--(NuclideIndex &); 	/// prefix operator --i
NuclideIndex operator++(NuclideIndex &, int); 	/// postfix operator i++
NuclideIndex operator--(NuclideIndex &, int); 	/// postfix operator i--
NuclideIndex operator+(NuclideIndex, NuclideIndex);
NuclideIndex operator-(NuclideIndex, NuclideIndex);
NuclideIndex operator+(NuclideIndex, int);
NuclideIndex operator-(NuclideIndex, int);
NuclideIndex operator!(NuclideIndex);		/// prefix operator !i



struct Nuclide {
    NuclideIndex id; 	/// Isotopic id
    //const char *S;    /// Symbol name
    int A;              /// Atomic number
    int Z;              /// Proton number
    int N;              /// Neutron number 
    double dm;          /// Mass excess [MeV]
};


struct Reaction {
    ReactionIndex id;
    int type;           /// TODO ReationType type;
    NuclideIndex in_major;
    NuclideIndex in_minor;
    NuclideIndex out_minor;
    NuclideIndex out_major;
    double reverse;
    double forward;

	//Reaction();
	//Reaction(ReactionIndex, int, //ReactionType, 
	//		 NuclideIndex, NuclideIndex, NuclideIndex, NuclideIndex, 
	//		 double, double);
	Reaction()
	: 
		id(ReactionIndexOverflow),
		type(-1),
		in_major(NuclideIndexOverflow),
		in_minor(NuclideIndexOverflow),
		out_minor(NuclideIndexOverflow),
		out_major(NuclideIndexOverflow),
		reverse(0),
		forward(0)
	{};

	Reaction(ReactionIndex index, int type, //ReactionType type, <<< TODO
			NuclideIndex in1, NuclideIndex in2, 
			NuclideIndex out1, NuclideIndex out2, 
			double reverse, double forward) 
	: 
		id(index),
		type(type),
		in_major(in1),
		in_minor(in2),
		out_minor(out1),
		out_major(out2),
		reverse(reverse),
		forward(forward)
	{};
};

typedef std::map<NuclideIndex,double> NuclideMap;
typedef std::map<NuclideIndex,Nuclide> NuclideList;
typedef std::map<ReactionIndex,double> ReactionMap;
typedef std::map<ReactionIndex,Reaction> ReactionList;

class ReactionNetwork {
	double T9;
	Reaction reaction[];
	double f[];
	double r[];
	int loop, inc, ip;
	double dt;
	double Y0[], Y[], dY_dt[];
	double H, rhob;

	int nucl();
	int linearize();
};

#endif  /// __REACTION_H__
