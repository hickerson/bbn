#include "include.h"

/*----------------------------------------------------
 * type: 0-10, each type has a unique (#n1,#n2,#n3,#n4) quartet
 * n1: incoming nuclide number 
 * n2: incoming light nuclide number
 * n3: outgoing light nuclide number
 * n4: outgoing nuclide number
 * rev: reverse reaction coefficient
 * q: energy release in reaction
 *----------------------------------------------------*/
void setup_reactions(double reacparam[][8]) 
{
	double _reacparam[NNUCREAC][8] =  // TODO do a real loop
	{
    /// reaction   type n1   n2  n3  n4   rev   Q [MeV]
        {np,        0,  n,   0,  0,  p,   0.0,   0.0},		// n <-> p
        {t_evh,     0,  t,   0,  0,  h,   0.0,   0.0},		// H3 -> e- + v + He3
        {Li8_evaa,  3,  Li8, 0,  0,  a,   0.0,   0.0},		// Li8 -> e- + v + 2He4
        {B12_evC12, 0,  B12, 0,  0,  C12, 0.0,   0.0},		// B12 -> e- + v + C12
        {C14_evN14, 0,  C14, 0,  0,  N14, 0.0,   0.0},		// C14 -> e- + v + N14
        {B8_evaa,   3,  B8,  0,  0,  a,   0.0,   0.0},		// B8 -> e+ + v + 2He4
        {C11_evB11, 0,  C11, 0,  0,  B11, 0.0,   0.0},		// C11 -> e+ + v + B11
        {N12_evC12, 0,  N12, 0,  0,  C12, 0.0,   0.0},		// N12 -> e+ + v + C12
        {N13_evC13, 0,  N13, 0,  0,  C13, 0.0,   0.0},		// N13 -> e+ + v + C13
        {O14_evN14, 0,  O14, 0,  0,  N14, 0.0,   0.0},		// O14 -> e+ + v + N14
        {O15_evN15, 0,  O15, 0,  0,  N15, 0.0,   0.0},		// O15 -> e+ + v + N15
        {pn_gd,     1,  p,   n,  0,  d,   0.471, 25.82},	// H + n -> g + H2
        {dn_gt,     1,  d,   n,  0,  t,   1.63,  72.62},	// H2 + n -> g + H3
        {hn_ga,     1,  h,   n,  0,  a,   2.61,  238.81},	// He3 + n -> g + He4
        {Li6n_gLi7, 1,  Li6, n,  0,  Li7, 1.19,  84.17},	// Li6 + n -> g + Li7
        {hn_pt,     2,  h,   n,  p,  t,   1.002, 8.863},	// He3 + n -> p + H3
        {Be7n_pLi7, 2,  Be7, n,  p,  Li7, 0.998, 19.081},	// Be7 + n -> p + Li7
        {Li6n_at,   2,  Li6, n,  t,  a,   1.070, 55.494},	// Li6 + n -> a + H3    /// TODO swap
        {Be7n_aa,   4,  Be7, n,  0,  a,   4.70,  220.39},	// Be7 + n -> 2He4
        {dp_gh,     1,  d,   p,  0,  h,   1.63,  63.750},	// H2 + p -> g + He3
        {tp_ga,     1,  t,   p,  0,  a,   2.61,  229.932},	// H3 + p -> g + He4
        {Li6p_gBe7, 1,  Li6, p,  0,  Be7, 1.19,  65.054},	// Li6 + p -> g + Be7
        {Li6p_ah,   2,  Li6, p,  h,  a,   1.07,  46.631},	// Li6 + p -> a + He3   /// TODO swap 
        {Li7p_aa,   4,  Li7, p,  0,  a,   4.69,  201.291},	// Li7 + p -> 2He4
        {da_gLi6,   1,  a,   d,  0,  Li6, 1.53,  17.118},	// H2 + a -> g + Li6    /// TODO swap id
        {ta_gLi7,   1,  a,   t,  0,  Li7, 1.11,  28.640},	// H3 + a -> g + Li7    /// TODO swap id
        {ha_gBe7,   1,  a,   h,  0,  Be7, 1.11,  18.423},	// He3 + a -> g + Be7    /// TODO swap id
        {dd_nh,     5,  d,   0,  n,  h,   1.73,  37.935},	// 2d -> n + He3
        {dd_pt,     5,  d,   0,  p,  t,   1.73,  46.798},	// 2d -> p + H3
        {td_na,     2,  t,   d,  n,  a,   5.54,  204.117},	// t + d -> n + He4
        {hd_pa,     2,  h,   d,  p,  a,   5.55,  212.980},	// He3 + d -> p + He4
        {hh_ppa,    10, h,   0,  p,  a,   3.39,  149.230},	// 2He3 -> 2p + He4
        {Li7d_naa,  8,  Li7, d,  n,  a,   9.95,  175.476},	// Li7 + d -> n + 2He4
        {Be7d_paa,  8,  Be7, d,  p,  a,   9.97,  194.557},	// Be7 + d -> p + 2He4
        {Li7n_gLi8, 1,  Li7, n,  0,  Li8, 1.31,  23.59},	// Li7 + n -> g + Li8
        {B10n_gB11, 1,  B10, n,  0,  B11, 3.04,  132.95},	// B10 + n -> g + B11
        {B11n_gB12, 1,  B11, n,  0,  B12, 2.34,  39.10},	// B11 + n -> g + B12
        {C11n_pB11, 2,  C11, n,  p,  B11, 1.002, 32.080},	// C11 + n -> p + B11
        {B10n_aLi7, 2,  B10, n,  a,  Li7, 0.758, 32.382},	// B10 + n -> a + Li7
        {Be7p_gB8,  1,  Be7, p,  0,  B8,  1.30,  1.595},	// Be7 + p -> g + B8
        {Be9p_gB10, 1,  Be9, p,  0,  B10, 0.973, 76.427},	// Be9 + p -> g + B10
        {B10p_gC11, 1,  B10, p,  0,  C11, 3.03,  100.840},	// B10 + p -> g + C11
        {B11p_gC12, 1,  B11, p,  0,  C12, 7.01,  185.173},	// B11 + p -> g + C12
        {C11p_gN12, 1,  C11, p,  0,  N12, 2.33,  6.975},	// C11 + p -> g + N12
        {B12p_nC12, 2,  B12, p,  n,  C12, 3.00,  146.08},	// B12 + p -> n + C12
        {Be9p_aLi6, 2,  Be9, p,  a,  Li6, 0.618, 24.674},	// Be9 + p -> a + Li6
        {B10p_aBe7, 2,  B10, p,  a,  Be7, 0.754, 13.301},	// B10 + p -> a + Be7
        {B12p_aBe9, 2,  B12, p,  a,  Be9, 0.292, 79.89},	// B12 + p -> a + Be9
        {Li6a_gB10, 1,  Li6, a,  0,  B10, 1.58,  51.753},	// Li6 + a -> g + B10
        {Li7a_gB11, 1,  Li7, a,  0,  B11, 4.02,  100.538},	// Li7 + a -> g + B11
        {Be7a_gC11, 1,  Be7, a,  0,  C11, 4.02,  87.539},	// Be7 + a -> g + C11
        {B8a_pC11,  2,  B8,  a,  p,  C11, 3.08,  86.00},	// B8 + a -> p + C11 
        {Li8a_nB11, 2,  Li8, a,  n,  B11, 3.07,  76.96},	// Li8 + a -> n + B11
        {Be9a_nC12, 2,  Be9, a,  n,  C12, 10.3,  66.160},	// Be9 + a -> n + C12
        {Be9d_nB10, 2,  Be9, d,  n,  B10, 2.07,  50.63},	// Be9 + d -> n + B10
        {B10d_pB11, 2,  B10, d,  p,  B11, 6.44,  107.13},	// B10 + d -> p + B11
        {B11d_nC12, 2,  B11, d,  n,  C12, 14.9,  159.36},	// B11 + d -> n + C12
        {aan_gBe9,  7,  a,   n,  0,  Be9, 0.584, 18.260},	// 2He4 + n -> g + Be9
        {aaa_gC12,  6,  a,   0,  0,  C12, 2.00,  84.420},	// 3He4 -> g + C12
        {Li8p_naa,  8,  Li8, p,  n,  a,   3.58,  177.73},	// Li8 + p -> n + 2He4
        {B8n_paa,   8,  B8,  n,  p,  a,   3.58,  218.82},	// B8 + n -> p + 2He4
        {Be9p_daa,  8,  Be9, p,  d,  a,   0.807, 7.555},	// Be9 + p -> d + 2He4
        {B11p_aaa,  9,  B11, p,  0,  a,   3.50,  100.753},	// B11 + p -> 3He4
        {C11n_aaa,  9,  C11, n,  0,  a,   3.49,  132.83},	// C11 + n -> 3He4
        {C12n_gC13, 1,  C12, n,  0,  C13, 0.886, 57.41},	// C12 + n -> g + C13
        {C13n_gC14, 1,  C13, n,  0,  C14, 3.58,  94.88},	// C13 + n -> g + C14
        {N14n_gN15, 1,  N14, n,  0,  N15, 2.71,  125.74},	// N14 + n -> g + N15
        {N13n_pC13, 2,  N13, n,  p,  C13, 1.002, 34.846},	// N13 + n -> p + C13
        {N14n_pC14, 2,  N14, n,  p,  C14, 3.003, 7.263},	// N14 + n -> p + C14
        {O15n_pN15, 2,  O15, n,  p,  N15, 1.002, 41.037},	// O15 + n -> p + N15
        {O15n_aC12, 2,  O15, n,  a,  C12, 0.709, 98.661},	// O15 + n -> a + C12
        {C12p_gN13, 1,  C12, p,  0,  N13, 0.884, 22.553},	// C12 + p -> g + N13
        {C13p_gN14, 1,  C13, p,  0,  N14, 1.19,  87.621},	// C13 + p -> g + N14
        {C14p_gN15, 1,  C14, p,  0,  N15, 0.900, 118.452},	// C14 + p -> g + N15
        {N13p_gO14, 1,  N13, p,  0,  O14, 3.57,  53.706},	// N13 + p -> g + O14
        {N14p_gO15, 1,  N14, p,  0,  O15, 2.70,  84.678},	// N14 + p -> g + O15
        {N15p_gO16, 2,  N15, p,  0,  O16, 3.62,  140.734},	// N15 + p -> g + O16
        {N15p_aC12, 2,  N15, p,  a,  C12, 0.706, 57.623},	// N15 + p -> a + C12
        {C12a_gO16, 1,  C12, a,  0,  O16, 5.13,  83.111},	// C12 + a -> g + O16
        {B10a_pC13, 2,  B10, a,  p,  C13, 9.36,  47.16},	// B10 + a -> p + C13
        {B11a_pC14, 2,  B11, a,  p,  C14, 11.0,  9.098},	// B11 + a -> p + C14
        {C11a_pN14, 2,  C11, a,  p,  N14, 3.68,  33.915},	// C11 + a -> p + N14
        {N12a_pO15, 3,  N12, a,  p,  O15, 4.26,  111.87},	// N12 + a -> p + O15
        {N13a_pO16, 3,  N13, a,  p,  O16, 5.81,  60.557},	// N13 + a -> p + O16
        {B10a_nN13, 2,  B10, a,  n,  N13, 9.34,  12.287},	// B10 + a -> n + N13
        {B11a_nN14, 2,  B11, a,  n,  N14, 3.67,  1.835},	// B11 + a -> n + N14
        {B12a_nN15, 2,  B12, a,  n,  N15, 4.25,  88.47},	// B12 + a -> n + N15
        {C13a_nO16, 2,  C13, a,  n,  O16, 5.79,  25.711}	// C13 + a -> n + O16
    };
    /*
    /// reaction   type n1 n2 n3 n4  rev   q
        {np,        0,  1.,0.,0.,2.,0.,0.},		// n -> p
        {t_evh,     0,  4.,0.,0.,5.,0.,0.},		// H3 -> e- + v + He3
        {Li8_evaa,  3,  10.,0.,0.,6.,0.,0.},		// Li8 -> e- + v + 2He4
        {B12_evC12, 0,  16.,0.,0.,17.,0.,0.},		// B12 -> e- + n + C12
        {C14_evN14, 0,  21.,0.,0.,22.,0.,0.},		// C14 -> e- + v + N14
        {B8_evaa,   3,  11.,0.,0.,6.,0.,0.},		// B8 -> e+ + v + 2He4
        {C11_evB11, 0,  15.,0.,0.,14.,0.,0.},		// C11 -> e+ + v + B11
        {N12_evC12, 0,  18.,0.,0.,17.,0.,0.},		// N12 -> e+ + v + C12
        {N13_evC13, 0,  20.,0.,0.,19.,0.,0.},		// N13 -> e+ + v + C13
        {O14_evN14, 0,  23.,0.,0.,22.,0.,0.},		// O14 -> e+ + v + N14
        {O15_evN15, 0,  25.,0.,0.,24.,0.,0.},		// O15 -> e+ + v + N15
        {pn_gd,     1,  2.,1.,0.,3.,0.471,25.82},	// H + n -> g + H2
        {dn_gt,     1,  3.,1.,0.,4.,1.63,72.62},	// H2 + n -> g + H3
        {hn_ga,     1,  5.,1.,0.,6.,2.61,238.81},	// He3 + n -> g + He4
        {Li6n_gLi7, 1,  7.,1.,0.,8.,1.19,84.17},	// Li6 + n -> g + Li7
        {hn_pt,     2,  5.,1.,2.,4.,1.002,8.863},	// He3 + n -> p + H3
        {Be7n_pLi7, 2,  9.,1.,2.,8.,0.998,19.081},	// Be7 + n -> p + Li7
        {Li6n_at,   2,  7.,1.,4.,6.,1.070,55.494},	// Li6 + n -> a + H3
        {Be7n_aa,   4,  9.,1.,0.,6.,4.70,220.39},	// Be7 + n -> a + He4
        {dp_gh,     1,  3.,2.,0.,5.,1.63,63.750},	// H2 + p -> g + He3
        {tp_ga,     1,  4.,2.,0.,6.,2.61,229.932},	// H3 + p -> g + He4
        {Li6p_gBe7, 1,  7.,2.,0.,9.,1.19,65.054},	// Li6 + p -> g + Be7
        {Li6p_ah,   2,  7.,2.,5.,6.,1.07,46.631},	// Li6 + p -> a + He3
        {Li7p_aa,   4,  8.,2.,0.,6.,4.69,201.291},	// Li7 + p -> a + He4
        {da_gLi6,   1,  6.,3.,0.,7.,1.53,17.118},	// H2 + a -> g + Li6
        {ta_gLi7,   1,  6.,4.,0.,8.,1.11,28.640},	// H3 + a -> g + Li7
        {ha_gBe7,   1,  6.,5.,0.,9.,1.11,18.423},	// He3 + a -> g + Be7
        {dd_nh,     5,  3.,0.,1.,5.,1.73,37.935},	// H2 + d -> n + He3
        {dd_pt,     5,  3.,0.,2.,4.,1.73,46.798},	// H2 + d -> p + H3
        {td_na,     2,  4.,3.,1.,6.,5.54,204.117},	// H3 + d -> n + He4
        {hd_pa,     2,  5.,3.,2.,6.,5.55,212.980},	// He3 + d -> p + He4
        {hh_ppa,    10, 5.,0.,2.,6.,3.39,149.230},	// He3 + He3 -> 2p + He4
        {Li7d_naa,  8,  8.,3.,1.,6.,9.95,175.476},	// Li7 + d -> n + a + He4
        {Be7d_paa,  8,  9.,3.,2.,6.,9.97,194.557},	// Be7 + d -> p + a + He4
        {Li7n_gLi8, 1,  8.,1.,0.,10.,1.31,23.59},	// Li7 + n -> g + Li8
        {B10n_gB11, 1,  13.,1.,0.,14.,3.04,132.95},	// B10 + n -> g + B11
        {B11n_gB12, 1,  14.,1.,0.,16.,2.34,39.10},	// B11 + n -> g + B12
        {C11n_pB11, 2,  15.,1.,2.,14.,1.002,32.080},	// C11 + n -> p + B11
        {B10n_aLi7, 2,  13.,1.,6.,8.,0.758,32.382},	// B10 + n -> a + Li7
        {Be7p_gB8,  1,  9.,2.,0.,11.,1.30,1.595},	// Be7 + p -> g + B8
        {Be9p_gB10, 1,  12.,2.,0.,13.,0.973,76.427},	// Be9 + p -> g + B10
        {B10p_gC11, 1,  13.,2.,0.,15.,3.03,100.840},	// B10 + p -> g + C11
        {B11p_gC12, 1,  14.,2.,0.,17.,7.01,185.173},	// B11 + p -> g + C12
        {C11p_gN12, 1,  15.,2.,0.,18.,2.33,6.975},	// C11 + p -> g + N12
        {B12p_nC12, 2,  16.,2.,1.,17.,3.00,146.08},	// B12 + p -> n + C12
        {Be9p_aLi6, 2,  12.,2.,6.,7.,0.618,24.674},	// Be9 + p -> a + Li6
        {B10p_aBe7, 2,  13.,2.,6.,9.,0.754,13.301},	// B10 + p -> a + Be7
        {B12p_aBe9, 2,  16.,2.,6.,12.,0.292,79.89},	// B12 + p -> a + Be9
        {Li6a_gB10, 1,  7.,6.,0.,13.,1.58,51.753},	// Li6 + a -> g + B10
        {Li7a_gB11, 1,  8.,6.,0.,14.,4.02,100.538},	// Li7 + a -> g + B11
        {Be7a_gC11, 1,  9.,6.,0.,15.,4.02,87.539},	// Be7 + a -> g + C11
        {B8a_pC11,  2,  11.,6.,2.,15.,3.08,86.00},	// B8 + a -> p + C11 
        {Li8a_nB11, 2,  10.,6.,1.,14.,3.07,76.96},	// Li8 + a -> n + B11
        {Be9a_nC12, 2,  12.,6.,1.,17.,10.3,66.160},	// Be9 + a -> n + C12
        {Be9d_nB10, 2,  12.,3.,1.,13.,2.07,50.63},	// Be9 + d -> n + B10
        {B10d_pB11, 2,  13.,3.,2.,14.,6.44,107.13},	// B10 + d -> p + B11
        {B11d_nC12, 2,  14.,3.,1.,17.,14.9,159.36},	// B11 + d -> n + C12
        {aan_gBe9,  7,  6.,1.,0.,12.,0.584,18.260},	// He4 + a + n -> g + Be9
        {aaa_gC12,  6,  6.,0.,0.,17.,2.00,84.420},	// He4 + 2a -> g + C12
        {Li8p_naa,  8,  10.,2.,1.,6.,3.58,177.73},	// Li8 + p -> n + a + He4
        {B8n_paa,   8,  11.,1.,2.,6.,3.58,218.82},	// B8 + n -> p + a + He4
        {Be9p_daa,  8,  12.,2.,3.,6.,0.807,7.555},	// Be9 + p -> d + a + He4
        {B11p_aaa,  9,  14.,2.,0.,6.,3.50,100.753},	// B11 + p -> 2a + He4
        {C11n_aaa,  9,  15.,1.,0.,6.,3.49,132.83},	// C11 + n -> 2a + He4
        {C12n_gC13, 1,  17.,1.,0.,19.,0.886,57.41},	// C12 + n -> g + C13
        {C13n_gC14, 1,  19.,1.,0.,21.,3.58,94.88},	// C13 + n -> g + C14
        {N14n_gN15, 1,  22.,1.,0.,24.,2.71,125.74},	// N14 + n -> g + N15
        {N13n_pC13, 2,  20.,1.,2.,19.,1.002,34.846},	// N13 + n -> p + C13
        {N14n_pC14, 2,  22.,1.,2.,21.,3.003,7.263},	// N14 + n -> p + C14
        {O15n_pN15, 2,  25.,1.,2.,24.,1.002,41.037},	// O15 + n -> p + N15
        {O15n_aC12, 2,  25.,1.,6.,17.,0.709,98.661},	// O15 + n -> a + C12
        {C12p_gN13, 1,  17.,2.,0.,20.,0.884,22.553},	// C12 + p -> g + N13
        {C13p_gN14, 1,  19.,2.,0.,22.,1.19,87.621},	// C13 + p -> g + N14
        {C14p_gN15, 1,  21.,2.,0.,24.,0.900,118.452},	// C14 + p -> g + N15
        {N13p_gO14, 1,  20.,2.,0.,23.,3.57,53.706},	// N13 + p -> g + O14
        {N14p_gO15, 1,  22.,2.,0.,25.,2.70,84.678},	// N14 + p -> g + O15
        {N15p_gO16, 2,  24.,2.,0.,26.,3.62,140.734},	// N15 + p -> g + O16
        {N15p_aC12, 2,  24.,2.,6.,17.,0.706,57.623},	// N15 + p -> a + C12
        {C12a_gO16, 1,  17.,6.,0.,26.,5.13,83.111},	// C12 + a -> g + O16
        {B10a_pC13, 2,  13.,6.,2.,19.,9.36,47.16},	// B10 + a -> p + C13
        {B11a_pC14, 2,  14.,6.,2.,21.,11.0,9.098},	// B11 + a -> p + C14
        {C11a_pN14, 2,  15.,6.,2.,22.,3.68,33.915},	// C11 + a -> p + N14
        {N12a_pO15, 3,  18.,6.,2.,25.,4.26,111.87},	// N12 + a -> p + O15
        {N13a_pO16, 3,  20.,6.,2.,26.,5.81,60.557},	// N13 + a -> p + O16
        {B10a_nN13, 2,  13.,6.,1.,20.,9.34,12.287},	// B10 + a -> n + N13
        {B11a_nN14, 2,  14.,6.,1.,22.,3.67,1.835},	// B11 + a -> n + N14
        {B12a_nN15, 2,  16.,6.,1.,24.,4.25,88.47},	// B12 + a -> n + N15
        {C13a_nO16, 2,  19.,6.,1.,26.,5.79,25.711}	// C13 + a -> n + O16
     */
    int reac,i,j;
    //TODO fix
    //ReactionIndex first = np;
    //ReactionIndex last = C13a_nO16;
    //ReactionIndex reac;
    //for (reac = 0; reac < NNUCREAC; reac++)
    for (j = 0; j < NNUCREAC; j++)
        for (i = 0; i < 8; i++)
        {
            int index = _reacparam[j][0];
            reacparam[index][i] = _reacparam[j][i];
        }
}

void setup_nuclides(Nuclide nuclide[]) {
//void setup_nuclides(int A[], int Z[], double D[]) {
	// Nuclides: 1=n, 2=p, 3=H2, 4=H3, 5=He3, 6=He4, 7=Li6, 8=Li7, 9=Be7, 10=Li8, 11=B8, 12=Be9, 13=B10, 14=B11, 15=C11, 16=B12, 17=C12, 18=N12, 19=C13, 20=N13, 21=C14, 22=N14, 23=O14, 24=N15, 25=O15, 26=O16 */
	
	//double Am[NNUC+1] = { 0., 1., 1., 2., 3., 3., 4., 6., 7., 7., 8., 8., 9., 10., 11., 11., 12., 12., 12., 13., 13., 14., 14., 14., 15., 15., 16.}; /// Mass number 

		
	//double Zm[NNUC+1] = {0., 0., 1., 1., 1., 2., 2., 3., 3., 4., 3., 5., 4., 5., 5., 6., 5., 6., 7., 6., 7., 6., 7., 8., 7., 8., 8.}; /// Atomic number Z
		
	//double Dm[NNUC+1] = {0., 8.071388, 7.289028, 13.135825, 14.949915, 14.931325, 2.424931, 14.0864, 14.9078, 15.7696, 20.9464, 22.9212, 11.34758, 12.05086, 8.6680, 10.6506, 13.3690, 0., 17.3382, 3.125036, 5.3455, 3.019916, 2.863440, 8.006521, 0.101439, 2.8554, -4.737036}; /// mass excess DeltaM
    /*
    A[n] = A[p] = 1;
    A[d] = 2;
    A[t] = A[h] = 3;
    A[a] = 4;
    A[Li6] = 6;
    */
    //double nuclide[NNUC_MAX][4];

    /*-----------------------------------
    /   Data from http://www.nndc.bnl.gov/masses/mass.mas03
    /
    /   S - isotopic symbol
    /   A - atomic mass number
    /   Z - atomic charge number
    /   N - neutron number
    /   M - mass excess 
    /------------------------------------
    /                 S   A   Z   N   D   */
	// Nuclides: 1=n, 2=p, 3=H2, 4=H3, 5=He3, 6=He4, 
	//double Dm[NNUC+1] = {0., 8.071388, 7.289028, 13.135825, 14.949915, 14.931325, 2.424931, 
    nuclide[0] =   { 0,  0,  0,  0,  0 };
    nuclide[n] =   { n,  1,  0,  1,  8.071388};
    nuclide[p] =   { p,  1,  1,  0,  7.289028};
    nuclide[d] =   { d,  1,  1,  1,  13.135825};
    nuclide[t] =   { t,  2,  1,  2,  14.949915};
    nuclide[h] =   { h,  3,  2,  1,  14.931325};
    nuclide[a] =   { a,  4,  2,  2,  2.424931};

	// 14.0864, 14.9078, 15.7696, 20.9464,
    nuclide[Li6] = {Li6, 6,  3,  3,  14.9078};
    nuclide[Li7] = {Li7, 7,  3,  4,  15.7696};
    nuclide[Be7] = {Be7, 7,  4,  3,  15.7696};
    nuclide[Li8] = {Li8, 8,  3,  5,  20.9464};

    // 22.9212, 11.34758, 12.05086, 8.6680,
    nuclide[B8] =  {B8,  8,  5,  3, 22.9212};
    nuclide[Be9] = {Be9, 9,  4,  5, 11.34758};
    nuclide[B10] = {B10, 10, 5,  5, 12.05086};
    nuclide[B11] = {B11, 11, 5,  6, 8.6680};

    // 10.6506, 13.3690, 0., 17.3382, 
    nuclide[C11] = {C11, 11, 6,  5, 10.6506};
    nuclide[B12] = {B12  12, 5,  7, 13.3690};
    nuclide[C12] = {C12, 12, 6,  6, 0};
    nuclide[N12] = {N12, 12, 7,  5, 17.3382};

    // 3.125036, 5.3455, 3.019916, 2.863440, 
    nuclide[C13] = {C13, 13, 6,  7, 3.125036};
    nuclide[N13] = {N13, 14, 7,  6, 5.3455};
    nuclide[C14] = {C14, 14, 6,  8, 3.019916};
    nuclide[N14] = {N14, 14, 7,  7, 2.863440};

    // 8.006521, 0.101439, 2.8554, -4.737036
    nuclide[O14] = {O14, 14, 8,  6, 8.006521};
    nuclide[N15] = {N15, 15, 7,  8, 0.101439};
    nuclide[O15] = {O15, 15, 8,  7, 2.8554};
    nuclide[O16] = {O16, 16, 8,  8, -4.737036};

    NuclideIndex i;
};

}
    
/*----------------------------------------------------*/


int linearize(double T9, double reacparam[][8], double f[], double r[], int loop, int inc, int ip, double dt, double Y0[], double Y[], double dY_dt[], double H, double rhob)
/* solves for new abundances using gaussian elimination with back substitution */
{
	/* Number of nuclides (#n1,#n2,#n3,#n4) for each of the 11 reaction types */
	double nn1[11]={1.,1.,1.,1.,1.,2.,3.,2.,1.,1.,2.};
	double nn2[11]={0.,1.,1.,0.,1.,0.,0.,1.,1.,1.,0.};
	double nn3[11]={0.,0.,1.,0.,0.,1.,0.,0.,1.,0.,2.};
	double nn4[11]={1.,1.,1.,2.,2.,1.,1.,1.,2.,3.,1.};

    ReactionIndex first = np;
    ReactionIndex last = C13a_nO16;
	
	//int i,j,k,l,n,i1,j1,ind;
	//double cn1,cn2,cn3,cn4,rn1,rn2,rn3,rn4,yY[NNUC];
	double cn1,cn2,cn3,cn4,yY[NNUC];
	cn1=cn2=cn3=cn4=0.;
	int fail;
	double bdln;
	int ierror;
	int c0 = 0;
	//double rev[last+1],q9[last+1];
	double a[NNUC+1][NNUC+1],b[NNUC+1],yx[NNUC+1];
	int icnvm;
	double x[NNUC+1], a0[NNUC+1][NNUC+1], cx, sum, xdy, t;
	int nord,test;
	
    /*
    ReactionIndex reac;
	for (reac = first; reac <= last; reac++) 
	{
        // TODO use struct instead of raw double array.
		rev[reac]=reacparam[reac][6];
		q9[reac]=reacparam[reac][7];
	}
    */
	
	NuclideIndex i,j,k,l;
	for(i = 1; i <= NNUC; i++) 
        for(j = 1; j <= NNUC; j++) 
            a[i][j] = 0;

    ReactionIndex n,i1,j1;
	for (n = first; n <= last; n++) 
	{
        int type = reacparam[n][1];
		i = reacparam[n][2];
		j = reacparam[n][3];
		k = reacparam[n][4];
		l = reacparam[n][5];
		double Rn = reacparam[n][6];
		double Q9 = reacparam[n][7];
        int rn1, rn2, rn3, rn4;
        rn1=nn1[type];
        rn2=nn2[type];
        rn3=nn3[type];
        rn4=nn4[type];
			
		if (i <= NNUC && l <= NNUC)
		{
			switch(type) 
			{
				case 0:	{ /* (1,0,0,1) type */
					cn1=f[n];
					cn2=0.;
					cn3=0.;
					cn4=r[n];
					break;}

				case 1: { /* (1,1,0,1) type */
					r[n]=Rn*1e10*pow(T9,1.5)*exp(-Q9/T9)*f[n];
					f[n]=rhob*f[n];
					cn1=Y[j]*f[n]/2.;
					cn2=Y[i]*f[n]/2.;
					cn3=0.;
					cn4=r[n];
					break;}

				case 2:	{ /* (1,1,1,1) type */
					f[n]=rhob*f[n];
					r[n]=Rn*exp(-Q9/T9)*f[n];
					cn1=Y[j]*f[n]/2.;
					cn2=Y[i]*f[n]/2.;
					cn3=Y[l]*r[n]/2.;
					cn4=Y[k]*r[n]/2.;
					break;}

				case 3:	{ /* (1,0,0,2) type */
					cn1=f[n];
					cn2=0.;
					cn3=0.;
					cn4=Y[l]*r[n]/2.;
					break;}

				case 4:	{ /* (1,1,0,2) type */
					f[n]=rhob*f[n];
					r[n]=Rn*exp(-Q9/T9)*f[n];
					cn1=Y[j]*f[n]/2.;
					cn2=Y[i]*f[n]/2.;
					cn3=0.;
					cn4=Y[l]*r[n]/2.;
					break;}

				case 5:	{ /* (2,0,1,1) type */
					f[n]=rhob*f[n];
					r[n]=Rn*exp(-Q9/T9)*f[n];
					cn1=Y[i]*f[n]/2.;
					cn2=0.;
					cn3=Y[l]*r[n]/2.;
					cn4=Y[k]*r[n]/2.;
					break;}

				case 6:	{ /* (3,0,0,1) type */
					r[n]=Rn*1.e20*pow(T9,1.5)*pow(T9,1.5)*exp(-Q9/T9)*f[n];
					f[n]=rhob*rhob*f[n];
					cn1=Y[i]*Y[i]*f[n]/6.;
					cn2=0.;
					cn3=0.;
					cn4=r[n];
					break;}
		
				case 7:	{ /* (2,1,0,1) type */
					r[n]=Rn*1.e20*pow(T9,1.5)*pow(T9,1.5)*exp(-Q9/T9)*f[n];
					f[n]=rhob*rhob*f[n];
					cn1=Y[j]*Y[i]*f[n]/3.;
					cn2=Y[i]*Y[i]*f[n]/6.;
					cn3=0.;
					cn4=r[n];
					break;}

				case 8:	{ /* (1,1,1,2) type */
					f[n]=rhob*f[n];
					r[n]=Rn*1.e-10*pow(T9,-1.5)*rhob*exp(-Q9/T9)*f[n];
					cn1=Y[j]*f[n]/2.;
					cn2=Y[i]*f[n]/2.;
					cn3=Y[l]*Y[l]*r[n]/6.;
					cn4=Y[k]*Y[l]*r[n]/3.;
					break;}

				case 9:	{ /* (1,1,0,3) type */
					f[n]=rhob*f[n];
					r[n]=Rn*1.e-10*pow(T9,-1.5)*rhob*exp(-Q9/T9)*f[n];
					cn1=Y[j]*f[n]/2.;
					cn2=Y[i]*f[n]/2.;
					cn3=0.;
					cn4=Y[l]*Y[l]*r[n]/6.;
					break;}

				case 10:{ /* (2,0,2,1) type */
					f[n]=rhob*f[n];
					r[n]=Rn*1.e-10*pow(T9,-1.5)*rhob*exp(-Q9/T9)*f[n];
					cn1=Y[i]*f[n]/2.;
					cn2=0.;
					cn3=Y[l]*Y[k]*r[n]/3.;
					cn4=Y[k]*Y[k]*r[n]/6.;}
			}

			i=NNUC+1-i;
			j=NNUC+1-j;
			k=NNUC+1-k;
			l=NNUC+1-l;
			
			if(j<=NNUC) 
                a[j][i]+=rn2*cn1;
			if(k<=NNUC) 
                a[k][i]-=rn3*cn1;
			a[i][i]+=rn1*cn1;
			a[l][i]-=rn4*cn1;
			
			if (j<=NNUC) 
			{
				a[j][j]+=rn2*cn2;
				if(k<=NNUC) 
                    a[k][j]-=rn3*cn2;
				a[i][j]+=rn1*cn2;
				a[l][j]-=rn4*cn2;
			}
			
			if (k<=NNUC)
			{
				if(j<=NNUC) 
                    a[j][k]-=rn2*cn3;
				a[k][k]+=rn3*cn3;
				a[i][k]-=rn1*cn3;
				a[l][k]+=rn4*cn3;
			}
			if(j<=NNUC) 
                a[j][l]-=rn2*cn4;
			if(k<=NNUC) 
                a[k][l]+=rn3*cn4;
			a[i][l]-=rn1*cn4;
			a[l][l]+=rn4*cn4;
		}
	}
	
	bdln=H*3.*1.e-5;
	
	for(i=1;i<=NNUC;i++)
	{
		i1=NNUC+1-i;
		for(j=1;j<=NNUC;j++)
		{
			j1=NNUC+1-j;
			if(fabs(a[j][i])<bdln*Y0[j1]/Y0[i1]) 
                a[j][i]=0.;
			else a[j][i]*=dt;
		}
		a[i][i]+=1.;
		b[i1]=Y0[i];
	}

	if(loop==1) 
        icnvm=ip; 
    else 
        icnvm=c0;
	
	nord=0;
	fail=0;
	for(i=1;i<=NNUC;i++)
	{
		x[i]=b[i];
		yx[i]=0.;
	}
	
	if(icnvm==inc) 
        for(i=1;i<=NNUC;i++) 
            for(j=1;j<=NNUC;j++) 
                a0[j][i]=a[j][i];

	for(i=1;i<=NNUC;i++)
	{
		if(a[i][i]==0.)
		{
			fail=i;
			return fail;
		}

		for(j=i+1;j<=NNUC;j++)
		{
			if(a[j][i]!=0.)
			{
				cx=a[j][i]/a[i][i];
				for(k=i+1;k<=NNUC;k++) 
                    a[j][k]-=cx*a[i][k];
				a[j][i]=cx;
				x[j]-=cx*x[i];
			}
		}
	}

	do
	{	x[NNUC]/=a[NNUC][NNUC];
		yx[NNUC]+=x[NNUC];
		
		for(i=NNUC-1;i>=1;i--)
		{
			sum=0.;
			for(j=i+1;j<=NNUC;j++) 
                sum+=a[i][j]*x[j];
			x[i]=(x[i]-sum)/a[i][i];
			yx[i]+=x[i];
		}

		test=1;
	
		if(icnvm==inc) 
		{
			for(i=1;i<=NNUC;i++) 
			{
				if(yx[i]!=0.) 
				{
					xdy=fabs(x[i]/yx[i]);
					if(xdy>2.e-4) 
					{
						if(nord<1) 
						{
							nord++;
							
							for(j=1;j<=NNUC;j++) 
							{
								t = 0.;
								for(k=1;k<=NNUC;k++) 
                                    t+=a0[j][k]*yx[k];
								x[j]=b[j]-t;
							}

							for(j=2;j<=NNUC;j++) 
                                for(k=j+1;k<=NNUC;k++) 
                                    x[k]-=a[k][j]*x[j];
							break;
						}
						else
						{
							fail=-1;
							ierror=i;
							return fail;
						}
					}
					else test=0;
				}
				else test=0;
			}
		}
		else test=0;
	}
	while(test);

	for(i=1;i<=NNUC;i++) 
	{
		yY[i]=yx[NNUC+1-i];
		dY_dt[i]=(yY[i]-Y0[i])/dt;
	}

#ifdef DEBUG
	if(fail!=0)
	{
		if(fail==-1) 
            printf("y(%d) failed to converge\n",ierror);
		if(fail>=1) 
            printf("%d th diagonal term equals zero\n",fail);
	}
#endif		
	return fail;
}

/*----------------------------------------------------*/

int nucl(int err, struct relicparam paramrelic, double ratioH[])
/* Main routine with computes the abundance ratios H2_H, ..., Be7_H as well as the baryon-to-photon ratio eta, using the parameters contained in paramrelic. The err parameter is a switch to choose if the central (err=0), high (err=1) or low (err=2) values of the nuclear rates is used. If (err) is negative, the lower value of only the nuclear rate number "-err" is used. If (err=4), the value of the nuclear rates is taken (gaussianly) randomly for a MC analysis. */
{
	int i;
    ReactionIndex first = np;
    ReactionIndex last = C13a_nO16;
	double f[last+1],r[last+1];
	for(i=0;i<=NNUC;i++) 
        ratioH[i]=0.;
	double sd;
	double rhod, sum_Y;
	double sum_dY_dt, sum_ZY, dsd_dT9, dphie_dT9, dlna3_dT9;
    double dphie_dlna3, dphie_dZY, sum_DeltaMdY_dt, sum_ZdY_dt;
	double cosh1, cosh2, cosh3, cosh4, cosh5, cosh6, cosh7;
    double sinh1, sinh2, sinh3, sinh4, sinh5, sinh6, sinh7;
	double T90,h_eta0,phie0;
	double dtl;
	int loop;
	double dh_dt, dphie_dt, dT9_dt, dlnT9_dt;
	double dT90_dt, dh_dt0, dphie_dt0;
	double dY_dt0[NNUC+1],dY_dt[NNUC+1],Y0[NNUC+1],Y[NNUC+1];
	double dtmin;
	double z;
	double H;
	dphie_dt0=dh_dt0=dT90_dt=phie0=h_eta0=T90=0.;

	/* Nuclides: 1=n, 2=p, 3=H2, 4=H3, 5=He3, 6=He4, 7=Li6, 8=Li7, 9=Be7, 10=Li8, 11=B8, 12=Be9, 13=B10, 14=B11, 15=C11, 16=B12, 17=C12, 18=N12, 19=C13, 20=N13, 21=C14, 22=N14, 23=O14, 24=N15, 25=O15, 26=O16 */
	
	double Am[NNUC+1] = {
        0., 1., 1., 2., 3., 3., 4., 6., 7., 7., 8., 8., 9., 10., 11., 11., 12., 12., 12., 13., 13., 14., 14., 14., 15., 15., 16.}; /// Mass number 
		
	double Zm[NNUC+1] = {0., 0., 1., 1., 1., 2., 2., 3., 3., 4., 3., 5., 4., 5., 5., 6., 5., 6., 7., 6., 7., 6., 7., 8., 7., 8., 8.}; /// Atomic number Z
		
	double Dm[NNUC+1] = {0., 8.071388, 7.289028, 13.135825, 14.949915, 14.931325, 2.424931, 14.0864, 14.9078, 15.7696, 20.9464, 22.9212, 11.34758, 12.05086, 8.6680, 10.6506, 13.3690, 0., 17.3382, 3.125036, 5.3455, 3.019916, 2.863440, 8.006521, 0.101439, 2.8554, -4.737036}; /// mass excess DeltaM

    int A[NNUC_SIZE];
    int Z[NNUC_SIZE];
    int D[NNUC_SIZE];
    setup_nuclides(A,Z,C);

	double reacparam[NNUCREAC+1][8];
    setup_reactions(reacparam);
		
    ReactionIndex n;
	for(n = first; n <= last; n++)
	{
		f[n] = 0;
		r[n] = 0;
	}
	
    /*
	double norm=1.;
	if(paramrelic.xinu1!=0. || b parametric)
	{
		double f_tmp[2],r_tmp[2];
		rate_pn(0,paramrelic,f_tmp,r_tmp,0.00001,0.00001);
		//norm=1./f_tmp[1]/paramrelic.life_neutron;
	}
    */

	double cy=0.1;
	double ct=0.01;
	double T9i=100.;
	double T9f=0.01;
	double Ytmin =1.e-30;
	int inc=50;
	double dt0=1.e-4;

	int ltime=0;
	int nitmax=1000;
	int is=1;
	int ip=inc;
	int it=0;
	int fail=0;

	if((err<0)||(err>100000))
	{
		cy=0.5;
		ct=0.1;
		dt0=1.e-2;
		nitmax=10;
	}

	double T9=T9i;
	double Tnu=T9;
	double t=1./pow(T9*0.09615,2.);
	double dt=dt0;

	if (15.011 / T9 > 58.)
	{
		Y[1] = 1.e-25;
		Y[2] = 1.;
	} 
	else if (15.011 / T9 < -58.)
	{
		Y[1] = 1.;
		Y[2] = 1.e-25;
	} 
	else 
	{
		Y[1] = 1. / (exp(15.011 / T9) + 1.);
		Y[2] = 1. / (exp(-15.011 / T9) + 1.);
	}

	Y0[1]=Y[1];
	Y0[2]=Y[2];

	z=5.929862032115561/T9;
	double h_eta=paramrelic.eta0*33683.*2.75;
	double phie=h_eta*(Y[2]*1.784e-5)/(pow(z,3.)*0.5*(Lbessel(z)-Lbessel(2.*z)*2.+Lbessel(3.*z)*3.-Lbessel(4.*z)*4.+Lbessel(5.*z)*5.-Lbessel(6.*z)*6.+Lbessel(7.*z)*7.));
	double rhob0=h_eta*pow(T9,3.);	
	double rho_gamma,drho_gamma_dT9,P_gamma,rho_epem;
	double drho_epem_dT9,drho_epem_dphie,P_epem,rho_neutrinos,rho_baryons;
	double dM_epem_dT9,dN_epem_dphie;
	
	Y[3] =Y[1]*Y[2]*rhob0*exp(25.82/T9)/(pow(T9,1.5)*4.71e9);
	
	Y0[3]=Y[3];
	for (i = 4; i <= NNUC; ++i) 
	{
		Y[i]=Ytmin;
		Y0[i]=Y[i];
	}
	rate_weak(err,f);



	while(ltime == 0)
	{
		for(loop=1;loop<=2;loop++)
		{
			rhod=dark_density(T9/11.605,paramrelic);		
			sd=dark_entropy(T9/11.605,paramrelic)/ 1.1605e10;
			dsd_dT9=dark_entropy_derivative(T9/11.605,paramrelic)/11.605/1.1605e10;

			z=5.929862032115561/T9;
			Tnu=pow(h_eta*T9*T9*T9/rhob0,1./3.)*T9i;
			if (phie<=17.)
			{
				cosh1=cosh(phie);
				cosh2=cosh(phie*2.);
				cosh3=cosh(phie*3.);
				cosh4=cosh(phie*4.);
				cosh5=cosh(phie*5.);
				cosh6=cosh(phie*6.);
				cosh7=cosh(phie*7.);
				sinh1=sinh(phie);
				sinh2=sinh(phie*2.);
				sinh3=sinh(phie*3.);
				sinh4=sinh(phie*4.);
				sinh5=sinh(phie*5.);
				sinh6=sinh(phie*6.);
				sinh7=sinh(phie*7.);
			}
			else
			{
				cosh1=0.;
				cosh2=0.;
				cosh3=0.;
				cosh4=0.;
				cosh5=0.;
				cosh6=0.;
				cosh7=0.;
				sinh1=0.;
				sinh2=0.;
				sinh3=0.;
				sinh4=0.;
				sinh5=0.;
				sinh6=0.;
				sinh7=0.;
			}

			rho_gamma=8.418*pow(T9,4.);
			
			drho_gamma_dT9=rho_gamma*4./T9;
			
			P_gamma=rho_gamma/3.;
				rho_epem=(Mbessel(z)*cosh1-Mbessel(2.*z)*cosh2+Mbessel(3.*z)*cosh3-Mbessel(4.*z)*cosh4+Mbessel(5.*z)*cosh5-Mbessel(6.*z)*cosh6+Mbessel(7.*z)*cosh7)*3206.; /* rho_e+ + rho_e- */
			
			drho_epem_dT9=z/T9*3206.*(Nbessel(z)*cosh1-Nbessel(2.*z)*2.*cosh2+Nbessel(3.*z)*3.*cosh3-Nbessel(4.*z)*4.*cosh4+Nbessel(5.*z)*5.*cosh5-Nbessel(6.*z)*6.*cosh6+Nbessel(7.*z)*7.*cosh7); /* d(rho_e+ + rho_e-)/d(T9) */
			
			drho_epem_dphie=(Mbessel(z)*sinh1-Mbessel(2.*z)*2.*sinh2+Mbessel(3.*z)*3.*sinh3-Mbessel(4.*z)*4.*sinh4+Mbessel(5.*z)*5.*sinh5-Mbessel(6.*z)*6.*sinh6+Mbessel(7.*z)*7.*sinh7)*3206.; /* d(rho_e+ + rho_e-)/d(phie) */
			
			P_epem=(Lbessel(z)*cosh1/z-Lbessel(2.*z)*cosh2/(z*2.)+Lbessel(3.*z)*cosh3/(z*3.)-Lbessel(4.*z)*cosh4/(z*4.)+Lbessel(5.*z)*cosh5/(z*5.)-Lbessel(6.*z)*cosh6/(z*6.)+Lbessel(7.*z)*cosh7/(z*7.))*3206.; /* P_e+ + P_e- */
			
			rho_neutrinos=12.79264*neutdens(Tnu,paramrelic);
										
			rho_baryons=h_eta*T9*T9*T9;
			dM_epem_dT9=-(z*z*z/T9)*(sinh1*(Lbessel(z)*3.-z*Mbessel(z))-sinh2*(Lbessel(2.*z)*3.-z*2.*Mbessel(2.*z))+sinh3*(Lbessel(3.*z)*3.-z*3.*Mbessel(3.*z))-sinh4*(Lbessel(4.*z)*3.-z*4.*Mbessel(4.*z))+sinh5*(Lbessel(5.*z)*3.-z*5.*Mbessel(5.*z))-sinh6*(Lbessel(6.*z)*3.-z*6.*Mbessel(6.*z))+sinh7*(Lbessel(7.*z)*3.-z*7.*Mbessel(7.*z))); /* d(pi^2 (hbar*c)^3 (ne- - ne+)*z^3 / 2(m c^2)^3) / d(T9) */
			
			dN_epem_dphie=z*z*z*(cosh1*Lbessel(z)-cosh2*2.*Lbessel(2.*z)+cosh3*3.*Lbessel(3.*z)-cosh4*4.*Lbessel(4.*z)+cosh5*5.*Lbessel(5.*z)-cosh6*6.*Lbessel(6.*z)+cosh7*7.*Lbessel(7.*z));
			if(dN_epem_dphie!=0.) dN_epem_dphie=1./dN_epem_dphie; /* d(pi^2/2 N_A (hbar*c/k)^3 h sum Z_i Y_i)/d(phie) */
			
			H=sqrt(Gn*8.*pi/3.*(rho_gamma+rho_epem+rho_neutrinos+rho_baryons+rhod));
			
			rate_pn(err,paramrelic,f,r,T9,Tnu);
            /*
			f[np]*=norm;
			r[np]*=norm;
		    */
            	
			rate_all(err,f,T9);
		
			fail=linearize(T9,reacparam,f,r,loop,inc,ip,dt,Y0,Y,dY_dt,H,rho_baryons);

			if( fail>0 ) 
                return 0;
			
			sum_Y=0.;
			sum_ZY=0.;
			sum_dY_dt=0.;
			sum_DeltaMdY_dt=0.;
			sum_ZdY_dt=0.;

			for (i=1;i<=NNUC;i++)
			{
				sum_Y+=Y[i];
				sum_ZY+=Zm[i]*Y[i];
				sum_dY_dt+=dY_dt[i];
				sum_DeltaMdY_dt+=Dm[i]*dY_dt[i];
				sum_ZdY_dt+=Zm[i]*dY_dt[i];
			}
		
			dphie_dT9=dN_epem_dphie*(-1.07e-4*h_eta*sum_ZY/T9-dM_epem_dT9);
			dphie_dlna3=-dN_epem_dphie*3.568e-5*h_eta*sum_ZY;
			dphie_dZY=dN_epem_dphie*3.568e-5*h_eta;
						dlna3_dT9=-(drho_gamma_dT9+drho_epem_dT9+drho_epem_dphie*dphie_dT9+rho_baryons*1.388e-4*sum_Y+T9*1e9*dsd_dT9)/(rho_gamma+P_gamma+rho_epem+P_epem+rho_baryons*(9.25e-5*T9*sum_Y+1.388e-4*T9*sum_dY_dt/(H*3.)+sum_DeltaMdY_dt/(H*3.))+T9*1.e9*sd+drho_epem_dphie*(dphie_dlna3+dphie_dZY*sum_ZdY_dt/(H*3.)));

			dT9_dt=3.*H/dlna3_dT9;
			dlnT9_dt=dT9_dt/T9;
			dh_dt=-h_eta*(H*3.+dlnT9_dt*3.);
			dphie_dt=dphie_dT9*dT9_dt+dphie_dlna3*(H*3.)+dphie_dZY*sum_ZdY_dt;
			
			if (T9 <= T9f || dt < fabs(1e-16 / dlnT9_dt) || ip == inc) 
			{
				it++;
				for (i=1;i<=NNUC;i++) 
                    ratioH[i]=Y[i]/Y[2];
			
				ratioH[2]=Y[2]*Am[2];
				ratioH[6]=Y[6]*Am[6];
				for(i=1;i<=9;i++) 
                    ratioH[10]+=ratioH[i];
				ratioH[10]-=1.;
				ratioH[0] = h_eta / 33683.;
				if((it==nitmax)||(ip<inc))
                    ltime = 1;
			}

			if(loop==1)
			{
				if(ip==inc) 
                    ip=0;
				ip++;
				is++;
				if(is>3)
				{
					dtmin=fabs(1./dlnT9_dt)*ct;
					for (i=1;i<=NNUC;i++)
					{
						if ((dY_dt[i]!=0.)&&(Y[i]>Ytmin)) 
						{
							dtl=(fabs(Y[i]/dY_dt[i]))*cy*(pow(log10(Y[i])/log10(Ytmin),2.)+1.);
							if (dtl<dtmin) dtmin=dtl;
						}
					}
					if (dtmin>dt*1.5) dtmin=dt*1.5;
					dt=dtmin;
				}
				t+=dt;
			
				T90=T9;
				h_eta0=h_eta;
				phie0=phie;
				
				dT90_dt=dT9_dt;
				dh_dt0=dh_dt;
				dphie_dt0=dphie_dt;
	
				T9=T90+dT90_dt*dt;
				h_eta=h_eta0+dh_dt0*dt;
				phie=phie0+dphie_dt0*dt;
				
				for (i=1;i<=NNUC;i++) 
				{
					Y0[i]=Y[i];
					dY_dt0[i]=dY_dt[i];
					Y[i]=Y0[i]+dY_dt0[i]*dt;
					if(Y[i]<Ytmin) Y[i]=Ytmin;
				}
			}
			else /* if(loop==2) */
			{
				T9=T90+(dT9_dt+dT90_dt)*0.5*dt;
				h_eta=h_eta0+(dh_dt+dh_dt0)*0.5*dt;
				phie=phie0+(dphie_dt+dphie_dt0)*0.5*dt;

				for (i=1;i<=NNUC;i++) 
				{
					Y[i]=Y0[i]+(dY_dt[i]+dY_dt0[i])*0.5*dt;
					if (Y[i]<Ytmin) Y[i]=Ytmin;
				}
			}
		}	
	}
			
	ratioH[8]+=ratioH[9];
	ratioH[5]+=ratioH[4];
	
	for (i=0;i<=NNUC;i++) 
        ratioH[i]=fabs(ratioH[i]);
	
	return fail;
}

/*----------------------------------------------------*/

int nucl_failsafe(int err, struct relicparam paramrelic, double ratioH[])
/* This routine is similar to nucl(...), the only difference is that it does not try to optimize the calculation time. */
{
    ReactionIndex first = np;
    ReactionIndex last = C13a_nO16;
	int i;
	for(i=0;i<=NNUC;i++) ratioH[i]=0.;
	double f[last+1],r[last+1];
	double sd;
	double rhod, sum_Y;
	double sum_dY_dt, sum_ZY, dsd_dT9, dphie_dT9, dlna3_dT9;
    double dphie_dlna3, dphie_dZY, sum_DeltaMdY_dt, sum_ZdY_dt;
	double cosh1, cosh2, cosh3, cosh4, cosh5, cosh6, cosh7;
    double sinh1, sinh2, sinh3, sinh4, sinh5, sinh6, sinh7;
	double T90,h_eta0,phie0;
	double dtl;
	int loop;
	double dh_dt, dphie_dt, dT9_dt, dlnT9_dt;
	double dT90_dt, dh_dt0, dphie_dt0;
	double dY_dt0[NNUC+1],dY_dt[NNUC+1],Y0[NNUC+1],Y[NNUC+1];
	double dtmin;
	double z;
	double H;
	dphie_dt0=dh_dt0=dT90_dt=phie0=h_eta0=T90=0.;

	/* Nuclides: 1=n, 2=p, 3=H2, 4=H3, 5=He3, 6=He4, 7=Li6, 8=Li7, 9=Be7, 10=Li8, 11=B8, 12=Be9, 13=B10, 14=B11, 15=C11, 16=B12, 17=C12, 18=N12, 19=C13, 20=N13, 21=C14, 22=N14, 23=O14, 24=N15, 25=O15, 26=O16 */

	/* Atomic mass number A */
	double Am[NNUC+1] = {0., 1., 1., 2., 3., 3., 4., 6., 7., 7., 8., 8., 9., 10., 11., 11., 12., 12., 12., 13., 13., 14., 14., 14., 15., 15., 16.}; 
	
    /* Atomic charge number Z */
	double Zm[NNUC+1] = {0., 0., 1., 1., 1., 2., 2., 3., 3., 4., 3., 5., 4., 5., 5., 6., 5., 6., 7., 6., 7., 6., 7., 8., 7., 8., 8.}; 
		
	double Dm[NNUC+1] = {0., 8.071388, 7.289028, 13.135825, 14.949915, 14.931325, 2.424931, 14.0864, 14.9078, 15.7696, 20.9464, 22.9212, 11.34758, 12.05086, 8.6680, 10.6506, 13.3690, 0., 17.3382, 3.125036, 5.3455, 3.019916, 2.863440, 8.006521, 0.101439, 2.8554, -4.737036}; /* mass excess DeltaM */	

	double reacparam[last+1][8];
    setup_reactions(reacparam);
		
	for(i = first; i <= last; i++)
	{
		f[i] = 0.;
		r[i] = 0.;
	}
	
    /* DONE in loop now
	double norm=1.;
	if(paramrelic.xinu1!=0.)
	{
		double f_tmp[2],r_tmp[2];
		rate_pn(0,paramrelic,f_tmp,r_tmp,0.00001,0.00001);
		norm=1./f_tmp[1]/paramrelic.life_neutron;
	}
    */

	double cy=0.1;
	double ct=0.01;
	double T9i=100.;
	double T9f=0.01;
	double Ytmin =1.e-30;
	int inc=50;
	double dt0=1.e-4;

	int ltime=0;
	int nitmax=1000;
	int is=1;
	int ip=inc;
	int it=0;
	int fail=0;

	double T9=T9i;
	double Tnu=T9;
	double t=1./pow(T9*0.09615,2.);
	double dt=dt0;

	if (15.011 / T9 > 58.)
	{
		Y[1] = 1.e-25;
		Y[2] = 1.;
	} 
	else if (15.011 / T9 < -58.)
	{
		Y[1] = 1.;
		Y[2] = 1.e-25;
	} 
	else 
	{
		Y[1] = 1. / (exp(15.011 / T9) + 1.);
		Y[2] = 1. / (exp(-15.011 / T9) + 1.);
	}

	Y0[1]=Y[1];
	Y0[2]=Y[2];

	z=5.929862032115561/T9;
	double h_eta=paramrelic.eta0*33683.*2.75;
	double phie=h_eta*(Y[2]*1.784e-5)/(pow(z,3.)*0.5*(Lbessel(z)-Lbessel(2.*z)*2.+Lbessel(3.*z)*3.-Lbessel(4.*z)*4.+Lbessel(5.*z)*5.-Lbessel(6.*z)*6.+Lbessel(7.*z)*7.));
	double rhob0=h_eta*pow(T9,3.);	
	double rho_gamma,drho_gamma_dT9,P_gamma,rho_epem;
	double drho_epem_dT9,drho_epem_dphie,P_epem,rho_neutrinos,rho_baryons;
	double dM_epem_dT9,dN_epem_dphie;
	
	Y[3] =Y[1]*Y[2]*rhob0*exp(25.82/T9)/(pow(T9,1.5)*4.71e9);
	
	Y0[3]=Y[3];
	for (i = 4; i <= NNUC; ++i) 
	{
		Y[i]=Ytmin;
		Y0[i]=Y[i];
	}
	rate_weak(err,f);



	while(ltime == 0)
	{
		for(loop=1;loop<=2;loop++)
		{
			rhod=dark_density(T9/11.605,paramrelic);		
			sd=dark_entropy(T9/11.605,paramrelic)/ 1.1605e10;
			dsd_dT9=dark_entropy_derivative(T9/11.605,paramrelic)/11.605/1.1605e10;

			z=5.929862032115561/T9;
			Tnu=pow(h_eta*T9*T9*T9/rhob0,1./3.)*T9i;
			if (phie<=17.)
			{
				cosh1=cosh(phie);
				cosh2=cosh(phie*2.);
				cosh3=cosh(phie*3.);
				cosh4=cosh(phie*4.);
				cosh5=cosh(phie*5.);
				cosh6=cosh(phie*6.);
				cosh7=cosh(phie*7.);
				sinh1=sinh(phie);
				sinh2=sinh(phie*2.);
				sinh3=sinh(phie*3.);
				sinh4=sinh(phie*4.);
				sinh5=sinh(phie*5.);
				sinh6=sinh(phie*6.);
				sinh7=sinh(phie*7.);
			}
			else
			{
				cosh1=0.;
				cosh2=0.;
				cosh3=0.;
				cosh4=0.;
				cosh5=0.;
				cosh6=0.;
				cosh7=0.;
				sinh1=0.;
				sinh2=0.;
				sinh3=0.;
				sinh4=0.;
				sinh5=0.;
				sinh6=0.;
				sinh7=0.;
			}

			rho_gamma=8.418*pow(T9,4.);
			
			drho_gamma_dT9=rho_gamma*4./T9;
			
			P_gamma=rho_gamma/3.;
				rho_epem=(Mbessel(z)*cosh1-Mbessel(2.*z)*cosh2+Mbessel(3.*z)*cosh3-Mbessel(4.*z)*cosh4+Mbessel(5.*z)*cosh5-Mbessel(6.*z)*cosh6+Mbessel(7.*z)*cosh7)*3206.; /* rho_e+ + rho_e- */
			
			drho_epem_dT9=z/T9*3206.*(Nbessel(z)*cosh1-Nbessel(2.*z)*2.*cosh2+Nbessel(3.*z)*3.*cosh3-Nbessel(4.*z)*4.*cosh4+Nbessel(5.*z)*5.*cosh5-Nbessel(6.*z)*6.*cosh6+Nbessel(7.*z)*7.*cosh7); /* d(rho_e+ + rho_e-)/d(T9) */
			
			drho_epem_dphie=(Mbessel(z)*sinh1-Mbessel(2.*z)*2.*sinh2+Mbessel(3.*z)*3.*sinh3-Mbessel(4.*z)*4.*sinh4+Mbessel(5.*z)*5.*sinh5-Mbessel(6.*z)*6.*sinh6+Mbessel(7.*z)*7.*sinh7)*3206.; /* d(rho_e+ + rho_e-)/d(phie) */
			
			P_epem=(Lbessel(z)*cosh1/z-Lbessel(2.*z)*cosh2/(z*2.)+Lbessel(3.*z)*cosh3/(z*3.)-Lbessel(4.*z)*cosh4/(z*4.)+Lbessel(5.*z)*cosh5/(z*5.)-Lbessel(6.*z)*cosh6/(z*6.)+Lbessel(7.*z)*cosh7/(z*7.))*3206.; /* P_e+ + P_e- */
			
			rho_neutrinos=12.79264*neutdens(Tnu,paramrelic);
										
			rho_baryons=h_eta*T9*T9*T9;
			dM_epem_dT9=-(z*z*z/T9)*(sinh1*(Lbessel(z)*3.-z*Mbessel(z))-sinh2*(Lbessel(2.*z)*3.-z*2.*Mbessel(2.*z))+sinh3*(Lbessel(3.*z)*3.-z*3.*Mbessel(3.*z))-sinh4*(Lbessel(4.*z)*3.-z*4.*Mbessel(4.*z))+sinh5*(Lbessel(5.*z)*3.-z*5.*Mbessel(5.*z))-sinh6*(Lbessel(6.*z)*3.-z*6.*Mbessel(6.*z))+sinh7*(Lbessel(7.*z)*3.-z*7.*Mbessel(7.*z))); /* d(pi^2 (hbar*c)^3 (ne- - ne+)*z^3 / 2(m c^2)^3) / d(T9) */
			
			dN_epem_dphie=z*z*z*(cosh1*Lbessel(z)-cosh2*2.*Lbessel(2.*z)+cosh3*3.*Lbessel(3.*z)-cosh4*4.*Lbessel(4.*z)+cosh5*5.*Lbessel(5.*z)-cosh6*6.*Lbessel(6.*z)+cosh7*7.*Lbessel(7.*z));
			if(dN_epem_dphie!=0.) dN_epem_dphie=1./dN_epem_dphie; /* d(pi^2/2 N_A (hbar*c/k)^3 h sum Z_i Y_i)/d(phie) */
			
			H=sqrt(Gn*8.*pi/3.*(rho_gamma+rho_epem+rho_neutrinos+rho_baryons+rhod));
			
			rate_pn(err,paramrelic,f,r,T9,Tnu);
			//f[1]*=norm;
			//r[1]*=norm;
			
			rate_all(err,f,T9);
		
			fail=linearize(T9,reacparam,f,r,loop,inc,ip,dt,Y0,Y,dY_dt,H,rho_baryons);

			if(fail>0) return 0;
			
			sum_Y=0.;
			sum_ZY=0.;
			sum_dY_dt=0.;
			sum_DeltaMdY_dt=0.;
			sum_ZdY_dt=0.;

			for (i=1;i<=NNUC;i++)
			{
				sum_Y+=Y[i];
				sum_ZY+=Zm[i]*Y[i];
				sum_dY_dt+=dY_dt[i];
				sum_DeltaMdY_dt+=Dm[i]*dY_dt[i];
				sum_ZdY_dt+=Zm[i]*dY_dt[i];
			}
		
			dphie_dT9=dN_epem_dphie*(-1.07e-4*h_eta*sum_ZY/T9-dM_epem_dT9);
			dphie_dlna3=-dN_epem_dphie*3.568e-5*h_eta*sum_ZY;
			dphie_dZY=dN_epem_dphie*3.568e-5*h_eta;
						dlna3_dT9=-(drho_gamma_dT9+drho_epem_dT9+drho_epem_dphie*dphie_dT9+rho_baryons*1.388e-4*sum_Y+T9*1e9*dsd_dT9)/(rho_gamma+P_gamma+rho_epem+P_epem+rho_baryons*(9.25e-5*T9*sum_Y+1.388e-4*T9*sum_dY_dt/(H*3.)+sum_DeltaMdY_dt/(H*3.))+T9*1.e9*sd+drho_epem_dphie*(dphie_dlna3+dphie_dZY*sum_ZdY_dt/(H*3.)));

			dT9_dt=3.*H/dlna3_dT9;
			dlnT9_dt=dT9_dt/T9;
			dh_dt=-h_eta*(H*3.+dlnT9_dt*3.);
			dphie_dt=dphie_dT9*dT9_dt+dphie_dlna3*(H*3.)+dphie_dZY*sum_ZdY_dt;
			
			if (T9 <= T9f || dt < fabs(1e-16 / dlnT9_dt) || ip == inc) 
			{
				it++;
				for (i=1;i<=NNUC;i++) ratioH[i]=Y[i]/Y[2];
			
				ratioH[2]=Y[2]*Am[2];
				ratioH[6]=Y[6]*Am[6];
				for(i=1;i<=9;i++) ratioH[10]+=ratioH[i];
				ratioH[10]-=1.;
				ratioH[0] = h_eta / 33683.;
				if((it==nitmax)||(ip<inc)) ltime = 1;
			}

			if(loop==1)
			{
				if(ip==inc) ip=0;
				ip++;
				is++;
				if(is>3)
				{
					dtmin=fabs(1./dlnT9_dt)*ct;
					for (i=1;i<=NNUC;i++)
					{
						if ((dY_dt[i]!=0.)&&(Y[i]>Ytmin)) 
						{
							dtl=(fabs(Y[i]/dY_dt[i]))*cy*(pow(log10(Y[i])/log10(Ytmin),2.)+1.);
							if (dtl<dtmin) dtmin=dtl;
						}
					}
					if (dtmin>dt*1.5) dtmin=dt*1.5;
					dt=dtmin;
				}
				t+=dt;
			
				T90=T9;
				h_eta0=h_eta;
				phie0=phie;
				
				dT90_dt=dT9_dt;
				dh_dt0=dh_dt;
				dphie_dt0=dphie_dt;
	
				T9=T90+dT90_dt*dt;
				h_eta=h_eta0+dh_dt0*dt;
				phie=phie0+dphie_dt0*dt;
				
				for (i=1;i<=NNUC;i++) 
				{
					Y0[i]=Y[i];
					dY_dt0[i]=dY_dt[i];
					Y[i]=Y0[i]+dY_dt0[i]*dt;
					if(Y[i]<Ytmin) Y[i]=Ytmin;
				}
			}
			else /* if(loop==2) */
			{
				T9=T90+(dT9_dt+dT90_dt)*0.5*dt;
				h_eta=h_eta0+(dh_dt+dh_dt0)*0.5*dt;
				phie=phie0+(dphie_dt+dphie_dt0)*0.5*dt;

				for (i=1;i<=NNUC;i++) 
				{
					Y[i]=Y0[i]+(dY_dt[i]+dY_dt0[i])*0.5*dt;
					if (Y[i]<Ytmin) Y[i]=Ytmin;
				}
			}
		}	
	}
			
	ratioH[8]+=ratioH[9];
	ratioH[5]+=ratioH[4];
	
	for (i=0;i<=NNUC;i++) ratioH[i]=fabs(ratioH[i]);
	
	return fail;
}

/*----------------------------------------------------*/

int nucl_witherrors(int err, struct relicparam paramrelic, double ratioH[], double sigma_ratioH[])
/* Routine which computes the abundance ratios (in ratioH[]) and their uncertainties (in sigma_ratioH[]), using the parameters contained in paramrelic. The err parameter is a switch to choose the evaluation error method (0=no error, 1=high values of the nuclear rates, 2=low values, 3=linear error calculation). */
{	
    ReactionIndex first = np;
    ReactionIndex last = C13a_nO16;
	int ie,je;
	for(ie=0;ie<=NNUC;ie++) 
        ratioH[ie]=sigma_ratioH[ie]=0.;

	if(err==0)
	{
		if(nucl(0,paramrelic, ratioH)>0) 
            return 0; 
		else 
            return 1;
	}
	else if(err==1||err==2)
	{
		if(nucl(err,paramrelic, sigma_ratioH)>0) 
            return 0;
		if(nucl(0,paramrelic, ratioH)>0) 
            return 0;
		for(je=0;je<=NNUC;je++) 
            sigma_ratioH[je]=fabs(sigma_ratioH[je]-ratioH[je]);
		return 1;
	}
	else if(err==3)
	{	
		if(nucl(0, paramrelic, ratioH)>0) 
            return 0;
		
		double ratioH_ref[NNUC+1];
		int optfail=0;
		
		if(nucl(-10000, paramrelic, ratioH_ref)>0) 
            optfail=1;	
		for(ie=0;ie<=NNUC;ie++) 
            optfail+=isnan(ratioH_ref[ie]);
		
		double ratioH_tmp[NNUC+1];
		
		for(ie = first; ie <= last; ie++)
		{
			if(optfail==0)
			{
				if(nucl(-ie, paramrelic, ratioH_tmp)>0)
				{
					optfail=1;
					break;break;
				}
			}
						
			for(je=0;je<=NNUC;je++) 
                optfail += isnan(ratioH_tmp[je]);
			
            for(je=0;je<=NNUC;je++) 
                sigma_ratioH[je] += pow(ratioH_tmp[je]-ratioH_ref[je],2.);
		}

		for(ie=0;ie<=NNUC;ie++) 
            sigma_ratioH[ie] = sqrt(sigma_ratioH[ie]);		
		
        for(ie=0;ie<=NNUC;ie++) 
            if(sigma_ratioH[ie]/ratioH[ie]<1.e-10) 
                optfail+=1;
		
		if(optfail>0)
		{
			printf("Sorry, more precise calculation required, please wait...\n");

			for(ie=0;ie<=NNUC;ie++) 
                ratioH_ref[ie]=ratioH[ie];
			for(ie=0;ie<=NNUC;ie++) 
                sigma_ratioH[ie]=0.;
		    for(ie = first; ie <= last; ie++)
			{
				if(nucl_failsafe(-ie, paramrelic, ratioH_tmp)>0) 
                    return 0;
						
				for(je=0;je<=NNUC;je++) 
                    sigma_ratioH[je] += pow(ratioH_tmp[je]-ratioH_ref[je],2.);
			}
			for(ie=0;ie<=NNUC;ie++) 
                sigma_ratioH[ie] = sqrt(sigma_ratioH[ie]);
		}


		for(ie=0;ie<=NNUC;ie++) 
            sigma_ratioH[ie] *= ratioH[ie]/ratioH_ref[ie];

		return 1;
	}
	return 0;
}

/*----------------------------------------------------*/

int bbn_excluded(int err, struct relicparam paramrelic)
/* "container" function which computes the abundances of the elements using the parameters of paramrelic and compares them to observational limits. The err parameter is a switch to choose if the central (err=0), high (err=1) or low (err=2) values of the nuclear rates is used. Returns 1 if the considered BBN scenario is allowed, 0 otherwise. */
{	 
	double H2_H,Yp,Li7_H,Be7_H,He3_H,Li6_H;
	double ratioH[NNUC+1];
		
	if(nucl(err,paramrelic,ratioH)==0)
	{
		H2_H=ratioH[3];
		Yp=ratioH[6];
		Li7_H=ratioH[8];
		Be7_H=ratioH[9];
		He3_H=ratioH[5];
		Li6_H=ratioH[7];
	
#ifdef DEBUG
		printf("Yp\t\t H2/H\t\t He3/H2\t\t Li7/H\t\t Li6/Li7\t Be7/H\n");
		printf("%.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H/H2_H,Li7_H,Li6_H/Li7_H,Be7_H);	
#endif		
		if(isnan(Yp)||isnan(H2_H)||isnan(He3_H/H2_H)||isnan(Li7_H)||isnan(Li6_H/Li7_H)||isnan(Be7_H)) return -1;
		if((Yp<0.258)&&((H2_H>1.2e-5)&&(H2_H<5.3e-5))&&(He3_H/H2_H<1.52)&&(Li7_H>0.85e-10)&&(Li6_H/Li7_H<0.66)) return 0; /* Conservative intervals from hep-ph/0604251 */ 
		else return 1;
	}
	else return -1;
}
