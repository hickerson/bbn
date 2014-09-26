#ifndef __REACTION_HH__
#define __REACTION_HH__

#include <string>
#include <iostream>
#include <ostream>

using namespace std;

/**
 * Reaction
 *
 * Author: Kevin Peter Hickerson
 * Created: Sep 13, 2014
 * 
 * Template parameters:
 *      real: the mathematical type to represent a real number.
 *      n:    maximum number of nuclide types in or out of the reaction.
 */
template <class real, int n> 
class Reaction
{
    const static real reaction_data[88][8];
    const static real reaction_type[4][11];

    public: 
        string name;    /// Reaction name.
        int type;       /// Reaction type.
        int in[n];      /// Incoming nuclide types.
        int out[n];     /// Outgoing nuclide types.
        real reverse;   /// Reverse reaction coefficient.
        real energy;    /// Endpoint energy of this reaction.

	public: 
    Reaction(): 
        type(0),
        reverse(0),
        energy(0)
	{
		//verbose("Default constructor for Reaction");
		for (int i = 0; i <= n; i++) {
            in[i] = 0;
            out[i] = 0;
        }
	}

	public: Reaction(int i) 
	{
        /*
        iform(i) = reacpr[i-1][2-1]; 		/// Reaction type.
        ii(i) = reacpr[i-1][3-1]; 			/// Incoming nuclide type.
        jj(i) = reacpr[i-1][4-1]; 			/// Incoming nuclide type.
        kk(i) = reacpr[i-1][5-1]; 			/// Outgoing nuclide type.
        ll(i) = reacpr[i-1][6-1]; 			/// Outgoing nuclide type.
        rev(i) = reacpr[i-1][7-1]; 			/// Reverse reaction coefficient.
        q9(i) = reacpr[i-1][8-1]; 			/// Energy released.
        */
        int index = reaction_data[i-1][0];  /// Reaction index.
        std::cout << "i: " << i << " r[o]: " << index << "\n";
        assert(i == index);                 /// Reaction index.

        type = reaction_data[i-1][1];       /// Reaction type.
        in[i] = reaction_data[i-1][2]; 	    /// Incoming nuclide type.
        in[i] = reaction_data[i-1][3]; 	    /// Incoming nuclide type.
        out[i] = reaction_data[i-1][4]; 	/// Outgoing nuclide type.
        out[i] = reaction_data[i-1][5]; 	/// Outgoing nuclide type.
        reverse = reaction_data[i-1][6];    /// Reverse reaction coefficient.
        energy = reaction_data[i-1][7];     /// Energy released.
	}

	public: Reaction(const Reaction & copy): 
			type(copy.type),
			reverse(copy.reverse),
			energy(copy.energy)
	{
		//verbose("Copy constructor for Reaction");
		for (int i = 0; i <= n; i++) {
            in[i] = copy.in[i];
            out[i] = copy.out[i];
        }
	}

	public:~Reaction()				/// Does nothing.
	{
		//verbose("Empty destructor for Reaction");
		printf("Empty destructor for Reaction\n");
	}
	
	public: Reaction & operator = (const Reaction & other)
    {
		//verbose("Copy operator for Reaction");
        if (this != &other) // protect against invalid self-assignment
        {
			name = other.name;
			type = other.type;
			reverse = other.reverse;
			energy = other.energy;
            for (int i = 0; i < n; i++) {
                in[i] = other.in[i];
                out[i] = other.out[i];
            }
        }
        return *this;
    }

    int getNuclideIndicies(int& i, int& j, int& k, int& l)
    {
		i = in[0]; 		    /// ID of incoming nuclide i.
		j = in[1]; 		    /// ID of incoming nuclide j.
		k = out[0]; 		/// ID of outgoing nuclide k.
		l = out[1]; 		/// ID of outgoing nuclide l.
        return type;
    }

    int getNuclideCounts(int& i, int& j, int& k, int& l)
    {
        i = reaction_type[0][type-1]; 	/// # of incoming nuclide i.
        j = reaction_type[1][type-1]; 	/// # of incoming nuclide j.
        k = reaction_type[2][type-1]; 	/// # of outgoing nuclide k.
        l = reaction_type[3][type-1]; 	/// # of outgoing nuclide l.
        return i + j + k + l;
    }

    public: const void output(std::ostream& os) const
    {
        os << "name: " << name << "\n";
        os << "type: " << type << "\n";
        os << "in: ";
        for (int i = 0; i < n; i++)
            os << "\t" << in[i];
        os << "\nout: ";
        for (int i = 0; i < n; i++)
            os << "\t" << out[i];
        os << "reverse: " << reverse << "\n";
        os << "energy: " << energy << "\n";
        os << "\n";
    }
};

template <class real, int n> 
std::ostream& operator<<(std::ostream& os, const Reaction<real,n> & ep)
{
    ep.output(os);
    return os;
}


template <class real, int n> 
const real Reaction<real,n>::reaction_type[4][11] = {
    {1, 1, 1, 1, 1, 2, 3, 2, 1, 1, 2},
    {0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0},
    {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2},
    {1, 1, 1, 2, 2, 1, 1, 1, 2, 3, 1}
};


template <class real, int n> 
const real Reaction<real,n>::reaction_data[88][8] = {
    {1,	    1,	1,	0,	0,	2,	0,  	0},
    {2,	    1,	4,	0,	0,	5,	0,  	0},
    {3,	    4,	10,	0,	0,	6,	0,  	0},
    {4,	    1,	16,	0,	0,	17,	0,  	0},
    {5,	    1,	21,	0,	0,	22,	0,  	0},
    {6,	    4,	11,	0,	0,	6,	0,  	0},
    {7,	    1,	15,	0,	0,	14,	0,  	0},
    {8,	    1,	18,	0,	0,	17,	0,  	0},
    {9,	    1,	20,	0,	0,	19,	0,  	0},
    {10,	1,	23,	0,	0,	22,	0,  	0},
    {11,	1,	25,	0,	0,	24,	0,	    0},
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

#endif
