#ifndef __REACTION_HH__
#define __REACTION_HH__

#include <string>
#include <iostream>
#include <ostream>


/**
 * Reaction
 *
 * Author: Kevin Peter Hickerson
 * Created: Sep 13, 2014
 * 
 * Template parameters:
 *      real: the mathematical type to represent a number.
 *      n:    maximum number of nuclide types in or out of the reaction.
 */
template <class real, int n> 
class Reaction
{
	public: 
    int type;       /// Reaction type.
    int in[n];      /// Incoming nuclide types.
    int out[n];     /// Outgoing nuclide types.
    real reverse;   /// Reverse reaction coefficient.
    real energy;    /// Endpoint energy of this reaction.


    private:

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
        type = reaction_data[i-1][1];       /// Reaction type.
        reverse = reaction_data[i-1][6];    /// Reverse reaction coefficient.
        energy = reaction_data[i-1][7];     /// Energy released.
        in[i] = reaction_data[i-1][2]; 	    /// Incoming nuclide type.
        in[i] = reaction_data[i-1][3]; 	    /// Incoming nuclide type.
        out[i] = reaction_data[i-1][4]; 	/// Outgoing nuclide type.
        out[i] = reaction_data[i-1][5]; 	/// Outgoing nuclide type.
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
			type = other.type;
			reverse = other.reverse;
			energy = other.energy;
            for (int i = 0; i <= n; i++) {
                in[i] = other.in[i];
                out[i] = other.out[i];
            }
        }
        return *this;
    }

	public: double & Y(const unsigned index)
    {
        if (index <= 0 or index > nnuc)
        {
            std::cout << "Error in Reaction abundances.\n";
            std::cout << "Index out of bounds. (index = " << index << ")\n";
            exit(1);
        }
        return _Y[index];
    }

	public: double & V(const unsigned index)
    {
        if (index > 0 and index <= 3)
            return _V[index];

        std::cout << "Error in Reaction temporaries.\n";
        std::cout << "Index out of bounds. (index = " << index << ")\n";
        exit(1);
    }

    public: const void output(std::ostream& os) const
    {
        os << "T9: " << T9 << "\n";
        os << "hv: " << hv << "\n";
        os << "phie: " << phie << "\n";
        os << "V address is: " << _Y << "\n";
        for (int i = 1; i <= 3; i++)
            os << "\t" << _V[i];
        os << "\n";
        os << "Y address is: " << _Y << "\n";
        for (int i = 1; i <= nnuc; i++)
            os << "\t" << _Y[i];
        os << "\n";
    }
};

template <class real, int n> 
std::ostream& operator<<(std::ostream& os, const Reaction<real,n> & ep)
{
    ep.output(os);
    return os;
}


#endif
