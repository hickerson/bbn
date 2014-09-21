#ifndef __NUCLEARPARAMETERS_HH__
#define __NUCLEARPARAMETERS_HH__

#include <string>
#include <iostream>
#include <ostream>


//static const double NOT_USED = -5555555.5555555;

//using namespace std;

/**
 * NuclearParameters
 *
 * Author: Kevin Peter Hickerson
 * Created: Sep 12, 2014
 */
template <class real, int n> 
class NuclearParameters
{
	//public: static const double Ydefault = 0;
	public: static const unsigned nnuc = n;
			static const unsigned nvar = n + 3;

	public: real T9;	            // temperature in GK = 10^9 Kelvin
			real hv;                    
			real phie;
			//real _Y[nnuc+1];
			//real const _v[nvar+1];
			//real* const V; 
			//real* const Y; 
    private:
			real _V[nvar-nnuc+1]; 
			real _Y[nnuc+1];        // abundances in mass fraction

	public: NuclearParameters(): 
			T9(0),
			hv(0),
			phie(0)
			//Y(&(v[4]))
	{
		//verbose("Default constructor for NuclearParameters");
		for (int i = 1; i <= 3; i++)
			;//_V[i] = 0;
		for (int i = 1; i <= nnuc; i++)
			;//_Y[i] = 0;							/// Default abundance.
		printf("Default: Y address is: %p\n", _Y);
	}

	public: NuclearParameters(const NuclearParameters & copy): 
			T9(copy.T9),
			hv(copy.hv),
			phie(copy.phie)
			//Y(&(v[4]))
	{
		//verbose("Copy constructor for NuclearParameters");
		for (int i = 1; i <= 3; i++)
			_V[i] = copy._V[i];
		for (int i = 1; i <= nnuc; i++)
			_Y[i] = copy._Y[i];					/// Default abundance.
		printf("Copy: Y address is: %p\n", _Y);
	}

	public:~NuclearParameters()				/// Does nothing.
	{
		//verbose("Empty destructor for NuclearParameters");
		printf("Empty destructor for NuclearParameters\n");
	}
	
	public: NuclearParameters & operator = (const NuclearParameters & other)
    {
		//verbose("Copy operator for NuclearParameters");
        if (this != &other) // protect against invalid self-assignment
        {
			T9 = other.T9;
			hv = other.hv;
			phie = other.phie;
			for (int i = 1; i <= 3; i++)
				_V[i] = other._V[i];
			for (int i = 1; i <= nnuc; i++)
				_Y[i] = other._Y[i];				/// Default abundance.
			printf("=: Y address is: %p\n", _Y);
        }
        return *this;
    }

	public: double & Y(const unsigned index)
    {
        if (index <= 0 or index > nnuc)
        {
            std::cout << "Error in NuclearParameters abundances.\n";
            std::cout << "Index out of bounds. (index = " << index << ")\n";
            exit(1);
        }
        return _Y[index];
    }

	public: double & V(const unsigned index)
    {
        if (index > 0 and index <= 3)
            return _V[index];

        std::cout << "Error in NuclearParameters temporaries.\n";
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
std::ostream& operator<<(std::ostream& os, const NuclearParameters<real,n> & ep)
{
    ep.output(os);
    return os;
}


#endif
