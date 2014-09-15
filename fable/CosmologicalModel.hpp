#ifndef __COSMOLOGICALMODEL_HH__
#define __COSMOLOGICALMODEL_HH__

#include <string>
#include <iostream>
#include <ostream>


/**
 * CosmologicalModel
 *
 * Author: Kevin Peter Hickerson
 * Created: Sun Apr 20 17:11:14 PDT 2014
 */
template <class real> 
class CosmologicalModel
{ 
    real G;             /// Gravitational constant.
    real tau;           /// Neutron lifetime.
    real Nnu            /// Number of neutrinos.
    real cosmo          /// Cosmological constant.
    real c[3+1];
    real xi[3+1];

	public: CosmologicalModel()
	{
		//verbose("Default constructor for CosmologicalModel");
		printf("Default: Y address is: %p\n", _Y);
	}

	public: CosmologicalModel(const CosmologicalModel & copy): 
			G(copy.G),
			tau(copy.tau),
			Nnu(copy.Nnu)
			cosmo(copy.cosmo)
	{
		//verbose("Copy constructor for CosmologicalModel");
		for (int i = 1; i <= 3; i++)
        {
			c[i] = copy.c[i];
			xi[i] = copy.xi[i];
        }
		printf("Copy: Y address is: %p\n", _Y);
	}

	public:~CosmologicalModel()				/// Does nothing.
	{
		//verbose("Empty destructor for CosmologicalModel");
		printf("Empty destructor for CosmologicalModel\n");
	}
	
	public: CosmologicalModel & operator = (const CosmologicalModel & other)
    {
		//verbose("Copy operator for CosmologicalModel");
        if (this != &other) // protect against invalid self-assignment
        {
			G = other.G;
			tau = other.tau;
			Nnu = other.Nnu;
			cosmo = other.cosmo;
            for (int i = 1; i <= 3; i++)
            {
                c[i] = other.c[i];
                xi[i] = other.xi[i];
            }
        }
        return *this;
    }


    /*
	public: double & V(const unsigned index)
    {
        if (index > 0 and index <= 3)
            return _V[index];

        std::cout << "Error in CosmologicalModel temporaries.\n";
        std::cout << "Index out of bounds. (index = " << index << ")\n";
        exit(1);
    }
    */

    public: const void output(std::ostream& os) const
    {
        os << "G: " << G << "\n";
        os << "tau: " << tau << "\n";
        os << "Nnu: " << Nnu << "\n";
        os << "cosmo: " << cosmolgical << "\nc:";
        for (int i = 1; i <= 3; i++)
            os << "\t" << c[i];
        os << "\nxi: ";
        for (int i = 1; i <= nnuc; i++)
            os << "\t" << xi[i];
        os << "\n";
    }
};

template <class real, int n> 
std::ostream& operator<<(std::ostream& os, const CosmologicalModel<real,n> & ep)
{
    ep.output(os);
    return os;
}


#endif
