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
    public:
    real G;             /// Gravitational constant.
    //real dG;            /// Variation of Gravitational constant.
    real tau;           /// Neutron lifetime.
    real cosmo;         /// Cosmological constant.
    real eta;           /// Baryon-to-photon ratio.
    real Nnu;           /// Number of neutrinos.
    real T9i;           /// Initial temperature.
    real T9f;           /// Final temperature.

    private:
    real xi[3+1];       /// Neutrino degeneracy parameter.
    real c[3+1];        /// Deprecated coefficients.

    /*
	public: CosmologicalModel():
			G(0),
			dG(0),
			tau(0),
			cosmo(0),
			eta(0),
			Nnu(0),
			Ti(0),
			Tf(0)
	{
		//verbose("Default constructor for CosmologicalModel");
		for (int i = 0; i <= 3; i++)
        {
			c[i] = 0;
			xi[i] = 0;
        }
	}
    */

	public: CosmologicalModel(const CosmologicalModel & copy): 
			G(copy.G),
			dG(copy.dG),
			tau(copy.tau),
			cosmo(copy.cosmo),
			eta(copy.eta),
			Nnu(copy.Nnu),
			Ti(copy.Ti),
			Tf(copy.Tf)
	{
		//verbose("Copy constructor for CosmologicalModel");
		for (int i = 0; i <= 3; i++)
        {
			c[i] = copy.c[i];
			xi[i] = copy.xi[i];
        }
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
