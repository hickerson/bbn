#ifndef __COSMOLOGICALMODEL_HH__
#define __COSMOLOGICALMODEL_HH__

#include <string>
#include <iostream>
#include <ostream>


/**
 * CosmologicalModel
 *
 * Author: Kevin Peter Hickerson
 * Created: Sep 10, 2014
 */
template <class real> 
class CosmologicalModel
{ 
    public:
    real G;             /// Gravitational constant.
    real dG;            /// Variation of Gravitational constant. // TODO
    real tau;           /// Neutron lifetime.
    real cosmo;         /// Cosmological constant.
    real eta;           /// Baryon-to-photon ratio.
    real Nnu;           /// Number of neutrinos.
    real T9i;           /// Initial temperature.
    real T9f;           /// Final temperature.

    real xi[3+1];       /// Neutrino degeneracy parameter.
    //private:            /// TODO make dependent on number of neutrinos
    real c[3+1];        /// Deprecated coefficients.

	public: CosmologicalModel():
			G(0),
			dG(0),
			tau(0),
			cosmo(0),
			eta(0),
			Nnu(0),
			T9i(0),
			T9f(0)
	{
		//verbose("Default constructor for CosmologicalModel");
		for (int i = 0; i <= 3; i++)
        {
			c[i] = 0;
			xi[i] = 0;
        }
	}

	public: CosmologicalModel(const CosmologicalModel & copy): 
			G(copy.G),
			dG(copy.dG),
			tau(copy.tau),
			cosmo(copy.cosmo),
			eta(copy.eta),
			Nnu(copy.Nnu),
			T9i(copy.T9i),
			T9f(copy.T9f)
	{
		//verbose("Copy constructor for CosmologicalModel");
		for (int i = 0; i <= 3; i++)
        {
			xi[i] = copy.xi[i];
			c[i] = copy.c[i];
        }
	}

	public:~CosmologicalModel()				/// Does nothing.
	{
		//verbose("Empty destructor for CosmologicalModel");
		//printf("Empty destructor for CosmologicalModel\n");
	}
	
	public: CosmologicalModel & operator = (const CosmologicalModel & other)
    {
		//verbose("Copy operator for CosmologicalModel");
        if (this != &other) // protect against invalid self-assignment
        {
			G = other.G;
			dG = other.dG;
			tau = other.tau;
			Nnu = other.Nnu;
			cosmo = other.cosmo;
            for (int i = 1; i <= 3; i++)
            {
                xi[i] = other.xi[i];
                c[i] = other.c[i];
            }
        }
        return *this;
    }

    public: const void output(std::ostream& os, string sep) const
    {
        os << "G: " << G << sep;
        os << "dG: " << dG << sep;
        os << "tau: " << tau << sep;
        os << "Nnu: " << Nnu << sep;
        os << "cosmo: " << cosmo << sep;
        os << "c: ";
        for (int i = 1; i <= 3; i++)
            os << "\t" << c[i];
        os << "\nxi: ";
        for (int i = 1; i <= 3; i++)
            os << "\t" << xi[i];
        os << "\n";
    }
};

template <class real, int n> 
std::ostream& operator<<(std::ostream& os, const CosmologicalModel<real> & ep)
{
    ep.output(os);
    return os;
}


#endif
