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
#include <string.h>
#include <strings.h>
#include "Reaction.h"
#include <vector>
//#include <complex>
//#include <complex.h>

/*--------------------------------------------------------------------*/

/*#define DEBUG*/

/*--------------------------------------------------------------------*/

#define pi      3.1415926535897932385
#define zeta3   1.2020569031595942855
#define sqrt12  3.4641016151377545870
#define sqrtpi  1.7724538509055160273

#define hbar    6.58211889e-25    /// Planck constant in GeV.s
#define Gn      6.67428e-8        /// Gravitational constant in m^3/g.s^2
#define me      0.510998928       /// electron mass in MeV/c^2
#define kB      0.0861733238      /// Boltzmann constant in MeV/GK
#define Mplanck 1.2209e19         /// Planck mass in GeV/c^2


/*--------------------------------------------------------------------*/

/**
 * structure containing the cosmological model parameters 
 */
class CosmologyModel {
public:
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

	NuclideList nuclides;
	ReactionList reactions;

public:
    CosmologyModel() { 
        Init_cosmomodel(); 
    }
    CosmologyModel(double eta, double nbnu, double life_neutron, 
                   double xinu1, double xinu2, double xinu3) {
        Init_cosmomodel(); 
	    Init_cosmomodel_param(eta, nbnu, life_neutron, xinu1, xinu2, xinu3);
    }
    CosmologyModel(double eta, double nbnu, double life_neutron) {
        Init_cosmomodel(); 
	    Init_cosmomodel_param(eta, nbnu, life_neutron, 0, 0, 0);
    }

public:
	void Init_modeleff(int model_eff);
	void Init_cosmomodel();
	void Init_cosmomodel_param(double eta, double nbnu, double life_neutron, double xinu1, double xinu2, double xinu3);
	void Init_dark_density(double dd0, double ndd, double T_end);
	void Init_dark_entropy(double sd0, double nsd, double T_end);
	void Init_dark_entropySigmaD(double Sigmad0, double nSigmad, double T_end);
	void Init_nonthermal(double nt0, double nnt, double T_end);
	void Init_fierz(double eta, double nbnu, double life_neutron, double fierz);
	void Init_nuclides();
	void Init_reactions();

public:
	double heff(double Temp) const;
	double sgStar(double Temp) const;
	double geff(double Temp) const;
	double dark_density(double T) const;
	double dark_entropy(double T) const;
	double dark_entropy_derivative(double T) const;
	double dark_entropy_Sigmad(double T) const;
	double nonthermal(double T) const;
	double neutrino_density(double Tnu) const;
};


/**
 * Strongly typed Flag for type of error computation to perform and display
 */
enum ErrorType {
    mean,
    high,
    low,
    gaussian
};


/*--------------------------------------------------------------------*/

/* general.cpp */
double max(double x, double y);
double min(double x, double y);
#if __cplusplus!=201103L
double complex polylog(int n, int m, double x);
double complex CLi2(double complex x);
double complex CLi3(double complex x);
#endif 
double Li2(double x);
double Li3(double x);
double Li4(double x);
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
float erfinv(float x);	// TODO change to double


using namespace std;

struct distribution {
	double mean;
	double min;
	double max;
	double error;
	int samples;
	vector<double> values;

	distribution(double value)
	  : mean(value),
		min(value),
		max(value),
		error(0),
		samples(1) {
		values.push_back(mean);
	}

	/** 
	 * gaussian distribution
	 * mean		the central value of the gussian
	 * error	the standard error of the gaussian
	 * sigma	the number of standard deviations from the mean to include
	 */
	distribution(double mean, double error, double sigma, int samples)
	  :	mean(mean),
		min(mean-sigma*error),
		max(mean+sigma*error),
		error(error),
		samples(samples)
	{
		double x;
		double value=mean;
		if (samples < 1) 
			samples = 1;
		double xmin = erf(min);
		double xmax = erf(max);
		double bin = (xmax-xmin)/samples;
		for (x=xmin; x<=xmax; x+=bin) 
		{
			value = mean + error*erfinv(x);
			values.push_back(value);
		}
		//max = value;
		//min = mean-2*value;
	}

	/** 
	 * uniform distribution 
	 */
	distribution(double min, double max, int samples) :
		mean((max-min)/2),
		min(min),
		max(max),
		error((max-min)/sqrt12),
		samples(samples<2?2:samples)
	{
		double bin = (max-min)/samples;
		if (bin <= 0)
		for (double x=min; x<=max; x+=bin) 
			values.push_back(x);
	}

	distribution & operator=(double value) {
		mean = value;
		min = value;
		max = value;
		error = 0;
		samples = 1;
		values.clear();
		values.push_back(mean);
		return *this;
	}

	distribution & operator=(distribution d) {
		mean = d.mean;
		min = d.min;
		max = d.max;
		error = d.error;
		samples = d.samples;
		values.clear();
		values = d.values;
		return *this;
	}
	/*
	public:init(const char* str) {
		scanf("%f", str)
	}
	*/
};


/* bbnrate.cpp */
void rate_weak(int err, ReactionMap & f);
void rate_pn(int err, const CosmologyModel & relic,
            ReactionMap & f, ReactionMap & r, double T9, double Tnu);
void rate_all(int err, ReactionMap & f, double T9);




/* bbn.cpp */
void setup_reactions(ReactionList & reaction);
void setup_nuclides(NuclideList & nuclides);

int linearize(double T9, const ReactionList & reactions, 
              ReactionMap & f, ReactionMap & r, int loop, int inc, int ip, 
	double dt, NuclideMap & y0, NuclideMap & y, NuclideMap & dydt, 
	double H, double rhob);
int nucl(int err, const CosmologyModel & relic, NuclideMap & ratioH);
int nucl_failsafe(int err, const CosmologyModel & relic, NuclideMap & ratioH);
int nucl_errors(int err, const CosmologyModel & relic, NuclideMap & ratioH, NuclideMap & sigma_ratioH);
//int bbn_excluded(int err, CosmologyModel relic);

#endif  /// __INCLUDE_H__
