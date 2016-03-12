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
#include "Reaction.h"

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

public:
	double heff(double Temp);
	double sgStar(double Temp);
	double geff(double Temp);
	double dark_density(double T);
	double dark_entropy(double T);
	double dark_entropy_derivative(double T);
	double dark_entropy_Sigmad(double T);
	double nonthermal(double T);
	double neutdens(double Tnu);
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




/* bbnrate.c */
void rate_weak(int err, ReactionMap & f);
void rate_pn(int err, CosmologyModel relic,
            ReactionMap & f, ReactionMap & r, double T9, double Tnu);
void rate_all(int err, ReactionMap & f, double T9);




/* bbn.c */
void setup_reactions(ReactionList & reaction);
void setup_nuclides(int A[], int Z[], double dm[]);
int linearize(double T9, ReactionList & reactions, 
              ReactionMap & f, ReactionMap & r, int loop, int inc, int ip, 
	double dt, NuclideMap & y0, NuclideMap & y, NuclideMap & dydt, 
	double H, double rhob);
int nucl(int err, CosmologyModel relic, NuclideMap & ratioH);
int nucl_failsafe(int err, CosmologyModel relic, NuclideMap & ratioH);
int nucl_witherrors(int err, CosmologyModel relic, NuclideMap & ratioH, NuclideMap & sigma_ratioH);
int bbn_excluded(int err, CosmologyModel relic);

#endif  /// __INCLUDE_H__
