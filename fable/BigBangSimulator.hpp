#ifndef __BIGBANGSIMULATOR_HPP__
#define __BIGBANGSIMULATOR_HPP__

#include <fem.hpp> // Fortran EMulation library of fable module
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>

#include "EvolutionParameters.hpp"
//#include "NuclearParameters.hpp"
#include "CosmologicalModel.hpp"
#include "Reaction.hpp"

namespace bbn {

	using namespace fem::major_types;
	//using namespace boost::numeric::ublas;

	static const double NOT_USED = -5555555.5555555;

	struct common_compr0
	{
		double cy0;
		double ct0;
		double T9i0;
		double T9f0;
		double ytmin0;
		int inc0;

		common_compr0() :
			cy0(0),
			ct0(0),
			T9i0(0),
			T9f0(0),
			ytmin0(0),
			inc0(0)
		{}
	};

	struct common_compr
	{
		double cy;
		double ct;
		double T9i;
		double T9f;
		double ytmin;
		int inc;

		common_compr() :
			cy(0),
			ct(0),
			T9i(0),
			T9f(0),
			ytmin(0),
			inc(0)
		{}
	};

	struct common_varpr0
	{
		double dt0;
		double eta0;

		common_varpr0() :
			dt0(0),
			eta0(0)
		{}
	};

	struct common_varpr
	{
		double dt1;
		double eta1;

		common_varpr() :
			dt1(0),
			eta1(0)
		{}
	};

/*
    class ModelParameters
    {
      public:
		double g;
		double tau;
		double xnu;

      private:
		double _c[3+1];
		double _cosmo;
		double _xi[3+1];
    }

    
	struct common_modpr0
	{
		//arr<double> c0;
		double c0[3+1];
		double cosmo0;
		//arr<double> xi0;
		double xi0[3+1];

		common_modpr0() :
			//c0(dimension(3), fem::fill0),
			cosmo0(0)//,
			//xi0(dimension(3), fem::fill0)
		{
			std::cout << "common_modpr0 constructor.\n";
			for (int i = 0; i < 3+1; i++) {
				c0[i] = 0;
				xi0[i] = 0;
			}
		}
	};

	struct common_modpr
	{
		//double g;
		//double tau;
		//double xnu;

		//arr<double> c;
		double c[3+1];
		double cosmo;
		//arr<double> xi;
		double xi[3+1];

		common_modpr() :
			g(0),
			tau(0),
			xnu(0),
			//c(dimension(3), fem::fill0),
			cosmo(0)//,
			//xi(dimension(3), fem::fill0)
		{}

	    // TODO public: double & Y(const unsigned index)
	};
    */

	/*
	   struct common_recpr0
	   {
	   static const int nrec = 88;

	   arr<double, 2> reacpr;

	   common_recpr0() :
	   reacpr(dimension(nrec, 8), fem::fill0)
	   {}
	   };

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_recpr0::nrec;
#endif
	 */

/*
	struct common_recpr
	{
		static const int nrec = 88;

		arr<int> iform;
		arr<int> ii;
		arr<int> jj;
		arr<int> kk;
		arr<int> ll;
		arr<double> rev;
		arr<double> q9;

		common_recpr() :
			iform(dimension(nrec), fem::fill0),
			ii(dimension(nrec), fem::fill0),
			jj(dimension(nrec), fem::fill0),
			kk(dimension(nrec), fem::fill0),
			ll(dimension(nrec), fem::fill0),
			rev(dimension(nrec), fem::fill0),
			q9(dimension(nrec), fem::fill0)
		{}
	};
    */

	/*
#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_recpr::nrec;
#endif
	 */

/*
	struct common_evolp1
	{
		static const int nnuc = 26;
		static const int nvar = 26+3;

		double T9;
		double hv;
		double phie;
		//arr<double> y;
		//double y[nnuc+1];
		
		// shared memory array
		//double U.v[nvar+1];
		//double* U.y;

		common_evolp1() :
			U.T9(0),
			U.hv(0),
			U.phie(0)//,
			//y(dimension(nnuc), fem::fill0)
		{
			std::cout << "common_evolp1 constructor\n";
			//y[0] = NOT_USED;
			//for (int i = 1; i < nnuc+1; i++)
			v[0] = NOT_USED;
			for (int i = 1; i < nvar+1; i++)
				v[i] = 0;
			y = &(v[4]);
		}
	};

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
	const int common_evolp1::nnuc;
#endif

	struct common_evolp2
	{
		static const int nnuc = 26;

		double dT9;
		double dhv;
		double dphie;
		//arr<double> dydt;
		double dydt[nnuc+1];

		common_evolp2() :
			dT9(0),
			dhv(0),
			dphie(0)//,
			//dydt(dimension(nnuc), fem::fill0)
		{
			std::cout << "common_evolp2 constructor\n";
			dydt[0] = NOT_USED;
			for (int i = 1; i < nnuc+1; i++)
				dydt[i] = 0;
		}
	};

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
	const int common_evolp2::nnuc;
#endif

	struct common_evolp3
	{
		static const int nnuc = 26;

		double T90;
		double hv0;
		double phie0;
		//arr<double> y0;
		double y0[nnuc+1];

		common_evolp3() :
			T90(0),
			hv0(0),
			phie0(0)//,
			//y0(dimension(nnuc), fem::fill0)
		{
			std::cout << "common_evolp3 constructor\n";
			y0[0] = NOT_USED;
			for (int i = 1; i < nnuc+1; i++)
				y0[i] = 0;
		}
	};


#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
	const int common_evolp3::nnuc;
#endif
*/


	struct common_ttime
	{
		double t;
		double dt;
		double dlT9dt;

		common_ttime() :
			t(0),
			dt(0),
			dlT9dt(0)
		{}
	};

	struct common_thermcb
	{
		arr<double> thm;
		double hubcst;

		common_thermcb() :
			thm(dimension(14), fem::fill0),
			hubcst(0)
		{}
	};

	struct common_endens
	{
		double rhone0;
		double rhob0;
		double rhob;
		double rnb;

		common_endens() :
			rhone0(0),
			rhob0(0),
			rhob(0),
			rnb(0)
		{}
	};

	/*
	   struct common_nucdat
	   {
	   static const int nnuc = 26;

	   arr<double> am;
	   arr<double> zm;
	   arr<double> dm;

	   common_nucdat() :
	   am(dimension(nnuc), fem::fill0),
	   zm(dimension(nnuc), fem::fill0),
	   dm(dimension(nnuc), fem::fill0)
	   {}
	   };

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_nucdat::nnuc;
#endif
	 */

	/*
	   struct common_besselcb
	   {
	   double bl1;
	   double bl2;
	   double bl3;
	   double bl4;
	   double bl5;
	   double bm1;
	   double bm2;
	   double bm3;
	   double bm4;
	   double bm5;
	   double bn1;
	   double bn2;
	   double bn3;
	   double bn4;
	   double bn5;

	   common_besselcb() :
	   bl1(0),
	   bl2(0),
	   bl3(0),
	   bl4(0),
	   bl5(0),
	   bm1(0),
	   bm2(0),
	   bm3(0),
	   bm4(0),
	   bm5(0),
	   bn1(0),
	   bn2(0),
	   bn3(0),
	   bn4(0),
	   bn5(0)
	   {}
	   };

	   struct common_kays
	   {
	   double bk0;
	   double bk1;
	   double bk2;
	   double bk3;
	   double bk4;

	   common_kays() :
	   bk0(0),
	   bk1(0),
	   bk2(0),
	   bk3(0),
	   bk4(0)
	   {}
	   };
	 */

	struct common_flags
	{
		int ltime;
		int is;
		int ip;
		int it;
		int mbad;

		common_flags() :
			ltime(0),
			is(0),
			ip(0),
			it(0),
			mbad(0)
		{}
	};

	struct common_checkcb
	{
		int itime;

		common_checkcb() :
			itime(0)
		{}
	};

	// TODO make a class that is arrayed.
	struct common_outdat
	{
		static const int itmax = 40;
		static const int nnuc = 26;

		arr<double, 2> xout;
		arr<double, 2> thmout;
		arr<double> T9out;
		arr<double> tout;
		arr<double> dtout;
		arr<double> etaout;
		arr<double> hubout;

		common_outdat() :
			xout(dimension(itmax, nnuc), fem::fill0),
			thmout(dimension(itmax, 6), fem::fill0),
			T9out(dimension(itmax), fem::fill0),
			tout(dimension(itmax), fem::fill0),
			dtout(dimension(itmax), fem::fill0),
			etaout(dimension(itmax), fem::fill0),
			hubout(dimension(itmax), fem::fill0)
		{}
	};

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
	const int common_outdat::itmax;
	const int common_outdat::nnuc;
#endif

	struct common_nupar
	{
		double T9mev;
		double tnmev;
		double Tnu;
		double cnorm;
		double rhonu;
		int nu;

		common_nupar() :
			T9mev(0),
			tnmev(0),
			Tnu(0),
			cnorm(0),
			rhonu(0),
			nu(0)
		{}
	};

	struct common_runopt
	{
		int irun;
		int isize;
		int jsize;

		common_runopt() :
			irun(0),
			isize(0),
			jsize(0)
		{}
	};

	struct common_outopt
	{
		int nout;
		bool outfile;

		common_outopt() :
			nout(0),
			outfile(false)
		{}
	};

	struct common_tcheck
	{
		int itime;

		common_tcheck() :
			itime(0)
		{}
	};

	struct common :
		fem::common,
		common_compr0,
		common_compr,
		common_varpr0,
		common_varpr,
		//common_modpr0,
		//common_modpr,
		//common_recpr0,
		//common_recpr,
		//common_evolp1,
		//common_evolp2,
		//common_evolp3,
		common_ttime,
		common_thermcb,
		common_endens,
		//common_nucdat,
		//common_besselcb,
		//common_kays,
		common_flags,
		//common_checkcb,
		common_outdat,
		common_nupar,
		common_runopt,
		common_outopt,
		common_tcheck
	{
		static const int nrec = 88;
		static const int nnuc = 26;
		int Nreactions = 88;
		int Nnuclide = 26;
        
		//static const int si[] = {1, 1, 1, 1, 1, 2, 3, 2, 1, 1, 2};
		//static const int sj[] = {0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0};
		//static const int sk[] = {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2};
		//static const int sl[] = {1, 1, 1, 2, 2, 1, 1, 1, 2, 3, 1};

		// TODO put in Nuclide
        //vector< Nuclide<double> > nuclides(nnuc+1);
		static const double am[nnuc+1];
		static const double zm[nnuc+1];
		static const double dm[nnuc+1];

		// TODO put in a subclass?
		double a[nnuc+1][nnuc+1];
		double b[nnuc+1];
		//double yx[nnuc+1];
		/*
        boost::numeric::ublas::matrix<double> a(nnuc,nnuc);
        boost::numeric::ublas::vector<double> b(nnu);
        boost::numeric::ublas::vector<double> y(nnuc);
        *//*
        using namespace boost::numeric::ublas {
            matrix<double> a(nnuc,nnuc);
            vector<double> b(nnu);
            vector<double> y(nnuc);
        }
		*/

		//arr<double> f;
		//arr<double> r;
		//arr<double,double> a;
		//arr<double> b;
		//arr<double> yx;
		double f[nrec+1];
		double r[nrec+1];
        Reaction<double,2> reactions[nrec+1];
        //vector< Reaction<double,2> > reactions(nrec+1);
        //vector<double> f(nrec+1);
        //vector<double> r(nrec+1);
		static const double reacpr[nrec][8];


		void qvary(int, double);        // TODO remove?

	private:
		EvolutionParameters<double, nnuc> U, U0, dU, dUdt, dUdt0; 
        CosmologicalModel<double> M, M0;//, dM;

	public: 
		common();
		//	int argc,
		//	char const* argv[]);

		void program_new123();
		void help();
		void setcom();
		void setmod();
		void check();
		void driver();
		void run();

		void therm();
		void nudens();
		void start();

		//double xintd(
		//		double const&,
		//		double const&,
		//		double(* func)(const double&),
		//		int const&);
		template <int n>
		double xintd(
				double const&,
				double const&,
				int const&);

		void eqslin(
				int const&,
				int&);

		void sol( int const&);
		void derivs( int const& );
		void accum();
		void output();

		template <int n>
		double integrand( double const& );

		void rate0();
		void rate1(double const& );
		void rate2();
		void rate3();
		void rate4();

	private:
		double getBesselL(double);
		double getBesselM(double);
		double getBesselN(double);


        bool isNuclideIndex(int index) {
            if (index > 0 and index <= nnuc)
                return true;
            exit(1);
            return false;
        }
		double& y(unsigned index) { 
            isNuclideIndex(index);
            return U.Y(index); 
        }
		double& y0(unsigned index) { 
            isNuclideIndex(index);
            return U0.Y(index); 
        }
		double& dy(unsigned index) { 
            isNuclideIndex(index);
            return dU.Y(index); 
        }
		double& dydt(unsigned index) { 
            isNuclideIndex(index);
            return dUdt.Y(index); 
        }
		double& dydt0(unsigned index) { 
            isNuclideIndex(index);
            return dUdt0.Y(index); 
        }

		double y(unsigned index) const { 
            isNuclideIndex(index);
            return U.Y(index); 
        }
		double y0(unsigned index) const { 
            isNuclideIndex(index);
            return U0.Y(index); 
        }
		double dy(unsigned index) const { 
            isNuclideIndex(index);
            return dU.Y(index); 
        }
		double dydt(unsigned index) const { 
            isNuclideIndex(index);
            return dUdt.Y(index); 
        }
		double dydt0(unsigned index) const { 
            isNuclideIndex(index);
            return dUdt0.Y(index); 
        }


        bool isVIndex(int index) {
            if (index > 0 and index <= nnuc + 3)
                return true;
            exit(1);
            return false;
        }
		double& v(unsigned index) { 
            isVIndex(index);
            return index < 4? U.V(index) : U.Y(index-3); 
        }
		double& v0(unsigned index) { 
            isVIndex(index);
            return index < 4? U0.V(index) : U0.Y(index-3); 
        }
		double& dv(unsigned index) { 
            isVIndex(index);
            return index < 4? dU.V(index) : dU.Y(index-3); 
        }
		double& dvdt(unsigned index) { 
            isVIndex(index);
            return index < 4? dUdt.V(index) : dUdt.Y(index-3); 
        }
		double& dvdt0(unsigned index) { 
            isVIndex(index);
            return index < 4? dUdt0.V(index) : dUdt0.Y(index-3); 
        }

        public: void output(std::ostream& os) const
        {
            os << "Current state (U):\n"; 
            os << U << "\n";
            os << "Initial state (U0):\n"; 
            os << U0 << "\n";
            os << "Differential (dU):\n"; 
            os << dU << "\n";
            os << "Time derivative (dUdt):\n"; 
            os << dUdt << "\n";
            os << "Initial time derivative (dUdt0):\n"; 
            os << dUdt0 << "\n";

            return os;
        }
	};


    std::ostream& operator<<(std::ostream& os, const common & BBN)
    {
        BBN.output(os);
        return os;
    }

	const double common::am[26+1] = {
		NOT_USED,
		1.f, 1.f, 2.f, 3.f, 3.f, 4.f, 6.f, 7.f, 7.f, 8.f, 8.f, 9.f,
		10.f, 11.f, 11.f, 12.f, 12.f, 12.f, 13.f, 13.f, 14.f, 14.f,
		14.f, 15.f, 15.f, 16.f
	};

	const double common::zm[26+1] = {
		NOT_USED,
		0.f, 1.f, 1.f, 1.f, 2.f, 2.f, 3.f, 3.f, 4.f, 3.f, 5.f, 4.f,
		5.f, 5.f, 6.f, 5.f, 6.f, 7.f, 6.f, 7.f, 6.f, 7.f, 8.f, 7.f,
		8.f, 8.f
	};

	const double common::dm[26+1] = {
		NOT_USED,
		.008665f, .007825f, .014102f, .016050f, .016030f, .002603f,
		.015125f, .016004f, .016929f, .022487f, .024609f, .012186f,
		.012939f, .009305f, .011432f, .014354f, .000000f, .018641f,
		.003354f, .005738f, .003242f, .003074f, .008597f, .000108f,
		.003070f, -.005085f
	};

/*
	const double common::s[4][11] = {
        {1, 1, 1, 1, 1, 2, 3, 2, 1, 1, 2},
	    {0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0},
	    {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2},
	    {1, 1, 1, 2, 2, 1, 1, 1, 2, 3, 1}
    };
*/
	const double common::reacpr[88][8] = {
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

} // namespace bbn


#endif // GUARD
