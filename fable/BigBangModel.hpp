#ifndef BBN_CMN_HPP
#define BBN_CMN_HPP

#include <fem.hpp> // Fortran EMulation library of fable module

namespace bbn {

using namespace fem::major_types;

struct common_compr0
{
  float cy0;
  float ct0;
  float t9i0;
  float t9f0;
  float ytmin0;
  int inc0;

  common_compr0() :
    cy0(0),
    ct0(0),
    t9i0(0),
    t9f0(0),
    ytmin0(0),
    inc0(0)
  {}
};

struct common_compr
{
  float cy;
  float ct;
  float t9i;
  float t9f;
  float ytmin;
  int inc;

  common_compr() :
    cy(0),
    ct(0),
    t9i(0),
    t9f(0),
    ytmin(0),
    inc(0)
  {}
};

struct common_varpr0
{
  float dt0;
  float eta0;

  common_varpr0() :
    dt0(0),
    eta0(0)
  {}
};

struct common_varpr
{
  float dt1;
  float eta1;

  common_varpr() :
    dt1(0),
    eta1(0)
  {}
};

struct common_modpr0
{
  arr<float> c0;
  float cosmo0;
  arr<float> xi0;

  common_modpr0() :
    c0(dimension(3), fem::fill0),
    cosmo0(0),
    xi0(dimension(3), fem::fill0)
  {}
};

struct common_modpr
{
  float g;
  float tau;
  float xnu;
  arr<float> c;
  float cosmo;
  arr<float> xi;

  common_modpr() :
    g(0),
    tau(0),
    xnu(0),
    c(dimension(3), fem::fill0),
    cosmo(0),
    xi(dimension(3), fem::fill0)
  {}
};

struct common_recpr0
{
  static const int nrec = 88;

  arr<float, 2> reacpr;

  common_recpr0() :
    reacpr(dimension(nrec, 8), fem::fill0)
  {}
};

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_recpr0::nrec;
#endif

struct common_recpr
{
  static const int nrec = 88;

  arr<int> iform;
  arr<int> ii;
  arr<int> jj;
  arr<int> kk;
  arr<int> ll;
  arr<float> rev;
  arr<float> q9;

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

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_recpr::nrec;
#endif

struct common_evolp1
{
  static const int nnuc = 26;

  float t9;
  float hv;
  float phie;
  arr<float> y;

  common_evolp1() :
    t9(0),
    hv(0),
    phie(0),
    y(dimension(nnuc), fem::fill0)
  {}
};

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_evolp1::nnuc;
#endif

struct common_evolp2
{
  static const int nnuc = 26;

  float dt9;
  float dhv;
  float dphie;
  arr<float> dydt;

  common_evolp2() :
    dt9(0),
    dhv(0),
    dphie(0),
    dydt(dimension(nnuc), fem::fill0)
  {}
};

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_evolp2::nnuc;
#endif

struct common_evolp3
{
  static const int nnuc = 26;

  float t90;
  float hv0;
  float phie0;
  arr<float> y0;

  common_evolp3() :
    t90(0),
    hv0(0),
    phie0(0),
    y0(dimension(nnuc), fem::fill0)
  {}
};

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_evolp3::nnuc;
#endif

struct common_ttime
{
  float t;
  float dt;
  float dlt9dt;

  common_ttime() :
    t(0),
    dt(0),
    dlt9dt(0)
  {}
};

struct common_thermcb
{
  arr<float> thm;
  float hubcst;

  common_thermcb() :
    thm(dimension(14), fem::fill0),
    hubcst(0)
  {}
};

struct common_endens
{
  float rhone0;
  float rhob0;
  float rhob;
  float rnb;

  common_endens() :
    rhone0(0),
    rhob0(0),
    rhob(0),
    rnb(0)
  {}
};

struct common_nucdat
{
  static const int nnuc = 26;

  arr<float> am;
  arr<float> zm;
  arr<float> dm;

  common_nucdat() :
    am(dimension(nnuc), fem::fill0),
    zm(dimension(nnuc), fem::fill0),
    dm(dimension(nnuc), fem::fill0)
  {}
};

#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN
const int common_nucdat::nnuc;
#endif

/*
struct common_besselcb
{
  float bl1;
  float bl2;
  float bl3;
  float bl4;
  float bl5;
  float bm1;
  float bm2;
  float bm3;
  float bm4;
  float bm5;
  float bn1;
  float bn2;
  float bn3;
  float bn4;
  float bn5;

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
  float bk0;
  float bk1;
  float bk2;
  float bk3;
  float bk4;

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

struct common_outdat
{
  static const int itmax = 40;
  static const int nnuc = 26;

  arr<float, 2> xout;
  arr<float, 2> thmout;
  arr<float> t9out;
  arr<float> tout;
  arr<float> dtout;
  arr<float> etaout;
  arr<float> hubout;

  common_outdat() :
    xout(dimension(itmax, nnuc), fem::fill0),
    thmout(dimension(itmax, 6), fem::fill0),
    t9out(dimension(itmax), fem::fill0),
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
  float t9mev;
  float tnmev;
  float tnu;
  float cnorm;
  float rhonu;
  int nu;

  common_nupar() :
    t9mev(0),
    tnmev(0),
    tnu(0),
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
    outfile(fem::bool0)
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
	common_modpr0,
	common_modpr,
	common_recpr0,
	common_recpr,
	common_evolp1,
	common_evolp2,
	common_evolp3,
	common_ttime,
	common_thermcb,
	common_endens,
	common_nucdat,
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
	/*
	fem::variant_core common_rates;
	fem::variant_core common_lncoef;
	fem::cmn_sve check_sve;
	fem::cmn_sve knux_sve;
	fem::cmn_sve rate0_sve;
	fem::cmn_sve xintd_sve;
	fem::cmn_sve rate1_sve;
	fem::cmn_sve start_sve;
	fem::cmn_sve eqslin_sve;
	fem::cmn_sve sol_sve;
	fem::cmn_sve rate2_sve;
	fem::cmn_sve rate3_sve;
	fem::cmn_sve rate4_sve;
	fem::cmn_sve run_sve;
	fem::cmn_sve blockdata_unnamed_sve;
	fem::cmn_sve program_new123_sve;
*/
	common(
	int argc,
	char const* argv[])
	:
	fem::common(argc, argv)
	{}

	//static const int nrec = 88;

	void setcom(common &);
	void setmod(common &);
	void check(common &);
	void driver(common &);
	void run(common &);
};




} // namespace bbn

#endif // GUARD
