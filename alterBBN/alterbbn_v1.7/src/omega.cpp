#include "include.h"

/**--------------------------------------------------------------
 * modifies the model of the QCD equations of state 
 */
void CosmologyModel::Init_modeleff(int model_eff)
{
	int ie,je;
	
	if(model_eff==1)
	{
		const double tableA[276][3]=
		{
#include "sgStar_heff/sgStar_heff_A.tab"
  		};
		for(ie=0;ie<=275;ie++) 
			for(je=0;je<=2;je++) 
				table_eff[ie][je]=tableA[ie][je];
	}
	else if(model_eff==2)
	{
		const double tableB[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) table_eff[ie][je]=tableB[ie][je];
	}
	else if(model_eff==3)
	{
	const double tableB2[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B2.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) table_eff[ie][je]=tableB2[ie][je];
		}
	else if(model_eff==4)
	{
		const double tableB3[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B3.tab"
  		};	
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) table_eff[ie][je]=tableB3[ie][je];
	}
	else if(model_eff==5)
	{
		const double tableC[276][3]=
		{
#include "sgStar_heff/sgStar_heff_C.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) table_eff[ie][je]=tableC[ie][je];
	}
	else
	{
		const double tableold[276][3]=
		{
#include "sgStar_heff/sgStar_heff_old.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) table_eff[ie][je]=tableold[ie][je];
	}
	return;
}

/**--------------------------------------------------------------
 * computes heff at the temperature Temp 
 */
double CosmologyModel::heff(double Temp) const
{
	int ie;
	
	if(Temp>= table_eff[0][0]) return table_eff[0][2];
	
	if(Temp<=0.) return table_eff[275][2];

	ie=1;
	while(Temp<table_eff[ie][0]) ie++;
	
	double heff1,heff2,T1,T2;
	heff1=table_eff[ie][2];
	heff2=table_eff[ie-1][2];
	T1=table_eff[ie][0];
	T2=table_eff[ie-1][0];
			
	return (heff2-heff1)/(T2-T1)*(Temp-T1)+heff1;
}

/**-------------------------------------------------------------
 * computes sgStar at the temperature Temp 
 */
double CosmologyModel::sgStar(double Temp) const
{
	int ie;
	
	if(Temp>= table_eff[0][0]) return table_eff[0][1];
	
	if(Temp<=0.) return table_eff[275][1];

	ie=1;
	while(Temp<table_eff[ie][0]) ie++;

	double sgStar1,sgStar2,T1,T2;
	sgStar1=table_eff[ie][1];
	sgStar2=table_eff[ie-1][1];
	T1=table_eff[ie][0];
	T2=table_eff[ie-1][0];
		
	return (sgStar2-sgStar1)/(T2-T1)*(Temp-T1)+sgStar1;
}

/**--------------------------------------------------------------
 * computes geff at the temperature Temp 
 */
double CosmologyModel::geff(double Temp) const
{
	double heff0=heff(Temp);
	
	return pow(heff0/sgStar(Temp)*(1.+(heff(Temp*1.001)-heff(Temp*0.999))/0.006/heff0),2.);
}

/**--------------------------------------------------------------
 * initializes the parameters contained in paramrelic 
 * TODO replace with constructor?
 */
void CosmologyModel::Init_cosmomodel()
{
	eta0 = 6.19e-10; 			/// baryon-to-photon ratio from WMAP
	life_neutron = 880.1; 		/// neutron lifetime PDG2012
	nbnu = 3; 					/// number of neutrinos
	
	xinu1 = xinu2 = xinu3 = 0; 	/// neutrino chemical potentials?

	fierz = 0;  				/// beta-decay Fierz interference term

	dd0 = ndd = Tdend = 0;
	sd0 = nsd = Tsend = 0;
	nt0 = nnt = Tnend = 0;
	Sigmad0 = nSigmad = TSigmaend = 0;

	Init_modeleff(2);
	
	return;
}

/**--------------------------------------------------------------
 * modifies the values of the baryon-to-photon ratio eta, 
 * the effective number of neutrinos nbnu and 
 * the neutron lifetime life_neutron 
 */
void CosmologyModel::Init_cosmomodel_param(double eta, double nbnu, double life_neutron, double xinu1, double xinu2, double xinu3)
{
	this->eta0=eta;
	this->nbnu=nbnu;
	this->life_neutron=life_neutron;
	this->xinu1=xinu1;
	this->xinu2=xinu2;
	this->xinu3=xinu3;
	return;
}

/**--------------------------------------------------------------
 * modifies the values of the baryon-to-photon ratio eta, 
 * the effective number of neutrinos nbnu and the neutron 
 * lifetime life_neutron and the beta-decay 
 * Fierz interference term
 */
void CosmologyModel::Init_fierz(double eta, double nbnu, double life_neutron, double fierz)
{
    Init_cosmomodel_param(eta, nbnu, life_neutron, 0, 0, 0);
	this->fierz=fierz;
	return;
}

/**--------------------------------------------------------------
 * modifies the parameters of the dark energy density which 
 * appears in the Friedmann equation 
 */
void CosmologyModel::Init_dark_density(double dd0, double ndd, double T_end)
{
	this->dd0=dd0;
	this->ndd=ndd;
	this->Tdend=T_end;
	
	return;
}

/**--------------------------------------------------------------
 * modifies the parameters of the dark entropy density which 
 * appears in the cosmological equations 
 */
void CosmologyModel::Init_dark_entropy(double sd0, double nsd, double T_end)
{
	this->sd0=sd0;
	this->nsd=nsd;
	this->Tsend=T_end;
	
	return;
}

/**--------------------------------------------------------------
 * modifies the parameters of the dark entropy generation which 
 * appears in the cosmological equations
 */
void CosmologyModel::Init_dark_entropySigmaD(double Sigmad0, double nSigmad, double T_end)
{
	this->Sigmad0=Sigmad0;
	this->nSigmad=nSigmad;
	this->TSigmaend=T_end;
	
	return;
}

/**--------------------------------------------------------------
 * modifies the parameters of the non-thermal relic particle 
 * production which appears in the cosmological equations 
 */
void CosmologyModel::Init_nonthermal(double nt0, double nnt, double T_end)
{
	this->nt0=nt0;
	this->nnt=nnt;
	this->Tnend=T_end;
	
	return;
}

/**--------------------------------------------------------------
 * computes the dark energy density at temperature T
 */
double CosmologyModel::dark_density(double T) const
{
	if((dd0==0.)||(T<Tdend)) return 0.;
	
	double rho_rad_1MeV=pi*pi/30.*geff(1.e-3)*1.e-12;
	
	return dd0*rho_rad_1MeV*pow(T/1.e-3,ndd);
}

/**--------------------------------------------------------------
 * computes the dark entropy density at temperature T 
 */
double CosmologyModel::dark_entropy(double T) const
{
	if((sd0==0.)&&(Sigmad0==0.)) return 0.;
	
	if((Sigmad0==0.)&&(T<Tsend)) return 0.;
	
	if(Sigmad0==0.)
	{
		double s_rad_1MeV=2.*pi*pi/45.*heff(1.e-3)*1.e-9;
	
		return sd0*s_rad_1MeV*pow(T/1.e-3,nsd);
	}
	else
	{
		double lnT,dlnT,Ttmp;
		int ie,nmax;
		double integ=0.;
		
		nmax=50;
		
		lnT=log(1.e-15);
		
		dlnT=(log(T)-lnT)/nmax;
		
		for(ie=1;ie<nmax;ie++) 
		{
			lnT+=dlnT;
			Ttmp=exp(lnT);
			integ+=sgStar(Ttmp)*dark_entropy_Sigmad(Ttmp)/sqrt(1.+dark_density(Ttmp)/(pi*pi/30.*geff(Ttmp)*pow(Ttmp,4.)))/pow(heff(Ttmp),2.)/pow(Ttmp,5.);
		}
		
integ+=sgStar(T)*dark_entropy_Sigmad(T)/sqrt(1.+dark_density(T)/(pi*pi/30.*geff(T)*pow(T,4.)))/pow(heff(T),2.)/pow(T,5.)/2.;
		
		integ*=dlnT;

		//double Mplanck=1.2209e19;	

		return 3.*Mplanck*sqrt(5./4./pi/pi/pi)*heff(T)*T*T*T*integ;	
	}
}

/**--------------------------------------------------------------/
 * computes the dark energy entropy derivative at temperature T 
 */
double CosmologyModel::dark_entropy_derivative(double T) const
{
	if((sd0==0.)&&(Sigmad0==0.)) return 0.;
	
	if((Sigmad0==0.)&&(T<Tsend)) return 0.;
	
	if(Sigmad0==0.)
	{
		return (dark_entropy(T*1.001)-dark_entropy(T*0.999))/0.002/T;
	}
	else
	{
	
		//double Mplanck=1.2209e19;	
		return 3.*sgStar(T)/T/heff(T)*(sqrt(geff(T))*dark_entropy(T)-sqrt(5.*Mplanck/4./pi/pi/pi)/T/T*dark_entropy_Sigmad(T)/sqrt(1.+dark_density(T)/(pi*pi/30.*geff(T)*pow(T,4.))));
	}
}

/**--------------------------------------------------------------
 * computes the dark entropy production at temperature T
 */
double CosmologyModel::dark_entropy_Sigmad(double T) const
{
	if((sd0==0.)&&(Sigmad0==0.)) return 0.;
	
	if((Sigmad0==0.)&&(T<Tsend)) return 0.;
	
	//double Mplanck=1.2209e19;
	
	if(Sigmad0==0.)
	{	
		return 1./Mplanck*(sqrt(24.*pi*(pi*pi/30.*geff(T)*pow(T,4.)+dark_density(T)))*dark_entropy(T)-sqrt(4.*pi*pi*pi/45.)*heff(T)/sgStar(T)*pow(T,3)*sqrt(1.+dark_density(T)/(pi*pi/30.*geff(T)*pow(T,4.)))*dark_entropy_derivative(T));
	}
	else
	{
		if(T<TSigmaend) return 0.;

		double s_rad_1MeV=2.*pi*pi/45.*heff(1.e-3)*1.e-9;
	
		double Sigma_rad_1MeV= 1./Mplanck*sqrt(4.*pi*pi*pi/5.*geff(1.e-3))*(1.e-6)*s_rad_1MeV;
	
		return Sigmad0*Sigma_rad_1MeV*pow(T/1.e-3,nSigmad);
	}
}

/**--------------------------------------------------------------
 * computes the non-thermally produced relic particle number 
 * density at temperature T 
 */
double CosmologyModel::nonthermal(double T) const
{
	if((nt0==0.)||(T<Tnend)) return 0.;
	
	return nt0*1.e-50*pow(T/1.e-3,nnt);
}

/**--------------------------------------------------------------
 * computes the neutrino density in case of degeneracies at 
 * temperature Tnu 
 */
double CosmologyModel::neutrino_density(double Tnu) const
{
	if((xinu1==0.)&&(xinu2==0.)&&(xinu3==0.)) return nbnu*2.*pi*pi/30.*pow(Tnu,4.)*7./8.;

	int ie,je,n;
	double rho=0.;
	double xinu[4];
	double max1,max2,int1,int2;
	double x;
	
	xinu[1]=xinu1;
	xinu[2]=xinu2;
	xinu[3]=xinu3;
	
	for(ie=1;ie<=3;ie++)
	{
		if(fabs(xinu[ie])<=0.03) rho+=2.*pi*pi/30.*pow(Tnu,4.)*(7./8.+(15./(4*pi*pi))*xinu[ie]*xinu[ie]+(15./(8.*pow(pi,4.)))*pow(xinu[ie],4.));
		
		else if(fabs(xinu[ie])>=30.) rho+=pow(Tnu,4.)/(8.*pi*pi)*pow(xinu[ie],4.)*(1.+12.*1.645/xinu[ie]/xinu[ie]);
		
		else
		{
			max1=(88.029+xinu[ie])*Tnu;
			int1=0.;
			n=50;
			for(je=1;je<=n-1;je++)
			{
				x=(double)je/(double)n*max1;
				int1+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu-xinu[ie]));
			}
			int1+=0.5*1./(2.*pi*pi)*pow(max1,3.)/(1.+exp(max1/Tnu-xinu[ie]));
			int1*=max1/(double)n;
			rho+=int1;
			
			max2=(88.029-xinu[ie])*Tnu;
			if(max2>0.)
			{
				int2=0.;
				n=50;
				for(je=1;je<=n-1;je++)
				{
					x=(double)je/(double)n*max2;
					int2+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu+xinu[ie]));
				}
				int2+=0.5/(2.*pi*pi)*pow(max2,3.)/(1.+exp(max2/Tnu+xinu[ie]));
				int2*=max2/(double)n;
				rho+=int2;
			}		
		}
	}
	return rho;
}


/**----------------------------------------------------
 * type: 0-10, each type has a unique 
 * (#n1,#n2,#n3,#n4) quartet
 * n1: incoming nuclide number 
 * n2: incoming light nuclide number
 * n3: outgoing light nuclide number
 * n4: outgoing nuclide number
 * rev: reverse reaction coefficient
 * q: energy release in reaction
 */
void CosmologyModel::Init_reactions() //ReactionList & reaction) 
{
	//map<ReactionIndex, Reaction> reaction;
	Reaction r[NNUCREAC] = {
    /// beta decay         type  i    j    k    l    rev    Q   rate        Reaction
		Reaction(n_p,        0,  Nu1, Nu0, Nu0, H1,  0,     0       ),	/// n <-> p
        Reaction(H3_evHe3,   0,  H3,  Nu0, Nu0, He3, 0,     0       ),	/// H3 -> e- + v + He3
        Reaction(Li8_ev2He4, 3,  Li8, Nu0, Nu0, He4, 0,     0       ),	/// Li8 -> e- + v + 2He4
        Reaction(B12_evC12,  0,  B12, Nu0, Nu0, C12, 0,     0       ),	/// B12 -> e- + v + C12
        Reaction(C14_evN14,  0,  C14, Nu0, Nu0, N14, 0,     0       ),	/// C14 -> e- + v + N14
        Reaction(B8_ev2He4,  3,  B8,  Nu0, Nu0, He4, 0,     0       ),	/// B8 -> e+ + v + 2He4
        Reaction(C11_evB11,  0,  C11, Nu0, Nu0, B11, 0,     0       ),	/// C11 -> e+ + v + B11
        Reaction(N12_evC12,  0,  N12, Nu0, Nu0, C12, 0,     0       ),	/// N12 -> e+ + v + C12
        Reaction(N13_evC13,  0,  N13, Nu0, Nu0, C13, 0,     0       ),	/// N13 -> e+ + v + C13
        Reaction(O14_evN14,  0,  O14, Nu0, Nu0, N14, 0,     0       ),	/// O14 -> e+ + v + N14
        Reaction(O15_evN15,  0,  O15, Nu0, Nu0, N15, 0,     0       ),	/// O15 -> e+ + v + N15

    /// reaction           type  i    j    k    l    rev    Q[MeV]          Reaction
        Reaction(H1n_gH2,    1,  H1,  Nu1, Nu0, H2,  0.471, 25.82   ),	/// H1 + n -> g + H2
        Reaction(H2n_gH3,    1,  H2,  Nu1, Nu0, H3,  1.63,  72.62   ),	/// H2 + n -> g + H3
        Reaction(He3n_gHe4,  1,  He3, Nu1, Nu0, He4, 2.61,  238.81  ),	/// He3 + n -> g + He4
        Reaction(Li6n_gLi7,  1,  Li6, Nu1, Nu0, Li7, 1.19,  84.17   ),	/// Li6 + n -> g + Li7
        Reaction(He3n_pH3,   2,  He3, Nu1, H1,  H3,  1.002, 8.863   ),	/// He3 + n -> p + H3
        Reaction(Be7n_pLi7,  2,  Be7, Nu1, H1,  Li7, 0.998, 19.081  ),	/// Be7 + n -> p + Li7
        Reaction(Li6n_tHe4,  2,  Li6, Nu1, H3,  He4, 1.070, 55.494  ),	/// Li6 + n -> d + He4
        Reaction(Be7n_2He4,  4,  Be7, Nu1, Nu0, He4, 4.70,  220.39  ),	/// Be7 + n -> 2He4
        Reaction(H2p_gHe3,   1,  H2,  H1,  Nu0, He3, 1.63,  63.750  ),	/// H2 + p -> g + He3
        Reaction(H3p_gHe4,   1,  H3,  H1,  Nu0, He4, 2.61,  229.932 ),	/// H3 + p -> g + He4
        Reaction(Li6p_gBe7,  1,  Li6, H1,  Nu0, Be7, 1.19,  65.054  ),	/// Li6 + p -> g + Be7
        Reaction(Li6p_hHe4,  2,  Li6, H1,  He3, He4, 1.07,  46.631  ),	/// Li6 + p -> h + He4
        Reaction(Li7p_2He4,  4,  Li7, H1,  Nu0, He4, 4.69,  201.291 ),	/// Li7 + p -> 2He4
        Reaction(He4d_gLi6,  1,  He4, H2,  Nu0, Li6, 1.53,  17.118  ),	/// He4 + d -> g + Li6
        Reaction(He4t_gLi7,  1,  He4, H3,  Nu0, Li7, 1.11,  28.640  ),	/// He4 + t -> g + Li7
        Reaction(He4h_gBe7,  1,  He4, He3, Nu0, Be7, 1.11,  18.423  ),	/// He4 + h -> g + Be7
        Reaction(H2d_nHe3,   5,  H2,  Nu0, Nu1, He3, 1.73,  37.935  ),	/// 2H2 -> n + He3
        Reaction(H2d_pH3,    5,  H2,  Nu0, H1,  H3,  1.73,  46.798  ),	/// 2H2 -> p + H3
        Reaction(H3d_nHe4,   2,  H3,  H2,  Nu1, He4, 5.54,  204.117 ),	/// H3 + H2 -> n + He4
        Reaction(He3d_pHe4,  2,  He3, H2,  H1,  He4, 5.55,  212.980 ),	/// He3 + H2 -> p + He4
        Reaction(He3h_2pHe4, 10, He3, Nu0, H1,  He4, 3.39,  149.230 ),	/// 2He3 -> 2p + He4
        Reaction(Li7d_n2He4, 8,  Li7, H2,  Nu1, He4, 9.95,  175.476 ),	/// Li7 + d -> n + 2He4
        Reaction(Be7d_p2He4, 8,  Be7, H2,  H1,  He4, 9.97,  194.557 ),	/// Be7 + d -> p + 2He4
        Reaction(Li7n_gLi8,  1,  Li7, Nu1, Nu0, Li8, 1.31,  23.59   ),	/// Li7 + n -> g + Li8
        Reaction(B10n_gB11,  1,  B10, Nu1, Nu0, B11, 3.04,  132.95  ),	/// B10 + n -> g + B11
        Reaction(B11n_gB12,  1,  B11, Nu1, Nu0, B12, 2.34,  39.10   ),	/// B11 + n -> g + B12
        Reaction(C11n_pB11,  2,  C11, Nu1, H1,  B11, 1.002, 32.080  ),	/// C11 + n -> p + B11
        Reaction(B10n_aLi7,  2,  B10, Nu1, He4, Li7, 0.758, 32.382  ),	/// B10 + n -> a + Li7
        Reaction(Be7p_gB8,   1,  Be7, H1,  Nu0, B8,  1.30,  1.595   ),	/// Be7 + p -> g + B8
        Reaction(Be9p_gB10,  1,  Be9, H1,  Nu0, B10, 0.973, 76.427  ),	/// Be9 + p -> g + B10
        Reaction(B10p_gC11,  1,  B10, H1,  Nu0, C11, 3.03,  100.840 ),	/// B10 + p -> g + C11
        Reaction(B11p_gC12,  1,  B11, H1,  Nu0, C12, 7.01,  185.173 ),	/// B11 + p -> g + C12
        Reaction(C11p_gN12,  1,  C11, H1,  Nu0, N12, 2.33,  6.975   ),	/// C11 + p -> g + N12
        Reaction(B12p_nC12,  2,  B12, H1,  Nu1, C12, 3.00,  146.08  ),	/// B12 + p -> n + C12
        Reaction(Be9p_aLi6,  2,  Be9, H1,  He4, Li6, 0.618, 24.674  ),	/// Be9 + p -> a + Li6
        Reaction(B10p_aBe7,  2,  B10, H1,  He4, Be7, 0.754, 13.301  ),	/// B10 + p -> a + Be7
        Reaction(B12p_aBe9,  2,  B12, H1,  He4, Be9, 0.292, 79.89   ),	/// B12 + p -> a + Be9
        Reaction(Li6a_gB10,  1,  Li6, He4, Nu0, B10, 1.58,  51.753  ),	/// Li6 + a -> g + B10
        Reaction(Li7a_gB11,  1,  Li7, He4, Nu0, B11, 4.02,  100.538 ),	/// Li7 + a -> g + B11
        Reaction(Be7a_gC11,  1,  Be7, He4, Nu0, C11, 4.02,  87.539  ),	/// Be7 + a -> g + C11
        Reaction(B8a_pC11,   2,  B8,  He4, H1,  C11, 3.08,  86.00   ),	/// B8 + a -> p + C11 
        Reaction(Li8a_nB11,  2,  Li8, He4, Nu1, B11, 3.07,  76.96   ),	/// Li8 + a -> n + B11
        Reaction(Be9a_nC12,  2,  Be9, He4, Nu1, C12, 10.3,  66.160  ),	/// Be9 + a -> n + C12
        Reaction(Be9d_nB10,  2,  Be9, H2,  Nu1, B10, 2.07,  50.63   ),	/// Be9 + d -> n + B10
        Reaction(B10d_pB11,  2,  B10, H2,  H1,  B11, 6.44,  107.13  ),	/// B10 + d -> p + B11
        Reaction(B11d_nC12,  2,  B11, H2,  Nu1, C12, 14.9,  159.36  ),	/// B11 + d -> n + C12
        Reaction(He4an_gBe9, 7,  He4, Nu1, Nu0, Be9, 0.584, 18.260  ),	/// 2He4 + n -> g + Be9
        Reaction(He4aa_gC12, 6,  He4, Nu0, Nu0, C12, 2.00,  84.420  ),	/// 3He4 -> g + C12
        Reaction(Li8p_n2He4, 8,  Li8, H1,  Nu1, He4, 3.58,  177.73  ),	/// Li8 + p -> n + 2He4
        Reaction(B8n_p2He4,  8,  B8,  Nu1, H1,  He4, 3.58,  218.82  ),	/// B8 + n -> p + 2He4
        Reaction(Be9p_d2He4, 8,  Be9, H1,  H2,  He4, 0.807, 7.555   ),	/// Be9 + p -> d + 2He4
        Reaction(B11p_3He4,  9,  B11, H1,  Nu0, He4, 3.50,  100.753 ),	/// B11 + p -> 3He4
        Reaction(C11n_3He4,  9,  C11, Nu1, Nu0, He4, 3.49,  132.83  ),	/// C11 + n -> 3He4
        Reaction(C12n_gC13,  1,  C12, Nu1, Nu0, C13, 0.886, 57.41   ),	/// C12 + n -> g + C13
        Reaction(C13n_gC14,  1,  C13, Nu1, Nu0, C14, 3.58,  94.88   ),	/// C13 + n -> g + C14
        Reaction(N14n_gN15,  1,  N14, Nu1, Nu0, N15, 2.71,  125.74  ),	/// N14 + n -> g + N15
        Reaction(N13n_pC13,  2,  N13, Nu1, H1,  C13, 1.002, 34.846  ),	/// N13 + n -> p + C13
        Reaction(N14n_pC14,  2,  N14, Nu1, H1,  C14, 3.003, 7.263   ),	/// N14 + n -> p + C14
        Reaction(O15n_pN15,  2,  O15, Nu1, H1,  N15, 1.002, 41.037  ),	/// O15 + n -> p + N15
        Reaction(O15n_aC12,  2,  O15, Nu1, He4, C12, 0.709, 98.661  ),	/// O15 + n -> a + C12
        Reaction(C12p_gN13,  1,  C12, H1,  Nu0, N13, 0.884, 22.553  ),	/// C12 + p -> g + N13
        Reaction(C13p_gN14,  1,  C13, H1,  Nu0, N14, 1.19,  87.621  ),	/// C13 + p -> g + N14
        Reaction(C14p_gN15,  1,  C14, H1,  Nu0, N15, 0.900, 118.452 ),	/// C14 + p -> g + N15
        Reaction(N13p_gO14,  1,  N13, H1,  Nu0, O14, 3.57,  53.706  ),	/// N13 + p -> g + O14
        Reaction(N14p_gO15,  1,  N14, H1,  Nu0, O15, 2.70,  84.678  ),	/// N14 + p -> g + O15
        Reaction(N15p_gO16,  2,  N15, H1,  Nu0, O16, 3.62,  140.734 ),	/// N15 + p -> g + O16
        Reaction(N15p_aC12,  2,  N15, H1,  He4, C12, 0.706, 57.623  ),	/// N15 + p -> a + C12
        Reaction(C12a_gO16,  1,  C12, He4, Nu0, O16, 5.13,  83.111  ),	/// C12 + a -> g + O16
        Reaction(B10a_pC13,  2,  B10, He4, H1,  C13, 9.36,  47.16   ),	/// B10 + a -> p + C13
        Reaction(B11a_pC14,  2,  B11, He4, H1,  C14, 11.0,  9.098   ),	/// B11 + a -> p + C14
        Reaction(C11a_pN14,  2,  C11, He4, H1,  N14, 3.68,  33.915  ),	/// C11 + a -> p + N14
        Reaction(N12a_pO15,  3,  N12, He4, H1,  O15, 4.26,  111.87  ),	/// N12 + a -> p + O15
        Reaction(N13a_pO16,  3,  N13, He4, H1,  O16, 5.81,  60.557  ), 	/// N13 + a -> p + O16
        Reaction(B10a_nN13,  2,  B10, He4, Nu1, N13, 9.34,  12.287  ),	/// B10 + a -> n + N13
        Reaction(B11a_nN14,  2,  B11, He4, Nu1, N14, 3.67,  1.835   ),	/// B11 + a -> n + N14
        Reaction(B12a_nN15,  2,  B12, He4, Nu1, N15, 4.25,  88.47   ),	/// B12 + a -> n + N15
        Reaction(C13a_nO16,  2,  C13, He4, Nu1, O16, 5.79,  25.711  )	/// C13 + a -> n + O16
    };

	// TODO replace with copy contructor above and remove below
    for (int n=0; n<NNUCREAC; n++)
    {
        ReactionIndex id = r[n].id;
        reactions[id].id = id;
        reactions[id].type = r[n].type;
        reactions[id].in_major = r[n].in_major;
        reactions[id].in_minor = r[n].in_minor;   // TODO don't set to zero.
        reactions[id].out_minor = r[n].out_minor; // TODO don't set to zero.
        reactions[id].out_major = r[n].out_major;
        reactions[id].reverse = r[n].reverse;
        reactions[id].forward = r[n].forward;
		assert(r[n].type < 11);
    }
	//printf("Reaction Index Overflow: %d\n", ReactionIndexOverflow);
	//printf("C13a_nO16: %d\n", C13a_nO16);
	//printf("NNUCREAC: %d\n", NNUCREAC);
}


/**----------------------------------------------------
 *  Data from http://www.nndc.bnl.gov/masses/mass.mas03
 *
 *  S - isotopic symbol
 *  A - atomic mass number
 *  Z - atomic charge number
 *  N - neutron number
 *  M - mass excess 
 *-----------------------------------------------------*/
void CosmologyModel::Init_nuclides() 
{
	///  				   S    name    A   Z   N   dm
	nuclides[Nu0] = Nuclide(Nu0, "none", 0,  0,  0,  0 );
	nuclides[Nu1] = Nuclide(Nu1, "n",    1,  0,  1,  8.071388);
	nuclides[H1]  = Nuclide(H1,  "H",    1,  1,  0,  7.289028);
	nuclides[H2]  = Nuclide(H2,  "D",    1,  1,  1,  13.135825);
	nuclides[H3]  = Nuclide(H3,  "T",    2,  1,  2,  14.949915);
	nuclides[He3] = Nuclide(He3, "He3",  3,  2,  1,  14.931325);
	nuclides[He4] = Nuclide(He4, "He4",  4,  2,  2,  2.424931);

	nuclides[Li6] = Nuclide(Li6, "Li6",  6,  3,  3,  14.9078);
	nuclides[Li7] = Nuclide(Li7, "Li7",  7,  3,  4,  15.7696);
	nuclides[Be7] = Nuclide(Be7, "Be7",  7,  4,  3,  15.7696);
	nuclides[Li8] = Nuclide(Li8, "Li8",  8,  3,  5,  20.9464);

	nuclides[B8]  = Nuclide(B8,  "B8",   8,  5,  3,  22.9212);
	nuclides[Be9] = Nuclide(Be9, "Be9",  9,  4,  5,  11.34758);
	nuclides[B10] = Nuclide(B10, "B10",  10, 5,  5,  12.05086);
	nuclides[B11] = Nuclide(B11, "B11",  11, 5,  6,  8.6680);

	nuclides[C11] = Nuclide(C11, "C11",  11, 6,  5,  10.6506);
	nuclides[B12] = Nuclide(B12, "B12",  12, 5,  7,  13.3690);
	nuclides[C12] = Nuclide(C12, "C12",  12, 6,  6,  0);
	nuclides[N12] = Nuclide(N12, "N12",  12, 7,  5,  17.3382);

	nuclides[C13] = Nuclide(C13, "C13",  13, 6,  7,  3.125036);
	nuclides[N13] = Nuclide(N13, "N13",  14, 7,  6,  5.3455);
	nuclides[C14] = Nuclide(C14, "C14",  14, 6,  8,  3.019916);
	nuclides[N14] = Nuclide(N14, "N14",  14, 7,  7,  2.863440);

	nuclides[O14] = Nuclide(O14, "O14",  14, 8,  6,  8.006521);
	nuclides[N15] = Nuclide(N15, "N15",  15, 7,  8,  0.101439);
	nuclides[O15] = Nuclide(O15, "O15",  15, 8,  7,  2.8554);
	nuclides[O16] = Nuclide(O16, "O16",  16, 8,  8, -4.737036);
}
