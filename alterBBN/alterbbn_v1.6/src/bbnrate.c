#include "include.h"
#include "assert.h"

void rate_weak(int err, double f[])
/* calculates the nuclear forward rates of weak interaction reactions */
/* err=0: central values; err=1: high values; err=2: low values; err>100000: random gaussian error; err<0: error only for process number (-err) */
{
    //double Nbeta = 11;
    int first = t_evh;
    int last = O15_evN15;
	double ferrlow[last+1]; 
    double ferrhigh[last+1]; 
    double ferr[last+1];
	
	f[t_evh] = 1.78141141239e-9;	/// H3 -> e- + v + He3
	f[Li8_evaa]  = 0.827;	 	    /// Li8 -> e- + v + 2He4
	f[B12_evC12] = 34.3;	 	    /// B12 -> e- + v + C12
	f[C14_evN14] = 3.834e-12;	 	/// C14 -> e- + v + N14
	f[B8_evaa]   = 0.9;	 	        /// B8 -> e+ + v + 2He4
	f[C11_evB11] = 5.668e-4;	 	/// C11 -> e+ + v + B11
	f[N12_evC12] = 63.01;	 	    /// N12 -> e+ + v + C12 
	f[N13_evC13] = 0.001159;	 	/// N13 -> e+ + v + C13
	f[O14_evN14] = 0.0098171;	    /// O14 -> e+ + v + N14 
	f[O15_evN15] = 0.0056704; 	    /// O15 -> e+ + v + N15 

	if((err > 0) || (err <= -first && err >= -last))
	{
		ferrhigh[first] = 5.e-3;
		ferrlow[first] = -5.e-3;
		ferr[first] = 5.e-3;
		
	    int reac;
		for(reac = first+1; reac <= last; reac++)
		{
			ferrhigh[reac] = 9.;
			ferrlow[reac] = -0.9;
			ferr[reac] = ferrhigh[reac];
		}
		
		if(err == 1) 
		    for(reac = first; reac <= last; reac++)
                f[reac] *= fabs(1+ferrhigh[reac]);
		if(err == 2) 
		    for(reac = first; reac <= last; reac++)
                f[reac] *= fabs(1+ferrlow[reac]);
				
		if(err>100000)
		{
			srand((unsigned int)(getpid()+err));
		    for(reac = first; reac <= last; reac++)
                f[reac] *= fabs(1+ferr[reac]*rand_gauss());
		}
		
		if(err<0) 
            f[-err] *= fabs(1+ferr[-err]);
	}

	return;
}


/*----------------------------------------------------*

double weak_phase_space(double Ex, double Ee, double Enu)
{
    if(x>1.)
    {
        (x+b)*pow(Enu,2.)*sqrt(x*x-1.)
            /(1.+exp(-Ee/z9))
            /(1.+exp(-Enu/znu-xi1));
    }
    return 0;
}
*/

/*----------------------------------------------------*/


void rate_pn(int err, struct relicparam paramrelic, double f[], double r[], double T9, double Tnu)
/* calculates the nuclear forward and reverse rates f[] and r[] of the reaction p <-> n at the temperature T9 */
/*err=0: central values; err=1: high values; err=2: low values; err>100000: random gaussian error; err<0: error only for process number (-err)*/
{
	double ferr,rerr;
	ferr = rerr = 0.;
    double b = paramrelic.fierz;            /// beta-decay Fierz interference term */
	double tau = paramrelic.life_neutron;   /// measured neutron lifetime at T=0 in s */
    double xi1 = paramrelic.xinu1;          /// neutrino chemical potential
	
	if(((xi1==0)&&(b==0))||(Tnu==0))
	{
		int i;
		double z9 = T9*kB/me;               /// inverse of old z def in v1.4
	
		double af[13] = {
            0.15735,0.46172e1,-0.40520e2,0.13875e3,
            -0.59898e2,0.66752e2,-0.16705e2,0.38071e1,
            -0.39140,0.23590e-1,-0.83696e-4,-0.42095e-4,0.17675e-5};
	
		double bf[10] = {
            0.22211e2,-0.72798e2,0.11571e3,-0.11763e2,
            0.45521e2,-0.37973e1,0.41266e0,-0.26210e-1,0.87934e-3,-0.12016e-4};
	
		f[1]=1.;
		for(i=1;i<=13;i++) 
            f[1]+=af[i-1]*pow(z9,i);
		f[1]*=exp(-0.33979*z9)/tau;     /// n->p
		
		//if(1/z9 < 5.10998997931) // almost 10*me, but why?
		//if(1/z9 < 10) // maybe 10 works instead!? 
		if(z9 > 0.1) // maybe 0.1 works instead!? 
		{
			r[1]=-0.62173;
			for(i=1;i<=10;i++)
                r[1]+=bf[i-1]*pow(z9,i);
			r[1]*=exp(-2.8602/z9)/tau;  /// n->p
		}
		else 
            r[1]=0.;                    /// p->n
	}
	else
	{		
		double T9mev = T9*kB;
		double Tnumev = Tnu*kB;
		double z9 = T9*kB/me;
		double znu = Tnu*kB/me;
		double q = 1.29333217/me;   /// q=(mn-mp)/me

		double int1 = 0;
 		double int2 = 0;
		double int3 = 0;
		double int4 = 0;
 		int n = 1000;
        double minexp = 10;
		double x;
		int je;

        /// Integral 0: Normalizer for n->p rate
		double zmax = max(50*z9,q);
        double norm = 0;
		for(je=1;je<=n-1;je++)
		{
			x = 1+(double)je/(double)n*(q-1);
			norm += (x+b)*pow(x-q,2.)*sqrt(x*x-1);
		}
		norm *= (q-1)/(double)n;    /// PDF normalizer
		norm *= tau/0.9805;         /// lifetime and Fermi function correction

        /// Integral 1: 1st part of n->p rate
		double max1 = max(zmax,fabs(znu*(50.+xi1)+q));
		for(je=1;je<=n-1;je++)
		{
			x = 1+(double)je/(double)n*(max1-1);
			if(x>1)
			{
				int1+=(x+b)*pow(x-q,2.)*sqrt(x*x-1)
                    /(1+exp(-x/z9))
                    /(1+exp((x-q)/znu-xi1));
			}
		}
		if(max1>1) 
            int1+=0.5*(max1+b)*pow(max1-q,2.)*sqrt(max1*max1-1)
                /(1+exp(-max1/z9))
                /(1+exp((max1-q)/znu-xi1));
		int1*=(max1-1)/(double)n;

        /// Integral 2: 2nd part of n->p rate
		double max2=max(zmax,fabs(znu*(50.-xi1)-q));
		for(je=1;je<=n-1;je++)
		{
			x=1+(double)je/(double)n*(max2-1);
			if(x>1)
			{
				int2+=(x-b)*pow(x+q,2.)*sqrt(x*x-1)
                    /(1+exp(x/z9))
                    /(1+exp(-(x+q)/znu-xi1));
			}
		}
		if(max2>1) 
            int2+=0.5*(max2-b)*pow(max2+q,2.)*sqrt(max2*max2-1)
                /(1+exp(max2/z9))
                /(1+exp(-(max2+q)/znu-xi1));
		int2*=(max2-1)/(double)n;

        /// Integral 3: 1st part of p->n rate
		double max3=max(zmax,fabs(znu*(50.-xi1)-q));
		for(je=1;je<=n-1;je++)
		{
			x=1+(double)je/(double)n*(max3-1);
			if(x>1)
			{
				int3+=(x-b)*pow(x+q,2.)*sqrt(x*x-1)
                    /(1+exp(-x/z9))
                    /(1+exp((x+q)/znu+xi1));
			}
		}
		if(max3>1) 
            int3+=0.5*(max3-b)*pow(max3+q,2.)*sqrt(max3*max3-1)
                /(1+exp(-max3/z9))
                /(1+exp((max3+q)/znu+xi1));
		int3*=(max3-1)/(double)n;
 
        /// Integral 4: 2nd part of p->n rate
		double max4=max(zmax,fabs(znu*(50.+xi1)+q));
		for(je=1;je<=n-1;je++)
		{
			x=1+(double)je/(double)n*(max4-1);
			if(x>1)
			{
				int4+=(x+b)*pow(x-q,2.)*sqrt(x*x-1)
                    /(1+exp(x/z9))
                    /(1+exp(-(x-q)/znu+xi1));
			}
		}
		if(max4>1) 
            int4+=0.5*(max4+b)*pow(max4-q,2.)*sqrt(max4*max4-1)
                /(1+exp(max4/z9))
                /(1+exp(-(max4-q)/znu+xi1));
		int4*=(max4-1)/(double)n;

		f[1]=(int1+int2)/norm;
        r[1]=(int3+int4)/norm;
	}

	if((err!=0&&err!=-10000))
	{
		ferr=0.0006;
		rerr=0.0006;	
	}
	
	if(err==1) 
	{
		f[1]*=fabs(1+ferr);
		r[1]*=fabs(1+rerr);
	}
	
	if(err==2) 
	{
		f[1]*=fabs(1-ferr);
		r[1]*=fabs(1-rerr);
	}
	
	if(err>100000)
	{
        printf("Error: rates are going to be random.");
			srand((unsigned int)(getpid()+err));
			f[1]*=fabs(1+ferr*rand_gauss());
			r[1]*=fabs(1+rerr*rand_gauss());
	}
		
	if(err==-1)
	{
		f[1]*=fabs(1-ferr);
		r[1]*=fabs(1-rerr);
	}		
	
	return;
}

/*----------------------------------------------------*/

void rate_all(int err, double f[], double T9)
/* calculates the nuclear forward rates f[] of the nuclear reactions given below at the temperature T9 */
/* err=0: central values; err=1: high values; err=2: low values; err>100000: random gaussian error; err<0: error only for process number (-err) */
{
    int first = pn_gd;
    int last = C13a_nO16; 
	double ferrlow[last+1];
    double ferrhigh[last+1];
    double ferr[last+1];
	int reac;

	if((err!=0)&&(err!=-10000)) 
        for(reac = first; reac <= last; reac++)
        {
            ferrhigh[reac]=9.;
            ferrlow[reac]=-0.9;
            ferr[reac]=ferrhigh[reac];
        }

/* H + n -> g + H2 */
	if(T9<1.5) 
        f[pn_gd] = 44060.*(1.-2.7503695564153143*T9-3.5220409333117897*T9*T9-0.2093513619089196*T9*T9*T9+0.10659679579058313*sqrt(T9)+ 4.62948586627009*pow(T9,3./2.) +1.3459574632779876*pow(T9,5./2.));
	else
        f[pn_gd] = (1.-sqrt(T9)*0.8504+T9*0.4895-pow(T9,1.5)*0.09623+T9*0.008471*T9-T9*2.8e-4*pow(T9,1.5))*47420.;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<1.5) 
            ferrlow[pn_gd] = -(67.95212182165885*sqrt(1.255364135304537+618.8966688079934*T9+21793.605102078516*T9*T9+97573.42092925594*T9*T9*T9+85539.62453206794*pow(T9,4.)+14687.125326371255*pow(T9,5.)+247.57751149254238*pow(T9,6.)-2855.9302137640752*pow(T9,11./2.)-25.82889427679147*pow(T9,1./2.)-5043.279674275595*pow(T9,3./2.)-57307.155278753264*pow(T9,5./2.)-111071.90152148822*pow(T9,7./2.)-44154.571125115704*pow(T9,9./2.)))/f[pn_gd];
		else 
            ferrlow[pn_gd]=-0.078;
		
		ferrhigh[pn_gd]=-ferrlow[pn_gd];
		ferr[pn_gd]=ferrlow[pn_gd];
	} 
	
/* H2 + n -> g + H3 */
	f[dn_gt]=(T9*18.9+1.)*66.2;
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[dn_gt]=0.3;
		ferrlow[dn_gt]=-0.3;
		ferr[dn_gt]=ferrlow[dn_gt];
	}
	
/* He3 + n -> g + He4 */
	f[hn_ga]=(T9*905.+1.)*6.62;
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[hn_ga]=1.;
		ferrlow[hn_ga]=-0.5;
		ferr[hn_ga]=ferrlow[hn_ga];
	}
	
/* Li6 + n -> g + Li7 */
	f[Li6n_gLi7]=5100.;
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Li6n_gLi7]=0.3;
		ferrlow[Li6n_gLi7]=-0.3;
		ferr[Li6n_gLi7]=ferrlow[Li6n_gLi7];
	}
	
/* He3 + n -> p + H3 */
	if(T9<2.5) 
        f[hn_pt]=7.064935e8+6.733213571736319e8*T9+1.7181155480346258e9*T9*T9-4.5367658146835446e8*T9*T9*T9-1.2216728981712557e8*pow(T9,4.)-4.92736677238425e8*sqrt(T9)-1.3659670893994067e9*pow(T9,3./2.)-6.629932739639357e8*pow(T9,5./2.)+4.834951929033479e8*pow(T9,7./2.);
	else 
        f[hn_pt]=4.81732e8;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<2.5) 
            ferrlow[hn_pt]=-sqrt(0.00187*0.00187+5.4*pow(sqrt(3.51e11+3.112097416989299e11*T9+7.901575411070865e10*T9*T9+2.04541951561715e10*T9*T9*T9+5.110885965380451e9*pow(T9,4.)+3.9016700171412725e9*pow(T9,5.)+1.2106464640648174e9*pow(T9,6.)+2.842691804858251e8*pow(T9,7.)+2.5025023636054292e8*pow(T9,8.)-1.0919522573895195e9*pow(T9,11./2.)-5.074476577064073e11*pow(T9,1./2.)-9.073561744271307e8*pow(T9,13./2.)-4.935780126698165e8*pow(T9,15./2.)-1.3272119856586942e11*pow(T9,3./2.)-3.982502921484235e10*pow(T9,5./2.)-1.4832025658250046e10*pow(T9,7./2.)-3.9093487936349277e9*pow(T9,9./2.))/f[hn_pt],2.));
		else 
            ferrlow[hn_pt]=-sqrt(0.00187*0.00187+5.4*pow(224626./f[hn_pt],2.));
		
		ferrhigh[hn_pt]=-ferrlow[hn_pt];
		ferr[hn_pt]=ferrlow[hn_pt];
	}
	
/* Be7 + n -> p + Li7 */
	if(T9<2.5) 
        f[Be7n_pLi7]=6.8423032e9+1.7674863e10*T9+2.6622006e9*T9*T9-3.3561608e8*T9*T9*T9-5.9309139e6*pow(T9,4.)-1.4987996e10*sqrt(T9)-1.0576906e10*pow(T9,3./2.)+2.7447598e8*pow(T9,5./2.)+7.6425157e7*pow(T9,7./2.)-2.282944e7*pow(T9,-3./2.)/exp(0.050351813/T9);
	else 
        f[Be7n_pLi7]=1.28039e9;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<2.5)
		{
			ferrhigh[Be7n_pLi7]=sqrt(0.02082*0.02082+1.2*pow((6.8637241e9+1.7654706e10*T9+2.6528695e9*T9*T9-3.334255e8*T9*T9*T9-5.8841835e6*pow(T9,4.)-1.4998087e10*pow(T9,1./2.)-1.0549873e10*pow(T9,3./2.)+2.7244147e8*pow(T9,5./2.)+7.5882824e7*pow(T9,7./2.)-(2.2946239e7*pow(T9,-3./2.))/exp(0.05042127/T9))/f[Be7n_pLi7]-1.,2.));
			ferrlow[Be7n_pLi7]=-sqrt(0.02082*0.02082+1.2*pow((5.3356377e9+1.2258644e10*T9+1.6991251e9*T9*T9-2.011188e8*T9*T9*T9-3.0052814e6*pow(T9,4.)-1.0649687e10*pow(T9,1./2.)-7.1221438e9*pow(T9,3./2.)+1.8848519e8*pow(T9,5./2.)+4.2410535e7*pow(T9,7./2.)+(6.0102571e7*pow(T9,-3./2.))/exp(0.2761375/T9))/f[Be7n_pLi7]-1.,2.));
		}
		else
		{
			ferrhigh[Be7n_pLi7]=sqrt(0.02082*0.02082+1.2*pow(1.28652e9/f[Be7n_pLi7]-1.,2.));
			ferrlow[Be7n_pLi7]=-sqrt(0.02082*0.02082+1.2*pow(1.27454e9/f[Be7n_pLi7]-1.,2.));
		}
		ferr[Be7n_pLi7]=ferrlow[Be7n_pLi7];
	}
	
/* Li6 + n -> a + H3 */	
	f[Li6n_at]=pow(T9,-1.5)*2.54e9*exp(-2.39/T9)+(1.-pow(T9/(T9*49.18+1.),1.5)*0.261/pow(T9,1.5))*1.68e8;
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Li6n_at]=0.1;
		ferrlow[Li6n_at]=-0.1;
		ferr[Li6n_at]=ferrlow[Li6n_at];
	}

/* Be7 + n -> a + He4 */
	f[Be7n_aa]=(T9*3760.+1.)*20500.;
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Be7n_aa]=0.9;
		ferrlow[Be7n_aa]=-0.9;
		ferr[Be7n_aa]=ferrlow[Be7n_aa];
	}


/* H2 + p -> g + He3 */	
	if (T9<4.) 
        f[dp_gh]=((-15.709674-721.9142*T9+173.23945*T9*T9-77.371692*T9*T9*T9+126.82087e0*pow(T9,1./3.)-206.50853*pow(T9,2./3.)+2120.7339*pow(T9,4./3.)-369.61306*pow(T9,5./3.)+127.8378*pow(T9,7./3.)+100.68769*pow(T9,8./3.))*pow(T9,-2./3.))/exp(1.29042942e0*pow(T9,-1./3.));
	else 
        f[dp_gh]=2049.72;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<4.)
		{
			ferrhigh[dp_gh]=sqrt(0.02235*0.02235+1.2*pow((((-20.607762-651.42468*T9+690.01011*T9*T9-28.933696*T9*T9*T9+134.27673*pow(T9,1./3.)-148.86271*pow(T9,2./3.)+1513.5632*pow(T9,4./3.)-668.14865*pow(T9,5./3.)-11.831976*pow(T9,7./3.)-1.7647978*pow(T9,8./3.))*pow(T9,-2./3.))/exp(0.82258711*pow(T9,-1./3.)))/f[dp_gh]-1.,2.));
			ferrlow[dp_gh]=-sqrt(0.02235*0.02235+1.2*pow((((-10.782957+32828.304*T9-65943.049*T9*T9+679.86033*T9*T9*T9+727.88161*pow(T9,1./3.)-7736.0309*pow(T9,2./3.)-64848.454*pow(T9,4./3.)+84983.958*pow(T9,5./3.)+30784.918*pow(T9,7./3.)-7555.1171*pow(T9,8./3.))*pow(T9,-2./3.))/exp(2.44792659*pow(T9,-1./3.)))/f[dp_gh]-1.,2.));
		}
		else
		{
			ferrhigh[dp_gh]=sqrt(0.02235*0.02235+1.2*pow(2104.68/f[dp_gh]-1.,2.));
			ferrlow[dp_gh]=-sqrt(0.02235*0.02235+1.2*pow(1994.63/f[dp_gh]-1.,2.));
		}
		ferr[dp_gh]=ferrlow[dp_gh];		
	}
	
/* H3 + p -> g + He4 */	
    f[tp_ga]=pow(T9,-2./3.)*2.2e4*exp(-3.869/pow(T9,1./3.))*(pow(T9,1./3.)*0.108+1.+pow(T9,2./3.)*1.68+T9*1.26+pow(T9,4./3.)*0.551+pow(T9,5./3.)*1.06);
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[tp_ga]=0.2;
		ferrlow[tp_ga]=-0.2;
		ferr[tp_ga]=ferrlow[tp_ga];
	}

/* Li6 + p -> g + Be7 */
	f[Li6p_gBe7] = 1.25e6*pow(T9,-2./3.)*exp(-8.415/pow(T9,1./3.))*(1.-0.252*T9+5.19e-2*T9*T9-2.92e-3*T9*T9*T9);
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Li6p_gBe7]=0.21;
		ferrlow[Li6p_gBe7]=-0.21;
		ferr[Li6p_gBe7]=ferrlow[Li6p_gBe7];
	}

/* Li6 + p -> a + He3 */
	if (T9<2.5) 
        f[Li6p_ah]=((-7.4966212e7-1.9411561e10*T9+1.6262854e10*T9*T9+2.0533495e7*pow(T9,1./3.)+3.9547491e9*pow(T9,2./3.)+3.7907358e10*pow(T9,4./3.)-3.4313768e10*pow(T9,5./3.)-3.9965228e9*pow(T9,7./3.)+4.0333873e8*pow(T9,8./3.))*pow(T9,-2./3.))/exp(4.62619323*pow(T9,-1./3.));
	else f[Li6p_ah]=3.05102e7;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<2.5)
		{
			ferrhigh[Li6p_ah]=sqrt(0.0931*0.0931+1.2*pow((((4.6794127e7-1.5077363e9*T9-2.8610381e9*T9*T9-4.1883216e8*pow(T9,1./3.)+1.3422134e9*pow(T9,2./3.)-9.4597359e8*pow(T9,4./3.)+3.6073249e9*pow(T9,5./3.)+9.4073567e8*pow(T9,7./3.)-1.1547155e8*pow(T9,8./3.))*pow(T9,-2./3.))/exp(2.47110932*pow(T9,-1./3.)))/f[Li6p_ah]-1.,2.));
			ferrlow[Li6p_ah]=-sqrt(0.0931*0.0931+1.2*pow((((-2.1541443e7-2.1453941e10*T9+1.5165239e10*T9*T9-4.5939493e8*pow(T9,1./3.)+5.5208615e9*pow(T9,2./3.)+3.8266784e10*pow(T9,4./3.)-3.3068204e10*pow(T9,5./3.)-3.6300979e9*pow(T9,7./3.)+3.5841419e8*pow(T9,8./3.))*pow(T9,-2./3.))/exp(4.53459377*pow(T9,-1./3.)))/f[Li6p_ah]-1.,2.));
		}
		else
		{
			ferrhigh[Li6p_ah]=sqrt(0.0931*0.0931+1.2*pow(3.09863e7/f[Li6p_ah]-1.,2.));
			ferrlow[Li6p_ah]=-sqrt(0.0931*0.0931+1.2*pow(3.00895e7/f[Li6p_ah]-1.,2.));
		}
		ferr[Li6p_ah]=ferrlow[Li6p_ah];
	}

/* Li7 + p -> a + He4 */
	if (T9<2.5)
	{
		f[Li7p_aa]=((-8.9654123e7-2.5851582e8*T9-2.6831252e7*T9*T9+3.8691673e8*pow(T9,1./3.)+4.9721269e8*pow(T9,2./3.)+2.6444808e7*pow(T9,4./3.)-1.2946419e6*pow(T9,5./3.)-1.0941088e8*pow(T9,7./3.)+9.9899564e7*pow(T9,8./3.))*pow(T9,-2./3.))/exp(7.73389632*pow(T9,-1./3.));
		f[Li7p_aa]+=exp(-1.137519e0*T9*T9-8.6256687*pow(T9,-1./3.))*(3.0014189e7-1.8366119e8*T9+1.7688138e9*T9*T9-8.4772261e9*T9*T9*T9+2.0237351e10*pow(T9,4.)-1.9650068e10*pow(T9,5.)+7.9452762e8*pow(T9,6.)+1.3132468e10*pow(T9,7.)-8.209351e9*pow(T9,8.)-9.1099236e8*pow(T9,9.)+2.7814079e9*pow(T9,10.)-1.0785293e9*pow(T9,11.)+1.3993392e8*pow(T9,12.))*pow(T9,-2./3.);
	}
	else
	{
		f[Li7p_aa]=1.53403e6;
		f[Li7p_aa]+=84516.7;
	}
	if((err!=0)&&(err!=-10000))
	{
		if(T9<2.5)
		{
			ferrhigh[Li7p_aa]=sqrt(0.08*0.08+pow((((1.6425644e7-7.682657e8*T9+1.2461811e9*T9*T9-1.1914365e8*pow(T9,1./3.)+3.3659333e8*pow(T9,2./3.)+1.8234158e9*pow(T9,4./3.)-1.9962683e9*pow(T9,5./3.)-5.4978741e8*pow(T9,7./3.)+1.4214466e8*pow(T9,8./3.))*pow(T9,-2./3.))/exp(6.34172901*pow(T9,-1./3.))+exp(-1.0418442*T9*T9-5.5570697*pow(T9,-1./3.))*(-25145.507+1.0787318e6*T9-1.5899728e7*T9*T9+1.7182625e8*T9*T9*T9-8.3103078e8*pow(T9,4.)+2.1243451e9*pow(T9,5.)-2.872313e9*pow(T9,6.)+2.0104043e9*pow(T9,7.)-4.3859588e8*pow(T9,8.)-3.529339e8*pow(T9,9.)+2.9815567e8*pow(T9,10.)-8.8920729e7*pow(T9,11.)+9.9850915e6*pow(T9,12.))*pow(T9,-2./3.))/f[Li7p_aa]-1.,2.));
			ferrlow[Li7p_aa]=-sqrt(0.08*0.08+pow((((-2.9979375e6-7.8110137e8*T9+1.1816185e9*T9*T9+5.700657e6*pow(T9,1./3.)+1.330785e8*pow(T9,2./3.)+1.733923e9*pow(T9,4./3.)-1.8284296e9*pow(T9,5./3.)-5.0363158e8*pow(T9,7./3.)+1.1026194e8*pow(T9,8./3.))*pow(T9,-2./3.))/exp(5.35732631*pow(T9,-1./3.))+exp(-1.0068557*T9*T9-5.2092464*pow(T9,-1./3.))*(-14997.544+665017.06*T9-1.0880148e7*T9*T9+1.1299875e8*T9*T9*T9-5.3097151e8*pow(T9,4.)+1.3288827e9*pow(T9,5.)-1.7652952e9*pow(T9,6.)+1.2196578e9*pow(T9,7.)-2.6871614e8*pow(T9,8.)-2.0119802e8*pow(T9,9.)+1.7032325e8*pow(T9,10.)-5.0416533e7*pow(T9,11.)+5.6188182e6*pow(T9,12.))*pow(T9,-2./3.))/f[Li7p_aa]-1.,2.));
		}
		else
		{
			ferrhigh[Li7p_aa]=sqrt(0.08*0.08+pow((1.57087e6+85552.6)/f[Li7p_aa]-1.,2.));
			ferrlow[Li7p_aa]=-sqrt(0.08*0.08+pow((1.57087e6+85552.6)/f[Li7p_aa]-1.,2.));
		}
		ferr[Li7p_aa]=ferrlow[Li7p_aa];
	}

/* H2 + a -> g + Li6 */
    f[da_gLi6]=1.482e1*pow(T9,-2./3.)*exp(-7.435/pow(T9,1./3.))*(1.+6.572*T9+7.6e-2*T9*T9+2.48e-2*T9*T9*T9)+8.28e1*pow(T9,-3./2.)*exp(-7.904/T9);
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[da_gLi6]=0.289+5.612*exp(-3.*T9)-2.63*exp(-2.*T9)+0.773*exp(-T9);
		ferrlow[da_gLi6]=-0.9813+0.355*pow(T9,1./2.)-0.0411*T9;
		ferr[da_gLi6]=ferrlow[da_gLi6];
	}

/* H3 + a -> g + Li7 */
	if (T9<2.5) 
        f[ta_gLi7]=((0.094614248-4.9273133*T9+99.358965*T9*T9-989.81236*T9*T9*T9+4368.45*pow(T9,4.)+931.93597*pow(T9,5.)-391.07855*pow(T9,6.)+159.23101*pow(T9,7.)-34.407594*pow(T9,8.)+3.3919004*pow(T9,9.)+0.017556217*pow(T9,10.)-0.036253427*pow(T9,11.)+0.0031118827*pow(T9,12.)-0.00008714468*pow(T9,13.))*pow(T9,-1./2.))/(exp(8.4e-7*T9)*pow(1.+1.78616593*T9,3.));
	else 
        f[ta_gLi7]=807.406;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<2.5)
		{
			ferrhigh[ta_gLi7]=sqrt(0.0871*0.0871+pow((((0.083877015-4.5408918*T9+96.316095*T9*T9-1016.5548*T9*T9*T9+4809.4834*pow(T9,4.)-168.10236*pow(T9,5.)+208.81839*pow(T9,6.)-64.618239*pow(T9,7.)+10.478926*pow(T9,8.)-0.41782376*pow(T9,9.)-0.06453532*pow(T9,10.)+0.004777625*pow(T9,11.)+0.00020027244*pow(T9,12.)-0.000017864206*pow(T9,13.))*pow(T9,-1./2.))/(exp(9.3e-7*T9)*pow(1.+1.60170507*T9,3.)))/f[ta_gLi7]-1.,2.));
			ferrlow[ta_gLi7]=-sqrt(0.0871*0.0871+pow((((0.066096606-3.5622862*T9+75.13824*T9*T9-788.24146*T9*T9*T9+3705.8889*pow(T9,4.)-106.98552*pow(T9,5.)+139.5561*pow(T9,6.)-7.8984539*pow(T9,7.)-1.6035703*pow(T9,8.)-0.17508886*pow(T9,9.)+0.046425912*pow(T9,10.)+0.0030233156*pow(T9,11.)-0.00081682606*pow(T9,12.)+0.000034545163*pow(T9,13.))*pow(T9,-1./2.))/(exp(0.01111331*T9)*pow(1.+1.63277688*T9,3.)))/f[ta_gLi7]-1.,2.));
		}
		else
		{
			ferrhigh[ta_gLi7]=sqrt(0.0871*0.0871+pow(916.513/f[ta_gLi7]-1.,2.));
			ferrlow[ta_gLi7]=-sqrt(0.0871*0.0871+pow(698.635/f[ta_gLi7]-1.,2.));
		}
		ferr[ta_gLi7]=ferrlow[ta_gLi7];
	}

/* He3 + a -> g + Be7 */
	if(T9<2.5) 
        f[ha_gBe7]=((0.000046165644-0.00046036111*T9-0.021600946*T9*T9+0.069627779*T9*T9*T9+7.346612*pow(T9,4.)-95.123199*pow(T9,5.)+391.13123*pow(T9,6.)-187.23717*pow(T9,7.)+86.111544*pow(T9,8.)-21.630169*pow(T9,9.)+3.6006922*pow(T9,10.)-0.34322836*pow(T9,11.)+0.018106742*pow(T9,12.)-0.00035681506*pow(T9,13.))*pow(T9,-1./2.))/(exp(0.48102949*T9)*pow(1.+1.17917554*T9,3.));
	else 
        f[ha_gBe7]=149.06;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<2.5)
		{
			ferrhigh[ha_gBe7]=sqrt(0.048*0.048+2.1*pow((((0.000050712746-0.00048202784*T9-0.023831596*T9*T9+0.056033679*T9*T9*T9+8.408972*pow(T9,4.)-106.22688*pow(T9,5.)+434.78964*pow(T9,6.)-238.48007*pow(T9,7.)+94.757251*pow(T9,8.)-23.705813*pow(T9,9.)+3.8007127*pow(T9,10.)-0.37029512*pow(T9,11.)+0.019933598*pow(T9,12.)-0.00045281691*pow(T9,13.))*pow(T9,-1./2.))/(exp(0.2282385*T9)*pow(1.+1.31654256*T9,3.)))/f[ha_gBe7]-1.,2.));
			ferrlow[ha_gBe7]=-sqrt(0.048*0.048+2.1*pow((((0.000049798665-0.00047801278*T9-0.023362423*T9*T9+0.061790614*T9*T9*T9+8.0589558*pow(T9,4.)-102.19595*pow(T9,5.)+418.23687*pow(T9,6.)-229.34858*pow(T9,7.)+92.638713*pow(T9,8.)-23.370595*pow(T9,9.)+3.7644261*pow(T9,10.)-0.36726621*pow(T9,11.)+0.019784483*pow(T9,12.)-0.00044951929*pow(T9,13.))*pow(T9,-1./2.))/(exp(0.25325444*T9)*pow(1.+1.28569931*T9,3.)))/f[ha_gBe7]-1.,2.));
		}
		else
		{
			ferrhigh[ha_gBe7]=sqrt(0.048*0.048+2.1*pow(150.371/f[ha_gBe7]-1.,2.));
			ferrlow[ha_gBe7]=-sqrt(0.048*0.048+2.1*pow(147.698/f[ha_gBe7]-1.,2.));
		}
		ferr[ha_gBe7]=ferrlow[ha_gBe7];
	}

/* H2 + H2 -> n + He3 */
	if(T9<4.) 
        f[dd_nh]=((-1.8436156e6 - 6.1150115e7*T9 - 2.7251853e7*T9*T9 - 2.2800422e6*T9*T9*T9 - 252433.58*pow(T9,4.) - 284357.41*pow(T9,10./3.) + 906146.25*pow(T9,11./3.) + 1.2270083e7*pow(T9,1./3.) - 1.3680884e7*pow(T9,2./3.) + 1.328894e8*pow(T9,4./3.) - 1.1916242e7*pow(T9,5./3.) + 8.3705218e6*pow(T9,7./3.) + 2.2357751e6*pow(T9,8./3.))*pow(T9,-2./3.))/exp(1.*pow(T9,-1./3.));
	else 
        f[dd_nh]=4.98099e7;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<4.)
		{
			ferrhigh[dd_nh]=sqrt(0.01*0.01+13.8*pow((((-1.7643388e6-6.0708618e7*T9-2.9398403e7*T9*T9-1.6035942e6*T9*T9*T9-188921.3*pow(T9,4.)-345931.36*pow(T9,10./3.)+684862.04*pow(T9,11./3.)+1.1654271e7*pow(T9,1./3.)-1.2269338e7*pow(T9,2./3.)+1.26615e8*pow(T9,4./3.)-3.5412041e6*pow(T9,5./3.)+6.1432859e6*pow(T9,7./3.)+2.8526453e6*pow(T9,8./3.))*pow(T9,-2./3.))/exp(pow(T9,-1./3.)))/f[dd_nh]-1.,2.));
			ferrlow[dd_nh]=-sqrt(0.01*0.01+13.8*pow((((-1.7643388e6-6.0708618e7*T9-2.9398403e7*T9*T9-1.6035942e6*T9*T9*T9-188921.3*pow(T9,4.)-345931.36*pow(T9,10./3.)+684862.04*pow(T9,11./3.)+1.1654271e7*pow(T9,1./3.)-1.2269338e7*pow(T9,2./3.)+1.26615e8*pow(T9,4./3.)-3.5412041e6*pow(T9,5./3.)+6.1432859e6*pow(T9,7./3.)+2.8526453e6*pow(T9,8./3.))*pow(T9,-2./3.))/exp(pow(T9,-1./3.)))/f[dd_nh]-1.,2.));
		}
		else
		{
			ferrhigh[dd_nh]=sqrt(0.01*0.01+13.8*pow(4.99298e7/f[dd_nh]-1.,2.));
			ferrlow[dd_nh]=-sqrt(0.01*0.01+13.8*pow(4.969e7/f[dd_nh]-1.,2.));
		}
		ferr[dd_nh]=ferrlow[dd_nh];
	}
	
/* H2 + H2 -> p + H3 */
	if(T9<4.) 
        f[dd_pt]=((-5.8523126e6+2.3222535e8*T9-9.877862e6*T9*T9+5.2331507e7*pow(T9,1./3.)-1.7022642e8*pow(T9,2./3.)-1.1875268e8*pow(T9,4./3.)+5.2922232e7*pow(T9,5./3.))*pow(T9,-2./3.))/exp(1.0676573*pow(T9,-1./3.));
	else 
        f[dd_pt]=4.021e7;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<4.)
		{
			ferrhigh[dd_pt]=sqrt(0.01*0.01+12.3*pow((((-5.7455947e6+2.2299893e8*T9-9.5242666e6*T9*T9+5.1106128e7*pow(T9,1./3.)-1.651035e8*pow(T9,2./3.)-1.1215042e8*pow(T9,4./3.)+5.0522037e7*pow(T9,5./3.))*pow(T9,-2./3.))/exp(1.04599646*pow(T9,-1./3.)))/f[dd_pt]-1.,2.));
			ferrlow[dd_pt]=-sqrt(0.01*0.01+12.3*pow((((-5.8975937e6+2.3572024e8*T9-1.0002957e7*T9*T9+5.2820931e7*pow(T9,1./3.)-1.7219803e8*pow(T9,2./3.)-1.2126621e8*pow(T9,4./3.)+5.3809211e7*pow(T9,5./3.))*pow(T9,-2./3.))/exp(1.07535853*pow(T9,-1./3.)))/f[dd_pt]-1.,2.));
		}
		else
		{
			ferrhigh[dd_pt]=sqrt(0.01*0.01+12.3*pow(4.02597e7/f[dd_pt]-1.,2.));
			ferrlow[dd_pt]=-sqrt(0.01*0.01+12.3*pow(4.01609e7/f[dd_pt]-1.,2.));
		}
		ferr[dd_pt]=ferrlow[dd_pt];
	}

/* H3 + H2 -> n + He4 */
	if(T9<2.5) f[td_na]=6.2265733e8/(exp(0.49711597/T9)*pow(T9,0.56785403)) + exp(-0.23309803*T9*T9 - 1.342742*pow(T9,-1./3.))*(-8.1144927e7 + 2.2315324e9*T9 - 2.9439669e9*T9*T9 + 1.8764462e9*T9*T9 - 6.0511612e8*pow(T9,4.) + 9.5196576e7*pow(T9,5.) - 5.2901086e6*pow(T9,6.))*pow(T9,-2./3.);
	else f[td_na]=3.40249e8;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<2.5)
		{
			ferrhigh[td_na]=sqrt(0.0126*0.0126+1.4*pow((6.200594e8/((exp(0.49495969/T9)*pow(T9,0.56078105))+exp(-0.23797125*T9*T9-1.3784792*pow(T9,-1./3.))*(-8.7018245e7+2.4114301e9*T9-3.2227206e9*T9*T9+2.0779852e9*T9*T9*T9-6.7739586e8*pow(T9,4.)+1.0762439e8*pow(T9,5.)-6.0348254e6*pow(T9,6.))*pow(T9,-2./3.)))/f[td_na]-1.,2.));
			ferrlow[td_na]=-sqrt(0.0126*0.0126+1.4*pow((6.3798186e8/((exp(0.49598246/T9)*pow(T9,0.58460934))+exp(-0.33273637*T9*T9-1.0508793*pow(T9,-1./3.))*(-4.0964097e7+1.064899e9*T9-7.152721e8*T9*T9-1.4155217e8*T9*T9*T9+3.9276243e8*pow(T9,4.)-1.5817375e8*pow(T9,5.)+2.128034e7*pow(T9,6.))*pow(T9,-2./3.)))/f[td_na]-1.,2.));
		}
		else
		{
			ferrhigh[td_na]=sqrt(0.0126*0.0126+1.4*pow(3.41798e8/f[td_na]-1.,2.));
			ferrlow[td_na]=-sqrt(0.0126*0.0126+1.4*pow(3.38424e8/f[td_na]-1.,2.));
		}
		ferr[td_na]=ferrhigh[td_na];
	}

/* He3 + H2 -> p + He4 */
	if(T9<2.5) 
        f[hd_pa]=3.1038385e8/(exp(1.6190981/T9)*pow(T9,0.12159455))+exp(-0.0062340825*T9*T9-1.4540617*pow(T9,-1./3.))*(-3.1335916e7-6.2051071e8*T9-1.8782248e9*T9*T9+6.5642773e8*T9*T9*T9+1.530887e8*pow(T9,4.)-4.9542138e8*pow(T9,10./3.)-1.770285e8*pow(T9,11./3.)+1.14185e8*pow(T9,1./3.)-2.516526e7*pow(T9,13./3.)+1.7500204e8*pow(T9,2./3.)-1.7513362e9*pow(T9,4./3.)+5.2792247e9*pow(T9,5./3.)-3.32382e9*pow(T9,7./3.)+2.0346284e9*pow(T9,8./3.))*pow(T9,-2./3.);
	else 
        f[hd_pa]=1.55167e8;
	if((err!=0)&&(err!=-10000))
	{
		if(T9<2.5)
		{
			ferrhigh[hd_pa]=sqrt(0.00299*0.00299+3.8*pow((2.7540767e8/((exp(1.7895761/T9)*pow(T9,0.42535964))+exp(-0.011584496*T9*T9-1.7647266*pow(T9,-1./3.))*(-4.0539244e7-6.8068775e8*T9+1.6755542e9*T9*T9+1.3327241e9*T9*T9*T9+2.5284074e8*pow(T9,4.)-8.0072489e8*pow(T9,10./3.)-3.2332801e8*pow(T9,11./3.)+1.3990258e8*pow(T9,1./3.)-4.0197501e7*pow(T9,13./3.)+2.4121225e8*pow(T9,2./3.)-2.3960064e9*pow(T9,4./3.)+5.3331297e9*pow(T9,5./3.)-7.7996883e9*pow(T9,7./3.)+3.3487409e9*pow(T9,8./3.))*pow(T9,-2./3.)))/f[hd_pa]-1.,2.));
			ferrlow[hd_pa]=-sqrt(0.00299*0.00299+3.8*pow((2.7552759e8/((exp(1.5970464/T9)*pow(T9,0.0070474065))+exp(-0.0067819916*T9*T9-2.0484693*pow(T9,-1./3.))*(-4.6389646e6-3.2264085e9*T9+8.3768817e10*T9*T9+4.6593422e10*T9*T9*T9-5.3027407e9*pow(T9,4.)-5.1730322e10*pow(T9,10./3.)+2.3630624e10*pow(T9,11./3.)-2.3467142e8*pow(T9,1./3.)+4.8023403e8*pow(T9,13./3.)+1.6138031e9*pow(T9,2./3.)+5.5023454e9*pow(T9,4./3.)-2.9668793e10*pow(T9,5./3.)-9.5677252e10*pow(T9,7./3.)+2.4498194e10*pow(T9,8./3.))*pow(T9,-2./3.)))/f[hd_pa]-1.,2.));	
		}
		else
		{
			ferrhigh[hd_pa]=sqrt(0.00299*0.00299+3.8*pow(1.55567e8/f[hd_pa]-1.,2.));
			ferrlow[hd_pa]=-sqrt(0.00299*0.00299+3.8*pow(1.54638e8/f[hd_pa]-1.,2.));
		}
		ferr[hd_pa]=ferrlow[hd_pa];
	}

/* He3 + He3 -> 2p + He4 */
	f[hh_ppa]=5.59e10*pow(T9,-2./3.)*exp(-12.277/pow(T9,1./3.))*(1.-0.135*T9+2.54e-2*T9*T9 - 1.29e-3*T9*T9*T9);
	if((err!=0)&&(err!=-10000))
	{
		ferrlow[hh_ppa]=-1.e-2-7.34e-2-0.102*exp(-T9)+3.23e-4*T9;
		ferrhigh[hh_ppa]=-ferrlow[hh_ppa];
		ferr[hh_ppa]=ferrlow[hh_ppa];
	}

/* Li7 + H2 -> n + a + He4 */
	if (T9<2.5) 
        f[Li7d_naa]=1.66e11*pow(T9,-2./3.)*exp(-10.254/pow(T9,1./3.))+ 1.71e6*pow(T9,-3./2.)*exp(-3.246/T9)  + 1.49e10*pow(T9,-3./2.)*exp(-4.0894/T9)* (2.57e-2/T9 + 2.6314*pow(T9,-2./3.)-4.1929/pow(T9,1./3.) -2.1241+4.1136*pow(T9,1./3.));
	else 
        f[Li7d_naa]=1.37518e9;
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Li7d_naa]=0.5;
		ferrlow[Li7d_naa]=-0.5;
		ferr[Li7d_naa]=ferrlow[Li7d_naa];
	}

/* Be7 + H2 -> p + a + He4 */
	f[Be7d_paa]=pow(T9,-2./3.)*1.07e12*exp(-12.428/pow(T9,1./3.));
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Be7d_paa]=0.9;
		ferrlow[Be7d_paa]=-0.9;
		ferr[Be7d_paa]=ferrlow[Be7d_paa];
	}

/* Li7 + n -> g + Li8 */
	f[Li7n_gLi8]=3.144e3+4.26e3*pow(T9,-3./2.)*exp(-2.576/T9);
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Li7n_gLi8]=0.3;
		ferrlow[Li7n_gLi8]=-0.3;
		ferr[Li7n_gLi8]=ferrlow[Li7n_gLi8];
	}

/* B10 + n -> g + B11 */
	f[B10n_gB11]=66200.;

/* B11 + n -> g + B12 */
	f[B11n_gB12]=pow(T9,-1.5)*2400.*exp(-0.223/T9)+729.;

/* C11 + n -> p + B11 */
	f[C11n_pB11]=(1.-sqrt(T9)*0.048+T9*0.01)*1.69e8;

/* B10 + n -> a + Li7 */
	f[B10n_aLi7]=2.2e7*(1.+1.064*T9);
	if((err!=0)&&(err!=-10000))
	{
		ferrlow[B10n_aLi7]=-0.10794+2.2003e-2*T9-3.51236e-3*T9*T9;
		ferrhigh[B10n_aLi7]=-ferrlow[B10n_aLi7];
		ferr[B10n_aLi7]=ferrlow[B10n_aLi7];
	}
	
/* Be7 + p -> g + B8 */
	f[Be7p_gB8]=2.61e5*pow(T9,-2./3.)*exp(-10.264/pow(T9,1./3.))*(1.-5.11e-2*T9+ 4.68e-2*T9*T9-6.6e-3*T9*T9*T9+3.12e-4*T9*T9*T9*T9)+ 2.05e3*pow(T9,-3./2.)*exp(-7.345/T9);
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Be7p_gB8]=0.15;
		ferrlow[Be7p_gB8]=-0.15;
		ferr[Be7p_gB8]=ferrlow[Be7p_gB8];
	}

/* Be9 + p -> g + B10 */
	f[Be9p_gB10]=pow(T9,-2./3.)*1.33e7*exp(-10.359/pow(T9,1./3.)-T9/0.846*T9/0.846)*(pow(T9,1./3.)*0.04+1.+pow(T9,2./3.)*1.52+T9*0.428+pow(T9,4./3.)*2.15+pow(T9,5./3.)*1.54)+pow(T9,-1.5)*96400.*exp(-3.445/T9)+pow(T9,-1.5)*2.72e6*exp(-10.62/T9);

/* B10 + p -> g + C11 */
	f[B10p_gC11]=pow(T9,-2./3.)*1.68e6*exp(-pow(T9,-1./3.)*12.064)/(pow(pow(T9,2./3.)-0.0273,2.)+4.69e-4)*(1.+0.977*T9+1.87*T9*T9-0.272*T9*T9*T9+0.013*T9*T9*T9*T9); 

/* B11 + p -> g + C12 */
	f[B11p_gC12]=pow(T9,-2./3.)*4.58e7*exp(-pow(T9,-1./3.)*12.097-pow(T9/0.6,2.))*(1.+0.353*T9-0.842*T9*T9)+pow(T9,-3./2.)*6.82e3*exp(-1.738/T9)+2.8e4*pow(T9,0.104)*exp(-3.892/T9); 
	
/* C11 + p -> g + N12 */
	f[C11p_gN12]=pow(T9,-2./3.)*42400.*exp(-13.658/pow(T9,1./3.)-T9/1.627*T9/1.627)*(pow(T9,1./3.)*0.031+1.+pow(T9,2./3.)*3.11+T9*0.665+pow(T9,4./3.)*4.61+pow(T9,5./3.)*2.5)+pow(T9,-1.5)*8840.*exp(-7.021/T9);
	
/* B12 + p -> n + C12 */
	f[B12p_nC12]=pow(T9,-2./3.)*4.02e11*exp(-12.12/pow(T9,1./3.));
	ferr[45]=0.3;
	
/* Be9 + p -> a + Li6 */
	f[Be9p_aLi6]=2.11e11*pow(T9,-2./3.)*exp(-10.361/pow(T9,1./3.)-pow(T9/0.4,2.))*(1.-0.189*T9+3.52e1*T9*T9)+5.24e8*pow(T9,-3./2.)*exp(-3.446/T9)+4.65e8*exp(-4.396/T9)/(pow(T9,0.293));
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Be9p_aLi6]=9.849e-2+0.9084*exp(-20.*T9)+1.21e-2*T9-5.987e-4*T9*T9;
		ferrlow[Be9p_aLi6]=-0.15236-0.676965*exp(-20.*T9)+0.13113*T9*exp(-T9);
		ferr[Be9p_aLi6]=ferrlow[Be9p_aLi6];
	}

/* B10 + p -> a + Be7 */
	if(T9>0.8) 
        f[B10p_aBe7]=1.01e10*pow(T9,-2./3.)*exp(-12.064/pow(T9,1./3.))*(-1.+15.8*T9-2.6*T9*T9+0.125*T9*T9*T9);
	else 
        f[B10p_aBe7]= 2.56e10*pow(T9,-2./3.)*exp(-12.064/pow(T9,1./3.))*(1.+5.95*T9+2.92e1*T9*T9-3.16e2*T9*T9*T9+9.14e2*pow(T9,4.)-1.085e3*pow(T9,5.)+4.65e2*pow(T9,6.))/(4.7e-4+pow(pow(T9,2./3.)-2.6e-2,2.));
	if((err!=0)&&(err!=-10000))
	{
		ferrlow[B10p_aBe7]=-0.101-0.1234*T9+1.418e-2*T9*T9-5.79e-4*T9*T9*T9;
		ferrhigh[B10p_aBe7]=-ferrlow[B10p_aBe7];
		ferr[B10p_aBe7]=ferrlow[B10p_aBe7];
	}

/* B12 + p -> a + Be9 */
	f[B12p_aBe9]=pow(T9,-2./3.)*2.01e11*exp(-12.12/pow(T9,1./3.));
	ferr[B12p_aBe9]=0.3;
	
/* Li6 + a -> g + B10 */
	f[Li6a_gB10]=pow(T9,-2./3.)*4.06e6*exp(-18.79/pow(T9,1./3.)-T9/1.326*T9/1.326)*(pow(T9,1./3.)*0.022+1.+pow(T9,2./3.)*1.54+T9*0.239+pow(T9,4./3.)*2.2+pow(T9,5./3.)*0.869)+pow(T9,-1.5)*1910.*exp(-3.484/T9)+10100./T9*exp(-7.269/T9);
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[Li6a_gB10]=1.;
		ferrlow[Li6a_gB10]=-0.5;
		ferr[Li6a_gB10]=ferrlow[Li6a_gB10];
	}

/* Li7 + a -> g + B11 */
	if(T9>1.21) f[50]=1.187e3*pow(T9,-3./2.)*exp(-2.959/T9)+ 7.945e3*(1.+0.1466*T9-1.273e-2*T9*T9)* exp(-4.922/T9)/(pow(T9,2.3e-2));
	else f[50] = 9.72e7*pow(T9,-2./3.)*exp(-19.163/pow(T9,1./3.)-pow(T9/0.4,2.))* (1.+2.84*T9-7.89*T9*T9)+ 3.35e2*pow(T9,-3./2.)*exp(-2.959/T9)+ 1.04e4*exp(-4.922/T9)/(pow(T9,2.3e-2));
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[50]=0.248+9.6575e-2*exp(-.5*pow(T9-6.,2.))-5.796e-2*exp(-5.*pow(T9-0.2,2.))+0.7442*exp(-20.*T9);
		ferrlow[50]=-0.2246+3.9114e-2*exp(-0.5*pow(T9-6.,2.))+3.16e-2*exp(-5.*pow(T9-0.2,2.))-0.50145*exp(-20.*T9);
		ferr[50]=ferrlow[50];
	}

/* Be7 + a -> g + C11 */
	if (T9>2.) f[51] = 1.41e3*exp(-3.015/T9)*pow(T9,0.636);
	else f[51] = 1.29e10*pow(T9,-2./3.)*exp(-23.214/pow(T9,1./3.)-pow(T9/0.8,2.))*(1.-6.47*T9+19.5*T9*T9-19.3*T9*T9*T9)+ 1.25e4*pow(T9,-3./2.)*exp(-6.498/T9)+1.44e5*pow(T9,-3./2.)* exp(-10.177/T9)+1.63e4*exp(-15.281/T9)*pow(T9,0.178);
	if((err!=0)&&(err!=-10000))
	{
		ferrhigh[51]=0.3298+0.127*T9-3.611e-2*T9*T9+5.3544e-3*T9*T9*T9-2.6134e-4*pow(T9,4.);
		ferrlow[51]=-0.35187+5.31e-2*T9-4.684e-2*T9*T9+6.771e-3*T9*T9*T9-3.021e-4*pow(T9,4.);
		ferr[51]=ferrlow[51];
	}
	
/* B8 + a -> p + C11 */
	f[52]=pow(T9,-2./3.)*1.08e15*exp(-27.36/pow(T9,1./3.));
	
/* Li8 + a -> n + B11 */
	f[53]=pow(T9/(T9/15.1+1.),5./6.)*8.62e13*pow(T9,-1.5)*exp(-19.461/pow(T9/(T9/15.1+1.),1./3.));
	
/* Be9 + a -> n + C12 */
	f[54]=pow(T9,-2./3.)*4.62e13*exp(-23.87/pow(T9,1./3.)-T9/0.049*T9/0.049)*(pow(T9,1./3.)*0.017+1.+pow(T9,2./3.)*8.57+T9*1.05+pow(T9,4./3.)*74.51+pow(T9,5./3.)*23.15)+pow(T9,-1.5)*7.34e-5*exp(-1.184/T9)+pow(T9,-1.5)*0.227*exp(-1.834/T9)+pow(T9,-1.5)*1.26e5*exp(-4.179/T9)+exp(-12.732/T9)*2.4e8;
	
/* Be9 + H2 -> n + B10 */
	f[55]=pow(T9,-2./3.)*7.16e8*exp(6.44-12.6/pow(T9,1./3.));
	
/* B10 + H2 -> p + B11 */
	f[56]=pow(T9,-2./3.)*9.53e8*exp(7.3-14.8/pow(T9,1./3.));

/* B11 + H2 -> n + C12 */
	f[57]=pow(T9,-2./3.)*1.41e9*exp(7.4-14.8/pow(T9,1./3.));

/* He4 + a + n -> g + Be9 */
	f[58]=2.59e-6/((T9*0.344+1.)*(T9*T9))*exp(-1.062/T9);
	
/* He4 + 2a -> g + C12 */
	f[59]=pow(T9,-1.5)*2.79e-8*pow(T9,-1.5)*exp(-4.4027/T9)+pow(T9,-1.5)*1.35e-8*exp(-24.811/T9);
	
/* Li8 + p -> n + a + He4 */
	f[60]=pow(T9,-2./3.)*8.65e9*exp(-8.52/pow(T9,1./3.)-T9/2.53*T9/2.53)+pow(T9,-1.5)*2.31e9*exp(-4.64/T9);
	
/* B8 + n -> p + a + He4 */
	f[61]=4.02e8;
	
/* Be9 + p -> d + a + He4 */
	f[62]=pow(T9,-2./3.)*2.11e11*exp(-10.359/pow(T9,1./3.)-T9/0.52*T9/0.52)*(pow(T9,1./3.)*0.04+1.+pow(T9,2./3.)*1.09+T9*0.307+pow(T9,4./3.)*3.21+pow(T9,5./3.)*2.3)+5.79e8/T9*exp(-3.046/T9)+pow(T9,-0.75)*8.5e8*exp(-5.8/T9);
	
/* B11 + p -> 2 a + He4 */
	f[63]=pow(T9,-2./3.)*2.2e12*exp(-12.095/pow(T9,1./3.)-T9/1.644*T9/1.644)*(pow(T9,1./3.)*0.034+1.+pow(T9,2./3.)*0.14+T9*0.034+pow(T9,4./3.)*0.19+pow(T9,5./3.)*0.116)+pow(T9,-1.5)*4.03e6*exp(-1.734/T9)+pow(T9,-1.5)*6.73e9*exp(-6.262/T9)+3.88e9/T9*exp(-14.154/T9);

/* C11 + n -> 2a + He4 */
	f[64]=1.58e8;

/* C12 + n -> g + C13 */
	f[65]=450.;

/* C13 + n -> g + C14 */
	f[66]=pow(T9,-1.5)*2.38e5*exp(-1.67/T9)+119.;

/* N14 + n -> g + N15 */
	f[67]=9940.;

/* N13 + n -> p + C13 */
	f[68]= 1.178e8*(1.+3.36e-1*T9-3.792e-2*T9*T9+2.02e-3*T9*T9*T9);
	f[68]*=(1.+1.131*exp(-1.2892e1/T9+1.9e-2*T9));
	ferr[68]=0.3;

/* N14 + n -> p + C14 */
	f[69]=(sqrt(T9)*0.361+1.+T9*0.502)*2.39e5+1.112e8/sqrt(T9)*exp(-4.983/T9);

/* O15 + n -> p + N15 */	
	f[70] = 1.158e8*(1.+2.19e-1*T9-2.9e-2*T9*T9+1.73e-3*T9*T9*T9);
	f[70]*=(1.+3.87e-1*exp(-26.171/T9+1.18e-1*T9));

/* O15 + n -> a + C12 */
	f[71]=(sqrt(T9)*0.188+1.+T9*0.015)*3.5e7;
	
/* C12 + p -> g + N13 */
	
f[72]=2.e7*pow(T9,-2./3.)*exp(-13.692*pow(T9,-1./3.)-pow(T9/0.46,2.))*(1.+9.89*T9-59.8*T9*T9+266.*T9*T9*T9)+1.e5*pow(T9,-3./2.)*exp(-4.913/T9)+4.24e5*pow(T9,-3./2.)*exp(-21.62/T9);
	
/* C13 + p -> g + N14 */
	f[73]=9.57e7*pow(T9,-2./3.)*exp(-13.720*pow(T9,-1./3.)-T9*T9)*(1.+3.56*T9)+1.5e6*pow(T9,-3./2.)*exp(-5.930/T9)+6.83e5*pow(T9,-8.64e-1)*exp(-12.057/T9);
	f[72]*=(1.-2.07*exp(-37.938/T9));
	
/* C14 + p -> g + N15 */
	f[74]=pow(T9,-2./3.)*6.8e6*exp(-13.741/pow(T9,1./3.)-T9/5.721*T9/5.721)*(pow(T9,1./3.)*0.03+1.+pow(T9,2./3.)*0.503+T9*0.107+pow(T9,4./3.)*0.213+pow(T9,5./3.)*0.115)+pow(T9,-1.5)*5360.*exp(-3.811/T9)+pow(T9,-1./3.)*98200.*exp(-4.739/T9);
	
/* N13 + p -> g + O14 */
	f[75]=pow(T9,-2./3.)*4.04e7*exp(-15.202/pow(T9,1./3.)-T9/1.191*T9/1.191)*(pow(T9,1./3.)*0.027+1.-pow(T9,2./3.)*0.803-T9*0.154+pow(T9,4./3.)*5.+pow(T9,5./3.)*2.44)+pow(T9,-1.5)*2.43e5*exp(-6.348/T9);
	
/* N14 + p -> g + O15 */
	f[76]=pow(T9,-2./3.)*4.9e7*exp(-15.228/pow(T9,1./3.)-T9/3.294*T9/3.294)*(pow(T9,1./3.)*0.027+1.-pow(T9,2./3.)*0.778-T9*0.149+pow(T9,4./3.)*0.261+pow(T9,5./3.)*0.127)+pow(T9,-1.5)*2370.*exp(-3.011/T9)+exp(-12.53/T9)*21900.;
	
/* N15 + p -> g + O16 */
	f[77]=pow(T9,-2./3.)*9.78e8*exp(-15.251/pow(T9,1./3.)-T9/0.45*T9/0.45)*(pow(T9,1./3.)*0.027+1.+pow(T9,2./3.)*0.219+T9*0.042+pow(T9,4./3.)*6.83+pow(T9,5./3.)*3.32)+pow(T9,-1.5)*11100.*exp(-3.328/T9)+pow(T9,-1.5)*14900.*exp(-4.665/T9)+pow(T9,-1.5)*3.8e6*exp(-11.048/T9);
	
/* N15 + p -> a + C12 */
	if (T9>2.5) f[78]= 4.17e7*pow(T9,0.917)*exp(-3.292/T9);
	else f[78]= 1.12e12*pow(T9,-2./3.)*exp(-15.253*pow(T9,-1./3.)-pow(T9/0.28,2.))*(1.+4.95*T9+143.*T9*T9)+1.01e8*pow(T9,-3./2.)* exp(-3.643/T9)+1.19e9*pow(T9,-3./2.)*exp(-7.406/T9);
	
/* C12 + a -> g + O16 */
	f[79]=1.04e8/(T9*T9)*exp(-32.12/pow(T9,1./3.)-T9/3.496*T9/3.496)/pow(pow(T9,-2./3.)*0.0489+1.,2.)+1.76e8/(T9*T9)/pow(pow(T9,-2./3.)*0.2654+1.,2.)*exp(-32.12/pow(T9,1./3.))+pow(T9,-1.5)*1250.*exp(-27.499/T9)+pow(T9,5.)*0.0143*exp(-15.541/T9);
	
/* B10 + a -> p + C13 */
	f[80]=pow(T9,-2./3.)*9.6e14*exp(-27.99/pow(T9,1./3.));
	
/* B11 + a -> p + C14 */
	f[81]=pow(T9,-2./3.)*5.37e11*exp(-28.234/pow(T9,1./3.)-T9/0.347*T9/0.347)*(pow(T9,1./3.)*0.015+1.+pow(T9,2./3.)*5.575+T9*0.576+pow(T9,4./3.)*15.888+pow(T9,5./3.)*4.174)+pow(T9,-1.5)*0.00544*exp(-2.827/T9)+pow(T9,-1.5)*336.*exp(-5.178/T9)+5.32e6/pow(T9,0.375)*exp(-11.617/T9);
	
/* C11 + a -> p + N14 */
	f[82]=pow(T9/(T9*0.0478+1.+pow(T9,5./3.)*0.00756/pow(T9*0.0478+1.,2./3.)),5./6.)*7.15e15*pow(T9,-1.5)*exp(-31.883/pow(T9/(T9*0.0478+1.+pow(T9,5./3.)*0.00756/pow(T9*0.0478+1.,2./3.)),1./3.));
	f[82]*= (1.+0.140*exp(-0.275/T9-0.210*T9));

/* N12 + a -> p + O15 */
	f[83]=pow(T9,-2./3.)*5.59e16*exp(-35.6/pow(T9,1./3.));
	
/* N13 + a -> p + O16 */
	f[84]=pow(T9/(T9*0.0776+1.+pow(T9,5./3.)*0.0264/pow(T9*0.0776+1.,2./3.)),5./6.)*3.23e17*pow(T9,-1.5)*exp(-35.829/pow(T9/(T9*0.0776+1.+pow(T9,5./3.)*0.0264/pow(T9*0.0776+1.,2./3.)),1./3.));
	
/* B10 + a -> n + N13 */
	f[85]=pow(T9,-2./3.)*1.2e13*exp(-27.989/pow(T9,1./3.)-T9/9.589*T9/9.589);

/* B11 + a -> n + N14 */
	f[86]=pow(T9,-2./3.)*6.97e12*exp(-28.234/pow(T9,1./3.)-T9/0.14*T9/0.14)*(pow(T9,1./3.)*0.015+1.+pow(T9,2./3.)*8.115+T9*0.838+pow(T9,4./3.)*39.804+pow(T9,5./3.)*10.456)+pow(T9,-1.5)*1.79*exp(-2.827/T9)+pow(T9,-1.5)*1710.*exp(-5.178/T9)+pow(T9,0.6)*4.49e6*exp(-8.596/T9);

/* B12 + a -> n + N15 */
	f[87]=pow(T9,-2./3.)*3.04e15*exp(-28.45/pow(T9,1./3.));

/* C13 + a -> n + O16 */
	if(T9>4.) 
        f[88]= 7.59e6*pow(T9,1.078)*exp(-12.056/T9);
	else 
        f[88]= 3.78e14/T9/T9*exp(-32.333*pow(T9,-1./3.)-pow(T9/0.71,2.))*(1.+4.68e1*T9-2.92e2*T9*T9+7.38e2*T9*T9*T9)+2.3e7*pow(T9,0.45)*exp(-13.03/T9);
	f[88]*=(1.+7.3318e1*exp(-58.176/T9-1.98e-1*T9));
	
	for(reac=first; reac <= last; reac++)
        f[reac] = fabs(f[reac]);
	
	if(err>0)
	{
		if(err>100000) srand((unsigned int)(getpid()+err));
		for(reac = first; reac <= last; reac++) 
		{
			if((fabs(1.+ferrhigh[reac])<1.)||(fabs(1.+ferrhigh[reac])>10.)) 
                ferrhigh[reac]=9.;
			if((fabs(1.+ferrlow[reac])>1.)||(fabs(1.+ferrlow[reac])<0.1)) 
                ferrlow[reac]=-0.9;
			if(err==1) 
                f[reac]*=fabs(1.+ferrhigh[reac]);
			if(err==2) 
                f[reac]*=fabs(1.+ferrlow[reac]);
			
			if(err>100000) f[reac]*=fabs(1.+ferr[reac]*rand_gauss());
		}
	}
	else if(err<=-12&&err>=-last) 
        f[-err]*=fabs(1.+ferr[-err]);
	
	return;
}
