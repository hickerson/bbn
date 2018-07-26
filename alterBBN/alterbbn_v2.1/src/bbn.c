#include "include.h"

//#define DEBUG
	
/*----------------------------------------------------*/

int linearize(double T, double reacparam[][10], double f[], double r[],
    int loop, int inc, int ip, double dt, double Y0[], double Y[],
    double dY_dt[], double H, double rhob)
/* solves for new abundances using gaussian elimination with
 * back substitution */
{
    /* Number of nuclides (#n1,#n2,#n3,#n4,#5,#6) for each of the 12 reaction types */
    double nn1[12]={1.,1.,1.,1.,1.,2.,3.,2.,1.,1.,2.,1.};
    double nn2[12]={0.,1.,1.,0.,1.,0.,0.,1.,1.,1.,0.,1.};
    double nn3[12]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double nn4[12]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.};
    double nn5[12]={0.,0.,1.,0.,0.,1.,0.,0.,1.,0.,2.,1.};
    double nn6[12]={1.,1.,1.,2.,2.,1.,1.,1.,2.,3.,1.,1.};
	
    int i,j,g,h,k,l,n,i1,j1,ind;
    double cn1,cn2,cn3,cn4,cn5,cn6,rn1,rn2,rn3,rn4,rn5,rn6,yY[NNUC+1];
    cn1=cn2=cn3=cn4=cn5=cn6=0.;
    int fail;
    double bdln;
#ifdef DEBUG
    int ierror;
#endif
    int c0 = 0;
    int type[NNUCREAC+1],n1[NNUCREAC+1],n2[NNUCREAC+1],n3[NNUCREAC+1],
            n4[NNUCREAC+1],n5[NNUCREAC+1],n6[NNUCREAC+1];
    double rev[NNUCREAC+1],q9[NNUCREAC+1];
    double a[NNUC+1][NNUC+1],b[NNUC+1],yx[NNUC+1];
    int icnvm;
    double x[NNUC+1], a0[NNUC+1][NNUC+1], cx, sum, xdy, t;
    int nord,test;

    for (i=1;i<=NNUCREAC;i++)
    {
        type[i]=(int)reacparam[i][1];
        n1[i]=(int)reacparam[i][2];
        n2[i]=(int)reacparam[i][3];
        n3[i]=(int)reacparam[i][4];
        n4[i]=(int)reacparam[i][5];
        n5[i]=(int)reacparam[i][6];
        n6[i]=(int)reacparam[i][7];
        rev[i]=reacparam[i][8];
        q9[i]=reacparam[i][9];
    }
	
    for(i=1;i<=NNUC;i++) for(j=1;j<=NNUC;j++) a[i][j]=0.;

    for (n=1;n<=NNUCREAC;n++)
    {
        ind=type[n];
        i=n1[n];
        j=n2[n];
        g=n3[n];
        h=n4[n];
        k=n5[n];
        l=n6[n];
        if (i <= NNUC && l <= NNUC)
        {
            rn1=nn1[ind];
            rn2=nn2[ind];
            rn3=nn3[ind];
            rn4=nn4[ind];
            rn5=nn5[ind];
            rn6=nn6[ind];
			
            switch(ind)
            {
                case 0:	{ /* (1,0,0,0,0,1) type */
                    cn1=f[n];
                    cn2=0.;
                    cn3=0.;
                    cn4=0.;
                    cn5=0.;
                    cn6=r[n];
                    break;}

                case 1: { /* (1,1,0,0,0,1) type */
                    r[n]=rev[n]*0.987e10*pow(T,1.5)*exp(-q9[n]/T)*f[n];
                    f[n]=rhob*f[n];
                    cn1=Y[j]*f[n]/2.;
                    cn2=Y[i]*f[n]/2.;
                    cn3=0.;
                    cn4=0.;
                    cn5=0.;
                    cn6=r[n];
                    break;}

                case 2:	{ /* (1,1,0,0,1,1) type */
                    f[n]=rhob*f[n];
                    r[n]=rev[n]*exp(-q9[n]/T)*f[n];
                    cn1=Y[j]*f[n]/2.;
                    cn2=Y[i]*f[n]/2.;
                    cn3=0.;
                    cn4=0.;
                    cn5=Y[l]*r[n]/2.;
                    cn6=Y[k]*r[n]/2.;
                    break;}

                case 3:	{ /* (1,0,0,0,0,2) type */
                    cn1=f[n];
                    cn2=0.;
                    cn3=0.;
                    cn4=0.;
                    cn5=0.;
                    cn6=Y[l]*r[n]/2.;
                    break;}

                case 4:	{ /* (1,1,0,0,0,2) type */
                    f[n]=rhob*f[n];
                    r[n]=rev[n]*exp(-q9[n]/T)*f[n];
                    cn1=Y[j]*f[n]/2.;
                    cn2=Y[i]*f[n]/2.;
                    cn3=0.;
                    cn4=0.;
                    cn5=0.;
                    cn6=Y[l]*r[n]/2.;
                    break;}

                case 5:	{ /* (2,0,0,0,1,1) type */
                    f[n]=rhob*f[n];
                    r[n]=rev[n]*exp(-q9[n]/T)*f[n];
                    cn1=Y[i]*f[n]/2.;
                    cn2=0.;
                    cn3=0.;
                    cn4=0.;
                    cn5=Y[l]*r[n]/2.;
                    cn6=Y[k]*r[n]/2.;
                    break;}

                case 6:	{ /* (3,0,0,0,0,1) type */
                    r[n]=rev[n]*0.974e20*pow(T,1.5)*pow(T,1.5)*
                            exp(-q9[n]/T)*f[n];
                    f[n]=rhob*rhob*f[n];
                    cn1=Y[i]*Y[i]*f[n]/6.;
                    cn2=0.;
                    cn3=0.;
                    cn4=0.;
                    cn5=0.;
                    cn6=r[n];
                    break;}
		
                case 7:	{ /* (2,1,0,0,0,1) type */
                    r[n]=rev[n]*0.974e20*pow(T,1.5)*pow(T,1.5)*
                            exp(-q9[n]/T)*f[n];
                    f[n]=rhob*rhob*f[n];
                    cn1=Y[j]*Y[i]*f[n]/3.;
                    cn2=Y[i]*Y[i]*f[n]/6.;
                    cn3=0.;
                    cn4=0.;
                    cn5=0.;
                    cn6=r[n];
                    break;}

                case 8:	{ /* (1,1,0,0,1,2) type */
                    f[n]=rhob*f[n];
                    r[n]=rev[n]*1.013e-10*pow(T,-1.5)*rhob*exp(-q9[n]/T)*f[n];
                    cn1=Y[j]*f[n]/2.;
                    cn2=Y[i]*f[n]/2.;
                    cn3=0.;
                    cn4=0.;
                    cn5=Y[l]*Y[l]*r[n]/6.;
                    cn6=Y[k]*Y[l]*r[n]/3.;
                    break;}

                case 9:	{ /* (1,1,0,0,0,3) type */
                    f[n]=rhob*f[n];
                    r[n]=rev[n]*1.013e-10*pow(T,-1.5)*rhob*exp(-q9[n]/T)*f[n];
                    cn1=Y[j]*f[n]/2.;
                    cn2=Y[i]*f[n]/2.;
                    cn3=0.;
                    cn4=0.;
                    cn5=0.;
                    cn6=Y[l]*Y[l]*r[n]/6.;
                    break;}

                case 10:{ /* (2,0,0,0,2,1) type */
                    f[n]=rhob*f[n];
                    r[n]=rev[n]*1.013e-10*pow(T,-1.5)*rhob*exp(-q9[n]/T)*f[n];
                    cn1=Y[i]*f[n]/2.;
                    cn2=0.;
                    cn3=0.;
                    cn4=0.;
                    cn5=Y[l]*Y[k]*r[n]/3.;
                    cn6=Y[k]*Y[k]*r[n]/6.;
                    break;}
                    
                case 11:{ /* (1,1,0,1,1,1) type */
                    f[n]=rhob*f[n];
                    r[n]=rev[n]*1.013e-10*pow(T,-1.5)*rhob*exp(-q9[n]/T)*f[n];
                    cn1=Y[j]*f[n]/2.;
                    cn2=Y[i]*f[n]/2.;
                    cn3=0.;
                    cn4=Y[k]*Y[l]*r[n]/3.;
                    cn5=Y[h]*Y[l]*r[n]/3.;
                    cn6=Y[h]*Y[k]*r[n]/3.; }
            }

            // Invert indexes
            i=NNUC+1-i;
            j=NNUC+1-j;
            g=NNUC+1-g;
            h=NNUC+1-h;
            k=NNUC+1-k;
            l=NNUC+1-l;

            // Fill i (n1) nuclide column
            a[i][i]+=rn1*cn1;
            if(j<=NNUC) a[j][i]+=rn2*cn1;
            if(g<=NNUC) a[g][i]+=rn3*cn1;
            if(h<=NNUC) a[h][i]-=rn4*cn1;
            if(k<=NNUC) a[k][i]-=rn5*cn1;
            a[l][i]-=rn6*cn1;
			
            // Fill j (n2) nuclide column
            if (j<=NNUC)
            {
                a[i][j]+=rn1*cn2;
                a[j][j]+=rn2*cn2;
                if(g<=NNUC) a[g][j]+=rn3*cn2;
                if(h<=NNUC) a[h][j]-=rn4*cn2;
                if(k<=NNUC) a[k][j]-=rn5*cn2;
                a[l][j]-=rn6*cn2;
            }
			
            // Fill g (n3) nuclide column
            if (g<=NNUC)
            {
                a[i][g]+=rn1*cn3;
                if(j<=NNUC) a[j][g]+=rn2*cn3;
                a[g][g]+=rn3*cn3;
                if(h<=NNUC) a[h][g]-=rn4*cn3;
                if(k<=NNUC) a[k][g]-=rn5*cn3;
                a[l][g]-=rn6*cn3;
            }

            // Fill h (n4) nuclide column
            if (h<=NNUC)
            {
                a[i][h]-=rn1*cn4;
                if(j<=NNUC) a[j][h]-=rn2*cn4;
                if(g<=NNUC) a[g][h]-=rn3*cn4;
                a[h][h]+=rn4*cn4;
                if(k<=NNUC) a[k][h]+=rn5*cn4;
                a[l][h]+=rn6*cn4;
            }

            // Fill k (n5) nuclide column
            if (k<=NNUC)
            {
                a[i][k]-=rn1*cn5;
                if(j<=NNUC) a[j][k]-=rn2*cn5;
                if(g<=NNUC) a[g][k]-=rn3*cn5;
                if(h<=NNUC) a[h][k]+=rn4*cn5;
                a[k][k]+=rn5*cn5;
                a[l][k]+=rn6*cn5;
            }

            // Fill l (n6) nuclide column
            a[i][l]-=rn1*cn6;
            if(j<=NNUC) a[j][l]-=rn2*cn6;
            if(g<=NNUC) a[g][l]-=rn3*cn6;
            if(h<=NNUC) a[h][l]+=rn4*cn6;
            if(k<=NNUC) a[k][l]+=rn5*cn6;
            a[l][l]+=rn6*cn6;

        }
    }

    // Finish the A matrix
    bdln=H*3.*1.e-5;
	
    for(i=1;i<=NNUC;i++)
    {
        i1=NNUC+1-i;                                            // Invert rows
        for(j=1;j<=NNUC;j++)
        {
            j1=NNUC+1-j;                                        // Invert columns
            if(fabs(a[j][i]) < bdln*Y0[j1]/Y0[i1]) a[j][i]=0.;  // Set 0 if tiny
            else a[j][i]*=dt;                                   // Bring dt over to the other side
        }
        a[i][i]+=1.;                                            // Add identity matrix
        b[i1]=Y0[i];                                            // Initial abundances
    }

    if(loop==1) icnvm=ip; else icnvm=c0;
	
    nord=0;
    fail=0;
    // Set RH and solution vectors to initial values
    for(i=1;i<=NNUC;i++)
    {
        x[i]=b[i];
        yx[i]=0.;
    }
    // Save matrix
    if(icnvm==inc) for(i=1;i<=NNUC;i++) for(j=1;j<=NNUC;j++) a0[j][i]=a[j][i];
    // If zeros at pivot points, terminate matrix evaluation
    for(i=1;i<=NNUC;i++)
    {
        if(a[i][i]==0.)
        {
            fail=i;
            return fail;
        }
        // Triangularize matrix
        for(j=i+1;j<=NNUC;j++)
        {
            if(a[j][i]!=0.)
            {
                cx=a[j][i]/a[i][i];
                for(k=i+1;k<=NNUC;k++) a[j][k]-=cx*a[i][k];
                a[j][i]=cx;
                x[j]-=cx*x[i];
            }
        }
    }

    // Back substitution
    do
    {	x[NNUC]/=a[NNUC][NNUC];
        yx[NNUC]+=x[NNUC];
		
        for(i=NNUC-1;i>=1;i--)
        {
            sum=0.;
            for(j=i+1;j<=NNUC;j++) sum+=a[i][j]*x[j];
            x[i]=(x[i]-sum)/a[i][i];
            yx[i]+=x[i];
        }

        test=1;
	
        if(icnvm==inc)
        {
            for(i=1;i<=NNUC;i++)
            {
                if(yx[i]!=0.)
                {
                    xdy=fabs(x[i]/yx[i]);
                    if(xdy>2.e-4)
                    {
                        if(nord<1)
                        {
                            nord++;
							
                            for(j=1;j<=NNUC;j++)
                            {
                                t = 0.;
                                for(k=1;k<=NNUC;k++) t+=a0[j][k]*yx[k];
                                x[j]=b[j]-t;
                            }

                            for(j=2;j<=NNUC;j++)
                                for(k=j+1;k<=NNUC;k++)
                                    x[k]-=a[k][j]*x[j];
                            break;
                        }
                        else
                        {
                            fail=-1;
#ifdef DEBUG
                            ierror=i;
#endif
                            return fail;
                        }
                    }
                    else test=0;
                }
                else test=0;
            }
        }
        else test=0;
    }
    while(test);
    // Derivatives of abundances
    for(i=1;i<=NNUC;i++)
    {
        yY[i]=yx[NNUC+1-i];
        dY_dt[i]=(yY[i]-Y0[i])/dt;
    }

#ifdef DEBUG
    if(fail!=0)
    {
        if(fail==-1) printf("y(%d) failed to converge\n",ierror);
        if(fail>=1) printf("%d th diagonal term equals zero\n",fail);
    }
#endif
    return fail;
}

/*----------------------------------------------------*/

int nucl_single(struct relicparam* paramrelic, double ratioH[NNUC+1], struct errorparam* paramerror)
/* Main routine with computes the abundance ratios H2_H, ..., Be7_H as well as
 * the baryon-to-photon ratio eta, using the parameters contained in paramrelic->
 * The err parameter is a switch to choose if the central (err=0), high (err=1)
 * or low (err=2) values of the nuclear rates is used. If (err=3),
 * the lower value of only the nuclear rate number " errnumber" is used. If (err=4),
 * the value of the nuclear rates is taken (gaussianly) randomly for a
 * MC analysis. */
{
	if(paramrelic->err==3) if((paramerror->errnumber<0)||(paramerror->errnumber>NNUCREAC+1)) return 1;			

	if((paramrelic->err==3)&&(paramerror->errnumber==NNUCREAC+1)) paramerror->life_neutron=paramrelic->life_neutron+paramrelic->life_neutron_error;	
	else if(paramrelic->err==4) paramerror->life_neutron=paramrelic->life_neutron+paramrelic->life_neutron_error*paramerror->random[NNUCREAC+1];
	else paramerror->life_neutron=paramrelic->life_neutron;
					
    int i;
    for(i=1;i<=NNUC;i++) ratioH[i]=0.;
    double f[NNUCREAC+1],r[NNUCREAC+1];
    double sd;
    double rhod,Pd,rhodprev,Sigmarad,sum_Y;
    double drhod_dT=0;
    double sum_dY_dt, sum_ZY, dsd_dT, dphie_dT, dlna3_dT, dphie_dlna3,
            dphie_dZY, sum_DeltaMdY_dt, sum_ZdY_dt;
    double cosh1, cosh2, cosh3, cosh4, cosh5, cosh6, cosh7, sinh1, sinh2, sinh3,
            sinh4, sinh5, sinh6, sinh7;
    double cosh1W, cosh2W, cosh3W, cosh4W, cosh5W, cosh6W, cosh7W;
    double sinh1W, sinh2W, sinh3W, sinh4W, sinh5W, sinh6W, sinh7W;
    double T0,h_eta0,phie0,phiW0;
    double dtl;
    int loop;
    double dh_dt, dphie_dt, dT_dt, dlnT_dt, dphiW_dt;
    double dT0_dt, dh_dt0, dphie_dt0,dphiW0_dt;
    double dlna3_dTnu,dTnu_dt,dlnTnu_dt,Tnu0,dTnu0_dt;
    double dY_dt0[NNUC+1],dY_dt[NNUC+1],Y0[NNUC+1],Y[NNUC+1];
    double dtmin;
    double z;               // Parameterized electron mass -> z = m_e*c^2 / (k_B*T)
    double H;               // Hubble parameter
    double zW;              // Parameterized WIMP mass -> zW = m_WIMP*c^2 / (k_B*T)
    double wimp_mass_ratio; // Ratio between WIMP mass and electron mass
    double zWd;             // Value of zW at neutrino decoupling temperature Tnud
    double phiW;            // Parameterized WIMP chemical potential

    dTnu_dt=wimp_mass_ratio=zW=dlna3_dT=dphie_dt0=dh_dt0=dT0_dt=phie0=h_eta0=T0=Tnu0=dTnu0_dt=phiW0=dphiW0_dt=0.;

    /* Nuclides: 1=n, 2=p, 3=H2, 4=H3, 5=He3, 6=He4, 7=Li6, 8=Li7, 9=Be7,
     * 10=Li8, 11=B8, 12=Be9, 13=B10, 14=B11, 15=C11, 16=B12, 17=C12, 18=N12,
     * 19=C13, 20=N13, 21=C14, 22=N14, 23=O14, 24=N15, 25=O15, 26=O16 */

    double Am[NNUC+1] = {0., 1., 1., 2., 3., 3., 4., 6., 7., 7., 8., 8., 9.,
                         10., 11., 11., 12., 12., 12., 13., 13., 14., 14., 14.,
                         15., 15., 16.}; /* Atomic number A */

    double Zm[NNUC+1] = {0., 0., 1., 1., 1., 2., 2., 3., 3., 4., 3., 5., 4.,
                         5., 5., 6., 5., 6., 7., 6., 7., 6., 7., 8., 7., 8.,
                         8.}; /* Charge number Z */

    double Dm[NNUC+1] = {0., 8.071388, 7.289028, 13.135825, 14.949915,
                         14.931325, 2.424931, 14.0864, 14.9078, 15.7696,
                         20.9464, 22.9212, 11.34758, 12.05086, 8.6680, 10.6506,
                         13.3690, 0., 17.3382, 3.125036, 5.3455, 3.019916,
                         2.863440, 8.006521, 0.101439, 2.8554,
                         -4.737036}; /* mass excess DeltaM in MeV */

    double reacparam[NNUCREAC+1][10] =
    {

// type: 0-10, each type has a unique (#n1,#n2,#n3,#n4) quartet
// n1: incoming nuclide number
// n2: incoming light nuclide number
// n3: incoming lightest nuclide number
// n4: outgoing lightest nuclide number
// n5: outgoing light nuclide number
// n6: outgoing nuclide number
// rev: reverse reaction coefficient
// q: energy release in reaction in 10**9 Kelvin (K = ev/k_B)

//   reac# type n1 n2 n3 n4 n5 n6 rev q

    {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},              // none
    {1.,0.,1.,0.,0.,0.,0.,2.,0.,0.},              // n -> p
    {2.,0.,4.,0.,0.,0.,0.,5.,0.,0.},              // H3 -> e- + v + He3
    {3.,3.,10.,0.,0.,0.,0.,6.,0.,0.},             // Li8 -> e- + v + 2He4
    {4.,0.,16.,0.,0.,0.,0.,17.,0.,0.},            // B12 -> e- + v + C12
    {5.,0.,21.,0.,0.,0.,0.,22.,0.,0.},            // C14 -> e- + v + N14
    {6.,3.,11.,0.,0.,0.,0.,6.,0.,0.},             // B8 -> e+ + v + 2He4
    {7.,0.,15.,0.,0.,0.,0.,14.,0.,0.},            // C11 -> e+ + v + B11
    {8.,0.,18.,0.,0.,0.,0.,17.,0.,0.},            // N12 -> e+ + v + C12
    {9.,0.,20.,0.,0.,0.,0.,19.,0.,0.},            // N13 -> e+ + v + C13
    {10.,0.,23.,0.,0.,0.,0.,22.,0.,0.},           // O14 -> e+ + v + N14
    {11.,0.,25.,0.,0.,0.,0.,24.,0.,0.},           // O15 -> e+ + v + N15
    {12.,1.,2.,1.,0.,0.,0.,3.,0.477,25.815},       // H + n -> g + H2
    {13.,1.,3.,1.,0.,0.,0.,4.,1.65,72.612},        // H2 + n -> g + H3
    {14.,1.,5.,1.,0.,0.,0.,6.,2.63,238.794},       // He3 + n -> g + He4
    {15.,1.,7.,1.,0.,0.,0.,8.,1.20,84.132},        // Li6 + n -> g + Li7
    {16.,2.,5.,1.,0.,0.,2.,4.,1.001,8.863},       // He3 + n -> p + H3
    {17.,2.,9.,1.,0.,0.,2.,8.,1.001,19.080},      // Be7 + n -> p + Li7
    {18.,2.,7.,1.,0.,0.,4.,6.,1.068,55.503},      // Li6 + n -> a + H3
    {19.,4.,9.,1.,0.,0.,0.,6.,4.68,220.382},       // Be7 + n -> a + He4
    {20.,1.,3.,2.,0.,0.,0.,5.,1.65,63.749},       // H2 + p -> g + He3
    {21.,1.,4.,2.,0.,0.,0.,6.,2.63,229.931},      // H3 + p -> g + He4
    {22.,1.,7.,2.,0.,0.,0.,9.,1.20,65.053},       // Li6 + p -> g + Be7
    {23.,2.,7.,2.,0.,0.,5.,6.,1.067,46.640},       // Li6 + p -> a + He3
    {24.,4.,8.,2.,0.,0.,0.,6.,4.68,201.302},      // Li7 + p -> a + He4
    {25.,1.,6.,3.,0.,0.,0.,7.,1.55,17.109},       // H2 + a -> g + Li6
    {26.,1.,6.,4.,0.,0.,0.,8.,1.13,28.629},       // H3 + a -> g + Li7
    {27.,1.,6.,5.,0.,0.,0.,9.,1.13,18.412},       // He3 + a -> g + Be7
    {28.,5.,3.,0.,0.,0.,1.,5.,1.73,37.934},       // H2 + d -> n + He3
    {29.,5.,3.,0.,0.,0.,2.,4.,1.73,46.798},       // H2 + d -> p + H3
    {30.,2.,4.,3.,0.,0.,1.,6.,5.51,204.116},      // H3 + d -> n + He4
    {31.,2.,5.,3.,0.,0.,2.,6.,5.51,212.979},      // He3 + d -> p + He4
    {32.,10.,5.,0.,0.,0.,2.,6.,3.35,149.229},     // He3 + He3 -> 2p + He4
    {33.,8.,8.,3.,0.,0.,1.,6.,9.81,175.487},      // Li7 + d -> n + a + He4
    {34.,8.,9.,3.,0.,0.,2.,6.,9.83,194.566},      // Be7 + d -> p + a + He4
    {35.,1.,5.,4.,0.,0.,0.,7.,2.47,183.290},      // He3 + H3 -> g + Li6
    {36.,2.,7.,3.,0.,0.,1.,9.,2.52,39.237},       // Li6 + d -> n + Be7
    {37.,2.,7.,3.,0.,0.,2.,8.,2.52,58.317},       // Li6 + d -> p + Li7
    {38.,2.,5.,4.,0.,0.,3.,6.,1.59,166.181},      // He3 + H3 -> d + He4
    {39.,10.,4.,4.,0.,0.,1.,6.,3.34,131.503},     // H3 + H3 -> 2n + He4
    {40.,11.,5.,4.,0.,1.,2.,6.,3.34,140.366},     // He3 + H3 -> n + p + He4
    {41.,2.,8.,4.,0.,0.,1.,12.,3.55,121.136},     // Li7 + H3 -> n + Be9
    {42.,2.,9.,4.,0.,0.,2.,12.,3.55,140.215},     // Be7 + H3 -> p + Be9
    {43.,2.,8.,5.,0.,0.,2.,12.,3.55,129.999},     // Li7 + He3 -> p + Be9
    {44.,1.,8.,1.,0.,0.,0.,10.,1.33,23.589},       // Li7 + n -> g + Li8
    {45.,1.,13.,1.,0.,0.,0.,14.,3.07,132.920},     // B10 + n -> g + B11
    {46.,1.,14.,1.,0.,0.,0.,16.,2.37,39.111},      // B11 + n -> g + B12
    {47.,2.,15.,1.,0.,0.,2.,14.,1.001,32.086},    // C11 + n -> p + B11
    {48.,2.,13.,1.,0.,0.,6.,8.,0.755,32.371},     // B10 + n -> a + Li7
    {49.,1.,9.,2.,0.,0.,0.,11.,1.32,1.595},       // Be7 + p -> g + B8
    {50.,1.,12.,2.,0.,0.,0.,13.,0.986,76.424},    // Be9 + p -> g + B10
    {51.,1.,13.,2.,0.,0.,0.,15.,3.07,100.834},    // B10 + p -> g + C11
    {52.,1.,14.,2.,0.,0.,0.,17.,7.10,185.173},    // B11 + p -> g + C12
    {53.,1.,15.,2.,0.,0.,0.,18.,2.37,6.979},      // C11 + p -> g + N12
    {54.,2.,16.,2.,0.,0.,1.,17.,3.00,146.061},     // B12 + p -> n + C12
    {55.,2.,12.,2.,0.,0.,6.,7.,0.618,24.663},     // Be9 + p -> a + Li6
    {56.,2.,13.,2.,0.,0.,6.,9.,0.754,13.291},     // B10 + p -> a + Be7
    {57.,2.,16.,2.,0.,0.,6.,12.,0.291,79.903},     // B12 + p -> a + Be9
    {58.,1.,7.,6.,0.,0.,0.,13.,1.60,51.761},      // Li6 + a -> g + B10
    {59.,1.,8.,6.,0.,0.,0.,14.,4.07,100.549},     // Li7 + a -> g + B11
    {60.,1.,9.,6.,0.,0.,0.,15.,4.07,87.543},      // Be7 + a -> g + C11
    {61.,2.,11.,6.,0.,0.,2.,15.,3.07,85.948},      // B8 + a -> p + C11
    {62.,2.,10.,6.,0.,0.,1.,14.,3.07,76.960},      // Li8 + a -> n + B11
    {63.,2.,12.,6.,0.,0.,1.,17.,10.28,66.158},     // Be9 + a -> n + C12
    {64.,2.,12.,3.,0.,0.,1.,13.,2.06,50.609},      // Be9 + d -> n + B10
    {65.,2.,13.,3.,0.,0.,2.,14.,6.42,107.105},     // B10 + d -> p + B11
    {66.,2.,14.,3.,0.,0.,1.,17.,14.85,159.357},     // B11 + d -> n + C12
    {67.,7.,6.,1.,0.,0.,0.,12.,0.600,18.262},     // He4 + a + n -> g + Be9
    {68.,6.,6.,0.,0.,0.,0.,17.,2.06,84.420},      // He4 + 2a -> g + C12
    {69.,8.,10.,2.,0.,0.,1.,6.,3.54,177.713},      // Li8 + p -> n + a + He4
    {70.,8.,11.,1.,0.,0.,2.,6.,3.55,218.787},      // B8 + n -> p + a + He4
    {71.,8.,12.,2.,0.,0.,3.,6.,0.796,7.554},      // Be9 + p -> d + a + He4
    {72.,9.,14.,2.,0.,0.,0.,6.,3.45,100.753},     // B11 + p -> 2a + He4
    {73.,9.,15.,1.,0.,0.,0.,6.,3.46,132.838},      // C11 + n -> 2a + He4
    {74.,1.,17.,1.,0.,0.,0.,19.,0.898,57.400},     // C12 + n -> g + C13
    {75.,1.,19.,1.,0.,0.,0.,21.,3.62,94.884},      // C13 + n -> g + C14
    {76.,1.,22.,1.,0.,0.,0.,24.,2.74,125.715},     // N14 + n -> g + N15
    {77.,2.,20.,1.,0.,0.,2.,19.,1.001,34.846},    // N13 + n -> p + C13
    {78.,2.,22.,1.,0.,0.,2.,21.,3.00,7.263},     // N14 + n -> p + C14
    {79.,2.,25.,1.,0.,0.,2.,24.,1.001,41.037},    // O15 + n -> p + N15
    {80.,2.,25.,1.,0.,0.,6.,17.,0.707,98.659},    // O15 + n -> a + C12
    {81.,1.,17.,2.,0.,0.,0.,20.,0.896,22.554},    // C12 + p -> g + N13
    {82.,1.,19.,2.,0.,0.,0.,22.,1.21,87.621},     // C13 + p -> g + N14
    {83.,1.,21.,2.,0.,0.,0.,24.,0.912,118.452},   // C14 + p -> g + N15
    {84.,1.,20.,2.,0.,0.,0.,23.,3.62,53.705},     // N13 + p -> g + O14
    {85.,1.,22.,2.,0.,0.,0.,25.,2.73,84.678},     // N14 + p -> g + O15
    {86.,2.,24.,2.,0.,0.,0.,26.,3.67,140.733},    // N15 + p -> g + O16
    {87.,2.,24.,2.,0.,0.,6.,17.,0.706,57.622},	  // N15 + p -> a + C12
    {88.,1.,17.,6.,0.,0.,0.,26.,5.20,83.111},     // C12 + a -> g + O16
    {89.,2.,13.,6.,0.,0.,2.,19.,9.35,47.134},      // B10 + a -> p + C13
    {90.,2.,14.,6.,0.,0.,2.,21.,11.03,9.098},      // B11 + a -> p + C14
    {91.,2.,15.,6.,0.,0.,2.,22.,3.68,33.921},     // C11 + a -> p + N14
    {92.,2.,18.,6.,0.,0.,2.,25.,4.25,111.620},     // N12 + a -> p + O15
    {93.,2.,20.,6.,0.,0.,2.,26.,5.80,60.557},     // N13 + a -> p + O16
    {94.,2.,13.,6.,0.,0.,1.,20.,9.34,12.288},     // B10 + a -> n + N13
    {95.,2.,14.,6.,0.,0.,1.,22.,3.67,1.835},      // B11 + a -> n + N14
    {96.,2.,16.,6.,0.,0.,1.,24.,4.25,88.439},      // B12 + a -> n + N15
    {97.,2.,19.,6.,0.,0.,1.,26.,5.79,25.711},     // C13 + a -> n + O16
    {98.,2.,14.,3.,0.,0.,2.,16.,4.96,13.296},     // B11 + d -> p + B12
    {99.,2.,17.,3.,0.,0.,2.,19.,1.88,31.585},     // C12 + d -> p + C13
    {100.,2.,19.,3.,0.,0.,2.,21.,7.58,69.069},    // C13 + d -> p + C14
    };


    for(i=0;i<=NNUCREAC;i++)
    {
        f[i] = 0.;
        r[i] = 0.;
    }

    double norm=1.;
    if((paramrelic->wimp)||(paramrelic->xinu1!=0.))
    {
        double f_tmp[2],r_tmp[2];
        int err=paramrelic->err;
        paramrelic->err=0;
        rate_pn(f_tmp,r_tmp,0.00001,0.00001,paramrelic,paramerror);
        paramrelic->err=err;
        norm=1./f_tmp[1]/paramerror->life_neutron;
    }
    /* Note that Neff0 and Neff are not used in the program after they are computed. However, since they are
     * observables in CMB measurements, they are important variables when considering the effect of light WIMPS. */
    double Neff0;                   // Part of Neff that is a function of the WIMP mass
    double Neff;                    // Effective number of neutrinos
    double Ti=paramrelic->Tinit*K_to_GeV;    // Initial temperature in GeV
    double Tf=0.01*K_to_GeV;                // Final temperature in GeV
    double Ytmin =1.e-30;
    int inc=50;
    int ltime=0;
    int is=1;
    int ip=inc;
    int it=0;
    int fail=0;
	double cy,ct,dt0;
	int nitmax;

    if(paramrelic->failsafe==0)
    {
        /* Conservative limiting iteration values for optimization of computational time.
         * If the changes in the abundances and/or temperature are too large, the program is redirected to
         * the method "failsafe", which does not try to optimize the computational time
         * (see method nucl). */
        cy=0.5; // Limiting value of dY/dt which, together with ct, controls the step-size
        ct=0.1; // Limiting value of dT/dt
        dt0=1.e-2*s_to_GeV;
        nitmax=10;
    }
    else
    {
		cy=pow(10.,-paramrelic->failsafe);                  // Limiting value of dY/dt which, together with ct, controls the step-size
		ct=pow(10.,-1-paramrelic->failsafe);                 // Limiting value of dT/dt
		dt0=1.e-10*s_to_GeV;
		nitmax=pow(10,1+paramrelic->failsafe);
	}

	if(paramrelic->phi_model&&paramrelic->rho_phi==0.) paramrelic->rho_phi=(paramrelic->rhot_phi0)*pow(pi,2.)/15.*pow(Ti,4.);
	

    double dt=dt0;

    /* Initialization of relevant temperatures.
     * T: photon/e+- temp. Also the WIMP temperature in the case of EM coupled WIMPs
     * Tnu: SM neutrino temp. Simply redshifted in the case of no WIMPs or EM coupled WIMPs.
     * Tnud: SM neutrino decoupling temp. Instantaneous decoupling at ~2 MeV assumed.
     * Tnu_eq: Temp. of the equivalent neutrinos. Same as Tnu except in the case of WIMPs that couple only to SM neutrinos.
     * Note that we do not need to separately account for the WIMP temp. since this is equal to either T or Tnu, depending
     * on their coupling to the SM particles.
     */
    double T=Ti;
    double Tprev=T;
    double Tnu=T;
    double Tnud=Ti;   // Neutrino decoupling temperature
    double Tnu_eq=Tnu; // Temperature of the equivalent neutrinos. Same as Tnu except in the case of WIMPs with coupling=1

    if (DMpn / T > 58.)
    {
        Y[1] = 1.e-25;
        Y[2] = 1.;
    }
    else if (DMpn / T < -58.)
    {
        Y[1] = 1.;
        Y[2] = 1.e-25;
    }
    else
    {
        Y[1] = 1. / (exp(DMpn / T) + 1.);
        Y[2] = 1. / (exp(-DMpn / T) + 1.);
    }

    Y0[1]=Y[1];
    Y0[2]=Y[2];


    z=m_e/T;
    if (paramrelic->wimp)
    {
        wimp_mass_ratio = paramrelic->m_chi / 1.e+3 / m_e;
        zW = wimp_mass_ratio*z;                         // Always, since T=Tnu initially
        zWd = wimp_mass_ratio*z*T/Tnud;                // In the case that the iteration starts earlier than Tnud
    }

    double rho_gamma,P_gamma,drho_gamma_dT;
    double rho_epem,P_epem,drho_epem_dT,drho_epem_dphie,dM_epem_dT,dN_epem_dphie;
    double rho_neutrinos,P_neutrinos,drho_neutrinos_dTnu,rho_neuteq,P_neuteq,drho_neuteq;
    double rho_baryons=0.;
    double rho_wimp,P_wimp,drho_wimp_dTvar,rho_wimp_Tnud,P_wimp_Tnud,n_wimp,dn_wimp_dTvar_phiWpart,dn_wimp_dTvar_zpart;
    double entropy_wimp_gamma;  // Entropy of WIMP normalized to the photon entropy -> s(m_WIMP)/s_gamma
    double phi_chid;            // Entropy of WIMP at Tnud normalized to entropy of massless WIMP -> s(m_WIMP)/s(0)
    double Tvar=0.;
 
    drho_neutrinos_dTnu=P_neutrinos=0.;

   //-----------------------------------------------------------------
    rho_gamma=pow(pi,2.)/15.*pow(T,4.);
    if (paramrelic->wimp)
    {
        if (paramrelic->fermion)
        {
            rho_wimp=(Mbessel(zW)-Mbessel(2.*zW)+Mbessel(3.*zW)-Mbessel(4.*zW)+Mbessel(5.*zW)-Mbessel(6.*zW)
                      +Mbessel(7.*zW))*pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi;
            P_wimp=(Lbessel(zW)/zW-Lbessel(2.*zW)/(zW*2.)+Lbessel(3.*zW)/(zW*3.)-Lbessel(4.*zW)/(zW*4.)
                    +Lbessel(5.*zW)/(zW*5.)-Lbessel(6.*zW)/(zW*6.)+Lbessel(7.*zW)/(zW*7.))
                    *pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi;

        }
        else // bosonic wimp
        {
            rho_wimp=(Mbessel(zW)+Mbessel(2.*zW)+Mbessel(3.*zW)+Mbessel(4.*zW)+Mbessel(5.*zW)+Mbessel(6.*zW)
                      +Mbessel(7.*zW))*pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi;
            P_wimp=(Lbessel(zW)/zW+Lbessel(2.*zW)/(zW*2.)+Lbessel(3.*zW)/(zW*3.)+Lbessel(4.*zW)/(zW*4.)
                    +Lbessel(5.*zW)/(zW*5.)+Lbessel(6.*zW)/(zW*6.)+Lbessel(7.*zW)/(zW*7.))
                    *pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi;

        }

        if (paramrelic->EM_coupled)
        {
            entropy_wimp_gamma = (rho_wimp + P_wimp) / (4./3. * rho_gamma);
        }
        else
        {
            // Neutrino coupled WIMPs do not contribute to the total entropy IF we start at neutrino decoupling
            entropy_wimp_gamma = 0.;
        }
    }
    else
    {
        rho_wimp = 0.;
        P_wimp = 0.;
        drho_wimp_dTvar = 0.;
        entropy_wimp_gamma = 0.;
    }

    /* Find s_e/s_gamma assuming zero phie at early times */
    rho_epem=(Mbessel(z)-Mbessel(2.*z)+Mbessel(3.*z)-Mbessel(4.*z)+Mbessel(5.*z)-Mbessel(6.*z)+Mbessel(7.*z))*2.*pow(m_e,4.)/pow(pi,2.);
    P_epem=(Lbessel(z)/z-Lbessel(2.*z)/(z*2.)+Lbessel(3.*z)/(z*3.)-Lbessel(4.*z)/(z*4.)+Lbessel(5.*z)/(z*5.)-
            Lbessel(6.*z)/(z*6.)+Lbessel(7.*z)/(z*7.))*2.*pow(m_e,4.)/pow(pi,2.);
    double entropy_epem_gamma=(rho_epem+P_epem) / (4./3.*rho_gamma);

    /* The late-time (CMB-measured) value of eta (eta0) given by the user, together with entropy conservation is
     * used to approximate an initial value of eta, here parameterized through h_eta. eta is not constant
     * prior to e+- annihilation and is therefore evolved together with the temperature. */
    double h_eta=paramrelic->eta0*M_u*g_to_GeV*2.*zeta3/pow(pi,2.)*(1. + entropy_epem_gamma + entropy_wimp_gamma);

    double phie=h_eta*Y[2]/(M_u*g_to_GeV)*pow(pi,2.)/(2.*pow(z,3.)*(Lbessel(z)-Lbessel(2.*z)*2.+Lbessel(3.*z)*3.-Lbessel(4.*z)*4.
                                                      +Lbessel(5.*z)*5.-Lbessel(6.*z)*6.+Lbessel(7.*z)*7.));
    double rhob0=h_eta*pow(T,3.);              // Initial baryon density

    /* ############ FIND INITIAL TIME ########### */
    /* Strictly valid in the limit T->infty, but valid assumption for the high initial temp. here. */
    double t=sqrt(12.*pi*G*sigma_SB)/pow(Ti,2.);

	if(paramrelic->phi_model&&paramrelic->rho_phi!=0.) if(1.-dt*(paramrelic->Gamma_phi)<0.) paramrelic->rho_phi=0.;
 
    if(paramrelic->wimp||rhod!=0.)
    {
        /* In the presence of a WIMP H is larger at a given time, thus t is smaller
         * at a given T. Assuming H propto t early on this leads to a correction factor H_SBBN/H_WIMP.
         * The same reasoning holds for any arbitrary dark density implemented. */
        double H_SBBN, H_WIMP;
        cosh1=cosh(phie);
        cosh2=cosh(phie*2.);
        cosh3=cosh(phie*3.);
        cosh4=cosh(phie*4.);
        cosh5=cosh(phie*5.);
        cosh6=cosh(phie*6.);
        cosh7=cosh(phie*7.);
               
        rho_neutrinos=neutdens(Tnu,paramrelic);
        rho_neuteq=paramrelic->dNnu*pow(Tnu_eq,4.);
        rho_epem=(Mbessel(z)*cosh1-Mbessel(2.*z)*cosh2+Mbessel(3.*z)*
                  cosh3-Mbessel(4.*z)*cosh4+Mbessel(5.*z)*cosh5-
                  Mbessel(6.*z)*cosh6+Mbessel(7.*z)*cosh7)*2.*pow(m_e,4.)/pow(pi,2.);
                  
        H_SBBN=sqrt(G*8.*pi/3.*(rho_gamma+rho_epem+rho_neutrinos+rho_neuteq+rho_baryons));
        H_WIMP=sqrt(G*8.*pi/3.*(rho_gamma+rho_epem+rho_wimp+rho_neutrinos+rho_neuteq+rho_baryons+rhod));
        
        t=t*H_SBBN/H_WIMP;
    }
	
    Y[3]=1/(reacparam[12][8])*1.e-10*Y[1]*Y[2]*(rhob0/g_to_GeV*pow(cm_to_GeV,3.))*exp(reacparam[12][9]*K_to_GeV/T)/pow(T/K_to_GeV,1.5); // Initial statistical equilibrium of deuterium formation pre-BBN
	
	Y0[3]=Y[3];
    
	for (i = 4; i <= NNUC; ++i)
    {
        Y[i]=Ytmin;
        Y0[i]=Y[i];
    }
    rate_weak(f,paramrelic,paramerror);

    /* ############## MAIN LOOP ############ */
    
    rhod=rhodprev=dark_density(T,paramrelic);
    
    while(ltime == 0)
    {				
        for(loop=1;loop<=2;loop++)
        {
            /* ########### DARK ENERGY DENSITY AND ENTROPY ########### */
			if(paramrelic->phi_model&&paramrelic->rho_phi!=0.) if(1.-dt*(paramrelic->Gamma_phi)<0.) paramrelic->rho_phi=0.;
			 
            rhod=dark_density(T,paramrelic);
            Pd=dark_density_pressure(T,paramrelic);
			if(fabs(T-Tprev)>1.e-14)
			{
				if(rhod==0.) drhod_dT=0.; else drhod_dT=(rhod-rhodprev)/(T-Tprev);
 				Tprev=T;
				rhodprev=rhod;
			}
            sd=dark_entropy(T,paramrelic);
            dsd_dT=dark_entropy_derivative(T,paramrelic);
                                             
            Sigmarad=entropy_Sigmarad(T,paramrelic);
            
            z=m_e/T;
            if (paramrelic->wimp)
            {
                if (paramrelic->EM_coupled) Tvar=T;
                else Tvar=Tnu;
                zW = wimp_mass_ratio*z*T/Tvar;
            }

            if (phie<=17.)
            {
                cosh1=cosh(phie);
                cosh2=cosh(phie*2.);
                cosh3=cosh(phie*3.);
                cosh4=cosh(phie*4.);
                cosh5=cosh(phie*5.);
                cosh6=cosh(phie*6.);
                cosh7=cosh(phie*7.);
                sinh1=sinh(phie);
                sinh2=sinh(phie*2.);
                sinh3=sinh(phie*3.);
                sinh4=sinh(phie*4.);
                sinh5=sinh(phie*5.);
                sinh6=sinh(phie*6.);
                sinh7=sinh(phie*7.);
            }
            else
            {
                cosh1=0.;
                cosh2=0.;
                cosh3=0.;
                cosh4=0.;
                cosh5=0.;
                cosh6=0.;
                cosh7=0.;
                sinh1=0.;
                sinh2=0.;
                sinh3=0.;
                sinh4=0.;
                sinh5=0.;
                sinh6=0.;
                sinh7=0.;
            }

            /* ############ PHOTON DENSITY AND PRESSURE ############# */
            rho_gamma=pow(pi,2.)/15.*pow(T,4.);

            drho_gamma_dT=rho_gamma*4./T;

            P_gamma=rho_gamma/3.;

            /* ############ ELECTRON AND POSITRON DENSITY AND PRESSURE ########## */
            rho_epem=(Mbessel(z)*cosh1-Mbessel(2.*z)*cosh2+Mbessel(3.*z)*
                      cosh3-Mbessel(4.*z)*cosh4+Mbessel(5.*z)*cosh5-
                      Mbessel(6.*z)*cosh6+Mbessel(7.*z)*cosh7)*2.*pow(m_e,4.)/pow(pi,2.);
                      /* rho_e+ + rho_e- */

            drho_epem_dT=z/T*(Nbessel(z)*cosh1-Nbessel(2.*z)*2.*
                                      cosh2+Nbessel(3.*z)*3.*cosh3-
                                      Nbessel(4.*z)*4.*cosh4+Nbessel(5.*z)*
                                      5.*cosh5-Nbessel(6.*z)*6.*cosh6+
                                      Nbessel(7.*z)*7.*cosh7)*2.*pow(m_e,4.)/pow(pi,2.);
                                      /* d(rho_e+ + rho_e-)/d(T) */

            drho_epem_dphie=(Mbessel(z)*sinh1-Mbessel(2.*z)*2.*sinh2+
                             Mbessel(3.*z)*3.*sinh3-Mbessel(4.*z)*4.*sinh4+
                             Mbessel(5.*z)*5.*sinh5-Mbessel(6.*z)*6.*sinh6+
                             Mbessel(7.*z)*7.*sinh7)*2.*pow(m_e,4.)/pow(pi,2.);
                             /* d(rho_e+ + rho_e-)/d(phie) */

            P_epem=(Lbessel(z)*cosh1/z-Lbessel(2.*z)*cosh2/(z*2.)+
                    Lbessel(3.*z)*cosh3/(z*3.)-Lbessel(4.*z)*cosh4/(z*4.)+
                    Lbessel(5.*z)*cosh5/(z*5.)-Lbessel(6.*z)*cosh6/(z*6.)+
                    Lbessel(7.*z)*cosh7/(z*7.))*2.*pow(m_e,4.)/pow(pi,2.); /* P_e+ + P_e- */

            /* ########## WIMP DENSITY AND PRESSURE ########### */
            if (paramrelic->wimp)
            {
                if (paramrelic->fermion) // fermionic wimp
                {
                    rho_wimp=(Mbessel(zW)-Mbessel(2.*zW)+Mbessel(3.*zW)-Mbessel(4.*zW)
                              +Mbessel(5.*zW)-Mbessel(6.*zW)+Mbessel(7.*zW))
                            *pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi;
                    P_wimp=(Lbessel(zW)/zW-Lbessel(2.*zW)/(zW*2.)+Lbessel(3.*zW)/(zW*3.)
                            -Lbessel(4.*zW)/(zW*4.)+Lbessel(5.*zW)/(zW*5.)-Lbessel(6.*zW)/(zW*6.)
                            +Lbessel(7.*zW)/(zW*7.))*pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi;
                    /* Note that for self conj. particles, the factor cosh(n*phiW) is replaced by exp(n*phiW).
                     * However, we assume phiW=0 for those particles, and the double counting of
                     * non-self conj. particles is baked into the definition of g_chi, so it makes no difference,
                     * since cosh(0)=exp(0)=1. If, in the future, a varying phiW is to be implemented, the
                     * difference must be taken into account */
                    drho_wimp_dTvar=(1./Tvar)*pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi*
                            (zW*Nbessel(zW)-
                             2*zW*Nbessel(2.*zW)+
                             3*zW*Nbessel(3.*zW)-
                             4*zW*Nbessel(4.*zW)+
                             5*zW*Nbessel(5.*zW)-
                             6*zW*Nbessel(6.*zW)+
                             7*zW*Nbessel(7.*zW));
                }
                else // bosonic wimp
                {
                    rho_wimp=(Mbessel(zW)+Mbessel(2.*zW)+Mbessel(3.*zW)+Mbessel(4.*zW)
                              +Mbessel(5.*zW)+Mbessel(6.*zW)+Mbessel(7.*zW))
                            *pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi;
                    P_wimp=(Lbessel(zW)/zW+Lbessel(2.*zW)/(zW*2.)+Lbessel(3.*zW)/(zW*3.)
                            +Lbessel(4.*zW)/(zW*4.)+Lbessel(5.*zW)/(zW*5.)+Lbessel(6.*zW)/(zW*6.)
                            +Lbessel(7.*zW)/(zW*7.))
                            *pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi;
                    drho_wimp_dTvar=(1./Tvar)*pow(paramrelic->m_chi/1.e+3,4.)*paramrelic->g_chi*
                            (zW*Nbessel(zW)+
                             2*zW*Nbessel(2.*zW)+
                             3*zW*Nbessel(3.*zW)+
                             4*zW*Nbessel(4.*zW)+
                             5*zW*Nbessel(5.*zW)+
                             6*zW*Nbessel(6.*zW)+
                             7*zW*Nbessel(7.*zW));
                }
            }

            /* ########### NEUTRINO DENSITY ############ */
            rho_neutrinos=neutdens(Tnu,paramrelic);
            
            if ((paramrelic->wimp)&&((paramrelic->neut_coupled)||(paramrelic->neuteq_coupled)))
            {
                P_neutrinos=rho_neutrinos/3.;
                drho_neutrinos_dTnu=neutdens_deriv(Tnu,paramrelic);
                if (paramrelic->neut_coupled)
                {
                    // Eq. neutrinos have their own temperature and are decoupled, thus derivative is zero
                    rho_neuteq=2.*pow(pi,2.)/30.*7./8.*paramrelic->dNnu*pow(Tnu_eq,4.);
                    drho_neuteq=0.;
                }
                else
                {
                    // Common SM and eq. neutrino temperature
                    rho_neuteq=2.*pow(pi,2.)/30.*7./8.*paramrelic->dNnu*pow(Tnu,4.);
                    drho_neuteq=7.*pow(pi,2.)/30.*paramrelic->dNnu*pow(Tnu,3.);
                }
                P_neuteq=rho_neuteq/3.;
            }
            else
            {
                // SM and eq. neutrinos share same temperature
                rho_neuteq=2.*pow(pi,2.)/30.*7./8.*paramrelic->dNnu*pow(Tnu,4.);
                P_neuteq=rho_neuteq/3.;
                drho_neuteq=0.;
            }

            /* ############### BARYON DENSITY ################ */
            rho_baryons=h_eta*pow(T,3.);

            dM_epem_dT=-(pow(z,3.)/T)*(sinh1*(Lbessel(z)*3.-z*Mbessel(z))-sinh2*
                                     (Lbessel(2.*z)*3.-z*2.*Mbessel(2.*z))+
                                     sinh3*(Lbessel(3.*z)*3.-z*3.*
                                            Mbessel(3.*z))-sinh4*
                                     (Lbessel(4.*z)*3.-z*4.*Mbessel(4.*z))+
                                     sinh5*(Lbessel(5.*z)*3.-z*5.*
                                            Mbessel(5.*z))-sinh6*
                                     (Lbessel(6.*z)*3.-z*6.*Mbessel(6.*z))+
                                     sinh7*(Lbessel(7.*z)*3.-z*7.*
                                            Mbessel(7.*z)));
                                     /* d(pi^2 (ne- - ne+)*z^3 /
                                      * 2 m^3) / d(T) */

            dN_epem_dphie=pow(z,3.)*(cosh1*Lbessel(z)-cosh2*2.*Lbessel(2.*z)+
                                 cosh3*3.*Lbessel(3.*z)-cosh4*4.*Lbessel(4.*z)+
                                 cosh5*5.*Lbessel(5.*z)-cosh6*6.*Lbessel(6.*z)+
                                 cosh7*7.*Lbessel(7.*z));
            if(dN_epem_dphie!=0.) dN_epem_dphie=1./dN_epem_dphie;
                        /* d(pi^2/2 1/M_u h sum Z_i Y_i)/d(phie) */

            // Summing up all energy densities from different sources
            H=sqrt(G*8.*pi/3.*(rho_gamma+rho_epem+rho_wimp+rho_neutrinos+rho_neuteq+rho_baryons+rhod));

            rate_pn(f,r,T/K_to_GeV,Tnu/K_to_GeV,paramrelic,paramerror); // conversion to CGS units
            f[1]*=norm;
            r[1]*=norm;

            rate_all(f,T/K_to_GeV,paramrelic,paramerror); // conversion to CGS units

            fail=linearize(T/K_to_GeV,reacparam,f,r,loop,inc,ip,dt/s_to_GeV,Y0,Y,dY_dt,H*s_to_GeV,rho_baryons/g_to_GeV*pow(cm_to_GeV,3.)); // conversion to CGS units
			
			for(i=1;i<=NNUC;i++)
			{
				dY_dt[i] = dY_dt[i]/s_to_GeV; // we get derivatives back in GeV
			}
			
            if(fail!=0)
			{
				paramrelic->failsafe=1;
				return 1;
			}
			
            sum_Y=0.;
            sum_ZY=0.;
            sum_dY_dt=0.;
            sum_DeltaMdY_dt=0.;
            sum_ZdY_dt=0.;

            for (i=1;i<=NNUC;i++)
            {
                sum_Y+=Y[i];
                sum_ZY+=Zm[i]*Y[i];
                sum_dY_dt+=dY_dt[i];
                sum_DeltaMdY_dt+=Dm[i]/1000.*dY_dt[i]/(M_u*g_to_GeV);
                sum_ZdY_dt+=Zm[i]*dY_dt[i];
            }

            dphie_dT=dN_epem_dphie*(-3.*pow(pi,2.)/(2.*M_u*g_to_GeV)*h_eta*sum_ZY/T-dM_epem_dT);
            dphie_dlna3=-dN_epem_dphie*pow(pi,2.)/(2.*M_u*g_to_GeV)*h_eta*sum_ZY;
            dphie_dZY=dN_epem_dphie*pow(pi,2.)/(2.*M_u*g_to_GeV)*h_eta;

            if((paramrelic->wimp)&&((paramrelic->neut_coupled)||(paramrelic->neuteq_coupled)))
            {
                // WIMPs are coupled to neutrinos, so need to dynamically vary Tnu
                if (paramrelic->neuteq_coupled)
                {
                    dlna3_dTnu=-(drho_neutrinos_dTnu+drho_neuteq+drho_wimp_dTvar)/(rho_neutrinos+P_neutrinos+rho_neuteq
                                                                                   +P_neuteq+rho_wimp+P_wimp);
                }
                else
                {
                    dlna3_dTnu=-(drho_neutrinos_dTnu+drho_wimp_dTvar)/(rho_neutrinos+P_neutrinos+rho_wimp+P_wimp);
                }

                dTnu_dt=3.*H/dlna3_dTnu;
                dlnTnu_dt=dTnu_dt/Tnu;
                // No WIMP contribution to dlna3_dT, since they are not EM coupled
                dlna3_dT=-(drho_gamma_dT+drho_epem_dT+drho_epem_dphie*dphie_dT+rho_baryons*zeta*sum_Y
                            +drhod_dT-T*dsd_dT-T*Sigmarad/(H*3.)*dlna3_dT)/
                        (rho_gamma+P_gamma+rho_epem+P_epem+rho_baryons*(2./3.*zeta*T*sum_Y+zeta*T*sum_dY_dt/(H*3.)
                                                                                        +sum_DeltaMdY_dt/(H*3.))
                         +rhod+Pd-T*sd+drho_epem_dphie*(dphie_dlna3+dphie_dZY*sum_ZdY_dt/(H*3.)));
            }
            else
            {
                // No WIMPs (relevant WIMP parameters set to 0 earlier), or EM coupled WIMPs
                dlna3_dT=-(drho_gamma_dT+drho_epem_dT+drho_epem_dphie*dphie_dT+(rho_baryons*K_to_GeV)*zeta*sum_Y
                            +drhod_dT-T*dsd_dT-T*Sigmarad/(H*3.)*dlna3_dT)/
                        (rho_gamma+P_gamma+rho_epem+P_epem+rho_baryons*(2./3.*zeta*T*sum_Y+zeta*T*sum_dY_dt/(H*3.)
                                                                                        +sum_DeltaMdY_dt/(H*3.))
                         +rhod+Pd-T*sd+drho_epem_dphie*(dphie_dlna3+dphie_dZY*sum_ZdY_dt/(H*3.)));
            }
									
            dT_dt=3.*H/dlna3_dT;
            dlnT_dt=dT_dt/T;
            dh_dt=-h_eta*(H*3.+dlnT_dt*3.);
            dphie_dt=dphie_dT*dT_dt+dphie_dlna3*(H*3.)+dphie_dZY*sum_ZdY_dt;
                        
            if(paramrelic->phi_model) paramrelic->rho_phi=max(0.,paramrelic->rho_phi-dt*(3.*H+paramrelic->Gamma_phi)*paramrelic->rho_phi);
            
            if (T <= Tf || dt < fabs(1.e-16 / dlnT_dt) || ip == inc)
            {
                it++;
                for (i=1;i<=NNUC;i++) ratioH[i]=Y[i]/Y[2];

                ratioH[2]=Y[2]*Am[2];
                ratioH[6]=Y[6]*Am[6];
                for(i=1;i<=9;i++) ratioH[10]+=ratioH[i];
                ratioH[10]-=1.;
                ratioH[0] = h_eta / (M_u*g_to_GeV*2.*zeta3/pow(pi,2.)); // This is the value of eta
                if((it==nitmax)||(ip<inc)) ltime = 1;
            }

            if(loop==1)
            {
                if(ip==inc) ip=0;
                ip++;
                is++;
                if(is>3)
                {
                    dtmin=fabs(1./dlnT_dt)*ct;
                    for (i=1;i<=NNUC;i++)
                    {
                        if ((dY_dt[i]!=0.)&&(Y[i]>Ytmin))
                        {
                            dtl=(fabs(Y[i]/dY_dt[i]))*cy*(pow(log10(Y[i])/log10(Ytmin),2.)+1.);
                            if (dtl<dtmin) dtmin=dtl;
                        }
                    }
                    if (dtmin>dt*1.5) dtmin=dt*1.5;
                    dt=dtmin;
                }

                t+=dt;

                T0=T;
                h_eta0=h_eta;
                phie0=phie;

                dT0_dt=dT_dt;
                dh_dt0=dh_dt;
                dphie_dt0=dphie_dt;

                T=T0+dT0_dt*dt;
							
                if(T<0||isnan(T)||isinf(T))
				{
					paramrelic->failsafe=1;
					return 1;
                }
				
                h_eta=h_eta0+dh_dt0*dt;
                phie=phie0+dphie_dt0*dt;

                if ((paramrelic->wimp)&&((paramrelic->neut_coupled)||(paramrelic->neuteq_coupled)))
                {
                    // Tnu as dynamic variable in case of neutrino coupled WIMPs
                    Tnu0=Tnu;
                    dTnu0_dt=dTnu_dt;
                    if (T>=Tnud) Tnu=T;
                    else Tnu=Tnu0+dTnu0_dt*dt;
                    if ((paramrelic->neut_coupled)&&(T<Tnud))
                    {
                        // If WIMPs are coupled only to SM neutrinos, Tnu_eq is proportional to a^-1 after neutrino decoupling
                        Tnu_eq=pow(h_eta*pow(T,3.)/rhob0,1./3.)*Ti;
                    }
                    else
                    {
                        // If WIMPs are coupled to both sets of neutrinos, they all share the same temperature.
                        Tnu_eq=Tnu;
                    }
                }
                else
                {
                    // If not neutrino coupled WIMPs, Tnu and Tnu_eq is proportional to a^-1 after neutrino decoupling
                    if (T>Tnud) Tnu=T;
                    else Tnu=pow(h_eta*pow(T,3.)/rhob0,1./3.)*Ti;
                    Tnu_eq=Tnu;
                }

                for (i=1;i<=NNUC;i++)
                {
                    Y0[i]=Y[i];
                    dY_dt0[i]=dY_dt[i];
                    Y[i]=Y0[i]+dY_dt0[i]*dt;
                    if(Y[i]<Ytmin) Y[i]=Ytmin;
                }
            }
            else /* if(loop==2) */
            {
                T=T0+(dT_dt+dT0_dt)*0.5*dt;
                h_eta=h_eta0+(dh_dt+dh_dt0)*0.5*dt;
                phie=phie0+(dphie_dt+dphie_dt0)*0.5*dt;

                if ((paramrelic->wimp)&&((paramrelic->neut_coupled)||(paramrelic->neuteq_coupled)))
                {
                    if (T>=Tnud) Tnu=T;
                    else Tnu=Tnu0+(dTnu_dt+dTnu0_dt)*0.5*dt;
                    if ((paramrelic->neut_coupled)&&(T<Tnud)) Tnu_eq=pow(h_eta*pow(T,3.)/rhob0,1./3.)*Ti;
                    else Tnu_eq=Tnu;
                }
                else
                {
                    if (T>Tnud) Tnu=T;
                    else Tnu=pow(h_eta*pow(T,3.)/rhob0,1./3.)*Ti;
                    Tnu_eq=Tnu;
                }

                for (i=1;i<=NNUC;i++)
                {
                    Y[i]=Y0[i]+(dY_dt[i]+dY_dt0[i])*0.5*dt;
                    if (Y[i]<Ytmin) Y[i]=Ytmin;
                }
            }
        }
    }

    ratioH[8]+=ratioH[9];           // Be7 -> Li7 post BBN
    ratioH[5]+=ratioH[4];           // H3 -> He3 post BBN

    for (i=1;i<=NNUC;i++) ratioH[i]=max(0.,ratioH[i]);

	if(paramrelic->phi_model) paramrelic->rho_phi=0.;

	if(ratioH[3]<1.e-10||ratioH[3]<1.e-10) /* H2_H & Yp */
	{
		paramrelic->failsafe=1;
		return 1;
	}
	
    return 0;
}

/*----------------------------------------------------*/

int nucl_err(struct relicparam* paramrelic, double ratioH[NNUC+1], double cov_ratioH[NNUC+1][NNUC+1])
/* Routine which computes the abundance ratios (in ratioH[]) and their
 * covariance matrice (in cov_ratioH[][]), using the parameters contained in
 * paramrelic->err. The err parameter is a switch to choose the evaluation error
 * method (0=no error, 1=high values of the nuclear rates, 2=low values,
 * 3=linear error calculation, 4=random Gaussian error calculation). */
{	
    int ie,je;    
    
    for(ie=1;ie<=NNUC;ie++) ratioH[ie]=0.;
    for(ie=1;ie<=NNUC;ie++) for(je=1;je<=NNUC;je++) cov_ratioH[ie][je]=0.;

    if(paramrelic->err==0)
    {
        if(nucl(paramrelic,ratioH)>0) return 0;
    }
    
    else if(paramrelic->err==1||paramrelic->err==2)
    {
		double ratioH_tmp[NNUC+1];
        if(nucl(paramrelic,ratioH_tmp)>0) return 0;
        paramrelic->err=0;
        if(nucl(paramrelic,ratioH)>0) return 0;
        for(ie=1;ie<=NNUC;ie++) cov_ratioH[ie][ie]=fabs(ratioH_tmp[ie]-ratioH[ie]);
    }
    
    else if(paramrelic->err==3)
    {
		int optfail=0;
		int checkzeros=0;
		
		paramrelic->err=0;
        if(nucl(paramrelic,ratioH)>0) optfail=1;
        for(ie=1;ie<=NNUC;ie++) optfail+=isnan(ratioH[ie]);
		paramrelic->err=3;
		
		double ratioH_all[NNUCREAC+2][NNUC+1];
        
#if defined(_OPENMP)	
#pragma omp parallel for
#endif
        for(ie=0;ie<=NNUCREAC+1;ie++)
        {
			int je;
			double ratioH_tmp[NNUC+1];
 
            if(optfail==0)
            {
				struct errorparam paramerror;
				paramerror.errnumber=ie;
				
				for(je=1;je<=NNUC;je++) ratioH_tmp[je]=0.;
																
				if(nucl_single(paramrelic,ratioH_tmp,&paramerror)>0) optfail=1;

				if(ratioH_tmp[3]*ratioH_tmp[6]==0.) checkzeros+=1;

				for(je=1;je<=NNUC;je++) ratioH_all[ie][je]=ratioH_tmp[je];

				for(je=1;je<=NNUC;je++) optfail+=isnan(ratioH_tmp[je]);
            }
        }

		if(checkzeros>10) optfail=1;

        if(optfail>0)
        {			
			 if(paramrelic->failsafe!=1)
			 {
#ifdef DEBUG
			 	printf("Sorry, more precise calculation required, please wait...\n");
#endif
			 	paramrelic->failsafe=1;
			 	return nucl_err(paramrelic,ratioH,cov_ratioH);
			 }
			 else return 0;
		 }
		 else
		 {
			 int je,ke;
		
			for(ie=0;ie<=NNUCREAC+1;ie++) for(je=1;je<=NNUC;je++) for(ke=1;ke<=NNUC;ke++) if(ratioH_all[ie][je]*ratioH_all[ie][ke]!=0.) cov_ratioH[je][ke]+=(ratioH_all[ie][je]-ratioH[je])*(ratioH_all[ie][ke]-ratioH[ke]);

			for(je=1;je<=NNUC;je++) optfail+=isnan(ratioH[je])+(cov_ratioH[je][je]<0.)+(sqrt(cov_ratioH[je][je])/ratioH[je]<1.e-10);
			for(je=1;je<=NNUC;je++) for(ke=1;ke<=NNUC;ke++) optfail+=isnan(cov_ratioH[je][ke])+(pow(cov_ratioH[je][ke],2.)>1.0001*fabs(cov_ratioH[je][je]*cov_ratioH[ke][ke]));
		 }
        
    }
    
    else if(paramrelic->err==4)
    {
		int optfail=0;
				
		int niter=1000;

		srand((unsigned int)(getpid()));
		
		double ratioH_all[niter][NNUC+1];
				
#if defined(_OPENMP)	
#pragma omp parallel for
#endif
		for(ie=0;ie<niter;ie++)
		{
			int je;
			struct errorparam paramerror;
			double ratioH_tmp[NNUC+1];
			for(je=1;je<=NNUC;je++) ratioH_tmp[je]=0.;
			
			if(optfail==0)
            {
				for(je=0;je<=NNUCREAC+1;je++) paramerror.random[je]=rand_gauss();
							
				if(nucl_single(paramrelic,ratioH_tmp,&paramerror)>0) optfail=1;
						
				for(je=1;je<=NNUC;je++) ratioH_all[ie][je]=ratioH_tmp[je];
				
				for(je=1;je<=NNUC;je++) optfail+=isnan(ratioH_all[ie][je]);
			}
		}
						
		int je,ke;

		for(ie=0;ie<niter;ie++) for(je=1;je<=NNUC;je++) ratioH[je]+=ratioH_all[ie][je];
		for(je=1;je<=NNUC;je++) ratioH[je]/=niter;
		for(ie=0;ie<niter;ie++) for(je=1;je<=NNUC;je++) for(ke=1;ke<=NNUC;ke++) cov_ratioH[je][ke]+=(ratioH_all[ie][je]*ratioH_all[ie][ke]-ratioH[je]*ratioH[ke]);
		for(je=1;je<=NNUC;je++) for(ke=1;ke<=NNUC;ke++) cov_ratioH[je][ke]/=niter;
	
		for(je=1;je<=NNUC;je++) optfail+=isnan(ratioH[je])+(cov_ratioH[je][je]<0.)+(sqrt(cov_ratioH[je][je])/ratioH[je]<1.e-10);
		for(je=1;je<=NNUC;je++) for(ke=1;ke<=NNUC;ke++) optfail+=isnan(cov_ratioH[je][ke])+(pow(cov_ratioH[je][ke],2.)>1.0001*fabs(cov_ratioH[je][je]*cov_ratioH[ke][ke]));
		
		if(optfail>0)
		{			
			if(paramrelic->failsafe!=1)
			{
#ifdef DEBUG
			 	printf("Sorry, more precise calculation required, please wait...\n");
#endif
				paramrelic->failsafe=1;
				return nucl_err(paramrelic,ratioH,cov_ratioH);
			}
			else return 0;
		}
	}
		  
    return 1;
}

/*----------------------------------------------------*/

int nucl(struct relicparam* paramrelic, double ratioH[NNUC+1])
{
	if(paramrelic->err>=3) return -1;
	
	struct errorparam paramerror;
	
	int test=nucl_single(paramrelic,ratioH,&paramerror);
	
	if(paramrelic->failsafe==0)	return test;
	else
	{
		paramrelic->failsafe=1;
		return nucl_single(paramrelic,ratioH,&paramerror);
	}
}

/*----------------------------------------------------*/

int bbn_excluded(struct relicparam* paramrelic)
/* "container" function which computes the abundances of the elements using
 * the parameters of paramrelic and compares them to observational limits.
 * The err parameter is a switch to choose if the central (err=0), high (err=1)
 * or low (err=2) values of the nuclear rates is used. Returns 1 if the
 * considered BBN scenario is allowed, 0 otherwise. */
{	 
    double H2_H,Yp,Li7_H,Be7_H,He3_H,Li6_H;
    double ratioH[NNUC+1];
			
    if(nucl(paramrelic,ratioH)==0)
    {
        H2_H=ratioH[3];
        Yp=ratioH[6];
        Li7_H=ratioH[8];
        Be7_H=ratioH[9];
        He3_H=ratioH[5];
        Li6_H=ratioH[7];

#ifdef DEBUG
		printf("Yp\t\t H2/H\t\t He3/H2\t\t Li7/H\t\t Li6/Li7\t Be7/H\n");
		printf("%.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H/H2_H,Li7_H,Li6_H/Li7_H,Be7_H);
#endif		
        if(isnan(Yp)||isnan(H2_H)||isnan(He3_H/H2_H)||isnan(Li7_H)||isnan(Li6_H/Li7_H)||isnan(Be7_H)) return -1;
        if((Yp<0.258)&&((H2_H>1.2e-5)&&(H2_H<5.3e-5))&&(He3_H/H2_H<1.52)&&(Li7_H>0.85e-10)&&(Li6_H/Li7_H<0.66)) return 0;
        /* Conservative intervals from hep-ph/0604251 */
        else return 1;
    }
    else return -1;
}

/*----------------------------------------------------*/

int bbn_excluded_chi2(struct relicparam* paramrelic)
{
	int ie,je;
    double ratioH[NNUC+1],cov_ratioH[NNUC+1][NNUC+1];
	
    int nobs=2;
    double **cov,**invcov,prediction[NNUC+1],observed[NNUC+1],sigmaobs[NNUC+1];
    int translate[nobs];
    
    cov=(double **) malloc(nobs*sizeof(double *));
	for(ie=0;ie<nobs;ie++) cov[ie]=(double *) malloc(nobs*sizeof(double));	

    invcov=(double **) malloc(nobs*sizeof(double *));
	for(ie=0;ie<nobs;ie++) invcov[ie]=(double *) malloc(nobs*sizeof(double));	

    
    if(paramrelic->err>=3)
	{
		if(nucl_err(paramrelic,ratioH,cov_ratioH)==0) return -1;
				
		/* H2_H=ratioH[3]
		   Yp=ratioH[6]
		   He3_H=ratioH[5]
		   Li7_H=ratioH[8]
		* */
		
		translate[0]=3;
		translate[1]=6;
		translate[2]=5;
		translate[3]=8;
		
		observed[0]=2.527e-5;
		observed[1]=0.2449;
		observed[2]=1.1e-5;
		observed[3]=1.6e-10;

		sigmaobs[0]=0.030e-5;
		sigmaobs[1]=0.0040;
		sigmaobs[2]=0.2e-5;
		sigmaobs[3]=0.3e-10;
		
		for(ie=0;ie<nobs;ie++) prediction[ie]=ratioH[translate[ie]];
		
		for(ie=0;ie<nobs;ie++) for(je=0;je<nobs;je++) cov[ie][je]=pow(sigmaobs[ie],2.)*(ie==je)+cov_ratioH[translate[ie]][translate[je]];
		
		if(invert_matrix(nobs,cov,invcov)==0) return -1;
		
		double chi2=0.;
		for(ie=0;ie<nobs;ie++) for(je=0;je<nobs;je++) chi2+=(prediction[ie]-observed[ie])*invcov[ie][je]*(prediction[je]-observed[je]);

		double CL;
		
		switch(nobs)
		{
			case 1: CL=4.; break;
			case 2: CL=6.18; break;
			case 3: CL=8.02; break;
			case 4: CL=9.72; break;
			case 5: CL=11.31; break;
		}

		paramrelic->chi2=chi2;
		paramrelic->nobs=nobs;

		if(chi2>CL) return 1; else return 0;		
	}
	else return -1;
}
