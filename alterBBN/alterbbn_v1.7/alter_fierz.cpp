#include "src/include.h"

void printratios(const char *lable, const NuclideIndex ni[], NuclideMap & nm) {
	printf("%s\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n", lable,
            nm[He4], nm[H2], nm[He3], nm[Li7], nm[Li6], nm[Be7]);
}  

/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/* Included is the Fierz interference term, which extends */
/* the electroweak sector to include scalar and tensor    */
/* currents that may result from SUSY, leptoquarks, etc   */
/*--------------------------------------------------------*/
int main(int argc,char** argv)
{ 
	double eta, nbnu, tau, fierz;

	if(argc<5) 
  	{ 
        printf(" This program needs 4 parameters:\n"
            "   eta     value of the baryon-to-photon ratio\n"
            "   nbnu    number of neutrinos\n"
            "   tau     neutron lifetime\n"
            "   fierz   beta-decay Fierz interference term\n");
        exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%lf",&eta);
  		sscanf(argv[2],"%lf",&nbnu);
  		sscanf(argv[3],"%lf",&tau);
  		sscanf(argv[4],"%lf",&fierz);
  	}
	
	CosmologyModel relic;      /// The parameters from the big bang relic before bbn.
	relic.Init_cosmomodel();	
	relic.Init_fierz(eta,nbnu,tau,fierz);
	
	NuclideMap ratioH, sigma_ratioH;
	NuclideIndex ratioIndex[6] = {He4,H2,He3,Li7,Li6,Be7};
	double H2_H,He3_H,Yp,Li7_H,Li6_H,Be7_H;
	double sigma_H2_H,sigma_He3_H,sigma_Yp,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H;

	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
	nucl(2,relic,ratioH);
//	H2_H=ratioH[H2];Yp=ratioH[He4];Li7_H=ratioH[Li7];Be7_H=ratioH[Be7];He3_H=ratioH[He3];Li6_H=ratioH[Li6];
//	printf("  low:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
	printratios("  low:", ratioIndex, ratioH);

	nucl(0,relic,ratioH);
//	H2_H=ratioH[H2];Yp=ratioH[He4];Li7_H=ratioH[Li7];Be7_H=ratioH[Be7];He3_H=ratioH[He3];Li6_H=ratioH[Li6];
	//printf(" cent:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
	printratios(" cent:", ratioIndex, ratioH);
	
	nucl(1,relic,ratioH);
	H2_H=ratioH[H2];Yp=ratioH[He4];Li7_H=ratioH[Li7];Be7_H=ratioH[Be7];He3_H=ratioH[He3];Li6_H=ratioH[Li6];
	printf(" high:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
			
	if(nucl_witherrors(3,relic,ratioH,sigma_ratioH))
	{
		printf("With uncertainties:\n");
	    H2_H=ratioH[H2]; Yp=ratioH[He4]; Li7_H=ratioH[Li7]; Be7_H=ratioH[Be7]; He3_H=ratioH[He3]; Li6_H=ratioH[Li6];
	    sigma_H2_H=sigma_ratioH[H2];sigma_Yp=sigma_ratioH[He4];sigma_Li7_H=sigma_ratioH[Li7];sigma_Be7_H=sigma_ratioH[Be7];sigma_He3_H=sigma_ratioH[He3];sigma_Li6_H=sigma_ratioH[Li6];
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		
		printf("value:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
		printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);
	}
    printf("b    :\t %.3e\n",relic.fierz);
	
	int compat=bbn_excluded(0,relic);

	if(compat==1)
        printf("Excluded by BBN constraints\n");
	else if(compat==0)
        printf("Compatible with BBN constraints\n");
	else
        printf("Computation failed\n");

	return 1;
}
