#include "src/bbnio.h"

/*
void printlables() {
	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
}

void printratios(const char *lable, const NuclideIndex ni[], NuclideMap & nm) {
	printf("%s\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n", lable,
            nm[He4], nm[H2], nm[He3], nm[Li7], nm[Li6], nm[Be7]);
} 
*/ 

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
	NuclideIndex ni[6] = {He4,H2,He3,Li7,Li6,Be7};

    return computeratios(relic, ni, ratioH, sigma_ratioH);

    /**
	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
	nucl(2,relic,ratioH);
	printratios("  low:", ratioIndex, ratioH);

	nucl(0,relic,ratioH);
	printratios(" cent:", ratioIndex, ratioH);
	
	nucl(1,relic,ratioH);
	printratios(" high:", ratioIndex, ratioH);
			
	if(nucl_witherrors(3,relic,ratioH,sigma_ratioH))
	{
		printf("With uncertainties:\n");
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
	    printratios("value:", ratioIndex, ratioH);
	    printratios("  +/-:", ratioIndex, sigma_ratioH);
	}
    printf("    b:\t %.3e\n",relic.fierz);
	
	int compat=bbn_excluded(0,relic);

	if(compat==1)
        printf("Excluded by BBN constraints\n");
	else if(compat==0)
        printf("Compatible with BBN constraints\n");
	else
        printf("Computation failed\n");
	return 1;
        */
}
