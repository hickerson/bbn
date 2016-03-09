#include "src/bbnio.h"

/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/* Included is the Fierz interference term, which extends */
/* the electroweak sector to include scalar and tensor    */
/* currents that may result from SUSY, leptoquarks, etc   */
/*--------------------------------------------------------*/
int main(int argc,char** argv)
{ 
	double eta;

	if(argc<2) 
  	{ 
        printf(" This program needs 1 parameter:\n"
               "   eta     value of the baryon-to-photon ratio\n");
        exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%lf",&eta);
  	}
	
	CosmologyModel relic;       /// The parameters from the big bang relic before bbn.
	relic.Init_cosmomodel();    /// Standard Model constructor.
	relic.eta0 = eta;           /// value of the baryon-to-photon ratio.
	
	NuclideMap ratioH, sigma_ratioH;
	NuclideIndex ni[6] = {He4,H2,He3,Li7,Li6,Be7};

    return computeratios(relic, ni, ratioH, sigma_ratioH);
}
