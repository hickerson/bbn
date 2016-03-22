#include "src/bbnio.h"

/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/* Included is the Fierz interference term, which extends */
/* the electroweak sector to include scalar and tensor    */
/* currents that may result from SUSY, leptoquarks, etc   */
/*--------------------------------------------------------*/
int main(int argc,char** argv)
{ 
	double baryon, neutrinos, lifetime, fierz;

	if(argc<5) 
  	{ 
        printf(" This program needs 4 parameters:\n"
            "   eta     value of the baryon-to-photon ratio\n"
            "   nnu     number of neutrinos\n"
            "   tau     neutron lifetime\n"
            "   fierz   beta-decay Fierz interference term\n");
        exit(1); 
  	} 
	else 
  	{
  		sscanf(argv[1],"%lf",&baryon);
  		sscanf(argv[2],"%lf",&neutrinos);
  		sscanf(argv[3],"%lf",&lifetime);
  		sscanf(argv[4],"%lf",&fierz);
  	}
	
	/*
	distribution d(3,1,5,10);
	printf("size: %d\n",d.samples);
	for (int i = 1; i<10; i++)
		printf("%d: %f\n",i,d.values[i]);
		*/

	CosmologyModel relic;      /// The parameters from the big bang relic before bbn.
	relic.Init_cosmomodel();	
	relic.Init_fierz(baryon, neutrinos, lifetime, fierz);
	
	NuclideMap ratioH, sigma_ratioH;
	NuclideIndex ni[6] = {He4,H2,He3,Li7,Li6,Be7};	/// only display this subset

	printf("Primary central values:\n");	
    compute_ratios(relic, ni, ratioH, sigma_ratioH);
    bbn_excluded(0, relic, ni, ratioH);
	
	/*
	printf("\nFormated bounded range:\n");
	print_lables_bounds("eta", ni);
	//print_lables_bounds("Fierz b", ni);
    for (double b = fierz - 0.2; b <= fierz + 0.2; b += 0.01) {
		for (double eta = 1e-10; eta <= 10e-10; eta*=pow(10,0.1)) {
			for (double tau = lifetime - 1.1; tau <= lifetime + 1.1; tau += 2*1.1/3) {
				double nnu = neutrinos;
				relic.Init_fierz(eta,nnu,tau,fierz);
				if(nucl_witherrors(3, relic, ratioH, sigma_ratioH))
				{
					print_ratios_error_bounds(eta, ni, ratioH, sigma_ratioH);
				}
			}
		}
	}
	*/
}
