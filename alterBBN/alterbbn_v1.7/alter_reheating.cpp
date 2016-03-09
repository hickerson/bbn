#include "src/bbnio.h"

/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/* In this BBN simulation, dark matter parameters can     */
/* varied to test what altering the comological expansion */
/* history does to the abundances.                        */
/*--------------------------------------------------------*/
int main(int argc,char** argv)
{
	double Trh;
	double dd0,ndd;
	double Sigmad0,nSigmad;

	if(argc<6)
  	{
    	printf( " This program needs 5 parameter:\n"
                "   dd0     dark energy proportion at BBN time\n"
                "   ndd     dark energy decrease exponent (preferentially >4)\n"
                "   Sigmad0 dark entropy production proportion at BBN time\n"
                "   nSigmad dark entropy production decrease exponent (-1 for standard reheating)\n"
                "   Trh     reheating temperature (in GeV)\n");
      		exit(1);
  	}
	else
  	{
   		sscanf(argv[1],"%lf",&dd0);     /// dd0     dark energy proportion at BBN time
     	sscanf(argv[2],"%lf",&ndd);     /// ndd     dark energy decrease exponent
     	sscanf(argv[3],"%lf",&Sigmad0); /// Sigmad0 dark entropy production proportion at BBN time
     	sscanf(argv[4],"%lf",&nSigmad); /// nSigmad dark entropy production decrease exponent 
    	sscanf(argv[5],"%lf",&Trh);     /// Trh     reheating temperature [GeV]
  	}
	
	CosmologyModel relic;                               /// Parameters from before BBN
	relic.Init_cosmomodel();                            /// Standard Model constructor
	relic.Init_dark_density(dd0,ndd,Trh);               /// Dark energy parameters
	relic.Init_dark_entropySigmaD(Sigmad0,nSigmad,Trh); /// Dark entropy parameters
	
	NuclideMap ratioH, sigma_ratioH;
	NuclideIndex ni[6] = {He4,H2,He3,Li7,Li6,Be7};

    return computeratios(relic, ni, ratioH, sigma_ratioH);
}
