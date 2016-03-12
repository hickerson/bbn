#include "bbnio.h"

void printlables() {
	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
}

// TODO, change to vector< >
void printratios(const char *lable, const NuclideIndex ni[], NuclideMap & nm) {
	printf("%s\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n", lable,
            nm[He4], nm[H2], nm[He3], nm[Li7], nm[Li6], nm[Be7]);
}

int compute_ratios(CosmologyModel relic, NuclideIndex ni[], 
                   NuclideMap & ratioH, NuclideMap & sigma_ratioH)
{ 
    printlables();
	nucl(1, relic, ratioH);
	printratios("  low:", ni, ratioH);

	nucl(0, relic, ratioH);
	printratios(" cent:", ni, ratioH);
	
	nucl(2, relic, ratioH);
	printratios(" high:", ni, ratioH);
			
	if(nucl_witherrors(3, relic, ratioH, sigma_ratioH))
	{
		printf("With uncertainties:\n");
		printlables();
	    printratios("value:", ni, ratioH);
	    printratios("  +/-:", ni, sigma_ratioH);
	}
	
	/*
	int compat=bbn_excluded(ratioH);

	if(compat==1)
        printf("Excluded by BBN constraints\n");
	else if(compat==0)
        printf("Compatible with BBN constraints\n");
	else
        printf("Computation failed\n");

    return compat;
	*/
}

//compute_ratios(relic, ni, ratioH);
//compute_constraints(relic, ni, ratioH);

/**
 * "container" function which computes the abundances 
 * of the elements using the parameters of relic and 
 * compares them to observational limits. 
 * ---------------------------------------------------
 * The err parameter is a switch to choose if the 
 * central (err=0), * high (err=1) * or low (err=2) 
 * values of the nuclear rates is used. Returns 1 if 
 * the considered BBN scenario is allowed, 0 otherwise.
 */
int bbn_excluded(int err, CosmologyModel relic, 
	NuclideIndex ni[], NuclideMap & ratioH)//, NuclideMap & observedHigh, NuclideMap & observedLow)
{	 
	int error = nucl(err, relic, ratioH);
	if(error != 0) {
        printf("Computation failed.\n");
		return -1;
	}
#ifdef DEBUG
	printlables();
	printratios(" mean:", ni, ratioH);
	printratios(" obsH:", ni, observedLow);
	printratios(" obsL:", ni, observedHigh);
#endif	

	for (int i = 0; i < 6; i++) {
		if(isnan(ratioH[ni[i]]))
			printf("Computation diverged.\n");
	}

	/// Conservative intervals from hep-ph/0604251 
	//if ((Yp<0.258)&&((H2_H>1.2e-5)&&(H2_H<5.3e-5))
	//&&(He3_H/H2_H<1.52)&&(Li7_H>0.85e-10)&&(Li6_H/Li7_H<0.66))
	if (ratioH[He4] < 0.258
	and ratioH[H2] >1.2e-5
	and ratioH[H2] < 5.3e-5
	and ratioH[He3]/ratioH[H2] < 1.52
	and ratioH[Li7] > 0.85e-10
	and ratioH[Li6]/ratioH[Li7] < 0.66) {
		printf("Compatible with BBN constraints\n");
		return 0;
	}
	else {
		printf("Excluded by BBN constraints\n");
		return 1;
	}

	/**	 TODO make more general to all constaints 
	for (int i = 0; i < 6; i++) {
		if (ratiosH[ni[i]] < observedLow[ni[i]] or 
		    ratiosH[ni[i]] > observedHigh[ni[i]]) {
        	printf("Excluded by BBN constraints\n");
			return 1;
		}
	}
	*/
}

/**
 * 0 < Yp < 0.258
 * 1.2e-5 < H2_H < 5.3e-5
 * 0 < He3_H/H2_H < 1.52
 * 0.85e-10 < Li7_H
 * Li6_H/Li7_H < 0.66
compute_constraints(relic, ni, ratioH) {
	NuclideMap observedMean;
	NuclideMap observedMean;
	NuclideMap observedMean;
}
 */
