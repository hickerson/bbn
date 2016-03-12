#include "bbnio.h"

void printlables() {
	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
}

// TODO, change to vector< >
void printratios(const char *lable, const NuclideIndex ni[], NuclideMap & nm) {
	printf("%s\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n", lable,
            nm[He4], nm[H2], nm[He3], nm[Li7], nm[Li6], nm[Be7]);
}

int computeratios(CosmologyModel relic, NuclideIndex ni[], 
                  NuclideMap & ratioH, NuclideMap & sigma_ratioH)
{ 
    printlables();
	//nucl(ErrorType::low, relic,ratioH);
	nucl(1, relic,ratioH);
	printratios("  low:", ni, ratioH);

	//nucl(ErrorType::mean, relic,ratioH);
	nucl(0, relic,ratioH);
	printratios(" cent:", ni, ratioH);
	
	//nucl(ErrorType::high, relic,ratioH); // not sure yet why this doesn't work in ++11
	nucl(2, relic,ratioH);
	printratios(" high:", ni, ratioH);
			
	//if(nucl_witherrors(ErrorType::gaussian, relic,ratioH, sigma_ratioH))
	if(nucl_witherrors(3, relic,ratioH, sigma_ratioH))
	{
		printf("With uncertainties:\n");
		printlables();
	    printratios("value:", ni, ratioH);
	    printratios("  +/-:", ni, sigma_ratioH);
	}
	
	//int compat=bbn_excluded(ErrorType::mean,relic);
	int compat=bbn_excluded(0,relic);

	if(compat==1)
        printf("Excluded by BBN constraints\n");
	else if(compat==0)
        printf("Compatible with BBN constraints\n");
	else
        printf("Computation failed\n");

    return compat;
}
