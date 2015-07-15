#include "src/include.h"

/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/*--------------------------------------------------------*/

int main()
{ 
	struct relicparam paramrelic;
	double ratioH[NNUC+1],sigma_ratioH[NNUC+1];
	double H2_H,He3_H,Yp,Li7_H,Li6_H,Be7_H;
	double sigma_H2_H,sigma_He3_H,sigma_Yp,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H;
	
	Init_cosmomodel(&paramrelic);
	
	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
	nucl(2,paramrelic,ratioH);
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	printf("  low:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

	nucl(0,paramrelic,ratioH);
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	printf(" cent:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
	
	nucl(1,paramrelic,ratioH);
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
	printf(" high:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
			
	if(nucl_witherrors(3,paramrelic,ratioH,sigma_ratioH))
	{
		printf("With uncertainties:\n");
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
		sigma_H2_H=sigma_ratioH[3];sigma_Yp=sigma_ratioH[6];sigma_Li7_H=sigma_ratioH[8];sigma_Be7_H=sigma_ratioH[9];sigma_He3_H=sigma_ratioH[5];sigma_Li6_H=sigma_ratioH[7];
		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		
		printf("value:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
		printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);
	}
	else printf("Uncertainty calculation failed\n\n");
		
	int compat=bbn_excluded(0,paramrelic);

	if(compat==1) printf("Excluded by BBN constraints\n");
	else if(compat==0) printf("Compatible with BBN constraints\n");
	else printf("Computation failed\n");

	return 1;
}
