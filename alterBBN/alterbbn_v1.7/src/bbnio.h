#ifndef __BBNIO_H__
#define __BBNIO_H__

#include "include.h"



/*--------------------------------------------------------------------*/
void print_lables();
void print_ratios(const char *lable, const NuclideIndex ni[], NuclideMap & nm);
void print_ratios_errors(double, NuclideIndex ni[], NuclideMap & nm, NuclideMap & snm);
int compute_ratios(CosmologyModel relic, NuclideIndex ni[], 
                  NuclideMap & ratioH, NuclideMap & sigma_ratioH);
int bbn_excluded(int err, CosmologyModel relic, 
	NuclideIndex ni[], NuclideMap & ratioH);//, NuclideMap & observedHigh, NuclideMap & observedLow);

#endif  /// __BBNIO_H__
