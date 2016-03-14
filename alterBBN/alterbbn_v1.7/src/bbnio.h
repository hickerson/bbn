#ifndef __BBNIO_H__
#define __BBNIO_H__

#include "include.h"



/*-- output --------------------------------------------------------*/
void print_lables();
void print_lables(const char *title, const NuclideIndex ni[]);
void print_lables_errors(const char *title, const NuclideIndex ni[]);
void print_ratios(const char *lable, const NuclideIndex ni[],
                                           NuclideMap & nm);
void print_ratios_error_bars(double, const NuclideIndex ni[], 
										   NuclideMap & ratioH, 
										   NuclideMap & sigma_ratioH);
void print_ratios_error_bounds(double, NuclideIndex ni[], 
									   NuclideMap & ratioH,
									   NuclideMap & sigma_ratioH);
void print_ratios_bounds(double, NuclideIndex ni[], 
                                 NuclideMap & low_ratioH, 
								 NuclideMap & high_ratioH);

/*-- compute -------------------------------------------------------*/
int compute_ratios(CosmologyModel relic, NuclideIndex ni[], 
                  						 NuclideMap & ratioH, 
										 NuclideMap & sigma_ratioH);
int bbn_excluded(int err, CosmologyModel relic, NuclideIndex ni[], 
												NuclideMap & ratioH);
											//, NuclideMap & observedHigh, 
											//  NuclideMap & observedLow);

#endif  /// __BBNIO_H__
