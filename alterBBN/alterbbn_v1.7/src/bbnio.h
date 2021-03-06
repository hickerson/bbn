#ifndef __BBNIO_H__
#define __BBNIO_H__

#include "include.h"
#include <vector>
//#include <boost/lambda/lambda.hpp>
//#include <boost/math/special_functions/erf.hpp>



/*-- input ---------------------------------------------------------*/


/*-- output --------------------------------------------------------*/
void print_lables(const char *title, NuclideIndex ni[]);
void print_lables_errors(const char *title, NuclideIndex ni[]);
void print_lables_bounds(const char *title, NuclideIndex ni[]);
void print_ratios(const char *lable, NuclideIndex ni[],
                                     NuclideMap & nm);
void print_ratios_error_bars(double, NuclideIndex ni[], 
									 NuclideMap & ratioH, 
								     NuclideMap & sigma_ratioH);
void print_ratios_error_bounds(double, NuclideIndex ni[], 
									   NuclideMap & ratioH,
									   NuclideMap & sigma_ratioH);
void print_ratios_bounds(double, NuclideIndex ni[], 
                                 NuclideMap & low_ratioH, 
								 NuclideMap & high_ratioH);
const char * get_nuclide_name(const NuclideIndex ni);
void get_ratio_name(const NuclideIndex ni, char buffer[]);

/*-- compute -------------------------------------------------------*/
int compute_ratios(CosmologyModel relic, NuclideIndex ni[], 
                  						 NuclideMap & ratioH, 
										 NuclideMap & sigma_ratioH);
int bbn_excluded(int err, CosmologyModel relic, NuclideIndex ni[], 
												NuclideMap & ratioH);
											//, NuclideMap & observedHigh, 
											//  NuclideMap & observedLow);

#endif  /// __BBNIO_H__
