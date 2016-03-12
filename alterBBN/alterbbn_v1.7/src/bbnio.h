#ifndef __BBNIO_H__
#define __BBNIO_H__

#include "include.h"



/*--------------------------------------------------------------------*/
void printlables();
void printratios(const char *lable, const NuclideIndex ni[], NuclideMap & nm);
int computeratios(CosmologyModel relic, NuclideIndex ni[], 
                  NuclideMap & ratioH, NuclideMap & sigma_ratioH);

#endif  /// __BBNIO_H__
