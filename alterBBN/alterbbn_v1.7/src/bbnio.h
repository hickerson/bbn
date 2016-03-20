#ifndef __BBNIO_H__
#define __BBNIO_H__

#include "include.h"
#include <vector>


using namespace std;

struct distribution {
	double mean;
	double min;
	double max;
	double error;
	int samples;
	vector<double> values;

	distribution(double value) :
		mean(value),
		min(value),
		max(value),
		error(0),
		samples(1) {
		values.push_back(mean);
	}

	/** gaussian distribution */
	distribution(double mean, double error, int samples) :
		mean(mean),
		min(mean-samples*error),
		max(mean+samples*error),
		error(error),
		samples(samples) {
		values.push_back(mean);
	}

	distribution(double mean, double min, double max, int samples) :
		mean(mean),
		min(min),
		max(max),
		error((max-min)/2),
		samples(samples) {
		values.push_back(mean);
	}

	distribution & operator=(double value) {
		mean = value;
		min = value;
		max = value;
		error = 0;
		samples = 1;
		values.clear();
		values.push_back(mean);
		return *this;
	}

	/*
	public:init(const char* str) {
		scanf("%f", str)
	} */
};

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
