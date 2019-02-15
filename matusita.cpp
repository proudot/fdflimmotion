// matusita.cpp
//

#include "matusita.h"
#include <iostream>
#define LZZ_INLINE inline
matusita_distance::matusita_distance (const CImg <float> & values, uint init_bins)
{
		value_min=values.min(); value_max=values.max();
		histogram=values.get_histogram(init_bins,value_min,value_max);
}

float matusita_distance::get_distance (float tau)
{
	float distance=0;
	int bin_threshold=(tau-value_min)*histogram.size()/(value_max-value_min);

	cfl hist_low(histogram.data(),bin_threshold),hist_high(histogram.data(bin_threshold),histogram.size()- bin_threshold);

	float dist_high= hist_high.is_empty() ? 0 : (hist_high.get_sqrt()-sqrt(hist_high.max())).sqr().sum();
	float dist_low= hist_low.is_empty() ? 0 : (hist_low.get_sqrt()-sqrt(hist_low.max())).sqr().sum();
	distance=dist_low+dist_high;
	return distance;
}

float matusita_distance::minimize_distance_exhaustive(){
	float min_tau=0; 
	float min_dist= get_distance(value_min);
	tau_dist.assign(2,histogram.size());
	for( float tau=value_min ; tau < value_max ; tau+=(value_max-value_min)/histogram.size()){
		float adist = get_distance(tau);
		if(min_dist>adist){
			min_dist=adist; min_tau=tau;
		};
		tau_dist.add_data(tau,adist);
	}
	// cout << " min distance : " << min_dist ;
	// cout << " min distance tau : " << min_tau ;
	return min_tau;
}

#undef LZZ_INLINE
