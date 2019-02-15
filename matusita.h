// matusita.h
//

#ifndef LZZ_matusita_h
#define LZZ_matusita_h
#include "CImg.h"
#include "evo_visualizer.hpp"
using namespace cimg_library;
using namespace std;


typedef CImg<float> cfl;

#define LZZ_INLINE inline
struct matusita_distance
{
  CImg <float> histogram;
  float value_max;
  float value_min;
	Mult_evo_visualizer tau_dist;
  matusita_distance (const CImg <float> & values, uint init_bins);
  float get_distance (float tau);
	float minimize_distance_exhaustive();
};
#undef LZZ_INLINE
#endif
