#ifndef _BG_STACK_
#define _BG_STACK_ 



#include <iostream>
#include "math.h"
#include <algorithm>

//My include
#include "CImg.h"
#include "sin_fit.hpp"

using namespace cimg_library;
using namespace std;


struct Bg_estim{
	float inhomogene;
	float sig_noise;
	float alpha;
	float alpha_neg;
	int region_size;
	CImg<float> lf_map;
	Evo_visualizer MAD_Evo;
	Bg_estim(float sig){
		alpha=3;
		alpha_neg=6;
		region_size=7;
		sig_noise=sig;
	}


	template<class T>
	CImg<T> get_rebuild_bg(CImg<T> & stack){
		float Pi=3.141592;	 //  Pi value
		float ws=2.*Pi/12.;
		Sin_fit fitter;
		CImg<float> rebuild_stack(stack,"xyz");
		CImg<float> region;
		lf_map.assign(stack,"xy").fill(0);
		fitter.set_threshold(0.01);
		fitter.set_max_iter(100);
		unsigned int pixel_aprox_number=(stack.width()-region_size+1)*(stack.height()-region_size+1);
		unsigned int count=0;
		unsigned int percentage=0;
		float perc_ratio = 100.f/pixel_aprox_number;

		cout << " Processing pixels (%) : ";
		cimg_for_insideXY(stack,x,y,region_size/2){
			region=stack.get_crop(x-region_size/2,y-region_size/2,
														x+region_size/2,y+region_size/2);
			get_phase_bg(region,sig_noise,fitter);
			percentage=count++*perc_ratio;
			printf("%02d", percentage); 
			fflush(stdout);	
			cimg_forZ(rebuild_stack,z){
				rebuild_stack(x,y,z)= fitter.get_C()+fitter.get_amp()*sin(ws*z+fitter.get_phase());
			}
			lf_map(x,y)=fitter.get_phase();
			cout << "\b\b";
		}

		return rebuild_stack;
	};

	float get_phase_bg(CImg<float> & values,float sig,Sin_fit & fitter){
	// Use of an assymetric cost function to suppress the influence of the particle in the bg
		fitter.assign(values,0);
		fitter.set_starting_value();
		// starting_value_bg(fitter,values);
		fitter.set_threshold(0.01);
		fitter.set_max_iter(100);
		fitter.run();
		// CImg<unsigned char> gr = fitter.display_fit(); 
		
		CImg<float> r = fitter.get_residue();
		// inhomogene=fabs(sqrt(r.get_abs().variance()));
		// cout << "inhomogene : " << inhomogene << endl;

		// float MAD = 1.4826*(r.get_abs() - r.get_abs().median()).get_abs().median();
		// float MAD = sqrt(r.variance(2));
		// MAD_Evo.add_data(MAD);
		
		// Leclerc_asym cf(MAD,alpha_neg*MAD,alpha);
		Leclerc_asym cf(sig,alpha_neg*sig,alpha);
		fitter.set_cost_function(&cf);
		fitter.irls_no_init();
		// CImg<unsigned char> gr1 = fitter.display_fit();
		// (gr,gr1).display();
		return fitter.get_beta(2);
	};
};

#endif /* _VAR_STAB_H_ */
