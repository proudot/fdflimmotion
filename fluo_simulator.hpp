//
// fluo_simulator.hpp
//
// Created by Philippe Roudot 2012-07-17
//
// Copyright (C) INRIA 2012 all right reserved
// 

#ifndef _FLUO_SIMULATOR_H_
#define _FLUO_SIMULATOR_H_

 #include <iostream>
#include <fstream>
#include "math.h"
#include <assert.h>
#include <list>
#include <stdlib.h>
#include "limits.h"
#include "float.h"

#include "CImg.h"
#include "simu_spot.hpp"

using namespace cimg_library;
using namespace std;

typedef CImg<float> cfl;

float inline frand() { return (float)rand() / (float)(RAND_MAX); }; //returns from 0 to 1
float inline frand(float Min, float Max) { return (frand() * (Max - Min)) + Min; }; //returns from Min to Max

#define NOISE_FREE 0
#define POISSON_NOISE 1
#define GAUSSIAN_NOISE 2
#define POISSON_GAUSSIAN_NOISE 3

typedef CImg<ushort> cus;

class fluo_simulator{
	// Produce a <frame_nb> stack simulating fluorescently tagged
	// protein in moving vesicle inside an observed in 
public:
	fluo_simulator(){
		x_size=600;
		y_size=500;
		bleach_rate=0;
		frame_nb=1;
		spot_number=10;
		spot_patch_size=21;
		spot_sigma=1.5;
		sigma_noise_spot_intensity=500;
		sigma_spot_move=1;
		stochastic_scale=false;
		noise_type=NOISE_FREE;
		SNR=0;
	}

	inline cfl get_simulation(){return simulation;}
	inline list<Simu_spot<ushort> >  get_spot_list(){return spot_list;}
	inline void clear(){simulation.assign();}

	fluo_simulator & set_background(cus & bg){
		background=bg;
		return *this;
	}

	fluo_simulator & set_size(float x, float y){
		x_size = x; y_size= y; return *this;
	}

	fluo_simulator & set_bleach_rate(float b){
		bleach_rate=b; return *this;
	}

	fluo_simulator &  set_frame_nb(int val ){
		frame_nb=val;return *this;
	}

	fluo_simulator &  set_spot_number(int val ){
		spot_number=val; return *this;
	}

	fluo_simulator &  set_sigma_noise_spot_intensity(float val ){
		sigma_noise_spot_intensity=val; return *this;
	}

	fluo_simulator &  set_sigma_spot_move(float val ){
		sigma_spot_move=val; return *this;
	}

	fluo_simulator & set_scaling_matrix(cfl & val ){
		scaling_matrix=val; return *this;
	}

	fluo_simulator & set_stochastic_scale(bool b){
		stochastic_scale=b; return *this;
	}

	fluo_simulator & set_gaussian_noise(){
		noise_type=GAUSSIAN_NOISE; return *this;
	}

	fluo_simulator & set_poisson_gaussian_noise(){
		noise_type=POISSON_GAUSSIAN_NOISE; return *this;
	}

	fluo_simulator & set_poisson_noise(){
		noise_type=POISSON_NOISE; return *this;
	}

	fluo_simulator & set_SNR(float val ){
		SNR=val; return *this;
	}
	
	inline void set_spot_sigma(float val ){spot_sigma=val;};


	// some preconfigured setup
	fluo_simulator & config_moving_spot_no_bg(){
		x_size=400;
		y_size=300;
		frame_nb=5;
		sigma_spot_move=3;
		return *this;
	}

	fluo_simulator & config_moving_spot(const char * bg_file){
		cus background(bg_file);
		background.slice(0);
		set_background(background);
		frame_nb=5;
		return *this;
	}
	fluo_simulator & config_scaling(const char * bg_file){
		cus background(bg_file);
		background.slice(0);
		set_background(background);
		frame_nb=10;
		scaling_matrix=cfl(2,2,1,1,
											 1.005,0.001,
											 0.001,1./1.005);
		frame_nb=5;
		return *this;
	}

	fluo_simulator & config_confocal_fluorescence(const char * bg_file){
		cus background(bg_file);
		background.slice(0);
		set_background(background);
		frame_nb=30;
		noise_type=POISSON_GAUSSIAN_NOISE;
		scaling_matrix=cfl(2,2,1,1,
											 1.005,-0.001,
											 0.001,1.005);
		stochastic_scale=true;
		SNR=5;
		sigma_spot_move=1.;
		bleach_rate=0.05;
		return *this;
	}

	// Simulation builder, main method.
	fluo_simulator & build_simulation(){
		// If background has not been set by the user, simulation
		// background is a 0 filled image of <x_size> and <y_size> is used
		// as background.
		if (background.is_empty()){
			background.assign(x_size,y_size).fill(0);
		}
		
		// Automatic threshold for spot support
		ushort cell_threshold=background.median()+sqrt(background.get_slice(0).variance());

		// border mask for spot drawing
		CImg<char> no_border_mask(background,"xy",1);
		cimg_for_borderXY(no_border_mask,x,y,spot_patch_size/2) no_border_mask(x,y)=0;
		
		// Spot list creation on the first frame (not drawn)
		CImg<char> spot_possible_locus(no_border_mask);
		spot_possible_locus.mul(background.get_threshold(cell_threshold,false,false));
		create_spots_list(spot_number,spot_possible_locus,background);  

		// For each frame redefine background (eg following scale_matrix),
		// draw and move spot.
		std::list<Simu_spot<ushort> >::iterator s=spot_list.begin();
		cus precedent_warping=background;
		CImg<float> img_center(1,2,1,1,background.width()/2,background.height()/2);
		cimg_for1(frame_nb,z){
			cus current=background;

			// Dealing with scale. Image is scaled then recentered.
			cfl global_scale, global_translation;
			if(!scaling_matrix.is_empty()||stochastic_scale){
				global_scale.assign(2,2,1,1,1,0,0,1);
				if(!scaling_matrix.is_empty()) global_scale*=scaling_matrix;
				cfl scaling_noise;
				if(stochastic_scale){
					scaling_noise.assign(2,2,1,1,
															 1+cimg::grand()*0.005,cimg::grand()*0.005,
															 cimg::grand()*0.005,1./(1.+cimg::grand()*0.005)
															 );
					global_scale*=scaling_noise;
				}
				global_translation=img_center-global_scale*img_center;
				precedent_warping=image_scaling(precedent_warping,global_scale,global_translation);
				current=precedent_warping;
			}

			// New support for spot after scaling
			CImg<char> available_mask=no_border_mask.get_mul(current.get_threshold(cell_threshold));
			int i;

			// For each spot, draw spot, move and define new intensity.
			for( s = spot_list.begin(),i=0; s != spot_list.end(); ++s,i++) {
				if(z<s->get_move_list().size()){
					ushort color = max(0.f,s->get_intensity(z)
														 - (float)current.linear_atXY(s->get_center()(0),s->get_center()(1)));
					s->draw_patch(&color);
					s->apply_patch(current);
					if(z!=frame_nb-1){
						if (global_scale.is_empty()){ 
							s->brownian_new_move(current,sigma_spot_move,available_mask);
						}else{
							s->brownian_new_move(current,sigma_spot_move,available_mask,
																	 global_scale,global_translation);
						}

						s->gaussian_new_int(sigma_noise_spot_intensity);
					}
				}
			}  

			// Bleach is imposed on the whole image 
			if (bleach_rate) current*=exp(-(float)z * bleach_rate);

			simulation.append(current,'z');

		}  
		if(SNR){
			switch(noise_type){
			case POISSON_NOISE: poisson_noise(simulation,SNR); break;
			case GAUSSIAN_NOISE: gaussian_noise(simulation,SNR); break;
			case POISSON_GAUSSIAN_NOISE: poisson_gaussian_noise(simulation,SNR); break;
			default:  cout <<"unknown noise type"<<endl;
			}
		}

		return *this;
	}
	

private:
	CImg<ushort> simulation;
	CImg<ushort> background;
	int x_size;
	int y_size;
	float bleach_rate;
	int frame_nb;
	list<Simu_spot<ushort> > 	spot_list;
	int spot_number;
	int spot_patch_size;
	float sigma_noise_spot_intensity;
	float sigma_spot_move;
	float spot_sigma;
	CImg<float> scaling_matrix;
	bool stochastic_scale;
	int noise_type;
	float SNR;

	
	cus image_scaling(cus & img, const cfl & scale,const cfl & shift_vector){
		cfl warping_vector(img.width(),img.height(),1,2);
		float *ptrx= warping_vector.data(0,0,0,0),*ptry=warping_vector.data(0,0,0,1);
		cfl invert_scale=scale.get_invert();
		cimg_forXY(warping_vector,x,y){
			*ptrx++=invert_scale(0,0)*x+invert_scale(1,0)*y;
			*ptry++=invert_scale(0,1)*x+invert_scale(1,1)*y;
		}
		img.warp(warping_vector,false,true,1);

		img.shift(shift_vector(0),shift_vector(1),0,0,1); 
		return img;
	}

	fluo_simulator & create_spots_list(int spot_number,CImg<char> & available_mask,cus & abackground){ 
		// Create spots object. Set starting point considering available mask (available_mask(x,y)==1).
		// set starting intensity given background.
		// srand ( time(NULL) );
		int x=0.;
		int y=0.;
		int dynamic=0;

		Simu_spot<ushort> sp;
		for (int i=0 ; i<spot_number ; i++){
			do{
				x=frand(0,available_mask.width());
				y=frand(0,available_mask.height());
			}while( !available_mask(x,y) )	;
			ushort color;
			if(abackground.max()!=abackground.min()){
				color = max(0,min((int)cimg::grand()*((USHRT_MAX   - abackground.max())/5)+ abackground.max(),USHRT_MAX));
			}else{
				color=(USHRT_MAX)/2;
			}
			
			dynamic=true;
			sp=Simu_spot<ushort>(x,y,spot_patch_size,dynamic,spot_sigma);
			sp.set_intensities(CImg<float>(1,1,1,1,color));
			spot_list.push_back(sp);
		}
		return *this;
	}

	fluo_simulator & gaussian_noise(cus & current , float SNR){

		// float noise_std=sqrt(current.variance()/pow((float)10,SNR/10));
		float noise_std=(current.mean())/(SNR);
		current.noise(noise_std,0); // gaussian noise
		return *this;
	}

	fluo_simulator & poisson_gaussian_noise(cus & current , float SNR){
		cfl current_fl=			current;
		current_fl*=SNR*SNR/current.mean();
		current_fl.noise(0,3);
		current_fl+=current_fl.get_fill(1000).noise(60,0);
		current=current_fl.get_cut(0,FLT_MAX);
		return *this;
	}
	fluo_simulator & poisson_noise(cus & current , float SNR){

		cfl current_fl=			current;
		current_fl*=SNR*SNR/current.mean();
		current_fl.noise(0,3);
		current=current_fl.get_cut(0,FLT_MAX); 
		return *this;
	}
public :

};


float process_phase(const float & stack_phase,const float & ref_phase, const float & tau_ref, const float & w){
	float phase_shift=ref_phase-atan(w*tau_ref);
	float lifetime=tan(stack_phase-phase_shift)/w;
	return lifetime;
}

void FD_FLIM_simu_motion(){
	float mv=2;
	float spot_sig=3;
	float sn=5;
	int size=100;
	int frame_nb=12;

	cfl phase_map(size,size,frame_nb); phase_map.fill(0);
	cfl total_non_modulated(size,size,frame_nb); total_non_modulated.fill(0);
	for(int sidx=0;sidx<sn;sidx++){
		fluo_simulator simu;
		simu.config_moving_spot_no_bg();
		simu.set_sigma_spot_move(mv);
		simu.set_spot_number(1);
		simu.set_frame_nb(frame_nb);
		simu.set_spot_sigma(spot_sig);
		// simu.set_SNR(snr);
		// simu.set_poisson_gaussian_noise();
		simu.set_size(size,size);
		simu.build_simulation();
		cfl non_modulated=simu.get_simulation();
		non_modulated.normalize(0,1);
		total_non_modulated.max(non_modulated);
		phase_map+=non_modulated.get_threshold(0.2)*(-1+0.4*(sidx));
	}
	total_non_modulated.save("non_modulated.tif");

	// Local max-based mixing
	cfl phase_modulated_bg=total_non_modulated.get_fill(10);
	cfl phase_modulated_spot=total_non_modulated.normalize(0,20);
	phase_modulated_bg.fill("i*(2+0.5*cos(z*2*pi/12))",false);
	phase_modulated_spot.fill("i*(2+0.5*cos(z*2*pi/12.+0.1))",false);
	cfl phase_modulated_max=phase_modulated_bg;
	phase_modulated_max.max(phase_modulated_spot);
	phase_modulated_max.save("simu-fdflim-max.tif");

	// Mask based mixing
	cfl total_non_modulated_off=total_non_modulated+10;
	phase_map.save("phase_map.tif");
	cfl phase_modulated_mask=total_non_modulated_off;
	cimg_forXYZ(phase_modulated_mask,x,y,z){
		phase_modulated_mask(x,y,z)*=(2+0.5*cos(z*2*cimg::PI/12+phase_map(x,y,z)));
	}
	phase_modulated_mask.save("simu-fdflim-mask.tif");

	cfl lf_map=phase_map.get_slice(0);
	float w=2*cimg::PI*0.04;	 // the Li-Flim software is usually set up to 40
	cimg_for(lf_map,ptr,float){
		*ptr=process_phase(*ptr,0,4.1,w);
	}
	lf_map.save("lf_gt.tif");
}






void test_fluo_simulator(const char * file_name)
// This function test <fluo_simulator> behavior, it also gives usage
// example
{
	{
		cout <<" basic no bg test" << endl;
		fluo_simulator simu;
		simu.config_moving_spot_no_bg();
		simu.build_simulation();
		simu.get_simulation().save("simu_no_bg.tiff");
	}

	{
		cout << " add <file_name> first slice as bg" << endl;
		fluo_simulator simu;
		simu.config_moving_spot(file_name);
		simu.build_simulation();
		simu.get_simulation().save("simu_bg.tiff");
	}

	{
		cout << " test warping no bleach" << endl;
		fluo_simulator simu;
		simu.config_scaling(file_name);
		simu.build_simulation();
		simu.get_simulation().save("simu_warp.tiff");
	}

	{
		cout << " test warping no bleach stochastic " << endl;
		fluo_simulator simu;
		cus background(file_name); background.slice(0);
		simu.set_background(background).set_frame_nb(10).set_stochastic_scale(true);
		simu.set_SNR(2);
		simu.set_poisson_noise();
		simu.set_sigma_spot_move(1.);
		simu.build_simulation();
		simu.get_simulation().save("simu_stochastic_warp.tiff");
	}

	{
		cout << " test warping no bleach stochastic + scaling tendancies " << endl;
		fluo_simulator simu;
		simu.config_confocal_fluorescence(file_name);
		simu.build_simulation();
		simu.get_simulation().save("simu_confocale_fluorescence.tiff");
	}
}

void testing_switch(int test_idx){
	switch(test_idx){
	case 1:
		{
			FD_FLIM_simu_motion();
		}

	}
}

#endif

