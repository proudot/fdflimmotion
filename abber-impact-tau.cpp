#include "CImg.h"
#include "CFigure.h"
#include "ICCD_noise_model.hpp"
#include "lifetime.hpp"

using namespace cimg_library;
using namespace std;

int main(int argc, char *argv[])
{
	int phase_nb=12;

	CImg<> simu_ref(200,200);
	float thousand=1000;
	simu_ref.draw_gaussian(simu_ref.width()/2,simu_ref.height()/2,50,&thousand);
	CImgList<> list(phase_nb,simu_ref);
	simu_ref=list.get_append('z');
	cimg_forXYZ(simu_ref,x,y,z){
		simu_ref(x,y,z)*=1+0.9*sin(z*2*cimg::PI/phase_nb+cimg::PI/2);
	}

	CImg<> simu_ref_gt=simu_ref;

	float 
		g0_e=10,
		g0_sig=2, 
		gccd=2, 
		agauss_e=1000, 
		agauss_sig=50;
	
	float
		x_0=300,
		y_0=300,
		sigma_x = 1.457419306943089339e+02,
		sigma_y = 1.409734354363447721e+02,
		o = 2.403725089884512223e-01;

	ICCD_noise_variance * iccd_noise_model=new ICCD_noise_variance;
	iccd_noise_model->set_parameters(g0_e,g0_sig,gccd,agauss_e,agauss_sig);
	iccd_noise_model->simulate_noise(simu_ref);

	phase_processor pp; 
	pp.assign(simu_ref.width(),simu_ref.height(),phase_nb);
	pp.noise_model=iccd_noise_model; // can be NULL
	CImg<> phase_model_iccd= pp.process_pixels(simu_ref,16,0);
	CImg<> phase_fourier= pp.process_pixels(simu_ref,1,0);
	phase_model_iccd.print();
	phase_fourier.print();
	
	const CImg<float> hist_iccd = phase_model_iccd.histogram(256,1,2);
	const CImg<float> hist_fourier = phase_fourier.histogram(256,1,2);
	hist_fourier.get_append(hist_iccd,'c').display_graph(0,3);
	
	// (simu_ref_gt,simu_ref).display();

	// ICCD_abber_noise_model noise_model;
	// noise_model.set_parameters(g0_e,g0_sig,gccd,agauss_e,agauss_sig,x_0,y_0,sigma_x,sigma_y,o);
	// noise_model.simulate_noise(simu_ref);
	

	delete iccd_noise_model;
	
	return 0;
}
