#include "CImg.h"
#include "CFigure.h"
#include "ICCD_noise_model.hpp"
#include "lifetime.hpp"

using namespace cimg_library;
using namespace std;

int main(int argc, char *argv[])
{
	int phase_nb=12;
	int expNumber=1;
	float photonCount=1000;
	
	float 
		g0_e=10,
		g0_sig=2, 
		gccd=2, 
		agauss_e=0, 
		agauss_sig=50;
	
	float
		x_0=100,
		y_0=100,
		sigma_x = 50,
		sigma_y = 50,
		o = 0.;


	CImg<> phaseSample=CImg<float>::sequence(50,0,2*cimg::PI);
	CImg<> ampSample=CImg<float>::sequence(50,0,3.);
	CImg<> CSample=CImg<float>::sequence(50,g0_e*photonCount*0.1,g0_e*photonCount*1.9);


	float trueAmp=0.9;

	int x=100,y=100;
	
	CImg<> energy(phaseSample.size(),ampSample.size(),CSample.size());
	

	CImg<> simu_ref(200,200);
	simu_ref=photonCount;
	CImgList<> list(phase_nb,simu_ref);
	simu_ref=list.get_append('z');
	cimg_forXYZ(simu_ref,x,y,z){
		simu_ref(x,y,z)*=(1+trueAmp*sin(z*2*cimg::PI/phase_nb+cimg::PI/2+1));
	}
	CImg<> simu_ref_gt=simu_ref;


	ICCD_abber_noise_variance * iccd_abber_noise_model=new ICCD_abber_noise_variance;
	iccd_abber_noise_model->set_parameters(g0_e,g0_sig,gccd,agauss_e,agauss_sig,
																				 x_0, y_0, sigma_x, sigma_y , o);
	iccd_abber_noise_model->simulate_noise(simu_ref);
	//(simu_ref_gt,simu_ref).display();


	Sin_fit_fact fitter;
	fitter.measure_nb=phase_nb;
	fitter.set_threshold(0.001);
	fitter.set_max_iter(100);
		
	CImg<float> Y(1,1,simu_ref.depth(),1);
	CImg<float> sig_weigth(simu_ref.depth());
	cimg_forZ(simu_ref,z){
		Y(z)=simu_ref(x,y,z);
		sig_weigth(z)=1./sqrt(iccd_abber_noise_model->get_variance(Y(z),x,y));
	};
	fitter.set_residue_weigth(sig_weigth);

	fitter.assign(Y,0);
	fitter.set_starting_value();
	//fitter.display_fit();
	fitter.print_state(cout);

	cimg_foroff(phaseSample,pi){ cimg_foroff(ampSample,ai) {cimg_foroff(CSample,ci){
				fitter.set_phase(phaseSample[pi]);fitter.set_amp(ampSample[ai]);fitter.set_C(CSample[ci]);
				energy(pi,ai,ci)=fitter.get_residue().magnitude();
			}}}
	//energy.display();
	
	
	energy.get_unroll('x').save("energy.dlm");

// fitter.run();
// if(irls){
// 	float MAD = sqrt(fitter.get_residue().variance(2)) ;
// 	fitter.get_cost_function()->set_sig_noise(2*MAD);
// 	fitter.irls_no_init();
// }

	
return 0;
}
