#include "CImg.h"
#include "CFigure.h"
#include "ICCD_noise_model.hpp"
#include "lifetime.hpp"
#include "fluo_simulator.hpp"
#include "sin_fit_fact.hpp"
using namespace cimg_library;
using namespace std;

int main(int argc, char *argv[])
{
	int phase_nb=12;
	int expNumber=1;
	CImg<> photonCountTested=CImg<float>::sequence(expNumber,1000,1000);
	photonCountTested.print("photonCount");
	CImg<float> fourierRMSE(1,expNumber,1,1);
	//CImg<float> ICCDBalancedRMSE(1,expNumber,1,1);
	CImg<float> ICCDAbberBalancedRMSE(1,expNumber,1,1);
	
	float 
		g0_e=10,
		g0_sig=2, 
		gccd=2, 
		agauss_e=1000, 
		agauss_sig=50;
	
	float
		x_0=100,
		y_0=100,
		sigma_x = 50,
		sigma_y = 50,
		o = 0.;

	ICCD_abber_noise_variance * iccd_abber_noise_model=new ICCD_abber_noise_variance;
	iccd_abber_noise_model->set_parameters(g0_e,g0_sig,gccd,agauss_e,agauss_sig,
																				 x_0, y_0, sigma_x, sigma_y , o);


	// Motion displacement
	float mv=2;
	float spot_sig=3;
	float sn=1;
	int size=200;
	int frame_nb=12;
	float phaseBG=1;
	float phaseSpot=cimg::PI/2+1;

	Sin_fit_fact fitter;
	fitter.measure_nb=phase_nb;
	fitter.set_threshold(0.001);
	fitter.set_max_iter(100);

	cfl Y(1,1,frame_nb,1);
	cfl sig_weigth(Y);
	for(int i=0;i<photonCountTested.size();i++){
		float photonCount=photonCountTested(i);
		cfl phaseFourier(sn);
		cfl phaseAbber(sn);
		for(int sidx=0;sidx<sn;sidx++){
			cfl phase_map(size,size,frame_nb); phase_map.fill(0);
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
			list<Simu_spot<ushort> > spotl=simu.get_spot_list();
			Simu_spot<ushort> spot=spotl.front();
			cfl non_modulated=simu.get_simulation();
			non_modulated.normalize(0,1);
			phase_map=non_modulated.get_threshold(0.2)*phaseSpot;
			phase_map.display();
			cfl non_modulated_off=non_modulated*(1+photonCount/2);
			non_modulated_off.display();
			
			cfl phase_modulated_mask=non_modulated_off;
			cimg_forXYZ(phase_modulated_mask,x,y,z){
				phase_modulated_mask(x,y,z)*=(1+0.9*sin(z*2*cimg::PI/frame_nb+phase_map(x,y,z)));
			}
			phase_modulated_mask.display();
			iccd_abber_noise_model->simulate_noise(phase_modulated_mask);
			spot.set_intensities_from_move(phase_modulated_mask,0);
			spot.print_all(cout);
			cimg_foroff(Y,z){
				Y(z)=spot.get_intensity(z);
				sig_weigth(z)=1./sqrt(iccd_abber_noise_model->get_variance(Y(z),spot.get_move(z)(0),spot.get_move(z)(1)));
			}
			fitter.assign(Y,0);
			fitter.set_starting_value();
			fitter.run();
			fitter.display_fit();
			phaseFourier(sidx)=fitter.get_phase();
			//fitter.set_residue_weigth(sig_weigth);
			fitter.set_starting_value();
			fitter.run();
			phaseAbber(sidx)=fitter.get_phase();
		}
		fourierRMSE(i)=sqrt((phaseFourier-phaseSpot).sqr().mean());
		ICCDAbberBalancedRMSE(i)=sqrt((phaseAbber-phaseSpot).sqr().mean());

		// Local max-based mixing
		// cfl phase_modulated_bg=total_non_modulated.get_fill(10);
		// cfl phase_modulated_spot=total_non_modulated.normalize(0,20);
		// phase_modulated_bg.fill("i*(2+0.5*cos(z*2*pi/12))",false);
		// phase_modulated_spot.fill("i*(2+0.5*cos(z*2*pi/12.+0.1))",false);
		// cfl phase_modulated_max=phase_modulated_bg;
		// phase_modulated_max.max(phase_modulated_spot);
		// phase_modulated_max.save("simu-fdflim-max.tif");

		// // Mask based mixing
		// cfl total_non_modulated_off=total_non_modulated+10;
		// phase_map.save("phase_map.tif");
		// cfl phase_modulated_mask=total_non_modulated_off;
		// cimg_forXYZ(phase_modulated_mask,x,y,z){
		// 	phase_modulated_mask(x,y,z)*=(2+0.5*cos(z*2*cimg::PI/12+phase_map(x,y,z)));
		// }
		// phase_modulated_mask.save("simu-fdflim-mask.tif");

	cout << "fourier: " << fourierRMSE(i) << endl;
	cout << "iccd abber balanced: " << ICCDAbberBalancedRMSE(i) << endl;
	}
	delete iccd_abber_noise_model;

	// noise_model.simulate_noise(simu_ref);
	photonCountTested.print("photonCountTested");
	fourierRMSE.print("fourier rmse");
	(photonCountTested).get_append(fourierRMSE,'x').save_dlm("ICCDfourierRMSE.dlm");
	(photonCountTested).get_append(ICCDAbberBalancedRMSE,'x').save_dlm("ICCDAbberWeightedRMSE.dlm");
	CFigure fig(600,400);
	fig.plot(photonCountTested,fourierRMSE,"r-");
	fig.plot(photonCountTested,ICCDAbberBalancedRMSE,"g-");
	fig.set_axis();
	fig.display();
	
	return 0;
}
