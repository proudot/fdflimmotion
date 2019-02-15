#include "CImg.h"
#include "CFigure.h"
#include "ICCD_noise_model.hpp"
#include "lifetime.hpp"

using namespace cimg_library;
using namespace std;

int main(int argc, char *argv[])
{
	int phase_nb=12;
	int expNumber=15;
	CImg<> photonCountTested=CImg<float>::sequence(expNumber,100,1000);
	photonCountTested.print("photonCount");
	CImg<float> fourierRMSE(1,expNumber,1,1);
	CImg<float> ICCDBalancedRMSE(1,expNumber,1,1);
	CImg<float> ICCDAbberBalancedRMSE(1,expNumber,1,1);
	
	for(int i=0;i<photonCountTested.size();i++){

		CImg<> simu_ref(200,200);
	float photoCount=photonCountTested(i);
	simu_ref=photoCount;
	CImgList<> list(phase_nb,simu_ref);
	simu_ref=list.get_append('z');
	cimg_forXYZ(simu_ref,x,y,z){
		simu_ref(x,y,z)*=(1+0.9*sin(z*2*cimg::PI/phase_nb+cimg::PI/2+1));
	}
	CImg<> simu_ref_gt=simu_ref;

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
	iccd_abber_noise_model->simulate_noise(simu_ref);
	//(simu_ref_gt,simu_ref).display();

	ICCD_noise_variance * iccd_noise_model=new ICCD_noise_variance;
	iccd_noise_model->set_parameters(g0_e,g0_sig,gccd,agauss_e,agauss_sig);


	phase_processor pp; 
	pp.assign(simu_ref.width(),simu_ref.height(),phase_nb);
	pp.noise_model=iccd_noise_model; // can be NULL
	CImg<> phase_model_iccd= pp.process_pixels(simu_ref,16,0);
	pp.noise_model=iccd_abber_noise_model; // can be NULL
	CImg<> phase_model_abber_iccd= pp.process_pixels(simu_ref,16,0);
	CImg<> phase_fourier= pp.process_pixels(simu_ref,1,0);
	// phase_model_iccd.print();
	// phase_fourier.print();
	
	const CImg<float> hist_iccd = phase_model_iccd.get_histogram(256,2,4);
	const CImg<float> hist_abber_iccd = phase_model_abber_iccd.get_histogram(256,2,4);
	const CImg<float> hist_fourier = phase_fourier.get_histogram(256,2,4);
	//hist_fourier.get_append(hist_iccd,'c').get_append(hist_abber_iccd,'c').display_graph(0,3);

	CImg<float> error_iccd = phase_model_iccd-cimg::PI/2-1;
	CImg<float> error_fourier = phase_fourier-cimg::PI/2-1;
	CImg<float> error_iccd_abber = phase_model_abber_iccd-cimg::PI/2-1;
	error_fourier.print("error fourier");
	
	// get rid of obviously diverged estimation
	// float mad_fourier=sqrt(error_fourier.variance(2));
	// float mad_iccd=sqrt(error_iccd.variance(2));
	// float obvious_err_thresh_iccd=error_iccd.median()+9*mad_iccd;
	// float obvious_err_thresh_fourier=error_fourier.median()+9*mad_fourier;
	// CImg<> error_fourier_inlier(error_fourier);
	// int count=0;
	// cimg_foroff(error_fourier,i){ if(error_fourier.get_abs()[i]<obvious_err_thresh_fourier){error_fourier_inlier[count++]=error_fourier[i];}};
	// error_fourier_inlier.assign(error_fourier_inlier.data(),count,1,1,1);
	// cout << "fourier outlier count :" << error_fourier.size()-count << endl;
	// CImg<> error_iccd_inlier(error_iccd);
	// count=0;
	// cimg_foroff(error_iccd,i){ if(error_iccd.get_abs()[i]<obvious_err_thresh_iccd){error_iccd_inlier[count++]=error_iccd[i];}};
	// error_iccd_inlier.assign(error_iccd_inlier.data(),count,1,1,1);
	// cout << "iccd outlier count :" << error_iccd.size()-count << endl;

	fourierRMSE(i)=sqrt(error_fourier.sqr().mean());
	ICCDBalancedRMSE(i)=sqrt(error_iccd.sqr().mean());
	ICCDAbberBalancedRMSE(i)=sqrt(error_iccd_abber.sqr().mean());
	cout << "fourier: " << fourierRMSE(i) << endl;
	cout << "iccd balanced: " << ICCDBalancedRMSE(i) << endl;
	cout << "iccd abber balanced: " << ICCDAbberBalancedRMSE(i) << endl;
	delete iccd_noise_model;
	delete iccd_abber_noise_model;
	}
	// (simu_ref_gt,simu_ref).display();

	// ICCD_abber_noise_model noise_model;
	// noise_model.set_parameters(g0_e,g0_sig,gccd,agauss_e,agauss_sig,x_0,y_0,sigma_x,sigma_y,o);
	// noise_model.simulate_noise(simu_ref);
	photonCountTested.print("photonCountTested");
	fourierRMSE.print("fourier rmse");
	(photonCountTested/10).get_append(fourierRMSE,'x').save_dlm("build/ICCDfourierRMSE.dlm");
	(photonCountTested/10).get_append(ICCDBalancedRMSE,'x').save_dlm("build/ICCDweightedRMSE.dlm");
	(photonCountTested/10).get_append(ICCDAbberBalancedRMSE,'x').save_dlm("build/ICCDAbberWeightedRMSE.dlm");
	CFigure fig(600,400);
	fig.plot(photonCountTested,fourierRMSE,"r-");
	fig.plot(photonCountTested,ICCDBalancedRMSE,"b-");
	fig.plot(photonCountTested,ICCDAbberBalancedRMSE,"g-");
	fig.set_axis();
	fig.display();
	

	
	
	return 0;
}
