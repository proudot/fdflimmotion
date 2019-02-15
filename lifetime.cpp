#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include "math.h"
#include <algorithm>
#include <time.h>

//My include
#include "lifetime.hpp"

using namespace cimg_library;
using namespace std;

template <class T>
static float estimateNoiseVar(CImg<T> & stack){
	CImg<float> value(stack.depth());
	cimg_forX(value,z){
		value(z)=stack.get_slice(z).noise_variance(3);
	}
	return value.max();
}

CImg<int> gauss_outlier_labeling(cfl & img,float mean,float std){
	cfl label = img;
	label-=mean;
	label.abs();
	label/=std;
	CImg<int> res = label;
	return res;
}

cfl cum_diff(const cfl & stack){
	cfl res(stack.get_slice(0),"xy",0);
	cfl tmp;
	for (int i = 1; i < stack.depth(); ++i) {
		tmp=stack.get_slice(0);
		tmp-=stack.get_slice(i);
		res+=tmp.get_abs();
	}
	return res;
}

cfl renorm_from_first(const cd & stack){
	float std0=sqrt(stack.get_slice(0).variance());
	float m0=(stack.get_slice(0).mean());
	float std = 0;
	float m = 0;
	cfl res(stack.get_slice(0));
	for (int i = 1; i < stack.depth(); ++i) {
		std=sqrt(stack.get_slice(i).variance());
		m=(stack.get_slice(i).mean());
		res.append((stack.get_slice(i)-m)*std0/std+m0,'z');
	}
	return res;
}

CImg<float> strip_zeros(CImg<float> & res){
	float non_null_ptr[res.size()];
	unsigned int count=0;
	unsigned int i = 0;
	while(i<res.size()){
		if(res[i]!=0){
			non_null_ptr[count++]=res[i];
		}
		i++;
	}
	CImg<float> non_null(non_null_ptr,count);
	return non_null;
}

int main(int argc,char** argv){

	// Handling options
	const char * file_name  = cimg_option("-tiff",(char*)NULL,"input stack (MANDATORY)");
	const char * ref_file_name  = cimg_option("-ref",(char*)NULL,"input ref stack (MANDATORY if no -ref_phase)");
	const char * ref_phase_file_name  = cimg_option("-ref_phase",(char*)NULL,"input ref phase stack (MANDATORY if no -ref)");
	const char * analysis_mask_name  = cimg_option("-plot_mask",(char*)NULL,"mask used for fit plotting");
	const char * mask_file_name= cimg_option("-m",(char*)NULL,"processing mask");
	const char * noise_calibr_stack_filename  = cimg_option("-noise",ref_file_name, "input noise stack for calibration (reference is taken if empty)");
	const int phase_nb = cimg_option("-phase_nb",12,"Number of signal use for phase modulation");
	int total_frame_count = cimg_option("-frame_nb",0,"Number of frame to be considered in the stack.");
	const int ref_phase_nb = cimg_option("-ref_phase_nb",12,"Number of signal use for phase modulation on reference");
	const int ROI = cimg_option("-ROI",1,"Select a rectangle ROI to process (1 for enabling).");
	const float threshold  = cimg_option("-threshold",0,"Only pixels above threshold are processed.");
	const int type  = cimg_option("-type",1,"Phase retrieval: 0: least square. 1: Fourier. 2: irls. 4: assym irls.");
	const int irls  = cimg_option("-irls",0,".");
	const int print_fit_stats = cimg_option("-presidue",0,"Print residue of sinusoidal fit");
	// const int robust_stats  = cimg_option("-rstats",0,"Compute robust mean and variance");
	// const float nsigma  = cimg_option("-nsigma",2,"Choosen sigma for outlier detection in robust stats.");

	if(file_name==NULL) return(0);

	// Reading reference and calibration stack (supposed of reasonable size)
	CImg<img_type> ref_stack;
	CImg<img_type> noise_calibr_stack;
	if(ref_file_name!=NULL){
		ref_stack.assign(ref_file_name);
		ref_stack.slices(0,ref_phase_nb-1);	}
	if(noise_calibr_stack_filename!=ref_file_name){
		noise_calibr_stack.assign(noise_calibr_stack_filename);
		// noise_calibr_stack.slices(0,ref_phase_nb-1);
	}else {
		noise_calibr_stack=ref_stack;
	}

	// Noise from reference or specific calibration stack
	INoise_variance * noise_model_ptr=NULL;
	CFigure noise_fit_fig;
	if((type==5)||(type==9)||(type==10)||(type==15)){
			cfl cov_noise;
			ICCD_noise_variance * iccd_noise_model=new ICCD_noise_variance;
			iccd_noise_model->ccd_model.b_size=20;
			iccd_noise_model->ccd_model.robust=true;
			iccd_noise_model->estimate_param(noise_calibr_stack);
			cov_noise=iccd_noise_model->ccd_model.estimator->ieng_robust_cov_matrix();
			CFigure fig(400,300);
			iccd_noise_model->ccd_model.estimator->build_1d_graph(fig);

			int font_size=12;
			fig.xlabel("E(I(x))",font_size);
			fig.ylabel("var(I(x))",font_size);
			fig.set_axis_font_size(font_size);
			// fig.set_axis(1500,18960,0,300000);

			// Computation using ICCD noise model (computed off line and hardcoded) to compare both aproach.
			ICCD_abber_noise_model * iccd_abber_noise_model=new ICCD_abber_noise_model;
			float a=1.936317977185310595e-04,
				b=-7.142794355502983528e-02,
				c=1.207309731908114918e+04,
				x_0=2.479598488795075468e+02,
				y_0=3.716442639085927340e+02,
				sigma_x = 1.457419306943089339e+02,
				sigma_y = 1.409734354363447721e+02,
				o = 2.403725089884512223e-01;

			iccd_abber_noise_model->set_parameters(a,b,c,x_0,y_0,sigma_x,sigma_y,o);


			// Computation using CCD noise model to compare both aproach.
			CCD_noise_variance * ccd_noise_model=new CCD_noise_variance;
			ccd_noise_model->b_size=20;
			ccd_noise_model->scheme=3;
			ccd_noise_model->robust=true;
			ccd_noise_model->estimate_param(ref_stack);
			cov_noise=ccd_noise_model->estimator->ieng_robust_cov_matrix();
			ccd_noise_model->estimator->plot_fit(fig,"g-");

			fig.set_axis(0,ccd_noise_model->estimator->get_X().max(),
									 0,ccd_noise_model->estimator->get_Y().max());
			fig.erase();
			fig.replot();
			fig.save("quadratic_noise_fit.png");


			// save .dat for gnuplot
			iccd_noise_model->ccd_model.estimator->get_X()
				.get_append(iccd_noise_model->ccd_model.estimator->get_Y(),'x').save_dlm("quadratic_fit.dat");

			iccd_noise_model->print();
			noise_model_ptr=iccd_noise_model;
			noise_fit_fig=fig;
			delete ccd_noise_model;

		}

	if((type==14)){
		// Computation using ICCD noise model confocale (computed off line and hardcoded) to compare both aproach.
		ICCD_abber_noise_variance * iccd_abber_noise_model=new ICCD_abber_noise_variance;

		float
			a=			 4.660194036013196945e-05,
			b=			 9.417017447162129073e+00,
			c=			-4.890463402256602421e+03,
			x_0=		 2.924863798666465868e+02,
			y_0=		 3.857255041914487492e+02,
			sigma_x= 1.000716942674229699e+02,
			sigma_y= 9.730745973387014658e+01,
			o =      1.697514865139863882e-01;

		iccd_abber_noise_model->set_parameters(a,b,c,x_0,y_0,sigma_x,sigma_y,o);
		iccd_abber_noise_model->print();
		noise_model_ptr=iccd_abber_noise_model;
	}


	if((type==17)){
		// Computation using ICCD noise model widefield (computed off line and hardcoded) to compare both aproach.
		ICCD_abber_noise_variance * iccd_abber_noise_model=new ICCD_abber_noise_variance;

		float
			a=1.936317977185310595e-04,
			b=-7.142794355502983528e-02,
			c=1.207309731908114918e+04,
			x_0=2.479598488795075468e+02,
			y_0=3.716442639085927340e+02,
			sigma_x = 1.457419306943089339e+02,
			sigma_y = 1.409734354363447721e+02,
			o = 2.403725089884512223e-01;

		iccd_abber_noise_model->set_parameters(a,b,c,x_0,y_0,sigma_x,sigma_y,o);
		iccd_abber_noise_model->print();
		noise_model_ptr=iccd_abber_noise_model;
	}
	if((type==16)){
		// Computation using ICCD noise model widefield (computed off line and hardcoded) to compare both aproach.
		ICCD_abber_noise_variance * iccd_abber_noise_model=new ICCD_abber_noise_variance;

		float
			a=1.936317977185310595e-04,
			b=-7.142794355502983528e-02,
			c=1.207309731908114918e+04,
			x_0=2.479598488795075468e+02,
			y_0=3.716442639085927340e+02,
			sigma_x = 1.457419306943089339e+02,
			sigma_y = 1.409734354363447721e+02,
			o = 2.403725089884512223e-01;

		iccd_abber_noise_model->set_parameters(a,b,c,x_0,y_0,sigma_x,sigma_y,o);
		iccd_abber_noise_model->print();
		noise_model_ptr=iccd_abber_noise_model;
	}



	// cropping
	cus first_image; first_image.load_tiff(file_name,0,1);
	CImg<img_type> ROI_coord;
	CImg<img_type> cropped_first_image(first_image);
	if(ROI){
		if (ROI==2){
			ROI_coord.assign(5);
			ROI_coord(0)=271 ;
			ROI_coord(1)= 273;
			ROI_coord(3)= 280;
			ROI_coord(4)= 287;
		}else{
			ROI_coord=first_image.get_slice(0).get_select("select a ROI.");
		}

		if(ref_file_name!=NULL) ref_stack.crop(ROI_coord(0),ROI_coord(1), 0,
																					 ROI_coord(3),ROI_coord(4),ref_stack.depth()-1);
		cropped_first_image.crop(ROI_coord(0),ROI_coord(1),ROI_coord(3),ROI_coord(4));
	}

	if (!ref_stack.is_empty()) ref_stack.save("crop_ref.tif");

	// Pixel to compute mask
	cfl cropped_first_image_float(cropped_first_image);
	cfl first_image_float(first_image);
	phase_processor pp; pp.assign(cropped_first_image.width(),cropped_first_image.height(),phase_nb);
	cfl mask_of_interest=cropped_first_image_float.get_threshold(threshold).get_dilate(30,30);
	if(mask_file_name){
		CImg<> explicit_mask(mask_file_name);
		mask_of_interest.mul(explicit_mask);
	}
	pp.pixel_to_compute.assign(mask_of_interest);



	if((type == 6)||(type==7)||(type==11)){
		cout << ":: noise variance computation using OLS residual::" << endl;
		OLS_Noise * noise_model_array=new OLS_Noise;
		noise_model_array->compute_residues(ref_stack,pp.pixel_to_compute);
		noise_model_array->compute_variance();
		noise_model_ptr=noise_model_array;
	}

	if(type == 8){
		cout << ":: noise variance computation using sample variance ::" << endl;
		Sample_Noise * noise_model_array=new Sample_Noise();
		noise_model_array->compute_empirical_variance(ref_stack,pp.pixel_to_compute);
		noise_model_ptr=noise_model_array;
	}


	pp.noise_model=noise_model_ptr; // can be NULL

	cfl ref_phases;
	if((ref_file_name != NULL)&&((ref_phase_file_name == NULL)))  {
		ref_phases = pp.process_pixels(ref_stack,type,0);
		ref_phases.save("ref_phases.tiff");
	}

	if(ref_phase_file_name){
		ref_phases.assign(ref_phase_file_name);
	}

	cfl sequence_lifetimes;
	cfl sequence_lifetimes_crop;
	cfl sequence_phases;
	cfl sequence_sub_stack;
	cus sub_orig_stack;
	cus sub_stack;

	// counting frames number if not specified
	TIFF *tif = TIFFOpen(file_name,"r");
	if (total_frame_count==0) do ++total_frame_count; while (TIFFReadDirectory(tif));
	TIFFClose(tif);

	for(int frame_count=0;
			frame_count<=total_frame_count-phase_nb;
			frame_count+=phase_nb){

		// read and crop substack
		sub_orig_stack.load_tiff(file_name,frame_count,frame_count+phase_nb-1);
		if(ROI){
			sub_stack=sub_orig_stack.get_crop(ROI_coord(0),ROI_coord(1),0
																				,ROI_coord(3),ROI_coord(4),sub_orig_stack.depth()-1);
		}else{
			sub_stack=sub_orig_stack;
		}

		// optly estimate noise variance
		float err_sig=0;
		if( (type==2)||(type==3)||(type==4))  {
			err_sig  = estimateNoiseVar(sub_stack);
			cout << "estimated noise variance : "  << err_sig << endl;
		}

		pp.cov_type=0;

		// optly estimate noise variance
		if((type == 6)||(type==7)){
			cout << ":: noise variance computation using OLS residual::" << endl;
			delete pp.noise_model;
			OLS_Noise * noise_model_array=new OLS_Noise;
			noise_model_array->compute_residues(sub_stack,pp.pixel_to_compute);
			noise_model_array->compute_variance();
			pp.noise_model=noise_model_array;
		}

		// optly estimate noise variance
		if(type == 8){
			cout << ":: noise variance computation using sample variance ::" << endl;
			delete pp.noise_model;
			Sample_Noise * noise_model_array=new Sample_Noise();
			noise_model_array->compute_empirical_variance(sub_stack,pp.pixel_to_compute);
			pp.noise_model=noise_model_array;
		}

		// Compute phase on main stack
		pp.random_fit=false; 				// plot fitting on random pixels for diagnosis
		if(analysis_mask_name) pp.pixel_to_analyse.assign(analysis_mask_name).slice(0); // plot fitting on a precise mask
		cfl stack_phases=pp.process_pixels(sub_stack,type,err_sig,irls);


		// Compute lifetime
		lifetime_processor lp;
		cfl lifetimes;
		if((ref_file_name != NULL)||(ref_phase_file_name != NULL)){
			lifetimes = lp.process_phase(stack_phases,ref_phases,4.1);}
		else{lifetimes = lp.process_phase(stack_phases); }

		if ((type == 5)||(type==6)) delete pp.noise_model;

		// Print stats on the estimator
		if((type!=1)&&(print_fit_stats)){
			cout << "mean cov :: " << endl;
			print_val(pp.cov_matrix);
			cout << ":: printing fit :: " << endl;
			CFigure residue_fig;
			CImg<float> x(pp.residues.at(0).size()); x.sequence(0,phase_nb-1); x.transpose();
			cfl X,R;
			R=pp.residues.get_append('x');
			X.assign(pp.residues.size(),pp.residues.at(0).size());
			cimg_forY(X,l) X.get_shared_line(l).fill(l);
			cimg_forY(R,l){ cout << "std dev residue " << l
													 << " : " << sqrt(R.get_shared_line(l).variance())
													 << " mean : " << R.get_shared_line(l).mean()
													 << endl;
			}

			// for(unsigned int i =0 ; i< pp.residues.size(); i++){
			// 	CImg<float> r = pp.residues.at(i);
			// 	X.append(x,'x');
			// 	R.append(r,'x');
			// }

			X.get_vector().append(R.get_vector(),'x').save_dlm("residue.dat");
			residue_fig.set_axis(X,R);
			residue_fig.plot(X,R,"bx");
			residue_fig.save("residue.png");
		}


		// Save various output at substack level
		{
			// unprocessed lifetime map
			sequence_lifetimes_crop.append(lifetimes,'z');


			// unprocessed lifetime map in spatial context
			if(ROI){
				CImg<float> lf;
				CImg<float> sp;
				lf=sub_orig_stack.get_slice(0); lf.fill(NO_VALUE);//////////
				sp=lf;
				lf.draw_image(ROI_coord(0),ROI_coord(1),lifetimes);
				sp.draw_image(ROI_coord(0),ROI_coord(1),stack_phases);
				sequence_phases.append(sp,'z');
				sequence_lifetimes.append(lf,'z');
				sequence_sub_stack.append(sub_stack,'z');


				lf.save("lifetime.tiff");
			}else
				{
					sequence_phases.append(stack_phases,'z');
					sequence_lifetimes.append(lifetimes,'z');
				};
		}
		lp.non_null.get_vector().save_dlm("lifetime.dat");
			// save row histogram
			int bins=100;
			int bin_min=1;
			int bin_max=4;
			CImg<float> raw_hist;
			raw_hist.assign(bins).get_sequence(bin_min,bin_max).get_transpose()
				.append(lp.non_null.get_histogram(bins,bin_min,bin_max).get_transpose(),'x')
				.save_dlm("lifetime_hist.dat");
		// lp.non_null.histogram(100).display_graph(0,3,1,0,lp.non_null.min(),lp.non_null.max());
		}

			// robustly labeled pixels ( 0 = inlier, increasing value mean increasing outlier probability)
			// CImg<int> label = gauss_outlier_labeling(res,estimated_mean,MAD);
			// label.save("label.tiff");
		// Save various output at sequence level.
		{
			sequence_lifetimes.save("lifetime.tiff");
			sequence_lifetimes_crop.save("lifetime_crop.tiff");
			sequence_phases.save("phase.tiff");
			if(ROI){
				sequence_sub_stack.save("cropped_stack.tiff");
			}
			// non-robustly labeled pixels ( 0 = inlier, increasing value mean increasing outlier probability)
			// CImg<int> label_no_estimate = gauss_outlier_labeling(res,mean,sigma);
			// label_no_estimate.save("label_no_estimate.tiff");

			// label difference
			// (label-label_no_estimate).save("label_diff.tiff");

			// robust 3 sigmas outlier outlining with inlier values (0 = outlier)
			// CImg<short > bad_lf_mask=  (1-label.get_threshold(nsigma));
			// res.get_mul(bad_lf_mask).save("robust_inlier.tiff");

			if(ref_file_name!= NULL) ref_stack.save("cropped_ref.tiff");

			// cout << "lifetime_processing" << endl;
			// Lifetime post-processing

			// store non-null lifetime
			CImg<float> non_null = strip_zeros(sequence_lifetimes_crop);

			// estimate variance and mean
			float mean = non_null.mean();
			float sigma =  sqrt(non_null.variance());
			float sigmad =  sqrt(non_null.variance(2));

			// Some shiny output
			cout << "mean :" << mean ;
			cout << " , median  : " << non_null.median();
			cout << " , sigma :" << sigma ;
			cout << " , sigmad :" << sigmad ;
			cout << " , var :" << sigma*sigma << endl;

		// 	// non-robustly labeled pixels ( 0 = inlier, increasing value mean increasing outlier probability)
		// 	CImg<int> label_no_estimate = gauss_outlier_labeling(sequence_lifetimes,mean,sigma);
		// 	label_no_estimate.save("label_no_estimate.tiff");

		// 	// normal 3 sigma  outlier outlining with inlier values (0 = outlier)
		// 	CImg<short> bad_lf_mask=  (1-label_no_estimate.get_threshold(nsigma));
		// 	(sequence_lifetimes.get_mul(bad_lf_mask)).save("normal_inlier.tiff");

		// // estimate variance and mean robustly
		// if(robust_stats){
		// 	mean_estimator me;
		// 	CImg<float> r ;
		// 	float MAD = 1.4826*(non_null - non_null.median()).get_abs().median();
		// 	CImg<float> beta(1,1,1,1,non_null.median());
		// 	// here maxiter/threshold 100/0.01 instead of 10/0.1 ... not tested
		// 	me.init_residue(non_null,beta);
		// 	Leclerc cf_leclerc(MAD,2);
		// 	me.set_cost_function(&cf_leclerc);
		// 	me.ils_no_init(r);
		// 	float estimated_mean = r(0);

		// 	// Some shiny output
		// 	cout << "estimated mean : " << estimated_mean << endl;
		// 	cout << "estimated sigma :" << MAD 	<< endl;
		// }


		}

	// }
	// else{
	// 	// Here we do not work on the whole frame but on a list of spot defining trajectories and
	// 	// feature size
	// 	list<Spot<img_type> > spots=restore_spots<Spot<img_type> >(spot_file_name);
	// 	Lf_switch<img_type> lf_switch;
	// 	std::list<Spot<img_type> >::iterator s=spots.begin();
	// 	for ( s = spots.begin(); s != spots.end(); ++s) {
	// 		Int_estimator<img_type>::intensities_from_move(&(*s),stack,1);
	// 		lf_switch.run(&(*s),stack,type);
	// 	}

		// if( type == 1) compute_lifetime_fourier(feature_l,2.42342,4.1,false);

		// 	switch(type){
		// 	case 0:
		// 		compute_lifetime(spots,2.42342,4.1,false);
		// 		break;
		// 	case 1:
		// 		compute_lifetime_fourier(spots,2.45956,4.1,false);
		// 		break;
		// 	case 2:
		// 		{
		// 		varEstimator<img_type> ve;
		// 		float err_sig = ve.estimate(stack);   // float err_sig = 450;
		// 		compute_lifetime_irls(spots,2.42342,4.1, err_sig ,3,&dcost_leclerc,false);
		// 		break;
		// 		}
		// 	default:
		// 		cout << " no fit " << endl;
		// 		break;
		// 	}
		// 	CImg<float> lifetimes = get_array(spots,&Spot<img_type>::get_lifetime);
		// 	cout << "var lifetime " << lifetimes.variance() << endl;
		// 	cout << "mean lifetime " << lifetimes.mean() << endl;
		// 	cout << "min lifetime " << lifetimes.min() << endl;
		// 	cout << "max lifetime " << lifetimes.max() << endl;

	// 	save_spots(spots,"results.spot");


	// 	if(truth !=NULL){
	// 		// Compaison with a truth spot file
	// 		list<True_Spot<img_type> > true_spots=restore_spots<True_Spot<img_type> >(truth);
	// 		match_spots<Spot<img_type> ,True_Spot<img_type> >(spots,true_spots);
	// 		erase_undetected<img_type>(true_spots);
	// 		print_match<img_type>(true_spots,"results.txt");
	// 	}
	// }
	return 0;


}
