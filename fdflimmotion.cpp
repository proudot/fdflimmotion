// Purpose : Challenge ISBI 2012

#include <iostream>
#include <fstream>
#include <list>
#include "math.h"
#include <algorithm>
#include <time.h>

//My include
#include "feature.hpp"
//#include "tracker.hpp"
#include "feature_high.hpp"
#include "CImg.h"
#include "feature_detector.hpp"
#include "lf_estimator.hpp"
#include "Matcher.hpp"
#include "evo_visualizer.hpp"
#include "lifetime.hpp"

using namespace cimg_library;
using namespace std;

// memory
// #define  ERROR_TH				0.001
// #define  MAX_ITER				100

#define  ERROR_TH				0.001
#define  MAX_ITER				10
#define  RESIDUE_MAX		1000
#define FEATURE_DIST	  10
#define WINDOW_SIZE			7
#define NB_FEATURE_MAX	30

int main(int argc,char** argv){

	const char * file_name  = cimg_option("-tiff",(char*)NULL,"input file");
	const char * ref_file_name  = cimg_option("-ref",(char*)NULL,"input reference stack");
	const char * ref_phase_file_name  = cimg_option("-ref_phase",(char*)NULL,"input phase of the reference stack");
	const char * mask_file_name= cimg_option("-m",(char*)NULL,"processing mask");
	const float ref_lifetime=cimg_option("-ref_lf",4.1,"fluorescence lifetime of the reference sample");
	const int plain_lf_type=cimg_option("-plainlf",1,"lifetime estimation algorithm for background (default: Fourier)");
	const char * noise_calibr_stack_filename  = cimg_option("-noise",ref_file_name, "input noise stack for calibration (reference is taken if empty)");
	const int nb_f_max = cimg_option("-feature_nb",NB_FEATURE_MAX,"nb feature max");
	const int frame_nb = cimg_option("-frame_nb",0,"nb of frame to process");
	const float psf_std = cimg_option("-psf_std",1.,"psd standard deviation");
	const int ref_frame_nb = cimg_option("-ref_frame_nb",frame_nb,"nb of frame to process");
  const int option  = cimg_option("-option",17,"an option for various testing purpose");
  const int lf_type  = cimg_option("-lf",13,"method for inloop lifetime retrieval");
  const int irls   = cimg_option("-irls",1,"irls or not.");
  const int window_size   = cimg_option("-ws",WINDOW_SIZE,"window_size");
  const int alpha_robust   = cimg_option("-alpha",18,"window_size");
	const float threshold  = cimg_option("-threshold",0,"Only pixels above threshold are processed.");
	const int inpainting = cimg_option("-inpainting",0,"Activate or deactivate inpainting");

	typedef float img_type;

	CImg<img_type> stack; stack.load_tiff(file_name);
	if(frame_nb>0) stack.slices(0,frame_nb-1);

	cfl ref_phases,ref_stack,noise_calibr_stack;
	if((ref_file_name != NULL)&&((ref_phase_file_name == NULL)))  {
		ref_stack.assign(ref_file_name).slices(0,ref_frame_nb-1);
	}


	if(noise_calibr_stack_filename!=ref_file_name){
		noise_calibr_stack.assign(noise_calibr_stack_filename);
	}else {
		noise_calibr_stack=ref_stack;
	}

	phase_processor pp; pp.assign(stack.width(),stack.height(),frame_nb);
	cfl float_first_frame=stack.get_slice(0);
	cfl mask_of_interest=float_first_frame.get_threshold(threshold).get_dilate(30,30);
	if(mask_file_name){
		CImg<> explicit_mask(mask_file_name);
		mask_of_interest.mul(explicit_mask);				
	}

	pp.pixel_to_compute.assign(mask_of_interest);

	// Use different noise model if specified by plain_lf_type
	INoise_variance * noise_model_ptr=NULL;
	CFigure noise_fit_fig;
	if((plain_lf_type==5)||(plain_lf_type==9)||(plain_lf_type==10)||(plain_lf_type==15)){
		cfl cov_noise;
		ICCD_noise_variance * iccd_noise_model=new ICCD_noise_variance;
		iccd_noise_model->ccd_model.b_size=20;
		iccd_noise_model->ccd_model.robust=true;
		iccd_noise_model->estimate_param(noise_calibr_stack);
		cov_noise=iccd_noise_model->ccd_model.estimator->ieng_robust_cov_matrix();
		CFigure fig(400,300);
		iccd_noise_model->ccd_model.estimator->build_1d_graph(fig);

		int font_size=12;
		fig.xlabel("E(I(x))",font_size).ylabel("var(I(x))",font_size).set_axis_font_size(font_size).save("quadratic_noise_fit.png");
		iccd_noise_model->print();
		noise_model_ptr=iccd_noise_model;
	}

	if((plain_lf_type==14)){
		// Computation using ICCD noise model (computed off line and hardcoded) to compare both aproach.
		ICCD_abber_noise_variance * iccd_abber_noise_model=new ICCD_abber_noise_variance;
		float a=	 4.660194036013196945e-05, b=			 9.417017447162129073e+00, c=			 -4.890463402256602421e+03,
			x_0=		 2.924863798666465868e+02, y_0=		 3.857255041914487492e+02,
			sigma_x =1.000716942674229699e+02, sigma_y =9.730745973387014658e+01,
			o =      1.697514865139863882e-01;
		iccd_abber_noise_model->set_parameters(a,b,c,x_0,y_0,sigma_x,sigma_y,o);
		iccd_abber_noise_model->print();
		noise_model_ptr=iccd_abber_noise_model;
	}
	pp.noise_model=noise_model_ptr; // can be NULL


	// Feature detection set up and first frame detection
	cout << ":: Feature detection :: " << endl;
	list<Feature<img_type> > feature_l; // global feature list
	feature_detector fd; fd.nb_feature_max=nb_f_max; fd.feature_dist=FEATURE_DIST; fd.feature_size=window_size;
	feat_highlight fh; fh.print_lost=false; fh.scale=3;
	feature_l=fd.feature_detection(stack.get_slice(0),psf_std);

	// Frame to frame tracking class setup
	int matcher_type;
	Match_model_ni_I<img_type> * matcher;
	matcher = new Match_model_ni_I<img_type>();
	matcher_type=TEMPLATE_MATCHER;
	matcher->set_threshold(ERROR_TH);
	matcher->set_max_iter(MAX_ITER);
	Leclerc cost_f; cost_f.set_alpha(alpha_robust); matcher->set_cost_function(&cost_f); // \alpha = 10 dans track_no_model

	// First tracking round using simple gaussian fitting
	cout << ":: First motion guess using gaussian fit ::"<< endl;
	CImg<img_type> current,next;
	gaussian_fit<img_type> * gf= new gaussian_fit<img_type>();
	gf->set_max_iter(10); gf->set_threshold(0.001);
	typedef std::list<Feature<img_type> >::iterator fit; fit s;
	for(int frame_idx=0;frame_idx<stack.depth();frame_idx++){
		CImg<img_type> current=stack.get_slice(frame_idx);
		gf->set_current(&current);
		for (fit s=	feature_l.begin(); s != feature_l.end(); ++s) {
			if(s->get_lost()) continue;
			gf->set_feature(&(*s));
			gf->init();
			gf->run();
			if(!gf->matching_diagnosis(RESIDUE_MAX,false)){
				if(frame_idx<2) s->set_lost(true);
			}else{ s->set_frame_move(); }
		}	// features
	} // frames
	CImg<unsigned char> fh_stack=fh.build_stack_color(stack,feature_l);
	fh.draw_tracks_on_img(fh_stack,feature_l).save_tiff("gf_tracks.tiff");


	// Iterative estimation of lifetime and motion
	cout << ":: Lifetime and motion joint estimation :: " << endl;
	// Variable watcher...
	int max_iter=8;
	Evo_visualizer scale_evo(stack.depth()*(nb_f_max+10));
	Evo_visualizer scale_model_evo(stack.depth()*(nb_f_max+10));
	cfl betas_ls;
	cfl betas_irls;
	list<cfl> f_betas_ls;
	list<cfl> f_betas_irls;
	Evo_visualizer * phase_evo[feature_l.size()]; cimg_for1(feature_l.size(),i) phase_evo[i] = new Evo_visualizer();
	{
		int end_frame_idx=stack.depth();

		int nb_lost=0;
		Lf_switch<img_type> lf_switch;
		int feature_idx=0;
		for (fit s=	feature_l.begin(); s != feature_l.end(); ++s,feature_idx++) {
			if(s->get_lost()) continue;
			cout << "feature ID " << feature_idx <<" : " <<flush;
			matcher->set_feature(&(*s));
			float old_phase;
			int nb_iter=0;
			do{
				if(s->get_lost()) break;
				old_phase=s->get_phase();
				char filename[40]; cimg::number_filename("fit_iter.dat",nb_iter,2,filename);
				lf_switch.run(&(*s),stack,lf_type,filename);
				cimg::number_filename("fit_iter.png",feature_idx,2,filename);
				cimg::number_filename(filename,nb_iter,2,filename);
				lf_switch.fig.save(filename);


				phase_evo[feature_idx]->add_data(s->get_phase());

				// Feature move list backup and erase.
				CImg<float> tmp0=s->get_move(0);
				s->erase_move_list();
				s->set_center(tmp0);
				s->set_new_move(s->get_center());

				// track over end_frame_idx frames
				for(int frame_idx=1;frame_idx<end_frame_idx;frame_idx++){
					if(s->get_lost()) break;
					current=stack.get_slice(frame_idx);
					cfl blurred_current=current.get_blur(2.f);
					matcher->set_current(&current);
					((Match_smodel<img_type> *)matcher)->set_frame_idx(frame_idx);

					matcher->init();
					matcher->run();
					// char filename[40]; cimg::number_filename("beta_feature.png",feature_idx,2,filename);
					// matcher->beta_evo.plot(); matcher->beta_evo.save(filename);

					if(irls){
						float scale_parameter=sqrt(matcher->get_residue().variance(2)); // MAD
						cost_f.set_sig_noise(scale_parameter);
						matcher->irls_no_init();
					}
					if(!matcher->matching_diagnosis(RESIDUE_MAX,false)){
						s->set_lost(true);
						cout << "lost on frame " << frame_idx << ", ";
					}else{
						s->set_frame_move();
					}
				} //frames

				nb_iter++;
				//}while((nb_iter<max_iter));
			}while((std::abs((s->get_phase()-old_phase)/old_phase) > 0.000000001)&&(nb_iter<max_iter));
			// keep the number of feature to nb_feature_max
			cout << "nb iter : " << nb_iter;
			cout << "phase : :" << s->get_phase();
			cout << endl;
		} // features

	}
	delete matcher ;
	delete gf;

	cout << ":: Feature Recording ::" << endl;
	// Save final spot
	save_spots(feature_l,"results.spot",false);
	challenge_write_spot_XML(feature_l,"results.xml");

	CImg<> stack_norm=stack;
	cimg_forZ(stack_norm,z) stack_norm.get_shared_plane(z).normalize(0,1000);
	fh_stack=fh.build_stack_color(stack_norm,feature_l);
	fh.draw_tracks_on_img(fh_stack,feature_l).save_tiff("tracks.tiff");

	// Printing phase evolution (DISABLED)
	// CFigure phase_evo_global;
	// for (unsigned int i = 0; i < feature_l.size(); ++i)
	// 	{
	// 		if(phase_evo[i]->data_idx>1){
	// 			phase_evo[i]->plot(phase_evo_global);
	// 			phase_evo[i]->plot();
	// 			char filename[40]; cimg::number_filename("phase_evo_f.png",i,2,filename);
	// 			phase_evo[i]->save(filename);
	// 			//phase_evo[i]->fig.save_data(strcat(filename,".dat"));
	// 		}
	// 	}
	// phase_evo_global.erase().set_axis().replot().save("phase_evo_global.png");


	if(ref_phase_file_name){
		ref_phases.assign(ref_phase_file_name);
	}

	if((ref_file_name != NULL)&&((ref_phase_file_name == NULL)))  {
		ref_phases = pp.process_pixels(ref_stack,1); // Warning testing with Fourier for speed !!!
		ref_phases.save("ref_phases.tiff");
	}

	// Final lifetime computation
	if(!ref_phases.is_empty()){
		lifetime_processor lp;
		float Pi=3.141592;
		float w=2*Pi*0.04;	 // the Li-Flim software is usually set up to 40
		CImg<> spot_lifetime(feature_l.size(),2);
		int i=0;
		cout << "size : " << feature_l.size() << endl;
		for (fit s=	feature_l.begin(); s != feature_l.end(); ++s,i++) {
			if(s->get_lost()) continue;
			cout << "spot : " << i << endl;
			float lifetime = lp.process_phase(s->get_phase(),ref_phases(s->get_x(),s->get_y()),ref_lifetime,w);
			s->set_lifetime(lifetime);
			spot_lifetime(i,0)=s->get_move_list().size();
			spot_lifetime(i,1)=lifetime;
		}
		spot_lifetime.save("spot_lifetime.dlm");
	}

	// Lifetime recontruction
	if(plain_lf_type){
		cout << ":: plain lifetime reconstruction ::" << endl;
		cfl stack_phases=pp.process_pixels(stack,plain_lf_type);
		stack_phases.save("plain_phase.tif");
		lifetime_processor lp;
		CImg<> plain_lf = lp.process_phase(stack_phases,ref_phases,ref_lifetime,mask_of_interest);
		plain_lf.save("plain_lf.tif");
		CImg<> bg_lf(plain_lf);
		cfl plain_lf_inpainted;
		cfl stack_phases_inpainted;
		if(inpainting){
			cout << ":: inpainting ::" << endl;
			plain_lf_inpainted=lifetime_repair(plain_lf,2,30,mask_of_interest);
			stack_phases_inpainted=lifetime_repair(stack_phases,2,30,mask_of_interest);
			plain_lf.mul(mask_of_interest); 
		}
		cout << ":: applying patches ::" << endl;
		float zero=0;
		for (fit s=	feature_l.begin(); s != feature_l.end(); ++s) {
			float lifetime=s->get_lifetime();
			float phase=s->get_phase();
			if(s->get_move_list().size()>6) plain_lf.draw_circle(s->get_move(0)(0),s->get_move(0)(1),5,&lifetime);
			if(s->get_move_list().size()>6) stack_phases.draw_circle(s->get_move(0)(0),s->get_move(0)(1),5,&phase);
			if(inpainting){
				if(s->get_move_list().size()>6) plain_lf_inpainted.draw_circle(s->get_move(0)(0),s->get_move(0)(1),5,&lifetime);
				if(s->get_move_list().size()>6) stack_phases_inpainted.draw_circle(s->get_move(0)(0),s->get_move(0)(1),5,&phase);
			}

			if(s->get_move_list().size()>1) bg_lf.draw_circle(s->get_move(0)(0),s->get_move(0)(1), 5,&zero);
		}
		plain_lf.print(); bg_lf.print();
		if(inpainting){
		plain_lf_inpainted.save("rebuild_lifetime_inpainted.tiff");
		stack_phases_inpainted.save("rebuild_phase_inpainted.tiff");
		}
		plain_lf.save("rebuild_lifetime.tiff");
		stack_phases.save("rebuild_phase.tiff");
		bg_lf.save("bg_lf.tiff");
		float *bg_lf_ptr=bg_lf.data(), *ptr=bg_lf.data(); cimg_foroff(bg_lf,off){if(*ptr){*bg_lf_ptr++=*ptr;} ptr++;}
		bg_lf.assign(bg_lf.data(),bg_lf_ptr-bg_lf.data()).save("bg_lf.dlm");

		float *plain_lf_ptr=plain_lf.data(); ptr=plain_lf.data(); cimg_foroff(plain_lf,off){if(*ptr){*plain_lf_ptr++=*ptr;} ptr++;}
		plain_lf.assign(plain_lf.data(),plain_lf_ptr-plain_lf.data()).save("plain_lf.dlm");
	}

		cimg_for1(feature_l.size(),i) delete phase_evo[i];
}

