#ifndef _LIFETIME_H_
#define _LIFETIME_H_

#include "CImg.h"
#include "inpaint.hpp"
#include "sin_fit_fact.hpp"
#include "bg_stack_estimator.hpp"
#include "sin_fit.hpp"
#include "ICCD_noise_model.hpp"
using namespace cimg_library;

float NO_VALUE=-10;


typedef unsigned short img_type;
typedef CImg<float> cfl;
typedef CImg<img_type> cd;
typedef CImg<img_type> cus;

static float get_phase_asym_var(CImg<float> & values,CImg<float> & weights,Sin_fit & fitter,Asym_cos_function * acf){
	static CImg<float> beta;
	fitter.assign(values,0);
	fitter.set_residue_weigth(weights);
	fitter.set_starting_value();
	fitter.run(beta);
	float MAD = sqrt(fitter.get_residue().variance(2)) ;
	acf->set_sig_pos(0.1*MAD);
	acf->set_sig_neg(3*MAD);			//
	fitter.irls_no_init();
	return fitter.get_phase();
}

static float get_phase_ls_var(CImg<float> & values,CImg<float> & weights,
															Sin_fit & fitter,float sig = 0,int irls=0){

	fitter.assign(values,0);
	fitter.set_residue_weigth(weights);
	fitter.set_starting_value();
	fitter.run();
	if(irls){
		float MAD = sqrt(fitter.get_residue().variance(2)) ;
		fitter.get_cost_function()->set_sig_noise(2*MAD);
		fitter.irls_no_init();
	}

	return fitter.get_phase();
}

static float get_phase_ls(CImg<float> & values,Sin_fit & fitter){
	static CImg<float> beta;
	fitter.assign(values,0);
	fitter.set_starting_value();
	fitter.run(beta);
	return fitter.get_phase();
}

static float get_phase_irls(CImg<float> & values, float sig,Sin_fit & fitter){
	static CImg<float> beta;
	fitter.assign(values,0);
	fitter.set_starting_value();
	fitter.run(beta);
	CImg<float> r = fitter.get_residue();
	float MAD = sqrt(fitter.get_residue().variance(2)) ;
	fitter.get_cost_function()->set_sig_noise(2*MAD);
	fitter.irls_no_init();

	return fitter.get_phase();
}

void starting_value_bg(Sin_fit & fitter,CImg<float> &val){
	fitter.set_C(val.median());
	fitter.set_amp(val.mean());
}

float get_phase_fourier(CImg<float> & values){
	int nb_frame = values.size();
	CImg<float> Cos(nb_frame);
  CImg<float> Sin(nb_frame);
	float phase;
	cimg_forX(Cos,k){
		Cos(k)=cos(2*cimg::PI*(k)/(nb_frame));
		Sin(k)=sin(2*cimg::PI*(k)/(nb_frame));
	}
	float fcos=0, fsin=0;
	cimg_forX(values ,x){
		fcos+=values(x)*Cos(x);
		fsin+=values(x)*Sin(x);
	}
	if (fcos){ phase=-atan(fsin/fcos)+cimg::PI/2;}
	else {phase=0;}

	return phase;
}

cfl emp= CImg<float>::empty();


struct INoise_variance{ virtual float get_variance(float x,float coordx=0, float coordy=0)=0; };

struct ICCD_noise_variance: public INoise_variance, public ICCD_noise_model{
	float get_variance(float x,float coordx=0,float coordy=0 ){return ICCD_noise_model::get_variance(x);}};

struct CCD_noise_variance: public INoise_variance, public CCD_noise_model{
	float get_variance(float x,float coordx=0,float coordy=0 ){return CCD_noise_model::get_variance(x);}};

struct ICCD_abber_noise_variance: public INoise_variance, public ICCD_abber_noise_model{
	float get_variance(float I,float coordx=0,float coordy=0 ){
		return ICCD_abber_noise_model::get_variance(I,coordx,coordy);}};


struct OLS_Noise:public INoise_variance{
	CImg<float> residues; 				// Nx12 residuals for N samples
	cfl variances;

	template<class T>
	void compute_residues(const CImg<T> &stack,  CImg<float> pixel_to_compute=CImg<float>::empty()){
		if (pixel_to_compute.is_empty()) pixel_to_compute.assign(stack,"xy",1);
		Sin_fit_fact fitter;
		cimg_forXY(stack,x,y){
			if(pixel_to_compute(x,y)){
				CImg<float> Y(1,1,stack.depth(),1);
				cimg_forZ(Y,z){Y(z)=stack(x,y,z);};
				fitter.assign(Y,0);
				fitter.set_starting_value();
				fitter.run();
				residues.append(fitter.get_residue(),'x');
			}
		}
	}

	void compute_variance(){
		variances.assign(residues.height());
		cimg_forY(residues,y){
			variances(y)=residues.get_shared_line(y).variance(2);
		}
	}

	float get_variance(float z,float x,float y ){return variances(z);};
};

struct Sample_Noise:public INoise_variance{
	CImg<float> residues; 				// Nx12 residuals for N samples
	cfl variances;

	template<class T>
	void compute_empirical_variance(const CImg<T> &stack,
																	CImg<float> pixel_to_compute=CImg<float>::empty()){
		variances.assign(stack.depth());
		if (pixel_to_compute.is_empty()) pixel_to_compute.assign(stack,"xy",1);
		CImg<float> Y(1,1,stack.depth(),1);
		cimg_forZ(stack,z){ cimg_forXY(stack,x,y){ if(pixel_to_compute(x,y)){
					variances(z)=stack.get_shared_plane(z).variance(3);
				}}}
	}
	float get_variance(float z,float x , float y ){return variances(z);};
};

struct phase_processor{
	INoise_variance * noise_model;
	CImg<float> cov_matrix;
	CImgList<float> cov_matrix_list;
	int cov_type;					// Covariance matrix type :
	// 0 : empirical
	// 1 : model based (linear square)
	CImgList<float> residues;
	CImg<float> phases;
	CImg<float> pixel_to_compute;
	CImg<float> pixel_to_analyse;
	int phase_nb;
	float random_fit;							// choose some random point to plot estimator.

	phase_processor(cus & stack){
		assign(stack);
	}

	phase_processor(){}

	void assign(cus & stack){
		assign(stack.width(),stack.height(),stack.depth());
	}
	void assign(int w, int h, int pnb){
		noise_model=NULL;
		cov_matrix.assign(3,3).fill(0);
		cov_type=0;
		pixel_to_compute.assign(w,h).fill(1);
		pixel_to_analyse.assign(w,h).fill(0);
		phase_nb=pnb;
		random_fit=false;
	}

	template <class T>
	CImg<float> process_pixels(CImg<T> & stack, int type=0,float sig=0,int irls=0){
		// float phase_shift=phase_ref-atan(w*tau_ref);

		float phase;

		int nb_div=0;
		phases.assign(stack.width(),stack.height());

		Sin_fit_fact fitter;
		fitter.measure_nb=phase_nb;
		fitter.set_threshold(0.001);
		fitter.set_max_iter(100);

		// dealing with cost function
		Cost_function * cf=NULL;
		if((irls)||(type==2)){
			cf = new Leclerc(0,1);
			fitter.set_cost_function(cf);
		}
		Asym_cos_function * acf=NULL;
		switch(type){
		case 13:
		case 10:
		case 9:
		case 7:
			{
				acf= new Leclerc_bis_asym(0,0,1);
				fitter.set_cost_function(acf);

				break;
			}
		}

		CImg<float> betas; 						// beta collection for covariance type 0

		cfl blurred_stack;
		if (type==5) {
			blurred_stack.assign(stack);
			cimg_forZ(blurred_stack,z) blurred_stack.get_shared_plane(z).blur(2.f);
		}
		int count=0;

		int plot_sample_nb=100;
		int sample_distance=phases.size()/plot_sample_nb;
		CFigure plot_all_fit;

		cimg_forXY(stack,x,y){
			if(pixel_to_compute(x,y)){
				switch(type){
				case 0:
					{
						CImg<float> Y(1,1,stack.depth(),1);
						cimg_forZ(Y,z){Y(z)=stack(x,y,z);};
						phase=get_phase_ls(Y,fitter);
						break;
					}
				case 1:
					{
						CImg<float> Y(stack.depth());
						cimg_forX(Y,z){
							Y(z)=stack(x,y,z);
						}
						phase=get_phase_fourier(Y);
						break;
					}
				case 2:
					{
						CImg<float> Y(1,1,stack.depth(),1);
						cimg_forZ(Y,z){Y(z)=stack(x,y,z);};
						phase=get_phase_irls(Y,sig,fitter);
						break;
					}
				case 3:
					{
						int region_size=5;
						Bg_estim bg_estim(sig);
						CImg<float> region;
						region=stack.get_crop(x-region_size/2,y-region_size/2,
																	x+region_size/2,y+region_size/2);
						phase=bg_estim.get_phase_bg(region,sig,fitter);
						break;
					}
				case 4:
					{
						CImg<float> Y(1,1,stack.depth(),1); cimg_forZ(Y,z){Y(z)=stack(x,y,z);};
						Bg_estim bg_estim(sig);
						bg_estim.alpha=3;
						bg_estim.alpha_neg=6;
						phase=bg_estim.get_phase_bg(Y,sig,fitter);
						break;
					}
				case 14:
				case 16:
				case 5: // Weights depend on intensity value
					{
						if (noise_model==NULL) cout << "ERROR : A noise model must be plug in phase processor."<<endl;
						CImg<float> Y(1,1,stack.depth(),1);
						CImg<float> sig_weigth(stack.depth());
						cimg_forZ(Y,z){
							Y(z)=stack(x,y,z);
							sig_weigth(z)=1./sqrt(noise_model[0].get_variance(Y(z),x,y));
						};

						phase=get_phase_ls_var(Y,sig_weigth,fitter,sig,irls);
						break;
					}
				case 6:
				case 8: // Weights depend on phase modulation
					{
						if (noise_model==NULL) cout << "ERROR on noise_model"<<endl;
						CImg<float> Y(1,1,stack.depth(),1);
						CImg<float> sig_weigth(stack.depth());
						cimg_forZ(Y,z){
							Y(z)=stack(x,y,z);
							sig_weigth(z)=1./sqrt(noise_model[0].get_variance(z));
						};

						phase=get_phase_ls_var(Y,sig_weigth,fitter,sig,irls);
						break;
					}
				case 7: // Weights assym depend on phase modulation (not sigma but variance)
					{
						if (noise_model==NULL) cout << "ERROR on noise_model"<<endl;
						CImg<float> Y(1,1,stack.depth(),1);
						CImg<float> sig_weigth(stack.depth());
						cimg_forZ(Y,z){
							Y(z)=stack(x,y,z);
							sig_weigth(z)=1./(noise_model[0].get_variance(z));
						};

						phase=get_phase_asym_var(Y,sig_weigth,fitter,acf);
						break;
					}
				case 9: // weighted + assymetric
					{
						if (noise_model==NULL) cout << "error on noise_model"<<endl;
						CImg<float> Y(1,1,stack.depth(),1);
						CImg<float> sig_weigth(stack.depth());
						cimg_forZ(Y,z){
							Y(z)=stack(x,y,z);
							sig_weigth(z)=1./(noise_model[0].get_variance(Y(z)));
						};

						phase=get_phase_asym_var(Y,sig_weigth,fitter,acf);
						break;
					}
				case 10: // local variance + assymetric
					{
						if (noise_model==NULL) cout << "error on noise_model"<<endl;
						int region_size=5;
						CImg<float> region=stack.get_crop(x-region_size/2,y-region_size/2,
																							x+region_size/2,y+region_size/2);
						CImg<float> sig_weigth(stack.depth()*region_size*region_size);
						cimg_foroff(region,i){
							sig_weigth[i]=1./sqrt(noise_model[0].get_variance(region[i]));
						};

						phase=get_phase_asym_var(region,sig_weigth,fitter,acf);
						break;
					}
				case 15:
				case 11: // Weights using  variance instead of std.
					{
						if (noise_model==NULL) cout << "ERROR on noise_model"<<endl;
						CImg<float> Y(1,1,stack.depth(),1);
						CImg<float> sig_weigth(stack.depth());
						cimg_forZ(Y,z){
							Y(z)=stack(x,y,z);
							sig_weigth(z)=1./(noise_model[0].get_variance(z));
						};

						phase=get_phase_ls_var(Y,sig_weigth,fitter);
						break;
					}
				case 12: // Weights using local pseudo residuals
					{
						static int ws=7;
						CImg<float> Y(1,1,stack.depth(),1);
						CImg<float> sig_weigth(stack.depth());
						cimg_forZ(Y,z){
							Y(z)=stack(x,y,z);
							sig_weigth(z)=1./sqrt(stack.get_crop(x-ws/2,y-ws/2,z,x+ws/2,y+ws/2,z,true)
																		.get_pseudo_residuals().variance(3));
						};
						phase=get_phase_ls_var(Y,sig_weigth,fitter,sig,irls);
						break;
					}
				case 13: // Weights using local pseudo residuals on assymetric results.
					{
						static int ws=7;
						CImg<float> Y(1,1,stack.depth(),1);
						CImg<float> sig_weigth(stack.depth());
						cimg_forZ(Y,z){
							Y(z)=stack(x,y,z);
							sig_weigth(z)=1./sqrt(stack.get_crop(x-ws/2,y-ws/2,z,x+ws/2,y+ws/2,z,true)
																		.get_pseudo_residuals().variance(3));
						};
						phase=get_phase_asym_var(Y,sig_weigth,fitter,acf);
						break;
					}
				}


				if(phase != phase) {		// Nan detection trick
					cout << "Nan phase processed" << endl;
					phase=NO_VALUE;
					nb_div++;
				}

				phases(x,y)=phase;

				// Analyse fit using plot and print if needed.
				if(((count%sample_distance ==0)&&(random_fit)||(pixel_to_analyse(x,y)))&&(type!=1)){
					char filename[40]; cimg::number_filename("plot_sin_fit.png",x,4,filename);
					cimg::number_filename(filename,y,4,filename);
					char wfilename[40]; cimg::number_filename("weighted_sin_fit.png",x,4,wfilename);
					cimg::number_filename(wfilename,y,4,wfilename);
					if(irls) cout << "pixel : ("
												<< x << "," << y
												<< "). noise std : " <<  fitter.get_cost_function()->sig_noise
												<< endl;
					fitter.build_1d_graph().save(filename);
					fitter.build_1d_weighted_graph().save(wfilename);
					if(irls) if ((x==1)&&(y==7)){
							fitter.get_cost_function()->get_cost_graph().save("cost_at_pixel.png");
							fitter.get_cost_function()->get_rho_graph().save("rho_at_pixel.png");
						}
					if (irls) if ((x==0)&&(y==14)){
							fitter.ieng_robust_cov_matrix().save_ascii("cov_0_14.txt");
						}
					if (irls) if ((x==0)&&(y==0)){
							fitter.ieng_robust_cov_matrix().save_ascii("cov_0_0.txt");
						}

					// fitter.build_1d_graph(plot_all_fit);
 				}

				if(type!=1){
					residues.push_back(fitter.get_residue());
					switch(cov_type){
					case 0:
						betas.append(fitter.get_beta(),'x');
						break;
					case 1:
						cov_matrix+=fitter.estimate_cov_matrix();
						break;
					case 2:
						cfl cov = fitter.ieng_robust_cov_matrix();
						cov_matrix+=cov;
						cov_matrix_list.push_back(cov);
						break;
					}
				}

				count++;
			}
			else{
				phases(x,y)=NO_VALUE;
			}
		}

		if(type!=1) switch(cov_type){
			case 0:
				cimg_forY(betas,i){
					betas.get_shared_line(i)-=betas.get_shared_line(i).mean();
				}
				cov_matrix=betas*betas.get_transpose()/count;
				break;
			case 1:
				cov_matrix/=count;
				break;
			}

		cout << "nb diverge : " <<  nb_div << endl;

		plot_all_fit.erase().set_axis().replot();
		plot_all_fit.save("plot_sin_fit_global.png");

		delete acf;
		delete cf;

		return phases;
	}
};


struct lifetime_processor{
	CImg<float> non_null;
	// Computing lifetime using a reference phase map and and a reference lifetime.

	CImg<float> process_phase(CImg<float> & stack_phase, CImg<float> & ref_phase,float tau_ref,CImg<> & mask=CImg<>::empty()){
		CImg<float> lifetime(stack_phase,"xy",0);
		float w=2*cimg::PI*0.04;	 // the Li-Flim software is usually set up to 40 Mhz
		non_null.assign(stack_phase,"xyzc");
		float * non_null_ptr=non_null.data();
		float * lifetime_ptr=lifetime.data();
    float * ref_phase_ptr=ref_phase.data();
    float * stack_phase_ptr=stack_phase.data();
		if(mask.is_empty()) mask.assign(stack_phase).fill(1);
    float * mask_ptr=mask.data();

		cimg_forXY(stack_phase,x,y){
			if((*stack_phase_ptr!=NO_VALUE)&&(*mask_ptr)){
				*lifetime_ptr=process_phase(*stack_phase_ptr,*ref_phase_ptr,tau_ref,w);
				*non_null_ptr=*lifetime_ptr; non_null_ptr++;
			}
			else{ *lifetime_ptr=0;}
			lifetime_ptr++;
			ref_phase_ptr++;
			stack_phase_ptr++;
			mask_ptr++;
		}
		non_null.assign(non_null.data(),non_null_ptr-non_null.data());

		return lifetime;
	}

	// Computing lifetime using a reference phase and and a reference lifetime.
	CImg<float> process_phase(CImg<float> & stack_phase,float ref_phase=2.41794, float tau_ref=4.1){
		CImg<float> lifetime(stack_phase,"xy",0);
		float w=2*cimg::PI*0.04;	 // the Li-Flim software is usually set up to 40
		non_null.assign(stack_phase,"xyzc");
		float * non_null_ptr=non_null.data();
		float * lifetime_ptr=lifetime.data();
    float * stack_phase_ptr=stack_phase.data();

		cimg_forXY(stack_phase,x,y){
			if(*stack_phase_ptr!=NO_VALUE){
				*lifetime_ptr=process_phase(*stack_phase_ptr,ref_phase,tau_ref,w);
				*non_null_ptr=*lifetime_ptr; non_null_ptr++;
			}
			else{ *lifetime_ptr=NO_VALUE;}
			lifetime_ptr++;
			stack_phase_ptr++;
		}
		non_null.assign(non_null.data(),non_null_ptr-non_null.data());

		return lifetime;
	}

	float process_phase(const float & stack_phase,const float & ref_phase, const float & tau_ref, const float & w){
		float phase_shift=ref_phase-atan(w*tau_ref);
		float lifetime=tan(stack_phase-phase_shift)/w;
		return lifetime;
	}

};

//! local variance
template<typename T>
CImg<T> blur_local_variance(const CImg<T> & img, int n, int variance_method=1){
	CImg<T> res(img.width(),img.height(),img.depth(),img.spectrum());
	T *ptrd = res._data;
	const int hl = n/2, hr = hl - 1 + n%2;
	cimg_forXYZC(res,x,y,z,c) { // 3d
		const int
			x0 = x - hl, y0 = y - hl, z0 = z-hl, x1 = x + hr, y1 = y + hr, z1 = z+hr,
			nx0 = x0<0?0:x0, ny0 = y0<0?0:y0, nz0 = z0<0?0:z0,
			nx1 = x1>=img.width()?img.width()-1:x1, ny1 = y1>=img.height()?img.height()-1:y1,
		  nz1 = z1>=img.depth()?img.depth()-1:z1;
		*(ptrd++) = img.get_crop(nx0,ny0,nz0,c,nx1,ny1,nz1,c).variance(variance_method);
	}
	return res;
}


CImg<unsigned char> get_lifetime_outliers_mask_mad( CImg<> & lf_map,float scale = 3, int patch_size=30){
	CImg<> MAD_map =  sqrt(blur_local_variance(lf_map,patch_size,2));
	CImg<> med=lf_map.get_blur_median(patch_size);
	CImg<unsigned char> mask(lf_map,"xyzc",0);
	mask=lf_map.operator-(lf_map.get_min(med+(scale*MAD_map))).get_threshold(0,false,true);
	mask+=lf_map.get_max(med-(scale*MAD_map)).operator-(lf_map).get_threshold(0,false,true);
	mask.threshold(0,false,true);
	return mask;
}

CImg<unsigned char> get_lifetime_outliers_mask_pseudo_mad( CImg<> & lf_map,float scale = 3, int patch_size=30){
	CImg<> med=lf_map.get_blur_median(patch_size);
	CImg<> MAD_map =  1.4826*(lf_map-med).get_abs().get_blur_median(patch_size);
	CImg<unsigned char> mask(lf_map,"xyzc",0);
	mask=lf_map.operator-(lf_map.get_min(med+(scale*MAD_map))).get_threshold(0,false,true);
	mask+=lf_map.get_max(med-(scale*MAD_map)).operator-(lf_map).get_threshold(0,false,true);
	mask.threshold(0,false,true);
	return mask;
}

CImg<> lifetime_repair(CImg<> & lf_map,float scale = 3, int patch_size=30, CImg<> process_mask=CImg<>::empty()){
	//CImg<unsigned char> mask=get_lifetime_outliers_mask_mad(lf_map,scale,patch_size);
	CImg<unsigned char> mask=get_lifetime_outliers_mask_pseudo_mad(lf_map,scale,patch_size);
	(mask,mask.get_mul(process_mask.get_threshold(0,false,true))).display();
	mask.mul(process_mask.get_threshold(0,false,true));
	mask.save("inpaint_mask.tif");
	return get_inpaint_interp_space(lf_map,mask);
}

int main_lifetime_repair(int argc, char *argv[])
{
	const char * lf_file  = cimg_option("-i",(char *)NULL, "plain lf file for lifetime reconstruction");
	const float scale  = cimg_option("-s",3, "");
	const int patch_size  = cimg_option("-N",30, "");
	CImg<> lf(lf_file);
	(lf,lifetime_repair(lf,scale,patch_size)).save("lf_inpaint.tif");
	return 0;
}



#endif /* _LIFETIME_H_ */
