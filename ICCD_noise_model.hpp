#ifndef _ICCD_NOISE_MODEL_H_
#define _ICCD_NOISE_MODEL_H_

//My include
#include "CImg.h"
#include "sin_fit.hpp"

using namespace cimg_library;
using namespace std;

typedef CImg<float>  cfl;

class quadratic_model : public Sin_fit{
public:
	quadratic_model(CImg<float> vars, CImg<float> means){
		init(vars,means);
	}
	quadratic_model(){};

	void init(CImg<float> & vars, CImg<float> & means){
		CImg<float> input_beta(3); input_beta.fill(0);
		Sin_fit::init(vars,means,input_beta);
	}

protected:
	float reg_funct( float x){
		return beta(0)+x*beta(1)+x*x*beta(2);
	}

	float reg_funct_deriv( int beta_idx,float x){
		float res;
		switch(beta_idx){
		case 0:
			res = 1;
			break;
		case 1:
			res=x ;
			break;
		case 2:
			res=x*x;
				break;
		default:
			cout << "FATAL ERROR in deriv_quadratic_funct" << endl;
		}

		return res;
	}

};

class linear_model : public Sin_fit{
public:
	linear_model(CImg<float> vars, CImg<float> means){
		init(vars,means);
	}

	linear_model(){};

	void init(CImg<float> & vars, CImg<float> & means){
		CImg<float> input_beta(2); input_beta.fill(0);
		Sin_fit::init(vars,means,input_beta);
	}

protected:
	float reg_funct( float x){
		return beta(0)+x*beta(1);
	}

	float reg_funct_deriv( int beta_idx,float x){
		float res=0;
		switch(beta_idx){
		case 0:
			res = 1;
			break;
		case 1:
			res=x ;
			break;
		default:
			cout << "FATAL ERROR in linear_derive_funct" << endl;
		}
		return res;
	}

};

struct CCD_noise_model{
	double g0,edc,gauss_e,gauss_sig;
	CImg<double> cov;
	int mean_method;
	int variance_method;
	bool visu;
	bool robust;
	Sin_fit * estimator;
	int scheme;
	int b_size;
	CCD_noise_model(){
		g0=0;
		b_size=7;
		edc=0;
		gauss_e=0;
		gauss_sig=0;
		mean_method=1;
		variance_method=3;
		visu=false;
		robust=true;
		scheme=2;
		estimator= new linear_model();
	};
	~CCD_noise_model(){delete estimator;
		if (robust) delete estimator->get_cost_function();}
	template<class T>
	void simulate_noise( CImg<T> & input){
		input.noise(0,3);
		input*=g0;
		input.noise(gauss_sig,0);
		input+=gauss_e;
	}
	template<class T>
	void estimate_param(const CImg<T> & input)
	// This method consider the model as a linear one, without
	// considering the variance of the the gain
	{
		CImg<float> vars,means,sizes;
		switch (scheme){
		case 0:  input.gat_compute_blocks_stats(input,means,vars,mean_method,variance_method,visu);
			break;
		case 2:  input.gat_compute_quadtree_stats(input,means,vars,mean_method,variance_method,visu) ;
			break;
		case 3:  compute_block_pixel(input,means,vars) ;
			break;
		// case 4:  input.gat_compute_quadtree_stats(input,means,vars,sizes,mean_method,variance_method,visu) ;
		// 	break;
		default: input.gat_compute_quadtree_stats(input,means,vars,mean_method,variance_method,visu);
			break;
		}
		estimator->init(means,vars);
		// if (scheme==4) estimator->set_residue_weigth(sizes);
		estimator->run();
		if(robust){
			float scale_parameter;
			// if(scheme==4){
			// 	cout << "compute scheme 4 scale" << endl;
			// 	CImg<double> tmp(sizes.sum());
			// 	int i =0;
			// 	const CImg<float> & r =estimator->get_residue();
			// 	cimg_foroff(sizes,s){
			// 		for(int k = 0 ; k < sizes(s);k++){
			// 			tmp(i++)=r(s);
			// 		}
			// 	}
			// 	scale_parameter=sqrt(tmp.variance(2));

			// }
			// else
			scale_parameter=sqrt(estimator->get_residue().variance(2));

			Leclerc * cost_leclerc= new Leclerc(scale_parameter,18);
			estimator->set_cost_function(cost_leclerc);
			estimator->irls_no_init();
		}
		compute_gain_and_bias();
	}


	void compute_gain_and_bias(){
		g0=estimator->get_beta()(1);
		edc = estimator->get_beta()(0);
	}

	float get_variance(float x){
		return edc+g0*x;
	}
	template < class T >
	void compute_block_pixel(const CImg<T>& input,CImg<float>& means,CImg<float> & vars){
		CImg<float> tmp_res;
		float v, m;
		CImg<float> var_array(input.size()/b_size);
		CImg<float> mean_array(input.size()/b_size);
		int count=0;
		cimg_forZ(input,z){
			const CImg<T> & current= input.get_slice(z);
			CImg<float> residuals=current.get_pseudo_residuals();
			char filename[40]; cimg::number_filename("pseudo_residual.tiff",z,2,filename);
			residuals.save(filename);
			for(int y = 0; y < current.height()-b_size; y+=b_size){
				for(int x = 0; x < current.width()-b_size; x+=b_size){
					tmp_res=residuals.get_crop(x,y,x+b_size-1,y+b_size-1);
					v=tmp_res.variance(variance_method);
					m=current.get_crop(x,y,x+b_size-1,y+b_size-1).median();
					if(v>0 && !isnan(v) && !isnan(m)){
						var_array[count]=v;
						mean_array[count]=m;
						count++;
					}
				}}
		}
		means.assign(mean_array.data(),count);
		vars.assign(var_array.data(),count);
	}

	void display(){estimator->display_fit();}


	template<class T>
	cfl apply_gat(CImg<T> & input){
		cfl gat_input(input,"xyzc");
		cimg_foroff(gat_input,i){
			gat_input[i]=gat((float)input[i]);
		}
		return gat_input;
	}

	float gat(float x){
		float res=2*std::sqrt(g0*x+3.f/8.f*g0*g0+edc)/g0;
		return res;
	}


	float inv_gat(float x){
		float res=(cimg::sqr(g0*x/2) -3.f/8.f*g0*g0-edc)/g0;
		return res;
	}
	template<class T>
	cfl apply_inv_gat(CImg<T> & input){
		cfl gat_input(input,"xyzc");
		cimg_foroff(gat_input,i){
			gat_input[i]=inv_gat((float)input[i]);
		}
		return gat_input;
	}


	template<class T>
	cfl apply_gat_local_approx(CImg<T> & input){
		cfl gat_input(input);
		cimg_foroff(gat_input,i){
			gat_input[i]=2*std::sqrt(g0*(input[i])+edc)/g0;
		}
		return gat_input;
	}
	void print(){
		cout << "estimated parameter : " << endl;
		cout << "gain  : " << g0<<endl;
		cout << "DC/RO gaussian noise : N(" << gauss_e<<","<<gauss_sig<<")"<<endl;
		cout << "constant compoment : " << edc << endl;
		// cout << "cov : " << cov(0,0) << " , " << cov(0,1) << " , " << cov(1,1)  << endl;
		cout << endl;
	}
};

struct ICCD_noise_model{
	double g0_e,gccd,g0_sig,gauss_e,gauss_sig;
	float C0,C1,C2; 							// coef ruling the quadratic_model
	CImg<double> cov;
	CCD_noise_model ccd_model;
	ICCD_noise_model(){
		delete ccd_model.estimator;
		ccd_model.estimator= new quadratic_model();
		ccd_model.scheme=3;
		ccd_model.robust=true;
	}

	template <class T>
	void estimate_param(const CImg<T> & input){
		ccd_model.estimate_param(input);
		compute_gain_and_bias();
	}

	void compute_gain_and_bias(){
		C0=ccd_model.estimator->get_beta(0);
		C1=ccd_model.estimator->get_beta(1);
		C2=ccd_model.estimator->get_beta(2);
	}
	template<class T>
	void simulate_noise( CImg<T> & input){
		input.noise(0,3);
		cfl giccd(input,"xyzc",g0_e); giccd.noise(g0_sig,0);
		input.mul(giccd);
		input.noise(gauss_sig,0);
		input+=gauss_e;
	}

	void set_parameters(float ag0_e,float ag0_sig,float gccd, float agauss_e,float agauss_sig){
		g0_e=ag0_e; g0_sig=ag0_sig; gauss_e=agauss_e; gauss_sig=agauss_sig;
		C2 = cimg::sqr(g0_sig/g0_e);
		C1 = gccd*(g0_e+cimg::sqr(g0_sig)/g0_e) - 2*gauss_e*C2;
		C0 = C2*gauss_e*gauss_e - gauss_e*gccd*(g0_e+cimg::sqr(g0_sig)/g0_e) +cimg::sqr(gauss_sig);
	}

	float gat_local_approx(float x){
		float t1 = (C2*(x*(C2*x+C1)+C0));
		float res = std::log(abs(2*std::sqrt(abs(t1)) +2*C2*x+C1))/std::sqrt(abs(C2));
		return res;
	}

	float inv_gat_local_approx(float x){
		float res=		(std::exp(-sqrt(C2)*x)*(-2*C1*std::exp(sqrt(C2)*x)
																					-4*C2*C0+std::exp(2*sqrt(C2)*x)+C1*C1))/(4*C2);
		return res;
	}

	float get_variance(float x){
		return C0+C1*x+C2*x*x;
	}

	template<class T>
	cfl apply_gat_local_approx(CImg<T> & input){
		cfl gat_input(input,"xyzc");
		cimg_foroff(gat_input,i){
			gat_input[i]=gat_local_approx((float)input[i]);
		}
		return gat_input;
	}
	template<class T>
	cfl apply_inv_gat_local_approx(CImg<T> & input){
		cfl gat_input(input,"xyzc");
		cimg_foroff(gat_input,i){
			gat_input[i]=inv_gat_local_approx((float)input[i]);
		}
		return gat_input;
	}

	void print(){
		cout << "estimated parameter : " << endl;
		// cout << "gain noise : N(" << g0_e<<","<<g0_sig<<")"<<endl;
		// cout << "DC/RO gaussian noise : N(" << gauss_e<<","<<gauss_sig<<")"<<endl;
		cout << "C0, C1, C2 " << C0 << " " << C1 << " " << C2 << endl;
	}
};

struct ICCD_abber_noise_model{
	float sigma_x,sigma_y,x_0,y_0,offset; // coef ruling spatial abbeatation
	ICCD_noise_model iccd_model;
	ICCD_abber_noise_model(){
	}

	template <class T>
	void estimate_param(const CImg<T> & input){
		cout << "WARNING :  estimate_param is not implemented for ICCD_abber_noise_model" << endl;
	}

	void compute_gain_and_bias(){
		// C0=ccd_model.estimator->get_beta(0);
		// C1=ccd_model.estimator->get_beta(1);
		// C2=ccd_model.estimator->get_beta(2);
	}

	float kernsig_to_noisescale(float kernsig, int ks){
		static CImg<> kernel; 
		float ka=1;
		int km=floor(ks/2);
		kernel.assign(ks,ks).draw_gaussian(km,km,kernsig,&ka);
		return (kernel/kernel.sum()).sqr().sum();
	}

	CImgList<> kernsigs_to_noisescales(int ks,int sampleNb){
		static CImg<> kernel; 
		int km=floor(ks/2);
		CImg<> noisescales(sampleNb);
		CImg<> kernsigs(sampleNb); kernsigs.sequence(0.0001,km+1);
		kernsigs.print("kernsigs");
		float *scale_ptr=noisescales.data();
		float *sig_ptr=kernsigs.data();
		cimg_foroff(kernsigs,i) *scale_ptr++=kernsig_to_noisescale(*sig_ptr++,ks);
		return (kernsigs,noisescales);
	}


	
	template<class T>
	void simulate_noise( CImg<T> & input){
		iccd_model.simulate_noise(input);
		CImg<T> abber_free(input);
		int ks=5;
		int km=floor(ks/2);
		float ka=1;
		CImgList<> kern_scale=kernsigs_to_noisescales(ks,10000);
		CImg<> scaleMap(input.width(),input.height());
		CImg<> noiseSigMap(scaleMap);
		cimg_for_insideXY(input,x,y,km){
			float noisesig=((exp(-(x-x_0)*(x-x_0)/(2*sigma_x*sigma_x)-(y-y_0)*(y-y_0)/(2*sigma_y*sigma_y))+offset));
			CImg<> minFun((kern_scale[1]-noisesig).abs());
			float &min_ref=minFun.min();
			float bestScaleIdx=(&min_ref-minFun.data());
			float blursig=kern_scale[0](bestScaleIdx);
			scaleMap(x,y)=blursig;
			noiseSigMap(x,y)=noisesig;
			CImg<> kernel; kernel.assign(ks,ks).draw_gaussian(km,km,blursig,&ka); kernel=kernel/kernel.sum();
			cimg_forZ(input,z){
				float filteredVal=0;
				for(int xi=-km;xi<(km+1);xi++)
					for(int yi=-km;yi<(km+1);yi++){
						filteredVal+=abber_free(x+xi,y+yi,z)*kernel(km+xi,km+yi);
					}
				input(x,y,z)=filteredVal;
			}
		}
		input.crop(2*km,2*km,0,input.width()-2*km-1,input.height()-2*km-1,input.depth()-1);

		// CFigure fig; fig.plot(kern_scale[0],kern_scale[1],"r-").set_axis();
		// fig.plot(scaleMap,noiseSigMap,"b+").display();
		// scaleMap.display();
		// (abber_free,input).display();
	}


	void set_parameters(float aC2, float aC1, float aC0  ,
											float ax_0, float ay_0,float asigma_x, float asigma_y, float aoffset){
		iccd_model.C2 = aC2; iccd_model.C1 = aC1; iccd_model.C0 = aC0;
		sigma_x=asigma_x;	sigma_y=asigma_y;
		x_0=ax_0; y_0=ay_0; 
		offset=aoffset;
	}

	void set_parameters(float ag0_e,float ag0_sig,float gccd, float agauss_e,float agauss_sig,
											float ax_0, float ay_0,float asigma_x, float asigma_y, float aoffset){
		iccd_model.set_parameters( ag0_e, ag0_sig, gccd,  agauss_e, agauss_sig);
		sigma_x=asigma_x;	sigma_y=asigma_y;
		x_0=ax_0; y_0=ay_0; 
		offset=aoffset;
	}

	float get_variance(float I,float x,float y){
		return iccd_model.get_variance(I)*
			(exp(-(x-x_0)*(x-x_0)/(2*sigma_x*sigma_x)-(y-y_0)*(y-y_0)/(2*sigma_y*sigma_y))+offset);
	}

	void print(){
		cout << "estimated parameter : " << endl;
		cout << "C0, C1, C2 " << iccd_model.C0 << " " << iccd_model.C1 << " " << iccd_model.C2 << endl;
		cout << "abberation center : " << x_0 << " , " << y_0 << endl;
		cout << "abberation scale : " << sigma_x << " , " << sigma_y << endl;
	}
};


void test_ccd_noise_model_bias()
{
	CCD_noise_model noise_model;
	noise_model.g0=1;
	noise_model.gauss_e=100;
	noise_model.gauss_sig=10;
	noise_model.edc=cimg::sqr(noise_model.gauss_sig)-  noise_model.g0*noise_model.gauss_e;

	float max_ramp = 600;
	int test_value_number = 600;
	int realisation_number=1000;
	cfl ramp(realisation_number,test_value_number);
	cimg_forY(ramp,y) ramp.get_shared_line(y).fill(y*max_ramp/test_value_number);
	cfl true_ramp=cfl::sequence(test_value_number,0,max_ramp);


	noise_model.simulate_noise(ramp);
	ramp.print("noisy ramp");

	cfl stabilized_ramp=noise_model.apply_gat(ramp);
	stabilized_ramp.print("stab");

	cfl mean_stab_ramp(test_value_number);
	cimg_forX(mean_stab_ramp,x)
		mean_stab_ramp(x)=stabilized_ramp.get_shared_line(x).mean();
	cfl inverse_ramp_expectation=noise_model.apply_inv_gat(mean_stab_ramp);

	cfl bias=true_ramp-inverse_ramp_expectation;

	bias.print("ccd bias");
	CFigure fig_bias;
	cfl X=cfl::sequence(test_value_number,0,max_ramp);
	fig_bias.set_axis(X,bias).plot(X,bias,"r-").save("ccd_bias.png");
	X(0)=1;
	cfl bias_coef=bias.get_div(X);
	bias_coef(0)=0;X(0)=0;
	fig_bias.clear().plot(X,bias_coef,"r-").ylabel("bias coeff").set_axis().replot().display();
	fig_bias.clear().plot(X,bias_coef,"r-").erase().ylabel("bias coeff").
		set_axis().replot().save("ccd_bias_coef.png");

};

void test_iccd_noise_model_bias()
{
	ICCD_noise_model noise_model;
	noise_model.set_parameters(10,3,1,100,10);

	float max_ramp = 600;
	int test_value_number = 600;
	int realisation_number=1000;
	cfl ramp(realisation_number,test_value_number);
	cimg_forY(ramp,y) ramp.get_shared_line(y).fill(y*max_ramp/test_value_number);
	cfl true_ramp=cfl::sequence(test_value_number,0,max_ramp);

	noise_model.simulate_noise(ramp);
	ramp.print("noisy ramp");

	cfl stabilized_ramp=noise_model.apply_gat_local_approx(ramp);
	stabilized_ramp.print("stab");

	cfl mean_stab_ramp(test_value_number);
	cimg_forX(mean_stab_ramp,x)
		mean_stab_ramp(x)=stabilized_ramp.get_shared_line(x).mean();
	cfl inverse_ramp_expectation=noise_model.apply_inv_gat_local_approx(mean_stab_ramp);

	cfl bias=true_ramp-inverse_ramp_expectation;

	bias.print("iccd bias");
	CFigure fig_bias;
	cfl X=cfl::sequence(test_value_number,0,max_ramp);
	fig_bias.set_axis(X,bias).plot(X,bias,"r-").save("iccd_bias.png");
	X(0)=1;
	cfl bias_coef=bias.get_div(X);
	bias_coef(0)=0;X(0)=0;
	fig_bias.clear().plot(X,bias_coef,"r-").ylabel("bias coeff").set_axis().replot().display();
	fig_bias.clear().plot(X,bias_coef,"r-").erase().ylabel("bias coeff").
		set_axis().replot().save("iccd_bias_coef.png");

};

#endif /* _ICCD_NOISE_MODEL_H_ */
