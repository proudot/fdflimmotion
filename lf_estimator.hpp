#ifndef _LF_ESTIMATOR_H_
#define _LF_ESTIMATOR_H_

#include "CImg.h"
#include "spot.hpp"
#include "sin_fit.hpp"
#include"sin_fit_bg_exp.hpp"
#include "gauss_kernel.hpp"
#include "sin_fit.hpp"
#include "int_estimator.hpp"
#include "bg_stack_estimator.hpp"
#include "ICCD_noise_model.hpp"

template<class T>
class lf_estimator
{
public:
	CImg<float> weigths;
	inline lf_estimator(){
		noise_model=NULL;
		phase_ref=2.42342;
		tau_ref=	4.1 ;};

	inline lf_estimator(Spot<T> * s, Sin_fit * est){
		noise_model=NULL;
		init(s,est);
	};

	virtual void init(Spot<T> * s, Sin_fit * est){
		noise_model=NULL;
		sp=s;
		estimator=est;
		phase_ref=2.42342;
		tau_ref=	4.1 ;
	};

	inline void record_spot(){
		float Pi=3.141592;	 //  Pi value
		float w=2*Pi*0.04;	 // the Li-Flim software is usually set up to 40
		float phase;

		sp->set_C(estimator->get_C());
		sp->set_amp(estimator->get_amp());
		phase=estimator->get_phase();
		sp->set_phase(phase);
		float phase_shift=phase_ref-atan(w*tau_ref);
		float lifetime=tan(-phase_shift+phase)/w;
		if(lifetime != lifetime) lifetime=0; 		// Nan detection trick
		sp->set_lifetime(fabs(lifetime));
	};
	Sin_fit * get_estimator(){return estimator;};
	Spot<T> * sp;
	Sin_fit * estimator;
	ICCD_noise_model * noise_model;
	float phase_ref;
	float tau_ref;

};

template<class T>
class lf_int : public lf_estimator<T>
{
public:
	inline lf_int(Spot<T> * s, Sin_fit * est){
		init(s,est);
	};

	inline void init(Spot<T> * s, Sin_fit * est){
		lf_estimator<T>::init(s,est);
		Sin_fit & fitter = *(this->estimator);
		CImg<float> Y;
		CImg<float> X;

		Y.assign(this->sp->get_intensities());
		X.assign(this->sp->get_intensities());
		cimg_forX(X,i){X(i)=i;};

		fitter.init(X,Y);
	}
};
template<class T>
class lf_max : public lf_estimator<T>
{
public:
	inline lf_max(Spot<T> * s, Sin_fit * est,CImg<T> & stack){
		init(s,est,stack);
	};
	lf_max(){};

	inline void init(Spot<T> * s, Sin_fit * est,CImg<T> &stack){
		lf_estimator<T>::init(s,est);
		Sin_fit & fitter = *(this->estimator);
		CImg<float> Y;
		CImg<float> X;

		Y.assign(s->get_move_list().size());
		X.assign(s->get_move_list().size());
		cimg_for1(s->get_move_list().size(),t){
			s->set_center(s->get_move(t));
			cfl feature_value=s->get_values(stack.get_slice(t));
			Y(t)=feature_value.max();
			X(t)=t;
		};

		fitter.init(X,Y);
	}
};
CImgList<float> g_display;
template<class T>
struct sig_estimator{
	CImgList<float>  gpatch;
	float sig_x,sig_y;
	inline float get_sig_x(){return sig_x;};
	sig_estimator(Spot<T> * f,CImg<T> & img){
		gpatch=CImgList<float>();
		f->gradient(img,gpatch,Gauss_kernel::k,Gauss_kernel::k_deriv);
	}
	inline CImg<float> run(){
		g_display.push_back(gpatch);
		float gmax, gmin;
		int xmax,ymax,xmin,ymin;

		max_idx(gpatch[0],gmax,xmax,ymax);
		CImg<float> minus_g=- gpatch[0];
		max_idx(minus_g,gmin,xmin,ymin);
		sig_x=fabs(xmax-xmin)/2;

		max_idx(gpatch[1],gmax,xmax,ymax);
		minus_g=-gpatch[1];
		max_idx(minus_g,gmin,xmin,ymin);
		sig_y=fabs(ymax-ymin)/2;

		CImg<float> res(2);
		res(0)=sig_x;
		res(1)=sig_y;
		return res;
	};
};

// float  test_sig_est[30*12];
// float  test_sig_fit[30*12];
// int test_sig_nb=0;

			// sig_estimator<T> se(s,c);
			// se.run();
			// test_sig_est[test_sig_nb]= se.get_sig_x();
			// test_sig_fit[test_sig_nb]= s->get_sig_x();
			// test_sig_nb++;

template<class T>
class lf_multiple : public lf_estimator<T>
{
public:
	inline lf_multiple(Spot<T> * s, Sin_fit * est,CImg<T> & stack){
		init(s,est,stack);
	};
	lf_multiple(){};

	inline void init(Spot<T> * aSpot, Sin_fit * est,CImg<T> & stack){
		lf_estimator<T>::init(aSpot,est);
		Sin_fit & fitter = *(this->estimator);
		Spot<T> * s = (this->sp);
		CImg<float> feature_value;
		CImg<float> Y;
		CImg<float> X;
		CImg<float> patch(s->get_width(),s->get_height());
		CImg<T> c;
		int idx;

		Y.assign(s->get_height()*s->get_width()*s->get_move_list().size());
		X.assign(s->get_height()*s->get_width()*s->get_move_list().size());
		idx=0;

		// float Pi=3.141592;
		// float wf=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
		// float w=2*Pi/12;
		// float threshold = 4000;
		// float C_bg=1;
		// float amp_bg=0.5;
		// float phase_bg=(atan(2.5*wf)+2.42342-atan(4.1*wf));
		// CImg<float> bg(Y);
		// cimg_forX(bg,t){
		// 	bg(t) = threshold*(C_bg+amp_bg*sin(w*t+phase_bg));
		// }

		for(unsigned int t = 0; t < s->get_move_list().size();t++){
			s->set_center(s->get_move(t));
			c=stack.get_slice(t);
			feature_value=get_f_values(c);
			T white =1;
			CImg<float> rounded_center=s->get_center().get_round();
			float dx = s->get_center()(0)-rounded_center(0);
			float dy = s->get_center()(1)-rounded_center(1);
			patch.draw_gaussian(patch.width()/2+dx,patch.height()/2+dy,s->get_sig_x(),&white);
			cimg_forXY(feature_value,x,y){
				Y(idx)=feature_value(x,y);
				// Y(idx)=(feature_value(x,y) - bg(t))/patch(x,y);
				X(idx)=t;
				idx++;
			}
		}

		fitter.init(X,Y);
	}
protected:

	inline virtual CImg<float> get_f_values(CImg<T> & img){
		return this->sp->get_values(img);
	}
};
template<class T>
class lf_multiple_no_interp_renorm : public lf_estimator<T>
{
public:
	inline lf_multiple_no_interp_renorm(Spot<T> * s, Sin_fit * est,CImg<T> & stack){
		init(s,est,stack);
		this->noise_model=NULL;
	};

	inline void init(Spot<T> * aSpot, Sin_fit * est,CImg<T> & stack){
		lf_estimator<T>::init(aSpot,est);
		Sin_fit & fitter = *(this->estimator);
		Spot<T> * s = (this->sp);
		CImg<float> feature_value;
		CImg<float> Y;
		CImg<float> X;
		CImg<float> patch(s->get_width(),s->get_height());
		CImg<T> c;
		int idx;

		Y.assign(s->get_height()*s->get_width()*s->get_move_list().size());
		X.assign(s->get_height()*s->get_width()*s->get_move_list().size());
		// this->weigths.assi gn(s->get_height()*s->get_width()*s->get_move_list().size());
		idx=0;

		for(unsigned int t = 0; t < s->get_move_list().size();t++){
			s->set_center(s->get_move(t));
			c=stack.get_slice(t);
			feature_value=get_f_values(c);
			T white =1;
			CImg<float> rounded_center=s->get_center().get_round();
			float dx = s->get_center()(0)-rounded_center(0);
			float dy = s->get_center()(1)-rounded_center(1);
			patch.draw_gaussian(patch.width()/2+dx,patch.height()/2+dy,s->get_sig_x(),&white);
			cimg_forXY(feature_value,x,y){
				Y(idx)=feature_value(x,y)/patch(x,y);
				// this->weigths(idx)=sqrt(this->noise_model->get_variance(feature_value(x,y))) /patch(x,y);
				X(idx)=t;
				idx++;
			}
		}


		fitter.init(X,Y);
	}
protected:

	inline virtual CImg<float> get_f_values(CImg<T> & img){
		return this->sp->get_values_no_interp(img);
	}
};
template<class T>
class lf_multiple_no_interp_renorm_bg : public lf_estimator<T>
{
public:
	inline lf_multiple_no_interp_renorm_bg(Spot<T> * s, Sin_fit * est,CImg<T> & stack,float asig){
		sig=asig;
		init(s,est,stack);
	};
	float sig;

	inline CImg<T> get_current_bg(CImg<T> & stack){
		CImg<T> region;
		CImg<T> rebuild;
		Bg_estim bge(sig);


		// int x_min,x_max,y_min,y_max;
		// x_min=s->get_move(0)(0);
		// x_max=s->get_move(0)(1);
		// y_min=s->get_move(0)(0);
		// y_max=s->get_move(0)(1);
		// for(unsigned int t = 0; t < s->get_move_list().size();t++){
		// 	if (s->get_move(t)(0) < x_min) x_min=s->get_move(t)(0);
		// 	if (s->get_move(t)(1) < y_min) y_min=s->get_move(t)(1);
		// 	if (s->get_move(t)(0) > x_max) x_max=s->get_move(t)(0);
		// 	if (s->get_move(t)(1) > y_max) y_max=s->get_move(t)(1);
		// }

		// region=stack.get_crop(x_min-region_size,);
		cimg_forZ(stack,z){
			region.append(this->sp->get_values_no_interp(stack.get_shared_plane(z)),'z');
		}
		rebuild=bge.get_rebuild_bg(region);

		return rebuild;
	};

	inline void init(Spot<T> * aSpot, Sin_fit * est,CImg<T> & stack){
		lf_estimator<T>::init(aSpot,est);
		Sin_fit & fitter = *(this->estimator);
		Spot<T> * s = (this->sp);
		CImg<float> feature_value;
		CImg<float> Y;
		CImg<float> X;
		CImg<float> patch(s->get_width(),s->get_height());
		CImg<T> c;
		int idx;
		CImg<T> bg;

		Y.assign(s->get_height()*s->get_width()*s->get_move_list().size());
		X.assign(s->get_height()*s->get_width()*s->get_move_list().size());
		idx=0;

		for(unsigned int t = 0; t < s->get_move_list().size();t++){
			s->set_center(s->get_move(t));
			c=stack.get_slice(t);
			feature_value=get_f_values(c);
			T white =1;
			CImg<float> rounded_center=s->get_center().get_round();
			float dx = s->get_center()(0)-rounded_center(0);
			float dy = s->get_center()(1)-rounded_center(1);
			float sigma = 1.5;
			patch.draw_gaussian(patch.width()/2+dx,patch.height()/2+dy,sigma,&white);
			bg=get_current_bg(stack).get_slice(t);
			cimg_forXY(feature_value,x,y){
				Y(idx)=(feature_value(x,y)-bg(x,y))/patch(x,y);
				X(idx)=t;
				idx++;
			}
		}

		fitter.init(X,Y);
	}
protected:

	inline virtual CImg<float> get_f_values(CImg<T> & img){
		return this->sp->get_values_no_interp(img);
	}
};

template<class T>
class lf_multiple_no_interp : public lf_multiple<T>{
public:
	inline lf_multiple_no_interp(Spot<T> * s, Sin_fit * est,CImg<T> & stack){
		this->init(s,est,stack);
	};

protected:
	inline CImg<float> get_f_values(CImg<T> & img){
		return this->sp->get_values_no_interp(img);
	}
};


template<class T>
class lf_multiple_deriv : public lf_estimator<T>
{
public:
	inline lf_multiple_deriv(Spot<T> * s, Sin_fit * est,CImg<T> & stack){
		init(s,est,stack);
	};
	lf_multiple_deriv(){};

	inline void init(Spot<T> * aSpot, Sin_fit * est,CImg<T> & stack){
		lf_estimator<T>::init(aSpot,est);
		Sin_fit & fitter = *(this->estimator);
		Spot<T> * s = (this->sp);
		CImg<float> feature_value;
		CImgList<float> g;
		CImg<float> Y;
		CImg<float> X;
		int idx;

		Y.assign(s->get_height()*s->get_width()*s->get_move_list().size()*2);
		X.assign(s->get_height()*s->get_width()*s->get_move_list().size()*2);
		idx=0;
		for(unsigned int t = 0; t < s->get_move_list().size();t++){
			s->set_center(s->get_move(t));
			CImg<T> c = stack.get_slice(t);
			g=get_gradient(c);
			cimg_forXY(g[0],x,y){
				Y(idx)=fabs(g[0](x,y));
				Y(idx+1)=fabs(g[1](x,y));
				X(idx)=t;
				X(idx+1)=t;
				idx+=2;
			}
		}

		fitter.init(X,Y);
	}
	inline virtual CImgList<float> get_gradient(CImg<T> & img){
		CImgList<float> g;
		this->sp->gradient(img,g,Gauss_kernel::k,Gauss_kernel::k_deriv);
		return g;
	}

};


template<class T>
class lf_multiple_deriv_no_interp : public lf_multiple_deriv<T>
{
public:
	inline lf_multiple_deriv_no_interp(Spot<T> * s, Sin_fit * est,CImg<T> & stack){
		this->init(s,est,stack);
	};
	inline CImgList<float> get_gradient(CImg<T> & img){
		CImgList<float> g;
		this->sp->get_gradient_no_interp(img,g,Gauss_kernel::k,Gauss_kernel::k_deriv);
		return g;
	}
};

// template<class T>
// struct max_idx {
// 	CImg<T> & img;
// 	T max;
// 	int x;
// 	int y;
// 	max_idx(CImg<T> & i):img(i){};
// 	inline float run(){
// 		T prov_max=img(0,0);
// 		cimg_forXY(img,i,j){
// 			if (img(i,j) > prov_max){
// 				prov_max=img(x,y);
// 				x=i;
// 				y=j;
// 			}
// 			max=prov_max;
// 			return max;
// 		}
// 	}
// }
template<class T>
class lf_multiple_max_deriv : public lf_multiple_deriv<T>
{
public:
	inline lf_multiple_max_deriv(Spot<T> * s, Sin_fit * est,CImg<T> & stack){
		init(s,est,stack);
	};
	lf_multiple_max_deriv(){};

	inline void init(Spot<T> * aSpot, Sin_fit * est,CImg<T> & stack){
		lf_estimator<T>::init(aSpot,est);
		Sin_fit & fitter = *(this->estimator);
		Spot<T> * s = (this->sp);
		CImg<float> feature_value;
		CImgList<float> g;
		CImg<float> Y;
		CImg<float> X;
		int idx;

		Y.assign(s->get_move_list().size()*4);
		X.assign(s->get_move_list().size()*4);
		idx=0;
		float sig = 1.5;
		float coef = sig/exp(-1/2);
		int x_m,y_m;
		float local_max;
		// guess sigma using max locus


		for(unsigned int t = 0; t < s->get_move_list().size();t++){
			s->set_center(s->get_move(t));
			CImg<T> c = stack.get_slice(t);
			g=this->get_gradient(c);
			max_idx(g[0],local_max,x_m,y_m);
			Y(idx)=g[0].max()*coef;
			Y(idx+1)=fabs(g[0].min()*coef);
			Y(idx+2)=g[1].max()*coef;
			Y(idx+3)=fabs(g[1].min()*coef);
			X(idx)=t;
			X(idx+1)=t;
			X(idx+2)=t;
			X(idx+3)=t;
			idx+=4;
		}

		fitter.init(X,Y);
	}
};

template<class T>
class lf_multiple_max_deriv_no_interp : public lf_multiple_max_deriv<T>
{
public:
	inline lf_multiple_max_deriv_no_interp(Spot<T> * s, Sin_fit * est,CImg<T> & stack){
		this->init(s,est,stack);
	};
	inline CImgList<float> get_gradient(CImg<float> img){
		CImgList<float> g;
		this->sp->get_gradient_no_interp(img,g,Gauss_kernel::k,Gauss_kernel::k_deriv);
		return g;
	}

};

// struct ICCD_noise_variance: public INoise_variance, public ICCD_noise_model{
// 	float get_variance(float x ){return ICCD_noise_model::get_variance(x);}};

template<class T>
struct Lf_switch{
	lf_estimator<T> * lf_est;
	ICCD_noise_model noise_model;
	CFigure fig;
	void run(Spot<T> * s, CImg<T> & stack, int type ,const char * name_dat=NULL){

	if(type !=1){
		// Sin_fit_C_amp estimator;

		int scale_parameter_type_lf=0;
		switch(type){
		case 0 : {
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_int<T>(&(*s),estimator);
			break;
		}
		case 2:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_int<T>(&(*s),estimator);
			scale_parameter_type_lf=1;
			break;
		}
		case 3:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple<T>(&(*s),estimator,stack);
			break;
		}
		case 4:{
			Sin_fit * estimator = new Sin_fit();
			lf_est = new lf_multiple<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=1;
			break;
		}
		case 5:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=2;
			break;
		}
		case 6:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_deriv<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=0;
			break;
		}
		case 7:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_max_deriv<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=0;
			break;
		}
		case 8:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_no_interp<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=0;
			break;
		}
		case 9:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_max_deriv_no_interp<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=0;
			break;
		}
		case 10:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_deriv_no_interp<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=0;
			break;
		}
		case 11:{
			Sin_fit_bg_exp<T> * estim_bg_exp = new Sin_fit_bg_exp<T>(s);
			lf_est = new lf_multiple<T>(&(*s),estim_bg_exp,stack);
			estim_bg_exp->set_starting_value(stack);
			break;
		}
		case 12:{
			Sin_fit_bg_exp<T> * estim_bg_exp = new Sin_fit_bg_exp<T>(s);
			lf_est = new lf_multiple_no_interp<T>(&(*s),estim_bg_exp,stack);
			estim_bg_exp->set_starting_value(stack);
			break;
		}
		case 13:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_no_interp_renorm<T>(&(*s),estimator,stack);
			// lf_est->noise_model = new ICCD_noise_model();
			// lf_est->noise_model->estimate_param(stack);
			break;
		}
		case 14:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_no_interp_renorm<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=2;
			break;
		}
		// case 15:{
		// 	Sin_fit * estimator = new Sin_fit();
		// 	varEstimator<T> ve;
		// 	float sig=ve.estimate(stack.get_crop(s->get_move(0)(0)-10,s->get_move(0)(1)-10,
		// 																			 s->get_move(0)(0)+10,s->get_move(0)(1)+10));
		// 	lf_est = new lf_multiple_no_interp_renorm_bg<T>(&(*s),estimator,stack,sig);
		// 	scale_parameter_type_lf=1;
		// 	break;
		// }
		// case 16:{
		// 	static varEstimator<T> ve;
		// 	static float sig = ve.estimate(stack);
		// 	static Bg_estim bge(sig);
		// 	static CImg<T> bg_stack = bge.get_rebuild_bg(stack);
		// 	static bool done = false;
		// 	if(!done ){ bg_stack.save("bg_simu_plain.tiff"); done = true; cout << "save_bg" << endl;}
		// 	CImg<T> stack_no_bg=stack-bg_stack;
		// 	Sin_fit * estimator = new Sin_fit();;
		// 	lf_est = new lf_multiple_no_interp_renorm<T>(&(*s),estimator,stack_no_bg);
		// 	break;
		// }
		case 17:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_no_interp_renorm<T>(&(*s),estimator,stack);
			lf_est->noise_model=&noise_model;
			estimator->set_residue_weigth(lf_est->weigths);
			break;
		}
		case 18:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_multiple_no_interp_renorm<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=2;
			break;
		}
		case 19:{
			Sin_fit * estimator = new Sin_fit();;
			lf_est = new lf_max<T>(&(*s),estimator,stack);
			scale_parameter_type_lf=0;
			break;
		}
		default:
			return;
		}
		Sin_fit * est = lf_est->get_estimator();
		est->set_max_iter(100);
		est->set_threshold(0.01);
		est->run();
		fig.assign();
		fig= est->build_1d_graph();
		if (name_dat!=NULL) est->save_X_Y_gnuplot(name_dat);

		// CImg<> phase_x;
		// CImg<> phase_y;
		// phase_x.assign(2,1,1,1,est->get_beta(2),est->get_beta(2));
		// phase_y.assign(2,1,1,1,est->get_Y().min(),est->get_Y().max());
		// fig.plot(phase_x,phase_y,"r-");
		// 	cout << " C ls : " << est->get_beta(0);
		// 	cout << " A ls : " << est->get_beta(1);
		// 	cout << " phase ls : " << est->get_beta(2);


		if(scale_parameter_type_lf){
			float err_sig=0;
			switch(scale_parameter_type_lf){
			case 1 : {
				// varEstimator<T> ve;
				// err_sig = ve.estimate(stack);
				err_sig = sqrt(stack.get_slice(0).get_pseudo_residuals().variance(3));
				break;
			}
			case 2: {
				CImg<float> r = est->get_residue();
				err_sig = sqrt(r.variance(2));
				break;
			}
			}
			Leclerc cf_leclerc(err_sig,18);
			est->set_cost_function(&cf_leclerc);
			est->irls_no_init();
			est->plot_fit(fig,"g-");
			// phase_x.assign(2,1,1,1,est->get_beta(2),est->get_beta(2));
			// phase_y.assign(2,1,1,1,est->get_Y().min(),est->get_Y().max());
			//fig.plot(phase_x,phase_y,"g-");
			//cout << " C irls : " << est->get_beta(0);
			// cout << " A irls : " << est->get_beta(1);
			// cout << " phase irls : " << est->get_beta(2);
		}
		// if (est->has_diverged()){
		// 	lf_est->sp->set_lost(true);
		// }else
		{ lf_est->record_spot(); }
		delete lf_est->get_estimator();
		// if(type==13) delete lf_est->noise_model;
		delete lf_est;
	}
};
};
#endif /* _LF_ESTIMATOR_H_ */
