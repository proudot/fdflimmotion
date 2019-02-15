#ifndef _MATCH_MODEL_H
#define _MATCH_MODEL_H


#include <iostream>
#include "math.h"
#include "Matcher.hpp"


using namespace cimg_library;
using namespace std;

template<class T>
class Match_smodel: public Matcher<T> {
public :
	// reuse the feature !
	Match_smodel(){};
	inline Match_smodel( Feature<T> * f,  CImg<T> *aCurrent, int aFrame_idx):
		Matcher<T>::Matcher(f,aCurrent){this->frame_idx=aFrame_idx;};
	inline void assign( Feature<T> * f,  CImg<T> *aCurrent, int aFrame_idx){
		this->feature=f; this->current=aCurrent; this->frame_idx=aFrame_idx;};
	virtual void fill_residue() ;
	virtual void fill_J();
	virtual void init();
	inline CImg<float> comp_patch(){ patch_build(); return this->patch;}
	virtual void set_starting_value();
	inline CImg<float> get_patch(){return patch;};
	inline CImg<float> get_res_patch(){return res_patch;};
	virtual void update_feature_state();
	// virtual int irls();						//
	inline void set_frame_idx(float val ){frame_idx=val;};

	
protected :
	CImg<float> patch;
	int frame_idx;
	CImg<float> res_patch;
	virtual void init_estimator();
	virtual void patch_build();
	virtual float modele(float x, float y);
};

// template <class T>
// int Match_smodel<T>::irls()
// {
// 	Feature<T> &feature = *(this->feature);
// 	static ofstream log("log_irls.log");
// 	CImg<float> b;
// 	feature.set_starting_center(feature.get_center());
// 	feature.set_sig_x(2);
// 	init_estimator();
// 	Sin_fit::irls(b, 0, sqrt(2),	&dcost_leclerc);
// 	this->beta=b;
// 	update_feature_state();
// 	feature.set_residue(this->residue.magnitude());
// 	feature.print_all(log);
// 	return this->get_nb_iter();
// }

template<class T>
void Match_smodel<T>::init(){
	// Feature<T> &feature = *(this->feature);
	// feature.set_sig_x(2);				// WARNING
	patch_build();
}

template<class T>
void Match_smodel<T>::init_estimator(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	residue.assign(feature.get_width()*feature.get_height());
	beta.assign(2).fill(0);
	this->init_residue(residue,beta);
}

// TODO optimize(to much recomp)
template<class T>
void Match_smodel<T>::patch_build(){
	Feature<T> &feature = *(this->feature);
	patch.assign(feature.get_width(),
									feature.get_height());
	res_patch.assign(patch);
	float  C= feature.get_C();
	float  amplitude= feature.get_amp();
	int  frame_idx= this->frame_idx;
	float Pi=3.141592;
	float w=2.*Pi/12.;                                   /* FIXME */
	float  phase = feature.get_phase();
	float color_current=C+amplitude*sin(frame_idx*w+phase);
	cimg_forXY(patch,x,y){
		patch(x,y) =	modele(x,y);
	}
	// patch/=patch.mean();
	// mask.assign(patch/patch.mean());
	// mask.fill(1);
	patch*=color_current;
}

template<class T>
float Match_smodel<T>::modele(float x, float y){
	Feature<T> &feature = *(this->feature);
	float x_c= feature.get_width()/2;
	float y_c= feature.get_height()/2;
	float sig = feature.get_sig_x()  ;
	return exp(-((x-x_c)*(x-x_c) + (y-y_c)*(y-y_c))/(2*sig*sig));
}

template<class T>
void Match_smodel<T>::set_starting_value ()
{
	// float Pi=3.141592;
	// float wt=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
	Feature<T> &feature = *(this->feature);
	feature.set_C(17999.5);								 // exact amplitude of simulated spot
	feature.set_amp(11999.7);								 // exact amplitude of simulated spot
	feature.set_phase(2.184);	/* FIXME */
}		/* -----  end of method Sin_fit::starting_value  ----- */



template<class T>
void Match_smodel<T>::update_feature_state(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & beta = this->beta;
	CImg<float> c = feature.get_starting_center()+beta;
	feature.set_center(c);
}

template<class T>
void Match_smodel<T>::fill_residue() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & residue = this->residue;
	update_feature_state();
	CImg<float> f_values_current;
	feature.interpolate(current,f_values_current);

	float i = 0.;

	cimg_forXY(f_values_current,x,y){
		res_patch(x,y)= f_values_current(x,y) - patch(x,y);
		residue(i++) = res_patch(x,y);
	}

}

template<class T>
void Match_smodel<T>::fill_J() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & k = this->k ;
	CImg<float> & J = this->J ;
	CImg<float> & k_deriv = this->k_deriv;
	CImgList<float> G;
	feature.gradient(current,G,k,k_deriv);
	int x= 0;
	int y=0;
	float sig = feature.get_sig_x();

	cimg_forY(J,i){
		x = i % feature.get_width();
		y = i / feature.get_height();
		J(0,i)=-G[0](x,y)/2+patch(x,y)*(x-feature.get_width()/2)/(2*sig*sig);
		J(1,i)=-G[1](x,y)/2+patch(x,y)*(y-feature.get_height()/2)/(2*sig*sig);
	}

};


template<class T>
class Match_model: public Match_smodel<T> {
public :
	Match_model(){};
	inline Match_model( Feature<T> * f,  CImg<T> *aCurrent, int aFrame_idx):
		Match_smodel<T>::Match_smodel(f,aCurrent,aFrame_idx){};
	virtual void update_feature_state();
	virtual void fill_residue() ;
	virtual void fill_J();

protected :
	virtual void init_estimator();
	void fill_spring();
};

template<class T>
void Match_model<T>::fill_spring ()
{
	// float sigma0=2;
	// float delta=0.01;
	// this->spring(2)=delta*(this->beta(2)-sigma0 );
}		/* -----  end of method Sin_fit::build_delta_beta  ----- */

template<class T>
void Match_model<T>::init_estimator(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	residue.assign(feature.get_width()*feature.get_height());
	beta.assign(3).fill(0);
	beta(2)=feature.get_sig_x();

	// Allocate memorie for fitting
	this->init_residue(residue,beta);
}

template<class T>
void Match_model<T>::update_feature_state(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & beta = this->beta;
	CImg<float> d(1,2);
	d(0)=beta(0);
	d(1)=beta(1);
	CImg<float> c = feature.get_starting_center()+d;
	feature.set_center(c);
	feature.set_sig_x(beta(2));
}

template<class T>
void Match_model<T>::fill_residue() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & residue = this->residue;
	update_feature_state();
	this->patch_build();

	CImg<float> f_values_current;
	feature.interpolate(current,f_values_current);

	int i = 0;
	cimg_forXY(f_values_current,x,y){
		residue(i++) =
			f_values_current(x,y) - this->patch(x,y);
	}
}

template<class T>
void Match_model<T>::fill_J() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & patch = this->patch ;
	CImg<float> & k = this->k ;
	float sig_x = feature.get_sig_x();
	CImg<float> & J = this->J ;
	CImg<float> & k_deriv = this->k_deriv;
	CImgList<float> G;
	feature.gradient(current,G,k,k_deriv);
	int x= 0;
	int y=0;

	float x_c= feature.get_width()/2;
	float y_c= feature.get_height()/2;
	cimg_forY(J,i){
		x = i % feature.get_width();
		y = i / feature.get_height();
		J(0,i)=-G[0](x,y)/2+(patch(x,y))*(x-x_c)/(2*sig_x*sig_x);
		J(1,i)=-G[1](x,y)/2+(patch(x,y))*(y-y_c)/(2*sig_x*sig_x);
		J(2,i)= patch(x,y)*((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c)) /(sig_x*sig_x*sig_x);
		// J(3,i)= modele_exp(x,y)*sin(w*frame_idx+beta(4));
		// J(4,i)= modele_exp(x,y)*beta(3)*cos(w*frame_idx+beta(4));
		// J(5,i)= patch(x,y)*2*(x-x_c)*(x-x_c)/(sig_x*sig_x*sig_x);
		// J(6,i)= patch(x,y)*2*(y-y_c)*(y-y_c)/(sig_y*sig_y*sig_y);
	}
};

template<class T>
class Match_model_no_interp: public Match_model<T> {
public :
	virtual void fill_residue() ;
	virtual void fill_J();
	virtual inline float get_I_max(){return I_max;};
	virtual void init();
	float I_max;
protected :
	virtual void patch_build();
};

template<class T>
void Match_model_no_interp<T>::init(){
	Feature<T> &feature = *(this->feature);
	float  C= feature.get_C();
	float  amplitude= feature.get_amp();
	int  frame_idx= this->frame_idx;
	float Pi=3.141592;
	float w=2.*Pi/12.;                                   /* FIXME */
	float  phase = feature.get_phase();

	I_max=C+amplitude*sin(frame_idx*w+phase);
}


template<class T>
void Match_model_no_interp<T>::patch_build(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & patch = this->patch;
	patch.assign(feature.get_width(),
									feature.get_height());
	float color_current=get_I_max();
	CImg<float> rounded_center=feature.get_center().get_round();
	float dx = feature.get_center()(0)-rounded_center(0);
	float dy = feature.get_center()(1)-rounded_center(1);
	// float one =1 ;
	// patch.draw_gaussian(patch.width()/2+dx,patch.height()/2+dy,feature.get_sig_x(),&one);
	// color_current/=patch.mean();
	patch.draw_gaussian(patch.width()/2+dx,patch.height()/2+dy,feature.get_sig_x(),&color_current);
	if(this->flag_gat){
		this->patch.gat(this->g,this->e);
	}
}


template<class T>
void Match_model_no_interp<T>::fill_residue() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & residue = this->residue;
	this->update_feature_state();
	this->patch_build();

	// Here is the main difference, we do not interpolate the values of f
	CImg<float> f_values_current;
	CImg<float> rounded_center=feature.get_center().get_round();
	float & fcx= (rounded_center(0));
	float & fcy= (rounded_center(1));
	f_values_current = current.get_crop(fcx-feature.get_width()/2,
																			fcy-feature.get_height()/2,
																			fcx+feature.get_width()/2,
																			fcy+feature.get_height()/2);
	int i = 0;
	cimg_forXY(f_values_current,x,y){
		residue(i++) = f_values_current(x,y) - this->patch(x,y);
	}
}

template<class T>
void Match_model_no_interp<T>::fill_J() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & patch = this->patch ;
	CImg<float> & k = this->k ;
	CImg<float> & J = this->J ;
	CImg<float> & k_deriv = this->k_deriv;
	CImgList<float> G;
	float sig_x = feature.get_sig_x();
	CImg<float> rounded_center=feature.get_center().get_round();
	CImg<float> center_bk=feature.get_center();
	feature.set_center(rounded_center); // "Don't try this at home" 3 lines
	feature.gradient(current,G,k,k_deriv);
	feature.set_center(center_bk);

	int x= 0;
	int y=0;

	float x_c= feature.get_width()/2;
	float y_c= feature.get_height()/2;
	float dx = feature.get_center()(0)-rounded_center(0);
	float dy = feature.get_center()(1)-rounded_center(1);
	x_c+=dx;
	y_c+=dy;

	cimg_forY(J,i){
		x = i % feature.get_width();
		y = i / feature.get_height();
		J(0,i)=-G[0](x,y)/2+(patch(x,y))*(x-x_c)/(2*sig_x*sig_x);
		J(1,i)=-G[1](x,y)/2+(patch(x,y))*(y-y_c)/(2*sig_x*sig_x);
		J(2,i)= patch(x,y)*((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c))/(sig_x*sig_x*sig_x);
	}
};

template<class T>
class Match_model_ni_I: public Match_model_no_interp<T> {
public :
	// Number 17
	void fill_J();
	virtual inline float get_I_max(){return this->beta(3);};

protected:
	virtual inline void fill_J_2(int i,int x, int y){
		Feature<T> &feature = *(this->feature);
		CImg<float> & patch = this->patch ;
		float sig_x = feature.get_sig_x();
		float x_c= feature.get_width()/2;
		float y_c= feature.get_height()/2;
		CImg<float> rounded_center=feature.get_center().get_round();
		float dx = feature.get_center()(0)-rounded_center(0);
		float dy = feature.get_center()(1)-rounded_center(1);
		x_c+=dx;
		y_c+=dy;
		this->J(2,i)= patch(x,y)*((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c)) /(sig_x*sig_x*sig_x);
	}
	virtual inline void fill_J_3(int i,int x, int y){this->J(3,i)=this->patch(x,y)/get_I_max();};
	virtual void init_estimator();
};

template<class T>
void Match_model_ni_I<T>::init_estimator(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	residue.assign(feature.get_width()*feature.get_height());
	beta.assign(4).fill(0);
	beta(2)=feature.get_sig_x();
	beta(3)=this->I_max;

	// Allocate memorie for fitting
	this->init_residue(residue,beta);
}


template<class T>
void Match_model_ni_I<T>::fill_J() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & patch = this->patch ;
	CImg<float> & k = this->k ;
	float sig_x = feature.get_sig_x();
	CImg<float> & J = this->J ;
	CImg<float> & k_deriv = this->k_deriv;
	CImgList<float> G;
	CImg<float> rounded_center=feature.get_center().get_round();
	CImg<float> center_bk=feature.get_center();
	feature.set_center(rounded_center); // "Don't try this at home" 3 lines
	feature.gradient(current,G,k,k_deriv);
	feature.set_center(center_bk);

	int x= 0;
	int y=0;

	float x_c= feature.get_width()/2;
	float y_c= feature.get_height()/2;
	float dx = feature.get_center()(0)-rounded_center(0);
	float dy = feature.get_center()(1)-rounded_center(1);
	x_c+=dx;
	y_c+=dy;
	cimg_forY(J,i){
		x = i % feature.get_width();
		y = i / feature.get_height();
		J(0,i)=-G[0](x,y)/2+(patch(x,y))*(x-x_c)/(2*sig_x*sig_x);
		J(1,i)=-G[1](x,y)/2+(patch(x,y))*(y-y_c)/(2*sig_x*sig_x);
		fill_J_2(i,x,y);
		fill_J_3(i,x,y);
	}
};


template<class T>
class Match_model_ni_I_alt_0: public Match_model_ni_I<T> {
	// Number 19
	virtual inline void fill_J_3(int i, int x, int y){this->J(3,i)=this->patch(x,y)/this->patch.max();};
};

template<class T>
class Match_model_ni_I_sigma_fixe: public Match_model_ni_I<T> {
	// Number 20
public:
	void init_estimator();
	virtual inline float get_I_max(){return this->beta(2);};
	protected:

	inline virtual void update_feature_state(){Match_smodel<T>::update_feature_state();};
	virtual inline void fill_J_3(int i, int x, int y){};
	virtual inline void fill_J_2(int i, int x, int y){this->J(2,i)=this->patch(x,y)/this->patch.max();};
};
template<class T>
void Match_model_ni_I_sigma_fixe<T>::init_estimator(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	residue.assign(feature.get_width()*feature.get_height());
	beta.assign(3).fill(0);
	beta(2)=this->I_max;
	// Allocate memorie for fitting
	this->init_residue(residue,beta);
}

template<class T>
class Match_model_deriv: public Match_model_no_interp<T> {
public :
	virtual void fill_residue() ;
	virtual void fill_J();
	virtual void init_estimator();
	CImgList<float> build_patch_x_y();
};

template<class T>
void Match_model_deriv<T>::init_estimator(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	residue.assign(feature.get_width()*feature.get_height()*2);
	beta.assign(3).fill(0);
	beta(2)=feature.get_sig_x();
	this->init_residue(residue,beta);
}

template<class T>
CImgList<float> Match_model_deriv<T>::build_patch_x_y(){
	CImgList<float> res(this->patch,this->patch);

	float sig_x = this->feature->get_sig_x();
	CImg<float> rounded_center=this->feature->get_center().get_round();
	float x_c=this->feature->get_width()/2;
	float y_c=this->feature->get_height()/2;
	float dx =this->feature->get_center()(0)-rounded_center(0);
	float dy =this->feature->get_center()(1)-rounded_center(1);
	int x_r,y_r;

	cimg_forXY(this->patch,x,y){
		x_r=(x-x_c-dx);
		y_r=(y-y_c-dy);
		res[0]*=x_r/(sig_x*sig_x);
		res[1]*=y_r/(sig_x*sig_x);
	}
	return res;
}

template<class T>
void Match_model_deriv<T>::fill_residue() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & residue = this->residue;
	this->update_feature_state();
	this->patch_build();
	CImgList<float> pxy = build_patch_x_y();

	// Here is the main difference, we do not interpolate the values of f
	CImgList<float> G;
	CImg<float> rounded_center=feature.get_center().get_round();
	CImg<float> center_bk=feature.get_center();
	feature.set_center(rounded_center); // "Don't try this at home" 3 lines
	feature.gradient(current,G,this->k,this->k_deriv);
	feature.set_center(center_bk);
	int i = 0;
	cimg_forXY(G[0],x,y){
		residue(i++) = G[0](x,y) + pxy[0](x,y);
		residue(i++) = G[1](x,y) + pxy[1](x,y);
	}
}

template<class T>
void Match_model_deriv<T>::fill_J() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & patch = this->patch ;
	CImg<float> & J = this->J ;
	CImg<float> & kernel = this->k ;
	CImg<float> & kernel_deriv = this->k_deriv;
	CImgList<float> G(2);
	CImgList<float> GG0(2);
	CImgList<float> GG1(2);
	float sig_x = feature.get_sig_x();

	CImg<float> rounded_center=feature.get_center().get_round();
	CImg<float> center_bk=feature.get_center();
	feature.set_center(rounded_center); // "Don't try this at home"  lines

  unsigned int border_size=max(kernel.width(),kernel_deriv.width())/2;
  int w= feature.get_width()+border_size*2;
  int h= feature.get_height()+border_size*2;
	float c0 = feature.get_center()(0);
	float c1 = feature.get_center()(1);
	CImg<float> normal_crop=current.get_crop(c0 - w/2, c1 - h/2,
																	 c0 + w/2, c1 + h/2);


  //Compute gradient
  G[0]=normal_crop ;
	G[0].convolve(kernel_deriv);
  G[0].convolve(kernel.get_transpose());
	GG0[0]=G[0];
	GG0[0].convolve(kernel_deriv);
  GG0[0].convolve(kernel.get_transpose());
	GG0[1]=G[0];
	GG0[1].convolve(kernel);
  GG0[1].convolve(kernel_deriv.get_transpose());

  GG0[0].crop( border_size,
							 border_size,
							 border_size+feature.get_width()-1,
							 border_size+feature.get_height()-1);
  GG0[1].crop( border_size,
							 border_size,
							 border_size+feature.get_width()-1,
							 border_size+feature.get_height()-1);

  G[1]=normal_crop ;
	G[1].convolve(kernel);
  G[1].convolve(kernel_deriv.get_transpose());
	GG1[0]=G[1];
	GG1[0].convolve(kernel_deriv);
  GG1[0].convolve(kernel.get_transpose());
	GG1[1]=G[1];
	GG1[1].convolve(kernel);
  GG1[1].convolve(kernel_deriv.get_transpose());

  GG1[0].crop( border_size,
							 border_size,
							 border_size+feature.get_width()-1,
							 border_size+feature.get_height()-1);
  GG1[1].crop( border_size,
							 border_size,
							 border_size+feature.get_width()-1,
							 border_size+feature.get_height()-1);

	feature.set_center(center_bk);

	int i=0;
	float x_c= feature.get_width()/2;
	float y_c= feature.get_height()/2;
	float sig_x_2=sig_x*sig_x;
	float x_r;
	float y_r;
	float x_r_2;
	float y_r_2;
	float dx = feature.get_center()(0)-rounded_center(0);
	float dy = feature.get_center()(1)-rounded_center(1);
	for ( int x = 0; x < feature.get_width();x++)
		for ( int y = 0; y < feature.get_height();y++) {
		x_r=(x-x_c-dx);
		y_r=(y-y_c-dy);
		x_r_2=x_r*x_r;
		y_r_2=y_r*y_r;
		J(0,i)=-GG0[0](x,y)/2+(patch(x,y))*(1-x_r_2/sig_x_2)/(2*sig_x_2);
		J(1,i)=-GG0[1](x,y)/2+(patch(x,y))*x_r*y_r/(2*sig_x_2*sig_x_2);
		J(2,i)= patch(x,y)*((x_r_2+y_r_2)/sig_x_2-1)*x_r/(sig_x*sig_x*sig_x);
		i++;
		J(0,i)=-GG1[0](x,y)/2+(patch(x,y))*x_r*y_r/(2*sig_x_2*sig_x_2);
		J(1,i)=-GG1[1](x,y)/2+(patch(x,y))*(1-y_r_2/sig_x_2)/(2*sig_x_2);
		J(2,i)= patch(x,y)*((x_r_2+y_r_2)/sig_x_2-1)*y_r/(sig_x*sig_x*sig_x);
		i++;
	}
};
#endif
