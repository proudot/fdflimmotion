#ifndef _GAUSSIAN_FIT_H_
#define _GAUSSIAN_FIT_H_


#include "Matcher.hpp"
#include "bg_estimator.hpp"


template<class T>
class gaussian_fit: public Matcher<T> {
public :
	gaussian_fit(){
		reset_I_max=true;
		patch.assign(1).fill(0);};
	virtual void fill_residue();
	virtual void fill_J();
	virtual void update_feature_state();
	virtual  void patch_build();
	virtual	inline float get_I_max(){return this->beta(3);};
	virtual	inline float get_init_I_max(){return this->I_max;};
	inline void set_I_max(float val ){I_max=val;};
	virtual inline void fill_deriv_I_max(int i, int x , int y){ this->J(3,i)= this->patch(x,y)/this->beta(3); };
	virtual inline void fill_deriv_sigma(int i, int x , int y){
		float sig_x=this->feature->get_sig_x();
		float x_c= this->feature->get_width()/2;
		float y_c= this->feature->get_height()/2;
		this->J(2,i)= this->patch(x,y)*((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c)) /(sig_x*sig_x*sig_x);};
	inline CImg<float> get_patch(){return patch;};


  virtual	void init();


	bool reset_I_max;
protected :
	CImg<float> patch;
	float I_max;
	virtual void init_estimator();
};


template<class T>
void gaussian_fit<T>::init(){
	Feature<T> &feature = *(this->feature);
	if(reset_I_max) I_max=(this->current)->cubic_atXY(feature.get_center()(0),feature.get_center()(1));
	// int width=feature.get_width();
	// int height=feature.get_height();
	// CImg<T> slice=this->current->crop(feature.get_center()(0)-width/2,feature.get_center()(1)-height/2,
	// 					 feature.get_center()(0)+width/2,feature.get_center()(1)+height/2);
	// I_max=slice.max();
}

template<class T>
void gaussian_fit<T>::init_estimator(){
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	Feature<T> &feature = *(this->feature);

	residue.assign(feature.get_width()*feature.get_height());
	beta.assign(4).fill(0);
	beta(2)=feature.get_sig_x();
	beta(3)=I_max;
	this->init_residue(residue,beta);
}

template<class T>
void gaussian_fit<T>::update_feature_state(){
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
void gaussian_fit<T>::patch_build(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & patch = this->patch;
	patch.assign(feature.get_width(),
									feature.get_height());
	float color_current=get_I_max();
	CImg<float> rounded_center=feature.get_center().get_round();
	float dx = feature.get_center()(0)-rounded_center(0);
	float dy = feature.get_center()(1)-rounded_center(1);
	patch.draw_gaussian(patch.width()/2+dx,patch.height()/2+dy,feature.get_sig_x(),&color_current);
}

template<class T>
void gaussian_fit<T>::fill_residue() {
	Feature<T> &feature = *(this->feature);
	CImg<T> &current = *(this->current);
	CImg<float> & residue = this->residue;
	this->update_feature_state();
	this->patch_build();

	// Here is a trick, we avoid  interpolating the values of the current image
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
		residue(i++) =
			f_values_current(x,y) - this->patch(x,y);
	}
}

template<class T>
void gaussian_fit<T>::fill_J() {
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
		fill_deriv_sigma(i,x,y);
		fill_deriv_I_max(i,x,y);
	}
}




template<class T>
class gaussian_bg_fixe_fit: public gaussian_fit<T> {
public :
	inline ~gaussian_bg_fixe_fit(){delete firs_bg_estimator;};
	virtual inline  float get_bg_value(){return bg_value;};
	virtual	inline  void set_bg_value(float val ){bg_value=val;};
	virtual inline void fill_deriv_I_max(int i, int x , int y){ this->J(3,i)= (this->patch(x,y)-get_bg_value())/this->get_I_max();};
	inline void fill_deriv_sigma(int i, int x , int y){
		float sig_x=this->feature->get_sig_x();
		float x_c= this->feature->get_width()/2;
		float y_c= this->feature->get_height()/2;
		this->J(2,i)= (this->patch(x,y)-get_bg_value())*((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c)) /(sig_x*sig_x*sig_x);};
	bg_estimator<T> * firs_bg_estimator;
	void patch_build();
	void init();
	virtual void store_f_values();
	virtual void restore_f_values();
protected :
	float bg_value;
};

template<class T>
void gaussian_bg_fixe_fit<T>::store_f_values(){
	float bg=get_bg_value();
	float I=this->get_I_max();
	this->feature->f_values.assign(2,1,1,1,bg,I);
}

template<class T>
void gaussian_bg_fixe_fit<T>::restore_f_values(){
	bg_value=this->feature->f_values(0);
	this->I_max=this->feature->f_values(1);
}
template<class T>
void gaussian_bg_fixe_fit<T>::init(){
	Feature<T> &feature = *(this->feature);
	bg_value=firs_bg_estimator->compute((*this->current),
																					feature.get_center()(0)-feature.get_width()/2,
																					feature.get_center()(1)-feature.get_height()/2,
																					feature.get_center()(0)+feature.get_width()/2,
																					feature.get_center()(1)+feature.get_height()/2
																					);
	if (this->reset_I_max){
		this->I_max=(this->current)->cubic_atXY(feature.get_center()(0),feature.get_center()(1)) ;
	  this->I_max-=bg_value;
	}
}


template<class T>
void gaussian_bg_fixe_fit<T>::patch_build(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & patch = this->patch;
	this->patch.assign(feature.get_width(),feature.get_height());
	float color_current=this->get_I_max();
	CImg<float> rounded_center=feature.get_center().get_round();
	float dx = feature.get_center()(0)-rounded_center(0);
	float dy = feature.get_center()(1)-rounded_center(1);

	patch.draw_gaussian(patch.width()/2+dx,patch.height()/2+dy,feature.get_sig_x(),&color_current);
	this->patch+=get_bg_value();
}

template<class T>
class gaussian_bg_fixe_I_fixe_fit: public gaussian_bg_fixe_fit<T> {
public :
	gaussian_bg_fixe_I_fixe_fit(){
		this->reset_I_max=true;
	};
	inline float get_I_max(){return this->I_max;};
	inline void fill_deriv_I_max(int i, int x , int y){};
protected :
	void init_estimator();
};


template<class T>
void gaussian_bg_fixe_I_fixe_fit<T>::init_estimator(){
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	Feature<T> &feature = *(this->feature);

	residue.assign(feature.get_width()*feature.get_height());
	beta.assign(3).fill(0);
	beta(2)=feature.get_sig_x();
	this->init_residue(residue,beta);
}


template<class T>
class gaussian_bg_fixe_sigma_fixe_fit: public gaussian_bg_fixe_fit<T> {
public :
	gaussian_bg_fixe_sigma_fixe_fit(){
		this->reset_I_max=true;
		this->patch.assign(1).fill(0);
	};
	virtual void update_feature_state();
	inline void fill_deriv_I_max(int i, int x , int y){
		this->J(2,i)= (this->patch(x,y)-this->get_bg_value())/this->get_I_max(); };
	inline void fill_deriv_sigma(int i, int x , int y){ };
	inline float get_I_max(){return this->beta(2);};
protected :
	void init_estimator();
};

template<class T>
void gaussian_bg_fixe_sigma_fixe_fit<T>::update_feature_state(){
	Feature<T> &feature = *(this->feature);
	CImg<float> & beta = this->beta;
	CImg<float> d(1,2);
	d(0)=beta(0);
	d(1)=beta(1);
	CImg<float> c = feature.get_starting_center()+d;
	feature.set_center(c);
}

template<class T>
void gaussian_bg_fixe_sigma_fixe_fit<T>::init_estimator(){
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	Feature<T> &feature = *(this->feature);

	residue.assign(feature.get_width()*feature.get_height());
	beta.assign(3).fill(0);
	beta(2)=this->I_max;
	this->init_residue(residue,beta);
}



template<class T>
class gaussian_bg_fit: public gaussian_bg_fixe_fit<T> {
public :
	gaussian_bg_fit(){
		this->reset_I_max=true;
		this->patch.assign(1).fill(0);
	};
	inline float get_bg_value(){return this->beta(4);};
protected :
	void init_estimator();
};

template<class T>
void gaussian_bg_fit<T>::init_estimator(){
	CImg<float> & residue = this->residue;
	CImg<float> & beta = this->beta;
	Feature<T> &feature = *(this->feature);

	residue.assign(feature.get_width()*feature.get_height());
	beta.assign(5).fill(0);
	beta(2)=feature.get_sig_x();
	beta(3)=this->I_max;
	beta(4)=this->bg_value;
	this->init_residue(residue,beta);
}
#endif /* _GAUSSIAN_FIT_H_ */
