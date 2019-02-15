#ifndef _FEATURE_H
#define _FEATURE_H

#include "CImg.h"
#include "spot.hpp"

#include <iostream>
#include "math.h"
#include <assert.h>
#include <vector>
#include "sin_fit.hpp"


using namespace cimg_library;
using namespace std;



struct Feature_data{
	// this class represents data that are relative to a feature and are deleted with a feature.
	// Though it does not characterized the feature, nor are vital for cycle of the feature.
	// Very handy to record data relative to the feature.
	CImg<float> betas_ls; 						// parameter estimate relative to the feature.
	CImg<float> betas_irls; 						// parameter estimate relative to the feature.
	
};

template<class T>
class Feature : public Spot<T>{

private:
  

	//starting center for each sequence of iteration.
  CImg<float> starting_center;
  double residue;
	bool lost;
	int lost_count;


public:
	int lost_idx;


  CImg<float> eigen_values;
  CImg<float> G;
  CImg<float> e;
  float laplacian;
  float hessian;
	T color_index;
	CImg<float> f_values;	// it can be usefull to save this values for various purposes
	CImgList<float> f_grad_values; // though it should not be automatic for obvious memory mgt reason
  Feature_data attached_data;

	



  Feature() ; 
  Feature(float x, float y, int w_s,bool) ; 
	Feature(const Feature<T> & f);
	Feature(const Spot<T> & f);

  inline void set_starting_center(const CImg<float>& img){ starting_center=img;}
  inline CImg<float>  get_starting_center(){ return starting_center;}
  inline void set_residue(float r){residue=r;}
  inline float get_min_eigen(){return eigen_values.min();}
	inline bool get_lost() const {return lost;}
	inline int get_lost_count(){return lost_count;}
	inline void incr_lost_count(){ lost_count++; }
	inline void reset_lost_count(){ lost_count=0; }
	inline void set_lost(bool value){ lost=value;}
  inline void move(const CImg<float>& aVector) { Spot<T>::center+=aVector; }
  inline float get_det() {return G.det();}
  inline float get_hessian() {return hessian;}
  double get_residue();

  void interpolate_window(const CImg<T>& current,const CImg<T>& next);
  bool is_out(const CImg<T> & img) ; 
  void set_frame_move() ; 
	CImg<float> track(  const CImg<T>& current,
											const CImg<T>& next,
											const CImg<float>& feature_values,
											const CImgList<float>& feature_grad,
											CImg<float>& kernel,
											CImg<float>& kernel_deriv,
											bool lighting_insensitive); 


	void get_interp_and_gradient(	const CImg<T>& img,
																CImg<float>& interp_window ,
																CImgList<float>& grad,
																CImg<float>& kernel,
																CImg<float>& kernel_deriv);
	void get_gradient_and_interp(	const CImg<T>& img,
																CImg<float>& interp_window ,
																CImgList<float>& grad,
																CImg<float>& kernel,
																CImg<float>& kernel_deriv);
	void get_interp(CImg<T> & img,CImg<float> & interp_window );


};
// END OF HEADER

// This function regroup all the needed interpolation computation.
// Indeed  we need : 
//  * the interpolated value of the feature window
//  * interpolated values in the neighborhood to compute gradient.
// 
// Then we compute the gradient immediatly to get rid of unuseful data and
// keeping only the feature in-window values.
// Gradient computation is done with two kernel of one dimension k,k_deriv:
// It's fitted for the derivative of the laplacien but also for classical
// gaussian.
template<class T>
CImg<float> Feature<T>::track( const CImg<T>& current,
															 const CImg<T>& next,
															 const CImg<float>& feature_values,
															 const CImgList<float>& feature_grad,
															 CImg<float>& kernel,
															 CImg<float>& kernel_deriv,
															 bool lighting_insensitive)
{

	CImg<float> matched_values;
  CImgList<float> grad;
  float diff;
	float grad_sum_x=0.;
	float grad_sum_y=0.;
	float alpha=1.;
	float alpha_grad=1.;
	float beta=0.;
	float var_cur=0;
	float var_next=0;
	float mean_cur=0;
	float mean_next=0;

  //members init
  G.fill(0.);
  e.fill(0.);

  get_gradient_and_interp(next,
													matched_values,
													grad,
													kernel,
													kernel_deriv);

// 	if(lighting_insensitive==1){
// 		mean_sq=matched_values.get_mul(matched_values).mean();
// 		mean_sq=sqrt(feature_values.get_mul(feature_values).mean()/mean_sq);
// 		alpha=mean_sq; 
// 		beta= feature_values.mean()-matched_values.mean()*mean_sq;
// 		alpha_grad=sqrt(feature_values.mean()/matched_values.mean());
// 	}
	if(lighting_insensitive){
		mean_next=matched_values.mean();
		mean_cur=feature_values.mean();
		cimg_forXY(matched_values,x,y){
			var_next+=(matched_values(x,y)-mean_next)*(matched_values(x,y)-mean_next);
			var_cur+=(feature_values(x,y)-mean_cur)*(feature_values(x,y)-mean_cur);
		}
		alpha=sqrt(var_cur/var_next); 
		beta= mean_cur - mean_next*alpha;
		alpha_grad=alpha;
	}else{
		// mean_sq=matched_values.get_mul(matched_values).mean();
		// mean_sq=sqrt(feature_values.get_mul(feature_values).mean()/mean_sq);
		// alpha=mean_sq; 
		// beta= feature_values.mean()-matched_values.mean()*mean_sq;
		// alpha_grad=sqrt(feature_values.mean()/matched_values.mean());
	}
		

	residue=0;
  cimg_forXY(matched_values,x,y){
		grad_sum_x=(feature_grad[0](x,y)+alpha_grad*grad[0](x,y)); 
		grad_sum_y=(feature_grad[1](x,y)+alpha_grad*grad[1](x,y));
		// Theory tell us to divide grad's sum by a factor 2, but
		// standford's implementation does not recommand it to "avoid
		// overshooting".
		
		G(0,0)+=grad_sum_x*grad_sum_x;
		G(1,0)+=grad_sum_x*grad_sum_y; 
		G(0,1)+=grad_sum_x*grad_sum_y;
		G(1,1)+=grad_sum_y*grad_sum_y;

		diff=feature_values(x,y)-alpha*matched_values(x,y)-beta;
		residue+=std::fabs(feature_values(x,y) - alpha*matched_values(x,y) - beta);
		e(0)+=grad_sum_x*diff;
		e(1)+=grad_sum_y*diff; 
  }
	
	CImg<float> d=e.solve(G);
  return d;
}




template<class T>
void Feature<T>::get_interp(CImg<T> & img,CImg<float> & interp_window ){
	float pos_x=0.;
	float pos_y=0.;
	interp_window.assign(this->width,this->height);
  cimg_forXY(interp_window,x,y){
		pos_x=Spot<T>::center(0)+float(x-this->width/2);
		pos_y=Spot<T>::center(1)+float(y-this->height/2);
		interp_window(x,y)=img.linear_atXY(pos_x,pos_y);
  }
}

template<class T>
void  Feature<T>::get_gradient_and_interp(const CImg<T>& img,
																					CImg<float>& interp_window ,
																					CImgList<float>& grad,
																					CImg<float>& kernel,
																					CImg<float>& kernel_deriv){
  CImg<float> normal_crop;
	CImgList<float> grad_tmp(2);
	float pos_x=0.;
	float pos_y=0.;
	int w_s = this->width;
	//could be done in an other function, let here to keep the simetry with the 
	//get_interp_and_gradient function
	interp_window.assign(this->width,this->height);
  cimg_forXY(interp_window,x,y){
		pos_x=Spot<T>::center(0)+float(x-w_s/2);
		pos_y=Spot<T>::center(1)+float(y-w_s/2);
		interp_window(x,y)=img.linear_atXY(pos_x,pos_y);
  }


  unsigned int border_size=max(kernel.width(),kernel_deriv.width())/2;
  int w= this->width+border_size*2;
  int h= this->height+border_size*2;
	normal_crop=img.get_crop((int)Spot<T>::center(0) - w/2, (int)Spot<T>::center(1) - h/2,
													 (int)Spot<T>::center(0) + w/2, (int)Spot<T>::center(1) + h/2);

	

  //Compute gradient
	//interpolated_crop.display();
  grad_tmp[0]=normal_crop ;
	grad_tmp[0].convolve(kernel_deriv);
  grad_tmp[0].convolve(kernel.get_transpose());

  grad_tmp[1]=normal_crop ;
	grad_tmp[1].convolve(kernel);
  grad_tmp[1].convolve(kernel_deriv.get_transpose());

  grad_tmp[0].crop( border_size,
										border_size,
										border_size+this->width,
										border_size+this->height);
  grad_tmp[1].crop( border_size,
										border_size,
										border_size+this->width,
										border_size+this->height);

	float float_diff_x=(int)Spot<T>::center(0)-Spot<T>::center(0);
	float float_diff_y=(int)Spot<T>::center(1)-Spot<T>::center(1);
	grad.assign(grad_tmp[0],grad_tmp[0]);
  cimg_forXY(grad[0],x,y){
		pos_x=float(x)+float_diff_x;
		pos_y=float(y)+float_diff_y;
		grad[0](x,y)=grad_tmp[0].linear_atXY(pos_x,pos_y);
		grad[1](x,y)=grad_tmp[1].linear_atXY(pos_x,pos_y);
  }
}


template<class T>
void  Feature<T>::get_interp_and_gradient(const CImg<T>& img,
																								CImg<float>& interp_window ,
																								CImgList<float>& grad,
																								CImg<float>& kernel,
																								CImg<float>& kernel_deriv){
	grad.assign(2);

  unsigned int border_size=max(kernel.width(),kernel_deriv.width())/2;
  CImg<float> interpolated_crop( this->width+border_size*2,
                                  this->height+border_size*2);
  
  int w=interpolated_crop.width();
  int h=interpolated_crop.height();
	float pos_x=0.;
	float pos_y=0.;

  cimg_forXY(interpolated_crop,x,y){
		pos_x=Spot<T>::center(0)+float(x-h/2);
		pos_y=Spot<T>::center(1)+float(y-w/2);
		interpolated_crop(x,y)=img.linear_atXY(pos_x,pos_y);
  }
	
	//crop interpolated value for other manipulation (diff, residue ...)
  interp_window.assign(interpolated_crop.get_crop(  border_size,
                                              border_size,
                                              border_size+this->width,
                                              border_size+this->height));
  //Compute gradient
	//interpolated_crop.display();
  grad[0]=interpolated_crop ;
	grad[0].convolve(kernel_deriv);
  grad[0].convolve(kernel.get_transpose());
  grad[0].crop( border_size,
                border_size,
                border_size+this->width,
                border_size+this->height);

  grad[1]=interpolated_crop ;
	grad[1].convolve(kernel);
  grad[1].convolve(kernel_deriv.get_transpose());
  grad[1].crop( border_size,
                border_size,
                border_size+this->width,
                border_size+this->height);


}

template<class T>
double Feature<T>::get_residue() { 
  return residue;
  } 

template<class T>
bool Feature<T>::is_out(const CImg<T> & img) { 
  return (Spot<T>::center(0)<this->width/2)
		||(Spot<T>::center(1)<this->height/2)
		||(Spot<T>::center(1)>=(img.height()-this->height/2))
		||(Spot<T>::center(0)>=(img.width()-this->width/2));
  } 

template<class T>
void Feature<T>::set_frame_move() { 
	this->move_list.push_back(this->center);
} 


// template<class T>
// Feature<T>::Feature(const Spot<T> & f):Spot<T>(f)  { 
//   G=CImg<float>(2,2,1,1,1);
//   e=CImg<float>(1,2,1,1,1);
// 	residue=0;
// 	lost=false;
// 	lost_count=0; 
// }

template<class T>
Feature<T>::Feature(const Feature<T> & f):Spot<T>(f)  { 
	starting_center=f.starting_center;
	residue=f.residue;
	lost=f.lost;
	G=f.G;
	e=f.e;
	eigen_values=f.eigen_values;
	color_index=f.color_index;
	lost_count=f.lost_count;
	f_values=f.f_values;
	// attached_data=f.attached_data;
}

template<class T>
Feature<T>::Feature(float x, float y, int w_s,bool set_frame_move) : Spot<T>(x,y,set_frame_move){ 
	this->height=w_s;
	this->width=w_s;
  G=CImg<float>(2,2,1,1,1);
  e=CImg<float>(1,2,1,1,1);
	color_index=((int)x*10)*(int)y*10 % 256;
	residue=0;
	lost=false;
	lost_count=0; 
	this->starting_frame_idx=0;
}

template<class T>
Feature<T>::Feature(){residue=0;}

/********** end of class **********/


// template<class T>
// class KLT: public Sin_fit {
// public :
// 	// reuse the feature !
// 	KLT();
// 	inline KLT( Feature<T> & f,  CImg<T> &aCurrent, CImg<T> & aNext) :
// 		current(aCurrent),next(aNext),feature(f) {} ;
// 	virtual void fill_residue() ;
// 	virtual void fill_J();
// 	virtual int run();
// 	virtual int run_other_starting_center();
// 	virtual int irls_other_starting_center();
// 	virtual int irls ( );
// 	virtual int irls_no_init(float th, int m_iter );
// 	inline void set_feature(Feature<T> &  val ){feature=val;};

// protected :
// 	CImg<T> & current;
// 	CImg<T> & next;
// 	Feature<T> & feature;
// 	CImg<float> f_values;
// 	CImgList<float> f_grad_values;
// 	float var_cur;
// 	float var_next;
// 	float mean_cur;
// 	float mean_next;
// 	float bias;
// 	virtual void init();
// 	static CImg<float> k;
// 	static CImg<float> k_deriv;
// 	static void build_kernel(	float sigma,
// 														CImg<float> & my_gaussian_kernel,
// 														CImg<float> & my_gaussian_kernel_deriv);
// 	static inline   CImg<float> get_k (){ 
// 		CImg<float> k1;
// 		CImg<float> k2;
// 		build_kernel(1.,k1,k2);
// 		return k1;
// 	}
// 	static inline   CImg<float> get_k_deriv (){ 
// 		CImg<float> k1;
// 		CImg<float> k2;
// 		build_kernel(1.,k1,k2);
// 		return k2;
// 	}
// };

// template<class T> CImg<float> KLT<T>::k=KLT<T>::get_k();
// template<class T> CImg<float> KLT<T>::k_deriv=KLT<T>::get_k_deriv();



// template<class T>
// void KLT<T>::init(){
// 	feature.get_gradient_and_interp(current,
// 																	f_values,
// 																	f_grad_values,
// 																	k,
// 																	k_deriv );
// 	mean_cur=f_values.mean();
// 	var_cur = 0 ;
// 	cimg_forXY(f_values,x,y){
// 		var_cur+=(f_values(x,y)-mean_cur)*(f_values(x,y)-mean_cur);
// 	}
// 	residue.assign(feature.get_width()*feature.get_height());
// 	beta.assign(2).fill(0);
// 	init_residue(residue,beta);
// }

// template <class T>
// int KLT<T>::irls (  )
// {
// 	CImg<float> d;
// 	feature.set_starting_center(feature.get_center());
// 	init();
// 	varEstimator<T> ve;
// 	float err_sig = ve.estimate(current);
// 	cout << "estimated sig : " << err_sig << endl;
	
// 	Sin_fit::irls(d, err_sig, sqrt(2),	&dcost_leclerc);
// 	CImg<float> c = feature.get_starting_center()+d;
// 	feature.set_center(c);
// 	feature.set_residue(residue.magnitude());
// 	return get_nb_iter();
// }

// template <class T>
// int KLT<T>::run()
// {
// 	CImg<float> d;
// 	CImg<float> c; 
// 	feature.set_starting_center(feature.get_center());
// 	init();
// 	Sin_fit::run(d);
// 	c = feature.get_starting_center()+d;
// 	feature.set_center(c);
// 	feature.set_residue(residue.magnitude());
// 	return get_nb_iter();
// }
// template <class T>
// int KLT<T>::irls_no_init(float th, int m_iter )
// {
// 	Feature<T> &feature = this->feature;
// 	static ofstream log("log_KLT_flim.log");
// 	CImg<float> b;
// 	feature.set_starting_center(feature.get_center());
// 	init();
// 	this->set_threshold(th);
// 	this->set_max_iter(m_iter); 
// 	Sin_fit::irls_no_init(b, sqrt(2),	&dcost_leclerc);
// 	this->beta=b; 
// 	CImg<float> c = feature.get_starting_center()+b;
// 	feature.set_center(c);
// 	feature.set_residue(this->residue.magnitude());
// 	feature.print_all(log); 
// 	return this->get_nb_iter();
// }
// template <class T>
// int KLT<T>::irls_other_starting_center()
// {
// 	CImg<float> d;
// 	CImg<float> c; 
// 	init();
// 	varEstimator<T> ve;
// 	float err_sig = ve.estimate(current);
// 	Sin_fit::irls(d, err_sig, 3,	&dcost_leclerc);
// 	c = feature.get_starting_center()+d;
// 	feature.set_center(c);
// 	feature.set_residue(residue.magnitude());
// 	return get_nb_iter();
// }
// template <class T>
// int KLT<T>::run_other_starting_center()
// {
// 	CImg<float> d;
// 	CImg<float> c; 
// 	init();
// 	Sin_fit::run(d);
// 	c = feature.get_starting_center()+d;
// 	feature.set_center(c);
// 	feature.set_residue(residue.magnitude());
// 	return get_nb_iter();
// }

// template<class T>
// void KLT<T>::fill_residue() { 
// 	bias = 1;
// 	float gain = 0;

// 	CImg<float> f_values_next;
// 	CImg<float> c = feature.get_starting_center()+beta;
// 	feature.set_center(c);
// 	feature.interpolate(next,f_values_next);
// 	mean_next=f_values_next.mean();
// 	var_next = 0 ;
// 	cimg_forXY(f_values_next,x,y){
// 		var_next+=(f_values_next(x,y)-mean_next)*(f_values_next(x,y)-mean_next);
// 	}

// 	bias=sqrt(var_cur/var_next);
// 	gain= mean_cur - mean_next*bias;
// 	float i = 0.;
// 	cimg_forXY(f_values,x,y){
// 		residue(i++) = f_values(x,y) - bias*f_values_next(x,y) - gain ;
// 	}
// };

// template<class T>
// void KLT<T>::fill_J() { 
// 	CImgList<float> G;
// 	feature.gradient(next,G,k,k_deriv);
// 	int x= 0;
// 	int y=0;
// 	cimg_forY(J,i){
// 		x = i % feature.get_width();
// 		y = i / feature.get_width();
// 		J(0,i)=(f_grad_values[0](x,y)+bias*G[0](x,y))/2.;
// 		J(1,i)=(f_grad_values[1](x,y)+bias*G[1](x,y))/2.;
// 	}
		
// };

// template<class T>
// void			KLT<T>::build_kernel(	float sigma,
// 																CImg<float> & my_gaussian_kernel,
// 																CImg<float> & my_gaussian_kernel_deriv){
// 	T black[1]={1};
// 	//approx of ln(10)
// 	float sqrtln10=1.5174;
// 	float sqrt2=1.4142;

// 	//compute the size to avoid value below 0.01
// 	int gaussian_kernel_size= 2*(int)(2*sigma*sqrtln10)+1; 
// 	int gaussian_kernel_deriv_size= 2*(int)(2*sqrt2*sigma*sqrtln10)+1; 
// 	my_gaussian_kernel.assign(gaussian_kernel_size);
// 	my_gaussian_kernel.draw_gaussian(gaussian_kernel_size/2,sigma,black);
// 	my_gaussian_kernel_deriv.assign(gaussian_kernel_deriv_size);
// 	my_gaussian_kernel_deriv.draw_gaussian(gaussian_kernel_deriv_size/2,sigma,black);

// 	float den=0.;
// 	cimg_forX(my_gaussian_kernel_deriv,x){
// 		my_gaussian_kernel_deriv(x)=(gaussian_kernel_deriv_size/2-x)*my_gaussian_kernel_deriv(x);
// 		den-=(x-gaussian_kernel_deriv_size/2 )*my_gaussian_kernel_deriv(x);
// 	}
// 	my_gaussian_kernel_deriv/=den;
// };



// template<class T>
// class KLT_flim: public KLT<T> {
// public :
// 	// reuse the feature !
// 	KLT_flim();
// 	KLT_flim( Feature<T> & f,  CImg<T> &aCurrent, CImg<T> & aNext,
// 						int aFrame_idx,
// 						float aAmplitude,
// 						float aPhase,
// 						float aW  
// 						) ;
// 	void fill_residue() ;
// 	void fill_J();
// 	int run();
// 	int  irls ( );
// 	void set_starting_value (  );

// protected :
// 	float amplitude;
// 	float w;
// 	float phase;
// 	float frame_idx;
// };

// template<class T>
// void KLT_flim<T>::set_starting_value (  )
// 	{
// 		float Pi=3.141592;
// 		float wt=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
// 		amplitude =2.5*4000;								 // exact amplitude of simulated spot
// 		w=2.*Pi/12.;                                   /* FIXME */
// 		phase=atan(2.5*wt)+2.42342-atan(4.1*wt);	/* FIXME */
// 	}		/* -----  end of method Sin_fit::starting_value  ----- */
// // 1.62303
// template<class T>
// KLT_flim<T>::KLT_flim(Feature<T> & f,
// 											CImg<T> &aCurrent,
// 											CImg<T> & aNext,
// 											int aFrame_idx,
// 											float aAmplitude,
// 											float aPhase=0.,
// 											float aW=0.):
// 	KLT<T>::KLT(f,aCurrent,aNext)
// {
// 	frame_idx=aFrame_idx;
// 	// build_kernel(1.,k,k_deriv);
// 	set_starting_value();
// 	if(aAmplitude != 0 ) amplitude = aAmplitude;
// 	if(aPhase != 0. ) phase = aPhase;
// 	if(aW != 0. ) w = aW;

// 	this->feature.get_gradient_and_interp(this->current,
//  																	KLT<T>::f_values,
// 																	KLT<T>::f_grad_values,
// 																	KLT<T>::k,
// 																	KLT<T>::k_deriv );
// 	KLT<T>::residue.assign(KLT<T>::feature.get_width()*KLT<T>::feature.get_height());
// 	KLT<T>::beta.assign(5).fill(0);
// 	KLT<T>::beta(2)=amplitude;
// 	KLT<T>::beta(3)=w;
// 	KLT<T>::beta(4)=phase;
// 	KLT<T>::init_residue(this->residue,this->beta);
// };

// template <class T>
// int KLT_flim<T>::run()
// {
// 	Feature<T> &feature = this->feature;
// 	ofstream log("log_KLT_flim.log");
// 	CImg<float> b;
// 	CImg<float> d(1,2);
// 	feature.set_starting_center(feature.get_center());
// 	Sin_fit::run(b);
// 	d(0)=b(0);
// 	d(1)=b(1);
// 	feature.set_phase(b(4));
	
// 	CImg<float> c = feature.get_starting_center()+d;
// 	feature.set_center(c);
// 	feature.set_residue(this->residue.magnitude());
// 	// log << "starting center\t: " 
// 	// 		<< feature.get_starting_center()(0) <<" "
// 	// 		<< feature.get_starting_center()(1) << endl;
// 	// log << "matched center \t: " 
// 	// 		<< feature.get_center()(0) << " " <<
// 	// 		feature.get_center()(1) << endl;
// 	// log << "amplitude \t:" << b(2) <<endl;
// 	// log << "frequency \t:" << b(3) <<endl;
// 	// log << "phase \t:"     << b(4) << endl<<endl;
// 	// log.close();
// 	return this->get_nb_iter();
// }
// template <class T>
// int KLT_flim<T>::irls()
// {
// 	Feature<T> &feature = this->feature;
// 	ofstream log("log_KLT_flim.log");
// 	CImg<float> b;
// 	CImg<float> d(1,2);
// 	feature.set_starting_center(feature.get_center());
// 	varEstimator<T> ve;
// 	float err_sig = ve.estimate(this->current);
// 	Sin_fit::irls(b, err_sig, 3,	&dcost_leclerc);
// 	d(0)=b(0);
// 	d(1)=b(1);
// 	feature.set_phase(b(4));
	
// 	CImg<float> c = feature.get_starting_center()+d;
// 	feature.set_center(c);
// 	feature.set_residue(this->residue.magnitude());
// 	log << "starting center\t: " 
// 			<< feature.get_starting_center()(0) <<" "
// 			<< feature.get_starting_center()(1) << endl;
// 	log << "matched center \t: " 
// 			<< feature.get_center()(0) << " " <<
// 			feature.get_center()(1) << endl;
// 	log << "amplitude \t:" << b(2) <<endl;
// 	log << "frequency \t:" << b(3) <<endl;
// 	log << "phase \t:"     << b(4) << endl<<endl;
// 	log.close();
// 	return this->get_nb_iter();
// }
// template<class T>
// void KLT_flim<T>::fill_residue() { 
// 	CImg<float> f_values_next;
// 	CImg<float> c = this->feature.get_starting_center()+this->beta;
// 	this->feature.set_center(c);
// 	this->feature.interpolate(this->next,f_values_next);
// 	CImg<float> &beta = this->beta;
// 	float i = 0.;
// 	cimg_forXY(this->f_values,x,y){
// 		this->residue(i++) = this->f_values(x,y) - f_values_next(x,y) - beta(2)*beta(3)*cos(beta(3)*frame_idx + beta(4)) ;
// 	}
// }
// template<class T>
// void KLT_flim<T>::fill_J() { 
// 	CImgList<float> G;
// 	CImg<float> &beta = this->beta;
// 	CImg<float> &J = this->J;
// 	CImgList<float> & f_grad_values = this->f_grad_values ;


// 	Feature<T> &feature = this->feature;
// 	feature.gradient(this->next,G,this->k,this->k_deriv);
// 	int x= 0;
// 	int y=0;
// 	cimg_forY(J,i){
// 		x = i % feature.get_width();
// 		y = i / feature.get_width();
// 		J(0,i)=(f_grad_values[0](x,y)+G[0](x,y))/2;
// 		J(1,i)=(f_grad_values[1](x,y)+G[1](x,y))/2;
// 		J(2,i)=beta(3)*cos(beta(3)*frame_idx + beta(4)); 
// 		J(3,i)=beta(2)*(cos(beta(3)*frame_idx + beta(4)) - beta(3)*frame_idx*sin(beta(3)*frame_idx+beta(4)));
// 		J(4,i)=-beta(3)*beta(2)*sin(beta(3)*frame_idx + beta(4)); 
// 	}
		
// };



// template<class T>
// class KLT_flim0: public KLT<T> {
// public :
// 	// reuse the feature !
// 	inline KLT_flim0( Feature<T> & f,  CImg<T> &aCurrent, CImg<T> & aNext, int aFrame_idx):
// 		KLT<T>::KLT(f,aCurrent,aNext),frame_idx(aFrame_idx){};
// 	virtual void fill_residue() ;
// 	virtual void fill_J();
// 	virtual void set_starting_value (  );

// protected :
// 	int frame_idx;
// 	CImg<float> local_var;
// 	float color_current;
// 	float color_next;
// 	virtual void init();
// 	float local_value(CImg<T> & img);
// 	void modele_variation();
// };

// template<class T> 
// void KLT_flim0<T>::modele_variation(){
// 	float sig = 1.5 ;
// 	local_var.assign(this->feature.get_width(),
// 									 this->feature.get_height());
// 	float x_c= this->feature.get_width()/2;
// 	float y_c= this->feature.get_height()/2;
// 	float  phase = this->feature.get_phase();
// 	float  C= this->feature.get_C();

// 	float Pi=3.141592;
// 	float w=2.*Pi/12;
// 	float  amplitude= this->feature.get_amp();
// 	color_current=C+amplitude*sin(frame_idx*w+phase);
// 	color_next=C+amplitude*sin((frame_idx+1)*w+phase);
// 	cimg_forXY(local_var,x,y){
// 		local_var(x,y) =	exp(-((x-x_c)*(x-x_c) + (y-y_c)*(y-y_c))/(2*sig*sig))*(color_current - color_next);
// 	}
// }

// template<class T> 
// float KLT_flim0<T>::local_value(CImg<T> & img){
// 	float cx= this->feature.get_center()(0);
// 	float cy= this->feature.get_center()(1);
// 	float local_value;
// 	CImg<float> local_values=img.get_crop(cx-2*this->feature.get_width(),
// 																									cy-2*this->feature.get_height(),
// 																									cx+2*this->feature.get_width(),
// 																									cy+2*this->feature.get_height());
	
// 	// local_values.sort();
// 	// local_values.assign(local_values.data(),2*this->feature.get_width()*this->feature.get_height());
// 	local_value=local_values.median();
// 	return local_value;
// }

// template<class T>
// void KLT_flim0<T>::set_starting_value (  )
// {
// 	float Pi=3.141592;
// 	float wt=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
// 	this->feature.set_C(1.5*4*4000);								 // exact amplitude of simulated spot
// 	this->feature.set_amp(3*4000);								 // exact amplitude of simulated spot
// 	this->feature.set_phase(atan(2.5*wt)+2.42342-atan(4.1*wt));	/* FIXME */
// }		/* -----  end of method Sin_fit::starting_value  ----- */
// template<class T>
// void KLT_flim0<T>::init(){
// 	Feature<T> &feature = this->feature;
// 	CImgList<float> & f_grad_values = this->f_grad_values ;
// 	CImg<float> & f_values = this->f_values ;
// 	CImg<T> & current = this->current ;
// 	// CImg<T> & next = this->next ;
// 	CImg<float> & k = this->k ;
// 	CImg<float> & k_deriv = this->k_deriv;
// 	CImg<float> & residue = this->residue;
// 	CImg<float> & beta = this->beta;

// 	modele_variation();

// 	feature.get_gradient_and_interp(current,
// 																	f_values,
// 																	f_grad_values,
// 																	k,
// 																	k_deriv );
// 	residue.assign(feature.get_width()*feature.get_height());
// 	beta.assign(2).fill(0);
// 	this->init_residue(residue,beta);
// }
// template<class T>
// void KLT_flim0<T>::fill_residue() { 
// 	Feature<T> &feature = this->feature;
// 	CImg<float> & f_values = this->f_values ;
// 	CImg<T> & next = this->next ;
// 	CImg<float> & residue = this->residue;
// 	CImg<float> & beta = this->beta;
	

// 	CImg<float> f_values_next;
// 	CImg<float> c = feature.get_starting_center()+beta;
// 	feature.set_center(c);
// 	feature.interpolate(next,f_values_next);

// 	float i = 0.;
	
// 	cimg_forXY(f_values,x,y){
// 		// residue(i++) = f_values(x,y) - f_values_next(x,y) - local_var(x,y);
// 		residue(i++) = f_values(x,y) - f_values_next(x,y)*color_current/color_next ;
// 	}
// };
// template<class T>
// void KLT_flim0<T>::fill_J() { 
// 	Feature<T> &feature = this->feature;
// 	CImgList<float> & f_grad_values = this->f_grad_values ;
// 	CImg<T> & next = this->next ;
// 	CImg<float> & k = this->k ;
// 	CImg<float> & J = this->J ;
// 	CImg<float> & k_deriv = this->k_deriv;
// 	CImgList<float> G;
// 	feature.gradient(next,G,k,k_deriv);
// 	int x= 0;
// 	int y=0;

// 	float sig = 2 ;
// 	cimg_forY(J,i){
// 		x = i % feature.get_width();
// 		y = i / feature.get_width();
// 		// J(0,i)=(f_grad_values[0](x,y)+G[0](x,y))/2 +local_var(x,y)*(x-feature.get_width()/2)/(2*sig*sig);
// 		// J(1,i)=(f_grad_values[1](x,y)+G[1](x,y))/2.+local_var(x,y)*(y-feature.get_height()/2)/(2*sig*sig);
// 		J(0,i)=(f_grad_values[0](x,y)+G[0](x,y)*color_current/color_next)/2 ;
// 		J(1,i)=(f_grad_values[1](x,y)+G[1](x,y)*color_current/color_next)/2;
// 	}
		
// };


#endif /* _FEATURE_H */

