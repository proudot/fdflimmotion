/*
 * =====================================================================================
 *
 *       Filename:  sin_fit_bg_exp.hpp
 *
 *    Description:  toolbox to fit parameter to a a sinusoid
 *
 *        Version:  1.0
 *        Created:  01/05/11 11:54:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  , 
 *
 * =====================================================================================
 */

#ifndef SIN_FIT_BG_EXP_H

#define SIN_FIT_BG_EXP_H



#include	"CImg.h"
#include "sin_fit.hpp"

using namespace cimg_library;
using namespace std;

template <class T>
class Sin_fit_bg_exp: public Sin_fit
{
public:

	/* ====================  LIFECYCLE     ======================================= */
	Sin_fit_bg_exp(Spot<T> * sp):spot(sp){};

	void init ( 	CImg<float> & input_X,
								CImg<float> & input_Y
								); 
	inline void set_starting_value (){};
	void set_starting_value (CImg<T> & stack );
	
	
	inline CImg<float> get_beta(){return beta;};
	inline float get_amp_bg(){return beta(4);};
	inline float get_phase_bg(){return beta(5);};
	inline float get_C_bg(){return beta(3);};
	inline void set_amp_bg(float a){ beta(4)=a;};
	inline void set_phase_bg(float p){ beta(2)=p;};
	inline void set_C_bg(float c){ beta(3)=c;};
	inline float get_bg(int idx){
		float Pi=3.141592;
		float w=2.*Pi/12.;
		return get_C_bg()+get_amp_bg()*(w*idx+get_phase_bg());
	}

	
	void fill_J ();
	void fill_residue ();
	CImg<float> get_deriv_patch(int beta_idx,int t);


protected:
	CImg<float> get_deriv_C_bg(int t);
	CImg<float> get_deriv_amp_bg(int t);
	CImg<float> get_deriv_phase_bg(int t);
	CImg<float> get_deriv_C(int t);
	CImg<float> get_deriv_amp(int t);
	CImg<float> get_deriv_phase(int t);

	/* ====================  DATA MEMBERS  ======================================= */
	// CImg<float> bg_values;				// bg_values consider as constant on the feature
	CImg<float> exp_patch;
	Spot<T> * spot;

};

// Sin_fit_bg_exp<T>::Sin_fit_bg_exp<T>(){};

template <class T>
void Sin_fit_bg_exp<T>::init ( 	CImg<float> & input_X,
															CImg<float> & input_Y
															)
{
	float white = 1;

	exp_patch.assign(spot->get_height(),spot->get_width());
	float sig = 1.5;
	exp_patch.draw_gaussian(exp_patch.width()/2,exp_patch.height()/2,sig,&white);

	beta.assign(6);

	Sin_fit::init(input_X,input_Y,beta);
	set_starting_value();
};		/* -----  end of method Sin_fit_bg_exp<T>::Sin_fit_bg_exp<T>  ----- */

template <class T>
void Sin_fit_bg_exp<T>::fill_J ()
{
	CImg<float> deriv_patch;
	cimg_forY(beta,j){
		int i = 0;
		for (unsigned int  t = 0 ; t < spot->get_move_list().size(); ++t)
			{
				deriv_patch=get_deriv_patch(j,t);
				cimg_forXY(deriv_patch,x,y){
					J(j,i)=deriv_patch(x,y);
					i++;
				}	
			}
	}

}

template <class T>
CImg<float> Sin_fit_bg_exp<T>::get_deriv_patch(int beta_idx,int t){
	CImg<float> res;
	switch(beta_idx){
	case 0:
		res=get_deriv_C(t);
		break;
	case 1:
		res=get_deriv_amp(t);
		break;
	case 2:
		res=get_deriv_phase(t);
		break;
	case 3:
		res=get_deriv_C_bg(t);
		break;
	case 4:
		res=get_deriv_amp_bg(t);
		break;
	case 5:
		res=get_deriv_phase_bg(t);
		break;
	default:
		cout <<"error in get_deriv_patch()" << endl;
		break;
	}
	return res;
};
	
template <class T>
CImg<float> Sin_fit_bg_exp<T>::get_deriv_C_bg(int t){
	return exp_patch.get_fill(1);
}
template <class T>
CImg<float> Sin_fit_bg_exp<T>::get_deriv_amp_bg(int t){
		float Pi=3.141592;
		float w=2.*Pi/12.;
		return exp_patch.get_fill( sin(w*t+get_phase_bg()));
}
template <class T>
CImg<float> Sin_fit_bg_exp<T>::get_deriv_phase_bg(int t){
		float Pi=3.141592;
		float w=2.*Pi/12.;
		return exp_patch.get_fill( get_amp_bg()*cos(w*t+get_phase_bg()));
}
template <class T>
CImg<float> Sin_fit_bg_exp<T>::get_deriv_C(int t){
	return exp_patch;
}
template <class T>
CImg<float> Sin_fit_bg_exp<T>::get_deriv_amp(int t){
		float Pi=3.141592;
		float w=2.*Pi/12.;
		return exp_patch*sin(w*t+get_phase());
}
template <class T>
CImg<float> Sin_fit_bg_exp<T>::get_deriv_phase(int t){
		float Pi=3.141592;
		float w=2.*Pi/12.;
		return exp_patch*get_amp()*cos(w*t+get_phase());
}

template <class T>
void Sin_fit_bg_exp<T>::fill_residue ()
{
		
	int i = 0;
	for (unsigned int  t = 0 ; t < spot->get_move_list().size(); ++t)
		{
			cimg_forXY(exp_patch,x,y){
				residue(i)=Y(i)-(exp_patch(x,y)*(reg_funct(X(i)))+get_bg(t));;
				i++;
				
			}
		}

}		/* -----  end of method Sin_fit_bg_exp<T>::build_delta_beta  ----- */


template <class T>
void Sin_fit_bg_exp<T>::set_starting_value ( CImg<T> & stack )
{
	float Pi=3.141592;
	float w=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
	CImg<T> I_max(spot->get_move_list().size());
	CImg<T> I_min(spot->get_move_list().size());
	CImg<T> val;
	for (unsigned int  t = 0 ; t < spot->get_move_list().size(); ++t)
		{
			spot->set_center(spot->get_move(t));
			val=spot->get_values(stack.get_slice(t));
			I_max(t)=val.max();
			I_min(t)=val.min();
		}
	set_C(I_max.mean());
	set_amp((I_max.max()-I_max.min())/2);	/* FIXME */
	set_phase(atan(2.5*w)+2.42342-atan(4.1*w));	/* FIXME */
	set_C_bg(I_min.mean());
	set_amp_bg((I_min.max()-I_min.min())/2);
	set_phase_bg(get_phase());
}		/* -----  end of method Sin_fit_bg_exp<T>::starting_value  ----- */


#endif /* end of include guard: SIN_FIT_BG_EXP_H */ 
