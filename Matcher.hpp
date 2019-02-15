#ifndef _IMATCHER_H_
#define _IMATCHER_H_
#include "sin_fit.hpp"
#include "feature.hpp"
#include "gauss_kernel.hpp"
#include "CFigure.h"

#define IMAGE_MATCHER 0
#define TEMPLATE_MATCHER 1

struct Matcher_data{
	// Very handy to record data relative to the feature.
	CImg<float> betas_ls;
	CImg<float> betas_irls;
	Evo_visualizer scale_evo;
};

template<class T>
class Matcher :  public Sin_fit
{
public:
	Matcher(){flag_gat=false;};
	inline Matcher( Feature<T> * f,  CImg<T> *aCurrent) :
		current(aCurrent),feature(f) {} ;
	virtual void fill_residue()=0 ;
	virtual void fill_J()=0;
	virtual int run();
	inline virtual void init(){};
	virtual int irls_no_init();
	virtual void update_feature_state()=0;
	inline void set_feature(Feature<T> *  val ){feature=val;};
	inline Feature<T> * get_feature( ){return feature;};
	inline void set_current(CImg<T> * val ){current=val;};
	inline CImg<T> * get_current(){return current;};
	CImgList<float> robust_weigth_analizer(Mult_evo_visualizer & evo );
	inline CImg<float> get_f_values(){return f_values;};
	inline CImgList<float> get_f_grad_values(){return f_grad_values;};
	inline void set_f_values(CImg<float> & val ){f_values=val;};
	inline void set_f_grad_values(CImgList<float> & val ){f_grad_values=val;};
	virtual void store_f_values();
	virtual void restore_f_values();



	bool flag_gat;
	float g,e;
	inline void exec_init_estimator(){init_estimator();};
	int matching_diagnosis(float residue_max, bool verbose=false);
	bool image_matcher;
	bool template_matcher;
	Matcher_data attached_data;

protected :
	inline void fill_spring (){};
	CImg<T> * current;
	Feature<T> * feature;
	CImg<float> f_values;
	CImgList<float> f_grad_values;
	virtual void init_estimator()=0;
	static CImg<float> k;
	static CImg<float> k_deriv;
};

template<class T> CImg<float> Matcher<T>::k=Gauss_kernel::get_k();
template<class T> CImg<float> Matcher<T>::k_deriv=Gauss_kernel::get_k_deriv();

template <class T>
int Matcher<T>::run()
{
	CImg<float> d;
	CImg<float> c;
	float sig_x_bck=feature->get_sig_x();
	feature->set_starting_center(feature->get_center());
	init_estimator();
	Sin_fit::run(d);
	update_feature_state();
	feature->set_residue(residue.magnitude());
	if (get_nb_iter() >= get_max_iter()) feature->set_sig_x(sig_x_bck);
	return get_nb_iter();
}
template <class T>
int Matcher<T>::irls_no_init()
{
	Feature<T> &feature = *(this->feature);
	CImg<float> b;
	float sig_x_bck=feature.get_sig_x();
	feature.set_starting_center(feature.get_center());
	init_estimator();
	Sin_fit::irls_no_init(b);
	update_feature_state();
	if (get_nb_iter() >= get_max_iter()) feature.set_sig_x(sig_x_bck);
	feature.set_residue(this->residue.magnitude());
	return this->get_nb_iter();
}

template<class T>
void Matcher<T>::store_f_values(){
	feature->f_values=get_f_values();
	feature->f_grad_values=get_f_grad_values();
}

template<class T>
void Matcher<T>::restore_f_values(){
	set_f_values(feature->f_values);
	set_f_grad_values(feature->f_grad_values);
}

template< class T>
int Matcher<T>::matching_diagnosis(float residue_max, bool verbose){
	int res = 0;
	if((this->get_residue().median()>residue_max)
		 ||(this->has_diverged() )
			||(feature->is_out(*current))){

		if (verbose){
			cout <<"lost due to  " ;

			if(this->has_diverged() )  			cout << "TOO MUCH ITER.";
			if(this->get_residue().median()>residue_max)
				cout<<"RESIDUE TOO LARGE";
			if(feature->is_out(*current))  								cout << "IS OUT";
		}
		res = 0;

	}else{
		if(verbose) cout << "matched";
		res=1;
	}
	return res;
}
template<class T>
CImgList<float> Matcher<T>::robust_weigth_analizer(Mult_evo_visualizer & evo ){
	CImgList<float> patch_result;

	cimg_for1(evo.nb_step(),s){
		cfl beta=evo.get_vector(s);
		this->set_beta(beta);
		cfl r = this->get_residue();
		cfl J = this->get_J();
		cfl rj(r.size());
		cfl wrj(r.size());
		cfl rr(r.size());
		cfl w(r.size());
		cfl p(r.size());

		cimg_foroff(r,i) rj[i]=r[i]*J(0,i);
		cimg_foroff(r,i) rr[i]=r[i]*r[i];
		if(this->get_cost_function()){
			cimg_foroff(r,i) w[i]=(this->get_cost_function()->get_cost(r[i]));
			cimg_foroff(r,i) wrj[i]=r[i]*J(0,i)*w[i];
			cimg_foroff(r,i) p[i]=this->get_cost_function()->get_rho(r[i]);
		}

		w.assign(w._data,this->get_feature()->get_width(),this->get_feature()->get_height());
		rj.assign(rj._data,this->get_feature()->get_width(),this->get_feature()->get_height());
		wrj.assign(wrj._data,this->get_feature()->get_width(),this->get_feature()->get_height());
		rr.assign(rr._data,this->get_feature()->get_width(),this->get_feature()->get_height());
		r.assign(r._data,this->get_feature()->get_width(),this->get_feature()->get_height());
		p.assign(p._data,this->get_feature()->get_width(),this->get_feature()->get_height());

		patch_result.push_back(r.
													 append(w.get_normalize(r.min(),r.max())).
													 append(rr.get_normalize(r.min(),r.max())).
													 append(p).
													 append(rj).
													 append(wrj)
													 );
	}

	return patch_result;

}
#endif /* _IMATCHER_H_ */
