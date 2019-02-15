#ifndef _INT_ESTIMATOR_H_
#define _INT_ESTIMATOR_H_

#include "CImg.h"
#include "spot.hpp"
#include "feature.hpp"
#include "gaussian_fit.hpp"
#include "match_model_test.hpp"

template<class T>
class Int_estimator
{
public:

	inline Int_estimator(Spot<T> * s, CImg<T> * st):sp(s),stack(st){};
	inline void run(){
			CImg<T> intensities(sp->get_move_list().size());
			cimg_forX(intensities,z){
				intensities(z)=get_int(z);
			}
			sp->set_intensities(intensities);
	};
	virtual float get_int(int idx)=0;
	Spot<T> * sp;
	CImg<T> * stack;


	static  void intensities_from_move(Spot<T> * f, CImg<T> & st, int type =2 );

};

template<class T>
class Int_bg:public Int_estimator<T>{
public:
	Int_bg(Spot<T> * s, CImg<T> * st):Int_estimator<T>(s,st){
		bg_int.assign(this->sp->get_move_list().size());
	};
	CImg<T>	 bg_int;
	inline virtual float get_int(int idx){
				this->sp->set_center(this->sp->get_move(idx));
				CImg<T> c = this->stack->get_slice(idx);
				CImg<T> patch=this->sp->get_values(c);
				CImg<float> gauss;
				CImg<float> A;
				T white =1;

				patch.vector();

				gauss.assign(this->sp->get_width(),this->sp->get_height());
				float sig = 1.5;
				gauss.draw_gaussian(gauss.width()/2,gauss.height()/2,sig,&white);
				// gauss.draw_gaussian(gauss.width()/2,gauss.height()/2,this->sp->get_sig_x(),&white);
				gauss.vector();

				A.assign(1,gauss.size());
				A.fill(1);
				A.append(gauss);
				CImg<float> t_A=A.get_transpose();

				CImg<float> res= (t_A*patch).solve(t_A*A);
				bg_int(idx)=res(0);
				return res(1);

	};

};

template<class T>
class Int_gf : public Int_estimator<T>
{
public:
	Int_gf(Spot<T> * s, CImg<T> * st):Int_estimator<T>(s,st){};
	inline virtual float get_int(int idx){
				this->sp->set_center(this->sp->get_move(idx));
				CImg<T> c = this->stack->get_slice(idx);
				Feature<T> f ( *(this->sp));
				gaussian_fit<T> gf;
				gf.set_current(&c); // Warning, taking the address of a temporary
				gf.set_feature(&f);
				gf.init();
				gf.run();
				return gf.get_I_max();
	};
};

template<class T>
class Int_gf_fixe : public Int_estimator<T>
{
public:
	Int_gf_fixe(Spot<T> * s, CImg<T> * st):Int_estimator<T>(s,st){};
	inline virtual float get_int(int idx){
		gaussian_bg_fixe_fit<T> gf;
		Feature<T> f ( *(this->sp));
		CImg<T> c = this->stack->get_slice(idx);
		f.set_center(f.get_move(idx));
		gf.set_current(&c); // Warning, taking the address of a temporary
		gf.set_feature(&f);
		gf.firs_bg_estimator = new bg_min<T>();
		gf.init();
		gf.run();
		return gf.get_I_max();

	};
};

template<class T>
class Int_from_patch : public Int_estimator<T>
{
public:
	Int_from_patch(Spot<T> * s, CImg<T> * st):Int_estimator<T>(s,st){};
	inline virtual float get_int(int idx){
		CImg<float> patch;
		CImg<float> values;
		Feature<T> f ( *(this->sp));
		float res;

		Match_model_no_interp<T> mm;

		CImg<T> c = this->stack->get_slice(idx);
		f.set_center(f.get_move(idx));
		mm.set_current(&c); // Warning, taking the address of a temporary
		mm.set_feature(&f);
		mm.set_frame_idx(idx);
		mm.init();
		patch=mm.comp_patch();
		patch/=mm.get_I_max();
		values=f.get_values_no_interp(c);
		res=values.sum()/patch.sum();

		return res;
	};
};

template<class T>
class Int_from_deriv : public Int_estimator<T>
{
public:
	Int_from_deriv(Spot<T> * s, CImg<T> * st):Int_estimator<T>(s,st){};
	inline virtual float get_int(int idx){
		CImgList<float> patch;
		CImg<float> values;
		Feature<T> f ( *(this->sp));
		CImgList<float> G;
		CImg<T> c = this->stack->get_slice(idx);
		float res;

		Match_model_deriv<T> mm;
		f.set_center(f.get_move(idx));
		mm.set_feature(&f);
		mm.set_current(&c); // Warning, taking the address of a temporary
		mm.set_frame_idx(idx);
		mm.init();
		mm.comp_patch();
		patch=mm.build_patch_x_y();
		patch[0]/=mm.get_I_max();
		patch[1]/=mm.get_I_max();
		CImg<float> rounded_center=f.get_center().get_round();
		CImg<float> center_bk=f.get_center();
		f.set_center(rounded_center); // "Don't try this at home" 3 lines
		f.gradient(c,G,Gauss_kernel::k,Gauss_kernel::k_deriv);
		f.set_center(center_bk);

		res=(G[0].sum()/patch[0].sum() + G[1].sum()/patch[1].sum())/2;
		return res;
	};
};

template<class T>
class Int_from_deriv_no_patch : public Int_estimator<T>
{
public:
	Int_from_deriv_no_patch(Spot<T> * s, CImg<T> * st):Int_estimator<T>(s,st){};
	inline virtual float get_int(int idx){
		CImgList<float> patch;
		CImg<float> values;
		Spot<T> & f = *(this->sp);
		CImgList<float> G;
		CImg<T> c = this->stack->get_slice(idx);
		float res;

		f.set_center(f.get_move(idx));
		CImg<float> rounded_center=f.get_center().get_round();
		CImg<float> center_bk=f.get_center();
		f.set_center(rounded_center); // "Don't try this at home" 3 lines
		f.gradient(c,G,Gauss_kernel::k,Gauss_kernel::k_deriv);
		f.set_center(center_bk);

		res=(G[0].max() + (-G[0]).max()+G[1].max() + (-G[1]).max());
		return res;
	};
};

template <class T>
void Int_estimator<T>::intensities_from_move(Spot<T> * f, CImg<T> & st, int type)
{
		switch(type){
		case 3 :
		{
			Int_gf<T> ie(f,&st);
			ie.run();
			break;
		}
		case 4:
		{
			Int_gf_fixe<T> ie(f,&st);
			ie.run();
			break;
		}
		case 5:
		{
			Int_from_patch<T> ie(f,&st);
			ie.run();
			break;
		}
		case 6:
		{
			Int_from_deriv<T> ie(f,&st);
			ie.run();
			break;
		}
		case 7:
		{
			Int_from_deriv_no_patch<T> ie(f,&st);
			ie.run();
			break;
		}
		case 8:
		{
			Int_bg<T> ie(f,&st);
			ie.run();
			break;
		}
	default:
		{
			f->set_intensities_from_move(st,type);
			break;
		}
	}
}


#endif /* _INT_ESTIMATOR_H_ */
