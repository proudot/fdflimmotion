/*
 * =====================================================================================
 *
 *       Filename:  sin_fit.hpp
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

#ifndef SIN_FIT_H
#define SIN_FIT_H

// #include <varEstimator.hpp>
#include	"CImg.h"
#include	<complex>
#include <fstream>
#include <iostream>
#include "evo_visualizer.hpp"
#include "CFigure.h"

using namespace cimg_library;
using namespace std;

template < class T >
void print_val ( CImg<T> M , ostream & of= cout)
{
	cimg_forY(M,y){
		cimg_forX(M,x){
			of << M(x,y) << "\t\t";
		}
		of <<endl;
	}
}		/* -----  end of template function print_val  ----- */


struct Cost_function {
	float alpha;
	float sig_noise;
	Cost_function(float s,float a):alpha(a),sig_noise(s){};
	Cost_function(){};
	virtual float get_cost(float r)=0;
	virtual float get_rho(float r)=0;
	inline  void set_sig_noise(float s){sig_noise=s;};
	inline  void set_alpha(float a){alpha=a;};
	inline virtual  CFigure get_cost_graph(){
		int sample_nb=1000;
		cfl samples(sample_nb),residues(sample_nb);
		float delta_r= 10*sqrt(alpha/2)*sig_noise/sample_nb;
		float r=0;int i=0;
		for(r = -5*sqrt(alpha/2)*sig_noise, i=0; i < sample_nb; r+= delta_r,i++){ // 
			samples(i)=get_cost(r);
			residues(i)=r;
		}
		CFigure fig; fig.set_axis_mode(1).plot(residues,samples,"b-");
		return fig;
	}

	inline virtual CFigure get_rho_graph(){
		int sample_nb=1000;
		cfl samples(sample_nb),residues(sample_nb);
		float delta_r= 10*sqrt(alpha/2)*sig_noise/sample_nb;
		float r=0;int i=0;
		for(r = -5*sqrt(alpha/2)*sig_noise, i=0; i < sample_nb; r+= delta_r,i++){ // 
			samples(i)=get_rho(r);
			residues(i)=r;
		}
		CFigure fig; fig.set_axis_mode(1).plot(residues,samples,"b-");
		return fig;
	}
};

struct Tukey : public Cost_function{

	Tukey(float sig,float alph):Cost_function(sig,alph){};
	Tukey(){};
	inline float get_cost(float r){
		return max(pow(1-(r/(sig_noise*alpha))*(r/(sig_noise*alpha)),(float)2),(float)0);
		// return (r/sig_noise)*max(pow(1-(r/(sig_noise*alpha))*(r/(sig_noise*alpha)),(float)2),(float)0);
	}
	inline float get_rho(float x){
		cout << "WARNING : tukey rho is not implemented" << endl;
		return x;
	}
};

struct Pseudotukey : public Cost_function{

	Pseudotukey(float sig,float alph):Cost_function(sig,alph){};
	inline float get_cost(float x){
		float res = 0;
		if(x>0){ res = get_cost(x,sig_noise,alpha);}
		else   { res = 1;}
		return res;
	}
	inline float get_rho(float x){
		cout << "WARNING : Pseudotukey rho is not implemented" << endl;
		return x;
	}

	static inline float get_cost(float x,float sig,float alp){
		return max(1-(x/(sig*alp))*(x/(sig*alp)),(float)0);
	}
};

struct Leclerc : public Cost_function{
	Leclerc(float sig,float alph):Cost_function(sig,alph){};
	Leclerc(){};
	inline float get_cost(float x){
		return exp(-x*x/(alpha*sig_noise*sig_noise))/(alpha*sig_noise*sig_noise);
		// return 2*x*exp(-x*x/(alpha*sig_noise*sig_noise))/(alpha*sig_noise*sig_noise);
	}
	inline float get_rho(float x){
		return 1-exp(-x*x/(alpha*sig_noise*sig_noise));
	}

};
struct Leclerc_bis : public Cost_function{
	Leclerc_bis(float sig,float alph):Cost_function(sig,alph){};
	Leclerc_bis(){};
	inline float get_cost(float x){
		return exp(-x*x/(alpha*sig_noise*sig_noise));
		// return 2*x*exp(-x*x/(alpha*sig_noise*sig_noise))/(alpha*sig_noise*sig_noise);
	}
	inline float get_rho(float x){
		return 1-exp(-x*x/(alpha*sig_noise*sig_noise));
	}

};
struct L1 : public Cost_function{
	L1(){};
	inline float get_cost(float x){
		return 1/fabs(x);
	}
	inline float get_rho(float x){
		return fabs(x);
	}

};

struct L1_approx : public Cost_function{

	L1_approx(){epsilon=0.0000001;};
	float epsilon;
	inline float get_cost(float x){
		return x/sqrt(epsilon+pow(x,2.f));
	}
	float get_rho(float x){
		return sqrt(epsilon+pow(x,2.f));
	}

};

struct German_McLure : public Cost_function{

	German_McLure(float sig,float alph):Cost_function(sig,alph){};
	German_McLure(){};
	inline float get_cost(float x){
		return 1/pow(alpha*sig_noise*sig_noise+x*x,float(2));
	}
	inline float get_rho(float x){
		cout << "WARNING : german_mclure rho is not implemented" << endl;
		return x;
	}

};

struct Asym_cos_function : public Cost_function{
	Cost_function * cf;
	float sig_neg;

	inline void set_sig_neg(float val ){sig_neg=val;};
	inline void set_sig_pos(float val ){sig_noise=val;};
	
	Asym_cos_function(float sig_pos,float sig_neg,float alph,	Cost_function * c)
		:Cost_function(sig_pos,alph),cf(c),sig_neg(sig_neg) { };

	inline float get_cost(float x){
		float res = 0 ;
		cf->set_alpha(alpha);
		if(x > 0){
			cf->set_sig_noise(sig_noise);
			res = cf->get_cost(x);
		}
		else{
			cf->set_sig_noise(sig_neg);
			res = cf->get_cost(x);
		}
		return res;
	}
	
	inline float get_rho(float x){
		float res = 0 ;
		cf->set_alpha(alpha);
		if(x > 0){
			cf->set_sig_noise(sig_noise);
			res = cf->get_rho(x);
		}
		else{
			cf->set_sig_noise(sig_neg);
			res = cf->get_rho(x);
		}
		return res;
	}
	inline CFigure get_cost_graph(){
		int sample_nb=1000;
		cfl samples(sample_nb),residues(sample_nb);
		float delta_r= 4*sqrt(alpha/2)*(sig_noise+sig_neg)/sample_nb;
		float r=0;int i=0;
		for(r = -4*sqrt(alpha/2)*(sig_noise+sig_neg)/2, i=0; i < sample_nb; r+= delta_r,i++){ // 
			samples(i)=get_cost(r);
			residues(i)=r;
		}
		CFigure fig; fig.set_axis_mode(1).plot(residues,samples,"b-");
		return fig;
	}

	inline CFigure get_rho_graph(){
		int sample_nb=1000;
		cfl samples(sample_nb),residues(sample_nb);
		float delta_r= 4*sqrt(alpha/2)*(sig_noise+sig_neg)/sample_nb;
		float r=0;int i=0;
		for(r = -4*sqrt(alpha/2)*(sig_noise+sig_neg)/2, i=0; i < sample_nb; r+= delta_r,i++){ // 
			samples(i)=get_rho(r);
			residues(i)=r;
		}
		CFigure fig; fig.set_axis_mode(1).plot(residues,samples,"b-");
		return fig;
	}
};

struct Tukey_asym: public Asym_cos_function{
	Tukey_asym(float sig_pos, float sig_neg, float alp)
		:Asym_cos_function(sig_pos,sig_neg,alp,NULL){
		cf = new Tukey(0,0);
	};
	~Tukey_asym(){delete cf;}
};

struct Leclerc_asym: public Asym_cos_function{
	Leclerc_asym(float sig_pos, float sig_neg, float alp)
		:Asym_cos_function(sig_pos,sig_neg,alp,NULL){
		cf = new Leclerc(0,0);
	};
	~Leclerc_asym(){delete cf;}
};
struct Leclerc_bis_asym: public Asym_cos_function{
	Leclerc_bis_asym(float sig_pos, float sig_neg, float alp)
		:Asym_cos_function(sig_pos,sig_neg,alp,NULL){
		cf = new Leclerc_bis(0,0);
	};
	~Leclerc_bis_asym(){delete cf;}
};
class Sin_fit
{
public:
	/* ====================  LIFECYCLE     ======================================= */
	Sin_fit ();
	Sin_fit ( 	CImg<float> & input_X,

							CImg<float> & input_Y
							);                             /* constructor */
	Sin_fit ( 	CImg<float> & input_Y);
	void assign( 	CImg<float> & stack, int option);
	void assign(  	CImg<float> & input_Y);
	void init ( 	CImg<float> & input_X,
								CImg<float> & input_Y,
								CImg<float> & input_beta
								);
	virtual void init (
						 CImg<float> & input_X,
						 CImg<float> & input_Y
						 );                             /* constructor */
	void test (  );
	void run ( CImg<float> & beta_output= CImg<float>::empty() );
	void w_fit ( CImg<float> & beta_output, CImg<float> & W );
	void irls ( CImg<float> & beta_output ,
							float err_sig ,
							float alpha,
							float (*pt2CostFunc)(float, float,float),
							CImg<float> *);
	void irls_no_init  ( CImg<float> & beta_output = CImg<float>::empty());
	inline CImg<float> get_residue() { fill_residue(); return residue;};
	inline CImg<float> get_J() { fill_J(); return J;};
	CImg<unsigned char> display_fit();
	virtual void set_starting_value ();
	inline void set_threshold(float val ){threshold=val;};
	inline void set_max_iter(int  val ){max_iter=val;};

	inline void set_scale_parameter (float val ){scale_parameter =val;cost_function->set_sig_noise(val);};

	void init_residue ( 	CImg<float> & aRESIDUE,
												CImg<float> & input_beta
												);

	// covariance matrix
	CImg<float> estimate_cov_matrix();
	float estimate_s_model ();
	CImg<float> cipra_robust_cov_matrix();
	CImg<float> ieng_robust_cov_matrix();


	// Getter & setters
	inline int get_nb_iter(){return nb_iter;};
	inline int get_max_iter(){return max_iter;};
	inline CImg<float> get_beta(){return beta;};
	inline void set_beta(CImg<float> & val ){ beta=val;};
	inline float get_beta(int i){return beta(i);};
	virtual inline float get_amp(){return beta(1);};
	virtual inline float get_phase(){return beta(2);};
	virtual inline float get_C(){return beta(0);};
	virtual inline void set_amp(float a){ beta(1)=a;};
	virtual inline void set_phase(float p){ beta(2)=p;};
	virtual inline void set_C(float c){ beta(0)=c;};
	inline void set_cost_function(Cost_function * val ){cost_function=val;};
	inline Cost_function * get_cost_function(){return cost_function;};
	inline void set_residue_weigth(CImg<float> & val ){residue_weigth=val.get_vector();};
	inline CImg<float> get_X(){return X;};
	inline CImg<float> get_Y(){return Y;};
	inline bool has_diverged(){return nb_iter>= max_iter; }
	virtual void fill_residue ();
	void generates_gnuplot_residual_plan(const char * file_name);



	Mult_evo_visualizer beta_evo;

	// Plot fit
	CFigure build_1d_graph(const char * marker="bx");
	void build_1d_graph(CFigure&,const char * marker="bx");
	void save_X_Y_gnuplot(const char * filename);

	void plot_fit(CFigure & fig,const char * );
	void plot_weighted_fit(CFigure & fig ,const char * color );
	CFigure build_1d_weighted_graph();
	void build_1d_weighted_graph(CFigure &fig);
	CFigure plot_residue();
	CFigure plot_residue(CFigure & fig);
	void print_state ( 	ostream & of);

	int measure_nb;



protected:
	float cos_func(float x );
	float cos_func_deriv ( int var_idx,float x);
	virtual float reg_funct(float x );
	virtual float reg_funct_deriv ( int var_idx,float x);
	virtual void fill_J ();
	void update_weights (  CImg<float> err_sig , float alpha, float (*pt2CostFunc)(float, float,float),
												 CImg<float> & W);
	void update_weights (  float err_sig , float alpha, float (*pt2CostFunc)(float, float,float),
												 CImg<float> & W);
	void recenter_beta ();
	virtual void fill_spring ();

	CImg<float>  Y;
	CImg<float > W;
	CImg<float> residue_weigth;

	CImg<float> X;
	Cost_function * cost_function;
	CImg<float>  beta;
	CImg<float> J;
	CImg<float>  residue;
	CImg<float> spring;
	CImg<float>  delta_beta;
	float threshold;
	float scale_parameter;
	float alpha_line;
	int max_iter;
	int nb_iter;
private:
	/* ====================  DATA MEMBERS  ======================================= */

	void cholesky_decomp_rec ( CImg<complex<float> > &A,
														 CImg<complex<float> > &id,
														 CImg<complex<float> > & A_copy,
														 CImg<complex<float> > & L,
														 int n,
														 CImg<complex<float> > & res);
	void cholesky_decomp ( CImg<float> & input,
												 CImg<float > & L,
												 CImg<float > & L_trans_conj);
};


cfl empirical_cov_matrix(cfl & betas,cfl & mean);

#endif /* end of include guard: SIN_FIT_H */

 
