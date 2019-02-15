#ifndef _SIN_FIT_FACT_H_
#define _SIN_FIT_FACT_H_

#include "sin_fit.hpp"


class Sin_fit_fact : public Sin_fit {
public:
		void set_starting_value();
protected:
	float reg_funct(float x){
		float Pi=3.141592;
		float w=2.*Pi/measure_nb;
		return get_C()*(1+get_amp()*sin(w*x+get_phase()));
	}
	float reg_funct_deriv(int var_idx, float x){
		float res=0;
		float pi=3.141592;
		float w=2.*pi/measure_nb;
	
		switch ( var_idx ) {
			case 0:	
				res = (1+get_amp()*sin(w*x+get_phase()));
				break;
				
			case 1:	
				res = get_C()*sin(w*x+get_phase());
				break;
				
				case 2:	
					res = get_C()*get_amp()*cos(w*x+get_phase());
				break;
			default:	
				printf("error : not enough variable\n");
				break;
	}				/* -----  end switch  ----- */
	return res;
	}
	
};

	void Sin_fit_fact::set_starting_value (  )
{
	float Pi=3.141592;
	float w=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
	set_C(Y.mean());
	set_amp((Y.max()-Y.min())/(2*get_C()));	/* FIXME */
	// beta(2)=2.*Pi/12.;                                   /* FIXME */
	set_phase(atan(2.5*w)+2.42342-atan(4.1*w));	/* FIXME */
}		/* -----  end of method Sin_fit::starting_value  ----- */

#endif /* _SIN_FIT_FACT_H_ */
