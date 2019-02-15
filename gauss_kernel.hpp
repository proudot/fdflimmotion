#ifndef _GAUSS_KERNEL_H_
#define _GAUSS_KERNEL_H_

#include "CImg.h"
using namespace cimg_library;

class Gauss_kernel
{
public:
	static CImg<float>  k;
	static CImg<float> k_deriv;
	static inline   CImg<float> get_k (){ 
		CImg<float> k1;
		CImg<float> k2;
		build_kernel(1.,k1,k2);
		return k1;
	}
	static inline   CImg<float> get_k_deriv (){ 
		CImg<float> k1;
		CImg<float> k2;
		build_kernel(1.,k1,k2);
		return k2;
	}
protected:
	static void build_kernel(	float sigma,
														CImg<float> & my_gaussian_kernel,
														CImg<float> & my_gaussian_kernel_deriv);
	
};

CImg<float> Gauss_kernel::k=Gauss_kernel::get_k();
CImg<float> Gauss_kernel::k_deriv=Gauss_kernel::get_k_deriv();

void	Gauss_kernel::build_kernel(	float sigma,
																			CImg<float> & my_gaussian_kernel,
																			CImg<float> & my_gaussian_kernel_deriv){
	float black[1]={1};
	//approx of ln(10)
	float sqrtln10=1.5174;
	float sqrt2=1.4142;

	//compute the size to avoid value below 0.01
	int gaussian_kernel_size= 2*(int)(2*sigma*sqrtln10)+1; 
	int gaussian_kernel_deriv_size= 2*(int)(2*sqrt2*sigma*sqrtln10)+1; 
	my_gaussian_kernel.assign(gaussian_kernel_size);
	my_gaussian_kernel.draw_gaussian(gaussian_kernel_size/2,sigma,black);
	my_gaussian_kernel_deriv.assign(gaussian_kernel_deriv_size);
	my_gaussian_kernel_deriv.draw_gaussian(gaussian_kernel_deriv_size/2,sigma,black);

	float den=0.;
	cimg_forX(my_gaussian_kernel_deriv,x){
		my_gaussian_kernel_deriv(x)=(gaussian_kernel_deriv_size/2-x)*my_gaussian_kernel_deriv(x)/(sigma*sigma);
		den-=(x-gaussian_kernel_deriv_size/2 )*my_gaussian_kernel_deriv(x)/(sigma*sigma); // density to avoid summing the kernel 
	}
	my_gaussian_kernel_deriv/=den;
	my_gaussian_kernel/=my_gaussian_kernel.sum();
}

#endif /* _GAUSS_KERNEL_H_ */
