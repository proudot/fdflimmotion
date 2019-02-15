#ifndef SIMU_SPOT_H
#define SIMU_SPOT_H

#include "CImg.h"

#include <iostream>
#include "math.h"
#include	<list>
#include <assert.h>
#include <vector>
#include "spot.hpp"

using namespace cimg_library;
using namespace std;

struct Param_simu {

	Param_simu(){

	}
};


template<class T>
class Simu_spot : public Spot<T>  {

private:
  


	void draw_patch_gaussian( const T * color);
	bool dynamic;
	float gaussian_sigma;
	int spot_patch_size;
	int move_square;
public:
	CImg<T> value;
	inline void set_dynamic(bool val ){dynamic=val;};
	inline void set_gaussian_sigma(float val ){gaussian_sigma=val;};
	inline void set_spot_patch_size(int val ){spot_patch_size=val;};
	
	Simu_spot();
	Simu_spot(int,int,int,bool dyn= false,float agaussian_sigma=1.5);
	void draw_patch(const  T* ,int type=0);
	void apply_patch(CImg<T>& img);
	void build_kernel(int,CImg<float> &);
	void scale_patch(const float scale);
	void move(CImg<T> &);
	Simu_spot<T>& brownian_new_move(CImg<T> & img,float sigma,const CImg<char> & available_mask,
																	const CImg<float> & scale,const CImg<float> & translation_matrix);
	Simu_spot<T>& brownian_new_move(CImg<T> & img,float sigma,const CImg<char> & available_mask);
	Simu_spot<T>& gaussian_new_int(float sigma);
	void rand_z_move();
		
};
// END OF HEADER

template<class T>
Simu_spot<T>::Simu_spot():Spot<T>(){
		gaussian_sigma				=				 1.5 ;
		spot_patch_size			  =				 51;
		move_square=2;
		dynamic=true;
};

template<class T>
Simu_spot<T>::Simu_spot(int x, int y, int aspot_patch_size,bool bDynamic,float agaussian_sigma) :
	Spot<T>(x,y,spot_patch_size)
{ 
	gaussian_sigma = agaussian_sigma ;
		move_square=2;
		spot_patch_size = aspot_patch_size;
		dynamic=bDynamic;
}

template<class T>
void Simu_spot<T>::rand_z_move(){
	float move = (float)(rand()%100+90)/100.;
	value*=move;
}
template<class T>
void Simu_spot<T>::scale_patch(const float scale){
	value*=scale;
}

template<class T>
void Simu_spot<T>::draw_patch(const T * color, int type){
	switch(type){
	case 0:
		draw_patch_gaussian(color);
		break;
	default:
		cout << " unknown type " << endl;
		break;
	}
}

template< class T >
void Simu_spot<T>::draw_patch_gaussian(const T * color){
	value.assign(spot_patch_size,spot_patch_size);
	// float sigma_err_x=(float)((rand() % 20 )-10)/10;
	// float sigma_err_y=(float)((rand() % 20 )-10)/10;
	// float sigma_err_xy=(float)((rand() % 20 )-10)/10;
	float sigma_err_x=0;
	float sigma_err_y=0;
	float sigma_err_xy=0;
	float scale=0.4;
	float dx = this->center(0)-floor(this->center(0));
	float dy = this->center(1)-floor(this->center(1));
	value.draw_gaussian(spot_patch_size/2+dx,spot_patch_size/2+dy,
											CImg<float>(2,2,1,1,
																	gaussian_sigma + scale*sigma_err_x,
																	scale*sigma_err_xy,
																	scale*sigma_err_xy,
																	gaussian_sigma + scale*sigma_err_y),
												color);

	Spot<T>::width=value.width();
	Spot<T>::height=value.height();
}

template< class T >
void Simu_spot<T>::apply_patch(CImg<T>& img){
	// DIRTY 
	// DIRTY 
	// DIRTY 
	int img_x,img_y;
	// T mean = 0;
	// cimg_forXY(value,x,y){
	// 	img_x=Spot<T>::center(0)+x-Spot<T>::width/2;
	// 	img_y=Spot<T>::center(1)+y-Spot<T>::height/2;
	// 	mean += img(img_x,img_y);
	// }
	// mean/=value.width()*value.height();
	// draw_patch_gaussian(&color);

	cimg_forXY(value,x,y){
		img_x=floor(Spot<T>::center(0))+x-Spot<T>::width/2;
		img_y=floor(Spot<T>::center(1))+y-Spot<T>::height/2;
		img(img_x,img_y)+=value(x,y);
	}
	// img_x=Spot<T>::center(0)-Spot<T>::width/2;
	// img_y=Spot<T>::center(1)-Spot<T>::height/2;
	// img.draw_image(img_x,img_y,0,0,value,value.get_normalize(0,1));
}


template<class T>
void Simu_spot<T>::move(CImg<T> & img){
	CImg<float> c(this->center);
	if(dynamic){
		do{
			this->center=c;
			this->center+=Spot<T>::center.get_rand(-move_square,move_square);
		}while(this->is_out(img));
	};
	set_new_move(Spot<T>::center);
}

template<class T>
Simu_spot<T>& Simu_spot<T>::brownian_new_move(CImg<T> & img,float sigma, const CImg<char> & available_mask){
	if(dynamic){
		float x,y;
		do{
			 x =(this->get_x()+ cimg::grand()*sigma);
			 y =(this->get_y()+ cimg::grand()*sigma);
			
		}while(!available_mask(x,y));
		this->set_x(x);this->set_y(y);
	};
	this->set_new_move(this->center);
	return *this;
}

template<class T>
Simu_spot<T>& Simu_spot<T>::brownian_new_move(CImg<T> & img,float sigma,
																							const CImg<char> & available_mask,
																							const CImg<float> & scale_matrix,
																							const CImg<float> & translation_matrix)
{

	//	CImg<unsigned char> img_disp(available_mask);
	//  img_disp.append(available_mask,'c').append(available_mask,'c');
	const	unsigned char red[3]={5,0,0};

	this->set_center(scale_matrix*this->get_center()+translation_matrix);
	//	img_disp.draw_point(this->get_x(),this->get_y(),red).display();
	int count = 0 ;
 	if(dynamic){
		bool available_pixel=false;
		do{
 			float x =(this->get_x()+ cimg::grand()*sigma);
			float y =(this->get_y()+ cimg::grand()*sigma);
			if(!this->is_out(available_mask)) available_pixel=available_mask(x,y);
			this->set_x(x);this->set_y(y);
			count++;
			if(count%5==0) { sigma*=2;}
		}while((!available_pixel)&&(count<20));	
	}
	if(count<20) 	this->set_new_move(this->center);
	return *this;
}

template<class T>
Simu_spot<T>& Simu_spot<T>::gaussian_new_int(float sigma){
	this->intensities.append(CImg<float>(1,1,1,1,cimg::grand()*sigma+this->get_intensity(0)),'x');
	return *this;	
}

#endif /* SIMU_SPOT_H */

