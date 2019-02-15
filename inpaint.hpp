//
// inpaint.hpp
//
// Created by Philippe Roudot 2013-06-27
//
// Copyright (C) INRIA 2013 all right reserved
// 



#ifndef _INPAINT_H_
#define _INPAINT_H_
#include "CImg.h"
#include <iostream>
using namespace cimg_library;
using namespace std;

//!	inpaint a 2D image on mask using basic anisotropic linear interpolation
void inpaint_interp(CImg<> &img, CImg<unsigned char> & mask){
	unsigned char* mptr = mask.data();
	CImg<> valw(img),valn(img),vale(img),vals(img);
	int x_start=0,y_start=0;
	cimg_forXY(mask,x,y){
		if((*mptr++)&&(x>0)&&(y>0)){
			valw(x,y,0)=valw(x-1,y,0); 
			valn(x,y)=valn(x,y-1); 
		}
	};
	mptr=mask.data(mask.size()-1);
	for(int y=mask.height()-1; y>-1; y--){
		for(int x=mask.width()-1; x>-1; x--){
			if(*mptr--){
				vale(x,y)=vale(x+1,y);
				vals(x,y)=vals(x,y+1);
			}
		}
	};
	img=(((vale+valw+valn+vals)/4)-img).get_mul(mask.get_threshold(0,false,true))+img;
};

CImg<> get_inpaint_interp(CImg<> & img,CImg<unsigned char> & mask){
	CImg<> res(img); inpaint_interp(res,mask); return res;
}

//!	inpaint a 2D image on mask using basic anisotropic linear interpolation.
//  Take space into account
void inpaint_interp_space(CImg<> &img, CImg<unsigned char> & mask){
	unsigned char* mptr = mask.data();
	CImg<> valw(img),valn(img),vale(img),vals(img);
	CImg<> distw(img,"xyzc",0),distn(img,"xyzc",0),diste(img,"xyzc",0),dists(img,"xyzc",0);
	int x_start=0,y_start=0;
	cimg_forXY(mask,x,y){
		if((*mptr++)&&(x>0)&&(y>0)){
			valw(x,y,0)=valw(x-1,y,0); distw(x,y)=distw(x-1,y)+1;
			valn(x,y)=valn(x,y-1); distn(x,y)=distn(x,y-1)+1;
		}
	};
	mptr=mask.data(mask.size()-1);
	for(int y=mask.height()-1; y>-1; y--){
		for(int x=mask.width()-1; x>-1; x--){
			if((*mptr--)&&(y<(mask.height()-1))&&(x<(mask.width()-1))){
				vale(x,y)=vale(x+1,y); diste(x,y)=diste(x+1,y)+1;
				vals(x,y)=vals(x,y+1); dists(x,y)=dists(x,y+1)+1;
			}
		}
	};

	CImg<> total_dist=diste+distw+dists+distn;
	img=( ( (vale.get_mul(total_dist-diste)+valw.get_mul(total_dist-distw)
					 + valn.get_mul(total_dist-distn)+vals.get_mul(total_dist-dists)
					 ).get_div(3*total_dist+(1-mask)) 
					)-img).get_mul(mask.get_threshold(0,false,true))+img;
};


CImg<> get_inpaint_interp_space(CImg<> & img,CImg<unsigned char> & mask){
	CImg<> res(img); inpaint_interp_space(res,mask); return res;
}


void test_inpaint(){
	CImg<unsigned char> test_image(400,400,1,1,0);
	((test_image.noise(30,2)*=100).blur(5)+=10);
	
	CImg<unsigned char> mask(test_image,"xyz",0);
	const unsigned char one=1;
	mask.draw_circle(200,200,20,&one);

	CImg<> masked_image=test_image.get_mul(1-mask);
	(test_image,
	 masked_image,
	 get_inpaint_interp(masked_image,mask),
	 get_inpaint_interp_space(masked_image,mask)
	 ).display();
}

#endif /* _INPAINT_H_ */
