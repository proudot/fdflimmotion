#ifndef _BG_ESTIMATOR_H_
#define _BG_ESTIMATOR_H_

#include "CImg.h"

template<class T>
class bg_estimator
{
public:
	inline bg_estimator(){};
	virtual float compute(CImg<T> img,int x0,int y0,int x1,int y1)=0;
	float bg_estimate;
};

template<class T>
class bg_min : public bg_estimator<T>
{
public:
	inline float compute(CImg<T> img,int x0,int y0,int x1,int y1){
		this->bg_estimate=img(x0,y0);
		for (int i = x0; i <= x1 ; ++i)
			for (int j = y0; j <= y1 ; ++j)
				{
					if( this->bg_estimate > img(i,j) ) this->bg_estimate= img(i,j);
				}
		return this->bg_estimate;
	}
};

template<class T>
class bg_LTS : public bg_estimator<T>
{
public:
	inline float compute(CImg<T> img,int x0,int y0,int x1,int y1){
		CImg<T> values=img.crop(x0,y0,x1,y1);
		values.unroll('x');
		values.sort();
		values=values.get_shared_points(0,values.size()/2);
		this->bg_estimate=mean_type(values);
		return this->bg_estimate ;
	}
protected:
	inline float mean_type(CImg<T> v){return v.median();};
};

template<class T>
class bg_LTS_mean: public bg_LTS<T>
{
protected:
	inline float mean_type(CImg<T> v){return v.mean();};
};


template<class T>
class bg_surround : public bg_estimator<T>
{
public:
	inline float compute(CImg<T> img,int x0,int y0,int x1,int y1){
		float bg_value;
		for (int i = x0; i <= x1; ++i) {		
			bg_value+=img(i,y0);
			bg_value+=img(i,y1);
		}

		for (int i = y0; i <= y1; ++i)
			{
				bg_value+=img(x0,i);
				bg_value+=img(x1,i);
			}
		bg_value/= 2*(x1-x0+1+(y1-y0+1)) - 4;
		this->bg_estimate = bg_value;
		return bg_value;
	}
};


#endif /* _BG_ESTIMATOR_H_ */
