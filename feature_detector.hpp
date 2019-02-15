#ifndef _FEATURE_DETECTOR_H_
#define _FEATURE_DETECTOR_H_

#include "CImg.h"
#include "matusita.h"
#include <iostream>

#include <ctime>

using namespace cimg_library;
using namespace std;



void	build_kernel(	float sigma,
										CImg<float> & my_gaussian_kernel,
										CImg<float> & my_gaussian_kernel_deriv=CImg<float>::empty()){
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

#define DEFAULT_WINDOW_SIZE 7
#define DEFAULT_KERNEL_SIGMA 1.

struct light_feature{
	int x,y;
	float metric;
	light_feature(int ax,int ay, float ametric):x(ax),y(ay),metric(ametric){};
};

bool compare_light_feature(const light_feature& a, const light_feature& b )
{
  return a.metric > b.metric;
}

struct feature_calc{
	int window_size;
	float kernel_sigma;
	list<light_feature> f_list;
	cfl non_null_data;

	feature_calc(){
		window_size=DEFAULT_WINDOW_SIZE;
		kernel_sigma=DEFAULT_KERNEL_SIGMA;
	}

template<class T>
CImg<float>		measure_all_feature(const CImg<T> & img){
	return measure_all_feature_optim(img, CImg<char>(img,"xyz",1));
}

template<class T>
CImg<float>		measure_all_feature(const CImg<T> & img,const CImg<char> & neighbor_free){
	CImgList<float> G(3,img.width(),img.height());
	CImgList<float> grads(2);
	CImg<float> min_eigen_values(img,"xyz");
	min_eigen_values.fill(0);
	CImg<float> grad_gaussian_kernel,grad_gaussian_kernel_deriv;

	build_kernel(kernel_sigma,grad_gaussian_kernel,grad_gaussian_kernel_deriv);

	//Tensor computation//

	//Compute gradient using gaussian derivative
	grads.assign(img,img);
	grads[0].convolve(grad_gaussian_kernel_deriv);
	grads[0].convolve(grad_gaussian_kernel.get_transpose());
	grads[1].convolve(grad_gaussian_kernel);
	grads[1].convolve(grad_gaussian_kernel_deriv.get_transpose());

	//computing local tensor
	G[0]= grads[0].get_pow(2);
	G[1]= grads[0].get_mul(grads[1]);
	G[2]= grads[1].get_pow(2);


	// Tensor integration on feature size
	CImg<float> eigenlists;
	CImg<float> eigenvalues;
	CImg<float> tensor(2,2);
	cimg_forZ(img,z){
		cimg_for_insideXY(img,x,y,window_size){
			if(neighbor_free(x,y,z)){
				tensor.fill(0);
				for(int i=-window_size/2;i<=window_size/2 ;i++){
					float *stop;
					float *ptrs;
					stop= G[0].data(x+window_size/2,y+i);
					for (ptrs = (G[0].data(x-window_size/2,y+i)); ptrs<=stop;ptrs++){
						tensor(0,0)+=*ptrs;
					}
					stop= G[1].data(x+window_size/2,y+i);
					for (ptrs = (G[1].data(x-window_size/2,y+i)); ptrs<=stop; ptrs++){
						tensor(0,1)+=*ptrs;
					}
					stop= G[2].data(x+window_size/2,y+i);
					for (ptrs = (G[2].data(x-window_size/2,y+i)); ptrs<=stop;ptrs++ ){
						tensor(1,1)+=*ptrs;
					}
				}
				tensor(1,0)=tensor(0,1);

				// Store minimum eigen value;
				tensor.symmetric_eigen(eigenvalues,eigenlists);
				float eigenmin=eigenvalues.min();
				min_eigen_values(x,y,z)=eigenmin;

				f_list.push_back(light_feature(x,y,eigenmin));
			}
		}
	}
	// std::cout<< "elapsed time in measure_all_feature "
	// 				 <<  ( ( std::clock() - start ) / (double)CLOCKS_PER_SEC ) <<endl;


	return min_eigen_values;
}

template<class T>
CImg<float>		measure_all_feature_optim(const CImg<T> & img,const CImg<char> & neighbor_free){
	CImgList<float> G(3,img.width(),img.height());
	CImgList<float> grads(2);
	CImg<float> min_eigen_values(img,"xyz",0);
	non_null_data.assign(img.size());	int	non_null_data_idx=0;

	CImg<float> grad_gaussian_kernel,grad_gaussian_kernel_deriv;
	f_list.clear();

	build_kernel(kernel_sigma,grad_gaussian_kernel,grad_gaussian_kernel_deriv);

	//Tensor computation//

	//Compute gradient using gaussian derivative
	grads.assign(img,img);
	grads[0].convolve(grad_gaussian_kernel_deriv);
	grads[0].convolve(grad_gaussian_kernel.get_transpose());
	grads[1].convolve(grad_gaussian_kernel);
	grads[1].convolve(grad_gaussian_kernel_deriv.get_transpose());

	//computing local tensor
	G[0]= grads[0].get_pow(2);
	G[1]= grads[0].get_mul(grads[1]);
	G[2]= grads[1].get_pow(2);


	// Tensor integration on feature size
	float *ptrG0=G[0].data(),*ptrG1=G[1].data(),*ptrG2=G[2].data();
	const char * ptrR=neighbor_free.data();
	int y_jump=G[0].width()-window_size;
	int offset=-window_size/2*(1+G[0].width());

	cimg_forZ(img,z){
		ptrG0+=G[0].width()*window_size; ptrG1+=G[1].width()*window_size; ptrG2+=G[2].width()*window_size;
		ptrR+=neighbor_free.width()*window_size;
		cimg_for_insideY(img,y,window_size){
			ptrG0+=window_size; ptrG1+=window_size; ptrG2+=window_size;
			ptrR+=window_size;
			cimg_for_insideX(img,x,window_size){
					float tul=0;						// tensor upper left
					float tdr=0;						// tensor down right
					float tdl=0;						// tensor (symmmetric) down left
					float *start1,*start2,*start3;
					start1= ptrG0+offset; start2= ptrG1+offset; start3= ptrG2+offset;
					for(int i=-window_size/2;i<=window_size/2 ;i++){
						for (int j = 0; j< window_size; j++){
							tul+=*start1++; tdl+=*start2++; tdr+=*start3++;
						}
						start1+=y_jump;
						start2+=y_jump;
						start3+=y_jump;
					}

					float eigenmin = (tul+tdr - std::sqrt(pow(tul-tdr,2.f) + 4*pow(tdl,2.f)))/2;
					non_null_data(non_null_data_idx++)=eigenmin;


					if(*ptrR){
						min_eigen_values(x,y,z)=eigenmin;
						f_list.push_back(light_feature(x,y,eigenmin));
					}
				ptrR++;
				ptrG0++; ptrG1++; ptrG2++;
			}
			ptrG0+=window_size; ptrG1+=window_size; ptrG2+=window_size;
			ptrR+=window_size;
		}
		ptrG0+=G[0].width()*window_size; ptrG1+=G[1].width()*window_size; ptrG2+=G[2].width()*window_size;
		ptrR+=neighbor_free.width()*window_size;
	}

	non_null_data.assign(non_null_data.data(),non_null_data_idx);

	return min_eigen_values;
}
};

struct feature_detector{
	float feature_dist;
	int feature_size;
	int nb_feature_max;
	int nb_feature_to_detect;
	int nb_feature_detected;
	float th_scale_param;
	CImgList<float> metric_save;
	CImgList<float> neighbor_free_before_selection_save;
	CImgList<float> detected_spot_display;
	Evo_visualizer selected_metric;
	float min_threshold;

	feature_detector(){
	 feature_dist=10;
	 feature_size=7;
	 nb_feature_max=100;
	 th_scale_param=9;
	 nb_feature_to_detect=0; // nb feature to detect to keep the same number of non lost feature.
	 nb_feature_detected=0;
	 min_threshold=100000;
	}

	template<class T>
	CImg<char> build_neighboor_free(const CImg<T> & image,list<Feature<T> > & detected_feature_l)
	// Build an image showing locus with no detected spot
	// Also count the number of spot to be detected to reach the total expect number of spot.
	{
		CImg<char> neighbor_free;
		neighbor_free.assign(image.width(),image.height(),1,1,1);
		int neighborhood_occ=0;
		typename list<Feature<T> >::iterator it=detected_feature_l.begin();
		nb_feature_to_detect=nb_feature_max;
		cimg_for1(detected_feature_l.size(),i){
			if (!it->get_lost()){
				nb_feature_to_detect--;
				neighbor_free.draw_circle(it->get_x(),it->get_y(),feature_dist,&neighborhood_occ);
			}
			++it;
		}
		return neighbor_free;
	}


	template<class T>
	list<Feature<T> > feature_redetection(const CImg<T> & image,
																				list<Feature<T> > & detected_feature_l,
																				int frame_idx){
		CImg<char> neighbor_free;
		neighbor_free=build_neighboor_free(image,detected_feature_l);
		neighbor_free_before_selection_save.push_back(neighbor_free);

		feature_calc fc;
		fc.window_size=feature_size;
		cfl metric = fc.measure_all_feature_optim(image,neighbor_free);
		metric_save.push_back(metric);

		list<light_feature> lf_list=fc.f_list;
		list<Feature<T> > new_feature=feature_selection(image,lf_list,nb_feature_to_detect);

		typename list<Feature<T> >::iterator it=new_feature.begin();
		for(;it !=new_feature.end();++it) it->set_starting_frame_idx(frame_idx);

		return new_feature;
	}

	template<class T>
	list<Feature<T> > feature_redetection_by_threshold(const CImg<T> & image,
																										 list<Feature<T> > & detected_feature_l,
																										 int frame_idx,float threshold){

		CImg<char> neighbor_free;

		feature_calc fc;
		fc.window_size=feature_size;
		cfl metric=fc.measure_all_feature_optim(image,neighbor_free);
		// float m; float threshold= fabs(fc.non_null_data.min_max(m)-m)/3;
		// cout << "metric.median() : " << metric.median() << endl;
		// cout << "sqrt metric.variance(2) : " << sqrt(metric.variance(2)) << endl;
		// float	threshold= fc.non_null_data.median()+th_scale_param*sqrt(fc.non_null_data.variance(2));
		// // float m; float threshold= fabs(fc.non_null_data.min_max(m)-m)/3;
		// cout << "threshold MAD : " << threshold << endl;

		metric_save.push_back(metric);

		list<light_feature> lf_list=fc.f_list;
		list<Feature<T> > new_feature=feature_selection(image,lf_list,threshold);
		typename list<Feature<T> >::iterator it=new_feature.begin();
		for(;it !=new_feature.end();++it) it->set_starting_frame_idx(frame_idx);
		return new_feature;
	}

	template<class T>
	list<Feature<T> > feature_redetection_by_number_then_histogram(const CImg<T> & image,
																										 list<Feature<T> > & detected_feature_l,
																										 int frame_idx,int a_nb_feature_to_detect){

		// Kepp considering already tracked vesicle for threshold estimation
		CImg<char> neighbor_free;
		neighbor_free.assign(image.width(),image.height(),1,1,1);

		feature_calc fc;
		fc.window_size=feature_size;
		cfl metric=fc.measure_all_feature_optim(image,neighbor_free);
		metric_save.push_back(metric);

		// sort extracted metric (attached to their position)
		fc.f_list.sort(compare_light_feature);

		// Keep <a_nb_feature_to_detect> best features that are spaced by
		// <feature_dist> pixels.
		char neighborhood_occ=0;
		CImg<char> neighbor_free_selection;
		neighbor_free_selection.assign(image.width(),image.height(),1,1,1);

		list<light_feature>::iterator f=fc.f_list.begin();
		list<light_feature> N_best_features;
		cfl N_best_metric(a_nb_feature_to_detect);float * N_best_metric_ptr=N_best_metric.data();
		while (N_best_features.size()<a_nb_feature_to_detect){
			if (neighbor_free_selection(f->x,f->y)){
				// A perimeter to avoid features to be selected in this area
				neighbor_free_selection.draw_circle(f->x,f->y, feature_dist ,&neighborhood_occ);
				*N_best_metric_ptr++=f->metric;
				N_best_features.push_back(*f);
			}
			++f;
		}

		// Outliers metrics are detected using the matusita distance on
		// the metric histogram
		matusita_distance dist( N_best_metric,a_nb_feature_to_detect/2);
		float threshold = dist.minimize_distance_exhaustive();
		if (threshold < min_threshold) min_threshold = threshold;
		threshold=min_threshold;
		cout << "threshold : " << threshold << " " ;

		//	Display and print selected metrics histograms stats
		CFigure fig;
		float blue[3]={0,0,255};


		cfl X; X.assign(dist.histogram.size()+1).sequence(dist.value_min,dist.value_max);
		cfl Y(1); Y= dist.histogram.get_append(Y,'x');
		fig.bar(X,Y,blue);					//
		fig.erase().set_axis().replot();

		CImgList<> tau_dist=dist.tau_dist.get_data();
		tau_dist[1].normalize(fig.get_ymin(),fig.get_ymax());
		fig.plot(tau_dist[0], tau_dist[1],"r-");
		char filename[40]; cimg::number_filename("matusita_dist_minimizing.png",frame_idx,2,filename);
		fig.save(filename);

		// Considering threshold and the already detectd spot we pick the
		// best spot out the <a_nb_feature_to_detect> best features
		CImg<char> neighbor_free_tracked_spot;
		neighbor_free_tracked_spot=build_neighboor_free(image,detected_feature_l);
		for(f=N_best_features.begin();f!=N_best_features.end();){
			if((f->metric <threshold)||(!neighbor_free_tracked_spot(f->x,f->y))) {
				f=N_best_features.erase(f); }
			else
				{f++;}
		}

		list<Feature<T> > new_feature;
		nb_feature_detected=0;
		for(f=N_best_features.begin();f!=N_best_features.end();f++){
				Feature<T> feature(f->x,f->y,feature_size,false);
				feature.set_starting_frame_idx(frame_idx);
				new_feature.push_back(feature);
				nb_feature_detected++;
		}

		return new_feature;
	}


	template<class T>
	list<Feature<T> > feature_detection(const CImg<T> & image,const float sig=1.){
			feature_calc fc;
			fc.window_size=feature_size;
			fc.kernel_sigma=sig;
			fc.measure_all_feature(image);
			list<light_feature> lf_list=fc.f_list;
			return feature_selection(image,lf_list,nb_feature_max);
	}

	template<class T>
	list<Feature<T> > feature_selection(const CImg<T> & image,list<light_feature> & lf_list,int a_nb_feature_to_detect,float threshold=0){

		// feature selection based on eigen values
		lf_list.sort(compare_light_feature);

		list<Feature<T> > feature_l;
		int nb_feature=0;
		char neighborhood_occ=0;
		CImg<char> neighbor_free;
		neighbor_free.assign(image.width(),image.height(),1,1,1);
		CImg<char> neighbor_free_display;
		neighbor_free_display.assign(image.width(),image.height(),1,1,1);
		selected_metric.assign(lf_list.size());

		list<light_feature>::iterator f=lf_list.begin();
		nb_feature_detected=0;
		while (nb_feature<a_nb_feature_to_detect){
			if (neighbor_free(f->x,f->y)){
				//Setting a perimeter to avoided other feature to be selected in this
				//area
				neighbor_free.draw_circle(f->x,f->y, feature_dist ,&neighborhood_occ);
				selected_metric.add_data(f->metric);
				// neighbor_free_display.draw_circle(f->x,f->y, 3 ,&neighborhood_occ);
				nb_feature++;
				Feature<T> feature(f->x,f->y,feature_size,false);
				feature_l.push_back(feature);
				nb_feature_detected++;
			}
			++f;
		} // end for each selected feature ...
		// detected_spot_display.push_back(neighbor_free_display);
		return feature_l;
	}

	template<class T>
	list<Feature<T> > feature_selection(const CImg<T> & image,list<light_feature> & lf_list,float threshold){
		cout << "feature selection by threshold" << endl;

		// feature selection based on eigen values
		lf_list.sort(compare_light_feature);

		// feature selection based on eigen values
		list<Feature<T> > feature_l;
		char neighborhood_occ=0;
		CImg<char> neighbor_free;
		neighbor_free.assign(image.width(),image.height(),1,1,1);

		CImg<char> neighbor_free_display;
		neighbor_free_display.assign(image.width(),image.height(),1,1,1);

		nb_feature_detected=0;
		list<light_feature>::iterator f=lf_list.begin();
		while (f->metric > threshold){
			if (neighbor_free(f->x,f->y)){
				// Avoiding selecting high metric is the neighborhood
				neighbor_free.draw_circle(f->x,f->y, feature_dist ,&neighborhood_occ);
				neighbor_free_display.draw_circle(f->x,f->y, 3 ,&neighborhood_occ);
				Feature<T> feature(f->x,f->y,feature_size,false);
				feature_l.push_back(feature);
				nb_feature_detected++;
			}
			++f;
		} // end for each selected feature ...
		detected_spot_display.push_back(neighbor_free_display);
		return feature_l;
	}
};

#endif /* _FEATURE_DETECTOR_H_ */

